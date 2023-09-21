#include "iftRobot.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Stream.h"


iftRobot *iftCreateRobot(const char *orig_path, const char *gt_path, int stop_iter, float stop_acur)
{
    iftRobot *bot = iftAlloc(1, sizeof (*bot));

    bot->orig = iftReadImageByExt(orig_path);
    bot->error = iftCreateImage(bot->orig->xsize, bot->orig->ysize, bot->orig->zsize);
    bot->segm = iftCreateImage(bot->orig->xsize, bot->orig->ysize, bot->orig->zsize);

    bot->gt_img = iftReadImageByExt(gt_path);
    int mid = iftMaximumValue(bot->gt_img) / 2;
    for (int p = 0; p < bot->gt_img->n; p++) {
        bot->gt_img->val[p] = (bot->gt_img->val[p] > mid);
    }

    bot->mrk_radius = iftCircular(2.0f);
    bot->limit_mrk_rad = iftCircular(3.0f);

    bot->stop_iter = stop_iter;
    bot->cur_iter = 0;

    bot->stop_acur = stop_acur;
    bot->cur_acur = 0.0f;

    bot->converged = false;
    bot->max_seed_per_iter = 1;

    return bot;
}


void iftResetRobotSegmentation(iftRobot *bot) {
    iftSetImage(bot->error, 0);
    iftSetImage(bot->segm, 0);
    bot->cur_iter = 0;
    bot->cur_acur = 0.0f;
    bot->converged = false;
    iftDestroyLabeledSet(&bot->seeds);
}


void iftDestroyRobot(iftRobot **bot)
{
    iftRobot *aux = *bot;
    if (aux != NULL) {
        iftDestroyImage(&aux->orig);
        iftDestroyImage(&aux->error);
        iftDestroyImage(&aux->gt_img);
        iftDestroyImage(&aux->segm);
        iftDestroyLabeledSet(&aux->seeds);
        iftDestroyAdjRel(&aux->mrk_radius);
        iftDestroyAdjRel(&aux->limit_mrk_rad);
    }
    aux = NULL;
}


void iftRobotUpdateError(iftRobot *bot)
{
    for (int p = 0; p < bot->segm->n; p++) {
        bot->error->val[p] = 0;
        if (bot->segm->val[p] != bot->gt_img->val[p]) {
            bot->error->val[p] = bot->gt_img->val[p] + 1; // increases for that bg errors have label 1}
        }
    }

    bot->cur_iter++;
    bot->cur_acur = iftDiceSimilarity(bot->segm, bot->gt_img);
}


void iftRobotPrintInfo(iftRobot *bot)
{
    printf("Iteration: %d\n", bot->cur_iter);
    printf("Accuracy: %f\n", bot->cur_acur);
}


bool iftRobotContinue(iftRobot *bot)
{
    bool out = true;
    if (bot->cur_acur > bot->stop_acur) {
        printf("LIDS Bot has passed the set accuracy:\nCriteria: %f\nCurrent: %f\n", bot->stop_acur, bot->cur_acur);
        out = false;
    }
    if (bot->cur_iter > bot->stop_iter) {
        printf("LIDS Bot has finished the desired %d iterations\n", bot->stop_iter);
        out = false;
    }
    if (bot->converged) {
        printf("LIDS Bot is not able to insert anymore markers\n");
        out = false;
    }
    return out;
}


void iftRobotFindSeedsMSSkel(iftRobot *bot) {

    iftImage *relabel_img = NULL, *dist_img;
    iftAdjRel *A8nb = iftCircular(sqrtf(2.0f));

    iftFImage *msskel = iftMSSkel2D(bot->error, A8nb, IFT_INTERIOR, &dist_img, &relabel_img);

    iftIntArray *labels = iftGetObjectLabels(relabel_img);
    iftIntArray *max_vals = iftCreateIntArray(labels->n + 1); // include bg
    iftIntArray *max_pos = iftCreateIntArray(labels->n + 1); // include bg
    iftIntArray *max_area = iftCreateIntArray(labels->n + 1);
    iftIntArray *max_dist = iftCreateIntArray(labels->n + 1);
    max_pos->val[0] = -1; // bg could not have skeleton, then there is no pixel for it

    /* New priority for seeds selection */
    iftDestroyIntArray(&labels);
    iftImage *area_img = iftCreateImageFromImage(bot->error);

    for (int p = 0; p < msskel->n; p++) {
        // bg pixels has object relabeled values after edt computing, then
        // we need to check its original value to guarantee the real label of bg pixels
        int label = (bot->gt_img->val[p] == 0) ? 0 : relabel_img->val[p];

        if(area_img->val[p] == 0) {
            iftObjectAreaFromPixel(relabel_img, p, area_img);
        }

        if ((area_img->val[p] >= max_area->val[label]) && (msskel->val[p] >= max_vals->val[label]) && (dist_img->val[p] > max_dist->val[label]) ) {
            max_vals->val[label] = msskel->val[p];
            max_pos->val[label]  = p;
            max_area->val[label] = area_img->val[p];
            max_dist->val[label] = dist_img->val[p];
        }
    }
    iftDestroyImage(&area_img);

    iftQuickSort(max_area->val, max_pos->val, 0, max_pos->n - 1, IFT_DECREASING);

    iftDestroyImage(&relabel_img);
    iftDestroyImage(&dist_img);
    iftDestroyIntArray(&max_vals);
    iftDestroyIntArray(&max_area);
    iftDestroyIntArray(&max_dist);

    iftLabeledSet *out = bot->seeds;
    int count_markers = 0;
    int count_pixels = 0;
    for (int o = 0; o < max_pos->n && count_markers < bot->max_seed_per_iter; o++) {
        int q = max_pos->val[o];
        if (q != -1 && bot->error->val[q] != 0) {
            iftVoxel u = iftGetVoxelCoord(bot->error, q);

            /* shift the seed if it's touching a border */
            for (int i = bot->limit_mrk_rad->n - 1; i > 0; i--){
                iftVoxel v = iftGetAdjacentVoxel(bot->limit_mrk_rad, u, i);
                if(iftValidVoxel(bot->error, v) && (u.x != v.x && u.y != v.y && u.z != v.z)){
                    u.x -= bot->limit_mrk_rad->dx[i];
                    u.y -= bot->limit_mrk_rad->dy[i];
                    u.z -= bot->limit_mrk_rad->dz[i];
                }
            }

            int in_border = 0;

            /* check if the seed is touching a border */
            for (int i = 0; i < bot->limit_mrk_rad->n; i++){
                iftVoxel v = iftGetAdjacentVoxel(bot->limit_mrk_rad, u, i);
                if(iftValidVoxel(bot->gt_img, v) && (bot->gt_img->val[q] != iftImgVoxelVal(bot->gt_img, v))){
                    in_border = 1;
                }
            }

            if (in_border == 0) {
                count_markers++;
                for (int i = 0; i < bot->mrk_radius->n; i++) {
                    iftVoxel v = iftGetAdjacentVoxel(bot->mrk_radius, u, i);

                    if (iftValidVoxel(bot->gt_img, v)){
                        /* It was adding the same seed multiple times */
                        int p = iftGetVoxelIndex(bot->gt_img, v);
                        if (!iftLabeledSetHasElement(out, p)) {
                            iftInsertLabeledSet(&out, p, bot->gt_img->val[q]);
                            count_pixels++;
                        }
                    }
                }
            }
        }
    }
    bot->seeds = out;

    if (count_pixels == 0)
        bot->converged = true;

    iftDestroyAdjRel(&A8nb);
    iftDestroyFImage(&msskel);
    iftDestroyIntArray(&max_pos);
}


void iftRobotFindSeedsCenterOfMass(iftRobot *bot)
{
    iftIntArray *labels = iftGetObjectLabels(bot->error);
    iftIntArray *max_pos = iftCreateIntArray(labels->n + 1);
    iftIntArray *max_area = iftCreateIntArray(labels->n + 1);// include bg
    max_pos->val[0] = -1; // bg could not have skeleton, then there is no pixel for it

    /* New priority for seeds selection */
    iftImage *area_img = iftCreateImageFromImage(bot->error);

    for (int p = 0; p < bot->error->n; p++) {
        // bg pixels has object relabeled values after edt computing, then
        // we need to check its original value to guarantee the real label of bg pixels
        int label = bot->error->val[p];

        if (max_pos->val[label] == 0) {
            if (area_img->val[p] == 0) {
                iftObjectAreaFromPixel(bot->error, p, area_img);
            }

            if (area_img->val[p] > max_area->val[label]) {

                double x = 0.0f;
                double y = 0.0f;

                int count = 0;
                for (int q = 0; q < bot->error->n; q++){
                    if (area_img->val[p] == area_img->val[q]) {
                        x += iftGetXCoord(area_img, q);
                        y += iftGetYCoord(area_img, q);
                        count++;
                    }
                }

                iftVoxel v;
                v.x = x / count;
                v.y = y / count;
                v.z = 0;
                max_pos->val[label] = iftGetVoxelIndex(area_img, v);
            }
        }
    }

    iftQuickSort(max_area->val, max_pos->val, 0, max_pos->n - 1, IFT_DECREASING);

    iftDestroyIntArray(&labels);
    iftDestroyIntArray(&max_area);
    iftDestroyImage(&area_img);

    iftLabeledSet *out = bot->seeds;
    int count_markers = 0;
    int count_pixels = 0;
    for (int o = 0; o < max_pos->n && count_markers < bot->max_seed_per_iter; o++) {
        int q = max_pos->val[o];
        if (q != -1) {
            iftVoxel u = iftGetVoxelCoord(bot->gt_img, q);

//            /* shift the seed if is touching a border */
            for (int i = bot->limit_mrk_rad->n-1; i > 0; i--) {
                iftVoxel v = iftGetAdjacentVoxel(bot->limit_mrk_rad, u, i);
                if (iftValidVoxel(bot->gt_img, v) && (u.x != v.x && u.y != v.y && u.z != v.z)) {
                    u.x -= bot->limit_mrk_rad->dx[i];
                    u.y -= bot->limit_mrk_rad->dy[i];
                    u.z -= bot->limit_mrk_rad->dz[i];
                }
            }

            int in_border = 0;

            /* check if the seed it's touching a border */
            for (int i = 0; i < bot->limit_mrk_rad->n; i++){
                iftVoxel v = iftGetAdjacentVoxel(bot->limit_mrk_rad, u, i);
                if(iftValidVoxel(bot->gt_img, v) && (bot->gt_img->val[q] != iftImgVoxelVal(bot->gt_img, v))){
                    in_border = 1;
                }
            }

            if (in_border == 0) {
                count_markers++;
                for (int i = 0; i < bot->mrk_radius->n; i++) {
                    iftVoxel v = iftGetAdjacentVoxel(bot->mrk_radius, u, i);

                    if (iftValidVoxel(bot->gt_img, v)) {
                        /* It was adding the same seed multiple times */
                        int p =iftGetVoxelIndex(bot->gt_img, v);
                        if (!iftLabeledSetHasElement(out, p)) {
                            iftInsertLabeledSet(&out, p, bot->gt_img->val[q]);
                            count_pixels++;
                        }
                    }
                }
            }
        }
    }
    if (count_pixels == 0)
        bot->converged = true;

    bot->seeds = out;

    iftDestroyIntArray(&max_pos);
}
