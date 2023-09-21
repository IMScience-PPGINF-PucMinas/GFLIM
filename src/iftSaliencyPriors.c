#include "iftSaliencyPriors.h"

iftSaliencyArray *iftFocusPrior(iftImage *superpixel_img, iftMImage *mimg, float variance){
    iftIntArray *region_sizes;
    iftDblArray *region_sizes_prob;
    region_sizes = iftCountLabelSpels(superpixel_img);
    region_sizes_prob = iftCreateDblArray(region_sizes->n-1);
    /*------------Find edges ------------*/
    iftImage *border_image = iftBorderProbImage(mimg);
    int threshold = iftOtsu(border_image);

    /*-------------- Compute Focusness ------------------*/
    iftSaliencyArray *focus = iftCreateSaliencyArray(region_sizes_prob->n);
    iftIntArray *border_counter = iftCreateIntArray(region_sizes_prob->n);

    iftAdjRel *A = iftRectangular(5, 5);

    for (int p=0; p < superpixel_img->n; p++) {
        int superpixel = superpixel_img->val[p] - 1;
        iftVoxel u = iftGetVoxelCoord(superpixel_img,p);
        for (int i=0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(superpixel_img,v)){
                int q = iftGetVoxelIndex(superpixel_img,v);
                if (superpixel_img->val[p] != superpixel_img->val[q]){
                    focus->val[superpixel] -= (float)(border_image->val[q] >= threshold);
                    border_counter->val[superpixel]+=1;
                }
            } else {
                border_counter->val[superpixel]+=1;
                break;
            }
        }
    }

    for(int s = 0; s < focus->n; s++)
        focus->val[s]= expf(-focus->val[s]/((float)border_counter->val[s]*variance));

    float max=0, min=IFT_INFINITY_FLT;
    for(int s=0; s < focus->n; s++){
        if(focus->val[s] > max)
            max = focus->val[s];
        if(focus->val[s] < min)
            min = focus->val[s];
    }

    for(int s = 0; s < focus->n; s++)
        focus->val[s] = 1 * (focus->val[s] - min) / (max - min);

    iftDestroyIntArray(&region_sizes);
    iftDestroyIntArray(&border_counter);
    iftDestroyDblArray(&region_sizes_prob);
    iftDestroyImage(&border_image);
    iftDestroyAdjRel(&A);

    focus->superpixel_image = superpixel_img;

    return focus;
}

iftSaliencyArray *iftFocusXRayPrior(iftImage *superpixel_img, iftMImage *mimg, float variance, float border_threshold){
    iftIntArray *region_sizes;
    iftDblArray *region_sizes_prob;
    region_sizes = iftCountLabelSpels(superpixel_img);
    region_sizes_prob = iftCreateDblArray(region_sizes->n-1);
    /*---------- Find edges -----------*/
    iftImage *border_image = iftBorderProbImage(mimg);
    iftImage *bin_border_image = iftThreshold(border_image, iftRound(iftOtsu(border_image)*border_threshold), iftMaximumValue(border_image), 255);

    /*----------- Compute Focusness --------*/
    iftSaliencyArray *focus = iftCreateSaliencyArray(region_sizes_prob->n);
    iftIntArray *border_counter = iftCreateIntArray(region_sizes_prob->n);

    iftAdjRel *A = iftRectangular(5, 5);

    for (int p=0; p < superpixel_img->n; p++) {
        int s = superpixel_img->val[p] - 1;
        focus->val[s] += (float)bin_border_image->val[p];
    }

    float max_focus = 0;
    for(int s = 0; s < focus->n; s++)
        if (focus->val[s] > max_focus)
            max_focus = focus->val[s];

    for(int s = 0; s < focus->n; s++)
        focus->val[s]= 1 - expf(-(float)focus->val[s]/(max_focus * variance));

    float max=0, min=IFT_INFINITY_FLT;
    for(int s=0; s < focus->n; s++){
        if(focus->val[s] > max)
            max = focus->val[s];
        if(focus->val[s] < min)
            min = focus->val[s];
    }

    for(int s = 0; s < focus->n; s++)
        focus->val[s] = 1 * (focus->val[s] - min) / (max - min);

    iftDestroyIntArray(&region_sizes);
    iftDestroyIntArray(&border_counter);
    iftDestroyDblArray(&region_sizes_prob);
    iftDestroyImage(&border_image);
    iftDestroyImage(&bin_border_image);
    iftDestroyAdjRel(&A);

    return focus;
}

iftSaliencyArray *iftImageCenterSaliencyPrior(iftImage *superpixel_img, float variance){
    iftIntArray *region_sizes;
    iftFloatArray *region_sizes_prob;
    iftSaliencyArray *region_to_center;

    region_sizes = iftCountLabelSpels(superpixel_img);
    region_sizes_prob = iftCreateFloatArray(region_sizes->n-1);

    region_to_center = iftCreateSaliencyArray((int)region_sizes_prob->n);

    for(int p = 0; p < superpixel_img->n; p++){
        iftVoxel u;
        u = iftGetVoxelCoord(superpixel_img, p);

        float normalized_x, normalized_y, normalized_x_center, normalized_y_center;
        normalized_x = (float) u.x / (float)(superpixel_img->xsize-1);
        normalized_y = (float) u.y / (float)(superpixel_img->ysize-1);
        normalized_x_center =  (((float)superpixel_img->xsize/2)-1) / (float)(superpixel_img->xsize-1);
        normalized_y_center =  (((float)superpixel_img->ysize/2)-1) / (float)(superpixel_img->ysize-1);
        region_to_center->val[superpixel_img->val[p]-1]+= sqrtf(square((normalized_x - normalized_x_center)) + square(normalized_y - normalized_y_center));
    }
    for(int i = 0; i < region_sizes_prob->n; i++) {
        region_to_center->val[i] = expf(-(region_to_center->val[i] / (float)region_sizes->val[i + 1])/variance);
    }

    float max_distance = 0;
    float min_distance = IFT_INFINITY_FLT;

    for(int rj = 0; rj < region_sizes_prob->n; rj++){ //Get Max
        if(region_to_center->val[rj] > max_distance)
            max_distance = region_to_center->val[rj];
        if(region_to_center->val[rj] < min_distance)
            min_distance = region_to_center->val[rj];
    }

    for(int rj = 0; rj < region_sizes_prob->n; rj++)
        region_to_center->val[rj] = (region_to_center->val[rj] - min_distance) /(max_distance-min_distance);

    iftDestroyIntArray(&region_sizes);
    iftDestroyFloatArray(&region_sizes_prob);

    region_to_center->superpixel_image = superpixel_img;

    return region_to_center;
}

iftSaliencyArray *iftForegroundDistanceSaliencyPrior(iftImage *superpixel_img, iftImage *saliency, float variance){
    iftIntArray *region_sizes;
    iftFloatArray *region_sizes_prob;

    region_sizes = iftCountLabelSpels(superpixel_img);
    region_sizes_prob = iftCreateFloatArray(region_sizes->n-1);

    iftSaliencyArray *region_to_foreground = iftCreateSaliencyArray(region_sizes_prob->n);
    iftImage *segmented = iftCreateImageFromImage(saliency);

    int otsu = iftOtsu(saliency);
    for(int p = 0; p < segmented->n; p++)
        if(saliency->val[p] >= otsu)
            segmented->val[p] = 1;
        else
            segmented->val[p] = 0;

    iftAdjRel *A = iftCircular(sqrtf(10));
    iftImage *mask = iftSelectImageDomain(superpixel_img->xsize, superpixel_img->ysize, superpixel_img->zsize);
    iftDestroyImage(&mask);
    iftDestroyAdjRel(&A);

    int number_foreground = 0;
    for(int p = 0; p < segmented->n; p++)
        if(segmented->val[p] == 1)
            number_foreground+=1;
    /*-------------- Foreground Scribbles ---------------*/
    iftVoxelArray *centroids = iftGeometricCentersFromLabelImage(superpixel_img);
    iftIntArray *query = iftCreateIntArray(centroids->n);
    for(int s = 0; s < region_to_foreground->n; s++)
        region_to_foreground->val[s] = IFT_INFINITY_FLT;
    if(number_foreground > 0) {

        for (int p = 0; p < superpixel_img->n; p++)
            if (segmented->val[p] == 1)
                query->val[superpixel_img->val[p] - 1] += 1;

        for(int s = 0; s < query->n; s++)
            if(((float)query->val[s]/(float)region_sizes->val[s + 1]) > 0.5 )
                query->val[s] = 1;
            else
                query->val[s] = 0;

        for (int p = 0; p < superpixel_img->n; p++) {
            iftVoxel u;
            u = iftGetVoxelCoord(superpixel_img, p);

            float normalized_x, normalized_y;
            normalized_x = (float) u.x / (float)(superpixel_img->xsize - 1);
            normalized_y = (float) u.y / (float)(superpixel_img->ysize - 1);
            for (int s = 0; s < centroids->n; s++) {
                if(query->val[s]) {
                    float normalized_x_scribble, normalized_y_scribble, distance;
                    normalized_x_scribble = (float) centroids->val[s].x / (float) (superpixel_img->xsize - 1);
                    normalized_y_scribble = (float) centroids->val[s].y / (float) (superpixel_img->ysize - 1);
                    distance = sqrtf(
                            square((normalized_x - normalized_x_scribble)) +
                            square(normalized_y - normalized_y_scribble));
                    region_to_foreground->val[superpixel_img->val[p] - 1] = iftMin(
                            region_to_foreground->val[superpixel_img->val[p] - 1], distance);
                }
            }

        }
    }else{
        iftDestroyIntArray(&query);
        iftDestroyVoxelArray(&centroids);
        iftDestroyImage(&segmented);
        iftDestroyIntArray(&region_sizes);
        iftDestroyFloatArray(&region_sizes_prob);
        return iftImageCenterSaliencyPrior(superpixel_img, variance);
    }
    iftDestroyVoxelArray(&centroids);

    for(int i = 0; i < region_sizes_prob->n; i++) {
//        printf("to_foreground %f \n", region_to_foreground->val[i]);
        region_to_foreground->val[i] = expf(-region_to_foreground->val[i] / variance);
    }

    iftDestroyImage(&segmented);
    float max_distance = 0;
    float min_distance = IFT_INFINITY_FLT;

    for(int rj = 0; rj < region_sizes_prob->n; rj++){ //Get Max
        if(region_to_foreground->val[rj] > max_distance)
            max_distance = region_to_foreground->val[rj];
        if(region_to_foreground->val[rj] < min_distance)
            min_distance = region_to_foreground->val[rj];
    }

    for(int rj = 0; rj < region_sizes_prob->n; rj++)
        region_to_foreground->val[rj] = (region_to_foreground->val[rj] - min_distance) /(max_distance-min_distance);

    iftDestroyIntArray(&region_sizes);
    iftDestroyFloatArray(&region_sizes_prob);

    region_to_foreground->superpixel_image = superpixel_img;

    return region_to_foreground;
}

iftSaliencyArray *iftForegroundDistanceSaliencyPriorNew(iftImage *superpixel_img, iftImage *saliency, float variance){
    iftIntArray *region_sizes;
    iftFloatArray *region_sizes_prob;

    region_sizes = iftCountLabelSpels(superpixel_img);
    region_sizes_prob = iftCreateFloatArray(region_sizes->n-1);

    iftSaliencyArray *region_to_foreground = iftCreateSaliencyArray(region_sizes_prob->n);
    iftImage *segmented = iftCreateImageFromImage(saliency);

    int otsu = iftOtsu(saliency);
    for(int p = 0; p < segmented->n; p++)
        if(saliency->val[p] >= otsu)
            segmented->val[p] = 1;
        else
            segmented->val[p] = 0;

    iftAdjRel *A = iftCircular(sqrtf(10));
    iftImage *mask = iftSelectImageDomain(superpixel_img->xsize, superpixel_img->ysize, superpixel_img->zsize);
    iftDestroyImage(&mask);
    iftDestroyAdjRel(&A);

    int number_foreground = 0;
    for(int p = 0; p < segmented->n; p++)
        if(segmented->val[p] == 1)
            number_foreground+=1;
    /*-------------- Foreground Scribbles ---------------*/
    iftVoxelArray *centroids = iftGeometricCentersFromLabelImage(superpixel_img);
    iftIntArray *query = iftCreateIntArray(centroids->n);
    for(int s = 0; s < region_to_foreground->n; s++)
        region_to_foreground->val[s] = IFT_INFINITY_FLT;
    if(number_foreground > 0) {

        for (int p = 0; p < superpixel_img->n; p++)
            if (segmented->val[p] == 1)
                query->val[superpixel_img->val[p] - 1] += 1;

        for(int s = 0; s < query->n; s++)
            if(((float)query->val[s]/(float)region_sizes->val[s + 1]) > 0.5)
                query->val[s] = 1;
            else
                query->val[s] = 0;

        for (int p = 0; p < superpixel_img->n; p++) {
            iftVoxel u;
            u = iftGetVoxelCoord(superpixel_img, p);

            float normalized_x, normalized_y;
            normalized_x = (float) u.x / (float)(superpixel_img->xsize - 1);
            normalized_y = (float) u.y / (float)(superpixel_img->ysize - 1);
            for (int s = 0; s < centroids->n; s++) {
                if(query->val[s]) {
                    float normalized_x_scribble, normalized_y_scribble, distance;
                    normalized_x_scribble = (float) centroids->val[s+1].x / (float) (superpixel_img->xsize - 1);
                    normalized_y_scribble = (float) centroids->val[s+1].y / (float) (superpixel_img->ysize - 1);
                    distance = sqrtf(
                            square((normalized_x - normalized_x_scribble)) +
                            square(normalized_y - normalized_y_scribble));
                    region_to_foreground->val[superpixel_img->val[p] - 1] = iftMin(
                            region_to_foreground->val[superpixel_img->val[p] - 1], distance);
                }
            }

        }
    }else{
        iftDestroyIntArray(&query);
        iftDestroyVoxelArray(&centroids);
        iftDestroyImage(&segmented);
        iftDestroyIntArray(&region_sizes);
        iftDestroyFloatArray(&region_sizes_prob);
        return iftImageCenterSaliencyPrior(superpixel_img, variance);
    }
    iftDestroyVoxelArray(&centroids);

    for(int s = 0; s < region_sizes_prob->n; s++) {
//        printf("to_foreground %f \n", region_to_foreground->val[i]);
        region_to_foreground->val[s] = expf(-region_to_foreground->val[s] / variance);
    }

    iftDestroyImage(&segmented);
    float max_distance = 0;
    float min_distance = IFT_INFINITY_FLT;

    for(int s = 0; s < region_sizes_prob->n; s++){ //Get Max
        if(region_to_foreground->val[s] > max_distance)
            max_distance = region_to_foreground->val[s];
        if(region_to_foreground->val[s] < min_distance)
            min_distance = region_to_foreground->val[s];
    }

    for(int s = 0; s < region_sizes_prob->n; s++)
        region_to_foreground->val[s] = (region_to_foreground->val[s] - min_distance) /(max_distance-min_distance);

    iftDestroyIntArray(&region_sizes);
    iftDestroyFloatArray(&region_sizes_prob);

    region_to_foreground->superpixel_image = superpixel_img;

    return region_to_foreground;
}

iftSaliencyArray *iftForegroundFeatureDistanceSaliencyPrior(iftImage *superpixel_img, iftImage *saliency, iftMImage *normalized_feats, float variance){
    iftIntArray *region_sizes;
    region_sizes = iftCountLabelSpels(superpixel_img);
    int number_superpixels = (int)region_sizes->n-1;

    iftSaliencyArray *feature_based_saliency = iftCreateSaliencyArray(number_superpixels);
    iftFloatArray *band_weight = iftCreateFloatArray((int)normalized_feats->m);

    iftFImage *normalized_saliency = iftCreateFImage(saliency->xsize, saliency->ysize, saliency->zsize);

    int min_saliency = IFT_INFINITY_INT, max_saliency = 0;

    for(int p = 0; p < saliency->n; p++){
        if(saliency->val[p] < min_saliency)
            min_saliency = saliency->val[p];
        if(saliency->val[p] > max_saliency)
            max_saliency = saliency->val[p];
    }

    for(int p = 0; p < saliency->n; p++)
        normalized_saliency->val[p] = (float)(saliency->val[p] - min_saliency) / (float)(max_saliency - min_saliency);

    int foreground_counter = 0;
    int threshold = iftOtsu(saliency);

    for(int p = 0; p < normalized_feats->n; p++)
        if(saliency->val[p] >= threshold){
            foreground_counter+=1;
            for (int b = 0; b < normalized_feats->m; b++){
                band_weight->val[b] += square(normalized_feats->val[p][b] - normalized_saliency->val[p]);
            }
        }

    for (int b = 0; b < normalized_feats->m; b++) {
        band_weight->val[b] = sqrtf(band_weight->val[b]) / (float) foreground_counter;
    }

    double max_weight = 0, min_weight = IFT_INFINITY_DBL;
    for(int b = 0; b < band_weight->n; b++){
        if(band_weight->val[b] > max_weight)
            max_weight = band_weight->val[b];
        if(band_weight->val[b] < min_weight)
            min_weight = band_weight->val[b];

    }

    for(int b = 0; b < band_weight->n; b++) {
        band_weight->val[b] = (float)(2 * (band_weight->val[b] - min_weight) / (max_weight - min_weight)) - 1;
//        printf("band weight[%d] %f\n", b+1, band_weight->val[b]);
    }

    for(int p = 0; p < saliency->n; p ++){
        int s = superpixel_img->val[p] - 1;
        for(int b = 0; b < normalized_feats->m; b++)
            feature_based_saliency->val[s] += normalized_feats->val[p][b] * band_weight->val[b];
    }

    for(int s = 0; s < number_superpixels; s++) {
        feature_based_saliency->val[s] = feature_based_saliency->val[s] / ((float) region_sizes->val[s + 1] * (float) normalized_feats->m);
    }

    float max_feature = IFT_INFINITY_FLT_NEG, min_feature = IFT_INFINITY_FLT;
    for(int s = 0; s < number_superpixels; s++){
        if(feature_based_saliency->val[s] > max_feature)
            max_feature = feature_based_saliency->val[s];
        if(feature_based_saliency->val[s] < min_feature)
            min_feature = feature_based_saliency->val[s];
    }


    for(int s = 0; s < number_superpixels; s++) {
        feature_based_saliency->val[s] = (feature_based_saliency->val[s] -min_feature) / (max_feature - min_feature);
    }

    for(int s = 0; s < number_superpixels; s++) {
//        printf("feature_based sal [%d] %f \n", s+1, feature_based_saliency->val[s]);
        feature_based_saliency->val[s] = expf(-feature_based_saliency->val[s] / variance);
//        printf("feature_based sal[%d] %f \n", s+1, feature_based_saliency->val[s]);
    }

    iftDestroyFloatArray(&band_weight);
    iftDestroyIntArray(&region_sizes);
    iftDestroyFImage(&normalized_saliency);

    return feature_based_saliency;
}

iftSaliencyArray *iftEllipseMatchingSaliencyPrior(iftImage *superpixel_img, float variance, int min_size, int max_size){
    iftVoxelArray *centroids = iftGeometricCentersFromLabelImage(superpixel_img);

    iftTensorScale *tensor_scale = iftSuperpixelToTensorScale(superpixel_img, 20, min_size, max_size);
    iftSaliencyArray *ellipseMatch = iftCreateSaliencyArray(centroids->n-1);
    iftImage *ellipse_match_image = iftCreateImageFromImage(superpixel_img);
    for(int p = 0; p < superpixel_img->n; p++) {
        iftVoxel u = iftGetVoxelCoord(superpixel_img, p);
        int s = superpixel_img->val[p]-1;
        float ellipse_eq = sqrtf(square((u.x - tensor_scale->pos_focus->val[s].x)) + (square((u.y - tensor_scale->pos_focus->val[s].y))))
                           + sqrtf(square((u.x - tensor_scale->neg_focus->val[s].x)) + (square((u.y - tensor_scale->neg_focus->val[s].y))));
        if (ellipse_eq <= 2*tensor_scale->major_axis->val[s])
            ellipseMatch->val[s]+=1;
    }

    iftDestroyImage(&ellipse_match_image);

    iftIntArray *region_sizes = iftCountLabelSpels(superpixel_img);

    float min_ellipse_score = IFT_INFINITY_FLT;
    for(int s = 0; s < ellipseMatch->n; s++){
        if(tensor_scale->area->val[s] >= (float)min_size || tensor_scale->area->val[s] <= (float)max_size) {
            ellipseMatch->val[s] /= (float)region_sizes->val[s + 1];
//            printf("Match = %f\n", ellipseMatch->val[s]);
            ellipseMatch->val[s] = expf(ellipseMatch->val[s] / variance);
            if (ellipseMatch->val[s] < min_ellipse_score)
                min_ellipse_score = ellipseMatch->val[s];
        }
    }

    for(int s = 0; s < ellipseMatch->n; s++){
        if(tensor_scale->area->val[s] < (float)min_size || tensor_scale->area->val[s] > (float)max_size)
            ellipseMatch->val[s] = min_ellipse_score;
    }

    float max = 0;
    float min = IFT_INFINITY_FLT;
    for(int rj = 0; rj < ellipseMatch->n; rj++){ //Get Max
        if(ellipseMatch->val[rj] > max)
            max = ellipseMatch->val[rj];
        if(ellipseMatch->val[rj] < min)
            min = ellipseMatch->val[rj];
    }

    for(int rj = 0; rj < ellipseMatch->n; rj++)
        ellipseMatch->val[rj] = (ellipseMatch->val[rj] - min) /(max-min);

    iftDestroyVoxelArray(&centroids);
    iftDestroyTensorScale(&tensor_scale);
    iftDestroyIntArray(&region_sizes);

    ellipseMatch->superpixel_image = superpixel_img;
    return ellipseMatch;
}

iftSaliencyArray *iftImageScribbleDistanceSaliencyPrior(iftImage *superpixel_img, iftLabeledSet *scribbles, float variance){
    iftIntArray *region_sizes;
    iftFloatArray *region_sizes_prob;
    iftSaliencyArray *region_to_scribble;
    iftImage *scribble_image = iftCreateImage(superpixel_img->xsize, superpixel_img->ysize, superpixel_img->zsize);
    iftLabeledSetToImage(scribbles, scribble_image, 1); //background = 1, foreground = 2

    region_sizes = iftCountLabelSpels(superpixel_img);
    region_sizes_prob = iftCreateFloatArray(region_sizes->n-1);
    region_to_scribble = iftCreateSaliencyArray((int) region_sizes_prob->n);

    int number_of_scribbled_pixels = 0;
    for(int p = 0; p < superpixel_img->n; p++)
        if(scribble_image->val[p] == 2)
            number_of_scribbled_pixels++;

    /*-------------- Foreground Scribbles ---------------*/
    if(number_of_scribbled_pixels > 0) {
        iftVoxelArray *scribbles_voxels = iftCreateVoxelArray(number_of_scribbled_pixels);

        int scribble_number = 0;
        for (int p = 0; p < superpixel_img->n; p++)
            if (scribble_image->val[p] == 2) {
                iftVoxel u;
                u = iftGetVoxelCoord(superpixel_img, p);
                scribbles_voxels->val[scribble_number].x = u.x;
                scribbles_voxels->val[scribble_number].y = u.y;

                scribble_number++;
            }
        for (int p = 0; p < superpixel_img->n; p++) {
            iftVoxel u;
            u = iftGetVoxelCoord(superpixel_img, p);

            float normalized_x, normalized_y, min_distance = IFT_INFINITY_FLT;
            normalized_x = (float) u.x / (float)(superpixel_img->xsize - 1);
            normalized_y = (float) u.y / (float)(superpixel_img->ysize - 1);
            for (int scribble = 0; scribble < scribbles_voxels->n; scribble++) {
                float normalized_x_scribble, normalized_y_scribble, distance;
                normalized_x_scribble = (float) scribbles_voxels->val[scribble].x / (float)(superpixel_img->xsize - 1);
                normalized_y_scribble = (float) scribbles_voxels->val[scribble].y / (float)(superpixel_img->ysize - 1);
                distance = sqrtf(
                        square((normalized_x - normalized_x_scribble)) + square(normalized_y - normalized_y_scribble));

                if (distance < min_distance)
                    min_distance = distance;
            }

            region_to_scribble->val[superpixel_img->val[p] - 1] += min_distance;
        }

        iftDestroyVoxelArray(&scribbles_voxels);
    }
    /*-------------- Background Scribbles ---------------*/

    if(0) {
        for (int p = 0; p < superpixel_img->n; p++)
            if (scribble_image->val[p] == 1)
                number_of_scribbled_pixels++;
        if (number_of_scribbled_pixels > 0) {
            iftVoxelArray *scribbles_voxels = iftCreateVoxelArray(number_of_scribbled_pixels);

            int scribble_number = 0;
            float scribble_sum_x = 0, scribble_sum_y = 0;
            for (int p = 0; p < superpixel_img->n; p++)
                if (scribble_image->val[p] == 1) {
                    iftVoxel u;
                    u = iftGetVoxelCoord(superpixel_img, p);
                    scribbles_voxels->val[scribble_number].x = u.x;
                    scribbles_voxels->val[scribble_number].y = u.y;

                    scribble_sum_x += (float)u.x;
                    scribble_sum_y += (float)u.y;

                    scribble_number++;
                }

            //    scribble_mean_point_x = scribble_sum_x/scribble_number;
            //    scribble_mean_point_y = scribble_sum_y/scribble_number;

            //
            for (int p = 0; p < superpixel_img->n; p++) {
                iftVoxel u;
                u = iftGetVoxelCoord(superpixel_img, p);

                float normalized_x, normalized_y, min_distance = IFT_INFINITY_FLT;
                normalized_x = (float) u.x / (float)(superpixel_img->xsize - 1);
                normalized_y = (float) u.y / (float)(superpixel_img->ysize - 1);
                for (int scribble = 0; scribble < scribbles_voxels->n; scribble++) {
                    float normalized_x_scribble, normalized_y_scribble, distance;
                    normalized_x_scribble = (float) scribbles_voxels->val[scribble].x / (float)(superpixel_img->xsize - 1);
                    normalized_y_scribble = (float) scribbles_voxels->val[scribble].y / (float)(superpixel_img->ysize - 1);
                    distance = sqrtf(
                            square((normalized_x - normalized_x_scribble)) +
                            square(normalized_y - normalized_y_scribble));

                    if (distance < min_distance)
                        min_distance = distance;
                }

                //        float normalized_x_scribble = scribble_mean_point_x/(superpixel_img->xsize-1);
                //        float normalized_y_scribble = scribble_mean_point_y/(superpixel_img->ysize-1);

                //        min_distance = sqrtf(square((normalized_x - normalized_x_scribble)) + square(normalized_y - normalized_y_scribble));
                region_to_scribble->val[superpixel_img->val[p] - 1] -= min_distance;
            }
        }
    }

    for(int i = 0; i < region_sizes_prob->n; i++)
        region_to_scribble->val[i] = expf(-(region_to_scribble->val[i] / (float)region_sizes->val[i + 1])/variance);

    float max_distance = 0;
    float min_distance = IFT_INFINITY_FLT;

    for(int rj = 0; rj < region_sizes_prob->n; rj++){ //Get Max
        if(region_to_scribble->val[rj] > max_distance)
            max_distance = region_to_scribble->val[rj];
        if(region_to_scribble->val[rj] < min_distance)
            min_distance = region_to_scribble->val[rj];
    }

    for(int rj = 0; rj < region_sizes_prob->n; rj++)
        region_to_scribble->val[rj] = (region_to_scribble->val[rj] - min_distance) /(max_distance-min_distance);

    iftDestroyIntArray(&region_sizes);
    iftDestroyFloatArray(&region_sizes_prob);

    region_to_scribble->superpixel_image = superpixel_img;

    return region_to_scribble;
}

void addAllPriorsToGraph(iftSaliencyGraph **saliency_graph_ref, iftImage *superpixel_img, iftMImage *features, iftImage *saliency, iftITSELFParameters *params){
    iftSaliencyGraph *saliency_graph = *saliency_graph_ref;
    float center_variance = params->prior_params.center_variance, focus_variance = params->prior_params.focus_variance, ellipse_variance = params->prior_params.ellipse_variance,
            foreground_distance_variance = params->prior_params.foreground_distance_variance, foreground_feature_variance = params->prior_params.foreground_feature_variance;
    if(foreground_feature_variance > 0) {
        iftSaliencyArray *foreground_feature_prior = iftForegroundFeatureDistanceSaliencyPrior(superpixel_img, saliency, features, foreground_feature_variance);
        iftAddPriorToGraph(&saliency_graph, foreground_feature_prior);
        //iftWriteSaliencyImage(superpixel_img, foreground_feature_prior, "priors/foreground_feature_variance.png", 0);
    }
    if(foreground_distance_variance > 0) {
        iftSaliencyArray *foreground_distance_prior = iftForegroundDistanceSaliencyPrior(superpixel_img, saliency, foreground_distance_variance);
        iftAddPriorToGraph(&saliency_graph, foreground_distance_prior);
        iftWriteSaliencyImage(superpixel_img, foreground_distance_prior, "priors/foreground_distance_prior.png", 0);
    }
    if(center_variance > 0){
        iftSaliencyArray *center_prior = iftImageCenterSaliencyPrior(superpixel_img, center_variance);
        iftAddPriorToGraph(&saliency_graph, center_prior);
//        iftWriteSaliencyImage(superpixel_img, center_prior, "priors/ellipse_prior.png", 0);
    }
    if(focus_variance > 0){
        iftSaliencyArray *focus_prior = iftFocusPrior(superpixel_img, features, focus_variance);
        iftAddPriorToGraph(&saliency_graph, focus_prior);
//        iftWriteSaliencyImage(superpixel_img, focus_prior, "priors/focus_prior.png", 0);
    }
    if(ellipse_variance > 0){
        int min_size = 400, max_size = 3000;
        iftSaliencyArray *ellipse_prior = iftEllipseMatchingSaliencyPrior(superpixel_img, ellipse_variance, min_size, max_size);
        iftAddPriorToGraph(&saliency_graph, ellipse_prior);
//        iftWriteSaliencyImage(superpixel_img, ellipse_prior, "priors/ellipse_prior.png", 0);
    }

}


void iftAddPriorToGraph(iftSaliencyGraph **ref_saliency_graph, iftSaliencyArray *new_prior){
    iftSaliencyGraph *saliency_graph = *ref_saliency_graph;

    iftAppendPriorToList(&(saliency_graph->prior_list), new_prior, saliency_graph->npriors);

    saliency_graph->npriors+=1;
}

void iftAppendPriorToList(iftSaliencyArray ***prior_list, iftSaliencyArray *new_prior, int npriors){
    iftSaliencyArray **new_prior_list = (iftSaliencyArray**) iftAlloc((size_t) npriors+1, sizeof(iftSaliencyArray*));

    iftSaliencyArray **old_prior_list = *prior_list;

    for(int i = 0; i<npriors; i++)
        new_prior_list[i] = old_prior_list[i];
    new_prior_list[npriors] = new_prior;
    iftFree(old_prior_list);

    *prior_list = new_prior_list;
}

void iftResetPriorList(iftSaliencyGraph **ref_saliency_graph){
    iftSaliencyGraph *saliency_graph = *ref_saliency_graph;
    iftSaliencyArray *aux;
    for(int i = 0; i < saliency_graph->npriors; i++){
        aux = saliency_graph->prior_list[i];
        iftDestroySaliencyArray(&aux);
    }
    iftFree(saliency_graph->prior_list);
    saliency_graph->npriors = 0;

    for(int p = 0; p < saliency_graph->combined_priors->n; p++)
        saliency_graph->combined_priors->val[p] = 1;
}

void iftDestroyPriorList(iftSaliencyGraph **ref_saliency_graph){
    iftSaliencyGraph *saliency_graph = *ref_saliency_graph;
    iftSaliencyArray *aux;
    iftImage *aux_img;
    if(saliency_graph->prior_list != NULL){
        for(int i = 0; i < saliency_graph->npriors; i++){
            aux = saliency_graph->prior_list[i];
            if(i %2 == 1){
                aux_img = aux->superpixel_image;
                iftDestroyImage(&aux_img);
            }
            iftDestroySaliencyArray(&aux);
        }
        iftFree(saliency_graph->prior_list);
    }
    saliency_graph->npriors = 0;
    saliency_graph->prior_list = NULL;

    if(saliency_graph->combined_priors != NULL){
        for(int p = 0; p < saliency_graph->combined_priors->n; p++)
            saliency_graph->combined_priors->val[p] = 1;
    }
}


/* Ellipse Matching Auxiliary Functions */
iftSet *iftBoundingBoxToSeedSet(iftBoundingBox bb, iftImage *image){
    iftSet *seeds = NULL;
    for(int p = 0; p < image->n; p++){
        iftVoxel t = iftGetVoxelCoord(image, p);
        if(t.x > bb.begin.x && t.x < bb.end.x && t.y > bb.begin.y && t.y < bb.end.y) {
            if (image->val[p] == 1)
                iftInsertSet(&seeds, p);
        }else if( ((t.x == bb.begin.x-1 || t.x == bb.end.x+1) && (t.y >= bb.begin.y-1 && t.y <= bb.end.y+1)) ||
                  ((t.y == bb.begin.y-1 || t.y == bb.end.y+1) && (t.x >= bb.begin.x-1 && t.x <= bb.end.x+1)) )
            iftInsertSet(&seeds, p);

    }
    return seeds;
}

iftImage *iftFillEllipse(iftImage *label_img, int label, iftVoxel *new_centroid, float *new_area){
    int adj_size = 20;
    iftImage *object = iftCreateImageFromImage(label_img);
    for(int p = 0; p < label_img->n; p++)
        if(label_img->val[p] == label)
            object->val[p] = 1;

    iftVoxel box_center;
    iftBoundingBox bb = iftMinBoundingBox(object, &box_center);
    iftSet *seeds = iftBoundingBoxToSeedSet(bb, object);

    bb.begin.x = iftMax(bb.begin.x - adj_size, 0);
    bb.begin.y = iftMax(bb.begin.y - adj_size, 0);
    bb.end.x =iftMin(adj_size + bb.end.x, label_img->xsize-1);
    bb.end.y =iftMin(adj_size + bb.end.y, label_img->ysize-1);

    iftImage *mask = iftCreateImageFromImage(label_img);
    iftFillBoundingBoxInImage(mask, bb, 1);

    iftImage *filled = iftCloseBasins(object, seeds, mask);
    iftDestroyImage(&mask);
    iftDestroySet(&seeds);

    iftVoxelArray *centroids = iftGeometricCentersFromLabelImage(filled);
    new_centroid->x = centroids->val[1].x;
    new_centroid->y = centroids->val[1].y;
    new_centroid->z = centroids->val[1].z;


    float area = 0;
    for(int p = 0; p < filled->n; p++){
        area+=(float)filled->val[p];
    }
    *new_area = area;

    iftDestroyImage(&object);
    iftDestroyVoxelArray(&centroids);

    return filled;


}

iftTensorScale *iftSuperpixelToTensorScale(iftImage *label_img, int m_pairs, int min_size, int max_size){
    iftImage *distance_image;
    iftImage *sqrt_distance_image;
    iftTensorScale *ts;
    iftIntArray *region_sizes = iftCountLabelSpels(label_img);

    iftVoxelArray *centroids = iftGeometricCentersFromLabelImage(label_img);

    iftAdjRel *A = iftCircular((float)1.0);
    distance_image = iftEuclDistTrans(label_img, A, IFT_INTERIOR, NULL, NULL, NULL);
    iftDestroyAdjRel(&A);

    sqrt_distance_image = iftCreateImageFromImage(label_img);
    for(int p = 0; p < distance_image->n; p++)
        sqrt_distance_image->val[p] = (int) sqrt(distance_image->val[p]);


    ts = (iftTensorScale *)malloc(sizeof(iftTensorScale));
    ts->orientation = iftCreateFloatArray(centroids->n - 1);
    ts->anisotropy = iftCreateFloatArray(centroids->n - 1);
    ts->major_axis = iftCreateFloatArray(centroids->n - 1);
    ts->minor_axis = iftCreateFloatArray(centroids->n - 1);
    ts->m_pairs     = m_pairs;
    ts->pos_focus = iftCreateVoxelArray(centroids->n-1);
    ts->neg_focus = iftCreateVoxelArray(centroids->n-1);
    ts->area = iftCreateFloatArray(centroids->n-1);
    ts->n = centroids->n-1;

    iftVector *tau = (iftVector *)malloc(sizeof(iftVector)*m_pairs);
    iftPoint *epsilon = (iftPoint *)malloc(sizeof(iftPoint)*m_pairs);

    float teta = (float)0.0;
    for(int i=0;i<m_pairs;i++){
        tau[i].x = cosf(teta);
        tau[i].y = sinf(teta);
        tau[i].z = (float)0.0;

        teta += ((float)PI/(float)m_pairs);
    }

    int min_dist_to_border, v, d,d1,d2;
    float gSxy, gSy2_x2, x_displacement,y_displacement,xc,yc,cosx, siny, aux, acc, wt;
    iftVoxel line;

    float e = (float)0.000001;
    iftImage *filled;
    iftImage *sqrt_distance_current;
    for(int l = 1; l < centroids->n; l++){
        if(region_sizes->val[l] >= min_size && region_sizes->val[l] <= max_size){
            iftVoxel new_centroid;
            float new_area;
            filled = iftFillEllipse(label_img, l, &new_centroid, &new_area);
            centroids->val[l] = new_centroid;
            ts->area->val[l-1] = new_area;
            iftDestroyImage(&distance_image);
            distance_image = iftEuclDistTrans(filled, A, IFT_INTERIOR, NULL, NULL, NULL);
            sqrt_distance_current = iftCreateImageFromImage(distance_image);
            for(int p = 0; p < distance_image->n; p++)
                sqrt_distance_current->val[p] = (int) sqrt(distance_image->val[p]);
        } else {
            filled = iftCreateImageFromImage(label_img);
            for (int p = 0; p < label_img->n; p++) {
                if (label_img->val[p] == l)
                    filled->val[p] = 1;
                else
                    filled->val[p] = 0;
            }
            sqrt_distance_current = sqrt_distance_image;
        }

        iftVoxel centroid = centroids->val[l];
        int p = iftGetVoxelIndex(filled, centroid);

        min_dist_to_border = sqrt_distance_current->val[p];
        gSxy = gSy2_x2 = (float)0.0;
        xc = (float)centroid.x+(float)0.5;
        yc = (float)centroid.y+(float)0.5;

        for(int k = 0; k < m_pairs; k++){
            cosx = tau[k].x;
            siny = tau[k].y;
            v = min_dist_to_border;
            d1 = d2 = 0;
            while(1){
                x_displacement = (float)v*cosx;
                y_displacement = (float)v*siny;

                if(d1==0){
                    line.x = (int)(x_displacement+xc);
                    line.y = (int)(y_displacement+yc);
                    line.z = 0;
                    d1 = sqrt_distance_current->val[iftGetVoxelIndex(filled, line)];
                    if(d1 == 0)
                        break;
                }

                if(d2==0){
                    line.x = (int)(xc-x_displacement);
                    line.y = (int)(yc-y_displacement);
                    line.z = 0;
                    d2 = sqrt_distance_current->val[iftGetVoxelIndex(filled, line)];
                    if(d2 == 0)
                        break;
                }

                d = (d1<d2)?d1:d2;
                d1 -= d;
                d2 -= d;
                v += d;
            }

            epsilon[k].x = x_displacement;
            epsilon[k].y = -y_displacement;
            gSxy -= x_displacement*y_displacement;            //gSxy += x*(-y);
            gSy2_x2 += (y_displacement+x_displacement)*(y_displacement-x_displacement); //(y*y-x*x);

        }
        /*--------- TETA ------------ */
        if(gSy2_x2==0.0){
            if(gSxy>0.0)
                teta=(float)(PI/2.0);
            else
                teta=(float)(-PI/2.0);
        }
        else{
            teta = atanf((gSxy+gSxy)/gSy2_x2);

            if(gSxy<0.0 && teta<0.0)
                teta+=(float)PI;
            else if(gSxy>0.0 && teta>0.0)
                teta-=(float)PI;
            else if(teta==0.0 && gSy2_x2>0.0)
                teta=(float)PI;
        }
        teta /= (float)2.0;
        //----------------- A & B ---------------------------------
        float minor_axis_squared = (float)distance_image->val[p] + e; //b2
        float minor_axis   = sqrtf((int)minor_axis_squared); //b1

        acc = wt = (float)0.0;
        float sin_teta = sinf(teta);
        float cos_teta = cosf(teta);
        float rotated_y_displacement, rotated_x_displacement, rotated_y_squared, rotated_x_squared, aa, w, major_axis_squared;
        for(int k = 0; k < m_pairs; k++){
            x_displacement = epsilon[k].x;
            y_displacement = epsilon[k].y;
            rotated_x_displacement = x_displacement*cos_teta - y_displacement*sin_teta;
            rotated_y_displacement = y_displacement*cos_teta + x_displacement*sin_teta;
            rotated_x_squared = square(rotated_x_displacement);
            rotated_y_squared = square(rotated_y_displacement);
            if(rotated_y_squared<minor_axis_squared){
                aa = minor_axis_squared*rotated_x_squared/(minor_axis_squared-rotated_y_squared);
//                printf("k = %d aa = %f rotated_y_squared = %f minor_axis_squared = %f\n", k, aa, rotated_y_squared, minor_axis_squared);
                if(aa>=minor_axis_squared){
                    w = (rotated_y_displacement<0.0)?(minor_axis+rotated_y_displacement):(minor_axis-rotated_y_displacement);
                    acc += w*aa;
                    wt  += w;
                }
            } else {
                aa = minor_axis_squared * rotated_x_squared / (minor_axis_squared - rotated_y_squared);
            }

            if(region_sizes->val[l] >= min_size && region_sizes->val[l] <= max_size)
                iftDestroyImage(&sqrt_distance_current);
        }


        if(wt>0.0)
            major_axis_squared = acc / wt;
        else
            major_axis_squared = minor_axis_squared;

        aux = (float)1.0-minor_axis_squared/major_axis_squared;

        if(aux<0.0)
            aux = (float)0.0;

        ts->anisotropy->val[l-1] = sqrtf(aux);
        ts->minor_axis->val[l-1] = minor_axis; //sqrtf(b2);
        ts->major_axis->val[l-1] = sqrtf(major_axis_squared);
        ts->area->val[l-1] = IFT_PI * ts->major_axis->val[l-1] * ts->minor_axis->val[l-1];
//            if(l == 185)
//                printf("Anis %f MIN %f MAJ %f AREA %d MAX_AREA %d\n ", ts->anisotropy->val[l-1], ts->minor_axis->val[l-1], ts->major_axis->val[l-1], ts->area->val[l-1], max_size);

        if(teta<0.0)
            teta+=(float)PI;
        if(teta>PI)
            teta = (float) PI;
        teta = (float)PI-teta;

        ts->orientation->val[l-1] = teta;
        iftDestroyImage(&filled);
        if(region_sizes->val[l] >= min_size && region_sizes->val[l] <= max_size)
            iftDestroyImage(&sqrt_distance_current);
    }

    for(int superpixel = 0; superpixel < centroids->n-1; superpixel++) {
        iftVoxel c = centroids->val[superpixel + 1];
        iftVoxel focus_pos, focus_neg; //,vertex, co_vertex;
        float focus_distance = sqrtf(
                square(ts->major_axis->val[superpixel]) - square(ts->minor_axis->val[superpixel]));
        double theta = ts->orientation->val[superpixel];
        focus_pos.x = c.x + (int) (focus_distance * cos(theta));
        focus_pos.y = c.y + (int) (focus_distance * -sin(theta));
        focus_neg.x = c.x - (int) (focus_distance * cos(theta));
        focus_neg.y = c.y - (int) (focus_distance * -sin(theta));
        ts->pos_focus->val[superpixel] = focus_pos;
        ts->neg_focus->val[superpixel] = focus_neg;
    }


    free(tau);
    free(epsilon);
    iftDestroyImage(&sqrt_distance_image);
    iftDestroyImage(&distance_image);
    iftDestroyIntArray(&region_sizes);

    iftDestroyVoxelArray(&centroids);



    return ts;
}

void iftDestroyTensorScale(iftTensorScale **tensor_scale_ref){
    iftTensorScale *ts = *tensor_scale_ref;
    if(ts != NULL){
        iftDestroyFloatArray(&ts->minor_axis);
        iftDestroyFloatArray(&ts->orientation);
        iftDestroyFloatArray(&ts->anisotropy);
        iftDestroyFloatArray(&ts->major_axis);
        iftDestroyVoxelArray(&ts->neg_focus);
        iftDestroyVoxelArray(&ts->pos_focus);
        iftDestroyFloatArray(&ts->area);
    }
    free(ts);
    *tensor_scale_ref = NULL;

}