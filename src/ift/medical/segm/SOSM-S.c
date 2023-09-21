#include "ift/medical/segm/SOSM-S.h"

#include "ift/core/dtypes/Dict.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Json.h"
#include "ift/core/io/NumPy.h"
#include "ift/medical/segm/AdaPro.h"
#include "iftParamOptimizationProblems.h"
#include "iftSegmentation.h"


/********************** PRIVATE FUNCTIONS *************************/
iftSOSMSObjModel *_iftCreateSOSMSObjModel(int label) {
    iftSOSMSObjModel *obj_model = (iftSOSMSObjModel*) calloc(1, sizeof(iftSOSMSObjModel));
    obj_model->label       = label;
    
    return obj_model;
}


void _iftDestroySOSMSObjModel(iftSOSMSObjModel **model) {
    if (model != NULL) {
        iftSOSMSObjModel *aux = *model;
        
        if (aux != NULL) {
            iftDestroyFImage(&aux->prob_atlas);
            iftFree(aux);
            *model = NULL;
        }
    }
}

iftSOSMSObjModel *_iftReadSOSMSObjModel(const char *obj_dir) {
    char *info_path = iftJoinPathnames(2, obj_dir, "info.json");
    iftDict *info   = iftReadJson(info_path);
    
    int label                   = iftGetLongValFromDict("label", info);
    iftSOSMSObjModel *obj_model = _iftCreateSOSMSObjModel(label);
    
    iftIntArray *shape              = iftGetIntArrayFromDict("template-shape", info);
    obj_model->template_shape.xsize = shape->val[0];
    obj_model->template_shape.ysize = shape->val[1];
    obj_model->template_shape.zsize = shape->val[2];
    iftDestroyIntArray(&shape);
    
    char *prob_atlas_path = iftJoinPathnames(2, obj_dir, "prob_atlas.npy");
    if (!iftFileExists(prob_atlas_path))
        iftError("Prob Atlas \"prob_atlas.npy\" not found for object model with label %d",
                 "_iftReadSOSMSObjModel", label);
   obj_model->prob_atlas = iftReadFImage(prob_atlas_path);
    
    iftIntArray *begin = iftGetIntArrayFromDict("prob-atlas-location:begin", info);
    obj_model->begin   = iftIntArrayToVoxel(begin);
    iftDestroyIntArray(&begin);
    
    iftIntArray *search_region_begin = iftGetIntArrayFromDict("prob-atlas-location:search-region:begin", info);
    obj_model->search_region.begin   = iftIntArrayToVoxel(search_region_begin);
    iftDestroyIntArray(&search_region_begin);
    
    iftIntArray *search_region_end = iftGetIntArrayFromDict("prob-atlas-location:search-region:end", info);
    obj_model->search_region.end   = iftIntArrayToVoxel(search_region_end);
    iftDestroyIntArray(&search_region_end);
    
    iftDestroyDict(&info);
    iftFree(info_path);
    iftFree(prob_atlas_path);
    
    return obj_model;
}


void _iftWriteSOSMSObjModel(iftSOSMSObjModel *obj_model, const char *out_dir) {
    if (!iftDirExists(out_dir))
        iftMakeDir(out_dir);
    
    // writes the displacement vector into Json Information file
    iftDict *info        = iftCreateDict();
    info->erase_elements = true;
    
    char key[128];
    strcpy(key, "label");
    iftInsertIntoDict(key, obj_model->label, info);
    
    iftIntArray *template_shape = iftCreateIntArray(3);
    template_shape->val[0] = obj_model->template_shape.xsize;
    template_shape->val[1] = obj_model->template_shape.ysize;
    template_shape->val[2] = obj_model->template_shape.zsize;
    iftInsertIntoDict("template-shape", template_shape, info);
    
    strcpy(key, "prob-atlas-location:begin");
    iftIntArray *begin = iftVoxelToIntArray(obj_model->begin);
    iftInsertIntoDict(key, begin, info);
    
    strcpy(key, "prob-atlas-location:search-region:begin");
    iftIntArray *search_region_begin = iftVoxelToIntArray(obj_model->search_region.begin);
    iftInsertIntoDict(key, search_region_begin, info);
    
    strcpy(key, "prob-atlas-location:search-region:end");
    iftIntArray *search_region_end = iftVoxelToIntArray(obj_model->search_region.end);
    iftInsertIntoDict(key, search_region_end, info);
    
    char *info_path = iftJoinPathnames(2, out_dir, "info.json");
    iftWriteJson(info, info_path);
    iftFree(info_path);
    
    char *prob_atlas_path = iftJoinPathnames(2, out_dir, "prob_atlas.npy");
//    iftWriteFImageAsNumPy(obj_model->prob_atlas, prob_atlas_path);
    iftFree(prob_atlas_path);
}


iftSOSMS *_iftCreateSOSMS(const iftIntArray *labels, const iftImage *template_img,
                          bool create_obj_models) {
    if (labels == NULL)
        iftError("Array of Labels is NULL", "_iftCreateSOSMS");
    if (template_img == NULL)
        iftError("Template Image is NULL", "_iftCreateSOSMS");
    
    iftSOSMS *sosms_s = (iftSOSMS*) calloc(1, sizeof(iftSOSMS));
    
    sosms_s->labels = iftCreateIntArray(labels->n);
    iftCopyIntArray(sosms_s->labels->val, labels->val, labels->n);
    
    sosms_s->template_img = iftCopyImage(template_img);
    
    sosms_s->obj_models = (iftSOSMSObjModel**) iftAlloc(labels->n, sizeof(iftSOSMSObjModel*));
    if (create_obj_models) {
        for (size_t o = 0; o < labels->n; o++) {
            sosms_s->obj_models[o] = _iftCreateSOSMSObjModel(labels->val[o]);
        }
    }
    
    return sosms_s;
}


iftSOSMS *_iftReadSOSMS(const char *path, char **tmp_dir) {
    *tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);
    iftUnzipFile(path, *tmp_dir);
    
    char *info_path = iftJoinPathnames(2, *tmp_dir, "info.json");
    if (!iftFileExists(info_path)) {
        iftRemoveDir(*tmp_dir);
        iftError("Info Json file \"info.json\" not found", "_iftReadSOSMS");
    }
    
    iftDict *info         = iftReadJson(info_path);
    info->erase_elements  = true;
    iftIntArray *labels = iftGetIntArrayFromDict("labels", info);
    
    char *template_path = iftJoinPathnames(2, *tmp_dir, "template.nii.gz");
    if (!iftFileExists(template_path)) {
        iftRemoveDir(*tmp_dir);
        iftError("Template \"template.nii.gz\" not found", "_iftReadSOSMS");
    }
    iftImage *template_img = iftReadImageByExt(template_path);
    iftFree(template_path);
    
    iftSOSMS *sosm_s  = _iftCreateSOSMS(labels, template_img, false);
    iftDestroyImage(&template_img);
    
    iftDestroyDict(&info);
    iftFree(info_path);
    
    char *obj_models_dir = iftJoinPathnames(2, *tmp_dir, "obj_models");
    
    for (int o = 0; o < sosm_s->labels->n; o++) {
        char label_str[32];
        sprintf(label_str, "%d", sosm_s->labels->val[o]);
        char *obj_dir = iftJoinPathnames(2, obj_models_dir, label_str);
        
        sosm_s->obj_models[o] = _iftReadSOSMSObjModel(obj_dir);
        
        iftFree(obj_dir);
    }
    iftFree(obj_models_dir);
    
    return sosm_s;
}


char *_iftWriteSOSMS(const iftSOSMS *sosm_s) {
    char *tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);
    
    // writing template
    char *template_path = iftJoinPathnames(2, tmp_dir, "template.nii.gz");
    iftWriteImageByExt(sosm_s->template_img, template_path);
    iftFree(template_path);
    
    // writing a Json with information
    char *info_path      = iftJoinPathnames(2, tmp_dir, "info.json");
    iftDict *info        = iftCreateDict();
    info->erase_elements = true;
    
    iftIntArray *labels = iftCreateIntArray(sosm_s->labels->n);
    iftCopyIntArray(labels->val, sosm_s->labels->val, sosm_s->labels->n);
    iftInsertIntoDict("labels", labels, info);
    
    iftWriteJson(info, info_path);
    iftDestroyDict(&info);
    iftFree(info_path);
    
    char *obj_models_dir = iftJoinPathnames(2, tmp_dir, "obj_models");
    iftMakeDir(obj_models_dir);
    
    for (int o = 0; o < sosm_s->labels->n; o++) {
        char obj_str[32];
        sprintf(obj_str, "%d", sosm_s->labels->val[o]);
        char *out_obj_dir = iftJoinPathnames(2, obj_models_dir, obj_str);
        _iftWriteSOSMSObjModel(sosm_s->obj_models[o], out_obj_dir);
        iftFree(out_obj_dir);
    }
    iftFree(obj_models_dir);
    
    return tmp_dir;
}



/********************** PUBLIC FUNCTIONS *************************/
void iftDestroySOSMS(iftSOSMS **sosm_s) {
    if (sosm_s != NULL) {
        iftSOSMS *aux = *sosm_s;
        
        if (aux != NULL) {
            int n_objs = aux->labels->n;
            
            iftDestroyIntArray(&aux->labels);
            iftDestroyImage(&aux->template_img);
            for (int o = 0; o < n_objs; o++) {
                _iftDestroySOSMSObjModel(&aux->obj_models[o]);
            }
            
            iftFree(aux);
            *sosm_s = NULL;
        }
    }
}


iftSOSMS *iftReadSOSMS(const char *path) {
    char *tmp_dir = NULL;
    
    iftSOSMS *sosm_s = _iftReadSOSMS(path, &tmp_dir);
    
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
    
    return sosm_s;
}


void iftWriteSOSMS(const iftSOSMS *sosm_s, const char *path) {
    char *tmp_dir = _iftWriteSOSMS(sosm_s);
    
    iftZipDirContent(tmp_dir, path);
    
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
}


iftSOSMS *iftTrainSOSMS(const iftFileSet *atlas_set, const iftImage *template_img,
                        const iftIntArray *labels) {
    iftSOSMS *sosm_s = _iftCreateSOSMS(labels, template_img, true);
    
    // mbb_cvs[o][i] = center voxel of the min. bounding boxel for the object labels->val[o]
    //                 from the atlas atlas_set->files[i]->path
    iftVoxel **mbb_cvs = iftAlloc(labels->n, sizeof(iftVoxel*));
    for (size_t o = 0; o < labels->n; o++) {
        mbb_cvs[o] = iftAlloc(atlas_set->n, sizeof(iftVoxel));
        sosm_s->obj_models[o]->prob_atlas = iftCreateFImage(template_img->xsize, template_img->ysize,
                                                            template_img->zsize);
    }


    #pragma omp parallel for
    for (size_t i = 0; i < atlas_set->n; i++) {
        iftImage *atlas     = iftReadImageByExt(atlas_set->files[i]->path);
        iftBoundingBox *mbb = iftMinLabelsBoundingBox(atlas, sosm_s->labels, NULL);
        
        for (size_t o = 0; o < labels->n; o++){
            mbb_cvs[o][i] = iftBoundingBoxCenterVoxel(mbb[o]);
            
            for (int z = mbb[o].begin.z; z <= mbb[o].end.z; z++)
                for (int y = mbb[o].begin.y; y <= mbb[o].end.y; y++)
                    for (int x = mbb[o].begin.x; x <= mbb[o].end.x; x++)
                            #pragma omp atomic
                            iftImgVal(sosm_s->obj_models[o]->prob_atlas, x, y, z) +=
                                    (iftImgVal(atlas, x, y, z) == labels->val[o]);
        }
        
        iftDestroyImage(&atlas);
        iftFree(mbb);
    }
    
    for (size_t o = 0; o < labels->n; o++) {
        iftSOSMSObjModel *obj_model = sosm_s->obj_models[o];
        
        iftBoundingBox prob_atlas_mbb = iftFMinBoundingBox(obj_model->prob_atlas, NULL);
        
        // begin voxel of the prob atlas on template image's coordinate space
        obj_model->begin.x = prob_atlas_mbb.begin.x;
        obj_model->begin.y = prob_atlas_mbb.begin.y;
        obj_model->begin.z = prob_atlas_mbb.begin.z;
        
        obj_model->template_shape = iftGetImageDomain(template_img);
        
        // for each object, the search region is the min bounding box that fits all center voxels from the
        // min bounding box of the atlases
        obj_model->search_region = iftMinBoundingBoxOfVoxels(mbb_cvs[o], atlas_set->n);
        iftFree(mbb_cvs[o]);
        
        iftFImage *cropped_prob_atlas = iftFExtractROI(obj_model->prob_atlas, prob_atlas_mbb);
        iftDestroyFImage(&obj_model->prob_atlas);
        
        // averaging
        #pragma omp parallel for
        for (int p = 0; p < cropped_prob_atlas->n; p++)
            cropped_prob_atlas->val[p] /= atlas_set->n;
        obj_model->prob_atlas = cropped_prob_atlas;
    }
    iftFree(mbb_cvs);
    
    return sosm_s;
}


iftImage *iftSegmentBySOSMS(const iftImage *img, const iftImage *grad_img_in, const iftSOSMS *sosm_s) {
    iftImage *grad_img  = NULL;
    if (grad_img_in == NULL)
        grad_img = iftComputeGradient(img, IFT_IMAGE_BASINS);
    else grad_img = iftCopyImage(grad_img_in);
    
    iftImage *seg_img       = iftCreateImageFromImage(img);
    iftFImage *prob_seg_img = iftCreateFImage(img->xsize, img->ysize, img->zsize);
    
    // for each object shape model
#pragma omp parallel for
    for (int o = 0; o < sosm_s->labels->n; o++) {
        printf("- Label: %d - Finding Seeds\n", sosm_s->labels->val[o]);
        iftSet *certain_obj_region = NULL;
        iftSet *forbidden          = NULL;
        iftLabeledSet *seeds = iftFindSOSMSSeeds(img, sosm_s->obj_models[o]);
        
        // The search region is a bounding box where the center of the min bounding box of the prob atlas can move
        // Since we first find the model seeds, we will move them instead of moving the prob. atlas
        // Then, the search region will be the maximum displacement from the center of the bounding box search region
        // to its extremities
        iftVector max_disp;
        max_disp.x = ceil((sosm_s->obj_models[o]->search_region.end.x - sosm_s->obj_models[o]->search_region.begin.x) / 2.0);
        max_disp.y = ceil((sosm_s->obj_models[o]->search_region.end.y - sosm_s->obj_models[o]->search_region.begin.y) / 2.0);
        max_disp.z = ceil((sosm_s->obj_models[o]->search_region.end.z - sosm_s->obj_models[o]->search_region.begin.z) / 2.0);
        
        printf("- Label: %d - Finding Best Seeds\n", sosm_s->labels->val[o]);
        iftVector best_disp = iftMSPSObjModel(grad_img, sosm_s->obj_models[o]->label, seeds, certain_obj_region, forbidden, 3, 4, max_disp);
        // iftVector best_disp = {0, 0, 0}; // to turn of object optimum location
        
        
        printf("- Label: %d - Translating Best Seeds\n", sosm_s->labels->val[o]);
        iftLabeledSet *final_seeds       = iftTranslateLabeledSet(seeds, grad_img, best_disp);
        // iftImage *seeds_img = iftCreateImageFromImage(img);
        // iftLabeledSetToImage(final_seeds, seeds_img, true);
        // iftWriteImageByExt(seeds_img, "out/seeds_obj_%d.nii.gz", sosm_s->labels->val[o]);
        
        iftSet *final_forbidden          = iftTranslateSet(forbidden, grad_img, best_disp);
        iftSet *final_certain_obj_region = iftTranslateSet(certain_obj_region, grad_img, best_disp);
        
        printf("- Label: %d - Delineating\n", sosm_s->labels->val[o]);
        iftImage *seg_obj_img = iftWatershed(grad_img, NULL, final_seeds, final_forbidden);
        iftSetToImage(final_certain_obj_region, seg_obj_img, sosm_s->obj_models[o]->label);
        
        printf("- Label: %d - Merging delineation\n", sosm_s->labels->val[o]);
        #pragma omp parallel for schedule(auto)
        for (int p = 0; p < seg_obj_img->n; p++) {
            if (seg_obj_img->val[p] != 0) {
                #pragma omp critical
                {
                    if (seg_img->val[p] == 0) {
                        seg_img->val[p] = seg_obj_img->val[p];
                    }
                        // tie break
                    else {
                        iftVoxel u = iftGetVoxelCoord(seg_obj_img, p);
                        iftVoxel v = iftVectorSub(u, sosm_s->obj_models[o]->begin); // voxel u on cropped prob atlas domain
                        float prob = (iftFValidVoxel(sosm_s->obj_models[o]->prob_atlas, v)) ? iftImgVoxelVal(sosm_s->obj_models[o]->prob_atlas, v) : 0.0f;
                        
                        if (prob > prob_seg_img->val[p]) {
                            seg_img->val[p] = seg_obj_img->val[p];
                            prob_seg_img->val[p] = prob;
                        }
                    }
                }
            }
        }
        
        // cleaning up
        iftDestroyLabeledSet(&seeds);
        iftDestroyLabeledSet(&final_seeds);
        iftDestroySet(&final_forbidden);
        iftDestroySet(&final_certain_obj_region);
        iftDestroyImage(&seg_obj_img);
    }
    iftDestroyImage(&grad_img);
    iftDestroyFImage(&prob_seg_img);
    
    return seg_img;
}


iftLabeledSet *iftFindSOSMSSeeds(const iftImage *test_img, const iftSOSMSObjModel *obj_model) {
    // it assumes that the test image is registered on the template image (or vice-versa), used to
    // build the object model
    iftFImage *prob_atlas = iftProbAtlasOnTemplateImageDomain(obj_model->prob_atlas, obj_model->template_shape, obj_model->begin);
    // iftWriteImageByExt(iftFImageToImage(prob_atlas, 255), "tmp/prob_atlas_%d.nii.gz", obj_model->label);
    
    // inner seeds = object's borders
    iftImage *obj_mask   = iftFThreshold(prob_atlas, 1, IFT_INFINITY_FLT, 1);
    iftSet *inner_seeds = iftObjectBorderSet(obj_mask, NULL);
    iftDestroyImage(&obj_mask);
    
    // outer seeds = border of the background
    iftImage *bg_mask   = iftFThreshold(prob_atlas, 0, 0, 1);
    iftSet *outer_seeds = iftObjectBorderSet(bg_mask, NULL);
    iftDestroyImage(&bg_mask);
    iftDestroyFImage(&prob_atlas);
    
    
    iftLabeledSet *all_seeds = NULL;
    iftInsertSetIntoLabeledSet(&inner_seeds, obj_model->label, &all_seeds);
    iftInsertSetIntoLabeledSet(&outer_seeds, 0, &all_seeds);
    
    return all_seeds;
}
