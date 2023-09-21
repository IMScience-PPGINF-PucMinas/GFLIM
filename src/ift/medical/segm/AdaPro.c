#include "ift/medical/segm/AdaPro.h"

#include "ift/core/dtypes/Dict.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Json.h"
#include "ift/core/io/NumPy.h"
#include "ift/core/tools/String.h"
#include "ift/medical/registration/Elastix.h"
#include "iftDataSet.h"
#include "iftMathMorph.h"
#include "iftRegistration.h"
#include "iftSegmentation.h"



/********************** PRIVATE FUNCTIONS *************************/
iftObjModel *_iftCreateObjModel(int label) {
    iftObjModel *obj_model    = (iftObjModel*) calloc(1, sizeof(iftObjModel));
    obj_model->label          = label;
    obj_model->erosion_radius = obj_model->dilation_radius = 0.0;
    obj_model->begin.x        = obj_model->begin.y = obj_model->begin.z = 0;
    
    return obj_model;
}


void _iftDestroyObjModel(iftObjModel **obj_model) {
    iftObjModel *aux = *obj_model;
    
    if (aux != NULL) {
        iftDestroyFImage(&aux->prob_atlas);
        iftFree(aux);
        *obj_model = NULL;
    }
}


iftAdaPro *_iftCreateAdaPro(const iftIntArray *labels, const iftImage *template_img) {
    iftAdaPro *adapro = ((iftAdaPro*) iftAlloc(1, sizeof(iftAdaPro)));
    adapro->template_img = iftCopyImage(template_img);
    
    adapro->labels = iftCreateIntArray(labels->n);
    iftCopyIntArray(adapro->labels->val, labels->val, labels->n);
    
    adapro->obj_models = (iftObjModel**) calloc(adapro->labels->n, sizeof(iftObjModel*));
    for (size_t o = 0; o < labels->n; o++)
        adapro->obj_models[o] = _iftCreateObjModel(labels->val[o]);
    
    return adapro;
}


iftImage *_iftAdaProTextureClassification(const iftImage *img, const iftAdaPro *adapro) {
    iftImage *rough_mask_dilated = iftDilateAdaProRoughSegmentation(adapro);
    iftBoundingBox bb_brain = iftMinBoundingBox(rough_mask_dilated, NULL);
    
    iftImage *mask = iftCreateImageFromImage(img);
    iftAdjRel *A = (iftIs3DImage(img)) ? iftSpheric(5.0) : iftCircular(5.0);
    int size = 10;

    #pragma omp parallel for
    for (int z = bb_brain.begin.z; z <= bb_brain.end.z; z += size)
        for (int y = bb_brain.begin.y; y <= bb_brain.end.y; y += size)
            for (int x = bb_brain.begin.x; x <= bb_brain.end.x; x += size) {
                iftBoundingBox bb = {.begin = {x, y, z}};
                bb.end.x = bb.begin.x + (size-1);
                bb.end.y = bb.begin.y + (size-1);
                bb.end.z = bb.begin.z + (size-1);
                
                // fit to the bb of the brain
                if (bb.end.x > bb_brain.end.x)
                    bb.end.x = bb_brain.end.x;
                if (bb.end.y > bb_brain.end.y)
                    bb.end.y = bb_brain.end.y;
                if (bb.end.z > bb_brain.end.z)
                    bb.end.z = bb_brain.end.z;
                
                
                iftDataSet *Zroi = iftImageROIToDataSetByAdjacencyFeats(img, rough_mask_dilated, bb, A);
                iftSetStatus(Zroi, IFT_TEST);
                iftSVMLinearClassifyOVA(adapro->svm, Zroi, IFT_TEST, NULL);

                #pragma omp parallel for
                for (int s = 0; s < Zroi->nsamples; s++)
                    mask->val[Zroi->sample[s].id] = Zroi->sample[s].label - 1; // background is 1 in Zroi
                
                iftDestroyDataSet(&Zroi);
            }
    iftDestroyImage(&rough_mask_dilated);
    iftDestroyAdjRel(&A);
    
    return mask;
}


/********************** PUBLIC FUNCTIONS *************************/
iftFImage *iftProbAtlasOnTemplateImageDomain(const iftFImage *prob_atlas, iftImageDomain template_shape,
                                             iftVoxel begin) {
    iftFImage *prob_atlas_on_test_space = iftCreateFImage(template_shape.xsize,
                                                          template_shape.ysize,
                                                          template_shape.zsize);
    iftFInsertROI(prob_atlas, prob_atlas_on_test_space, begin);
    
    return prob_atlas_on_test_space;
}


void iftDestroyAdaPro(iftAdaPro **adapro) {
    iftAdaPro *aux = *adapro;
    
    if (aux != NULL) {
        for (int o = 0; o < aux->labels->n; o++)
            _iftDestroyObjModel(&aux->obj_models[o]);
        
        iftDestroyIntArray(&aux->labels);
        iftDestroyImage(&aux->template_img);
        iftDestroyImage(&aux->rough_mask);
        iftDestroyImage(&aux->train_voxels_img);
        iftDestroySVM(aux->svm);
        aux->svm = NULL;
        iftDestroyStrArray(&aux->elastix_files);
        
        iftFree(aux);
        *adapro = NULL;
    }
}


iftAdaPro *iftReadAdaPro(const char *format, ...) {
    va_list args;
    char path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(path, format, args);
    va_end(args);
    
    char *tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);
    iftUnzipFile(path, tmp_dir);
    
    char *info_path = iftJoinPathnames(2, tmp_dir, "info.json");
    iftDict *info = iftReadJson(info_path);
    
    iftIntArray *labels = iftGetIntArrayFromDict("labels", info);
    iftImage *template_img = iftReadImageByExt(iftJoinPathnames(2, tmp_dir, "template.nii.gz"));
    iftAdaPro *adapro = _iftCreateAdaPro(labels, template_img);
    
    char *train_voxels_img_path = iftJoinPathnames(2, tmp_dir, "training_voxels.nii.gz");
    if (iftFileExists(train_voxels_img_path))
        adapro->train_voxels_img = iftReadImageByExt(train_voxels_img_path);
    iftFree(train_voxels_img_path);
    
    char *rough_mask_path = iftJoinPathnames(2, tmp_dir, "rough_mask.nii.gz");
    if (iftFileExists(rough_mask_path))
        adapro->rough_mask = iftReadImageByExt(rough_mask_path);
    iftFree(rough_mask_path);
    
    char *svm_path = iftJoinPathnames(2, tmp_dir, "svm.zip");
    if (iftFileExists(svm_path))
        adapro->svm = iftReadSVM(svm_path);
    iftFree(svm_path);
    
    char *elastix_dir = iftJoinPathnames(2, tmp_dir, "elastix_registration_files");
    if (iftDirExists(elastix_dir)) {
        iftFileSet *elastix_fset = iftLoadFileSetFromDirBySuffix(elastix_dir, ".txt", 0);
        iftSetElastixFilesToAdaPro(adapro, elastix_fset);
        iftDestroyFileSet(&elastix_fset);
    }
    iftFree(elastix_dir);
    
    char *obj_models_dir = iftJoinPathnames(2, tmp_dir, "obj_models");
    for (int o = 0; o < adapro->labels->n; o++) {
        iftObjModel *obj_model = _iftCreateObjModel(adapro->labels->val[o]);
        
        char obj_str[32];
        sprintf(obj_str, "%d", adapro->labels->val[o]);
        char *obj_dir = iftJoinPathnames(2, obj_models_dir, obj_str);
        
        char *prob_atlas_path = iftJoinPathnames(3, obj_models_dir, obj_str, "prob_atlas.npy");
        obj_model->prob_atlas = iftReadFImage(prob_atlas_path);
        iftFree(prob_atlas_path);
        iftFree(obj_dir);
        
        iftDict *obj_model_dict = iftGetDictFromDict(iftConcatStrings(2, "object-models:", obj_str), info);
        obj_model->label = iftGetLongValFromDict("label", obj_model_dict);
        
        iftIntArray *template_shape = iftGetIntArrayFromDict("template-shape", obj_model_dict);
        obj_model->template_shape.xsize = template_shape->val[0];
        obj_model->template_shape.ysize = template_shape->val[1];
        obj_model->template_shape.zsize = template_shape->val[2];
        
        iftIntArray *begin = iftGetIntArrayFromDict("begin", obj_model_dict);
        obj_model->begin = iftIntArrayToVoxel(begin);
        
        obj_model->erosion_radius = iftGetDblValFromDict("erosion-radius", obj_model_dict);
        obj_model->dilation_radius = iftGetDblValFromDict("dilation-radius", obj_model_dict);
        
        adapro->obj_models[o] = obj_model;
    }
    
    iftDestroyIntArray(&labels);
    iftDestroyImage(&template_img);
    
    iftRemoveDir(tmp_dir);
    iftFree(obj_models_dir);
    iftFree(tmp_dir);
    
    return adapro;
}


void iftWriteAdaPro(const iftAdaPro *adapro, const char *format, ...) {
    va_list args;
    char path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(path, format, args);
    va_end(args);
    
    if (!iftEndsWith(path, ".zip"))
        iftError("Invalid extension for AdaPro: %s... Try *.zip", "iftWriteAdaPro", path);
    
    char *parent_dir = iftParentDir(path);
    if (!iftDirExists(parent_dir))
        iftMakeDir(parent_dir);
    iftFree(parent_dir);
    
    char *tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);
    
    // writing template
    char *template_path = iftJoinPathnames(2, tmp_dir, "template.nii.gz");
    iftWriteImageByExt(adapro->template_img, template_path);
    iftFree(template_path);
    
    // writing rough segmentation mask
    char *rough_mask_path = iftJoinPathnames(2, tmp_dir, "rough_mask.nii.gz");
    iftWriteImageByExt(adapro->rough_mask, rough_mask_path);
    iftFree(rough_mask_path);
    
    // writing markers
    if (adapro->train_voxels_img != NULL) {
        char *train_voxels_img_path = iftJoinPathnames(2, tmp_dir, "training_voxels.nii.gz");
        iftWriteImageByExt(adapro->train_voxels_img, train_voxels_img_path);
        iftFree(train_voxels_img_path);
    }
    
    // writing linear SVM texture classifier
    if (adapro->svm != NULL) {
        char *svm_path = iftJoinPathnames(2, tmp_dir, "svm.zip");
        iftWriteSVM(adapro->svm, svm_path);
        iftFree(svm_path);
    }
    
    // writing elastix registration files
    if (adapro->elastix_files) {
        char *elastix_dir = iftJoinPathnames(2, tmp_dir, "elastix_registration_files");
        iftMakeDir(elastix_dir);
        
        for (long i = 0; i < adapro->elastix_files->n; i++) {
            char out_filename[512];
            sprintf(out_filename, "ElastixRegistrationFile.%02ld.txt", i);
            char *out_path = iftJoinPathnames(2, elastix_dir, out_filename);
            
            iftWriteStringToFile(adapro->elastix_files->val[i], out_path);
            
            iftFree(out_path);
        }
        iftFree(elastix_dir);
    }
    
    // writing a Json with information
    char *info_path = iftJoinPathnames(2, tmp_dir, "info.json");
    iftDict *info   = iftCreateDict();
    
    iftInsertIntoDict("labels", adapro->labels, info);
    
    char *obj_models_dir = iftJoinPathnames(2, tmp_dir, "obj_models");
    iftMakeDir(obj_models_dir);
    
    for (int o = 0; o < adapro->labels->n; o++) {
        iftObjModel *obj_model = adapro->obj_models[o];
        
        char obj_str[32];
        sprintf(obj_str, "%d", adapro->labels->val[o]);
        char *obj_dir = iftJoinPathnames(2, obj_models_dir, obj_str);
        iftMakeDir(obj_dir);
        
        char *prob_atlas_path = iftJoinPathnames(2, obj_dir, "prob_atlas.npy");
        iftWriteFImage(obj_model->prob_atlas, prob_atlas_path);
        iftFree(prob_atlas_path);
        iftFree(obj_dir);
        
        iftDict *obj_model_dict = iftCreateDict();
        iftInsertIntoDict("label", obj_model->label, obj_model_dict);
        
        iftIntArray *template_shape = iftCreateIntArray(3);
        template_shape->val[0] = obj_model->template_shape.xsize;
        template_shape->val[1] = obj_model->template_shape.ysize;
        template_shape->val[2] = obj_model->template_shape.zsize;
        iftInsertIntoDict("template-shape", template_shape, obj_model_dict);
        
        iftIntArray *begin = iftVoxelToIntArray(obj_model->begin);
        iftInsertIntoDict("begin", begin, obj_model_dict);
        
        iftInsertIntoDict("erosion-radius", obj_model->erosion_radius, obj_model_dict);
        iftInsertIntoDict("dilation-radius", obj_model->dilation_radius, obj_model_dict);
        
        iftInsertIntoDict(iftConcatStrings(2, "object-models:", obj_str), obj_model_dict, info);
    }
    iftFree(obj_models_dir);
    
    iftWriteJson(info, info_path);
    iftZipDirContent(tmp_dir, path);
    
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
    iftFree(info_path);
    iftDestroyDict(&info);
}


iftAdaPro *iftTrainAdaPro(const iftFileSet *atlas_set, const iftImage *template_img, const iftImage *train_voxels_img,
                          const iftIntArray *labels, const iftFloatArray *e_radius_arr,
                          const iftFloatArray *d_radius_arr, double C, const iftFileSet *elastix_fset) {
    iftAdaPro *adapro = _iftCreateAdaPro(labels, template_img);
    
    adapro->train_voxels_img = iftCopyImage(train_voxels_img);
    
    for (long o = 0; o < labels->n; o++)
        adapro->obj_models[o]->prob_atlas = iftCreateFImage(template_img->xsize, template_img->ysize,
                                                            template_img->zsize);
    
    // summing all label images to build the prob. atlases
    #pragma omp parallel for
    for (size_t i = 0; i < atlas_set->n; i++) {
        iftImage *atlas     = iftReadImageByExt(atlas_set->files[i]->path);
        iftBoundingBox *mbb = iftMinLabelsBoundingBox(atlas, adapro->labels, NULL);
        
        for (size_t o = 0; o < labels->n; o++)
            for (int z = mbb[o].begin.z; z <= mbb[o].end.z; z++)
                for (int y = mbb[o].begin.y; y <= mbb[o].end.y; y++)
                    for (int x = mbb[o].begin.x; x <= mbb[o].end.x; x++)
                            #pragma omp atomic
                            iftImgVal(adapro->obj_models[o]->prob_atlas, x, y, z) +=
                                    (iftImgVal(atlas, x, y, z) == labels->val[o]);
        
        iftDestroyImage(&atlas);
        iftFree(mbb);
    }
    
    
    // creating a rough brain mask
    adapro->rough_mask = iftCreateImageFromImage(template_img);
    for (long o = 0; o < labels->n; o++)
            #pragma omp parallel for
            for (int p = 0; p < adapro->rough_mask->n; p++)
                if (!iftAlmostZero(adapro->obj_models[o]->prob_atlas->val[p]))
                    adapro->rough_mask->val[p] = 1;
    
    iftImage *aux = adapro->rough_mask;
    adapro->rough_mask = iftCloseBin(adapro->rough_mask, 8);
    iftDestroyImage(&aux);
    
    // building object models
    for (long o = 0; o < labels->n; o++) {
        iftObjModel *obj_model = adapro->obj_models[o];
        
        obj_model->erosion_radius = e_radius_arr->val[o];
        obj_model->dilation_radius = d_radius_arr->val[o];
        
        iftBoundingBox prob_atlas_mbb = iftFMinBoundingBox(obj_model->prob_atlas, NULL);
        
        // begin voxel of the prob atlas on template image's coordinate space
        obj_model->begin.x = prob_atlas_mbb.begin.x;
        obj_model->begin.y = prob_atlas_mbb.begin.y;
        obj_model->begin.z = prob_atlas_mbb.begin.z;
        
        obj_model->template_shape = iftGetImageDomain(adapro->template_img);
        
        iftFImage *cropped_prob_atlas = iftFExtractROI(obj_model->prob_atlas, prob_atlas_mbb);
        iftDestroyFImage(&obj_model->prob_atlas);
        
        // averaging
        #pragma omp parallel for
        for (int p = 0; p < cropped_prob_atlas->n; p++)
            cropped_prob_atlas->val[p] /= atlas_set->n;
        obj_model->prob_atlas = cropped_prob_atlas;
    }
    
    // training texture classifier
    if (train_voxels_img) {
        // training the linear SVM texture classifier
        iftAdjRel *A = (iftIs3DImage(template_img)) ? iftSpheric(5.0) : iftCircular(5.0);
        iftDataSet *Ztrain = iftImageToDataSetByAdjacencyFeats(template_img, train_voxels_img, A);
        iftSetStatus(Ztrain, IFT_TRAIN);
        adapro->svm = iftCreateSVM(IFT_LINEAR, IFT_OVA, C, 1.0 / Ztrain->nfeats); // create SVM
        iftSVMTrainOVA(adapro->svm, Ztrain, Ztrain);
        iftDestroyAdjRel(&A);
    }
    
    
    // saving the Elastix Registration Files
    if (elastix_fset)
        iftSetElastixFilesToAdaPro(adapro, elastix_fset);
    
    return adapro;
}


iftImage *iftSegmentByAdaPro(const iftImage *img, const iftAdaPro *adapro, const char *aux_basename) {
    iftImage *seg_img       = iftCreateImageFromImage(img);
    iftFImage *prob_seg_img = iftCreateFImage(img->xsize, img->ysize, img->zsize); // used to treat ties
    
    // texture classification by linear SVM
    iftImage *mask = NULL;
    if (adapro->svm != NULL) {
        puts("  - Texture Classification");
        mask = _iftAdaProTextureClassification(img, adapro);
        
        if (aux_basename)
            iftWriteImageByExt(mask, "%s_classification_mask.nii.gz", aux_basename);
    }
    
    iftImage *grad_img = iftImageBasins(img, NULL);
    if (aux_basename)
        iftWriteImageByExt(grad_img, "%s_gradient.nii.gz", aux_basename);

    #pragma omp parallel for
    for (int o = 0; o < adapro->labels->n; o++) {
        int label = adapro->labels->val[o];
        printf("  -- Label: %d - Finding Seeds\n", label);
        iftSet *forbidden = NULL;
        iftLabeledSet *seeds = iftFindObjModelSeeds(adapro->obj_models[o], mask, &forbidden);
        
        if (aux_basename) {
            iftImage *seeds_img = iftCreateImageFromImage(img);
            iftLabeledSetToImage(seeds, seeds_img, true);
            iftWriteImageByExt(seeds_img, "%s_seeds_%d.nii.gz", aux_basename, label);
            iftSetToImage(forbidden, seeds_img, 5);
            iftLabeledSetToImage(seeds, seeds_img, true);
            iftWriteImageByExt(seeds_img, "%s_seeds_and_forbidden_%d.nii.gz", aux_basename, label);
            iftDestroyImage(&seeds_img);
        }
        
        printf("  -- Label: %d - Delineating\n", label);
        iftImage *seg_obj_img = iftWatershed(grad_img, NULL, seeds, forbidden);

        #pragma omp critical
        {
            printf("  -- Label: %d - Merging delineation\n", adapro->labels->val[o]);
            for (int p = 0; p < seg_obj_img->n; p++) {
                if (seg_obj_img->val[p] != 0) {
                    iftVoxel u = iftGetVoxelCoord(seg_obj_img, p);
                    iftVoxel v = iftVectorSub(u, adapro->obj_models[o]->begin); // voxel u on cropped prob atlas domain
                    float prob = (iftFValidVoxel(adapro->obj_models[o]->prob_atlas, v)) ? iftImgVoxelVal(adapro->obj_models[o]->prob_atlas, v) : 0.0f;
                    
                    if ((seg_img->val[p] == 0) || (prob > prob_seg_img->val[p])) {
                        seg_img->val[p] = seg_obj_img->val[p];
                        prob_seg_img->val[p] = prob;
                    }
                }
            }
        } // #pragma critical
        
        // cleaning up
        iftDestroySet(&forbidden);
        iftDestroyLabeledSet(&seeds);
        iftDestroyImage(&seg_obj_img);
    }
    
    puts("  - Smoothing Segmentation");
    iftImage *final_seg_img = iftSmoothRegionsByDiffusion(seg_img, img, 0.0, 1);
    // iftImage *final_seg_img = iftSmoothRegionsByDiffusion(seg_img, img, 0.5, 5);
    
    iftDestroyFImage(&prob_seg_img);
    iftDestroyImage(&mask);
    iftDestroyImage(&grad_img);
    iftDestroyImage(&seg_img);
    
    return final_seg_img;
}


iftFileSet *iftWriteAdaProElastixParamFiles(const iftAdaPro *adapro, const char *basename) {
    if (adapro->elastix_files == NULL)
        iftWarning("AdaPro model does not have Elastix Parameter Files. NULL is returned", "iftWriteAdaProElastixParamFiles");
    
    char *final_basename = (basename) ? iftCopyString(basename) : iftMakeTempPathname(NULL, NULL, NULL);
    
    if (final_basename) {
        char *parent_dir = iftParentDir(final_basename);

        if (!iftDirExists(parent_dir))
            iftMakeDir(parent_dir);
        iftFree(parent_dir);
    }
    
    iftFileSet *elastix_files = iftCreateFileSet(adapro->elastix_files->n);
    for (long i = 0; i < elastix_files->n; i++) {
        char *out_path = iftCopyString("%s.%ld.txt", final_basename, i);
        elastix_files->files[i] = iftCreateFile(out_path);
        iftWriteStringToFile(adapro->elastix_files->val[i], elastix_files->files[i]->path);
        iftFree(out_path);
    }
    
    iftFree(final_basename);
    
    return elastix_files;
}


void iftRegisterAdaProOnTestImageByElastix(iftAdaPro *adapro, const iftImage *test_img, const iftFileSet *elastix_files,
                                           iftFileSet **def_fields_out) {
    iftImage *template_old = adapro->template_img;
    iftFileSet *def_fields = NULL;
    iftImage *template_reg = iftRegisterImageByElastix(template_old, test_img, NULL, NULL,
                                                       elastix_files, NULL, &def_fields,NULL);
    
    adapro->template_img = iftLinearStretch(template_reg, iftMinimumValue(template_reg),
                                            iftMaximumValue(template_reg), 0, 4095);
    iftDestroyImage(&template_reg);
    
    // mapping the brain mask
    iftImage *aux = adapro->rough_mask;
    adapro->rough_mask = iftTransformImageByTransformix(adapro->rough_mask, def_fields->files[def_fields->n - 1]->path);
    iftDestroyImage(&aux);
    
    iftImageDomain test_img_dom = iftGetImageDomain(test_img);

    #pragma omp parallel for
    for (int o = 0; o < adapro->labels->n; o++) {
        iftImage *prob_atlas_img         = iftFImageToImage(adapro->obj_models[o]->prob_atlas, 4095);
        iftImage *prob_atlas_on_template = iftCreateImageFromImage(template_old);
        iftInsertROI(prob_atlas_img, prob_atlas_on_template, adapro->obj_models[o]->begin);
        iftDestroyImage(&prob_atlas_img);
        
        iftImage *mapped_prob_atlas_img = iftTransformImageByTransformix(prob_atlas_on_template,
                                                                         def_fields->files[def_fields->n - 1]->path);
        iftDestroyImage(&prob_atlas_on_template);
        
        iftBoundingBox mbb               = iftMinBoundingBox(mapped_prob_atlas_img, NULL);
        iftImage *cropped_prob_atlas_img = iftExtractROI(mapped_prob_atlas_img, mbb);
        iftDestroyImage(&mapped_prob_atlas_img);
        
        iftDestroyFImage(&adapro->obj_models[o]->prob_atlas);
        adapro->obj_models[o]->prob_atlas = iftImageToFImageMaxVal(cropped_prob_atlas_img, 1.0);
        iftDestroyImage(&cropped_prob_atlas_img);
        
        adapro->obj_models[o]->begin          = mbb.begin;
        adapro->obj_models[o]->template_shape = iftGetImageDomain(test_img);
        
        adapro->obj_models[o]->template_shape = test_img_dom;
    }
    iftDestroyImage(&template_old);
    
    if (def_fields_out)
        *def_fields_out = def_fields;
    else {
        iftRemoveFileSet(def_fields);
        iftDestroyFileSet(&def_fields);
    }
}



iftLabeledSet *iftFindObjModelSeeds(const iftObjModel *obj_model, const iftImage *clf_mask, iftSet **forbidden) {
    // it assumes that the test image is registered on the template image (or vice-versa), used to
    // build the object model
    iftFImage *prob_atlas = iftProbAtlasOnTemplateImageDomain(obj_model->prob_atlas, obj_model->template_shape,
                                                              obj_model->begin);
    // iftWriteImageByExt(iftFImageToImage(prob_atlas, 255), "tmp/prob_atlas_%d.nii.gz", obj_model->label);
    
    // outer seeds = border of the background
    iftImage *prob_atlas_bin = iftFThreshold(prob_atlas, 0.000001, 1, 1);
    iftSet *outer_seeds = NULL;
    iftImage *dilation  = iftDilateBin(prob_atlas_bin, &outer_seeds, obj_model->dilation_radius);
    iftImage *rough_mask_dilated = iftCloseBin(dilation, 8);
    iftDestroyImage(&dilation);
    iftDestroyImage(&prob_atlas_bin);
    
    // inner seeds = object's borders
    iftImage *obj_mask = iftFThreshold(prob_atlas, 1, IFT_INFINITY_FLT, 1);
    iftSet *inner_seeds = NULL;
    iftImage *erosion  = iftErodeBin(obj_mask, &inner_seeds, obj_model->erosion_radius);
    
    // filtering the voxels classified as background
    if (clf_mask != NULL) {
        // removing the object seeds classified as bg
        iftSet *filt_inner_seeds = NULL;
        while (inner_seeds != NULL) {
            int p = iftRemoveSet(&inner_seeds);
            if (clf_mask->val[p])
                iftInsertSet(&filt_inner_seeds, p);
        }
        inner_seeds = filt_inner_seeds;
        
        // forbidding all object voxels in rough segmentation classied as bg
        for (int p = 0; p < rough_mask_dilated->n; p++) {
            if (rough_mask_dilated->val[p] && !clf_mask->val[p])
                iftInsertSet(forbidden, p);
        }
    }
    
    iftDestroyImage(&obj_mask);
    iftDestroyImage(&erosion);
    iftDestroyFImage(&prob_atlas);
    iftDestroyImage(&rough_mask_dilated);
    
    iftLabeledSet *all_seeds = NULL;
    iftInsertSetIntoLabeledSet(&inner_seeds, obj_model->label, &all_seeds);
    iftInsertSetIntoLabeledSet(&outer_seeds, 0, &all_seeds);
    
    return all_seeds;
}


void iftSetElastixFilesToAdaPro(iftAdaPro *adapro, const iftFileSet *elastix_fset) {
    iftDestroyStrArray(&adapro->elastix_files);
    
    adapro->elastix_files = iftCreateStrArray(elastix_fset->n);
    
    for (long i = 0; i < elastix_fset->n; i++) {
        iftFree(adapro->elastix_files->val[i]); // its creation function allocates a string here
        adapro->elastix_files->val[i] = iftReadFileAsString(elastix_fset->files[i]->path);
    }
}


iftImage *iftDilateAdaProRoughSegmentation(const iftAdaPro *adapro) {
    float max_dilation = IFT_INFINITY_FLT_NEG;
    for (int o = 0; o < adapro->labels->n; o++)
        if (adapro->obj_models[o]->dilation_radius > max_dilation)
            max_dilation = adapro->obj_models[o]->dilation_radius;
    
    iftSet *S = NULL;
    iftImage *rough_mask_dilated = iftDilateBin(adapro->rough_mask, &S, max_dilation);
    iftDestroySet(&S);
    
    return rough_mask_dilated;
}



