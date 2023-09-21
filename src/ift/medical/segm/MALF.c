#include "ift/medical/segm/MALF.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/IntMatrix.h"
#include "ift/core/io/Dir.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/tools/OS.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftSegmentation.h"
#include "iftSimilarity.h"


iftFileSet *iftAtlasSelectionByNMI(const iftImage *test_img, const iftFileSet *train_img_set, int n_atlases,
                                   const iftImage *test_mask, const iftFileSet *train_mask_set,
                                   iftIntArray **selected_idxs_out) {
    if (n_atlases > train_img_set->n)
        iftError("Number of Selected Atlases is > than number of Training Images: %d > %d",
                 "iftAtlasSelectionByNMI", n_atlases, train_img_set->n);
    else if (n_atlases > train_img_set->n) {
        iftWarning("Number of Selected Atlases == number of Training Images\n" \
                   "A copy of the training image will be returned",
                   "iftAtlasSelectionByNMI");
        return iftCopyFileSet(train_img_set);
    }
    
    if (train_mask_set && train_mask_set->n != train_img_set->n)
        iftError("Number of train masks is different from the number of train images: %ld != %ld",
                 "iftAtlasSelectionByNMI", train_mask_set->n, train_img_set->n);
    
    iftDblArray *scores = iftCreateDblArray(train_img_set->n);

    iftImage *final_test_img = (test_mask) ? iftMask(test_img, test_mask) : iftCopyImage(test_img);
    
    #pragma omp parallel for
    for (int f = 0; f < train_img_set->n; f++) {
        iftImage *train_img = iftReadImageByExt(train_img_set->files[f]->path);
        iftImage *final_train_img = (train_mask_set) ? iftReadImageByExt(train_mask_set->files[f]->path) : train_img;
        
        scores->val[f] = iftNormalizedMutualInformation(final_train_img, final_test_img);
        
        if (train_img != final_train_img)
            iftDestroyImage(&final_train_img);
        iftDestroyImage(&train_img);
    }
    
    iftIntArray *indexes = iftIntRange(0, train_img_set->n, 1);
    iftDQuickSort(scores->val, indexes->val, 0, scores->n, IFT_DECREASING);
    
    iftFileSet *chosen_set = iftCreateFileSet(n_atlases);
    iftIntArray *selected_idxs = iftCreateIntArray(n_atlases);

   #pragma omp parallel for
    for (int i = 0; i < n_atlases; i++) {
        selected_idxs->val[i] = indexes->val[i];
        chosen_set->files[i] = iftCreateFile(train_img_set->files[selected_idxs->val[i]]->path);
    }
    iftDestroyIntArray(&indexes);
    iftDestroyDblArray(&scores);
    
    if (selected_idxs_out)
        *selected_idxs_out = selected_idxs;
    else iftDestroyIntArray(&selected_idxs);
    
    iftDestroyImage(&final_test_img);
    
    return chosen_set;
}


iftImage *iftSegmentByClassicalMALF(const iftImage *test_img, const iftFileSet *train_label_set) {
    int max_label = IFT_INFINITY_INT_NEG;

    #pragma omp parallel for schedule(auto)
    for (int f = 0; f < train_label_set->n; f++) {
        iftImage *label_img = iftReadImageByExt(train_label_set->files[f]->path);
        int max_img_label   = iftMaximumValue(label_img);
        
        if (max_label < max_img_label)
            max_label = max_img_label;
        iftDestroyImage(&label_img);
    }
    
    iftIntMatrix *label_occurs = iftCreateIntMatrix(max_label+1, test_img->n);

    #pragma omp parallel for schedule(auto)
    for (int f = 0; f < train_label_set->n; f++) {
        iftImage *label_img = iftReadImageByExt(train_label_set->files[f]->path);
        
        // counts label occurrences on each voxel
        for (int p = 0; p < label_img->n; p++) {
            int label = label_img->val[p];
            iftMatrixElem(label_occurs, label, p)++;
        }
        
        iftDestroyImage(&label_img);
    }
    
    iftImage *final_seg_img = iftMajorityVoting(test_img, label_occurs);
    iftDestroyIntMatrix(&label_occurs);
    
    return final_seg_img;
}



iftImage *iftSegmentByMALFSTAPLE(const iftImage *test_img, const iftFileSet *train_atlas_set) {
    char *tmpdir = iftMakeTempDir("tmpdir_staple_", NULL, NULL);
    char *out_img_path = iftJoinPathnames(2, tmpdir, "segmentation.nii.gz");
    
    char *args = iftConcatStrings(2, "-o ", out_img_path);
    
    for (int f = 0; f < train_atlas_set->n; f++) {
        char *aux = args;
        args = iftConcatStrings(3, args, " ", train_atlas_set->files[f]->path);
        iftFree(aux);
    }
    iftRunProgram("crlSTAPLE", args);
    iftFree(args);
    
    args = iftConcatStrings(3, out_img_path, " ", out_img_path);
    iftRunProgram("crlIndexOfMaxComponent", args);
    
    iftImage *seg_img = iftReadImageByExt(out_img_path);
    
    iftFree(args);
    iftRemoveDir(tmpdir);
    iftFree(tmpdir);
    iftFree(out_img_path);
    
    return seg_img;
}

