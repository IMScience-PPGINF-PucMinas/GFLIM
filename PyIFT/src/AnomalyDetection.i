%include "iftDataSet.i"
%include "iftImage.i"
%include "iftMImage.i"

%newobject iftRegErrorMagnitude;
%feature("autodoc", "2");
iftImage *iftRegErrorMagnitude(const iftImage *img, const iftImage *template_img, const iftImage *bias);

%newobject iftMeanRegErrorMagnitude;
%feature("autodoc", "2");
iftImage *iftMeanRegErrorMagnitude(const iftFileSet *img_set, const iftImage *template_img, bool add_stdev_error);

%newobject iftISFOnRegErrors;
%feature("autodoc", "2");
iftImage *iftISFOnRegErrors(const iftImage *img, const iftImage *reg_error_mag, const iftImage *template_img,
                            const iftImage *label_img, const iftDblArray *alphas, const iftDblArray *betas,
                            const iftDblArray *thres_factors, const iftDblArray *min_dist_to_borders,
                            const iftIntArray *n_seeds_on_correct_region);

%newobject iftISFOnRegErrorsFast;
%feature("autodoc", "2");
iftImage *iftISFOnRegErrorsFast(const iftImage *img, const iftImage *reg_error_mag, const iftImage *template_img,
                                const iftImage *label_img, const iftDblArray *alphas, const iftDblArray *betas,
                                const iftDblArray *thres_factors, const iftDblArray *min_dist_to_borders,
                                const iftIntArray *n_seeds_on_correct_region);

%newobject iftGridSamplingOnDomes;
%feature("autodoc", "2");
    iftIntArray *iftGridSamplingOnDomes(const iftImage *img, const iftImage *bin_mask, int n_samples_on_flat_region,
                                        float thres_factor, iftImage **domes_out);

%newobject iftExtractSupervoxelHistRegErrorsFeats;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelHistRegErrorsFeats(const iftImage *test_reg_error_mag, const iftImage *test_svoxels_img,
                                                    const iftFileSet *train_set, int n_bins, int *n_svoxels_out);

%newobject iftExtractSupervoxelHistRegErrorsFeatsDilation;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelHistRegErrorsFeatsDilation(const iftImage *test_img, const iftImage *test_svoxels_img,
                                                            const iftFileSet *train_set, const iftImage *template_img,
                                                            int n_bins, float radius, const iftImage *bias, int *n_svoxels_out);

%newobject iftExtractSupervoxelBandHistFeats;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelBandHistFeats(const iftMImage *test_filt_img, const iftImage *test_svoxels_img,
                                               const iftFileSet *train_set, int n_bins, int *n_svoxels_out);

%newobject iftWriteSupervoxelDataSets;
%feature("autodoc", "2");
void iftWriteSupervoxelDataSets(const iftDataSet **Zarr, int n_supervoxels, const char *out_dir);

%newobject iftComputeLinearAttenuationWeightsByEDT;
%feature("autodoc", "2");
iftFImage *iftComputeLinearAttenuationWeightsByEDT(const iftImage *label_img, float max_attenuation_factor);

%newobject iftComputeExponentialAttenuationWeightsByEDT;
%feature("autodoc", "2");
iftFImage *iftComputeExponentialAttenuationWeightsByEDT(const iftImage *label_img, float max_attenuation_factor, float exponent);

%newobject iftWeightedRegErrorMagnitude;
%feature("autodoc", "2");
iftImage *iftWeightedRegErrorMagnitude(const iftImage *img, const iftImage *template_img, const iftFImage *weights);

%newobject iftRemoveSVoxelsByVolAndMeanRegError;
%feature("autodoc", "2");
iftImage *iftRemoveSVoxelsByVolAndMeanRegError(const iftImage *svoxels_img, const iftImage *reg_error_mag,
                                               int min_vol, float min_mean_reg_error_on_svoxel);

%newobject iftISFOnAttentionMap;
%feature("autodoc", "2");
iftImage *iftISFOnAttentionMap(const iftImage *img, const iftImage *attention_map, const iftImage *target_img,
                               const iftImage *label_img, const iftDblArray *alphas, const iftDblArray *betas,
                               const iftDblArray *thres_factors, const iftDblArray *min_dist_to_borders,
                               const iftIntArray *n_seeds_on_correct_region);

%newobject iftExtractSupervoxelBICFeats;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelBICFeats(const iftImage *test_img, const iftImage *test_svoxels_img,
                                          const iftFileSet *train_set, int n_bins_per_channel,
                                          int *n_svoxels_out);

%newobject iftExtractSupervoxelAttentionMapFeats;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelAttentionMapFeats(const iftImage *test_attention_map, const iftImage *test_svoxels_img,
                                                   const iftFileSet *train_set, int n_bins, int *n_svoxels_out);

%newobject iftExtractSupervoxelLBPFeats;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelLBPFeats(const iftImage *test_img, const iftImage *test_svoxels_img,
                                          const iftFileSet *train_set, int n_bins, int *n_svoxels_out);

%newobject iftExtractSupervoxelVLBPFeats;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelVLBPFeats(const iftImage *test_img, const iftImage *test_svoxels_img,
                                           const iftFileSet *train_set, int n_bins, int *n_svoxels_out);

%newobject iftExtractSupervoxelTextureFeats;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelTextureFeats(const iftImage *test_img, const iftImage *test_attention_map,
                                              const iftImage *test_svoxels_img, const iftFileSet *train_set,
                                              const iftFileSet *train_attention_map_set,
                                              int n_bins, int *n_svoxels_out);

