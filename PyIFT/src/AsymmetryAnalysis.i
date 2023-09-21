%include "iftDataSet.i"
%include "iftImage.i"
%include "iftMImage.i"



typedef enum {
    IFT_RIGHT_BRAIN_SIDE, IFT_LEFT_BRAIN_SIDE
} iftBrainSide;

%newobject iftBrainAsymMap;
%feature("autodoc", "2");
iftImage *iftBrainAsymMap(const iftImage *img, const iftImage *bias);

%newobject iftMeanBrainAsymMap;
%feature("autodoc", "2");
iftImage *iftMeanBrainAsymMap(const iftFileSet *img_set, bool add_stdev_asymmetries);

%newobject iftMeanBrainDiffMap;
%feature("autodoc", "2");
iftImage *iftMeanBrainDiffMap(const iftFileSet *img_set, const iftImage *template_img, bool add_stdev_asymmetries);

%newobject iftGridSamplingByBrainAsymmetries;
%feature("autodoc", "2");
iftIntArray * iftGridSamplingByBrainAsymmetries(const iftImage *asym_map, const iftImage *bin_mask,
                                                int min_samples_on_symmetries, float thres_factor);

%newobject iftBuildBrainHemisphereMImage;
%feature("autodoc", "2");
iftMImage *iftBuildBrainHemisphereMImage(const iftImage *img);

%newobject iftBuildBrainHemisphereMImageAsym;
%feature("autodoc", "2");
iftMImage *iftBuildBrainHemisphereMImageAsym(const iftImage *img, const iftImage *asym_map);

%newobject iftSymmISF;
%feature("autodoc", "2");
iftImage *iftSymmISF(const iftImage *img, const iftImage *bin_mask, float alpha, float beta, float thres_factor,
                     float min_dist_to_border, int n_seeds_on_symmetric_regions, const iftImage *normal_asymmap);

%newobject iftSymmOISF;
%feature("autodoc", "2");
iftImage *iftSymmOISF(const iftImage *img, const iftImage *bin_mask, int n_supervoxels, float alpha, float beta,
                      float gamma, const iftImage *normal_asymmap, float thres_factor);

%newobject iftExtractBrainSide;
%feature("autodoc", "2");
iftImage *iftExtractBrainSide(const iftImage *img, iftBrainSide side);

%newobject iftExtractSupervoxelHAAFeats;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelHAAFeats(const iftImage *test_img, const iftImage *test_sym_svoxels,
                                          const iftFileSet *train_set, int n_bins, const iftImage *normal_asym_map,
                                          int *n_svoxels_out);

%newobject iftExtractSupervoxelBICAsymmFeats;
%feature("autodoc", "2");
iftDataSet **iftExtractSupervoxelBICAsymmFeats(const iftImage *test_img, const iftImage *test_sym_svoxels,
                                               const iftFileSet *train_set, int n_bins_per_channel, const iftImage *normal_asym_map,
                                               int *n_svoxels_out);

%newobject iftBuildRefDataSupervoxelDataSets;
%feature("autodoc", "2");
void iftBuildRefDataSupervoxelDataSets(iftDataSet **Zarr, int n_supervoxels, const char *test_img_path,
                                       const char *test_supervoxels_path, const iftFileSet *train_set);

