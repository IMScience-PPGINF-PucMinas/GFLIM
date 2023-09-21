%include "iftMatrix.i"
%include "iftDataSet.i"
%include "iftMImage.i"
%include "iftAdjacency.i"
%include "iftIterativeOPF.i"

%newobject iftReadFLIMArch;
%feature("autodoc", "2");
  iftFLIMArch *iftReadFLIMArch(char *filename);

%newobject iftWriteFLIMArch;
%feature("autodoc", "2");
  void iftWriteFLIMArch(iftFLIMArch *arch, char *filename);

%feature("autodoc", "2");
  void iftFLIMLearnModel(char *orig_dir, char *markers_dir, char *param_dir, iftFLIMArch *arch);

%feature("autodoc", "2");
void iftFLIMLearnLayer(char *activ_dir, char *markers_dir, char *param_dir, int param_index, iftFLIMArch *arch, char *output_dir);

%feature("autodoc", "2");
  void iftFLIMExtractFeatures(char *orig_dir, char *image_list, iftFLIMArch *arch, char *param_dir, char *feat_dir, char *object_dir, int device);

%feature("autodoc", "2");
  void iftFLIMExtractFeaturesFromLayer(char *orig_dir, char *image_list, iftFLIMArch *arch, char *param_dir, int layer_index,
                                     char *feat_dir, char *object_dir, int device);

%feature("autodoc", "2");
iftAdjRel *iftFLIMAdjRelFromKernel(iftFLIMLayer layer, bool dim3D);

%feature("autodoc", "2");
//iftMImage     *iftFLIMAveragePooling(iftMImage *img, int width, int height, int depth, int stride);

%feature("autodoc", "2");
//iftMImage     *iftFLIMMaxPooling(iftMImage *mimg, int width, int height, int depth, int stride);

%feature("autodoc", "2");
iftMatrix *iftFLIMSelectKernelsManual(char *kernel_bank_path, char *selected_kernels_path);

