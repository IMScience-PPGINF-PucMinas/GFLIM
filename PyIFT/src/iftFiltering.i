%include "iftCommon.i"
%include "iftImage.i"
%include "iftFImage.i"
%include "iftMImage.i"
%include "iftAdjacency.i"
%include "iftSort.i"
%include "iftMatrix.i"
%include "iftKernel.i"

%newobject iftLinearFilter;
%feature("autodoc", "2");
iftImage  *iftLinearFilter(const iftImage *img, iftKernel *K);

%newobject iftLinearFilterInRegion;
%feature("autodoc", "2");
iftImage  *iftLinearFilterInRegion(iftImage *img, iftImage *mask, iftKernel *K);

%newobject iftFastLinearFilter;
%feature("autodoc", "2");
iftImage  *iftFastLinearFilter(iftImage *img, iftKernel *K, char crop);

%newobject iftFLinearFilter;
%feature("autodoc", "2");
iftFImage *iftFLinearFilter(iftFImage *img, iftKernel *K);

%newobject iftFLinearFilterInRegion;
%feature("autodoc", "2");
iftFImage  *iftFLinearFilterInRegion(iftFImage *img, iftImage *mask, iftKernel *K);

%newobject iftFastFLinearFilter;
%feature("autodoc", "2");
iftFImage *iftFastFLinearFilter(iftFImage *img, iftKernel *K, char crop);

%newobject iftMLinearFilter;
%feature("autodoc", "2");
iftMImage *iftMLinearFilter(iftMImage *img, iftMKernel *K);

%newobject iftMMLinearFilter;
%feature("autodoc", "2");
iftMImage *iftMMLinearFilter(iftMImage *img, iftMMKernel *k_bank);

%newobject iftMedianFilter;
%feature("autodoc", "2");
iftImage *iftMedianFilter(const iftImage *img, iftAdjRel *Ain);

%newobject iftMeanFilter;
%feature("autodoc", "2");
iftImage *iftMeanFilter(const iftImage *img, const iftAdjRel *Ain);

%newobject iftMMedianFilter;
%feature("autodoc", "2");
iftMImage *iftMMedianFilter(iftMImage *img, iftAdjRel *A);

%newobject iftModaFilter;
%feature("autodoc", "2");
iftImage  *iftModaFilter(iftImage *img, iftAdjRel *A);

%newobject iftSobelGradientMagnitude;
%feature("autodoc", "2");
iftImage  *iftSobelGradientMagnitude(iftImage *img);

%newobject iftImageToMatrix;
%feature("autodoc", "2");
iftMatrix *iftImageToMatrix(iftImage *img, iftFastAdjRel *F, char crop);

%newobject iftImageAdjToMatrix;
%feature("autodoc", "2");
iftMatrix *iftImageAdjToMatrix(iftImage *img, iftAdjRel *A);

%newobject iftMatrixToImage;
%feature("autodoc", "2");
iftImage  *iftMatrixToImage(iftMatrix *M, int xsize, int ysize, int zsize);

%newobject iftFImageToMatrix;
%feature("autodoc", "2");
iftMatrix *iftFImageToMatrix(iftFImage *img, iftFastAdjRel *F, char crop);

%newobject iftMatrixToFImage;
%feature("autodoc", "2");
iftFImage  *iftMatrixToFImage(iftMatrix *M, int xsize, int ysize, int zsize);

%newobject iftMImageToMatrix;
%feature("autodoc", "2");
iftMatrix *iftMImageToMatrix(iftMImage *img, iftFastAdjRel *F);

%newobject iftMatrixToMImage;
%feature("autodoc", "2");
iftMImage *iftMatrixToMImage(iftMatrix *M, int xsize, int ysize, int zsize, int nbands, char band_orientation);

%newobject iftMatrixToMImageArray;
%feature("autodoc", "2");
iftMImageArray *iftMatrixToMImageArray(iftMatrix *M, int nImgs, int xsize, int ysize, int zsize, int nbands, char band_orientation);

%newobject iftMImageToFeatureMatrix;
%feature("autodoc", "2");
iftMatrix *iftMImageToFeatureMatrix(iftMImage *mimg, iftAdjRel *A, float *fill_band_with);

%newobject iftMImageArrayToFeatureMatrix;
%feature("autodoc", "2");
iftMatrix *iftMImageArrayToFeatureMatrix(iftMImageArray *mimgArray, iftAdjRel *A);

%newobject iftKernelToMatrix;
%feature("autodoc", "2");
iftMatrix *iftKernelToMatrix(iftKernel *K);

%newobject iftMKernelToMatrix;
%feature("autodoc", "2");
iftMatrix *iftMKernelToMatrix(iftMKernel *K);

%newobject iftMMKernelToMatrix;
%feature("autodoc", "2");
iftMatrix *iftMMKernelToMatrix(iftMMKernel *k_bank);

%newobject iftSmoothImage;
%feature("autodoc", "2");
iftImage  *iftSmoothImage(iftImage *img, iftAdjRel *A, float sigma);

%newobject iftSmoothImageInRegion;
%feature("autodoc", "2");
iftImage  *iftSmoothImageInRegion(iftImage *img, iftImage *mask, iftAdjRel *A, float sigma);

%newobject iftFastBilateralFilter2D;
%feature("autodoc", "2");
iftImage *iftFastBilateralFilter2D(iftImage *img, int s_sigma, float r_sigma);

%newobject iftFastMBilateralFilter2D;
%feature("autodoc", "2");
iftMImage *iftFastMBilateralFilter2D(iftMImage *img, int s_sigma, float r_sigma);

%newobject iftNormalizeImage;
%feature("autodoc", "2");
iftImage *iftNormalizeImage(iftImage *img, iftAdjRel *A, int Imax);

%newobject iftAlphaPooling;
%feature("autodoc", "2");
iftImage *iftAlphaPooling(iftImage *img, iftAdjRel *A, int stride, float alpha);

%newobject iftRandomKernelBankAsMatrix;
%feature("autodoc", "2");
iftMatrix *iftRandomKernelBankAsMatrix(int size, int nbands, int nKernels);

%newobject iftMConvolution;
%feature("autodoc", "2");
iftMImage *iftMConvolution(iftMImage *mimg, iftMatrix *kernels, int stride);

%newobject iftMConvolutionArray;
%feature("autodoc", "2");
iftMImageArray *iftMConvolutionArray(iftMImageArray *mimgArray, iftMatrix *kernels, int stride);

%newobject iftMMaxPooling;
%feature("autodoc", "2");
iftMImage *iftMMaxPooling(iftMImage *mimg, float radius, int stride);

%newobject iftMMaxPoolingArray;
%feature("autodoc", "2");
iftMImageArray *iftMMaxPoolingArray(iftMImageArray *mimgArray, float radius, int stride);

%newobject iftMMaxPoolingRoi;
%feature("autodoc", "2");
iftMImage *iftMMaxPoolingRoi(iftMImage *mimg, iftRoiArray *roiArray, int stride);

%newobject iftMMaxPoolingRoiArray;
%feature("autodoc", "2");
iftMImageArray *iftMMaxPoolingRoiArray(iftMImageArray *mimgArray, iftRoiArray *roiArray, int stride);

%newobject iftMMinPooling;
%feature("autodoc", "2");
iftMImage *iftMMinPooling(iftMImage *mimg, float radius, int stride);

%newobject iftMMinPoolingArray;
%feature("autodoc", "2");
iftMImageArray *iftMMinPoolingArray(iftMImageArray *mimgArray, float radius, int stride);

%newobject iftMReLU;
%feature("autodoc", "2");
iftMImage *iftMReLU(iftMImage *mimg);

%newobject iftMReLUArray;
%feature("autodoc", "2");
iftMImageArray *iftMReLUArray(iftMImageArray *mimg);

