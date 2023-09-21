%include "iftCommon.i"
%include "iftImage.i"
%include "iftFImage.i"
%include "iftAdjacency.i"
%include "iftMatrix.i"



typedef struct {
    iftMatrix      *data; /* matrix of pixel features */
    float          **val; /* fast access for each pixel feature (data row) ex. mimg->val[pixel][band]*/
    int            xsize,ysize,zsize; /* image dimensions of each band */
    float          dx,dy,dz;  /* voxel size */
    int            *tby, *tbz; /* LUT to speed up index to/from coordinate conversions */
    unsigned long  n,m; /* number of voxels and number of bands */
} iftMImage;

%extend iftMImage {

	iftImage* __getitem__(int band){
	    iftMImage* mimg = ($self);
	    iftImage* img;
	
	    img = iftMImageToImage(mimg, 255, band);
	    return img;
	}
	
	void __setitem__(int band, iftImage* img){
	    iftMImage* mimg = ($self);
	
	    if(band >= mimg->m){
	        SWIG_Error(12, "Index out of MImage bands");
	    }
	
	    for(int i = 0; i < img->n; i++){
	        mimg->val[i][band] = (float) img->val[i];
	    }
	}
	
	/* bugado
	void AddImage(iftImage* img){
	    iftMImage* mimg = ($self);
	    mimg->m++;
	    iftRealloc(mimg->band, sizeof(iftBand) * mimg->m);
	
	    for(int i = 0; i < img->n; i++){
	        mimg->band[mimg->m - 1].val[i] = (float) img->val[i];
	    }
	}
	*/
	
	PyObject* AsNumPy(void) {
	
	    iftMImage *mimg = ($self);
	
	    int n_dims = 4;
	    npy_intp *dims = calloc(n_dims, sizeof (*dims));
	    dims[0] = mimg->zsize;
	    dims[1] = mimg->ysize;
	    dims[2] = mimg->xsize;
	    dims[3] = mimg->m;
	
	    int n_bands = mimg->m;
	
	    PyArrayObject *result = PyArray_SimpleNew(n_dims, dims, NPY_FLOAT);
	    float *data = array_data(result);
	
	    for (int p = 0; p < mimg->n; p++) {
	        int start_voxel_idx = p * n_bands;
	        for (int b = 0; b < n_bands; b++) {
	            data[start_voxel_idx + b] = mimg->val[p][b];
	        }
	    }
	
	    free(dims);
	    return PyArray_Return(result);
	}
	
	void WriteAsNumPy(const char *filename) {
	    iftMImage *mimg = ($self);
	    iftWriteMImageAsNumPy(mimg, filename);
	}
	

	~iftMImage() {
		iftMImage* ptr = ($self);
		iftDestroyMImage(&ptr);
	}
};



typedef struct {
  iftMImage **val;
  int n;
} iftMImageArray;

%extend iftMImageArray {

	iftMImage * __getitem__(int index)
	{
	    iftMImageArray *arr = ($self);
	
	    if (index < 0 || index >= arr->n) {
	        SWIG_Error(0, "Index out of bonds.");
	        return NULL;
	    }
	
	    return iftCopyMImage(arr->val[index]);
	}

	~iftMImageArray() {
		iftMImageArray* ptr = ($self);
		iftDestroyMImageArray(&ptr);
	}
};

%newobject iftMGetVoxelCoord;
%feature("autodoc", "2");
iftVoxel    iftMGetVoxelCoord(const iftMImage *img, int p);

%newobject iftCreateMImage;
%feature("autodoc", "2");
iftMImage  *iftCreateMImage(int xsize,int ysize,int zsize, int nbands);

%newobject iftMValidVoxel;
%feature("autodoc", "2");
char iftMValidVoxel(const iftMImage *img, iftVoxel v);

%newobject iftImageToMImage;
%feature("autodoc", "2");
iftMImage  *iftImageToMImage(const iftImage *img, char color_space);

%newobject iftMImageToImage;
%feature("autodoc", "2");
  iftImage   *iftMImageToImage(const iftMImage *img, int Imax, int band);

%feature("autodoc", "2");
static inline bool iftIs3DMImage(const iftMImage *img) ;

%newobject iftReadMImage;
%feature("autodoc", "2");
iftMImage   *iftReadMImage(const char *filename);

%newobject iftReadMImageFromNumPy;
%feature("autodoc", "2");
iftMImage *iftReadMImageFromNumPy(const char *npy_path, ...);

%feature("autodoc", "2");
void  	     iftWriteMImage(iftMImage *mimg, const char *filename);

%feature("autodoc", "2");
void iftWriteMImageAsNumPy(const iftMImage *mimg, const char *npy_path, ...);

%feature("autodoc", "2");
void        iftWriteMImageBands(iftMImage *mimg, char *base_filename);

%newobject iftMImageBasins;
%feature("autodoc", "2");
iftImage   *iftMImageBasins(const iftMImage *img, iftAdjRel *A);

%newobject iftMImageGradient;
%feature("autodoc", "2");
iftImage   *iftMImageGradient(const iftMImage *img, iftAdjRel *A, int Imax);

%newobject iftGridSampling;
%feature("autodoc", "2");
iftImage *iftGridSampling(iftMImage *img, iftImage *mask, int nsamples);

%newobject iftAltMixedSampling;
%feature("autodoc", "2");
iftImage *iftAltMixedSampling(iftMImage *img, iftImage *mask, int nsamples);

%newobject iftSelectNSamplesFromMask;
%feature("autodoc", "2");
iftImage *iftSelectNSamplesFromMask(iftMImage *img, iftImage *mask1, int nsamples);

%newobject iftExtractImageFeatures;
%feature("autodoc", "2");
iftMImage *iftExtractImageFeatures(const iftImage *img, const iftImage *label_img, const iftAdjRel *A, bool voxel_coord);

%newobject iftKernelTransformMImage;
%feature("autodoc", "2");
iftMImage *iftKernelTransformMImage(const iftMImage *input, const iftDoubleMatrix *train_data,
                                    const iftDoubleMatrix *omega, iftKernelFunction *K,
                                    double param1, double param2, double scale);

%newobject iftColorQuantization;
%feature("autodoc", "2");
iftMImage *iftColorQuantization(const iftMImage *input, float lower_bound, float upper_bound, int bins_per_band);

%newobject iftMImageCentralize;
%feature("autodoc", "2");
iftMImage *iftMImageCentralize(const iftMImage *mimg);

%newobject iftMGaussianFilter;
%feature("autodoc", "2");
iftMImage *iftMGaussianFilter(const iftMImage *mimg, float radius, float stddev);

%newobject iftBinarySetsNCA;
%feature("autodoc", "2");
iftMImage *iftBinarySetsNCA(const iftMImage *mimg, const iftSet *obj, const iftSet *bkg,
                            int d_out, int n_iters, double learn_rate);

%feature("autodoc", "2");
void iftRescaleMImage(iftMImage *mimg, bool single_max);

%newobject iftStackGrayImages;
%feature("autodoc", "2");
iftMImage *iftStackGrayImages(int n, ...);

%newobject iftMGetVoxelIndex_pyift;
%feature("autodoc", "2");
int iftMGetVoxelIndex_pyift(iftMImage *s, iftVoxel v);

%inline %{

iftMImage *CreateMImageFromNumPy(PyObject *input) {
    int is_new = 0;
    PyArrayObject *ary = obj_to_array_allow_conversion(input, NPY_FLOAT, &is_new);

    if (ary == NULL)
        return NULL;

    if (!require_contiguous(ary)) {
        SWIG_Error(0, "Input numpy array is not contiguous");
        return NULL;
    }

    int n_dims = array_numdims(ary);

    // (ysize, xsyze, nbands) or
    // (zsize, ysize, xsyze, nbands)
    int xsize, ysize, zsize, n_bands; 

    if (n_dims == 3) {
        zsize = 1;
        ysize = array_size(ary, 0);
        xsize = array_size(ary, 1);
        n_bands = array_size(ary, 2);
    }
    else if (n_dims == 4) {
        zsize = array_size(ary, 0);
        ysize = array_size(ary, 1);
        xsize = array_size(ary, 2);
        n_bands = array_size(ary, 3);
    }
    else {
        char error[256];
        sprintf(error, "Image Numpy Array must have 3 or 4 dimensions, %d found.", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    iftMImage *mimg = iftCreateMImage(xsize, ysize, zsize, n_bands);

    float *ptr = array_data(ary);

    for (int p = 0; p < mimg->n; p++) {
        int start_voxel_idx = p * n_bands;

        for (int b = 0; b < n_bands; b++) {
            mimg->val[p][b] = ptr[start_voxel_idx + b];
        }
    }
    
    return mimg;
}


%}

%newobject CreateMImageFromNumPy;
%feature("autodoc", "2");
iftMImage *CreateMImageFromNumPy(PyObject *input) ;

