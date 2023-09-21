%include "iftCommon.i"
%include "iftImage.i"
%include "iftMImage.i"
%include "iftAdjacency.i"
%include "iftFImage.i"
%include "iftRadiometric.i"
%include "iftMatrix.i"
%include "iftSort.i"
%include "iftKernel.i"
%include "iftDescriptors.i"



typedef enum {
    IFT_UNKNOWN         = 0x01,  // undefined status (default)
    IFT_TEST            = 0x02,  // testing sample
    IFT_TRAIN           = 0x04,  // training sample
    IFT_OUTLIER         = 0x08,  // a sample that does not fit to the training set
    IFT_ERROR           = 0x10,  // misclassified sample
    IFT_ARTIFICIAL      = 0x20,  // a sample created by some data augmentation technique
    IFT_SUPERVISED      = 0x40,  // it is a training sample whose true label has been confirmed by a specialist
    IFT_LABELPROPAGATED = 0x80,  // it is a training sample whose true label might not be correct, since it has been propagated by some procedure
    IFT_PROTOTYPE       = 0x100, // it is a representative of a group or class 
    IFT_CENTROID        = 0x200, // it is the centroid of a group
    IFT_ALL             = IFT_UNKNOWN | IFT_TEST | IFT_TRAIN | IFT_OUTLIER | IFT_ERROR |  IFT_ARTIFICIAL | IFT_SUPERVISED | IFT_LABELPROPAGATED | IFT_PROTOTYPE | IFT_CENTROID
} iftSampleStatus;



typedef enum {
    IFT_WEIGHT = 0,
    IFT_LABEL = 1,
    IFT_CLASS = 2,
    IFT_GROUP = 3,
    IFT_POINT = 4,
} iftFeatureLabel;



typedef enum {
    IFT_POSITIVE_CLASS = 1,
    IFT_NEGATIVE_CLASS = 2
} iftBinaryProblemLabel;



typedef struct {
    char **status;
    int niters;
    int nsamples;
} iftSampler;

%extend iftSampler {

	~iftSampler() {
		iftSampler* ptr = ($self);
		iftDestroySampler(&ptr);
	}
};



typedef struct {
    float *feat;    // feature values - it should be always a pointer to data in iftDataSet struct
    double *projection;
    int truelabel;   // 1,2,...,nclasses
    int label;   // 1,2,....,nclasses
    int id;      // identification which may be either a voxel or the sample's index in the file set
    int group; // 1,2,...,nClusters (used for OPF_unsup)
    // related image or the position of an image in the
    // directory with a related list of images.
    bool isSupervised;
    bool isLabelPropagated;
    int numberTimesChecked;

    float weight;  // probability density, path cost, etc
    uint status;   // Status of the sample: TEST, IFT_TRAIN, IFT_ERROR, IFT_OUTLIER.
} iftSample;

%extend iftSample {

	

};



typedef enum {
  IFT_REF_DATA_NULL, IFT_REF_DATA_IMAGE, IFT_REF_DATA_FIMAGE, IFT_REF_DATA_MIMAGE, IFT_REF_DATA_FILESET, IFT_REF_DATA_MMKERNEL,
  IFT_REF_DATA_CSV
} iftRefDataType;



typedef struct {
    iftSample *sample;   // sample
    int capacity; //max number of samples (pre allocates this space)
    int nfeats;   // number of features
    int nsamples; // number of samples
    int ngroups;  // number of clusters
    int nclasses; // number of ground truth classes
    int ntrainsamples; // number of training samples
    iftArcWeightFun iftArcWeight; // distance function
    int function_number; // This field is used to store the ID of the distance function set with iftSetDistanceFunction. It allows the
    // proper setting of the distance function when using iftWriteOPFDataset and iftReadOPFDataset
    float *alpha;   // coefficients used for arc weight computation
    iftRefDataType ref_data_type; // type of the Reference Data (mainly used to read and write datasets)
    void *ref_data; // it might be an image, for voxel datasets, a text file with image filenames, for image datasets, a region graph, for supervoxel datasets, etc.
    iftFeatSpaceParam fsp; // parameters of the feature scape transformation
    iftMatrix* data;
    iftDoubleMatrix* projection;
} iftDataSet;

%extend iftDataSet {

	~iftDataSet() {
		iftDataSet* ptr = ($self);
		iftDestroyDataSet(&ptr);
	}
	PyObject *GetData(void) {
	    iftDataSet *Z = ($self);
	
	    if (Z->data == NULL)
	        return Py_BuildValue(""); // returning a PyObject as None
	    else {
	        npy_intp dims[2] = {Z->nsamples, Z->nfeats};
	        PyArrayObject *ref = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, (float*) Z->data->val);
	        PyArrayObject *npy = PyArray_SimpleNew(2, dims, NPY_FLOAT32);
	        PyArray_CopyInto(npy, ref); // (dst, src)
	
	        return PyArray_Return(npy);
	    }
	}
	
	PyObject *GetProjection(void) {
	    iftDataSet *Z = ($self);
	
	    if (Z->projection == NULL)
	        return Py_BuildValue(""); // returning a PyObject as None
	    else {
	        npy_intp dims[2] = {Z->projection->nrows, Z->projection->ncols};
	        PyArrayObject *ref = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT64, (double*) Z->projection->val);
	        PyArrayObject *npy = PyArray_SimpleNew(2, dims, NPY_FLOAT64);
	        PyArray_CopyInto(npy, ref); // (dst, src)
	
	        return PyArray_Return(npy);
	
	    }
	}
	
	
	PyObject *GetTrueLabels(void) {
	    iftDataSet *Z = ($self);
	
	    npy_intp dims[1] = {Z->nsamples};
	    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_INT32);
	
	    int *ptr = (int *) array_data(npy);
	    for (int s = 0; s < Z->nsamples; s++)
	        ptr[s] = Z->sample[s].truelabel;
	
	    return PyArray_Return(npy);
	}
	
	PyObject *GetLabels(void) {
	    iftDataSet *Z = ($self);
	
	    npy_intp dims[1] = {Z->nsamples};
	    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_INT32);
	
	    int *ptr = (int *) array_data(npy);
	    for (int s = 0; s < Z->nsamples; s++)
	        ptr[s] = Z->sample[s].label;
	
	    return PyArray_Return(npy);
	}
	
	PyObject *GetIds(void) {
	    iftDataSet *Z = ($self);
	
	    npy_intp dims[1] = {Z->nsamples};
	    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_INT32);
	
	    int *ptr = (int *) array_data(npy);
	    for (int s = 0; s < Z->nsamples; s++)
	        ptr[s] = Z->sample[s].id;
	
	    return PyArray_Return(npy);
	}
	
	PyObject *GetGroups(void) {
	    iftDataSet *Z = ($self);
	
	    npy_intp dims[1] = {Z->nsamples};
	    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_INT32);
	
	    int *ptr = (int *) array_data(npy);
	    for (int s = 0; s < Z->nsamples; s++)
	        ptr[s] = Z->sample[s].group;
	
	    return PyArray_Return(npy);
	}
	
	PyObject *GetIsSupervised(void) {
	    iftDataSet *Z = ($self);
	
	    npy_intp dims[1] = {Z->nsamples};
	    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_BOOL);
	
	    bool *ptr = (bool *) array_data(npy);
	    for (int s = 0; s < Z->nsamples; s++)
	        ptr[s] = Z->sample[s].isSupervised;
	
	    return PyArray_Return(npy);
	}
	
	PyObject *GetIsLabelPropagated(void) {
	    iftDataSet *Z = ($self);
	
	    npy_intp dims[1] = {Z->nsamples};
	    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_BOOL);
	
	    bool *ptr = (bool *) array_data(npy);
	    for (int s = 0; s < Z->nsamples; s++)
	        ptr[s] = Z->sample[s].isLabelPropagated;
	
	    return PyArray_Return(npy);
	}
	
	PyObject *GetWeights(void) {
	    iftDataSet *Z = ($self);
	
	    npy_intp dims[1] = {Z->nsamples};
	    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	
	    float *ptr = (float*) array_data(npy);
	    for (int s = 0; s < Z->nsamples; s++)
	        ptr[s] = Z->sample[s].weight;
	
	    return PyArray_Return(npy);
	}
	
	PyObject *GetStatus(void) {
	    iftDataSet *Z = ($self);
	
	    npy_intp dims[1] = {Z->nsamples};
	    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_UINT32);
	
	    unsigned int *ptr = (unsigned int *) array_data(npy);
	    for (int s = 0; s < Z->nsamples; s++)
	        ptr[s] = Z->sample[s].status;
	
	    return PyArray_Return(npy);
	}
	
	PyObject *GetAlphas(void) {
	    iftDataSet *Z = ($self);
	
	    if (Z->alpha == NULL)
	        return Py_BuildValue(""); // returning a PyObject as None
	    else {
	        npy_intp dims[1] = {Z->nfeats};
	        PyArrayObject *ref = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT32, (float*) Z->alpha);
	        PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	        PyArray_CopyInto(npy, ref); // (dst, src)
	
	        return PyArray_Return(npy);
	    }
	}
	
	PyObject *GetRefData(void) {
	    iftDataSet *Z = ($self);
	
	    if (Z->ref_data == NULL) {
	        return PyList_New(0); // return an empty List
	    }
	    else {
	        if (Z->ref_data_type == IFT_REF_DATA_FILESET) {
	            PyObject *list = PyList_New(Z->nsamples);
	
	            iftFileSet *ref_data = (iftFileSet *) Z->ref_data;
	
	            for (int i = 0; i < ref_data->n; i++)
	                PyList_SetItem(list, i, PyString_FromString(ref_data->files[i]->path));
	            return list;
	        }
	        if (Z->ref_data_type == IFT_REF_DATA_CSV) {
	            iftCSV *ref_data = (iftCSV *) Z->ref_data;
	
	            PyObject *csv = PyList_New(ref_data->nrows);
	
	            for (int r = 0; r < ref_data->nrows; r++) {
	                PyObject *row = PyList_New(ref_data->ncols);
	                
	                for (int c = 0; c < ref_data->ncols; c++)
	                    PyList_SetItem(row, c, PyString_FromString(ref_data->data[r][c]));
	                PyList_SetItem(csv, r, row);
	            }
	            
	            return csv;
	        }
	        else return PyList_New(0); // return an empty List
	    }
	}
	
	
	void SetRefData(PyObject *list) {
	
	    /// Example of list: ref_data = [['00000001_000001.png'], ['00000001_000002.png'], ['00000002_000001.png'], ['00000002_000002.png']]
	
	    iftDataSet *Z = ($self);
	
	//    if (Z->ref_data != NULL)
	//        SWIG_Error(0, "Existing Reference Data is not NULL");
	 
	    if (!PyList_CheckExact(list)){
	        SWIG_Error(0, "Input \"list\" must be a list.");
	        return;
	    }
	
	    long list_size = PyList_Size(list);
	
	    if (list_size != Z->nsamples) {
	        char error[512];
	        sprintf(error, "Size of the reference data (list of paths) is != from the number of samples: %ld != %ld", list_size, Z->nsamples);
	        SWIG_Error(0, error);
	        return;
	    }
	    
	    
	    PyObject *item0 = PyList_GetItem(list, 0);
	
	    if (!PyList_CheckExact(item0)){
	        SWIG_Error(0, "Input \"item0\" must be a list.");
	        return;
	    }
	
	    long n_cols = PyList_Size(item0);
	
	    if (n_cols == 1) {
	        Z->ref_data_type = IFT_REF_DATA_FILESET;
	        
	        iftFileSet *ref_data = iftCreateFileSet(Z->nsamples);
	        
	        for (int i = 0; i < Z->nsamples; i++) {
	            PyObject *row = PyList_GetItem(list, i);
	
	            if (!PyList_CheckExact(row)){
	               SWIG_Error(0, "Input \"row\" must be a list.");
	               return;
	            }
	
	            char *path = PyUnicodeToString(PyList_GetItem(row, 0));
	
	            ref_data->files[i] = iftCreateFile(path);
	            Z->sample[i].id = i;
	        }
	        Z->ref_data = ref_data;
	    }
	    else {
	        Z->ref_data_type = IFT_REF_DATA_CSV;
	
	        iftCSV *ref_data = iftCreateCSV(Z->nsamples, n_cols);
	
	        for (int r = 0; r < Z->nsamples; r++) {
	            PyObject *row = PyList_GetItem(list, r);
	            
	            for (int c = 0; c < n_cols; c++) {
	                strcpy(ref_data->data[r][c], PyUnicodeToString(PyList_GetItem(row, c)));
	                Z->sample[r].id = r;
	            }
	        }
	        
	        Z->ref_data = ref_data;
	    }
	}
	
	
	void SetTrueLabels(PyObject *in_npy) {
	    char error[512];
	    iftDataSet *Z = ($self);
	
	    PyArrayObject *npy = NULL;
	    if (in_npy != Py_None) {
	        if (PyArray_TYPE(in_npy) != NPY_INT32) {
	            SWIG_Error(0, "Numpy array must be INT32");
	            return;
	        }
	
	        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);
	
	        int n_dims = array_numdims(npy);
	        if (n_dims != 1) {
	            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        npy_intp *shape = array_dimensions(npy); // (nsamples,)
	        if (shape[0] != Z->nsamples) {
	            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        int *ptr = (int *) array_data(npy);
	        for (int s = 0; s < Z->nsamples; s++)
	            Z->sample[s].truelabel = ptr[s];
	
	        Z->nclasses = iftCountNumberOfClassesDataSet(Z);
	    }
	
	}
	
	
	void SetId(PyObject *in_npy) {
	    char error[512];
	    iftDataSet *Z = ($self);
	
	    PyArrayObject *npy = NULL;
	    if (in_npy != Py_None) {
	        if (PyArray_TYPE(in_npy) != NPY_INT32) {
	            SWIG_Error(0, "Numpy array must be INT32");
	            return;
	        }
	
	        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);
	
	        int n_dims = array_numdims(npy);
	        if (n_dims != 1) {
	            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        npy_intp *shape = array_dimensions(npy); // (nsamples,)
	        if (shape[0] != Z->nsamples) {
	            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        int *ptr = (int *) array_data(npy);
	        for (int s = 0; s < Z->nsamples; s++)
	            Z->sample[s].id = ptr[s];
	   
	    }
	}
	
	
	void SetLabels(PyObject *in_npy) {
	    char error[512];
	    iftDataSet *Z = ($self);
	
	    PyArrayObject *npy = NULL;
	    if (in_npy != Py_None) {
	        if (PyArray_TYPE(in_npy) != NPY_INT32) {
	            SWIG_Error(0, "Numpy array must be INT32");
	            return;
	        }
	
	        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);
	
	        int n_dims = array_numdims(npy);
	        if (n_dims != 1) {
	            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        npy_intp *shape = array_dimensions(npy); // (nsamples,)
	        if (shape[0] != Z->nsamples) {
	            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        int *ptr = (int *) array_data(npy);
	        for (int s = 0; s < Z->nsamples; s++)
	            Z->sample[s].label = ptr[s];
	
		int tmp = iftCountNumberOfClassesDataSet(Z);
		if (tmp > Z->nclasses)
		    Z->nclasses = tmp;     
	    }
	}
	
	void SetWeight(PyObject *input)
	{
	    char error[512];
	    PyArrayObject *ary = obj_to_array_no_conversion(input, NPY_FLOAT);
	
	    if (ary == NULL) return;
	
	    if (!require_contiguous(ary))
	        SWIG_Error(0, "Input numpy array is not contiguous");
	
	    iftDataSet *Z = ($self);
	    float* ptr = array_data(ary);
	
	
	    int n_dims = array_numdims(ary);
	    if (n_dims != 1) {
	            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
	            SWIG_Error(0, error);
	            return;
	        }
	
	    npy_intp *shape = array_dimensions(ary); // (nsamples,)
	    if (shape[0] != Z->nsamples) {
	            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
	            SWIG_Error(0, error);
	            return;
	    }
	
	    for (int s = 0; s < Z->nsamples; s++)
	        Z->sample[s].weight = ptr[s];
	
	}
	
	void SetGroups(PyObject *in_npy) {
	    char error[512];
	    iftDataSet *Z = ($self);
	
	    PyArrayObject *npy = NULL;
	    if (in_npy != Py_None) {
	        if (PyArray_TYPE(in_npy) != NPY_INT32) {
	            SWIG_Error(0, "Numpy array must be INT32");
	            return;
	        }
	
	        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);
	
	        int n_dims = array_numdims(npy);
	        if (n_dims != 1) {
	            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        npy_intp *shape = array_dimensions(npy); // (nsamples,)
	        if (shape[0] != Z->nsamples) {
	            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        int *ptr = (int *) array_data(npy);
	        for (int s = 0; s < Z->nsamples; s++)
	            Z->sample[s].group = ptr[s];
	
	        Z->ngroups = iftCountNumberOfGroupsDataSet(Z);
	    }
	}
	
	void SetGroupsFromImage(const iftImage *img, bool increment)
	{
	    iftDataSet *Z = ($self);
	    if (Z->nsamples != img->n)
	        return;
	
	    for (int i = 0; i < img->n; i++)
	        Z->sample[i].group = img->val[i] + increment;
	
	    Z->ngroups = iftCountNumberOfGroupsDataSet(Z);
	}
	
	
	void SetNClasses(int n){
	  
	    iftDataSet *Z = ($self);
	    Z->nclasses = n;
	}
	
	
	void SetData(PyObject *input)
	{
	    PyArrayObject *ary = obj_to_array_no_conversion(input, NPY_FLOAT);
	
	    if (ary == NULL) return;
	
	    if (!require_contiguous(ary))
	        SWIG_Error(0, "Input numpy array is not contiguous");
	
	    iftDataSet *Z = ($self);
	    float* ptr = array_data(ary);
	
	    int n_dims = array_numdims(ary);
	    if (n_dims != 2)
	        SWIG_Error(0, "Input must be a 2-D numpy array");
	
	    int ary_xsize = array_size(ary, 1);
	    int ary_ysize = array_size(ary, 0);
	
	    if (Z->nsamples != ary_ysize || Z->nfeats != ary_xsize)
	        SWIG_Error(0, "Input must be the same size as the previously set data");
	
	    for (int i = 0; i < Z->nsamples; i++) {
	        for (int j = 0; j < Z->nfeats; j++) {
	            Z->sample[i].feat[j] = ptr[i * Z->nfeats + j];
	        }
	    }
	}
	
	void SetStatus(PyObject *in_npy) {
	    char error[512];
	    iftDataSet *Z = ($self);
	
	    PyArrayObject *npy = NULL;
	    if (in_npy != Py_None) {
	        if (PyArray_TYPE(in_npy) != NPY_INT32) {
	            SWIG_Error(0, "Numpy array must be INT32");
	            return;
	        }
	
	        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);
	
	        int n_dims = array_numdims(npy);
	        if (n_dims != 1) {
	            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        npy_intp *shape = array_dimensions(npy); // (nsamples,)
	        if (shape[0] != Z->nsamples) {
	            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
	            SWIG_Error(0, error);
	            return;
	        }
	
	        int *ptr = (int *) array_data(npy);
	        for (int s = 0; s < Z->nsamples; s++)
	            Z->sample[s].status = ptr[s];
	   
	    }
	}
	
	void SetNTrainSamples(int n){
	  
	    iftDataSet *Z = ($self);
	    Z->ntrainsamples = n;
	}
	
	
	

};

%newobject iftKFold;
%feature("autodoc", "2");
iftSampler *iftKFold(int nsamples, int k);

%newobject iftInverseKFold;
%feature("autodoc", "2");
iftSampler *iftInverseKFold(int nsamples, int k);

%newobject iftRandomSubsampling;
%feature("autodoc", "2");
iftSampler *iftRandomSubsampling(int nsamples, int n, int ntrain);

%newobject iftBalancedRandomSubsampling;
%feature("autodoc", "2");
iftSampler *iftBalancedRandomSubsampling(const int *labels, int n_samples, int n_iters, int n_train);

%feature("autodoc", "2");
void iftDataSetSampling(iftDataSet *Z, const iftSampler *sampler, int iteration);

%newobject iftGetDataSetTrueLabels;
%feature("autodoc", "2");
iftIntArray *iftGetDataSetTrueLabels(const iftDataSet *dataset);

%newobject iftCreateDataSet;
%feature("autodoc", "2");
iftDataSet *iftCreateDataSet(int nsamples, int nfeats);

%newobject iftCopyDataSet;
%feature("autodoc", "2");
iftDataSet *iftCopyDataSet(const iftDataSet *Z, bool copy_feats);

%newobject iftDataSetObjectMap;
%feature("autodoc", "2");
iftImage *iftDataSetObjectMap(const iftDataSet *Z, const iftImage *comp, int max_val, int label);

%feature("autodoc", "2");
void iftSetRefData(iftDataSet *Z, const void *ref_data, iftRefDataType ref_data_type);

%newobject iftExtractGroup;
%feature("autodoc", "2");
iftDataSet *iftExtractGroup(iftDataSet *Z, int group);

%newobject iftExtractSamples;
%feature("autodoc", "2");
iftDataSet *iftExtractSamples(const iftDataSet *Z, iftSampleStatus status);

%newobject iftDataSetToLabelImage;
%feature("autodoc", "2");
iftImage *iftDataSetToLabelImage(const iftDataSet *Z, const iftImage *comp, bool decrement_labels, iftFeatureLabel label_type);

%newobject iftDataSetClusterInformationToLabelImage;
%feature("autodoc", "2");
iftImage *iftDataSetClusterInformationToLabelImage(const iftDataSet *Z, const iftImage *comp, bool decrement_groups);

%feature("autodoc", "2");
void iftSetStatus(iftDataSet *Z, iftSampleStatus status);

%feature("autodoc", "2");
void iftAddStatus(iftDataSet *Z, iftSampleStatus status);

%feature("autodoc", "2");
void iftRemoveStatus(iftDataSet *Z, iftSampleStatus status);

%newobject iftMImageToDataSet;
%feature("autodoc", "2");
iftDataSet *iftMImageToDataSet(const iftMImage *mimg, const iftImage *label_img, float coord_scale);

%newobject iftSplitDataSetAt;
%feature("autodoc", "2");
void iftSplitDataSetAt(const iftDataSet *Z, int sample_idx, iftDataSet **Z1, iftDataSet **Z2);

%feature("autodoc", "2");
void iftWriteDataSet(const iftDataSet *Z, const char *pathname, ...);

%newobject iftReadDataSet;
%feature("autodoc", "2");
iftDataSet *iftReadDataSet(const char *pathname, ...);

%feature("autodoc", "2");
void iftLabelDataSetFromSeeds(iftDataSet *Z, iftLabeledSet *seeds, iftImage* region);

%feature("autodoc", "2");
void iftSelectDataSetClassesByCluster(iftDataSet *Z, float threshold);

%newobject iftBoundingBoxArrayToDataSet;
%feature("autodoc", "2");
iftDataSet *iftBoundingBoxArrayToDataSet(const iftImage *img, const iftBoundingBoxArray *bb_ary, int n_feats);

%newobject iftPatchesFromSuperpixels;
%feature("autodoc", "2");
iftDataSet *iftPatchesFromSuperpixels(char *input_dir, char *output, int nsuperpixels, int patch_size);

%newobject iftDataSetFromAllSeeds;
%feature("autodoc", "2");
  iftDataSet *iftDataSetFromAllSeeds(char *markers_dir, char *mimages_dir, iftAdjRel *A);

%inline %{

iftDataSet *CreateDataSetFromNumPy(PyObject *in_ary, PyObject *in_truelabel) {
    
    PyArrayObject* ary = obj_to_array_no_conversion(in_ary, NPY_FLOAT32);
    PyArrayObject* truelabel = NULL;
    if (in_truelabel != Py_None) {
        truelabel = obj_to_array_no_conversion(in_truelabel, NPY_INT32);
    }
    
    if (ary == NULL) return NULL;
    
    int n_dims = array_numdims(ary);
    if (n_dims != 2) {
        SWIG_Error(0, "Feature Matrix (data) must be 2-dimensional");
        return NULL;
    }

    npy_intp *shape = array_dimensions(ary); // (nsamples, nfeats)

    iftDataSet *Z = iftCreateDataSet(shape[0], shape[1]);

    if (truelabel != NULL) {
        int n_dims_tl = array_numdims(truelabel);
        if (n_dims_tl != 1) {
            SWIG_Error(0, "True Label Array (truelabel) must be 1-dimensional");
            return NULL;
        }

        npy_intp *shape_tl = array_dimensions(truelabel);

        if (shape[0] != shape_tl[0]) {
            char msg[512];
            sprintf(msg, "Different dimensions for Feature Matrix rows and True Label size: %ld != %ld", shape[0], shape_tl[0]);
            SWIG_Error(0, msg);
            return NULL;
        }

        // copying the true labels from numpy array to iftDataSet 
        int *truelabel_ptr = (int*) array_data(truelabel);
        
        #pragma omp parallel for
        for (int s = 0; s < Z->nsamples; s++) 
            Z->sample[s].truelabel = truelabel_ptr[s];

	Z->nclasses = iftCountNumberOfClassesDataSet(Z);

    }

    

    // copying the feature matrix from numpy array to iftDataSet 
    float *data_ptr = (float*) array_data(ary);
    memcpy(Z->data->val, data_ptr, sizeof (*data_ptr) * Z->nsamples * Z->nfeats);

    return Z;
}

float (*iftWeightedL2NormPtr(void))(float *f1, float *f2, float *alpha, int n) {
    return iftDistance1;
}


%}

%newobject CreateDataSetFromNumPy;
%feature("autodoc", "2");
iftDataSet *CreateDataSetFromNumPy(PyObject *in_ary, PyObject *in_truelabel) ;

float (*iftWeightedL2NormPtr(void))(float *f1, float *f2, float *alpha, int n) ;

