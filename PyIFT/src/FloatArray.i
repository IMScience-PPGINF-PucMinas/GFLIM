


typedef struct {
    /** Number of elements. */
    long n;
    /** Array of float values. */
    float *val;
} iftFloatArray;

%extend iftFloatArray {

	PyObject* AsNumPy(void)
	{
	    iftFloatArray *a = ($self);
	    npy_intp dims[1] = {a->n};
	    PyArrayObject* result = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	    
	    float *data = array_data(result);
	    memcpy(data, a->val, a->n * sizeof (*data));
	    
	    return PyArray_Return(result);
	}
	
	void FromNumPy(PyObject* input) {
	    
	    PyArrayObject *ary = obj_to_array_no_conversion(input, NPY_FLOAT32);
	    
	    if (ary == NULL) return;
	    
	    if (!require_contiguous(ary))
	        SWIG_Error(0, "Input numpy array is not contiguous");
	    
	    iftFloatArray *a = ($self);
	    float *ptr = array_data(ary);
	    
	    int n_dims = array_numdims(ary);
	    if (n_dims != 1)
	        SWIG_Error(0, "Input is not a 1D array");
	    
	    int ary_size = array_size(ary, 0);
	    
	    if (ary_size != a->n)
	    {
	        a->n = ary_size;
	        a->val = iftRealloc(a->val, a->n * sizeof *a->val);
	        if (a->val == NULL)
	            iftError(MSG_MEMORY_ALLOC_ERROR, "FromNumpy");
	    }
	    memcpy(a->val, ptr, a->n * sizeof *a->val);
	}
	
	

	~iftFloatArray() {
		iftFloatArray* ptr = ($self);
		iftDestroyFloatArray(&ptr);
	}
};

%newobject iftFloatArrayAsString;
%feature("autodoc", "2");
char *iftFloatArrayAsString(const iftFloatArray *farr);

%inline %{

iftFloatArray *CreateFloatArrayFromNumPy(PyObject *input) {
    PyArrayObject *arr = obj_to_array_no_conversion(input, NPY_FLOAT32);
    if (!require_contiguous(arr))
        SWIG_Error(0, "Input numpy array is not contiguous");
    
    int n_dims = array_numdims(input);
    if (n_dims != 1) {
        SWIG_Error(0, "Input Array must be 1-dimensional");
        return NULL;
    }
    
    npy_intp *shape = array_dimensions(input);
    
    iftFloatArray *farr = iftCreateFloatArray(shape[0]);
    
    float *ptr = array_data(input);
    memcpy(farr->val, ptr, sizeof (*ptr) * farr->n);
    
    return farr;
}

%}

%newobject CreateFloatArrayFromNumPy;
%feature("autodoc", "2");
iftFloatArray *CreateFloatArrayFromNumPy(PyObject *input) ;

