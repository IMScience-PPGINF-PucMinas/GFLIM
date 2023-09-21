


typedef struct {
    /** Number of elements. */
    long n;
    /** Array of double values. */
    double *val;
} iftDblArray;

%extend iftDblArray {

	PyObject* AsNumPy(void)
	{
	    iftDblArray *a = ($self);
	    npy_intp dims[1] = {a->n};
	    PyArrayObject* result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
	    
	    double *data = array_data(result);
	    memcpy(data, a->val, a->n * sizeof (*data));
	    
	    return PyArray_Return(result);
	}
	
	void FromNumPy(PyObject* input) {
	    
	    PyArrayObject *ary = obj_to_array_no_conversion(input, NPY_DOUBLE);
	    
	    if (ary == NULL) return;
	    
	    if (!require_contiguous(ary))
	        SWIG_Error(0, "Input numpy array is not contiguous");
	    
	    iftDblArray *a = ($self);
	    double *ptr = array_data(ary);
	    
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
	
	

	~iftDblArray() {
		iftDblArray* ptr = ($self);
		iftDestroyDblArray(&ptr);
	}
};

%inline %{

iftDblArray *CreateDblArrayFromNumPy(PyObject *input) {
    PyArrayObject *arr = obj_to_array_no_conversion(input, NPY_DOUBLE);
    if (!require_contiguous(arr))
        SWIG_Error(0, "Input numpy array is not contiguous");
    
    int n_dims = array_numdims(input);
    if (n_dims != 1) {
        SWIG_Error(0, "Input Array must be 1-dimensional");
        return NULL;
    }
    
    npy_intp *shape = array_dimensions(input);
    
    iftDblArray *darr = iftCreateDblArray(shape[0]);
    
    double *ptr = array_data(input);
    memcpy(darr->val, ptr, sizeof (*ptr) * darr->n);
    
    return darr;
}

%}

%newobject CreateDblArrayFromNumPy;
%feature("autodoc", "2");
iftDblArray *CreateDblArrayFromNumPy(PyObject *input) ;

