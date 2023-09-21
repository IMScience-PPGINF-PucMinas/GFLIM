//! swig(newobject, stable)
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

//! swig(stable)
float (*iftWeightedL2NormPtr(void))(float *f1, float *f2, float *alpha, int n) {
    return iftDistance1;
}
