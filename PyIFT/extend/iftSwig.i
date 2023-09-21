%module pyift

%{
#define SWIG_FILE_WITH_INIT
#include "ift.h"
%}

#ifndef PYIFT_DEBUG
// #define PYIFT_DEBUG
#endif

%include "numpy.i"
%include "typemaps.i"
%fragment("NumPy_Fragments");

%init %{
    import_array();
%}

%rename("%(strip:[ift])s") "";

#include <string.h>
#include <Python.h>


/************** FRAGMENTS ******************/

%fragment("PyUnicodeToString", "header") {
    char* PyUnicodeToString(PyObject *obj){
        PyObject *temp_bytes = PyUnicode_AsEncodedString(obj, "ASCII", "strict");
        return PyBytes_AS_STRING(temp_bytes);
    }
}

/************** TYPEMAPS ******************/

%typemap(in, fragment="PyUnicodeToString") (const char *filename) {
    $1 = PyUnicodeToString($input);
};

%typemap(in, fragment="PyUnicodeToString") (const char *filepath) {
    $1 = PyUnicodeToString($input);
};

%typemap(in, fragment="PyUnicodeToString") (const char *format) {
    $1 = PyUnicodeToString($input);
};

%typemap(out) (const char *) {
   $result = PyString_FromString($1);
};

// reference https://stackoverflow.com/questions/36184402/how-to-apply-a-swig-typemap-for-a-double-pointer-struct-argument
%typemap(in) iftSet ** (iftSet *temp)
{
    // Alternatively, check if $input is a 0 integer `PyObject`...
    if (!SWIG_IsOK(SWIG_ConvertPtr($input, (void **) &temp, $*1_descriptor, SWIG_POINTER_DISOWN)))
        temp = NULL;
    $1 = &temp;
}

%typemap(argout) iftSet ** {
    PyObject* temp = NULL;
    if ($result == Py_None) {
        $result = PyList_New(0);
    } else if (!PyList_Check($result)) {
        temp = $result;
        $result = PyList_New(1);
        PyList_SetItem($result, 0, temp);
    }
    if (SWIG_as_voidptr($1) != NULL && SWIG_as_voidptr(*$1) !=  NULL) {
        temp = SWIG_NewPointerObj(SWIG_as_voidptr(*$1), SWIGTYPE_p_iftSet, SWIG_POINTER_OWN);
        PyList_Append($result, temp);
        Py_DECREF(temp);
    }
}

%typemap(in) iftImage ** (iftImage *temp)
{
    // Alternatively, check if $input is a 0 integer `PyObject`...
    if (!SWIG_IsOK(SWIG_ConvertPtr($input, (void **) &temp, $*1_descriptor, SWIG_POINTER_DISOWN)))
        temp = NULL;
    $1 = &temp;
}

%typemap(argout) iftImage ** {
    PyObject* temp = NULL;
    if ($result == Py_None) {
        $result = PyList_New(0);
    } else if (!PyList_Check($result)) {
        temp = $result;
        $result = PyList_New(1);
        PyList_SetItem($result, 0, temp);
    }
    if (SWIG_as_voidptr($1) != NULL && SWIG_as_voidptr(*$1) !=  NULL) {
        temp = SWIG_NewPointerObj(SWIG_as_voidptr(*$1), SWIGTYPE_p_iftImage, SWIG_POINTER_OWN);
        PyList_Append($result, temp);
        Py_DECREF(temp);
    }
}

%typemap(in) iftFImage ** (iftFImage *temp)
{
    // Alternatively, check if $input is a 0 integer `PyObject`...
    if (!SWIG_IsOK(SWIG_ConvertPtr($input, (void **) &temp, $*1_descriptor, SWIG_POINTER_DISOWN)))
        temp = NULL;
    $1 = &temp;
}

%typemap(argout) iftFImage ** {
    PyObject* temp = NULL;
    if ($result == Py_None) {
        $result = PyList_New(0);
    } else if (!PyList_Check($result)) {
        temp = $result;
        $result = PyList_New(1);
        PyList_SetItem($result, 0, temp);
    }
    if (SWIG_as_voidptr($1) != NULL && SWIG_as_voidptr(*$1) !=  NULL) {
        temp = SWIG_NewPointerObj(SWIG_as_voidptr(*$1), SWIGTYPE_p_iftFImage, SWIG_POINTER_OWN);
        PyList_Append($result, temp);
        Py_DECREF(temp);
    }
}

/************** OTHERS ******************/

typedef float (*iftKnnGraphCutFun)(iftKnnGraph *graph);

typedef float (*iftArcWeightFun)(float *f1, float *f2, float *alpha, int n);

typedef double (iftKernelFunction)(const double *x1, const double *x2, int d, double param1, double param2);
