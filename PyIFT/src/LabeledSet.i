%include "iftCommon.i"
%include "iftImage.i"



typedef struct {
  int elem;
  int label;
  int marker;
  int handicap;
  struct ift_labeledset *next;
} iftLabeledSet;

%extend iftLabeledSet {

	iftLabeledSet* __add__(iftLabeledSet* lset){
	    iftLabeledSet *copy = iftCopyOrderedLabeledSet(($self));
	    iftConcatLabeledSet(&copy,&lset);
	    return copy;
	}
	
	iftLabeledSet* __sub__(iftLabeledSet* lset){
	    iftLabeledSet *copy = iftCopyOrderedLabeledSet(($self));
	    iftRemoveSubsetLabeledSet(&copy,&lset);
	    return copy;
	}
	
	PyObject* AsDict(){
	    iftLabeledSet* s = ($self);
	
	    PyObject *dict = PyDict_New();
	    while(s != NULL){
	        PyDict_SetItem(dict, PyInt_FromLong(s->elem), PyInt_FromLong(s->label));
	        s = s->next;
	    }
	
	    return dict;
	}
	
	void Print(){
	    iftLabeledSet *ptr = ($self);
	
	    while(ptr != NULL){
	        printf("%d, %d\n",ptr->elem, ptr->label);
	        ptr = ptr->next;
	    }
	}

	~iftLabeledSet() {
		iftLabeledSet* ptr = ($self);
		iftDestroyLabeledSet(&ptr);
	}
};

%inline %{

iftLabeledSet *CreateLabeledSetFromCoord(PyObject *_index, PyObject *_labels)
{
    PyArrayObject* index = obj_to_array_no_conversion(_index, NPY_INT32);
    PyArrayObject* labels = obj_to_array_no_conversion(_labels, NPY_INT32);

    if (!index || !labels) return NULL;

    if (array_numdims(index) != 1 || array_numdims(labels) != 1) {
        SWIG_Error(0, "Input must be a 1-dimensional array");
        return NULL;
    }

    if (!require_contiguous(labels) || !require_contiguous(index)) {
        SWIG_Error(0, "Inputs array is not contiguous");
        return NULL;
    }

    if (array_size(index, 0) != array_size(labels, 0)) {
        SWIG_Error(0, "Inputs must have the same length");
        return NULL;
    }

    int *d_index = array_data(index);
    int *d_labels = array_data(labels);

    iftLabeledSet *s = NULL;

    int size = array_size(index, 0);
    for (int i = 0; i < size; i++) {
        iftInsertLabeledSet(&s, d_index[i], d_labels[i]);
    }

    return s;
}


%}

%newobject CreateLabeledSetFromCoord;
%feature("autodoc", "2");
iftLabeledSet *CreateLabeledSetFromCoord(PyObject *_index, PyObject *_labels)
;

