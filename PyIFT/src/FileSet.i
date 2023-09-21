


typedef struct {
    /** Number of Files */
    long n;
    /** Array of Files */
    iftFile **files;
} iftFileSet;

%extend iftFileSet {

	~iftFileSet() {
		iftFileSet* ptr = ($self);
		iftDestroyFileSet(&ptr);
	}
	PyObject *GetPaths(void) {
	    iftFileSet *fset = ($self);
	
	    if (fset == NULL) {
	        return PyList_New(0); // return an empty List
	    }
	    else {
	        PyObject *list = PyList_New(fset->n);
	
	        for (int i = 0; i < fset->n; i++)
	            PyList_SetItem(list, i, PyString_FromString(fset->files[i]->path));
	        return list;
	    }
	}
	

};

%newobject iftLoadFileSetFromDirBySuffix;
%feature("autodoc", "2");
iftFileSet *iftLoadFileSetFromDirBySuffix(const char *dir_pathname, const char *ext, int depth);

%newobject iftLoadFileSetFromDirOrCSV;
%feature("autodoc", "2");
iftFileSet *iftLoadFileSetFromDirOrCSV(const char *file_entry, long hier_levels, bool sort_pathnames);

