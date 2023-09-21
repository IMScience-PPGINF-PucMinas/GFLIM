%include "iftCommon.i"



typedef enum {
    YCbCr_CSPACE,
    YCbCrNorm_CSPACE,
    RGB_CSPACE,
    RGBNorm_CSPACE,
    GRAY_CSPACE,
    GRAYNorm_CSPACE,
    WEIGHTED_YCbCr_CSPACE,
    LAB_CSPACE,
    LABNorm_CSPACE,
	LABNorm2_CSPACE,
    HSV_CSPACE,
    UNDEFINED_CSPACE
} iftColorSpace;



typedef struct {
  int val[3];
	float alpha;
} iftColor;

%extend iftColor {

	void __setitem__(int i, int value){
	    if(i >= 0 && i < 3){
	        ($self)->val[i] = value;
	    }
	}
	
	void SetValues(PyObject *tuple){
	    iftColor *color = ($self);
	
	    color->val[0] = (int)PyInt_AsLong(PyTuple_GetItem(tuple, 0));
	    color->val[1] = (int)PyInt_AsLong(PyTuple_GetItem(tuple, 1));
	    color->val[2] = (int)PyInt_AsLong(PyTuple_GetItem(tuple, 2));
	}
	
	PyObject* GetValues(){
	    iftColor *color = ($self);
	    PyObject *tuple = PyTuple_New(3);
	
	    PyTuple_SetItem(tuple, 0, PyInt_FromLong(color->val[0]));
	    PyTuple_SetItem(tuple, 1, PyInt_FromLong(color->val[1]));
	    PyTuple_SetItem(tuple, 2, PyInt_FromLong(color->val[2]));
	
	    return tuple;
	}

};



typedef struct {
  iftColor *color;
  int ncolors;
} iftColorTable;

%extend iftColorTable {

	~iftColorTable() {
		iftColorTable* ptr = ($self);
		iftDestroyColorTable(&ptr);
	}
	
	PyObject *AsNumPy(void) {
	    iftColorTable *ctb = ($self);
	
	    npy_intp dims[2] = {ctb->ncolors, 3};
	    PyArrayObject *npy = PyArray_SimpleNew(2, dims, NPY_INT32);
	
	    int *ptr = (int *) array_data(npy);
	
	    for (int c = 0; c < ctb->ncolors; c++) {
	        int i = c * 3;
	        ptr[i]   = ctb->color[c].val[0]; // R
	        ptr[i+1] = ctb->color[c].val[1]; // G
	        ptr[i+2] = ctb->color[c].val[2]; // B
	    }
	
	    return PyArray_Return(npy);
	}
	
	void __setitem__(int idx, const char *hex)
	{
	    iftColorTable *ctb = ($self);
	    
	    if (idx >= ctb->ncolors)
	        SWIG_Error(12, "Index out of iftColorTable range");
	    else
	        ctb->color[idx] = iftRGBtoYCbCr(iftHexaColorToRGB(hex), 255);
	}
	

};

%newobject iftCreateColorTable;
%feature("autodoc", "2");
iftColorTable *iftCreateColorTable(int ncolors);

%newobject iftCreateFColorMatrix;
%feature("autodoc", "2");
iftFColorMatrix *iftCreateFColorMatrix(int ncolors);

%newobject iftConvertRGBColorTableToYCbCrColorTable;
%feature("autodoc", "2");
void iftConvertRGBColorTableToYCbCrColorTable(iftColorTable *colorTable_rgb, int normalization_value);

%newobject iftConvertYCbCrColorTableToRGBColorTable;
%feature("autodoc", "2");
iftColorTable *iftConvertYCbCrColorTableToRGBColorTable(const iftColorTable *ctb, int normalization_value);

%newobject iftBlueToRedColorTable;
%feature("autodoc", "2");
iftColorTable *iftBlueToRedColorTable(int ncolors);

%newobject iftCreateRandomColorTable;
%feature("autodoc", "2");
iftColorTable *iftCreateRandomColorTable(int n_colors);

%feature("autodoc", "2");
iftColor iftRGBtoYCbCr(iftColor cin, int normalization_value);

%feature("autodoc", "2");
iftColor iftYCbCrtoRGB(iftColor cin, int normalization_value);

%feature("autodoc", "2");
iftColor iftHexaColorToRGB(const char *hexa_color);

%newobject iftCreateGradientColorTable;
%feature("autodoc", "2");
iftColorTable *iftCreateGradientColorTable(const iftStrArray *hexa_colors, int n_colors);

%newobject iftGrayColorTable;
%feature("autodoc", "2");
iftColorTable *iftGrayColorTable(int n_colors);

%newobject iftIronColorTable;
%feature("autodoc", "2");
iftColorTable *iftIronColorTable(int n_colors);

%newobject iftHotColorTable;
%feature("autodoc", "2");
iftColorTable *iftHotColorTable(int n_colors);

%newobject iftRainbowColorTable;
%feature("autodoc", "2");
iftColorTable *iftRainbowColorTable(int n_colors);

%newobject iftCategoricalColorTable;
%feature("autodoc", "2");
iftColorTable *iftCategoricalColorTable(int n_colors);

%newobject iftReverseRainbowColorTable;
%feature("autodoc", "2");
iftColorTable *iftReverseRainbowColorTable(int n_colors);

%newobject iftHeatMapColorTable;
%feature("autodoc", "2");
iftColorTable *iftHeatMapColorTable(int n_colors);

%newobject iftRedHotColorTable;
%feature("autodoc", "2");
iftColorTable *iftRedHotColorTable(int n_colors);

%newobject iftGreenHotColorTable;
%feature("autodoc", "2");
iftColorTable *iftGreenHotColorTable(int n_colors);

%newobject iftBlueHotColorTable;
%feature("autodoc", "2");
iftColorTable *iftBlueHotColorTable(int n_colors);

%newobject iftRedYellowHotColorTable;
%feature("autodoc", "2");
iftColorTable *iftRedYellowHotColorTable(int n_colors);

%newobject iftRedToBlueColorTable;
%feature("autodoc", "2");
iftColorTable *iftRedToBlueColorTable(int n_colors);

%newobject iftViridisColorTable;
%feature("autodoc", "2");
iftColorTable *iftViridisColorTable(int n_colors);

%newobject iftPlasmaColorTable;
%feature("autodoc", "2");
iftColorTable *iftPlasmaColorTable(int n_colors);

%newobject iftCategory21ColorTable;
%feature("autodoc", "2");
iftColorTable *iftCategory21ColorTable(int n_colors);

%newobject iftCategory10ColorTable;
%feature("autodoc", "2");
iftColorTable *iftCategory10ColorTable(int n_colors);

%inline %{

iftColorTable *CreateColorTableFromNumPy(PyObject *input) {
    PyArrayObject *arr = obj_to_array_no_conversion(input, NPY_INT32);
    if (!require_contiguous(arr))
        SWIG_Error(0, "Input numpy array is not contiguous");

    // (n_colors, 3) ==> 3 channels: [:, 0] = Y, [:, 0] = Cb, [:, 0] = Cr
    int n_dims = array_numdims(arr);
    if (n_dims != 2) {
        char error[512];
        sprintf(error, "Color Numpy Array must have 2 dimensions. Found: %d", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    npy_intp *shape = array_dimensions(arr); // (n_colors, 3)
    if (shape[1] != 3) {
        char error[512];
        sprintf(error, "Number of Channels in the Color Number Array must have 3. Found: %d", shape[1]);
        SWIG_Error(0, error);
        return NULL;
    }

    int *ptr = array_data(arr);

    iftColorTable *cmap = iftCreateColorTable(shape[0]);
    
    for (int c = 0; c < cmap->ncolors; c++) {
        int i = c * 3;
        cmap->color[c].val[0] = ptr[i];   // Y
        cmap->color[c].val[1] = ptr[i+1]; // Cb
        cmap->color[c].val[2] = ptr[i+2]; // Cr
    }

    return cmap;
}




%}

%newobject CreateColorTableFromNumPy;
%feature("autodoc", "2");
iftColorTable *CreateColorTableFromNumPy(PyObject *input) ;

