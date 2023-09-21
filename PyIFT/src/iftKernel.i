%include "iftCommon.i"
%include "iftAdjacency.i"
%include "iftMatrix.i"
%include "iftSpectrum.i"



typedef struct {
  iftAdjRel *A;
  float     *weight;
} iftKernel;

%extend iftKernel {

	~iftKernel() {
		iftKernel* ptr = ($self);
		iftDestroyKernel(&ptr);
	}
};



typedef struct {
  iftAdjRel *A;
  iftBand   *weight;
  int        nbands;
} iftMKernel;

%extend iftMKernel {

	~iftMKernel() {
		iftMKernel* ptr = ($self);
		iftDestroyMKernel(&ptr);
	}
};



typedef struct {
  iftAdjRel *A;
  iftBand  **weight; /* one kernel per row and one band per column */
  float	    *bias; /* bias provided from SVM -> one bias per Kernel Multiband */
  int        nkernels;
  int        nbands;
  iftMatrix *W; /* Whitening transform */
  float     *mean; /* for centralization before whitening (size is W->ncols) */ 
} iftMMKernel;

%extend iftMMKernel {

	~iftMMKernel() {
		iftMMKernel* ptr = ($self);
		iftDestroyMMKernel(&ptr);
	}
};

%newobject iftRandomMMKernel;
%feature("autodoc", "2");
iftMMKernel *iftRandomMMKernel(iftAdjRel *A, int nbands, int nkernels);

%newobject iftGaussianKernel;
%feature("autodoc", "2");
iftKernel   *iftGaussianKernel(float radius, float stdev);

%newobject iftGaussianKernel2D;
%feature("autodoc", "2");
iftKernel   *iftGaussianKernel2D(float radius, float stdev);

%newobject iftSobelXKernel;
%feature("autodoc", "2");
iftKernel   *iftSobelXKernel(void);

%newobject iftSobelYKernel;
%feature("autodoc", "2");
iftKernel   *iftSobelYKernel(void);

%newobject iftSobelZKernel;
%feature("autodoc", "2");
iftKernel   *iftSobelZKernel(void);

%newobject iftSobelXKernel2D;
%feature("autodoc", "2");
iftKernel   *iftSobelXKernel2D(void);

%newobject iftSobelYKernel2D;
%feature("autodoc", "2");
iftKernel   *iftSobelYKernel2D(void);

%newobject iftGabor2D;
%feature("autodoc", "2");
iftKernel   *iftGabor2D(float gw,float gh,float gx0,float gy0,float wfreq,float worient,float wphase,iftAdjRel* A);

