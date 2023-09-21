%include "iftCommon.i"



typedef struct {
    int *dx, *dy, *dz, *dt;
    /* displacements to achieve the n adjacent voxels. */
    int n; /* number of adjacent voxels. */
} iftAdjRel;

%extend iftAdjRel {

	PyObject* __getitem__(int i){
	    PyObject *displacement = PyTuple_New(3);
	
	    PyTuple_SetItem(displacement, 0, PyInt_FromLong(($self)->dx[i]));
	    PyTuple_SetItem(displacement, 1, PyInt_FromLong(($self)->dy[i]));
	    PyTuple_SetItem(displacement, 2, PyInt_FromLong(($self)->dz[i]));
	
	    return displacement;
	}

	~iftAdjRel() {
		iftAdjRel* ptr = ($self);
		iftDestroyAdjRel(&ptr);
	}
};

%newobject iftCreateAdjRel;
%feature("autodoc", "2");
iftAdjRel *iftCreateAdjRel(int n);

%newobject iftSpheric;
%feature("autodoc", "2");
iftAdjRel *iftSpheric(float r);

%newobject iftHyperSpheric;
%feature("autodoc", "2");
iftAdjRel *iftHyperSpheric(float r);

%newobject iftSemiHyperSpheric;
%feature("autodoc", "2");
iftAdjRel *iftSemiHyperSpheric(float r, bool positive);

%newobject iftHemispheric;
%feature("autodoc", "2");
iftAdjRel *iftHemispheric(float r, char axis, int direction);

%newobject iftSphericEdges;
%feature("autodoc", "2");
iftAdjRel *iftSphericEdges(float r);

%newobject iftCircular;
%feature("autodoc", "2");
iftAdjRel *iftCircular(float r);

%newobject iftCircularEdges;
%feature("autodoc", "2");
iftAdjRel *iftCircularEdges(float r);

%newobject iftClockCircular;
%feature("autodoc", "2");
iftAdjRel *iftClockCircular(float r);

%newobject iftRightSide;
%feature("autodoc", "2");
iftAdjRel *iftRightSide(iftAdjRel *A, float r);

%newobject iftLeftSide;
%feature("autodoc", "2");
iftAdjRel *iftLeftSide(iftAdjRel *A, float r);

%newobject iftRectangular;
%feature("autodoc", "2");
iftAdjRel *iftRectangular(int xsize, int ysize);

%newobject iftRectangularWithDilation;
%feature("autodoc", "2");
iftAdjRel *iftRectangularWithDilation(int xsize, int ysize, int sx, int sy);

%newobject iftCross;
%feature("autodoc", "2");
iftAdjRel *iftCross(int xsize, int ysize);

%newobject iftCuboid;
%feature("autodoc", "2");
iftAdjRel *iftCuboid(int xsize, int ysize, int zsize);

%newobject iftHyperCuboid;
%feature("autodoc", "2");
iftAdjRel *iftHyperCuboid(int xsize, int ysize, int zsize, int tsize);

%newobject iftCuboidWithDilation;
%feature("autodoc", "2");
iftAdjRel *iftCuboidWithDilation(int xsize, int ysize, int zsize, int sx, int sy, int sz);

%newobject iftCopyAdjacency;
%feature("autodoc", "2");
iftAdjRel *iftCopyAdjacency(const iftAdjRel *A);

%feature("autodoc", "2");
void iftMaxAdjShifts(const iftAdjRel *A, int *dx, int *dy, int *dz);

%newobject iftGetAdjacentVoxel;
%feature("autodoc", "2");
iftVoxel iftGetAdjacentVoxel(const iftAdjRel *A, iftVoxel u, int adj);

