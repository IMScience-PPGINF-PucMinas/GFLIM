


typedef struct {
    /** The leftmost and topmost point of the Bounding Box */
    iftVoxel begin;
    /** The rightmost and bottommost point of the Bounding Box */
    iftVoxel end;
} iftBoundingBox;

%extend iftBoundingBox {

};



typedef struct {
    iftBoundingBox *val;
    long n;
} iftBoundingBoxArray;

%extend iftBoundingBoxArray {

	~iftBoundingBoxArray() {
		iftBoundingBoxArray* ptr = ($self);
		iftDestroyBoundingBoxArray(&ptr);
	}
};

%feature("autodoc", "2");
long iftBoundingBoxBoxVolume(iftBoundingBox bb);

%feature("autodoc", "2");
void iftPrintBoundingBox(iftBoundingBox bb);

%newobject iftCentralizeBoundingBox;
%feature("autodoc", "2");
iftBoundingBox iftCentralizeBoundingBox(iftBoundingBox bb, iftVoxel new_center);

