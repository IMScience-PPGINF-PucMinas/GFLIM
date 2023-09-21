
%feature("autodoc", "2");
typedef unsigned short ushort;



typedef struct {
    float *val;
} iftBand;

%extend iftBand {

};



typedef struct {
    int x, y, z, t;
} iftVoxel;

%extend iftVoxel {

};



typedef struct {
    int xsize;
    int ysize;
    int zsize;
    int tsize;
} iftImageDomain;

%extend iftImageDomain {

};

