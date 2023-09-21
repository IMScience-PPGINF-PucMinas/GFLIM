%include "iftMemory.i"
%include "iftDataTransference.i"



typedef enum {
    IFT_INTERIOR=0, IFT_EXTERIOR=1, IFT_BOTH=2
} iftSide;



typedef enum {
    IFT_AXIS_X, IFT_AXIS_Y, IFT_AXIS_Z
} iftAxis;

%feature("autodoc", "2");
char *iftFormattedTime(float runtime);

%feature("autodoc", "2");
void iftRandomSeed(unsigned int);

%feature("autodoc", "2");
long iftNormalizationValue(long maxval);

