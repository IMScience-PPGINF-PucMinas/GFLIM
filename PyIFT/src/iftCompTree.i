%include "iftCommon.i"
%include "iftImage.i"
%include "iftAdjacency.i"

%newobject iftCreateMaxTree;
%feature("autodoc", "2");
iftCompTree *iftCreateMaxTree(const iftImage *img);

%newobject iftCreateMinTree;
%feature("autodoc", "2");
iftCompTree *iftCreateMinTree(const iftImage *img);

