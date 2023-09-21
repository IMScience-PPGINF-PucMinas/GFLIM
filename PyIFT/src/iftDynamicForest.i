%include "iftFImage.i"
%include "iftMImage.i"

%newobject iftDynamicTreesObject;
%feature("autodoc", "2");
iftImage *iftDynamicTreesObject(iftMImage *features, iftLabeledSet *seeds, iftAdjRel *A);

%newobject iftDynamicTreesRoot;
%feature("autodoc", "2");
iftImage *iftDynamicTreesRoot(iftMImage *features, iftLabeledSet *seeds, iftAdjRel *A);

