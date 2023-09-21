
%feature("autodoc", "2");
int iftLocalMaximum(const iftImage *weights, int p, const iftAdjRel *disk);

%feature("autodoc", "2");
int iftFurthestInError(const iftImage *source, const iftImage *target, const iftImage *mask, int p);

%newobject iftQueryForAnchorsPosition;
%feature("autodoc", "2");
iftSet *iftQueryForAnchorsPosition(const iftImage *gt_contour, const iftImage *source_mask,
                                   const iftImage *gt_mask, const iftSet *anchors);

%newobject iftFurtherThanThreshold;
%feature("autodoc", "2");
iftSet *iftFurtherThanThreshold(const iftVoxelArray *anchors,
                                const iftImage *mask,
                                float threshold);

