%include "iftAdjacency.i"
%include "iftCommon.i"
%include "iftDataSet.i"
%include "iftImage.i"
%include "iftMatrix.i"
%include "iftMathMorph.i"
%include "iftMImage.i"
%include "iftRegion.i"

%newobject iftExtractRemovalMarkers;
%feature("autodoc", "2");
iftSet 	      *iftExtractRemovalMarkers(iftLabeledSet **s);

%newobject iftLabelObjBorderSet;
%feature("autodoc", "2");
iftLabeledSet *iftLabelObjBorderSet(iftImage *bin, iftAdjRel *A);

%newobject iftImageBorderLabeledSet;
%feature("autodoc", "2");
iftLabeledSet *iftImageBorderLabeledSet(iftImage *img);

%newobject iftLabelCompSet;
%feature("autodoc", "2");
iftLabeledSet *iftLabelCompSet(iftImage *bin, iftAdjRel *A);

%newobject iftFuzzyModelToLabeledSet;
%feature("autodoc", "2");
iftLabeledSet *iftFuzzyModelToLabeledSet(iftImage *model);

%newobject iftMAdjustSeedCoordinates;
%feature("autodoc", "2");
iftLabeledSet *iftMAdjustSeedCoordinates(iftLabeledSet *S, iftMImage *input, iftMImage *output);

%newobject iftAdjustSeedCoordinates;
%feature("autodoc", "2");
iftLabeledSet *iftAdjustSeedCoordinates(iftLabeledSet *Sin, iftImage *orig, iftMImage *output);

%newobject iftImageBorderSet;
%feature("autodoc", "2");
iftSet        *iftImageBorderSet(const iftImage *img);

%newobject iftMultiObjectBorderLabeledSet;
%feature("autodoc", "2");
iftLabeledSet *iftMultiObjectBorderLabeledSet(iftImage *img, iftAdjRel *A);

%newobject iftObjectToSet;
%feature("autodoc", "2");
iftSet *iftObjectToSet(const iftImage *label_img);

%newobject iftObjectToLabeledSet;
%feature("autodoc", "2");
iftLabeledSet *iftObjectToLabeledSet(const iftImage *label_img);

%newobject iftObjectsToList;
%feature("autodoc", "2");
iftList *iftObjectsToList(const iftImage *label_img);

%newobject iftObjectBorderSet;
%feature("autodoc", "2");
iftSet *iftObjectBorderSet(const iftImage *label_img, iftAdjRel *Ain);

%newobject iftBackgroundBorderSet;
%feature("autodoc", "2");
iftSet *iftBackgroundBorderSet(const iftImage *label_img, iftAdjRel *Ain);

%newobject iftObjectBorderBitMap;
%feature("autodoc", "2");
iftBMap *iftObjectBorderBitMap(const iftImage *label, const iftAdjRel *A, int *n_border_spels);

%feature("autodoc", "2");
void iftMaskImageToSet(const iftImage *mask, iftSet **S);

%newobject iftMaskImageFromSet;
%feature("autodoc", "2");
  iftImage *iftMaskImageFromSet(iftSet *S, int xsize, int ysize, int zsize);

%feature("autodoc", "2");
void iftSetToImage(const iftSet *S, iftImage *img, int obj);

%feature("autodoc", "2");
void iftListToImage(const iftList *L, iftImage *img, int obj);

%newobject iftMaskToList;
%feature("autodoc", "2");
iftList *iftMaskToList(const iftImage *mask);

%feature("autodoc", "2");
void iftLabeledSetToImage(const iftLabeledSet *S, iftImage *img, bool increment_label);

%feature("autodoc", "2");
void iftIntArrayToImage(const iftIntArray *iarr, iftImage *img, int obj);

%newobject iftMaskToIntArray;
%feature("autodoc", "2");
iftIntArray *iftMaskToIntArray(const iftImage *mask);

%newobject iftRegionsToLabeledSet;
%feature("autodoc", "2");
iftLabeledSet *iftRegionsToLabeledSet(const iftImage *label, bool decrement_label);

%newobject iftEndPoints;
%feature("autodoc", "2");
iftSet        *iftEndPoints(iftImage *skel, iftAdjRel *A);

%newobject iftFindPathOnSkeleton;
%feature("autodoc", "2");
iftSet        *iftFindPathOnSkeleton(iftImage *skel, iftAdjRel *A, int src, int dst);

%newobject iftSkeletonPoints;
%feature("autodoc", "2");
iftSet        *iftSkeletonPoints(iftImage *skel);

%newobject iftObjectBorders;
%feature("autodoc", "2");
iftImage *iftObjectBorders(const iftImage *label_img, const iftAdjRel *Ain, bool keep_border_labels,
                           bool include_image_frame);

%newobject iftFindAndLabelObjectBorders;
%feature("autodoc", "2");
iftImage *iftFindAndLabelObjectBorders(const iftImage *label_img, const iftAdjRel *Ain);

%newobject iftEasyLabelComp;
%feature("autodoc", "2");
iftImage      *iftEasyLabelComp(iftImage *bin, iftAdjRel *A);

%newobject iftLabelComp;
%feature("autodoc", "2");
iftImage      *iftLabelComp(iftImage *bin, iftAdjRel *A);

%newobject iftSelectLargestComp;
%feature("autodoc", "2");
iftImage *iftSelectLargestComp(const iftImage *bin, const iftAdjRel *Ain);

%newobject iftSelectSmallestComp;
%feature("autodoc", "2");
iftImage      *iftSelectSmallestComp(iftImage *bin, iftAdjRel *A);

%newobject iftSelectKLargestComp;
%feature("autodoc", "2");
iftImage      *iftSelectKLargestComp(iftImage *bin, iftAdjRel *A, int K);

%newobject iftSelectKSmallestComp;
%feature("autodoc", "2");
iftImage      *iftSelectKSmallestComp(iftImage *bin, iftAdjRel *A, int K);

%newobject iftComponentArea;
%feature("autodoc", "2");
iftImage      *iftComponentArea(iftImage *bin, iftAdjRel *A);

%newobject iftSelectCompAboveArea;
%feature("autodoc", "2");
iftImage      *iftSelectCompAboveArea(iftImage *bin, iftAdjRel *A, int thres);

%newobject iftSelectCompBelowArea;
%feature("autodoc", "2");
iftImage      *iftSelectCompBelowArea(iftImage *bin, iftAdjRel *A, int thres);

%newobject iftSelectCompInAreaInterval;
%feature("autodoc", "2");
iftImage      *iftSelectCompInAreaInterval(iftImage *bin, iftAdjRel *A, int thres_min, int thres_max);

%newobject iftRegionalMaxima;
%feature("autodoc", "2");
iftImage      *iftRegionalMaxima(iftImage *img);

%newobject iftRegionalMinima;
%feature("autodoc", "2");
iftImage      *iftRegionalMinima(iftImage *img);

%newobject iftRegionalMaximaInRegion;
%feature("autodoc", "2");
iftImage      *iftRegionalMaximaInRegion(iftImage *img, iftImage *mask);

%newobject iftRegionalMinimaInRegion;
%feature("autodoc", "2");
iftImage      *iftRegionalMinimaInRegion(iftImage *img, iftImage *mask);

%newobject iftRootVoxels;
%feature("autodoc", "2");
iftImage      *iftRootVoxels(iftImage *pred);

%newobject iftLeafVoxels;
%feature("autodoc", "2");
iftImage      *iftLeafVoxels(iftImage *pred, iftAdjRel *A);

%newobject iftLeafVoxelsOnMask;
%feature("autodoc", "2");
iftImage      *iftLeafVoxelsOnMask(iftImage *pred, iftAdjRel *A, iftImage *mask);

%newobject iftRegionArea;
%feature("autodoc", "2");
iftImage      *iftRegionArea(iftImage *label);

%newobject iftSelectLargestRegion;
%feature("autodoc", "2");
iftImage *iftSelectLargestRegion(const iftImage *label_img);

%newobject iftSelectSmallestRegion;
%feature("autodoc", "2");
iftImage      *iftSelectSmallestRegion(iftImage *label);

%newobject iftSelectKLargestRegions;
%feature("autodoc", "2");
iftImage      *iftSelectKLargestRegions(iftImage *label, int K);

%newobject iftSelectKLargestRegionsAndPropagateTheirLabels;
%feature("autodoc", "2");
iftImage      *iftSelectKLargestRegionsAndPropagateTheirLabels(iftImage *label, iftAdjRel *A, int K);

%newobject iftSelectKSmallestRegions;
%feature("autodoc", "2");
iftImage      *iftSelectKSmallestRegions(iftImage *label, int K);

%newobject iftSelectRegionsAboveArea;
%feature("autodoc", "2");
iftImage      *iftSelectRegionsAboveArea(iftImage *label, int thres);

%newobject iftSelectRegionsAboveAreaAndPropagateTheirLabels;
%feature("autodoc", "2");
iftImage      *iftSelectRegionsAboveAreaAndPropagateTheirLabels(iftImage *label, int thres);

%newobject iftSelectRegionsBelowArea;
%feature("autodoc", "2");
iftImage      *iftSelectRegionsBelowArea(iftImage *label, int thres);

%newobject iftSelectRegionsInAreaInterval;
%feature("autodoc", "2");
iftImage      *iftSelectRegionsInAreaInterval(iftImage *label, int min_area, int max_area);

%newobject iftSelectRegionsInAreaIntervalInplace;
%feature("autodoc", "2");
void           iftSelectRegionsInAreaIntervalInplace(iftImage *label, int min_area, int max_area);

%newobject iftLabelContPixel;
%feature("autodoc", "2");
iftImage  *iftLabelContPixel(iftImage *bin);

%newobject iftReadSeeds;
%feature("autodoc", "2");
iftLabeledSet *iftReadSeeds(const iftImage *img, const char *filename, ...);

%newobject iftMReadSeeds;
%feature("autodoc", "2");
iftLabeledSet *iftMReadSeeds(iftMImage *img, char *filename);

%feature("autodoc", "2");
void iftMWriteSeeds(iftLabeledSet *seed, iftMImage *image, char *filename);

%feature("autodoc", "2");
void iftWriteSeeds(iftLabeledSet *seed, const iftImage *image, const char *filename, ...);

%newobject iftSeedImageFromLabeledSet;
%feature("autodoc", "2");
iftImage	  *iftSeedImageFromLabeledSet(iftLabeledSet* labeled_set, iftImage *image);

%newobject iftLabeledSetFromSeedImage;
%feature("autodoc", "2");
iftLabeledSet *iftLabeledSetFromSeedImage(iftImage* seed_image, bool decrement);

%newobject iftLabeledSetFromSeedImageMarkersAndHandicap;
%feature("autodoc", "2");
iftLabeledSet *iftLabeledSetFromSeedImageMarkersAndHandicap(iftImage* seed_image, iftImage *marker, iftImage *handicap);

%newobject iftGetSeeds;
%feature("autodoc", "2");
iftLabeledSet* iftGetSeeds(iftLabeledSet* S, int nelem, int label);

%newobject iftGetMisclassifiedSeeds;
%feature("autodoc", "2");
iftLabeledSet* iftGetMisclassifiedSeeds(iftLabeledSet* S, int nelem, int label, iftImage* gt_image, iftImage* cl_image);

%newobject iftFastLabelComp;
%feature("autodoc", "2");
iftImage *iftFastLabelComp(const iftImage *bin, const iftAdjRel *Ain);

%newobject iftBinaryMaskToSet;
%feature("autodoc", "2");
iftSet *iftBinaryMaskToSet(iftImage *mask);

%newobject iftHBasins;
%feature("autodoc", "2");
iftImage *iftHBasins(iftImage *img, int H);

%newobject iftHDomes;
%feature("autodoc", "2");
iftImage *iftHDomes(iftImage *img, int H);

%newobject iftGridSamplingOnMask;
%feature("autodoc", "2");
iftIntArray *iftGridSamplingOnMask(const iftImage *bin_mask, float radius,
                                   int initial_obj_voxel_idx, long n_samples);

%newobject iftGridSamplingForPatchExtraction;
%feature("autodoc", "2");
iftIntArray *iftGridSamplingForPatchExtraction(iftImageDomain img_dom, int patch_xsize, int patch_ysize,
                                               int patch_zsize, int stride_x, int stride_y, int stride_z);

%newobject iftBoundingBoxesAroundVoxels;
%feature("autodoc", "2");
iftBoundingBoxArray *iftBoundingBoxesAroundVoxels(const iftImage *img, const iftIntArray *voxel_indices, int size);

%newobject iftLabelComponentsBySeeds;
%feature("autodoc", "2");
iftImage *iftLabelComponentsBySeeds(iftImage *comp, iftLabeledSet *seeds, bool incr);

%newobject iftDataSetToLabeledSeeds;
%feature("autodoc", "2");
iftLabeledSet *iftDataSetToLabeledSeeds(iftDataSet *Z, iftImage *comp);

%newobject iftConnectObjectSeeds;
%feature("autodoc", "2");
iftLabeledSet *iftConnectObjectSeeds(iftImage *img, iftLabeledSet *seeds);

%newobject iftPropagateSeedsToCluster;
%feature("autodoc", "2");
iftLabeledSet *iftPropagateSeedsToCluster(const iftDataSet *Z, int truelabel, float purity);

%newobject iftSelectSeedsForEnhancement;
%feature("autodoc", "2");
iftLabeledSet *iftSelectSeedsForEnhancement(iftLabeledSet *seeds, iftImage *groups, int label, float threshold);

%newobject iftSeedsFeatures;
%feature("autodoc", "2");
iftDoubleMatrix *iftSeedsFeatures(iftLabeledSet *seeds, const iftMImage *mimg, int label);

%feature("autodoc", "2");
void iftSeedsSuperpixelBins(iftLabeledSet *seeds, const iftMImage *mimg, const iftImage *superpixel,
                            iftDoubleMatrix *data, iftIntArray *label);

%newobject iftDrawDilatedSeeds;
%feature("autodoc", "2");
iftImage *iftDrawDilatedSeeds(const iftImage *image, const iftLabeledSet *seeds,
                              const iftAdjRel *A, const iftColorTable *ctb_rgb);

%newobject iftRemoveFrameComp;
%feature("autodoc", "2");
  iftImage *iftRemoveFrameComp(iftImage *bin);

