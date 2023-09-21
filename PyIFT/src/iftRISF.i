%include "iftMatrix.i"
%include "iftDataSet.i"
%include "iftImage.i"
%include "iftAdjacency.i"



typedef enum {
  IFT_ISF_UNDEFINED_SAMPLING = 0,
  IFT_ISF_GRID_SAMPLING = 1,
  IFT_ISF_MIXED_SAMPLING = 2,
  IFT_ISF_RANDOM_SAMPLING = 3,
  IFT_ISF_GEODESIC_SAMPLING = 4,
  IFT_ISF_LATTICE_SAMPLING = 5,
  IFT_ISF_HEXGRID_SAMPLING = 6
} iftISFSampling;



typedef enum {
  /** Invalid fallback for NULL values. */
  IFT_ISF_UNDEFINED_PATHCOST = 0,
  IFT_ISF_ADD_ROOT_PATHCOST = 1,
  IFT_ISF_ADD_MEAN_PATHCOST = 2,
  IFT_ISF_ADD_ARC_PATHCOST = 3,
  IFT_ISF_MAX_MEAN_PATHCOST = 4
} iftISFPathCost;



typedef enum {
  /** Invalid fallback for NULL values. */
  IFT_ISF_UNDEFINED_SEEDRECOMP = 0,
  /** Medoid according to some pre-specified feature distance. */
  IFT_ISF_MEDOID_SEEDRECOMP = 1,
  /** Geometric medoid. */
  IFT_ISF_CENTROID_SEEDRECOMP = 2
} iftISFSeedRecomp;



typedef struct {
  int nSuperpixels;
  int nIters;
  int nSmoothIters;
  /** Superpixels with size below threshold are removed in the end. */
  float smallSpThreshold; /*< Value relative to mean sp size. */
  iftISFSampling sampling;
  iftISFPathCost pathCost; /*< See iftISFPathCost. */
  iftISFSeedRecomp seedRecomp; /*< See iftISFSeedRecomp. */
} iftISFOptions;

%extend iftISFOptions {

	~iftISFOptions() {
		iftISFOptions* ptr = ($self);
		iftDestroyISFOptions(&ptr);
	}
};

%newobject iftInitISFOptions;
%feature("autodoc", "2");
iftISFOptions * iftInitISFOptions(int nSuperpixels, int nIters, int nSmoothIters, float smallSpThreshold, iftISFSampling sampling, iftISFPathCost pathCost, iftISFSeedRecomp seedRecomp);



typedef struct {
  /** Number of nodes. */
  int nNodes;
  /** Node adjacencies for explicit graphs. */
  iftSet **nodeAdj;
  /** Parametric space features (e.g. LAB color). */
  iftMatrix *feats; /*< Node based */
  iftMatrix *seedFeats; /*< Superpixel based */
  /** Function used to calculate distances in parametric space. */
  iftArcWeightFun distFun; /*< float(float, float, float, int) */
  /** Alpha for selecting of weighting feats in distFun. */
  float *distAlpha;
  /** Geometric space features (e.g. centroid). */
  iftMatrix *pos; /*< Node based */
  iftMatrix *seedPos; /*< Superpixel based */
  /** Path cost parameters. */
  float alpha;
  float beta;
  /** Look up table from image to nodes & spatial domain info. */
  iftImage * refImg; /*< Generally a superpixel label map. */
  long long int * spSize; /*< w.r.t. refImg/nNodes. */
  /** General Forest annotation info indexed by nodes. */
  iftForestAnnotationInfo * ann;
  /** Adjacency relation for implicit graphs. */
  iftAdjRel *A;
  /** Flag for COMPLETE, IMPLICIT or EXPLICIT graphs. */
  char type;
} iftSuperpixelIGraph;

%extend iftSuperpixelIGraph {

	~iftSuperpixelIGraph() {
		iftSuperpixelIGraph* ptr = ($self);
		iftDestroySuperpixelIGraph(&ptr);
	}
};

%newobject iftInitSuperpixelIGraph;
%feature("autodoc", "2");
iftSuperpixelIGraph *iftInitSuperpixelIGraph(const iftImage *refImg);

%feature("autodoc", "2");
void iftSetSuperpixelIGraphImplicitAdjacency(iftSuperpixelIGraph *igraph, iftAdjRel *A);

%feature("autodoc", "2");
void iftSetSuperpixelIGraphExplicitRAG(iftSuperpixelIGraph *igraph, const iftAdjRel *A);

%feature("autodoc", "2");
void iftSetSuperpixelIGraphFeatures(iftSuperpixelIGraph *igraph, const iftMatrix *superpixelFeatures, iftArcWeightFun distFun, float alpha, float beta);

%newobject iftSuperpixelSegmentationByRISF;
%feature("autodoc", "2");
iftImage * iftSuperpixelSegmentationByRISF(iftSuperpixelIGraph *igraph, const iftISFOptions *opt, const iftImage *img);

%newobject iftComputeSuperpixelFeaturesByGeometricCenter;
%feature("autodoc", "2");
iftMatrix * iftComputeSuperpixelFeaturesByGeometricCenter(const iftImage *superpixelLabelMap);

%newobject iftComputeSuperpixelFeaturesByColorSpaceMean;
%feature("autodoc", "2");
iftMatrix * iftComputeSuperpixelFeaturesByColorSpaceMean(const iftImage *superpixelLabelMap, const iftImage *img, iftColorSpace colorSpace);

%newobject iftComputeSuperpixelFeaturesByMImageMean;
%feature("autodoc", "2");
iftMatrix * iftComputeSuperpixelFeaturesByMImageMean(const iftMImage *mimg, const iftImage *superpixelLabelMap);

