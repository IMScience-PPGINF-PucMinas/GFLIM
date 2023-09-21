%include "iftAdjacency.i"
%include "iftCommon.i"
%include "iftCompTree.i"
%include "iftImage.i"
%include "iftGraphics.i"
%include "iftMImage.i"
%include "iftSegmentation.i"



typedef struct {
  iftINode   *node;    /* list of graph nodes */
  int         nnodes;  /* number of graph nodes */
  int         nfeats;  /* number of image features */
  iftImage   *index;   /* node index */
  float     **feat;    /* image features */
  int        *label, *marker, *root, *pred; /* forest annotation */
  double      *pvalue;  /* forest annotation */
  iftAdjRel  *A;       /* adjacency relation (implicit graphs only) */
  char        type;    /* COMPLETE, IMPLICIT, EXPLICIT */
} iftIGraph;

%extend iftIGraph {

	~iftIGraph() {
		iftIGraph* ptr = ($self);
		iftDestroyIGraph(&ptr);
	}
};

%newobject iftImplicitIGraph;
%feature("autodoc", "2");
  iftIGraph *iftImplicitIGraph(iftMImage *img, const iftImage *mask, iftAdjRel *A);

%newobject iftExplicitIGraph;
%feature("autodoc", "2");
  iftIGraph *iftExplicitIGraph(const iftMImage *img, const iftImage *mask, const iftImage *label, iftAdjRel *A);

%newobject iftSpatialKnnIGraph;
%feature("autodoc", "2");
  iftIGraph *iftSpatialKnnIGraph(iftMImage *img, iftImage *mask, iftAdjRel *A, int K);

%newobject iftSpatialIGraph;
%feature("autodoc", "2");
iftIGraph *iftSpatialIGraph(iftMImage *img, iftImage *mask, iftAdjRel *A, float df);

%newobject iftKnnIGraph;
%feature("autodoc", "2");
  iftIGraph *iftKnnIGraph(iftMImage *img, iftImage *mask, int K);

%newobject iftIGraphPathValue;
%feature("autodoc", "2");
  iftImage *iftIGraphPathValue(iftIGraph *igraph);

%newobject iftIGraphWeight;
%feature("autodoc", "2");
  iftFImage *iftIGraphWeight(iftIGraph *igraph);

%newobject iftIGraphLabel;
%feature("autodoc", "2");
  iftImage *iftIGraphLabel(iftIGraph *igraph);

%newobject iftIGraphRoot;
%feature("autodoc", "2");
  iftImage *iftIGraphRoot(iftIGraph *igraph);

%feature("autodoc", "2");
  int iftIGraphISF_Mean(iftIGraph *igraph, iftImage *seeds, double alpha, double beta0, int niters);

%feature("autodoc", "2");
  int iftIGraphISF_Root(iftIGraph *igraph, iftImage *seeds, double alpha, double beta0, int niters);

%newobject iftExtract_ISF_MIX_ROOT_Superpixels;
%feature("autodoc", "2");
iftImage *iftExtract_ISF_MIX_ROOT_Superpixels(iftImage *img, int nsuperpixels, float alpha, float beta, int niters, int smooth_niters);

%newobject iftExtract_ISF_MIX_MEAN_Superpixels;
%feature("autodoc", "2");
iftImage *iftExtract_ISF_MIX_MEAN_Superpixels(iftImage *img, int nsuperpixels, float alpha, float beta, int niters, int smooth_niters);

%newobject iftSuperpixelCenterSetFromIGraph;
%feature("autodoc", "2");
iftSet *iftSuperpixelCenterSetFromIGraph(iftIGraph *igraph);

