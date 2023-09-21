%include "iftCommon.i"
%include "iftDataSet.i"
%include "iftSort.i"
%include "iftRegion.i"



typedef struct {
  iftMSTNode *node;     // node
  int         nnodes;   // number of nodes
  float       maxarcw;  // maximum arc weight in the tree
  float       minarcw;  // minimum arc weight in the tree
  iftDataSet *Z;        // Each graph node is a training sample in Z
} iftMST;

%extend iftMST {

	~iftMST() {
		iftMST* ptr = ($self);
		iftDestroyMST(&ptr);
	}
};

%newobject iftCreateMST;
%feature("autodoc", "2");
iftMST *iftCreateMST(iftDataSet *Z);

