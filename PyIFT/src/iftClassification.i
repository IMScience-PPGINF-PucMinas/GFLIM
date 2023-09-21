%include "iftCommon.i"
%include "iftClustering.i"
%include "iftDataSet.i"
%include "iftMetrics.i"
%include "iftMST.i"
%include "iftDeepLearning.i"
%include "iftSVM.i"
%include "iftArgumentList.i"
%include "iftDistanceFunctions.i"
%include "iftGenericVector.i"



typedef struct {
	iftCplNode *node;          // list of nodes
	float      *pathval;       // path value for each node
	int        *ordered_nodes; // list of nodes ordered by path value
	int         nnodes;        // number of nodes
	iftFHeap   *Q;             // priority queue
	iftDataSet *Z;             // Each graph node is a training sample in Z
} iftCplGraph;

%extend iftCplGraph {

	PyObject *AsDict(void)
	{
	    iftCplGraph *graph = ($self);
	
	    PyObject *dict = PyDict_New();
	    for (int i = 0; i < graph->nnodes; i++)
	    {
	        PyObject *list = PyList_New(1);
	        PyList_SetItem(list, 0, PyLong_FromLong((long) graph->node[i].pred));
	        PyDict_SetItem(dict, PyLong_FromLong((long) i), list);
	    }
	
	    return dict;
	}
	

	~iftCplGraph() {
		iftCplGraph* ptr = ($self);
		iftDestroyCplGraph(&ptr);
	}
};

%feature("autodoc", "2");
void iftWriteCplGraph(const iftCplGraph *graph, const char *pathname);

%feature("autodoc", "2");
iftCplGraph* iftReadCplGraph(const char* pathname);

%newobject iftCreateCplGraph;
%feature("autodoc", "2");
iftCplGraph *iftCreateCplGraph(iftDataSet *Z);

%feature("autodoc", "2");
void         iftSupTrain(iftCplGraph *graph);

%feature("autodoc", "2");
int iftClassify(const iftCplGraph *graph, iftDataSet *Ztest);

%newobject iftClassifyByOneClassOPF;
%feature("autodoc", "2");
void iftClassifyByOneClassOPF(iftDataSet *Ztest, const iftKnnGraph *graph,
                              float quantile_of_k);

%feature("autodoc", "2");
int          iftClassifyWithCertaintyValues(const iftCplGraph *graph, iftDataSet *Ztest);

%newobject iftSemiSupTrain;
%feature("autodoc", "2");
iftCplGraph* iftSemiSupTrain(iftDataSet *Ztrain);

