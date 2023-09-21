%include "iftCommon.i"
%include "iftDataSet.i"
%include "iftAdjSet.i"
%include "iftMSPS.i"
%include "iftMST.i"
%include "iftSort.i"
%include "iftMetrics.i"



typedef struct {
    /** List of nodes in the graph */
    iftKnnNode *node;
    /** List of path value of the nodes */
    float      *pathval;
    /** List of nodes ordered by its path value */
    int        *ordered_nodes;
    /** Number of nodes of the graph */
    int         nnodes;
    /** Maximum arc weight for each value of k */
    float      *maxarcw;
    /** Maximum number of neighbors possible for each node */
    int         kmax;
    /** Best number of neighbors for each node */
    int         k;
    /** Priority queue */
    iftFHeap   *Q;
    /** Corresponding dataset */
    iftDataSet *Z;
} iftKnnGraph;

%extend iftKnnGraph {

	
	PyObject *AsDict(void)
	{
	    iftKnnGraph *graph = ($self);
	
	    PyObject *dict = PyDict_New();
	    for (int i = 0; i < graph->nnodes; i++)
	    {
	        PyObject *list = PyList_New(graph->kmax);
	        int j = 0;
	        for (iftAdjSet *adj = graph->node[i].adj; adj != NULL; adj = adj->next)
	        {
	            PyList_SetItem(list, j, PyInt_FromLong((long) adj->node));
	            ++j;
	        }
	        PyDict_SetItem(dict, PyInt_FromLong((long) i), list);
	    }
	
	    return dict;
	}
	

	~iftKnnGraph() {
		iftKnnGraph* ptr = ($self);
		iftDestroyKnnGraph(&ptr);
	}
};

%newobject iftCreateKnnGraph;
%feature("autodoc", "2");
iftKnnGraph *iftCreateKnnGraph(iftDataSet *Z, int kmax);

%feature("autodoc", "2");
iftKnnGraph *iftReadKnnGraph(const char *format, ...);

%feature("autodoc", "2");
void iftWriteKnnGraph(const iftKnnGraph *graph, const char *path, ...);

%feature("autodoc", "2");
int  iftUnsupTrain(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);

%newobject iftNormalizedCut;
%feature("autodoc", "2");
float  iftNormalizedCut(iftKnnGraph *graph);

%newobject iftExtractCentroidsFromDataSetAsDataSet;
%feature("autodoc", "2");
iftDataSet* iftExtractCentroidsFromDataSetAsDataSet(iftDataSet *Z, bool returnRealCent, bool usePrototypes);

%newobject iftMSTtoKnnGraph;
%feature("autodoc", "2");
iftKnnGraph *iftMSTtoKnnGraph(iftMST* mst, int number_neigh);

%inline %{

float (*iftNormalizedCutPtr(void))(iftKnnGraph *graph) {
    return iftNormalizedCut;
}


%}

float (*iftNormalizedCutPtr(void))(iftKnnGraph *graph) ;

