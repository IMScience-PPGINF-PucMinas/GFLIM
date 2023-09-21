%include "iftImage.i"
%include "iftMImage.i"
%include "iftDataSet.i"
%include "iftMetrics.i"



typedef struct {
    double *mean;
    iftSet *begin;
    iftSet *end;
    int dim;
    int size;
} iftDynamicSet;

%extend iftDynamicSet {

	iftSet *getSet(void)
	{
	    return ($self)->begin;
	}

	~iftDynamicSet() {
		iftDynamicSet* ptr = ($self);
		iftDestroyDynamicSet(&ptr);
	}
};



typedef struct {
    iftImage *label;
    iftFImage *cost;
    iftImage *root;
    iftImage *pred;
    iftImage *order;
    iftImage *delay;
    iftDynamicSet **dyn_trees;
} iftDTForest;

%extend iftDTForest {

	
	iftImage *GetLabel(void)
	{
	    return iftCopyImage(($self)->label);
	}
	
	
	iftImage *GetRoot(void)
	{
	    return iftCopyImage(($self)->root);
	}
	
	
	iftImage *GetPred(void)
	{
	    return iftCopyImage(($self)->pred);
	}
	
	
	iftImage *GetOrder(void)
	{
	    return iftCopyImage(($self)->order);
	}
	
	
	iftImage *GetDelay(void)
	{
	    return iftCopyImage(($self)->delay);
	}
	
	
	iftFImage *GetCost(bool sqrt = false)
	{
	    iftFImage *aux = ($self)->cost;
	    iftFImage *copy = iftCreateFImage(aux->xsize, aux->ysize, aux->zsize);
	
	    for (int i = 0; i < aux->n; i++)
	        copy->val[i] = aux->val[i];
	
	    if (sqrt) {
	       for (int i = 0; i < copy->n; i++) {
	           copy->val[i] = sqrtf(copy->val[i]);
	       }
	    }
	    return copy;
	}
	
	
	PyObject *GetSets(void)
	{
	    iftDTForest *forest = ($self);
	
	    PyObject *dict = PyDict_New();
	    for (int i = 0; i < forest->label->n; i++) {
	        if (forest->dyn_trees[i]) {
	            PyObject *val = SWIG_NewPointerObj(SWIG_as_voidptr(forest->dyn_trees[i]->begin),
	                                               SWIGTYPE_p_iftSet, 0);
	            PyDict_SetItem(dict, PyInt_FromLong(i), val);
	        }
	    }
	
	    return dict;
	}
	
	
	

	~iftDTForest() {
		iftDTForest* ptr = ($self);
		iftDestroyDTForest(&ptr);
	}
};

%newobject iftCreateDTForest;
%feature("autodoc", "2");
iftDTForest *iftCreateDTForest(const iftMImage *mimg, const iftAdjRel *A, const iftLabeledSet *seeds,
                               float delta, float gamma);

%newobject iftDynamicSetObjectPolicy;
%feature("autodoc", "2");
iftImage *iftDynamicSetObjectPolicy(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, bool use_dist);

%newobject iftDynamicSetRootPolicy;
%feature("autodoc", "2");
iftImage *iftDynamicSetRootPolicy(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist);

%newobject iftDynamicSetRootPolicyInMask;
%feature("autodoc", "2");
iftImage *iftDynamicSetRootPolicyInMask(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, iftImage *mask);

%newobject iftDynamicSetMinRootPolicy;
%feature("autodoc", "2");
iftImage *iftDynamicSetMinRootPolicy(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist);

%newobject iftDynamicSetWithCluster;
%feature("autodoc", "2");
iftImage *iftDynamicSetWithCluster(iftMImage *mimg, iftImage *cluster, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist);

%newobject iftDynamicSetRootEnhanced;
%feature("autodoc", "2");
iftImage *iftDynamicSetRootEnhanced(iftMImage *mimg, iftImage *objmap, iftAdjRel *A, iftLabeledSet *seeds, int h, float alpha, bool use_dist);

%newobject iftDynamicSetMinRootEnhanced;
%feature("autodoc", "2");
iftImage *iftDynamicSetMinRootEnhanced(iftMImage *mimg, iftImage *objmap, iftAdjRel *A, iftLabeledSet *seeds, int h, float alpha, bool use_dist);

%newobject iftDynTreeRoot;
%feature("autodoc", "2");
iftImage *iftDynTreeRoot(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int delta,
                         float gamma, iftImage *objmap, float alpha);

%newobject iftDynTreeClosestRoot;
%feature("autodoc", "2");
iftImage *iftDynTreeClosestRoot(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int delta,
                                float gamma, iftImage *objmap, float alpha);

%newobject iftDTRootWeightsMap;
%feature("autodoc", "2");
iftMImage *iftDTRootWeightsMap(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist);

