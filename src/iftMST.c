#include "iftMST.h"

#include "ift/core/io/Stream.h"


iftMST *iftCreateMST(iftDataSet *Z)
{
	iftMST *mstree=(iftMST *)iftAlloc(1,sizeof(iftMST));
	int u,v,s,t;
	iftFHeap *Q;
	float *cost;
	int   *pred;
	float  arcw;

	if (Z->ntrainsamples == 0){
        iftError("No samples for MST generation", "iftCreateMST");
	}
	mstree->node = (iftMSTNode *)iftAlloc(Z->ntrainsamples,sizeof(iftMSTNode));
	if (mstree->node == NULL)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateMST");

	mstree->nnodes  = Z->ntrainsamples;
	mstree->maxarcw = IFT_INFINITY_FLT_NEG;
	mstree->minarcw = IFT_INFINITY_FLT;
	mstree->Z       = Z;

	for (u=0; u < mstree->nnodes; u++) {
		mstree->node[u].adj = NULL;
		mstree->node[u].maxarcadj = IFT_NIL;
		mstree->node[u].color = IFT_WHITE;
	}

	int j = 0;
	for (int s=0; s < Z->nsamples; s++){
	  if (iftHasSampleStatus(Z->sample[s], IFT_TRAIN)){
	    Z->sample[j].weight = 0.0;
	    mstree->node[j].sample = s;
	    j++;
	  }
	}

	pred    = iftAllocIntArray(mstree->nnodes);
	cost    = iftAllocFloatArray(mstree->nnodes);
	Q       = iftCreateFHeap(mstree->nnodes,cost);

	for (u=0; u < mstree->nnodes; u++) {
	  cost[u]= IFT_INFINITY_FLT; pred[u]= IFT_NIL;
	}
	cost[0]=0;
	iftInsertFHeap(Q,0);
	mstree->maxarcw = 0.0;

	while ( !iftEmptyFHeap(Q) ) {

		u = iftRemoveFHeap(Q);
		s = mstree->node[u].sample;
		v = pred[u];

		if (v != IFT_NIL){ /* Create the arcs of the MST and compute sample weights */
			t    = mstree->node[v].sample;

			if(iftDist == NULL)
				arcw = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
			else
				arcw = iftDist->distance_table[s][t];

			iftInsertSet(&mstree->node[u].adj,v);
			iftInsertSet(&mstree->node[v].adj,u);
			if (arcw > mstree->maxarcw)
				mstree->maxarcw = arcw;
			if (arcw < mstree->minarcw)
				mstree->minarcw = arcw;
			if (arcw > Z->sample[s].weight)
				Z->sample[s].weight = arcw;
			if (arcw > Z->sample[t].weight)
				Z->sample[t].weight = arcw;
		}

		/* Propagate paths */

		for (v=0; v < mstree->nnodes; v++){
			if (Q->color[v] != IFT_BLACK){
				t = mstree->node[v].sample;

				if(iftDist == NULL)
					arcw = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
				else
					arcw = iftDist->distance_table[s][t];

				if ( arcw < cost[v] ) {
					pred[v] = u;
					cost[v] = arcw;
					if(Q->color[v] == IFT_WHITE)
						iftInsertFHeap(Q, v);
					else
						iftGoUpFHeap(Q, Q->pos[v]);
				}
			}
		}
	}
	iftDestroyFHeap(&Q);

	/* update maxarcadj */
	for (u=0; u < mstree->nnodes; u++) {
		iftSet *set = mstree->node[u].adj;
		if (set != NULL) {
			mstree->node[u].maxarcadj = set->elem;
			set = set->next;
			while (set != NULL) {
				if (cost[mstree->node[u].maxarcadj] < cost[set->elem]) {
					mstree->node[u].maxarcadj = set->elem;
				}
				set = set->next;
			}
		}
	}

	iftFree(cost);
	iftFree(pred);

	return(mstree);
}

//Creates MST from the complete graph defined by the given samples
iftMST *iftCreateMSTFromSet(iftDataSet *Z, iftSet* samples){
	iftMST *mstree=(iftMST *)iftAlloc(1,sizeof(iftMST));

	//Number of nodes in the MST
	int size = iftSetSize(samples);

	//Allocates space for mstree
	mstree->node = (iftMSTNode *)iftAlloc(size,sizeof(iftMSTNode));
	mstree->nnodes  = size;
	mstree->maxarcw = IFT_INFINITY_FLT_NEG;
	mstree->minarcw = IFT_INFINITY_FLT;
	mstree->Z       = Z;

	//Prepares the data structures
	int i;
	for (i = 0; i < mstree->nnodes; i++) {
		mstree->node[i].adj = NULL;
		mstree->node[i].maxarcadj = IFT_NIL;
		mstree->node[i].color = IFT_WHITE;
	}

	//Creates a correspondence between samples in the dataset and nodes in the MST
	iftSet *s = samples;
	i = 0;
	while(s){
		int sample = s->elem;
		Z->sample[sample].weight = 0.0;
		mstree->node[i].sample = sample;

		s = s->next;
		i++;
	}

	int *pred = iftAllocIntArray(mstree->nnodes);
	float *cost = iftAllocFloatArray(mstree->nnodes);
	iftFHeap *Q = iftCreateFHeap(mstree->nnodes,cost);

	for (i = 0; i < mstree->nnodes; i++){
		cost[i] = IFT_INFINITY_FLT;
		pred[i] = IFT_NIL;
	}

	cost[0] = 0;
	iftInsertFHeap(Q,0);
	mstree->maxarcw = 0.0;

	while ( !iftEmptyFHeap(Q) ) {
		//u is the current node
		int u = iftRemoveFHeap(Q);
		//s is the corresponding sample in the dataset
		int s = mstree->node[u].sample;

		int v = pred[u];

		if (v != IFT_NIL){ /* Create the arcs of the MST and compute sample weights */
			//t is the predecessor's corresponding sample in the dataset
			int t    = mstree->node[v].sample;
			float arcw = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
			iftInsertSet(&mstree->node[u].adj,v);
			iftInsertSet(&mstree->node[v].adj,u);
			if (arcw > mstree->maxarcw)
				mstree->maxarcw = arcw;
			if (arcw < mstree->minarcw)
				mstree->minarcw = arcw;
			if (arcw > Z->sample[s].weight)
				Z->sample[s].weight = arcw;
			if (arcw > Z->sample[t].weight)
				Z->sample[t].weight = arcw;
		}

		/* Propagate paths */

		for (v=0; v < mstree->nnodes; v++){
			if (Q->color[v] != IFT_BLACK){
				int t = mstree->node[v].sample;
				float arcw = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
				if ( arcw < cost[v] ) {
					pred[v] = u;
					cost[v] = arcw;
					if(Q->color[v] == IFT_WHITE)
						iftInsertFHeap(Q, v);
					else
						iftGoUpFHeap(Q, Q->pos[v]);
				}
			}
		}
	}
	iftDestroyFHeap(&Q);

	/* update maxarcadj */
	int u;
	for (u= 0; u < mstree->nnodes; u++) {
		iftSet *set = mstree->node[u].adj;
		if (set != NULL) {
			mstree->node[u].maxarcadj = set->elem;
			set = set->next;
			while (set != NULL) {
				if (cost[mstree->node[u].maxarcadj] < cost[set->elem]) {
					mstree->node[u].maxarcadj = set->elem;
				}
				set = set->next;
			}
		}
	}

	iftFree(cost);
	iftFree(pred);

	return(mstree);

}

void iftDestroyMST(iftMST **mstree)
{
	iftMST *aux=*mstree;
	int u;

	if (aux != NULL) {
		for (u=0; u < aux->nnodes; u++)
			iftDestroySet(&aux->node[u].adj);
		iftFree(aux->node);
		iftFree(aux);
		*mstree=NULL;
	}
}

void iftNormalizeSampleWeightMST(iftMST *mstree)
{
	int s;
	iftDataSet *Z = mstree->Z;

	if ((mstree->maxarcw - mstree->minarcw) > 0.0)
		for (s=0; s < Z->nsamples; s++) {
			Z->sample[s].weight = ((Z->sample[s].weight-mstree->minarcw)/
								   (mstree->maxarcw-mstree->minarcw));
			Z->sample[s].weight = ((IFT_MAXWEIGHT - 1.0) * Z->sample[s].weight) + 1.0;
		}
}

int iftSelectSupTrainSamplesByMST(iftDataSet *Z, float train_perc)
{
	iftMST *mstree;
	int     c, u, v, s, t, nsamples, *class_size;
	iftSet *adj;
	int    *index;
	float  *weight;

	if (Z->nclasses == 0)
        iftError("There are no classes", "iftSelectSupTrainSamplesByMST");


	iftResetDataSet(Z);

	/* Verify the trivial case */

	if (train_perc == 1.0) {
		for (s=0; s < Z->nsamples; s++)
			iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
		Z->ntrainsamples = Z->nsamples;
		return(Z->ntrainsamples);
	}

	/* Compute the size of each class */

	class_size = iftAllocIntArray(Z->nclasses+1);
	for (s=0; s < Z->nsamples; s++)
		class_size[Z->sample[s].truelabel] ++;

	/* Verify number of available samples per class and update
       class_size with it */

	Z->ntrainsamples = 0;
	for (c=1; c <= Z->nclasses; c++) {
		nsamples = (int)(train_perc*class_size[c]);
		if (nsamples >= class_size[c])
			nsamples = class_size[c];
		if (nsamples <= 0){
			fprintf(stderr,"For class %d\n",c);
            iftError("No available samples", "iftSelectSupTrainSamplesByMST");
		}
		class_size[c]=nsamples;
		Z->ntrainsamples += nsamples;
	}

	/* Compute MST */

	mstree = iftCreateMST(Z);

	/* For each class, first select the samples with adjacent ones from
       distinct classes. */

	nsamples = 0;
	for (u=0; (u < mstree->nnodes)&&(nsamples < Z->ntrainsamples); u++) {
		s = mstree->node[u].sample;
		if ((iftHasSampleStatus(Z->sample[s], IFT_TEST))&&(class_size[Z->sample[s].truelabel]>0)) {
			for (adj = mstree->node[u].adj; adj != NULL; adj = adj->next){
				v = adj->elem;
				t = mstree->node[v].sample;
				if (Z->sample[s].truelabel != Z->sample[t].truelabel){
					iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
					class_size[Z->sample[s].truelabel]--;
					nsamples += 1;
					break;
				}
			}
		}
	}

	iftDestroyMST(&mstree);

	/* Second, select the remaining samples by following their
       increasing order of weight */

	if (nsamples < Z->ntrainsamples) {     /* Sort samples by weight */

		index  = iftAllocIntArray(Z->nsamples);
		weight = iftAllocFloatArray(Z->nsamples);

		for (s=0; s < Z->nsamples; s++) {
			index[s]  = s;
			weight[s] = Z->sample[s].weight;
		}

		iftFHeapSort(weight, index, Z->nsamples, IFT_INCREASING);

		/* Select the remaining samples */

		for (s=0; (s < Z->nsamples)&&(nsamples < Z->ntrainsamples); s++) {
			if ((iftHasSampleStatus(Z->sample[index[s]], IFT_TEST))&&
				(class_size[Z->sample[index[s]].truelabel]>0)){
				iftSetSampleStatus(&Z->sample[index[s]], IFT_TRAIN);
				class_size[Z->sample[index[s]].truelabel]--;
				nsamples += 1;
			}
		}

		iftFree(weight);
		iftFree(index);
	}

	iftFree(class_size);

	return(Z->ntrainsamples);
}

typedef struct iftSortMSTNode {
	int origIndex;
	float weight;
	iftMSTNode node;
} iftSortMSTNode;

static int cmpSortMSTNodeAscending(const void *p1, const void *p2) {
	const float weight1 = ((iftSortMSTNode *) p1)->weight;
	const float weight2 = ((iftSortMSTNode *) p2)->weight;
	if (weight1 < weight2) {
		return -1;
	} else if (weight1 > weight2) {
		return 1;
	}
	return 0;
}

static int cmpSortMSTNodeDescending(const void *p1, const void *p2) {
	const float weight1 = ((iftSortMSTNode *) p1)->weight;
	const float weight2 = ((iftSortMSTNode *) p2)->weight;
	if (weight1 < weight2) {
		return 1;
	} else if (weight1 > weight2) {
		return -1;
	}
	return 0;
}

void iftSortNodesByWeightMST(iftMST *mstree, int order) {
	// create an array of iftSortMSTNode to be sorted by qsort
	iftSortMSTNode *sortNodes = (iftSortMSTNode *) iftAlloc(mstree->nnodes, sizeof(iftSortMSTNode));
	if (sortNodes == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftSortNodesByWeightMST");
	}
	int i;
	for (i = 0; i < mstree->nnodes; i++) {
		sortNodes[i].origIndex = i;
		sortNodes[i].node = mstree->node[i];
		sortNodes[i].weight = mstree->Z->sample[mstree->node[i].sample].weight;
	}

	// sort nodes by weight
	if (order == IFT_DECREASING) {
		qsort(sortNodes, mstree->nnodes, sizeof(iftSortMSTNode), cmpSortMSTNodeDescending);
	} else {
		qsort(sortNodes, mstree->nnodes, sizeof(iftSortMSTNode), cmpSortMSTNodeAscending);
	}

	// create index conversion table
	int *orig2SortIdxTable = iftAllocIntArray(mstree->nnodes);
	for (i = 0; i < mstree->nnodes; i++) {
		orig2SortIdxTable[sortNodes[i].origIndex] = i;
	}

	for (i = 0; i < mstree->nnodes; i++) {
		// copy sorted node
		mstree->node[i] = sortNodes[i].node;
		// update sorted indices
		mstree->node[i].maxarcadj = orig2SortIdxTable[mstree->node[i].maxarcadj];
		iftSet *adj = mstree->node[i].adj;
		while (adj != NULL) {
			adj->elem = orig2SortIdxTable[adj->elem];
			adj = adj->next;
		}
	}

	iftFree(orig2SortIdxTable);
	iftFree(sortNodes);
}

iftSet *iftSelectSamplesForUserLabelingMST(iftMST *mstree, int n) {
	const iftSample *sample = mstree->Z->sample;
	iftSet *selectedNodes = NULL;
	int u, setSize = 0;
	for (u = 0; u < mstree->nnodes; u++) {
		const int s = mstree->node[u].sample;
		const int t = mstree->node[mstree->node[u].maxarcadj].sample;
		if (sample[s].label != sample[t].label) {
			if (sample[s].truelabel == IFT_NIL) { // if not labeled by the user
				setSize += iftUnionSetElem(&selectedNodes, s);
				if (setSize == n) break;
			}
			if (sample[t].truelabel == IFT_NIL) { // if not labeled by the user
				setSize += iftUnionSetElem(&selectedNodes, t);
				if (setSize == n) break;
			}
		}
	}
	return selectedNodes;
}


void iftSetNodeColors(iftMST *mst, iftSet *samples, char color){
	int u;
	for(u = 0; u < mst->nnodes; u++){
		if(iftSetHasElement(samples,mst->node[u].sample)){
			mst->node[u].color = color;
		}
	}
}

iftMST *iftCreateSuperpixelActiveLearningMST(iftRegionGraph* graph){
	iftMST *mstree = (iftMST*)iftAlloc(1,sizeof(iftMST));

	int size = graph->nnodes;
	iftDataSet *Z = graph->dataset;

	//Allocates space for mstree
	mstree->node = (iftMSTNode *)iftAlloc(size,sizeof(iftMSTNode));
	mstree->nnodes  = size;
	mstree->maxarcw = IFT_INFINITY_FLT_NEG;
	mstree->minarcw = IFT_INFINITY_FLT;
	mstree->Z       = graph->dataset;

	//Prepares the data structures
	int *pred = iftAllocIntArray(mstree->nnodes);
	float *cost = iftAllocFloatArray(mstree->nnodes);
	iftFHeap *Q = iftCreateFHeap(mstree->nnodes,cost);

	int i;
	for (i = 0; i < mstree->nnodes; i++) {
		mstree->node[i].adj = NULL;
		mstree->node[i].maxarcadj = IFT_NIL;
		mstree->node[i].color = IFT_WHITE;

		//Creates a correspondence between samples in the dataset and nodes in the MST
		mstree->node[i].sample = i;
		mstree->Z->sample[i].weight = 0.0;

		cost[i] = IFT_INFINITY_FLT;
		pred[i] = IFT_NIL;
	}

	for(i = 0; i < mstree->nnodes; i++){
		if(Q->color[i] == IFT_WHITE){
			cost[i] = 0;
			iftInsertFHeap(Q,i);
			mstree->maxarcw = 0.0;

			while ( !iftEmptyFHeap(Q) ) {
				//u is the current node
				int u = iftRemoveFHeap(Q);
				int v = pred[u];

				if (v != IFT_NIL){ /* Create the arcs of the MST  */
					float arcw = cost[u];
					iftInsertSet(&mstree->node[u].adj,v);
					iftInsertSet(&mstree->node[v].adj,u);
					if (arcw > mstree->maxarcw)
						mstree->maxarcw = arcw;
					if (arcw < mstree->minarcw)
						mstree->minarcw = arcw;
					if (arcw > Z->sample[u].weight)
						Z->sample[u].weight = arcw;
					if (arcw > Z->sample[v].weight)
						Z->sample[v].weight = arcw;
				}

				/* Propagate paths */

				iftSet *adj = graph->node[u].adjacent;
				while(adj){
					v = adj->elem;
					if (Z->sample[u].label != Z->sample[v].label){
						if(Q->color[v] != IFT_BLACK){
							float arcw = Z->iftArcWeight(Z->sample[u].feat,Z->sample[v].feat,Z->alpha,Z->nfeats);
							if ( arcw < cost[v] ) {
								pred[v] = u;
								cost[v] = arcw;
								if(Q->color[v] == IFT_WHITE)
									iftInsertFHeap(Q, v);
								else
									iftGoUpFHeap(Q, Q->pos[v]);
							}
						}
					}

					adj = adj->next;
				}
			}
		}
	}

	iftDestroyFHeap(&Q);

	/* update maxarcadj */
	int u;
	for (u= 0; u < mstree->nnodes; u++) {
		iftSet *set = mstree->node[u].adj;
		if (set != NULL) {
			mstree->node[u].maxarcadj = set->elem;
			set = set->next;
			while (set != NULL) {
				if (cost[mstree->node[u].maxarcadj] < cost[set->elem]) {
					mstree->node[u].maxarcadj = set->elem;
				}
				set = set->next;
			}
		}
	}

	iftFree(cost);
	iftFree(pred);

	return(mstree);
}
