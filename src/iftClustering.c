#include "iftClustering.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Json.h"
#include "ift/core/io/Stream.h"

/**
@file
@brief A description is missing here
*/
/* --------------------------- Private Methods ---------------------------- */

typedef struct ift_unsupfeatselprob {
    iftDataSet *Z;
    int kmax;
} iftUnsupFeatSelProb;

typedef struct ift_BestkProblem {
    iftKnnGraph *graph;
    iftKnnGraphCutFun iftGraphCutFun;
} iftBestkProblem;

iftUnsupFeatSelProb *iftCreateUnsupFeatSelProb(iftDataSet *Z, int kmax);
void                 iftDestroyUnsupFeatSelProb(iftUnsupFeatSelProb **prob);
void                 iftUnsupFeatSelecMSDeltas(iftMSPS *msps, float gamma);
float                iftUnsupFeatSelecMSPSFitness(void *problem, float *theta);
void                 iftBestkMSDeltas(iftMSPS *msps, float gamma);
float                iftBestkFitness(void *problem, float *theta);
int                  iftBestkByKnnGraphCutMSPS(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);
int                  iftRootNode(iftKnnGraph *graph, int u);

int iftUnsupOPFRootMap(iftKnnGraph *graph);
int iftUnsupTrainRootMap(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun, char fast_computation_k);

iftUnsupFeatSelProb *iftCreateUnsupFeatSelProb(iftDataSet *Z, int kmax)
{
    iftUnsupFeatSelProb *prob=(iftUnsupFeatSelProb *)iftAlloc(1,sizeof(iftUnsupFeatSelProb));
    prob->Z    = Z;
    prob->kmax = kmax;

    iftSelectUnsupTrainSamples(Z,1.0,100);

    return(prob);
}

void iftDestroyUnsupFeatSelProb(iftUnsupFeatSelProb **prob)
{
    if (*prob != NULL){
        iftFree(*prob);
        *prob     = NULL;
    }
}

void iftUnsupFeatSelecMSDeltas(iftMSPS *msps, float gamma)
{
    int i,j;

    for (i=0; i < msps->n; i++) {
        msps->delta->val[i] = 0.005;
        for (j=1; j < msps->m; j++) {
            msps->delta->val[iftGetMatrixIndex(msps->delta,i,j)] = msps->delta->val[iftGetMatrixIndex(msps->delta,i,j-1)]*gamma;
        }
    }
}

float iftUnsupFeatSelecMSPSFitness(void *problem, float *theta)
{
    float value=0.0;
    int i;
    iftKnnGraph *graph;
    iftUnsupFeatSelProb *prob= (iftUnsupFeatSelProb *) problem;
    iftDataSet *Z=prob->Z;

    /* Check the limits: constraints of the problem */

    for (i=0; i < Z->nfeats; i++) {
        if (theta[i] < 0.0) theta[i]=0.0;
        if (theta[i] > 1.0) theta[i]=1.0;
    }

    for (i=0; i < Z->nfeats; i++)
        prob->Z->alpha[i]=theta[i];

    /* Compute fitness value */

    graph = iftCreateKnnGraph(prob->Z,prob->kmax);
    iftUnsupTrain(graph, iftNormalizedCut);
    value = iftNormalizedCut(graph);
    iftDestroyKnnGraph(&graph);

    return(value);
}


/* Compute clusters on a KNN-Graph */

int iftUnsupOPF(iftKnnGraph *graph)
{
    int s,t,l,i,j,u,v;
    float tmp;
    iftAdjSet *adj;
    iftSet    *adjplat;
    iftDataSet *Z=graph->Z;

    // Initialization
    iftResetFHeap(graph->Q);
    iftSetRemovalPolicyFHeap(graph->Q, IFT_MAXVALUE);
    for(s = 0; s < graph->Z->nsamples; s++)
        iftRemoveSampleStatus(&Z->sample[s], IFT_PROTOTYPE);

    for (u = 0; u < graph->nnodes; u++){
        s = graph->node[u].sample;
        graph->pathval[u]     = Z->sample[s].weight-1.0;
        graph->node[u].root   = u;
        graph->node[u].pred   = IFT_NIL;
        graph->node[u].predWeight = 0;
        Z->sample[s].group    = 0;
        iftInsertFHeap(graph->Q, u);
    }

    // Optimum-Path Forest Computation
    l = 1; j = 0;
    while (!iftEmptyFHeap(graph->Q)){
        u=iftRemoveFHeap(graph->Q);
        graph->ordered_nodes[j]=u; j++;
        s = graph->node[u].sample;

        if (graph->node[u].root == u){ // root node
            graph->pathval[u]    = Z->sample[s].weight;
            graph->node[u].pred = IFT_NIL;
            graph->node[u].predWeight = Z->sample[s].weight;
            Z->sample[s].group = l;
            iftAddSampleStatus(&Z->sample[s], IFT_PROTOTYPE);
            l++;
        }

        // extend optimum paths
        for (adj=graph->node[u].adj, i=1; i <= graph->k; i++, adj = adj->next) {
            v = adj->node;
            t = graph->node[v].sample;

            if (graph->Q->color[v] != IFT_BLACK ) {
                tmp = iftMin(graph->pathval[u], Z->sample[t].weight);
                if (tmp > graph->pathval[v]){
                    graph->node[v].root  = graph->node[u].root;
                    graph->node[v].pred = u;
                    graph->node[v].predWeight = tmp;
                    Z->sample[t].group = Z->sample[s].group;
                    graph->pathval[v]    = tmp;
                    iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
                }
            }
        }

        adjplat = graph->node[u].adjplat;  // extend optimum paths on plateaus
        while (adjplat != NULL){
            v = adjplat->elem;
            t = graph->node[v].sample;
            if (graph->Q->color[v] != IFT_BLACK ) {
                tmp = iftMin(graph->pathval[u], Z->sample[t].weight);
                if (tmp > graph->pathval[v]){
                    graph->node[v].root  = graph->node[u].root;
                    graph->node[v].pred = u;
                    graph->node[v].predWeight = tmp;
                    Z->sample[t].group = Z->sample[s].group;
                    graph->pathval[v]    = tmp;
                    iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
                }
            }
            adjplat = adjplat->next;
        }
    }
    iftResetFHeap(graph->Q);
    Z->ngroups = l-1;

    /* mark prototypes as centroids */
    for(s = 0; s < graph->Z->nsamples; s++)
        if(iftHasSampleStatus(Z->sample[s], IFT_PROTOTYPE))
            iftAddSampleStatus(&Z->sample[s], IFT_CENTROID);
    
    return(Z->ngroups);
}

/* Compute displacements for best k by MSPS optimization */

void iftBestkMSDeltas(iftMSPS *msps, float gamma)
{
    int i,j;

    for (i=0; i < msps->n; i++) {
        msps->delta->val[i] = 1;
        //    printf("%lf ",msps->error->val[0][i]);
        for (j=1; j < msps->m; j++) {
            msps->delta->val[iftGetMatrixIndex(msps->delta,i,j)] = msps->delta->val[iftGetMatrixIndex(msps->delta,i,j-1)]*gamma;
            //      printf("%lf ",msps->error->val[j][i]);
        }
        //    printf("\n");
    }

}


/* MSPS fitness function for the best k problem */

float iftBestkFitness(void *problem, float *theta)
{
    float value=0.0;
    iftBestkProblem *prob = (iftBestkProblem *)problem;
    iftKnnGraph *graph = prob->graph;

    /* Check the limits: constraints of the problem */

    if (theta[0] <  1) theta[0]=1;
    if (theta[0] > graph->kmax) theta[0]=graph->kmax;

    /* Compute fitness value */

    graph->k = (int)theta[0];
    //  printf("For k=%d\n",graph->k);
    iftPDFByRange(graph);
    iftUnsupOPF(graph);
    value = prob->iftGraphCutFun(graph);
    //  printf("clustering done\n");
    //  printf("Fitness %lf for k %d\n",value,graph->k);

    return(value);
}

int iftBestkByKnnGraphCutMSPS(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun)
{
    iftMSPS *msps;
    iftBestkProblem *prob = (iftBestkProblem *)iftAlloc(1,sizeof(iftBestkProblem));

    prob->graph = graph;
    prob->iftGraphCutFun = iftGraphCutFun;
    msps  = iftCreateMSPS(1,5,iftBestkFitness,prob);
    iftBestkMSDeltas(msps,3.0);
    iftMSPSMin(msps);
    iftDestroyMSPS(&msps);
    iftFree(prob);

    return(graph->k);
}

int iftRootNode(iftKnnGraph *graph, int u)
{
    if (graph->node[u].root == u)
        return(u);
    else
        return(graph->node[u].root = iftRootNode(graph,graph->node[u].root));
}

/* Compute clusters on a KNN-Graph and propagates the root indices instead of the cluster labels*/
int iftUnsupOPFRootMap(iftKnnGraph *graph)
{
    int s,t,i,j,u,v;
    float tmp;
    iftAdjSet *adj;
    iftSet    *adjplat;
    iftDataSet *Z=graph->Z;

    // Initialization
    iftResetFHeap(graph->Q);
    iftSetRemovalPolicyFHeap(graph->Q, IFT_MAXVALUE);

    for (u = 0; u < graph->nnodes; u++){
        s = graph->node[u].sample;
        graph->pathval[u]     = Z->sample[s].weight-1.0;
        graph->node[u].root   = u;
        Z->sample[s].group    = -1;
        iftInsertFHeap(graph->Q, u);
    }

    // Optimum-Path Forest Computation
    int cluster_count=0;
    j = 0;
    while (!iftEmptyFHeap(graph->Q)){
        u=iftRemoveFHeap(graph->Q);
        graph->ordered_nodes[j]=u; j++;
        s = graph->node[u].sample;

        if (graph->node[u].root == u){ // root node
            graph->pathval[u]    = Z->sample[s].weight;
            /* in the sample label we are putting the root index instead of the cluster label */
            Z->sample[s].group   = s;
            cluster_count++;
        }

        // extend optimum paths
        for (adj=graph->node[u].adj, i=1; i <= graph->k; i++, adj = adj->next) {
            v = adj->node;
            t = graph->node[v].sample;

            if (graph->Q->color[v] != IFT_BLACK ) {
                tmp = iftMin(graph->pathval[u], Z->sample[t].weight);
                if (tmp > graph->pathval[v]){
                    graph->node[v].root  = graph->node[u].root;
                    Z->sample[t].group   = Z->sample[s].group;
                    graph->pathval[v]    = tmp;
                    iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
                }
            }
        }

        adjplat = graph->node[u].adjplat;  // extend optimum paths on plateaus
        while (adjplat != NULL){
            v = adjplat->elem;
            t = graph->node[v].sample;
            if (graph->Q->color[v] != IFT_BLACK ) {
                tmp = iftMin(graph->pathval[u], Z->sample[t].weight);
                if (tmp > graph->pathval[v]){
                    graph->node[v].root  = graph->node[u].root;
                    Z->sample[t].group   = Z->sample[s].group;
                    graph->pathval[v]    = tmp;
                    iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
                }
            }
            adjplat = adjplat->next;
        }

    }

    iftResetFHeap(graph->Q);

    Z->ngroups = cluster_count;

    return(Z->ngroups);
}

/*
 * Clusters the samples of a KNN Graph by the normal OPF method as done by the iftUnsupTrain function, but this function propagates the root index instead of the root label. This function is used in the divide and conquer methods. It returns the number of groups.
 */
int iftUnsupTrainRootMap(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun, char fast_computation_k)
{
    int ngroups;

    if (fast_computation_k)
        iftFastBestKByKnnGraphCut(graph,iftGraphCutFun);
    else
        iftBestkByKnnGraphCut(graph,iftGraphCutFun);

    iftPDFByRange(graph);
    ngroups=iftUnsupOPFRootMap(graph);

    return(ngroups);
}


/* --------------------- Public Methods --------------------------------------*/

void iftPDFByRange(iftKnnGraph *graph)
{ // NOLINT

    iftDataSet   *Z = graph->Z;
    float K         = 2.0*graph->maxarcw[graph->k]*graph->maxarcw[graph->k]/9.0;
    float maximum   = IFT_INFINITY_FLT_NEG, minimum= IFT_INFINITY_FLT;

    if (graph->kmax == 0)
        iftError("This is not a Knn-Graph\n", "iftPDFByRange");

    // Compute probability density function
    #pragma omp parallel for // NOLINT
    for (int u=0; u < graph->nnodes; u++) {
        int s = graph->node[u].sample;

        if (graph->node[u].adjplat != NULL)
            iftDestroySet(&(graph->node[u].adjplat));

        Z->sample[s].weight = 0;

        iftAdjSet *adj_u = NULL;
        int i;
        for (adj_u = graph->node[u].adj, i = 1; i <= graph->k; i++, adj_u = adj_u->next) {
            Z->sample[s].weight += (1.0 - expf((-adj_u->arcw*adj_u->arcw)/K));

            if (i == iftMax(graph->k / 2, 1))
	      graph->node[u].maxarcw = adj_u->arcw; // update radius for the maximum distance to the k closest neighbors
        }
    }

    for (int u=0; u < graph->nnodes; u++) {
        int s = graph->node[u].sample;
        if (Z->sample[s].weight > maximum)
            maximum = Z->sample[s].weight;
        if (Z->sample[s].weight < minimum)
            minimum = Z->sample[s].weight;
    }

    if (maximum > minimum){ /* it is mandatory to keep normalization */
        #pragma omp parallel for
        for (int u=0; u < graph->nnodes; u++){
            int s = graph->node[u].sample;
            Z->sample[s].weight = ((IFT_MAXWEIGHT - 1.0) * (maximum - Z->sample[s].weight) / (maximum - minimum)) + 1.0;
        }
    }


    // Add adjacent nodes on density plateaus if one is a neighbor of the other and the contrary doesn't happen
    #pragma omp parallel for shared(graph)
    for (int u=0; u < graph->nnodes; u++){
        iftAdjSet   *adj_u,*adj_v=NULL;
        int i,j;
        int s = graph->node[u].sample;
        for (adj_u=graph->node[u].adj,i=1;i <= graph->k;i++,adj_u=adj_u->next) {
            int v = adj_u->node;
            int t = graph->node[v].sample;

            if (iftAlmostZero(Z->sample[t].weight - Z->sample[s].weight))
            {
                for (adj_v=graph->node[v].adj,j=1;j <= graph->k; j++,adj_v=adj_v->next) {
                    if (u == adj_v->node) break;
                }
                if (j > graph->k) {
                    #pragma omp critical
                    iftInsertSet(&(graph->node[v].adjplat),u);
                }
            }
        }
    }
}



int iftBestkByKnnGraphCut(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun)
{
    int bestk=1;
    float gcm,min_gcm= IFT_INFINITY_FLT;

    /* for each value of K in the interval [1,KMAX]*/
    for (int k = graph->kmax; (k >= 1)&&(min_gcm != 0.0); k--) {
        graph->k = k;                   /* test the new value of K*/
        iftPDFByRange(graph);           /* compute the PDF with the new value of K*/
        iftUnsupOPF(graph);             /* cluster the samples*/
        gcm = iftGraphCutFun(graph);    /* evaluate the function cut and update the min value*/
        if (gcm < min_gcm){
            bestk   = k;
            min_gcm = gcm;
        }
//        printf("%d %.3f\n",k,gcm);
    }
    graph->k = bestk;                   /* get the value of K that minimizes the cut function*/

    return(graph->k);
}

int iftFastBestKByKnnGraphCut(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun)
{
    int bestk=1;
    float gcm,min_gcm= IFT_INFINITY_FLT;

    for (int k = graph->kmax; (k >= 1)&&(min_gcm != 0.0); k--) {
        graph->k = k;
        iftPDFByRange(graph);
        iftUnsupOPF(graph);
        gcm = iftGraphCutFun(graph);
        if (gcm < min_gcm){
            bestk   = k;
            min_gcm = gcm;
        }
        else{
            break;
        }
    }
//    printf("best k -> %d\n",bestk);
    graph->k = bestk;

    return(graph->k);
}

iftKnnGraph *iftCreateKnnGraph(iftDataSet *Z, int kmax)
 {

    iftKnnGraph *graph=(iftKnnGraph *)iftAlloc(1,sizeof(iftKnnGraph));
    int          nnodes=Z->ntrainsamples;
    iftAdjSet   *adj;

    if (nnodes == 0){
        iftError("No samples for training", "iftCreateKnnGraph");
    }
    if ((kmax >= nnodes)||(kmax < 0))
        iftError("Invalid number kmax of arcs %d", "iftCreateKnnGraph", kmax);

    if (kmax==0)
        kmax=1;

    graph->nnodes = nnodes;
    graph->node   = (iftKnnNode *)iftAlloc(nnodes,sizeof(iftKnnNode));
    if (graph->node == NULL){
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateKnnGraph");
    }

    graph->pathval       = iftAllocFloatArray(nnodes);
    graph->ordered_nodes = iftAllocIntArray(nnodes);
    graph->Q        = iftCreateFHeap(nnodes,graph->pathval);
    graph->maxarcw  = iftAllocFloatArray(kmax+1);
    graph->kmax     = kmax;
    graph->k        = kmax;
    graph->Z        = Z;
    
    for (int u=0; u < graph->nnodes; u++){
        graph->node[u].adj      = NULL;
        graph->node[u].adjplat  = NULL;
        graph->node[u].sample   = IFT_NIL;
        graph->node[u].maxarcw  = 0.0;
        graph->node[u].root     = u;
    }

    int j = 0;
    for (int s=0; s < Z->nsamples; s++){
      if (iftHasSampleStatus(Z->sample[s], IFT_TRAIN)){
	  graph->node[j].sample = s;
	  j++;
        }
    }


#pragma omp parallel for
    for (int u=0; u < graph->nnodes; u++){

        float *d=iftAllocFloatArray(kmax+2);
        int *nn=iftAllocIntArray(kmax+2);

        int s = graph->node[u].sample;

        for (int k=1; k <= kmax; k++){
            d[k]= IFT_INFINITY_FLT;
            nn[k]=-1;
        }

        if (!Z->iftArcWeight)
            iftError("Arc weight function must be defined", "iftCreateKnnGraph");

        // Compute the k closest training nodes of node u
        int  i;
        int k  = 2;
        for (int v=0; v < graph->nnodes; v++){
            if (u != v){
                int t = graph->node[v].sample;
                if (iftDist == NULL) {
                    d[k] = Z->iftArcWeight(Z->sample[s].feat, Z->sample[t].feat, Z->alpha, Z->nfeats);
                } else {
                    d[k] = iftDist->distance_table[s][t];
                }
                nn[k]  = v;
                i      = k;
                while ((i > 1) && (d[i]<d[i-1])){ // sort in the increasing
                    // order of distance
                    float dist   = d[i];
                    int n = nn[i];
                    d[i]   = d[i-1]; nn[i]   = nn[i-1];
                    d[i-1] = dist;   nn[i-1] = n;
                    i--;
                }
                if (k<=kmax) k++;
            }
        }


        /* Set an initial maximum arc weight for each node and insert k-nearest
     neighbors in the adjacency list, taking into account that insertion must keep their increasing order of distance.
        */
        graph->node[u].maxarcw = d[iftMax(kmax / 2, 1)]; // Median value makes it more robut to outliers. iftPDFByRange updates this parameter with the best maximum arc weight for each node.

        for (k=kmax; k >= 1; k--){ // Insertion in AdjSet is LIFO
            iftInsertAdjSet(&(graph->node[u].adj),nn[k],d[k]);
        }

        iftFree(d);
        iftFree(nn);
    }

    // Compute the maximum arc weight in the graph for each value of k.
    for (int k=1; k <= kmax; k++){
        graph->maxarcw[k] = IFT_INFINITY_FLT_NEG;
    }

    int k;
    for (int u=0; u < graph->nnodes; u++){
        for (adj = graph->node[u].adj,k=1; k <= kmax; k++, adj = adj->next){
            if (adj->arcw > graph->maxarcw[k]){
                graph->maxarcw[k] = adj->arcw;
            }
        }
    }

    return(graph);
}

iftKnnGraph *iftCreateKnnGraphInt1D(iftDataSet *Z, int kmax)
{
    iftKnnGraph *graph=(iftKnnGraph *)iftAlloc(1,sizeof(iftKnnGraph));
    int          nnodes=Z->ntrainsamples;
    iftAdjSet   *adj;

    if (Z->ref_data_type != IFT_REF_DATA_IMAGE){
        iftError("Reference Data Type must be an Integer Image", "iftCreateKnnGraphInt1D");
    }

    if (Z->nfeats > 1){
        iftError("Number of features must be 1", "iftCreateKnnGraphInt1D");
    }

    if (nnodes == 0){
        iftError("No samples for training", "iftCreateKnnGraphInt1D");
    }
    if ((kmax >= nnodes)||(kmax < 0))
        iftError("Invalid number kmax of arcs %d", "iftCreateKnnGraphInt1D", kmax);

    if (kmax==0)
        kmax=1;

    graph->nnodes = nnodes;
    graph->node   = (iftKnnNode *)iftAlloc(nnodes,sizeof(iftKnnNode));
    if (graph->node == NULL){
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateKnnGraphInt1D");
    }

    graph->pathval       = iftAllocFloatArray(nnodes);
    graph->ordered_nodes = iftAllocIntArray(nnodes);
    graph->Q        = iftCreateFHeap(nnodes,graph->pathval);
    graph->maxarcw  = iftAllocFloatArray(kmax+1);
    graph->kmax     = kmax;
    graph->k        = kmax;
    graph->Z        = Z;

    for (int u=0; u < graph->nnodes; u++){
        graph->node[u].adj      = NULL;
        graph->node[u].adjplat  = NULL;
        graph->node[u].sample   = IFT_NIL;
        graph->node[u].maxarcw  = 0.0;
        graph->node[u].root     = u;
    }

    int j = 0;
    for (int s=0; s < Z->nsamples; s++){
        if (iftHasSampleStatus(Z->sample[s], IFT_TRAIN)){
            graph->node[j].sample = s;
            j++;
        }
    }

    /* sort the nodes by their 1D value */
    int *ordered_value = iftAllocIntArray(Z->nsamples);
    int *ordered_index = iftAllocIntArray(Z->nsamples);

    for(int u=0; u < graph->nnodes; u++) {
      int s = graph->node[u].sample;
      ordered_value[u] = Z->sample[s].feat[0];
      ordered_index[u] = u;
    }
    iftBucketSort(ordered_value, ordered_index, graph->nnodes, IFT_INCREASING);

    /* choose the k closest for every sample (considering that the samples are now ordered by the 1D value) */
#pragma omp parallel for
    for(int u=0; u < graph->nnodes; u++) {
      int s = graph->node[u].sample;
      int val = Z->sample[s].feat[0];

      /* find the position of graph->node[u] inside the ordered array (binary search) */
      int low = 0;
      int high = graph->nnodes -1;
      int k=0;

      bool b=false;
      int mid;
      while(low <= high) {
        if(!b) {
          mid = (float)graph->nnodes/(float)val;
          b=true;
        }
        else {
          mid = (float)(low+high)/2.0;
        }

        if(val < ordered_value[mid])
          high = mid - 1;
        else if(val > ordered_value[mid])
          low = mid + 1;
        else  {
          k = mid;
          break;
        }
      }

      /* determine the position of the first and last samples for the k neighbors */
      int k1 = (k - ((float)kmax/2.0) <= 0) ? 0 : k - ((float)kmax/2.0);
      int k2 = k1 + (kmax-1);

      if(k2 > (graph->nnodes - 1)) {
        k1 -= (k2 - (graph->nnodes - 1));
        k2 = graph->nnodes - 1;
      }

      /* pick the k samples between k1 and k2 (in reverse order of proximity, because insertion in Adjet is LIFO) */
      for(int k=k1; k<=k2; k++) {
        iftInsertAdjSet(&(graph->node[u].adj), ordered_index[k], fabsf(ordered_value[k] - Z->sample[s].feat[0]));
      }

      /* Set an initial maximum arc weight for each node and insert k-nearest */
      graph->node[u].maxarcw = ordered_value[(int)iftMax((float)(k2-k1)/2.0, 1)]; // Median value makes it more robust to outliers. iftPDFByRange updates this parameter with the best maximum arc weight for each node.
    }

    /* Compute the maximum arc weight in the graph for each value of k */
    for (int k=1; k <= kmax; k++){
        graph->maxarcw[k] = IFT_INFINITY_FLT_NEG;
    }

    int k;
    for (int u=0; u < graph->nnodes; u++) {
        for (adj = graph->node[u].adj, k=1; k <= kmax; k++, adj = adj->next){
            if (adj->arcw > graph->maxarcw[k]){
                graph->maxarcw[k] = adj->arcw;
            }
        }
    }

    iftFree(ordered_value);
    iftFree(ordered_index);

    return(graph);
}

void iftDestroyKnnGraph(iftKnnGraph **graph)
{
    int u;
    iftKnnGraph *aux=(*graph);

    if (aux!=NULL){
        for (u=0; u < aux->nnodes; u++){
            if (aux->node[u].adj != NULL){
                iftDestroyAdjSet(&(aux->node[u].adj));
            }
            if (aux->node[u].adjplat != NULL){
                iftDestroySet(&(aux->node[u].adjplat));
            }
        }
        iftFree(aux->maxarcw);
        iftFree(aux->node);
        iftFree(aux->pathval);
        iftFree(aux->ordered_nodes);
        iftDestroyFHeap(&(aux->Q));
        iftFree(aux);
        (*graph) = NULL;
    }
}


iftKnnGraph *iftReadKnnGraph(const char *format, ...) {
    va_list args;
    char path[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(path, format, args);
    va_end(args);

    char *tmpdir = iftMakeTempDir("knngraph_", NULL, NULL);

    iftUnzipFile(path, tmpdir);

    char *info_path = iftJoinPathnames(2, tmpdir, "info.json");
    iftDict *info = iftReadJson(info_path);

    iftKnnGraph *graph = (iftKnnGraph*) iftAlloc(1, sizeof(iftKnnGraph));

    // reading the traiing dataset
    char *train_dataset_set = iftJoinPathnames(2, tmpdir, "dataset.zip");
    graph->Z = iftReadDataSet(train_dataset_set);
    iftFree(train_dataset_set);

    graph->nnodes = iftGetLongValFromDict("nnodes", info);
    graph->kmax = iftGetLongValFromDict("kmax", info);
    graph->k = iftGetLongValFromDict("k", info);

    // reading pathval
    iftDblArray *pathval = iftGetDblArrayFromDict("pathval", info);
    graph->pathval = iftAllocFloatArray(pathval->n);
    iftCopyDblArrayToFloatArray(graph->pathval, pathval->val, pathval->n);
    iftDestroyDblArray(&pathval);

    // reading maxarcw
    iftDblArray *maxarcw = iftGetDblArrayFromDict("maxarcw", info);
    graph->maxarcw = iftAllocFloatArray(maxarcw->n);
    iftCopyDblArrayToFloatArray(graph->maxarcw, maxarcw->val, maxarcw->n);
    iftDestroyDblArray(&maxarcw);

    // writing ordered_nodes
    iftIntArray *ordered_nodes = iftGetIntArrayFromDict("ordered_nodes", info);
    graph->ordered_nodes = ordered_nodes->val;
    ordered_nodes->val = NULL;
    iftDestroyIntArray(&ordered_nodes);


    graph->node = (iftKnnNode*) iftAlloc(graph->nnodes, sizeof(iftKnnNode));
    for (int i = 0; i < graph->nnodes; i++) {
        char node_key[512];
        sprintf(node_key, "nodes:%d", i);
        iftDict *node_info = iftGetDictFromDict(node_key, info);

        graph->node[i].root       = iftGetLongValFromDict("root", node_info);
        graph->node[i].sample     = iftGetLongValFromDict("sample", node_info);
        graph->node[i].pred       = iftGetLongValFromDict("pred", node_info);
        graph->node[i].maxarcw    = iftGetDblValFromDict("maxarcw", node_info);
        graph->node[i].predWeight = iftGetDblValFromDict("predWeight", node_info);


        iftIntArray *adjplat = iftGetIntArrayFromDict("adjplat", node_info);
        graph->node[i].adjplat = iftIntArrayToSet(adjplat);

        iftMatrix *adj = iftGetDblMatrixFromDict("adj", node_info);
        graph->node[i].adj = iftMatrixToAdjSet(adj);
    }

    iftFree(info_path);
    iftDestroyDict(&info);


    iftRemoveDir(tmpdir);
    iftFree(tmpdir);

    return graph;
}


void iftWriteKnnGraph(const iftKnnGraph *graph, const char *format, ...) {
    va_list args;
    char path[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(path, format, args);
    va_end(args);


    char *parent_dir = iftParentDir(path);
    if (!iftDirExists(parent_dir))
        iftMakeDir(parent_dir);
    iftFree(parent_dir);


    char *tmpdir = iftMakeTempDir("knngraph_", NULL, NULL);

    // writing the traiing dataset
    char *train_dataset_set = iftJoinPathnames(2, tmpdir, "dataset.zip");
    iftWriteDataSet(graph->Z, train_dataset_set);
    iftFree(train_dataset_set);


    iftDict *info = iftCreateDict();
    iftInsertIntoDict("nnodes", graph->nnodes, info);
    iftInsertIntoDict("k", graph->k, info);
    iftInsertIntoDict("kmax", graph->kmax, info);

    // writing pathval
    iftDblArray *pathval = iftCreateDblArray(graph->nnodes);
    iftCopyFloatArrayToDblArray(pathval->val, graph->pathval, graph->nnodes);
    iftInsertIntoDict("pathval", pathval, info);
    
    // writing maxarcw
    iftDblArray *maxarcw = iftCreateDblArray(graph->nnodes);
    iftCopyFloatArrayToDblArray(maxarcw->val, graph->maxarcw, graph->nnodes);
    // for some reason or another, there could be some nan in the maxarcw array
    for (int i = 0; i < maxarcw->n; i++) {
        if (isnan(maxarcw->val[i]))
            maxarcw->val[i] = 0.0;
    }
    iftInsertIntoDict("maxarcw", maxarcw, info);

    // writing ordered_nodes
    iftIntArray *ordered_nodes = iftCreateIntArray(graph->nnodes);
    iftCopyIntArray(ordered_nodes->val, graph->ordered_nodes, graph->nnodes);
    iftInsertIntoDict("ordered_nodes", ordered_nodes, info);

    for (int i = 0; i < graph->nnodes; i++) {
        iftKnnNode node = graph->node[i];

        iftDict *node_info = iftCreateDict();
        iftInsertIntoDict("maxarcw", node.maxarcw, node_info);
        iftInsertIntoDict("root", node.root, node_info);
        iftInsertIntoDict("sample", node.sample, node_info);
        iftInsertIntoDict("pred", node.pred, node_info);
        iftInsertIntoDict("predWeight", node.predWeight, node_info);

        iftIntArray *adjplat = iftSetToArray(node.adjplat);
        iftInsertIntoDict("adjplat", adjplat, node_info);

        iftMatrix *adj = iftAdjSetToMatrix(node.adj);
        iftInsertIntoDict("adj", adj, node_info);

        char node_key[512];
        sprintf(node_key, "nodes:%d", i);
        iftInsertIntoDict(node_key, node_info, info);
    }

    char *info_path = iftJoinPathnames(2, tmpdir, "info.json");
    iftWriteJson(info, info_path);
    iftFree(info_path);
    iftDestroyDict(&info);

    iftZipDirContent(tmpdir, path);

    iftRemoveDir(tmpdir);
    iftFree(tmpdir);
}


int iftUnsupTrain(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun)
{
    int ngroups;

    iftBestkByKnnGraphCut(graph, iftGraphCutFun);
    iftPDFByRange(graph);
    ngroups=iftUnsupOPF(graph);

    return(ngroups);
}

int iftFastUnsupTrain(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun)
{
    int ngroups;

    iftFastBestKByKnnGraphCut(graph, iftGraphCutFun);
    iftPDFByRange(graph);
    ngroups=iftUnsupOPF(graph);

    return(ngroups);
}

int iftUnsupTrainWithCClusters(iftKnnGraph *graph, int c)
{
    int ngroups=0;
    iftDataSet *Z=graph->Z;

    if (c >= graph->nnodes)
      c = 1;
    
    /* Find a number higher or equal to c clusters */

    for (graph->k=graph->kmax; (graph->k >= 1) && (ngroups < c); graph->k--){
        iftPDFByRange(graph);
        ngroups=iftUnsupOPF(graph);
    }
    
    if (ngroups > c) {
      
      /* Find the c most relevant roots */
      iftSet *R=NULL;
      int nroots=0;
      for (int i=0; (i < graph->nnodes)&&(nroots < c); i++){
	int u = graph->ordered_nodes[i];
	if (graph->node[u].root == u) {
	  iftInsertSet(&R,u);
	  nroots++;
	}
      }
      
      /* Make the remaining roots to point to their closest relevant root */
      for (int u = 0; u < graph->nnodes; u++) {
	if (graph->node[u].root == u){
	  float dist, min_dist = IFT_INFINITY_FLT;
	  int closest_root     = IFT_NIL;
	  iftSet *S=R;
	  int s = graph->node[u].sample;
	  while(S != NULL) {
	    int v = S->elem;
	    int t = graph->node[v].sample;
	    if (iftDist == NULL)
	      dist = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
	    else
	      dist = iftDist->distance_table[s][t];
	    if (dist < min_dist){
	      min_dist     = dist;
	      closest_root = v;
	    }
	    S = S->next;
	  }
	  if (closest_root != u)
	    iftRemoveSampleStatus(&Z->sample[graph->node[u].sample], \
				  IFT_PROTOTYPE);
	  graph->node[u].root = closest_root;
	}
      }
      
      /* Relabel the root set */
      
      int l=0;
      while (R != NULL) {
	int u = iftRemoveSet(&R);
	int s = graph->node[u].sample;
	l++; Z->sample[s].group = l; 
      }
      Z->ngroups = l;
      
      
      /* Propagate labels to the remaining nodes */
      
      for (int u = 0; u < graph->nnodes; u++) {
	int r = iftRootNode(graph,u);
	int root_sample = graph->node[r].sample;
	int s = graph->node[u].sample;
	Z->sample[s].group = Z->sample[root_sample].group;
      }	  
    }

    return(Z->ngroups);
}

int iftFastUnsupTrainWithCClustersBinarySearch(iftKnnGraph *graph, int c)
{
    int ngroups=0;
    iftDataSet *Z=graph->Z;

    if (c > graph->nnodes)
        iftError("Number of clusters is higher than the number of nodes", "iftUnsupTrainWithCClusters");

    /* Find a number higher or equal to c clusters by binary search */

    int k_left = 1, k_right = graph->kmax, k_mid = (float)(k_left + k_right) / 2.0;
    int ngroups_left, ngroups_right;
    while(k_left != k_right && ngroups < c) {
        graph->k = k_left;
        iftPDFByRange(graph);
        ngroups_left=iftUnsupOPF(graph);

        graph->k = k_right;
        iftPDFByRange(graph);
        ngroups_right=iftUnsupOPF(graph);

        /* we choose the side whose number of groups is closer (and greater) to c */
        if((ngroups_left - c) < (ngroups_right - c)) {
            k_right = k_mid - 1;
            ngroups = ngroups_left;
        }
        else {
            k_left = k_mid + 1;
            ngroups = ngroups_right;
        }

        k_mid = (float)(k_left + k_right) / 2.0;
    }

    if (ngroups > c) {

        /* Find the c most relevant roots */
        iftSet *R=NULL;
        int nroots=0;
        for (int i=0; (i < graph->nnodes)&&(nroots < c); i++){
            int u = graph->ordered_nodes[i];
            if (graph->node[u].root == u) {
                iftInsertSet(&R,u);
                nroots++;
            }
        }

        /* Make the remaining roots to point to their closest relevant root */
        for (int u = 0; u < graph->nnodes; u++) {
            if (graph->node[u].root == u){
                float dist, min_dist= IFT_INFINITY_FLT;
                int closest_root    = IFT_NIL;
                iftSet *S=R;
                int s = graph->node[u].sample;
                while(S != NULL) {
                    int v = S->elem;
                    int t = graph->node[v].sample;
                    if (iftDist == NULL)
                        dist = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
                    else
                        dist = iftDist->distance_table[s][t];
                    if (dist < min_dist){
                        min_dist     = dist;
                        closest_root = v;
                    }
                    S = S->next;
                }
                graph->node[u].root = closest_root;
            }
        }

        /* Relabel the root set */

        int l=0;
        while (R != NULL) {
            int u = iftRemoveSet(&R);
            int s = graph->node[u].sample;
            l++; Z->sample[s].group = l;
        }
        Z->ngroups = l;


        /* Propagate labels to the remaining nodes */

        for (int u = 0; u < graph->nnodes; u++) {
            int r = iftRootNode(graph,u);
            int root_sample = graph->node[r].sample;
            int s = graph->node[u].sample;
            Z->sample[s].group = Z->sample[root_sample].group;
        }
    }

    return(ngroups);
}

int iftUnsupClassify(iftKnnGraph *graph, iftDataSet *Z)
{

    iftDataSet *Z1=graph->Z;
    float K=2.0*graph->maxarcw[graph->k]*graph->maxarcw[graph->k]/9.0;

    if (Z1->nfeats != Z->nfeats)
        iftError("Incompatible datasets", "iftUnsupClassify");

    // Propagate group labels
    int noutliers=0;

#pragma omp parallel for shared(Z,Z1,graph,K,noutliers)
    for (int t=0; t < Z->nsamples; t++)
    {
      if ((Z==Z1)&&(iftHasSampleStatus(Z->sample[t], IFT_TRAIN) || iftHasSampleStatus(Z->sample[t], IFT_PROTOTYPE)))
            continue;

        /* Adan Echemendia modified this code, each label will be has default label -1, the previous code put this label
         * to 0 but this carry problems with methods of propagating the ids of the cluster representatives, so in this
         * methods a sample can be classified with label 0 because it was conquered by the representative with id 0 of a
         * cluster. This doesn't produce any consequence to the calls to this method */
        Z->sample[t].group   = -1;
        Z->sample[t].weight  = 0.0;
        float dist,min_dist = IFT_INFINITY_FLT;
        int closest_node = IFT_NIL;

        for (int i = 0; i < graph->nnodes; i++){
            int u    = graph->ordered_nodes[i];
            int s    = graph->node[u].sample;
            if (iftDist==NULL)
                dist = Z1->iftArcWeight(Z1->sample[s].feat,Z->sample[t].feat,Z1->alpha,Z1->nfeats);
            else
                dist = iftDist->distance_table[s][t];
            if (dist <= graph->node[u].maxarcw){
                Z->sample[t].group    = Z1->sample[s].group;
                Z->sample[t].weight   = graph->pathval[u]*expf(-dist/K);
                break;
            }else{ // compute closest node for the case of outliers
                if (dist < min_dist) {
                    min_dist       = dist;
                    closest_node   = u;
                }
            }
        }

        if (Z->sample[t].group==-1){ /* t is an outlier */
#pragma omp atomic
            noutliers++;
            iftAddSampleStatus(&Z->sample[t], IFT_OUTLIER);
            /* it is important for background removal */
            Z->sample[t].weight = graph->pathval[closest_node]*exp(-(min_dist)/K);
            int s = graph->node[closest_node].sample;
            Z->sample[t].group  = Z1->sample[s].group;
        }
    }

    for (int t=0; t < Z->nsamples; t++)
        if (Z->sample[t].group > Z->ngroups )
            Z->ngroups = Z->sample[t].group;

    return(noutliers);
}

int iftSupOPFKnnClassify(iftKnnGraph *graph, iftDataSet *Z)
{
  
  iftDataSet *Z1=graph->Z;
  float K=2.0*graph->maxarcw[graph->k]*graph->maxarcw[graph->k]/9.0;

  if (Z1->nclasses == 0)
      iftError("There are no classes", "iftSupKnnClassify");
  
  if (Z1->nfeats != Z->nfeats)
      iftError("Incompatible datasets", "iftSupKnnClassify");
  

  int noutliers=0;
  
#pragma omp parallel for shared(Z,Z1,graph,K,noutliers)
  for (int t=0; t < Z->nsamples; t++)
    {
      
      
      if ((Z==Z1)&&(iftHasSampleStatus(Z->sample[t], IFT_TRAIN))){
	continue;
      }
      
      Z->sample[t].label   = 0;
      Z->sample[t].weight  = 0.0;
      float dist,min_dist  = IFT_INFINITY_FLT;
      int closest_node     = IFT_NIL;

      for (int i = 0; i < graph->nnodes; i++){
	int u    = graph->ordered_nodes[i];
	int s    = graph->node[u].sample;
	if (iftDist==NULL)
	  dist = Z1->iftArcWeight(Z1->sample[s].feat,Z->sample[t].feat,Z1->alpha,Z1->nfeats);
	else
	  dist = iftDist->distance_table[s][t];

	if (dist <= graph->node[u].maxarcw){
	  Z->sample[t].label    = Z1->sample[s].truelabel;
	  Z->sample[t].weight   = graph->pathval[u]*expf(-dist/K);
	  break;
	}else{ // compute closest node for the case of outliers
	  if (dist < min_dist) {
	    min_dist       = dist;
	    closest_node   = u;
	  }
	}
      }

      if (Z->sample[t].label==0){ /* t is an outlier */
#pragma omp atomic
	noutliers++;
	iftAddSampleStatus(&Z->sample[t], IFT_OUTLIER);
	/* it is important for background removal */
	Z->sample[t].weight = graph->pathval[closest_node]*exp(-(min_dist)/K);
	int s = graph->node[closest_node].sample;
	Z->sample[t].label  = Z1->sample[s].truelabel;
      }
    }

    if (Z->nclasses < Z1->nclasses)
        Z->nclasses = Z1->nclasses;
  
    return(noutliers);
}



int iftUnsupClassifyDataSetLessOutliers(iftKnnGraph *graph, iftDataSet *Z)
{
    iftDataSet *Z1=graph->Z;
    float K=2.0*graph->maxarcw[graph->k]*graph->maxarcw[graph->k]/9.0;

    if (Z1->nfeats != Z->nfeats)
        iftError("Incompatible datasets", "iftUnsupClassify");

    // Propagate group labels
    int noutliers=0;

#pragma omp parallel for shared(Z,Z1,graph,K,noutliers)
    for (int t=0; t < Z->nsamples; t++)
    {
      if ((Z==Z1)&&(iftHasSampleStatus(Z->sample[t], IFT_TRAIN)))
            continue;

        Z->sample[t].group   = -1; /* see the previous methods to check this*/
        Z->sample[t].weight  = 0.0;
        float dist;

        for (int i = 0; i < graph->nnodes; i++){
            int u    = graph->ordered_nodes[i];
            int s    = graph->node[u].sample;
            if (iftDist==NULL)
                dist = Z1->iftArcWeight(Z1->sample[s].feat,Z->sample[t].feat,Z1->alpha,Z1->nfeats);
            else
                dist = iftDist->distance_table[s][t];
            if (dist <= graph->node[u].maxarcw){
                Z->sample[t].group    = Z1->sample[s].group;
                Z->sample[t].weight   = graph->pathval[u]*expf(-dist/K);
                break;
            }
        }

        if (Z->sample[t].group==-1){ /* t is an outlier */
#pragma omp atomic
            noutliers++;
            iftAddSampleStatus(&Z->sample[t], IFT_OUTLIER);
        }

    }

    /*TODO Check this, I think I can remove it*/
    for (int t=0; t < Z->nsamples; t++)
        if (Z->sample[t].group > Z->ngroups )
            Z->ngroups = Z->sample[t].group;

    return(noutliers);
}

/* In order to improve this function, we should provide higher
   priority to select outliers closer to the roots. This will make the
   algorithm to avoid real outliers, such as border voxels when Z
   comes from an image, in the training set.
 */

iftKnnGraph *iftUnsupLearn(iftDataSet *Z, float kmax_perc, iftKnnGraphCutFun iftGraphCutFun, int niterations, char debug)
{
    int kmax,  s, idx, noutliers, new_samples_from_test=0, new_samples_from_outliers=0;
    iftKnnGraph *graph=NULL;

    if (Z->ntrainsamples == 0)
        iftError("No samples for training", "iftUnsupLearn");

    noutliers = 0;
    idx = 0;

    do {

        if (noutliers > 0) {

            // Compute the number of new samples from outliers
            new_samples_from_outliers = noutliers * 0.05;

            if (new_samples_from_outliers==0) break;

            new_samples_from_test = iftMin(new_samples_from_outliers, Z->nsamples - Z->ntrainsamples - noutliers);

            /* Select new samples from test */
            if (debug)
                printf("Getting %d new samples from test set\n",new_samples_from_test);

            while(new_samples_from_test > 0) {
                do {
                    s = iftRandomInteger(0,Z->nsamples-1);
                } while (!iftHasSampleStatus(Z->sample[s], IFT_TEST));
                iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
                new_samples_from_test--;
                Z->ntrainsamples++;
            }

            /* Select new samples from outliers */
            if (debug)
                printf("Getting %d new samples from outlier set\n",new_samples_from_outliers);

            while(new_samples_from_outliers > 0) {
                do {
                    s = iftRandomInteger(0,Z->nsamples-1);
                } while (!iftHasSampleStatus(Z->sample[s], IFT_OUTLIER));
                iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
                new_samples_from_outliers--;
                Z->ntrainsamples++;
            }

            // Reset the remaining outliers to the IFT_TEST state and put prototypes to train status.
            for (s=0; s < Z->nsamples; s++) {
                if (iftHasSampleStatus(Z->sample[s], IFT_PROTOTYPE))
                    iftRemoveSampleStatus(&Z->sample[s], IFT_PROTOTYPE);
                else
                    if (iftHasSampleStatus(Z->sample[s], IFT_OUTLIER)){
                        iftRemoveSampleStatus(&Z->sample[s], IFT_OUTLIER);
                    }
            }
        }

        kmax  = (int) iftMax((kmax_perc * Z->ntrainsamples), 1);

        // Retraining

        if (graph != NULL)
            iftDestroyKnnGraph(&graph);
        graph = iftCreateKnnGraph(Z,kmax);
        iftUnsupTrain(graph, iftGraphCutFun);
        noutliers = iftUnsupClassifyDataSetLessOutliers(graph, Z);
        idx++;

        if (debug){
            fprintf(stdout,"Iter. %d: kmax %d, Train. samples %d, Outliers %d\n", idx, kmax, Z->ntrainsamples ,noutliers);fflush(stdout);
        }


    } while ((idx < niterations)&&(noutliers > Z->ntrainsamples * 0.05));

    if (debug)
    printf("The final number of outliers detected in the learning proccess was %d\n",noutliers);

    return(graph);
}

/* Trying to improve the origina iftUnsupLearn differentiating between real outliers and false positive outliers.
 */
iftKnnGraph *iftUnsupLearn2(iftDataSet *Z, float kmax_perc, iftKnnGraphCutFun iftGraphCutFun, int niterations, char
debug)
{
    iftKnnGraph *graph=NULL;
    int kmax;

    if (Z->ntrainsamples == 0)
        iftError("No samples for training", "iftUnsupLearn2");

    int noutliers = 0;
    int idx = 0;

    do {

        if (noutliers > 0) {

            /* creates a new dataset with the outliers of the previous iteration */
            iftDataSet *outlier_datset= iftCreateDataSet(noutliers, Z->nfeats);
            int *outliers_pos=iftAllocIntArray(noutliers);
            int q=0;

            for(int s=0;s<Z->nsamples;s++){
                if (iftHasSampleStatus(Z->sample[s], IFT_OUTLIER)){
                    iftCopySample(&Z->sample[s],&outlier_datset->sample[q],Z->nfeats,true);
                    outliers_pos[q]=s;
                    q++;
                }
            }
            if (debug){
                printf("Created a new outlier dataset with %d samples\n",noutliers);
            }

            /* train the outlier dataset */
            if (graph != NULL)
                iftDestroyKnnGraph(&graph);

            iftSetStatus(outlier_datset,IFT_TRAIN);
            kmax  = (int) iftMax((0.1 * outlier_datset->ntrainsamples), 1);
            if (debug)
                printf("The k-max of the outlier dataset is %d\n",kmax);

            graph = iftCreateKnnGraph(outlier_datset,kmax);
            iftUnsupTrain(graph, iftGraphCutFun);

            if (debug){
                printf("The outlier dataset was clustered in %d groups\n",outlier_datset->ngroups);
            }

            /* differentiate between the true positive outliers to the false positive outliers. We count the number of outliers that have each group, in a uniform distribution there must be approximately (total_outliers/cluster_number) outliers in each group so if a group has less than 0.5*(total_outliers/cluster_number) outliers we would consider all samples belonging to them as true
             * positive outliers and we wouldn't add them to the training set of Z. All remaining samples would be added as train samples of the original dataset Z. */
            int *count_samples_by_cluster=iftAllocIntArray(outlier_datset->ngroups);
            int l;
            for(int s=0;s<outlier_datset->nsamples;s++){
                l=outlier_datset->sample[s].group;
                count_samples_by_cluster[l-1]++;
            }

            char *is_group_outlier = iftAllocCharArray(outlier_datset->ngroups);
            for (int i=0;i<outlier_datset->ngroups;i++){
                if ((float)count_samples_by_cluster[i]< 0.5*(noutliers/outlier_datset->ngroups)){
                    is_group_outlier[i]=1;
                    if (debug)
                        printf("Cluster %d is composed by %d true positive outliers and aren't added to the training "
                                       "set\n",i,count_samples_by_cluster[i]);
                }
                else{
                    is_group_outlier[i]=0;
                    if (debug)
                        printf("Cluster %d with %d samples are false positive outliers and will be added to the "
                                       "training set\n",i,count_samples_by_cluster[i]);
                }
            }

            int orig_s;
            for (int s=0;s<outlier_datset->nsamples;s++){
                l=outlier_datset->sample[s].group;
                orig_s=outliers_pos[s];
                if (is_group_outlier[l-1])
                    /* this is a real outlier so it's not added to the training set*/
                    iftSetSampleStatus(&Z->sample[orig_s], IFT_TEST);
                else{
                    /* this is a false positive outlier so we added to the training set*/
                    iftSetSampleStatus(&Z->sample[orig_s], IFT_TRAIN);
                    Z->ntrainsamples++;
                }

            }

            iftDestroyDataSet(&outlier_datset);
            iftFree(outliers_pos);
            iftFree(is_group_outlier);
            iftFree(count_samples_by_cluster);

        }

        kmax  = (int) iftMax((kmax_perc * Z->ntrainsamples), 1);

        // Retraining
        if (graph != NULL)
            iftDestroyKnnGraph(&graph);
        graph = iftCreateKnnGraph(Z,kmax);

        iftFastUnsupTrain(graph, iftGraphCutFun);
        noutliers = iftUnsupClassify(graph, Z);
        idx++;

        if (debug)
            printf("Iter. %d: kmax %d, Train. samples %d, Outliers %d\n", idx, kmax, Z->ntrainsamples ,noutliers);

    } while ((idx < niterations)&&(noutliers > Z->ntrainsamples * 0.05));

    if (debug)
        printf("The final number of outliers detected in the learning proccess was %d\n",noutliers);

    return(graph);
}

int iftSemiSupClassify(iftKnnGraph *graph, iftDataSet *Z)
{
    int     u,s,t,i,closest_node,noutliers;
    float   dist, min_dist;
    iftDataSet *Z1=graph->Z;
    float K=2.0*graph->maxarcw[graph->k]*graph->maxarcw[graph->k]/9.0;


    if (!iftIsSupervisedDataSet(Z1))
        iftError("It requires a supervised training set", "iftSemiSupClassify");

    if (Z1->nfeats != Z->nfeats)
        iftError("Incompatible datasets", "iftSemiSupClassify");

    // Replace root labels by root classes in the graph

    for (u=0; u < graph->nnodes; u++){
        if (u == graph->node[u].root){
            s = graph->node[u].sample;
            Z1->sample[s].label = Z1->sample[s].truelabel;
        }
    }

    // Propagate root class in the graph

    for (u=0; u < graph->nnodes; u++)
        if (u != graph->node[u].root){
            s = graph->node[u].sample;
            t = graph->node[graph->node[u].root].sample;
            Z1->sample[s].label = Z1->sample[t].label;
        }


    // Propagate root classes to the remaining samples in the dataset

    noutliers=0;
    for (t=0; t < Z->nsamples; t++)
    {

        if ((Z==Z1)&&(iftHasSampleStatus(Z->sample[t], IFT_TRAIN))){
            continue;
        }

        Z->sample[t].label   = 0;
        Z->sample[t].weight  = 0.0;
        min_dist = IFT_INFINITY_FLT; closest_node = IFT_NIL;

        for (i = 0; i < graph->nnodes; i++){
            u    = graph->ordered_nodes[i];
            s    = graph->node[u].sample;

            if (iftDist==NULL)
                dist = Z1->iftArcWeight(Z1->sample[s].feat,Z->sample[t].feat,Z1->alpha,Z1->nfeats);
            else
                dist = iftDist->distance_table[s][t];

            if (dist <= graph->node[u].maxarcw){
                Z->sample[t].label    = Z1->sample[s].label;
                Z->sample[t].weight   = graph->pathval[u]*exp(-dist/K);
                break;
            }else{ // compute closest node for the case of outliers
                if (dist < min_dist) {
                    min_dist       = dist;
                    closest_node   = u;
                }
            }
        }

        if (Z->sample[t].label==0){ /* t is an outlier */
            noutliers++;
            iftSetSampleStatus(&Z->sample[t], IFT_OUTLIER);
            /* it is better to identify outliers */
            Z->sample[t].weight = 0;
            //	Z->sample[t].weight = graph->pathval[closest_node]*exp(-(min_dist)/K);
            s                   = graph->node[closest_node].sample;
            Z->sample[t].label  = Z1->sample[s].label;
        }
    }

    return(noutliers);
}

void iftUnsupFeatSelecByMSPS(iftDataSet *Z, int kmax)
{
    iftUnsupFeatSelProb *prob;
    iftMSPS *msps;
    int i;

    prob = iftCreateUnsupFeatSelProb(Z,kmax);
    msps = iftCreateMSPS(Z->nfeats,3,iftUnsupFeatSelecMSPSFitness,prob);
    for (i=0; i < Z->nfeats; i++)
        msps->theta[i]=Z->alpha[i];
    iftUnsupFeatSelecMSDeltas(msps,3.0);
    iftMSPSMin(msps);
    iftDestroyUnsupFeatSelProb(&prob);
    for (i=0; i < Z->nfeats; i++){
        if (msps->theta[i] < 0.0) msps->theta[i]=0.0;
        if (msps->theta[i] > 1.0) msps->theta[i]=1.0;
        Z->alpha[i]=msps->theta[i];
    }
    iftDestroyMSPS(&msps);
}


float iftNormalizedCut(iftKnnGraph *graph)
{
    float  normalized_cut;
    float *acumI; //acumulated complementar weights inside each group
    float *acumB; //acumulated complementar weights on the boundary between groups
    iftDataSet *Z=graph->Z;

    acumI = iftAllocFloatArray(Z->ngroups+1);
    acumB = iftAllocFloatArray(Z->ngroups+1);

#pragma omp parallel for shared(graph,acumI,acumB,Z)
    for (int u = 0; u < graph->nnodes; u++){
        iftAdjSet *adj;
        int s = graph->node[u].sample;
        int i;
        for (adj=graph->node[u].adj,i=1; i <= graph->k; i++,adj=adj->next) {
            int v = adj->node;
            int t = graph->node[v].sample;
            if (Z->sample[s].group == Z->sample[t].group){
#pragma omp atomic
                acumI[Z->sample[s].group] += (graph->maxarcw[graph->k]-adj->arcw);
            }else{
#pragma omp atomic
                acumB[Z->sample[s].group] += (graph->maxarcw[graph->k]-adj->arcw);
            }
        }
    }

    normalized_cut  = 0.0;
#pragma omp parallel for reduction(+:normalized_cut)
    for (int l=1; l <= Z->ngroups; l++){
        if ((acumI[l]+acumB[l])>0.0)
            normalized_cut +=  acumB[l]/(acumI[l]+acumB[l]);
    }

    iftFree(acumI);
    iftFree(acumB);

    return(normalized_cut);
}

iftDataSet *iftGraphBoundaryReduction(iftKnnGraph *graph, iftDataSet *Z) {
    const iftSample *sample = Z->sample;

    uchar *zSelection = iftAllocUCharArray(Z->nsamples);

    int s, t;

    if (graph->nnodes == Z->nsamples && graph->Z == Z) {
        int u;
        for (u = 0; u < graph->nnodes; u++) {
            iftAdjSet *adjSet = graph->node[u].adj;
            s = graph->node[u].sample;
            while (adjSet != NULL) {
                t = graph->node[adjSet->node].sample;
                if (sample[s].label != sample[t].label) {
                    zSelection[s] = 1;
                    zSelection[t] = 1;
                }
                adjSet = adjSet->next;
            }
        }
    } else {
        float maxarcw = 0.0f;
        for (s = 0; s < graph->nnodes; s++) {
            maxarcw += graph->node[s].maxarcw;
        }
        maxarcw /= graph->nnodes;

#pragma omp parallel for schedule(dynamic) shared(Z,sample,zSelection,maxarcw)
        for (s = 0; s < Z->nsamples; s++) {
            int t;
            for (t = s + 1; t < Z->nsamples; t++) {
                float weight;
                if (iftDist==NULL)
                    weight = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
                else
                    weight = iftDist->distance_table[s][t];

                if (weight <= maxarcw) {
                    if (sample[s].label != sample[t].label) {
                        zSelection[s] = 1;
                        zSelection[t] = 1;
                    }
                }
            }
        }
    }

    // count number of samples in the boundary
    int nsamples = 0;
    for (s = 0; s < Z->nsamples; s++) {
        if (zSelection[s] == 1) {
            nsamples++;
        }
    }

    iftDataSet *reducedZ = iftCreateDataSet(nsamples, Z->nfeats);
    reducedZ->nclasses = Z->nclasses;
    reducedZ->ngroups = Z->ngroups;

    for (s = 0, t = 0; s < Z->nsamples; s++) {
        if (zSelection[s] == 1) {
            float *feat = reducedZ->sample[t].feat;
            memcpy(feat, Z->sample[s].feat, sizeof(float) * Z->nfeats);
            reducedZ->sample[t] = Z->sample[s];
            reducedZ->sample[t].feat = feat;
            reducedZ->sample[t].truelabel = IFT_NIL;
            reducedZ->sample[t].label = IFT_NIL;
            t++;
        }
    }
    assert(t == nsamples);

    iftFree(zSelection);

    return reducedZ;
}

iftSet *iftGetKnnBoundarySamples(iftKnnGraph *graph) {
    iftDataSet *Z = graph->Z;
    const iftSample *sample = Z->sample;
    uchar *zSelection = iftAllocUCharArray(Z->nsamples);
    int s, t, u;
    iftSet *boundary=NULL;

    // Select boundary samples
    if (graph->nnodes == Z->nsamples && graph->Z == Z) {
        for (u = 0; u < graph->nnodes; u++) {
            iftAdjSet *adjSet = graph->node[u].adj;
            s = graph->node[u].sample;
            while (adjSet != NULL) {
                t = graph->node[adjSet->node].sample;
                if (sample[s].label != sample[t].label) {
                    zSelection[s] = 1;
                    zSelection[t] = 1;
                }
                adjSet = adjSet->next;
            }
        }
    } else {
        float maxarcw=0.0;
        for (u = 0; u < graph->nnodes; u++) {
            maxarcw += graph->node[u].maxarcw;
        }
        maxarcw /= graph->nnodes;

#pragma omp parallel for schedule(dynamic) shared(Z,sample,zSelection,maxarcw)
        for (s = 0; s < Z->nsamples; s++) {
            int t;
            for (t = s + 1; t < Z->nsamples; t++) {
                float weight;
                if (iftDist==NULL)
                    weight = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
                else
                    weight = iftDist->distance_table[s][t];
                if (weight <= maxarcw) {
                    if (sample[s].label != sample[t].label) {
                        zSelection[s] = 1;
                        zSelection[t] = 1;
                    }
                }
            }
        }
    }

    // Get boundary samples

    for (s = 0; s < Z->nsamples; s++) {
        if (zSelection[s] == 1) {
            iftInsertSet(&boundary,s);
        }
    }

    iftFree(zSelection);

    return boundary;
}

iftSet *iftGetKnnRootSamples(iftKnnGraph *graph) {
    iftSet *roots=NULL;
    int u;

    for(u = 0; u < graph->nnodes;u++){
        if( u == graph->node[u].root ){
            iftInsertSet(&roots,graph->node[u].sample);
        }
    }

    return roots;
}

/**
 * Compare two iftKnnNode by maxarcw in increasing order.
 * @return -1 if a<b, 1 if a>b, or 0 if a=b.
 */
static int iftCompareKnnNodeByMaxArcWInc(const void *a, const void *b)
{
    iftKnnNode *nodeA = (iftKnnNode *) a;
    iftKnnNode *nodeB = (iftKnnNode *) b;
    if (nodeA->maxarcw < nodeB->maxarcw) {
        return -1;
    } else if (nodeA->maxarcw > nodeB->maxarcw) {
        return 1;
    } else {
        return 0;
    }
}

/**
 * Return an array of size Z->ngroups+1 containing in each position a list of sample
 * indexes of each cluster, sorted by their distance to the root node (increasing distance).
 * Thus, the root node is the first element of each list.
 * This code assumes that each cluster contains only one root node and labels start at 1.
 * @param graph Input clustered graph.
 * @param nNodes Output array of size Z->ngroups+1 with the number of nodes in each cluster.
 * @return Array of size Z->ngroups+1 of iftFIFO*.
 */
int **iftGetNodesPerClusterOrderedByDistanceToRoot(iftKnnGraph *graph, int **nNodes)
{
    int i, j;
    const iftDataSet *Z = graph->Z;
    const iftKnnNode *node = graph->node;

    // sanity tests
    if (graph == NULL) {
        iftError("Input graph is null.", "iftGetNodesPerClusterOrderedByDistanceToRoot");
    }
    if (nNodes == NULL) {
        iftError("nNodes is null.", "iftGetNodesPerClusterOrderedByDistanceToRoot");
    }
    // check if each cluster have only one root node and count number of nodes
    int *nNode = (int *) iftAlloc(Z->ngroups + 1, sizeof(int));
    int *nRoot = (int *) iftAlloc(Z->ngroups + 1, sizeof(int));
    for (i = 0; i < graph->nnodes; i++) {
        const int sampleIdx = node[i].sample;
        const int sampleLabel = Z->sample[sampleIdx].label;
        nNode[sampleLabel]++;
        if (i == node[i].root) { // if node is root
            nRoot[sampleLabel]++;
        }
    }
    char msg[128];
    for (i = 1; i <= Z->ngroups; i++) {
        int error = 0;
        if (nNode[i] == 0) {
            sprintf(msg, "Cluster label %d is empty.", i);
            error = 1;
        }
        if (nRoot[i] == 0) {
            sprintf(msg, "Cluster label %d - no root node.", i);
            error = 1;
        }
        if (nRoot[i] > 1) {
            sprintf(msg, "Cluster label %d - %d root nodes - must have only one root node.",
                    i, nRoot[i]);
            error = 1;
        }
        if (error) {
            iftFree(nNode);
            iftFree(nRoot);
            iftError(msg, "iftGetNodesPerClusterOrderedByDistanceToRoot");
        }
    }
    iftFree(nRoot);

    // create an array of nodes for each cluster
    iftKnnNode **cluster = (iftKnnNode **) iftAlloc((Z->ngroups + 1), sizeof(iftKnnNode *));
    int *iCluster = (int *) iftAlloc((Z->ngroups + 1), sizeof(int));
    for (i = 1; i <= Z->ngroups; i++) {
        cluster[i] = (iftKnnNode *) iftAlloc(nNode[i], sizeof(iftKnnNode));
        iCluster[i] = 1; // index 0 for root node
    }
    for (i = 0; i < graph->nnodes; i++) {
        int sampleIdx = graph->node[i].sample;
        int sampleLabel = Z->sample[sampleIdx].label;
        if (i == graph->node[i].root) { // if node is root
            cluster[sampleLabel][0].sample = sampleIdx;
            cluster[sampleLabel][0].maxarcw = 0.0f;
        } else {
            cluster[sampleLabel][iCluster[sampleLabel]++].sample = sampleIdx;
        }
    }
    iftFree(iCluster);

    // compute distances
    for (i = 1; i <= Z->ngroups; i++) {
        // index of the root node of cluster i
        int rootIdx = cluster[i][0].sample;
        for (j = 1; j < nNode[  i]; j++) {
            // index of the node j of cluster i
            int nodeIdx = cluster[i][j].sample;
            // compute distance from root to node
            if (iftDist == NULL)
                cluster[i][j].maxarcw = Z->iftArcWeight(Z->sample[rootIdx].feat, Z->sample[nodeIdx].feat, Z->alpha, Z->nfeats);
            else
                cluster[i][j].maxarcw = iftDist->distance_table[rootIdx][nodeIdx];
        }

        // sort nodes of cluster i by increasing distance
        if (nNode[i] > 1) {
            // sort nodes 1 to nNode[i]
            qsort(&(cluster[i][1]), nNode[i] - 1, sizeof(iftKnnNode), iftCompareKnnNodeByMaxArcWInc);
        }
    }

    // convert array of nodes to array of indexes
    int **ret = (int **) iftAlloc(sizeof(int *), (Z->ngroups + 1));
    ret[0] = NULL;
    for (i = 1; i <= Z->ngroups; i++) {
        ret[i] = (int *) iftAlloc(sizeof(int), nNode[i]);
        for (j = 0; j < nNode[i]; j++) {
            ret[i][j] = cluster[i][j].sample;
        }
        iftFree(cluster[i]);
    }
    iftFree(cluster);

    *nNodes = nNode;

    return ret;
}

/**
 * Return an array of size Z->ngroups+1 containing in each position a list of sample
 * indexes of each cluster, sorted by their decreasing weight.
 * Thus, the root node is the first element of each list.
 * This code assumes that each cluster contains only one root node and labels start at 1.
 * @param graph Input clustered graph.
 * @param nNodes Output array of size Z->ngroups+1 with the number of nodes in each cluster.
 * @return Array of size Z->ngroups+1 of iftFIFO*.
 */
int **iftGetNodesPerClusterOrderedByWeight(iftKnnGraph *graph, int **nNodes)
{
    int i, j;
    const iftDataSet *Z = graph->Z;
    const iftKnnNode *node = graph->node;

    // sanity tests
    if (graph == NULL) {
        iftError("Input graph is null.", "iftGetNodesPerClusterOrderedByWeight");
    }
    if (nNodes == NULL) {
        iftError("nNodes is null.", "iftGetNodesPerClusterOrderedByWeight");
    }
    // check if each cluster have only one root node and count number of nodes
    int *nNode = (int *) iftAlloc(Z->ngroups + 1, sizeof(int));
    int *nRoot = (int *) iftAlloc(Z->ngroups + 1, sizeof(int));
    for (i = 0; i < graph->nnodes; i++) {
        const int sampleIdx = node[i].sample;
        const int sampleLabel = Z->sample[sampleIdx].label;
        nNode[sampleLabel]++;
        if (i == node[i].root) { // if node is root
            nRoot[sampleLabel]++;
        }
    }
    char msg[128];
    for (i = 1; i <= Z->ngroups; i++) {
        int error = 0;
        if (nNode[i] == 0) {
            sprintf(msg, "Cluster label %d is empty.", i);
            error = 1;
        }
        if (nRoot[i] == 0) {
            sprintf(msg, "Cluster label %d - no root node.", i);
            error = 1;
        }
        if (nRoot[i] > 1) {
            sprintf(msg, "Cluster label %d - %d root nodes - must have only one root node.",
                    i, nRoot[i]);
            error = 1;
        }
        if (error) {
            iftFree(nNode);
            iftFree(nRoot);
            iftError(msg, "iftGetNodesPerClusterOrderedByWeight");
        }
    }
    iftFree(nRoot);

    // create an array of nodes for each cluster
    iftKnnNode **cluster = (iftKnnNode **) iftAlloc((Z->ngroups + 1), sizeof(iftKnnNode *));
    int *iCluster = (int *) iftAlloc((Z->ngroups + 1), sizeof(int));
    for (i = 1; i <= Z->ngroups; i++) {
        cluster[i] = (iftKnnNode *) iftAlloc(nNode[i], sizeof(iftKnnNode));
        iCluster[i] = 1; // index 0 for root node
    }
    for (i = 0; i < graph->nnodes; i++) {
        int sampleIdx = graph->node[i].sample;
        int sampleLabel = Z->sample[sampleIdx].label;
        if (i == graph->node[i].root) { // if node is root
            cluster[sampleLabel][0].sample = sampleIdx;
            cluster[sampleLabel][0].maxarcw = 0.0f;
        } else {
            cluster[sampleLabel][iCluster[sampleLabel]++].sample = sampleIdx;
        }
    }
    iftFree(iCluster);

    // compute maxarcw of each node
    for (i = 1; i <= Z->ngroups; i++) {
        for (j = 1; j < nNode[i]; j++) {
            // index of the node j of cluster i
            int nodeIdx = cluster[i][j].sample;
            // assign the complement of node weight to sort in decreasing order
            cluster[i][j].maxarcw = IFT_MAXWEIGHT - Z->sample[nodeIdx].weight;
//      cluster[i][j].maxarcw = Z->sample[nodeIdx].weight;
        }

        // sort nodes of cluster i by increasing distance
        if (nNode[i] > 1) {
            // sort nodes 1 to nNode[i]
            qsort(&(cluster[i][1]), nNode[i] - 1, sizeof(iftKnnNode), iftCompareKnnNodeByMaxArcWInc);
        }
    }

    // convert array of nodes to array of indexes
    int **ret = (int **) iftAlloc(sizeof(int *), (Z->ngroups + 1));
    ret[0] = NULL;
    for (i = 1; i <= Z->ngroups; i++) {
        ret[i] = (int *) iftAlloc(sizeof(int), nNode[i]);
        for (j = 0; j < nNode[i]; j++) {
            ret[i][j] = cluster[i][j].sample;
        }
        iftFree(cluster[i]);
    }
    iftFree(cluster);

    *nNodes = nNode;

    return ret;
}

void iftPropagateClusterTrueLabels(iftKnnGraph *graph)
{
    iftDataSet *Z=graph->Z;

    if (Z->nclasses == 0)
      iftError("Dataset with undefined classes","iftPropagateClusterTrueLabels");
    
    
    for (int i=0; i < graph->nnodes; i++)
        if (graph->node[i].root == i){
            int s = graph->node[i].sample;
            Z->sample[s].label = Z->sample[s].truelabel;
        }

    for (int i=0;  i < graph->nnodes; i++){
        int s = graph->node[i].sample;
        int r = graph->node[i].root;
        Z->sample[s].label = Z->sample[graph->node[r].sample].label;
    }
}

void  iftPropagateClusterTrueLabels2(iftDataSet *Z){

    if (Z->nclasses < 1)
        iftError("The dataset must have the true labels associated", "iftPropagateClusterTrueLabels2");

    int *representatives = iftAllocIntArray(Z->ngroups+1);
    for (int s=0;s<Z->nsamples;s++){
        /* check if the sample is the representative of its group*/
        if (iftHasSampleStatus(Z->sample[s], IFT_PROTOTYPE)){
            /* save the representative position */
            representatives[Z->sample[s].group]=s;
            Z->sample[s].label=Z->sample[s].truelabel;
        }
    }

    /* Propagate the result for the other samples*/
    int group;
    for (int s=0;s<Z->nsamples;s++){
        if (!iftHasSampleStatus(Z->sample[s], IFT_PROTOTYPE)){
            group =Z->sample[s].group;
            Z->sample[s].label=Z->sample[representatives[group]].label;
        }
    }

    iftFree(representatives);
}

void  iftPropagateClusterTrueLabels3(iftDataSet *Z){
    if (Z->nclasses < 1)
        iftError("The dataset must have the true labels associated", "iftPropagateClusterTrueLabels3");

    int *matrix=iftAllocIntArray(Z->ngroups*Z->nclasses);

    for (int s=0;s<Z->nsamples;s++){
        matrix[(Z->sample[s].group-1)*Z->nclasses + (Z->sample[s].truelabel -1)]++;
    }

    int *new_labels=iftAllocIntArray(Z->ngroups+1);
    int best_true_label;
    int higher_count_true_label;

    for (int i=0;i<Z->ngroups;i++){
        best_true_label=-1;
        higher_count_true_label=0;
        for (int j=0;j<Z->nclasses;j++){
            if (matrix[i*Z->nclasses+j] > higher_count_true_label){
                best_true_label=j+1;
                higher_count_true_label=matrix[i*Z->nclasses+j];
            }
        }
        new_labels[i+1]=best_true_label;
    }

    /* Propagate the true label with the most occurrences to the rest of the cluster*/
#pragma omp parallel for
    for (int s=0;s<Z->nsamples;s++){
        Z->sample[s].label=new_labels[Z->sample[s].group];
    }

    iftFree(matrix);
    iftFree(new_labels);
}

bool iftSelectSupTrainSamplesWithNoOutliers(iftDataSet *Z, int max_ntrainsamples)
{
    /* limit the number of training samples by max_ntainsamples */

  if (Z->nsamples > max_ntrainsamples){
    iftSampler* sampler = iftRandomSubsampling(Z->nsamples, 1, max_ntrainsamples);
	  iftDataSetSampling(Z, sampler, 0);
    iftDestroySampler(&sampler);
  }
  else iftSetStatus(Z, IFT_TRAIN);

  /* compute the groups in a knn-graph with low value
     k=MAX(5,0.05max_ntrainsamples), and mark in the training set as
     outliers the samples whose truelabels are inconsistent with the
     truelabels of their groups. */

  //iftKnnGraph *kgraph = iftCreateKnnGraph(Z,1);
  iftKnnGraph *kgraph = iftCreateKnnGraph(Z,iftMax(5,0.05*Z->ntrainsamples));
  iftUnsupTrain(kgraph, iftNormalizedCut);
  for (int u=0; u < kgraph->nnodes; u++) {
    int s = kgraph->node[u].sample;
    int t = kgraph->node[kgraph->node[u].root].sample;
    if (Z->sample[s].truelabel != Z->sample[t].truelabel) {
      iftSetSampleStatus(&Z->sample[s], IFT_OUTLIER);
      Z->ntrainsamples--;
    }
  }

  /* the training set is considered invalid when there are no
     remaining training samples */

  if (Z->ntrainsamples<=0){
    iftDestroyKnnGraph(&kgraph);
    return(false);
  }

  /* the training set is considered invalid when there are no training
     samples for any given truelabel */

  int *nsamples_per_class = iftAllocIntArray(Z->nclasses);
  for (int s=0; s < Z->nsamples; s++)
    if (iftHasSampleStatus(Z->sample[s], IFT_TRAIN))
      nsamples_per_class[Z->sample[s].truelabel-1]++;

  for (int i=0; i < Z->nclasses; i++)
    if (nsamples_per_class[i]==0){
            iftDestroyKnnGraph(&kgraph);
            iftFree(nsamples_per_class);
            return(false);
    }

  iftDestroyKnnGraph(&kgraph);
  iftFree(nsamples_per_class);

  /* the training set is valid */

  return(true);
}

/* The divide and conquer solution takes the samples with a random strategy and only allows two levels. */
void iftDivideAndConquerUnsupOPF2Levels(iftDataSet *orig_dataset, int nb_partitions, int
nb_train_samples, iftKnnGraphCutFun iftGraphCutFun,
                                        float kmax1, float kmax2, char debug) {

    iftDataSet **first_level_datasets;

    int first_level_dataset_count= iftSplitDatasetIntoRandomlyDataSetArray(orig_dataset,nb_partitions,&first_level_datasets, false);

    if (debug){
        printf("They were created %d datasets\n\n",first_level_dataset_count);
    }

    /* the total number of clusters in the first level */
    int first_level_cluster_total=0;
    /* saves the number of clusters of each dataset of the first level */
    int *first_level_cluster_count= iftAllocIntArray(first_level_dataset_count);

    if (debug){
        printf("Number of samples in each dataset\n");
        for (int i=0;i<first_level_dataset_count;i++)
            printf("%d ",first_level_datasets[i]->nsamples);
        printf("\n\n");
    }

    if (nb_train_samples > 0){
        for (int i=0;i<first_level_dataset_count;i++)
            if (nb_train_samples > first_level_datasets[i]->nsamples)
                iftError("The number of train samples is greather than a partition in the dataset",
                         "iftDivideAndConquerUnsupOPF2Levels");
    }

    if (nb_train_samples > 0){
        /* cluster the datasets of the first level with the required number of training samples*/
#pragma omp parallel for reduction(+:first_level_cluster_total)
        for (int i=0;i<first_level_dataset_count;i++){

            iftSampler* sampler = iftRandomSubsampling(first_level_datasets[i]->nsamples, 1, nb_train_samples);
            iftDataSetSampling(first_level_datasets[i], sampler, 0);
            iftDestroySampler(&sampler);

            if (first_level_datasets[i]->ntrainsamples == 0)
                iftError("No samples for training in dataset %d", "iftDivideAndConquerUnsupOPF2Levels",
                         i);

            int kmax;
            if (kmax1 >= 1.0)
                kmax=kmax1;
            else
                kmax  = (int) iftMax((kmax1 * first_level_datasets[i]->ntrainsamples), 1);

            iftKnnGraph *graph = iftCreateKnnGraph(first_level_datasets[i],kmax);
            first_level_cluster_count[i]= iftUnsupTrainRootMap(graph, iftGraphCutFun,1);
            iftUnsupClassify(graph,first_level_datasets[i]);
            first_level_cluster_total+=first_level_cluster_count[i];

            /* destroys the KNN Graph */
            iftDestroyKnnGraph(&graph);
        }
    }
    else{
        /* cluster the datasets of the first level using all samples as training samples*/
#pragma omp parallel for reduction(+:first_level_cluster_total)
        for (int i=0;i<first_level_dataset_count;i++){

            iftSetStatus(first_level_datasets[i],IFT_TRAIN);

            if (first_level_datasets[i]->ntrainsamples == 0)
                iftError("No samples for training in dataset %d", "iftDivideAndConquerUnsupOPF2Levels",
                         i);

            int kmax;
            if (kmax1 >= 1.0)
                kmax=kmax1;
            else
                kmax  = (int) iftMax((kmax1 * first_level_datasets[i]->ntrainsamples), 1);

            iftKnnGraph *graph = iftCreateKnnGraph(first_level_datasets[i],kmax);
            first_level_cluster_count[i]= iftUnsupTrainRootMap(graph, iftGraphCutFun,1);
            first_level_cluster_total+=first_level_cluster_count[i];

            /* destroys the KNN Graph */
            iftDestroyKnnGraph(&graph);
        }
    }

    if (debug){
        printf("Total number of groups resulting in the first level -> %d\n\n",first_level_cluster_total);
        printf("Number of groups by partition in the first level\n");
        for (int i=0;i<first_level_dataset_count;i++){
            printf("%d ",first_level_cluster_count[i]);
        }
        printf("\n\n");
    }

    /* creates a new dataset with the representatives of the first level clusters. This will be the second level of the hierarchy*/
    iftDataSet *second_level_dataset= iftCreateDataSetWithoutFeatsVector(first_level_cluster_total, orig_dataset->nfeats);

    iftCopyRefData(second_level_dataset, orig_dataset->ref_data, orig_dataset->ref_data_type);
    // second_level_dataset->ref_data_type=orig_dataset->ref_data_type;
    // second_level_dataset->ref_data=orig_dataset->ref_data;
    second_level_dataset->nclasses=orig_dataset->nclasses;

#pragma omp parallel for
    for (int i=0;i< first_level_dataset_count;i++){
        int q=0;
        int initial_pos=0;
        for (int k=0;k < i;k++)
            initial_pos+=first_level_cluster_count[k];

        for (int s=0;s <first_level_datasets[i]->nsamples;s++){
            /* check if the sample is a cluster root in its dataset */
            if (s == first_level_datasets[i]->sample[s].group){
                iftCopySample(&first_level_datasets[i]->sample[s],&second_level_dataset->sample[initial_pos+q],
                              orig_dataset->nfeats,false);
                q++;
                /* check if we already found all cluster roots in the current dataset */
                if (q == first_level_cluster_count[i]){
                    break;
                }

            }
        }
    }
    /* Mark all samples as train samples in the second level */
    iftSetStatus(second_level_dataset,IFT_TRAIN);

    int kmax;
    if (kmax2 >= 1.0)
        kmax=kmax2;
    else
        kmax= (int) iftMax((kmax2 * second_level_dataset->ntrainsamples), 1);

    iftKnnGraph *graph=iftCreateKnnGraph(second_level_dataset,kmax);

    /* cluster the samples of the the dataset in the second level propagating the regularly cluster labels*/
//    iftUnsupTrain(graph, iftGraphCutFun);
    iftFastUnsupTrain(graph,iftGraphCutFun);

    orig_dataset->ngroups=second_level_dataset->ngroups;

    /* destroys the second level graph */
    iftDestroyKnnGraph(&graph);

    /* Propagates down the groups of the dataset of the second level */
#pragma omp parallel for
    for (int s=0;s<second_level_dataset->nsamples;s++){
        int sample_id=second_level_dataset->sample[s].id;
        orig_dataset->sample[sample_id].group=second_level_dataset->sample[s].group;
        iftSetSampleStatus(&orig_dataset->sample[sample_id], second_level_dataset->sample[s].status);
    }

    /* Destroys the second level dataset */
    iftDestroyDataSetWithoutFeatVectors(&second_level_dataset);

    /* Traverse the datasets of the first level of the hierarchy. For each dataset i, all samples of i contain as label the root index of the clustering result in this dataset. This roots already have its real label in the original dataset, because these labels were propagated from the second level of the hierarchy. So, for each sample s of dataset i, we find the label of its root in the original dataset and assign it as the label of s in the original dataset*/
#pragma omp parallel for
    for (int i=0;i< first_level_dataset_count;i++){
        for (int s=0;s< first_level_datasets[i]->nsamples;s++){
            int sample_id=first_level_datasets[i]->sample[s].id;
            int sample_root=first_level_datasets[i]->sample[s].group;
            int root_id=first_level_datasets[i]->sample[sample_root].id;
            orig_dataset->sample[sample_id].group=orig_dataset->sample[root_id].group;
            orig_dataset->sample[sample_id].weight=first_level_datasets[i]->sample[s].weight;
            if (!iftHasSampleStatus(orig_dataset->sample[sample_id], IFT_PROTOTYPE)){
                iftSetSampleStatus(&orig_dataset->sample[sample_id], first_level_datasets[i]->sample[s].status);
            }
        }
    }

    /* destroys all datasets of the first level of the hierarchy*/
    iftDestroyDataSetArray(&first_level_datasets,first_level_dataset_count,false);
    iftFree(first_level_datasets);
    iftFree(first_level_cluster_count);
}


iftKnnGraph *iftClusterImageByDivideAndConquerUnsupOPF2Levels(iftDataSet *orig_dataset, int num_blocks,
                                                              float train_percent,
                                                              iftKnnGraphCutFun iftGraphCutFun,
                                                              float kmax_first_level, float kmax_second_level,
                                                              char debug){

    if (orig_dataset->ref_data_type != IFT_REF_DATA_IMAGE && orig_dataset->ref_data_type != IFT_REF_DATA_MIMAGE && orig_dataset->ref_data_type != IFT_REF_DATA_FIMAGE){
        iftError("The referenced data type of the dataset is not an image ",
                 "iftClusterImageByDivideAndConquerUnsupOPF2Levels");
    }

    iftDataSet **first_level_datasets;
    iftMImage *mimg=(iftMImage *)orig_dataset->ref_data;

    int first_level_dataset_count=iftMImageTilesToDataSetArray(orig_dataset,num_blocks,&first_level_datasets,false);

    if (debug){
        printf("Total number of samples -> %d\n",orig_dataset->nsamples);
        printf("Count of datasets in first level -> %d\n",first_level_dataset_count);
    }

    /* the total number of clusters in the first level */
    int first_level_cluster_total=0;
    /* saves the number of clusters of each dataset of the first level */
    int *first_level_cluster_count= iftAllocIntArray(first_level_dataset_count);

#pragma omp parallel for reduction(+:first_level_cluster_total)
    for (int i=0;i<first_level_dataset_count;i++){

        iftBoundingBox *bb=(iftBoundingBox *)first_level_datasets[i]->ref_data;

        iftImage *mask1 = iftSelectRegionOfInterest(mimg->xsize,mimg->ysize,mimg->zsize,bb->begin,bb->end);
        iftImage *mask_sampling = iftGridSampling(mimg, mask1, (int)train_percent);

        int train_samples_nb=0;
#pragma omp parallel for reduction(+:train_samples_nb)
        for (int s=0;s<first_level_datasets[i]->nsamples;s++){
            int voxel=first_level_datasets[i]->sample[s].id;
            if (mask_sampling->val[voxel]){
                iftSetSampleStatus(&first_level_datasets[i]->sample[s], IFT_TRAIN);
                train_samples_nb++;
            }
        }
        first_level_datasets[i]->ntrainsamples=train_samples_nb;

        if (debug){
            printf("Founded %d train samples in dataset %d\n",train_samples_nb,i);
        }

        iftDestroyImage(&mask1);
        iftDestroyImage(&mask_sampling);

        if (first_level_datasets[i]->ntrainsamples == 0)
            iftError("No samples for training in dataset %d of level %d",
                     "iftClusterImageByDivideAndConquerUnsupOPF2Levels",
                     i, 0);

        int kmax;
        if (kmax_first_level < 1.0)
            kmax = (int) iftMax((kmax_first_level * first_level_datasets[i]->ntrainsamples), 1);
        else
            kmax=kmax_first_level;

        iftKnnGraph *graph = iftCreateKnnGraph(first_level_datasets[i],kmax);
        first_level_cluster_count[i]= iftUnsupTrainRootMap(graph, iftGraphCutFun, 1);
        iftUnsupClassify(graph, first_level_datasets[i]);

        first_level_cluster_total+=first_level_cluster_count[i];

        // iftFree(first_level_datasets[i]->ref_data);
        /* destroys the KNN Graph */
        iftDestroyKnnGraph(&graph);
    }

    /* creates a new dataset with the representatives of the first level clusters. This will be the second level of the hierarchy*/
    iftDataSet *second_level_dataset= iftCreateDataSetWithoutFeatsVector(first_level_cluster_total,
                                                                         orig_dataset->nfeats);

    if (debug){
        printf("Count of samples in the second level ->%d\n",second_level_dataset->nsamples);
    }

    iftCopyRefData(second_level_dataset, orig_dataset->ref_data, orig_dataset->ref_data_type);
    // second_level_dataset->ref_data_type=orig_dataset->ref_data_type;
    // second_level_dataset->ref_data=orig_dataset->ref_data;

#pragma omp parallel for
    for (int i=0;i< first_level_dataset_count;i++){
        int q=0;
        int initial_pos=0;
        for (int k=0;k < i;k++)
            initial_pos+=first_level_cluster_count[k];

        for (int s=0;s <first_level_datasets[i]->nsamples;s++){

            /* check if the sample is a cluster root in its dataset */
            if (s == first_level_datasets[i]->sample[s].group){
                iftCopySample(&first_level_datasets[i]->sample[s],&second_level_dataset->sample[initial_pos+q],
                              orig_dataset->nfeats,false);
                q++;
                /* check if we already found all cluster roots in the current dataset */
                if (q == first_level_cluster_count[i]){
                    break;
                }

            }
        }
    }

    /* Mark all samples as train samples in the second level */
    iftSetStatus(second_level_dataset,IFT_TRAIN);

    int kmax;
    if (kmax_second_level < 1.0)
        kmax= (int) iftMax((kmax_second_level * second_level_dataset->ntrainsamples), 1);
    else
        kmax=kmax_second_level;

    iftKnnGraph *graph=iftCreateKnnGraph(second_level_dataset,kmax);

    /* cluster the samples of the the second level dataset propagating the regularly cluster labels*/
//  iftFastUnsupTrain(graph, iftGraphCutFun);
    iftUnsupTrain(graph, iftGraphCutFun);

    if (debug){
        printf("The resulting k in the second level dataset is %d\n",graph->k);
        printf("The resulting number of groups in the second level dataet is %d\n",second_level_dataset->ngroups);
    }

    /* Propagates down the groups of the dataset of the second level */
#pragma omp parallel for
    for (int s=0;s<second_level_dataset->nsamples;s++){
        int sample_id=second_level_dataset->sample[s].id;
        orig_dataset->sample[sample_id].group=second_level_dataset->sample[s].group;
        iftSetSampleStatus(&orig_dataset->sample[sample_id], second_level_dataset->sample[s].status);
    }
    orig_dataset->ngroups=second_level_dataset->ngroups;

    if (debug){
        int count_prototypes=0;
        for (int s=0;s<orig_dataset->nsamples;s++)
            if (iftHasSampleStatus(orig_dataset->sample[s], IFT_PROTOTYPE))
                count_prototypes++;
        printf("Number of prototypes after second level propagation is %d\n",count_prototypes);
    }

    /* Traverse the datasets of the first level of the hierarchy. For each dataset i, all samples of i contain as group the root index of the clustering result in this dataset. This roots already have its real group in the original dataset, because these group were propagated from the second level of the hierarchy. So, for each sample s of dataset i, we find the group of its root in the original dataset and assign it as the group of s in the original dataset*/
#pragma omp parallel for
    for (int i=0;i< first_level_dataset_count;i++){
        for (int s=0;s< first_level_datasets[i]->nsamples;s++){
            int sample_id=first_level_datasets[i]->sample[s].id;
            int sample_root=first_level_datasets[i]->sample[s].group;
            int root_id=first_level_datasets[i]->sample[sample_root].id;
            orig_dataset->sample[sample_id].group=orig_dataset->sample[root_id].group;
            orig_dataset->sample[sample_id].weight=first_level_datasets[i]->sample[s].weight;
            /* signalize the outliers*/
            if (!iftHasSampleStatus(orig_dataset->sample[sample_id], IFT_PROTOTYPE)){
                iftSetSampleStatus(&orig_dataset->sample[sample_id], first_level_datasets[i]->sample[s].status);
            }
        }
    }

    if (debug){
        int count_prototypes=0;
        for (int s=0;s<orig_dataset->nsamples;s++)
            if (iftHasSampleStatus(orig_dataset->sample[s], IFT_PROTOTYPE))
                count_prototypes++;
        printf("Number of prototypes after first level propagation is %d\n",count_prototypes);
    }

    /* destroys all datasets of the first level of the hierarchy*/
    iftDestroyDataSetArray(&first_level_datasets,first_level_dataset_count,false);
    iftFree(first_level_datasets);
    iftFree(first_level_cluster_count);

    /* we are returning the knn graph of the second level, we need to destroy it in the caller function, we need also to destroy the related dataset that doesn't have proper feature vectors  */
    return graph;
}


/* This function assumes that the dataset contains an image. This is a divide and conquer approach but here I'm staying at the first level. It is designed to divide a image in few partitions because in each partition we are using some samples as train samples and we are classifying the rest of them */
int iftClusterImageByDivideAndConquerUnsupOPF1Level(iftDataSet *orig_dataset, int num_blocks, float train_percent,
                                                    iftKnnGraphCutFun iftGraphCutFun, float k_max, char debug) {

    if (orig_dataset->ref_data_type != IFT_REF_DATA_IMAGE && orig_dataset->ref_data_type != IFT_REF_DATA_MIMAGE && orig_dataset->ref_data_type != IFT_REF_DATA_FIMAGE){
        iftError("The referenced data type of the dataset is not an image ", "iftImageDivideAndConquerByBlocksUnsupOPF");
    }

    iftDataSet **first_level_datasets;
    iftMImage *mimg=(iftMImage *)orig_dataset->ref_data;

    int first_level_dataset_count=iftMImageTilesToDataSetArray(orig_dataset,num_blocks,&first_level_datasets,false);

    if (debug){
        printf("There were created %d datasets",first_level_dataset_count);
    }

    /* the total number of clusters in the first level */
    int first_level_cluster_total=0;
    /* saves the number of clusters of each dataset of the first level */
    int *first_level_cluster_count= iftAllocIntArray(first_level_dataset_count);

#pragma omp parallel for reduction(+:first_level_cluster_total)
    for (int i=0;i<first_level_dataset_count;i++){

        iftBoundingBox *bb=(iftBoundingBox *)first_level_datasets[i]->ref_data;

        iftImage *mask1 = iftSelectRegionOfInterest(mimg->xsize,mimg->ysize,mimg->zsize,bb->begin,bb->end);
        iftImage *mask_sampling = iftGridSampling(mimg, mask1, (int)train_percent);

        int train_samples_nb=0;
#pragma omp parallel for reduction(+:train_samples_nb)
        for (int s=0;s<first_level_datasets[i]->nsamples;s++){
            int voxel=first_level_datasets[i]->sample[s].id;
            if (mask_sampling->val[voxel]){
                iftSetSampleStatus(&first_level_datasets[i]->sample[s], IFT_TRAIN);
                train_samples_nb++;
            }
        }
        first_level_datasets[i]->ntrainsamples=train_samples_nb;

        iftDestroyImage(&mask1);
        iftDestroyImage(&mask_sampling);

        if (first_level_datasets[i]->ntrainsamples == 0)
            iftError("No samples for training in dataset %d of level %d",
                     "iftClusterImageByDivideAndConquerUnsupOPF1Level",
                     i, 0);
        int kmax;
        if (k_max < 1.0)
            kmax  = (int) iftMax((k_max * first_level_datasets[i]->ntrainsamples), 1);
        else
            kmax=k_max;

        iftKnnGraph *graph = iftCreateKnnGraph(first_level_datasets[i],kmax);
        iftFastUnsupTrain(graph,iftGraphCutFun);
        first_level_cluster_count[i]=first_level_datasets[i]->ngroups;
        first_level_cluster_total+=first_level_cluster_count[i];

        iftUnsupClassify(graph,first_level_datasets[i]);

        /* destroys the KNN Graph */
        // iftFree(first_level_datasets[i]->ref_data);
        iftDestroyKnnGraph(&graph);
    }

    /* Traverse the datasets of the first level and propagate their labels to the original dataset. Note that we can't mix labels between differents datasets*/
#pragma omp parallel for
    for (int i=0;i< first_level_dataset_count;i++){
        int initial_label=0;
        for(int k=0;k<i;k++)
            initial_label+=first_level_cluster_count[k];

        for (int s=0;s< first_level_datasets[i]->nsamples;s++){
            int sample_id=first_level_datasets[i]->sample[s].id;
            orig_dataset->sample[sample_id].group=first_level_datasets[i]->sample[s].group+initial_label;
            orig_dataset->sample[sample_id].weight=first_level_datasets[i]->sample[s].weight;
            /* signalize the outliers*/
            iftSetSampleStatus(&orig_dataset->sample[sample_id], first_level_datasets[i]->sample[s].status);
        }
    }
    orig_dataset->ngroups=first_level_cluster_total;

    /* destroys all datasets of the first level of the hierarchy*/
    iftDestroyDataSetArray(&first_level_datasets,first_level_dataset_count,false);
    iftFree(first_level_datasets);
    iftFree(first_level_cluster_count);

    return first_level_cluster_total;
}

iftDataSet* iftExtractCentroidsFromDataSetAsDataSet(iftDataSet *Z, bool returnRealCent, bool usePrototypes)
{
    if (Z->ngroups <= 0)
        iftError("The provided Dataset does not contain any groups", "iftExtractCentroidsFromDataSetAsDataSet");

    iftDataSet *cent = NULL;

    /* if the prototypes will be used as centroids, mark them as centroids and return them */
    if(usePrototypes) {
        cent = iftCreateDataSet(Z->ngroups, Z->nfeats);
        iftIntArray *protPerGroup = iftCreateIntArray(Z->ngroups);

        int prot = 0;
        for(int p=0; p<Z->nsamples; p++) {
            if(iftHasSampleStatus(Z->sample[p], IFT_PROTOTYPE)) {
                iftAddSampleStatus(&Z->sample[p], IFT_CENTROID);
                int g = Z->sample[p].group - 1;
                if(!protPerGroup->val[g]) {
                    iftCopySample(&Z->sample[p], &cent->sample[g], Z->nfeats, true);
                    protPerGroup->val[g] = 1;
                }
                else {
                    if(Z->sample[p].weight > cent->sample[g].weight)
                        iftCopySample(&Z->sample[p], &cent->sample[g], Z->nfeats, true);
                }
                prot++;
            }
        }

        if(prot == 0)
            iftError("The dataset does not contain any prototype", "iftExtractCentroidsFromDataSetAsDataSet");
    }
    /* if not, compute the centroids, mark them as prototypes and return them */
    else {
        float **centroids = (float**) iftAlloc(Z->ngroups+1, sizeof(float*));
        int *nsamples = iftAllocIntArray(Z->ngroups+1);
        int *minDistIdx = iftAllocIntArray(Z->ngroups+1);

        /* compute the mean vectors (centroids) for each cluster */
        for(int g=1; g<Z->ngroups+1; g++) {
            centroids[g] = iftAllocFloatArray(Z->nfeats);
            minDistIdx[g] = IFT_NIL;
        }

        for(int p=0; p<Z->nsamples; p++) {
            if(!iftHasSampleStatus(Z->sample[p], IFT_OUTLIER)) {
                iftSumFloatArrays(centroids[Z->sample[p].group], Z->sample[p].feat, centroids[Z->sample[p].group], Z->nfeats);
                nsamples[Z->sample[p].group]++;
            }
        }

        for(int g=1; g<Z->ngroups+1; g++) {
            for(int f=0; f<Z->nfeats; f++) {
                centroids[g][f] = centroids[g][f]/iftMax((float)nsamples[g],1.0);
            }
        }

        /* find the closest samples to the centroids */
        if(returnRealCent == 1) {
            for(int g=1; g<Z->ngroups+1; g++) {
                float minDist = IFT_INFINITY_FLT;
                for(int p=0; p<Z->nsamples; p++) {
                    if(!iftHasSampleStatus(Z->sample[p], IFT_OUTLIER)) {
                        float dist = Z->iftArcWeight(centroids[g],Z->sample[p].feat, Z->alpha, Z->nfeats);
                        if(dist < minDist) {
                            minDist = dist;
                            minDistIdx[g] = p;
                        }
                    }
                }
                iftAddSampleStatus(&Z->sample[minDistIdx[g]], IFT_CENTROID);
            }
        }

        /* create the dataset to be returned */
        cent = iftCreateDataSet(Z->ngroups, Z->nfeats);
        if(returnRealCent == 1) {
            for(int g = 1; g < Z->ngroups+1; g++)
                iftCopySample(&Z->sample[minDistIdx[g]], &cent->sample[g-1], Z->nfeats, true);
        }
        else {
            for(int g = 1; g < Z->ngroups+1; g++)
                for (int f = 0; f < Z->nfeats; f++)
                    cent->sample[g-1].feat[f] = centroids[g][f];
        }

        /* free memory */
        for(int g=1; g<Z->ngroups+1; g++) {
            iftFree(centroids[g]);
        }
        iftFree(centroids);
        iftFree(nsamples);
        iftFree(minDistIdx);
    }
    cent->nclasses = Z->nclasses;
    cent->ngroups = Z->ngroups;
    iftCopyRefData(cent, Z->ref_data, Z->ref_data_type);

    return cent;
}

iftMatrix* iftExtractCentroidsFromDataSetAsMatrix(iftDataSet *Z, bool returnRealCent, bool usePrototypes)
{
    if (Z->ngroups <= 0)
        iftError("The provided Dataset does not contain any groups", "iftExtractCentroidsFromDataSetAsMatrix");

    /* compute the centroids as a dataset and then convert them to a matrix */
    iftDataSet *centDataSet = iftExtractCentroidsFromDataSetAsDataSet(Z, returnRealCent, usePrototypes);
    iftMatrix *centMatrix = iftDataSetToFeatureMatrix(centDataSet);

    iftDestroyDataSet(&centDataSet);

    return centMatrix;
}

iftKnnGraph *iftMSTtoKnnGraph(iftMST* mst, int number_neigh)
{
    int kmax   = number_neigh;
    if (number_neigh >= mst->nnodes)
        iftError("Number of neighbors cannot be greater than or equal to the number of nodes", "iftMSTtoKnnGraph");

    float *cost     = iftAlloc(mst->nnodes, sizeof *cost);
    iftFHeap *H     = iftCreateFHeap(mst->nnodes, cost);
    iftDataSet *Z   = mst->Z;

    /* Initializing KNN Graph */

    iftKnnGraph *graph = iftAlloc(1, sizeof *graph);

    graph->nnodes        = mst->nnodes;
    graph->node          = iftAlloc(graph->nnodes, sizeof *graph->node);

    graph->pathval       = iftAllocFloatArray(graph->nnodes);
    graph->ordered_nodes = iftAllocIntArray(graph->nnodes);
    graph->Q             = iftCreateFHeap(graph->nnodes, graph->pathval);
    graph->maxarcw       = iftAllocFloatArray(kmax+1);
    graph->kmax          = kmax;
    graph->k             = kmax;
    graph->Z             = Z;

    if (graph->node == NULL || graph->pathval == NULL ||
        graph->ordered_nodes == NULL || graph->maxarcw == NULL)
    {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftMSTtoKnnGraph");
    }

    for (int u = 0; u < graph->nnodes; u++)
    {
        graph->node[u].adj     = NULL;
        graph->node[u].adjplat = NULL;
        graph->node[u].sample  = IFT_NIL;
        graph->node[u].maxarcw = 0.0f;
        graph->node[u].root    = u;
        graph->node[u].sample  = mst->node[u].sample;
    }

    for (int u = 0; u < mst->nnodes; u++)
    {
        cost[u]     = IFT_INFINITY_FLT;
    }

    for (int u = 0; u < mst->nnodes; u++)
    {
        iftSet *dirty = NULL;
        int     s     = mst->node[u].sample;

        for (iftSet *S = mst->node[u].adj; S != NULL; S = S->next)
        {
            int v = S->elem;
            int t = mst->node[v].sample;

            iftInsertSet(&dirty, v);
            cost[v]     = Z->iftArcWeight(Z->sample[s].feat, Z->sample[t].feat, Z->alpha, Z->nfeats);
            iftInsertFHeap(H, v);
        }

        iftSet *adj=NULL;

        int c = 0;
        while (c < number_neigh)
        {
            int v = iftRemoveFHeap(H);
            int t = mst->node[v].sample;

            iftInsertSet(&adj, v);

            for (iftSet *S = mst->node[v].adj; S != NULL; S = S->next)
            {
                int w = S->elem;
                int r = mst->node[w].sample;

                if (H->color[w] == IFT_WHITE && w != u)
                {
                    cost[w] = Z->iftArcWeight(Z->sample[t].feat, Z->sample[r].feat, Z->alpha, Z->nfeats);
                    iftInsertFHeap(H, w);
                    iftInsertSet(&dirty, w);
                }
            }

            c++;
        }

        for (int k = 0; k < graph->k; k++)
        {
            int v = iftRemoveSet(&adj);
            iftInsertAdjSet(&graph->node[u].adj, v, cost[v]);

            if (cost[v] > graph->node[u].maxarcw)
                graph->node[u].maxarcw = cost[v];

        }

        for (iftSet *S = dirty; S != NULL; S = S->next)
        {
            int v           = S->elem;
            cost[v]         = IFT_INFINITY_FLT;
            H->color[v]     = IFT_WHITE;
            H->pos[v]       = -1;
            H->node[v]      = -1;
        }
        H->last     = -1;

        iftDestroySet(&dirty);
    }

    /* compute maximum arc-weight for each k */

    for (int k = 1; k <= kmax; k++)
    {
        graph->maxarcw[k] = IFT_INFINITY_FLT_NEG;
    }

    for (int u = 0; u < graph->nnodes; u++)
    {
        int k = 1;
        for (iftAdjSet *adj = graph->node[u].adj; k <= kmax; adj = adj->next, k++)
        {
            if (adj->arcw > graph->maxarcw[k]) {
                graph->maxarcw[k] = adj->arcw;
            }
        }
    }

    iftFree(cost);
    iftDestroyFHeap(&H);

    return graph;
}

iftMatrix *iftKnnNodeToMatrix(iftKnnGraph *graph, int node_index)
{
    iftMatrix *A = NULL;
    iftAdjSet *adj = NULL;
    int nfeats = graph->Z->nfeats;
    int i,f;

    A = iftCreateMatrix(graph->Z->nfeats,graph->kmax);

    adj = graph->node[node_index].adj;

    for (i = 0; adj != NULL; i++){ //col
        int adj_index = adj->node;
        for (f = 0 ; f < nfeats; f++) //row
            A->val[iftGetMatrixIndex(A,f,i)] = graph->Z->sample[adj_index].feat[f];
        adj = adj->next;
    }

    return A;
}

iftIntArray *iftGetKnnNeighborSamplesIndex(iftKnnGraph *graph, int node_index)
{
    iftIntArray *arr = iftCreateIntArray(graph->kmax);
    iftAdjSet *adj = NULL;

    iftSetIntArray(arr->val,graph->kmax,-1);

    adj = graph->node[node_index].adj;
    for (int i = 0; adj != NULL; i++){
        int adj_index = adj->node;
        arr->val[i] = adj_index;
        adj = adj->next;
    }

    return arr;
}


//void iftJoinClusterByImageAdjacency(iftDataSet *Z, iftImage *comp, int k)
//{
//    int n_groups = Z->ngroups + 1;
//    iftMatrix *dist = iftCreateMatrix(n_groups, n_groups);
//
//    iftSetMatrix(dist, IFT_INFINITY_FLT);
//
////    iftAdjRel *A = iftCircular(1.0);
////
////    for (int i = 0; i < comp->n; i++)
////    {
////        int comp_i = comp->val[i];
////        int grp_i = Z->sample[comp_i].group;
////        iftVoxel u = iftGetVoxelCoord(comp, i);
////
////        for (int j = 1; j < A->n; j++)
////        {
////            iftVoxel v = iftGetAdjacentVoxel(A, u, j);
////            if (iftValidVoxel(comp, v))
////            {
////                int comp_j = comp->val[iftGetVoxelIndex(comp, v)];
////                int grp_j = Z->sample[comp_j].group;
////
////                float tmp = Z->iftArcWeight(Z->sample[comp_i].feat, Z->sample[comp_j].feat,
////                                            Z->alpha, Z->nfeats);
////
////                if (grp_i < grp_j && iftMatrixElem(dist, grp_j, grp_i) > tmp) {
////                    iftMatrixElem(dist, grp_j, grp_i) = tmp;
////                } else if (grp_j < grp_i && iftMatrixElem(dist, grp_i, grp_j) > tmp) {
////                    iftMatrixElem(dist, grp_i, grp_j) = tmp;
////                }
////            }
////        }
////
////    }
//
//    for (int i = 0; i < Z->nsamples - 1; i++)
//    {
//        int grp_i = Z->sample[i].group;
//        for (int j = i + 1; j < Z->nsamples; j++)
//        {
//            int grp_j = Z->sample[j].group;
//            float tmp = Z->iftArcWeight(Z->sample[i].feat, Z->sample[j].feat, Z->alpha, Z->nfeats);
//
//            if (grp_i < grp_j && iftMatrixElem(dist, grp_j, grp_i) > tmp) {
//                iftMatrixElem(dist, grp_j, grp_i) = tmp;
//            } else if (grp_j < grp_i && iftMatrixElem(dist, grp_i, grp_j) > tmp) {
//                iftMatrixElem(dist, grp_i, grp_j) = tmp;
//            }
//        }
//    }
//
//    iftFHeap *H = iftCreateFHeap(dist->n, dist->val);
//
//    for (int i = 0; i < dist->n; i++)
//    {
//        iftInsertFHeap(H, i);
//    }
//
//    int cur_k = n_groups;
//
//    /* k + 1 because the group 0 is empty and won't be removed */
//    while (cur_k > k + 1)
//    {
//        int p = iftRemoveFHeap(H);
//
//        int grp_i = iftGetMatrixRow(dist, p);
//        int grp_j = iftGetMatrixCol(dist, p);
//
//        for (int i = 0; i < Z->nsamples; i++)
//        {
//            if (Z->sample[i].group == grp_j) {
//                Z->sample[i].group = grp_i;
//            }
//        }
//
//        for (int i = 0; i < n_groups; i++)
//        {
//            if (grp_j < i) {
//
//                if (iftMatrixElem(dist, i, grp_j) < iftMatrixElem(dist, i, grp_i)) {
//                    iftMatrixElem(dist, i, grp_i) = iftMatrixElem(dist, i, grp_j);
//                    iftGoUpFHeap(H, H->pos[iftGetMatrixIndex(dist, i, grp_i)]);
//                }
//                iftMatrixElem(dist, i, grp_j) = IFT_INFINITY_FLT;
//                iftGoDownFHeap(H, H->pos[iftGetMatrixIndex(dist, i, grp_j)]);
//
//            } else if (grp_j > i) {
//                if (grp_i < i) {
//                    if (iftMatrixElem(dist, grp_j, i) < iftMatrixElem(dist, i, grp_i)) {
//                        iftMatrixElem(dist, i, grp_i) = iftMatrixElem(dist, grp_j, i);
//                        iftGoUpFHeap(H, H->pos[iftGetMatrixIndex(dist, i, grp_i)]);
//                    }
//                } else {
//                    if (iftMatrixElem(dist, grp_j, i) < iftMatrixElem(dist, grp_i, i)) {
//                        iftMatrixElem(dist, grp_i, i) = iftMatrixElem(dist, grp_j, i);
//                        iftGoUpFHeap(H, H->pos[iftGetMatrixIndex(dist, grp_i, i)]);
//                    }
//                }
//                iftMatrixElem(dist, grp_j, i) = IFT_INFINITY_FLT;
//                iftGoDownFHeap(H, H->pos[iftGetMatrixIndex(dist, grp_j, i)]);
//            }
//        }
//
//        --cur_k;
//    }
//
////    iftDestroyFHeap(&H);
////    iftDestroyMatrix(&dist);
//}

iftFloatArray *iftClusterPurity(iftKnnGraph   *graph, bool compute_OPF)
{
  iftDataSet *Z = graph->Z;
  
  if (Z->nclasses == 0)
    iftError("It requires an annotated dataset","iftClusterPurity");

  if (Z->ngroups == 0)
    compute_OPF=true;
  
  iftFloatArray *purity   =  iftCreateFloatArray(Z->nclasses+1);
  int           *nsamples =  iftCountSamplesPerClassDataSet(Z);
						     
  if (compute_OPF){
    iftPDFByRange(graph);
    iftUnsupOPF(graph);
  }
    
  for (int u=0; u < graph->nnodes; u++){
    int s = graph->node[u].sample;
    int r = graph->node[graph->node[u].root].sample;
    if (Z->sample[s].truelabel == Z->sample[r].truelabel){
      purity->val[Z->sample[r].truelabel]++;
    }
  }

  for (int c=1; c <= Z->nclasses; c++){
    purity->val[c] /= nsamples[c];
    purity->val[0] += purity->val[c];
  }
  purity->val[0] /= Z->nclasses;
  
  return(purity);
}
