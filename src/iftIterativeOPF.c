#include "ift.h"
#include "iftIterativeOPF.h"

void iftIterativeOPF(iftGraph *graph, int n_clusters, int max_iter, iftConnFuncType conn_func_type, iftCentUpdateType cent_update_type) {

  if(graph->Z == NULL) {
    iftError("Dataset is empty", "iftIterativeOPF");
  }
    
  if (graph->Z->nsamples == 0) {
    iftError("No samples for clustering", "iftIterativeOPF");
  }

  if (cent_update_type == IFT_XY_COORD && graph->Z->ref_data_type != IFT_REF_DATA_MIMAGE) {
    iftError("Center update type IFT_IMG_COORD can't be used with given data type", "iftIterativeOPF");
  }
  
  graph->Z->ngroups = n_clusters;

  if (graph->centroids == NULL)
    graph->centroids = iftRandomIntegers(0, graph->nnodes - 1, n_clusters);

  graph->center_hist    = iftAlloc(1, sizeof(*graph->center_hist));
  graph->center_hist->n = 0;
  
  int iter           = 0;

  iftDynamicSet **S = iftAlloc(n_clusters, sizeof (*S));

  while (!iftCheckConvergence(graph->centroids, graph->center_hist->val, n_clusters, graph->center_hist->n) && iter < max_iter) {
    // Inserting current centroid to centroid history
    iftInsertCenterHist(graph->center_hist, graph->centroids, n_clusters);

    iftDynamicOptimumPathForest(graph, S, conn_func_type);

    if (cent_update_type == IFT_DATA_FEAT) {
      iftUpdateCentersByFeatures(graph, S, graph->centroids, n_clusters);
    } else if (cent_update_type == IFT_XY_COORD) {
      iftUpdateCentersByXYCoord(graph, S, graph->centroids, n_clusters);
    } else {
      iftError("Unknown option for center update type", "iftIteratedWatersheds");
    }

    iftResetFHeap(graph->Q); 

    for (int c = 0; c < n_clusters; c++)
      iftDestroyDynamicSet(&S[c]);
  
    iter += 1;
  }

  iftFree(S);
}


void iftIteratedWatersheds(iftGraph *graph, int n_clusters, int max_iter, iftConnFuncType conn_func_type, iftCentUpdateType cent_update_type, bool use_arcw) {
  if(graph->Z == NULL) {
    iftError("Dataset is empty", "iftIteratedWatersheds");
  }
    
  if (graph->Z->nsamples == 0) {
    iftWarning("No samples for clustering", "iftIteratedWatersheds");
    return;
  }

  if (cent_update_type == IFT_XY_COORD && graph->Z->ref_data_type != IFT_REF_DATA_MIMAGE) {
    iftError("Center update type IFT_XY_COORD can't be used with given data type", "iftIterativeOPF");
  }

  graph->Z->ngroups = n_clusters;
  if (graph->centroids == NULL)
    graph->centroids = iftRandomIntegers(0, graph->nnodes - 1, n_clusters);

  graph->center_hist    = iftAlloc(1, sizeof(*graph->center_hist));
  graph->center_hist->n = 0;

  int iter = 0;
  
  iftDynamicSet **S = iftAlloc(n_clusters, sizeof (*S));
  
  while (!iftCheckConvergence(graph->centroids, graph->center_hist->val, n_clusters, graph->center_hist->n) && iter < max_iter) {
    // Inserting current centroid to centroid history
    iftInsertCenterHist(graph->center_hist, graph->centroids, n_clusters);
    
    iftOptimumPathForest(graph, S, conn_func_type, use_arcw);
    
    if (cent_update_type == IFT_DATA_FEAT) {
      iftUpdateCentersByFeatures(graph, S, graph->centroids, n_clusters);
    } else if (cent_update_type == IFT_XY_COORD) {
      iftUpdateCentersByXYCoord(graph, S, graph->centroids, n_clusters);
    } else {
      iftError("Unknown option for center update type", "iftIteratedWatersheds");
    }
    iftResetFHeap(graph->Q);
    
    for (int c = 0; c < n_clusters; c++)
      iftDestroyDynamicSet(&S[c]);

    iter += 1;
  }

  iftFree(S);
}


void iftInsertCenterHist(iftCenterHist *center_hist, int *centroids, int n_clusters) {
  center_hist->val = iftRealloc(center_hist->val, (center_hist->n + 1) * sizeof(*center_hist->val));
  center_hist->val[center_hist->n] = iftAllocIntArray(n_clusters);
  iftCopyIntArray(center_hist->val[center_hist->n], centroids, n_clusters);
  center_hist->n += 1;
}


void iftPrintCenterHist(iftCenterHist *center_hist, int n_clusters) {
  for (int i = 0; i < center_hist->n; i++) {
    iftPrintIntArray(center_hist->val[i], n_clusters);
  }
}


void iftInsertSetDynamicSet(iftGraph *graph, iftDynamicSet *S, int u) {
    if (S->size) {
        iftSet *a = (iftSet*) iftAlloc(1, sizeof *a);
        if (!a)
            iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertDynamicSet");
        a->elem = u;
        S->end->next = a;
        S->end = a;
    } else {
        S->begin = (iftSet*) iftAlloc(1, sizeof *S->begin);
        S->begin->elem = u;
        S->end = S->begin;
    }
    
    S->size += 1;
    
    int s = graph->node[u].sample;

    for (int i = 0; i < S->dim; i++) { 
        S->mean[i] += (graph->Z->sample[s].feat[i] - S->mean[i]) / S->size;
    }
}


void iftUpdateCentersByFeatures(iftGraph *graph, iftDynamicSet **S, int *centroids, int n_clusters) {
  float min_dist, dist;
  int min_root, u, s;
  iftSet *aux = NULL;
  float *mean = iftAlloc(graph->Z->nfeats, sizeof(*mean));

  for (int i = 0; i < n_clusters; i++) {
    min_dist = IFT_INFINITY_FLT;
    aux = S[i]->begin;
    min_root = IFT_NIL;
    while (aux != NULL) {
      u = aux->elem;
      s = graph->node[u].sample;
      //iftCopyDblArrayToFloatArray(mean, S[i]->mean, graph->Z->nfeats);
      //dist = iftEuclideanDistance(mean, graph->Z->sample[s].feat, graph->Z->nfeats);
      dist = iftEuclideanDistanceDblFlt(S[i]->mean, graph->Z->sample[s].feat, graph->Z->nfeats);
      if (dist < min_dist) {
          min_dist = dist;
          min_root = u;
      }
      aux = aux->next;
    }
	  if (min_root != IFT_NIL)
	    centroids[i] = min_root;
  }  

  iftFree(mean);
}


void iftUpdateCentersByXYCoord(iftGraph *graph, iftDynamicSet **S, int *centroids, int k) {
  float min_dist, dist;
  int min_root, size, s;
  iftSet *aux = NULL;
  iftMImage *mimg = graph->Z->ref_data;
  float *coord = iftAllocFloatArray(3);
  iftVoxel u;
  
  // computing the means
  float** means = iftAlloc(k, sizeof(float *));
  for (int i = 0; i < k; i++) {
    means[i] = iftAllocFloatArray(2);
    aux      = S[i]->begin;
    size     = 0;
    // computing the individual means
    while (aux != NULL) {
      s = aux->elem;
      u = iftMGetVoxelCoord(mimg, s);

      size += 1;
      means[i][0] += (u.x - means[i][0]) / size;
      means[i][1] += (u.y - means[i][1]) / size;
      aux = aux->next;
    }
    
    // computing the closest node to the centroid
    min_dist = IFT_INFINITY_FLT;
    aux = S[i]->begin;
    min_root = IFT_NIL;
    while (aux != NULL) {
      s = aux->elem;
      u = iftMGetVoxelCoord(mimg, s);
      coord[0] = u.x;
      coord[1] = u.y;
      dist = iftEuclideanDistance(means[i], coord, 2);
      if (dist < min_dist) {
        min_dist = dist;
        min_root = s;
      }
      aux = aux->next;
    }
    centroids[i] = min_root;
  }
}


/*void iftUpdateCentroids(iftGraph *graph, iftDynamicSet **S, int *centroids, int k) {
  float min_dist, dist;
  int min_root, u, s;
  iftSet *aux = NULL;
  float *mean = iftAlloc(graph->Z->nfeats, sizeof(*mean));

  for (int i = 0; i < k; i++) {
    min_dist = IFT_INFINITY_FLT;
    aux = S[i]->begin;
    min_root = IFT_NIL;
    while (aux != NULL) {
      u = aux->elem;
      s = graph->node[u].sample;
      iftCopyDblArrayToFloatArray(mean, S[i]->mean, graph->Z->nfeats);
      dist = iftEuclideanDistance(mean, graph->Z->sample[s].feat, graph->Z->nfeats);
      if (dist < min_dist) {
          min_dist = dist;
          min_root = u;
      }
      aux = aux->next;
    }
	  if (min_root != IFT_NIL)
	    centroids[i] = min_root;
  }  
}


int iftUpdateSingleCentroid(iftGraph *graph, iftDynamicSet *S) {
  float min_dist, dist;
  int min_root, u, s;
  iftSet *aux = NULL;
  float *mean = iftAlloc(graph->Z->nfeats, sizeof(*mean));

  min_dist = IFT_INFINITY_FLT;
  aux = S->begin;
  min_root = IFT_NIL;
  while (aux != NULL) {
    u = aux->elem;
    s = graph->node[u].sample;
    iftCopyDblArrayToFloatArray(mean, S->mean, graph->Z->nfeats);
    dist = iftEuclideanDistance(mean, graph->Z->sample[s].feat, graph->Z->nfeats);
    if (dist < min_dist) {
        min_dist = dist;
        min_root = u;
    }
    aux = aux->next;
  }

  return min_root;
}


void iftMUpdateCentroids(iftGraph *graph, iftDynamicSet **S, int *centroids, int k) {
  float min_dist, dist;
  int min_root, size, s;
  iftSet *aux = NULL;
  iftMImage *mimg = graph->Z->ref_data;
  float *coord = iftAllocFloatArray(3);
  iftVoxel u;
  
  // computing the means
  float** means = iftAlloc(k, sizeof(float *));
  for (int i = 0; i < k; i++) {
    means[i] = iftAllocFloatArray(2);
    aux      = S[i]->begin;
    size     = 0;
    // computing the individual means
    while (aux != NULL) {
      s = aux->elem;
      u = iftMGetVoxelCoord(mimg, s);

      size += 1;
      means[i][0] += (u.x - means[i][0]) / size;
      means[i][1] += (u.y - means[i][1]) / size;
      aux = aux->next;
    }
    
    // computing the closest node to the centroid
    min_dist = IFT_INFINITY_FLT;
    aux = S[i]->begin;
    min_root = IFT_NIL;
    while (aux != NULL) {
      s = aux->elem;
      u = iftMGetVoxelCoord(mimg, s);
      coord[0] = u.x;
      coord[1] = u.y;
      dist = iftEuclideanDistance(means[i], coord, 2);
      if (dist < min_dist) {
        min_dist = dist;
        min_root = s;
      }
      aux = aux->next;
    }
    centroids[i] = min_root;
  }
}


void iftReadDelaunayTriangulation(const char *pathname, iftGraph *graph) {
    if (!iftFileExists(pathname))
        iftError("The file pathname \"%s\" does not exists!", "iftReadDelaunayTriangulation", pathname);
    
    long nrows;
    int u, v;
    char strSeparator[2] = {' ', '\0'};

    FILE *fp = fopen(pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadDelaunayTriangulation", pathname);

    char *line = iftGetLine(fp);
    char *node;
    nrows = atoi(line);

    iftSList *sl = NULL;
    iftList *l = NULL;
    iftIntArray *iarr = NULL;

    for (int i = 0; i < nrows; i++) {
        line = iftGetLine(fp);
        sl = iftSplitString(line, strSeparator);
        l = iftCreateList();

        node = iftRemoveSListHead(sl);
        while (strcmp(node, "") != 0) {
            iftInsertListIntoTail(l, atoi(node));
            node = iftRemoveSListHead(sl);
        }
        
        iarr = iftListToIntArray(l);

        for (int j = 0; j < iarr->n - 1; j++) {
            u = iarr->val[j];
            v = iarr->val[j + 1];
            if (iftAdjSetHasElement(graph->node[u].adj, v) == 0)
                iftInsertAdjSet(&graph->node[u].adj, v, 0.0);

            if (iftAdjSetHasElement(graph->node[v].adj, u) == 0)
                iftInsertAdjSet(&graph->node[v].adj, u, 0.0);
        }

        iftDestroySList(&sl);
        iftDestroyIntArray(&iarr);
        iftDestroyList(&l);
        iftFree(line);
    }

    fclose(fp);
}*/


iftGraph *iftCreateGraph(iftDataSet *Z) {
  iftGraph *graph = (iftGraph *) iftAlloc(1, sizeof(iftGraph));
  int nnodes = Z->ntrainsamples;

  if (nnodes == 0) {
    iftError("No samples for training", "iftCreateGraph");
  }

  graph->nnodes = nnodes;
  graph->node   = (iftGraphNode *) iftAlloc (nnodes, sizeof(iftGraphNode));

  if (graph->node == NULL){
    iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateKnnGraph");
  }

  graph->is_complete   = false;
  graph->centroids     = NULL;
  graph->pathval       = iftAllocFloatArray(nnodes);
  graph->ordered_nodes = iftAllocIntArray(nnodes);
  graph->Q             = iftCreateFHeap(nnodes, graph->pathval);
  graph->Z             = Z;

#pragma omp parallel for
  for (int u = 0; u < graph->nnodes; u++){
    graph->node[u].adj     = NULL;
    graph->node[u].adjplat = NULL;
    graph->node[u].sample  = IFT_NIL;
    graph->node[u].maxarcw = 0.0;
    graph->node[u].root    = u;
    graph->node[u].pred    = IFT_NIL;
  }

  int u = 0;
  for (int s=0; s < Z->nsamples; s++){
    if (iftHasSampleStatus(Z->sample[s], IFT_TRAIN)){
      graph->node[u].sample = s;
      u++;
    }
  }

  return(graph);
}


void iftDestroyGraph(iftGraph **graph) {
  int u;
  iftGraph *aux = (*graph);

  if (aux != NULL) {
    for (u = 0; u < aux->nnodes; u++) {
      if (aux->node[u].adj != NULL) {
        iftDestroyAdjSet(&(aux->node[u].adj));
      }
      if (aux->node[u].adjplat != NULL){
        iftDestroySet(&(aux->node[u].adjplat));
      }
    }
    iftFree(aux->node);
    iftFree(aux->pathval);
    iftFree(aux->ordered_nodes);
    iftDestroyFHeap(&(aux->Q));
    iftFree(aux);
    (*graph) = NULL;
  }
}


void iftSetMGraphAdjacencySets(iftGraph *graph, iftMImage *mimg, iftAdjRel *A) {
  int p, q, j;

  for (p = 0; p < mimg->n; p++) {
    iftVoxel u = iftMGetVoxelCoord(mimg, p);
    for (j = 1; j < A->n; j++) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, j);
        if (iftMValidVoxel(mimg, v)) {
            q = iftMGetVoxelIndex(mimg, v);
            // Setting the adjacency relation for graph
            iftInsertAdjSet(&graph->node[p].adj, q, 0.0);
        }
    }
  }
}


void iftUpdateOptimumPathCost(iftGraph *graph, iftIntArray *opt_labels, int *opt_centroids, float *opt_cost) {
  float cost;
  int s;
  cost = iftSumFloatArray(graph->pathval, graph->nnodes);
  if (cost < *opt_cost) {
    *opt_cost = cost;
    for (int u = 0; u < graph->nnodes; u++) {
      s = graph->node[u].sample;
      opt_labels->val[u] = graph->Z->sample[s].group;
    }
    iftCopyIntArray(opt_centroids, graph->centroids, graph->Z->ngroups);
  }
}


void iftUpdateOptimumPathCostKMeans(iftGraph *graph, iftIntArray *opt_labels, float **mean_centroids, float **opt_centroids, float *opt_cost) {
  float cost;
  int s;
  cost = iftSumFloatArray(graph->pathval, graph->nnodes);
  if (cost < *opt_cost) {
    *opt_cost = cost;
    for (int u = 0; u < graph->nnodes; u++) {
      s = graph->node[u].sample;
      opt_labels->val[u] = graph->Z->sample[s].group;
    }

    for (int c = 0; c < graph->Z->ngroups; c++)
      for (int i = 0; i < graph->Z->nfeats; i++)
        opt_centroids[c][i] = mean_centroids[c][i];
  }
}

void iftDynamicOptimumPathForest(iftGraph *graph, iftDynamicSet **S, iftConnFuncType conn_func_type) {
  int u, v, s, t, c, g_id;
  float w, tmp = 0.0;
  float *cluster_mean = iftAlloc(graph->Z->nfeats, sizeof(*cluster_mean));
  
  for (u = 0; u < graph->nnodes; u++)
    graph->pathval[u]   = IFT_INFINITY_FLT;
  
  for (c = 0; c < graph->Z->ngroups; c++) {
    S[c]                      = iftCreateDynamicSet(graph->Z->nfeats);
    u                         = graph->centroids[c];
    s                         = graph->node[u].sample;
    graph->Z->sample[s].group = c + 1;
    graph->pathval[u]         = 0;
    graph->node[u].pred       = u;
  }
  
  for (u = 0; u < graph->nnodes; u++)
    iftInsertFHeap(graph->Q, u);
  
  while (!iftEmptyFHeap(graph->Q)) {
    u = iftRemoveFHeap(graph->Q);
    s = graph->node[u].sample;
    g_id = graph->Z->sample[s].group - 1;

    iftInsertSetDynamicSet(graph, S[g_id], u);

    // If the graph is complete
    if (graph->is_complete) { 

      for (int v = 0; v < graph->nnodes; v++) {
	      t = graph->node[v].sample;
        if (graph->Q->color[v] != IFT_BLACK && iftHasSampleStatus(graph->Z->sample[t], IFT_TRAIN)) {

          w = iftEuclideanDistanceDblFlt(S[g_id]->mean, graph->Z->sample[t].feat, graph->Z->nfeats);
          
          if (conn_func_type == IFT_MAX) {
            tmp = iftMax(graph->pathval[u], w);
          } else if (conn_func_type == IFT_SUM) {
            tmp = graph->pathval[u] + w;
          } else {
            iftError("Unknown option for connectivity function type", "iftDynamicOptimumPathForest");
          }

          if (tmp < graph->pathval[v]) {
            graph->node[v].pred       = u;
            graph->node[v].root       = graph->node[u].root;
            graph->pathval[v]         = tmp;
            graph->Z->sample[t].group = graph->Z->sample[s].group;
            iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
          } 
        }
      }
    // Otherwise 
    } else {
      
      iftAdjSet *aux = graph->node[u].adj;
      while (aux != NULL) {
        v = aux->node;
        t = graph->node[v].sample;
	
        if (graph->Q->color[v] != IFT_BLACK && iftHasSampleStatus(graph->Z->sample[t], IFT_TRAIN)) {

          w = iftEuclideanDistanceDblFlt(S[g_id]->mean, graph->Z->sample[t].feat, graph->Z->nfeats);

          if (conn_func_type == IFT_MAX) {
            tmp = iftMax(graph->pathval[u], w);
          } else if (conn_func_type == IFT_SUM) {
            tmp = graph->pathval[u] + w;
          } else {
            iftError("Unknown option for connectivity function type", "iftDynamicOptimumPathForest");
          }

          if (tmp < graph->pathval[v]) {
            graph->node[v].pred       = u;
            graph->node[v].root       = graph->node[u].root;
            graph->pathval[v]         = tmp;
            graph->Z->sample[t].group = graph->Z->sample[s].group;
            iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
          }
        }
        aux = aux->next;        
      }
    }
  }
}


void iftOptimumPathForest(iftGraph *graph, iftDynamicSet **S, iftConnFuncType conn_func_type, bool use_arcw) {
  int u, v, s, t, c;
  float w, tmp = 0.0;
  
  for (u = 0; u < graph->nnodes; u++)
      graph->pathval[u] = IFT_INFINITY_FLT;
  
  for (c = 0; c < graph->Z->ngroups; c++) {
    S[c]                                        = iftCreateDynamicSet(graph->Z->nfeats);
    u                                           = graph->centroids[c];
    s                                           = graph->node[u].sample;
    graph->Z->sample[s].group                   = c + 1;
    graph->pathval[u]                           = 0;
    graph->node[u].pred                         = IFT_NIL;
  }
  
  for (u = 0; u < graph->nnodes; u++)
    iftInsertFHeap(graph->Q, u);
  
  while (!iftEmptyFHeap(graph->Q)) {
    u = iftRemoveFHeap(graph->Q);
    s = graph->node[u].sample;
    
    iftInsertSetDynamicSet(graph, S[graph->Z->sample[s].group - 1], u);
    
    // If the graph is complete
    if (graph->is_complete) { 
      for (int v = 0; v < graph->nnodes; v++) {
	      t = graph->node[v].sample;
        if (graph->Q->color[v] != IFT_BLACK && iftHasSampleStatus(graph->Z->sample[t], IFT_TRAIN)) {
          
          w = iftEuclideanDistance(graph->Z->sample[s].feat, graph->Z->sample[t].feat, graph->Z->nfeats);

          if (conn_func_type == IFT_MAX) {
            tmp = iftMax(graph->pathval[u], w);
          } else if (conn_func_type == IFT_SUM) {
            tmp = graph->pathval[u] + w;
          } else {
            iftError("Unknown option for connectivity function type", "iftOptimumPathForest");
          }
          
          if (tmp < graph->pathval[v]) {
            graph->node[v].pred       = u;
            graph->node[v].root       = graph->node[u].root;
            graph->pathval[v]         = tmp;
            graph->Z->sample[t].group = graph->Z->sample[s].group;
            iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
          } 
        }
      }
    // Otherwise
    } else {
      iftAdjSet *aux = graph->node[u].adj;
      while (aux != NULL) {
        v = aux->node;
	      t = graph->node[v].sample;

        if (graph->Q->color[v] != IFT_BLACK && iftHasSampleStatus(graph->Z->sample[t], IFT_TRAIN)) {
          
          w = use_arcw ? aux->arcw : iftSquaredEuclideanDistance(graph->Z->sample[s].feat, graph->Z->sample[t].feat, graph->Z->nfeats);

          if (conn_func_type == IFT_MAX) {
            tmp = iftMax(graph->pathval[u], w);
          } else if (conn_func_type == IFT_SUM) {
            tmp = graph->pathval[u] + w;
          } else {
            iftError("Unknown option for connectivity function type", "iftOptimumPathForest");
          }

          if (tmp < graph->pathval[v]) {
            graph->node[v].pred       = u;
            graph->node[v].root       = graph->node[u].root;
            graph->pathval[v]         = tmp;
            graph->Z->sample[t].group = graph->Z->sample[s].group;
            iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
          }
        }
        aux = aux->next;
      }
    }
  }
}

iftGraph *iftCreateConnectedKnnGraph(iftDataSet *Z, int k, bool set_centers, int n_clusters, const char *filename) {
  iftKnnGraph *kgraph = iftCreateKnnGraph(Z, k);
  iftGraph    *graph  = iftCreateGraph(Z);
  graph->is_complete  = false;

  int   u, v;

  // keeping symmetry on the graph
  for (u = 0; u < kgraph->nnodes; u++) {
    iftAdjSet *aux = kgraph->node[u].adj;        
    while (aux != NULL) {
      v = aux->node;
      if (iftAdjSetHasElement(kgraph->node[v].adj, u) == 0)
        iftInsertAdjSet(&kgraph->node[v].adj, u, aux->arcw);
      aux = aux->next;
    }
  }

  // copying adjacency from kgraph
  for (u = 0; u < graph->nnodes; u++) {
      iftAdjSet *aux = kgraph->node[u].adj;
      while (aux != NULL) {
        v = aux->node;
        iftInsertAdjSet(&graph->node[u].adj, v, aux->arcw);
        aux = aux->next;
      }
  }

  int n_comp = iftLabelConnectedComponents(graph);
  int n_comp_r = n_comp;
  int *conn_comp_l = iftAllocIntArray(n_comp + 1);
  int *conn_comp_lr = iftAllocIntArray(n_comp + 1);

  for (int i = 0; i <= n_comp; i++) {
    conn_comp_l[i] = i;
    conn_comp_lr[i] = i;
  }

  iftDestroyKnnGraph(&kgraph);

  if (set_centers && n_comp < n_clusters) {
    char buffer[IFT_STR_DEFAULT_SIZE];
    sprintf(buffer, "Number of clusters %d is greater than number of connected components %d", n_clusters, n_comp);
    iftError(buffer, "iftCreateConnectedKnnGraph");
  }

  if (n_comp == 1)
    return graph;

  iftEdges *edges = iftGetOrderedMSTEdges(Z, IFT_INCREASING);

  /* Make the graph with k-nearest neighbors a connected graph */
  for (u = 0; u < edges->n_edges; u++) {
    int ix = edges->ind[u];
    int uu = edges->edge[ix].u;
    int vv = edges->edge[ix].v;
    if(conn_comp_l[graph->node[uu].root] != conn_comp_l[graph->node[vv].root]) {
      float arcw = iftEuclideanDistance(graph->Z->sample[uu].feat, graph->Z->sample[vv].feat, graph->Z->nfeats);          
      iftInsertAdjSet(&graph->node[uu].adj, vv, arcw);
      iftInsertAdjSet(&graph->node[vv].adj, uu, arcw);

      // propagating label between connected components
      int target_label = conn_comp_l[graph->node[vv].root];
      for (int i = 1; i < n_comp + 1; i++) {
        if (conn_comp_l[i] == target_label)
          conn_comp_l[i] = conn_comp_l[graph->node[uu].root];
      }
      
      if (n_comp_r > n_clusters) {
        iftCopyIntArray(conn_comp_lr, conn_comp_l, n_comp + 1);
        n_comp_r -= 1;
      }
    }
  }
  
  /*char result_filename[IFT_STR_DEFAULT_SIZE];
  iftCSV *csv = iftCreateCSV(graph->nnodes, graph->Z->nfeats + 1);
  char buffer[IFT_STR_DEFAULT_SIZE];
  for (int u = 0; u < graph->nnodes; u++) {
      int s = graph->node[u].sample;
      for (int j = 0; j < graph->Z->nfeats; ++j) {
          sprintf(buffer, "%f", graph->Z->sample[s].feat[j]);
          strcpy(csv->data[s][j], iftCopyString(buffer));
      }
      sprintf(buffer, "%d", conn_comp_lr[graph->node[u].root]);
      strcpy(csv->data[s][graph->Z->nfeats], iftCopyString(buffer));
  }

  strcpy(result_filename, "");
  strcat(result_filename, "./result_csv/");
  strcat(result_filename, filename);
  strcat(result_filename, "_conn_comp_rest.csv");
  iftWriteCSV(csv, result_filename, ',');
  iftDestroyCSV(&csv);

  for (u = 0; u < graph->nnodes; u++) {
      if (u == 0) {
          graph->pathval[u] = 0;
          graph->node[u].pred = IFT_NIL;
      } else {
          graph->pathval[u] = IFT_INFINITY_FLT;
      }
      iftInsertFHeap(graph->Q, u);
  }
  
  while (!iftEmptyFHeap(graph->Q)) {
    u = iftRemoveFHeap(graph->Q);
    s = graph->node[u].sample;
    if (graph->node[u].pred != IFT_NIL) {
      pred_u = graph->node[u].pred;
      t = graph->node[pred_u].sample;
      if (iftAdjSetHasElement(graph->node[u].adj, pred_u) == 0 && conn_comp[graph->node[u].root] != conn_comp[graph->node[pred_u].root]) {
        float arcw = iftEuclideanDistance(graph->Z->sample[t].feat, graph->Z->sample[s].feat, graph->Z->nfeats);          
        iftInsertAdjSet(&graph->node[u].adj, pred_u, arcw);
        iftInsertAdjSet(&graph->node[pred_u].adj, u, arcw);
        conn_comp[graph->node[u].root] = conn_comp[graph->node[pred_u].root];
        if (n_rest > n_clusters) {
          conn_comp_rest[graph->node[u].root] = conn_comp_rest[graph->node[pred_u].root];
          n_rest -= 1;
        }
      }
    }
    
    for (v = 0; v < graph->nnodes; v++) {
      if (graph->Q->color[v] != IFT_BLACK) {
        t = graph->node[v].sample;
        tmp = iftEuclideanDistance(graph->Z->sample[t].feat, graph->Z->sample[s].feat, graph->Z->nfeats);
        if (tmp < graph->pathval[v]) {
          graph->pathval[v]   = tmp;
          graph->node[v].pred = u;
          iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
        }
      }
    }
  }*/

  if (set_centers) {
    iftList *conn_comp = iftAlloc(n_comp + 1, sizeof(*conn_comp));
    int *conn_comp_size = iftAllocIntArray(n_comp + 1);
    graph->centroids = iftAllocIntArray(n_clusters);

    for (int u = 0; u < graph->nnodes; u++) {
      iftInsertListIntoTail(&conn_comp[conn_comp_lr[graph->node[u].root]], u);
      conn_comp_size[conn_comp_lr[graph->node[u].root]] += 1;
    }

    int j = 0;
    for (int i = 1; i <= n_comp + 1; i++) {
      float min = IFT_INFINITY_FLT;
      if (conn_comp_size[i] != 0) {
        //int ix = iftRandomInteger(0, conn_comp_size[i] - 1);
        int indx = 0;
        iftIntArray * iarray = iftListToIntArray(&conn_comp[i]);
        for (int u = 0; u < iarray->n; u++) {
          if (graph->Z->sample[graph->node[iarray->val[u]].sample].weight < min) {
            min = graph->Z->sample[graph->node[iarray->val[u]].sample].weight;
            indx = iarray->val[u];
          }
        }
        //graph->centroids[j] = iarray->val[ix];
        graph->centroids[j] = indx;
        j += 1;
      }
    }
  }

  return graph;
}


int iftLabelConnectedComponents(iftGraph *graph) {
  int u, v, w, l = 1;

  for (u = 0; u < graph->nnodes; u++) {
    graph->node[u].root = 0;
  }

  iftResetFHeap(graph->Q);

  for (u = 0; u < graph->nnodes; u++) {
    if (graph->node[u].root == 0) {
      graph->node[u].root = l;
      iftInsertFHeap(graph->Q, u);
      while (!iftEmptyFHeap(graph->Q)) {
        v = iftRemoveFHeap(graph->Q);
        iftAdjSet *aux = graph->node[v].adj;
        while (aux != NULL) {
          w = aux->node;
          if (graph->node[w].root == 0) {
              graph->node[w].root = graph->node[u].root;
              iftInsertFHeap(graph->Q, w);
          }
          aux = aux->next;
        }
      }
      l += 1;
    }
  }

  return l - 1;
}


int iftCheckGraphSymmetry(iftGraph *graph) {
  for (int u = 0; u < graph->nnodes; u++) {
    iftAdjSet *aux = graph->node[u].adj;
    while (aux != NULL) {
      int v = aux->node;  
      if (iftAdjSetHasElement(graph->node[v].adj, u) == 0)
        return 0;
      aux = aux->next;
    }
  }
  return 1;
}


int *iftSelectInitialCentersPercGraph(iftGraph *graph, int n_clusters, iftConnFuncType conn_func_type, float percentage, bool use_arcw) {
  
  int *indices = NULL;
  int *out_seeds = iftAllocIntArray(n_clusters);
  int index = percentage * (graph->nnodes - 1);

  for (int i = 0; i < n_clusters; i++) {
    int nseeds = i > 0 ? i : i + 1;

    iftDynamicSet **S = iftAlloc(nseeds, sizeof(*S));
    graph->centroids = iftAllocIntArray(nseeds);

    if (i == 0) {
      graph->centroids[i] = iftRandomInteger(0, graph->nnodes - 1);
    } else {
      for (int j = 0; j < i; j++)
        graph->centroids[j] = out_seeds[j];
    }
    
    graph->Z->ngroups = nseeds;

    iftOptimumPathForest(graph, S, conn_func_type, use_arcw);
    
    out_seeds[i] = iftMaxIndexFloatArray(graph->pathval, graph->nnodes);

    indices = iftAllocIntArray(graph->nnodes);
    for (int j = 0; j < graph->nnodes; j++) {
        indices[j] = j;
    }
    
    iftFQuickSort(graph->pathval, indices, 0, graph->nnodes - 1, IFT_INCREASING);
    
    out_seeds[i] = indices[index];

    iftFree(indices);
    /*for (int c = 0; c < n_clusters; c++)
      if (S[c] == NULL) {
        printf("brunao\n");
      } else {
        printf("chiquisao\n");
      }*/
      
      //iftDestroyDynamicSet(S[c]);
  }

  return out_seeds;
}


int *iftSelectInitialCenters(iftDataSet *Z, int n_clusters, iftConnFuncType conn_func_type) {
  
  int *out_seeds = iftAllocIntArray(n_clusters);

  iftGraph *graph;

  for (int i = 0; i < n_clusters; i++) {
    int nseeds = i > 0 ? i : i + 1;

    iftDynamicSet **S = iftAlloc(nseeds, sizeof(*S));
    graph = iftCreateGraph(Z);
    graph->centroids = iftAllocIntArray(nseeds);

    if (i == 0) {
      graph->centroids[i] = iftRandomInteger(0, graph->nnodes - 1);
    } else {
      for (int j = 0; j < i; j++)
        graph->centroids[j] = out_seeds[j];
    }
    
    graph->is_complete = true;
    graph->Z->ngroups = nseeds;

    iftOptimumPathForest(graph, S, conn_func_type, false);

    out_seeds[i] = iftMaxIndexFloatArray(graph->pathval, graph->nnodes);
    
    iftDestroyGraph(&graph);
  }

  return out_seeds;
}


/*void iftSetInitialCenters(iftGraph *graph, int k) {
  float min, tmp;
  int i = 1, j, u, v, s, t, next_centroid;
  graph->centroids = iftAllocIntArray(k);
  // selecting first center
  graph->centroids[0] = iftRandomInteger(0, graph->nnodes - 1);
  float *closest_centroid = iftAllocFloatArray(graph->nnodes);

  // selecting the rest of centers
  while (i < k) {
    for (u = 0; u < graph->nnodes; u++) {
      s   = graph->node[u].sample;
      min = IFT_INFINITY_FLT;
      for (j = 0; j < i; j++) {
        v = graph->centroids[j];
        t = graph->node[v].sample;
        tmp = iftEuclideanDistance(graph->Z->sample[s].feat, graph->Z->sample[t].feat, graph->Z->nfeats);
        min = tmp > min ? min : tmp;
      }
      closest_centroid[u] = min;
    }
    next_centroid = iftMaxIndexFloatArray(closest_centroid, graph->nnodes);
    graph->centroids[i] = next_centroid;
    i += 1;
  }

  iftFree(closest_centroid);
}*/


iftDataSet *iftKernelTrick(iftDataSet *Z, int n_clusters, iftConnFuncType conn_func_type, int *init_seeds) {
  iftDataSet *Zout = iftCreateDataSet(Z->nsamples, n_clusters);
  iftGraph *graph;

  Zout->iftArcWeight    = Z->iftArcWeight;
  Zout->function_number = Z->function_number;
  iftCopyRefData(Zout, Z->ref_data, Z->ref_data_type);
  Zout->ntrainsamples   = Z->ntrainsamples;

  for (int s = 0; s < Z->nsamples; s++) {
    Zout->sample[s].id        = Z->sample[s].id;
    Zout->sample[s].truelabel = Z->sample[s].truelabel;
    Zout->sample[s].label     = Z->sample[s].label;
    Zout->sample[s].group     = Z->sample[s].group;
    iftSetSampleStatus(&Zout->sample[s], Z->sample[s].status);
  }

  if (init_seeds == NULL) {
    init_seeds = iftAllocIntArray(n_clusters);
    init_seeds = iftRandomIntegers(0, Z->nsamples - 1, n_clusters);
  }

  for (int i = 0; i < n_clusters; i++) {
    graph = iftCreateGraph(Z);
    graph->centroids = iftAllocIntArray(1);
    graph->centroids[0] = init_seeds[i];
    graph->is_complete = true;
    graph->Z->ngroups = 1;

    iftDynamicSet **S = iftAlloc(1, sizeof(*S));

    iftOptimumPathForest(graph, S, conn_func_type, false);

    for (int s = 0; s < Z->nsamples; s++) {
      Zout->sample[s].feat[i] = graph->pathval[s];
    }

    iftDestroyGraph(&graph);
  }

  return Zout;
}


iftEdges *iftGetOrderedMSTEdges(iftDataSet *Z, uchar order) {
  iftEdges *edges = iftAlloc(1, sizeof *edges);
  edges->n_edges  = 0;

  iftMST *mst = iftCreateMST(Z);

  for (int u = 0; u < mst->nnodes; u++) {
    iftSet *adj = mst->node[u].adj;
    while (adj != NULL) {
      edges->edge = iftRealloc(edges->edge, (edges->n_edges + 1) * sizeof(*edges->edge));
      edges->arcw = iftRealloc(edges->arcw, (edges->n_edges + 1) * sizeof(*edges->arcw));
      edges->ind = iftRealloc(edges->ind, (edges->n_edges + 1) * sizeof(*edges->ind));
      edges->edge[edges->n_edges].u = u;
      edges->edge[edges->n_edges].v = adj->elem;
      edges->arcw[edges->n_edges] = iftEuclideanDistance(Z->sample[mst->node[u].sample].feat, Z->sample[mst->node[adj->elem].sample].feat, Z->nfeats);
      edges->ind[edges->n_edges] = edges->n_edges;
      edges->n_edges += 1;
      adj = adj->next;
    }
  }

  iftFQuickSort(edges->arcw, edges->ind, 0, edges->n_edges - 1, order);

  return edges;
}


bool iftDynamicSetHasElement(iftDynamicSet *S, int v) {
  for (iftSet *set = S->begin; set != NULL; set = set->next)
    if (set->elem == v)
      return true;
  return false;
}


float iftEuclideanDistanceDblFlt(double *d, float *f, int n) {
  float dist = 0.0f;
  for (int i = 0; i < n; i++) {
      dist += (d[i]-f[i]) * (d[i]-f[i]);
  }

  return(sqrtf(dist));
}


bool iftCheckConvergence(int *curr_cent, int **center_hist, int n_clusters, int n_hist) {
  for (int i = 0; i < n_hist; i++)
    if (iftCompareIntArrays(curr_cent, center_hist[i], n_clusters))
      return true;
  return false;
}


void iftFindOptimumCentersDynOPF(iftGraph *graph, iftCenterHist *center_hist, int n_clusters, iftConnFuncType conn_func_type) {
  float cost = FLT_MAX;
  int *ocenters = iftAllocIntArray(n_clusters); 
  iftDynamicSet **S;
  for (int i = 0; i < center_hist->n; i++) {
    iftCopyIntArray(graph->centroids, center_hist->val[i], n_clusters);
    S = iftAlloc(n_clusters, sizeof (*S));   
    iftDynamicOptimumPathForest(graph, S, conn_func_type);
    float tmp = iftSumFloatArray(graph->pathval, graph->nnodes);
    if (tmp < cost) {
      cost = tmp;
      iftCopyIntArray(ocenters, graph->centroids, n_clusters);
    }
    for (int c = 0; c < n_clusters; c++)
      iftDestroyDynamicSet(&S[c]);
  }
  iftCopyIntArray(graph->centroids, ocenters, n_clusters);
  S = iftAlloc(n_clusters, sizeof (*S));   
  iftDynamicOptimumPathForest(graph, S, conn_func_type);
  for (int c = 0; c < n_clusters; c++)
      iftDestroyDynamicSet(&S[c]);
  //iftFree(S);
  //iftFree(ocenters);
}


void iftFindOptimumCentersOPF(iftGraph *graph, iftCenterHist *center_hist, int n_clusters, iftConnFuncType conn_func_type, bool use_arcw) {
  float cost = FLT_MAX;
  int *ocenters = iftAllocIntArray(n_clusters); 
  iftDynamicSet **S;
  for (int i = 0; i < center_hist->n; i++) {
    iftCopyIntArray(graph->centroids, center_hist->val[i], n_clusters);
    S = iftAlloc(n_clusters, sizeof (*S));   
    iftOptimumPathForest(graph, S, conn_func_type, use_arcw);
    float tmp = iftSumFloatArray(graph->pathval, graph->nnodes);
    if (tmp < cost) {
      cost = tmp;
      iftCopyIntArray(ocenters, graph->centroids, n_clusters);
    }
    for (int c = 0; c < n_clusters; c++)
      iftDestroyDynamicSet(&S[c]);
    iftFree(S);
  }
  iftCopyIntArray(graph->centroids, ocenters, n_clusters);
  S = iftAlloc(n_clusters, sizeof (*S));   
  iftOptimumPathForest(graph, S, conn_func_type, use_arcw);
  for (int c = 0; c < n_clusters; c++)
      iftDestroyDynamicSet(&S[c]);
  iftFree(S);
  iftFree(ocenters);
}

/*iftImage *iftSuperpixelSegmentationByIDT(iftMImage *mimg, iftAdjRel *A, int num_superpixels, int max_iterations, iftSet **roots) {
    iftDataSet *Z = iftMImageToDataSet(mimg, NULL, 0);
    iftImage *label_img;

    Z->ntrainsamples = Z->nsamples;
    iftSetStatus(Z, IFT_TRAIN);
    
    iftGraph *graph = iftCreateGraph(Z);    
    iftSetMGraphAdjacencySets(graph, mimg, A);
    iftIterativeOPF(graph, num_superpixels, max_iterations, IFT_MAX, IFT_XY_COORD);
    iftFindOptimumCentersDynOPF(graph, graph->center_hist, num_superpixels, IFT_MAX);    

    for (int i = 0; i < num_superpixels; i++)
        iftInsertSet(roots, graph->centroids[i]);

    label_img = iftDataSetClusterInformationToLabelImage(graph->Z, NULL, false);

    return label_img;
}*/

iftIGraph *iftSuperpixelSegmentationByIDT(iftMImage *mimg, iftImage *seeds, iftAdjRel *A, int num_superpixels, int max_iterations) {
  int j, p, s, *seed, *optseed, **shist, nhist, nseeds, iter = 0;
  double totpathcost, minpathcost = IFT_INFINITY_DBL;
  bool converged = false;
  iftIGraph *igraph;

  igraph = iftImplicitIGraph(mimg, NULL, A);

  double *pvalue = iftAllocDoubleArray(igraph->nnodes);  

  iftDynamicSet **S = iftAlloc(num_superpixels, sizeof(*S));
  
  if (seeds == NULL) {
    nseeds = num_superpixels;
    seed = iftRandomIntegers(0, igraph->nnodes - 1, nseeds);
  } else {
    nseeds = 0;
    seed = iftAllocIntArray(nseeds);
    for (s = 0; s < igraph->nnodes; s++) {  
      p = igraph->node[s].voxel;
      if (seeds->val[p] != 0) {
        seed[nseeds] = s;
        nseeds += 1;
      }
    }
  }

  optseed = iftAllocIntArray(nseeds);

  shist = iftAlloc(1, sizeof(*shist));
  nhist = 0;

  while (!converged && iter < max_iterations) {
    shist = iftRealloc(shist, (nhist + 1) * sizeof(*shist));
    shist[nhist] = iftAllocIntArray(nseeds);
    iftCopyIntArray(shist[nhist], seed, nseeds);
    nhist += 1;

    iftIDynamicIFT(igraph, mimg, S, seed, nseeds, pvalue);
    
    totpathcost = iftSumDoubleArray(pvalue, igraph->nnodes);
    if (totpathcost < minpathcost) {
      minpathcost = totpathcost;
      iftCopyIntArray(optseed, seed, nseeds);
    }

    iftIRecomputeSeeds(igraph, S, seed, nseeds);       

    for (j = 0; j < nseeds; j++) {
      iftDestroyDynamicSet(&S[j]);
    }

    iter += 1;

    for (j = 0; j < nhist; j++) {
      if (iftCompareIntArrays(seed, shist[j], nseeds)) {
        converged = true;
        break;
      }
    }
  }

  iftIDynamicIFT(igraph, mimg, S, optseed, nseeds, pvalue);

  for (j = 0; j < nseeds; j++) {
    iftDestroyDynamicSet(&S[j]);
  }
  
  iftFree(S);

  return igraph;
}

void iftIDynamicIFT(iftIGraph *igraph, iftMImage *mimg, iftDynamicSet **S, int *seed, int nseeds, double *pvalue) {
  double w, tmp = 0.0;
  int i, j, p, q, s, t;
  iftVoxel u, v;
  iftDHeap *Q;

  for (s = 0; s < igraph->nnodes; s++) {
    p                 = igraph->node[s].voxel;
    igraph->pvalue[p] = IFT_INFINITY_DBL;
    pvalue[s]         = igraph->pvalue[p];
    igraph->pred[p]   = IFT_NIL;
  }

  i = 0;
  for (j = 0; j < nseeds; j++) {
    S[j]              = iftCreateDynamicSet(igraph->nfeats);
    s                 = seed[j];
    p                 = igraph->node[s].voxel;
    igraph->pvalue[p] = 0;
    pvalue[s]         = igraph->pvalue[p];
    igraph->pred[p]   = IFT_NIL;
    igraph->label[p]  = i + 1;
    igraph->root[p]   = p;
    i += 1;
  } 

  Q = iftCreateDHeap(igraph->nnodes, pvalue);   

  for (s = 0; s < igraph->nnodes; s++) {
    iftInsertDHeap(Q, s);
  }

  while (!iftEmptyDHeap(Q)) {
    s = iftRemoveDHeap(Q);
    p = igraph->node[s].voxel;
    igraph->pvalue[p] = pvalue[s];
    u = iftGetVoxelCoord(igraph->index, p); 

    iftInsertDynamicSet(S[igraph->label[p] - 1], mimg, p);

    for (i = 1; i < igraph->A->n; i++) {
      v = iftGetAdjacentVoxel(igraph->A, u, i);
      if (iftValidVoxel(igraph->index, v)) {
        q = iftGetVoxelIndex(igraph->index, v);
        t = igraph->index->val[q];

        if (t != IFT_NIL && Q->color[t] != IFT_BLACK) {
          w = (double) iftEuclideanDistanceDblFlt(S[igraph->label[p] - 1]->mean, 
                                                  igraph->feat[q], 
                                                  igraph->nfeats);
          tmp = iftMax(igraph->pvalue[p], w);
          
          if (tmp < pvalue[t]) {
            pvalue[t]        = tmp;
            igraph->root[q]  = igraph->root[p];
            igraph->label[q] = igraph->label[p];
            igraph->pred[q]  = p;
            iftGoUpDHeap(Q, Q->pos[t]);
          }
        }
      }
    }
  }
}

void iftIRecomputeSeeds(iftIGraph *igraph, iftDynamicSet **S, int *seed, int nseeds) {
  int j, p, s, dsize, next_seed;
  float **meanp, min_dist, dist, *pcoor;
  iftVoxel u;
  iftSet *aux;

  meanp = iftAlloc(nseeds, sizeof(*meanp));
  for (j = 0; j < nseeds; j++) {
    // computing the mean pixel for each dynamic set
    meanp[j] = iftAllocFloatArray(2);
    aux      = S[j]->begin;
    dsize     = 0;
    while (aux != NULL) {
      s = aux->elem;
      p = igraph->node[s].voxel;
      u = iftGetVoxelCoord(igraph->index, p);

      dsize += 1;

      meanp[j][0] += (u.x - meanp[j][0]) / dsize;
      meanp[j][1] += (u.y - meanp[j][1]) / dsize;
      
      aux = aux->next;
    }

    // seed recomputation
    min_dist = IFT_INFINITY_FLT;
    next_seed = IFT_NIL;
    aux = S[j]->begin;
    pcoor = iftAllocFloatArray(2);
    while (aux != NULL) {
      s = aux->elem;
      p = igraph->node[s].voxel;
      u = iftGetVoxelCoord(igraph->index, p);

      pcoor[0] = u.x;
      pcoor[1] = u.y;

      dist = iftFeatDistance(meanp[j], pcoor, 2);
      
      if (dist < min_dist) {
        min_dist = dist;
        next_seed = s;
      }

      aux = aux->next;        
    }

    seed[j] = next_seed; 
  }
}