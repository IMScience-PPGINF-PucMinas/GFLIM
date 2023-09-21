#include "iftHierarchicCluster.h"

#include "ift/core/io/Stream.h"

/**
@file
@brief Clustering methods using iftClusterHierarchy structure over iftDataSet.
*/

/* TODO
 *  - Add sampleID initialization
 *  - Test pretty much everything
 */

/* Private Functions */

iftClusterHierarchyCluster **iftInitClusterArray(int n);
iftClusterHierarchyCluster *iftInitCluster();
iftClusterSample *iftInitSampleArray(int n);
void iftResizeSampleArray(iftClusterSample **SS, int size);
iftClusterHierarchyCluster *iftDataSetToCluster(iftDataSet *Z);
void iftApplyClusterPartition(iftClusterHierarchy *H, int level, int cluster, int nPartitions, iftClusterHierarchyCluster ***CC);
void iftDestroyCluster(iftClusterHierarchyCluster **C);
void iftAddNewLevelToClusterHierarchy(iftClusterHierarchy *H);
inline float iftDistanceKmeansCp(float *f1, float *f2, int dim)
{
  float dist = 0.0f;

  for (int i = 0; i < dim; i++)
    dist += (f1[i] - f2[i]) * (f1[i] - f2[i]);

  return dist;
}

// for qsort and bsearch
int iftHierarchyCmp (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}

/*-------------------*/

/**
 * @brief Initialize an empty hierarchy.
 *
 * @param featureSize Size of samples feature vectors.
 * @param path Optional writable path for disk storage operations.
 * @return Pointer to created hierarchy.
 */
iftClusterHierarchy *iftCreateClusterHierarchy(int featureSize, char *path)
{
  iftClusterHierarchy *H = (iftClusterHierarchy*) iftAlloc(1, sizeof *H);

  if (H == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateClusterHierarchy");

  if (path != NULL) {
    H->path = iftAllocCharArray(strlen(path) + 1);
    strcpy(H->path, path);
  } else {
      iftError("Hierarchy needs a working path.", "iftCreateClusterHierarchy");
  }
  H->totalSamples = 0;
  H->featureSize = featureSize;
  H->levelCount = 0;
  iftAddNewLevelToClusterHierarchy(H);

  return H;
}

iftClusterHierarchyCluster **iftInitClusterArray(int n)
{
  if (n <= 0)
    return NULL;

  iftClusterHierarchyCluster **C = (iftClusterHierarchyCluster**) iftAlloc(n, sizeof *C);

  if (C == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftInitClusterArray");

  for (int cc = 0; cc < n; cc++)
    C[cc] = iftInitCluster();

  return C;
}

iftClusterHierarchyCluster *iftInitCluster()
{
  iftClusterHierarchyCluster *C = (iftClusterHierarchyCluster*) iftAlloc(1, sizeof *C);

  if (C == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftInitCluster");

  C->sampleCount = 0;
  C->status = CLUSTER_HIERARCHY_IN_MEMORY;
  C->sample = NULL;

  return C;
}

iftClusterSample *iftInitSampleArray(int n)
{
  iftClusterSample *S = (iftClusterSample*) iftAlloc(n, sizeof *S);

  if (S == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftInitSampleArray");

  return S;
}

void iftResizeSampleArray(iftClusterSample **SS, int size)
{
  *SS = (iftClusterSample*) iftRealloc(*SS, size * sizeof *SS);

  if (*SS == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftResizeSampleArray");
}

void iftAddNewLevelToClusterHierarchy(iftClusterHierarchy *H)
{
  H->levelCount += 1;
  H->level = (iftClusterHierarchyLevel**) iftRealloc(H->level, H->levelCount * sizeof *(H->level));

  if (H->level == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftAddNewLevelToClusterHierarchy");

  H->level[H->levelCount - 1] = (iftClusterHierarchyLevel*) iftAlloc(1, sizeof *(H->level[H->levelCount - 1]));

  if (H->level[H->levelCount - 1] == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftAddNewLevelToClusterHierarchy");

  H->level[H->levelCount - 1]->size = 0;
  H->level[H->levelCount - 1]->clusterCount = 0;
  H->level[H->levelCount - 1]->cluster = NULL;
}

iftClusterHierarchyCluster *iftDataSetToCluster(iftDataSet *Z) 
{
  iftClusterHierarchyCluster *C = iftInitCluster();

  C->sampleCount = Z->nsamples;
  C->sample = (iftClusterSample*) iftAlloc(Z->nsamples, sizeof *(C->sample));

  if (C->sample == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftDataSetToCluster");

  for (int s = 0; s < Z->nsamples; s++) {
    C->sample[s].dsSample = &(Z->sample[s]);
    C->sample[s].isPromoted = 0;
    C->sample[s].parentCluster = -1;
  }

  return C;
}


void iftAddDataSetToClusterHierarchy(iftClusterHierarchy *H, int level, iftDataSet *Z)
{
  if (H->featureSize != Z->nfeats)
      iftError("DataSet and Hierarchy features don't match", "iftAddDataSetToClusterHierarchy");

  iftClusterHierarchyCluster *C = iftDataSetToCluster(Z);

  H->totalSamples += Z->nsamples;
  H->level[level]->size += Z->nsamples;
  H->level[level]->clusterCount += 1;
  H->level[level]->cluster = (iftClusterHierarchyCluster**) iftRealloc(H->level[level]->cluster, H->level[level]->clusterCount * sizeof *(H->level[level]->cluster));

  if (H->level[level]->cluster == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftAddDataSetToClusterHierarchy");

  H->level[level]->cluster[H->level[level]->clusterCount-1] = C;
}

// Creates a temporary dataset over a cluster to apply existing
//  algorithms over iftDataset.
iftDataSet *iftClusterToDataSet(iftClusterHierarchy *H, int level, int cluster) 
{
  iftClusterHierarchyCluster *C = H->level[level]->cluster[cluster];
  iftDataSet *Z = iftCreateDataSet(C->sampleCount, H->featureSize);

  // Currently copying only features, more sample info may be necessary
  //  for some cluster partitions.
  for (int s = 0; s < C->sampleCount; s++) {
    for (int i = 0; i < H->featureSize; i++)
      Z->sample[s].feat[i] = C->sample[s].dsSample->feat[i];
  }

  return Z;
}

void iftApplyKNNOnCluster(iftClusterHierarchy *H, int level, int cluster, int k, int maxIterations, float minImprovement) {
  iftDataSet *Z = iftClusterToDataSet(H, level, cluster);
  iftDataSet *Zk = iftKmeansInitCentroidsFromSamples(Z,k);

  iftKmeansRun(0, Z, &Zk, maxIterations, minImprovement);

  // Create a cluster for each label
  iftClusterHierarchyCluster **C = iftInitClusterArray(Z->ngroups);

  // Count number of samples per cluster and calculates distance from centroids
  float *centroidDistance = iftAllocFloatArray(Z->nsamples);
  int *distanceIndex = iftAllocIntArray(Z->nsamples);
  for (int s = 0; s < Z->nsamples; s++) {
    int label = Z->sample[s].label - 1;
    C[label]->sampleCount += 1;
    centroidDistance[s] = iftDistanceKmeansCp(Zk->sample[label].feat, Z->sample[s].feat, Z->nfeats);
    distanceIndex[s] = s;
  }

  iftFHeapSort(centroidDistance, distanceIndex, Z->nsamples, IFT_INCREASING);

  // Allocate sample data for each cluster
  for (int c = 0; c < Z->ngroups; c++)
    C[c]->sample = iftInitSampleArray(C[c]->sampleCount);

  // Assign ordered samples to new Clusters
  int *clusterPos = iftAllocIntArray(Z->ngroups);
  for (int s = 0; s < Z->nsamples; s++) {
    int label = Z->sample[distanceIndex[s]].label - 1;
    C[label]->sample[clusterPos[label]] = H->level[level]->cluster[cluster]->sample[distanceIndex[s]];
    clusterPos[label] += 1;
  }

  iftFree(centroidDistance);
  iftFree(distanceIndex);
  iftFree(clusterPos);
  iftDestroyDataSet(&Zk);

  // Assign new clusters to the Hierarchy
  iftApplyClusterPartition(H, level, cluster, Z->ngroups, &C);

  iftDestroyDataSet(&Z);
}

void iftApplyOPFOnCluster(iftClusterHierarchy *H, int level, int cluster, int kmax) {
  iftDataSet *Z = iftClusterToDataSet(H, level, cluster);
  iftSetStatus(Z,IFT_TRAIN);

  iftKnnGraph *graph = iftCreateKnnGraph(Z,kmax);
  iftUnsupTrain(graph, iftNormalizedCut);

  // Create a cluster for each label
  iftClusterHierarchyCluster **C = iftInitClusterArray(Z->ngroups);

  // Count number of samples per Cluster
  for (int s = 0; s < Z->nsamples; s++) {
    int label = Z->sample[s].label - 1;
    C[label]->sampleCount += 1;
  }

  // Allocate sample data for each cluster
  for (int c = 0; c < Z->ngroups; c++)
    C[c]->sample = iftInitSampleArray(C[c]->sampleCount);

  // Assign ordered graph samples to clusters
  int *clusterPos = iftAllocIntArray(Z->ngroups);
  for (int n = 0; n < graph->nnodes; n++) {
    int sample = graph->node[graph->ordered_nodes[n]].sample;
    int label = Z->sample[sample].label - 1;
    C[label]->sample[clusterPos[label]] = H->level[level]->cluster[cluster]->sample[sample];
    clusterPos[label] += 1;
  }

  iftDestroyKnnGraph(&graph);
  iftFree(clusterPos);

  // Assign new clusters to the Hierarchy
  iftApplyClusterPartition(H, level, cluster, Z->ngroups, &C);

  iftDestroyDataSet(&Z);
}

void iftApplyClusterPartition(iftClusterHierarchy *H, int level, int cluster, int nPartitions, iftClusterHierarchyCluster ***CC)
{
  iftClusterHierarchyCluster **C = *CC;
  iftClusterHierarchyLevel *L = H->level[level];
  iftClusterHierarchyCluster *oldCluster = L->cluster[cluster];
  int clusterTail = L->clusterCount;
  L->clusterCount += nPartitions - 1; // already have 1 from old cluster
  L->cluster = (iftClusterHierarchyCluster**) iftRealloc(L->cluster, L->clusterCount * sizeof *(L->cluster));
  if (L->cluster == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftApplyClusterPartition");

  // Assign first new cluster where the old one was and the rest on the tail
  L->cluster[cluster] = C[0];
  for (int cc = 1; cc < nPartitions; cc++) {
    L->cluster[clusterTail] = C[cc];
    clusterTail += 1;
  }

  iftDestroyCluster(&oldCluster);
  iftFree(C);
  *CC = NULL;
}

void iftApplyClassSplitOnCluster(iftClusterHierarchy *H, int level, int cluster)
{
  iftClusterHierarchyLevel *L = H->level[level];

  // Avoid degenerate case
  if (L->cluster[cluster]->sampleCount <= 0)
    return;

  // Count how many classes the cluster have
  iftSet *classSet = NULL;
  iftInsertSet(&classSet, L->cluster[cluster]->sample[0].dsSample->truelabel);
  for (int i = 1; i < L->cluster[cluster]->sampleCount; ++i)
    iftUnionSetElem(&classSet, L->cluster[cluster]->sample[i].dsSample->truelabel);
  int classCount = iftSetSize(classSet);

  // Map existing classes to a contiguous index
  int *classIndex = iftAllocIntArray(classCount);
  iftSet *setAux = classSet;
  for (int i = 0; i < classCount; ++i, setAux = setAux->next)
    classIndex[i] = setAux->elem;
  qsort(classIndex, classCount, sizeof *classIndex, iftHierarchyCmp);
  iftDestroySet(&classSet);

  // Create cluster array where each cluster corresponds to a class
  iftClusterHierarchyCluster **C = iftInitClusterArray(classCount);

  // Count number of samples per Cluster
  for (int s = 0; s < L->cluster[cluster]->sampleCount; s++) {
    int classVal = L->cluster[cluster]->sample[s].dsSample->truelabel;
    int *indexPtr = (int*) bsearch(&classVal, classIndex, classCount, sizeof *classIndex, iftHierarchyCmp);
    int index = (int)(indexPtr - classIndex);
    C[index]->sampleCount += 1;
  }

  // Allocate sample data for each cluster
  for (int c = 0; c < classCount; c++)
    C[c]->sample = iftInitSampleArray(C[c]->sampleCount);

  // Assing samples to new clusters
  int *clusterPos = iftAllocIntArray(classCount);
  for (int s = 0; s < L->cluster[cluster]->sampleCount; s++) {
    int classVal = L->cluster[cluster]->sample[s].dsSample->truelabel;
    int *indexPtr = (int*) bsearch(&classVal, classIndex, classCount, sizeof *classIndex, iftHierarchyCmp);
    int index = (int)(indexPtr - classIndex);
    C[index]->sample[clusterPos[index]++] = L->cluster[cluster]->sample[s];
  }

  iftFree(classIndex);
  iftFree(clusterPos);
  
  // Apply partition in the hierarchy
  iftApplyClusterPartition(H, level, cluster, classCount, &C);
} 

// Moves the first cluster to the end of the cluster array and
//  puts an empty cluster in its place.
void iftMoveReferenceCluster(iftClusterHierarchy *H, int level)
{
  int size = H->level[level]->clusterCount += 1;
  iftClusterHierarchyCluster **C = H->level[level]->cluster;
  C = (iftClusterHierarchyCluster**) iftRealloc(C, size * sizeof *C);

  if (C == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftMoveReferenceCluster");

  C[size-1] = C[0];
  C[0] = iftInitCluster();
}

/* Promoted samples always go to first cluster of next level.
 * Use iftMoveReferenceCluster to free it if necessary. */
void iftPromoteClusterSample(iftClusterHierarchy *H, int level, int cluster, int sample) {
  // Allocate new level if it does not exist yet
  if (level + 1 == H->levelCount)
    iftAddNewLevelToClusterHierarchy(H);

  // Create cluster 0 if necessary to fill with promoted samples
  if (H->level[level+1]->clusterCount == 0) {
    H->level[level+1]->clusterCount += 1;
    H->level[level+1]->cluster = iftInitClusterArray(1);
  }

  iftClusterHierarchyCluster *C = H->level[level+1]->cluster[0];
  iftClusterSample *S = &(H->level[level]->cluster[cluster]->sample[sample]);

  // Add sample to last position of the sample array
  H->level[level+1]->size += 1;
  C->sampleCount += 1;
  iftResizeSampleArray(&C->sample, C->sampleCount);
  C->sample[C->sampleCount - 1].dsSample = S->dsSample;
  C->sample[C->sampleCount - 1].parentCluster = cluster;
  C->sample[C->sampleCount - 1].isPromoted = 0;
  S->isPromoted = 1;
}

void iftDestroyCluster(iftClusterHierarchyCluster **C)
{
  iftFree((*C)->sample);
  iftFree(*C);
  *C = NULL;
}

void iftDestroyHierarchy(iftClusterHierarchy **H)
{
  for (int l = 0; l < (*H)->levelCount; l++) {
    for (int c = 0; c < (*H)->level[l]->clusterCount; c++) {
      iftFree((*H)->level[l]->cluster[c]->sample);
      iftFree((*H)->level[l]->cluster[c]);
    }
    iftFree((*H)->level[l]->cluster);
    iftFree((*H)->level[l]);
  }
  iftFree((*H)->level);
  iftFree((*H)->path);
  iftFree((*H));
  *H = NULL;
}

void iftDumpHierarchyMetadata(iftClusterHierarchy *H, int showClusters, int showSamples, int showFeatures)
{
  // Main cluster struct
  printf("Working path: %s\n", H->path);
  printf("Total samples: %d\n", H->totalSamples);
  printf("Feature size: %d\n", H->featureSize);
  printf("Number of levels: %d\n", H->levelCount);
  for (int i = 0; i < H->levelCount; i++) {
    // Level struct
    iftClusterHierarchyLevel *L = H->level[i];
    printf("\n- Hierarchy level %d -\n", i);
    printf("Level size: %d\n", L->size);
    printf("Cluster Count: %d\n", L->clusterCount);
    
    if (!showClusters)
      continue;

    for (int j = 0; j < L->clusterCount; j++) {
      // Cluster struct
      iftClusterHierarchyCluster *C = L->cluster[j];
      printf("-- Cluster %d --\n", j);
      printf("Sample Count: %d\n", C->sampleCount);
      printf("Cluster Load Status: %d\n", C->status);

      if(!showSamples)
        continue;

      for (int s = 0; s < C->sampleCount; s++) {
        // Sample struct
        printf("--- Sample %d ---\n", s);
        printf("Parent cluster: %d\n", C->sample[s].parentCluster);
        printf("Label: %d\n", C->sample[s].dsSample->label);
        printf("Is Promoted: %d\n", C->sample[s].isPromoted);

        if(!showFeatures)
          continue;

        printf("Features:");
        for (int dim = 0; dim < H->featureSize; dim++)
          printf(" %f", C->sample[s].dsSample->feat[dim]);
        printf("\n");
      }
    }
  }
}

