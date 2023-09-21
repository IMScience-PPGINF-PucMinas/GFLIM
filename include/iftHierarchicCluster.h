#ifndef _IFT_HIERARCHIC_CLUSTER_H_
#define _IFT_HIERARCHIC_CLUSTER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftDataSet.h"
#include "iftClustering.h"
#include "iftKmeans.h"
#include "iftGraphics.h"

  typedef struct ift_cluster_hierarchy iftClusterHierarchy; 
  typedef struct ift_cluster_hierarchy_level iftClusterHierarchyLevel;
  typedef struct ift_cluster iftClusterHierarchyCluster;
  typedef struct ift_cluster_sample iftClusterSample;

  typedef enum {
    CLUSTER_HIERARCHY_IN_MEMORY,
    CLUSTER_HIERARCHY_IN_DISK,
    CLUSTER_HIERARCHY_INVALID_STATE
  } iftClusterLoadStatus;

  struct ift_cluster_hierarchy {
    char *path;     /* Where samples can be stored in disk */
    int totalSamples;
    int featureSize;
    int levelCount;
    iftClusterHierarchyLevel **level;
  };

  struct ift_cluster_hierarchy_level {
    int size;
    int clusterCount;
    iftClusterHierarchyCluster **cluster;
  };

  struct ift_cluster {
    int sampleCount;
    iftClusterLoadStatus status;
    iftClusterSample *sample; // When possible, ordered by relevance (e.g. sample[0] is centroid)
  };

  struct ift_cluster_sample {
    const iftSample *dsSample; // Sample from an iftDataset in memory
    int sampleID; // Identifies sample globally (including in disk)
    int isPromoted;
    int parentCluster;
  };

  iftClusterHierarchy *iftCreateClusterHierarchy(int featureSize, char *path);
  void                 iftAddDataSetToClusterHierarchy(iftClusterHierarchy *H, int level, iftDataSet *Z);
  iftDataSet          *iftClusterToDataSet(iftClusterHierarchy *H, int level, int cluster);
  void                 iftApplyKNNOnCluster(iftClusterHierarchy *H, int level, int cluster, int K, int maxIterations, float minImprovement);
  void                 iftApplyOPFOnCluster(iftClusterHierarchy *H, int level, int cluster, int kmax);
  void                 iftApplyRandomSplitOnCluster(iftClusterHierarchy *H, int level, int cluster); //empty TODO
  void                 iftApplyClassSplitOnCluster(iftClusterHierarchy *H, int level, int cluster); 
  void                 iftMoveReferenceCluster(iftClusterHierarchy *H, int level);
  void                 iftPromoteClusterSample(iftClusterHierarchy *H, int level, int cluster, int sample);
  void                 iftDestroyHierarchy(iftClusterHierarchy **H);

  /* Standard hierarchy builders
   *   Fixed hierarchy architectures for ease of use.
   *   TODO Define some
   */

  /* Debug only */
  void iftDumpHierarchyMetadata(iftClusterHierarchy *H, int showClusters, int showSamples, int showFeatures);

#ifdef __cplusplus
}
#endif

#endif /*_IFT_HIERARCHIC_CLUSTER_H_*/
