#ifndef IFT_KMEANS_H_
#define IFT_KMEANS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftDataSet.h"



/**
 * @brief Cluster a dataset by KMeans method
 * @author Adan Echemendia
 * @param Z A dataset
 * @param k Desired number of clusters
 * @param maxIterations Maximum number of iterations
 * @param minImprovement Threshold to stop the mean computations
 * @param kmedoids 1 if we want to do kmedoids, 0 if we want the original algorithm
 * @char grid_sampling 1 if the centers are initially chosen by grid sampling in case that the dataset represents an
 * image
 * @char debug 1 if we want to print debug information, 0 otherwise
 */
void iftClusterDataSetByKMeans(iftDataSet *Z, int nb_cluster, int maxIterations, float minImprovement, char kmedoids,
															 char grid_sampling, char debug);


void  iftKmeansRun(int bKmedoids,iftDataSet* Z,iftDataSet** pZk,int maxIterations,float minImprovement);
void  iftSphericalKmeansRun(iftDataSet* Z,iftDataSet** pZk,int maxIterations);

void iftSimpleKmeansRun(iftDataSet* Z, iftDataSet** pZk, int maxIterations,float minImprovement);

/**
 * @brief Initializes k centroid getting random samples from dataset Z
 * @param Z A dataset
 * @param k Number of centroids
 * @return A new dataset with the centroids
 */
iftDataSet* iftKmeansInitCentroidsFromSamples(iftDataSet* Z, int k);

iftDataSet* iftKmeansInitCentroidsRandomNormal(iftDataSet* Z, int k);

iftDataSet *iftNearestSamplesFromCentroids(iftDataSet* Zorig, iftDataSet *Zc, iftDataSet* Zk, int m, int include_centroids);
iftDataSet *iftBuildDataSetAsNearestSamplesFromCentroids(iftDataSet* Zorig, iftDataSet *Zc,
			iftDataSet* Zk, int include_centroids, char *method);

#ifdef __cplusplus
}
#endif

#endif
