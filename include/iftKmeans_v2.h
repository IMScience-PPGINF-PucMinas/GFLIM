#ifndef IFT_KMEANS_V2_H_
#define IFT_KMEANS_V2_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftDataSet.h"
#include "ift/segm/DynamicTrees.h"  


/**
 * @brief Distance function between two feature (float) vectors.
 * @date Dez, 2021
 * @author Alexandre Falcao.
 */
typedef float (*iftKmeansDistance)(float *f1, float *f2, int n);
  
typedef struct ift_centroid {
    /** Centroid id **/
    int   id;
    /** Number of elements in hist_feat **/
    int   n_hist;
    /** Centroid features **/
    float *feat;
    /** Centroid history features **/
    float **hist_feat;
} iftCentroid;

typedef struct ift_centroid_hist {
    /** Number of clusters **/
    int          n_clusters;
    /** Number of iterations **/
    int          n_iter;
    /** Number of features **/
    int          n_feats;
    /** Centroids **/
    iftCentroid  **centroid;
} iftCentroidsHist;

/** 
 * Clusters id are stored in the group attribute of each sample. The most representative elements 
 * (i.e. samples closest to the arithmetic centroids or medoids) are labeled as IFT_PROTOTYPE once 
 * the algorithm finished its execution. In order to get the arithmetic centroids, it is required to 
 * pass a structure iftCentroidsHist as input parameter for iftKMeans. An example of this procedure 
 * is in the demos, in the file ift/demo/Clustering/iftKMeansDemo.c 
 *
 * Examples of Kmeans calling:
 * 
 * iftDataSet *Z = ...;
 * int n_clusters = 2;
 * int max_iter = 20;
 * float tol = 0.001;
 * 
 * iftKMeans(Z, n_clusters, max_iter, tol, NULL, NULL);
 * 
 * If using iftCentroidsHist
 * 
 * iftCentroidsHist *cent_hist = iftAlloc(1, sizeof(*cent_hist));
 * 
 * iftKMeans(Z, n_clusters, max_iter, tol, NULL, cent_hist);
 * 
 * If centroids are chosen a priori
 * 
 * int *init_centroids = iftAllocIntArray(n_clusters);
 * init_centroids[0] = 10; // it can be any other sample index in Z
 * init_centroids[1] = 20; // it can be any other sample index in Z
 * 
 * iftKMeans(Z, n_clusters, max_iter, tol, init_centroids, NULL);
 * 
 * or also 
 * 
 * iftKMeans(Z, n_clusters, max_iter, tol, init_centroids, cent_hist);
**/

/**
 * @brief K-means algorithm
 * @param iftDataSet *Z --- input dataset
 * @param n_clusters --- number of clusters
 * @param max_iter --- maximum number of iterations
 * @param tol --- (e.g. 0.001) tolerance of the difference of cluster centers of two consecutive iterations
 * @param init_centroids --- (e.g. NULL) (optional) set of initial centroids
 * @param centroids_hist --- (e.g. NULL) (optional) variable to store the centroids' history
 * @param distance function in the feature space (default is
 * iftEuclideanDistance). 
 * @author David Aparco Cardenas  
 * @date Jun 2nd, 2020
 */

void iftKMeans(iftDataSet *Z, int n_clusters, int max_iter, float tol, int *init_centroids, iftCentroidsHist *centroids_hist, iftKmeansDistance distance);

/**
 * @brief Set status of clusters' representatives (samples closest to centroids) to IFT_PROTOTYPE.
 * @param iftDataSet *Z --- input dataset
 * @param S --- dynamic sets representing the clusters
 * @author David Aparco Cardenas 
 * @date Jun 2nd, 2020
 */

void iftSetDataSetPrototypeStatus(iftDataSet *Z, iftDynamicSet **S);

/**
 * @brief Insert centroid history.
 * @param centroid --- centroid
 * @param centroid_feat --- new centroid history record to insert
 * @param n_feats --- number of features of each sample
 * @author David Aparco Cardenas 
 * @date Jun 2nd, 2020
 */

void iftInsertCentroidHistory(iftCentroid *centroid, float *centroid_feat, int n_feats);

/**
 * @brief Get dataset prototypes.
 * @param Z --- input dataset
 * @return array of id's of samples whose status is IFT_PROTOTYPE
 * @author David Aparco Cardenas 
 * @date Jun 2nd, 2020
 */

int *iftGetDataSetPrototypes(iftDataSet *Z);

/**
 * @brief Insert sample to dynamic set.
 * @param Z --- input dataset
 * @param S --- dynamic sets representing the clusters
 * @param s --- id of sample to be inserted
 * @author David Aparco Cardenas 
 * @date Jun 2nd, 2020
 */

void iftInsertSampleDynamicSet(iftDataSet *Z, iftDynamicSet *S, int s);

/**
 * @brief K-medoids algorithm
 * @param iftDataSet *Z --- input dataset
 * @param n_clusters --- number of clusters
 * @param max_iter --- maximum number of iterations
 * @param distance function in the feature space (default is
 * iftEuclideanDistance). 
 * @author David Aparco Cardenas  
 * @date Dez, 2020
 */

void iftKMedoids(iftDataSet *Z, int n_clusters, int max_iter, iftKmeansDistance distance);


#ifdef __cplusplus
}
#endif

#endif
