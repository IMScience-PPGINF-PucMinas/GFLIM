/**
 * @file iftClustering.h
 * @brief Definitions and functions about clustering methods.
 * @author Adan Echemendia
 * @date May 17, 2016
 *
  */
#ifndef IFT_CLUSTERING_H_
#define IFT_CLUSTERING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/FHeap.h"

#include "iftCommon.h"
#include "iftDataSet.h"
#include "iftAdjSet.h"
#include "iftMSPS.h"
#include "iftMST.h"
#include "iftSort.h"
#include "iftMetrics.h"

/**
 * @brief Node of a KNN Graph
 * @author Adan Echemendia
 * @date May 17, 2016
  */
typedef struct ift_knnnode {
    /**  Maximum arc weight from the node to its neighbors */
    float maxarcw;
    /** Corresponding root node in the graph */
    int   root;
    /** Corresponding training sample in the original dataset */
    int   sample;
    /** predecessor node*/
    int pred;
    /** predecessor node weight*/
    float predWeight;
    /** List of adjacent nodes */
    iftAdjSet *adj;
    /** List of adjacent nodes on plateaus of density */
    iftSet    *adjplat;
} iftKnnNode;

/**
 * @brief KNN Graph composed by the training samples of a dataset
 * @author Adan Echemendia
 * @date May 17, 2016
 *
 */
//! swig(extend = iftKnnGraphExt.i, destroyer = iftDestroyKnnGraph)
typedef struct ift_knngraph {
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

/**
 * @brief Generic type to define a cut function
 * @details [long description]
 *
 * @param graph Input clustered graph
 * @return  A graph cut function value
 */
typedef float (*iftKnnGraphCutFun)(iftKnnGraph *graph);

/**
 * @brief Creates a KNN Graph from the training samples of a dataset
 * @author Adan Echemendia
 * @date May 17, 2016
 * @param Z Input dataset
 * @param kmax Maximum number of neighbors possible for each node
 * @return A KNN Graph
 *
 */
//! swig(newobject, stable)
iftKnnGraph *iftCreateKnnGraph(iftDataSet *Z, int kmax);

/**
 * @brief Fast method to create a KNN Graph from the training samples of a dataset ONLY when they are in 1D and of type integer. It first sorts the 1D data and then choose the k nearest neighbors in the ordered data
 * @author Cesar Castelo
 * @date Nov 28, 2017
 * @param Z Input dataset
 * @param kmax Maximum number of neighbors possible for each node
 * @return A KNN Graph
 *
 */
iftKnnGraph *iftCreateKnnGraphInt1D(iftDataSet *Z, int kmax);

/**
 * @brief Destroys a KNN Graph
 * @param graph KNN Graph to destroy
 */
void iftDestroyKnnGraph(iftKnnGraph **graph);


/**
 * @brief Read a KNN Graph.
 * @author Samuka Martins
 * @date May 15, 2018
 */
 //! swig()
iftKnnGraph *iftReadKnnGraph(const char *format, ...);


/**
 * @brief Write a KNN Graph.
 * @author Samuka Martins
 * @date May 15, 2018
 */
//! swig()
void iftWriteKnnGraph(const iftKnnGraph *graph, const char *path, ...);


/**
 * @brief It selects a maximum given number of labeled training
 * samples and removes potential outliers from that training
 * set. These outliers are identified as samples that fall into
 * clusters, as computed for a pdf with many domes, whose roots are
 * from different classes. If all samples are removed from any given
 * class, then the set is considered invalid for training. Otherwise,
 * the remaining training samples are considered valid.
 *
 * @author Alexandre Falcao
 * @date Jul 2, 2016
 * @param Input data set for training sample selection.
 * @param Maximum number of training samples.
 * @return Boolean decision about the vality of the data set for training.
 *
 * @note This function needs to be improved in that the Kmax should be
 * selected according to the problem at hand. If it is too high (e.g.,
 * 5 for some problems) it might eliminate too many informative
 * samples.
 */
bool iftSelectSupTrainSamplesWithNoOutliers(iftDataSet *Z, int max_ntrainsamples);

/**
 * @brief Clusters the samples of a KNN Graph by OPF
 * @details First, it finds the best k value (adjacency parameter) in the KNN Graph minimizing the graph cut function, given as a
 * parameter, in the interval [1, graph->kmax]. After that, computes the clusters with OPF.
 * @author Adan Echemendia
 * @date May 17, 2016
 * @param graph Input KNN Graph
 * @param iftGraphCutFun Graph cut function
 * @return Resulting number of clusters
 */
//! swig(stable)
int  iftUnsupTrain(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);

/**
 * @brief Clusters the samples of a KNN Graph by OPF. The difference with the function iftUnsupTrain is that this function does a local search, instead an exhaustive search, to find the adjacency parameter K
 * @author Adan Echemendia
 * @date May 17, 2016
 * @param graph Input KNN Graph
 * @param iftGraphCutFun Graph cut function
 * @return Resulting number of clusters
 */
int iftFastUnsupTrain(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);

/**
 * @brief Clusters the samples of a KNN Graph in a specified number of groups by OPF.
 * @details The returned value is the number of clusters given as result by OPF. This number is then adjusted to
 * remain with the specified number of groups.
 * @author Adan Echemendia
 * @date May 17, 2016
 * @param graph Input KNN Graph
 * @param c Desired number of clusters
 * @return Original number of groups given by OPF.
 */
int  iftUnsupTrainWithCClusters(iftKnnGraph *graph, int c);

/**
 * @brief Clusters the samples of a KNN Graph in a specified number of groups by OPF.
 * This version searchs for the best k using a binary search strategy.
 * @details The returned value is the number of clusters given as result by OPF. This number is then adjusted to
 * remain with the specified number of groups.
 * @author Cesar Castelo
 * @date Feb 20, 2018
 * @param graph Input KNN Graph
 * @param c Desired number of clusters
 * @return Original number of groups given by OPF.
 */
int iftFastUnsupTrainWithCClustersBinarySearch(iftKnnGraph *graph, int c);

/**
 * @brief Propagates group labels in the samples of a dataset given a clustered graph
 * @details [long description]
 *
 * @param graph Input clustered graph
 * @param Z Input ataset
 *
 * @return Number of outliers in the dataset
 */
int iftUnsupClassify(iftKnnGraph *graph, iftDataSet *Z);


 /**
 * @brief Propagates the truelabel of the group representatives to another dataset 
 * @details [long description]
 *
 * @param graph Input clustered graph
 * @param Z Input dataset
 *
 * @return Number of outliers in the dataset
 */
 
  int iftSupOPFKnnClassify(iftKnnGraph *graph, iftDataSet *Z);

  
/**
 * @brief Propagates labels to a dataset given a clustered graph. If a sample of <b>Z</b> is recognized as outlier it gets the label 0 in the dataset. This outlier treatment is the difference of this function with the function <b>iftUnsupClassify</b>
 * @details [long description]
 * @author Adan Echemendia
 * @param graph Input clustered graph
 * @param Z Input ataset
 *
 * @return Number of outliers in the dataset
 */
int iftUnsupClassifyDataSetLessOutliers(iftKnnGraph *graph, iftDataSet *Z);


/**
 * @brief Propagates labels from a clustered graph (representative sample set) to a large dataset
 * @details At first, the function propagates the labels using a clustered graph that contains a little initial set.
 * If the process results in a high number of outliers, that means that the "representative sample set" is not good enough and in the next iteration, some of these outliers will be included in the "representative sample set". The process continues until reach a maximum number of iterations or detect a low number of outliers. Note that at the end, the method is trying to learn the PDF of a large dataset from a little sample set.
 *
 * @param Z Large input dataset
 * @param kmax_perc The maximum number of neighbors possible for each node will be this percent value of the number of training samples in the dataset
 * @param iftGraphCutFun Graph cut function
 * @param niterations Maximum number of iterations
 * @return A clustered graph, composed by a representative sample set of a large dataset
 */
iftKnnGraph *iftUnsupLearn(iftDataSet *Z, float kmax_perc, iftKnnGraphCutFun iftGraphCutFun, int niterations,char debug);

/**
 * @brief Propagates labels from a clustered graph (representative sample set) to a large dataset. This function tries to improve the
 * <b>iftUnsupLearn</b> function differentiating between real outliers and false positive outliers.
 * @details After the clustering and propagation phases, it's created a dataset with the resulting outliers. Now, instead of adding to the training set any outlier, we are trying to differentiate among false positive and real outliers. We are interested in adding only the last ones to the next training set. To try to make this distinction we are going to cluster the dataset of outliers. Assuming that there must be approximately (total_outliers/cluster_number) outliers in each group, if a cluster has less than 0.5*(total_outliers/cluster_number) outliers, we consider all samples belonging to it as real outliers. These samples aren't added to the training set. All remaining samples are considered false positive outliers and are added to the training set for the next iteration
 * @param Z Large input dataset
 * @param kmax_perc The maximum number of neighbors possible for each node will be this percent value of the number of training samples in the dataset
 * @param iftGraphCutFun Graph cut function
 * @param niterations Maximum number of iterations
 * @return A clustered graph, composed by a representative sample set of a large dataset
 */
iftKnnGraph *iftUnsupLearn2(iftDataSet *Z, float kmax_perc, iftKnnGraphCutFun iftGraphCutFun, int niterations, char
debug);

/**
 * @brief Clusters a big dataset using a divide and conquer strategy with 2 levels and the unsupervised OPF method.
 * @note The samples of each partition are selected randomly without repetition.
 * @param dataset A big dataset to be clustered
 * @param nb_partitions Number of partitions to divide the big dataset
 * @param nb_train_samples Number of train samples for each partition. If this number is 0, we use all samples in the
 * partion as training samples
 * @param iftGraphCutFun A graph cut function
 * @param kmax_perc1 k max percent in the first level
 * @param kmax_perc2 k max percent in the second level
 * @param debug 1 to print debug information 0 otherwise
 */
void iftDivideAndConquerUnsupOPF2Levels(iftDataSet *dataset, int nb_partitions, int nb_train_samples, iftKnnGraphCutFun
                                        iftGraphCutFun,
                                        float kmax_perc1, float kmax_perc2, char debug);


/**
 * @brief Clusters a big image using a divide and conquer strategy, tile division, and the unsupervised OPF method. The hierarchy has two levels. The first level is formed by partitioning the image into blocks and clustering each one by OPF. The second level is composed by the cluster representatives in the first level. The samples of the second level are clustered, and the result is
 * propagated downwards to all samples of the original dataset.
 * @param orig_dataset A big dataset (that references an image) to be clustered
 * @param num_blocks Number of tiles to partition the image
 * @param train_percent Percent of train samples for each partition
 * @param iftGraphCutFun A graph cut function
 * @param kmax_perc1 KMax percent in the first level
 * @param kmax_perc2 KMax percent in the second level
 * @param debug 1 if we want debug information, 0 otherwise
 * @returns The KnnGraph resulting of the clustering of the second level
 *
 */
iftKnnGraph *iftClusterImageByDivideAndConquerUnsupOPF2Levels(iftDataSet *orig_dataset, int num_blocks,
                                                              float train_percent,
                                                              iftKnnGraphCutFun iftGraphCutFun, float kmax_perc1,
                                                              float kmax_perc2, char debug);

/**
 * @brief Clusters a big image using a divide and conquer strategy, tile division and the unsupervised OPF method. The hierarchy has only one level. This function differs with the <b>iftDivideAndConquerRandomlyUnsupOPF2Levels</b> function that here we are clustering each partition independently and there is not a second level to join cluster prototypes. Two adjacent clusters, found in different blocks, are joined if they have similar color information.
 * @param orig_dataset A big dataset (that references an image) to be clustered
 * @param num_blocks Number of tiles to partition the image
 * @param train_percent Percent of train samples for each partition
 * @param iftGraphCutFun A graph cut function
 * @param k_max_percent The maximum number of neighbors possible for each node will be this percent value of the number of training samples in the each partitioned dataset
 * @param debug 1 if we want debug information, 0 otherwise
 * @returns The number of clusters found
 */
int iftClusterImageByDivideAndConquerUnsupOPF1Level(iftDataSet *orig_dataset, int num_blocks, float train_percent,
                                                    iftKnnGraphCutFun iftGraphCutFun, float k_max_percent, char debug);

/**
 * @brief Propagates labels to the samples of a dataset <b>Z</b> using a semi-supervised strategy
 * @param graph Supervised clustered graph
 * @param Z Input dataset
 * @return Number of outliers in the dataset Z
 */
int  iftSemiSupClassify(iftKnnGraph *graph, iftDataSet *Z);

/**
 * @brief Selects the best features in a dataset by multi-scale parameter search
 * @details [long description]
 *
 * @param Z Input dataset
 * @param kmax Maximum number of neighbors possible for each node
 */
void  iftUnsupFeatSelecByMSPS(iftDataSet *Z, int kmax);

/**
 * @brief Computes a normalized cut in a clustered graph
 * @details [long description]
 *
 * @param graph Input clustered graph
 * @return The normalized cut value
 */
//! swig(newobject)
float  iftNormalizedCut(iftKnnGraph *graph);

/**
 * @brief Finds the best <b>k</b> value (adjacency parameter) for a graph of samples given a cut function
 * @details For each <b>k</b> in the interval [1, graph->kmax], the function computes the PDF of all nodes, groups
 * them by OPF, and determines the cut function value. It returns the <b>k</b> value that minimizes the aforementioned
 * cut function value.
 * @param graph Graph containing the samples to cluster
 * @param iftGraphCutFun Cut Function
 * @return Optimum k value for the clustered graph
 */
int  iftBestkByKnnGraphCut(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);

/**
 * @brief Finds the best <br>k</br> value (adjacency parameter) for a graph of samples given a cut function. In this function we are not doing an exhaustive search into the interval [0,kmax] to find <b>k</b>. The idea is do a local search from [kmax..1] and stopping when the optimization function stops growing. We assume that this value of <b>k</br> is good enough, so if we need a closer
 * view of the dataset we can reduce the <b>kmax</b> parameter
 * iftBestkByKnnGraphCut
 * @param graph KNN Graph
 * @param iftGraphCutFun Cut Function
 * @return A local optimum for the adjacency parameter <b>k</b>
 */
int iftFastBestKByKnnGraphCut(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);

/**
 * @brief Builds a new dataset with the boundary samples of an input clustered dataset
 * @details [long description]
 *
 * @param graph Input clustered graph
 * @param Z Input clustered dataset
 *
 * @return A new dataset with the boundary samples of the clustered dataset Z
 */
iftDataSet  *iftGraphBoundaryReduction(iftKnnGraph *graph, iftDataSet *Z);

/**
 * @brief Gets the boundary samples from the dataset contained in a clustered graph
 * @details [long description]
 *
 * @param graph Input clustered graph
 * @return Set of boundary samples from the dataset contained in a clustered graph
 */
iftSet  *iftGetKnnBoundarySamples(iftKnnGraph *graph);

/**
 * @brief Gets the samples corresponding to the root nodes of a clustered graph
 * @details [long description]
 *
 * @param graph Input clustered graph
 * @return The sample set corresponding to the root nodes of a clustered graph
 */
iftSet  *iftGetKnnRootSamples(iftKnnGraph *graph);

/**
 * @brief Computes an array containing, in each position, a list of sample indexes of each cluster sorted by their distance to the root node.
 * @details The array has size Z->nlabels+1 containing, in each position, a list of sample indexes of each cluster sorted by their distance
 * to the root node (increasing distance). The root node is the first element of each list. This function assumes that each cluster
 * contains only one root node and labels start at 1.
 *
 * @param graph Input clustered graph
 * @param nNodes Output array of size Z->nlabels+1 with the number of nodes in each cluster
 *
 * @return Array that contains in each position a list of sample indexes of each cluster, sorted by their distance to the root node
 */
int  **iftGetNodesPerClusterOrderedByDistanceToRoot(iftKnnGraph *graph, int **nNodes);

/**
 * @brief Computes an array containing, in each position, a list of sample indexes of each cluster sorted by their decreasing weight
 * @details The array has size Z->nlabels+1 containing, in each position, a list of sample indexes of each cluster sorted by their
 * decreasing weight. The root node is the first element of each list. This function assumes that each cluster contains only one root
 * node and labels start at 1.
 *
 * @param graph Input clustered graph
 * @param nNodes Output array of size Z->nlabels+1 with the number of nodes in each cluster
 *
 * @return Array that contains in each position a list of sample indexes of each cluster, sorted by their decreasing weight
 */
int  **iftGetNodesPerClusterOrderedByWeight(iftKnnGraph *graph, int **nNodes);

/**
 * @brief Propagates, for each cluster of the graph, the true label of the corresponding root
 * @details [long description]
 *
 * @param graph Input clustered graph
 */
void  iftPropagateClusterTrueLabels(iftKnnGraph *graph);

/**
 * @brief Propagates, for each cluster of the dataset, the true label of the corresponding prototype. The prototypes
 * are identified because they have a status value equal to IFT_PROTOTYPE.
 * @author Adan Echemendia
 */
void  iftPropagateClusterTrueLabels2(iftDataSet *Z);

/**
 * @brief Propagates, for each cluster of the dataset, the true label more common in it.
 * @author Adan Echemendia
 */
void  iftPropagateClusterTrueLabels3(iftDataSet *Z);

void  iftPDFByRange(iftKnnGraph *graph);

int   iftUnsupOPF(iftKnnGraph *graph);

/**
 * @brief Extract the centroids from a previously clustered iftDataSet. If realCentroids is TRUE, it returns the closest
 * elements to the computed centroids (mean vectors), if FALSE it returns the computed centroids (mean vectors of each group)
 * @param Z Dataset that was previously clustered
 * @param returnRealCent Whether to return the real centroids or the mean vector of the group instead
 * @param usePrototypes Whether to return the prototypes of each group as the centroids instead of compute them (this SHOULD be
 * used as TRUE only if the dataset was clustered with OPF)
 * @return A dataset containing the centroids
 * @author Cesar Castelo
 * @date Abr 23, 2018
 */
//! swig(newobject)
iftDataSet* iftExtractCentroidsFromDataSetAsDataSet(iftDataSet *Z, bool returnRealCent, bool usePrototypes);

/**
 * @brief Extracts the centroids from a previously clustered iftDataSet. If realCentroids is TRUE, it returns the closest
 * elements to the computed centroids (mean vectors), if FALSE it returns the computed centroids (mean vectors of each group)
 * @param Z Dataset that was previously clustered
 * @param returnRealCent Whether to return the real centroids or the mean vector of the group instead
 * @param usePrototypes Whether to return the prototypes of each group as the centroids instead of compute them (this SHOULD be
 * used as TRUE only if the dataset was clustered with OPF)
 * @return A matrix containing the centroids
 * @author Cesar Castelo
 * @date Abr 23, 2018
 */
iftMatrix* iftExtractCentroidsFromDataSetAsMatrix(iftDataSet *Z, bool returnRealCent, bool usePrototypes);


/**
 * @brief Creates a KNN graph measuring the distance according to the arc-weights on the MST graph
 * @param mst           minimum spanning tree
 * @param number_neigh  number of neighbors on graph
 * @return Knn graph for clustering
 * @author Jordão Bragantini
 * @date Aug, 2018
 */

//! swig(newobject)
iftKnnGraph *iftMSTtoKnnGraph(iftMST* mst, int number_neigh);

/**
 * @brief Transforms the nearest neighbor's features of a given node into a matrix
 * @param graph        Knn graph
 * @param node_index   Node whose neighbors will be considered
 * @return matrix with the neighbors as rows and features as columns
 * @author Azael Sousa
 * @date Dec 03, 2019
 */
iftMatrix *iftKnnNodeToMatrix(iftKnnGraph *graph, int node_index);

/**
 * @brief For a given node, computes its neighbor's index in the dataset
 * @param graph        Knn graph
 * @param node_index   Node whose neighbors will be considered
 * @return An IntArray with the sample indexes of the respective neighbors
 * @author Azael Sousa
 * @date Dec 03, 2019
 */
iftIntArray *iftGetKnnNeighborSamplesIndex(iftKnnGraph *graph, int node_index);

//void iftJoinClusterByImageAdjacency(iftDataSet *Z, iftImage* suppxl, int k);

/**
 * @brief Computes the class purity of the clusters for each class. It
 * can also compute the clusters based on the graph->k nearest
 * neighbors. It requires an annotated dataset in graph->Z. 
 * @param graph --- input knn graph 
 * @param compute_OPF --- input option to either compute or not the clusters. 
 * @return An FloatArray with the average purity followed by the purity per class. 
 * @author Alexandre Falcão 
 * @date Dec 31st, 2019
 */
  
iftFloatArray *iftClusterPurity(iftKnnGraph   *graph, bool compute_OPF);

#ifdef __cplusplus
}
#endif

#endif
