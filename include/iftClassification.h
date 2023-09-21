/**
 * @file iftClassification.h
 * @brief Definitions and functions about Classification algorithms.
 */

#ifndef IFT_CLASSIFICATION_H
#define IFT_CLASSIFICATION_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/FHeap.h"

#include "iftCommon.h"
#include "iftClustering.h"
#include "iftDataSet.h"
#include "iftMetrics.h"
#include "iftMST.h"
#include "iftDeepLearning.h"
#include "iftSVM.h"

#include "iftArgumentList.h"
#include "iftDistanceFunctions.h"
#include "iftGenericVector.h"



//void srand(iftRandomLinearCongruentialGenerator randomLCG, unsigned int seed)
//{
//    randomLCG.X = seed;
//}

typedef enum ift_Classifier_type {
	CLASSIFIER_CUSTOM,
	CLASSIFIER_SVM,
	CLASSIFIER_OPF_SUP,
    CLASSIFIER_KNN,
	CLASSIFIER_RANDOM
} iftClassifierType;



typedef void (*iftTrainFunction) (void* classifierObj, iftDataSet* dataSet);
typedef void (*iftPredictFunction) (void* classifierObj, iftDataSet* dataSet);
// typedef void (*FreeFunction) (void* object);

typedef struct _ift_genericClassifier {
	void * classifier;
	FreeFunction freeFunctionClassifier;
	iftTrainFunction trainFunction;
	iftPredictFunction predictFunction;
	const char* classifierName;
}iftGenericClassifier;

iftGenericClassifier* iftCreateGenericClassifier(iftClassifierType classifierType);
void iftTrainGenericClassifier(iftGenericClassifier *genericClassifier, iftDataSet *dataSet);
void iftPredictGenericClassifier(iftGenericClassifier* genericClassifier, iftDataSet* dataSet);
void iftDestroyGenericClassifier(iftGenericClassifier** pGenericClassifier );


/////////////SVM classifier
iftSVM* iftCreateSVMClassifier();
void iftTrainSvmVoidPointer(void *svmClassifier, iftDataSet *dataSet);
void iftPredictSvmVoidPointer(void *svmClassifier, iftDataSet *dataSet);
void iftDestroyClassifierVoidPointerSvm(void* psvmClassifier);
////////////////

/**
 * Defines a node in a complete graph
 */
typedef struct ift_cplnode {
	int   pred;     // predecessor node
	int   sample;   // corresponding training sample in Z
	int	  length;
} iftCplNode;

/**
 * Defines a complete graph
 */
//! swig(extend = iftCplGraphExt.i, destroyer = iftDestroyCplGraph)
typedef struct ift_cplgraph {
	iftCplNode *node;          // list of nodes
	float      *pathval;       // path value for each node
	int        *ordered_nodes; // list of nodes ordered by path value
	int         nnodes;        // number of nodes
	iftFHeap   *Q;             // priority queue
	iftDataSet *Z;             // Each graph node is a training sample in Z
} iftCplGraph;


/////////////////////////////////OPF-SUP classifier
typedef struct _ift_OPF_SUP {
	iftCplGraph* completeGraph;
	iftDistanceTable* distanceTable;
	iftSet* protoTypes;
    bool recomputeDistanceTable;
}iftOPFSUPClassifier;

iftOPFSUPClassifier* iftCreateOPFSUPClassifier();
void iftDestroyClassifierOpfsup(iftOPFSUPClassifier** pOpfsupClassifier);
void iftDestroyClassifierVoidPointerOpfsup(void* pOpfsupClassifier);
void iftTrainOpfsup(iftOPFSUPClassifier* opfsup, iftDataSet* dataSet);
void iftTrainOpfsupVoidPointer(void *opfsupClassifier, iftDataSet *dataSet);
void iftPredictOpfsup(iftOPFSUPClassifier* opfsupClassifier, iftDataSet* dataSet);
void iftPredictOpfsupVoidPointer(void *opfsupClassifier, iftDataSet *dataSet);
////////////////////////////////////

typedef struct _ift_OPF_SEMI_SUP {
    iftDataSet *datasetTrain;
    iftMST *minimumSpanningTree;
    iftCplGraph *completeGraph  ;
    iftSet *adj ;
}iftOpfSemiSupClassifier;


/////////////////////////////////random classifier
typedef struct _ift_Random_Classifier {
	double* classProbabilities;
	size_t numberOfClasses;
	//int seed;
}iftRandomClassifier;

iftRandomClassifier* iftCreateRandomClassifier();
void iftDestroyRandomCliftActiveLearningStateiftActiveLearningStateassifier(iftRandomClassifier** pRandomClassifier);
void iftDestroyClassifierVoidPointerRandom(void* pRandomClassifier);
void iftTrainRandomClassifier(iftRandomClassifier* randomClassifier, iftDataSet* dataSet);
void iftTrainRandomClassifierGeneric(void* randomClassifier, iftDataSet* dataSet);
void iftPredictRandomClassifier(iftRandomClassifier* randomClassifier, iftDataSet* dataSet);
void iftPredictRandomClassifierGeneric(void* randomClassifier, iftDataSet* dataSet);
////////////////////////////////////


/////////////////////////////////knn classifier
typedef struct _ift_knn_Classifier {
	iftDataSet* datasetTrain;
	iftDistanceFunction distanceFunction;
    bool copyfeatsOfTraiSamples;
	iftArgumentList* distanceFunctionParameters;
	size_t kNearest;
}iftKnnClassifier;

iftKnnClassifier* iftCreateKnnClassifier();
void iftDestroyKnnClassifier(iftKnnClassifier** pKnnClassifier);
void iftDestroyClassifierVoidPointerKnn(void* pKnnClassifier);
void iftTrainKnnClassifier(iftKnnClassifier* knnClassifier, iftDataSet* dataSet);
void iftTrainKnnClassifierVoidPointer(void* knnClassifier, iftDataSet* dataSet);
void iftPredictKnnClassifier(iftKnnClassifier* knnClassifier, iftDataSet* dataSet);
void iftPredictKnnClassifierVoidPointer(void* knnClassifier, iftDataSet* dataSet);
////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////

/**
 * Defines a logistic regression classifier
 */
typedef struct ift_logreg {
	int nclasses; /* number of classes */
	int nfeats; 	/* number of features*/
	float error;	/* error*/
	iftMatrix* coef; /* coefficient matrix*/
} iftLogReg;

//LOGREG

////////////

/**
 * @brief Trains a logistic regressor over the IFT_TRAIN samples in the dataset Z.
 * @author Peixinho
 * @date March 2016
 */
iftLogReg* iftLogRegTrain(iftDataSet* Z);
/**
 * @brief Classifies the IFT_TEST samples in the dataset Z, using the logistic regressor.
 * @author Peixinho
 * @date March 2016
 */
void iftLogRegClassify(iftLogReg* logreg, iftDataSet* Z);

/**
 * @brief Destroys logistic regressor classifier.
 * @author Peixinho
 * @date March 2016
 */
void iftDestroyLogReg(iftLogReg** logreg);

/**
 * @brief Reads a logistic regressor classifier from a file
 * @param pathname The pathname of the file
 * @return A logistic regressor classifier
 */
iftLogReg *iftReadLogReg(const char *pathname);

/**
 * @brief Writes a logistic regressor classifer to a file
 * @param logreg A logistic regressor classifier
 * @param pathname The pathname of the file
 */
void iftWriteLogReg(const iftLogReg *logreg, const char *pathname);

/**
 * @brief Writes a complete graph to a file
 * @param graph A complete graph
 * @param pathname The pathname of the file
 */
//! swig()
void iftWriteCplGraph(const iftCplGraph *graph, const char *pathname);

/**
 * @brief Reads a complete graph from a file
 * @param pathname Pathname of the file
 * @return A complete graph
 */
 //! swig()
iftCplGraph* iftReadCplGraph(const char* pathname);

/**
 * @brief Creates a complete graph from the training samples of the dataset <b>Z</b>
 * @param Z A dataset
 * @return The complete graph
 */
//! swig(newobject, stable)
iftCplGraph *iftCreateCplGraph(iftDataSet *Z);

/**
 * @brief Destroys a complete graph
 * @param graph A complete graph
 */
void         iftDestroyCplGraph(iftCplGraph **graph);
void         iftSupTrain2(iftCplGraph *graph);
//! swig(stable)
void         iftSupTrain(iftCplGraph *graph);
iftCplGraph *iftSupLearn(iftDataSet *Z);


/**
 * @brief Classifies testing samples (status = IFT_TEST) from the dataset Ztest by using
 * Supervised OPF.
 * @author Samuel Martins
 * @date Jun 14, 2016
 *
 * @note The min cost of the testing sample after classification is stored in its field <b>weight</b>
 * in the dataset.
 * @note If a testing sample is classified incorrectly, its status will be IFT_ERROR.
 * 
 * @param  graph Supervised OPF (Complete Graph).
 * @param  Ztest Dataset with the testing samples to be classified.
 * @return       Number of classification errors.
 */
//! swig(stable)
int iftClassify(const iftCplGraph *graph, iftDataSet *Ztest);


/**
 * @brief Classify by One Class OPF.
 *
 * If a testing sample is inside the adjacency radius of any training sample, it belongs this training sample's cluster
 * and its label is IFT_POSITIVE_CLASS
 * Otherwise, it is an outlier (status is IFT_OUTLIER) and its predicted label is IFT_NEGATIVE_CLASS.
 * 
 * Each training sample s has its own adjacency radius which corresponds to its distance to its i-th nearest
 * neighbor (i in [1, graph->k]) which is defined as:
 * 
 * i = ceil(quantile_of_k * graph->k)
 * 
 * where graph->k is the number of neighbors of each training sample, and 
 * quantile_of_k is the quantile used to define the rank of i.
 * 
 * quantile_of_k = 0.0  ==> we consider the nearest neighbor,
 * quantile_of_k = 0.5  ==> the median value,
 * quantile_of_k = 1.0  ==> the k-th nearest neighbor (the farthest one).
 * 
 * @note Only the samples with status IFT_TEST in Ztest will be classified.
 * 
 * @param Ztest Dataset for classification (only the samples with status IFT_TEST will be classified) 
 * @param graph Trained KnnGraph (by iftUnsupTrain) with only the training samples of the target class.
 * @param quantile_of_k Quantile used to find the rank of the i-th nearest neighbor that defines the adjacency radius
 * of each training sample. It must be in [0, 1]
 * 
 * @author Samuka Martins
 * @date Dec 3, 2019
 */
//! swig(newobject)
void iftClassifyByOneClassOPF(iftDataSet *Ztest, const iftKnnGraph *graph,
                              float quantile_of_k);



/**
 * Classifies testing samples (status = IFT_TEST) from the dataset Ztest by using
 * Supervised OPF and computes a certainty rate for each testing sample. 
 * @author Samuel Martins
 * @date Jun 14, 2016
 *
 * Similar to @ref iftClassify(), this function also computes the
 * certainty value that the test samples belong to their predicted
 * labels. \n For that, the test sample is initially "conquered" by
 * the train. sample that offered the minimum cost
 * <b>min_cost1</b>. Suppose that its label is <b>A</b>. \n Then, the
 * function searches another class/truelabel different from <b>A</b> whose training
 * sample offers the minimum cost <b>min_cost2</b> to the testing
 * sample. \n Finally, the certainty value is computed by: \n
 * <code><b>certainty = (min_cost2) / (min_cost1 +
 * min_cost2)</b></code>
 *
 * @note The certainty value of the testing sample after classification is stored as its <b>weight</b>
 * in the dataset.
 * @note If a testing sample is classified incorrectly, its status will be IFT_ERROR.
 * 
 * @param  graph Supervised OPF (Complete Graph).
 * @param  Ztest Dataset with the testing samples to be classified.
 * @return Number of classification errors.
 */
 //! swig(stable)
int          iftClassifyWithCertaintyValues(const iftCplGraph *graph, iftDataSet *Ztest);

void         iftBorderClassify(iftCplGraph *graph, iftDataSet *Z,int truelabel);
iftSet      *iftPrototypesByMST(iftCplGraph *graph);
void         iftSupFeatSelection(iftDataSet *Z, float train_perc, int nteams);
int iftClassifyByKNN(iftDataSet *Z1, iftDataSet *Z, int k);
iftSet      *iftGetMSTBoundarySamples(iftCplGraph *graph, iftMST *mst, int nsamples);
float        iftGetKFoldCrossValidationMeanOPFAccuracy(iftDataSet *Z, const int k);
iftSet      *iftGetRootDistanceBasedSamples(iftCplGraph *trainGraph, iftDataSet *Z,
											int **sampleLists, int *nSamples, uchar status, int statusCond,
											int nSelect);
float iftFindMinimumCostAndClassifySample(const iftCplGraph *graph, iftDataSet *Ztest, int t);

/** @brief This function performs OPF's semi-supervised algorithm directly on a complete graph. 
  *
  * @author Thiago V. Spina and A.X. Falcao
  * @date November, 14th 2017
*/
//! swig(newobject, stable)
iftCplGraph* iftSemiSupTrain(iftDataSet *Ztrain);



/**
 * @brief   Trains an Image (Texture/Brightness) OPF Classifier from the voxel brightness from a
 *          a set of markers over the training image.
 * @author  Samuel Martins
 * @date    Aug 21, 2016
 * @ingroup Classification
 *
 * @param train_img     The training image.
 * @param train_markers The training markers used to get the voxel brightness from the Training Image.
 * @param smooth_factor A smooth factor. Use 0.0 to ignore it.
 *
 * @return An Image (Texture/Brightness) OPF Classifier.
 */
iftCplGraph *iftTrainImageClassifierByOPF(const iftImage *train_img, const iftLabeledSet *train_markers,
										  float smooth_factor);



/**
 * @brief Classifies an Image by a Texture OPF Classifier.
 *
 * The classifier can be multi-label, however, the output fuzzy image is an only with resulting of the
 * "merge" of all fuzzinesses from the non-background voxels (voxel != 0).
 * 
 * @param  test_img      Image to be classified.
 * @param  clf           OPF Classifier.
 * @param  smooth_factor Factor to smooth the input image. Use 0.0 to ignore the smoothness.
 * @param  max_range     Maximum range (value) of the images. Ex: For 8-bit image --> 255
 * @param  fuzzy_map     Return (if != NULL) the resulting fuzzy map.
 * @return               Classified Image.
 * 
 * @author Samuel Martins
 * @date Sep 8, 2016
 * @ingroup Classification
 */
iftImage *iftClassifyImageByOPF(const iftImage *test_img, const iftCplGraph *clf, float smooth_factor,
								int max_range, iftImage **fuzzy_map);



#ifdef __cplusplus
}
#endif

#endif
