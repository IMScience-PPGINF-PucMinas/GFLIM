//
// Created by deangeli on 11/08/17.
//

#ifndef IFT_IFTACTIVELEARNING_H
#define IFT_IFTACTIVELEARNING_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftGenericVector.h"
#include "iftDataSet.h"
#include "iftClassification.h"
#include "iftClustering.h"

/////////////////////////////////SELECTOR//////////////////////////////////////////////////
typedef void (*iftPreProcessingSelectorFunction) (void* selectorObject);
typedef iftIntArray* (*iftSelecSamplesSelectorFunction) (void* selectorObject);

typedef struct _iftGenericSelector {
    void* selector;
    const char* selectorName;
    FreeFunction freeFunctionSelector;

    iftPreProcessingSelectorFunction preProcessingSelectorFunction;
    iftSelecSamplesSelectorFunction selecSamplesSelectorFunction;
    iftIntArray* selectedSamplesIndex;
    iftGenericVector* preProcessingOuput;
    iftGenericVector* selectSamplesOutput;

}iftGenericSelector;

iftGenericSelector* iftCreateGenericSelector();
void iftDestroyGenericSelector(iftGenericSelector** pGenericSelector);
void iftPreProcessGenericSelector(iftGenericSelector* genericSelector);
iftIntArray* iftSelectSamplesGenericSelector(iftGenericSelector* genericSelector);



//////////RDS SELECTOR
typedef struct _iftSelectorRDS iftSelectorRDS;
typedef void (*iftFindReference2SampleRDSFunction) (iftSample* sample, iftSelectorRDS* rds);

struct _iftSelectorRDS {
    iftDataSet *dataSet;
    bool destroyDataset;
    iftKnnGraph* opfClusteringGrapgh;
    bool createAndComputeGraph;
    float** samplesDistance2Root;
    bool destroyGraph;
    float OPFUNSUP_kmaxPercentage;
    int numberClusters;
    int numberOfSelectedSamplesPerCluster;
    int currentLearningCycle;
    char* path2SamplesFolder;
    char currentFileName[100];
    char fileExtension[10];

    iftFindReference2SampleRDSFunction findReference2SampleRDSFunction;

    iftGenericVector* clusterOfSamplesIndex; //vector <  vector<int>  >  || all clusters -> clusterIndex -> sampleIndex
    iftGenericVector* selectedSamples; //vector<int>  ||  sampleIndex
    iftGenericVector* selectedSamplesPerCycle; //vector <  vector<int>  >  || all lists -> list at iteration -> sampleIndex
    bool autoSetStatus;
    iftKnnGraphCutFun iftGraphCutFun;

    iftIntArray* selectedSamplesIndex_typeSafe;
    //void* classifier;
    bool destroyClassifier;
    iftGenericClassifier* classifier;
};

iftSelectorRDS* iftCreateSelectorRDS(iftDataSet *dataSet);
void iftDestroySelectorRDS(iftSelectorRDS** selectorRDS);
void iftDestroySelectorRDSVoidPointer(void* selectorRDSvp);

void iftPreProcessSelectorRDS(iftSelectorRDS* selectorRDS);
void iftPreProcessSelectorRDSVoidPointer(void *selectorRDSVp);

void iftSelecSamplesByRDS(iftSelectorRDS* selectorRDS);
iftIntArray* iftSelecSamplesByRDSVoidPointer(void *selectorRDSvp);

void iftPrintSelectorInfoRDS(iftSelectorRDS *selectorRDS, bool printSampleName);
void iftPrintSelectedSamplesInfoRDS(iftSelectorRDS* selectorRDS);

void iftSuperviseSelectedSamplesByRDS(iftSelectorRDS* selectorRDS);
int  iftUnsupTrainForRDS(iftKnnGraph *graph,int nclasses);

////////////////

//////////RANDOM SELECTOR
typedef struct _iftSelectorRandom {
    iftDataSet *dataSet;
    iftGenericVector* selectedSamples; // vector<int>
    int* shuffledIndices;
    int shuffledIndicesSize;
    int numberOfSelectedSamples;
    bool clearDatasetWhenDestroyed;
    iftIntArray* selectedSamplesIndex_typeSafe;
}iftSelectorRandom;

iftSelectorRandom* iftCreateSelectorRandom(iftDataSet *dataSet);
void iftDestroySelectorRandom(iftSelectorRandom** selectorRandom);
void iftDestroySelectorRandomVoidPointer(void* selectorRandomvp);
void iftSelecSamplesByRandom(iftSelectorRandom* selectorRandom);
iftIntArray* iftSelecSamplesByRandomVoidPointer(void* selectorRandomvp);
////////////////////////

//////////RWS SELECTOR
typedef struct _iftSelectorRWS {
    iftDataSet *dataSet;
    bool destroyDataset;

    iftKnnGraph* opfClusteringGrapgh;
    bool createAndComputeGraph;
    bool destroyGraph;
    iftKnnGraphCutFun iftGraphCutFun;


    bool destroyClassifier;
    iftGenericClassifier* classifier;
    iftClassifierType classifierType;
    bool autoSetStatus;
    float OPF_kmaxPercentage;
    int numberClusters;
    iftGenericVector* clusterOfSamplesIndex; //vector <  vector<int>  >  || all clusters -> clusterIndex -> sampleIndex
    int numberOfSelectedSamples;
    int numberOfSelectedSamplesPerCluster;
    iftIntArray* selectedSamplesIndex_typeSafe;

}iftSelectorRWS;

iftSelectorRWS* iftCreateSelectorRWS(iftDataSet *dataSet);
void iftPreProcessSelectorRWS(iftSelectorRWS* selectorRWS);
void iftSelecSamplesByRWS(iftSelectorRWS* selectorRWS);
iftIntArray* iftSelecSamplesByRWSVoidPointer(void *selectorRDSvp);
///////////////////////////////////


///////////////////////ACTIVE LEARNING//////////////////////////////
//typedef iftActiveLearningState (*iftActiveLearningStateFunction) (void* activeLearningObject);
//
//typedef struct _iftGenericActiveLearning {
//    void* activeLeaner;
//    const char* activeLearnerName;
//    FreeFunction freeFunctionActiveLeaner;
//    iftActiveLearningState currentState;
//
//
//    iftActiveLearningStateFunction preProcessingActiveLearningFunction;
//    iftActiveLearningStateFunction selecSamplesActiveLearningFunction;
//    iftActiveLearningStateFunction supervisingActiveLearningFunction;
//    iftActiveLearningStateFunction trainActiveLearningFunction;
//    iftActiveLearningStateFunction posProcessingActiveLearningFunction;
//}iftGenericActiveLearning;
//
//iftGenericActiveLearning* iftCreateGenericActiveLearner();
//void iftDestroyGenericActiveLearning(iftGenericActiveLearning** pActiveLearning);
//void iftActiveLearningNextState(iftGenericActiveLearning*  genericActiveLearning);

///////////COMMOM active learning
typedef enum ift_ActiveLearningState {
    ACTIVELEARNING_BEGIN,
    ACTIVELEARNING_PREPROCESSING_STATE,
    ACTIVELEARNING_SELECT_SAMPLES_STATE,
    ACTIVELEARNING_SUPERVISING_STATE,
    ACTIVELEARNING_CLASSIFIER_TRAIN_STATE,
    ACTIVELEARNING_POSPROCESSING_STATE,
    ACTIVELEARNING_ENDING,
} iftActiveLearningState;

typedef struct _iftCommonActiveLearning iftCommonActiveLearning;
typedef void (*iftActiveLearningStateFunction) (iftCommonActiveLearning* activeLearningObject);
typedef void (*iftActiveLearningSynchronizeFunction) (iftCommonActiveLearning* activeLearning,
                                                      iftGenericSelector* selector, iftGenericClassifier* classifier);

struct _iftCommonActiveLearning {
    iftGenericSelector* selector;
    iftGenericClassifier* classifier;
    size_t currentLearningCycle;
    iftActiveLearningState nextState;
    bool activeLearningFinished;


    //trasition state functions. this is useful to control the active learning flow, and to check variables
    iftActiveLearningStateFunction beginActiveLearningFunction;
    iftActiveLearningStateFunction preProcessingActiveLearningFunction;
    iftActiveLearningStateFunction selecSamplesActiveLearningFunction;
    iftActiveLearningStateFunction supervisingActiveLearningFunction;
    iftActiveLearningStateFunction trainActiveLearningFunction;
    iftActiveLearningStateFunction posProcessingActiveLearningFunction;
    iftActiveLearningStateFunction endActiveLearningFunction;

    //implementation for each state.
    iftActiveLearningSynchronizeFunction functionUsed2PreProcessing;
    iftActiveLearningSynchronizeFunction functionUsed2SelectSamples;
    iftActiveLearningSynchronizeFunction functionUsed2SuperviseSamples;
    iftActiveLearningSynchronizeFunction functionUsed2TrainClassifier;
    iftActiveLearningSynchronizeFunction functionUsed2PosProcessing;


    const char* path2Samples;
};

iftCommonActiveLearning* iftCreateCommonActiveLearning();
void iftDestroyCommonActiveLearning(iftCommonActiveLearning** activeLearningObject);

void iftInitializeCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning);
void iftPreprocessCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning);
void iftSelecSamplesCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning);
void iftSuperviseSamplesgCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning);
void iftTrainClassifierCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning);
void iftPosProcessCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning);

void iftSelecSamplesSelectorRDS(iftCommonActiveLearning* activeLearning,
                                iftGenericSelector* selectorRDS,
                                iftGenericClassifier* classifier);

void iftSuperviseSamplesSelectorRDS(iftCommonActiveLearning* activeLearning,
                                    iftGenericSelector* selectorRDS,
                                    iftGenericClassifier* classifier);





void iftCommonActiveLearningGoToNextState(iftCommonActiveLearning* activeLearning);
///////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////MACHINE LEARNING PRE_PROCESSING////////////

////PCA
typedef struct _iftPCA {
    double numberComponents;
    iftMatrix* rotationMatrix;
    iftMatrix* rotationMatrix_original;
    iftMatrix* inverseRotationMatrix_original;
    iftMatrix* variance;
    iftMatrix* variance_frac;
    iftMatrix* mean;
    iftMatrix* inputData;
    iftMatrix* outputData;
    iftDataSet* inputDataset;
    iftDataSet* outputDataset;
    bool destroyInputData;
    bool destroyOutData;
}iftPrincipalComponentAnalysis;

iftPrincipalComponentAnalysis* iftCreatePrincipalComponentAnalysis(iftMatrix* matrix);
iftPrincipalComponentAnalysis* iftCreatePrincipalComponentAnalysisGivenDataset(iftDataSet* dataSet);
void iftDestroyPrincipalComponentAnalysis(iftPrincipalComponentAnalysis** pca_pointer);
void iftComputePrincipalComponentAnalysis(iftPrincipalComponentAnalysis *pca);
iftDataSet* iftComputePrincipalComponentAnalysisGivenDataset(iftDataSet *dataSet, double nComponents);
void iftApplyPricipalComponentAnalysis(iftPrincipalComponentAnalysis* pca,iftMatrix* matrix);
void iftApplyPricipalComponentAnalysisGivenDataset(iftPrincipalComponentAnalysis* pca,iftDataSet* dataSet);
void iftComputeKnnGraphStatistics(iftDataSet* dataSet, iftDoubleMatrix** clusterMeanp,iftDoubleMatrix** clusterStandardDeviationp);
//////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

#endif //IFT_IFTACTIVELEARNING_H
