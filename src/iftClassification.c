#include "iftClassification.h"

#include "ift/core/io/Dir.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"



/********************** PRIVATE FUNCTIONS *************************/


iftGenericClassifier* iftCreateGenericClassifier(iftClassifierType classifierType){

    iftGenericClassifier* genericClassifier = (iftGenericClassifier*)iftAlloc(1,sizeof(iftGenericClassifier));
    genericClassifier->classifier = NULL;
    genericClassifier->trainFunction = NULL;
    genericClassifier->predictFunction = NULL;
    genericClassifier->freeFunctionClassifier = NULL;

    switch (classifierType){
        case CLASSIFIER_CUSTOM:
            return genericClassifier;
        case CLASSIFIER_SVM:
            genericClassifier->classifierName = "SVM";
            genericClassifier->classifier = (void*)iftCreateSVMClassifier();
            genericClassifier->trainFunction = iftTrainSvmVoidPointer;
            genericClassifier->predictFunction = iftPredictSvmVoidPointer;
            genericClassifier->freeFunctionClassifier = iftDestroyClassifierVoidPointerSvm;
            break;
        case CLASSIFIER_OPF_SUP:
            genericClassifier->classifierName = "OPF-SUP";
            genericClassifier->classifier = (void*)iftCreateOPFSUPClassifier();
            genericClassifier->trainFunction = iftTrainOpfsupVoidPointer;
            genericClassifier->predictFunction = iftPredictOpfsupVoidPointer;
            genericClassifier->freeFunctionClassifier = iftDestroyClassifierVoidPointerOpfsup;
            break;
        case CLASSIFIER_KNN:
            genericClassifier->classifier = "KNN";
            genericClassifier->classifier = (void*)iftCreateKnnClassifier();
            genericClassifier->trainFunction = iftTrainKnnClassifierVoidPointer;
            genericClassifier->predictFunction = iftPredictKnnClassifierVoidPointer;
            genericClassifier->freeFunctionClassifier = iftDestroyClassifierVoidPointerKnn;
            break;
        case CLASSIFIER_RANDOM:
            genericClassifier->classifierName = "Random";
            genericClassifier->classifier = (void*)iftCreateRandomClassifier();
            genericClassifier->trainFunction = iftTrainRandomClassifierGeneric;
            genericClassifier->predictFunction = iftPredictRandomClassifierGeneric;
            genericClassifier->freeFunctionClassifier = iftDestroyClassifierVoidPointerRandom;
            break;
        default:
            return genericClassifier;
    }

    return genericClassifier;
}


void iftTrainGenericClassifier(iftGenericClassifier *genericClassifier, iftDataSet *dataSet){
    if(genericClassifier == NULL){
        iftError("Generic classifier is null\n","iftFit_genericClassifier");
        return;
    }
    if(genericClassifier->classifier == NULL){
        iftError("Inner classifier is empty \n","iftFit_genericClassifier");
        return;
    }
    if(genericClassifier->trainFunction == NULL){
        iftError("Fit function not defined \n","iftFit_genericClassifier");
        return;
    }

    genericClassifier->trainFunction(genericClassifier->classifier,dataSet);
}

void iftPredictGenericClassifier(iftGenericClassifier *genericClassifier, iftDataSet *dataSet){
    if(genericClassifier == NULL){
        iftError("Generic classifier is null\n","iftFit_genericClassifier");
        return;
    }
    if(genericClassifier->classifier == NULL){
        iftError("Inner classifier is empty \n","iftFit_genericClassifier");
        return;
    }
    if(genericClassifier->predictFunction == NULL){
        iftError("Predict function not defined \n","iftFit_genericClassifier");
        return;
    }

    genericClassifier->predictFunction(genericClassifier->classifier,dataSet);
}

void iftDestroyGenericClassifier(iftGenericClassifier** pGenericClassifier ){
    if(pGenericClassifier == NULL || *pGenericClassifier == NULL){
        return;
    }
    iftGenericClassifier* genericClassifier = *pGenericClassifier;
    if(genericClassifier->freeFunctionClassifier){
        genericClassifier->freeFunctionClassifier(genericClassifier->classifier);
    }
    iftFree(genericClassifier);
    *pGenericClassifier = NULL;
}


/////////SVM
iftSVM* iftCreateSVMClassifier(){
    iftSVM* svmClassifier =  iftAlloc(1,sizeof(iftSVM));
    svmClassifier->multiClass = IFT_OVO;

    svmParameter *param = (svmParameter *) iftAlloc(1, sizeof(svmParameter));

    // Default libsvm parameter values
    param->svm_type = C_SVC;
    param->kernel_type = LINEAR;
    param->degree = 3;
    param->gamma = 0;  // 1/num_features
    param->coef0 = 0;
    param->nu = 0.5;
    param->cache_size = 100;
    param->C = 1;
    param->eps = 1e-5;
    param->p = 0.1;
    param->shrinking = 1;
    param->probability = 0;
    param->nr_weight = 0;
    param->weight_label = NULL;
    param->weight = NULL;
    svmClassifier->params = param;

    svmClassifier->params->kernel_type = LINEAR;
    svmClassifier->params->C = 1;
    svmClassifier->params->gamma = 1;
    return svmClassifier;
}
void iftDestroyClassifierVoidPointerSvm(void* psvmClassifier){
    iftSVM* aux = ((iftSVM*)psvmClassifier);
    iftDestroySVMClassifier(&aux);
}
void iftTrainSvmVoidPointer(void *svmClassifier, iftDataSet *dataSet){
    iftSVM* svm = (iftSVM*)svmClassifier;
    iftSVMTrain(svm,dataSet);
    //fit_svm(svm, dataSet,sample_status);
}

void iftPredictSvmVoidPointer(void *svmClassifier, iftDataSet *dataSet){
    iftSVM* svm = (iftSVM*)svmClassifier;
    iftSVMClassify(svm,dataSet,IFT_TEST);
}
////////////

//////////OPF-SUP
iftOPFSUPClassifier* iftCreateOPFSUPClassifier(){
    iftOPFSUPClassifier* opfsupClassifier = (iftOPFSUPClassifier*)iftAlloc(1,sizeof(iftOPFSUPClassifier));
    opfsupClassifier->completeGraph = NULL;
    opfsupClassifier->distanceTable = NULL;
    opfsupClassifier->recomputeDistanceTable = true;
    return opfsupClassifier;
}

void iftDestroyClassifierOpfsup(iftOPFSUPClassifier** pOpfsupClassifier){
    if(pOpfsupClassifier == NULL || *pOpfsupClassifier == NULL){
        return;
    }
    iftOPFSUPClassifier* aux = *pOpfsupClassifier;
    if(aux->completeGraph != NULL){
        iftDestroyCplGraph(&(aux->completeGraph));
        (*pOpfsupClassifier)->completeGraph = NULL;
    }
    if(aux->distanceTable != NULL){
        iftDestroyDistanceTable(&(aux->distanceTable));
        (*pOpfsupClassifier)->distanceTable = NULL;
    }
    iftFree((*pOpfsupClassifier));
    *pOpfsupClassifier = NULL;
}

void iftDestroyClassifierVoidPointerOpfsup(void* pOpfsupClassifier){
    iftOPFSUPClassifier* opfsupClassifier = ((iftOPFSUPClassifier*)pOpfsupClassifier);
    iftDestroyClassifierOpfsup(&opfsupClassifier);

}

void iftTrainOpfsupVoidPointer(void *opfsupClassifier, iftDataSet *dataSet){
    iftOPFSUPClassifier* opfsup = (iftOPFSUPClassifier*)opfsupClassifier;
    iftTrainOpfsup(opfsup,dataSet);
}

void iftTrainOpfsup(iftOPFSUPClassifier* opfsup, iftDataSet* dataSet){
    if(dataSet == NULL){
        iftError("Dataset is empty","iftTrainOpfsup");
    }
    //iftOPFSUPClassifier* opfsup = (iftOPFSUPClassifier*)opfsupClassifier;

    //////creating graph
    iftCplGraph *graph=(iftCplGraph *)iftAlloc(1,sizeof(iftCplGraph));
    int s, u, nnodes = iftCountSamples(dataSet,IFT_TRAIN);
    dataSet->ntrainsamples = nnodes;
    if (nnodes == 0){
        iftWarning("No samples for training", "iftCreateCplGraph");
        return;
    }

    graph->nnodes = nnodes;
    graph->node   = (iftCplNode *)iftAlloc(nnodes,sizeof(iftCplNode));

    if (graph->node == NULL){
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateCplGraph");
    }

    graph->pathval       = iftAllocFloatArray(nnodes);
    graph->ordered_nodes = iftAllocIntArray(nnodes);
    graph->Q             = iftCreateFHeap(nnodes,graph->pathval);
    graph->Z             = dataSet;
    u=0;
    for (s=0; s < dataSet->nsamples; s++) {
        if (iftHasSampleStatus(dataSet->sample[s], IFT_TRAIN)){
            graph->node[u].pred     = IFT_NIL;
            graph->node[u].sample   = s;
            u++;

        }
    }
    if(opfsup->completeGraph != NULL){
        iftDestroyCplGraph(&(opfsup->completeGraph));
    }
    opfsup->completeGraph = graph;

    /////////////////////////////////

    /////////////computing distance table

    if(opfsup->distanceTable != NULL){
        if(opfsup->recomputeDistanceTable == true){
            iftDestroyDistanceTable(&(opfsup->distanceTable));
            if(dataSet->iftArcWeight){
                opfsup->distanceTable = iftCompEuclDistanceTable(dataSet, dataSet);
            }else{
                opfsup->distanceTable = iftCompDistanceTable(dataSet, dataSet);
            }
        }
    }else{
        if(dataSet->iftArcWeight){
            opfsup->distanceTable = iftCompEuclDistanceTable(dataSet, dataSet);
        }else{
            opfsup->distanceTable = iftCompDistanceTable(dataSet, dataSet);
        }
    }


    //////////////////////////

    /////computing prototypes by mst
    iftSet  *S = NULL;
    int  t, v;
    iftDataSet *Z=graph->Z;
    float dist;

    // initialization
    for (u = 0; u < graph->nnodes; u++) {
        graph->pathval[u]  = IFT_INFINITY_FLT;
        s = graph->node[u].sample;
        Z->sample[s].label = 0;
        graph->node[u].pred  = IFT_NIL;
    }

    graph->pathval[0]   = 0;
    iftInsertFHeap(graph->Q, 0); // IFT_GRAY to 0

    // Minimum Spanning Tree by Degenerated IFT
    while ( !iftEmptyFHeap(graph->Q) ) {

        u = iftRemoveFHeap(graph->Q); // BLACK to u
        s = graph->node[u].sample;
        v = graph->node[u].pred;

        if (v != IFT_NIL){
            t = graph->node[v].sample;
            if (Z->sample[s].truelabel != Z->sample[t].truelabel){
                if (Z->sample[s].label == 0){
                    Z->sample[s].label = 1;
                    iftInsertSet(&S,u);
                }
                if (Z->sample[t].label == 0){
                    Z->sample[t].label = 1;
                    iftInsertSet(&S,v);
                }
            }
        }

        for (v=0; v < graph->nnodes; v++){
            if (graph->Q->color[v] != IFT_BLACK){
                t = graph->node[v].sample;
                dist = opfsup->distanceTable->distance_table[s][t];

                if ( dist < graph->pathval[v] ) {
                    graph->node[v].pred = u;
                    graph->pathval[v]   = dist;
                    if(graph->Q->color[v] == IFT_WHITE)
                        iftInsertFHeap(graph->Q, v); // IFT_GRAY TO v
                    else
                        iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
                }
            }
        }
    }

    if(S == NULL){
        Z->sample[0].label = 1;
        iftInsertSet(&S,0);
    }

    iftResetFHeap(graph->Q);
    ///////

//    if(opfsup->protoTypes != NULL){
//        iftDestroySet(&(opfsup->protoTypes));
//    }
    opfsup->protoTypes = S;

    //Train
    int i;
    float tmp;
    // initialization

    for (u = 0; u < graph->nnodes; u++) {
        graph->pathval[u] = IFT_INFINITY_FLT;
        s = graph->node[u].sample;
        Z->sample[s].label  = 0;
        graph->node[u].pred = IFT_NIL;
    }
    while (S != NULL) {
        u = iftRemoveSet(&S);
        s = graph->node[u].sample;
        graph->pathval[u]     = 0;
        Z->sample[s].label    = Z->sample[s].truelabel;
        iftInsertFHeap(graph->Q, u); // IFT_GRAY to u
    }

    // IFT with fmax
    i=0;
    while ( !iftEmptyFHeap(graph->Q) ) {
        u = iftRemoveFHeap(graph->Q);
        graph->ordered_nodes[i]=u; i++;
        s = graph->node[u].sample;
        Z->sample[s].weight = graph->pathval[u];

        for (v=0; v < graph->nnodes; v++){
            if (graph->Q->color[v] != IFT_BLACK){
                t = graph->node[v].sample;
                if (graph->pathval[u] < graph->pathval[v]){
                    dist = opfsup->distanceTable->distance_table[s][t];
                    tmp  = iftMax(graph->pathval[u], dist);
                    if ( tmp < graph->pathval[v] ) {
                        graph->node[v].pred  = u;
                        Z->sample[t].label   = Z->sample[s].label;
                        graph->pathval[v]    = tmp;
                        if(graph->Q->color[v] == IFT_WHITE)
                            iftInsertFHeap(graph->Q, v); // IFT_GRAY to u
                        else
                            iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
                    }
                }
            }
        }
    }


    iftResetFHeap(graph->Q);
    Z->ngroups = Z->nclasses;
    ////////////////
}

void iftPredictOpfsup(iftOPFSUPClassifier* opfsup, iftDataSet* dataSet){
    //iftOPFSUPClassifier* opfsup = (iftOPFSUPClassifier*)opfsupClassifier;
    iftDataSet *Ztrain = opfsup->completeGraph->Z;

    // checkers
    if (Ztrain->nfeats != dataSet->nfeats) {
        iftError("Train. and Test. Datasets with different Number of Feats: %d != %d",
                 "iftClassify", Ztrain->nfeats, dataSet->nfeats);
    }


#pragma omp parallel for schedule(auto)
    for (int t = 0; t < dataSet->nsamples; t++) {
        if(iftHasSampleStatus(dataSet->sample[t], IFT_TEST)) {
            iftFindMinimumCostAndClassifySample(opfsup->completeGraph, dataSet, t);
        }
    }
}

void iftPredictOpfsupVoidPointer(void *opfsupClassifier, iftDataSet *dataSet){
    iftOPFSUPClassifier* opfsup = (iftOPFSUPClassifier*)opfsupClassifier;
    iftPredictOpfsup(opfsup,dataSet);
}

////////////////////////


/////////////RANDOM classifier
iftRandomClassifier* iftCreateRandomClassifier(){
    iftRandomClassifier* randomClassifier = iftAlloc(1,sizeof(iftRandomClassifier));
    randomClassifier->classProbabilities = NULL;
    randomClassifier->numberOfClasses = 0;
    //randomClassifier->seed = 0;
    return randomClassifier;
}

void iftDestroyRandomClassifier(iftRandomClassifier** pRandomClassifier){
    if(pRandomClassifier == NULL || *pRandomClassifier == NULL){
        return;
    }
    iftRandomClassifier* aux = *pRandomClassifier;
    if(aux->classProbabilities != NULL){
        free( (*pRandomClassifier)->classProbabilities );
    }
    iftFree((*pRandomClassifier));
    (*pRandomClassifier) = NULL;
}

void iftDestroyClassifierVoidPointerRandom(void* pRandomClassifier){
    iftRandomClassifier* randomClassifier = ((iftRandomClassifier*)pRandomClassifier);
    iftDestroyRandomClassifier(&randomClassifier);
}

void iftTrainRandomClassifierGeneric(void* randomClassifier, iftDataSet* dataSet){
    iftRandomClassifier* randClassifier = (iftRandomClassifier*)randomClassifier;
    iftTrainRandomClassifier(randClassifier,dataSet);
}

void iftTrainRandomClassifier(iftRandomClassifier* randomClassifier, iftDataSet* dataSet){
    //iftRandomClassifier* randClassifier = (iftRandomClassifier*)randomClassifier;
    if(dataSet->nclasses == 0){
        iftCountNumberOfClassesDataSet(dataSet);
    }
    randomClassifier->numberOfClasses = dataSet->nclasses;
    randomClassifier->classProbabilities = iftAlloc(randomClassifier->numberOfClasses,sizeof(double));
    for (int sampleIndex = 0; sampleIndex < dataSet->nsamples; ++sampleIndex) {

        if(dataSet->sample[sampleIndex].truelabel == 0){
            continue;
        }
        if(iftHasSampleStatus(dataSet->sample[sampleIndex], IFT_TRAIN)){
            randomClassifier->classProbabilities[dataSet->sample[sampleIndex].truelabel - 1]++;
        }
    }
    double total = 0;
    for (int classIndex = 0; classIndex < randomClassifier->numberOfClasses; ++classIndex){
        total += randomClassifier->classProbabilities[classIndex];
    }

    for (int classIndex = 0; classIndex < randomClassifier->numberOfClasses; ++classIndex){
        randomClassifier->classProbabilities[classIndex] /= (total + 0.00000001);
    }
}

void iftPredictRandomClassifierGeneric(void* randomClassifier, iftDataSet* dataSet){
    iftRandomClassifier* randClassifier = (iftRandomClassifier*)randomClassifier;
    iftPredictRandomClassifier(randClassifier,dataSet);
}

void iftPredictRandomClassifier(iftRandomClassifier* randomClassifier, iftDataSet* dataSet){
    //iftRandomClassifier* randClassifier = (iftRandomClassifier*)randomClassifier;
    int labelIndex;
    for (int i = 0; i < dataSet->nsamples; ++i) {
        if( !iftHasSampleStatus(dataSet->sample[i], IFT_TEST) ){
            continue;
        }
        double prob = iftRandomNormalDouble();
        double probAcc = 0;
        for (labelIndex = 0; labelIndex < randomClassifier->numberOfClasses; ++labelIndex) {
            probAcc += randomClassifier->classProbabilities[labelIndex];
            if(probAcc > prob){
                break;
            }
        }
        dataSet->sample[i].label = labelIndex+1;
    }
}

////////////////////////////////////

///////////Knn classifier
iftKnnClassifier* iftCreateKnnClassifier(){
    iftKnnClassifier* knnClassifier = iftAlloc(1,sizeof(iftKnnClassifier));
    knnClassifier->datasetTrain = NULL;
    knnClassifier->distanceFunction = iftComputeL2Norm;
    knnClassifier->distanceFunctionParameters = NULL;
    knnClassifier->kNearest = 1;
    knnClassifier->copyfeatsOfTraiSamples = false;
    return knnClassifier;
}

void iftDestroyKnnClassifier(iftKnnClassifier** pKnnClassifier){
    if(pKnnClassifier == NULL || *pKnnClassifier == NULL){
        return;
    }
    iftKnnClassifier* aux = *pKnnClassifier;
    if(aux->datasetTrain != NULL){
        iftDestroyDataSet( &((*pKnnClassifier)->datasetTrain) );
    }
    if(aux->distanceFunctionParameters != NULL){
        iftDestroyArgumentList( &((*pKnnClassifier)->distanceFunctionParameters) );
    }
    free((*pKnnClassifier));
    (*pKnnClassifier) = NULL;
}

void iftDestroyClassifierVoidPointerKnn(void* pKnnClassifier){
    iftKnnClassifier* knnClassifier = ((iftKnnClassifier*)pKnnClassifier);
    iftDestroyKnnClassifier(&knnClassifier);
}

void iftTrainKnnClassifierVoidPointer(void* knnClassifier, iftDataSet* dataSet){
    iftKnnClassifier* knn = (iftKnnClassifier*)knnClassifier;
    iftTrainKnnClassifier(knn,dataSet);
}

void iftTrainKnnClassifier(iftKnnClassifier* knnClassifier, iftDataSet* dataSet){
    int numberSamples = iftCountSamples(dataSet, IFT_TRAIN);
    if(numberSamples <= 0){
        iftWarning("There are not train samples","iftTrainKnnClassifier");
        return;
    }
    iftDataSet* datasetTrain = iftCreateDataSetAux(numberSamples,dataSet->nfeats,knnClassifier->copyfeatsOfTraiSamples);
    int currentSampleTrainIndex = 0;
    for (int sampleIndex = 0; sampleIndex < dataSet->nsamples; ++sampleIndex) {
        if(iftHasSampleStatus(dataSet->sample[sampleIndex], IFT_TRAIN)){
            iftCopySample( &(dataSet->sample[sampleIndex]),
                           &(datasetTrain->sample[currentSampleTrainIndex]),
                           dataSet->nfeats,
                           knnClassifier->copyfeatsOfTraiSamples);
            currentSampleTrainIndex++;
        }
    }
    datasetTrain->nclasses = dataSet->nclasses;
    if(datasetTrain->nclasses <= 0){
        iftCountNumberOfClassesDataSet(datasetTrain);
    }
    knnClassifier->datasetTrain = datasetTrain;
}

void iftPredictKnnClassifierVoidPointer(void* knnClassifier, iftDataSet* dataSet){
    iftKnnClassifier* knn = (iftKnnClassifier*)knnClassifier;
    iftPredictKnnClassifier(knn,dataSet);
}

void iftPredictKnnClassifier(iftKnnClassifier* knnClassifier, iftDataSet* dataSetTest){
    if(knnClassifier->datasetTrain == NULL){
        iftError("Classifier is not trained","iftPredictKnnClassifier");
        return;
    }
    if(dataSetTest->nfeats != knnClassifier->datasetTrain->nfeats){
        iftError("Feature vector mismatch","iftPredictKnnClassifier");
        return;
    }
    iftArgumentList* argumentList = knnClassifier->distanceFunctionParameters;
    iftDataSet* trainDataset = knnClassifier->datasetTrain;
    int featureVectorSize = trainDataset->nfeats;
    int *indices = (int *)calloc(knnClassifier->datasetTrain->nsamples,sizeof(int));
    double* distances = (double *)calloc(knnClassifier->datasetTrain->nsamples,sizeof(size_t));
    size_t *count = (size_t *)calloc(knnClassifier->datasetTrain->nclasses+1,sizeof(size_t));
    int sampleIndex;
    int trueLabel;
    int maxCount;
    int maxLabelIndex = 0;
    for (size_t testSampleIndex = 0; testSampleIndex < dataSetTest->nsamples; ++testSampleIndex) {
        if(iftHasSampleStatus(dataSetTest->sample[testSampleIndex], IFT_TEST)){
            for (size_t i = 0; i <= knnClassifier->datasetTrain->nclasses; ++i) {
                count[i] = 0;
            }
            for (int trainSampleindex = 0; trainSampleindex < knnClassifier->datasetTrain->nsamples; ++trainSampleindex) {
                indices[trainSampleindex] = trainSampleindex;
                distances[trainSampleindex] = knnClassifier->distanceFunction(trainDataset->sample[trainSampleindex].feat,
                                                                              dataSetTest->sample[testSampleIndex].feat,
                                                                              featureVectorSize,
                                                                              argumentList);
            }
            iftDQuickSort(distances, indices, 0, knnClassifier->datasetTrain->nsamples-1, IFT_INCREASING);
            for (size_t k = 0; k < knnClassifier->kNearest; ++k) {
                sampleIndex = indices[k];
                trueLabel = knnClassifier->datasetTrain->sample[sampleIndex].truelabel;
                count[trueLabel]++;
            }
            maxLabelIndex = 0;
            maxCount = 0;
            for (size_t i = 0; i <= knnClassifier->datasetTrain->nclasses; ++i) {
                if(count[i] > maxCount){
                    maxCount = count[i];
                    maxLabelIndex = i;
                }
            }
            dataSetTest->sample[testSampleIndex].label = maxLabelIndex;
        }
    }
    free(count);
    free(indices);
    free(distances);
}

/* ///////////////////////////////////////////////// */

/**
 * @brief Supervised OPF Classification
 *
 * For a given test sample from the dataset Ztest, it returns the
 * minimum cost and assigns a label to the sample.
 * 
 * @author Alexandre Falcao.
 * @param graph           Supervised OPF.
 * @param Ztest           Dataset with the testing sample to be classified.
 * @param t               Test sample in Ztest.
 * @return                the minimum cost to the test sample.
 */

float iftFindMinimumCostAndClassifySample(const iftCplGraph *graph, iftDataSet *Ztest, int t)
{

    int length = 0;
    int i = 0; // index to get the ith node, ordered by the path cost, in the ordered nodes array
    int u = graph->ordered_nodes[i]; // gets the index of the node
    int s = graph->node[u].sample; // gets the index of the train. sample related to the node u
    iftDataSet *Ztrain=graph->Z; // gets training set

    // gets the distance between the test sample and a training sample that follows the order of cost
    float weight = 0.0;
    if (iftDist == NULL)
        weight = Ztrain->iftArcWeight(Ztrain->sample[s].feat, Ztest->sample[t].feat, Ztrain->alpha, Ztrain->nfeats);
    else
        weight = iftDist->distance_table[s][t];

    // minimize the cost map
    float min_cost    = iftMax(graph->pathval[u], weight);
    int   label       = Ztrain->sample[s].label;

    // if the current training sample is not the last one and the
    // first minimum cost has not been found yet

    while ((i < graph->nnodes - 1) && (min_cost > graph->pathval[graph->ordered_nodes[i + 1]])) {
        u = graph->ordered_nodes[i + 1];
        s = graph->node[u].sample;

        if (iftDist == NULL)
            weight = Ztrain->iftArcWeight(Ztrain->sample[s].feat, Ztest->sample[t].feat, Ztrain->alpha, Ztrain->nfeats);
        else
            weight = iftDist->distance_table[s][t];

        // if the current training sample offers a lower cost than the current min cost, then "conquer" the test sample
        float tmp = iftMax(graph->pathval[u], weight);

        if (tmp < min_cost) {
            length = graph->node[u].length + 1;
            min_cost = tmp;
            label    = Ztrain->sample[s].label;
        }
        i++;
    }

    // classify
    Ztest->sample[t].label  = label;
    Ztest->sample[t].weight = length;

    return(min_cost);
}



/**
 * @brief Supervised OPF Classification with Certainty Values
 *
 * For a given test sample from the dataset Ztest, it returns the
 * minimum cost with respect to a second closest class (i.e., a class
 * different from the class that offered the minimum cost to that test
 * sample).
 * 
 * @author Alexandre Falcao 
 * @param graph           Supervised OPF.
 * @param Ztest           Dataset with the testing sample to be classified.
 * @param t               Test sample in Ztest.
 * @return                the minimum cost from the second closest class to the test sample.
 */

float iftFindMinimumCostFromAnotherClass(const iftCplGraph *graph, iftDataSet *Ztest, int t) {

    int i, u, s;
    iftDataSet *Ztrain=graph->Z;


    // find the first training sample from a class different than the
    // one used to classify t. We will assume that this sample can
    // always be found, since the data set cannot have a single class.
    for (i=0, u=graph->ordered_nodes[i], s=graph->node[u].sample;       \
       (i < graph->nnodes)&&(Ztrain->sample[s].label==Ztest->sample[t].label); i++) {
        u=graph->ordered_nodes[i];
        s=graph->node[u].sample;
    }

    // gets the distance between the test sample and that training sample.
    float weight = 0.0;
    if (iftDist == NULL)
        weight = Ztrain->iftArcWeight(Ztrain->sample[s].feat, Ztest->sample[t].feat, Ztrain->alpha, Ztrain->nfeats);
    else
        weight = iftDist->distance_table[s][t];

    // find the second minimum cost with respect to a class different
    // than the one used to classify t.

    float min_cost    = iftMax(graph->pathval[u], weight);

    // if the current training sample is not the last one and the
    // second minimum cost has not been found yet


    for (; i < graph->nnodes - 1; i++)
    {
        u = graph->ordered_nodes[i + 1];
        s = graph->node[u].sample;

        if (Ztrain->sample[s].label!=Ztest->sample[t].label){

            if (min_cost <= graph->pathval[u]) // we have found the second minimum cost
                return(min_cost);

            if (iftDist == NULL)
                weight = Ztrain->iftArcWeight(Ztrain->sample[s].feat, Ztest->sample[t].feat, Ztrain->alpha, Ztrain->nfeats);
            else
                weight = iftDist->distance_table[s][t];

            // if the current training sample offers a lower cost than the current min cost, then replace it
            float tmp = iftMax(graph->pathval[u], weight);
            if (tmp < min_cost) {
                min_cost = tmp;
            }
        }
    }

    return(min_cost);
}


/**
 * @brief      Builds a DataSet of brightness (1 feat per sample) from a test_img.
 *
 * @param[in]  test_img       The test image used to build the dataset.
 * @param[in]  smooth_factor  A smooth factor. Use 0.0 to ignore it.
 *
 * @return     Testing dataset of brightness (1 feat per sample).
 */
iftDataSet *iftBuildTestImageDataSet(const iftImage *test_img, float smooth_factor) {
    if (test_img == NULL)
        iftError("Testing Image is NULL", "iftClassifyImage");

    // gets a smooth image, if the smooth factor is > 0.0
    iftImage *smooth_img = NULL;
    if (smooth_factor > 0.0) {
        /* Smoothing image */
        iftKernel *K = NULL;
        if (iftIs3DImage(test_img))
            K = iftGaussianKernel(1.5, smooth_factor * iftMaximumValue(test_img));
        else
            K = iftGaussianKernel2D(1.5, smooth_factor * iftMaximumValue(test_img));

        smooth_img = iftLinearFilter(test_img, K);
        test_img   = smooth_img;

        iftDestroyKernel(&K);
    }


    iftMImage *mimg = NULL;
    if (!iftIsColorImage(test_img)) {
        mimg = iftImageToMImage(test_img, GRAY_CSPACE);
    } else {
        mimg = iftImageToMImage(test_img, LAB_CSPACE);
    }


    // gets the testing dataset
    iftDataSet *Z = iftMImageToDataSet(mimg, NULL,0);


    // DESTROYERS
    if (smooth_img != NULL) {
        iftDestroyImage(&smooth_img);
    }

    return Z;
}



/**
 * @brief Return the adjacency radius for each training sample of a Knn Graph
 * as its distance to its i-th nearest neighbor defined by the rank.
 * 
 * E.g: rank == 1 considers the distance to the nearest neighbor as adjacency radius.
 * rank == graph->k considers the distance to the k-th nearest neighbor (the farthest one).
 * 
 * @param graph Knn Graph.
 * @param rank Rank of the considered nearest neighbor that is used to define the adjacency radius. 
 * @return iftDblArray* (graph->nnodes,) Array with the adjacency radius of each training sample in graph.
 * 
 * @author Samuel Martins
 * @date Dec 3, 2019
 */
iftDblArray *iftAdjacencyRadiusKnnGraphByRank(const iftKnnGraph *graph, int rank) {
    iftDblArray *adj_radius = iftCreateDblArray(graph->nnodes);

    #pragma omp parallel for
    for (int u = 0; u < graph->nnodes; u++) {
        const iftAdjSet *adj = graph->node[u].adj;

        for (int k = 1; k < rank; k++) {
            adj = adj->next;
        }
        adj_radius->val[u] = adj->arcw;        
    }

    return adj_radius;
}
/******************************************************************/


/********************** PUBLIC FUNCTIONS *************************/
iftLogReg* iftCreateLogReg() {
    return (iftLogReg*) iftAlloc(1, sizeof(iftLogReg));
}

iftCplGraph *iftCreateCplGraph(iftDataSet *Z)
{
    iftCplGraph *graph=(iftCplGraph *)iftAlloc(1,sizeof(iftCplGraph));
    int s, u, nnodes = Z->ntrainsamples;

    if (nnodes == 0){
        iftError("No samples for training", "iftCreateCplGraph");
    }

    graph->nnodes = nnodes;
    graph->node   = (iftCplNode *)iftAlloc(nnodes,sizeof(iftCplNode));

    if (graph->node == NULL){
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateCplGraph");
    }

    graph->pathval       = iftAllocFloatArray(nnodes);
    graph->ordered_nodes = iftAllocIntArray(nnodes);
    graph->Q             = iftCreateFHeap(nnodes,graph->pathval);
    graph->Z             = Z;
    u=0;
    for (s=0; s < Z->nsamples; s++) {
      if (iftHasSampleStatus(Z->sample[s], IFT_TRAIN)){
            graph->node[u].pred     = IFT_NIL;
            graph->node[u].sample   = s;
            u++;
        }
    }

    return(graph);
}

iftLogReg* iftLogRegTrain(iftDataSet* Z) {
    iftMatrix* X = iftDataSetToFeatureMatrixHomogCoord(Z, IFT_TRAIN);
    iftMatrix* Y = iftDataSetToLabelsMatrix(Z, IFT_TRAIN);

    iftLogReg* logreg = iftCreateLogReg();
    logreg->nclasses = Z->nclasses;
    logreg->nfeats = Z->nfeats;

    logreg->coef = iftLeastSquares(X, Y);

    iftDestroyMatrix(&X);
    iftDestroyMatrix(&Y);

    return logreg;
}

void iftLogRegClassify(iftLogReg* logreg, iftDataSet* Z) {

    if(logreg->nfeats!=Z->nfeats) {
        iftError("The Classifier was trained with a different number of features.", "iftLogRegClassify");
    }

    iftMatrix* X = iftDataSetToFeatureMatrixHomogCoord(Z, IFT_TEST);
    iftMatrix* Y = iftMultMatrices(X, logreg->coef);

    int idx = 0;
    for (int i = 0; i < Z->nsamples; ++i) {
        if(iftHasSampleStatus(Z->sample[i], IFT_TEST)){
            Z->sample[i].weight = -IFT_INFINITY_FLT;
            for (int c = 0; c < logreg->nclasses; ++c) {
                float prob = iftSigm(Y->val[iftGetMatrixIndex(Y, c, idx)]);
                if(prob>Z->sample[i].weight) {
                    Z->sample[i].weight = prob;
                    Z->sample[i].label = c + 1;
                }
            }
            idx++;
        }
    }

    iftDestroyMatrix(&X);
    iftDestroyMatrix(&Y);
}

void iftDestroyLogReg(iftLogReg** logreg) {
    iftDestroyMatrix(&(*logreg)->coef);
    iftFree(*logreg);
    *logreg = NULL;
}


iftLogReg *iftReadLogReg(const char *pathname) {
    char *tmp_dir = iftMakeTempDir("tempdir_", NULL, NULL);

    iftUnzipFile(pathname, tmp_dir);

    char *info_path = iftJoinPathnames(2, tmp_dir, "info.txt");
    if (!iftFileExists(info_path))
        iftError("Info file not found", "iftReadLogReg");

    iftLogReg *logreg = iftCreateLogReg();

    FILE *fp = fopen(info_path, "rb");
    if (fscanf(fp, "%d %d %f\n", &logreg->nclasses, &logreg->nfeats, &logreg->error) != 3)
        iftError("Error when reading info file", "iftReadLogReg");
    fclose(fp);
    iftFree(info_path);

    char *mat_path = iftJoinPathnames(2, tmp_dir, "coef.bin");
    if (!iftFileExists(mat_path))
        iftError("Coef file not found", "iftReadLogReg");
    logreg->coef = iftReadRawMatrix(mat_path);
    iftFree(mat_path);

    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);

    return logreg;
}


void iftWriteLogReg(const iftLogReg *logreg, const char *pathname) {
    if (logreg == NULL)
        iftError("Log. Reg. Classifier is NULL", "iftWriteLogReg");
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftWriteLogReg");
    if (!iftEndsWith(pathname, ".zip"))
        iftError("Invalid Extension for the Log Reg Classifier: \"%s\"... ... Use *.zip",
                 "iftWriteLogReg", pathname);

    char *tmp_dir = iftMakeTempDir("tempdir_", NULL, NULL);

    char *info_path = iftJoinPathnames(2, tmp_dir, "info.txt");
    FILE *fp = fopen(info_path, "wb");
    fprintf(fp, "%d %d %f\n", logreg->nclasses, logreg->nfeats, logreg->error);
    fclose(fp);
    iftFree(info_path);

    char *mat_path = iftJoinPathnames(2, tmp_dir, "coef.bin");
    iftWriteRawMatrix(logreg->coef, mat_path);
    iftFree(mat_path);

    iftZipDirContent(tmp_dir, pathname);

    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
}



iftCplGraph *iftReadCplGraph(const char* pathname) {
    FILE* fp = NULL;
    char *path = NULL;
    char *tmp_dir = NULL;
    iftCplGraph *graph = NULL;
    iftDataSet* Z = NULL;

    if(pathname == NULL)
        iftError("Input Pathname must not be NULL", "iftReadCplGraph");

    if(!iftFileExists(pathname))
        iftError("Input Pathname does not refer to a valid file", "iftReadCplGraph", pathname);

    if(!iftCompareStrings(iftFileExt(pathname), ".zip"))
        iftError("Input Pathname for the complete graph must end with a \".zip\" file extension", "iftReadCplGraph");

    tmp_dir = iftMakeTempDir("tmp_complgraph", NULL, NULL);

    if(!iftDirExists(tmp_dir))
        iftError("Error when creating temporary dir %s for unzipping the complete graph file content",
                 "iftReadCplGraph",
                 tmp_dir);

    iftUnzipFile(pathname, tmp_dir);

    path = iftJoinPathnames(2, tmp_dir, "dataset.zip");

    if (!iftFileExists(path))
        iftError("Cannot open file: \"%s\". File \"%s\" does not refer to a valid complete graph", "iftReadCplGraph",
                 path, pathname);

    /* Loading data set */
    Z = iftReadOPFDataSet(path);
    iftFree(path);

    graph = iftCreateCplGraph(Z);

    /* Reading info */
    path = iftJoinPathnames(2, tmp_dir, "info.data");

    if (!iftFileExists(path))
        iftError("Cannot open file: \"%s\". File \"%s\" does not refer to a valid complete graph", "iftReadCplGraph",
                 path, pathname);

    fp = fopen(path, "rb");

    iftFree(path);

    // Writing an array of iftCplNode
    if(fread(graph->node, sizeof(iftCplNode), graph->nnodes, fp) != graph->nnodes)
        iftError("Error when reading the nodes of the complete graph", "iftReadCplGraph");

    if(fread(graph->pathval, sizeof(float), graph->nnodes, fp) != graph->nnodes)
        iftError("Error when reading the path values of the complete graph", "iftReadCplGraph");

    if(fread(graph->ordered_nodes, sizeof(int), graph->nnodes, fp) != graph->nnodes)
        iftError("Error when reading the order of nodes of the complete graph", "iftReadCplGraph");

    fclose(fp);


    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);

    return graph;
}


void iftWriteCplGraph(const iftCplGraph *graph, const char *pathname) {
    FILE* fp = NULL;
    char *path = NULL;
    char *tmp_dir = NULL;

    if (graph == NULL)
        iftError("Graph is NULL", "iftWriteCplGraph");
    if (pathname == NULL)
        iftError("Out Pathname is NULL", "iftWriteCplGraph");
    if (graph->Z == NULL)
        iftError("Complete graph\'s reference dataset is NULL", "iftWriteCplGraph");

    if(!iftCompareStrings(iftFileExt(pathname), ".zip"))
        iftError("The output Pathname for the complete graph must end with extension \".zip\"", "iftWriteCplGraph");

    /* Creating temporary dir */
    tmp_dir = iftMakeTempDir("tmp_cmplgraph", NULL, NULL);

    /* Opening file that will contain the basic information about the graph */
    path = iftJoinPathnames(2, tmp_dir, "info.data");
    fp = fopen(path, "wb");

    if (fp == NULL) {
        iftError("Cannot open %s file", "iftWriteCplGraph", path);
    }
    iftFree(path);

    /* IMPORTANT: we do not save nnodes since it will be set by iftCreateCplGraph when recreating the graph from the data
       set
     */
    // Writing an array of iftCplNode
    if(fwrite(graph->node, sizeof(iftCplNode), graph->nnodes, fp) != graph->nnodes)
        iftError("Error when writing the nodes of the complete graph", "iftWriteCplGraph");

    if(fwrite(graph->pathval, sizeof(float), graph->nnodes, fp) != graph->nnodes)
        iftError("Error when writing the path values of the complete graph", "iftWriteCplGraph");

    if(fwrite(graph->ordered_nodes, sizeof(int), graph->nnodes, fp) != graph->nnodes)
        iftError("Error when writing the order of nodes of the complete graph", "iftWriteCplGraph");

    fclose(fp);

    // Writing the data graph's set
    path = iftJoinPathnames(2, tmp_dir, "dataset.zip");

    iftWriteOPFDataSet(graph->Z, path);

    iftFree(path);

    // Zippping the folder content
    iftZipDirContent(tmp_dir, pathname);

    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
}


void iftDestroyCplGraph(iftCplGraph **graph)
{
    iftCplGraph *aux = *graph;

    if (aux != NULL) {
        iftFree(aux->node);
        iftFree(aux->pathval);
        iftFree(aux->ordered_nodes);
        iftDestroyFHeap(&aux->Q);
        iftFree(aux);
        *graph = NULL;
    }
}


iftSet *iftPrototypesByMST(iftCplGraph *graph)
{
    iftSet  *S=NULL;
    int     s, t, u, v;
    iftDataSet *Z=graph->Z;
    float dist;


    // initialization

    for (u = 0; u < graph->nnodes; u++) {
        graph->pathval[u]  = IFT_INFINITY_FLT;
        s = graph->node[u].sample;
        Z->sample[s].label = 0;
        graph->node[u].pred  = IFT_NIL;
    }

    graph->pathval[0]   = 0;
    iftInsertFHeap(graph->Q, 0); // IFT_GRAY to 0

    // Minimum Spanning Tree by Degenerated IFT

    while ( !iftEmptyFHeap(graph->Q) ) {

        u = iftRemoveFHeap(graph->Q); // BLACK to u
        s = graph->node[u].sample;
        v = graph->node[u].pred;

        if (v != IFT_NIL){
            t = graph->node[v].sample;
            if (Z->sample[s].truelabel != Z->sample[t].truelabel){
                if (Z->sample[s].label == 0){
                    Z->sample[s].label = 1;
                    iftInsertSet(&S,u);
		    iftAddSampleStatus(&graph->Z->sample[s], IFT_PROTOTYPE);
                }
                if (Z->sample[t].label == 0){
                    Z->sample[t].label = 1;
                    iftInsertSet(&S,v);
		    iftAddSampleStatus(&graph->Z->sample[t], IFT_PROTOTYPE);
                }
            }
        }

        for (v=0; v < graph->nnodes; v++){
            if (graph->Q->color[v] != IFT_BLACK){
                t = graph->node[v].sample;
                if (iftDist == NULL)
                    dist = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
                else
                    dist = iftDist->distance_table[s][t];

                if ( dist < graph->pathval[v] ) {
                    graph->node[v].pred = u;
                    graph->pathval[v]   = dist;
                    if(graph->Q->color[v] == IFT_WHITE)
                        iftInsertFHeap(graph->Q, v); // IFT_GRAY TO v
                    else
                        iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
                }
            }
        }
    }

    iftResetFHeap(graph->Q);

    return(S);
}

iftSet* iftBipartitePrototypes(iftCplGraph* graph, int k) {

    iftDataSet* Z = graph->Z;

    iftSet* S = NULL;

    iftMatrix* adjacencyMatrix = iftCreateMatrix(graph->nnodes, graph->nnodes);
    iftSetFloatArray(adjacencyMatrix->val, adjacencyMatrix->n, IFT_INFINITY_FLT);

	#pragma omp parallel for
    for (int u = 0; u < graph->nnodes; ++u) {
      int s1 = graph->node[u].sample;
      for (int v = u+1; v < graph->nnodes; ++v) {
	int s2 = graph->node[v].sample;
	if(Z->sample[s1].truelabel != Z->sample[s2].truelabel) {
	  iftMatrixElem(adjacencyMatrix, u, v) = iftMatrixElem(adjacencyMatrix, v, u) = Z->iftArcWeight(Z->sample[s1].feat, Z->sample[s2].feat, Z->alpha, Z->nfeats);
	}
      }
    }

    
	int* vmin = iftAllocIntArray(k);
	float* min = iftAllocFloatArray(k);

    for (int u = 0; u < graph->nnodes; ++u) {
		iftSetFloatArray(min, k, IFT_INFINITY_FLT);
        for (int v = 0; v < graph->nnodes; ++v) {
            if(iftMatrixElem(adjacencyMatrix, u, v) < min[k-1]) {
                min[k-1] = iftMatrixElem(adjacencyMatrix, u, v);
                vmin[k-1] = v;

				iftFQuickSort(min, vmin, 0, k-1, IFT_INCREASING);
            }
        }

//		iftPrintFloatArray(min, k);

        for (int v = 0; v < graph->nnodes; ++v) {
            if (iftFindIntArrayElementIndex(vmin, k, v) == IFT_NIL) {
                iftMatrixElem(adjacencyMatrix, u, v) = IFT_INFINITY_FLT;
            }
        }
    }

	char* added = iftAllocCharArray(graph->nnodes);

    for (int u = 0; u < graph->nnodes; ++u) {
		for (int v = u + 1; v < graph->nnodes; ++v) {
			if (iftMatrixElem(adjacencyMatrix, u, v) != IFT_INFINITY_FLT
				&& iftMatrixElem(adjacencyMatrix, u, v) == iftMatrixElem(adjacencyMatrix, v, u)) {
				if(!added[u])
					iftInsertSet(&S, u);
				if(!added[v])
					iftInsertSet(&S, v);
				added[u] = added[v] = true;
			}
		}
	}

	iftFree(added);

    return S;
}


void  iftSupTrain2(iftCplGraph *graph)
{
//    iftSet *S = iftPrototypesByMST(graph);
	iftSet* S = iftBipartitePrototypes(graph, 10);
	int count = 0;

	iftSet* it = S;
	while(it!=NULL) {
		count+=1;
//		int sample = graph->node[it->elem].sample;
//		printf("%d (%d), ", it->elem, graph->Z->sample[sample].truelabel);
		it = it->next;
	}

	printf("# prototypes: %d\n", count);

	int s,t,i,u,v;
	float tmp, dist;
	iftDataSet *Z=graph->Z;

	// initialization

	for (u = 0; u < graph->nnodes; u++) {
		graph->pathval[u] = IFT_INFINITY_FLT;
		s = graph->node[u].sample;
		Z->sample[s].label  = 0;
		graph->node[u].pred = IFT_NIL;
	}
	while (S != NULL) {
		u = iftRemoveSet(&S);
		s = graph->node[u].sample;
		graph->pathval[u]     = 0;
		Z->sample[s].label    = Z->sample[s].truelabel;
		iftInsertFHeap(graph->Q, u); // IFT_GRAY to u
        graph->node[u].length = 0;
	}



	// IFT with fmax
	i=0;
	while ( !iftEmptyFHeap(graph->Q) ) {
		u = iftRemoveFHeap(graph->Q);
		graph->ordered_nodes[i]=u; i++;
		s = graph->node[u].sample;
		Z->sample[s].weight = graph->pathval[u];

		for (v=0; v < graph->nnodes; v++) {
			t = graph->node[v].sample;
			if (graph->Q->color[v] != IFT_BLACK && Z->sample[t].truelabel==Z->sample[s].truelabel) {
				if (graph->pathval[u] < graph->pathval[v]){
					if (iftDist == NULL){
						dist  = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
					}else{
						dist = iftDist->distance_table[s][t];
					}
					tmp  = iftMax(graph->pathval[u], dist);
					if ( tmp < graph->pathval[v] ) {
						graph->node[v].pred  = u;
                        graph->node[v].length = graph->node[u].length + 1;
						Z->sample[t].label   = Z->sample[s].label;
						graph->pathval[v]    = tmp;
						if(graph->Q->color[v] == IFT_WHITE)
							iftInsertFHeap(graph->Q, v); // IFT_GRAY to u
						else
							iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
					}
				}
			}
		}
	}

	iftResetFHeap(graph->Q);
	Z->ngroups = Z->nclasses;
}

void  iftSupTrain(iftCplGraph *graph) {
    iftSet *S = iftPrototypesByMST(graph);

//    iftSet* S = iftBipartitePrototypes(graph);

    int s, t, i, u, v;
    float tmp, dist;
    iftDataSet *Z = graph->Z;

    // initialization

    for (u = 0; u < graph->nnodes; u++) {
        graph->pathval[u] = IFT_INFINITY_FLT;
        s = graph->node[u].sample;
        Z->sample[s].label = 0;
        graph->node[u].pred = IFT_NIL;
    }
    while (S != NULL) {
        u = iftRemoveSet(&S);
        s = graph->node[u].sample;
        graph->pathval[u] = 0;
        Z->sample[s].label = Z->sample[s].truelabel;
        iftInsertFHeap(graph->Q, u); // IFT_GRAY to u
        graph->node[u].length = 0;
    }



    // IFT with fmax
    i = 0;
    while (!iftEmptyFHeap(graph->Q)) {
        u = iftRemoveFHeap(graph->Q);
        graph->ordered_nodes[i] = u;
        i++;
        s = graph->node[u].sample;
        Z->sample[s].weight = graph->pathval[u];

        for (v = 0; v < graph->nnodes; v++) {
            t = graph->node[v].sample;
            if (graph->Q->color[v] != IFT_BLACK && Z->sample[s].truelabel == Z->sample[t].truelabel) {
                if (graph->pathval[u] < graph->pathval[v]) {
                    if (iftDist == NULL) {
                        dist = Z->iftArcWeight(Z->sample[s].feat, Z->sample[t].feat, Z->alpha, Z->nfeats);
                    } else {
                        dist = iftDist->distance_table[s][t];
                    }
                    tmp = iftMax(graph->pathval[u], dist);
                    if (tmp < graph->pathval[v]) {
                        graph->node[v].pred = u;
                        graph->node[v].length = graph->node[u].length + 1;
                        Z->sample[t].label = Z->sample[s].label;
                        graph->pathval[v] = tmp;
                        if (graph->Q->color[v] == IFT_WHITE)
                            iftInsertFHeap(graph->Q, v); // IFT_GRAY to u
                        else
                            iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
                    }
                }
            }
        }
    }

    iftResetFHeap(graph->Q);
    Z->ngroups = Z->nclasses;
}

int iftClassifyByKNN(iftDataSet *Z1, iftDataSet *Z, int k)
{
    int s,t,n,nerrors=0,label,i,j,nn[k+1],nclass[Z1->nclasses+1];
    float d[k+1], dist;


    if ((k < 1)||(k > Z1->ntrainsamples))
        iftError("The value of k-nearest neighbors is invalid", "iftClassifyKNN");

    if (Z1->nfeats != Z->nfeats)
        iftError("Incompatible datasets", "iftClassifyByKNN");

    for (t = 0; t < Z->nsamples; t++){
        if (iftHasSampleStatus(Z->sample[t], IFT_TEST)){
            for (i=0; i <= k; i++) {
                d[i]  = IFT_INFINITY_FLT;
                nn[i] = IFT_NIL;
            }
            for (s = 0; s < Z1->nsamples; s++){
                if (iftHasSampleStatus(Z1->sample[s], IFT_TRAIN)){
                    if ((iftDist == NULL)||(Z1 != Z))
                        d[k] = Z1->iftArcWeight(Z1->sample[s].feat,Z->sample[t].feat, Z1->alpha,Z1->nfeats);
                    else
                        d[k] = iftDist->distance_table[s][t];

                    nn[k] = s;
                    j     = k;
                    while ((j > 0)&&(d[j]<d[j-1])){ /* sort in the increasing
                         order of distance */
                        dist   = d[j];   n       = nn[j];
                        d[j]   = d[j-1]; nn[j]   = nn[j-1];
                        d[j-1] = dist;   nn[j-1] = n;
                        j--;
                    }
                }
            }
            for (i=1; i <= Z->nclasses; i++)
                nclass[i]=0;
            for (i=0; i < k; i++)
                nclass[Z1->sample[nn[i]].truelabel]++;
            for (i=2,label=1; i <= Z->nclasses; i++) {
                if (nclass[i] > nclass[label])
                    label = i;
            }
            Z->sample[t].label = label;
            if (Z->sample[t].label != Z->sample[t].truelabel)
                nerrors ++;
        }
    }
    Z->ngroups = Z1->nclasses;

    return(nerrors);
}


int iftClassify(const iftCplGraph *graph, iftDataSet *Ztest) {
    iftDataSet *Ztrain = graph->Z;

    // checkers
    if (Ztrain->nfeats != Ztest->nfeats) {
        iftError("Train. and Test. Datasets with different Number of Feats: %d != %d",
                 "iftClassify", Ztrain->nfeats, Ztest->nfeats);
    }

    int n_errors = 0;

#pragma omp parallel for schedule(auto)
    for (int t = 0; t < Ztest->nsamples; t++) {
        if(iftHasSampleStatus(Ztest->sample[t], IFT_TEST)) {
            iftFindMinimumCostAndClassifySample(graph, Ztest, t);
        }
    }

#pragma omp parallel for schedule(auto)
    for (int t = 0; t < Ztest->nsamples; t++) {
        if (Ztest->sample[t].label != Ztest->sample[t].truelabel) {
            n_errors++;
            iftAddSampleStatus(&Ztest->sample[t], IFT_ERROR);
        }
    }

    Ztest->nclasses = Ztrain->nclasses;

    return n_errors;
}


void iftClassifyByOneClassOPF(iftDataSet *Ztest, const iftKnnGraph *graph,
                              float quantile_of_k) {
    if (quantile_of_k < 0.0 || quantile_of_k > 1.0) {
        iftError("Quantile value %d not in [0, 1]", "iftClassifyByOneClassOPF", quantile_of_k);
    }

    int rank = ceil(graph->k * quantile_of_k);
    if (rank == 0) { rank = 1; }
    if (rank > graph->k) { rank = graph->k; }
    // printf("rank %d/%d\n", rank, graph->k);
    
    iftDblArray *adj_radius = iftAdjacencyRadiusKnnGraphByRank(graph, rank);
    
    const iftDataSet *Ztrain = graph->Z;
    float square_sigma = (graph->maxarcw[graph->k] * graph->maxarcw[graph->k]) / 9.0;

    #pragma omp parallel for
    for (int t = 0; t < Ztest->nsamples; t++) {
        if (iftHasSampleStatus(Ztest->sample[t], IFT_TEST)) {
            Ztest->sample[t].group = -1;
            Ztest->sample[t].weight = 0.0;
            float min_dist = IFT_INFINITY_FLT;
            int closest_node = IFT_NIL;

            for (int i = 0; i < graph->nnodes; i++) {
                int u = graph->ordered_nodes[i];
                int s = graph->node[u].sample;

                float dist = 0.0;
                
                if (iftDist == NULL) {
                    dist = Ztrain->iftArcWeight(Ztrain->sample[s].feat,
                                                Ztest->sample[t].feat,
                                                Ztrain->alpha,
                                                Ztrain->nfeats);
                }
                else { dist = iftDist->distance_table[s][t]; }

                if (dist <= adj_radius->val[u]) {
                    Ztest->sample[t].group = Ztrain->sample[s].group;
                    Ztest->sample[t].weight = graph->pathval[u] * expf((-dist) / (2.0 * square_sigma));
                    Ztest->sample[t].label = IFT_POSITIVE_CLASS;
                    break;
                }
                // compute the closest node for the case of outliers
                else {
                    if (dist < min_dist) {
                        min_dist = dist;
                        closest_node = u;
                    }
                }
            }

            // t is an outlier
            if (Ztest->sample[t].group == -1) {
                iftAddSampleStatus(&Ztest->sample[t], IFT_OUTLIER);
                Ztest->sample[t].weight = graph->pathval[closest_node] * exp(-(min_dist) / (2.0 * square_sigma));

                int s = graph->node[closest_node].sample;
                Ztest->sample[t].group = Ztrain->sample[s].group;

                Ztest->sample[t].label = IFT_NEGATIVE_CLASS;
            }
        }
    }

    // compute the number of groups of the dataset
    for (int t = 0; t < Ztest->nsamples; t++) {
        if (Ztest->sample[t].group > Ztest->ngroups) {
            Ztest->ngroups = Ztest->sample[t].group;
        }
    }

    iftDestroyDblArray(&adj_radius);
}


int iftClassifyWithCertaintyValues(const iftCplGraph *graph, iftDataSet *Ztest) {
    iftDataSet *Ztrain = graph->Z;

    // checkers
    if (Ztrain->nfeats != Ztest->nfeats) {
        iftError("Train. and Test. Datasets with different Number of Feats: %d != %d",
                 "iftClassify", Ztrain->nfeats, Ztest->nfeats);
    }

    int n_errors = 0;

#pragma omp parallel for schedule(auto)
    for (int t = 0; t < Ztest->nsamples; t++) {
        // By default, we set training samples with weight 1.0
        Ztest->sample[t].weight = 1.0;

        // Otherwise, we compute the weight
        if(iftHasSampleStatus(Ztest->sample[t], IFT_TEST)) {
            float min_cost1 = iftFindMinimumCostAndClassifySample(graph, Ztest, t);
            float min_cost2 = iftFindMinimumCostFromAnotherClass(graph, Ztest, t);

            // assign certainty value
            if (!iftAlmostZero(min_cost1 + min_cost2))
                Ztest->sample[t].weight = min_cost2 / (min_cost1 + min_cost2);
            else
                Ztest->sample[t].weight = 0.5;
        }
    }

    for (int t = 0; t < Ztest->nsamples; t++) {
        if (Ztest->sample[t].label != Ztest->sample[t].truelabel) {
            n_errors++;
            iftAddSampleStatus(&Ztest->sample[t], IFT_ERROR);
        }
    }

    Ztest->nclasses = Ztrain->nclasses;

    return n_errors;
}


void iftClassifySample(iftCplGraph *graph, iftSample *sample){
    int u, s, i, label = 0;
    float tmp, weight, minCost;
    iftDataSet *Z1=graph->Z;

    i       = 0;
    u       = graph->ordered_nodes[i];
    s       = graph->node[u].sample;

    weight  = Z1->iftArcWeight(Z1->sample[s].feat,sample->feat,Z1->alpha,Z1->nfeats);
    minCost = iftMax(graph->pathval[u], weight);
    label   = Z1->sample[s].label;

    while((i < graph->nnodes-1)&&
          (minCost > graph->pathval[graph->ordered_nodes[i+1]])){

        u  = graph->ordered_nodes[i+1];
        s  = graph->node[u].sample;

        weight = Z1->iftArcWeight(Z1->sample[s].feat,sample->feat,Z1->alpha,Z1->nfeats);

        tmp = iftMax(graph->pathval[u], weight);
        if(tmp < minCost){
            minCost = tmp;
            label   = Z1->sample[s].label;
        }
        i++;
    }
    sample->label   = label;
    sample->weight  = minCost;
}

void  iftBorderClassify(iftCplGraph *graph, iftDataSet *Z, int truelabel)
{
    int u, s, t, i, indClass, indOtherClass;
    float tmp, weight, minCost[2];
    iftDataSet *Z1=graph->Z;

    if (Z1->nfeats != Z->nfeats)
        iftError("Incompatible datasets", "iftBinMembership");

    if (Z1->nclasses != 2)
        iftError("Number of classes is not 2", "iftBinMembership");

    for (t = 0; t < Z->nsamples; t++)
    {

        if ((Z==Z1)&&(iftHasSampleStatus(Z->sample[t], IFT_TRAIN)))
            continue;

        minCost[0] = IFT_INFINITY_FLT;
        minCost[1] = IFT_INFINITY_FLT;

        for (i=0; i < graph->nnodes; i++){

            u  = graph->ordered_nodes[i];
            s  = graph->node[u].sample;

            if ((iftDist==NULL)||(Z1 != Z))
                weight = Z1->iftArcWeight(Z1->sample[s].feat,Z->sample[t].feat,Z1->alpha,Z1->nfeats);
            else
                weight = iftDist->distance_table[s][t];

            if (Z1->sample[s].label==1) {
                tmp = iftMax(graph->pathval[u], weight);
                if(tmp < minCost[0]){
                    minCost[0] = tmp;
                }
            }else{
                tmp = iftMax(graph->pathval[u], weight);
                if(tmp < minCost[1]){
                    minCost[1] = tmp;
                }
            }

        }

        // index Class
        if ( truelabel == 1 ) {
            indClass = 0;
            indOtherClass = 1;
        }else{
            indClass = 1;
            indOtherClass = 0;
        }
        // Membership of the sample to the class
        Z->sample[t].weight = minCost[indOtherClass]/(minCost[indClass]+minCost[indOtherClass]+0.1);
        if (minCost[0] < minCost[1]) {
            Z->sample[t].label  = 1;
        }else{
            Z->sample[t].label  = 2;
        }

    }

    Z->ngroups = 2;
}


iftCplGraph *iftSupLearn(iftDataSet *Z)
{
    int          s, t, u, i, idx, nerrors=1;
    iftCplGraph *graph=NULL;


    if (Z->ntrainsamples == 0)
        iftError("No samples for training", "iftSupLearn");


    for (idx=1; (idx <= 10)&&(nerrors!=0); idx++) {
        if (graph != NULL) iftDestroyCplGraph(&graph);
        graph = iftCreateCplGraph(Z);
        iftSupTrain(graph);

        nerrors=iftClassify(graph,Z);
        for (s=0; (s < Z->nsamples)&&(nerrors!=0); s++) {
            if (iftHasSampleStatus(Z->sample[s], IFT_ERROR)){
                i = graph->nnodes-1;
                while (i >= 0){
                    u    = graph->ordered_nodes[i];
                    if (graph->node[u].pred != IFT_NIL) {
                        t    = graph->node[u].sample;
                        if (iftHasSampleStatus(Z->sample[t], IFT_TRAIN)){
                            if (Z->sample[t].truelabel == Z->sample[s].truelabel){
                                iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
                                iftSetSampleStatus(&Z->sample[t], IFT_TEST);
                                break;
                            }
                        }
                    }
                    i--;
                }
            }
        }
        //    printf("number of errors %d\n",nerrors);
    }

    if (graph != NULL) iftDestroyCplGraph(&graph);
    graph = iftCreateCplGraph(Z);
    iftSupTrain(graph);

    return(graph);
}





// Classify black MST nodes until /nsamples/ nodes are obtained whose edges belong to different labels
iftSet *iftGetMSTBoundarySamples(iftCplGraph *graph, iftMST *mst, int nsamples){
    iftSample *sample = mst->Z->sample;

    iftSet *selectedNodes = 0;

//  printf("iftGetMSTBoundarySamples\n");

    int u, setSize = 0;
    for (u = 0; (u < mst->nnodes) && (setSize < nsamples); u++) {
        int v = mst->node[u].maxarcadj;

        int s = mst->node[u].sample;
        int t = mst->node[v].sample;

        if(mst->node[u].color == IFT_WHITE || mst->node[v].color == IFT_WHITE){
            iftClassifySample(graph,&graph->Z->sample[s]);
            iftClassifySample(graph,&graph->Z->sample[t]);

            if (sample[s].label != sample[t].label){
//              printf("dist: %f\n",mst->Z->iftArcWeight(mst->Z->sample[s].feat,mst->Z->sample[t].feat,mst->Z->alpha,mst->Z->nfeats));
                if(mst->node[u].color == IFT_WHITE){
                    iftUnionSetElem(&selectedNodes,s);
                    mst->node[u].color = IFT_BLACK;
                    setSize++;
                }
                if(mst->node[v].color == IFT_WHITE && setSize < nsamples){
                    iftUnionSetElem(&selectedNodes,t);
                    mst->node[v].color = IFT_BLACK;
                    setSize++;
                }
            }
        }
    }

    // if the number of selected nodes is insufficient,
    // scan the list again and get white nodes
    for (u = 0; (setSize < nsamples) && (u < mst->nnodes); u++) {
        if (mst->node[u].color == IFT_WHITE) {
            iftUnionSetElem(&selectedNodes, mst->node[u].sample);
            mst->node[u].color = IFT_BLACK;
            setSize++;
        }

        const int v = mst->node[u].maxarcadj;
        if ((setSize < nsamples) && (mst->node[v].color == IFT_WHITE)) {
            iftUnionSetElem(&selectedNodes, mst->node[v].sample);
            mst->node[v].color = IFT_BLACK;
            setSize++;
        }
    }

    return selectedNodes;
}

/**
 * Compute the mean OPF accuracy of a k-fold cross validation using all samples
 * of the input dataset.
 * @param Z Input dataset.
 * @param k Number of partitions.
 * @return Mean OPF accuracy.
 */
float iftGetKFoldCrossValidationMeanOPFAccuracy(iftDataSet *Z, const int k) {
    // count the number of samples per class
    int nSamplesPerCls[Z->nclasses + 1];
    memset(nSamplesPerCls, 0, sizeof(int) * (Z->nclasses + 1));
    for (int i = 0; i < Z->nsamples; i++) {
        // sample class must be [1,Z->nclasses]
        if (Z->sample[i].truelabel <= 0 || Z->sample[i].truelabel > Z->nclasses) {
            iftError("Sample with invalid class!",
                     "iftGetKFoldCrossValidationMeanOPFAccuracy");
        }
        nSamplesPerCls[Z->sample[i].truelabel]++;
    }

    iftDataSet *Zfold[k];

    int nSamplesPerFold[k][Z->nclasses + 1];
    int nSamplesPerFoldSum[Z->nclasses + 1];
    memset(nSamplesPerFoldSum, 0, sizeof(int) * (Z->nclasses + 1));
    // create folds 1 to (k - 1)
    for (int i = 1; i < k; i++) {
        int nSamples = 0;
        for (int j = 1; j <= Z->nclasses; j++) {
            if (nSamplesPerFoldSum[j] == nSamplesPerCls[j]) {
                nSamplesPerFold[i][j] = 0;
            } else {
                int n = iftMax((nSamplesPerCls[j] / (float ) k) + 0.5f, 1.0f);
                if (n > (nSamplesPerCls[j] - nSamplesPerFoldSum[j])) {
                    n = nSamplesPerCls[j] - nSamplesPerFoldSum[j];
                }
                nSamplesPerFold[i][j] = n;
                nSamplesPerFoldSum[j] += n;
                nSamples += n;
            }
        }
        Zfold[i] = iftCreateDataSet(nSamples, Z->nfeats);
        Zfold[i]->ntrainsamples = 0;
        Zfold[i]->nclasses = Z->nclasses;
        Zfold[i]->iftArcWeight = Z->iftArcWeight;
        Zfold[i]->function_number = Z->function_number;
        iftCopyRefData(Zfold[i], Z->ref_data, Z->ref_data_type);
        // Zfold[i]->ref_data = Z->ref_data;
        Zfold[i]->fsp = iftCopyFeatSpaceParam(Z);
        memcpy(Zfold[i]->alpha, Z->alpha, sizeof(float) * Z->nfeats);
    }
    // create fold 0 with the remaining samples
    int nSamples = 0;
    for (int j = 1; j <= Z->nclasses; j++) {
        int n = nSamplesPerCls[j] - nSamplesPerFoldSum[j];
        if (n > 0) {
            nSamplesPerFold[0][j] = n;
            nSamplesPerFoldSum[j] += n;
            nSamples += n;
        } else {
            nSamplesPerFold[0][j] = 0;
        }
    }
    Zfold[0] = iftCreateDataSet(nSamples, Z->nfeats);
    Zfold[0]->ntrainsamples = 0;
    Zfold[0]->nclasses = Z->nclasses;
    Zfold[0]->iftArcWeight = Z->iftArcWeight;
    Zfold[0]->function_number = Z->function_number;
    iftCopyRefData(Zfold[0], Z->ref_data, Z->ref_data_type);
    // Zfold[0]->ref_data = Z->ref_data;
    Zfold[0]->fsp = iftCopyFeatSpaceParam(Z);
    memcpy(Zfold[0]->alpha, Z->alpha, sizeof(float) * Z->nfeats);

    // copy each training sample to a random fold
    for (int i = 0; i < Z->nsamples; i++) {
        int sampleCls = Z->sample[i].truelabel;
        bool found = false;
        do {
            int fold = iftRandomInteger(0, k - 1);
            if (nSamplesPerFold[fold][sampleCls] > 0) {
                int idx = Zfold[fold]->ntrainsamples++;
                float *tmp = Zfold[fold]->sample[idx].feat;
                Zfold[fold]->sample[idx] = Z->sample[i];
                Zfold[fold]->sample[idx].feat = tmp;
                memcpy(tmp, Z->sample[i].feat, sizeof(float) * Z->nfeats);
                nSamplesPerFold[fold][sampleCls]--;
                found = true;
            }
        } while (!found);
    }

    // do cross validation
    float acumAccuracy = 0;
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < k; i++) {
        iftDataSet *test = iftCopyDataSet(Zfold[i], true);
        for (int j = 0; j < test->nsamples; j++) {
            iftSetSampleStatus(&test->sample[j], IFT_TEST);
        }
        test->ntrainsamples = 0;

        iftDataSet *train = iftCreateDataSet(Z->nsamples - test->nsamples, Z->nfeats);
        train->nclasses = Z->nclasses;
        train->nfeats = Z->nfeats;
        train->iftArcWeight = Z->iftArcWeight;
        train->function_number = Z->function_number;
        iftCopyRefData(train, Z->ref_data, Z->ref_data_type);
        // train->ref_data = Z->ref_data;
        train->fsp = iftCopyFeatSpaceParam(Z);
        memcpy(train->alpha, Z->alpha, sizeof(float) * Z->nfeats);

        train->ntrainsamples = train->nsamples;

        // copy samples from all folds to the train dataset, skipping the test fold
        for (int j = 0, l = 0; j < k; j++) {
            // skip the test fold
            if (j == i)
                continue;

            for (int m = 0; m < Zfold[j]->nsamples; m++, l++) {
                float *tmp = train->sample[l].feat;
                train->sample[l] = Zfold[j]->sample[m];
                iftSetSampleStatus(&train->sample[l], IFT_TRAIN);
                train->sample[l].feat = tmp;
                memcpy(tmp, Zfold[j]->sample[m].feat, sizeof(float) * Z->nfeats);
            }
        }

        iftCplGraph *trainGraph = iftCreateCplGraph(train);
        iftSupTrain(trainGraph);

        iftClassify(trainGraph, test);

//    float accuracy = iftClassifAccuracy(test);
        float accuracy = iftOPFAccuracy(test, false);

#pragma omp critical
        {
            acumAccuracy += accuracy;
        }
        //printf("Accuracy fold %d: %f\n", i, accuracy);

        iftDestroyCplGraph(&trainGraph);
        iftDestroyDataSet(&train);
        iftDestroyDataSet(&test);
    }

    acumAccuracy /= k;

    for (int i = 0; i < k; i++) {
        iftDestroyDataSet(&Zfold[i]);
    }

    return acumAccuracy;
}

/**
 * Select at most nSelect samples from the input dataset, considering the
 * given status and using the following criteria:
 * 1) Get samples from cluster Ci (sampleLists[i]), closer to root
 *    ri (sampleLists[i][0]), whose labels given by the classifier (trainGraph)
 *    are different from the class of ri;
 * 2) If no samples are available in (1), select the most distant samples
 *    in Ci to root ri.
 * @param trainGraph Input classifier.
 * @param Z Input dataset.
 * @param sampleLists Array of size Z->ngroups+1 containing in each position
 *                    a list of samples of each label, sorted by their
 *                    increasing distance to the root node.
 * @param nSamples Array of size Z->ngroups+1 with the number of samples in
 *                 each list.
 * @param status Status of selected samples (equal or different).
 * @param statusCond Use 0 to select samples with the given status, or 1
 *                   to select samples without the given status.
 * @param nSelect Maximum number of samples to be selected.
 * @return Set with indexes of the selected samples with the given status.
 */
iftSet *iftGetRootDistanceBasedSamples(iftCplGraph *trainGraph, iftDataSet *Z,
                                       int **sampleLists, int *nSamples, uchar status,
                                       int statusCond, int nSelect) {
    int i, j;

    int *iSampleListsFwd = (int *) iftAlloc((Z->ngroups + 1), sizeof(int));
    int *iSampleListsBck = (int *) iftAlloc((Z->ngroups + 1), sizeof(int));
    int **selectedSamplesPerLabel = (int **) iftAlloc((Z->ngroups + 1), sizeof(int *));
    int *nSelectedSamplesPerLabel = (int *) iftAlloc(Z->ngroups + 1, sizeof(int));
    int nSelectedSamples = 0;

    for (i = 1; i <= Z->ngroups; i++) {
        iSampleListsFwd[i] = 1;
        iSampleListsBck[i] = nSamples[i] - 1;
        selectedSamplesPerLabel[i] = (int *) iftAlloc(nSelect, sizeof(int));
    }

    bool nodeFound;
    do {
        // for each label
        nodeFound = false;
#pragma omp parallel for schedule(dynamic) default(none) shared(trainGraph,Z,sampleLists,nSamples,statusCond,status,iSampleListsFwd,iSampleListsBck,selectedSamplesPerLabel,nSelectedSamplesPerLabel,nSelectedSamples,nodeFound) private(i,j)
        for (i = 1; i <= Z->ngroups; i++)   {
            // get the next available sample with label different from root
            int rootIdx = sampleLists[i][0];
            int nodeIdx = -1;
            int rootClass = Z->sample[rootIdx].truelabel;
            bool nodeFound2 = false;
            while (iSampleListsFwd[i] < nSamples[i]) {
                nodeIdx = sampleLists[i][iSampleListsFwd[i]];
                ++iSampleListsFwd[i];

                iftSample *sample = &(Z->sample[nodeIdx]);
                if (((statusCond == 0) && (iftHasSampleStatus(*sample, status)))
                    || ((statusCond == 1) && (!iftHasSampleStatus(*sample, status)))) {
                    // classify the sample and compare the label
                    iftClassifySample(trainGraph, sample);
                    if (sample->label != rootClass) {
                        nodeFound = nodeFound2 = true;
                        break;
                    }
                }
            }
            // if no sample found, get the most distant sample
            if (!nodeFound2) {
                while (iSampleListsBck[i] > 0) {
                    nodeIdx = sampleLists[i][iSampleListsBck[i]];
                    --iSampleListsBck[i];

                    iftSample *sample = &(Z->sample[nodeIdx]);
                    if (((statusCond == 0) && (iftHasSampleStatus(*sample, status)))
                        || ((statusCond == 1) && (!iftHasSampleStatus(*sample, status)))) {
                        // check if the sample was already selected
                        bool alreadySelected = false;
                        for (j = 0; j < nSelectedSamplesPerLabel[i]; j++) {
                            if (selectedSamplesPerLabel[i][j] == nodeIdx) {
                                alreadySelected = true;
                                break;
                            }
                        }
                        if (!alreadySelected) {
                            nodeFound = nodeFound2 = true;
                            break;
                        }
                    }
                }
            }
            // if all nodes are training samples, then go to next label
            if (!nodeFound2) {
                continue;
            }
            // add new sample
            selectedSamplesPerLabel[i][nSelectedSamplesPerLabel[i]] = nodeIdx;
            ++nSelectedSamplesPerLabel[i];
#pragma omp critical
            {
                ++nSelectedSamples;
            }
        }
    } while (nodeFound && (nSelectedSamples < nSelect));

    iftFree(iSampleListsFwd);
    iftFree(iSampleListsBck);

    iftSet *selectedSamples = NULL;
    nSelectedSamples = 0;
    for (i = 0; (i < nSelect) && (nSelectedSamples < nSelect); i++) {
        for (j = 1; (j <= Z->ngroups) && (nSelectedSamples < nSelect); j++) {
            if (i < nSelectedSamplesPerLabel[j]) {
                iftInsertSet(&selectedSamples, selectedSamplesPerLabel[j][i]);
                ++nSelectedSamples;
            }
        }
    }

    iftFree(nSelectedSamplesPerLabel);
    for (i = 1; i <= Z->ngroups; i++) {
        iftFree(selectedSamplesPerLabel[i]);
    }
    iftFree(selectedSamplesPerLabel);

    return selectedSamples;
}


/* This function performs OPF's semi-supervised algorithm directly on
   a dataset with supervised and unsupervised samples. The computation
   of a MST is not necessary. */

iftCplGraph *iftSemiSupTrain(iftDataSet *Ztrain) {
    iftCplGraph *graph = NULL;
    int i, u, v, s, t;
    float tmp, dist;


    graph = iftCreateCplGraph(Ztrain);

    // initialization
    for (u = 0; u < graph->nnodes; u++) {
      s = graph->node[u].sample;
      if (iftHasSampleStatus(graph->Z->sample[s], IFT_SUPERVISED)||(iftHasSampleStatus(graph->Z->sample[s], IFT_LABELPROPAGATED))){ 
	graph->pathval[u] = 0;
        graph->Z->sample[s].label = graph->Z->sample[s].truelabel;
      } else {	
	graph->pathval[u] = IFT_INFINITY_FLT;
      }
      iftInsertFHeap(graph->Q, u);
    }

    // IFT with fmax on a complete graph 
    i = 0;
    while (!iftEmptyFHeap(graph->Q)) {
        u = iftRemoveFHeap(graph->Q);
	s = graph->node[u].sample;
        graph->ordered_nodes[i] = u;
        i++;
        Ztrain->sample[s].weight = graph->pathval[u];

        for (v = 0;  v < graph->nnodes; v++){
	  if (graph->pathval[u] < graph->pathval[v]) {
	    t = graph->node[v].sample;
	    if (iftDist == NULL) {
	      dist = graph->Z->iftArcWeight(graph->Z->sample[s].feat, graph->Z->sample[t].feat, graph->Z->alpha, graph->Z->nfeats);
	    } else {
	      dist = iftDist->distance_table[s][t];
	    }
	    tmp = iftMax(graph->pathval[u], dist);
	    if (tmp < graph->pathval[v]) {
	      graph->node[v].pred = u;
	      graph->Z->sample[t].label = graph->Z->sample[s].label;
	      graph->pathval[v] = tmp;
	      iftGoUpFHeap(graph->Q, graph->Q->pos[v]);
	    }
	  }
        }
    }

    return(graph);
}


iftCplGraph *iftTrainImageClassifierByOPF(const iftImage *train_img, const iftLabeledSet *train_markers,
                                          float smooth_factor) {
    if (train_img == NULL)
        iftError("Training Image is NULL", "iftTrainImageClassifierByOPF");
    if (train_markers == NULL)
        iftError("Training Marker Set is NULL", "iftTrainImageClassifierByOPF");

    // gets a smooth image, if the smooth factor is > 0.0
    iftImage *smooth_img = NULL;
    if (smooth_factor > 0.0) {
        /* Smoothing image */
        iftKernel *K = NULL;
        if (iftIs3DImage(train_img))
            K = iftGaussianKernel(1.5, smooth_factor * iftMaximumValue(train_img));
	else
            K = iftGaussianKernel2D(1.5, smooth_factor * iftMaximumValue(train_img));

        smooth_img = iftLinearFilter(train_img, K);
        train_img  = smooth_img;

        iftDestroyKernel(&K);
    }

    iftMImage *mimg = NULL;
    if (!iftIsColorImage(train_img)) {
        mimg = iftImageToMImage(train_img, GRAY_CSPACE);
    } else {
        mimg = iftImageToMImage(train_img, LAB_CSPACE);
    }

    // gets the training dataset
    iftDataSet *Z = iftMImageSeedsToDataSet(mimg, train_markers);

    if (!iftSelectSupTrainSamplesWithNoOutliers(Z, 500))
        iftError("Seeds do not represent a valid training set", "main");

    iftCplGraph *graph = iftCreateCplGraph(Z);
    iftSupTrain(graph);

    // DESTROYERS
    if (smooth_img != NULL) {
        iftDestroyImage(&smooth_img);
    }

    return graph;
}


iftImage *iftClassifyImageByOPF(const iftImage *test_img, const iftCplGraph *clf, float smooth_factor,
                                int max_range, iftImage **fuzzy_map) {
    iftDataSet *Ztest = iftBuildTestImageDataSet(test_img, smooth_factor);

    iftClassifyWithCertaintyValues(clf, Ztest);

    iftImage *seg_img = iftDataSetToLabelImage(Ztest, NULL, true, IFT_CLASS);

    if (fuzzy_map != NULL) {
        *fuzzy_map = iftDataSetObjectMap(Ztest, NULL, max_range, -1); // -1 consider all labels
    }

    // DESTROYERS
    iftDestroyDataSet(&Ztest);

    return seg_img;
}
/******************************************************************/





















