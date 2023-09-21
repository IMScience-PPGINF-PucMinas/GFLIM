//
// Created by deangeli on 11/08/17.
//

#include "iftActiveLearning.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Stream.h"

//////////////////////////////////////SELECTOR//////////////////////////////////////////////////
iftGenericSelector* iftCreateGenericSelector(){
    iftGenericSelector* genericSelector = (iftGenericSelector*)iftAlloc(1,sizeof(iftGenericSelector));
    genericSelector->selector = NULL;
    genericSelector->freeFunctionSelector = NULL;

    genericSelector->preProcessingSelectorFunction = NULL;
    genericSelector->preProcessingOuput = NULL;
    genericSelector->selectedSamplesIndex = NULL;
    genericSelector->selecSamplesSelectorFunction = NULL;
    genericSelector->selectSamplesOutput = NULL;

    return  genericSelector;
}

void iftDestroyGenericSelector(iftGenericSelector** pGenericSelector){
    if(pGenericSelector == NULL || *pGenericSelector == NULL){
        return;
    }
    iftGenericSelector* aux= *pGenericSelector;
    iftDestroyGenericVector(&(aux->selectSamplesOutput));
    iftDestroyGenericVector(&(aux->preProcessingOuput));

    if(aux->freeFunctionSelector != NULL){
        aux->freeFunctionSelector(aux->selector);
    }
}

void iftPreProcessGenericSelector(iftGenericSelector* genericSelector){
    if(genericSelector->selector == NULL){
        iftError("Selector not defined","iftPreProcessSelector");
        return;
    }

    if(genericSelector->preProcessingSelectorFunction == NULL){
        return;
    }
    genericSelector->preProcessingSelectorFunction(genericSelector->selector);
}

iftIntArray* iftSelectSamplesGenericSelector(iftGenericSelector* genericSelector){
    if(genericSelector->selector == NULL){
        iftError("Selector not defined","iftPreProcessSelector");
        return NULL;
    }
    if(genericSelector->selecSamplesSelectorFunction == NULL){
        iftError("Selection function is not defined","iftSelectSamplesGenericSelector");
        return NULL;
    }
    genericSelector->selectedSamplesIndex = genericSelector->selecSamplesSelectorFunction(genericSelector->selector);
    return genericSelector->selectedSamplesIndex;
}


//////////RDS
iftSelectorRDS* iftCreateSelectorRDS(iftDataSet *dataSet){
    iftSelectorRDS* selectorRDS = (iftSelectorRDS*)calloc(1,sizeof(iftSelectorRDS));
    selectorRDS->dataSet = dataSet;
    selectorRDS->clusterOfSamplesIndex = NULL;
    selectorRDS->opfClusteringGrapgh = NULL;
    selectorRDS->samplesDistance2Root = NULL;
    selectorRDS->autoSetStatus = true;
    selectorRDS->numberClusters = 0;
    selectorRDS->OPFUNSUP_kmaxPercentage = 0.1;
    selectorRDS->iftGraphCutFun = iftNormalizedCut;
    selectorRDS->currentLearningCycle = 0;
    selectorRDS->selectedSamplesPerCycle = NULL;
    selectorRDS->numberOfSelectedSamplesPerCluster = 1;
    selectorRDS->destroyDataset = false;
    selectorRDS->destroyGraph = false;
    selectorRDS->destroyClassifier = false;
    selectorRDS->path2SamplesFolder = NULL;
    selectorRDS->createAndComputeGraph = true;
    selectorRDS->selectedSamplesIndex_typeSafe = NULL;
    return selectorRDS;
}

void iftPreProcessSelectorRDS(iftSelectorRDS* selectorRDS ){
    //float OPFUNSUP_kmaxPercentage = selectorRDS->OPFUNSUP_kmaxPercentage;
    //int maximumNumberIterations = selectorRDS->maximumNumberIterations;
    //iftKnnGraphCutFun iftGraphCutFun = selectorRDS->iftGraphCutFun;
    iftDataSet* dataSet = selectorRDS->dataSet;

    if(selectorRDS->autoSetStatus){
        iftSetStatus(dataSet,IFT_TRAIN);
        dataSet->ntrainsamples = dataSet->nsamples;
    }

    if(dataSet->iftArcWeight == NULL){
        dataSet->iftArcWeight = iftDistance1;
        dataSet->alpha = iftAlloc(dataSet->nfeats,sizeof(float));
        for (int i = 0; i < dataSet->nfeats; ++i) {
            dataSet->alpha[i] = 1;
        }
    }

    if(selectorRDS->createAndComputeGraph == true){
        selectorRDS->opfClusteringGrapgh = iftCreateKnnGraph(dataSet,selectorRDS->OPFUNSUP_kmaxPercentage*dataSet->nsamples);
        if(selectorRDS->dataSet->nclasses == 0){
            iftCountNumberOfClassesDataSet(selectorRDS->dataSet);
        }
        if(selectorRDS->dataSet->nclasses == 0){
            selectorRDS->numberClusters = iftUnsupTrainForRDS(selectorRDS->opfClusteringGrapgh, 20);
        }else{
            selectorRDS->numberClusters = iftUnsupTrainForRDS(selectorRDS->opfClusteringGrapgh, dataSet->nclasses);
        }
    }
    int * clusterSizes = iftAllocIntArray(selectorRDS->numberClusters);
    for (int sampleIndex = 0; sampleIndex < selectorRDS->dataSet->nsamples; ++sampleIndex) {
        clusterSizes[ (dataSet->sample[sampleIndex].group) - 1 ]++;
    }
    iftGenericVector* clusterOfSamplesIndex = iftCreateGenericNullVector(selectorRDS->numberClusters,sizeof(iftGenericVector*));
    float** clusterDistance = (float**)calloc(selectorRDS->numberClusters, sizeof(float *));
    for (int clusterIndex = 0; clusterIndex < selectorRDS->numberClusters; ++clusterIndex) {
        iftVectorAt(iftGenericVector*,clusterOfSamplesIndex,clusterIndex) = iftCreateGenericNullVector(clusterSizes[clusterIndex],sizeof(int));
        clusterDistance[clusterIndex] = (float*)calloc(clusterSizes[clusterIndex],sizeof(float));
    }
    // Fill sample indexes and their distance to the root
    int * currentSampleIndexInCluster = iftAllocIntArray(selectorRDS->numberClusters);
    for(int nodeIndex=0; nodeIndex < selectorRDS->opfClusteringGrapgh->nnodes; nodeIndex++){
        int sampleIndex = selectorRDS->opfClusteringGrapgh->node[selectorRDS->opfClusteringGrapgh->ordered_nodes[nodeIndex]].sample;
        int clusterIndex = dataSet->sample[sampleIndex].group-1;
        int sampleIndexInCluster = currentSampleIndexInCluster[clusterIndex];
        currentSampleIndexInCluster[clusterIndex]++;
        iftGenericVector* cluster = iftVectorAt(iftGenericVector*,clusterOfSamplesIndex,clusterIndex);
        iftVectorAt(int,cluster,sampleIndexInCluster) = sampleIndex;
        if (sampleIndexInCluster == 0) {
            //root
            clusterDistance[clusterIndex][sampleIndexInCluster] = 0.0;
        }else{
            int rootIndex = iftVectorAt(int,cluster,0);
            if (iftDist == NULL) {
                clusterDistance[clusterIndex][sampleIndexInCluster] = dataSet->iftArcWeight(dataSet->sample[sampleIndex].feat,
                                                                                            dataSet->sample[rootIndex].feat, dataSet->alpha, dataSet->nfeats);
            } else {
                clusterDistance[clusterIndex][sampleIndexInCluster] = iftDist->distance_table[sampleIndex][rootIndex];
            }
        }
    }
    //Sort cluster indexes by distance
    for (int clusterIndex = 0; clusterIndex < selectorRDS->numberClusters; ++clusterIndex){
        iftGenericVector* cluster = iftVectorAt(iftGenericVector*,clusterOfSamplesIndex,clusterIndex);
        iftFQuickSort(clusterDistance[clusterIndex], (int*) cluster->data, 1, clusterSizes[clusterIndex]-1, IFT_INCREASING);
    }
    selectorRDS->samplesDistance2Root = clusterDistance;
    selectorRDS->clusterOfSamplesIndex = clusterOfSamplesIndex;
    selectorRDS->numberOfSelectedSamplesPerCluster = 1;

    free(clusterSizes);
    free(currentSampleIndexInCluster);
    return;
}

void iftPreProcessSelectorRDSVoidPointer(void *selectorRDSVp){
    iftSelectorRDS* selectorRDS  = (iftSelectorRDS*)selectorRDSVp;
    iftPreProcessSelectorRDS(selectorRDS);
}

void iftSelecSamplesByRDS(iftSelectorRDS* selectorRDS){
    size_t numberOfSamples2Select = selectorRDS->numberOfSelectedSamplesPerCluster*selectorRDS->numberClusters;
    iftGenericVector* selectedSamples = iftCreateGenericVector(numberOfSamples2Select, selectorRDS->clusterOfSamplesIndex->elementSize);
    selectorRDS->selectedSamples = selectedSamples;
    int sampleIndex;
    iftDataSet* dataSet = selectorRDS->dataSet;

    if(selectorRDS->classifier == NULL){
        iftError("Classifier object not defined","iftSelecSamplesByRDS");
        return;
    }
    int numberOfSelectedSamples;
    iftDataSet* sampleDataset = iftCreateDataSetAux(1,dataSet->nfeats,false);
    for (int clusterId = 0; clusterId < selectorRDS->numberClusters; ++clusterId) {
        numberOfSelectedSamples = 0;
        iftGenericVector* cluster = iftVectorAt(iftGenericVector*,
                                                                 selectorRDS->clusterOfSamplesIndex,
                                                                 clusterId);

        while(numberOfSelectedSamples < selectorRDS->numberOfSelectedSamplesPerCluster){
            int rootSampleIndex = iftVectorAt(int,cluster,0);
            iftSample* rootSample = &(dataSet->sample[rootSampleIndex]);
            //you must supervise the root sample in cluster in order to supervise the remaining samples in cluster.
            if(rootSample->isSupervised == false){
                iftGenericVectorPushBack(int,selectedSamples,rootSampleIndex);
                break;
            }

            for (int i = 0; i < cluster->size; ++i) {
                sampleIndex = iftVectorAt(int,cluster,i);
                if(dataSet->sample[sampleIndex].isSupervised){
                    continue;
                }
                sampleDataset->sample[0].feat = dataSet->sample[sampleIndex].feat;
                iftPredictGenericClassifier(selectorRDS->classifier,sampleDataset);
                dataSet->sample[sampleIndex].label = sampleDataset->sample[0].label;
                if(rootSample->truelabel != sampleDataset->sample[0].label){
                    iftGenericVectorPushBack(int,selectedSamples,sampleIndex);
                    numberOfSelectedSamples++;
                }
                if(numberOfSelectedSamples == selectorRDS->numberOfSelectedSamplesPerCluster){
                    break;
                }
            }

            //if there are not samples in cluster
            if(numberOfSelectedSamples < selectorRDS->numberOfSelectedSamplesPerCluster){
                for (int i = cluster->size-1; i >= 0; --i) {
                    sampleIndex = iftVectorAt(int,cluster,i);
                    if(dataSet->sample[sampleIndex].isSupervised){
                        continue;
                    }
                    iftGenericVectorPushBack(int,selectedSamples,sampleIndex);
                    numberOfSelectedSamples++;
                    if(numberOfSelectedSamples == selectorRDS->numberOfSelectedSamplesPerCluster){
                        break;
                    }
                }
            }
            break;
        }



    }

    if(selectorRDS->selectedSamplesIndex_typeSafe != NULL){
        free(selectorRDS->selectedSamplesIndex_typeSafe->val);
        free(selectorRDS->selectedSamplesIndex_typeSafe);
    }
    
    selectorRDS->selectedSamplesIndex_typeSafe = iftAlloc(selectedSamples->size,sizeof(int));
    selectorRDS->selectedSamplesIndex_typeSafe->val = (int*)selectedSamples->data;
    selectorRDS->selectedSamplesIndex_typeSafe->n = selectedSamples->size;
}

iftIntArray* iftSelecSamplesByRDSVoidPointer(void *selectorRDSvp){
    iftSelectorRDS* selectorRDS  = (iftSelectorRDS*)selectorRDSvp;
    iftSelecSamplesByRDS(selectorRDS);
    return selectorRDS->selectedSamplesIndex_typeSafe;
}

void iftPrintSelectorInfoRDS(iftSelectorRDS *selectorRDS, bool printSampleName) {
    iftDataSet* dataSet = selectorRDS->dataSet;
    printf(" Knn-graph kmax: %f (%d out %d)\n",selectorRDS->OPFUNSUP_kmaxPercentage,
           (int)(selectorRDS->OPFUNSUP_kmaxPercentage*dataSet->nsamples),
           dataSet->nsamples);
    printf("Number of clusters: %d\n",selectorRDS->numberClusters);
    printf("Number of samples per cluster: %d\n",selectorRDS->numberOfSelectedSamplesPerCluster);
    int *countTrueLabelsInCluster = (int*)iftAlloc(dataSet->nclasses+1,sizeof(int));
    float maxDistance = 0;
    for (int i = 0; i < selectorRDS->clusterOfSamplesIndex->size; ++i) {
        maxDistance = 0;
        iftGenericVector* cluster = iftVectorAt(iftGenericVector*,selectorRDS->clusterOfSamplesIndex,i);
        for (int j = 0; j < dataSet->nclasses+1; ++j) {
            countTrueLabelsInCluster[j] = 0;
        }
        for (int k = 0; k < cluster->size; ++k){
            int sampleIndex = iftVectorAt(int,cluster,k);
            countTrueLabelsInCluster[dataSet->sample[sampleIndex].truelabel]++;
            if(selectorRDS->samplesDistance2Root[i][k] > maxDistance){
                maxDistance = selectorRDS->samplesDistance2Root[i][k];
            }
        }
        printf("cluster %d size: %lu => {",i,cluster->size);
        for (int j = 1; j < dataSet->nclasses+1; ++j){
            printf("%d ",countTrueLabelsInCluster[j]);
        }
        printf("} ");
        int root = iftVectorAt(int,cluster,0);
        printf(" root is sample %d (trueLabel: %d id: %d",root,dataSet->sample[root].truelabel,dataSet->sample[root].id);
        if(dataSet->ref_data_type == IFT_REF_DATA_FILESET){
            if(dataSet->ref_data != NULL){
                iftFileSet* fileSet = (iftFileSet*)dataSet->ref_data;
                int id = dataSet->sample[root].id;
                if(printSampleName == true){
                    printf(" name: %s",fileSet->files[id]->path);
                }
            }
        }
        printf(") ");
        printf( "max_Dist:%f ",maxDistance);
        printf("\n");
    }
    free(countTrueLabelsInCluster);
}

void iftPrintSelectedSamplesInfoRDS(iftSelectorRDS* selectorRDS){
    iftDataSet* dataSet = selectorRDS->dataSet;
    iftIntArray* indices = selectorRDS->selectedSamplesIndex_typeSafe;
    if(indices == NULL){
        return;
    }
    printf("Number of selected samples:%lu \n",indices->n);
    if(indices->val == NULL || indices->n == 0){
        return;
    }
    int c_clusterIndex = 0;
    int c_sample = 0;
    bool found = false;
    int root;
    while (selectorRDS->clusterOfSamplesIndex->size > c_clusterIndex && indices->n > c_sample){
        iftGenericVector* cluster = iftVectorAt(iftGenericVector*,selectorRDS->clusterOfSamplesIndex,c_clusterIndex);
        root = iftVectorAt(int,cluster,0);
        for (int i = 0; i < cluster->size; ++i) {
            found = false;

            if(indices->val[c_sample] == iftVectorAt(int,cluster,i)){
                int sampleIndex = indices->val[c_sample];
                float dist = selectorRDS->samplesDistance2Root[c_clusterIndex][i];
                printf("sampleIndex:%d| trueLabel:%d| id:%d| label:%d| cluster:%d| dist:%f| rootTrueLabel:%d| rootLabel:%d| ",sampleIndex,dataSet->sample[sampleIndex].truelabel,dataSet->sample[sampleIndex].id,
                       dataSet->sample[sampleIndex].label,c_clusterIndex,dist,dataSet->sample[root].truelabel,dataSet->sample[root].label);
                if(dataSet->ref_data_type == IFT_REF_DATA_FILESET){
                    if(dataSet->ref_data != NULL){
                        iftFileSet* fileSet = (iftFileSet*)dataSet->ref_data;
                        int id = dataSet->sample[sampleIndex].id;
                        printf(" name: %s",fileSet->files[id]->path);
                    }
                }
                printf(")\n");
                c_sample++;
                found = true;
                break;
            }
        }
        if(found == false){
            c_clusterIndex++;
        }
    }
}

void iftDestroySelectorRDS(iftSelectorRDS** selectorRDS){
    if(selectorRDS == NULL || *selectorRDS == NULL){
        return;
    }
    iftSelectorRDS* aux = *selectorRDS;
    if(aux->selectedSamplesPerCycle != NULL){
        iftDestroyGenericVector(&(aux->selectedSamplesPerCycle));
        (*selectorRDS)->selectedSamplesPerCycle = NULL;
    }
    if(aux->clusterOfSamplesIndex != NULL){
        iftDestroyGenericVector(&(aux->clusterOfSamplesIndex));
        (*selectorRDS)->clusterOfSamplesIndex = NULL;
    }
    if(aux->destroyDataset){
        iftDestroyDataSet(&(aux->dataSet));
        (*selectorRDS)->dataSet = NULL;
    }
    if(aux->destroyGraph){
        iftDestroyKnnGraph(&(aux->opfClusteringGrapgh));
        for (int i = 0; i < aux->numberClusters; ++i) {
            free(aux->samplesDistance2Root[i]);
        }
        free(aux->samplesDistance2Root);
        (*selectorRDS)->samplesDistance2Root = NULL;
    }
    if(aux->destroyClassifier){
        iftDestroyGenericClassifier(&(aux->classifier));
        (*selectorRDS)->classifier = NULL;
    }
    iftFree((*selectorRDS));
    *selectorRDS = NULL;
}

void iftDestroySelectorRDSVoidPointer(void* selectorRDSvp){
    iftSelectorRDS* selectorRDS  = (iftSelectorRDS*)selectorRDSvp;
    iftDestroySelectorRDS(&selectorRDS);
}

int iftUnsupTrainForRDS(iftKnnGraph *graph,int nclasses){
    int ngroups = 1;
    /* for each value of K in the interval [1,KMAX]*/
    int iter = 0;
    bool foundGoodClusters = false;
    while(foundGoodClusters == false){

        for (int k = graph->kmax; (k >= 1)&&(ngroups < 2*nclasses); k--) {
            graph->k = k;                   /* test the new value of K*/
            iftPDFByRange(graph);           /* compute the PDF with the new value of K*/
            ngroups = iftUnsupOPF(graph);             /* cluster the samples*/
            iter++;
            //printf("%d %d %d\n",ngroups,iter,k);
        }

        foundGoodClusters = true;
    }


    return(ngroups);
}

void iftComputeKnnGraphStatistics(iftDataSet* dataSet, iftDoubleMatrix** clusterMeanp,iftDoubleMatrix** clusterStandardDeviationp){
    if (dataSet == NULL){
        return;
    }
    int nclusters = dataSet->ngroups + 1;
    *clusterMeanp = iftCreateDoubleMatrix(dataSet->nfeats,nclusters);
    *clusterStandardDeviationp = iftCreateDoubleMatrix(dataSet->nfeats,nclusters);
    iftDoubleMatrix* clusterMean = *clusterMeanp;
    iftDoubleMatrix* clusterStandardDeviation = *clusterStandardDeviationp;

    int *countSamplesPerCluster = iftAlloc(nclusters,sizeof(int));
    int groupIndex = 0;

    //mean
    for (int sampleIndex = 0; sampleIndex < dataSet->nsamples; ++sampleIndex) {
        groupIndex = dataSet->sample[sampleIndex].group;
        for (int j = 0; j < dataSet->nfeats; ++j) {
            clusterMean->val[clusterMean->tbrow[groupIndex] + j] += dataSet->sample[sampleIndex].feat[j];
        }
        countSamplesPerCluster[groupIndex] += 1;
    }

    for (int clusterIndex = 0; clusterIndex < nclusters; ++clusterIndex) {
        for (int j = 0; j < dataSet->nfeats; ++j) {
            clusterMean->val[clusterMean->tbrow[clusterIndex] + j] /= (countSamplesPerCluster[clusterIndex] + 0.00000000001);
        }
    }
    //

    //std
    for (int sampleIndex = 0; sampleIndex < dataSet->nsamples; ++sampleIndex) {
        groupIndex = dataSet->sample[sampleIndex].group;
        for (int j = 0; j < dataSet->nfeats; ++j) {
            double var = clusterMean->val[clusterMean->tbrow[groupIndex] + j] - dataSet->sample[sampleIndex].feat[j];
            var = var*var;
            clusterStandardDeviation->val[clusterStandardDeviation->tbrow[groupIndex] + j] += var;
        }
    }
    for (int clusterIndex = 0; clusterIndex < nclusters; ++clusterIndex) {
        for (int j = 0; j < dataSet->nfeats; ++j) {
            double value = clusterStandardDeviation->val[clusterStandardDeviation->tbrow[clusterIndex] + j];
            value /= (countSamplesPerCluster[clusterIndex] + 0.0000000001);
            value = sqrt(value);
            clusterStandardDeviation->val[clusterStandardDeviation->tbrow[clusterIndex] + j] = value;
        }
    }
    printf("\n");
    printf("mean\n");
    iftPrintDoubleMatrix(clusterMean);
    printf("\n");

    printf("std\n");
    iftPrintDoubleMatrix(clusterStandardDeviation);
    printf("\n");

}

/////////////////

///RANDOM
iftSelectorRandom* iftCreateSelectorRandom(iftDataSet *dataSet){
    iftSelectorRandom* randomSelector = iftAlloc(1,sizeof(iftSelectorRandom));
    randomSelector->dataSet = dataSet;
    randomSelector->shuffledIndices = iftAlloc(dataSet->nsamples,sizeof(int));
    randomSelector->shuffledIndicesSize = dataSet->nsamples;
    randomSelector->selectedSamples = NULL;
    randomSelector->clearDatasetWhenDestroyed = false;
    randomSelector->selectedSamplesIndex_typeSafe = NULL;
    if(dataSet->nclasses > 0){
        randomSelector->numberOfSelectedSamples = dataSet->nclasses*2;
    }else{
        randomSelector->numberOfSelectedSamples = 4;
    }
    for (int i = 0; i < dataSet->nsamples; ++i) {
        randomSelector->shuffledIndices[i] = i;
    }
    iftShuffleIntArray(randomSelector->shuffledIndices,dataSet->nsamples);
    return randomSelector;
}

void iftSelecSamplesByRandom(iftSelectorRandom* selectorRandom){
    selectorRandom->selectedSamples = iftCreateGenericVector(selectorRandom->numberOfSelectedSamples,sizeof(int));
    iftDataSet* dataSet = selectorRandom->dataSet;
    int sampleIndex;
    for (int i = 0; i < selectorRandom->shuffledIndicesSize; ++i) {
        sampleIndex = selectorRandom->shuffledIndices[i];
        if(dataSet->sample[sampleIndex].isSupervised == false){
            iftGenericVectorPushBack(int,selectorRandom->selectedSamples,sampleIndex);
        }
        if(selectorRandom->selectedSamples->size >= selectorRandom->numberOfSelectedSamples){
            break;
        }
    }
    if(selectorRandom->selectedSamplesIndex_typeSafe != NULL){
        free(selectorRandom->selectedSamplesIndex_typeSafe->val);
        free(selectorRandom->selectedSamplesIndex_typeSafe);
    }
    selectorRandom->selectedSamplesIndex_typeSafe = iftAlloc(selectorRandom->selectedSamples->size,sizeof(int));
    selectorRandom->selectedSamplesIndex_typeSafe->val = (int*)selectorRandom->selectedSamples->data;
    selectorRandom->selectedSamplesIndex_typeSafe->n = selectorRandom->selectedSamples->size;
}

iftIntArray* iftSelecSamplesByRandomVoidPointer(void* selectorRandomvp){
    iftSelectorRandom* selectorRandom = (iftSelectorRandom*)selectorRandomvp;
    iftSelecSamplesByRandom(selectorRandom);
    return selectorRandom->selectedSamplesIndex_typeSafe;
}

void iftDestroySelectorRandomVoidPointer(void* selectorRandomvp){
    iftSelectorRandom* selectorRandom = (iftSelectorRandom*)selectorRandomvp;
    iftDestroySelectorRandom(&selectorRandom);
}

void iftDestroySelectorRandom(iftSelectorRandom** selectorRandom){
    if(selectorRandom == NULL || *selectorRandom == NULL){
        return;
    }
    iftSelectorRandom* aux = *selectorRandom;
    free(aux->shuffledIndices);
    iftDestroyGenericVector(&(aux->selectedSamples));
    if(aux->clearDatasetWhenDestroyed == true){
        iftDestroyDataSet(&(aux->dataSet));
    }
    free(*selectorRandom);
    *selectorRandom = NULL;
}

////////////

/////RWS
iftSelectorRWS* iftCreateSelectorRWS(iftDataSet *dataSet){
    iftSelectorRWS* selectorRWS = (iftSelectorRWS*)calloc(1,sizeof(iftSelectorRWS));
    selectorRWS->dataSet = dataSet;
    selectorRWS->createAndComputeGraph = true;
    selectorRWS->destroyDataset = false;
    selectorRWS->destroyGraph = false;
    selectorRWS->opfClusteringGrapgh = NULL;
    selectorRWS->iftGraphCutFun = iftNormalizedCut;
    selectorRWS->autoSetStatus = true;
    selectorRWS->OPF_kmaxPercentage = 0.05;
    selectorRWS->numberOfSelectedSamplesPerCluster = 1;
    return selectorRWS;
}

void iftPreProcessSelectorRWS(iftSelectorRWS* selectorRWS){
    iftKnnGraphCutFun iftGraphCutFun = selectorRWS->iftGraphCutFun;
    iftDataSet* dataSet = selectorRWS->dataSet;
    if(selectorRWS->autoSetStatus){
        iftSetStatus(dataSet,IFT_TRAIN);
        dataSet->ntrainsamples = dataSet->nsamples;
    }
    if(selectorRWS->createAndComputeGraph == true){
        selectorRWS->opfClusteringGrapgh = iftCreateKnnGraph(dataSet,selectorRWS->OPF_kmaxPercentage*dataSet->nsamples);
        selectorRWS->numberClusters = iftUnsupTrain(selectorRWS->opfClusteringGrapgh, iftGraphCutFun);
    }
    int * clusterSizes = iftAllocIntArray(selectorRWS->numberClusters);
    for (int sampleIndex = 0; sampleIndex < selectorRWS->dataSet->nsamples; ++sampleIndex) {
        clusterSizes[ (dataSet->sample[sampleIndex].group) - 1 ]++;
    }
    iftGenericVector* clusterOfSamplesIndex = iftCreateGenericNullVector(selectorRWS->numberClusters,sizeof(iftGenericVector*));
    for (int clusterIndex = 0; clusterIndex < selectorRWS->numberClusters; ++clusterIndex) {
        iftVectorAt(iftGenericVector*,clusterOfSamplesIndex,clusterIndex) = iftCreateGenericNullVector(clusterSizes[clusterIndex],sizeof(int));
    }
    // Fill in sample indexes and their distance to the root
    int * currentSampleIndexInCluster = iftAllocIntArray(selectorRWS->numberClusters);
    for(int nodeIndex=0; nodeIndex < selectorRWS->opfClusteringGrapgh->nnodes; nodeIndex++){
        int sampleIndex = selectorRWS->opfClusteringGrapgh->node[selectorRWS->opfClusteringGrapgh->ordered_nodes[nodeIndex]].sample;
        int clusterIndex = dataSet->sample[sampleIndex].group-1;
        int sampleIndexInCluster = currentSampleIndexInCluster[clusterIndex];
        currentSampleIndexInCluster[clusterIndex]++;
        iftGenericVector* cluster = iftVectorAt(iftGenericVector*,clusterOfSamplesIndex,clusterIndex);
        iftVectorAt(int,cluster,sampleIndexInCluster) = sampleIndex;
    }
    selectorRWS->clusterOfSamplesIndex = clusterOfSamplesIndex;
    free(clusterSizes);
    free(currentSampleIndexInCluster);
}

iftIntArray* iftSelecSamplesByRWSVoidPointer(void *selectorRWSvp){
    iftSelectorRWS* selectorRWS = (iftSelectorRWS*)selectorRWSvp;
    iftSelecSamplesByRWS(selectorRWS);
    return selectorRWS->selectedSamplesIndex_typeSafe;
}

void iftSelecSamplesByRWS(iftSelectorRWS* selectorRWS) {
    iftDataSet *dataSet = selectorRWS->dataSet;
    if (selectorRWS->classifier == NULL) {
        iftError("Classifier object not defined", "iftSelecSamplesByRDS");
        return;
    }
    size_t numberOfSamples2Select = selectorRWS->numberClusters;
    iftGenericVector *selectedSamples = iftCreateGenericVector(numberOfSamples2Select,
                                                               selectorRWS->clusterOfSamplesIndex->elementSize);
    //selectorRWS->selectedSamples = selectedSamples;
    int numberOfSelectedSamples;
    int sampleIndex;
    iftDataSet *sampleDataset = iftCreateDataSetAux(1, dataSet->nfeats, false);
    for (int clusterId = 0; clusterId < selectorRWS->numberClusters; ++clusterId) {
        numberOfSelectedSamples = 0;
        iftGenericVector *cluster = iftVectorAt(iftGenericVector*,
                                                                 selectorRWS->clusterOfSamplesIndex,
                                                                 clusterId);

        while (numberOfSelectedSamples < selectorRWS->numberOfSelectedSamplesPerCluster) {
            int rootSampleIndex = iftVectorAt(int, cluster, 0);
            iftSample *rootSample = &(dataSet->sample[rootSampleIndex]);
            if (rootSample->isSupervised == false) {
                iftGenericVectorPushBack(int, selectedSamples, rootSampleIndex);
                break;
            }
            for (int i = 0; i < cluster->size; ++i) {
                sampleIndex = iftVectorAt(int, cluster, i);
                if (dataSet->sample[sampleIndex].isSupervised) {
                    continue;
                }
                sampleDataset->sample[0].feat = dataSet->sample[sampleIndex].feat;
                iftPredictGenericClassifier(selectorRWS->classifier, sampleDataset);
                dataSet->sample[sampleIndex].label = sampleDataset->sample[0].label;
                if (rootSample->truelabel != sampleDataset->sample[0].label) {
                    iftGenericVectorPushBack(int, selectedSamples, sampleIndex);
                    numberOfSelectedSamples++;
                }
                if (numberOfSelectedSamples == selectorRWS->numberOfSelectedSamplesPerCluster) {
                    break;
                }
            }
        }
    }
    if(selectorRWS->selectedSamplesIndex_typeSafe != NULL){
        free(selectorRWS->selectedSamplesIndex_typeSafe->val);
        free(selectorRWS->selectedSamplesIndex_typeSafe);
        selectorRWS->selectedSamplesIndex_typeSafe = NULL;
    }
    selectorRWS->selectedSamplesIndex_typeSafe = iftAlloc(selectedSamples->size,sizeof(int));
    selectorRWS->selectedSamplesIndex_typeSafe->val = (int*)selectedSamples->data;
    selectorRWS->selectedSamplesIndex_typeSafe->n = selectedSamples->size;
    free(selectedSamples);
}
///////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////ACTIVE LEARNING//////////////////////////////////////////////////



iftCommonActiveLearning* iftCreateCommonActiveLearning(){
    iftCommonActiveLearning* activeLearning = iftAlloc(1,sizeof(iftCommonActiveLearning));
    activeLearning->classifier = NULL;
    activeLearning->selector = NULL;
    activeLearning->currentLearningCycle = 0;
    activeLearning->nextState = ACTIVELEARNING_BEGIN;
    activeLearning->functionUsed2PreProcessing = NULL;
    activeLearning->functionUsed2PosProcessing = NULL;

    activeLearning->beginActiveLearningFunction = iftInitializeCommonActiveLearning;
    activeLearning->preProcessingActiveLearningFunction = iftPreprocessCommonActiveLearning;
    activeLearning->selecSamplesActiveLearningFunction = iftSelecSamplesCommonActiveLearning;
    activeLearning->supervisingActiveLearningFunction = iftSuperviseSamplesgCommonActiveLearning;
    activeLearning->trainActiveLearningFunction = iftTrainClassifierCommonActiveLearning;
    activeLearning->posProcessingActiveLearningFunction = iftPosProcessCommonActiveLearning;
    activeLearning->endActiveLearningFunction = NULL;
    return activeLearning;
}

void iftDestroyCommonActiveLearning(iftCommonActiveLearning** activeLearningObject){
    if(activeLearningObject == NULL || *activeLearningObject == NULL){
        return;
    }
    iftCommonActiveLearning* aux = *activeLearningObject;
    if(aux->selector != NULL){
        iftDestroyGenericSelector(&(aux->selector));
        (*activeLearningObject)->selector = NULL;
    }
    if(aux->classifier != NULL){
        iftDestroyGenericClassifier(&(aux->classifier));
        (*activeLearningObject)->classifier = NULL;
    }
    iftFree((*activeLearningObject));
    (*activeLearningObject) = NULL;
}

void iftInitializeCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning){
    commonActiveLearning->nextState = ACTIVELEARNING_PREPROCESSING_STATE;
    commonActiveLearning->currentLearningCycle = 0;
    if(commonActiveLearning->preProcessingActiveLearningFunction == NULL){
        commonActiveLearning->nextState = ACTIVELEARNING_SELECT_SAMPLES_STATE;
    }

    if(commonActiveLearning->selector == NULL){
        iftError("selector is not defined","iftInitCommonActiveLearning");
        commonActiveLearning->nextState = ACTIVELEARNING_ENDING;
    }
    if(commonActiveLearning->classifier == NULL){
        iftError("classifier is not defined","iftInitCommonActiveLearning");
        commonActiveLearning->nextState = ACTIVELEARNING_ENDING;
    }
}

void iftPreprocessCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning){
    commonActiveLearning->nextState = ACTIVELEARNING_SELECT_SAMPLES_STATE;
    if(commonActiveLearning->selector == NULL){
        iftError("selector is not defined","iftPreprocessingCommonActiveLearning");
        commonActiveLearning->nextState = ACTIVELEARNING_ENDING;
    }else{
        if(commonActiveLearning->functionUsed2PreProcessing != NULL){
            commonActiveLearning->functionUsed2PreProcessing(commonActiveLearning,
                                                             commonActiveLearning->selector,
                                                             commonActiveLearning->classifier);
        }
    }
}

void iftSelecSamplesCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning){
    commonActiveLearning->nextState = ACTIVELEARNING_SUPERVISING_STATE;
    if(commonActiveLearning->selector == NULL){
        iftError("selector is not defined","iftPreprocessingCommonActiveLearning");
    }else{
        if(commonActiveLearning->functionUsed2SelectSamples != NULL){
            commonActiveLearning->functionUsed2SelectSamples(commonActiveLearning,
                                                             commonActiveLearning->selector,
                                                             commonActiveLearning->classifier);
        }else{
            iftError("selection behaviour is not defined (functionUsed2PreProcessing)",
                     "iftSelecSamplesCommonActiveLearning");
        }
    }
}

void iftSuperviseSamplesgCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning){
    commonActiveLearning->nextState = ACTIVELEARNING_CLASSIFIER_TRAIN_STATE;
    if(commonActiveLearning->functionUsed2SuperviseSamples == NULL){
        iftError("function to supervise samples is not defined","iftSupervisingCommonActiveLearning");
        commonActiveLearning->nextState = ACTIVELEARNING_ENDING;
    }else{
        commonActiveLearning->functionUsed2SuperviseSamples(commonActiveLearning,
                                                            commonActiveLearning->selector,
                                                            commonActiveLearning->classifier);
    }

}

void iftTrainClassifierCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning){
    commonActiveLearning->nextState = ACTIVELEARNING_POSPROCESSING_STATE;
    if(commonActiveLearning->posProcessingActiveLearningFunction == NULL){
        commonActiveLearning->nextState = ACTIVELEARNING_SELECT_SAMPLES_STATE;
    }
    if(commonActiveLearning->classifier == NULL){
        iftError("classifier is not defined","iftTrainClassifierCommonActiveLearning");
    }else{
        if(commonActiveLearning->functionUsed2TrainClassifier == NULL){
            iftError("function to train classifier is not defined","iftTrainClassifierCommonActiveLearning");
        }else{
            commonActiveLearning->functionUsed2TrainClassifier(commonActiveLearning,
                                                               commonActiveLearning->selector,
                                                               commonActiveLearning->classifier);
        }
    }
}

void iftPosProcessCommonActiveLearning(iftCommonActiveLearning* commonActiveLearning){
    commonActiveLearning->nextState = ACTIVELEARNING_SELECT_SAMPLES_STATE;
    if(commonActiveLearning->functionUsed2PosProcessing != NULL){
        commonActiveLearning->functionUsed2PosProcessing(commonActiveLearning,
                                                         commonActiveLearning->selector,
                                                         commonActiveLearning->classifier);
    }
}





void iftSelecSamplesSelectorRDS(iftCommonActiveLearning* activeLearning,
                                iftGenericSelector* selectorRDS,
                                iftGenericClassifier* classifier){
    iftUnusedParameter(classifier);
    if(selectorRDS == NULL){
        iftError("generic selector is NULL","iftPreprocessSelectorRDS");
        activeLearning->nextState = ACTIVELEARNING_ENDING;
    }
    if(selectorRDS->selector == NULL){
        iftError("RDS selector is NULL","iftPreprocessSelectorRDS");
        activeLearning->nextState = ACTIVELEARNING_ENDING;
    }
    iftSelectorRDS* rds = (iftSelectorRDS*)selectorRDS->selector;
    iftSelecSamplesByRDS(rds);
}

void iftCreateFileSampleNameRDS(iftSelectorRDS* rds,iftSample* sample){
    sprintf(rds->currentFileName,"%s/%06d_%08d.%s",rds->path2SamplesFolder,sample->truelabel,sample->id,rds->fileExtension);
    //printf("%s\n",rds->currentFileName);
}

void iftCreateTextFileInfoForSelectedSamplesRDS(iftSelectorRDS* rds){
    FILE *selectedSamplesFileInfo;
    char fileName[50];
    printf("supervising samples (RDS)... creating text files\n");
    sprintf(fileName,"selectedSamples_%d.txt",rds->currentLearningCycle);
    selectedSamplesFileInfo = fopen(fileName, "w");
    //bool isSampleFile = true;
    for (size_t i = 0; i < rds->selectedSamples->size; ++i) {
        int sampleIndex = iftVectorAt(int,rds->selectedSamples,i);
        iftSample* sample = &(rds->dataSet->sample[sampleIndex]);
        if(true){
            iftCreateFileSampleNameRDS(rds,sample);
        }
        fprintf(selectedSamplesFileInfo, "%03lu sampleIndex: %08d  sampleId: %08d trueLabel: %06d Label: %06d ref: %s\n",
                i, sampleIndex,sample->id,sample->truelabel,sample->label,rds->currentFileName);
    }
    fclose(selectedSamplesFileInfo);
}

void iftSuperviseSamplesSelectorRDS(iftCommonActiveLearning* activeLearning,
                                    iftGenericSelector* selectorRDS,
                                    iftGenericClassifier* classifier){
    iftSelectorRDS* rds = selectorRDS->selector;
//    bool isSampleFile = true;
//    bool getSampleFile = true;
    iftCreateTextFileInfoForSelectedSamplesRDS(rds);
}



void iftCommonActiveLearningGoToNextState(iftCommonActiveLearning* activeLearning){
    switch (activeLearning->nextState){
        case ACTIVELEARNING_BEGIN:
            if(activeLearning->beginActiveLearningFunction == NULL){
                iftWarning("begin function is not defined","iftCommonActiveLearningGoToNextState");
                activeLearning->nextState =  ACTIVELEARNING_PREPROCESSING_STATE;
            }else{
                activeLearning->beginActiveLearningFunction(activeLearning);
            }
            break;
        case ACTIVELEARNING_PREPROCESSING_STATE:
            if(activeLearning->preProcessingActiveLearningFunction == NULL){
                iftWarning("pre process function is not defined","iftCommonActiveLearningGoToNextState");
                activeLearning->nextState =  ACTIVELEARNING_SELECT_SAMPLES_STATE;
            }else{
                activeLearning->preProcessingActiveLearningFunction(activeLearning);
            }
            break;
        case ACTIVELEARNING_SELECT_SAMPLES_STATE:
            if(activeLearning->selecSamplesActiveLearningFunction == NULL){
                iftError("select samples function is not defined","iftCommonActiveLearningGoToNextState");
                activeLearning->nextState =  ACTIVELEARNING_ENDING;
            }else{
                activeLearning->selecSamplesActiveLearningFunction(activeLearning);
            }
            break;
        case ACTIVELEARNING_SUPERVISING_STATE:
            if(activeLearning->supervisingActiveLearningFunction == NULL){
                iftError("supervise samples function is not defined","iftCommonActiveLearningGoToNextState");
                activeLearning->nextState =  ACTIVELEARNING_ENDING;
            }else{
                activeLearning->supervisingActiveLearningFunction(activeLearning);
            }
            break;
        case ACTIVELEARNING_CLASSIFIER_TRAIN_STATE:
            if(activeLearning->trainActiveLearningFunction == NULL){
                iftError("train classifier function is not defined","iftCommonActiveLearningGoToNextState");
                activeLearning->nextState =  ACTIVELEARNING_ENDING;
            }else{
                activeLearning->trainActiveLearningFunction(activeLearning);
            }
            break;
        case ACTIVELEARNING_POSPROCESSING_STATE:
            if(activeLearning->posProcessingActiveLearningFunction == NULL){
                iftWarning("pos process function is not defined","iftCommonActiveLearningGoToNextState");
                activeLearning->nextState =  ACTIVELEARNING_SELECT_SAMPLES_STATE;
            }else{
                activeLearning->posProcessingActiveLearningFunction(activeLearning);
            }
            break;
        case ACTIVELEARNING_ENDING:
            if(activeLearning->endActiveLearningFunction == NULL){
                iftWarning("end function is not defined","iftCommonActiveLearningGoToNextState");
            }else{
                activeLearning->endActiveLearningFunction(activeLearning);
            }
            activeLearning->activeLearningFinished = true;
            break;
        default:
            iftError("unkwon state","iftCommonActiveLearningGoToNextState");
            activeLearning->nextState = ACTIVELEARNING_ENDING;
            break;
    }

}
///////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////MACHINE LEARNING PRE_PROCESSING////////////

////PCA
iftPrincipalComponentAnalysis* iftCreatePrincipalComponentAnalysis(iftMatrix* matrix){
    iftPrincipalComponentAnalysis* pca = iftAlloc(1,sizeof(iftPrincipalComponentAnalysis));
    pca->inputData = matrix;
    pca->mean = NULL;
    pca->rotationMatrix = NULL;
    pca->rotationMatrix_original = NULL;
    pca->rotationMatrix = NULL;
    pca->variance = NULL;
    pca->variance_frac = NULL;
    pca->outputData = NULL;
    pca->numberComponents = 0.99;
    pca->destroyInputData = false;
    pca->destroyOutData = false;
    pca->inputDataset = NULL;
    pca->outputDataset = NULL;
    return pca;
}

iftPrincipalComponentAnalysis* iftCreatePrincipalComponentAnalysisGivenDataset(iftDataSet* dataSet){
    iftPrincipalComponentAnalysis* pca =  iftCreatePrincipalComponentAnalysis(dataSet->data);
    pca->inputDataset = dataSet;
    return pca;
}

void iftDestroyPrincipalComponentAnalysis(iftPrincipalComponentAnalysis** pca_pointer){
    if(pca_pointer==NULL || *pca_pointer == NULL){
        return;
    }
    iftPrincipalComponentAnalysis* aux = *pca_pointer;
    if(aux->destroyInputData == true){
        if(aux->inputData != NULL){
            iftDestroyMatrix(&(aux->inputData));
        }
    }
    if(aux->mean != NULL){
        iftDestroyMatrix(&(aux->mean));
    }
    if(aux->rotationMatrix != NULL){
        iftDestroyMatrix(&(aux->rotationMatrix));
    }
    if(aux->variance != NULL){
        iftDestroyMatrix(&(aux->variance));
    }
    if(aux->variance_frac != NULL){
        iftDestroyMatrix(&(aux->variance_frac));
    }
    if(aux->rotationMatrix_original != NULL){
        iftDestroyMatrix(&(aux->rotationMatrix_original));
    }
    if(aux->inverseRotationMatrix_original != NULL){
        iftDestroyMatrix(&(aux->inverseRotationMatrix_original));
    }
    if(aux->destroyOutData == true){
        if(aux->outputData != NULL){
            iftDestroyMatrix(&(aux->outputData));
        }
    }
    iftFree(*pca_pointer);
    *pca_pointer = NULL;
}

void iftComputePrincipalComponentAnalysis(iftPrincipalComponentAnalysis *pca){
    if(pca == NULL){
        iftError("PCA object is null","iftComputePrincipalComponentAnalysis");
    }
    if(pca->inputData == NULL){
        iftError("inputData is null","iftComputePrincipalComponentAnalysis");
    }
    if(pca->inputData->nrows <= 0 || pca->inputData->ncols <= 0){
        iftError("inputData has invalid shape","iftComputePrincipalComponentAnalysis");
    }
    iftMatrix* inputMatrix = pca->inputData;
    iftMatrix* Y = iftCreateMatrix(inputMatrix->ncols,inputMatrix->nrows);
    int numberDimensions = inputMatrix->ncols;

    //data centralization
    double *mean=iftAlloc(numberDimensions,sizeof(double));
    if(pca->mean != NULL){
        iftDestroyMatrix(&(pca->mean));
    }
    pca->mean = iftCreateMatrix(numberDimensions,1);

    int k=0;
    for (int row=0; row < inputMatrix->nrows; row++) {
        for (int col = 0; col < inputMatrix->ncols; col++) {
            //truste me...this works
            Y->val[k] = inputMatrix->val[k];
            mean[col] += inputMatrix->val[k];
            k++;
        }
    }

    for (int col = 0; col < inputMatrix->ncols; col++) {
        mean[col] /= inputMatrix->nrows;
        pca->mean->val[col] = mean[col];
    }

    k=0;
    for (int row=0; row < inputMatrix->nrows; row++) {
        for (int col = 0; col < inputMatrix->ncols; col++) {
            Y->val[k] = (Y->val[k] - mean[col]);
            k++;
        }
    }

    float factor = 1.0 / sqrt(inputMatrix->ncols-1);
    iftMultMatrixByScalar(Y,factor);

    iftMatrix *eigenVectors_U,*eigenValues,*eigenVectors_Vt,* eigenVectors_V, *variance,*variance_frac;
    iftSingleValueDecomp(Y, &eigenVectors_U, &eigenValues, &eigenVectors_Vt);
    eigenVectors_V = iftTransposeMatrix(eigenVectors_Vt);

    variance = iftCreateMatrix(eigenValues->ncols,eigenValues->nrows);
    variance_frac = iftCreateMatrix(eigenValues->ncols,eigenValues->nrows);
    double sum = 0;
    for (int i = 0; i < eigenValues->n; ++i) {
        variance->val[i] = eigenValues->val[i]*eigenValues->val[i];
        sum += variance->val[i];
    }
    for (int i = 0; i < eigenValues->n; ++i) {
        variance_frac->val[i] = variance->val[i]/sum;
    }

    int ncols_reduced = 0;
    if(pca->numberComponents <= 1.0) {
        sum = 0.0;
        for (int i = 0; i < variance_frac->n; ++i) {
            ncols_reduced++;
            sum += variance_frac->val[i];
            if(sum>=pca->numberComponents) {
                break;
            }
        }
    }
    else {
        ncols_reduced = iftMin(pca->numberComponents, variance_frac->n);
    }
    iftMatrix* pca_rotation = iftSubMatrix(eigenVectors_V, 0, eigenVectors_V->nrows, 0, ncols_reduced);
    pca->rotationMatrix = pca_rotation;
    pca->rotationMatrix_original = eigenVectors_V;
    pca->inverseRotationMatrix_original = eigenVectors_Vt;
    pca->variance = variance;
    pca->variance_frac = variance_frac;

    if(pca->inputDataset == NULL){
        iftApplyPricipalComponentAnalysis(pca,inputMatrix);
    }else{
        iftApplyPricipalComponentAnalysisGivenDataset(pca,pca->inputDataset);
    }

    iftFree(mean);
    iftDestroyMatrix(&eigenVectors_U);
    iftDestroyMatrix(&eigenValues);
    iftDestroyMatrix(&Y);
}


void iftApplyPricipalComponentAnalysis(iftPrincipalComponentAnalysis* pca,iftMatrix* matrix){
    if(pca == NULL || matrix == NULL){
        iftError("Input object(s) is/are null","iftApplyPricipalComponentAnalysis");
    }
    if(matrix->ncols != pca->mean->ncols){
        iftError("Dimesion mismatch","iftApplyPricipalComponentAnalysis");
    }
    int k =0;
    for (int row=0; row < matrix->nrows; row++) {
        for (int col = 0; col < matrix->ncols; col++) {
            matrix->val[k] = (matrix->val[k] - pca->mean->val[col]);
            k++;
        }
    }
    pca->outputData = iftMultMatrices(matrix,pca->rotationMatrix);
    if(pca->inputDataset != NULL){
        iftDataSet* inputDataset = pca->inputDataset;
        iftDataSet* outputDataset = iftCreateDataSetAux(pca->outputData->nrows,pca->outputData->ncols,false);
        outputDataset->data = pca->outputData;
        outputDataset->nclasses = inputDataset->nclasses;
        outputDataset->nsamples = inputDataset->nsamples;
        outputDataset->ntrainsamples = inputDataset->ntrainsamples;
        outputDataset->ngroups = inputDataset->ngroups;
        outputDataset->projection = inputDataset->projection;
        for (int i = 0; i < outputDataset->nsamples; ++i) {
            outputDataset->sample[i].truelabel = inputDataset->sample[i].truelabel;
            outputDataset->sample[i].label = inputDataset->sample[i].label;
            outputDataset->sample[i].isSupervised = inputDataset->sample[i].isSupervised;
            iftSetSampleStatus(&outputDataset->sample[i], inputDataset->sample[i].status);
            outputDataset->sample[i].group = inputDataset->sample[i].group;
            outputDataset->sample[i].id = inputDataset->sample[i].id;
            outputDataset->sample[i].weight = inputDataset->sample[i].weight;
            outputDataset->sample[i].feat = &(iftMatrixElem(pca->outputData,0,i));
        }
        pca->outputDataset = outputDataset;
    }
}

void iftApplyPricipalComponentAnalysisGivenDataset(iftPrincipalComponentAnalysis* pca,iftDataSet* dataSet){
    iftMatrix* matrix = dataSet->data;
    if(pca == NULL || matrix == NULL){
        iftError("Input object(s) is/are null","iftApplyPricipalComponentAnalysis");
    }
    if(matrix->ncols != pca->mean->ncols){
        iftError("Dimesion mismatch","iftApplyPricipalComponentAnalysis");
    }
    int k =0;
    for (int row=0; row < matrix->nrows; row++) {
        for (int col = 0; col < matrix->ncols; col++) {
            matrix->val[k] = (matrix->val[k] - pca->mean->val[col]);
            k++;
        }
    }
    pca->outputData = iftMultMatrices(matrix,pca->rotationMatrix);
    iftDataSet* inputDataset = pca->inputDataset;
    iftDataSet* outputDataset = iftCreateDataSetAux(pca->outputData->nrows,pca->outputData->ncols,false);
    outputDataset->data = pca->outputData;
    outputDataset->nclasses = inputDataset->nclasses;
    outputDataset->nsamples = inputDataset->nsamples;
    outputDataset->ntrainsamples = inputDataset->ntrainsamples;
    outputDataset->ngroups = inputDataset->ngroups;
    outputDataset->projection = inputDataset->projection;
    iftCopyRefData(outputDataset, inputDataset->ref_data, inputDataset->ref_data_type);
    // outputDataset->ref_data_type = inputDataset->ref_data_type;
    // outputDataset->ref_data = inputDataset->ref_data;//be careful
    for (int i = 0; i < outputDataset->nsamples; ++i) {
        outputDataset->sample[i].feat = &(iftMatrixElem(pca->outputData,0,i));
        outputDataset->sample[i].truelabel = dataSet->sample[i].truelabel;
        outputDataset->sample[i].label = dataSet->sample[i].label;
        outputDataset->sample[i].isSupervised = dataSet->sample[i].isSupervised;
        iftSetSampleStatus(&outputDataset->sample[i], dataSet->sample[i].status);
        outputDataset->sample[i].group = dataSet->sample[i].group;
        outputDataset->sample[i].id = dataSet->sample[i].id;
        outputDataset->sample[i].weight = dataSet->sample[i].weight;
    }
    pca->outputDataset = outputDataset;
}

iftDataSet* iftComputePrincipalComponentAnalysisGivenDataset(iftDataSet *dataSet, double nComponents){
    iftPrincipalComponentAnalysis* pca = iftCreatePrincipalComponentAnalysisGivenDataset(dataSet);
    pca->destroyOutData = false;
    pca->destroyInputData = false;
    pca->numberComponents = nComponents;
    iftComputePrincipalComponentAnalysis(pca);
    iftDataSet* output = pca->outputDataset;
    iftDestroyPrincipalComponentAnalysis(&pca);
    return output;
}





//////////////////////////////////////////////////////////

