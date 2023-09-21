#include "iftDataSet.h"
#include "iftIterativeOPF.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include <iftRegion.h>
#include <iftMetrics.h>
#include "iftFiltering.h"
#include "iftGenericVectorUtil.h"
#include "iftClustering.h"


/**
@file
@brief A description is missing here
*/
iftDistanceTable   *iftDist               = NULL;
iftMinkowskiTable  *iftMinkowski          = NULL;

/********************** PRIVATE FUNCTIONS *************************/
void iftIsThereTrueLabelOrGroupZero(const iftDataSet *Z, const char *function_name) {
    for (int s = 0; s < Z->nsamples; s++)
        if (Z->sample[s].truelabel == 0) {
            iftWarning("DataSet with some sample(s) with true label with Zero value", function_name);
            break;
        }

    for (int s = 0; s < Z->nsamples; s++)
        if (Z->sample[s].group == 0) {
            iftWarning("DataSet with some sample(s) with group with Zero value", function_name);
            break;
        }
}


void iftSetDataSetIDsincrementally(iftDataSet *Z) {
    for (int s = 0; s < Z->nsamples; s++)
        Z->sample[s].id = s;
}


/**
* @brief General function to destroy a DataSet
* @param Z [description]
*/
void iftGeneralDestroyDataSet(iftDataSet **Z, bool destroy_feats) {
    if (Z != NULL) {
        iftDataSet *aux = *Z;

        if(aux != NULL) {
            if (aux->fsp.mean != NULL) {
                iftFree(aux->fsp.mean);
            }
            if (aux->fsp.stdev != NULL) {
                iftFree(aux->fsp.stdev);
            }
            if (aux->fsp.w != NULL) {
                iftFree(aux->fsp.w);
            }
            if (aux->fsp.R != NULL) {
                iftDestroyMatrix(&(aux->fsp.R));
            }
            if (aux->fsp.W != NULL) {
                iftDestroyMatrix(&(aux->fsp.W));
            }
            if (aux->alpha != NULL)
                iftFree(aux->alpha);

            if (aux->sample != NULL) {
                if (destroy_feats) {
                    for (int s = 0; s < aux->nsamples; s++) {
                        if (aux->sample[s].feat != NULL)
                            iftFree(aux->sample[s].feat);
                    }
                }
                iftFree(aux->sample);
            }
            if(aux->projection != NULL){
                iftDestroyDoubleMatrix(&(aux->projection));
            }

            iftFree(aux);
            *Z = NULL;
        }
    }
}

/**
* @brief Get the index of the first Non-NULL DataSet from an array.
* @author Samuel Martins
* @dataset Jun 6, 2016
*
* @param  Z_array    Array of DataSets which could have NULL datasets
* @param  n_datasets Number of Datasets (array size)
* @return            The index of the first Non-NULL DataSet from an array or -1 if all datasets are NULL.
*/
int iftGetFirstDataSetIndexFromArray(iftDataSet **Z_array, int n_datasets) {
    int i0 = -1;

    // since it could have NULL datasets, let's find the first non-NULL dataset
    for (int i = 0; i < n_datasets; i++) {
        if (Z_array[i] != NULL) {
            i0 = i;
            break;
        }
    }

    return i0;
}


int iftGetTotalNumberOfSamplesFromDataSets(iftDataSet **Z_array, int n_datasets) {
    int nt_samples = 0;
    int i0 = iftGetFirstDataSetIndexFromArray(Z_array, n_datasets);
    int n_feats    = Z_array[i0]->nfeats;

    // gets the total number of samples and checks if all datasets have the same number of feats
    for (int i = 0; i < n_datasets; i++) {
        if (Z_array[i] != NULL) {
            nt_samples += Z_array[i]->nsamples;

            if (Z_array[i]->nfeats != n_feats) {
                iftError("Datasets with different number of feats: %d, %d", "iftGetTotalNumberOfSamplesFromDataSets",
                         n_feats, Z_array[i]->nfeats);
            }
        }
    }

    return nt_samples;
}

iftSampler* iftCreateSampler(int nsamples, int niters) {
    iftSampler* sampler = (iftSampler*) iftAlloc(1, sizeof(iftSampler));
    sampler->niters     = niters;
    sampler->nsamples   = nsamples;

    sampler->status = (char**) iftAlloc(niters, sizeof(char*));
    for (int i = 0; i < niters; ++i) {
        sampler->status[i] = iftAllocCharArray(nsamples);
    }

    return sampler;
}

void iftDestroySampler(iftSampler** sampler) {
    iftSampler* s = *sampler;

    if (s != NULL) {
        for (int i = 0; i < s->niters; ++i) {
            iftFree(s->status[i]);
        }
        iftFree(s->status);
        iftFree(s);
        *sampler = NULL;
    }
}

void iftDataSetSampling(iftDataSet *Z, const iftSampler *sampler, int iteration) {

    if(iteration<0 || iteration>=sampler->niters) {
        iftError("Invalid iteration.", "iftDataSetSampling");
    }

    if(sampler->nsamples!=Z->nsamples) {
        iftError("Invalid number of samples.", "iftDataSetSampling");
    }

    Z->ntrainsamples = 0;

    for (int i = 0; i < Z->nsamples; ++i) {
      if (sampler->status[iteration][i]==IFT_TRAIN){
	iftRemoveSampleStatus(&Z->sample[i],IFT_TEST);      
	iftAddSampleStatus(&Z->sample[i],IFT_TRAIN);      
	Z->ntrainsamples++;
      } else { /* It has status of test sample */
	iftRemoveSampleStatus(&Z->sample[i],IFT_TRAIN);      
	iftAddSampleStatus(&Z->sample[i],IFT_TEST);            
      }
    }
}

iftSampler* iftKFold(int nsamples, int k) {
    int *folds = iftAllocIntArray(nsamples);

    for (int i = 0; i < nsamples; ++i) {
        folds[i] = i%k;
    }

    iftShuffleIntArray(folds, nsamples);

    iftSampler* sampler = iftCreateSampler(nsamples, k);

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < nsamples; ++j) {
            sampler->status[i][j] = folds[j]==i? IFT_TEST : IFT_TRAIN;
        }
    }

    iftFree(folds);

    return sampler;
}

iftSampler* iftInverseKFold(int nsamples, int k) {
    int *folds = iftAllocIntArray(nsamples);

    for (int i = 0; i < nsamples; ++i) {
        folds[i] = i%k;
    }

    iftShuffleIntArray(folds, nsamples);

    iftSampler* sampler = iftCreateSampler(nsamples, k);

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < nsamples; ++j) {
            sampler->status[i][j] = folds[j]==i? IFT_TRAIN : IFT_TEST;
        }
    }

    iftFree(folds);

    return sampler;
}

iftSampler* iftNKFold(int nsamples, int n, int k) {
    int *folds = iftAllocIntArray(nsamples);

    for (int i = 0; i < nsamples; ++i) {
        folds[i] = i%k;
    }

    iftSampler* sampler = iftCreateSampler(nsamples, n*k);

    int it = 0;
    for (int l = 0; l < n; ++l) {//repeat kfold n times

        iftShuffleIntArray(folds, nsamples);

        for (int i = 0; i < k; ++i) {//k fold
            for (int j = 0; j < nsamples; ++j) {
                sampler->status[it][j] = folds[j]==i? IFT_TEST : IFT_TRAIN;
            }
            it++;
        }
    }

    iftFree(folds);

    return sampler;
}

iftSampler* iftLeaveOneOut(int nsamples) {

    iftSampler* sampler = iftCreateSampler(nsamples, nsamples);

    for (int i = 0; i < nsamples; ++i) {
        memset(sampler->status, IFT_TRAIN, nsamples*sizeof(char));
    }

    for (int i = 0; i < nsamples; ++i) {
        sampler->status[i][i] = IFT_TEST;
    }

    return sampler;
}

iftSampler* iftRandomSubsampling(int nsamples, int n, int ntrain) {

    if(ntrain>nsamples) {
        iftError("Invalid number of training samples.", "iftRandomSubsampling");
    }

    iftSampler* sampler = iftCreateSampler(nsamples, n);

    int* status = iftAllocIntArray(nsamples);

    for (int i = 0; i < ntrain; ++i)
        status[i] = IFT_TRAIN;
    for (int i = ntrain; i < nsamples; ++i)
        status[i] = IFT_TEST;

    for (int i = 0; i < n; ++i) {
        iftShuffleIntArray(status, nsamples);
        for (int j = 0; j < nsamples; ++j) {
            sampler->status[i][j] = status[j];
        }
    }

    iftFree(status);

    return sampler;
}


iftSampler *iftStratifiedRandomSubsampling(const int *labels, int n_samples, int n_iters, int n_train) {
    if (labels == NULL)
        iftError("Label Array is NULL", "iftStratifiedRandomSubsampling");
    if (n_train > n_samples) {
        iftError("Invalid number of training samples to be sampled/selected: %d >= %d (number of samples).",
                 "iftStratifiedRandomSubsampling", n_train, n_samples);
    }

    ///////////
    // since that available labels not necessarily goes from 1 to n, we have to figure out the available labels (BUG this decision was a big mistake --- someone has to fix this and avoid samples with label 0 here)
    
    iftIntArray *classes = iftIntArrayUnique(labels, n_samples);
    int n_classes        = classes->n;

    // stores the pair (class, class idx)
    // ex: classes = [2, 5, 255] --> the class_idxs_dict for the class 255 is 2
    iftDict *class_idxs_dict = iftCreateDictWithApproxSize(n_classes);
    for (int cidx = 0; cidx < n_classes; cidx++) {
        int class_label = classes->val[cidx];
        iftInsertIntoDict(class_label, cidx, class_idxs_dict);
    }
    iftDestroyIntArray(&classes);
    ///////////

    ///////////
    // computes the number of samples of each class
    int *class_sizes = iftAllocIntArray(n_classes);
    for (int s = 0; s < n_samples; s++) {
        int class_label = labels[s];
        int cidx        = iftGetLongValFromDict(class_label, class_idxs_dict); // index of the class <class_label>
        class_sizes[cidx]++;
    }
    ///////////

    // Counts the number of train samples per class
    int *n_train_per_class = iftAllocIntArray(n_classes);
    int n_total_train      = 0;

    for (int cidx = 0; cidx < n_classes; cidx++) {
        n_train_per_class[cidx] = n_train * (class_sizes[cidx] / (1.0 * n_samples));
        n_total_train          += n_train_per_class[cidx];
    }

    // TODO: the better way is to add the remainder samples from the class with lowest size to the highest one
    // make sure we get <n_train> samples, it might sacrifice a perfect preservation of a priori proportion
    for (int cidx = 0; cidx < n_classes; cidx = (cidx+1) % n_classes) {
        if (n_total_train >= n_train)
            break;
        else {
            n_train_per_class[cidx]++;
            n_total_train++;
        }
    }
    ///////////

    ///////////
    // Stores the indices of the samples from each class
    int **class_sample_idxs = (int**) iftAlloc(n_classes, sizeof(int*));
    for (int cidx = 0; cidx < n_classes; cidx++) {
        class_sample_idxs[cidx] = iftAllocIntArray(class_sizes[cidx]); // the size of the queue is the number of samples of the class
    }

    int *counts = iftAllocIntArray(n_classes);
    for (int s = 0; s < n_samples; s++) {
        int class_label = labels[s];
        int cidx        = iftGetLongValFromDict(class_label, class_idxs_dict); // index of the class <class_label>
        class_sample_idxs[cidx][counts[cidx]++] = s;
    }
    ///////////

    iftSampler *sampler = iftCreateSampler(n_samples, n_iters);

    // For each iteration, sampling the training samples per class by shuffling their indices
    for (int it = 0; it < n_iters; it++) {
        for (int cidx = 0; cidx < n_classes; cidx++) {
            iftShuffleIntArray(class_sample_idxs[cidx], class_sizes[cidx]);

            // assigns the sampled/selected training samples
            for (int i = 0; i < n_train_per_class[cidx]; i++) {
                int s                  = class_sample_idxs[cidx][i]; // sample index
                sampler->status[it][s] = IFT_TRAIN;
            }

            // assigns the remainder samples as testing samples
            for (int i = n_train_per_class[cidx]; i < class_sizes[cidx]; i++) {
                int s                  = class_sample_idxs[cidx][i]; // sample index
                sampler->status[it][s] = IFT_TEST;
            }
        }
    }
    ///////////

    // DESTROYERS
    iftDestroyDict(&class_idxs_dict);
    iftFree(class_sizes);
    iftFree(n_train_per_class);
    for (int cidx = 0; cidx < n_classes; cidx++)
        iftFree(class_sample_idxs[cidx]);
    iftFree(class_sample_idxs);

    return sampler;
}


iftSampler *iftBalancedRandomSubsampling(const int *labels, int n_samples, int n_iters, int n_train) {
    if (labels == NULL)
        iftError("Label Array is NULL", "iftBalancedRandomSubsampling");
    if (n_train > n_samples) {
        iftError("Invalid number of training samples to be sampled/selected: %d >= %d (number of samples).",
                 "iftBalancedRandomSubsampling", n_train, n_samples);
    }

    ///////////
    // since that available labels not necessarily goes from 1 to n, we have to figure out the available labels
    iftIntArray *classes = iftIntArrayUnique(labels, n_samples);
    int n_classes        = classes->n;

    // stores the pair (class, class idx)
    // ex: classes = [2, 5, 255] --> the class_idxs_dict for the class 255 is 2
    iftDict *class_idxs_dict = iftCreateDictWithApproxSize(n_classes);
    for (int cidx = 0; cidx < n_classes; cidx++) {
        int class_label = classes->val[cidx];
        iftInsertIntoDict(class_label, cidx, class_idxs_dict);
    }
    iftDestroyIntArray(&classes);
    ///////////

    ///////////
    // computes the number of samples of each class
    int *class_sizes = iftAllocIntArray(n_classes);
    for (int s = 0; s < n_samples; s++) {
        int class_label = labels[s];
        int cidx        = iftGetLongValFromDict(class_label, class_idxs_dict); // index of the class <class_label>
        class_sizes[cidx]++;
    }
    ///////////

    // gets the min number of samples from the classes
    int min_n_classes = IFT_INFINITY_INT;
    for (int cidx = 0; cidx < n_classes; cidx++) {
        if (class_sizes[cidx] < min_n_classes) {
            min_n_classes = class_sizes[cidx];
        }
    }

    // checks the number of samples after balancing
    int n_total_balanced_samples = n_classes * min_n_classes;
    if (n_total_balanced_samples < n_train) {
        iftWarning("Number of Samples after balancing is < than the required num of train. samples: %d < %d\n" \
                 "Then, the number of train. samples will be %d\n", "iftBalancedRandomSubsampling",
                   n_total_balanced_samples, n_train, n_total_balanced_samples);
        n_train = n_total_balanced_samples;
    }

    // Counts the number of train samples per class
    int *n_train_per_class           = iftAllocIntArray(n_classes);
    int n_balanced_samples_per_class = n_train / n_classes;
    int n_total_train                = 0;

    for (int cidx = 0; cidx < n_classes; cidx++) {
        n_train_per_class[cidx] = n_balanced_samples_per_class;
        n_total_train          += n_train_per_class[cidx];
    }

    // make sure we get <n_train> samples, it might sacrifice a perfect class balance
    for (int cidx = 0; cidx < n_classes; cidx = (cidx+1) % n_classes) {
        if (n_total_train >= n_train)
            break;
        else {
            n_train_per_class[cidx]++;
            n_total_train++;
        }
    }
    ///////////

    ///////////
    // Stores the indices of the samples from each class
    // We used a queue just to indices the elements more easily
    int **class_sample_idxs = (int**) iftAlloc(n_classes, sizeof(int*));
    for (int cidx = 0; cidx < n_classes; cidx++) {
        class_sample_idxs[cidx] = iftAllocIntArray(class_sizes[cidx]); // the size of the queue is the number of samples of the class
    }

    int *counts = iftAllocIntArray(n_classes);
    for (int s = 0; s < n_samples; s++) {
        int class_label = labels[s];
        int cidx        = iftGetLongValFromDict(class_label, class_idxs_dict); // index of the class <class_label>
        class_sample_idxs[cidx][counts[cidx]++] = s;
    }
    ///////////

    iftSampler *sampler = iftCreateSampler(n_samples, n_iters);

    // For each iteration, sampling the training samples per class by shuffling their indices
    for (int it = 0; it < n_iters; it++) {
        for (int cidx = 0; cidx < n_classes; cidx++) {
            iftShuffleIntArray(class_sample_idxs[cidx], class_sizes[cidx]);

            // assigns the sampled/selected training samples
            for (int i = 0; i < n_train_per_class[cidx]; i++) {
                int s                  = class_sample_idxs[cidx][i]; // sample index
                sampler->status[it][s] = IFT_TRAIN;
            }

            // assigns the remainder samples as testing samples
            for (int i = n_train_per_class[cidx]; i < class_sizes[cidx]; i++) {
                int s                  = class_sample_idxs[cidx][i]; // sample index
                sampler->status[it][s] = IFT_TEST;
            }
        }
    }
    ///////////

    // DESTROYERS
    iftDestroyDict(&class_idxs_dict);
    iftFree(class_sizes);
    iftFree(n_train_per_class);
    for (int cidx = 0; cidx < n_classes; cidx++)
        iftFree(class_sample_idxs[cidx]);
    iftFree(class_sample_idxs);

    return sampler;
}


bool iftIsSupervisedDataSet(const iftDataSet *Z)
{
    if (Z->nclasses > 1)
        return(1);
    else
        return(0);
}

bool iftIsTrainingDataSet(const iftDataSet *Z)
{
    if (Z->ntrainsamples == Z->nsamples)
        return(1);
    else
        return(0);
}

bool iftIsTestingDataSet(const iftDataSet *Z)
{
    if (Z->ntrainsamples == 0)
        return(1);
    else
        return(0);
}

bool iftIsNormalizedDataSet(const iftDataSet *Z)
{

    if ((Z->fsp.mean != NULL)&&(Z->fsp.stdev != NULL)){
        return(1);
    }else{
        return(0);
    }

}

bool iftIsCentralizedDataSet(const iftDataSet *Z)
{
    if (Z->fsp.mean != NULL)
        return(1);
    else
        return(0);
}

char        iftIsTransformedDataSetByPCA(iftDataSet *Z)
{
    if (Z->fsp.R != NULL)
        return(1);
    else
        return(0);
}

bool iftIsWhitenedDataSet(const iftDataSet *Z)
{
    if (Z->fsp.W != NULL)
        return(1);
    else
        return(0);
}


iftDataSet* iftCreateDataSetAux(int nsamples, int nfeats, bool allocData){
    if (nsamples < 0)
        iftError("Invalid Number of Samples: %d < 0\n", "iftCreateDataSet", nsamples);
    if (nfeats < 0)
        iftError("Invalid Number of Feats: %d < 0\n", "iftCreateDataSet", nfeats);


    iftDataSet *Z = (iftDataSet*) iftAlloc(1, sizeof(iftDataSet));
    Z->ntrainsamples = 0;
    Z->capacity      = nsamples;
    Z->nsamples      = nsamples;
    Z->nfeats        = nfeats;
    Z->ngroups       = 0;
    Z->nclasses      = 0;
    Z->sample = (iftSample*) iftAlloc(Z->capacity, sizeof(iftSample));
    
    if (Z->sample == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDataSet");
    }

    iftSetDistanceFunction(Z, 1);
    Z->alpha = iftAllocFloatArray(nfeats);
    for (int f = 0; f < nfeats; f++) {
        Z->alpha[f] = 1.0;
    }

    Z->fsp.mean   = NULL;
    Z->fsp.stdev  = NULL;
    Z->fsp.R      = NULL;
    Z->fsp.W      = NULL;
    Z->fsp.w      = NULL;
    Z->fsp.nfeats = 0;
    Z->fsp.ncomps = 0;

    Z->ref_data_type = IFT_REF_DATA_NULL;
    Z->ref_data      = NULL;
    Z->projection =    NULL;

    if(allocData){
        Z->data = iftCreateMatrix(Z->nfeats, Z->capacity);
        for (int s = 0; s < nsamples; s++) {
            Z->sample[s].id        = IFT_NIL;
            Z->sample[s].truelabel = 0;
            Z->sample[s].label     = 0;
            Z->sample[s].group     = 0;
            Z->sample[s].feat = &(iftMatrixElem(Z->data, 0, s));
            iftSetSampleStatus(&Z->sample[s], IFT_UNKNOWN);
            Z->sample[s].weight    = 0.0;
        }
    }else{
        Z->data = NULL;
        for (int s = 0; s < nsamples; s++) {
            Z->sample[s].id        = IFT_NIL;
            Z->sample[s].truelabel = 0;
            Z->sample[s].label     = 0;
            Z->sample[s].group     = 0;
            Z->sample[s].feat = NULL;
            iftSetSampleStatus(&Z->sample[s], IFT_UNKNOWN);
            Z->sample[s].weight    = 0.0;
        }
    }

    return Z;
}


void iftGrowDataset(iftDataSet* dataset, int capacity) {
    if(dataset->capacity >= capacity) {
        return;
    }

    dataset->sample = (iftSample*) iftRealloc(dataset->sample, capacity*sizeof(iftSample));
    dataset->data->val = (float*) iftRealloc(dataset->data->val, capacity*dataset->nfeats*sizeof(float));
    dataset->capacity = capacity;
}

void iftAddSampleDataset(iftDataSet* dataset, iftSample sample) {

    if(dataset->nsamples >= dataset->capacity) {
        iftGrowDataset(dataset, dataset->capacity * 2);
    }

    bool copyFeats = sample.feat!=NULL;

    dataset->sample[dataset->nsamples].feat = &iftMatrixElem(dataset->data, 0, dataset->nsamples);
    iftCopySample(&sample, &(dataset->sample[dataset->nsamples]), dataset->nfeats, copyFeats);

    if(!copyFeats)
        dataset->sample[dataset->nsamples].feat = &iftMatrixElem(dataset->data, 0, dataset->nsamples);

    dataset->nsamples++;
    dataset->data->nrows++;
}

iftDataSet *iftCreateDataSet(int nsamples, int nfeats) {
    return iftCreateDataSetAux(nsamples,nfeats,true);
}


/*
* Create an iftDataSet with no Feature allocation - it used for address copy between iftDataSets
*/
iftDataSet *iftCreateDataSetWithoutFeatsVector(int nsamples, int nfeats) {
    if (nsamples < 0)
        iftError("Invalid Number of Samples: %d < 0\n", "iftCreateDataSetWithoutFeatsVector", nsamples);
    if (nfeats < 0)
        iftError("Invalid Number of Feats: %d < 0\n", "iftCreateDataSetWithoutFeatsVector", nfeats);

    iftDataSet *Z = (iftDataSet*) iftAlloc(1, sizeof(iftDataSet));
    Z->ntrainsamples = 0;
    Z->nsamples      = nsamples;
    Z->nfeats        = nfeats;
    Z->ngroups       = 0;
    Z->nclasses      = 0;
    iftSetDistanceFunction(Z, 1);

    Z->sample = (iftSample*) iftAlloc(nsamples, sizeof(iftSample));
    if (Z->sample == NULL)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDataSetWithoutFeatsVector");

    for (int s = 0; s < nsamples; s++) {
        Z->sample[s].id        = IFT_NIL;
        Z->sample[s].truelabel = 0;
        Z->sample[s].label     = 0;
        // Z->sample[s].feat      = iftAllocFloatArray(nfeats);
        iftSetSampleStatus(&Z->sample[s], IFT_UNKNOWN);
        Z->sample[s].weight    = 0.0;
    }

    Z->alpha = iftAllocFloatArray(nfeats);
    for (int f = 0; f < nfeats; f++)
        Z->alpha[f] = 1.0;

    Z->fsp.mean   = NULL;
    Z->fsp.stdev  = NULL;
    Z->fsp.R      = NULL;
    Z->fsp.W      = NULL;
    Z->fsp.w      = NULL;
    Z->fsp.nfeats = 0;
    Z->fsp.ncomps = 0;

    Z->ref_data_type = IFT_REF_DATA_NULL;
    Z->ref_data      = NULL;

    return Z;
}


iftDataSet **iftCreateDataSetArray(size_t n_datasets) {
    return ((iftDataSet**) iftAlloc(n_datasets, sizeof(iftDataSet*)));
}


void iftDestroyRefData(iftDataSet *Z) {
    if (Z->ref_data != NULL) {
        if (Z->ref_data_type == IFT_REF_DATA_IMAGE)
            iftDestroyImage(((iftImage **) &Z->ref_data));
        else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE)
            iftDestroyFImage(((iftFImage **) &Z->ref_data));
        else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE)
            iftDestroyMImage(((iftMImage **) &Z->ref_data));
        else if (Z->ref_data_type == IFT_REF_DATA_FILESET)
            iftDestroyFileSet(((iftFileSet **) &Z->ref_data));
        else if (Z->ref_data_type == IFT_REF_DATA_MMKERNEL)
            iftDestroyMMKernel(((iftMMKernel **) &Z->ref_data));
        else if (Z->ref_data_type == IFT_REF_DATA_CSV)
            iftDestroyCSV(((iftCSV **) &Z->ref_data));
        else iftError("Destroyer for the Ref. Data Type not implemented yet", "iftDestroyRefData");

        Z->ref_data = NULL;
    }
}


void iftDestroyDataSet(iftDataSet **Z) {

    if (Z != NULL) {
        iftDataSet *aux = *Z;

        if(aux != NULL) {
            if (aux->fsp.mean != NULL) {
                iftFree(aux->fsp.mean);
            }

            if (aux->fsp.stdev != NULL) {
                iftFree(aux->fsp.stdev);
            }

            if (aux->fsp.w != NULL) {
                iftFree(aux->fsp.w);
            }

            if (aux->fsp.R != NULL) {
                iftDestroyMatrix(&(aux->fsp.R));
            }

            if (aux->fsp.W != NULL) {
                iftDestroyMatrix(&(aux->fsp.W));
            }

            if (aux->alpha != NULL)
                iftFree(aux->alpha);

            if (aux->sample != NULL) {
                iftFree(aux->sample);
            }
	    
            if(aux->projection != NULL){
	        	iftDestroyDoubleMatrix(&(aux->projection));
            }

            if(aux->data != NULL){
	        	iftDestroyMatrix(& (aux->data));
            }

            iftDestroyRefData(aux);

            iftFree(aux);
            *Z = NULL;
        }
    }
}


void iftDestroyDataSetWithoutFeatVectors(iftDataSet **Z) {
    iftGeneralDestroyDataSet(Z, false);
}


void iftDestroyDataSetArray(iftDataSet ***Z_array, size_t n_datasets, bool destroy_feats) {
    if (Z_array != NULL) {
        iftDataSet **aux = *Z_array;

        if (destroy_feats) {
            for (size_t i = 0; i < n_datasets; i++)
                iftDestroyDataSet(&aux[i]);
        }
        else {
            for (size_t i = 0; i < n_datasets; i++){
                iftDestroyDataSetWithoutFeatVectors(&aux[i]);
            }
        }

        iftFree(aux);
        *Z_array = NULL;
    }
}


iftImage *iftDataSetObjectMap(const iftDataSet *Z, const iftImage *comp, int max_val, int label) {
    if (Z == NULL)
        iftError("DataSet is NULL", "iftDataSetObjectMap");
    if (Z->ref_data == NULL)
        iftError("Reference Data is NULL", "iftDataSetObjectMap");
    if ((Z->ref_data_type != IFT_REF_DATA_IMAGE)  &&
        (Z->ref_data_type != IFT_REF_DATA_FIMAGE) &&
        (Z->ref_data_type != IFT_REF_DATA_MIMAGE))
        iftError("Reference Data Type should be an Image (or FImage or MImage)", "iftDataSetObjectMap");


    iftImage *weight_img = NULL;

    if (Z->ref_data_type == IFT_REF_DATA_IMAGE){
        iftImage *img = (iftImage*) Z->ref_data;
        weight_img    = iftCreateImage(img->xsize, img->ysize, img->zsize);
        iftCopyVoxelSize(img, weight_img);
    } else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE) {
        iftFImage *fimg = (iftFImage*) Z->ref_data;
        weight_img      = iftCreateImage(fimg->xsize, fimg->ysize, fimg->zsize);
        iftCopyVoxelSize(fimg, weight_img);
    } else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE) {
        iftMImage *mimg = (iftMImage*) Z->ref_data;
        weight_img      = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
        iftCopyVoxelSize(mimg, weight_img);
    }

    float min_weight = 0, max_weight = 0;
    iftGetMinMaxWeightInDataSet(Z, &min_weight, &max_weight);

    if (!iftAlmostZero(max_weight-min_weight)){
        // consider all labels != 1 (background)
        if (label == -1) {
#pragma omp parallel for schedule(auto)
            for (int s = 0; s < Z->nsamples; s++) {
                if (Z->sample[s].label != 1) {
                    weight_img->val[Z->sample[s].id] = (int) (((Z->sample[s].weight - min_weight) / (max_weight - min_weight)) * max_val);
                }
            }
        }
            // gets the fuzzy map only for the required label
        else {
#pragma omp parallel for schedule(auto)
            for (int s = 0; s < Z->nsamples; s++) {
                if (Z->sample[s].label == label) {
                    weight_img->val[Z->sample[s].id] = (int) (((Z->sample[s].weight - min_weight) / (max_weight - min_weight)) * max_val);
                }
            }
        }
    }

    if (comp != NULL)
    {
        if (iftMaximumValue(comp) != Z->nsamples) {
            iftError("Number of compoments different from number of dataset samples", "iftDataSetObjectMap");
        }
        for (int p = 0; p < weight_img->n; p++)
        {
            int c = comp->val[p] - 1;
            weight_img->val[p] = weight_img->val[Z->sample[c].id];
        }
    }

    return weight_img;
}


iftDataSet *iftImageToDataSet(iftImage *img)
{
    iftDataSet *Z;
    int p;

    if (iftIsColorImage(img)){
        Z=iftCreateDataSet(img->n,3);
        for (p=0; p < Z->nsamples; p++) {
            Z->sample[p].feat[0] = img->val[p];
            Z->sample[p].feat[1] = img->Cb[p];
            Z->sample[p].feat[2] = img->Cr[p];
            Z->sample[p].id   = p;
        }
    }else{ /* Gray scale image */
        Z        = iftCreateDataSet(img->n,1);
        for (p=0; p < Z->nsamples; p++) {
            Z->sample[p].feat[0] = img->val[p];
            Z->sample[p].id      = p;
        }
    }
    iftSetStatus(Z,IFT_TEST);
    iftCopyRefData(Z, img, IFT_REF_DATA_IMAGE);
    // Z->ref_data      = img;
    // Z->ref_data_type = IFT_REF_DATA_IMAGE;

    return(Z);
}

iftDataSet *iftImageBoundingBoxToDataSet(iftImage *img, iftBoundingBox bb )
{
    if (img == NULL)
        iftError("Image is NULL", "iftImageBoundingBoxToDataSet");

    iftVoxel upper_left_voxel = bb.begin;
    iftVoxel bottom_right_voxel = bb.end;

    if (!iftValidVoxel(img, upper_left_voxel))
        iftError(
                "Initial Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for image with size (%d, %d, %d)",
                "iftImageBoundingBoxToDataSet", upper_left_voxel.x, upper_left_voxel.y, upper_left_voxel.z, img->xsize,
                img->ysize, img->zsize);
    if (!iftValidVoxel(img, bottom_right_voxel)) {
        iftWarning("Final Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for image with size (%d, %d, %d).\n It will copied/inserted what it is possible", "iftImageBoundingBoxToDataSet",bottom_right_voxel.x, bottom_right_voxel.y, bottom_right_voxel.z, img->xsize, img->ysize, img->zsize);
        // gets the valid ending voxel
        bottom_right_voxel.x = iftMin(bottom_right_voxel.x, img->xsize - 1);
        bottom_right_voxel.y = iftMin(bottom_right_voxel.y, img->ysize - 1);
        bottom_right_voxel.z = iftMin(bottom_right_voxel.z, img->zsize - 1);
    }

    iftDataSet *Z;
    int bb_voxel_number= (bottom_right_voxel.x - upper_left_voxel.x +1)*(bottom_right_voxel.y-upper_left_voxel.y +1)*(bottom_right_voxel.z-upper_left_voxel.z+1);

    if (iftIsColorImage(img)) {
        Z=iftCreateDataSet(bb_voxel_number,3);
        int q = 0;
        iftVoxel u;

        for (u.z = upper_left_voxel.z; u.z <= bottom_right_voxel.z; u.z++)
            for (u.y = upper_left_voxel.y; u.y <= bottom_right_voxel.y; u.y++)
                for (u.x = upper_left_voxel.x; u.x <= bottom_right_voxel.x; u.x++) {
                    int p = iftGetVoxelIndex(img, u);
                    Z->sample[q].feat[0] = img->val[p];
                    Z->sample[q].feat[1] = img->Cb[p];
                    Z->sample[q].feat[2] = img->Cr[p];
                    Z->sample[q].id   = p;
                    q++;
                }
    }
    else {
        Z=iftCreateDataSet(bb_voxel_number,1);
        int q = 0;
        iftVoxel u;

        for (u.z = upper_left_voxel.z; u.z <= bottom_right_voxel.z; u.z++)
            for (u.y = upper_left_voxel.y; u.y <= bottom_right_voxel.y; u.y++)
                for (u.x = upper_left_voxel.x; u.x <= bottom_right_voxel.x; u.x++) {
                    int p = iftGetVoxelIndex(img, u);
                    Z->sample[q].feat[0] = img->val[p];
                    Z->sample[q].id   = p;
                    q++;
                }
    }

    iftSetStatus(Z,IFT_TEST);
    iftCopyRefData(Z, img, IFT_REF_DATA_IMAGE);
    // Z->ref_data_type = IFT_REF_DATA_IMAGE;
    // Z->ref_data      = img;

    return (Z);
}

iftDataSet *iftMImageBoundingBoxToDataSet(iftDataSet *orig_dataset, iftBoundingBox bb, bool copy_feats )
{
    if (orig_dataset->ref_data_type != IFT_REF_DATA_MIMAGE)
        iftError("The referenced data type of the dataset is not an mimage ", "iftMImageBoundingBoxToDatasetArray");

    iftMImage *mimg=(iftMImage *)orig_dataset->ref_data;

    if (mimg == NULL)
        iftError("Image is NULL", "iftMImageBoundingBoxToDataSet");

    iftVoxel upper_left_voxel = bb.begin;
    iftVoxel bottom_right_voxel = bb.end;

    if (!iftValidVoxel(mimg, upper_left_voxel))
        iftError(
                "Initial Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for image with size (%d, %d, %d)",
                "iftImageBoundingBoxToDataSet", upper_left_voxel.x, upper_left_voxel.y, upper_left_voxel.z, mimg->xsize,
                mimg->ysize, mimg->zsize);
    if (!iftValidVoxel(mimg, bottom_right_voxel)) {
        iftWarning("Final Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for image with size (%d, %d, %d).\n It will copied/inserted what it is possible", "iftImageBoundingBoxToDataSet", bottom_right_voxel.x, bottom_right_voxel.y, bottom_right_voxel.z, mimg->xsize, mimg->ysize, mimg->zsize);
        // gets the valid ending voxel
        bottom_right_voxel.x = iftMin(bottom_right_voxel.x, mimg->xsize - 1);
        bottom_right_voxel.y = iftMin(bottom_right_voxel.y, mimg->ysize - 1);
        bottom_right_voxel.z = iftMin(bottom_right_voxel.z, mimg->zsize - 1);
    }

    iftDataSet *Z;
    int bb_voxel_number= (bottom_right_voxel.x - upper_left_voxel.x +1)*(bottom_right_voxel.y-upper_left_voxel.y +1)*(bottom_right_voxel.z-upper_left_voxel.z+1);

    if (copy_feats)
        Z=iftCreateDataSet(bb_voxel_number,mimg->m);
    else
        Z= iftCreateDataSetWithoutFeatsVector(bb_voxel_number, mimg->m);

    int q,x_size,y_size,p,i,j,k;
    iftVoxel u;
    x_size=bottom_right_voxel.x-upper_left_voxel.x+1;
    y_size=bottom_right_voxel.y-upper_left_voxel.y+1;

    for (int z = upper_left_voxel.z; z <= bottom_right_voxel.z; z++)
        for (int y = upper_left_voxel.y; y <= bottom_right_voxel.y; y++)
#pragma omp parallel for private(p,q,u,i,j,k)
            for (int x = upper_left_voxel.x; x <= bottom_right_voxel.x; x++) {
                k=z-upper_left_voxel.z;
                j=y-upper_left_voxel.y;
                i=x-upper_left_voxel.x;
                u.x=x;u.y=y;u.z=z;
                p=iftGetVoxelIndex(mimg, u);
                q= i + j*x_size + k*x_size*y_size;
                iftCopySample(&orig_dataset->sample[p],&Z->sample[q],orig_dataset->nfeats,copy_feats);
            }

    Z->nclasses=orig_dataset->nclasses;
    Z->ngroups=orig_dataset->ngroups;
    Z->ref_data_type = IFT_REF_DATA_MIMAGE;
    /* actually the ref data will be a pointer to the current bounding box*/
    iftBoundingBox *ref_bb=(iftBoundingBox *)iftAlloc(1,sizeof(iftBoundingBox));
    ref_bb->begin=bb.begin;
    ref_bb->end=bb.end;
    Z->ref_data=ref_bb;
    iftSetStatus(Z,IFT_TEST);

    return (Z);
}

int iftMImageTilesToDataSetArray(iftDataSet *orig_dataset, int n_tiles,iftDataSet *** out_datasets, bool copy_feats){

    if (n_tiles <= 0)
        iftError("Number of tiles %d <= 0", "iftMImageTilesToDataSetArray", n_tiles);
    if (orig_dataset->ref_data_type != IFT_REF_DATA_MIMAGE)
        iftError("The referenced data type of the dataset is not an mimage ", "iftImageTilesToDatasetArray");

    iftMImage *mimg=(iftMImage *)orig_dataset->ref_data;

    if (mimg == NULL)
        iftError("Image is NULL", "iftImageTilesToDataSetArray");

    // tries to find regular tiles for the given number of tiles
    int block_vol = mimg->n / (1.0 * n_tiles);

    int n_blocks_x, n_blocks_y, n_blocks_z;
    iftNumberOfEquallyDimensionedTilesWithGivenMaximumSize(mimg->xsize, mimg->ysize, mimg->zsize,block_vol, &n_blocks_x, &n_blocks_y, &n_blocks_z);

    iftImageTiles *tiles = iftMImageTilesByEquallyDividingAxes(mimg, n_blocks_x, n_blocks_y, n_blocks_z);

    iftDataSet **aux=iftCreateDataSetArray(tiles->ntiles);

#pragma omp parallel for
    for (int b = 0; b < tiles->ntiles; b++) {
        iftBoundingBox bb     = tiles->tile_coords[b];
        aux[b]=iftMImageBoundingBoxToDataSet(orig_dataset,bb,copy_feats);
    }

    *out_datasets=aux;
    int dataset_count=tiles->ntiles;
    iftDestroyImageTiles(&tiles);


    return dataset_count;
}

int iftSplitDatasetIntoRandomlyDataSetArray(iftDataSet *Z, int nb_partitions, iftDataSet ***datasets, bool
copy_features){

    if (nb_partitions <= 0)
        iftError("Number of partitions for dataset %d <= 0", "iftSplitDatasetIntoRandomlyDataSetArray",
                 nb_partitions);

    if (Z == NULL)
        iftError("Dataset to be splitted is NULL", "iftSplitDatasetIntoRandomlyDataSetArray");

    int nb_samples_for_dataset=Z->nsamples/nb_partitions;
    int nb_samples_last_dataset = Z->nsamples - nb_samples_for_dataset *(nb_partitions - 1);

    /* if the sample number of the last dataset are less than the half of the number of elements of dataset, add these samples to the previous dataset. I'm done this to forbid a dataset with a little sample number */
    iftDataSet **aux;
    if (nb_samples_last_dataset < nb_samples_for_dataset/2 && nb_partitions > 1){
        nb_partitions--;
        nb_samples_last_dataset = nb_samples_last_dataset + nb_samples_for_dataset;
    }

    aux=iftCreateDataSetArray(nb_partitions);

    for (int i=0;i<nb_partitions-1;i++){
        if (copy_features)
            aux[i]=iftCreateDataSet(nb_samples_for_dataset, Z->nfeats);
        else
            aux[i]= iftCreateDataSetWithoutFeatsVector(nb_samples_for_dataset, Z->nfeats);
        iftCopyRefData(aux[i], Z->ref_data, Z->ref_data_type);
        // aux[i]->ref_data_type = Z->ref_data_type;
        // aux[i]->ref_data      = Z->ref_data;
        aux[i]->nclasses=Z->nclasses;
        aux[i]->ngroups=Z->ngroups;
    }
    if (copy_features)
        aux[nb_partitions-1]=iftCreateDataSet(nb_samples_last_dataset, Z->nfeats);
    else
        aux[nb_partitions-1]= iftCreateDataSetWithoutFeatsVector(nb_samples_last_dataset, Z->nfeats);
    iftCopyRefData(aux[nb_partitions-1], Z->ref_data, Z->ref_data_type);
    // aux[nb_partitions-1]->ref_data_type = Z->ref_data_type;
    // aux[nb_partitions-1]->ref_data      = Z->ref_data;
    aux[nb_partitions-1]->nclasses=Z->nclasses;
    aux[nb_partitions-1]->ngroups=Z->ngroups;

    int *shuffled_indexes= iftAllocIntArray(Z->nsamples);
    for (int i=0; i < Z->nsamples; i++)
        shuffled_indexes[i]=i;
    iftShuffleIntArray(shuffled_indexes, Z->nsamples);

#pragma omp parallel for shared(Z,shuffled_indexes,aux)
    for (int i=0;i<nb_partitions-1;i++){
        int q=0;
        for (int j=i*nb_samples_for_dataset; j < (i+1)*nb_samples_for_dataset;j++){
            int s=shuffled_indexes[j];
            iftCopySample(&Z->sample[s], &aux[i]->sample[q], Z->nfeats, copy_features);
            q++;
        }
    }

    /*copy the elements of the last dataset */
    int q=0;
    for (int j=(nb_partitions-1)*nb_samples_for_dataset;j< Z->nsamples; j++){
        int s=shuffled_indexes[j];
        iftCopySample(&Z->sample[s], &aux[nb_partitions - 1]->sample[q], Z->nfeats, copy_features);
        q++;
    }

    *datasets=aux;
    iftFree(shuffled_indexes);
    return nb_partitions;
}


iftDataSet *iftSupervoxelsToDataSet(const iftImage *img, const iftImage *supervoxel,
                                    const iftIntArray *true_labels) {
    iftDataSet *Z;
    int r,s,p,nregions=iftMaximumValue(supervoxel);
    int *size=iftAllocIntArray(nregions);

    if (iftMinimumValue(supervoxel)<=0)
        iftError("Minimum label value must be 1", "iftImageCompsToDataSet");

    iftVerifyImageDomains(img,supervoxel,"iftSupervoxelsToDataset");

    if (iftIsColorImage(img)){
        Z=iftCreateDataSet(nregions,3);
        for (p=0; p < img->n; p++) {
            s = supervoxel->val[p] - 1;
            size[s]++;
            Z->sample[s].feat[0] += img->val[p];
            Z->sample[s].feat[1] += img->Cb[p];
            Z->sample[s].feat[2] += img->Cr[p];
        }

        for(r = 0; r < nregions; r++){
            Z->sample[r].feat[0] = Z->sample[r].feat[0]/size[r];
            Z->sample[r].feat[1] = Z->sample[r].feat[1]/size[r];
            Z->sample[r].feat[2] = Z->sample[r].feat[2]/size[r];
        }

    }else{ /* Gray scale image */
        Z        = iftCreateDataSet(nregions,1);
        for (p=0; p < img->n; p++) {
            s = supervoxel->val[p] - 1;
            size[s]++;
            Z->sample[s].feat[0] += img->val[p];
        }

        for(r = 0; r < nregions; r++){
            Z->sample[r].feat[0] = Z->sample[r].feat[0]/size[r];
        }
    }

    iftFree(size);

    //This may seem redundant, but it is useful for splitting the dataset
    for(r = 0; r < nregions; r++) {
        Z->sample[r].id = r + 1;
    }

    iftSetDistanceFunction(Z, 1);
    iftCopyRefData(Z, supervoxel, IFT_REF_DATA_IMAGE);
    // Z->ref_data      = (void*) supervoxel;
    // Z->ref_data_type = IFT_REF_DATA_IMAGE;

    // assigns the true labels
    if (true_labels == NULL) {
        for (int r = 0; r < nregions; r++)
            Z->sample[r].truelabel = r+1; // label from the supervoxel
        Z->nclasses = nregions;
    }
    else {
        for (int r = 0; r < nregions; r++) {
            Z->sample[r].truelabel = true_labels->val[r]; // label from the supervoxel
        }
        Z->nclasses = iftCountUniqueIntElems(true_labels->val, true_labels->n);
        iftIsThereTrueLabelOrGroupZero(Z, "iftSupervoxelsToDataSet");
    }

    return(Z);
}


//nbins should be n^3 for some natural number n
iftDataSet *iftSupervoxelsToHistogramDataSet(const iftImage *img, const iftImage *label, int nbins, iftIntArray *true_labels){
    int nregions = iftMaximumValue(label);
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));
    int nbits = 1;

    while ((1 << nbits) <= normalization_value) { /* number of bits per voxel */
        nbits++;
    }

    iftVerifyImageDomains(img,label,"iftSupervoxelsToHistogramDataSet");
    if (iftMinimumValue(label)<=0)
        iftError("Minimum label value must be 1", "iftSupervoxelsToHistogramDataSet");

    iftDataSet* Z = iftCreateDataSet(nregions,nbins);
    int *region_size = iftAllocIntArray(nregions);
    if (iftIsColorImage(img)){
        int xsize   = iftRound(pow(nbins, 1.0 / 3.0));
        int xysize  = xsize*xsize;

        if (xsize < 1)
            iftError("Insufficient number of bins", "iftSupervoxelsToHistogramDataSet");
        if (xsize > normalization_value)
            iftError("Excessive number of bins", "iftSupervoxelsToHistogramDataSet");

        float binsize = ((float)pow(2, nbits))/xsize;

        iftColor RGB, YCbCr;

        int p,b,r;
        for (p=0; p < img->n; p++){
            r = label->val[p] - 1;
            region_size[r]++;

            YCbCr.val[0] = img->val[p];
            YCbCr.val[1] = img->Cb[p];
            YCbCr.val[2] = img->Cr[p];
            RGB = iftYCbCrtoRGB(YCbCr,normalization_value);

            b = (int)( ((float)RGB.val[0]) / binsize )
                + xsize*(int)( ((float)RGB.val[1]) / binsize )
                + xysize*(int)( ((float)RGB.val[2]) / binsize );

            Z->sample[r].feat[b]++;
        }

        for(r = 0; r < nregions; r++){
            for (b=0; b < nbins; b++) {
                Z->sample[r].feat[b] /= region_size[r];
            }
        }
    }else{ /* Gray scale image */

        float binsize = ((float)pow(2, nbits))/nbins;

        int p,r;
        for (p=0; p < img->n; p++){
            r = label->val[p] - 1;
            region_size[r]++;
            Z->sample[r].feat[(int)(img->val[p]/binsize)]++;
        }

        int b;
        for(r = 0; r < nregions; r++){
            for (b=0; b < nbins; b++) {
                Z->sample[r].feat[b] /= region_size[r];
            }
        }
    }

    iftFree(region_size);

    //This may seem redundant, but it is useful for splitting the dataset
    int r;
    for(r = 0; r < nregions; r++){
        Z->sample[r].id = r + 1;
    }


    iftSetDistanceFunction(Z, 1);
    iftCopyRefData(Z, label, IFT_REF_DATA_IMAGE);
    // Z->ref_data      = (void*) label;
    // Z->ref_data_type = IFT_REF_DATA_IMAGE;

    // assigns the true labels
    if (true_labels != NULL) {
        for (int r = 0; r < nregions; r++) {
            Z->sample[r].truelabel = true_labels->val[r]; // label from the supervoxel
        }
        Z->nclasses = iftCountUniqueIntElems(true_labels->val, true_labels->n);
        iftIsThereTrueLabelOrGroupZero(Z, "iftSupervoxelsToHistogramDataSet");
    }

    return(Z);
}

iftDataSet *iftSupervoxelsToLabOrGrayHistogramDataSet(const iftImage *img, const iftImage *label, int nbins, int *region_size){

    int nregions = iftMaximumValue(label);
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));
    int nbits = 1;

    while ((1 << nbits) <= normalization_value) { /* number of bits per voxel */
        nbits++;
    }

    iftVerifyImageDomains(img,label,"iftSupervoxelsToLabOrGrayHistogramDataSet");
    if (iftMinimumValue(label)<=0)
        iftError("Minimum label value must be 1", "iftSupervoxelsToLabOrGrayHistogramDataSet");

    iftDataSet* Z = iftCreateDataSet(nregions,nbins);

    if (iftIsColorImage(img)){
        int xsize   = iftRound(pow(nbins, 1.0 / 3.0));
        int xysize  = xsize*xsize;


        if (xsize < 1)
            iftError("Insufficient number of bins", "iftSupervoxelsToLabOrGrayHistogramDataSet");
        if (xsize > normalization_value)
            iftError("Excessive number of bins", "iftSupervoxelsToLabOrGrayHistogramDataSet");

        float binsize = ((float)pow(2, nbits))/xsize;

#pragma omp parallel for
        for (int p=0; p < img->n; p++){
            int r = label->val[p] - 1;
            region_size[r]++;
            iftColor YCbCr;
            iftFColor Lab;

            YCbCr.val[0] = img->val[p];
            YCbCr.val[1] = img->Cb[p];
            YCbCr.val[2] = img->Cr[p];
            Lab=iftYCbCrtoLabNorm2(YCbCr,normalization_value);

            int b = (int)( (Lab.val[0]*normalization_value) / binsize )
                + xsize*(int)( (Lab.val[1]*normalization_value) / binsize )
                + xysize*(int)( (Lab.val[2]*normalization_value) / binsize );

            Z->sample[r].feat[b]++;
        }

    }else{ /* Gray scale image */

        float binsize = ((float)pow(2, nbits))/nbins;

        int p,r;
        for (p=0; p < img->n; p++){
            r = label->val[p] - 1;
            region_size[r]++;
            Z->sample[r].feat[(int)(img->val[p]/binsize)]++;
        }
    }

    iftSetDistanceFunction(Z,12);
    iftCopyRefData(Z, label, IFT_REF_DATA_IMAGE);
    // Z->ref_data      = (void*) label;
    // Z->ref_data_type = IFT_REF_DATA_IMAGE;

    return(Z);
}



typedef void (*iftExtractAdjacencyFeats)(const iftImage *img, int p, const iftAdjRel *A, iftSample *sample);

void iftExtractGrayAdjacencyFeats(const iftImage *img, int p, const iftAdjRel *A, iftSample *sample) {
    iftVoxel u = iftGetVoxelCoord(img, p);

    #pragma omp parallel for
    for (int i = 0; i < A->n; i++) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, i);

        if (iftValidVoxel(img, v))
            (*sample).feat[i] = iftImgVoxelVal(img, v);
    }
}

void iftExtractColorAdjacencyFeats(const iftImage *img, int p, const iftAdjRel *A, iftSample *sample) {
    iftVoxel u = iftGetVoxelCoord(img, p);

    for (int i = 0; i < A->n; i++) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, i);

        if (iftValidVoxel(img, v)) {
            (*sample).feat[i]              = iftImgVoxelVal(img, v);
            (*sample).feat[i + A->n]       = iftImgVoxelCb(img, v);
            (*sample).feat[i + (2 * A->n)] = iftImgVoxelCr(img, v);
        }
    }
}


iftDataSet *iftImageToDataSetByAdjacencyFeats(const iftImage *img, const iftImage *label_img,
                                              const iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = (iftIs3DImage(img)) ? iftSpheric(sqrtf(3.0)) : iftCircular(sqrtf(2.0));
    else A = iftCopyAdjacency(Ain);

    iftDataSet *Z = NULL;
    iftExtractAdjacencyFeats extractAdjacencyFeatsFunc = iftExtractGrayAdjacencyFeats;
    int nfeats = A->n;
    if (iftIsColorImage(img)) {
        nfeats = 3 * A->n;
        extractAdjacencyFeatsFunc = iftExtractColorAdjacencyFeats;
    }

    if (label_img) {
        iftVerifyImages(img, label_img, "iftImageToDataSetByAdjacencyFeats");
        iftIntArray *labels = iftGetObjectLabels(label_img);

        int n_obj_spels = img->n - iftCountBackgroundSpels(label_img);
        Z = iftCreateDataSet(n_obj_spels, nfeats);
        Z->nclasses = labels->n;
        iftDestroyIntArray(&labels);

        int s = 0;
        for (int p = 0; p < img->n; p++) {
            if (label_img->val[p] != 0) {
                extractAdjacencyFeatsFunc(img, p, A, &Z->sample[s]);
                Z->sample[s].truelabel = label_img->val[p];
                Z->sample[s].id        = p;
                s++;
            }
        }
    }
    else {
        Z = iftCreateDataSet(img->n, nfeats);

        #pragma omp parallel for
        for (int p = 0; p < img->n; p++) {
            extractAdjacencyFeatsFunc(img, p, A, &Z->sample[p]);
            Z->sample[p].id = p;
        }
    }

    iftCopyRefData(Z, img, IFT_REF_DATA_IMAGE);

    if (Ain != NULL)
        iftDestroyAdjRel(&A);

    return Z;
}


iftDataSet *iftImageROIToDataSetByAdjacencyFeats(const iftImage *img, const iftImage *label_img, iftBoundingBox bb, const iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = (iftIs3DImage(img)) ? iftSpheric(sqrtf(3.0)) : iftCircular(sqrtf(2.0));
    else A = iftCopyAdjacency(Ain);

    iftExtractAdjacencyFeats extractAdjacencyFeatsFunc = iftExtractGrayAdjacencyFeats;
    int nfeats = A->n;
    if (iftIsColorImage(img)) {
        nfeats = 3 * A->n;
        extractAdjacencyFeatsFunc = iftExtractColorAdjacencyFeats;
    }

    iftDataSet *Z = iftCreateDataSet(iftBoundingBoxBoxVolume(bb), nfeats);

    int s = 0;
    for (int z = bb.begin.z; z <= bb.end.z; z++)
        for (int y = bb.begin.y; y <= bb.end.y; y++)
            for (int x = bb.begin.x; x <= bb.end.x; x++) {
                int p = iftGetVoxelIndex(img, ((iftVoxel) {x, y, z}));
                Z->sample[s].id = p;
                if (label_img->val[p] != 0)
                    extractAdjacencyFeatsFunc(img, p, A, &Z->sample[s]);
                s++;
            }

    iftCopyRefData(Z, img, IFT_REF_DATA_IMAGE);

    if (Ain != NULL)
        iftDestroyAdjRel(&A);

    return Z;
}


iftDataSet *iftObjectToDataSet(iftImage *label)
{
    iftDataSet *Z;
    int s, p, nsamples;

    nsamples = 0;
    for (p=0; p < label->n; p++)
        if (label->val[p]!=0)
            nsamples++;

    if (iftIs3DImage(label)) {
        Z        = iftCreateDataSet(nsamples,3);
        for (p=0,s=0; p < label->n; p++){
            if (label->val[p]!=0){
                Z->sample[s].feat[0] = (float)iftGetXCoord(label,p);
                Z->sample[s].feat[1] = (float)iftGetYCoord(label,p);
                Z->sample[s].feat[2] = (float)iftGetZCoord(label,p);
                Z->sample[s].id      = p;
                s++;
            }
        }
    } else { /* 2D image */
        Z        = iftCreateDataSet(nsamples,2);
        for (p=0,s=0; p < label->n; p++){
            if (label->val[p]!=0){
                Z->sample[s].feat[0] = (float)iftGetXCoord(label,p);
                Z->sample[s].feat[1] = (float)iftGetYCoord(label,p);
                Z->sample[s].id   = p;
                s++;
            }
        }
    }

    iftCopyRefData(Z, label, IFT_REF_DATA_IMAGE);
    // Z->ref_data      = label;
    // Z->ref_data_type = IFT_REF_DATA_IMAGE;

    return(Z);
}


iftDataSet *iftImageSeedsToDataSet(iftImage *img, iftLabeledSet *S)
{
    iftDataSet *Z;
    int s,p,nsamples,nclasses,minclass,maxclass,inc=0;
    iftLabeledSet *Saux;

    Saux=S; nsamples = 0; minclass= IFT_INFINITY_INT; maxclass= IFT_INFINITY_INT_NEG;
    while (Saux != NULL) {
        nsamples++;
        if (Saux->label < minclass)
            minclass = Saux->label;
        if (Saux->label > maxclass)
            maxclass = Saux->label;
        Saux = Saux->next;
    }

    if (minclass == 0)
        inc = 1;

    nclasses = maxclass+inc;

    if (iftIsColorImage(img)){
        Z=iftCreateDataSet(nsamples,3);
        Saux=S; s = 0;
        while (Saux != NULL) {
            p = Saux->elem;
            Z->sample[s].feat[0] = img->val[p];
            Z->sample[s].feat[1] = img->Cb[p];
            Z->sample[s].feat[2] = img->Cr[p];
	    iftAddSampleStatus(&Z->sample[s], IFT_TRAIN);
	    iftAddSampleStatus(&Z->sample[s], IFT_SUPERVISED);
            Z->sample[s].truelabel   = Saux->label+inc;
            Z->sample[s].id   = p;
            s++;
            Saux = Saux->next;
        }
    }else{ /* Gray scale image */
        Z=iftCreateDataSet(nsamples,1);
        Saux=S; s = 0;
        while (Saux != NULL) {
            p = Saux->elem;
            Z->sample[s].feat[0] = img->val[p];
            Z->sample[s].truelabel   = Saux->label+inc;
	    iftAddSampleStatus(&Z->sample[s], IFT_TRAIN);
	    iftAddSampleStatus(&Z->sample[s], IFT_SUPERVISED);
            Z->sample[s].id = p;
            s++;
            Saux = Saux->next;
        }
    }

    iftCopyRefData(Z, img, IFT_REF_DATA_IMAGE);
    // Z->ref_data      = img;
    // Z->ref_data_type = IFT_REF_DATA_IMAGE;
    Z->nclasses      = nclasses;

    return(Z);
}

iftDataSet *iftMImageSeedsToDataSet(iftMImage *img, const iftLabeledSet *S)
{
    iftDataSet *Z;
    int s,p,nsamples,nclasses,minclass,maxclass,inc=0;
    const iftLabeledSet *Saux;

    Saux=S; nsamples = 0; minclass= IFT_INFINITY_INT; maxclass= IFT_INFINITY_INT_NEG;
    while (Saux != NULL) {
        nsamples++;
        if (Saux->label < minclass)
            minclass = Saux->label;
        if (Saux->label > maxclass)
            maxclass = Saux->label;
        Saux = Saux->next;
    }

    if (minclass == 0)
        inc = 1;

    nclasses = maxclass+inc;

    Z=iftCreateDataSet(nsamples,img->m);
    Saux=S; s = 0;
    while (Saux != NULL) {
        p = Saux->elem;
        for (int b=0; b < img->m; b++)
            Z->sample[s].feat[b] = img->val[p][b];
        Z->sample[s].truelabel   = Saux->label+inc;
	iftAddSampleStatus(&Z->sample[s], IFT_TRAIN);
	iftAddSampleStatus(&Z->sample[s], IFT_SUPERVISED);
        Z->sample[s].id   = p;
        s++;
        Saux = Saux->next;
    }

    iftCopyRefData(Z, img, IFT_REF_DATA_MIMAGE);
    // Z->ref_data      = img;
    // Z->ref_data_type = IFT_REF_DATA_MIMAGE;
    Z->nclasses = nclasses;

    return(Z);
}


iftDataSet *iftImageRegionToDataSet(iftImage *img, iftImage *region)
{
    iftDataSet *Z;
    int p,s,nsamples;

    iftVerifyImageDomains(img,region,"iftImageRegionToDataSet");

    nsamples = 0;
    for (p=0; p < img->n; p++)
        if (region->val[p]!=0)
            nsamples++;

    s=0;
    if (iftIsColorImage(img)){
        Z=iftCreateDataSet(nsamples,3);
        for (p=0; p < img->n; p++) {
            if (region->val[p]!=0) {
                Z->sample[s].feat[0] = img->val[p];
                Z->sample[s].feat[1] = img->Cb[p];
                Z->sample[s].feat[2] = img->Cr[p];
		iftAddSampleStatus(&Z->sample[s], IFT_TRAIN);
		iftAddSampleStatus(&Z->sample[s], IFT_SUPERVISED);
                Z->sample[s].id   = p;
                s++;
            }
        }
    }else{
        Z=iftCreateDataSet(nsamples,1);
        for (p=0; p < img->n; p++) {
            if (region->val[p]!=0){
                Z->sample[s].feat[0] = img->val[p];
                Z->sample[s].id   = p;
		iftAddSampleStatus(&Z->sample[s], IFT_TRAIN);
		iftAddSampleStatus(&Z->sample[s], IFT_SUPERVISED);
                s++;
            }
        }
    }

    iftCopyRefData(Z, img, IFT_REF_DATA_IMAGE);
    // Z->ref_data      = img;
    // Z->ref_data_type = IFT_REF_DATA_IMAGE;

    return(Z);
}

void iftResetDataSet(iftDataSet *Z)
{
    int s;

    Z->ntrainsamples=0;
    for (s=0; s < Z->nsamples; s++) {
        iftSetSampleStatus(&Z->sample[s], IFT_UNKNOWN);
        Z->sample[s].label  = 0;
        Z->sample[s].weight = 0.0;
    }
}

int iftSelectUnsupTrainSamples(iftDataSet *Z, float train_perc, int marked_times)
{

    int s, *sample, i, t, high;
    int    *count;

    // Verify if it is the trivial case of selecting all samples.

    if (train_perc == 1.0) {
        for (s=0; s < Z->nsamples; s++)
            iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
        Z->ntrainsamples = Z->nsamples;
        return(Z->nsamples);
    }

    if ((train_perc <= 0.0))
        iftError("Invalid percentage of training samples", "iftSelectUnsupTrainSamples");

    // Reset status

    iftSetStatus(Z, IFT_TEST);

    // Compute the number of training samples

    /* check if was passed a value of a percent in train_perc parameter*/
    if (train_perc < 1.0){
        Z->ntrainsamples = (int) (train_perc*Z->nsamples);
    }
    else
        Z->ntrainsamples = (int)train_perc;


    if (Z->ntrainsamples == 0)
        iftError("Percentage is too low for this dataset", "iftSelectUnsupTrainSamples");

    // Prepare samples for selection

    sample = iftAllocIntArray(Z->nsamples);
    count  = iftAllocIntArray(Z->nsamples);
    t=0;
    for (s=0; s < Z->nsamples; s++) {
        sample[t]=s;
        count[t]=marked_times;
        t++;
    }

    // Randomly select samples

    t = 0; high = Z->nsamples-1;
    while (t < Z->ntrainsamples) {
        i = iftRandomInteger(0,high);
        s = sample[i];
        if (count[i]==0){
            iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
            iftSwap(sample[i], sample[high]);
            iftSwap(count[i], count[high]);
            t++; high--;
        }else{
            count[i]--;
        }
    }
    iftFree(count);
    iftFree(sample);

    return(Z->ntrainsamples);
}

int iftSelectUnsupTrainSamplesByWeight(iftDataSet *Z, float train_perc)
{
    int s, *sample, i, t, high;
    int    *count;

    // Verify if it is the trivial case of selecting all samples.

    if (train_perc == 1.0) {
        for (s=0; s < Z->nsamples; s++)
            iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
        Z->ntrainsamples = Z->nsamples;
        return(Z->nsamples);
    }

    if ((train_perc <= 0.0)||(train_perc > 1.0))
        iftError("Invalid percentage of training samples", "iftSelectUnsupTrainSamples");

    // Reset status

    iftSetStatus(Z, IFT_TEST);

    // Compute the number of training samples

    Z->ntrainsamples = (int) (train_perc*Z->nsamples);
    if (Z->ntrainsamples == 0)
        iftError("Percentage is too low for this dataset", "iftSelectUnsupTrainSamples");

    // Prepare samples for selection

    sample = iftAllocIntArray(Z->nsamples);
    count  = iftAllocIntArray(Z->nsamples);
    t=0;
    for (s=0; s < Z->nsamples; s++) {
        sample[t]=s;
        count[t]= 100 * (IFT_MAXWEIGHT - Z->sample[s].weight) / IFT_MAXWEIGHT;
        t++;
    }

    // Randomly select samples


    t = 0; high = Z->nsamples-1;
    while (t < Z->ntrainsamples) {
        i = iftRandomInteger(0,high);
        s = sample[i];
        if (count[i]==0){
            iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
            iftSwap(sample[i], sample[high]);
            iftSwap(count[i], count[high]);
            t++; high--;
        }else{
            count[i]--;
        }
    }
    iftFree(count);
    iftFree(sample);

    return(Z->ntrainsamples);
}

void iftSetStatus(iftDataSet *Z, iftSampleStatus status) {
    for (int s = 0; s < Z->nsamples; s++)
        iftSetSampleStatus(&Z->sample[s], status);

    Z->ntrainsamples = (status & IFT_TRAIN) ? Z->nsamples : 0;
}


void iftAddStatus(iftDataSet *Z, iftSampleStatus status) {
    for (int s = 0; s < Z->nsamples; s++)
        iftAddSampleStatus(&Z->sample[s], status);

    // update the number of training samples
    if (status & IFT_TRAIN)
      Z->ntrainsamples = Z->nsamples;
}

void iftRemoveStatus(iftDataSet *Z, iftSampleStatus status) {
    // uint status_comp = ~status; // status complement (not operator)

    for (int s = 0; s < Z->nsamples; s++)
        iftRemoveSampleStatus(&Z->sample[s], status);

    // update the number of training samples
    if (status & IFT_TRAIN)
        Z->ntrainsamples = 0;
}

void iftSetStatusForSet(iftDataSet *Z, iftSet *set, iftSampleStatus status){
    iftSet *s = set;
    while(s){
        iftSetSampleStatus(&Z->sample[s->elem], status);
        s = s->next;
    }

    int ntrain = 0;
    int i;
    for (i=0; i < Z->nsamples; i++)
        if(iftHasSampleStatus(Z->sample[i], IFT_TRAIN))
            ntrain++;

    Z->ntrainsamples = ntrain;
}

int iftSelectSupTrainSamples(iftDataSet *Z, float train_perc)
{

    iftDeprecated("iftSelectSupTrainSamples", "iftDataSetSampling/iftStratifiedRandomSubsampling", "Use iftDataSetSampling with a iftStratifiedRandomSubsampling strategy.");

    int s, *sample=NULL, *count=NULL, i;
    int t, high, nsamples, c, *class_size=NULL;

    if (Z->nclasses == 0)
        iftError("There are no classes", "iftSelectSupTrainSamples");

    if ((train_perc <= 0.0)||(train_perc > 1.0))
        iftError("Invalid percentage of training samples", "iftSelectSupTrainSamples");

    // Reset status

    iftSetStatus(Z, IFT_TEST);

    // Verify if it is the trivial case of selecting all samples.

    if (iftAlmostZero(train_perc - 1.0)) {
        for (s=0; s < Z->nsamples; s++)
            iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
        Z->ntrainsamples = Z->nsamples;
        return(Z->nsamples);
    }

    // Compute the number of training samples

    Z->ntrainsamples = (int) (train_perc*Z->nsamples);
    if (Z->ntrainsamples == 0)
        iftError("Percentage is too low for this dataset", "iftSelectSupTrainSamples");


    // Count number of samples per class


    class_size = iftAllocIntArray(Z->nclasses+1);
    for (s=0; s < Z->nsamples; s++){
        class_size[Z->sample[s].truelabel]++;
    }

    // Verify percentage and number of training samples per class

    Z->ntrainsamples = 0;
    for (c=1; c <= Z->nclasses; c++) {
        nsamples = (int)(train_perc*class_size[c]);
        if (nsamples > class_size[c] || nsamples == 0)
            nsamples = class_size[c];
        Z->ntrainsamples += nsamples;
        if (nsamples <= 0){
            fprintf(stderr,"For class %d\n",c);
            iftError("No available samples", "iftSelectSupTrainSamples");
        }
    }

    // Randomly select samples


    for (c=1; c <= Z->nclasses; c++) {
        nsamples = (int)(train_perc*class_size[c]);
        if (nsamples > class_size[c])
            nsamples = class_size[c];
        sample = iftAllocIntArray(class_size[c]);
        count  = iftAllocIntArray(class_size[c]);
        t=0;
        for (s=0; s < Z->nsamples; s++)
            if (Z->sample[s].truelabel==c){
                sample[t]=s;
                count[t]=100;
                t++;
            }
        t = 0; high = class_size[c]-1;
        while (t < nsamples) {
            i = iftRandomInteger(0,high);
            s = sample[i];
            if (count[i]==0){
                iftSetSampleStatus(&Z->sample[s], IFT_TRAIN);
                iftSwap(sample[i], sample[high]);
                iftSwap(count[i], count[high]);
                t++; high--;
            }else{
                count[i]--;
            }
        }
        iftFree(count);
        iftFree(sample);
    }

    iftFree(class_size);

    return(Z->ntrainsamples);
}


iftDataSet *iftReadXYDataSet(char *filename)
{
    FILE *fp=fopen(filename,"r");
    int   s,nsamples, x, y;
    iftDataSet *Z;

    if (fscanf(fp,"%d ",&nsamples)!=1) iftError("Reading error", "iftReadXYDataSet");
    Z = iftCreateDataSet(nsamples,2);
    Z->nclasses = 0;
    for (s=0; s < Z->nsamples; s++){
        if (fscanf(fp,"%d %d %d",&x,&y,&Z->sample[s].truelabel)!=3) iftError("Reading error", "iftReadXYDataSet");
        if (Z->sample[s].truelabel > Z->nclasses)
            Z->nclasses = Z->sample[s].truelabel;
        Z->sample[s].feat[0] = x;
        Z->sample[s].feat[1] = y;
    }

    fclose(fp);
    return(Z);
}

/**
* Transforms the train dataset in a matrix with ones in the last column as an independent term.
*/
iftMatrix*iftDataSetToFeatureMatrixHomogCoord(iftDataSet *Z, char status) {

    int nsamples = 0;
    for (int i = 0; i < Z->nsamples; ++i) {
        if(iftHasSampleStatus(Z->sample[i], status))
            nsamples++;
    }

    iftMatrix* m = iftCreateMatrix(Z->nfeats+1, nsamples);
    int idx = 0;
    for (int i = 0; i < Z->nsamples; ++i) {
        if(!iftHasSampleStatus(Z->sample[i], status)) continue;
        for (int j = 0; j < Z->nfeats; ++j) {
            m->val[iftGetMatrixIndex(m, j, idx)] = Z->sample[i].feat[j];
        }
        m->val[iftGetMatrixIndex(m, Z->nfeats, idx)] = 1;
        idx++;
    }
    return m;
}

iftMatrix* iftDataSetToLabelsMatrix(iftDataSet *Z, char status) {

    int nsamples = 0;
    for (int i = 0; i < Z->nsamples; ++i) {
        if(iftHasSampleStatus(Z->sample[i], status))
            nsamples++;
    }

    int idx=0;
    iftMatrix* m = iftCreateMatrix(Z->nclasses, nsamples);
    for (int i = 0; i < Z->nsamples; ++i) {
        if (!iftHasSampleStatus(Z->sample[i], status)) continue;
        for (int c = 0; c < Z->nclasses; ++c) {
            m->val[iftGetMatrixIndex(m, c, idx)] = (Z->sample[i].truelabel == c+1);
        }
        idx++;
    }
    return m;
}

iftMatrix* iftDataSetMeanFeatures(const iftDataSet* Z) {
    iftMatrix* mean = iftCreateMatrix(Z->nfeats, 1);

    for (int i = 0; i < Z->nsamples; ++i) {
        for (int j = 0; j < Z->nfeats; ++j) {
            mean->val[j] += Z->sample[i].feat[j]/(float)Z->nsamples;
        }
    }

    return mean;
}

iftMatrix* iftDataSetStdDevFeatures(const iftDataSet* Z) {
    iftMatrix* mean = iftDataSetMeanFeatures(Z);
    iftMatrix* std = iftCreateMatrix(Z->nfeats, 1);

    double *temp = iftAllocDoubleArray(Z->nfeats);

    for (int i = 0; i < Z->nsamples; ++i) {
        for (int j = 0; j < Z->nfeats; ++j) {
            temp[j] += iftFastPow((double)(Z->sample[i].feat[j] - mean->val[j]), 2)/(double)Z->nsamples;
        }
    }

    for (int j = 0; j < Z->nfeats; ++j) {
        std->val[j] += (float) sqrt(temp[j]);
    }

    iftDestroyMatrix(&mean);
    iftFree(temp);

    return std;
}

iftMatrix* iftMatrixMultTranspose(const iftMatrix* M) {
    iftMatrix* T = iftTransposeMatrix(M);

    iftMatrix* mult = iftMultMatrices(M, T);
    iftDestroyMatrix(&T);

    return mult;
}

iftMatrix* iftTransposeMultMatrix(const iftMatrix* M) {
    iftMatrix* T = iftTransposeMatrix(M);

    iftMatrix* mult =  iftMultMatrices(T, M);
    iftDestroyMatrix(&T);

    return mult;
}

iftMatrix* iftLinearDiscriminantAnalysis(const iftDataSet *Z) {

    int nclasses = Z->nclasses;
    int nfeats = Z->nfeats;

    if(!iftIsSupervisedDataSet(Z)) {
        iftError("Linear Discriminant Analysis is a supervised method.", "iftLinearDiscriminantAnalysis");
    }

    if(nfeats<nclasses) {
        iftError("Number of features must be greater than the number of classes.", "iftLinearDiscriminantAnalysis");
    }

    float fact = 1.0/(float)(Z->nsamples - Z->nclasses);
    float sqrt_fact = sqrtf(fact);

    iftDataSet* Zt = iftExtractSamples(Z, IFT_TRAIN);
    iftMatrix* globalMean = iftDataSetMeanFeatures(Zt);

    iftMatrix* classesMean = iftCreateMatrix(nfeats, nclasses);
    iftIntArray* nsamplesPerClass = iftCreateIntArray(nclasses);

    for (int c = 0; c < Z->nclasses; ++c) {
        iftDataSet* Zc = iftExtractSamplesFromClass(Zt, c+1);
        iftMatrix* mean = iftDataSetMeanFeatures(Zc);

        nsamplesPerClass->val[c] = Zc->nsamples;

        for (int j = 0; j < nfeats; ++j) {
            iftMatrixElem(classesMean, j, c) = mean->val[j];
        }
        iftDestroyMatrix(&mean);
        iftDestroyDataSet(&Zc);
    }

    for(int s=0; s<Zt->nsamples; s++) {
        int c = Zt->sample[s].truelabel - 1;
        for (int j = 0; j < Zt->nfeats; ++j) {
            Zt->sample[s].feat[j] -= iftMatrixElem(classesMean, j, c);
        }
    }

    iftMatrix* std = iftDataSetStdDevFeatures(Zt);

    for (int i = 0; i < std->n; ++i) {
        if(iftAlmostZero(std->val[i]))
            std->val[i] = 1.0;
    }

    for (int s = 0; s < Zt->nsamples; ++s) {
        for (int j = 0; j < nfeats; ++j) {
            Zt->sample[s].feat[j] = (sqrt_fact * Zt->sample[s].feat[j]) / (std->val[j]);
        }
    }

    iftMatrix* X = iftDataSetToFeatureMatrix(Zt);

    //within variance
    iftMatrix *U, *S, *V;
    iftSingleValueDecomp(X, &U, &S, &V);

    iftDestroyMatrix(&X);

    iftMatrix* scaledFeats = iftCreateMatrix(V->ncols, V->nrows);

    fflush(stdout);

    //Scaling of within covariance is: V' 1/S
    for (int i = 0; i < V->nrows; ++i) {
        if(S->val[i] < 1e-12)
            break;
        for (int j = 0; j < V->ncols; ++j) {
            iftMatrixElem(scaledFeats, j, i) = (iftMatrixElem(V, j, i)/ std->val[j]) / S->val[i];
        }
    }

    iftDestroyMatrix(&U);
    iftDestroyMatrix(&S);
    iftDestroyMatrix(&V);

    iftTransposeMatrixInPlace(scaledFeats);
    iftMatrix* classCenteredMean = iftCreateMatrix(Zt->nfeats, Zt->nclasses);

    for (int c = 0; c < nclasses; ++c) {
        for (int j = 0; j < nfeats; ++j) {
            iftMatrixElem(classCenteredMean, j, c) = sqrtf(fact * nsamplesPerClass->val[c]) * (iftMatrixElem(classesMean, j, c) - globalMean->val[j]);
        }
    }

    iftMatrix* SwInvSb = iftMultMatrices(classCenteredMean, scaledFeats);

    iftSingleValueDecomp(SwInvSb, &U, &S, &V);

    iftTransposeMatrixInPlace(V);
    iftMatrix* scaling = iftMultMatrices(scaledFeats, V);

    iftMatrix* T = iftSubMatrix(scaling, 0, scaling->nrows, 0, Zt->nclasses-1);

    iftDestroyMatrix(&scaling);
    iftDestroyMatrix(&U);
    iftDestroyMatrix(&V);
    iftDestroyMatrix(&S);
    iftDestroyMatrix(&SwInvSb);
    iftDestroyMatrix(&scaledFeats);
    iftDestroyMatrix(&scaling);
    iftDestroyMatrix(&globalMean);
    iftDestroyMatrix(&classesMean);
    iftDestroyMatrix(&X);
    iftDestroyDataSet(&Zt);
    iftDestroyMatrix(&std);
    iftDestroyMatrix(&classCenteredMean);

    return T;
}

void iftWriteCSVDataSet(iftDataSet* dataset, const char* filename, char *fileset) {

  iftCSV* csv;

  if ((fileset!=NULL)&& (dataset->ref_data_type==IFT_REF_DATA_FILESET))
    csv = iftCreateCSV(dataset->nsamples, dataset->nfeats + 2);
  else
    csv = iftCreateCSV(dataset->nsamples, dataset->nfeats + 1);
    
  char buffer[IFT_STR_DEFAULT_SIZE];

    for (int s = 0; s < dataset->nsamples; ++s) {

      for (int j = 0; j < dataset->nfeats; ++j) {
	sprintf(buffer, "%f", dataset->sample[s].feat[j]);
	strcpy(csv->data[s][j], iftCopyString(buffer));
      }

      sprintf(buffer, "%d", dataset->sample[s].truelabel);
      strcpy(csv->data[s][dataset->nfeats], iftCopyString(buffer));

    }
    
    if((fileset != NULL) && (dataset->ref_data_type==IFT_REF_DATA_FILESET)){
      iftWriteFileSetAsCSV((iftFileSet*) dataset->ref_data, fileset);
      for (int s = 0; s < dataset->nsamples; ++s) {
	sprintf(buffer, "%d", dataset->sample[s].id);
	strcpy(csv->data[s][dataset->nfeats+1], iftCopyString(buffer));
      }
    }

    iftWriteCSV(csv, filename, ',');
    iftDestroyCSV(&csv);
}


iftDataSet *iftReadOPFDataSet(const char *pathname) {
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftReadOPFDataSet");
    if (iftEndsWith(pathname, ".data")) {
        iftError("Dataset \"%s\" has old format\n" \
               "Convert it using the program demo/Classification/iftConvertOldOPFDataSet.c",
                 "iftReadOPFDataSet", pathname);
    }

    char *tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);
    iftUnzipFile(pathname, tmp_dir);

    //////////////// READING DATASET HEADER AND SAMPLE INFORMATION ///////////////////
    char *info_path = iftJoinPathnames(2, tmp_dir, "info.data");
    FILE *fp = fopen(info_path, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadOPFDataSet", info_path);
    iftFree(info_path);


    int nsamples, nfeats, nclasses, ngroups;
    int function_number = IFT_NIL;
    if (fread(&nsamples, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when Reading Number of Samples", "iftReadOPFDataSet");
    }
    if (fread(&nfeats, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when Reading Number of Classes", "iftReadOPFDataSet");
    }
    if (fread(&nclasses, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when Reading Number of Classes", "iftReadOPFDataSet");
    }
    if (fread(&ngroups, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when Reading Number of Labels", "iftReadOPFDataSet");
    }
    if (fread(&function_number, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when Reading Function Number", "iftReadOPFDataSet");
    }


    iftDataSet *Z      = iftCreateDataSet(nsamples, nfeats);
    Z->nclasses        = nclasses;
    Z->ngroups         = ngroups;
    Z->function_number = function_number;
    Z->ntrainsamples   = 0;
    iftSetDistanceFunction(Z, function_number);

    //// READING SAMPLES
    for (int s = 0; s < Z->nsamples; s++) {

        if (fread(&Z->sample[s].id, sizeof(int), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when reading Id from Sample %d", "iftReadOPFDataSet", s);
        }
        if (fread(&Z->sample[s].truelabel, sizeof(int), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when reading True Label (class) from Sample %d", "iftReadOPFDataSet", s);
        }
        if (fread(&Z->sample[s].label, sizeof(int), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when reading Predicted Label from Sample %d", "iftReadOPFDataSet", s);
        }
        if (fread(&Z->sample[s].weight, sizeof(float), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when reading Weight from Sample %d", "iftReadOPFDataSet", s);
        }
        if (fread(&Z->sample[s].status, sizeof(uchar), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when reading Status from Sample %d", "iftReadOPFDataSet", s);
        }
        if (fread(Z->sample[s].feat, sizeof(float), Z->nfeats, fp) != Z->nfeats) {
            iftRemoveDir(tmp_dir);
            iftError("Error when reading Feat. Vector from Sample %d", "iftReadOPFDataSet", s);
        }

        if(!iftHasSampleStatus(Z->sample[s], IFT_ALL))
            iftSetSampleStatus(&Z->sample[s], IFT_UNKNOWN);

        Z->ntrainsamples += iftHasSampleStatus(Z->sample[s], IFT_TRAIN);
    }

    if (fread(Z->alpha, sizeof(float), Z->nfeats, fp) != Z->nfeats) {
        iftRemoveDir(tmp_dir);
        iftError("Error when reading Alpha Array", "iftReadOPFDataSet");
    }

    if (fread(&Z->ref_data_type, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when reading the Ref. Data Type", "iftWriteOPFDataSet");
    }

    fclose(fp);
    /////////////////////////////////////////////////////////////////////////////

    //////////////// READING FEAT. SPACE PARAM ///////////////////
    char *fsp_path = iftJoinPathnames(2, tmp_dir, "feat_space.fsp");
    Z->fsp = iftReadFeatSpaceParam(fsp_path);
    iftFree(fsp_path);
    /////////////////////////////////////////////////////////////

    //////////////// WRITING REF. DATA ///////////////////
    char *ref_data_path = NULL;
    if (Z->ref_data_type == IFT_REF_DATA_IMAGE) {
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.scn");
        if (!iftFileExists(ref_data_path)) { // tests if a 3D image exists
            iftFree(ref_data_path);
            ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.png");

            if (!iftFileExists(ref_data_path)) { // tests if a 2D image exists
                iftError("There is no Reference Data Image with extension *.[scn,png] in Dataset",
                         "iftWriteOPFDataset");
            }
        }

        Z->ref_data = iftReadImageByExt(ref_data_path);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE) {
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.npy");
        if (!iftFileExists(ref_data_path))
            iftError("There is no Reference Data FImage with extension *.fscn in Dataset",
                     "iftWriteOPFDataset");

        Z->ref_data = iftReadFImage(ref_data_path);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE) {
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.mimg");
        if (!iftFileExists(ref_data_path))
            iftError("There is no Reference Data MImage with extension *.mimg in Dataset",
                     "iftWriteOPFDataset");

        Z->ref_data = iftReadMImage(ref_data_path);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_FILESET) {
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.csv");
        if (!iftFileExists(ref_data_path))
            iftError("There is no Reference Data FILE \"%s\" in Dataset",
                     "iftWriteOPFDataset", ref_data_path);

        Z->ref_data = iftLoadFileSetFromCSV(ref_data_path, false);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_MMKERNEL) {
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.mmk");
        if (!iftFileExists(ref_data_path))
            iftError("There is no Reference Data MMKernel with extension *.mmk in Dataset",
                     "iftWriteOPFDataset");

        Z->ref_data = iftReadMMKernel(ref_data_path);
    }

    iftFree(ref_data_path);
    /////////////////////////////////////////////////////////////


    // DESTROYERS
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);

    return Z;
}


void iftWriteOPFDataSet(const iftDataSet *Z, const char *pathname) {

  if (Z == NULL)
      iftError("Dataset is NULL", "iftWriteOPFDataSet");

    if (pathname == NULL)
        iftError("Pathname is NULL", "iftWriteOPFDataSet");

    if (!iftEndsWith(pathname, ".zip"))
        iftError("Extension from Dataset's pathname \"%s\" should be *.zip",
                 "iftWriteOPFDataSet", pathname);
    if (Z->function_number <= 0)
        iftError("Function number is %d, which cannot be identified as any standard distance function for OPF.\n" \
               "That most likely indicates that a custom distance function is being used. " \
               "Please set the function number to a standard one before saving the OPF dataset to a file, "\
                "and then reload the proper function manually when necessary.",
                 "iftWriteOPFDataSetFilePointer", Z->function_number);

    char *tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);


    //////////////// WRITING DATASET HEADER AND SAMPLE INFORMATION ///////////////////
    char *info_path = iftJoinPathnames(2, tmp_dir, "info.data");
    FILE *fp = fopen(info_path, "wb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteOPFDataSet", info_path);
    iftFree(info_path);

    if (fwrite(&Z->nsamples, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when writing Number of Samples", "iftWriteOPFDataSet");
    }
    if (fwrite(&Z->nfeats, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when writing Number of Features", "iftWriteOPFDataSet");
    }
    if (fwrite(&Z->nclasses, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when writing Number of Classes", "iftWriteOPFDataSet");
    }
    if (fwrite(&Z->ngroups, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when writing Number of Groups", "iftWriteOPFDataSet");
    }
    if (fwrite(&Z->function_number, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when writing Function Number", "iftWriteOPFDataSet");
    }

    //// WRITING SAMPLES
    for (int s = 0; s < Z->nsamples; s++) {
        if (fwrite(&Z->sample[s].id, sizeof(int), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when writing Id from Sample %d", "iftWriteOPFDataSet", s);
        }
        if (fwrite(&Z->sample[s].truelabel, sizeof(int), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when writing True Label (class) from Sample %d", "iftWriteOPFDataSet", s);
        }
        if (fwrite(&Z->sample[s].label, sizeof(int), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when writing Predicted Label from Sample %d", "iftWriteOPFDataSet", s);
        }
        if (fwrite(&Z->sample[s].weight, sizeof(float), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when writing Weight from Sample %d", "iftWriteOPFDataSet", s);
        }
        if (fwrite(&Z->sample[s].status, sizeof(uchar), 1, fp) != 1) {
            iftRemoveDir(tmp_dir);
            iftError("Error when writing Status from Sample %d", "iftWriteOPFDataSet", s);
        }
        if (fwrite(Z->sample[s].feat, sizeof(float), Z->nfeats, fp) != Z->nfeats) {
            iftRemoveDir(tmp_dir);
            iftError("Error when writing Feat. Vector from Sample %d", "iftWriteOPFDataSet", s);
        }
    }

    if (fwrite(Z->alpha, sizeof(float), Z->nfeats, fp) != Z->nfeats) {
        iftRemoveDir(tmp_dir);
        iftError("Error when writing Alpha Array", "iftWriteOPFDataSet");
    }

    // to prevent that the ref_data is NULL and its type is != IFT_REF_DATA_NULL
    iftRefDataType ref_data_type = IFT_REF_DATA_NULL;
    if (Z->ref_data != NULL)
        ref_data_type = Z->ref_data_type;

    if (fwrite(&ref_data_type, sizeof(int), 1, fp) != 1) {
        iftRemoveDir(tmp_dir);
        iftError("Error when writing the Ref. Data Type", "iftWriteOPFDataSet");
    }
    fclose(fp);
    /////////////////////////////////////////////////////////////////////////////

    //////////////// WRITING FEAT. SPACE PARAM ///////////////////
    char *fsp_path = iftJoinPathnames(2, tmp_dir, "feat_space.fsp");
    iftWriteFeatSpaceParam(&Z->fsp, fsp_path);
    iftFree(fsp_path);
    /////////////////////////////////////////////////////////////

    //////////////// WRITING REF. DATA ///////////////////
    char *ref_data_path = NULL;
    if (Z->ref_data_type == IFT_REF_DATA_IMAGE) {
        iftImage *img = (iftImage*) Z->ref_data;

        if (iftIs3DImage(img)) {
            ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.scn");
        }
        else { // 2D
            ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.png"); // saves gray and color images
        }
        iftWriteImageByExt(img, ref_data_path);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE) {
        iftFImage *fimg = (iftFImage*) Z->ref_data;
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.npy");
        iftWriteFImage(fimg, ref_data_path);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE) {
        iftMImage *mimg = (iftMImage*) Z->ref_data;
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.mimg");
        iftWriteMImage(mimg, ref_data_path);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_FILESET) {
        iftFileSet *fset = (iftFileSet*) Z->ref_data;
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.csv");
        iftCSV *csv   = iftFileSetToCSV(fset);
        iftWriteCSV(csv, ref_data_path, ',');
        iftDestroyCSV(&csv);
    } else if (Z->ref_data_type == IFT_REF_DATA_MMKERNEL) {
        iftMMKernel *K = (iftMMKernel *) Z->ref_data;
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.mmk");
        iftWriteMMKernel(K,ref_data_path);
    } else if (Z->ref_data_type == IFT_REF_DATA_CSV) {
        iftCSV *csv = (iftCSV *) Z->ref_data;
        ref_data_path = iftJoinPathnames(2, tmp_dir, "ref_data.csv");
        iftWriteCSV(csv, ref_data_path, ',');
    }

    iftFree(ref_data_path);
    /////////////////////////////////////////////////////////////

    // ZIP ALL FILES
    iftZipDirContent(tmp_dir, pathname);


    // DESTROYERS
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
}



void iftSetRefData(iftDataSet *Z, const void *ref_data, iftRefDataType ref_data_type) {
    if (Z->ref_data != NULL) {
        iftWarning("Deallocating an existent reference data", "iftSetRefData");
        iftDestroyRefData(Z);
    }

    Z->ref_data_type = ref_data_type;
    
    if (Z->ref_data_type == IFT_REF_DATA_NULL) {
        return;
    } 
    if (Z->ref_data_type == IFT_REF_DATA_IMAGE) {
        iftImage *img = (iftImage*) ref_data;
        if (img->n != Z->nsamples)
            iftError("Number of Pixels/Voxels from the Ref. Image is != from the Number of Samples: %d != %d",
                     "iftSetRefData", img->n, Z->nsamples);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE) {
        iftFImage *fimg = (iftFImage*) ref_data;
        if (fimg->n != Z->nsamples)
            iftError("Number of Pixels/Voxels from the Ref. FImage is != from the Number of Samples: %d != %d",
                     "iftSetRefData", fimg->n, Z->nsamples);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE) {
        iftMImage *mimg = (iftMImage*) ref_data;
        if (mimg->n != Z->nsamples)
            iftError("Number of Pixels/Voxels from the Ref. MImage is != from the Number of Samples: %d != %d",
                     "iftSetRefData", mimg->n, Z->nsamples);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_FILESET) {
        iftFileSet *fset = (iftFileSet*) ref_data;
        if (fset->n != Z->nsamples)
            iftError("Number of Files from the File Set is != from the Number of Samples: %d != %d",
                     "iftSetRefData", fset->n, Z->nsamples);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_MMKERNEL) {
        iftMMKernel *K = (iftMMKernel*) ref_data;
        if (K->nkernels != Z->nsamples)
            iftError("Number of kernels from the MMKernel is != from the Number of Samples: %d != %d",
                     "iftSetRefData", K->nkernels, Z->nsamples);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_CSV) {
        iftCSV *csv = (iftCSV*) ref_data;
        if (csv->nrows != Z->nsamples)
            iftError("Number of rows from the CSV is != from the Number of Samples: %d != %d",
                     "iftSetRefData", csv->nrows, Z->nsamples);
    }
    else iftError("Setting this Ref. Data Type not implemented yet", "iftSetRefData");


    // sets the ids from 0 to Z->nsamples-1 --> sample[s].id = s
    iftSetDataSetIDsincrementally(Z);
    iftCopyRefData(Z, ref_data, ref_data_type);
}



void iftCopyRefData(iftDataSet *Z, const void *ref_data, iftRefDataType ref_data_type) {
  
    if (Z->ref_data != NULL) {
        iftWarning("Deallocating an existent reference data", "iftCopyRefData");
        iftDestroyRefData(Z);
    }

    if(ref_data != NULL) {
        Z->ref_data_type = ref_data_type;
        
        if (Z->ref_data_type == IFT_REF_DATA_IMAGE) {
            iftImage *img = (iftImage*) ref_data;
            Z->ref_data = iftCopyImage(img);
        }
        else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE) {
            iftFImage *fimg = (iftFImage*) ref_data;
            Z->ref_data = iftFCopyImage(fimg);
        }
        else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE) {
            iftMImage *mimg = (iftMImage*) ref_data;
            Z->ref_data = iftCopyMImage(mimg);
        }
        else if (Z->ref_data_type == IFT_REF_DATA_FILESET) {
            iftFileSet *fset = (iftFileSet*) ref_data;
            Z->ref_data = iftCopyFileSet(fset);
        }
        else if (Z->ref_data_type == IFT_REF_DATA_MMKERNEL) {
            iftMMKernel *K = (iftMMKernel*) ref_data;
            Z->ref_data = iftCopyMMKernel(K);
        }
        else if (Z->ref_data_type == IFT_REF_DATA_CSV) {
            iftCSV *csv = (iftCSV*) ref_data;
            Z->ref_data = iftCopyCSV(csv);
        }
        else iftError("Setting this Ref. Data Type not implemented yet", "iftCopyRefData");
    }
}

void iftWriteDataSet(const iftDataSet *dataset, const char *format, ...) {
    va_list args;
    char pathname[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(pathname, format, args);
    va_end(args);
    
    char *parent_dir = iftParentDir(pathname);

    if (!iftDirExists(parent_dir))
        iftMakeDir(parent_dir);
    iftFree(parent_dir);

    iftWriteDatasetGivenVersion(dataset,pathname,IFT_CURRENT_DATASET_VERSION);

}

void iftWriteDatasetGivenVersion(const iftDataSet *dataSet, const char *pathname, iftDatasetVersion version){

  if(version == OPFDATASET_0){
    iftWriteOPFDataSet(dataSet,pathname);
  } else{
    iftWriteDataSetFishStandard(dataSet,pathname,version);
  }
}

bool iftCompareSamples(const iftSample* sample1, const iftSample* sample2, int nfeats) {

    if(sample1->id!=sample2->id) {
        return false;
    }
    if(sample1->truelabel!=sample2->truelabel) {
        return false;
    }
    if(sample1->label!=sample2->label) {
        return false;
    }
    if(sample1->status!=sample2->status) {
        return false;
    }
    if(sample1->weight!=sample2->weight) {
        return false;
    }
    if(sample1->group!=sample2->weight) {
        return false;
    }
    if(sample1->isLabelPropagated!=sample2->isLabelPropagated) {
        return false;
    }
    if(sample1->isSupervised!=sample2->isSupervised) {
        return false;
    }
    if(sample1->numberTimesChecked!=sample2->numberTimesChecked) {
        return false;
    }

    for(int f = 0; f < nfeats; ++f) {
        float x = fabs((sample1->feat[f])-(sample2->feat[f]));
        printf("%f %f --> %f\n", (sample1->feat[f]), (sample2->feat[f]), x);
        if(x > IFT_EPSILON)
            return false;
    }

    return true;
}

bool iftCompareDataSets(const iftDataSet* dataset1, const iftDataSet* dataset2, bool compareFeatures, bool verbose) {

    if(dataset1->nsamples!=dataset2->nsamples) {
        if(verbose) {
            printf("Different number of samples => (%d,%d).\n", dataset1->nsamples, dataset2->nsamples);
        }
        return false;
    }
    if(dataset1->ntrainsamples!=dataset2->ntrainsamples) {
        if(verbose) {
            printf("Different number of training samples => (%d,%d).\n", dataset1->ntrainsamples, dataset2->ntrainsamples);
        }
        return false;
    }
    if(dataset1->nclasses!=dataset2->nclasses) {
        if(verbose) {
            printf("Different number of classes => (%d,%d).\n", dataset1->nclasses, dataset2->nclasses);
        }
        return false;
    }

    if(dataset1->function_number!=dataset2->function_number) {
        if(verbose) {
            printf("Different distance function => (%d,%d).\n", dataset1->function_number, dataset2->function_number);
        }
        return false;
    }

    for (int s = 0; s < dataset1->nsamples; ++s) {
        if(!iftCompareSamples(&(dataset1->sample[s]), &(dataset2->sample[s]), dataset1->nfeats)) {
            if(verbose) {
                printf("Difference in sample #%d.\n", s);
            }
            return false;
        }
    }

    return true;
}

void iftWriteDataSetFishStandard(const iftDataSet *dataSet, const char *pathname, iftDatasetVersion version){
    if (dataSet == NULL) {
        iftError("Dataset is NULL", "iftWriteDataSet");
    }
    if (pathname == NULL) {
        iftError("Pathname is NULL", "iftWriteDataSet");
    }

    char *tempDirPath = iftMakeTempDir("tmpdir_", NULL, NULL);
    char *buffer = NULL;
    FILE *fp = NULL;

    // struct stat st = {0};
    // if (stat(tempDirPath, &st) == -1) {
    //     mkdir(tempDirPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // }
    //////////////// WRITING DATASET HEADER AND SAMPLE INFORMATION ///////////////////
    int* dimensionsFeat = iftAlloc(2,sizeof(int));
    int* dimensionsProjection = iftAlloc(2,sizeof(int));
    int* dimensionsSampleInfo = iftAlloc(2,sizeof(int));
    int* dimensionsGeneralInfo = iftAlloc(2,sizeof(int));
    int* dimensionsAlphas = iftAlloc(2,sizeof(int));
    dimensionsFeat[0] = dataSet->nsamples;
    dimensionsFeat[1] = dataSet->nfeats;
    size_t totalElementsFeat = (size_t)dimensionsFeat[0]*(size_t)dimensionsFeat[1];
    dimensionsProjection[0] = dataSet->nsamples;
    dimensionsProjection[1] = 2;
    size_t totalElementsProjection = dimensionsProjection[0]*dimensionsProjection[1];
    dimensionsSampleInfo[0] = dataSet->nsamples;
    dimensionsSampleInfo[1] = 1;
    size_t totalElementsSampleInfo = dimensionsSampleInfo[0]*dimensionsSampleInfo[1];
    dimensionsGeneralInfo[0] =1;
    dimensionsGeneralInfo[1] =1;
    size_t totalElementsGeneralInfo = dimensionsGeneralInfo[0]*dimensionsGeneralInfo[1];
    dimensionsAlphas[0] = dataSet->nfeats;
    dimensionsAlphas[1] = 1;

    //Transfering sample info to arrays
    int* truelabels = NULL;
    int* labels = NULL;
    int *id = NULL;
    int *numberTimesChecked = NULL;
    int *cluster_id = NULL;
    uchar* status = NULL;
    uchar* isSupervised = NULL;
    uchar* isLabelPropagated = NULL;
    float* weight = NULL;
//    float* alpha = iftAlloc(dimensionsFeat[1],sizeof(float));
    if(dataSet->sample != NULL) {
        truelabels = iftAlloc(dimensionsSampleInfo[0],sizeof(int));
        labels = iftAlloc(dimensionsSampleInfo[0],sizeof(int));
        id = iftAlloc(dimensionsSampleInfo[0],sizeof(int));
        cluster_id = iftAlloc(dimensionsSampleInfo[0],sizeof(int));
        status = iftAlloc(dimensionsSampleInfo[0],sizeof(uchar));
        isSupervised = iftAlloc(dimensionsSampleInfo[0],sizeof(uchar));
        isLabelPropagated = iftAlloc(dimensionsSampleInfo[0],sizeof(uchar));
        weight = iftAlloc(dimensionsSampleInfo[0],sizeof(float));
        numberTimesChecked = iftAlloc(dimensionsSampleInfo[0],sizeof(int));
        for (int sampleIndex = 0; sampleIndex < dataSet->nsamples; ++sampleIndex) {
            truelabels[sampleIndex] = dataSet->sample[sampleIndex].truelabel;
            labels[sampleIndex] = dataSet->sample[sampleIndex].label;
            id[sampleIndex] = dataSet->sample[sampleIndex].id;
            cluster_id[sampleIndex] = dataSet->sample[sampleIndex].group;
            status[sampleIndex] = dataSet->sample[sampleIndex].status;
            isSupervised[sampleIndex] = dataSet->sample[sampleIndex].isSupervised;
            isLabelPropagated[sampleIndex] = dataSet->sample[sampleIndex].isLabelPropagated;
            weight[sampleIndex] = dataSet->sample[sampleIndex].weight;
            numberTimesChecked[sampleIndex] = dataSet->sample[sampleIndex].numberTimesChecked;
        }
    }
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"version.data");

    buffer = iftJoinPathnames(2, tempDirPath, "version.data");
    int ver = version;
    fp = fopen(buffer, "wb");
    if (fp == NULL){
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
    }
    writeGenericVectorPointerAsText(int, "%d", ',', '\n', (&ver), 1, fp);
    fclose(fp);
    iftFree(buffer);
    ////
    //writing feats
    if(dataSet->data != NULL) {
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"feat.data");
        buffer = iftJoinPathnames(2, tempDirPath, "feat.data");
        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',', '\n', dimensionsFeat, 2, fp);
        writeGenericVectorPointerAsBinary(float, dataSet->data->val, totalElementsFeat, fp);
        fclose(fp);
        iftFree(buffer);
    }else if(dataSet->sample != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"feat.data");
        buffer = iftJoinPathnames(2, tempDirPath, "feat.data");
        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        float *feats = (float*)iftAlloc(dataSet->nsamples*dataSet->nfeats,sizeof(float));
        long k=0;
        for (int i = 0; i < dataSet->nsamples; ++i) {
            for (int j = 0; j < dataSet->nfeats; ++j) {
                feats[k] = dataSet->sample[i].feat[j];
                k++;
            }
        }
        writeGenericVectorPointerAsText(int, "%d", ',', '\n', dimensionsFeat, 2, fp);
        writeGenericVectorPointerAsBinary(float, feats, totalElementsFeat, fp);

        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing projection
    if(dataSet->projection != NULL) {
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"projection.data");
        buffer = iftJoinPathnames(2, tempDirPath, "projection.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',', '\n', dimensionsProjection, 2, fp);
        writeGenericVectorPointerAsBinary(double, dataSet->projection->val, totalElementsProjection, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing alpha
    if(dataSet->alpha != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"alpha.data");
        buffer = iftJoinPathnames(2, tempDirPath, "alpha.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',', '\n', dimensionsAlphas, 2, fp);
        writeGenericVectorPointerAsBinary(float, dataSet->alpha, dataSet->nfeats, fp);
        fclose(fp);
        iftFree(buffer);
    } else {
      iftWarning("This is a bug","iftWriteDataSet");
    }
    ////
    //writing truelabels
    if(truelabels != NULL) {
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"true_label.data");
        buffer = iftJoinPathnames(2, tempDirPath, "true_label.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL) {
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',', '\n', dimensionsSampleInfo, 2, fp);
        writeGenericVectorPointerAsBinary(int, truelabels, totalElementsSampleInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing labels
    if(labels != NULL) {
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"label.data");
        buffer = iftJoinPathnames(2, tempDirPath, "label.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsSampleInfo,2,fp);
        writeGenericVectorPointerAsBinary(int, labels, totalElementsSampleInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing ids
    if(id != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"id.data");
        buffer = iftJoinPathnames(2, tempDirPath, "id.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsSampleInfo,2,fp);
        writeGenericVectorPointerAsBinary(int, id, totalElementsSampleInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing groups
    if(cluster_id != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"group.data");
        buffer = iftJoinPathnames(2, tempDirPath, "group.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsSampleInfo,2,fp);
        writeGenericVectorPointerAsBinary(int, cluster_id, totalElementsSampleInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing status
    if(status != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"status.data");
        buffer = iftJoinPathnames(2, tempDirPath, "status.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsSampleInfo,2,fp);
        writeGenericVectorPointerAsBinary(uchar, status, totalElementsSampleInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////

    //writing numberOfChecked
    if(status != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"n_checked.data");
        buffer = iftJoinPathnames(2, tempDirPath, "n_checked.data");
        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsSampleInfo,2,fp);
        writeGenericVectorPointerAsBinary(int, numberTimesChecked, totalElementsSampleInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////

    //writing isSupervised
    if(isSupervised != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"is_supervised.data");
        buffer = iftJoinPathnames(2, tempDirPath, "is_supervised.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsSampleInfo,2,fp);
        writeGenericVectorPointerAsBinary(uchar, isSupervised, totalElementsSampleInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }

    ////
    //writing isLabelPropagated
    if(isLabelPropagated != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"is_label_propagated.data");
        buffer = iftJoinPathnames(2, tempDirPath, "is_label_propagated.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsSampleInfo,2,fp);
        writeGenericVectorPointerAsBinary(uchar, isLabelPropagated, totalElementsSampleInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }

    ////
    //writing weight
    if(weight != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"weight.data");
        buffer = iftJoinPathnames(2, tempDirPath, "weight.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsSampleInfo,2,fp);
        writeGenericVectorPointerAsBinary(float, weight, totalElementsSampleInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing nclasses
    if(dataSet != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"n_classes.data");
        buffer = iftJoinPathnames(2, tempDirPath, "n_classes.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsGeneralInfo,2,fp);
        writeGenericVectorPointerAsBinary(int, (&(dataSet->nclasses)), totalElementsGeneralInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing functionNumber
    if(dataSet != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"function_number.data");
        buffer = iftJoinPathnames(2, tempDirPath, "function_number.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsGeneralInfo,2,fp);
        writeGenericVectorPointerAsBinary(int, (&(dataSet->function_number)), totalElementsGeneralInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing ngroups
    if(dataSet != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"n_labels.data");
        buffer = iftJoinPathnames(2, tempDirPath, "n_labels.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsGeneralInfo,2,fp);
        writeGenericVectorPointerAsBinary(int, (&(dataSet->ngroups)), totalElementsGeneralInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing ntrainSamples
    if(dataSet != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"n_train_samples.data");
        buffer = iftJoinPathnames(2, tempDirPath, "n_train_samples.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsGeneralInfo,2,fp);
        writeGenericVectorPointerAsBinary(int, (&(dataSet->ntrainsamples)), totalElementsGeneralInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }
    ////
    //writing iftRefDataType
    if(dataSet != NULL){
        // memset(buffer,0,sizeof(buffer));
        // strcat(buffer,tempDirPath);
        // strcat(buffer,"ref_data_type.data");
        buffer = iftJoinPathnames(2, tempDirPath, "ref_data_type.data");

        fp = fopen(buffer, "wb");
        if (fp == NULL){
            iftError(MSG_FILE_OPEN_ERROR, "iftWriteDataSet", buffer);
        }
        writeGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensionsGeneralInfo,2,fp);
        writeGenericVectorPointerAsBinary(iftRefDataType, (&(dataSet->ref_data_type)), totalElementsGeneralInfo, fp);
        fclose(fp);
        iftFree(buffer);
    }

    ////
    //////////////// WRITING REF. DATA ///////////////////
    if(dataSet != NULL) {
        char *ref_data_path = NULL;
        if (dataSet->ref_data_type == IFT_REF_DATA_IMAGE) {
            iftImage *img = (iftImage *) dataSet->ref_data;
            if (iftIs3DImage(img)) {
                ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.scn");
            } else { // 2D
                ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.png"); // saves gray and color images
            }
            iftWriteImageByExt(img, ref_data_path);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_FIMAGE) {
            iftFImage *fimg = (iftFImage *) dataSet->ref_data;
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.npy");
            iftWriteFImage(fimg, ref_data_path);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_MIMAGE) {
            iftMImage *mimg = (iftMImage *) dataSet->ref_data;
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.mimg");
            iftWriteMImage(mimg, ref_data_path);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_FILESET) {
            iftFileSet *fset = (iftFileSet *) dataSet->ref_data;
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.csv");
            iftCSV *csv = iftFileSetToCSV(fset);
            iftWriteCSV(csv, ref_data_path, ',');
            iftDestroyCSV(&csv);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_MMKERNEL) {
            iftMMKernel *K = (iftMMKernel *) dataSet->ref_data;
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.mmk");
            iftWriteMMKernel(K, ref_data_path);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_CSV) {
            iftCSV *csv = (iftCSV *) dataSet->ref_data;
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.csv");
            iftWriteCSV(csv, ref_data_path, ',');
        }
        iftFree(ref_data_path);
    }
    /////////////////////////////////////////////////////////////
    iftFree(truelabels);
    iftFree(labels);
    iftFree(id);
    iftFree(cluster_id);
    iftFree(status);
    iftFree(isSupervised);
    iftFree(isLabelPropagated);
    iftFree(weight);
    iftFree(numberTimesChecked);
    iftFree(dimensionsAlphas);
    iftFree(dimensionsFeat);
    iftFree(dimensionsGeneralInfo);
    iftFree(dimensionsProjection);
    iftFree(dimensionsSampleInfo);
    // ZIP ALL FILES
    iftZipDirContent(tempDirPath, pathname);
    // DESTROYERS
    iftRemoveDir(tempDirPath);
    iftFree(tempDirPath);
}



iftDataSet *iftReadDataSet(const char *format, ...) {
    if (format == NULL)
        iftError("Pathname is NULL", "iftReadDataSet");

    va_list args;
    char pathname[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(pathname, format, args);
    va_end(args);

    int randomNumber = iftRandomInteger(0,1073741824);
    char Numberstr[20];
    sprintf(Numberstr, "%d", randomNumber);
    char *tempDirPath = iftMakeTempDir("tmpdir_", NULL, NULL);
    char *buffer = NULL;
    iftUnzipFile(pathname, tempDirPath);
    FILE *fp;
    size_t totalSize;
    float* feats = NULL;
    float* alphas = NULL;
    double* projection = NULL;
    int* truelabels = NULL;
    int* labels = NULL;
    int *id = NULL;
    int *nChecked = NULL;
    int *cluster_id = NULL;
    uchar* status = NULL;
    uchar* isSupervised = NULL;
    uchar* isLabelPropagated = NULL;
    float* weight = NULL;
    iftDataSet* dataSet = NULL;
    // struct stat st = {0};
    // if (stat(tempDirPath, &st) == -1) {
    //     mkdir(tempDirPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // }

    ////
    buffer = iftJoinPathnames(2, tempDirPath, "version.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"version.data");
    fp = fopen(buffer, "rb");
    //old format
    if (fp == NULL){
        iftRemoveDir(tempDirPath);
        iftFree(tempDirPath);
        return iftReadOPFDataSet(pathname);
    }

    size_t * dimensions = iftAlloc(2,sizeof(size_t));
    int size;
    int nSamples = -1;
    int nFeats = -1;
    int projectionDimension = -1;
    int nclasses = 0;
    int functionNumber = 0;
    int ngroups = 0;
    int ntrainsamples = 0;
    iftRefDataType refDataType = IFT_REF_DATA_NULL;

    //reading feats
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "feat.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"feat.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(size_t, "%lu", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        feats = iftAlloc(totalSize,sizeof(float));
        readGenericVectorPointerAsBinary(float, feats, totalSize, fp);
        nSamples = dimensions[0];
        nFeats = dimensions[1];
        fclose(fp);
    }
    ////
    //reading Projection
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "projection.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"projection.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        projection = iftAlloc(totalSize,sizeof(double));
        readGenericVectorPointerAsBinary(double, projection, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        projectionDimension = dimensions[1];
        fclose(fp);
    }
    ////
    //reading alpha
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "alpha.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"alpha.data");
    fp = fopen(buffer, "rb");
    if(fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = dimensions[0];
        alphas = iftAlloc(totalSize,sizeof(float));
        readGenericVectorPointerAsBinary(float, alphas, totalSize, fp);
        fclose(fp);
    } 
    ////
    //reading truelabel
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "true_label.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"true_label.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        truelabels = iftAlloc(totalSize,sizeof(int));
        readGenericVectorPointerAsBinary(int, truelabels, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        fclose(fp);
    }
    ////
    //reading label
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "label.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"label.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        labels = iftAlloc(totalSize,sizeof(int));
        readGenericVectorPointerAsBinary(int, labels, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        fclose(fp);
    }
    ////
    //reading id
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "id.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"id.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        id = iftAlloc(totalSize,sizeof(int));
        readGenericVectorPointerAsBinary(int, id, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        fclose(fp);
    }
    ////
    ////
    //n checked
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "n_checked.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"n_checked.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        nChecked = iftAlloc(totalSize,sizeof(int));
        readGenericVectorPointerAsBinary(int, nChecked, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        fclose(fp);
    }
    ////
    //reading group
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "group.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"group.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        cluster_id = iftAlloc(totalSize,sizeof(int));
        readGenericVectorPointerAsBinary(int, cluster_id, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        fclose(fp);
    }
    ////
    //reading status
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "status.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"status.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        status = iftAlloc(totalSize,sizeof(uchar));
        readGenericVectorPointerAsBinary(uchar, status, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        fclose(fp);
    }
    ////
    //reading isSupervised
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "is_supervised.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"is_supervised.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        isSupervised = iftAlloc(totalSize,sizeof(uchar));
        readGenericVectorPointerAsBinary(uchar, isSupervised, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        fclose(fp);
    }

    ////
    //reading isLabelPropagated
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "is_label_propagated.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"is_label_propagated.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        isLabelPropagated = iftAlloc(totalSize,sizeof(uchar));
        readGenericVectorPointerAsBinary(uchar, isLabelPropagated, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        fclose(fp);
    }

    ////
    //reading weight
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "weight.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"weight.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' , '\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        weight = iftAlloc(totalSize,sizeof(float));
        readGenericVectorPointerAsBinary(float, weight, totalSize, fp);
        if(nSamples < 0){
            nSamples = dimensions[0];
        }
        fclose(fp);
    }
    ////
    //reading nclasses
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "n_classes.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"n_classes.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' ,'\n',dimensions,size,fp);
        totalSize = dimensions[0];
        readGenericVectorPointerAsBinary(int, (&nclasses), totalSize, fp);
        fclose(fp);
    }
    ////
    //reading functionNumber
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "function_number.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"function_number.data");
    fp = fopen(buffer, "rb");
    if(fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' ,'\n',dimensions,size,fp);
        totalSize = dimensions[0];
        readGenericVectorPointerAsBinary(int, (&functionNumber), totalSize, fp);
        fclose(fp);
    }
    ////
    //reading ngroups
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "n_labels.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"n_labels.data");
    fp = fopen(buffer, "rb");
    if(fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' ,'\n',dimensions,size,fp);
        totalSize = dimensions[0];
        readGenericVectorPointerAsBinary(int, (&ngroups), totalSize, fp);
        fclose(fp);
    }
    ////
    //reading ntrainSamples
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "n_train_samples.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"n_train_samples.data");
    fp = fopen(buffer, "rb");
    if(fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' ,'\n',dimensions,size,fp);
        totalSize = dimensions[0];
        readGenericVectorPointerAsBinary(int, (&ntrainsamples), totalSize, fp);
        fclose(fp);
    }
    ////
    //reading RefdataType
    iftFree(buffer);
    buffer = iftJoinPathnames(2, tempDirPath, "ref_data_type.data");
    // memset(buffer,0,sizeof(buffer));
    // strcat(buffer,tempDirPath);
    // strcat(buffer,"ref_data_type.data");
    fp = fopen(buffer, "rb");
    if (fp != NULL){
        size = 2;
        readGenericVectorPointerAsText(int, "%d", ',' ,'\n',dimensions,size,fp);
        totalSize = (size_t)dimensions[0]*(size_t)dimensions[1];
        readGenericVectorPointerAsBinary(iftRefDataType, (&refDataType), totalSize, fp);
        fclose(fp);
    }
    ////
    if(nSamples > 0 && nFeats > 0){
        dataSet = iftCreateDataSetAux(nSamples,nFeats,false);
        if(nFeats > 0 && feats != NULL){
            dataSet->nfeats = nFeats;
            size_t  shift = 0;
            for (int i = 0; i < nSamples; ++i) {
                dataSet->sample[i].feat = &(feats[shift]);
                shift += nFeats;
            }
            dataSet->data = (iftMatrix*)iftAlloc(1,sizeof(iftMatrix));
            dataSet->data->val = feats;
            dataSet->data->nrows = nSamples;
            dataSet->data->ncols = nFeats;
            dataSet->data-> n = nSamples*nFeats;
            dataSet->data->allocated = true;
        }
        if(projectionDimension > 0 && projection != NULL){
            dataSet->projection = iftCreateDoubleMatrix(projectionDimension,nSamples);
            dataSet->projection->val = projection;

            size_t  shift = 0;
            for (int i = 0; i < nSamples; ++i) {
                dataSet->sample[i].projection = &(projection[shift]);
                shift += projectionDimension;
            }
        }
        if(alphas != NULL){
            iftFree(dataSet->alpha);
            dataSet->alpha = alphas;
        } 
        for (int i = 0; i < nSamples; ++i) {
            if(truelabels != NULL){
                dataSet->sample[i].truelabel = truelabels[i];
            }
            if(labels != NULL){
                dataSet->sample[i].label = labels[i];
            }
            if(id != NULL){
                dataSet->sample[i].id = id[i];
            }
            if(cluster_id != NULL){
                dataSet->sample[i].group = cluster_id[i];
            }
            if(status != NULL){
                iftSetSampleStatus(&dataSet->sample[i], status[i]);
            }
            if(isSupervised != NULL){
                dataSet->sample[i].isSupervised = isSupervised[i];
            }
            if(isLabelPropagated != NULL){
                dataSet->sample[i].isLabelPropagated = isLabelPropagated[i];
            }
            if(weight  != NULL){
                dataSet->sample[i].weight = weight[i];
            }
            if(nChecked != NULL){
                dataSet->sample[i].numberTimesChecked = nChecked[i];
            }

        }
        dataSet->nclasses = nclasses;
        dataSet->ngroups = ngroups;
        dataSet->function_number = functionNumber;
        iftSetDistanceFunction(dataSet,functionNumber);
        dataSet->ntrainsamples = ntrainsamples;
        dataSet->ref_data_type = refDataType;
    }

    //////////////// Readinng REF. DATA ///////////////////
    if(dataSet != NULL) {
        char *ref_data_path = NULL;
        if (dataSet->ref_data_type == IFT_REF_DATA_IMAGE) {
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.scn");
            if (!iftFileExists(ref_data_path)) { // tests if a 3D image exists
                iftFree(ref_data_path);
                ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.png");

                if (!iftFileExists(ref_data_path)) { // tests if a 2D image exists
                    iftError("There is no Reference Data Image with extension *.[scn,png] in Dataset",
                             "iftWriteOPFDataset");
                }
            }

            dataSet->ref_data = iftReadImageByExt(ref_data_path);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_FIMAGE) {
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.npy");
            if (!iftFileExists(ref_data_path))
                iftError("There is no Reference Data FImage with extension *.fscn in Dataset",
                         "iftWriteOPFDataset");

            dataSet->ref_data = iftReadFImage(ref_data_path);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_MIMAGE) {
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.mimg");
            if (!iftFileExists(ref_data_path))
                iftError("There is no Reference Data MImage with extension *.mimg in Dataset",
                         "iftWriteOPFDataset");

            dataSet->ref_data = iftReadMImage(ref_data_path);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_FILESET) {
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.csv");
            if (!iftFileExists(ref_data_path))
                iftError("There is no Reference Data FILE \"%s\" in Dataset",
                         "iftWriteOPFDataset", ref_data_path);

            dataSet->ref_data = iftLoadFileSetFromCSV(ref_data_path, false);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_MMKERNEL) {
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.mmk");
            if (!iftFileExists(ref_data_path))
                iftError("There is no Reference Data MMKernel with extension *.mmk in Dataset",
                         "iftWriteOPFDataset");

            dataSet->ref_data = iftReadMMKernel(ref_data_path);
        } else if (dataSet->ref_data_type == IFT_REF_DATA_CSV) {
            ref_data_path = iftJoinPathnames(2, tempDirPath, "ref_data.csv");
            if (!iftFileExists(ref_data_path))
                iftError("There is no Reference Data FILE \"%s\" in Dataset",
                         "iftWriteOPFDataset", ref_data_path);

            dataSet->ref_data = iftReadCSV(ref_data_path, ',');
        }
        iftFree(ref_data_path);
    }

    iftFree(buffer);
    iftRemoveDir(tempDirPath);
    iftFree(tempDirPath);
    iftFree(truelabels);
    iftFree(labels);
    iftFree(id);
    iftFree(cluster_id);
    iftFree(status);
    iftFree(isSupervised);
    iftFree(weight);
    iftFree(dimensions);
    return dataSet;
}

iftDataSet* iftUpdateFileSetBaseNameFromDataset(iftDataSet *dataset, const char *newBaseName, const char *fileExtension){
    if(dataset == NULL){
        iftError("Dataset is empty","iftUpdateFileSetBaseNameFromDataset");
        return NULL;
    }
    if(dataset->ref_data_type != IFT_REF_DATA_FILESET || dataset->ref_data == NULL){
        iftError("Invalid reference data","iftUpdateFileSetBaseNameFromDataset");
        return NULL;
    }
    iftDataSet* outputDataset = iftCopyDataSet(dataset,true);
    iftFileSet* fileSetOriginal = (iftFileSet*)dataset->ref_data;
    iftFileSet* fileSetCopy = iftCopyFileSet(fileSetOriginal);
    iftCopyRefData(outputDataset, fileSetCopy, dataset->ref_data_type);
    // outputDataset->ref_data_type   = dataset->ref_data_type;
    // outputDataset->ref_data = (void*)fileSetCopy;

    for (int i = 0; i < fileSetOriginal->n; ++i) {
        char * SampleDirectory = iftParentDir(fileSetCopy->files[i]->path);
        char * SampleFileName = iftFilename(fileSetCopy->files[i]->path,fileExtension);
        char* fullPath = iftConcatStrings(3,newBaseName,SampleFileName,fileExtension);
        iftFree(fileSetCopy->files[i]->path);
        fileSetCopy->files[i]->path = fullPath;
        iftFree(SampleDirectory);
        iftFree(SampleFileName);
    }

//    for (int j = 0; j < outputDataset->nsamples; ++j) {
//        int sampleId = outputDataset->sample[j].id;
//        printf("%s %d\n",fileSetCopy->files[sampleId]->path,outputDataset->sample[j].truelabel);
//    }

    return outputDataset;
}

iftFeatSpaceParam iftReadFeatSpaceParam(const char *path) {
    if (path == NULL)
        iftError("Pathname is NULL", "iftReadFeatSpaceParam");

    iftFeatSpaceParam fsp;
    fsp.nfeats = 0;
    fsp.ncomps = 0;
    fsp.mean   = NULL;
    fsp.stdev  = NULL;
    fsp.w      = NULL;
    fsp.R      = NULL;
    fsp.W      = NULL;

    FILE *fp = fopen(path, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadFeatSpaceParam", path);

    if (fscanf(fp, "nfeats: %d\n", &fsp.nfeats) != 1)
        iftError("Error when reading Number of Feats", "iftReadFeatSpaceParam");

    int has_mean = 0;
    if (fscanf(fp, "has_mean: %d\n", &has_mean) != 1)
        iftError("Error when reading the Flag that indicates if the Mean Array is stored",
                 "iftReadFeatSpaceParam");
    if (has_mean) {
        fsp.mean = iftAllocFloatArray(fsp.nfeats);
        if (fread(fsp.mean, sizeof(float), fsp.nfeats, fp) != fsp.nfeats)
            iftError("Error when reading Mean Array of Feat. Space Param", "iftReadFeatSpaceParam");
    }

    int has_stdev = 0;
    if (fscanf(fp, "\nhas_stdev: %d\n", &has_stdev) != 1)
        iftError("Error when reading the Flag that indicates if the Stdev Array is stored",
                 "iftReadFeatSpaceParam");
    if (has_stdev) {
        fsp.stdev = iftAllocFloatArray(fsp.nfeats);
        if (fread(fsp.stdev, sizeof(float), fsp.nfeats, fp) != fsp.nfeats)
            iftError("Error when reading Stdev Array of Feat. Space Param", "iftReadFeatSpaceParam");
    }


    if (fscanf(fp, "\nncomps: %d\n", &fsp.ncomps) != 1)
        iftError("Error when reading the Number of Components", "iftReadFeatSpaceParam");

    int has_w = 0;
    if (fscanf(fp, "has_w: %d\n", &has_w) != 1)
        iftError("Error when reading the Flag that indicates if the w Array is stored",
                 "iftReadFeatSpaceParam");
    if (has_w) {
        fsp.w = iftAllocCharArray(fsp.ncomps);
        if (fread(fsp.w, sizeof(char), fsp.ncomps, fp) != fsp.ncomps)
            iftError("Error when reading W Array of Feat. Space Param", "iftReadFeatSpaceParam");
    }

    int has_rotation_matrix = 0;
    if (fscanf(fp, "\nhas_rotation_matrix: %d\n", &has_rotation_matrix) != 1)
        iftError("Error when reading the Flag that indicates if the Rotation Matrix is stored",
                 "iftReadFeatSpaceParam");
    if (has_rotation_matrix) {
        fsp.R = iftReadRawMatrixFILE(fp);
    }

    int has_whitening_matrix = 0;
    if (fscanf(fp, "\nhas_whitening_matrix: %d\n", &has_whitening_matrix) != 1)
        iftError("Error when reading the Flag that indicates if the Whitening Matrix is stored",
                 "iftReadFeatSpaceParam");
    if (has_whitening_matrix) {
        fsp.W = iftReadRawMatrixFILE(fp);
    }

    fclose(fp);

    return fsp;
}


void iftWriteFeatSpaceParam(const iftFeatSpaceParam *fsp, const char *path) {
    if (fsp == NULL)
        iftError("Feat. Space is NULL", "iftWriteFeatSpaceParam");
    if (path == NULL)
        iftError("Pathname is NULL", "iftWriteFeatSpaceParam");

    FILE *fp = fopen(path, "wb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteFeatSpaceParam", path);

    fprintf(fp, "nfeats: %d\n", (*fsp).nfeats);
    fprintf(fp, "has_mean: %d\n", (*fsp).mean != NULL);
    if ((*fsp).mean != NULL) {
        if (fwrite((*fsp).mean, sizeof(float), (*fsp).nfeats, fp) != (*fsp).nfeats)
            iftError("Error when writing Mean Array of Feat. Space Param", "iftWriteFeatSpaceParam");
    }

    fprintf(fp, "\nhas_stdev: %d\n", (*fsp).stdev != NULL);
    if ((*fsp).stdev != NULL) {
        if (fwrite((*fsp).stdev, sizeof(float), (*fsp).nfeats, fp) != (*fsp).nfeats)
            iftError("Error when writing Stdev Array of Feat. Space Param", "iftWriteFeatSpaceParam");
    }

    fprintf(fp, "\nncomps: %d\n", (*fsp).ncomps);
    fprintf(fp, "has_w: %d\n", (*fsp).w != NULL);
    if ((*fsp).w != NULL) {
        if (fwrite((*fsp).w, sizeof(char), (*fsp).ncomps, fp) != (*fsp).ncomps)
            iftError("Error when writing W Array of Feat. Space Param", "iftWriteFeatSpaceParam");
    }

    fprintf(fp, "\nhas_rotation_matrix: %d\n", (*fsp).R != NULL);
    if ((*fsp).R != NULL) {
        iftWriteRawMatrixFILE((*fsp).R, fp);
    }
    fprintf(fp, "\nhas_whitening_matrix: %d\n", (*fsp).W != NULL);
    if ((*fsp).W != NULL) {
        iftWriteRawMatrixFILE((*fsp).W, fp);
    }

    fclose(fp);
}


iftDataSet *iftExtractSamples(const iftDataSet *Z, iftSampleStatus status)
{
    int         i, s, t, nsamples;
    iftDataSet *Z1;

    nsamples=0;
    for (s=0; s < Z->nsamples; s++)
        if (iftHasSampleStatus(Z->sample[s], status)){
            nsamples++;
        }

    

    if (nsamples == 0)
        iftError("There are no samples from the desired status", "iftExtractSamples");

    Z1                  = iftCreateDataSet(nsamples,Z->nfeats);

    Z1->nclasses        = Z->nclasses;
    Z1->ngroups         = Z->ngroups;
    Z1->iftArcWeight    = Z->iftArcWeight;
    Z1->function_number = Z->function_number;
    iftCopyRefData(Z1, Z->ref_data, Z->ref_data_type);

    t=0;
    for (s=0; s < Z->nsamples; s++){
      if (iftHasSampleStatus(Z->sample[s], status)){
	if (Z->ref_data != NULL)
	  Z1->sample[t].id     = Z->sample[s].id;
	else
	  Z1->sample[t].id     = t;

	Z1->sample[t].isSupervised  = Z->sample[s].isSupervised;
	Z1->sample[t].isLabelPropagated  = Z->sample[s].isLabelPropagated;
	Z1->sample[t].numberTimesChecked  = Z->sample[s].numberTimesChecked;
	Z1->sample[t].truelabel  = Z->sample[s].truelabel;
	Z1->sample[t].label   = Z->sample[s].label;
	Z1->sample[t].group   = Z->sample[s].group;
	Z1->sample[t].weight  = Z->sample[s].weight;
	iftSetSampleStatus(&Z1->sample[t],Z->sample[s].status);
	
	for (i=0; i < Z->nfeats; i++)
	  Z1->sample[t].feat[i] = Z->sample[s].feat[i];
	t++;
      }
    }
    
    if (status & IFT_TRAIN){
        Z1->ntrainsamples = Z1->nsamples;
    }

    
    Z1->fsp = iftCopyFeatSpaceParam(Z);

    for (i=0; i < Z->nfeats; i++)
        Z1->alpha[i] = Z->alpha[i];

    return(Z1);
}

iftDataSet *iftExtractClass(iftDataSet *Z, int truelabel)
{
    int i, s, t, class_size=0;
    iftDataSet *Z1;

    if (!iftIsSupervisedDataSet(Z))
        iftError("It requires class information", "iftExtractClass");

    for (s=0; s < Z->nsamples; s++)
        if (Z->sample[s].truelabel == truelabel)
            class_size++;

    if (class_size==0){
        char msg[100];
        sprintf(msg,"There are no samples from the desired class %d",truelabel);
        iftWarning(msg,"iftExtractClass");
        return(NULL);
    }

    Z1 = iftCreateDataSet(class_size,Z->nfeats);
    Z1->nclasses = truelabel;
    Z1->ngroups = Z->ngroups;
    Z1->iftArcWeight = Z->iftArcWeight;
    Z1->function_number = Z->function_number;
    iftCopyRefData(Z1, Z->ref_data, Z->ref_data_type);

    t=0;
    for (s=0; s < Z->nsamples; s++)
        if (Z->sample[s].truelabel == truelabel){
	  if (Z->ref_data != NULL)
	    Z1->sample[t].id     = Z->sample[s].id;
	  else
	    Z1->sample[t].id     = t;
	  Z1->sample[t].weight = Z->sample[s].weight;
	  Z1->sample[t].truelabel = Z->sample[s].truelabel;
	  Z1->sample[t].group = Z->sample[s].group;
	  Z1->sample[t].label = Z->sample[s].label;
	  iftSetSampleStatus(&Z1->sample[t], IFT_TRAIN);
	  iftAddSampleStatus(&Z1->sample[t], IFT_SUPERVISED);
	  
	  for (i=0; i < Z->nfeats; i++)
	    Z1->sample[t].feat[i] = Z->sample[s].feat[i];
	  t++;
        }

    Z1->fsp = iftCopyFeatSpaceParam(Z);

    for (i=0; i < Z->nfeats; i++)
        Z1->alpha[i] = Z->alpha[i];

    return(Z1);
}

iftDataSet *iftExtractGroup(iftDataSet *Z, int group)
{
    int i, s, t, group_size=0;
    iftDataSet *Z1;

    if (Z->ngroups <= 0)
        iftError("It requires group information", "iftExtractGroup");

    for (s=0; s < Z->nsamples; s++)
        if (Z->sample[s].group == group)
            group_size++;

    if (group_size==0){
        char msg[100];
        sprintf(msg,"There are no samples from the desired group %d", group);
        iftWarning(msg,"iftExtractGroup");
        return(NULL);
    }

    Z1 = iftCreateDataSet(group_size,Z->nfeats);
    Z1->ngroups = group;
    Z1->nclasses = Z->nclasses;
    Z1->iftArcWeight = Z->iftArcWeight;
    Z1->function_number = Z->function_number;
    iftCopyRefData(Z1, Z->ref_data, Z->ref_data_type);

    t=0;
    for (s=0; s < Z->nsamples; s++)
        if (Z->sample[s].group == group){
	  if (Z->ref_data != NULL)
	    Z1->sample[t].id     = Z->sample[s].id;
	  else
	    Z1->sample[t].id     = t;
	  Z1->sample[t].weight = Z->sample[s].weight;
	  Z1->sample[t].truelabel = Z->sample[s].truelabel;
	  Z1->sample[t].group = Z->sample[s].group;
	  Z1->sample[t].label = Z->sample[s].label;
	  iftSetSampleStatus(&Z1->sample[t], IFT_TRAIN);
	  
            for (i=0; i < Z->nfeats; i++)
                Z1->sample[t].feat[i] = Z->sample[s].feat[i];
            t++;
        }

    Z1->fsp = iftCopyFeatSpaceParam(Z);

    for (i=0; i < Z->nfeats; i++)
        Z1->alpha[i] = Z->alpha[i];

    return(Z1);
}

iftDataSet *iftExtractObjectClasses(iftDataSet *Z)
{
    int i, s, t, class_size=0;
    iftDataSet *Z1;

    if (!iftIsSupervisedDataSet(Z))
        iftError("It requires class information", "iftExtractObjectClasses");

    for (s=0; s < Z->nsamples; s++)
        if (Z->sample[s].truelabel > 1)
            class_size++;

    if (class_size==0){
        char msg[100];
        sprintf(msg,"There are no samples from object classes");
        iftWarning(msg,"iftExtractObjectClasses");
        return(NULL);
    }

    Z1 = iftCreateDataSet(class_size,Z->nfeats);
    Z1->nclasses = 2;
    Z1->iftArcWeight = Z->iftArcWeight;
    Z1->function_number = Z->function_number;
    iftCopyRefData(Z1, Z->ref_data, Z->ref_data_type);

    t=0;
    for (s=0; s < Z->nsamples; s++)
        if (Z->sample[s].truelabel > 1){
	  if (Z->ref_data != NULL)
	    Z1->sample[t].id     = Z->sample[s].id;
	  else
	    Z1->sample[t].id     = t;
	  Z1->sample[t].label = Z->sample[s].label;
	  Z1->sample[t].weight = Z->sample[s].weight;
	  Z1->sample[t].truelabel = 2;
	  iftSetSampleStatus(&Z1->sample[t], IFT_TRAIN);
	  for (i=0; i < Z->nfeats; i++)
	    Z1->sample[t].feat[i] = Z->sample[s].feat[i];
	  t++;
        }

    Z1->fsp = iftCopyFeatSpaceParam(Z);

    for (i=0; i < Z->nfeats; i++)
        Z1->alpha[i] = Z->alpha[i];

    return(Z1);
}

iftImage *iftDataSetToLabelImage(const iftDataSet *Z, const iftImage *comp, bool decrement_labels, iftFeatureLabel label_type) {
    if (Z == NULL)
        iftError("Dataset is NULL", "iftDataSetToLabelImage");
    if (Z->ref_data == NULL)
        iftError("ReferData is NULL", "iftDataSetToLabelImage");
    if ((Z->ref_data_type != IFT_REF_DATA_IMAGE) &&
        (Z->ref_data_type != IFT_REF_DATA_FIMAGE) &&
        (Z->ref_data_type != IFT_REF_DATA_MIMAGE))
        iftError("Reference Data Type should be an Image", "iftDataSetToLabelImage");

    iftImage *label_img = NULL;
    if (Z->ref_data_type == IFT_REF_DATA_IMAGE) {
        iftImage *img = (iftImage*) Z->ref_data;
        label_img     = iftCreateImage(img->xsize, img->ysize, img->zsize);
        iftCopyVoxelSize(img, label_img);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE) {
        iftFImage *fimg = (iftFImage*) Z->ref_data;
        label_img       = iftCreateImage(fimg->xsize, fimg->ysize, fimg->zsize);
        iftCopyVoxelSize(fimg, label_img);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE) {
        iftMImage *mimg = (iftMImage*) Z->ref_data;
        label_img       = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
        iftCopyVoxelSize(mimg, label_img);
    }

    // decrement_labels is bool (0 for false, 1 for true)
    // thus, if true, the samples' labels will be decremented

    switch (label_type) {

    case IFT_LABEL:

      if (Z->nclasses == 0)
          iftError("DataSet has no labeled samples", "iftDataSetToLabelImage");

      for (int p = 0; p < Z->nsamples; p++) {
        label_img->val[Z->sample[p].id] = Z->sample[p].label - decrement_labels;
      }
      break;
    case IFT_CLASS:

      if (Z->nclasses == 0)
          iftError("DataSet has no labeled samples", "iftDataSetToLabelImage");

      for (int p = 0; p < Z->nsamples; p++) {
        label_img->val[Z->sample[p].id] = Z->sample[p].truelabel - decrement_labels;
      }
      break;
    case IFT_GROUP:

      if (Z->ngroups == 0)
          iftError("DataSet has no grouped samples", "iftDataSetToLabelImage");

      for (int p = 0; p < Z->nsamples; p++) {
        label_img->val[Z->sample[p].id] = Z->sample[p].group - decrement_labels;
      }
      break;
    default:
        iftError("Label type must be IFT_LABEL, IFT_CLASS or IFT_GROUP", "iftDataSetToLabelImage");
    }
    
    if (comp != NULL)
    {
        if (iftMaximumValue(comp) != Z->nsamples) {
            iftError("Number of compoments different from number of dataset samples", "iftDataSetToLabelImage");
        }
        for (int p = 0; p < label_img->n; p++)
        {
            int c = comp->val[p] - 1;
            label_img->val[p] = label_img->val[Z->sample[c].id];
        }
    }

    return label_img;
}

iftImage *iftDataSetClusterInformationToLabelImage(const iftDataSet *Z, const iftImage *comp, bool decrement_groups) {
    if (Z == NULL)
        iftError("Dataset is NULL", "iftDataSetClusterInformationToLabelImage");
    if (Z->ref_data == NULL)
        iftError("ReferData is NULL", "iftDataSetClusterInformationToLabelImage");
    if ((Z->ref_data_type != IFT_REF_DATA_IMAGE) &&
        (Z->ref_data_type != IFT_REF_DATA_FIMAGE) &&
        (Z->ref_data_type != IFT_REF_DATA_MIMAGE))
        iftError("Reference Data Type should be an Image", "iftDataSetClusterInformationToLabelImage");

    iftImage *label_img = NULL;
    if (Z->ref_data_type == IFT_REF_DATA_IMAGE) {
        iftImage *img = (iftImage*) Z->ref_data;
        label_img     = iftCreateImage(img->xsize, img->ysize, img->zsize);
        iftCopyVoxelSize(img, label_img);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE) {
        iftFImage *fimg = (iftFImage*) Z->ref_data;
        label_img       = iftCreateImage(fimg->xsize, fimg->ysize, fimg->zsize);
        iftCopyVoxelSize(fimg, label_img);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE) {
        iftMImage *mimg = (iftMImage*) Z->ref_data;
        label_img       = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
        iftCopyVoxelSize(mimg, label_img);
    }

    // decrement_labels is bool (0 for false, 1 for true)
    // thus, if true, the samples' groups will be decremented
    for (int p = 0; p < Z->nsamples; p++) {
        label_img->val[Z->sample[p].id] = Z->sample[p].group - decrement_groups;
    }

    if (comp != NULL)
    {
        if (iftMaximumValue(comp) != Z->nsamples) {
            iftError("Number of compoments different from number of dataset samples", "iftDataSetToLabelImage");
        }
        for (int p = 0; p < label_img->n; p++)
        {
            int c = comp->val[p] - 1;
            label_img->val[p] = label_img->val[Z->sample[c].id];
        }
    }

    return label_img;
}

iftImage *iftDataSetClustersToQuantizedImage(const iftDataSet *Z, bool decrement_groups) {
    if (Z == NULL)
        iftError("Dataset is NULL", "iftDataSetClustersToQuantizedImage");
    if (Z->ref_data == NULL)
        iftError("ReferData is NULL", "iftDataSetClustersToQuantizedImage");
    if ((Z->ref_data_type != IFT_REF_DATA_IMAGE) &&
        (Z->ref_data_type != IFT_REF_DATA_FIMAGE) &&
        (Z->ref_data_type != IFT_REF_DATA_MIMAGE))
        iftError("Reference Data Type should be an Image", "iftDataSetClustersToQuantizedImage");

    iftImage *quant_img = NULL;
    if (Z->ref_data_type == IFT_REF_DATA_IMAGE) {
        iftImage *img = (iftImage*) Z->ref_data;
        quant_img     = iftCreateColorImage(img->xsize, img->ysize, img->zsize, iftImageDepth(img));
        iftCopyVoxelSize(img, quant_img);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE) {
        iftFImage *fimg = (iftFImage*) Z->ref_data;
        quant_img       = iftCreateColorImage(fimg->xsize, fimg->ysize, fimg->zsize, iftFImageDepth(fimg));
        iftCopyVoxelSize(fimg, quant_img);
    }
    else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE) {
        iftMImage *mimg = (iftMImage*) Z->ref_data;
        quant_img       = iftCreateColorImage(mimg->xsize, mimg->ysize, mimg->zsize, iftMImageDepth(mimg));
        iftCopyVoxelSize(mimg, quant_img);
    }

    /* compute quantity of pixels for each class and sum of feature values */
    iftFloatArray **sum_val = iftAlloc(Z->nfeats, sizeof(iftFloatArray*));
    iftIntArray *count_val = iftCreateIntArray(Z->ngroups);

    for (int f = 0; f < Z->nfeats; f++) {
        sum_val[f] = iftCreateFloatArray(Z->ngroups);
    }
    #pragma omp parallel for
        for (int p = 0; p < Z->nsamples; p++) {
            count_val->val[Z->sample[p].group - decrement_groups] += 1;
            for (int f = 0; f < Z->nfeats; f++) {
                sum_val[f]->val[Z->sample[p].group - decrement_groups] += Z->sample[p].feat[f];
            }
        }

    /* set the mean color for each pixel */
    #pragma omp parallel for
        for (int p = 0; p < Z->nsamples; p++) {
            quant_img->val[Z->sample[p].id] = (int)(sum_val[0]->val[Z->sample[p].group - decrement_groups] / (float)(count_val->val[Z->sample[p].group - decrement_groups]));
            quant_img->Cb[Z->sample[p].id] = (ushort)(sum_val[1]->val[Z->sample[p].group - decrement_groups] / (float)(count_val->val[Z->sample[p].group - decrement_groups]));
            quant_img->Cr[Z->sample[p].id] = (ushort)(sum_val[2]->val[Z->sample[p].group - decrement_groups] / (float)(count_val->val[Z->sample[p].group - decrement_groups]));
        }

    return quant_img;
}

iftImage *iftDataSetWeight(iftDataSet *Z)
{
    iftImage *weight=NULL;
    int p;

    if (Z->ref_data_type != IFT_REF_DATA_IMAGE &&
        Z->ref_data_type != IFT_REF_DATA_MIMAGE &&
        Z->ref_data_type != IFT_REF_DATA_FIMAGE){
        iftError("Reference data must be an image", "iftDataSetWeight");
    }

    if (Z->ref_data_type == IFT_REF_DATA_IMAGE){
        iftImage *img = (iftImage *)Z->ref_data;
        weight = iftCreateImage(img->xsize,img->ysize,img->zsize);
    }else{
        if (Z->ref_data_type == IFT_REF_DATA_FIMAGE){
            iftFImage *img = (iftFImage *)Z->ref_data;
            weight = iftCreateImage(img->xsize,img->ysize,img->zsize);
        }else{
            if (Z->ref_data_type == IFT_REF_DATA_MIMAGE){
                iftMImage *img = (iftMImage *)Z->ref_data;
                weight = iftCreateImage(img->xsize,img->ysize,img->zsize);
            }
        }
    }

    if (weight != NULL){
        for (p=0; p < Z->nsamples; p++){
            weight->val[Z->sample[p].id] = iftRound(Z->sample[p].weight);
        }
    }else{
        iftError("It requires an associated image", "iftDataSetWeight");
    }

    return(weight);
}


iftMatrix* iftComputeConfusionMatrix(iftDataSet *dataset, bool normalized) {

    iftMatrix* confusionMatrix= iftCreateMatrix(dataset->nclasses,dataset->nclasses);
    int* countPerClass = iftAllocIntArray(dataset->nclasses);
    int col=0;
    int row=0;
    for (int sampleIndex = 0; sampleIndex < dataset->nsamples; ++sampleIndex) {

        if(iftHasSampleStatus(dataset->sample[sampleIndex], IFT_TEST)) {
            col = dataset->sample[sampleIndex].label - 1;
            row = dataset->sample[sampleIndex].truelabel - 1;
            iftMatrixElem(confusionMatrix, col, row) += 1;
            countPerClass[row]++;
        }
    }

    if(normalized){
        for (int i = 0; i < confusionMatrix->nrows; ++i) {
            if( countPerClass[i] == 0){
                continue;
            }
            for (int j = 0; j < confusionMatrix->ncols; ++j) {
                iftMatrixElem(confusionMatrix,j,i) /= (countPerClass[i]+0.0000001);
            }
        }
    }

    iftFree(countPerClass);
    return confusionMatrix;
}

iftMatrix *iftDataSetToFeatureMatrix(const iftDataSet *Z) {
    iftMatrix *X;

    /* Create the feature matrix */
    X = iftCreateMatrix(Z->nfeats, Z->nsamples);

    for (int i = 0; i < X->nrows; i++)
        for (int j = 0; j < X->ncols; j++)
            X->val[iftGetMatrixIndex(X, j, i)] = Z->sample[i].feat[j];

    return X;
}



iftDataSet *iftFeatureMatrixToDataSet(iftMatrix *X)
{
    iftDataSet *Z;
    int s,i,index,row,col;


    /* Create the DataSet */

    Z  = iftCreateDataSet(X->nrows, X->ncols);
    for (row=0, s=0, index=0; row < X->nrows; row++,s++) {
        for (col=0,i=0; col < X->ncols; col++,i++) {
            Z->sample[s].feat[i] = X->val[index];
            index++;
        }
    }

    return(Z);
}


void iftPrintDataSetInfo(iftDataSet *Z, bool printStatisticsInfo)
{
  int nclasses = iftCountNumberOfClassesDataSet(Z);
  int ngroups  = iftCountNumberOfGroupsDataSet(Z);
  int nsamples = 0;

    printf("Number of Samples: %d\n", Z->nsamples);

    nsamples = iftCountSamples(Z, IFT_TRAIN);
    if(nsamples>0)
        printf("\tTRAIN: %d\n", nsamples);
    
    nsamples = iftCountSamples(Z, IFT_SUPERVISED);
    if(nsamples>0)
        printf("\tSUPERVISED: %d\n", nsamples);

    nsamples = iftCountSamples(Z, IFT_LABELPROPAGATED);
    if(nsamples>0)
        printf("\tLABELPROPAGATED: %d\n", nsamples);
 
    nsamples = iftCountSamples(Z, IFT_TEST);
    if(nsamples>0)
        printf("\tTEST: %d\n", nsamples);

    nsamples = iftCountSamples(Z, IFT_OUTLIER);
    if(nsamples>0)
        printf("\tOUTLIER: %d\n", nsamples);

    nsamples = iftCountSamples(Z, IFT_ERROR);
    if(nsamples>0)
        printf("\tERROR: %d\n", nsamples);

    nsamples = iftCountSamples(Z, IFT_UNKNOWN);
    if(nsamples>0)
        printf("\tUNKNOWN: %d\n", nsamples);

    printf("\nNumber of Classes: %d\n", nclasses);
    printf("\nNumber of Groups: %d\n", ngroups);

    int       *truelabel = iftAllocIntArray(Z->nsamples);
    int       *group     = iftAllocIntArray(Z->nsamples);

    for (int s=0; s < Z->nsamples; s++) { 
      truelabel[s]=Z->sample[s].truelabel;
      group[s]=Z->sample[s].group;
    }

    int *elem, *quant;
    int ntruelabels = iftCountElemsPerUniqueIntElem(truelabel, Z->nsamples, &elem, &quant);
    
    for (int l=0; l < ntruelabels; l++) {
      printf("truelabel %d: %d samples\n",elem[l],quant[l]);
    }
    iftFree(elem); iftFree(quant);

   ngroups = iftCountElemsPerUniqueIntElem(group, Z->nsamples, &elem, &quant);
    
    for (int l=0; l < ngroups; l++) {
      printf("group %d: %d samples\n",elem[l],quant[l]);
    }
    iftFree(elem); iftFree(quant);
    
    iftFree(truelabel); iftFree(group);
    
    printf("\nFeature Vector Size: %d\n", Z->nfeats);

    if(printStatisticsInfo == false){
        return;
    }

    iftFloatArray* mean = iftCreateFloatArray(Z->nfeats);
    iftFloatArray* std =  iftCreateFloatArray(Z->nfeats);
    iftFloatArray* max = iftCreateFloatArray(Z->nfeats);
    iftFloatArray* min = iftCreateFloatArray(Z->nfeats);

    iftSetFloatArray(min->val, min->n, IFT_INFINITY_FLT);
    iftSetFloatArray(max->val, max->n, IFT_INFINITY_FLT_NEG);

    for (int i = 0; i < Z->nsamples; ++i) {
        for (int j = 0; j < Z->nfeats; ++j) {
            min->val[j]   = iftMin(min->val[j], Z->sample[i].feat[j]);
            max->val[j]   = iftMax(max->val[j], Z->sample[i].feat[j]);
            mean->val[j] += Z->sample[i].feat[j]/Z->nsamples;
        }
    }

    for (int i = 0; i < Z->nsamples; ++i) {
        for (int j = 0; j < Z->nfeats; ++j) {
            std->val[j] += iftFastPow(Z->sample[i].feat[j] - mean->val[j], 2)/(Z->nsamples-1);
        }
    }

    for (int j = 0; j < Z->nfeats; ++j) {
        std->val[j] = sqrt(std->val[j]);
    }

    printf("\nFeatures Info\n");
    printf("\tMean: "); iftPrintFormattedFloatArray(mean->val, mean->n, "\t{", "}\n", ",\t");
    printf("\tStd: "); iftPrintFormattedFloatArray(std->val, mean->n, "\t{", "}\n", ",\t");
    printf("\tMin: "); iftPrintFormattedFloatArray(min->val, mean->n, "\t{", "}\n", ",\t");
    printf("\tMax: "); iftPrintFormattedFloatArray(max->val, mean->n, "\t{", "}\n", ",\t");

    for (int c = 1; c <= Z->nclasses; ++c) {

        iftSetFloatArray(mean->val, mean->n, 0.0);
        iftSetFloatArray(max->val, max->n, IFT_INFINITY_FLT_NEG);
        iftSetFloatArray(min->val, min->n, IFT_INFINITY_FLT);
        nsamples = iftCountSamplesTrueLabel(Z, c, IFT_ALL);

        for (int i = 0; i < Z->nsamples; ++i) {
            if(Z->sample[i].truelabel == c) {
                for (int j = 0; j < Z->nfeats; ++j) {
                    min->val[j] = iftMin(min->val[j], Z->sample[i].feat[j]);
                    max->val[j] = iftMax(max->val[j], Z->sample[i].feat[j]);
                    mean->val[j] += Z->sample[i].feat[j] / nsamples;
                }
            }
        }

        for (int i = 0; i < Z->nsamples; ++i) {
            if(Z->sample[i].truelabel == c) {
                for (int j = 0; j < Z->nfeats; ++j) {
                    std->val[j] += iftFastPow(Z->sample[i].feat[j] - mean->val[j], 2) / (nsamples - 1);
                }
            }
        }

        for (int j = 0; j < Z->nfeats; ++j) {
            std->val[j] = sqrt(std->val[j]);
        }

        printf("\t#%d\n", c);
        printf("\t\tMean: "); iftPrintFormattedFloatArray(mean->val, mean->n, "\t{", "}\n", ",\t");
        printf("\t\tStd: "); iftPrintFormattedFloatArray(std->val, mean->n, "\t{", "}\n", ",\t");
        printf("\t\tMin: "); iftPrintFormattedFloatArray(min->val, mean->n, "\t{", "}\n", ",\t");
        printf("\t\tMax: "); iftPrintFormattedFloatArray(max->val, mean->n, "\t{", "}\n", ",\t");

    }

    iftDestroyFloatArray(&mean);
    iftDestroyFloatArray(&std);
    iftDestroyFloatArray(&max);
    iftDestroyFloatArray(&min);
}

iftMatrix *iftDatasetCovarianceMatrix(iftDataSet *Z) {
    iftMatrix *cov, *X, *Xt;
    float K = 1.0 / ((float) Z->nsamples - 1);

    if (!iftIsCentralizedDataSet(Z))
        iftError("Dataset must be centralized", "iftDatasetCovarianceMatrix");

    X = iftDataSetToFeatureMatrix(Z);
    Xt = iftTransposeMatrix(X);
    iftMultMatrixByScalar(X, K);
    cov = iftMultMatrices(Xt, X);

    iftDestroyMatrix(&X);
    iftDestroyMatrix(&Xt);

    return (cov);
}


/*
* The datasets must have the same number classes, samples and features.
* Use iftValidateDataSets for that.
*/
iftFeatSpaceParam iftComputeOverallMeanAndStdevFromDatasets(iftDataSet **Zs, int num_Z) {
    iftFeatSpaceParam fsp;
    int nsamples = 0, nfeats = 0, ntrainsamples = 0;

    fsp.mean = fsp.stdev = NULL;
    fsp.R = NULL;
    fsp.W = NULL;
    fsp.w = NULL;
    fsp.nfeats = 0;
    fsp.ncomps = 0;

    if (num_Z >= 1) {
        nfeats   = Zs[0]->nfeats;
        nsamples = Zs[0]->nsamples;
    } else
        iftError("Number of Datasets < 1", "iftComputeOverallMeanAndStdevFromDatasets");

    for (int it = 0; it < num_Z; it++)
        if ((!iftIsTrainingDataSet(Zs[it])) || iftIsNormalizedDataSet(Zs[it]))
            iftError("It requires a non-normalized Training set", "iftComputeOverallMeanAndStdevFromDatasets");

    fsp.mean   = iftAllocFloatArray(nfeats);
    fsp.stdev  = iftAllocFloatArray(nfeats);
    fsp.nfeats = nfeats;

    // Compute Overall Average
    puts("# Computing Overall Feature Average #\n");
    for (int it = 0; it < num_Z; it++) {
        for (int s = 0; s < nsamples; s++) {
            for (int i = 0; i < nfeats; i++) {
                fsp.mean[i] += Zs[it]->sample[s].feat[i];
            }
        }
        ntrainsamples += Zs[it]->ntrainsamples;
    }

    for (int i = 0; i < nfeats; i++)
        fsp.mean[i] /= ntrainsamples;

    // Compute Standard Deviation
    puts("# Computing Overall Feature Standard Deviation #\n");
    for (int it = 0; it < num_Z; it++) {
        for (int s = 0; s < nsamples; s++) {
            for (int i = 0; i < nfeats; i++) {
                fsp.stdev[i] += (Zs[it]->sample[s].feat[i] - fsp.mean[i])
                                * (Zs[it]->sample[s].feat[i] - fsp.mean[i]);
            }
        }
    }

    for (int i = 0; i < nfeats; i++)
        fsp.stdev[i] = sqrtf(fsp.stdev[i] / (ntrainsamples - 1));

    return fsp;
}

void iftNormalizeDataSetByZScoreInPlace(iftDataSet *Z, iftFeatSpaceParam *fsp, float stdev_factor)
{
    int s,i;
    float *mean,*stdev;

    if (iftIsNormalizedDataSet(Z))
        iftError("It requires a non-normalized training set", "iftNormalizeDataSetByZScoreInPlace");

    if (fsp == NULL || fsp->mean == NULL){
        /* we need to compute the mean and the standard deviation here*/
        mean=iftAllocFloatArray(Z->nfeats);
        stdev=iftAllocFloatArray(Z->nfeats);

        for (s=0; s < Z->nsamples; s++) {
            if(iftHasSampleStatus(Z->sample[s], IFT_TRAIN)) {
                for (i = 0; i < Z->nfeats; i++) {
                    mean[i] += Z->sample[s].feat[i];
                }
            }
        }
        for (i=0; i < Z->nfeats; i++){
            mean[i] /= Z->ntrainsamples;
        }

        for (s=0; s < Z->nsamples; s++) {
            if(iftHasSampleStatus(Z->sample[s], IFT_TRAIN)) {
                for (i = 0; i < Z->nfeats; i++) {
                    stdev[i] += (Z->sample[s].feat[i] - mean[i]) * (Z->sample[s].feat[i] - mean[i]);
                }
            }
        }
        for (i=0; i < Z->nfeats; i++) {
            stdev[i] = sqrtf( stdev[i] / Z->ntrainsamples ) + stdev_factor;
        }
    }
    else{
        if (fsp->nfeats != Z->nfeats) {
            iftError("Number is different between Dataset and FeatSpaceParam - (%d - %d)", "iftNormalizeDataSetByZScoreInPlace",
                fsp->nfeats, Z->nfeats);
        }
        mean=fsp->mean;
        stdev=fsp->stdev;
    }

    for (s=0; s < Z->nsamples; s++) {
        for (i=0; i < Z->nfeats; i++) {
	  Z->sample[s].feat[i] = (Z->sample[s].feat[i]-mean[i])/(stdev[i]);
	}	
    }
    
    if (iftIsCentralizedDataSet(Z)){
        iftFree(Z->fsp.mean);
    }

    Z->fsp.mean   = mean;
    Z->fsp.stdev  = stdev;
    Z->fsp.nfeats = Z->nfeats;
}

iftDataSet *iftNormalizeDataSetByZScore(iftDataSet *Z, iftFeatSpaceParam *fsp, float stdev_factor)
{
    if (iftIsNormalizedDataSet(Z))
        iftError("It requires a non-normalized training set", "iftNormalizeDataSetByZScore");

    iftDataSet *Zn = iftCopyDataSet(Z, true);
    iftNormalizeDataSetByZScoreInPlace(Zn, fsp, stdev_factor);
    
    return Zn;
}

iftDataSet *iftNormOneDataSet(iftDataSet *Z)
{
    int s,i;
    iftDataSet *Zn = iftCopyDataSet(Z, true);
    double norm;

    for (s=0; s < Zn->nsamples; s++) {
        norm = 0.0;
        for (i=0; i < Zn->nfeats; i++) {
            norm += Zn->sample[s].feat[i]*Zn->sample[s].feat[i];
        }
        if (norm > IFT_EPSILON) {
            norm = sqrt(norm);
            for (i=0; i < Zn->nfeats; i++) {
                Zn->sample[s].feat[i] /= norm;
            }
        } else {
            iftError("DataSet contains 0-feature vectors", "iftNormOneDataSet");
        }
    }

    return(Zn);
}

void iftNormalizeSamples(iftDataSet *Z)
{
    int s,i;
    float norm;
    float *mean=iftAllocFloatArray(Z->nfeats);


    if (!iftIsCentralizedDataSet(Z)){
        for (s=0; s < Z->nsamples; s++) {
            for (i=0; i < Z->nfeats; i++) {
                mean[i] += Z->sample[s].feat[i];
            }
        }

        for (i=0; i < Z->nfeats; i++){
            mean[i] /= Z->nsamples;
        }
    }

    for (s=0; s < Z->nsamples; s++) {
        norm = 0;
        for (i=0; i < Z->nfeats; i++) {
            Z->sample[s].feat[i] = (Z->sample[s].feat[i]-mean[i]);
            norm += Z->sample[s].feat[i]*Z->sample[s].feat[i];
        }
        norm = sqrtf(norm);
        if (norm > IFT_EPSILON){
            for (i=0; i < Z->nfeats; i++) {
                Z->sample[s].feat[i] = Z->sample[s].feat[i]/norm;
            }
        }
    }

    iftFree(mean);
}

void iftPowerNormalizeDataSetInPlace(iftDataSet *Z, float alpha)
{
    for (int s=0; s < Z->nsamples; s++) {
        for (int i=0; i < Z->nfeats; i++) {
            float sign = Z->sample[s].feat[i] > 0 ? 1 : -1;
            Z->sample[s].feat[i] = powf(Z->sample[s].feat[i], alpha) * sign;
        }
    }
}

void iftNormalizeDataSetByL2NormInPlace(iftDataSet *Z)
{
    for (int s=0; s < Z->nsamples; s++) {
        float sum = 0;
        for (int i=0; i < Z->nfeats; i++) {
            sum += powf(Z->sample[s].feat[i], 2.0);
        }
        for (int i=0; i < Z->nfeats; i++) {
            Z->sample[s].feat[i] /= sqrtf(sum);
        }
    }
}

iftDataSet* iftReadCSVDataSet(const char* filepath, char separator, int NumberColumnsAfterFeatures, char *fileset) {

  iftCSV* csv = iftReadCSV(filepath, separator);


  iftDataSet* dataset = iftCreateDataSet(csv->nrows, csv->ncols - NumberColumnsAfterFeatures);

    for (int i=0; i < csv->nrows; ++i) {
      dataset->sample[i].id = IFT_NIL;
      for(int j=0; j<csv->ncols-NumberColumnsAfterFeatures; ++j) {
	
	float val = 0.0f;
	sscanf(csv->data[i][j], "%f", &val);
	dataset->sample[i].feat[j] = val;
      }
    }

    if (NumberColumnsAfterFeatures == 1) { 
      for (int i=0; i < csv->nrows; ++i) {
	int val = 0;
	sscanf(csv->data[i][csv->ncols-1], "%d", &val);
	dataset->sample[i].truelabel = val;
      }
    } else {
      if (NumberColumnsAfterFeatures == 2) { 
	for (int i=0; i < csv->nrows; ++i) {
	  int val = 0;
	  sscanf(csv->data[i][csv->ncols-2], "%d", &val);
	  dataset->sample[i].truelabel = val;	  
	  sscanf(csv->data[i][csv->ncols-1], "%d", &val);
	  dataset->sample[i].id = val;
	}
      }
    }

    if (fileset != NULL) {
      dataset->ref_data_type = IFT_REF_DATA_FILESET;
      if (!iftFileExists(fileset))
	iftError("There is no Reference Data FILE \"%s\" in Dataset",
		 "iftReadCSVDataSet", fileset);      
      dataset->ref_data = iftLoadFileSetFromCSV(fileset, false);
    }

    iftCountNumberOfClassesDataSet(dataset);
    iftDestroyCSV(&csv);

    return dataset;
}

iftDataSet* iftReadCSVImageDataSet(const char* filepath, char separator,
                                   bool isSupervised, int truelabelColumn,
                                   const char* imagesPath) {
    iftCSV* csv = iftReadCSV(filepath, separator);
    iftDataSet* dataset = NULL;
    FILE* fileImagesPath = fopen(imagesPath, "r"); /* should check the result */
    char **matrixPaths = (char **)calloc(csv->nrows,sizeof(char *));
    char line[256];

    if(isSupervised){
        if(truelabelColumn < 0) {
            truelabelColumn = csv->ncols + truelabelColumn;
        }
        dataset = iftCreateDataSet(csv->nrows, csv->ncols - isSupervised);

        for (int i=0; i < csv->nrows; ++i) {
            size_t featIndex = 0;
            if (fgets(line, sizeof(line), fileImagesPath) == NULL)
                iftError("Error when reading a line", "iftReadCSVImageDataSet");
            matrixPaths[i] = calloc(256,sizeof(char));
            for(int j=0; j<csv->ncols; ++j) {
                float val = 0.0f;
                sscanf(csv->data[i][j], "%f", &val);

                //dataset->sample[i].id = i;
                if(truelabelColumn == j) {
                    dataset->sample[i].truelabel = (int) val;
                }
                else{
                    dataset->sample[i].feat[featIndex] = val;
                    featIndex++;
                }

            }
        }
    }
    iftCountNumberOfClassesDataSet(dataset);
    fclose(fileImagesPath);
    iftDestroyCSV(&csv);
    return dataset;
}




iftDataSet* iftDataSetTransformBasis(const iftDataSet *Z, const iftMatrix *T) {

    if(Z->nfeats!=T->nrows) {
        iftError("Transformation matrix and dataset have incompatible dimensions. Z = (%d,%d), T = (%d, %d)",
                 "iftDataSetTransformBasis", Z->nsamples, Z->nfeats, T->nrows, T->ncols);
    }

    iftDataSet* ZT = iftCreateDataSet(Z->nsamples, T->ncols);

    iftMatrix* feats = iftCreateMatrix(Z->nfeats, 1);
    iftMatrix* featsT = iftCreateMatrix(ZT->nfeats, 1);

    for (int i = 0; i < Z->nsamples; ++i) {

        for (int j = 0; j < Z->nfeats; ++j)
            feats->val[j] = Z->sample[i].feat[j];

        iftMultMatricesInPlace(feats, T, false, false, &featsT);

        for (int j = 0; j < ZT->nfeats; ++j)
            ZT->sample[i].feat[j] = featsT->val[j];

        ZT->sample[i].id     = Z->sample[i].id;
        ZT->sample[i].truelabel  = Z->sample[i].truelabel;
        ZT->sample[i].label  = Z->sample[i].label;
        iftSetSampleStatus(&ZT->sample[i], Z->sample[i].status);
        ZT->sample[i].weight = Z->sample[i].weight;
    }

    ZT->nclasses        = Z->nclasses;
    ZT->ngroups         = Z->ngroups;
    ZT->iftArcWeight    = Z->iftArcWeight;
    ZT->function_number = Z->function_number;
    iftCopyRefData(ZT, Z->ref_data, Z->ref_data_type);
    // ZT->ref_data        = Z->ref_data;
    // ZT->ref_data_type   = Z->ref_data_type;
    ZT->ntrainsamples = Z->ntrainsamples;

    ZT->fsp = iftCopyFeatSpaceParam(Z);
    
    //maybe we should multiply to concatenate multiple rotation?
    ZT->fsp.R = iftCopyMatrix(T);
    ZT->fsp.ncomps = T->ncols;

    iftDestroyMatrix(&feats);
    iftDestroyMatrix(&featsT);

    return ZT;
}

iftDataSet *iftWhiteningTransform(iftDataSet *Z)
{
    int           s,i,row,col;
    iftDataSet   *Zw = NULL, *Zn=NULL,*Zc=NULL;
    iftMatrix    *Cov,*U,*Uinv,*S,*Vt,*X,*Sinv,*W,*SV,*Y;

    if ((!iftIsTrainingDataSet(Z))||(iftIsWhitenedDataSet(Z)))
        iftError("It requires a non-whitened training set", "iftWhiteningTransform");

    if (!iftIsCentralizedDataSet(Z)){ /* translate dataset */
        Zc = iftNormalizeContrastDataSet(Z);
        Zn = iftCentralizeDataSet(Zc);
        iftDestroyDataSet(&Zc);
    }else{
        Zn = Z;
    }

    /* Compute Principal Components */

    Cov  = iftDatasetCovarianceMatrix(Zn);
    iftSingleValueDecomp(Cov,&U,&S,&Vt);
    iftDestroyMatrix(&Cov);

    /* Compute S^(-1/2) for whitening */

    Sinv = iftCreateMatrix(Zn->nfeats,Zn->nfeats);
    float Seps = powf(10.,(log(Zn->nfeats)/log(2)-2)/2.); // 8 => 0.1, 16 => 0.01
    for (i=0; i < Zn->nfeats; i++)
        Sinv->val[iftGetMatrixIndex(Sinv,i,i)] = 1.0/(sqrt(S->val[i]+Seps));
//    Sinv->val[iftGetMatrixIndex(Sinv,i,i)] = 1.0/(sqrt(S->val[i]+IFT_EPSILON));
//    Sinv->val[iftGetMatrixIndex(Sinv,i,i)] = 1.0/(sqrt(S->val[i]+0.1)); // Menotti on 20130917 - good for 8x8 patches
//Sinv->val[iftGetMatrixIndex(Sinv,i,i)] = 1.0/(sqrt(S->val[i]+0.01)); // Menotti on 20130917 - good for 16x16 patches

    iftDestroyMatrix(&S);

    /* Compute Whitening Transform W = U Sinv Uinv */

    Uinv = iftTransposeMatrix(U);
    SV   = iftMultMatrices(Sinv,U);
    W    = iftMultMatrices(Uinv,SV);

    iftDestroyMatrix(&SV);
    iftDestroyMatrix(&Vt);
    iftDestroyMatrix(&U);
    iftDestroyMatrix(&Uinv);
    iftDestroyMatrix(&Sinv);

    /* Apply whitening transformation: Y = X W , where X is feature
     matrix of Zn */

    X = iftDataSetToFeatureMatrix(Zn);
    Y = iftMultMatrices(X,W);

    iftDestroyMatrix(&X);

    Zw = iftCreateDataSet(Zn->nsamples,Zn->nfeats);
    Zw->iftArcWeight = Zn->iftArcWeight;
    Zw->function_number = Zn->function_number;
    Zw->nclasses        = Zn->nclasses;
    iftCopyRefData(Zw, Zn->ref_data, Zn->ref_data_type);
    // Zw->ref_data        = Zn->ref_data;
    // Zw->ref_data_type   = Zn->ref_data_type;
    for (s=0, row=0; s < Zw->nsamples; s++, row++) {
        Zw->sample[s].truelabel  = Zn->sample[s].truelabel;
        iftSetSampleStatus(&Zw->sample[s], Zn->sample[s].status);
        Zw->sample[s].id     = Zn->sample[s].id;

        /* Feature Matrix to DataSet, correcting initial translation */

        for (i=0, col=0; col < Zw->nfeats; col++, i++)
            Zw->sample[s].feat[i] = Y->val[iftGetMatrixIndex(Y,col,row)] ;
    }

    iftDestroyMatrix(&Y);

    Zw->ntrainsamples  = Zn->ntrainsamples;
    for (i=0; i < Zw->nfeats; i++){
        Zw->alpha[i]=1.0;
    }

    /* Copy feature space parameters */

    Zw->fsp    = iftCopyFeatSpaceParam(Zn);
    Zw->fsp.W  = iftCopyMatrix(W);

    if (Zn != Z){
        iftDestroyDataSet(&Zn);
    }


    iftDestroyMatrix(&W);

    return(Zw);

}

iftDataSet *iftNormalizeContrastDataSet(iftDataSet *Z)
{
    int s,i;
    float mean, stdev;
    iftDataSet *Zn = iftCopyDataSet(Z, true);


    if (!iftIsTrainingDataSet(Z))
        iftError("It requires a training set", "iftNormalizeContrastDataSet");

    /* This is a feature normalization for kernel datasets in
     convolution neural networks. Since it is done independently
     sample by sample, it does not count as a normalized dataset for
     classification purpose */

    for (s=0; s < Zn->nsamples; s++) {
        mean = 0.0;
        for (i=0; i < Zn->nfeats; i++) {
            mean += Zn->sample[s].feat[i];
        }
        mean /= Zn->nfeats;

        stdev = 0.0;
        for (i=0; i < Zn->nfeats; i++) {
            Zn->sample[s].feat[i] = Zn->sample[s].feat[i] - mean;
            stdev += Zn->sample[s].feat[i]*Zn->sample[s].feat[i];
        }
        stdev = sqrtf( (stdev + 0.1) / (Zn->nfeats-1) );

//    if (stdev > IFT_EPSILON)
        for (i=0; i < Zn->nfeats; i++) {
            Zn->sample[s].feat[i] /= stdev;
        }
    }

    return(Zn);

}

void iftCentralizeDataSetInPlace(iftDataSet *Z)
{
    int s,i;
    float *mean=iftAllocFloatArray(Z->nfeats);

    if(Z->ntrainsamples<=0) {
        iftFree(mean);
        iftError("It requires at least one training sample.", "iftCentralizeDataSetInPlace");
    }

    if (iftIsCentralizedDataSet(Z)||iftIsNormalizedDataSet(Z)) {
        iftFree(mean);
        iftError("It requires a non-centralized data set", "iftCentralizeDataSetInPlace");
    }

    for (s=0; s < Z->nsamples; s++) {
        if(iftHasSampleStatus(Z->sample[s], IFT_TRAIN)) {
            for (i = 0; i < Z->nfeats; i++) {
                mean[i] += Z->sample[s].feat[i];
            }
        }
    }

    for (i=0; i < Z->nfeats; i++) {
        mean[i] /= Z->ntrainsamples;
    }

    for (s=0; s < Z->nsamples; s++) {
        for (i = 0; i < Z->nfeats; i++) {
            Z->sample[s].feat[i] = (Z->sample[s].feat[i] - mean[i]);
        }
    }

    Z->fsp.mean = mean;
    Z->fsp.nfeats = Z->nfeats;
}

iftDataSet *iftCentralizeDataSet(iftDataSet *Z)
{
    iftDataSet *Zc = iftCopyDataSet(Z, true);
    iftCentralizeDataSetInPlace(Zc);

    return(Zc);
}

iftDataSet* iftUnitNormDataSet(iftDataSet* Z) {
    int n = Z->nsamples;
    int dim = Z->nfeats;
    iftDataSet *Zun = iftCopyDataSet(Z, true);

    for (int s=0; s < n ; s++) {
        float norm = 0;
        for (int f = 0; f < dim; f++ )
            norm += Zun->sample[s].feat[f] * Zun->sample[s].feat[f];
        norm = sqrt(norm + IFT_EPSILON);
        for( int f = 0; f < dim; f++ )
            Zun->sample[s].feat[f] /= norm;
    }

    return (Zun);
}

void iftCopyClassifResult(iftDataSet *Z1, iftDataSet *Z2)
{
    int s;

    if ( (Z1->nsamples != Z2->nsamples) ||
         ((Z1->nclasses > 0)&&(Z1->nclasses != Z2->nclasses)) )
        iftError("Incompatible datasets", "iftCopyClassifResult");


    Z2->ngroups     = Z1->ngroups;
    Z2->ntrainsamples = Z1->ntrainsamples;
    for (s=0; s < Z1->nsamples; s++) {
        Z2->sample[s].weight = Z1->sample[s].weight;
        Z2->sample[s].label  = Z1->sample[s].label;
        iftSetSampleStatus(&Z2->sample[s], Z1->sample[s].status);
    }

}

iftDataSet *iftImageBorderDataSet(iftDataSet *Z1, int size, int nsamples)
{
    iftDataSet *Z;
    iftImage *img = (iftImage *)Z1->ref_data;
    iftVoxel u;
    int  p, i, j, s, t, dx, dy, dz, n;
    int *sample = iftAllocIntArray(Z1->nsamples);

    if (img == NULL) {
        iftError("This is not an image dataset", "iftBorderFromImageDataSet");
    }

    /* Get candidates on the image's border */

    i = 0;
    if (iftIs3DImage(img)){
        for (s=0; s < Z1->nsamples; s++) {
            p  = Z1->sample[s].id;
            u  = iftGetVoxelCoord(img,p);
            dx = iftMin(u.x, img->xsize - 1 - u.x);
            dy = iftMin(u.y, img->ysize - 1 - u.y);
            dz = iftMin(u.z, img->zsize - 1 - u.z);
            if ((dx < size)||(dy < size)||(dz < size)){ // p is on the image's border
                sample[i] = s; i++;
            }
        }
    }else{ /* 2D image */
        for (s=0; s < Z1->nsamples; s++) {
            p  = Z1->sample[s].id;
            u  = iftGetVoxelCoord(img,p);
            dx = iftMin(u.x, img->xsize - 1 - u.x);
            dy = iftMin(u.y, img->ysize - 1 - u.y);
            if ((dx < size)||(dy < size)){ // p is on the image's border
                sample[i] = s; i++;
            }
        }
    }
    n=i;

    if (n < nsamples){
        iftWarning("Adjusting number of samples","iftBorderFromImageDataSet");
        nsamples = n;
    }


    Z        = iftCreateDataSet(nsamples,Z1->nfeats);
    iftCopyRefData(Z, Z1->ref_data, Z1->ref_data_type);
    // Z->ref_data        = Z1->ref_data;
    // Z->ref_data_type   = Z1->ref_data_type;

    /* Select random samples from the image's border */
    t = 0;
    while (nsamples > 0) {
        i = iftRandomInteger(0,n-1);
        s = sample[i];
        p = Z1->sample[s].id;
        for (j=0; j < Z->nfeats; j++)
            Z->sample[t].feat[j] = Z1->sample[s].feat[j];
        Z->sample[t].id     = p;
        t++;
        sample[i]            = sample[n-1];
        n--; nsamples--;
    }

    iftFree(sample);

    return(Z);
}


iftDataSet *iftEliminateAmbiguousSamples(iftDataSet *Z)
{
    iftDataSet *Z1;
    int   i, s, t, nsamples=Z->nsamples;
    char *sample = iftAllocCharArray(Z->nsamples);
    float dist;

    if (!iftIsSupervisedDataSet(Z))
        iftError("There is no class information", "iftEliminateAmbiguousSamples");

    if (!iftIsTestingDataSet(Z))
        iftError("It requires a testing dataset", "iftEliminateAmbiguousSamples");

    if (iftIsCentralizedDataSet(Z))
        iftError("Features have been processed", "iftEliminateAmbiguousSamples");

    for (s=0; s < Z->nsamples-1; s++)
        for (t=s+1; t < Z->nsamples; t++){
            if (Z->sample[s].truelabel != Z->sample[t].truelabel){
                if (iftDist == NULL) {
                    dist = Z->iftArcWeight(Z->sample[s].feat,Z->sample[t].feat,Z->alpha,Z->nfeats);
                }else{
                    dist = iftDist->distance_table[s][t];
                }

                if (dist==0.0){
                    sample[s]=1; sample[t]=1;
                    nsamples-=2;
                }
            }
        }

    Z1 = iftCreateDataSet(nsamples,Z->nfeats);
    Z1->iftArcWeight    = Z->iftArcWeight;
    Z1->function_number = Z->function_number;
    Z1->ngroups         = 0;
    Z1->nclasses        = 0;
    iftCopyRefData(Z1, Z->ref_data, Z->ref_data_type);
    // Z1->ref_data        = Z->ref_data;
    // Z1->ref_data_type   = Z->ref_data_type;

    t = 0;
    for (s=0; s < Z->nsamples; s++) {
        if (sample[s]==0){
            Z1->sample[t].truelabel     = Z->sample[s].truelabel;
            Z1->sample[t].label     = Z->sample[s].label;
            Z1->sample[t].id        = Z->sample[s].id;
            if (Z1->sample[t].truelabel > Z1->nclasses)
                Z1->nclasses = Z1->sample[t].truelabel;
            if (Z1->sample[t].label > Z1->ngroups)
                Z1->ngroups = Z1->sample[t].label;
            for (i=0; i < Z1->nfeats; i++)
                Z1->sample[t].feat[i]=Z->sample[s].feat[i];
            iftSetSampleStatus(&Z1->sample[t], IFT_UNKNOWN);
            Z1->sample[t].weight    = Z->sample[s].weight;
            t++;
        }
    }

    for (i=0; i < Z1->nfeats; i++)
        Z1->alpha[i]=Z->alpha[i];

    Z1->fsp = iftCopyFeatSpaceParam(Z);

    iftFree(sample);

    if (Z1->nclasses < Z->nclasses)
        iftWarning("Classes have been eliminated","iftEliminateAmbiguousSamples");

    return(Z1);
}

iftFeatSpaceParam iftCopyFeatSpaceParam(const iftDataSet *Z)
{
    iftFeatSpaceParam fsp;
    int i,j;

    fsp.mean = fsp.stdev = NULL;
    fsp.R = NULL;
    fsp.W = NULL;
    fsp.w = NULL;
    fsp.nfeats = 0;
    fsp.ncomps = 0;

    if (Z->fsp.mean != NULL){
        fsp.mean = iftAllocFloatArray(Z->fsp.nfeats);
        for (i=0; i < Z->fsp.nfeats; i++)
            fsp.mean[i] = Z->fsp.mean[i];
        fsp.nfeats = Z->fsp.nfeats;
    }

    if (Z->fsp.stdev != NULL){
        fsp.stdev = iftAllocFloatArray(Z->fsp.nfeats);
        for (i=0; i < Z->fsp.nfeats; i++)
            fsp.stdev[i] = Z->fsp.stdev[i];
    }

    if (Z->fsp.R != NULL){
        fsp.R = iftCreateMatrix(Z->fsp.R->ncols,Z->fsp.R->nrows);
        fsp.ncomps = Z->fsp.ncomps;
        for (j=0; j < Z->fsp.R->nrows; j++)
            for (i=0; i < Z->fsp.R->ncols; i++)
                fsp.R->val[iftGetMatrixIndex(fsp.R,i,j)] = Z->fsp.R->val[iftGetMatrixIndex(Z->fsp.R,i,j)];
    }

    if (Z->fsp.W != NULL){
        fsp.W = iftCreateMatrix(Z->fsp.W->ncols,Z->fsp.W->nrows);
        fsp.ncomps = Z->fsp.ncomps;
        for (j=0; j < Z->fsp.W->nrows; j++)
            for (i=0; i < Z->fsp.W->ncols; i++)
                fsp.W->val[iftGetMatrixIndex(fsp.W,i,j)] = Z->fsp.W->val[iftGetMatrixIndex(Z->fsp.W,i,j)];
    }

    if (Z->fsp.w != NULL){
        fsp.w = iftAllocCharArray(Z->fsp.ncomps);
        for (i=0; i < Z->fsp.ncomps; i++)
            fsp.w[i] = Z->fsp.w[i];
    }

    return(fsp);
}


iftDataSet *iftCopyDataSet(const iftDataSet *Z, bool copy_feats) {
    iftDataSet *Zcopy = iftCreateDataSetAux(Z->nsamples, Z->nfeats,false);

    Zcopy->nclasses        = Z->nclasses;
    Zcopy->ngroups         = Z->ngroups;
    Zcopy->nfeats          = Z->nfeats;
    Zcopy->iftArcWeight    = Z->iftArcWeight;
    Zcopy->function_number = Z->function_number;
    iftCopyRefData(Zcopy, Z->ref_data, Z->ref_data_type);
    Zcopy->ntrainsamples   = Z->ntrainsamples;

    if(copy_feats){
        Zcopy->data = iftCopyMatrix(Z->data);
    }else{
        Zcopy->data = iftCreateMatrixPointer(Z->data->val,Z->nfeats,Z->nsamples);
    }

    for (int s = 0; s < Z->nsamples; s++) {
        Zcopy->sample[s].id        = Z->sample[s].id;
        Zcopy->sample[s].truelabel = Z->sample[s].truelabel;
        Zcopy->sample[s].isSupervised = Z->sample[s].isSupervised;
        Zcopy->sample[s].label     = Z->sample[s].label;
        Zcopy->sample[s].group     = Z->sample[s].group;
        iftSetSampleStatus(&Zcopy->sample[s], Z->sample[s].status);
        Zcopy->sample[s].weight = Z->sample[s].weight;
        Zcopy->sample[s].feat = &(Zcopy->data->val[s*Zcopy->nfeats]);
        Zcopy->sample[s].isLabelPropagated = Z->sample[s].isLabelPropagated;
        Zcopy->sample[s].isSupervised = Z->sample[s].isSupervised;
    }
    Zcopy->ntrainsamples = Z->ntrainsamples;
    Zcopy->fsp = iftCopyFeatSpaceParam(Z);

    for (int i = 0; i < Zcopy->nfeats; i++)
        Zcopy->alpha[i] = Z->alpha[i];

    if(Z->projection != NULL){
        if(Zcopy->projection == NULL){
            Zcopy->projection = iftCreateDoubleMatrix(Z->projection->ncols,Z->projection->nrows);
        }
        for (int k = 0; k < Z->projection->n; ++k) {
            Zcopy->projection->val[k] = Z->projection->val[k];
        }
    }


    return(Zcopy);
}


// Check if all DataSets from a list have the same number of classes, samples and feats
void iftValidateDataSets(iftDataSet **Zs, int num_Z) {
    int nclasses = 0, nsamples = 0, nfeats = 0;
    char msg[200];

    // First iteration out of loop to get the (nclasses, nsamples, nfeats) for dataset validation
    if (num_Z >= 1) {
        nclasses = Zs[0]->nclasses;
        nfeats   = Zs[0]->nfeats;
        nsamples = Zs[0]->nsamples;
    } else
        iftError("Number of Datasets < 1", "iftValidateDataSets");


    // Dataset Validation Loop
    for (int i = 1; i < num_Z; i++) {
        if (nclasses != Zs[i]->nclasses) {
            sprintf(msg, "IFT_ERROR: Number of Classes is different between datasets: Zs[%d] --- Zs[%d]\n" \
            "%d != %d", 0, i, nclasses, Zs[i]->nclasses);
            iftError(msg, "validateDatasets");
        } else if (nsamples != Zs[i]->nsamples) {
            sprintf(msg, "IFT_ERROR: Number of Samples is different between datasets: Zs[%d] --- Zs[%d]\n" \
            "%d != %d", 0, i, nsamples, Zs[i]->nsamples);
            iftError(msg, "validateDatasets");
        } else if (nfeats != Zs[i]->nfeats) {
            sprintf(msg, "IFT_ERROR: Number of Features is different between datasets: Zs[%d] --- Zs[%d]\n" \
            "%d != %d", 0, i, nfeats, Zs[i]->nfeats);
            iftError(msg, "validateDatasets");
        }
    }
}


/*
* In the output rotation matrix, the eigenvectors are the columns
*/
iftMatrix *iftRotationMatrixByPCA(iftDataSet *Z) {
    iftMatrix *A, *U, *S, *Vt, *V;

    iftDataSet* Ztrain = iftExtractSamples(Z, IFT_TRAIN);

    if(Ztrain->nsamples<=0) {
        iftError("Dataset must have training samples", "iftRotationMatrixByPCA");
    }

    if (!iftIsCentralizedDataSet(Ztrain))
        iftError("Dataset must be centralized", "iftRotationMatrixByPCA");

    /* Compute Principal Components */
    if (Z->nsamples < Z->nfeats){
        A = iftDataSetToFeatureMatrix(Ztrain);
    }
    else{
        A = iftDatasetCovarianceMatrix(Ztrain);
    }

    iftDestroyDataSet(&Ztrain);

    iftSingleValueDecomp(A, &U, &S, &Vt);
    iftDestroyMatrix(&A);
    iftDestroyMatrix(&S);
    iftDestroyMatrix(&U);

    V = iftTransposeMatrix(Vt);
    iftDestroyMatrix(&Vt);

    return V;
}

iftMatrix *iftRotationMatrixAndSingularValuesByPCA(iftDataSet *Z, iftMatrix **S) {
    iftMatrix *A, *U, *S_aux, *Vt, *V;

    if (!iftIsCentralizedDataSet(Z))
        iftError("Dataset must be centralized", "iftRotationMatrixByPCA");

    /* Compute Principal Components */
    if (Z->nsamples < Z->nfeats)
        A = iftDataSetToFeatureMatrix(Z);
    else
        A = iftDatasetCovarianceMatrix(Z);

    iftSingleValueDecomp(A, &U, &S_aux, &Vt);
    if (S == NULL)
        iftDestroyMatrix(&S_aux);
    else *S = S_aux;

    iftDestroyMatrix(&A);
    iftDestroyMatrix(&U);

    V = iftTransposeMatrix(Vt);
    iftDestroyMatrix(&Vt);

    return V;
}

iftDataSet *iftAlignDataSetByPCA(iftDataSet *Z)
{
    int         s,i,j;
    iftDataSet *Z1, *Z2;
    iftMatrix  *V;

    /* Centralize DataSet */

    iftSetStatus(Z,IFT_TRAIN);
    Z1 = iftCentralizeDataSet(Z);

    /* Compute rotation matrix according to PCA */

    V  = iftRotationMatrixByPCA(Z1);

    /* Rotate Dataset */

    Z2 = iftCreateDataSet(Z1->nsamples,Z1->nfeats);
    Z2->iftArcWeight      = Z1->iftArcWeight;
    Z2->function_number   = Z1->function_number;
    Z2->nclasses          = Z1->nclasses;
    iftCopyRefData(Z2, Z1->ref_data, Z1->ref_data_type);
    // Z2->ref_data          = Z1->ref_data;
    // Z2->ref_data_type     = Z1->ref_data_type;

    for (s=0; s < Z2->nsamples; s++) {
        Z2->sample[s].truelabel  = Z1->sample[s].truelabel;
        iftSetSampleStatus(&Z2->sample[s], Z1->sample[s].status);
        Z2->sample[s].id         = Z1->sample[s].id;

        /* Rotate Feature Space */

        for (j=0; j < Z2->nfeats; j++) {
            Z2->sample[s].feat[j]=0;
            for (i=0; i < V->nrows; i++)
                Z2->sample[s].feat[j] += Z1->sample[s].feat[i]*V->val[iftGetMatrixIndex(V,j,i)];

        }
    }

    Z2->ntrainsamples  = Z1->ntrainsamples;

    /* Copy feature space parameters */

    Z2->fsp    = iftCopyFeatSpaceParam(Z1);
    Z2->fsp.R  = iftCopyMatrix(V);
    Z2->fsp.ncomps = Z1->nfeats;
    for (i=0; i < Z2->nfeats; i++){
        Z2->alpha[i]=1.0;
    }

    iftDestroyMatrix(&V);
    iftDestroyDataSet(&Z1);

    return(Z2);
}

iftDataSet *iftNormalizeTestDataSet(iftDataSet *Z1, iftDataSet *Z2)
{
    iftDataSet *Zn;
    int s,i;

    if (!iftIsNormalizedDataSet(Z1))
        iftError("It requires a normalized data set as reference", "iftNormalizeTestDataSet");

    if (iftIsCentralizedDataSet(Z2))
        iftError("It requires an unprocessed data set", "iftNormalizeTestDataSet");

    Zn   = iftCopyDataSet(Z2, true);

    for (s=0; s < Zn->nsamples; s++) { /* normalize testing set */
        for (i=0; i < Zn->nfeats; i++) {
            Zn->sample[s].feat[i] = Zn->sample[s].feat[i]-Z1->fsp.mean[i];
            if (Z1->fsp.stdev[i] > IFT_EPSILON) {
                Zn->sample[s].feat[i] /= Z1->fsp.stdev[i];
            }
        }
    }

    Zn->fsp = iftCopyFeatSpaceParam(Z1);

    return(Zn);
}


iftDataSet *iftNormalizeTestDataSet2(iftDataSet *Z, iftFeatSpaceParam fsp) {
    iftDataSet *Zn;
    char msg[200];

    if ((iftIsCentralizedDataSet(Z)) || (!iftIsTestingDataSet(Z)))
        iftError("It requires an unprocessed testing set", "iftNormalizeTestDataSet2");

    if (fsp.nfeats != Z->nfeats) {
        sprintf(msg, "Number is different between Dataset and FeatSpaceParam - (%d - %d)",
                fsp.nfeats, Z->nfeats);
        iftError(msg, "iftNormalizeTestDataSet2");
    }

    Zn = iftCopyDataSet(Z, true);

    for (int s = 0; s < Zn->nsamples; s++) { /* normalize testing set */
        for (int i = 0; i < Zn->nfeats; i++) {
            Zn->sample[s].feat[i] = Zn->sample[s].feat[i] - fsp.mean[i];
            if (fsp.stdev[i] > IFT_EPSILON)
                Zn->sample[s].feat[i] /= fsp.stdev[i];
        }
    }
    Zn->fsp=fsp;


    return (Zn);
}


iftDataSet *iftCentralizeTestDataSet(iftDataSet *Z1, iftDataSet *Z2)
{
    iftDataSet *Zn;
    int s,i;

    if (!iftIsCentralizedDataSet(Z1))
        iftError("It requires a centralized data set reference", "iftCentralizeTestDataSet");

    if (iftIsCentralizedDataSet(Z2))
        iftError("It requires an unprocessed testing set", "iftCentralizeTestDataSet");

    if (Z2->nfeats != Z1->fsp.nfeats)
        iftError("Number of Feats is incompatible between Z1.fsp and Z2", "iftCentralizeTestDataSet");

    Zn   = iftCopyDataSet(Z2, true);

    for (s=0; s < Zn->nsamples; s++) { /* normalize testing set */
        for (i=0; i < Zn->nfeats; i++) {
            Zn->sample[s].feat[i] = Zn->sample[s].feat[i]-Z1->fsp.mean[i];
        }
    }

    Zn->fsp = iftCopyFeatSpaceParam(Z1);

    return(Zn);
}



iftFeatSpaceParam iftCreateFeatSpaceParam() {
    iftFeatSpaceParam fsp;

    fsp.mean = fsp.stdev = NULL;
    fsp.R = NULL;
    fsp.W = NULL;
    fsp.w = NULL;
    fsp.nfeats = fsp.ncomps = 0;

    return fsp;
}


/*
* Read only the Rotation Matrix, nfeats, ncomps and mean
*/
iftFeatSpaceParam iftReadPCAMatrix(char *filename) {
    iftFeatSpaceParam fsp = iftCreateFeatSpaceParam();
    int ncols, nrows;
    FILE *fp = NULL;

    fp = fopen(filename, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadPCAMatrix", filename);

    if (fread(&fsp.nfeats, sizeof(int), 1, fp) != 1)
        iftError("Reading Error", "iftReadPCAMatrix");
    if (fread(&(fsp.ncomps), sizeof(int), 1, fp) != 1)
        iftError("Reading error", "iftReadPCAMatrix");

    fsp.mean = iftAllocFloatArray(fsp.nfeats);
    if (fread(fsp.mean, sizeof(float), fsp.nfeats, fp) != fsp.nfeats)
        iftError("Reading error", "iftReadPCAMatrix");

    // Read Rotation Matrix
    if (fread(&ncols, sizeof(int), 1, fp) != 1)
        iftError("Reading error", "iftReadPCAMatrix");
    if (fread(&nrows, sizeof(int), 1, fp) != 1)
        iftError("Reading error", "iftReadPCAMatrix");

    fsp.R = iftCreateMatrix(ncols, nrows);
    if (fread(fsp.R->val, sizeof(double), fsp.R->n, fp) != fsp.R->n)
        iftError("Reading error", "iftReadPCAMatrix");

    fclose(fp);

    return fsp;
}

/*
* Write only the Rotation Matrix, nfeats, ncomps and mean
*/
void iftWritePCAMatrix(iftFeatSpaceParam fsp, char *filename) {
    FILE *fp = NULL;

    fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWritePCAMatrix", filename);

    if (fwrite(&fsp.nfeats, sizeof(int), 1, fp) != 1)
        iftError("Writing Error", "iftWritePCAMatrix");
    if (fwrite(&(fsp.ncomps), sizeof(int), 1, fp) != 1)
        iftError("Writing error", "iftWritePCAMatrix");
    if (fwrite(fsp.mean, sizeof(float), fsp.nfeats, fp) != fsp.nfeats)
        iftError("Writing error", "iftWritePCAMatrix");

    // Write Rotation Matrix
    if (fwrite(&fsp.R->ncols, sizeof(int), 1, fp) != 1)
        iftError("Writing error", "iftWritePCAMatrix");
    if (fwrite(&fsp.R->nrows, sizeof(int), 1, fp) != 1)
        iftError("Writing error", "iftWritePCAMatrix");
    if (fwrite(fsp.R->val, sizeof(double), fsp.R->n, fp) != fsp.R->n)
        iftError("Writing error", "iftWritePCAMatrix");

    fclose(fp);
}



iftDataSet *iftTransformTestDataSetByPCA(iftDataSet *Z1, iftDataSet *Z2) {
    iftDataSet *Zp;

    if ((!iftIsTransformedDataSetByPCA(Z1)) || (!iftIsTrainingDataSet(Z1)))
        iftError("It requires a transformed training set by PCA", "iftTransformTestDataSetByPCA");

    if ((!iftIsCentralizedDataSet(Z2)) || (!iftIsTestingDataSet(Z2)))
        iftError("It requires a centralized testing set", "iftTransformTestDataSetByPCA");

    if (Z1->fsp.R->nrows != Z2->nfeats) // Reminder: The eigen-vectors are the columns in Z1->fsp.R
        iftError("Incompatible projection", "iftTransformTestDataSetByPCA");

    /* Create projected testing set */
    Zp = iftCreateDataSet(Z2->nsamples, Z1->nfeats);
    Zp->nclasses      = Z2->nclasses;
    Zp->fsp           = iftCopyFeatSpaceParam(Z1);
    iftCopyRefData(Zp, Z2->ref_data, Z2->ref_data_type);
    // Zp->ref_data      = Z2->ref_data;
    // Zp->ref_data_type = Z2->ref_data_type;

    /* compute projection */
#pragma omp parallel for
    for (int s = 0; s < Zp->nsamples; s++) {
        float *feat = iftAllocFloatArray(Z1->fsp.ncomps);

        Zp->sample[s].truelabel = Z2->sample[s].truelabel;

        for (int j = 0; j < Z1->fsp.ncomps; j++) {
            feat[j] = 0.0;
            for (int i = 0; i < Z1->fsp.R->nrows; i++)
                feat[j] += Z2->sample[s].feat[i]
                           * Z1->fsp.R->val[iftGetMatrixIndex(Z1->fsp.R,j,i)];
        }

        if (Z1->fsp.w != NULL) { /* Sup PCA */
            int j = 0;
            for (int i = 0; i < Z1->fsp.ncomps; i++) {
                if (Z1->fsp.w[i]) {
                    Zp->sample[s].feat[j] = feat[i];
                    j++;
                }
            }
        } else { /* PCA */
            for (int j = 0; j < Z1->fsp.ncomps; j++) {
                Zp->sample[s].feat[j] = feat[j];
            }
        }

        iftSetSampleStatus(&Zp->sample[s], IFT_TEST);
        Zp->sample[s].id = Z2->sample[s].id;
        iftFree(feat);
    }

    return (Zp);
}


/*
* m = number of samples
* n = number of feats
* D = min(m, n)
* c = number_of_components
*
* Z1 = (m, n)
* Z1->fsp.R = (n, D)
* Z2 = (m, c)
*
* Rt = (D, n)
*
* Zout = <Z1, Rt> + Z1->fsp.mean
*
*/
iftDataSet *iftInverseTransformTestDataSetByPCA(iftDataSet *Z1, iftDataSet *Z2) {
    iftDataSet *Zinv;

    if ((!iftIsTransformedDataSetByPCA(Z1)) || (!iftIsTrainingDataSet(Z1)))
        iftError("It requires a transformed training set by PCA", "iftInverseTransformTestDataSetByPCA");

    iftMatrix *Rt = iftTransposeMatrix(Z1->fsp.R); // The eigen-vectors became the rows of the matrix Rt

    if ((Z1->fsp.ncomps != Z2->nfeats)) // Reminder: The eigen-vectors became the rows of the matrix Rt
        iftError("Incompatible projection", "iftInverseTransformTestDataSetByPCA");

    /* Create projected testing set */
    Zinv                = iftCreateDataSet(Z2->nsamples, Rt->ncols);
    Zinv->nclasses      = Z2->nclasses;
    Zinv->fsp           = iftCopyFeatSpaceParam(Z1);
    iftCopyRefData(Zinv, Z2->ref_data, Z2->ref_data_type);
    // Zinv->ref_data      = Z2->ref_data;
    // Zinv->ref_data_type = Z2->ref_data_type;

    /* compute projection */
#pragma omp parallel for
    for (int s = 0; s < Zinv->nsamples; s++) {
        float *feat = iftAllocFloatArray(Rt->ncols);

        Zinv->sample[s].truelabel = Z2->sample[s].truelabel;

        for (int j = 0; j < Rt->ncols; j++) {
            feat[j] = 0.0;
            for (int i = 0; i < Z1->fsp.ncomps; i++)
                feat[j] += Z2->sample[s].feat[i] * Rt->val[iftGetMatrixIndex(Rt, j, i)];
        }

        if (Z1->fsp.w != NULL) { /* Sup PCA - TODO: CHECK THE CORRECTNESS */
            int j = 0;
            for (int i = 0; i < Rt->ncols; i++) {
                if (Z1->fsp.w[i]) {
                    Zinv->sample[s].feat[j] = feat[i] + Z1->fsp.mean[i];
                    j++;
                }
            }
        } else { /* PCA */
            for (int j = 0; j < Rt->ncols; j++) {
                // Inverse Transformation requires to sum up the feature mean of the Rotation matrix
                Zinv->sample[s].feat[j] = feat[j] + Z1->fsp.mean[j];
            }
        }

        iftSetSampleStatus(&Zinv->sample[s], IFT_TEST);
        Zinv->sample[s].id = Z2->sample[s].id;
        iftFree(feat);
    }

    return Zinv;
}


void iftMultDataSetByScalar(iftDataSet *Z, float scalar)
{
    int nsamples = Z->nsamples;
    int nfeats = Z->nfeats;

    for (int s=0; s < nsamples; s++)
        for (int f=0; f < nfeats; f++)
            Z->sample[s].feat[f] *= scalar;
}


iftDataSet *iftTransFeatSpaceByPCA(iftDataSet *Z, int num_of_comps)
{
    int i;
    iftDataSet *Zpca;
    iftMatrix  *V;

    if ((Z->ntrainsamples<=0)||(!iftIsCentralizedDataSet(Z)))
        iftError("It requires a centralized training set", "iftTransFeatSpaceByPCA");


    if ((num_of_comps > Z->nfeats)||(num_of_comps > Z->nsamples))
        iftError("There is no need to reduce feature space", "iftTransFeatSpaceByPCA");

    /* Compute rotation matrix according to PCA */

    V  = iftRotationMatrixByPCA(Z);

    if(num_of_comps<V->ncols) {
        iftMatrix* sub = iftSubMatrix(V, 0, V->nrows, 0, num_of_comps);
        iftSwap(V, sub);
        iftDestroyMatrix(&sub);
    }

    /* Rotate dataset and select the principal componentes */
    Zpca = iftDataSetTransformBasis(Z, V);

    /* Copy feature space parameters */
    Zpca->fsp.ncomps = num_of_comps;
    for (i=0; i < Zpca->nfeats; i++){
        Zpca->alpha[i]=1.0;
    }

    iftDestroyMatrix(&V);

    return(Zpca);
}


iftDataSet *iftTransFeatSpaceBySupPCA(iftDataSet *Z, int num_of_comps)
{
    iftDataSet *Zt,*Zp;
    int s,t,c,a,b,i,j,*index,ncomps;
    typedef struct ift_class {
        float *fmin, *fmax;
        int size;
    } iftClass;
    iftClass *classes;
    float *tot_overlap, overlap;
    float fmean_a, fmean_b;


    if ((Z->ntrainsamples<=0)||(!iftIsCentralizedDataSet(Z)))
        iftError("It requires a centralized training set", "iftTransFeatSpaceBySupPCA");


    if ((num_of_comps > Z->nfeats)||(num_of_comps > Z->nsamples)) {
        iftError("There is no need to reduce feature space", "iftTransFeatSpaceBySupPCA");
    }

    /* Compute PCA with higher number of components */

    ncomps = iftMin(iftMin(Z->nfeats, Z->nsamples), num_of_comps * 2);

    Zp = iftTransFeatSpaceByPCA(Z,ncomps);

    /* Find minimum and maximum feature values of each class */

    classes = (iftClass *) iftAlloc(Zp->nclasses,sizeof(iftClass));
    for (c=1; c <= Zp->nclasses; c++) {
        classes[c-1].fmin = iftAllocFloatArray(Zp->nfeats);
        classes[c-1].fmax = iftAllocFloatArray(Zp->nfeats);
        classes[c-1].size = 0;
        for (i=0; i < Zp->nfeats; i++) {
            classes[c-1].fmin[i]    = IFT_INFINITY_FLT;
            classes[c-1].fmax[i]    = IFT_INFINITY_FLT_NEG;
        }
    }

    for (s=0; s < Zp->nsamples; s++) {
        if(iftHasSampleStatus(Zp->sample[s], IFT_TRAIN)) {
            for (i = 0; i < Zp->nfeats; i++) {
                if (Zp->sample[s].feat[i] < classes[Zp->sample[s].truelabel - 1].fmin[i])
                    classes[Zp->sample[s].truelabel - 1].fmin[i] = Zp->sample[s].feat[i];
                if (Zp->sample[s].feat[i] > classes[Zp->sample[s].truelabel - 1].fmax[i])
                    classes[Zp->sample[s].truelabel - 1].fmax[i] = Zp->sample[s].feat[i];
            }
            classes[Zp->sample[s].truelabel - 1].size++;
        }
    }

    /* Compute the num_of_comps (features) that best separate distinct classes */
    tot_overlap = iftAllocFloatArray(Zp->nfeats);
    index       = iftAllocIntArray(Zp->nfeats);

    for (i=0; i < Zp->nfeats; i++) {
        tot_overlap[i]=0;
        index[i]=i;
    }

    for (s=0; s < Zp->nsamples; s++) {
        if (iftHasSampleStatus(Zp->sample[s], IFT_TRAIN)) {
            a = Zp->sample[s].truelabel;
            for (t = s + 1; t < Zp->nsamples; t++) {
                if (iftHasSampleStatus(Zp->sample[t], IFT_TRAIN)) {
                    b = Zp->sample[t].truelabel;
                    if (a != b) {
                        for (i = 0; i < Zp->nfeats; i++) {
                            fmean_a = (classes[a - 1].fmin[i] + classes[a - 1].fmax[i]) / 2.0;
                            fmean_b = (classes[b - 1].fmin[i] + classes[b - 1].fmax[i]) / 2.0;
                            overlap = 0;
                            if (fmean_a <= fmean_b) {
                                if ((Zp->sample[s].feat[i] >= classes[b - 1].fmin[i]) &&
                                    (Zp->sample[s].feat[i] <= classes[a - 1].fmax[i]))
                                    overlap += 1;//(class[a-1].fmax[i] - Zp->sample[s].feat[i])/(class[a-1].fmax[i] - class[b-1].fmin[i])*1.0/class[a-1].size;
                                if ((Zp->sample[t].feat[i] >= classes[b - 1].fmin[i]) &&
                                    (Zp->sample[t].feat[i] <= classes[a - 1].fmax[i]))
                                    overlap += 1;//(-class[b-1].fmin[i] + Zp->sample[t].feat[i])/(class[a-1].fmax[i] - class[b-1].fmin[i])*1.0/class[b-1].size;
                            } else {
                                if ((Zp->sample[s].feat[i] >= classes[a - 1].fmin[i]) &&
                                    (Zp->sample[s].feat[i] <= classes[b - 1].fmax[i]))
                                    overlap += 1;//(Zp->sample[s].feat[i] - class[a-1].fmin[i])/(class[b-1].fmax[i] - class[a-1].fmin[i])*1.0/class[a-1].size;
                                if ((Zp->sample[t].feat[i] >= classes[a - 1].fmin[i]) &&
                                    (Zp->sample[t].feat[i] <= classes[b - 1].fmax[i]))
                                    overlap += 1;//(-Zp->sample[t].feat[i] + class[b-1].fmax[i])/(class[b-1].fmax[i] - class[a-1].fmin[i])*1.0/class[b-1].size;
                            }
                            tot_overlap[i] += overlap;
                        }
                    }
                }
            }
        }
    }

    iftFQuickSort(tot_overlap, index, 0, Zp->nfeats-1, IFT_INCREASING);

    /* Free memory */
    for (c=1; c <= Zp->nclasses; c++) {
        iftFree(classes[c-1].fmin);
        iftFree(classes[c-1].fmax);
    }
    iftFree(classes);
    iftFree(tot_overlap);

    /* Create new dataset with the reduced feature space */

    Zt = iftCreateDataSet(Zp->nsamples,num_of_comps);

    Zt->fsp    = iftCopyFeatSpaceParam(Zp);
    Zt->fsp.w  = iftAllocCharArray(Zp->fsp.ncomps);


    for (i=0; i < Zt->nfeats; i++){
        Zt->fsp.w[index[i]]=1; // select num_of_comps
    }

    for (i=0; i < Zt->nfeats; i++){
        Zt->alpha[i]=1.0;
    }

    iftFree(index);

    Zt->iftArcWeight    = Zp->iftArcWeight;
    Zt->function_number = Zp->function_number;
    Zt->ngroups         = Zp->ngroups;
    Zt->nclasses        = Zp->nclasses;
    iftCopyRefData(Zt, Zp->ref_data, Zp->ref_data_type);
    // Zt->ref_data        = Zp->ref_data;
    // Zt->ref_data_type   = Zp->ref_data_type;

    for (s=0; s < Zp->nsamples; s++) {
        Zt->sample[s].truelabel  = Zp->sample[s].truelabel;
        for (i=0, j=0; i < Zt->fsp.ncomps; i++)
            if (Zt->fsp.w[i]==1){
                Zt->sample[s].feat[j]=Zp->sample[s].feat[i];
                j++;
            }
        iftSetSampleStatus(&Zt->sample[s], Zp->sample[s].status);
        Zt->sample[s].id        = Zp->sample[s].id;
    }
    Zt->ntrainsamples  = Zp->ntrainsamples;

    iftDestroyDataSet(&Zp);

    return(Zt);
}

iftVector iftPrincipalAxis(iftImage *label)
{
    iftMatrix  *V;
    iftDataSet *Z,*Zc;
    iftVector   paxis;

    Z   = iftObjectToDataSet(label);
    iftSetStatus(Z,IFT_TRAIN);
    Zc  = iftCentralizeDataSet(Z);
    V   = iftRotationMatrixByPCA(Zc);
    paxis.x = V->val[0];
    paxis.y = V->val[1*V->ncols];
    paxis.z = V->val[2*V->ncols];
    iftDestroyMatrix(&V);
    iftDestroyDataSet(&Z);
    iftDestroyDataSet(&Zc);

    return(paxis);
}


// iftArcWeight can be for example iftDistance1.

float iftDistance1(float *f1, float *f2, float *alpha, int n){
    float dist=0.0f;
    for (int i=0; i < n; i++){
        dist += (f1[i]-f2[i])*(f1[i]-f2[i])*alpha[i];
    }

    return(sqrtf(dist));
}

float iftDistance2(float *f1, float *f2, float *alpha, int n){
    return(((float) IFT_MAXWEIGHT * log(iftDistance1(f1, f2, alpha, n) + 1)));
}

float iftDistance3(float *f1, float *f2, float *alpha, int n){
    int i;
    float dist=0.0f;

    for (i=0; i < n; i++)
        dist += powf(fabs(f1[i]-f2[i]),alpha[i]);

    return(dist);
}

float iftDistance4(float *f1, float *f2, float *alpha, int n){
    return(((float) IFT_MAXWEIGHT * log(iftDistance3(f1, f2, alpha, n) + 1)));
}

float iftDistance5(float *f1, float *f2, float *alpha, int n) {
    float weight = iftDistance1(f1,f2,alpha,n);

    return(1.0 - exp(-weight/K_distance5));
}

/* It assumes that the feature vectors have norm one. Use
 iftNormOneDataSet. */

float iftDistance6(float *f1, float *f2, float *alpha, int n)
{
    int i;
    float inner_prod=0;

    for (i=0; i < n; i++) {
        inner_prod += f1[i]*f2[i]*alpha[i];
    }

    return(0.5 - inner_prod/2.0);
}

/* It assumes that the feature vectors have norm one. Use
 iftNormOneDataSet. */

float iftDistance7(float *f1, float *f2, float *alpha, int n)
{
    float dist=0.0f;
    int i, j;

    for (i=0; i < n; i++){
        j = (int)(iftMinkowski->mult_factor * (f2[i] - f1[i])) +
            iftMinkowski->mult_factor;
        dist += iftMinkowski->dist[j]*alpha[i];
    }

    return(dist);
}


float iftDistance8(float *f1, float *f2, float *alpha, int n) {
    float cov = iftCov(f1, f2, n);
    float var_f1 = iftVar(f1, n);
    float var_f2 = iftVar(f2, n);

    float corr_coef = (float) cov/sqrtf(var_f1*var_f2);
    return (float) (1 - fabs(corr_coef));
}

float iftDistance9(float *f1, float *f2, float *alpha, int n) {
    float cov = iftCov(f1, f2, n);
    float var_f1 = iftVar(f1, n);
    float var_f2 = iftVar(f2, n);

    float corr_coef = (float) cov/sqrtf(var_f1*var_f2);
    float mici = var_f1 + var_f2 - sqrtf(((var_f1+var_f2)*(var_f1+var_f2)) - 4*var_f1*var_f2*(1 - (corr_coef*corr_coef)));

    return mici;
}

float iftDistance10(float *f1, float *f2, float *alpha, int n) {
    int i;
    float dist=0.0f, sf1 = 0.0f, sf2 = 0.0f, num;

    for (i = 0; i < n; i++){
        sf1+=f1[i];
        sf2+=f2[i];
    }

    for (i=0; i < n; i++){
        num = f1[i] - f2[i];
        dist += (num*num)/(f1[i]+f2[i]+IFT_EPSILON);
    }
    return dist;
//    return(sqrtf(dist));
}

float iftDistance11(float *f1, float *f2, float *alpha, int n) {
    float dist = 0.0f;
    float num;

    for(n--; n>=0; n--)
    {
        num = sqrt(f1[n])-sqrt(f2[n]);
        dist += num*num;
    }

    // 0.707106781 == 1/sqrt(2)
    return 0.707106781*sqrt(dist);
}

float iftDistance12(float *f1, float *f2, float *alpha, int n) {
    float dist = 0.0f;

    for(n--; n>=0; n--)
    {
        dist += sqrtf(f1[n]*f2[n]);
    }

    return dist;
}

float iftEuclideanDistance(float *f1, float *f2, int n) {
    float dist = 0.0f;
    int i;

    for (i = 0; i < n; i++) {
        dist += (f1[i] - f2[i]) * (f1[i] - f2[i]);
    }

    return (sqrt(dist));
}


float iftSquaredEuclideanDistance(float *f1, float *f2, int n) {
    float dist = 0.0f;
    int i;

    for (i = 0; i < n; i++) {
        dist += (f1[i] - f2[i]) * (f1[i] - f2[i]);
    }

    return (dist);
}

float iftL1Norm(float *f1, float *f2, float *alpha, int n){
    double diff=0.0;
    double distance = 0.0;
    for (int i = 0; i < n; ++i) {
        diff = f1[i]  - f2[i];
        if(diff < 0){
            diff = -diff;
        }
        distance += diff;
    }
    if(n > 0){
        distance /= (n);
    }
    return ((float)distance);
}


/* Look-up table for Minkowski distance computation (iftDistance7) */
void iftSetMinkowskiTable(int mult_factor, float exp_factor)
{
    iftMinkowski         = (iftMinkowskiTable *)
            iftAlloc(1,sizeof(iftMinkowskiTable));
    iftMinkowski->nelems = 2*mult_factor+1;
    iftMinkowski->dist   = (float *) iftAlloc(iftMinkowski->nelems, sizeof(float));
    iftMinkowski->mult_factor = mult_factor;
    if (exp_factor <= IFT_EPSILON)
        iftError("Invalid exponential factor", "iftSetMinkowskiTable");
    iftMinkowski->exp_factor = exp_factor;

    for (int i=0; i < iftMinkowski->nelems; i++) {
        float x = fabs((float)(i-iftMinkowski->mult_factor)/(float)iftMinkowski->mult_factor);
        iftMinkowski->dist[i] = powf(x,iftMinkowski->exp_factor);
    }
}

void  iftDestroyMinkowskiTable()
{
    if (iftMinkowski != NULL) {
        iftFree(iftMinkowski->dist);
        iftFree(iftMinkowski);
    }
}



/* --------------------- Distance Table ------------------- */


iftDistanceTable *iftCreateDistanceTable(int nsamples1, int nsamples2)
{
    int s;
    iftDistanceTable *dt=(iftDistanceTable *)iftAlloc(1,sizeof(iftDistanceTable));

    dt->distance_table = (float **)iftAlloc(nsamples1,sizeof(float *));
    for (s=0; s < nsamples1; s++)
        dt->distance_table[s] = iftAllocFloatArray(nsamples2);
    dt->nsamples1 = nsamples1;
    dt->nsamples2 = nsamples2;

    return(dt);
}

void  iftDestroyDistanceTable(iftDistanceTable **dt)
{
    int s;
    iftDistanceTable *aux=*dt;

    if (aux != NULL) {
        for (s=0; s < aux->nsamples1; s++) {
            iftFree(aux->distance_table[s]);
        }
        iftFree(aux->distance_table);
        iftFree(aux);
        *dt = NULL;
    }
}

iftDistanceTable *iftReadDistanceTable(char *filename)
{
    int   s, nsamples1, nsamples2;
    FILE  *fp=fopen(filename,"rb");
    iftDistanceTable *dt;

    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadDistanceTable", filename);

    if (fread(&nsamples1,sizeof(int),1,fp)!=1) iftError("Reading error", "iftReadDistanceTable");
    if (fread(&nsamples2,sizeof(int),1,fp)!=1) iftError("Reading error", "iftReadDistanceTable");

    dt = iftCreateDistanceTable(nsamples1, nsamples2);
    for (s=0; s < nsamples1; s++) {
        if(fread(dt->distance_table[s],sizeof(float),nsamples2,fp)!=nsamples2)
            iftError("Reading error", "iftReadDistanceTable");
    }
    fclose(fp);

    return(dt);
}

iftDistanceTable *iftCompDistanceTable(const iftDataSet *Z1, const iftDataSet* Z2)
{
    if(Z1->nfeats!=Z2->nfeats) {
        iftError("Different number of features.", "iftCompDistanceTable");
    }

    iftDistanceTable *dt=iftCreateDistanceTable(Z1->nsamples, Z2->nsamples);

    if(Z1==Z2) /*Trick to compute faster the dataset to itself*/ {
        #pragma omp parallel for shared(Z1,Z2,dt)
        for (int s = 0; s < Z1->nsamples-1; s++) {
            for (int t = s+1; t < Z2->nsamples; t++) {
                dt->distance_table[t][s] = dt->distance_table[s][t] =
                        Z1->iftArcWeight(Z1->sample[s].feat,Z2->sample[t].feat,Z1->alpha,Z1->nfeats);
            }
        }
    }
    else {
        #pragma omp parallel for shared(Z1,Z2,dt) schedule(auto)
        for (int s = 0; s < Z1->nsamples; s++) {
            for (int t = 0; t < Z2->nsamples; t++) {
                dt->distance_table[s][t] = Z1->iftArcWeight(Z1->sample[s].feat, Z2->sample[t].feat, Z1->alpha, Z1->nfeats);
            }
        }
    }

    return(dt);
}

void iftPrintDistanceTable(iftDistanceTable *dt){

    fprintf(stdout, "\n");
    for (int r = 0; r < dt->nsamples1; r++) {
        for (int c = 0; c < dt->nsamples2; c++) {
            fprintf(stdout, "%6.5lf ", dt->distance_table[r][c]);
        }
        fprintf(stdout, "\n");
    }
}

void iftWriteDistanceTable(iftDistanceTable *dt, char *filename)
{
    FILE *fp=fopen(filename,"wb");
    int   s;

    fwrite(&dt->nsamples1,sizeof(int),1,fp);
    fwrite(&dt->nsamples2,sizeof(int),1,fp);
    for (s=0; s < dt->nsamples1; s++)
        fwrite(dt->distance_table[s],sizeof(float),dt->nsamples2,fp);

    fclose(fp);
}

void iftSetTrainingSupervoxelsFromSeeds(iftDataSet *dataset, iftImage *label_image, iftLabeledSet* seed){
    dataset->nclasses = 0;

    int region_label,pixel;

    while(seed != NULL) {
        pixel = seed->elem;
        region_label = label_image->val[pixel] - 1;
        dataset->sample[region_label].truelabel  = seed->label;
        iftSetSampleStatus(&dataset->sample[region_label], IFT_TRAIN);

        // Counting classes
        if (seed->label >= dataset->nclasses)
            dataset->nclasses = dataset->sample[region_label].truelabel + 1;

        seed = seed->next;
    }

    dataset->ntrainsamples = 0;
    int i;
    for (i=0; i < dataset->nsamples; i++)
        if (iftHasSampleStatus(dataset->sample[i], IFT_TRAIN))
            dataset->ntrainsamples++;
}



iftDataSet *iftMImageToDataSet(const iftMImage *mimg, const iftImage *label_img, float coord_scale)
{
    if ((label_img != NULL) && !iftIsDomainEqual(mimg, label_img))
        iftError("Multiband Image and Label Image have different domains: (%d, %d, %d) != (%d, %d, %d", "iftMImageToDataSet",
                  mimg->xsize, mimg->ysize, mimg->zsize, label_img->xsize, label_img->ysize, label_img->zsize);

    iftDataSet *Z = NULL;
    int nsamples  = mimg->n;
    int ndim = (coord_scale <= 0 ? 0 : (mimg->zsize > 1 ? 3 : 2));
    int nfeats    = mimg->m + ndim;

    if (label_img == NULL)
    {
        Z = iftCreateDataSet(nsamples, nfeats);
        for (int p = 0; p < mimg->n; p++) {
            for (int b = 0; b < mimg->m; b++)
                Z->sample[p].feat[b] = mimg->val[p][b];
            if (ndim > 0)
            {
                iftVoxel v = iftMGetVoxelCoord(mimg, p);
                Z->sample[p].feat[mimg->m] = v.x * coord_scale;
                Z->sample[p].feat[mimg->m + 1] = v.y * coord_scale;
                if (ndim > 2)
                    Z->sample[p].feat[mimg->m + 2] = v.z * coord_scale;
            }
            Z->sample[p].id = p;
        }
    }
    else {
        iftIntArray *labels = iftGetObjectLabels(label_img);

        if (labels->n == 0)
	       iftError("Empty label image", "iftMImageToDataSet");
        // binary label image (only one object)
        else if (labels->n == 1) {
            nsamples = iftCountTotalObjectSpels(label_img);
            Z        = iftCreateDataSet(nsamples, nfeats);
	
            int s = 0;
            for (int p = 0; p < mimg->n; p++) {
                if (label_img->val[p] != 0) {
                    for (int b = 0; b < mimg->m; b++) {
                        Z->sample[s].feat[b] = mimg->val[p][b];
                    }
                    if (ndim > 0)
                    {
                        iftVoxel v = iftMGetVoxelCoord(mimg, p);
                        Z->sample[s].feat[mimg->m] = v.x * coord_scale;
                        Z->sample[s].feat[mimg->m + 1] = v.y * coord_scale;
                        if (ndim > 2)
                            Z->sample[s].feat[mimg->m + 2] = v.z * coord_scale;
                    }
                    Z->sample[s].id = p;
                    s++;
                }
            }
        }
        // multi-label image
        else {
            nsamples = labels->n; // number of labels
            Z        = iftCreateDataSet(nsamples, nfeats);
	
            // the label values are the dict's keys and their indices are the corresponding dict's values
            iftDict *labels_dict = iftIntArrayToDict(labels->val, labels->n);
            iftPoint *gcs        = iftAlloc(labels->n, sizeof(iftPoint)); // one geometric center per object, zero-initialization
            iftIntArray *sizes   = iftCreateIntArray(labels->n); // one per object, zero-initialization

            for (int p = 0; p < mimg->n; p++) {
                if (label_img->val[p] != 0) {
                    int label = label_img->val[p];
                    int s     = iftGetLongValFromDict(label, labels_dict); // the object's index

                    for (int b = 0; b < mimg->m; b++)
                        Z->sample[s].feat[b] += mimg->val[p][b];
                    
                    iftVoxel u = iftMGetVoxelCoord(mimg, p);
                    gcs[s].x += u.x;
                    gcs[s].y += u.y;
                    gcs[s].z += u.z;
                    sizes->val[s]++;
                }
            }
      
            // for each object, average its feature vector and geometric center
            for (int s = 0; s < Z->nsamples; s++) {
                int label = labels->val[s];
    	  
                for (int b = 0; b < Z->nfeats; b++)
                    Z->sample[s].feat[b] /= sizes->val[s];
    	  
                iftVoxel gc;
                gc.x = iftRound(gcs[s].x / sizes->val[s]);
                gc.y = iftRound(gcs[s].y / sizes->val[s]);
                gc.z = iftRound(gcs[s].z / sizes->val[s]);
    	  
                // the sample's id is the geometric center of the object <label>
                // But, the computed gc might not belong to this object, then the next voxel from
                // the object will be the sample's id
                int id = iftGetVoxelIndex(label_img, gc);
                if (label_img->val[id] != label) { // gc is out of its superpixel
                    id = IFT_NIL;
                    int radii = sqrtf(sizes->val[s] / 2.0) + 1;
                    iftAdjRel *A = (iftIs3DMImage(mimg) ? iftSpheric(radii) : iftCircular(radii));
                    for (int i = 1; i < A->n; i++) {
                        iftVoxel v = iftGetAdjacentVoxel(A, gc, i);
                        if (!iftValidVoxel(label_img, v))
                            continue;
                        id = iftGetVoxelIndex(label_img, v);
                        if (label_img->val[id] == label)
                            break;
                    }
                    iftDestroyAdjRel(&A);
                }
                if (id == IFT_NIL) {
                    char *st = iftAlloc(512, sizeof *st);
                    snprintf(st, 512, "Sample id not found for label %d, this function only guarantee to find the centroids "
                                     "for connected components for efficiency purposes, TO NOT CHANGE THIS", label);
                    iftWarning(st, "iftMImageToDataSet");
                    iftFree(st);
                }
                Z->sample[s].id = id;
                if (ndim > 0)
                {
                    iftVoxel v = iftGetVoxelCoord(label_img, id);
                    Z->sample[s].feat[mimg->m] = v.x * coord_scale;
                    Z->sample[s].feat[mimg->m + 1] = v.y * coord_scale;
                    if (ndim > 2)
                        Z->sample[s].feat[mimg->m + 2] = v.z * coord_scale;
                }
            }
            iftDestroyDict(&labels_dict);
            iftFree(gcs);
            iftDestroyIntArray(&sizes);
        }
      
        iftDestroyIntArray(&labels);
    }

    iftCopyRefData(Z, mimg, IFT_REF_DATA_MIMAGE);
    iftSetStatus(Z,IFT_TEST);
    
    return Z;
}


iftDataSet *iftMImageToDataSetUsingAdjacency(iftMImage *mimg, iftAdjRel *A)
{
    iftDataSet *Z;

    Z=iftCreateDataSet(mimg->n, mimg->m * A->n);

#pragma omp parallel for shared(Z,mimg,A)
    for (int p=0; p < mimg->n; p++) {
        iftVoxel u = iftMGetVoxelCoord(mimg, p);
        for(int b=0; b < mimg->m; b++) {
            for (int i = 0; i < A->n; i += 1) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);
                if (iftMValidVoxel(mimg, v)) {
                    int q = iftMGetVoxelIndex(mimg, v);
                    Z->sample[p].feat[i+b*A->n] = mimg->val[q][b];
                }
            }
        }
        Z->sample[p].id              = p;
    }
    iftCopyRefData(Z, mimg, IFT_REF_DATA_MIMAGE);
    // Z->ref_data       = mimg;
    // Z->ref_data_type  = IFT_REF_DATA_MIMAGE;

    iftSetStatus(Z,IFT_TEST);
    return(Z);
}

iftDataSet* iftMImageToDataSetInRegion(iftMImage* mimg, iftImage *mask)
{
    iftDataSet* Z;
    int p,b,nsamples=0,s;

    for(p = 0; p < mask->n; p++)
        if (mask->val[p])
            nsamples++;

    Z = iftCreateDataSet(nsamples, mimg->m);
    iftCopyRefData(Z, mimg, IFT_REF_DATA_MIMAGE);
    // Z->ref_data       = mimg;
    // Z->ref_data_type  = IFT_REF_DATA_MIMAGE;

    for(p = 0, s = 0; p < mimg->n; p++)
    {
        if (mask->val[p]){
            for(b = 0; b < mimg->m ; b++)
                Z->sample[s].feat[b] = mimg->val[p][b];
            Z->sample[s].id = p;
            s++;
        }
    }

    iftSetStatus(Z,IFT_TEST);
    return Z;
}


void iftImageGTToDataSet(iftImage* imgGT,iftDataSet* Z)
{
    if (Z->nsamples != imgGT->n) {
        char msg[200];
        sprintf(msg,"Number of samples (%d) in dataset and voxels (%d) in GT image are different\n\n",Z->nsamples,imgGT->n);
        iftError(msg, "iftGTToDataset");
    }

    /*todo we can remove this piece of code if we do the call to iftRelabelGrayscaleImage function before calling this
    * function*/
    iftDict *dict=iftCreateDict();
    int current_val=0;

    for(int p = 0; p < imgGT->n; p++) {
        if (!iftDictContainKey(imgGT->val[p],dict,NULL)){
            iftInsertIntoDict(imgGT->val[p],++current_val,dict);
        }
    }

#pragma omp parallel for shared(Z,dict)
    for(int p = 0; p < imgGT->n; p++) {
        Z->sample[p].truelabel = (int)iftGetLongValFromDict(imgGT->val[p],dict);
    }

    Z->nclasses = current_val;

    iftDestroyDict(&dict);
}

iftDataSet *iftMSupervoxelsToDataSet(iftMImage *mimg, iftImage *label)
{
    iftDataSet *Z;
    int r,s,p,b,nregions=iftMaximumValue(label);
    int *size=iftAllocIntArray(nregions);

    if (iftMinimumValue(label)<=0)
        iftError("Minimum label value must be 1", "iftImageCompsToDataSet");

    //iftVerifyImageDomains(mimg,label,"iftMSupervoxelsToDataSet");

    Z=iftCreateDataSet(nregions,mimg->m);
    for (p=0; p < mimg->n; p++) {
        s = label->val[p] - 1;
        size[s]++;
        for(b=0; b < mimg->m;b++)
            Z->sample[s].feat[b] += mimg->val[p][b];
    }

    for(r = 0; r < nregions; r++){
        for(b=0; b < mimg->m; b++)
            Z->sample[r].feat[b] = Z->sample[r].feat[b]/size[r];
    }

    iftFree(size);

    //This may seem redundant, but it is useful for splitting the dataset
    for(r = 0; r < nregions; r++){
        Z->sample[r].id = r + 1;
    }


    iftSetDistanceFunction(Z, 1);
    iftSetStatus(Z,IFT_TEST);

    return(Z);
}

//Features: sum band 1, sum band 2, sum band 3, relative area, npixels
iftDataSet* iftSupervoxelsToMeanSizeDataSet(iftImage* image, iftImage* supervoxel, iftColorSpace colorspace){

    int normalization_value = iftNormalizationValue(iftMaximumValue(image));

    if(!iftIsColorImage(image))
        iftError("Color image expected.", "iftSupervoxelsToMeanSizeDataSet");
    if(colorspace != LAB_CSPACE && colorspace != RGB_CSPACE && colorspace != YCbCr_CSPACE)
        iftError("Invalid colorspace.", "iftSupervoxelsToMeanSizeDataSet");

    iftVerifyImageDomains(image, supervoxel,"iftSupervoxelsToMeanSizeDataSet");

    int nregions = iftMaximumValue(supervoxel);
    int *size = iftAllocIntArray(nregions);

    iftDataSet* dataset = iftCreateDataSet(nregions, 5);

    int p;
    for(p = 0; p < image->n; p++){
        int r = supervoxel->val[p] - 1;
        size[r]++;

        iftColor color;
        color.val[0] = image->val[p];
        color.val[1] = image->Cb[p];
        color.val[2] = image->Cr[p];

        if(colorspace == LAB_CSPACE){
            color = iftYCbCrtoRGB(color,normalization_value);
            iftFColor fcolor = iftRGBtoLab(color,normalization_value);

            dataset->sample[r].feat[0] += fcolor.val[0];
            dataset->sample[r].feat[1] += fcolor.val[1];
            dataset->sample[r].feat[2] += fcolor.val[2];
        }else {
            if (colorspace == RGB_CSPACE)
                color = iftYCbCrtoRGB(color,normalization_value);

            dataset->sample[r].feat[0] += color.val[0];
            dataset->sample[r].feat[1] += color.val[1];
            dataset->sample[r].feat[2] += color.val[2];
        }
    }

    int r;
    for(r = 0; r < nregions; r++){
        if(size[r] == 0)
            iftError("Empty region", "iftRegionsToMeanStdSizeDataset");

        //Relative size
        dataset->sample[r].feat[3] = ((float)size[r])/image->n;

        //Number of pixels
        dataset->sample[r].feat[4] = size[r];

        dataset->sample[r].id = r + 1;
    }

    dataset->iftArcWeight    = iftDistMeanSizeSupervoxel;
    dataset->function_number = IFT_NIL;
    iftCopyRefData(dataset, supervoxel, IFT_REF_DATA_IMAGE);
    // dataset->ref_data        = (void*)supervoxel;
    // dataset->ref_data_type   = IFT_REF_DATA_IMAGE;

    iftFree(size);
    iftSetStatus(dataset,IFT_TEST);

    return dataset;
}

//Features (9): (normalized) mean L, mean a, mean b, std L, std a, std b, skewness L, skewness a, skewness b
iftDataSet* iftSupervoxelsToLabColorMomentsDataset(iftImage* image, iftImage* label_image){

    int normalization_value = iftNormalizationValue(iftMaximumValue(image));

    if(!iftIsColorImage(image))
        iftError("Color image expected.", "iftSupervoxelsToColorMomentsDataset");

    iftVerifyImageDomains(image, label_image,"iftSupervoxelsToColorMomentsDataset");

    int nregions = iftMaximumValue(label_image);
    int *size = iftAllocIntArray(nregions);

    iftDataSet* dataset = iftCreateDataSet(nregions, 9);

    iftMImage *lab_img = iftCreateMImage(image->xsize, image->ysize, image->zsize, 3);

    int p;
    for(p = 0; p < image->n; p++){
        int r = label_image->val[p] - 1;
        if(r < 0)
            continue;

        size[r]++;

        iftColor color;
        color.val[0] = image->val[p];
        color.val[1] = image->Cb[p];
        color.val[2] = image->Cr[p];

        color = iftYCbCrtoRGB(color,normalization_value);
        iftFColor fcolor = iftRGBtoLab(color,normalization_value);

        // L in [0, 99.998337]
        // a in [-86.182236, 98.258614]
        // b in [-107.867744, 94.481682]
        lab_img->val[p][0] = (fcolor.val[0] - 0.)/(99.998337 - 0.);
        lab_img->val[p][1] = (fcolor.val[1] + 86.182236)/(98.258614 + 86.182236);
        lab_img->val[p][2] = (fcolor.val[2] + 107.867744)/(94.481682 + 107.867744);

        dataset->sample[r].feat[0] += lab_img->val[p][0];
        dataset->sample[r].feat[1] += lab_img->val[p][1];
        dataset->sample[r].feat[2] += lab_img->val[p][2];
    }

    int r;
    for(r = 0; r < nregions; r++){
        if(size[r] == 0)
            iftError("Empty region", "iftSupervoxelsToColorMomentsDataset");

        dataset->sample[r].id = r+1;

        //Means
        dataset->sample[r].feat[0] /= size[r];
        dataset->sample[r].feat[1] /= size[r];
        dataset->sample[r].feat[2] /= size[r];
    }

    for(p = 0; p < image->n; p++){
        int r = label_image->val[p] - 1;
        if(r < 0)
            continue;

        dataset->sample[r].feat[3] += pow(dataset->sample[r].feat[0] - lab_img->val[p][0], 2);
        dataset->sample[r].feat[4] += pow(dataset->sample[r].feat[1] - lab_img->val[p][1], 2);
        dataset->sample[r].feat[5] += pow(dataset->sample[r].feat[2] - lab_img->val[p][2], 2);

        dataset->sample[r].feat[6] += pow(dataset->sample[r].feat[0] - lab_img->val[p][0], 3);
        dataset->sample[r].feat[7] += pow(dataset->sample[r].feat[1] - lab_img->val[p][1], 3);
        dataset->sample[r].feat[8] += pow(dataset->sample[r].feat[2] - lab_img->val[p][2], 3);

    }

    for(r = 0; r < nregions; r++){
        //Standard deviation
        dataset->sample[r].feat[3] = sqrt(dataset->sample[r].feat[3]/size[r]);
        dataset->sample[r].feat[4] = sqrt(dataset->sample[r].feat[4]/size[r]);
        dataset->sample[r].feat[5] = sqrt(dataset->sample[r].feat[5]/size[r]);
    }

    for(r = 0; r < nregions; r++){
        //Skewness
        if(iftAlmostZero(pow(dataset->sample[r].feat[3], 3)))
            dataset->sample[r].feat[6] = 0;
        else
            dataset->sample[r].feat[6] = (dataset->sample[r].feat[6]/size[r])/pow(dataset->sample[r].feat[3], 3);

        if(iftAlmostZero(pow(dataset->sample[r].feat[4], 3)))
            dataset->sample[r].feat[7] = 0;
        else
            dataset->sample[r].feat[7] = (dataset->sample[r].feat[7]/size[r])/pow(dataset->sample[r].feat[4], 3);

        if(iftAlmostZero(pow(dataset->sample[r].feat[5], 3)))
            dataset->sample[r].feat[8] = 0;
        else
            dataset->sample[r].feat[8] = (dataset->sample[r].feat[8]/size[r])/pow(dataset->sample[r].feat[5], 3);
    }


    iftSetDistanceFunction(dataset, 1);
    iftCopyRefData(dataset, label_image, IFT_REF_DATA_IMAGE);
    // dataset->ref_data      = (void*)label_image;
    // dataset->ref_data_type = IFT_REF_DATA_IMAGE;

    iftFree(size);
    iftDestroyMImage(&lab_img);

    return dataset;
}

//Features (3*bins_per_band): (normalized histogram) L1 ... L(bins_per_band), a1 ... a(bins_per_band), b1 ... b(bins_per_band)
//Default bins_per_band is 4
iftDataSet* iftSupervoxelsToLabHistogramDataset(iftImage* image, iftImage* label_image, int bins_per_band){

    int normalization_value = iftNormalizationValue(iftMaximumValue(image));

    if(!iftIsColorImage(image))
        iftError("Color image expected.", "iftSupervoxelsToLabHistogramDataset");
    iftVerifyImageDomains(image, label_image,"iftSupervoxelsToLabHistogramDataset");

    int nregions = iftMaximumValue(label_image);
    int *size = iftAllocIntArray(nregions);

    iftDataSet* dataset = iftCreateDataSet(nregions, bins_per_band*3);

    iftMImage *lab_img = iftCreateMImage(image->xsize, image->ysize, image->zsize, 3);

    int p,r;
    for(p = 0; p < image->n; p++){
        int r = label_image->val[p] - 1;
        if(r < 0)
            continue;

        size[r]++;

        iftColor color;
        color.val[0] = image->val[p];
        color.val[1] = image->Cb[p];
        color.val[2] = image->Cr[p];

        color = iftYCbCrtoRGB(color,normalization_value);
        iftFColor fcolor = iftRGBtoLab(color,normalization_value);

        // L in [0, 99.998337]
        // a in [-86.182236, 98.258614]
        // b in [-107.867744, 94.481682]
        lab_img->val[p][0] = (fcolor.val[0] - 0.)/(99.998337 - 0.);
        lab_img->val[p][1] = (fcolor.val[1] + 86.182236)/(98.258614 + 86.182236);
        lab_img->val[p][2] = (fcolor.val[2] + 107.867744)/(94.481682 + 107.867744);

        dataset->sample[r].feat[iftMin((int)(lab_img->val[p][0] * bins_per_band), bins_per_band - 1)]++;
        dataset->sample[r].feat[iftMin((int)(lab_img->val[p][1] * bins_per_band), bins_per_band - 1) + bins_per_band]++;
        dataset->sample[r].feat[iftMin((int)(lab_img->val[p][2] * bins_per_band), bins_per_band - 1) + 2 * bins_per_band]++;
    }

    int i;
    for(r = 0; r < nregions; r++){
        if(size[r] == 0)
            iftError("Empty region", "iftSupervoxelsToLabHistogramDataset");

        dataset->sample[r].id = r+1;

        for(i = 0; i < dataset->nfeats; i++)
            dataset->sample[r].feat[i] /= size[r];
    }

    iftSetDistanceFunction(dataset, 1);
    iftCopyRefData(dataset, label_image, IFT_REF_DATA_IMAGE);
    // dataset->ref_data      = (void*)label_image;
    // dataset->ref_data_type = IFT_REF_DATA_IMAGE;

    iftFree(size);
    iftDestroyMImage(&lab_img);

    return dataset;
}

//Features:
//color: (2*bins_per_band**3) INT1 ... INT(bins_per_band**3), EXT1 ... EXT(bins_per_band**3)
//grayscale: (2*bins_per_band) INT1 ... INT(bins_per_band), EXT1 ... EXT(bins_per_band)
//Default bins_per_band is 4
iftDataSet *iftSupervoxelsToBICDataset(iftImage *img, iftImage *super_img, int bins_per_band,
                                       int max_range, iftIntArray *true_labels) {
    // checkers
    if (img == NULL)
        iftError("Input Image is NULL", "iftSupervoxelsToBICDataset");
    if (super_img == NULL)
        iftError("Supervoxel Image is NULL", "iftSupervoxelsToBICDataset");
    iftVerifyImageDomains(img, super_img, "iftSupervoxelsToBICDataset");
    if (bins_per_band <= 0)
        iftError("Bins per Band %d <= 0", "iftSupervoxelsToBICDataset", bins_per_band);

    int n_regions = iftMaximumValue(super_img);

    if ((true_labels != NULL) && (true_labels->n < n_regions))
        iftError("Number of Input True Labels %lu is not enough for the number of supervoxels %lu\n",
                 "iftSupervoxelsToBICDataset", true_labels->n, n_regions);

    int norm_value;
    if (max_range <=0 )
        norm_value = iftNormalizationValue(iftMaximumValue(img));
    else norm_value = max_range;

    int bins = bins_per_band;

    if (iftIsColorImage(img))
        bins = bins_per_band * bins_per_band * bins_per_band;

    int *size = iftAllocIntArray(n_regions);

    iftImage *quant_img = iftQuantize(img, bins_per_band, norm_value);
    iftDataSet *dataset = iftCreateDataSet(n_regions, bins * 2);

    iftAdjRel *A = NULL;
    if (iftIs3DImage(img))
        A = iftSpheric(1.0);
    else
        A = iftCircular(1.0);

    for (int p = 0; p < img->n; p++) {
        int r = super_img->val[p] - 1;
        if (r < 0)
            continue;
        size[r]++;

        iftVoxel u = iftGetVoxelCoord(img, p);

        int border = 0;
        int i;
        for (i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftValidVoxel(img, v)) {
                int q = iftGetVoxelIndex(img, v);

                if (quant_img->val[p] != quant_img->val[q]) {
                    border = 1;
                    break;
                }
            }
        }

        if (border) {
            dataset->sample[r].feat[quant_img->val[p]]++;
        } else {
            dataset->sample[r].feat[quant_img->val[p] + bins]++;
        }
    }

    //This normalization could be done independently for each histogram, but this choice is deliberate.
    for (int r = 0; r < dataset->nsamples; r++)
        for (int i = 0; i < dataset->nfeats; i++)
            dataset->sample[r].feat[i] = dataset->sample[r].feat[i] / size[r];

    iftFree(size);
    iftDestroyImage(&quant_img);
    iftDestroyAdjRel(&A);

    // assigns the true labels
    if (true_labels == NULL) {
        for (int r = 0; r < dataset->nsamples; r++)
            dataset->sample[r].truelabel = r+1; // label from the supervoxel
    }
    else {
        for (int r = 0; r < dataset->nsamples; r++)
            dataset->sample[r].truelabel = true_labels->val[r]; // label from the supervoxel
    }

    //This is not the usual distance for computing BIC distances. However, this function was implemented with binary classification in mind.
    iftSetDistanceFunction(dataset, 2);
    iftCopyRefData(dataset, super_img, IFT_REF_DATA_IMAGE);
    // dataset->ref_data = (void *) super_img;
    // dataset->ref_data_type  = IFT_REF_DATA_IMAGE;

    return dataset;
}

//This version of LBP is grayscale and rotation invariant
//Features: (59) LBPH1, ..., LBPH59
iftDataSet* iftSupervoxelsToUniformLBP2D(iftImage* image, iftImage* label_image){
    iftVerifyImageDomains(image, label_image,"iftSupervoxelsToUniformLBP");
    if(iftIs3DImage(image))
        iftError("Two dimensional image expected", "iftSupervoxelsToUniformLBP");

    int nregions = iftMaximumValue(label_image);
    int *size = iftAllocIntArray(nregions);

    int nn = 8; //number of neighbors
    iftAdjRel* A = iftCircular(1.5);

    iftDataSet* dataset = iftCreateDataSet(nregions, 3 + (nn * (nn - 1)));

    int p;
    for(p = 0; p < image->n; p++){
        int r = label_image->val[p] - 1;
        if(r < 0)
            continue;
        size[r]++;

        iftVoxel u = iftGetVoxelCoord(image,p);

        int signed_texture[8] = {0};

        int i;
        for(i = 1; i < A->n; i++){
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);
            if(iftValidVoxel(image,v)){
                int q = iftGetVoxelIndex(image,v);

                if(image->val[q] >= image->val[p])
                    signed_texture[i - 1] = 1;
            }
        }

        int changes = 0;
        for(i = 0; i < nn - 1; i++)
            changes += abs(signed_texture[i] - signed_texture[i + 1]);

        int lbp = 0;
        if(changes <= 2){ //Uniform pattern
            int nones = 0;
            int first_one = -1;
            int first_zero = -1;

            for(i = 0; i < nn; i++){
                if(signed_texture[i] == 1){
                    nones++;
                    if(first_one < 0)
                        first_one = i;
                }
                else{
                    if(first_zero < 0)
                        first_zero = i;
                }

            }

            if(nones == 0)
                lbp = 0;
            else if (nones == nn)
                lbp = nn * (nn - 1) + 1; //57
            else{
                int rot_index = 0;
                if(first_one == 0)
                    rot_index = nones - first_zero;
                else
                    rot_index = nn - first_one;

                lbp = 1 + (nones - 1) *nn + rot_index;
            }
        }else{ //Non-uniform pattern
            lbp = nn * (nn - 1) + 2; //58
        }

        if(lbp < 0 || lbp >= dataset->nfeats)
            iftError("Invalid LBP code", "iftSupervoxelsToUniformLBP");

        dataset->sample[r].feat[lbp]++;
    }

    int r,i;
    for(r = 0; r < dataset->nsamples; r++)
        for(i = 0; i < dataset->nfeats; i++)
            dataset->sample[r].feat[i] = dataset->sample[r].feat[i]/size[r];

    iftFree(size);
    iftDestroyAdjRel(&A);


    iftSetDistanceFunction(dataset, 2);
    iftCopyRefData(dataset, label_image, IFT_REF_DATA_IMAGE);
    // dataset->ref_data      = (void*)label_image;
    // dataset->ref_data_type = IFT_REF_DATA_IMAGE;

    return dataset;
}

//This is not the same HOG usually used for image classification.
//This is a simplified version to be used on superpixels instead of fixed windows.
//Default nbins is 16
iftDataSet* iftSupervoxelsToSimplifiedHOG2D(iftImage* image, iftImage* label_image, int nbins){
    iftVerifyImageDomains(image, label_image,"iftSupervoxelsToSimplifiedHOG");
    if(iftIs3DImage(image))
        iftError("Two dimensional image expected", "iftSupervoxelsToSimplifiedHOG");

    int nregions = iftMaximumValue(label_image);
    int *size = iftAllocIntArray(nregions);

    iftKernel* sobel_x = iftSobelXKernel2D();
    iftKernel* sobel_y = iftSobelYKernel2D();

    iftImage *image_x = iftLinearFilter(image, sobel_x);
    iftImage *image_y = iftLinearFilter(image, sobel_y);

    float *thetas = iftAllocFloatArray(nbins);

    iftDataSet *dataset = iftCreateDataSet(nregions, nbins);

    int i;
    float binsize = IFT_PI / nbins;

    for(i = 1; i < nbins; i++)
        thetas[i] = i*binsize;

    int p;
    for(p = 0; p < image->n; p++){
        int r = label_image->val[p] - 1;
        if(r < 0)
            continue;
        size[r]++;

        //Unsigned binning
        float dx = image_x->val[p];
        float dy = image_y->val[p];

        float radius = sqrt(dx*dx + dy*dy);
        if(dy < 0){
            dy = -dy; dx = -dx;
        }
        float theta = IFT_PI / 2;
        if(!(iftAlmostZero(dx) && iftAlmostZero(dy)))
            theta = atan2(dy,dx);

        //Boundary guarantee
        theta = iftMax(0, iftMin(theta, IFT_PI));

        //Find interval containing theta. Accumulate the radius proportional to the distance to each vector.
        for(i = 1; i < nbins; i++){
            if( thetas[i] > theta){
                //Divide between bin i and bin i-1
                dataset->sample[r].feat[i] += (1 - (thetas[i] - theta)/binsize)*radius;
                dataset->sample[r].feat[i - 1] += ((thetas[i] - theta)/binsize)*radius;

                break;
            }
        }
        if(i == nbins){
            //Divide between bin nbins -1 and bin 0
            dataset->sample[r].feat[0] += (1 - (IFT_PI - theta) / binsize) * radius;
            dataset->sample[r].feat[nbins - 1] += ((IFT_PI - theta) / binsize) * radius;
        }
    }

    int r;
    for(r = 0; r < dataset->nsamples; r++){
        if(size[r] == 0)
            iftError("Empty region", "iftSupervoxelsToSimplifiedHOG");

        float binsum = 0;
        for(i = 0; i < dataset->nfeats; i++)
            binsum += dataset->sample[r].feat[i];
        for(i = 0; i < dataset->nfeats; i++)
            dataset->sample[r].feat[i] = dataset->sample[r].feat[i]/binsum;
    }

    iftDestroyKernel(&sobel_x);
    iftDestroyKernel(&sobel_y);

    iftDestroyImage(&image_x);
    iftDestroyImage(&image_y);

    iftFree(thetas);
    iftFree(size);


    iftSetDistanceFunction(dataset, 1);
    iftCopyRefData(dataset, label_image, IFT_REF_DATA_IMAGE);
    // dataset->ref_data      = (void*)label_image;
    // dataset->ref_data_type = IFT_REF_DATA_IMAGE;

    return dataset;
}

iftDataSet* iftConcatDatasetFeatures(iftDataSet** datasets, int ndatasets){
    if(ndatasets < 1)
        return 0;

    int nsamples = datasets[0]->nsamples;

    int nfeats = 0;
    int i,j,k;
    for(i = 0; i < ndatasets; i++){
        nfeats += datasets[i]->nfeats;

        if(nsamples != datasets[i]->nsamples)
            iftError("Incompatible datasets.", "iftConcatDatasetFeatures");
    }

    iftDataSet *dataset = iftCreateDataSet(nsamples, nfeats);

    int f;
    for(i = 0; i < nsamples; i++){
        f = 0;
        for(j = 0; j < ndatasets; j++){
            for(k = 0; k < datasets[j]->nfeats; k++){
                dataset->sample[i].feat[f] = datasets[j]->sample[i].feat[k];
                f++;
            }
        }
    }

    return dataset;
}

float iftDistMeanSizeSupervoxel(float *f1, float *f2, float *alpha, int n){
    float dist = 0.0;

    int i;
    for(i = 0; i < 3; i++){
        dist += powf(f1[i]/f1[4] - f2[i]/f2[4],2)*alpha[i];
    }

    dist += powf(f1[3] + f2[3], 2)*alpha[3];

    return dist;
}

float* iftMergeMeanSizeSupervoxel(float *f1, float*f2, float*alpha, int n){
    float *output = iftAllocFloatArray(5);

    int i;
    for(i = 0; i < 5; i++){
        output[i] = f1[i] + f2[i];
    }

    return output;
}

float* iftMergeVector(float *f1, float*f2, float*alpha, int n){
    float *output = iftAllocFloatArray(n);

    int i;
    for(i = 0; i < n; i++){
        output[i] = (f1[i] + f2[i])/2;
    }
    return output;
}

//Features: sum band1, sum band2, sum band3, relative area, npixels,
//          xmax, ymax, xmin, ymin, band1 bins
//          (bins_per_band), band2 bins (bins_per_band),
//          band 3 bins (bins_per_band), lbp bins
//          (normalization_value), bins_per_band

iftDataSet* iftSupervoxelsToSelectiveSearchDataset(iftImage *image, iftImage* supervoxels, int bins_per_band, iftColorSpace colorspace){
    int normalization_value = iftNormalizationValue(iftMaximumValue(image));

    if(!iftIsColorImage(image))
        iftError("Color image expected", "iftSupervoxelsToSelectiveSearchDataset");
    if(colorspace != LAB_CSPACE && colorspace != RGB_CSPACE && colorspace != YCbCr_CSPACE)
        iftError("Invalid colorspace.", "iftSupervoxelsToMeanSizeDataSet");
    if(bins_per_band < 1)
        iftError("Invalid number of bins", "iftSupervoxelsToSelectiveSearchDataset");

    iftDataSet *mss_dataset = iftSupervoxelsToMeanSizeDataSet(image, supervoxels, colorspace);
    int nregions = mss_dataset->nsamples;

    int *xmax = iftAllocIntArray(nregions);
    int *ymax = iftAllocIntArray(nregions);
    int *xmin = iftAllocIntArray(nregions);
    int *ymin = iftAllocIntArray(nregions);

    int **colorh = (int**)iftAlloc(nregions, sizeof(int*));
    int **lbph = (int**)iftAlloc(nregions,sizeof(int*));

    int i;
    for(i = 0; i < nregions; i++){
        colorh[i] = iftAllocIntArray(bins_per_band * 3);
        lbph[i] = iftAllocIntArray(normalization_value+1);

        xmin[i] = IFT_INFINITY_INT;
        ymin[i] = IFT_INFINITY_INT;
    }

    float bin_size = ((float)normalization_value+1)/bins_per_band;

    iftAdjRel *adj;

    if (iftIs3DImage(image))
        adj = iftSpheric(sqrtf(3.0));
    else
        adj = iftCircular(sqrtf(2.0));

    int p;
    for(p = 0; p < image->n; p++){
        int region = supervoxels->val[p] - 1;

        iftColor color;
        color.val[0] = image->val[p];
        color.val[1] = image->Cb[p];
        color.val[2] = image->Cr[p];

        int binb1 = 0, binb2 = 0, binb3 = 0;

        if(colorspace == LAB_CSPACE){
            color = iftYCbCrtoRGB(color,normalization_value);
            iftFColor fcolor = iftRGBtoLab(color,normalization_value);
            color = iftLabtoQLab(fcolor,normalization_value);
        }else
        if (colorspace == RGB_CSPACE)
            color = iftYCbCrtoRGB(color,normalization_value);

        binb1 = (int)(color.val[0]/bin_size);
        binb2 = (int)(color.val[1]/bin_size);
        binb3 = (int)(color.val[2]/bin_size);

        if(binb1 < 0 || binb2 < 0 || binb3 < 0 || binb1 >= bins_per_band || binb2 >= bins_per_band || binb3 >= bins_per_band)
            iftError("Invalid color space", "iftSupervoxelsToSelectiveSearchDataset");

        colorh[region][binb1]++;
        colorh[region][binb2 + bins_per_band]++;
        colorh[region][binb3 + 2*bins_per_band]++;

        int lbp_code = 0;
        iftVoxel u = iftGetVoxelCoord(image, p);

        for(i = 1; i < adj->n; i++){
            lbp_code = lbp_code << 1;

            iftVoxel v = iftGetAdjacentVoxel(adj, u, i);

            if(iftValidVoxel(image, v)){
                int q = iftGetVoxelIndex(image, v);

                if(image->val[p] > image->val[q])
                    lbp_code = lbp_code + 1;
            }
        }

        if(lbp_code < 0 || lbp_code > normalization_value)
            iftError("Invalid lbp code", "iftSupervoxelsToSelectiveSearchDataset");

        lbph[region][lbp_code]++;

        xmax[region] = iftMax(xmax[region], u.x);
        ymax[region] = iftMax(ymax[region], u.y);

        xmin[region] = iftMin(xmin[region], u.x);
        ymin[region] = iftMin(ymin[region], u.y);

    }

    iftDataSet *dataset = iftCreateDataSet(nregions, mss_dataset->nfeats + 4 + 3*bins_per_band + normalization_value + 1);

    int j;
    for(i = 0; i < mss_dataset->nsamples; i++){
        for(j = 0; j < mss_dataset->nfeats; j++){
            dataset->sample[i].feat[j] = mss_dataset->sample[i].feat[j];
        }

        dataset->sample[i].feat[j] = xmax[i];
        dataset->sample[i].feat[j+1] = ymax[i];
        dataset->sample[i].feat[j+2] = xmin[i];
        dataset->sample[i].feat[j+3] = ymin[i];

        int colorh_index = j+4;

        for(j = 0; j < bins_per_band; j++){
            dataset->sample[i].feat[j + colorh_index] = colorh[i][j];
            dataset->sample[i].feat[j + colorh_index + bins_per_band] = colorh[i][j + bins_per_band];
            dataset->sample[i].feat[j + colorh_index + 2*bins_per_band] = colorh[i][j + 2*bins_per_band];
        }

        int lbph_index = j + colorh_index + 2*bins_per_band;
        for(j = 0; j <= normalization_value; j++)
            dataset->sample[i].feat[j + lbph_index] = lbph[i][j];

        dataset->sample[i].feat[j + lbph_index] = bins_per_band;

        iftFree(colorh[i]);
        iftFree(lbph[i]);
    }

    dataset->iftArcWeight = iftDistSelectiveSearchSupervoxel;
    dataset->function_number = IFT_NIL;
    iftCopyRefData(dataset, supervoxels, IFT_REF_DATA_IMAGE);
    // dataset->ref_data      = (void*)supervoxels;
    // dataset->ref_data_type = IFT_REF_DATA_IMAGE;

    iftFree(colorh);
    iftFree(lbph);

    iftFree(xmax);
    iftFree(ymax);
    iftFree(xmin);
    iftFree(ymin);

    iftDestroyAdjRel(&adj);
    iftDestroyDataSet(&mss_dataset);

    return dataset;
}

float* iftMergeSelectiveSearchSupervoxel(float *f1, float *f2, float *alpha, int n){
    float *output_mss = iftMergeMeanSizeSupervoxel(f1, f2, alpha, 5);

    float *output = iftAllocFloatArray(n);
    int j;
    for(j = 0; j < 5; j++)
        output[j] = output_mss[j];

    int bins_per_band = f1[n - 1];

    output[5] = iftMax(f1[5], f2[5]);
    output[6] = iftMax(f1[6], f2[6]);
    output[7] = iftMin(f1[7], f2[7]);
    output[8] = iftMin(f1[8], f2[8]);

    int colorh_index = 9;
    for(j = 0; j < bins_per_band; j++){
        output[j + colorh_index] = f1[j + colorh_index] + f2[j + colorh_index];
        output[j + colorh_index + bins_per_band] = f1[j + colorh_index + bins_per_band] + f2[j + colorh_index + bins_per_band];
        output[j + colorh_index + 2*bins_per_band] = f1[j + colorh_index + 2*bins_per_band] + f2[j + colorh_index + 2*bins_per_band];
    }

    int lbph_index = j + colorh_index + 2*bins_per_band;
    int normalization_value = n - lbph_index;
    for(j = 0; j <= normalization_value; j++)
        output[j + lbph_index] = f1[j + lbph_index] + f2[j + lbph_index];

    output[j + lbph_index] = f1[n-1];

    iftFree(output_mss);

    return output;
}

//The 7 alphas weight: mean band1, mean band2, mean band3, relative area, color histogram intersection, lbp histogram, rectangularity
float iftDistSelectiveSearchSupervoxel(float *f1, float *f2, float *alpha, int n){
    float mss_dist = 0;

    int i;
    for(i = 0; i < 3; i++){
        mss_dist += powf(f1[i]/f1[4] - f2[i]/f2[4],2)*alpha[i];
    }

    mss_dist += powf(f1[3] + f2[3], 2)*alpha[3];

    int bins_per_band = f1[n-1];

    int size1 = f1[4];
    int size2 = f2[4];

    int colorh_index = 9;

    float inv_color_intersection = 0;
    int j;
    for(j = 0; j < bins_per_band; j++){
        inv_color_intersection += iftMin(f1[j + colorh_index] / size1, f2[j + colorh_index] / size2);
        inv_color_intersection += iftMin(f1[j + colorh_index + bins_per_band] / size1, f2[j + colorh_index + bins_per_band] / size2);
        inv_color_intersection += iftMin(f1[j + colorh_index + 2 * bins_per_band] / size1, f2[j + colorh_index + 2 * bins_per_band] / size2);
    }
    inv_color_intersection = alpha[4]*(1 - inv_color_intersection/3);

    float inv_lbp_intersection = 0;
    int lbph_index = j + colorh_index + 2*bins_per_band;
    int normalization_value = n - lbph_index;
    for(j = 0; j <= normalization_value; j++){
        inv_lbp_intersection += iftMin(f1[j + lbph_index] / size1, f2[j + lbph_index] / size2);
    }
    inv_lbp_intersection = alpha[5]*(1 - inv_lbp_intersection);

    int new_xmax = iftMax(f1[5], f2[5]);
    int new_ymax = iftMax(f1[6], f2[6]);
    int new_xmin = iftMin(f1[7], f2[7]);
    int new_ymin = iftMin(f1[8], f2[8]);

    float n_pixels_image = (f1[4] + f2[4])/(f1[3] + f2[3]);
    int size_bb = (new_xmax - new_xmin + 1)*(new_ymax - new_ymin + 1);
    float inv_rectangularity = alpha[6]*((size_bb/n_pixels_image) - f1[3] - f2[3]);

    return sqrtf(mss_dist) + inv_color_intersection + inv_lbp_intersection + inv_rectangularity;
}

iftDataSet* iftMImageToEdgesDataSet(iftMImage* mimg,iftAdjRel* A)
{
    fprintf(stderr,"edges: %ld, pixelsize: %ld\n",(A->n-1)*mimg->n,mimg->m);

    iftDataSet* Zedges = iftCreateDataSet((A->n-1)*mimg->n, mimg->m);

    // generating vertical edges
    int s,p,q,b;
    for(s = 0,p = 0; p < mimg->n; p++)
    {
        iftVoxel pixel = iftMGetVoxelCoord(mimg,p);

        int i;
        for(i = 1; i < A->n; i++,s++)
        {
            iftVoxel neighbor = iftGetAdjacentVoxel(A,pixel,i);

            if (iftMValidVoxel(mimg,neighbor))
                q = iftMGetVoxelIndex(mimg,neighbor);
            else
                q = p;

            for(b = 0; b < mimg->m ; b++)
                Zedges->sample[s].feat[b] = sqrt(mimg->val[p][b] * mimg->val[q][b]);
            Zedges->sample[s].id = s;
        }
    }

    return Zedges;
}

iftDataSet* iftMImageToLabeledEdgesDataSet(iftMImage* mimg,iftImage* imgGT, iftImage* imglabels,iftAdjRel* A)
{ // 2d iftMImage
    fprintf(stderr,"edges: %ld, pixelsize: %ld\n",(A->n-1)*mimg->n,mimg->m);

    iftVerifyImageDomains(imgGT,imglabels,"iftMImageToLabeledEdgesDataSet");

    iftDataSet* Zedges = iftCreateDataSet((A->n-1)*mimg->n, mimg->m);
    Zedges->nclasses = 2;

    // generating vertical edges
    int s,p,q,b;
    for(s = 0, p = 0; p < mimg->n; p++)
    {
        iftVoxel pixel   = iftMGetVoxelCoord(mimg ,p);

        int i;
        for(i = 1; i < A->n; i++)
        {
            iftVoxel neighbor   = iftGetAdjacentVoxel(A,pixel  ,i);

            if (iftMValidVoxel(mimg,neighbor))
                q   = iftMGetVoxelIndex(mimg ,neighbor);

            else
                q   = p;

            if ( (imgGT->val[p] != imgGT->val[q]) && (imglabels->val[p] != imglabels->val[q]) )
            {
                Zedges->sample[s].id = s;
                Zedges->sample[s].truelabel = 2;  // edges
                for(b = 0; b < mimg->m ; b++)
                    Zedges->sample[s].feat[b] = sqrt(mimg->val[p][b] * mimg->val[q][b]);
                s++;
            }
            else
            {
                if ( (imglabels->val[p] == imglabels->val[q]) && (imgGT->val[p] == imgGT->val[q]) )
                {
                    Zedges->sample[s].id = s;
                    Zedges->sample[s].truelabel = 1; // non edges
                    for(b = 0; b < mimg->m ; b++)
                        Zedges->sample[s].feat[b] = sqrt(mimg->val[p][b] * mimg->val[q][b]);
                    s++;
                }
            }
        }
    }
    Zedges->nsamples = s; // gambiarra

    return Zedges;
}

iftFImage* iftEdgesDataSetToFImage(iftDataSet* dataset,iftMImage* mimg,iftAdjRel* A)
{
    fprintf(stderr,"(%ld x %d):\n",mimg->n,dataset->nsamples);
    if ( (A->n-1)*mimg->n != dataset->nsamples)
        iftError("incompatible dimensions ", "iftEdgesDataSetToImage");

    iftFImage* output = iftCreateFImage(mimg->xsize,mimg->ysize,mimg->zsize);

    int s,p,q;
    for(s = 0;s < dataset->nsamples; s++)
    {
        if (dataset->sample[s].label == 2) // edges
        {
            p = s / (A->n-1);
            iftVoxel pixel = iftMGetVoxelCoord(mimg,p);
            q = s - p * (A->n-1);

            iftVoxel neighbor = iftGetAdjacentVoxel(A,pixel,q+1);
            if (iftMValidVoxel(mimg,neighbor))
            {
                q = iftMGetVoxelIndex(mimg,neighbor);
                if (dataset->sample[s].weight > 0)
                {
                    output->val[p] += dataset->sample[s].weight;
                    output->val[q] += dataset->sample[s].weight;
                }
            }
        }
    }

//  float beta = 1.0; // to be adjusted/learned
//  for(p = 0 ; p < output->n; p++)
//      output->val[p] = powf(output->val[p], beta);

    return output;
}

iftImage* iftPixelDataSetToImage(iftDataSet* dataset,iftImage* img)
{
    //  fprintf(stderr,"(%d x %d):\n",img->n,dataset->nsamples);
    if ( img->n != dataset->nsamples) {
        iftError("incompatible dimensions ", "iftEdgesDataSetToImage");
    }

    iftImage* output = iftCreateImage(img->xsize,img->ysize,img->zsize);

    for(int s = 0;s < dataset->nsamples; s++) {
        if ( (dataset->sample[s].label == 2) && (dataset->sample[s].truelabel == 2) )
            output->val[s] = 255;
        else if ( (dataset->sample[s].label == 2) && (dataset->sample[s].truelabel == 1) )
            output->val[s] = 170;
        else if ( (dataset->sample[s].label == 1) && (dataset->sample[s].truelabel == 1) )
            output->val[s] =   0;
        else if ( (dataset->sample[s].label == 1) && (dataset->sample[s].truelabel == 2) )
            output->val[s] =  85;
    }

    return output;
}

iftFImage* iftPixelDataSetToFImage(iftDataSet* dataset,iftImage* img)
{
    //  fprintf(stderr,"(%d x %d):\n",img->n,dataset->nsamples);
    if ( img->n != dataset->nsamples) {
        iftError("incompatible dimensions ", "iftEdgesDataSetToImage");
    }

    iftFImage* output = iftCreateFImage(img->xsize,img->ysize,img->zsize);

    for(int s = 0;s < dataset->nsamples; s++) {
        if (dataset->sample[s].label == 2) // Border pixels
            output->val[s] =  dataset->sample[s].weight;
        else
            output->val[s] = -dataset->sample[s].weight;
    }

    return output;
}



iftDataSet *iftMMKernelToDataSet(iftMMKernel *K)
{
    iftDataSet *Z=iftCreateDataSet(K->nkernels, K->nbands*K->A->n);

    for (int k=0, s=0; k < K->nkernels; k++, s++) {
        for (int b=0; b < K->nbands; b++) {
            for (int i=0; i < K->A->n; i++) {
                Z->sample[s].feat[i+b*K->A->n] = K->weight[k][b].val[i];
            }
        }
        Z->sample[s].id = k;
    }
    iftSetStatus(Z,IFT_TRAIN);
    iftCopyRefData(Z, K, IFT_REF_DATA_MMKERNEL);
    // Z->ref_data      = K;
    // Z->ref_data_type = IFT_REF_DATA_MMKERNEL;

    return(Z);
}

iftMMKernel *iftDataSetToMMKernel(iftDataSet *Z)
{
    iftMMKernel *K1, *K2;

    if (Z->ref_data == NULL)
        iftError("It requires a MMKernel dataset", "iftDataSetToMMKernel");

    K1 = (iftMMKernel *) Z->ref_data;
    K2 = iftCreateMMKernel(K1->A,K1->nbands,K1->nkernels);

    for (int k=0, s=0; k < K1->nkernels; k++, s++) {
        for (int b=0; b < K1->nbands; b++) {
            for (int i=0; i < K1->A->n; i++) {
                K2->weight[k][b].val[i] = Z->sample[s].feat[i+b*K1->A->n];
            }
        }
    }

    return(K2);
}


void iftExtractSamplesOfMImages(iftMImage *img, int truelabel, iftAdjRel *A, int sample_stride, FILE **fp) {
    iftFastAdjRel *F = iftCreateFastAdjRel(A, img->tby, img->tbz);

    int x, y, z;
    // Used to calculate the index of the FILEs (number of the dataset/region sampled)
    int tby = (img->xsize - 2*F->bx)/sample_stride; // amount of splits in axis x
    int tbz = ((img->ysize - 2*F->by)/sample_stride) * tby; // amount of splits in axis (x, y)

#pragma omp parallel for private(x, y, z) shared(img, truelabel, A, sample_stride, fp, tby, tbz, F)
    for (y = F->by; y < img->ysize - F->by; y += sample_stride)
        for (x = F->bx; x < img->xsize - F->bx; x += sample_stride)
            for (z = F->bz; z < img->zsize - F->bz; z += sample_stride) {
                int p = x + img->tby[y] + img->tbz[z];
                int index = (x - F->bx)/sample_stride + ((y - F->by)*tby)/sample_stride +
                            ((z - F->bz)*tbz)/sample_stride;

                iftFeatures *features = iftCreateFeatures(A->n * img->m);
//              printf("****** p = %d\n", p);
//              printf("------ index = %d\n\n", index);

                int b;
                for (b = 0; b < img->m; b++) {
                    int i;
                    for (i = 0; i < A->n; i++) {
                        int q = p + F->dq[i];
                        features->val[i + b*A->n] = img->val[q][b];
                    }
                }

#pragma omp critical
                iftWriteFeaturesInFile(features, truelabel, fp[index]);

                iftDestroyFeatures(&features);
            }
}

void iftSetDistanceFunction(iftDataSet *Z, int function_number)
{
    Z->function_number = function_number;
    switch(function_number) {
        case 1:
            Z->iftArcWeight = iftDistance1;
            break;
        case 2:
            Z->iftArcWeight = iftDistance2;
            break;
        case 3:
            Z->iftArcWeight = iftDistance3;
            break;
        case 4:
            Z->iftArcWeight = iftDistance4;
            break;
        case 5:
            Z->iftArcWeight = iftDistance5;
            break;
        case 6:
            Z->iftArcWeight = iftDistance6;
            break;
        case 7:
            Z->iftArcWeight = iftDistance7;
            break;
        case 8:
            Z->iftArcWeight = iftDistance8;
            break;
        case 9:
            Z->iftArcWeight = iftDistance9;
            break;
        case 10:
            Z->iftArcWeight = iftDistance10;
            break;
        case 11:
            Z->iftArcWeight = iftDistance11;
            break;
        case 12:
            Z->iftArcWeight = iftDistance12;
        break;
        default:
            Z->iftArcWeight = iftDistance1;
            Z->function_number = 1; // Since the function_number is not one of the previous cases, we reset it to 1
            iftWarning("Function number %02d is invalid. Distance function %02d will be used instead.", "iftSetDistanceFunction",
                       function_number, Z->function_number);

    }
}

/* Dataset for image segmentation based on region merging: the dataset
 is created with multi-band image attributes for each region. */

iftDataSet* iftRegionMergingDataSet(iftMImage *image, iftImage* label_image){

    iftDataSet *dataset    = NULL;
    int         nregions   = iftMaximumValue(label_image);
    int        *size       = iftAllocIntArray(nregions);
    int         b, p, r;

    dataset       = iftCreateDataSet(nregions, image->m+1);
    for(p = 0; p < image->n; p++){
        r = label_image->val[p] - 1;
        size[r]++; // compute region size
        // compute sum of attribute values in each region
        for (b=0; b < image->m; b++)
            dataset->sample[r].feat[b] += image->val[p][b];
    }
    for(r = 0; r < nregions; r++){
        if(size[r] == 0){
            char msg[200];
            sprintf(msg,"Missing region %d from %d to %d labels\n",r+1,1,nregions);
            iftError(msg, "iftRegionGraphDataSet");
        }
        // region size
        dataset->sample[r].feat[image->m] = (float)size[r];
        // sample id is the region label
        dataset->sample[r].id = r + 1;
    }

    dataset->iftArcWeight    = iftRegionMergingDist;
    dataset->function_number = IFT_NIL;
    iftCopyRefData(dataset, label_image, IFT_REF_DATA_IMAGE);
    // dataset->ref_data        = (void *)label_image;
    iftFree(size);

    return(dataset);
}

float iftRegionMergingDist(float *f1, float *f2, float *alpha, int n)
{
    float dist=0.0;

    for (int i=0; i < n-1; i++) {
        dist += pow((f1[i]/f1[n-1]-f2[i]/f2[n-1]),2.0)*alpha[i];
    }

    return(dist);
}

float *iftRegionMergingFeat(float *f1, float*f2, float*alpha, int n){
    float *output = iftAllocFloatArray(n);

    int i;
    for(i = 0; i < n; i++){
        output[i] = (f1[i] + f2[i]);
    }
    return (output);
}


void iftSwitchSamples(iftSample *s1, iftSample *s2, int nfeats) {
    iftSample aux;

    aux.feat   = iftAllocFloatArray(nfeats);
    aux.id     = (*s1).id;
    aux.truelabel  = (*s1).truelabel;
    aux.label  = (*s1).label;
    iftSetSampleStatus(&aux, (*s1).status);
    aux.weight = (*s1).weight;
    for (int i = 0; i < nfeats; i++)
        aux.feat[i] = (*s1).feat[i];

    (*s1).id     = (*s2).id;
    (*s1).truelabel  = (*s2).truelabel;
    (*s1).label  = (*s2).label;
    iftSetSampleStatus(s1, (*s2).status);
    (*s1).weight = (*s2).weight;
    for (int i = 0; i < nfeats; i++)
        (*s1).feat[i] = (*s2).feat[i];

    (*s2).id     = aux.id;
    (*s2).truelabel  = aux.truelabel;
    (*s2).label  = aux.label;
    iftSetSampleStatus(s2, aux.status);
    (*s2).weight = aux.weight;
    for (int i = 0; i < nfeats; i++)
        (*s2).feat[i] = aux.feat[i];
}

// The feature vector of sout must be allocated previously
void iftCopySample(iftSample *src, iftSample *dst, int n_feats, bool copy_feats) {
    (*dst).id     = (*src).id;
    (*dst).truelabel  = (*src).truelabel;
    (*dst).label  = (*src).label;
    (*dst).group  = (*src).group;
    iftSetSampleStatus(dst, (*src).status);
    (*dst).weight = (*src).weight;
    (*dst).isSupervised = (*src).isSupervised;
    (*dst).isLabelPropagated = (*src).isLabelPropagated;
    if (copy_feats) {
        for (int i = 0; i < n_feats; i++){
            (*dst).feat[i] = (*src).feat[i];
        }
    }
    else (*dst).feat = (*src).feat;
}


/*
* Switch the Not SV Samples from Ztrain per Errors Samples from Ztest.
* not_SVs = array that stores the indices of the NOT SVs Samples from Ztrain
* error_samples = array that stores the indices of the Error Samples from Ztest
*/
void iftSwitchNotSVsPerErrorSamples(iftDataSet *Ztrain, iftDataSet *Ztest, int *not_SVs,
                                    int num_not_SVs, int *error_samples, int num_errors) {

    if (Ztrain->nfeats != Ztest->nfeats)
        iftError("Ztrain->nfeats != Ztest->nfeats", "iftSwitchNotSVsPerErrors");

    int *count  = iftAllocIntArray(num_not_SVs);
    for (int s = 0; s < num_not_SVs; s++) {
        count[s]  = 100;
    }

    // Randomly select samples
    int t = 0, high = num_not_SVs-1, i;
    int idx_not_SV, idx_error;
    while ((t < num_errors) && (t < num_not_SVs)) {
        i = iftRandomInteger(0, high);

        if (count[i] == 0) {
            idx_not_SV = not_SVs[i];
            idx_error  = error_samples[t];
//          printf("# not_SVs[%d]: %d <-----> error_sample[%d]: %d\n", i, not_SVs[i], t, error_samples[t]);
            iftSwap(not_SVs[i], not_SVs[high]);
            iftSwap(count[i],  count[high]);
            iftSwitchSamples(&Ztrain->sample[idx_not_SV], &Ztest->sample[idx_error], Ztrain->nfeats);
            t++;
            high--;
        } else {
            count[i]--;
        }
    }

    iftFree(count);
}


iftDataSet *iftSelectNegativeSamples(iftDataSet *Z, int positive_class) {
    int nsamples = 0;
    iftDataSet *Zout = NULL;

    // Count the number of samples of the selected class
    for (int s = 0; s < Z->nsamples; s++)
        if (Z->sample[s].truelabel != positive_class)
            nsamples++;

    if (nsamples > 0) {
        Zout = iftCreateDataSet(nsamples, Z->nfeats);
        Zout->nclasses = Z->nclasses - 1;
        int idx = 0;

        for (int s = 0; s < Z->nsamples; s++)
            if (Z->sample[s].truelabel != positive_class) {
                iftCopySample(&Z->sample[s], &Zout->sample[idx], Z->nfeats, true);
                idx++;
            }
        printf("\n********************************\n");
        printf("- Zout: (nclasses, nsamples, nfeats): (%d, %d, %d)\n", Zout->nclasses, Zout->nsamples, Zout->nfeats);
        printf("********************************\n");
    }

    return Zout;
}


iftDataSet *iftExtractSamplesFromClass(const iftDataSet *Z, int labelClass) {

    int nsamples = 0;

    for (int i = 0; i < Z->nsamples; ++i) {
        int truelabel = Z->sample[i].truelabel;
        if(truelabel == labelClass) nsamples++;
    }

    iftDataSet* Zc = iftCreateDataSet(nsamples, Z->nfeats);

    Zc->nclasses = 1;
    for (int i = Z->nsamples-1; i >= 0; --i) {
        if(Z->sample[i].truelabel == labelClass)
            iftCopySample(&(Z->sample[i]), &(Zc->sample[--nsamples]), Z->nfeats, true);
    }

    return Zc;
}

iftDataSet *iftExtractSamplesFromClass_omp(const iftDataSet *Z, int labelClass, int opt) {

    int nsamples = 0;

    #pragma omp parallel for if(opt == 6) reduction(+:nsamples) shared(Z, labelClass) // OK
    for (int i = 0; i < Z->nsamples; ++i) {
        if(Z->sample[i].truelabel == labelClass) nsamples++;
    }

    iftDataSet* Zc = iftCreateDataSet(nsamples, Z->nfeats);
    Zc->nclasses = 1;

    //#pragma omp parallel for shared(Z, Zc, labelClass) // NOP
    for (int i = Z->nsamples-1; i >= 0; --i) {
        if(Z->sample[i].truelabel == labelClass)
            iftCopySample(&(Z->sample[i]), &(Zc->sample[--nsamples]), Z->nfeats, true);
    }

    return Zc;
}

int iftCountNumberOfClassesDataSet(iftDataSet *Z) {

    int *labels = iftAllocIntArray(Z->nsamples);

    for (int i = 0; i < Z->nsamples; ++i) {
        labels[i] = Z->sample[i].truelabel;
    }

    iftIntArray *unique = iftIntArrayUnique(labels, Z->nsamples);
    Z->nclasses = iftMaxValIntArray(unique->val, unique->n, NULL);
    iftFree(labels);
    
    return Z->nclasses;
}

int iftCountNumberOfGroupsDataSet(iftDataSet *Z) {

    int *groups = iftAllocIntArray(Z->nsamples);

    for (int i = 0; i < Z->nsamples; ++i) {
        groups[i] = Z->sample[i].group;
    }

    iftIntArray *unique = iftIntArrayUnique(groups, Z->nsamples);
    Z->ngroups = iftMaxValIntArray(unique->val, unique->n, NULL);
    iftFree(groups);
    
    return Z->ngroups;
}

void iftUpdateNumbTrainSamples(iftDataSet *Z)
{
    Z->ntrainsamples = iftCountSamples(Z, IFT_TRAIN);
}


int* iftCountSamplesPerClassDataSet(iftDataSet *Z)
{
    int *sampPerClass = iftAllocIntArray(Z->nclasses+1);

    for(int p=0; p<Z->nsamples; p++) {
      if (Z->sample[p].truelabel>0)
	sampPerClass[Z->sample[p].truelabel]++;
    }

    return sampPerClass;
}

int* iftCountSamplesPerGroupDataSet(iftDataSet *Z)
{
    int *sampPerGroup = iftAllocIntArray(Z->ngroups+1);

    for(int p=0; p<Z->nsamples; p++) {
      if (Z->sample[p].group>0)
        sampPerGroup[Z->sample[p].group]++;
    }

    return sampPerGroup;
}

void iftGetMinMaxWeightInDataSet(const iftDataSet *Z, float *min_weight, float *max_weight) {
    if (Z == NULL)
        iftError("Dataset is NULL", "iftGetMinMaxWeightInDataSet");
    if (Z->nsamples == 0)
        iftError("Dataset is Empty (no samples)", "iftGetMinMaxWeightInDataSet");
    if (min_weight == NULL)
        iftError("Reference for Min Weight is NULL", "iftGetMinMaxWeightInDataSet");
    if (max_weight == NULL)
        iftError("Reference for Max Weight is NULL", "iftGetMinMaxWeightInDataSet");

    *min_weight = IFT_INFINITY_FLT;
    *max_weight = IFT_INFINITY_FLT_NEG;

    for (int s = 0; s < Z->nsamples; s++) {
        if (Z->sample[s].weight < *min_weight)
            *min_weight = Z->sample[s].weight;
        if (Z->sample[s].weight > *max_weight)
            *max_weight = Z->sample[s].weight;
    }
}

void iftCopyDataSetToDataSet(iftDataSet *Zsrc, iftDataSet *Zdst, int begin) {
    if (Zdst == NULL)
        iftError("The destination dataset is NULL", "iftCopyDataSetToDataSet");

    if (Zdst->data == NULL) {
        Zdst->data = iftCreateMatrix(Zsrc->nfeats, Zsrc->capacity);
    }

    if (begin + Zsrc->nsamples > Zdst->nsamples)
        iftError("The destination dataset does not have the capacity to receive the source dataset", "iftCopyDataSetToDataSet");

    if(Zsrc->nfeats != Zdst->nfeats)
        iftError("The source and destination datasets must have the same number of features", "iftCopyDataSetToDataSet");

    if(Zsrc->nclasses != Zdst->nclasses)
        iftError("The source and destination datasets must have the same number of classes", "iftCopyDataSetToDataSet");

    if(Zsrc->ngroups != Zdst->ngroups)
        iftError("The source and destination datasets must have the same number of groups", "iftCopyDataSetToDataSet");


    if(Zsrc->ref_data_type == IFT_REF_DATA_FILESET && Zdst->ref_data == NULL)
        Zdst->ref_data = iftCreateFileSet(Zdst->nsamples);

    Zdst->ntrainsamples += Zsrc->ntrainsamples;
    for (int i = 0; i < Zsrc->nsamples; i++) {
        iftCopySample(&Zsrc->sample[i], &Zdst->sample[i+begin], Zsrc->nfeats, true);

        /* copy files to the new fileset */
        if(Zsrc->ref_data_type == IFT_REF_DATA_FILESET) {
            int id = Zsrc->sample[i].id;
            iftFileSet *fsAux1 = (iftFileSet*)Zsrc->ref_data;
            iftFileSet *fsAux2 = (iftFileSet*)Zdst->ref_data;
            fsAux2->files[i+begin] = iftCreateFile(fsAux1->files[id]->path);
            iftUpdateDataSetFileInfo(fsAux2->files[i+begin]);
        }
    }

    for (int s = 0; s < Zdst->nsamples; s++)
        Zdst->sample[s].id = s;
}

iftDataSet *iftMergeDataSetArray(iftDataSet **Z_array, int n_datasets, bool copy_feats) {
    if (Z_array == NULL)
        iftError("Array of Datasets is NULL", "iftMergeDataSetArray");
    if (n_datasets <= 0)
        iftError("Number of Datasets %d is <= 0", "iftMergeDataSetArray", n_datasets);

    // some dataset(s) could have NULL
    int i0 = iftGetFirstDataSetIndexFromArray(Z_array, n_datasets);

    // all datasets are NULL
    if (i0 == -1)
        iftError("All datasets from array are NULL", "iftMergeDataSetArray");

    int n_feats         = Z_array[i0]->nfeats;
    int n_total_samples = iftGetTotalNumberOfSamplesFromDataSets(Z_array, n_datasets);

    iftDataSet *Z = NULL;
    if (copy_feats)
        Z = iftCreateDataSet(n_total_samples, n_feats);
    else
        Z = iftCreateDataSetWithoutFeatsVector(n_total_samples, n_feats);

    iftFileSet *newFileSet = NULL;
    if(Z_array[i0]->ref_data_type == IFT_REF_DATA_FILESET)
        newFileSet = iftCreateFileSet(n_total_samples);

    int s = 0;
    int n_classes = 0;
    int n_groups = 0;

    for (int i = i0; i < n_datasets; i++) {
        if (Z_array[i] != NULL) {
            n_classes = iftMax(n_classes, Z_array[i]->nclasses);
            n_groups = iftMax(n_groups, Z_array[i]->ngroups);
            Z->ntrainsamples += Z_array[i]->ntrainsamples;

            for (int t = 0; t < Z_array[i]->nsamples; t++) {
                // copies samples and their feature vectors in merged dataset
                iftCopySample(&Z_array[i]->sample[t], &Z->sample[s], n_feats, copy_feats);
                Z->sample[s].id = s; // ID re-numbering

                /* copy files to the new fileset */
                if(Z_array[i0]->ref_data_type == IFT_REF_DATA_FILESET) {
                    int id = Z_array[i]->sample[t].id;
                    iftFileSet *fileSetAux = (iftFileSet*)Z_array[i]->ref_data;
                    newFileSet->files[s] = iftCreateFile(fileSetAux->files[id]->path);
                    iftUpdateDataSetFileInfo(newFileSet->files[s]);
                }
                s++;
            }
        }
    }
    Z->nclasses = n_classes;
    Z->ngroups  = n_groups;
    //iftIsThereTrueLabelOrGroupZero(Z, "iftMergeDataSetArray");

    if(Z_array[i0]->ref_data_type == IFT_REF_DATA_FILESET) {
        iftCopyRefData(Z, newFileSet, IFT_REF_DATA_FILESET);
    }

    return Z;
}


iftDataSet *iftMergeDataSetArrayBySampling(iftDataSet **Z_array, int n_datasets, iftSampler **samplers,
                                      bool copy_feats) {
    if (Z_array == NULL)
        iftError("Array of Datasets is NULL", "iftMergeDataSetArrayBySampling");
    if (n_datasets <= 0)
        iftError("Number of Datasets %d is <= 0", "iftMergeDataSetArrayBySampling", n_datasets);
    if (samplers == NULL)
        iftError("Array of Samplers is NULL", "iftMergeDataSetArrayBySampling");

    // counts the total number of selected samples
    int n_total_samples = 0;
    for (int i = 0; i < n_datasets; i++) {
        if (Z_array[i]->nsamples != samplers[i]->nsamples)
            iftError("DataSet %d has different number of samples from its sampler: %d != %d",
                     "iftMergeDataSetArrayBySampling", Z_array[i]->nsamples, samplers[i]->nsamples);

#pragma omp parallel for reduction(+:n_total_samples)
        for (int t = 0; t < samplers[i]->nsamples; t++)
            if (samplers[i]->status[0][t] == IFT_TRAIN) {
                n_total_samples++;
            }
    }

    int n_feats   = Z_array[0]->nfeats; // it is expected that all datasets have the same number of feats
    iftDataSet *Z = iftCreateDataSet(n_total_samples, n_feats);

    // TODO: we could parallelize this code if we save the start sample index in the joined dataset for
    // each dataset
    // copies the samples
    int *true_labels = iftAllocIntArray(n_total_samples);
    int *labels      = iftAllocIntArray(n_total_samples);
    int s = 0;
    for (int i = 0; i < n_datasets; i++) {
        for (int t = 0; t < samplers[i]->nsamples; t++) {
            // just the first sampler's iteration is considered
            if (samplers[i]->status[0][t] == IFT_TRAIN) {

                iftCopySample(&Z_array[i]->sample[t], &Z->sample[s], n_feats, copy_feats);
                true_labels[s] = Z->sample[s].truelabel;
                labels[s]      = Z->sample[s].label;
                s++;
            }
        }
    }
    Z->nclasses = iftCountUniqueIntElems(true_labels, n_total_samples);
    Z->ngroups  = iftCountUniqueIntElems(labels, n_total_samples);
    iftIsThereTrueLabelOrGroupZero(Z, "iftMergeDataSetArrayBySampling");


    // DESTROYERS
    iftFree(true_labels);
    iftFree(labels);

    return Z;
}

iftIntArray * iftGetDataSetTrueLabels(const iftDataSet *dataset) {
    iftIntArray* labels = iftCreateIntArray(dataset->nsamples);

#pragma omp parallel for
    for (int i = 0; i < labels->n; ++i) {
        labels->val[i] = dataset->sample[i].truelabel;
    }

    return labels;
}

void iftSampleFileSet(iftFileSet *Z, const iftSampler *sampler, int iteration) {
    if(iteration<0 || iteration>=sampler->niters) {
        iftError("Invalid iteration.", "iftSampleFileSet");
    }

    if(sampler->nsamples!=Z->n) {
        iftError("Invalid number of samples.", "iftSampleFileSet");
    }

    for (int i = 0; i < Z->n; ++i) {
        Z->files[i]->status = sampler->status[iteration][i];
    }

}


iftFileSet *iftExtractFileSamples(const iftFileSet *fileSet, iftSampleStatus status) {
    int nfiles = 0;

    for (int i = 0; i < fileSet->n; ++i) {
        if (iftHasStatus(fileSet->files[i]->status, status)) { nfiles++; }
    }

    iftFileSet* extracted = (iftFileSet*) iftAlloc(1, sizeof(iftFileSet));
    extracted->files = (iftFile**) iftAlloc(nfiles, sizeof(iftFile*));
    extracted->n = nfiles;

    int extractedIdx = 0;
    for (int i = 0; i < fileSet->n; ++i) {
        if (iftHasStatus(fileSet->files[i]->status, status)) {
            extracted->files[extractedIdx++] = iftCopyFile(fileSet->files[i]);
        }
    }

    return extracted;
}

void iftSplitFileSetInTwo(iftFileSet *fileset, int nSamplesTrain, int sampMet, iftFileSet **trainFileset, iftFileSet **testFileset)
{
    iftSampler *sampler = NULL;
    int *labels = NULL;
    
    switch (sampMet) {
        case 0:
            sampler = iftRandomSubsampling(fileset->n, 1, nSamplesTrain);
            break;

        case 1:
            labels = iftFileSetLabels(fileset);
            sampler = iftStratifiedRandomSubsampling(labels, fileset->n, 1, nSamplesTrain);
            iftFree(labels);
            break;

        default:
            iftError("Invalid sampling method chosen", "iftSplitFileSetInTwo");
    }

    iftSampleFileSet(fileset, sampler, 0);

    if(*trainFileset != NULL) iftDestroyFileSet(trainFileset);
    *trainFileset = iftExtractFileSamples(fileset, IFT_TRAIN);

    if(*testFileset != NULL) iftDestroyFileSet(testFileset);
    *testFileset = iftExtractFileSamples(fileset, IFT_TEST);

    iftDestroySampler(&sampler);
}

int iftCountSamples(const iftDataSet *Z, iftSampleStatus status) {
    int s, nsamples;

    nsamples=0;
    for (s=0; s < Z->nsamples; s++)
        if (iftHasSampleStatus(Z->sample[s], status)){
            nsamples++;
        }

    return nsamples;
}

int iftCountSamplesTrueLabel(const iftDataSet* Z, int truelabel, iftSampleStatus status) {
    int nsamples = 0;
    for (int i = 0; i < Z->nsamples; ++i) {
        iftSampleStatus s = Z->sample[i].status;
        if((s & status) && Z->sample[i].truelabel==truelabel)
            nsamples++;
    }
    return nsamples;
}

int iftCountSamplesLabel(const iftDataSet* Z, int label, iftSampleStatus status) {
    int nsamples = 0;
    for (int i = 0; i < Z->nsamples; ++i) {
        if((iftHasSampleStatus(Z->sample[i], status)) && Z->sample[i].label==label)
            nsamples++;
    }
    return nsamples;
}

iftDistanceTable* iftCompEuclDistanceTable(const iftDataSet *Z1, const iftDataSet* Z2) {

  iftMatrix *X, *Y;
  iftMatrix *XX, *YY, *XY, *YT;

  X   = iftDataSetToFeatureMatrix(Z1);

  if(Z1==Z2)/* Trick to improve the computation of dataset to itself */
    Y = X;
  else
    Y = iftDataSetToFeatureMatrix(Z2);

  YT = iftTransposeMatrix(Y);
  XY = iftMultMatrices(X, YT);
  iftDestroyMatrix(&YT);

  iftPointWiseMatrixPowInPlace(X, 2);
  XX = iftMatrixSumRow(X);
  iftDestroyMatrix(&X);

  if(Z1==Z2) { /* Trick to improve the computation of dataset to itself */
    YY = XX;
  }else {
    iftPointWiseMatrixPowInPlace(Y, 2);
    YY = iftMatrixSumRow(Y);
    iftDestroyMatrix(&Y);
  }

  iftDistanceTable* dist = iftCreateDistanceTable(Z1->nsamples, Z2->nsamples);

  if(Z1==Z2) { /* Trick to improve the computation of dataset to itself */
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < Z1->nsamples-1; ++i) {
      for (int j = i + 1; j < Z2->nsamples; ++j) {
    double distance = XX->val[i] + YY->val[j] - 2 * iftMatrixElem(XY, j, i);
    dist->distance_table[i][j] = dist->distance_table[j][i] = (float) sqrt(distance);
      }
    }
  }
  else {
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < Z1->nsamples; ++i) {
      for (int j = 0; j < Z2->nsamples; ++j) {
    double distance = XX->val[i] + YY->val[j] - 2 * iftMatrixElem(XY, j, i);
      dist->distance_table[i][j] = (float) sqrt(distance);
      }
    }
  }


  iftDestroyMatrix(&XX);
  iftDestroyMatrix(&XY);

  if(Z1!=Z2) {
    iftDestroyMatrix(&YY);
  }

  return dist;
}


iftDataSet *iftMergeDataSets(const iftDataSet *Z1, const iftDataSet *Z2) {
    if ((Z1 == NULL) && (Z2 == NULL))
        iftError("Z1 and Z2 are NULL", "iftMergeDataSets");
    if ((Z1 != NULL) && (Z2 == NULL))
        return(iftCopyDataSet(Z1, true));
    if ((Z1 == NULL) && (Z2 != NULL))
        return(iftCopyDataSet(Z2, true));

    if (Z1->nfeats != Z2->nfeats)
        iftError("Number of features is different between Z1 and Z2", "iftMergeDataSets");

//  if(Z1->ngroups != 0 || Z2->ngroups != 0)
//      iftError("Data sets must be unlabeled (ngroups must be zero)", "iftMergeDataSets");

    if (Z1->iftArcWeight != Z2->iftArcWeight)
        iftError("iftArcWeight is different between Z1 and Z2", "iftMergeDataSets");

    if((iftIsCentralizedDataSet(Z1) && !iftIsCentralizedDataSet(Z2)) ||
        (!iftIsCentralizedDataSet(Z1) && iftIsCentralizedDataSet(Z2)))
        iftError("One of the datasets is standardized and the other is not, merging the datasets is not possible",
            "iftMergeDataSets");

    if(iftIsCentralizedDataSet(Z1) && iftIsCentralizedDataSet(Z1)) {
        if(!iftCompareFloatArrays(Z1->fsp.mean, Z2->fsp.mean, Z1->fsp.nfeats) ||
            !iftCompareFloatArrays(Z1->fsp.stdev, Z2->fsp.stdev, Z1->fsp.nfeats))
            iftError("Both datasets are standardized, but their FSP params are different (mean/stdev), merging the datasets is not possible",
                "iftMergeDataSets");
    }

    if ((Z1->ref_data != Z2->ref_data) && (Z1->ref_data_type != IFT_REF_DATA_FILESET) && (Z1->ref_data_type != IFT_REF_DATA_CSV))
        iftError("Data sets must have the same reference data or both be NULL or both be IFT_REF_DATA_FILESET or IFT_REF_DATA_CSV",
            "iftMergeDataSets");

    int nsamples = Z1->nsamples + Z2->nsamples;
    int nfeats = Z1->nfeats;

    iftDataSet *Z3    = iftCreateDataSet(nsamples, nfeats);
    Z3->nclasses      = iftMax(Z1->nclasses, Z2->nclasses);
    Z3->iftArcWeight  = Z1->iftArcWeight;
    Z3->ntrainsamples = Z1->ntrainsamples + Z2->ntrainsamples;
    Z3->ngroups       = iftMax(Z1->ngroups, Z2->ngroups);
    Z3->function_number = Z1->function_number;
    Z3->fsp             = iftCopyFeatSpaceParam(Z1);

    for(int f = 0; f < Z1->nfeats; f++)
        Z3->alpha[f] = Z1->alpha[f];

    iftIntArray *true_labels = iftCreateIntArray(Z3->nsamples);

    // Copy samples from Z1 to Z3
    #pragma omp parallel for
    for (int s = 0; s < Z1->nsamples; s++) {
        iftCopySample(&Z1->sample[s], &Z3->sample[s], Z1->nfeats, true);
        true_labels->val[s] = Z3->sample[s].truelabel;
    }

    // Copy samples from Z2 to Z3
    #pragma omp parallel for
    for (int s = Z1->nsamples; s < Z3->nsamples; s++) {
        iftCopySample(&Z2->sample[s - Z1->nsamples], &Z3->sample[s], Z1->nfeats, true);
        true_labels->val[s] = Z3->sample[s].truelabel;
    }

    iftIntArray *unique_true_labels = iftIntArrayUnique(true_labels->val, true_labels->n);
    Z3->nclasses = unique_true_labels->n;
    // exclude 0 truelabel
    if (unique_true_labels->val[0] == 0)
        Z3->nclasses--;
    iftDestroyIntArray(&true_labels);
    iftDestroyIntArray(&unique_true_labels);

    Z3->ref_data_type = Z1->ref_data_type;
    // merge the file sets
    if (Z1->ref_data_type == IFT_REF_DATA_FILESET) {
        Z3->ref_data = iftMergeFileSet(Z1->ref_data, Z2->ref_data);

        #pragma omp parallel for
        for (int s = 0; s < Z3->nsamples; s++)
            Z3->sample[s].id = s;
    }
    // merge CSVs
    else if (Z1->ref_data_type == IFT_REF_DATA_CSV) {
        Z3->ref_data = iftMergeCSVs(Z1->ref_data, Z2->ref_data);

        #pragma omp parallel for
        for (int s = 0; s < Z3->nsamples; s++)
            Z3->sample[s].id = s;
    }
    else {
        iftSetRefData(Z3, Z1->ref_data, Z1->ref_data_type);
    }


    return Z3;
}

void iftMergeDataSetsInPlace(iftDataSet **Z1, iftDataSet *Z2)
{
    if(*Z1 == NULL) {
        *Z1 = iftCopyDataSet(Z2, 1);
    }
    else {
        iftDataSet *dataSetAux = iftMergeDataSets(*Z1, Z2);
        iftDestroyDataSet(Z1);
        *Z1 = dataSetAux;
    }
}


iftDataSet *iftExtractDataSet(const iftDataSet *Z, int begin, int end) {
    if ((begin < 0) || (end >= Z->nsamples) || (begin > end))
        iftError("Invalid range [%d, %d]. Try a range in the interval [0, %d]",
                 "iftCopyDataSetInRange:", begin, end, Z->nsamples-1);

    int nsamples_out = end - begin + 1;
    iftDataSet *Zout = iftCreateDataSet(nsamples_out, Z->nfeats);
    Zout->iftArcWeight  = Z->iftArcWeight;
    Zout->ntrainsamples = 0;
    Zout->function_number = Z->function_number;
    Zout->fsp = iftCopyFeatSpaceParam(Z);

    // copying the first part of the dataset
    iftIntArray *true_labels_Zout = iftCreateIntArray(Zout->nsamples);
    iftIntArray *groups_Zout = iftCreateIntArray(Zout->nsamples);

    int ntrainsamples_out = 0;
    #pragma omp parallel for reduction(+:ntrainsamples_out)
    for (int s = begin; s <= end; s++) {
        int t = s - begin;
        iftCopySample(&Z->sample[s], &Zout->sample[t], Z->nfeats, true);
        true_labels_Zout->val[t] = Zout->sample[t].truelabel;
        groups_Zout->val[t] = Zout->sample[t].group;
        ntrainsamples_out += iftHasSampleStatus(Zout->sample[t], IFT_TRAIN);
    }
    Zout->ntrainsamples = ntrainsamples_out;

    iftIntArray *unique_true_labels_Zout = iftIntArrayUnique(true_labels_Zout->val, true_labels_Zout->n);
    iftIntArray *unique_groups_Zout = iftIntArrayUnique(groups_Zout->val, groups_Zout->n);
    Zout->nclasses = unique_true_labels_Zout->n;
    Zout->ngroups = unique_groups_Zout->n;

    #pragma omp parallel for
    for (int f = 0; f < Zout->nfeats; f++)
        Zout->alpha[f] = Z->alpha[f];

    Zout->ref_data_type = Z->ref_data_type;
    if (Z->ref_data_type == IFT_REF_DATA_FILESET) {
        iftFileSet *ref_data = (iftFileSet *) Z->ref_data;
        iftFileSet *ref_data_out = iftCreateFileSet(Zout->nsamples);
        

        #pragma omp parallel for
        for (int s = begin; s <= end; s++) {
            int t = s - begin;
            ref_data_out->files[t] = iftCreateFile(ref_data->files[s]->path);
            Zout->sample[t].id = t;
        }

        Zout->ref_data = ref_data_out;
    }
    else if (Z->ref_data_type == IFT_REF_DATA_CSV) {
        iftCSV *ref_data = (iftCSV *) Z->ref_data;
        iftCSV *ref_data_out = iftCreateCSV(Zout->nsamples, ref_data->ncols);
        

        #pragma omp parallel for
        for (int s = begin; s <= end; s++) {
            int t = s - begin;
            Zout->sample[t].id = t;

            for (int c = 0; c < ref_data->ncols; c++)
                strcpy(ref_data_out->data[t][c], ref_data->data[s][c]);
        }

        Zout->ref_data = ref_data_out;
    }
    else
        Zout->ref_data = Z->ref_data;

    iftDestroyIntArray(&true_labels_Zout);
    iftDestroyIntArray(&unique_true_labels_Zout);
    iftDestroyIntArray(&groups_Zout);
    iftDestroyIntArray(&unique_groups_Zout);

    return Zout;
}


void iftSplitDataSetAt(const iftDataSet *Z, int sample_idx, iftDataSet **Z1, iftDataSet **Z2) {
    if ((sample_idx <= 0) || (sample_idx >= Z->nsamples))
        iftError("Invalid sample idx for splitting: %d... Try one in [1, nsamples-1]",
                 "iftSplitDataSetAt", sample_idx);
    if (Z1 == NULL)
        iftError("First split dataset's referece Z1 is NULL", "iftSplitDataSetAt");
    if (Z2 == NULL)
        iftError("First split dataset's referece Z2 is NULL", "iftSplitDataSetAt");

    *Z1 = iftExtractDataSet(Z, 0, sample_idx-1);
    *Z2 = iftExtractDataSet(Z, sample_idx, Z->nsamples-1);
}


void iftConqueringOutliersByIFTInMImageDataset(iftDataSet *Z){

    if (Z == NULL)
        iftError("Dataset is NULL", "iftImageDatasetConqueringOutliersByIFT");
    if (Z->ref_data == NULL)
        iftError("Reference Data is NULL", "iftImageDatasetConqueringOutliersByIFT");
    if (Z->ref_data_type != IFT_REF_DATA_MIMAGE)
        iftError("Reference Data Type should be an MImage", "iftImageDatasetConqueringOutliersByIFT");

    iftFHeap  *Q = NULL;
    iftAdjRel *A = NULL;
    float *path_cost_array;
    int      i, p, q, tmp;
    iftVoxel    u, v;

    iftMImage *mimg = (iftMImage*) Z->ref_data;
    if (mimg->zsize >1)
        A = iftSpheric(1.0);
    else
        A = iftCircular(sqrtf(2.0));

    // Initialization
    path_cost_array = iftAllocFloatArray(Z->nsamples);
    Q     = iftCreateFHeap(Z->nsamples, path_cost_array);

    for (p = 0; p < Z->nsamples; p++)
    {
        if (Z->sample[p].label != -1) {
            path_cost_array[p] = 0.0;
            iftInsertFHeap(Q, p);
        } else {
            path_cost_array[p] = IFT_INFINITY_FLT;
        }
    }

    // Image Foresting Transform
    float dist;
    while (!iftEmptyFHeap(Q))
    {
        p = iftRemoveFHeap(Q);
        u = iftMGetVoxelCoord(mimg, p);

        for (i = 1; i < A->n; i++)
        {
            v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                q = iftMGetVoxelIndex(mimg, v);
                if (Q->color[q]!= IFT_BLACK) {
                    if (path_cost_array[q] > path_cost_array[p]) {
                        if (iftDist==NULL)
                            /* note here that we're calculating the arc weight with all the feature vector, if the dataset has
               * features with image coordinates we must set the corresponding alpha values in the dataset before
               * call this function*/
                            dist = Z->iftArcWeight(Z->sample[p].feat,Z->sample[q].feat,Z->alpha,Z->nfeats);
                        else
                            dist = iftDist->distance_table[p][q];
                        tmp = iftMax(path_cost_array[p], dist);
                        if (tmp < path_cost_array[q])  // For this path-value function,
                        {
                            Z->sample[q].label = Z->sample[p].label;
                            path_cost_array[q] = tmp;

                            if(Q->color[q] == IFT_WHITE)
                                // this implies that q has never
                                // been inserted in Q.
                                iftInsertFHeap(Q, q);
                            else
                                iftGoUpFHeap(Q, Q->pos[q]);
                        }
                    }
                }
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftDestroyFHeap(&Q);
    iftFree(path_cost_array);

}


void iftMinMaxFeatureScale(iftDataSet *dataset){

    float *minimus = (float*)iftAlloc(dataset->nfeats,sizeof(float));
    float *maximus = (float*)iftAlloc(dataset->nfeats,sizeof(float));

    for (int j = 0; j < dataset->nfeats; ++j) {
        minimus[j] = dataset->sample[0].feat[j];
        maximus[j] = dataset->sample[0].feat[j];
    }

    for (int j = 0; j < dataset->nfeats; ++j) {
        for (int i = 0; i < dataset->nsamples; ++i) {

            if(minimus[j] > dataset->sample[i].feat[j]){
                minimus[j] = dataset->sample[i].feat[j];
            }

            if(maximus[j] < dataset->sample[i].feat[j]){
                maximus[j] = dataset->sample[i].feat[j];
            }
        }
    }

    for (int j = 0; j < dataset->nfeats; ++j) {
//        printf("min:%f max:%f\n",minimus[j],maximus[j]);
        for (int i = 0; i < dataset->nsamples; ++i) {
            dataset->sample[i].feat[j] = (dataset->sample[i].feat[j] - minimus[j] + 1e-5)/(maximus[j]-minimus[j]+1e-5);
        }
    }

    iftFree(minimus);
    iftFree(maximus);

}

iftDataSet* iftSupervoxelsToLabColorMeanStdAndSizeDataset(iftImage* image, iftImage* label_image){

    int normalization_value = iftNormalizationValue(iftMaximumValue(image));

    if(!iftIsColorImage(image))
        iftError("Color image expected.", "iftSupervoxelsToLabColorMeanStdAndSizeDataset");

    iftVerifyImageDomains(image, label_image,"iftSupervoxelsToLabColorMeanStdAndSizeDataset");

    int nregions = iftMaximumValue(label_image);
    int *size = iftAllocIntArray(nregions);

    iftDataSet* dataset = iftCreateDataSet(nregions, 7);

    iftMImage *lab_img = iftCreateMImage(image->xsize, image->ysize, image->zsize, 3);

    iftColor color;
    iftFColor fcolor;
    for(int p = 0; p < image->n; p++){
        int r = label_image->val[p] - 1;
        if(r < 0)
            continue;

        size[r]++;
        color.val[0] = image->val[p];
        color.val[1] = image->Cb[p];
        color.val[2] = image->Cr[p];

        color = iftYCbCrtoRGB(color,normalization_value);
        fcolor = iftRGBtoLab(color,normalization_value);

        // L in [0, 99.998337]
        // a in [-86.182236, 98.258614]
        // b in [-107.867744, 94.481682]
        lab_img->val[p][0] = (fcolor.val[0] - 0.)/(99.998337 - 0.);
        lab_img->val[p][1] = (fcolor.val[1] + 86.182236)/(98.258614 + 86.182236);
        lab_img->val[p][2] = (fcolor.val[2] + 107.867744)/(94.481682 + 107.867744);

        dataset->sample[r].feat[0] += lab_img->val[p][0];
        dataset->sample[r].feat[1] += lab_img->val[p][1];
        dataset->sample[r].feat[2] += lab_img->val[p][2];
    }

#pragma omp parallel for
    for(int r = 0; r < nregions; r++){
        if(size[r] == 0)
            iftError("Empty region", "iftSupervoxelsToLabColorMeanStdAndSizeDataset");

        dataset->sample[r].id = r+1;

        //Means
        dataset->sample[r].feat[0] /= size[r];
        dataset->sample[r].feat[1] /= size[r];
        dataset->sample[r].feat[2] /= size[r];
        dataset->sample[r].feat[6] = size[r] / (float) image->n; // added
    }

#pragma omp parallel for
    for(int p = 0; p < image->n; p++){
        int r = label_image->val[p] - 1;
        if(r < 0)
            continue;

        dataset->sample[r].feat[3] += pow(dataset->sample[r].feat[0] - lab_img->val[p][0], 2);
        dataset->sample[r].feat[4] += pow(dataset->sample[r].feat[1] - lab_img->val[p][1], 2);
        dataset->sample[r].feat[5] += pow(dataset->sample[r].feat[2] - lab_img->val[p][2], 2);

    }

#pragma omp parallel for
    for(int r = 0; r < nregions; r++){
        //Standard deviation

        dataset->sample[r].feat[3] = sqrt(dataset->sample[r].feat[3]/size[r]);
        dataset->sample[r].feat[4] = sqrt(dataset->sample[r].feat[4]/size[r]);
        dataset->sample[r].feat[5] = sqrt(dataset->sample[r].feat[5]/size[r]);
    }

    dataset->iftArcWeight = iftDistance1;
    iftCopyRefData(dataset, label_image, IFT_REF_DATA_IMAGE);
    // dataset->ref_data = (void*)label_image;

    iftFree(size);
    iftDestroyMImage(&lab_img);

    return dataset;
}

iftDataSet* iftSupervoxelsToYCbCrMeanStdAndSizeDataset(iftImage* image, iftImage* label_image){

    int normalization_value = iftNormalizationValue(iftMaximumValue(image));

    if(!iftIsColorImage(image))
        iftError("Color image expected.", "iftSupervoxelsToLabColorMeanStdAndSizeDataset");

    iftVerifyImageDomains(image, label_image,"iftSupervoxelsToLabColorMeanStdAndSizeDataset");

    int nregions = iftMaximumValue(label_image);
    int *size = iftAllocIntArray(nregions);

    iftDataSet* dataset = iftCreateDataSet(nregions, 7);

    iftMImage *lab_img = iftCreateMImage(image->xsize, image->ysize, image->zsize, 3);

    for(int p = 0; p < image->n; p++){
        int r = label_image->val[p] - 1;
        if(r < 0)
            continue;

        size[r]++;
        lab_img->val[p][0] = ((float)image->val[p])/(float)normalization_value;
        lab_img->val[p][1] = ((float)image->Cb[p])/(float)normalization_value;
        lab_img->val[p][2] = ((float)image->Cr[p])/(float)normalization_value;

        dataset->sample[r].feat[0] += lab_img->val[p][0];
        dataset->sample[r].feat[1] += lab_img->val[p][1];
        dataset->sample[r].feat[2] += lab_img->val[p][2];
    }

#pragma omp parallel for
    for(int r = 0; r < nregions; r++){
        if(size[r] == 0)
            iftError("Empty region", "iftSupervoxelsToLabColorMeanStdAndSizeDataset");

        dataset->sample[r].id = r+1;

        //Means
        dataset->sample[r].feat[0] /= size[r];
        dataset->sample[r].feat[1] /= size[r];
        dataset->sample[r].feat[2] /= size[r];
        dataset->sample[r].feat[6] = size[r] / (float) image->n; // added
    }

#pragma omp parallel for
    for(int p = 0; p < image->n; p++){
        int r = label_image->val[p] - 1;
        if(r < 0)
            continue;

        dataset->sample[r].feat[3] += pow(dataset->sample[r].feat[0] - lab_img->val[p][0], 2);
        dataset->sample[r].feat[4] += pow(dataset->sample[r].feat[1] - lab_img->val[p][1], 2);
        dataset->sample[r].feat[5] += pow(dataset->sample[r].feat[2] - lab_img->val[p][2], 2);
    }

#pragma omp parallel for
    for(int r = 0; r < nregions; r++){
        //Standard deviation
        dataset->sample[r].feat[3] = sqrt(dataset->sample[r].feat[3]/size[r]);
        dataset->sample[r].feat[4] = sqrt(dataset->sample[r].feat[4]/size[r]);
        dataset->sample[r].feat[5] = sqrt(dataset->sample[r].feat[5]/size[r]);
    }

    dataset->iftArcWeight = iftDistance1;
    dataset->ref_data = (void*)label_image;

    iftFree(size);
    iftDestroyMImage(&lab_img);

    return dataset;
}

void iftAttachBinaryLabelImageToDataSet(iftDataSet *Z, iftImage *label)
{
    iftImage *grayImg = NULL;
    if(iftIsColorImage(label)) {
        grayImg = iftImageGray(label);
        iftDestroyImage(&label);
        label = grayImg;
    }

    if(Z->ref_data_type != IFT_REF_DATA_IMAGE && Z->ref_data_type != IFT_REF_DATA_FIMAGE && Z->ref_data_type != IFT_REF_DATA_MIMAGE) {
        iftError("The provided dataset's ref_data_type must be an image", "iftAttachBinaryLabelImageToDataSet");
    }

    bool isBinary = iftIsBinaryImage(label);
    if(!isBinary) {
        iftWarning("The provided label image has more than two labels, a threshold will be applied", "iftAttachBinaryLabelImageToDataSet");
    }

    int maxValue = iftMaximumValue(label);
    iftImage *bin = NULL;
    if (maxValue > 1 || !isBinary){
        bin = iftThreshold(label,128,255,1);
        iftDestroyImage(&label);
        label = bin;
    }

    for(int p = 0; p<Z->nsamples; p++) {
        Z->sample[p].truelabel = label->val[Z->sample[p].id];
    }
    Z->nclasses = 2;
}

iftDataSet *iftRemoveClassFromDataSet(iftDataSet *Z, int classId)
{
    /* create the new dataset */
    int *sampPerClass = iftCountSamplesPerClassDataSet(Z);
    iftDataSet *Z1 = iftCreateDataSet(Z->nsamples - sampPerClass[classId], Z->nfeats);

    Z1->nclasses = Z->nclasses - 1;
    Z1->ngroups = 0; // resets the number of groups
    Z1->function_number = Z->function_number;
    iftCopyRefData(Z1, Z->ref_data, Z->ref_data_type);
    // Z1->ref_data = Z->ref_data;
    // Z1->ref_data_type = Z->ref_data_type;
    iftSetDistanceFunction(Z1, Z->function_number);
    iftCopyFloatArray(Z1->alpha, Z->alpha, Z->nfeats);

    int t = 0;
    for(int s = 0; s < Z->nsamples; s++) {
        if(Z->sample[s].truelabel != classId) {
            iftCopySample(&Z->sample[s], &Z1->sample[t], Z->nfeats, true);
            t++;
        }
    }
    iftUpdateNumbTrainSamples(Z1);

    return Z1;
}


void iftLabelDataSetFromSeeds(iftDataSet *Z, iftLabeledSet *seeds, iftImage* region)
{
    iftImage    *truelabel      = NULL; // to store seeds as a true labeled image

    iftSetStatus(Z,IFT_TRAIN);
    
    if (Z->ntrainsamples != Z->nsamples)
        iftError("It requires that a full training set", "iftLabelDataSetFromSeeds");

    if (Z->ref_data_type != IFT_REF_DATA_IMAGE && Z->ref_data_type != IFT_REF_DATA_FIMAGE && Z->ref_data_type != IFT_REF_DATA_MIMAGE)
        iftError("Reference data must be an iftImage, iftMImage or iftFImage", "iftLabelDataSetFromSeeds");

    if (Z->ref_data_type == IFT_REF_DATA_IMAGE)
        truelabel = iftCreateImageFromImage((iftImage *)Z->ref_data);
    else if (Z->ref_data_type == IFT_REF_DATA_MIMAGE) {
        iftMImage *mimg = (iftMImage *)Z->ref_data;
        truelabel = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    } else if (Z->ref_data_type == IFT_REF_DATA_FIMAGE){
        iftFImage *fimg = (iftFImage *) Z->ref_data;
        truelabel = iftCreateImage(fimg->xsize, fimg->ysize, fimg->zsize);
    }

    iftLabeledSetToImage(seeds, truelabel, true);

    if (region == NULL) /* It only makes sense for pixel datasets that store all image pixels */
    {
        if (Z->nsamples != truelabel->n)
            iftError("No region image has been specified and the dataset does not use all pixels as samples.",
                     "iftLabelDataSetFromSeeds");

        for (int s = 0; s < Z->nsamples; s++) {
            int p = Z->sample[s].id;
            if (truelabel->val[p] > 0){
                Z->sample[s].truelabel = truelabel->val[p];
                iftAddSampleStatus(&Z->sample[s], IFT_SUPERVISED);
            }
        }
        Z->nclasses = iftCountNumberOfClassesDataSet(Z);
        iftDestroyImage(&truelabel);

        return;
    }

    /* Mark regions selected by the seeds */

    int n_labels = iftMaximumValue(truelabel) + 1;
    int n_regions = iftMaximumValue(region);
    int **marked_region = iftAllocIntMatrix(iftMaximumValue(region) + 1, n_labels); // to mark regions selected by the seeds

    iftVerifyImageDomains(region, truelabel, "iftLabelDataSetFromSeeds");
    if (n_regions != Z->nsamples) {
        iftError("Dataset is not compatible with the region image. They must have the same number of elements.",
                 "iftLabelDataSetFromSeeds");
    }

    for (int p = 0; p < truelabel->n; p++) {
        int l = truelabel->val[p];
        if (l > 0){
            marked_region[l][region->val[p]] += 1;
        }
    }

#pragma omp parallel for
    for (int s = 0; s < Z->nsamples; s++) {
        int p = Z->sample[s].id;
        int max_l = -1;
        int max = 0;
        for (int l = 1; l < n_labels; l++) {
            if (marked_region[l][region->val[p]] > max) {
                max = marked_region[l][region->val[p]];
                max_l = l;
            }
        }
        if (max_l > 0) {
            Z->sample[s].truelabel = max_l;
            iftAddSampleStatus(&Z->sample[s], IFT_SUPERVISED);
        }
    }

    Z->nclasses = iftCountNumberOfClassesDataSet(Z);

    iftFreeMatrix(marked_region, n_labels);
    iftDestroyImage(&truelabel);
}


void iftSelectDataSetClassesByCluster(iftDataSet *Z, float threshold)
{
    if (Z->ngroups == 0 || Z->nclasses == 0)
        iftError("Dataset must have labeled classes and groups", "iftSelectDataSetClassesByCluster");

    float **perc = iftAllocFloatMatrix(Z->ngroups + 1, Z->nclasses + 1);
    int *size = iftAllocIntArray(Z->ngroups + 1);

    for (int i = 0; i < Z->nsamples; i++)
    {
        int grp = Z->sample[i].group;
        int lb = Z->sample[i].truelabel;
        size[grp] += 1;
        perc[lb][grp] += 1;
    }

#pragma omp parallel for
    for (int i = 1; i <= Z->ngroups; i++)
    {
        for (int lb = 0; lb <= Z->nclasses; lb++) {
            perc[lb][i] /= size[i];
        }
    }

    for (int i = 0; i < Z->nsamples; i++)
    {
        int grp = Z->sample[i].group;
        int lb = Z->sample[i].truelabel;

        if ( (perc[lb][grp] + perc[0][grp]) < threshold) {
            Z->sample[i].truelabel = 0;
            iftRemoveSampleStatus(&Z->sample[i], IFT_TRAIN);

            if (iftHasSampleStatus(Z->sample[i], IFT_LABELPROPAGATED))
                iftRemoveSampleStatus(&Z->sample[i], IFT_LABELPROPAGATED);
        }
    }

    Z->nclasses = iftCountNumberOfClassesDataSet(Z);

    iftFreeMatrix(perc, Z->nclasses + 1);
    iftFree(size);
}


iftDataSet *iftBoundingBoxArrayToDataSet(const iftImage *img, const iftBoundingBoxArray *bb_ary, int n_feats)
{
    iftDataSet *Z = iftCreateDataSet(bb_ary->n, n_feats);

#pragma omp parallel for
    for (int i = 0; i < bb_ary->n; i++) {
        float *feat = Z->sample->feat;
        iftBoundingBox *bb = &bb_ary->val[i];
        int j = 0;
        for (int x = bb->begin.x; x < bb->end.x; x++) {
            for (int y = bb->begin.y; y < bb->end.y; y++) {
                for (int z = bb->begin.z; z < bb->end.z; z++) {
                    iftVoxel v = {.x = x, .y = y, .z = z};
                    int p = iftGetVoxelIndex(img, v);
                    feat[j] = img->val[p];
                    j++;
                }
            }
        }
    }

    return Z;
}


iftMatrix *iftHighDimensionalDistanceEstimationBasedOnPCA(iftDataSet *Z, int k, int n)
{
    iftKnnGraph *graph          = NULL;
    iftMatrix *U                = NULL;
    iftMatrix *S                = NULL;
    iftMatrix *Vt               = NULL;
    iftMatrix *dist             = NULL;
    iftMatrix *samples          = NULL;
    iftMatrix *samplesT         = NULL;
    iftMatrix *vector_basis     = NULL;
    iftMatrix *projection_coef  = NULL;
    iftMatrix *sum              = NULL;
    iftIntArray *neighbor_index = NULL;
    iftSample ref_sample;
    int f,i,j,l;
    int extended_k = k*3;
    /* We initially selected three times the number of neighbors because
     * Z is still in the high dimensional space> In this space we can not
     * be sure whether a sample is adjacent to another since they are all
     * equidistant. */

    if ((k < 0) || (extended_k > Z->nsamples))
        k = Z->nsamples;

    // initializing
    graph = iftCreateKnnGraph(Z,extended_k);
    dist  = iftCreateMatrix(Z->nsamples,Z->nsamples);

    //#pragma omp parallel for private(samples,samplesT,neighbor_index,ref_sample,U,S,Vt,vector_basis,projection_coef,sum,f,i,j,l)
    for (i = 0; i < Z->nsamples; i++){

        // Transforming the node's nearest neighbors into a matrix
        samples = iftKnnNodeToMatrix(graph,i);

        // Saving reference sample
        ref_sample.feat = (float *)calloc(Z->nfeats,sizeof(float));
        iftCopySample(&(Z->sample[i]),&ref_sample,Z->nfeats,true);

        // Centralizing matrix
        for (j = 0; j < samples->ncols; j++)
            for (l = 0; l < samples->nrows; l++)
                samples->val[iftGetMatrixIndex(samples, j, l)] -= ref_sample.feat[j];

        samplesT = iftTransposeMatrix(samples);

        // Defining hyperplane based on PCA
        iftSingleValueDecomp(samples,&U,&S,&Vt);
        iftDestroyMatrix(&samples);
        iftDestroyMatrix(&U);
        iftDestroyMatrix(&S);

        // Selecting eigenvectors from highest eigenvalues
        vector_basis = iftCreateMatrix(Vt->ncols,n);
        for (f = 0; f < n; f++)
            for (l = 0; l < Vt->ncols; l++)
                vector_basis->val[iftGetMatrixIndex(vector_basis,l,f)] = Vt->val[iftGetMatrixIndex(Vt,l,f)];
        iftDestroyMatrix(&Vt);

        // Computing linear dependency coefficients
        projection_coef = iftMultMatrices(vector_basis,samplesT);

        // Computing distance matrix
        // -- square matrix
        for (f = 0; f < projection_coef->ncols; f++)
            for (l = 0; l < projection_coef->nrows; l++)
                projection_coef->val[iftGetMatrixIndex(projection_coef,f,l)] *= projection_coef->val[iftGetMatrixIndex(projection_coef,f,l)];
        // -- sum for each sample
        sum = iftMatrixSumColumn(projection_coef);

        // Updating nearest neighbors and selecting the k nearest neighbors
        // -- Getting neighbor nodes indexes
        neighbor_index = iftGetKnnNeighborSamplesIndex(graph,i);
        iftMatrix *sorted_sum = iftCopyMatrix(sum);
        iftFQuickSort(sorted_sum->val,neighbor_index->val,0,neighbor_index->n-1,IFT_INCREASING);
        iftDestroyMatrix(&sum);

        // Creating sparse distance matrix
        for (f = 0; f < k; f++){
            dist->val[iftGetMatrixIndex(dist,neighbor_index->val[f],i)] = sorted_sum->val[f];
        }

        // Deallocating memory
        free(ref_sample.feat);
        iftDestroyIntArray(&neighbor_index);
        iftDestroyMatrix(&vector_basis);
        iftDestroyMatrix(&projection_coef);
        iftDestroyMatrix(&sorted_sum);
        iftDestroyMatrix(&samplesT);
    }
    iftDestroyKnnGraph(&graph);

    return dist;
}

int iftSampleNearestNeighbor(iftDataSet *Z, iftSample *s){

    int nearest_sample_id = 0;

    /* Compute distance to all samples of Z */
    float *dist = (float *)calloc(sizeof(float),Z->nsamples);
    for (int si = 0; si < Z->nsamples; si++){
        dist[si] = iftFeatDistance(Z->sample[si].feat, s->feat,Z->nfeats);
    }

    float min = iftMaxFloatArray(dist,Z->nsamples);
    for (int s = 0; s < Z->nsamples; s++){
        if (dist[s] < min){
            min = dist[s];
            nearest_sample_id = s+1;
        }
    }

    free(dist);

    return nearest_sample_id;
}

iftFloatArray *iftSampleDistanceToGroups(iftDataSet *Z, iftSample *s, iftIntArray **group_idx, int cluster_representative, bool sort)
{
    *group_idx = iftCreateIntArray(Z->ngroups);
    iftFloatArray *dists = iftCreateFloatArray(Z->ngroups);

    /* Compute distance to all samples of Z */
    float *dist = NULL;

    if (cluster_representative == 1) {
        dist = (float *) calloc(sizeof(float), Z->nsamples);

        for (int si = 0; si < Z->nsamples; si++){
            dist[si] = iftFeatDistance(Z->sample[si].feat, s->feat,Z->nfeats);
        }

        for (int g = 1; g <= Z->ngroups; g++) {
            float min = iftMaxFloatArray(dist,Z->nsamples);
            (*group_idx)->val[g] = g;
            for (int si = 0; si < Z->nsamples; si++) {
                if ((dist[si] < min) && (Z->sample[si].group == g)) {
                    min = dist[si];
                }
            }
            dists->val[g] = min;
        }

    } else if (cluster_representative == 2) {
        int pos = 0;
        for (int si = 0; si < Z->nsamples; si++){
            if (iftHasSampleStatus(Z->sample[si],IFT_PROTOTYPE)) {
                dists->val[pos++] = iftFeatDistance(Z->sample[si].feat, s->feat, Z->nfeats);
                (*group_idx)->val[Z->sample[si].group] = Z->sample[si].group;
            }
        }
    } else if (cluster_representative == 3) {
        iftDataSet *Zaux = iftExtractCentroidsFromDataSetAsDataSet(Z,FALSE,FALSE);
        int pos = 0;
        for (int si = 0; si < Zaux->nsamples; si++){
            dists->val[pos++] = iftFeatDistance(Zaux->sample[si].feat, s->feat, Zaux->nfeats);
            (*group_idx)->val[Z->sample[si].group] = Zaux->sample[si].group;
        }
    }

    if (sort)
        iftFQuickSort(dists->val,(*group_idx)->val,0,Z->ngroups,IFT_INCREASING);

    return dists;
}

iftDataSet *iftPatchesFromSuperpixels(char *input_dir, char *output, int nsuperpixels, int patch_size)
{
    timer *tstart;
    tstart = iftTic();

    /* Load input image/activation set */

    iftFileSet *fs_input = iftLoadFileSetFromDirBySuffix(input_dir, ".png", 1);
    int nimages = fs_input->n;

    iftMImage **input;
    int first = 0;
    int last = fs_input->n - 1;

    input = (iftMImage **)calloc(fs_input->n, sizeof(iftMImage *));

    for (int i = first; i <= last; i++)
    {
        // printf("Reading image: %d of %d\r", i + 1, (int)fs->n);
        // fflush(stdout);
      input[i] = iftImageToMImage(iftReadImageByExt(fs_input->files[i]->path), RGB_CSPACE);
    }
    // printf("\n");

    int nsupervoxels = nsuperpixels;
    int nsamples = nsuperpixels * nimages;
    int nclasses = iftFileSetLabelsNumber(fs_input);
    iftAdjRel *Q = NULL;
    char *basename1, *basename2;
    basename2 = iftFilename(output, ".zip");
    iftMakeDir(basename2);
    char filename[200];
    sprintf(filename, "%s.csv", basename2);
    FILE *fp = fopen(filename, "w");

    Q = iftRectangular(patch_size, patch_size);

    int nfeats = input[0]->m * Q->n;

    /* Create supervoxel dataset */

    iftDataSet *Z = iftCreateDataSet(nsamples, nfeats);
    printf("%d samples, %d feats\n", nsamples, nfeats);
    iftFileSet *fsRef = iftCreateFileSet(nsamples);

    iftSetStatus(Z, IFT_TRAIN);
    Z->nclasses = nclasses;

    int s = 0;
    // iterates through each image in input fileset
    for (int i = 0; i < nimages; i++)
    {
        basename1 = iftFilename(fs_input->files[i]->path, ".png");

        // printf("Processing file %s: %d of %d files\r", basename1, i + 1, nimages);
        // fflush(stdout);

        // true label of current file
        int truelabel = fs_input->files[i]->label;

        // It is mandatory to include the following line of code in order to select random initial centers.
        iftRandomSeed(time(NULL));

        iftSet *S;

        sprintf(filename, "%s/%s.%s", basename2, basename1, "png");

        iftImage *mask1, *seeds;
        iftAdjRel *A_;
        iftIGraph *igraph;

        /* Compute ISF superpixels */
        if (iftIs3DMImage(input[i]))
        {
            A_ = iftSpheric(1.0);
        }
        else
        {
            A_ = iftCircular(1.0);
        }

        mask1 = iftSelectImageDomain(input[i]->xsize, input[i]->ysize, input[i]->zsize);

        /* minima of a basins manifold in that domain */
        igraph = iftImplicitIGraph(input[i], mask1, A_);

        /* seed sampling for ISF */
        seeds = iftGridSampling(input[i], mask1, nsuperpixels);

        iftNumberOfElements(seeds);

        iftDestroyImage(&mask1);

        iftIGraphISF_Root(igraph, seeds, 0.1, 12, 2);

        for (int t = 0; t < igraph->nnodes; t++)
        {
            int p = igraph->node[t].voxel;
            int r = igraph->root[p];
            if (p == r)
            {
                iftInsertSet(&S, p);
            }
        }

        // printf("ISF proc time im ms: %f\n", iftCompTime(t1,t2));

        iftDestroyImage(&seeds);
        iftDestroyIGraph(&igraph);
        iftDestroyAdjRel(&A_);

        int j = 0;
        int invalid_voxel = 0;

        // for each supervoxel
        while ((j < nsupervoxels) && (S != NULL))
        {
            int p = iftRemoveSet(&S);
            iftVoxel u = iftMGetVoxelCoord(input[i], p);
            invalid_voxel = 0;
            int f = 0;
            // for each pixel in square patch
            for (int k = 0; k < Q->n; k++)
            {
                iftVoxel v = iftGetAdjacentVoxel(Q, u, k);
                if (iftMValidVoxel(input[i], v) == 1)
                {
                    int q = iftMGetVoxelIndex(input[i], v);
                    for (int b = 0; b < input[i]->m; b++)
                    {
                        Z->sample[s].feat[f] = input[i]->val[q][b];
                        f++;
                    }
                }
                else
                {
                    invalid_voxel = 1;
                    break;
                }
            }
            if (invalid_voxel == 0)
            {
                sprintf(filename, "%s-%03d-%03d-%03d-%06d", basename1, u.x, u.y, patch_size, s);
                fprintf(fp, "%s\n", filename);
                Z->sample[s].id = s;
                Z->sample[s].truelabel = truelabel;
                sprintf(filename, "%s/%s.%s", input_dir, basename1, "png");
                fsRef->files[s] = iftCreateFile(iftAbsPathname(filename));
                s++;
            }
            else
            { // do not include invalid voxels
                fsRef->n = s;
                Z->nsamples = s;
            }
            j++;
        }
    }

    for (int s = 0; s < Z->nsamples; s++)
        Z->sample[s].group = s + 1;
    Z->ngroups = Z->nsamples;

    // iftDestroyAdjRel(&A);
    // iftDestroyAdjRel(&Q);
    fclose(fp);
    iftAddStatus(Z, IFT_SUPERVISED);

    printf("\nDataset with %d samples\n", Z->nsamples);

    
    printf("%ld\n",fsRef->n);

    iftSetRefData(Z, fsRef, IFT_REF_DATA_FILESET);

    printf("Writing dataset...\n");
    // iftWriteDataSet(Z, output);

    printf("Destroying images\n");
    // destroy mimage set
    // iftMImage **aux = input;

    // for (int i = 0; i < nimages; i++)
    //     iftDestroyMImage(&aux[i]);
    // free(aux);
    // *input = NULL;

    // iftDestroyFileSet(&fs_input);
    // iftDestroyFileSet(&fsRef);
    // iftFree(basename2);

    printf("Done ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));

    return Z;
}

iftDataSet *iftDataSetFromAllSeeds(char *markers_dir, char *mimages_dir, iftAdjRel *A){
  iftFileSet *fs_seeds = iftLoadFileSetFromDirBySuffix(markers_dir, "-seeds.txt", 1);
  int nimages          = fs_seeds->n; 
  char filename[200];
  int nfeats           = 0, nsamples = 0;
  
  for (int i = 0; i < nimages; i++) {
    char *basename   = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
    sprintf(filename, "%s/%s.mimg", mimages_dir, basename);
    iftMImage *input = iftReadMImage(filename);
    if ((nfeats != 0)&&(nfeats != A->n * input->m))
      iftError("All images should have the same number of channels","iftDataSetFromAllSeeds");    
    nfeats           = A->n * input->m;
    iftLabeledSet *S = iftMReadSeeds(input, fs_seeds->files[i]->path);
    nsamples        += iftLabeledSetSize(S);
    iftFree(basename);
    iftDestroyMImage(&input);
    iftDestroyLabeledSet(&S);
  }

  /* create a dataset with patches from seeds selected in all
     training images */
	
  iftDataSet *Z = iftCreateDataSet(nsamples,nfeats);
  int s         = 0, incr = 0;
	
  for (int i = 0; i < nimages; i++) {
    char *basename   = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
    sprintf(filename, "%s/%s.mimg", mimages_dir, basename);
    iftMImage *input = iftReadMImage(filename);
    iftLabeledSet *S = iftMReadSeeds(input, fs_seeds->files[i]->path);
    while (S != NULL) {
      int label;
      int p = iftRemoveLabeledSet(&S, &label);
      if (label==0)
	incr=1;
      Z->sample[s].id = s;
      Z->sample[s].truelabel = label;
      if (Z->sample[s].truelabel > Z->nclasses)
	Z->nclasses = Z->sample[s].truelabel;
      iftVoxel u = iftMGetVoxelCoord(input, p);
      int j      = 0;
      for (int k = 0; k < A->n; k++) {
	iftVoxel v = iftGetAdjacentVoxel(A, u, k);
	if (iftMValidVoxel(input, v)) {
	  int q = iftMGetVoxelIndex(input, v);
	  for (int b = 0; b < input->m; b++) {
	    Z->sample[s].feat[j] = input->val[q][b];
	    j++;
	  }
	} else {
	  for (int b = 0; b < input->m; b++) {
	    Z->sample[s].feat[j] = 0;
	    j++;
	  }
	}
      }
      s++;
    }
    iftFree(basename);
    iftDestroyMImage(&input);
  }
  
  if (incr == 1){/* correct true labels and number of classes */
    for (s=0; s < Z->nsamples; s++)
      Z->sample[s].truelabel += 1;
    Z->nclasses += 1;
  }
  
  iftSetStatus(Z, IFT_TRAIN);
  iftAddStatus(Z, IFT_SUPERVISED);
  iftDestroyFileSet(&fs_seeds);

  return(Z);
}
