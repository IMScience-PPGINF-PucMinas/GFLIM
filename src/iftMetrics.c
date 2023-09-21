#include "iftMetrics.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/tools/String.h"


/****************************PRIVATE FUNCTIONS**********************/
float iftTruePositiveRateGivenErrors(iftErrorClassification* error)
{
    return ((float)error->tp)/((float)error->tp+error->fn);
}

float iftFalsePositiveRateGivenErrors(iftErrorClassification* error)
{
    return ((float)error->fp)/((float)error->fp+error->tn);
}

/********************** PUBLIC FUNCTIONS **********************/
float iftPrecisionGivenErrors(iftErrorClassification* error) {
    return ((float) error->tp) / (error->tp + error->fp);
}


float iftRecallGivenErrors(iftErrorClassification* error) {
    return ((float) error->tp) / (error->tp + error->fn);
}

float iftAccuracyGivenErrors(iftErrorClassification* error) {
    return ((float) error->tp + error->tn)
            / (error->tp + error->fp + error->tn + error->fn);
}

float iftProbabilityRatioGivenErrors(iftErrorClassification* error)
{
    return iftTruePositiveRateGivenErrors(error)/iftFalsePositiveRateGivenErrors(error);
}

float iftOddsRatioGivenErrors(iftErrorClassification* error)
{
    return ((float)error->tp*error->tn)/((float)error->fp*error->fn);
}

float iftFScoreGivenErrors(iftErrorClassification* error) {
    float precision = iftPrecisionGivenErrors(error);
    float recall = iftRecallGivenErrors(error);

    return (2 * precision * recall) / (precision + recall);
}

double iftSumDoubleArray (double *array, int length) {
    int i = 0;
    double sum = 0;
    for (i = 0; i < length; i++)
        sum += array[i];
    return (double) sum;
}

float iftSumFloatArray (float *array, int length) {
    int i = 0;
    float sum = 0;
    for (i = 0; i < length; i++)
        sum += array[i];
    return (float) sum;
}

int iftSumIntArray (int *array, int length) {
    int i = 0;
    int sum = 0;
    for (i = 0; i < length; i++)
        sum += array[i];
    return sum;
}

int iftMaxIndexFloatArray(float *array, int length) {
  float max = array[0];
  float max_index = 0;

  int i;
  for (i = 1; i < length; ++i) {
      if (array[i] > max) {
        max = array[i];
        max_index = i;
      }
  }

  return max_index;
}

int iftMinIndexFloatArray(float *array, int length) {
  float min = array[0];
  float min_index = 0;
  
  int i;
  for (i = 1; i<length; ++i) {
    if(array[i]<min) {
        min = array[i];
        min_index = i;
    }
  }

  return min_index;
}

float iftMaxFloatArray(float* array, int length) {
    float max = array[0];

    int i;
    for (i = 1; i<length; ++i) {
        if(array[i]>max) max = array[i];
    }

    return max;
}

float iftMinFloatArray(float* array, int length) {
    float min = array[0];

    int i;
    for (i = 1; i<length; ++i) {
        if(array[i]<min) min = array[i];
    }

    return min;
}

int iftMedianIntArray(int* array, int length) {
    int* sorted = iftAllocIntArray(length);
    iftIntArray* index = iftIntRange(0, length, 1);

    iftCopyIntArray(sorted, array, length);

    iftBucketSort(sorted, index->val, length, IFT_INCREASING);

    int median;
    if (length%2==0) {
        float m = 0;
        m = ((float)sorted[length/2] + (float)sorted[(length/2)-1])/2.0;
        median = (int) roundf(m);
    }
    else {
        median = sorted[length/2];
    }

    iftFree(sorted);
    iftDestroyIntArray(&index);

    return  median;
}

float iftWeightedMedianFloatArray(float* array, int length, float* weight) {
    float* sorted = iftAllocFloatArray(length);
    iftIntArray* index = iftIntRange(0, length, 1);

    iftCopyFloatArray(sorted, array, length);

    iftFHeapSort(sorted, index->val, length, IFT_INCREASING);

    float sum = iftSumFloatArray(weight, length);

    float acc = 0.0;
    int i;

    for (i = 0; i < length; ++i) {
        acc += weight[index->val[i]]/sum;
        if(acc>=0.5)
            break;
    }

    float median = sorted[i];

    iftFree(sorted);
    iftDestroyIntArray(&index);

    return  median;
}

float iftMedianFloatArray(float* array, int length) {
    float* sorted = iftAllocFloatArray(length);
    iftIntArray* index = iftIntRange(0, length, 1);

    iftCopyFloatArray(sorted, array, length);

    iftFHeapSort(sorted, index->val, length, IFT_INCREASING);

    float median = 0;
    if (length%2==0) {
        median = (sorted[length/2] + sorted[(length/2)-1])/2.0;
    }
    else {
        median = sorted[length/2];
    }

    iftFree(sorted);
    iftDestroyIntArray(&index);

    return  median;
}

float iftWeightedMeanFloatArray (float *array, int length, float* weight) {
    float mean = 0.0;
    float sum = iftSumFloatArray(weight, length);

    for (int i = 0; i < length; ++i) {
        mean += array[i]*weight[i]/sum;
    }

    return mean;
}

float iftMeanFloatArray (float *array, int length) {
    return (float) iftSumFloatArray(array, length) * 1.0 / length * 1.0;
}

void iftMultScalarFloatArray(float* array, int length, float scalar) {
    for (int i = 0; i < length; ++i) {
        array[i] *= scalar;
    }
}

/**
@brief Compute the average of a float array disregarding positions with value zero.\n

Given an float array of size N, this function sums all values different than 0 and divides by the number
of values different than zero.\n
 
@param  array                 Float *.    Float array;
@param  length                Integer.    Size of the float array

@return The average of all values different than zero.
*/
float iftMeanFloatArrayDiffZero(float *array, int length)
{
    int i = 0, counter=0;
    float sum = 0.0;
    for (i = 0; i < length; i++)
    {
        if (array[i] != 0)
        {
            sum += array[i];
            counter++;
        }
    }
    if (counter == 0)
        return 0.0;
    return (float) (sum/counter) * 1.0;
}

float iftVarianceFloatArray (float *array, int length) {
    float average = iftMeanFloatArray(array, length);
    float accumulator = 0;
    float variance = 0;
    int i = 0;
    for (i = 0; i < length; i++)
        accumulator += pow((array[i] - average), 2);
    variance = (float) (accumulator * 1.0) / (length * 1.0);
    return variance;
}

float iftStddevFloatArray (float *array, int length) {
    return sqrt(iftVarianceFloatArray(array, length));
}

bool iftCompareFloatArrays (float *array1, float *array2, int length) {
    int i = 0;
    bool equal = true;
    for (i = 0; i < length && equal; i++)
        equal = (array1[i] == array2[i]);
        
    return equal;
}

bool iftCompareIntArrays (int *array1, int *array2, int length) {
    int i = 0;
    bool equal = true;
    for (i = 0; i < length && equal; i++)
        equal = (array1[i] == array2[i]);
        
    return equal;
}

float  iftComputeMetricInDataSetByName(iftDataSet *Z, char *metric) {
    float metricVal = 0.0f;

    if(iftCompareStrings(metric, "fscore"))
    {
        if(Z->nclasses == 2)
            metricVal = iftFscore(Z, 1);
        else
            metricVal = iftFscoreMulticlass(Z);
    } else if (iftCompareStrings(metric, "true_positives"))
    {
        metricVal = iftTruePositives(Z);
    } else if (iftCompareStrings(metric, "kappa"))
    {
        metricVal = iftCohenKappaScore(Z);
    } else if (iftCompareStrings(metric, "skewed_true_positives"))
    {
        metricVal = iftSkewedTruePositives(Z);
    } else if (iftCompareStrings(metric, "normalized_accuracy"))
    {
        metricVal = iftNormAccuracy(Z);
    } else if (iftCompareStrings(metric, "weighted_accuracy"))
    {
        metricVal = iftClassifAccuracy(Z);
    } else if (iftCompareStrings(metric, "opf_accuracy"))
    {
        metricVal = iftOPFAccuracy(Z, false);
    }

    return metricVal;
}

float iftFscore(iftDataSet *Z, int PositiveClass) {
    float *ncorrect = NULL, TP, FN, FP, Rec, Prec, Fscore = 0.0f;
    int i, *nclass, ntestsamples = 0, NegativeClass;

    if (Z->nclasses != 2)
        iftError("This dataset does not have 2 classes", "iftFscore");

    ncorrect = iftAllocFloatArray(Z->nclasses + 1);
    nclass = iftAllocIntArray(Z->nclasses + 1);

    for (i = 0; i < Z->nsamples; i++) {
        if (iftHasSampleStatus(Z->sample[i], IFT_TEST)) {
            nclass[Z->sample[i].truelabel]++;
            ntestsamples++;
            if (Z->sample[i].truelabel == Z->sample[i].label) {
                ncorrect[Z->sample[i].truelabel]++;
            }
        }
    }

    if (PositiveClass == 1)
        NegativeClass = 2;
    else
        NegativeClass = 1;

    TP = ncorrect[PositiveClass];
    printf("- TP: %f\n", TP);
    // TN = ncorrect[NegativeClass];
    FN = (nclass[PositiveClass] - ncorrect[PositiveClass]);
    printf("- FN: %f\n", FN);
    FP = (nclass[NegativeClass] - ncorrect[NegativeClass]);
    printf("- FP: %f\n", FP);

    Rec = TP / (TP + FN);
    Prec = TP / (TP + FP);
    Fscore = 2 * (Rec * Prec) / (Rec + Prec);

    for (i = 1; i <= Z->nclasses; i++) {
        printf("cl %d = %f (%d/%d)\n", i, ncorrect[i] / (float) nclass[i],
                (int) ncorrect[i], nclass[i]);
    }
    printf("rec: %.2f, prec: %.2f\n", Rec, Prec);

    iftFree(ncorrect);
    iftFree(nclass);

    return (Fscore);
}

float iftFscoreMulticlass(iftDataSet *Z) {
    iftMatrix* confusion = iftComputeConfusionMatrix(Z, false);

    float fScore = 0.0f;
    /* compute the F-Score for each class */
    for(int cla = 1; cla <= Z->nclasses; cla++) {
        int TP = iftMatrixElem(confusion, cla-1, cla-1);
        int FP = 0, FN = 0, TN = 0;
        for(int col = 1; col <= Z->nclasses; col++) {
            if(col != cla) {
                FP += iftMatrixElem(confusion, col-1, cla-1);
                FN += iftMatrixElem(confusion, cla-1, col-1);
                for(int row = 1; row <= Z->nclasses; row++) {
                    if(row != cla) {
                        TN += iftMatrixElem(confusion, col-1, row-1);
                    }
                }
            }
        }
        float precision = (float)(TP)/(float)(TP + FP);
        float recall = (float)(TP)/(float)(TP + FN);
        fScore += (2.0 * (precision*recall) / (precision + recall));
    }
    fScore /= (float)Z->nclasses;

    return fScore;
}


float iftTruePositives(iftDataSet *Z) {
    int i, ncorrectsamples = 0, ntestsamples = 0;

    for (i = 0; i < Z->nsamples; i++) {
      if (iftHasSampleStatus(Z->sample[i], IFT_TEST)) {
        ntestsamples++;
        if (Z->sample[i].truelabel == Z->sample[i].label) {
          ncorrectsamples++;
        }
      }
    }

    return ((float) ncorrectsamples / (float) ntestsamples);
}

float iftCohenKappaScore(iftDataSet* Z) {

    iftMatrix* confusion = iftComputeConfusionMatrix(Z, false);

    iftWriteMatrix(confusion,"confusion.csv");
    
    double sum = iftMatrixSum(confusion);

    iftMultScalarFloatArray(confusion->val, confusion->n, 1.0/(sum+0.0000001));

    float observed = 0.0;

    for (int r = 0; r < confusion->nrows; ++r) {
        observed += iftMatrixElem(confusion, r, r);
    }

    iftMatrix* col = iftMatrixSumColumn(confusion);
    iftMatrix* row = iftMatrixSumRow(confusion);

    float expected = 0.0;

    for (int r = 0; r < confusion->nrows; ++r) {
        expected += col->val[r] * row->val[r];
    }

    iftDestroyMatrix(&row);
    iftDestroyMatrix(&col);
    iftDestroyMatrix(&confusion);

    return (observed - expected) / (1 - expected + 0.0000001);
}


float iftSkewedTruePositives(iftDataSet *Z) {
    float *TP = NULL, Acc = 0.0f, tot_weight = 0.0f;
    int i, *nclass, ntestsamples = 0;

    TP = iftAllocFloatArray(Z->nclasses + 1);
    nclass = iftAllocIntArray(Z->nclasses + 1);

    for (i = 0; i < Z->nsamples; i++) {
      if (!iftHasSampleStatus(Z->sample[i], IFT_TRAIN)) {
        nclass[Z->sample[i].truelabel]++;
        ntestsamples++;
        if (Z->sample[i].truelabel == Z->sample[i].label) {
          TP[Z->sample[i].truelabel]++;
        }
      }
    }

    for (i = 1; i <= Z->nclasses; i++) {
      /* printf("cl %d = %f (%d/%d)\n", i, TP[i] / (float) nclass[i], */
      /*     (int) TP[i], nclass[i]); */
      Acc += (TP[i] / (float) nclass[i])
        * (1.0 - ((float) nclass[i] / ntestsamples));
      tot_weight += (1.0 - ((float) nclass[i] / ntestsamples));
    }

    Acc /= tot_weight;
    iftFree(TP);
    iftFree(nclass);

    return (Acc);
}

float iftNormAccuracy(iftDataSet* Z) {
    iftFloatArray* truePositive = iftCreateFloatArray(Z->nclasses);
    iftIntArray* nclass = iftCreateIntArray(Z->nclasses);
    for (int i = 0; i < Z->nsamples; ++i) {
        if(iftHasSampleStatus(Z->sample[i], IFT_TEST)) {
            int trueLabel = Z->sample[i].truelabel;
            nclass->val[trueLabel-1]++;
            if(Z->sample[i].label == trueLabel)
                truePositive->val[trueLabel-1]+=1.0;
        }
    }
    for (int i = 0; i < truePositive->n; ++i) {
        if(nclass->val[i] == 0){
            truePositive->val[i] = 0.0;
        }else{
            truePositive->val[i]/=nclass->val[i];
        }

    }
	float normAcc = iftMeanFloatArray(truePositive->val, truePositive->n);
    iftDestroyIntArray(&nclass);
    iftDestroyFloatArray(&truePositive);

    return normAcc;
}


float iftClassifAccuracy(iftDataSet *Z) {
    float Acc = 0.0f, **error_matrix = NULL;
    float FP = 0.0f, FN = 0.0f;
    int i, *nclass = NULL, ntestsamples = 0;

    error_matrix = (float **) iftAlloc(Z->nclasses + 1, sizeof(float *));
    if (error_matrix == NULL)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftClassifAccuracy");
    for (i = 0; i <= Z->nclasses; i++)
        error_matrix[i] = iftAllocFloatArray(2);

    nclass = iftAllocIntArray(Z->nclasses + 1);

    for (i = 0; i < Z->nsamples; i++) {
        if (!iftHasSampleStatus(Z->sample[i], IFT_TRAIN)) {
            nclass[Z->sample[i].truelabel]++;
            ntestsamples++;
        }
    }

    for (i = 0; i < Z->nsamples; i++) {
        if (!iftHasSampleStatus(Z->sample[i], IFT_TRAIN)) {
            if (Z->sample[i].truelabel != Z->sample[i].label) {
                error_matrix[Z->sample[i].truelabel][0]++; // FN
                error_matrix[Z->sample[i].label][1]++; // FP
            }
        }
    }

    for (i = 1; i <= Z->nclasses; i++) {
        if (nclass[i] != 0) {
            error_matrix[i][0] /= nclass[i];
        }
        if ((ntestsamples - nclass[i]) != 0) {
            error_matrix[i][1] /= (float) (ntestsamples - nclass[i]);
        }
        FN += error_matrix[i][0];
        FP += error_matrix[i][1];
    }

    FP = (FP / Z->nclasses);
    FN = (FN / Z->nclasses);

    Acc = 1.0 - iftMax(FP, FN);

    for (i = 0; i <= Z->nclasses; i++)
        iftFree(error_matrix[i]);
    iftFree(error_matrix);
    iftFree(nclass);

    return (Acc);
}

/**
 * Compute the OPF accuracy of the given dataset, considering only
 * samples with status != IFT_TRAIN.
 * See http://www.ic.unicamp.br/~afalcao/libopf/node3.html
 * @param Z Input dataset.
 * @param printAccByClass Whether or not to print the accuracy of each class
 * @return OPF accuracy.
 */
float iftOPFAccuracy(iftDataSet *Z, bool printAccByClass) {
    float **error_matrix = NULL;
    int i, *nclass = NULL, ntestsamples = 0;

    error_matrix = (float **) iftAlloc(Z->nclasses + 1, sizeof(float *));
    if (error_matrix == NULL)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftOPFAccuracy");
    for (i = 0; i <= Z->nclasses; i++)
        error_matrix[i] = iftAllocFloatArray(2);

    nclass = iftAllocIntArray(Z->nclasses + 1);

    for (i = 0; i < Z->nsamples; i++) {
        if (!iftHasSampleStatus(Z->sample[i], IFT_TRAIN)) {
            nclass[Z->sample[i].truelabel]++;
            ntestsamples++;
        }
    }

    for (i = 0; i < Z->nsamples; i++) {
        if (!iftHasSampleStatus(Z->sample[i], IFT_TRAIN)) {
            if (Z->sample[i].truelabel != Z->sample[i].label) {
                error_matrix[Z->sample[i].truelabel][1]++; // FN
                error_matrix[Z->sample[i].label][0]++; // FP
            }
        }
    }

    int nlabels = 0;
    float error = 0.0f;
    for (i = 1; i <= Z->nclasses; i++) {

        if (nclass[i] != 0) {
            error_matrix[i][1] /= (float) nclass[i];
            nlabels++;
        }
        if ((ntestsamples - nclass[i]) != 0) {
            error_matrix[i][0] /= (float) (ntestsamples - nclass[i]);
        }
        if (nclass[i] != 0) {
            error += (error_matrix[i][0] + error_matrix[i][1]);
        }
        if(printAccByClass)
            printf("Class %d acc: %f\n", i, (1.0 - (error_matrix[i][0] + error_matrix[i][1])/2.0));
    }

    float Acc = 1.0 - (error / (2.0f * nlabels));

    for (i = 0; i <= Z->nclasses; i++)
        iftFree(error_matrix[i]);
    iftFree(error_matrix);
    iftFree(nclass);

    return (Acc);
}

iftFloatArray *iftTruePositivesByClass(iftDataSet *Z)
{
    iftMatrix* confusion = iftComputeConfusionMatrix(Z, false);
    iftIntArray *nclass = iftCreateIntArray(Z->nclasses + 1);
    iftFloatArray *acc = iftCreateFloatArray(Z->nclasses + 1);

    for (int i = 0; i < Z->nsamples; i++)
        if (!iftHasSampleStatus(Z->sample[i], IFT_TRAIN))
            nclass->val[Z->sample[i].truelabel]++;

    for (int c = 1; c < Z->nclasses+1; c++)
        acc->val[c] = iftMatrixElem(confusion, c-1, c-1) / (float)nclass->val[c];

    return acc;
}


float iftMeasureERR(float *negatives, int total_neg, float *positives, int total_pos, float threshold) 
{
  int FRcount = 0; /* false rejection count */
  int FAcount = 0; /* false acceptance count */
  float FRR, FAR;

    // Compute the FAR (False Acceptance Rate)
    for (int i = 0; i < total_neg; i++)
        if (negatives[i] >= threshold)
            FAcount++;
    if (total_neg > 0)
        FAR = (float) FAcount/total_neg;
    else
        FAR = 0;

    // Compute the FRR (False Rejection Rate)
    for (int i = 0; i < total_pos; i++)
        if (positives[i] < threshold)
            FRcount++;
    if (total_pos > 0)
        FRR = (float) FRcount/total_pos;
    else
        FRR = 0;

    printf("- FAR: %f, FRR: %f\n", FAR, FRR);
    float ERR = (float) (FAR+FRR)/2;

    return ERR;
}


void iftMeasureFarFrr(iftDataSet *Z, float *FAR, float *FRR) {
    int FAcount = 0;
    int FRcount = 0;
    int total_pos = 0;
    int total_neg = 0;

    for (int s = 0; s < Z->nsamples; s++)
        if (Z->sample[s].truelabel == 1) {
            total_pos++;
            if (Z->sample[s].truelabel != Z->sample[s].label)
                FRcount++;
        } else if (Z->sample[s].truelabel == 2) {
            total_neg++;
            if (Z->sample[s].truelabel != Z->sample[s].label)
                FAcount++;
        }

    *FAR = (float) FAcount / total_neg;
    *FRR = (float) FRcount / total_pos;

    printf("\n***** FRcount = %d, total_pos: %d, FRR: %f\n", FRcount, total_pos, *FRR);
    printf("***** FAcount = %d, total_neg: %d, FAR: %f\n\n", FAcount, total_neg, *FAR);
}


/*
 * Compute the Confusion Matrix for the Ztest.
 * IMPORTANT: The dataset must be Binary with class 1 = positive class, class 2 = negative class
 */
iftErrorClassification iftMeasureErrorClassification(iftDataSet *Ztest) {
    if (Ztest->nclasses != 2)
        iftError("It is not a Binary Dataset", "iftComputeConfusionMatrix");

    int num_errors[2] = {0, 0}; // [0] = class 1 (positive); [1] = class 2 (negative)
    int num_hits[2]   = {0, 0}; // [0] = class 1 (positive); [1] = class 2 (negative)
    int POS_IDX = 0, NEG_IDX = 1;

    for (int s = 0; s < Ztest->nsamples; s++) {
        if (Ztest->sample[s].truelabel == Ztest->sample[s].label) { // hit
            num_hits[Ztest->sample[s].truelabel - 1]++;
        }
        else { // error
            num_errors[Ztest->sample[s].truelabel - 1]++;
        }
    }

    iftErrorClassification cm;
    cm.tp = num_hits[POS_IDX];
    cm.tn = num_hits[NEG_IDX];
    cm.fn = num_errors[POS_IDX];
    cm.fp = num_errors[NEG_IDX];

    printf("TP: %d/%d, FN: %d/%d", cm.tp, (cm.tp+cm.fn), cm.fn, (cm.tp+cm.fn));
    printf("TN: %d/%d, FP: %d/%d", cm.tn, (cm.tn+cm.fp), cm.fp, (cm.tn+cm.fp));

    return cm;
}


/*
 * Variation of the accuracy founded in Papa, J. P. - A discrete approach for supervised pattern recognition.
 */
float iftWeightedAccuracyBinaryClasses(iftErrorClassification cm) {
    int num_neg  = cm.tn + cm.fp;
    int num_pos  = cm.tp + cm.fn;

    float e1 = (float) cm.fp/num_neg;
    float e2 = (float) cm.fn/num_pos;

    float acc = 1 - (e1 + e2)/2;

    puts("**********************************");
    printf("- TP: %d, FN: %d, FP: %d, TN: %d\n", cm.tp, cm.fn, cm.fp, cm.tn);
    printf("- e1: %f\n", e1);
    printf("- e2: %f\n", e2);
    printf("- acc: %f\n", acc);
    puts("**********************************\n");

    return acc;
}

/*
 * Compute the Matthews Correlation Coefficient.
 * It returns a value between −1 and +1. A coefficient of +1 represents a perfect prediction,
 * 0 no better than random prediction and −1 indicates total disagreement between prediction and observation.
 */
float iftMeasureMCC(iftErrorClassification cm) {
    int TP = cm.tp, FN = cm.fn, FP = cm.fp, TN = cm.tn;
//  printf("* TP: %d, TN: %d, FP: %d, FN: %d\n", TP, TN, FP, FN);

    float numerator = (float) (TP*TN) - (FP*FN);
    float denominator = (float) (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN);
    float acc = (float) numerator/sqrt(denominator);

//  puts("**********************************");
//  printf("- mcc = %.0f/%f = %f\n", numerator, denominator, acc);
//  puts("**********************************\n");

    return acc;
}


/**
 * @brief Distance of Cosines.
 *
 * This function computes the distance of cosines of two float vectors.
 * It assumes that the feature vectors have norm one. Use iftNormOneDataSet.\n
 *
 * @param f1 Float Vector 1.
 * @param f2 Float Vector 2.
 * @param n Size of the vectors
 * @return Int value from 0 (totally similar) to 1 (totally dissimilar).
*/
float iftCosineDistance(float *f1, float *f2, int n) {
    float inner_prod = 0;

    for (int i = 0; i < n; i++)
        inner_prod += f1[i]*f2[i];

    return(0.5 - inner_prod/2.0);
}


/**
 * @brief Distance of Cosines.
 *
 * This function computes the distance of cosines of two float vectors.
 *
 * @param f1 Float Vector 1.
 * @param f2 Float Vector 2.
 * @param n Size of the vectors
 * @return Int value from 0 (totally similar) to 1 (totally dissimilar).
*/
float iftCosineDistance2(float *f1, float *f2, int n) {
    float v1_norm = 0.0;
    float v2_norm = 0.0;
    float v1_dot_v2 = 0.0;

    for (int i = 0; i < n; i++) {
        v1_dot_v2 += f1[i]*f2[i];
        v1_norm   += f1[i]*f1[i];
        v2_norm   += f2[i]*f2[i];
    }
    v1_norm = sqrtf(v1_norm);
    v2_norm = sqrtf(v2_norm);

    if (v1_norm == 0)
        iftError("V1 norm = 0", "iftCosineDistance2");
    if (v2_norm == 0)
        iftError("V2 norm = 0", "iftCosineDistance2");

    float cos_sim = (float) v1_dot_v2 / (v1_norm*v2_norm);

    return (float) 0.5 - cos_sim/2;
}


/**
 * @brief Compute the Correlation Coefficient of two vectors.
 *
 * This function computes the Correlation Coefficient value -1 <= p(x, y) <= 1.
 * If x and y are linearly dependent, p(x, y) = -1 or 1. If p(x, y) = 0, x and y are totally uncorrelated.
 * Reference: http://www.math.uiuc.edu/~hildebr/461/variance.pdf
 *
 * @param x Float Vector 1.
 * @param y Float Vector 2.
 * @param n Size of the vectors
 * @return Correlation Coefficient
*/
float iftCorrCoef(float *x, float *y, int n) {
    float cov = iftCov(x, y, n);
    float var_x = iftVar(x, n);
    float var_y = iftVar(y, n);

    return (float) cov/sqrtf(var_x*var_y);
}

/**
 * @brief Compute the Distance based on the Correlation Coefficient of two vectors.
 *
 * This function computes the Distance based on Correlation Coefficient value 0 <= dist(x, y) <= 1.
 * If x and y are linearly dependent, dist(x, y) = 0.
 * Reference: http://dl.acm.org/citation.cfm?id=507477
 *
 * @param x Float Vector 1.
 * @param y Float Vector 2.
 * @param n Size of the vectors
 * @return Distance based on the Correlation Coefficient
*/
float iftCorrCoefDist(float *x, float *y, int n) {
    return (float) (1 - fabs(iftCorrCoef(x, y, n)));
}


/**
 * @brief Compute the Mutual Information Compression Index of two vectors.
 *
 * This function computes the Mutual Information Compression Index value. 0 <= mici(x,y) <= 0.5*(var(x) + var(y)).
 * If mici(x,y) = 0, x and y are linearly related.
 * Reference: http://dl.acm.org/citation.cfm?id=507477
 *
 * @param x Float Vector 1.
 * @param y Float Vector 2.
 * @param n Size of the vectors
 * @return Mutual Information Compression Index
 *
*/
float iftMICI(float *x, float *y, int n) {
    float var_x = iftVar(x, n);
    float var_y = iftVar(y, n);
    float corr_coef = iftCorrCoef(x, y, n);

    float mici = var_x + var_y - sqrtf(((var_x+var_y)*(var_x+var_y)) - 4*var_x*var_y*(1 - (corr_coef*corr_coef)));

    return mici;
}



int iftBestSampleScoreSum(const iftMatrix *score_matrix, iftMetricType metric_type, double *out_best_score) {
    int most_similar  = -1; // idx from the most similar sample
    double best_score = 0.0;

    // the code was replicated due to optimization issues
    if (metric_type == IFT_SIMILARITY_FUNCTION) {
        best_score = 0.0;

        // finds the most similar image from the Score Matrix
        for (int i = 0; i < score_matrix->nrows; i++) {
            double aux = 0.0;

            #pragma omp parallel for reduction(+: aux)
            for (int j = 0; j < score_matrix->ncols; j++)
                aux += score_matrix->val[iftGetMatrixIndex(score_matrix, j, i)];

            printf("aux = %f\n", aux);
            printf("best_score = %f\n", best_score);

            if (aux > best_score) {
                most_similar = i;
                best_score   = aux;
            }
        }
    }
    else if (metric_type == IFT_DISTANCE_FUNCTION) {
        best_score = IFT_INFINITY_DBL;

        // finds the most similar image from the Score Matrix
        for (int i = 0; i < score_matrix->nrows; i++) {
            double aux = 0.0;

            #pragma omp parallel for reduction(+: aux)
            for (int j = 0; j < score_matrix->ncols; j++)
                aux += score_matrix->val[iftGetMatrixIndex(score_matrix, j, i)];

            printf("aux = %f\n", aux);
            printf("best_score = %f\n", best_score);

            if (aux < best_score) {
                most_similar = i;
                best_score   = aux;
            }
        }
    }
    else
        iftError("Invalid Image Similarity Function Type... Try: IFT_SIMILARITY_FUNCTION or IFT_DISTANCE_FUNCTION",
                 "iftBestSampleScoreSum");

    printf("* best_score = %f\n", best_score);

    
    if (out_best_score != NULL)
        *out_best_score = best_score;

    return most_similar;
}




float iftKullbackLeiblerDivergence(const iftFloatArray *P, const iftFloatArray *Q) {
    if (P->n != Q->n)
        iftError("Distributions with different sizes: %ld != %ld", "iftKullbackLeiblerDivergence",
                  P->n, Q->n);

    // normalizing the distribuitions
    float sum_P = iftSumFloatArrayElems(P->val, P->n);
    iftFloatArray *Pnorm = iftCreateFloatArray(P->n);
    iftScaleFloatArray(P->val, 1.0 / sum_P, Pnorm->val, P->n);

    float sum_Q = iftSumFloatArrayElems(Q->val, Q->n);
    iftFloatArray *Qnorm = iftCreateFloatArray(Q->n);
    iftScaleFloatArray(Q->val, 1.0 / sum_Q, Qnorm->val, Q->n);

    float KL = 0.0;

    #pragma omp parallel for reduction(+:KL)
    for (int i = 0; i < Pnorm->n; i++)
        if ((Pnorm->val[i] > 0))
            KL += (Pnorm->val[i] * (log(Pnorm->val[i]) - log(Qnorm->val[i])));

    iftDestroyFloatArray(&Pnorm);
    iftDestroyFloatArray(&Qnorm);

    return KL;
}

double iftMImageMahalanobis(const iftMImage *mimg, const iftMatrix *M, int p, int q)
{
    double *diff = iftAlloc(mimg->m, sizeof *diff);
    double *aux = iftAlloc(mimg->m, sizeof *aux);

    for (int i = 0; i < mimg->m; i++) {
        aux[i] = 0;
        diff[i] = mimg->val[p][i] - mimg->val[q][i];
    }

    for (int i = 0; i < mimg->m; i++) {
        for (int j = 0; j < mimg->m; j++) {
            aux[i] += diff[j] * M->val[j * mimg->m + i];
        }
    }

    double out = 0;
    for (int i = 0; i < mimg->m; i++) {
        out += aux[i] * diff[i];
    }

    iftFree(diff);
    iftFree(aux);

    return out;
}

double iftMImageMahalanobisDouble(const iftMImage *mimg, const double *M, int p, int q)
{
    double *diff = iftAlloc(mimg->m, sizeof *diff);
    double *aux = iftAlloc(mimg->m, sizeof *aux);

    for (int i = 0; i < mimg->m; i++) {
        aux[i] = 0;
        diff[i] = mimg->val[p][i] - mimg->val[q][i];
    }

    for (int i = 0; i < mimg->m; i++) {
        for (int j = 0; j < mimg->m; j++) {
            aux[i] += diff[j] * M[j * mimg->m + i];
        }
    }

    double out = 0;
    for (int i = 0; i < mimg->m; i++) {
        out += aux[i] * diff[i];
    }

    iftFree(diff);
    iftFree(aux);

    return out;
}

