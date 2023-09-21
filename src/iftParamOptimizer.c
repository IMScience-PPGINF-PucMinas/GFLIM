//
// Created by peixinho on 9/25/15.
//

#include "iftParamOptimizer.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Stream.h"


int iftParamsSize(iftDict* opt, int* paramsSize) {
    int i=0;
    iftIntArray* intArray = NULL;
    iftDblArray* dblArray = NULL;
    iftStrArray* strArray = NULL;
    iftDict* dict = NULL;

    for (iftKeyVal k = iftDictFirst(opt); iftDictValidKey(k); k = iftDictNext(opt, k)) {

        switch (k.val.type) {
            case IFT_INT_ARRAY_TYPE:
                intArray = iftGetIntArrayVal(k.val);
                paramsSize[i] = intArray->n;
                break;
            case IFT_DBL_ARRAY_TYPE:
                dblArray = iftGetDblArrayVal(k.val);
                paramsSize[i] = dblArray->n;
                break;
            case IFT_STR_ARRAY_TYPE:
                strArray = iftGetStrArrayVal(k.val);
                paramsSize[i] = strArray->n;
                break;
            case IFT_DICT_TYPE:
                dict = iftGetDictVal(k.val);
                i += iftParamsSize(dict, paramsSize+i) - 1;
                break;
            default:
                iftError("%s dictionary parameter is not a string/double/integer array", "iftParamsSize",
                         iftGValToString(k.key));
                break;
        }
        ++i;
    }

    return i;
}

int iftParamsNumberDict(iftDict* opt) {
    int nparams = 0;

    for(iftKeyVal k = iftDictFirst(opt); iftDictValidKey(k); k = iftDictNext(opt, k)) {
        if(k.val.type == IFT_DICT_TYPE) {
            iftDict* dict = iftGetDictVal(k.val);
            nparams += iftParamsNumberDict(dict);
        }
        else {
            nparams++;
        }
    }
    return nparams;
}

int iftParamsDict(iftDict* opt, int* idx, iftDict* params) {
    int i=0;
    iftIntArray* intArray = NULL;
    iftDblArray* dblArray = NULL;
    iftStrArray* strArray = NULL;
    iftDict* dict = NULL;

    for (iftKeyVal k = iftDictFirst(opt); iftDictValidKey(k); k = iftDictNext(opt, k)) {

        switch (k.val.type) {
            case IFT_INT_ARRAY_TYPE:
                intArray = iftGetIntArrayVal(k.val);
                iftInsertKeyValIntoDict(iftCopyGVal(k.key), iftCreateGVal(intArray->val[idx[i]]), params);
                break;
            case IFT_DBL_ARRAY_TYPE:
                dblArray = iftGetDblArrayVal(k.val);
                iftInsertKeyValIntoDict(iftCopyGVal(k.key), iftCreateGVal(dblArray->val[idx[i]]), params);
                break;
            case IFT_STR_ARRAY_TYPE:
                strArray = iftGetStrArrayVal(k.val);
                iftInsertKeyValIntoDict(iftCopyGVal(k.key), iftCreateGVal(strArray->val[idx[i]]), params);
                break;
            case IFT_DICT_TYPE:
                dict = iftCreateDict();
                iftInsertKeyValIntoDict(iftCopyGVal(k.key), iftCreateGVal(dict), params);
                i += iftParamsDict(iftGetDictVal(k.val), &idx[i], dict) - 1;
                break;
            default:
                iftError("%s dictionary parameter is not a string/double/integer array", "iftGridSearch",
                         iftGValToString(k.key));
                break;
        }
        ++i;
    }
    return i;
}

iftDict* iftGridSearch(iftDict *opt, iftDictFitnessFunc f, iftDict *problem) {

    int nparams = iftParamsNumberDict(opt);
    int* idx = iftAllocIntArray(nparams);
    iftDict *params = iftCreateDict();

    double max = IFT_INFINITY_FLT_NEG;
    int *maxIdx = iftAllocIntArray(nparams);
    double meanEvalTime = 0.0;

    int* paramsSize = iftAllocIntArray(nparams);

    iftParamsSize(opt, paramsSize);
    ////iftPrintIntArray(paramsSize, nparams);

    int count = 0;
    for(;;) {
        timer *t1 = iftTic();
        iftParamsDict(opt, idx, params);

        ////iftPrintIntArray(idx, nparams);

        double eval = f(problem, params);

        if(eval>max) {
            max = eval;
            iftCopyIntArray(maxIdx, idx, nparams);
        }

        //increment index
        for (int i = nparams - 1; i >= 0; --i) {
            idx[i]++;
            if (idx[i] >= paramsSize[i])
                idx[i] = 0;
            else break;
        }
        meanEvalTime += iftCompTime(t1, iftToc()) / 1000.0;
        //iftFree(t1);
        count++;

        //has gone through all the values
        if(iftSumIntArray(idx, nparams) <= 0) break;
    }
    meanEvalTime /= (double)count;

    iftParamsDict(opt, maxIdx, params);
    /*insert the best function objetive value*/
    iftInsertIntoDict("best_func_val", max, params);
    iftInsertIntoDict("mean_eval_time", meanEvalTime, params);

    iftFree(idx);
    iftFree(maxIdx);
    iftFree(paramsSize);

    return params;
}


/**
 * @brief An adapter to transform a iftClassifierFitnessFunc into a iftFitnessFunc to be used in the optimization methods.
 */
double iftClassifierFitnessFuncAdapter(iftDict *problem, iftDict *params) {

    iftDataSet* dataset = (iftDataSet*) iftGetPtrValFromDict("dataset", problem);
    iftSampler* sampler = (iftSampler*) iftGetPtrValFromDict("sampler", problem);
    iftClassifierFitnessFunc func = (iftClassifierFitnessFunc) (size_t) iftGetPtrValFromDict("fitness", problem);
    iftDict* p = iftGetDictFromDict("problem", problem);

    double acc = 0.0;
    for (int it = 0; it < sampler->niters; ++it) {
		iftDataSetSampling(dataset, sampler, it);
        acc += func(p, dataset, params) / sampler->niters;
    }

    return acc;
}

/**
 * @brief An adapter to transform a iftClassifierFitnessFunc into a iftFitnessFunc to be used in the optimization methods.
 */
double iftDescriptorFitnessFuncAdapter(iftDict *problem, iftDict *params) {

    iftDescriptorFitnessFunc f = (iftDescriptorFitnessFunc)(size_t)iftGetPtrValFromDict("fitness", problem);
    iftDict* p = iftGetDictFromDict("problem", problem);
    iftFileSet* fileSet = (iftFileSet*)iftGetPtrValFromDict("fileset", problem);

    double acc = f(p, fileSet, params);

    return acc;
}

iftDict * iftGridSearchClassifier(iftDict *opt, iftClassifierFitnessFunc f, iftDataSet *Z, iftSampler *sampler,
                                  iftDict *problem) {

    iftDict* classifierProblem = iftCreateDict();
    iftInsertIntoDict("sampler", (void*)sampler, classifierProblem);
    iftInsertIntoDict("problem", problem, classifierProblem);
    iftInsertIntoDict("fitness", (void*)(size_t)f, classifierProblem);//ansi forbids function pointer to void* conversion
    iftInsertIntoDict("dataset", (void*)Z, classifierProblem);

    iftDict* best = iftGridSearch(opt, iftClassifierFitnessFuncAdapter, classifierProblem);

    iftDestroyDict(&classifierProblem);

    return best;
}

iftDict * iftGridSearchDescriptor(iftDict *opt, iftDescriptorFitnessFunc f, iftFileSet *fileSet, iftDict *problem) {

    iftDict* descriptorProblem = iftCreateDict();

    iftInsertIntoDict("fileset", (void*)fileSet, descriptorProblem);
    iftInsertIntoDict("problem", problem, descriptorProblem);
    iftInsertIntoDict("fitness", (void*)(size_t)f, descriptorProblem);//ansi forbids function pointer to void* conversion

    iftDict* best = iftGridSearch(opt, iftDescriptorFitnessFuncAdapter, descriptorProblem);

    iftDestroyDict(&descriptorProblem);

    return best;
}

iftDict* iftRandomSearch(iftDict *opt, iftDictFitnessFunc f, int ntrials, iftDict *problem) {

    int nparams = iftParamsNumberDict(opt);
    int* idx = iftAllocIntArray(nparams);
    int* maxIdx = iftAllocIntArray(nparams);
    iftDict* params = iftCreateDict();
    double max = IFT_INFINITY_DBL_NEG;
    double meanEvalTime = 0.0;

    int* paramsSize = iftAllocIntArray(nparams);
    iftParamsSize(opt, paramsSize);

    ////iftPrintIntArray(paramsSize, nparams);
    for (int trial = 0; trial < ntrials; ++trial) {
        timer *t1 = iftTic();
        for (int i = 0; i < nparams; ++i) {
            idx[i] = iftRandomInteger(0, paramsSize[i] - 1);
        }

        iftParamsDict(opt, idx, params);

        double eval = f(problem, params);

        if(eval>max) {
            max = eval;
            iftCopyIntArray(maxIdx, idx, nparams);
        }
        meanEvalTime += iftCompTime(t1, iftToc()) / 1000.0;
        iftFree(t1);
    }
    meanEvalTime /= (double)ntrials;

    iftParamsDict(opt, maxIdx, params);

    /*insert the best function objetive value and the mean execution time */
    iftInsertIntoDict("best_func_val", max, params);
    iftInsertIntoDict("mean_eval_time", meanEvalTime, params);

    iftFree(idx);
    iftFree(maxIdx);

    return params;
}

iftDict * iftRandomSearchClassifier(iftDict *optimizer, iftClassifierFitnessFunc func, iftDataSet *dataset,
                                    iftSampler *sampler, int ntrials, iftDict *problem) {

    iftDict* classifierProblem = iftCreateDict();
    iftInsertIntoDict("sampler", (void*)sampler, classifierProblem);
    iftInsertIntoDict("problem", problem, classifierProblem);
    iftInsertIntoDict("fitness", (void*)(size_t)func, classifierProblem);//ansi forbids function pointer to void* conversion
    iftInsertIntoDict("dataset", (void*)dataset, classifierProblem);

    iftDict* best = iftRandomSearch(optimizer, iftClassifierFitnessFuncAdapter, ntrials, classifierProblem);

    iftDestroyDict(&classifierProblem);

    return best;
}

iftDict * iftRandomSearchDescriptor(iftDict *opt, iftDescriptorFitnessFunc f, iftFileSet *fileSet, int ntrials,
                                    iftDict *problem) {

    iftDict* descriptorProblem = iftCreateDict();

    iftInsertIntoDict("fileset", (void*)fileSet, descriptorProblem);
    iftInsertIntoDict("problem", problem, descriptorProblem);
    iftInsertIntoDict("fitness", (void*)(size_t)f, descriptorProblem);//ansi forbids function pointer to void* conversion

    iftDict* best = iftRandomSearch(opt, iftDescriptorFitnessFuncAdapter, ntrials, descriptorProblem);

    iftDestroyDict(&descriptorProblem);

    return best;
}

iftIntArray *iftUniformIntSearchSpace(int begin, int end, int step)
{
    if(begin > end)
        iftError("Begin must be lower than end", "iftUniformIntSearchSpace");

    if(step <= 0)
        iftError("Step must be positive", "iftUniformIntSearchSpace");

    int val; int i;
    iftIntArray *array = iftCreateIntArray(floor((float)(end-begin)/(float)step) + 1);
    for(val = begin, i = 0; val <= end; val += step, i++) {
        array->val[i] = val;
    }
    return array;
}

iftDblArray *iftUniformDblSearchSpace(double begin, double end, double step)
{
    if(begin > end)
        iftError("Begin must be lower than end", "iftUniformDblSearchSpace");

    if(step <= 0)
        iftError("Step must be positive", "iftUniformDblSearchSpace");

    double val; int i;
    iftDblArray *array = iftCreateDblArray(floor((end-begin)/step) + 1);
    for(val = begin, i = 0; val <= end; val += step, i++) {
        array->val[i] = val;
    }
    return array;
}
