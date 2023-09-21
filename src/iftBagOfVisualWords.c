//
// Created by Cesar Castelo on Jan 04, 2018
//
#include "iftBagOfVisualWords.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/IntMatrix.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/FileSet.h"
#include "ift/core/io/File.h"
#include "ift/core/io/Json.h"
#include "iftDeepLearning.h"


/* ---------------------------------------------------------------------------------------------------------------------*/
/*                                                Initialization functions                                              */
/* ---------------------------------------------------------------------------------------------------------------------*/

iftBagOfVisualWords *iftCreateBovw()
{

    iftBagOfVisualWords* bovw = (iftBagOfVisualWords *) iftAlloc(1, sizeof(iftBagOfVisualWords));

    bovw->localFeats = NULL;
    bovw->dict = NULL;

    bovw->localFeatsOrderHrchy = NULL;
    bovw->orderHrchies = NULL;

    bovw->dictLearnFileset = NULL;
    bovw->dictLearnMasksFileset = NULL;

    bovw->intPointDetector = NULL;
    bovw->intPointDetectorId = 0;

    bovw->localFeatExtractor = NULL;
    bovw->localFeatExtractorId = 0;

    bovw->dictEstimator = NULL;
    bovw->dictEstimatorId = 0;

    bovw->codFunc = NULL;
    bovw->codFuncId = 0;

    bovw->arcWeightFunc = NULL;
    bovw->distFuncId = 0;

    bovw->funcParams = NULL;

    return bovw;
}

void iftBovwSetFunctionPointers(iftBagOfVisualWords *bovw, iftBovwIntPointDetectorId intPointDetector,
                                iftBovwLocalFeatExtractorId localFeatExtractor, iftBovwDictEstimatorId dictEstimator,
                                iftBovwCodFuncId codFunc, iftBovwDistFuncId distFunc, iftDict *funcParams)
{
    /* set interest point detector function pointer */
    bovw->intPointDetectorId = intPointDetector;
    switch(intPointDetector) {
        case BOVW_INT_POINT_DETECTOR_RANDOM:
            bovw->intPointDetector = iftBovwRandomIntPointDetector; break;
        case BOVW_INT_POINT_DETECTOR_GRID:
            bovw->intPointDetector = iftBovwGridIntPointDetector; break;
        case BOVW_INT_POINT_DETECTOR_UNSUP_SPIX_ISF:
            bovw->intPointDetector = iftBovwUnsupSuperpixelISFIntPointDetector; break;
        case BOVW_INT_POINT_DETECTOR_SUP_SPIX_ISF:
            bovw->intPointDetector = iftBovwSupSuperpixelISFIntPointDetector; break;
        default:
            bovw->intPointDetector = NULL;
    }

    /* set feature extractor function pointer */
    bovw->localFeatExtractorId = localFeatExtractor;
    switch(localFeatExtractor) {
        case BOVW_LOCAL_FEAT_EXTRACTOR_RAW:
            bovw->localFeatExtractor = iftBovwRawFeatExtractor; break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_BIC:
            bovw->localFeatExtractor = iftBovwBICFeatExtractor; break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_LBP:
            bovw->localFeatExtractor = iftBovwLBPFeatExtractor; break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_BRIEF:
            bovw->localFeatExtractor = iftBovwBriefFeatExtractor; break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_CONV:
            bovw->localFeatExtractor = iftBovwConvolutionalFeatExtractor; break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_CONV_MULTI_LAYER:
            bovw->localFeatExtractor = iftBovwConvolutionalMultiLayerFeatExtractor; break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_DEEP_FEATS_MIMG:
            bovw->localFeatExtractor = iftBovwDeepFeaturesMImageFeatExtractor; break;
        default:
            bovw->localFeatExtractor = NULL;
    }

    /* determine the number of order hierarchies */
    int nOrderHrchies;
    switch(dictEstimator) {
        case BOVW_DICT_ESTIMATOR_UNSUP_KMEANS:
        case BOVW_DICT_ESTIMATOR_UNSUP_OPF:
            nOrderHrchies = 0; break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_PIX_LABEL:
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMG_CLASS:
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMAGE:
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_POSITION:
        case BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL:
        case BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS:
        case BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE:
        case BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION:
        case BOVW_DICT_ESTIMATOR_SUP_MANUAL_IMG_CLASS:
        case BOVW_DICT_ESTIMATOR_SUP_MANUAL_POSITION:
            nOrderHrchies = 1; break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMG_CLASS_AND_POSITION:
        case BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS_AND_POSITION:
        case BOVW_DICT_ESTIMATOR_SUP_MANUAL_IMG_CLASS_AND_POSITION:
            nOrderHrchies = 2; break;
        default:
            nOrderHrchies = 0;
    }
    if(nOrderHrchies > 0) {
        bovw->orderHrchies = iftCreateIntArray(nOrderHrchies);
    }

    /* set dictionary estimator function pointer */
    bovw->dictEstimatorId = dictEstimator;
    switch(dictEstimator) {
        case BOVW_DICT_ESTIMATOR_UNSUP_KMEANS:
            bovw->dictEstimator = iftBovwUnsupKMeansDictEstimator;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_PIX_LABEL:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyKMeansDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_PIX_LABEL;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMG_CLASS:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyKMeansDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMAGE:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyKMeansDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_IMAGE;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_POSITION:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyKMeansDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_POSITION;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMG_CLASS_AND_POSITION:
            bovw->dictEstimator = iftBovwSupTwoLevelOrderHrchyKMeansDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS;
            bovw->orderHrchies->val[1] = BOVW_LOCAL_FEATS_ORDERING_POSITION;
            break;
        case BOVW_DICT_ESTIMATOR_UNSUP_OPF:
            bovw->dictEstimator = iftBovwUnsupOPFDictEstimator;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyOPFDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_PIX_LABEL;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyOPFDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyOPFDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_IMAGE;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyOPFDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_POSITION;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS_AND_POSITION:
            bovw->dictEstimator = iftBovwSupTwoLevelOrderHrchyOPFDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS;
            bovw->orderHrchies->val[1] = BOVW_LOCAL_FEATS_ORDERING_POSITION;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_MANUAL_IMG_CLASS:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyManualDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_MANUAL_POSITION:
            bovw->dictEstimator = iftBovwSupOneLevelOrderHrchyManualDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_POSITION;
            break;
        case BOVW_DICT_ESTIMATOR_SUP_MANUAL_IMG_CLASS_AND_POSITION:
            bovw->dictEstimator = iftBovwSupTwoLevelOrderHrchyManualDictEstimator;
            bovw->orderHrchies->val[0] = BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS;
            bovw->orderHrchies->val[1] = BOVW_LOCAL_FEATS_ORDERING_POSITION;
            break;
        default:
            bovw->dictEstimator = NULL;
    }

    /* set coding function pointer */
    bovw->codFuncId = codFunc;
    switch(codFunc) {
        case BOVW_COD_FUNC_HARD_ASGMT:
            bovw->codFunc = iftBovwHardAsgmtCodFunc; break;
        case BOVW_COD_FUNC_SOFT_ASGMT:
            bovw->codFunc = iftBovwSoftAsgmtCodFunc; break;
        case BOVW_COD_FUNC_HARD_ASGMT_BATCH:
            bovw->codFunc = iftBovwHardAsgmtBatchCodFunc; break;
        case BOVW_COD_FUNC_SOFT_ASGMT_BATCH:
            bovw->codFunc = iftBovwSoftAsgmtBatchCodFunc; break;
        case BOVW_COD_FUNC_2ND_ORDER_STATS_FISHER_VECTORS:
            bovw->codFunc = iftBovw2ndOrderStatsFisherVectorsCodFunc; break;
        default:
            bovw->codFunc = NULL;
    }

    /* set distance function pointer */
    bovw->distFuncId = distFunc;
    switch(distFunc) {
        case BOVW_DIST_FUNC_1:
            bovw->arcWeightFunc = iftDistance1; break;
        case BOVW_DIST_FUNC_2:
            bovw->arcWeightFunc = iftDistance2; break;
        case BOVW_DIST_FUNC_3:
            bovw->arcWeightFunc = iftDistance3; break;
        case BOVW_DIST_FUNC_4:
            bovw->arcWeightFunc = iftDistance4; break;
        case BOVW_DIST_FUNC_5:
            bovw->arcWeightFunc = iftDistance5; break;
        case BOVW_DIST_FUNC_6:
            bovw->arcWeightFunc = iftDistance6; break;
        case BOVW_DIST_FUNC_7:
            bovw->arcWeightFunc = iftDistance7; break;
        case BOVW_DIST_FUNC_8:
            bovw->arcWeightFunc = iftDistance8; break;
        case BOVW_DIST_FUNC_9:
            bovw->arcWeightFunc = iftDistance9; break;
        case BOVW_DIST_FUNC_10:
            bovw->arcWeightFunc = iftDistance10; break;
        case BOVW_DIST_FUNC_11:
            bovw->arcWeightFunc = iftDistance11; break;
        case BOVW_DIST_FUNC_12:
            bovw->arcWeightFunc = iftDistance12; break;
        default:
            bovw->arcWeightFunc = iftDistance1;
    }

    /* set the function pointers' params */
    bovw->funcParams = iftCopyDict(funcParams);
}

void iftDestroyBovw(iftBagOfVisualWords **bovw)
{
    if((*bovw)->dict != NULL) {
        if(iftBovwIsFisherVectorsEncoding((*bovw)->codFuncId)) {
            if((*bovw)->dictEstimatorId == BOVW_DICT_ESTIMATOR_UNSUP_OPF) {
                iftBovwOPFClustModelFV *clustModel = (iftBovwOPFClustModelFV*)((*bovw)->dict);
                iftDestroyBovwOPFClustModelFV(&clustModel);
            }

            if( (*bovw)->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL ||
                (*bovw)->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS ||
                (*bovw)->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE ||
                (*bovw)->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION) {
                
                iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)((*bovw)->dict);
                int nClasses = iftCountUniqueIntElems(iftMatrixRowPointer((*bovw)->localFeatsOrderHrchy, 0), (*bovw)->localFeats->nsamples);
                for(int c = 1; c < nClasses+1; c++)
                    iftDestroyBovwOPFClustModelFV(&(clustModels[c]));
                iftFree(clustModels);
            }
        }
        else {
            iftDataSet *dict = (iftDataSet*)(*bovw)->dict;
            iftDestroyDataSet(&dict);
        }
    }

    if((*bovw)->dictLearnFileset != NULL)
        iftDestroyFileSet(&((*bovw)->dictLearnFileset));

    if((*bovw)->dictLearnMasksFileset != NULL)
        iftDestroyFileSet(&((*bovw)->dictLearnMasksFileset));

    if((*bovw)->localFeats != NULL)
        iftDestroyDataSet(&((*bovw)->localFeats));

    if((*bovw)->localFeatsOrderHrchy != NULL)
        iftDestroyIntMatrix(&((*bovw)->localFeatsOrderHrchy));

    if((*bovw)->orderHrchies != NULL)
        iftDestroyIntArray(&((*bovw)->orderHrchies));

    if((*bovw)->funcParams != NULL)
        iftDestroyDict(&((*bovw)->funcParams));

    iftFree(*bovw);
}

void iftBovwWrite(iftBagOfVisualWords *bovw, char *filename)
{
    if(!iftCompareStrings(iftFileExt(filename), ".bovw"))
        iftError("Invalid file extension", "iftBovwWrite");

    char *tmpDir = iftMakeTempDir("tmpdir_", "bovw", "./");
    char auxFilename[2048];

    /* create JSON with dict params */
    iftDict *dictParams = iftCreateDict();
    iftInsertIntoDict("num_local_feats", bovw->localFeats->nsamples, dictParams);
    iftInsertIntoDict("dim_local_feats", bovw->localFeats->nfeats, dictParams);
    iftInsertIntoDict("int_point_detector_id", bovw->intPointDetectorId, dictParams);
    iftInsertIntoDict("local_feat_extractor_id", bovw->localFeatExtractorId, dictParams);
    iftInsertIntoDict("dict_estimator_id", bovw->dictEstimatorId, dictParams);
    iftInsertIntoDict("cod_func_id", bovw->codFuncId, dictParams);
    iftInsertIntoDict("dist_func_id", bovw->distFuncId, dictParams);
    iftInsertIntoDict("func_params", bovw->funcParams, dictParams);
    
    sprintf(auxFilename, "%s/dict_params.json", tmpDir);
    iftWriteJson(dictParams, auxFilename);

    /* dictionary files */
    if(iftBovwIsFisherVectorsEncoding(bovw->codFuncId)) {
        if(bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_UNSUP_OPF) {
            iftBovwOPFClustModelFV *clustModel = (iftBovwOPFClustModelFV*)(bovw->dict);
            sprintf(auxFilename, "%s/dict.model", tmpDir);
            iftWriteBovwOPFClustModelFV(clustModel, auxFilename);
        }
        if( bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION) {
            
            iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)(bovw->dict);
            int nClasses = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 0), bovw->localFeats->nsamples);
            for(int c = 1; c < nClasses+1; c++) {
                sprintf(auxFilename, "%s/dict_%d.model", tmpDir, c);
                iftWriteBovwOPFClustModelFV(clustModels[c], auxFilename);
            }
        }

    }
    else {
        iftDataSet *dict = (iftDataSet*)bovw->dict;
        sprintf(auxFilename, "%s/dict.zip", tmpDir);
        iftWriteDataSet(dict, auxFilename);
    }

    /* write individual files */
    if(bovw->localFeatsOrderHrchy) {
        sprintf(auxFilename, "%s/local_feats_order_hrchy.npy", tmpDir);
        iftWriteIntMatrix(bovw->localFeatsOrderHrchy, auxFilename);
    }

    if(bovw->orderHrchies) {
        sprintf(auxFilename, "%s/order_hrchies.npy", tmpDir);
        iftWriteIntArray(bovw->orderHrchies, auxFilename);
    }

    sprintf(auxFilename, "%s/dict_learn_fileset.csv", tmpDir);
    iftWriteFileSetAsCSV(bovw->dictLearnFileset, auxFilename);

    if(bovw->dictLearnMasksFileset) {
        sprintf(auxFilename, "%s/train_masks_fileset.csv", tmpDir);
        iftWriteFileSetAsCSV(bovw->dictLearnMasksFileset, auxFilename);
    }

    /* create the zip file */
    iftZipDirContent(tmpDir, filename);
    iftRemoveDir(tmpDir);
}

iftBagOfVisualWords *iftBovwRead(const char *filename)
{
    if(!iftCompareStrings(iftFileExt(filename), ".bovw"))
        iftError("Invalid file extension", "iftBovwRead");

    char auxFilename[2048];
    char *tmpDir = iftMakeTempDir("tmpdir_", "bovw", "./");
    iftUnzipFile(filename, tmpDir);

    iftBagOfVisualWords *bovw = iftCreateBovw();

    /* read json file with dict params */
    sprintf(auxFilename, "%s/dict_params.json", tmpDir);
    iftDict *dictParams = iftReadJson(auxFilename);

    bovw->intPointDetectorId = iftGetLongValFromDict("int_point_detector_id", dictParams);
    bovw->localFeatExtractorId = iftGetLongValFromDict("local_feat_extractor_id", dictParams);
    bovw->dictEstimatorId = iftGetLongValFromDict("dict_estimator_id", dictParams);
    bovw->codFuncId = iftGetLongValFromDict("cod_func_id", dictParams);
    bovw->distFuncId = iftGetLongValFromDict("dist_func_id", dictParams);
    iftDict *funcParams = iftGetDictFromDict("func_params", dictParams);
    
    iftBovwSetFunctionPointers(bovw, bovw->intPointDetectorId, bovw->localFeatExtractorId, bovw->dictEstimatorId, bovw->codFuncId, bovw->distFuncId, funcParams);

    /* dictionary files */
    if(iftBovwIsFisherVectorsEncoding(bovw->codFuncId)) {
        if(bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_UNSUP_OPF) {
            sprintf(auxFilename, "%s/dict.model", tmpDir);
            bovw->dict = iftReadBovwOPFClustModelFV(auxFilename);
        }
        if( bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION) {
            
            int nClasses = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 0), bovw->localFeats->nsamples);
            bovw->dict = (iftBovwOPFClustModelFV**)iftAlloc(nClasses+1, sizeof(iftBovwOPFClustModelFV*));
            iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)bovw->dict;
            for(int c = 1; c < nClasses+1; c++) {
                sprintf(auxFilename, "%s/dict_%d.model", tmpDir, c);
                clustModels[c] = iftReadBovwOPFClustModelFV(auxFilename);
            }
        }
    }
    else {
        sprintf(auxFilename, "%s/dict.zip", tmpDir);
        bovw->dict = iftReadDataSet(auxFilename);
    }

    /* read individual files */
    if(bovw->localFeatsOrderHrchy) {
        sprintf(auxFilename, "%s/local_feats_order_hrchy.npy", tmpDir);
        bovw->localFeatsOrderHrchy = iftReadIntMatrix(auxFilename);
    }

    if(bovw->orderHrchies) {
        sprintf(auxFilename, "%s/order_hrchies.npy", tmpDir);
        bovw->orderHrchies = iftReadIntArray(auxFilename);
    }

    sprintf(auxFilename, "%s/dict_learn_fileset.csv", tmpDir);
    bovw->dictLearnFileset = iftLoadFileSetFromCSV(auxFilename, true);

    if(bovw->dictLearnMasksFileset) {
        sprintf(auxFilename, "%s/train_masks_fileset.csv", tmpDir);
        bovw->dictLearnMasksFileset = iftLoadFileSetFromCSV(auxFilename, true);
    }

    iftRemoveDir(tmpDir);

    return bovw;
}

iftBovwOPFClustModelFV *iftCreateBovwOPFClustModelFV(iftKnnGraph *graph, iftMatrix *means, iftMatrix **covariances, iftMatrix **covariancesInv, iftFloatArray *covariancesDet, iftFloatArray *weights)
{
    iftBovwOPFClustModelFV *model = (iftBovwOPFClustModelFV*)iftAlloc(1, sizeof(iftBovwOPFClustModelFV));
    model->graph = graph;
    model->means = means;
    model->covariances = covariances;
    model->covariancesInv = covariancesInv;
    model->covariancesDet = covariancesDet;
    model->weights = weights;
    model->nFeatsPerWord = means->ncols;

    return model;
}

void iftDestroyBovwOPFClustModelFV(iftBovwOPFClustModelFV **model)
{
    if((*model)->graph != NULL)
        iftDestroyKnnGraph(&((*model)->graph));

    if((*model)->means != NULL)
        iftDestroyMatrix(&((*model)->means));

    if((*model)->covariancesDet != NULL)
        iftDestroyFloatArray(&((*model)->covariancesDet));

    if((*model)->covariances != NULL) {
        for(int g = 0; g < (*model)->weights->n; g++)
            iftDestroyMatrix(&((*model)->covariances[g]));
        iftFree((*model)->covariances);
    }

    if((*model)->covariancesInv != NULL) {
        for(int g = 0; g < (*model)->weights->n; g++)
            iftDestroyMatrix(&((*model)->covariancesInv[g]));
        iftFree((*model)->covariancesInv);
    }

    if((*model)->weights != NULL)
        iftDestroyFloatArray(&((*model)->weights));
}

void iftWriteBovwOPFClustModelFV(iftBovwOPFClustModelFV *model, char *filename)
{
    if(!iftCompareStrings(iftFileExt(filename), ".model"))
        iftError("Invalid file extension", "iftWriteBovwOPFClustModelFV");

    char *tmpDir = iftMakeTempDir("tmpdir_", "clust_model", "./");
    char auxFilename[2048];

    if(model->graph) {
        sprintf(auxFilename, "%s/knn-graph.graph", tmpDir);
        iftWriteKnnGraph(model->graph, auxFilename);
    }

    if(model->means) {
        sprintf(auxFilename, "%s/means.npy", tmpDir);
        iftWriteMatrix(model->means, auxFilename);
    }

    if(model->weights) {
        sprintf(auxFilename, "%s/weights.npy", tmpDir);
        iftWriteFloatArray(model->weights, auxFilename);
    }

    if(model->covariances) {
        for(int g = 0; g < model->weights->n; g++) {
            sprintf(auxFilename, "%s/covariances_%d.npy", tmpDir, g+1);
            iftWriteMatrix(model->covariances[g], auxFilename);
        }
    }

    if(model->covariancesInv) {
        for(int g = 0; g < model->weights->n; g++) {
            sprintf(auxFilename, "%s/covariances_inv_%d.npy", tmpDir, g+1);
            iftWriteMatrix(model->covariancesInv[g], auxFilename);
        }
    }

    if(model->covariancesDet) {
        sprintf(auxFilename, "%s/covariances_det.npy", tmpDir);
        iftWriteFloatArray(model->covariancesDet, auxFilename);
    }

    iftZipDirContent(tmpDir, filename);
    iftRemoveDir(tmpDir);
}

iftBovwOPFClustModelFV *iftReadBovwOPFClustModelFV(char *filename)
{
    if(!iftCompareStrings(iftFileExt(filename), ".model"))
        iftError("Invalid file extension", "iftReadBovwOPFClustModelFV");

    char auxFilename[2048];
    char *tmpDir = iftMakeTempDir("tmpdir_", "clust_model", "./");
    iftUnzipFile(filename, tmpDir);

    sprintf(auxFilename, "%s/knn-graph.graph", tmpDir);
    iftKnnGraph *graph = iftReadKnnGraph(auxFilename);

    sprintf(auxFilename, "%s/means.npy", tmpDir);
    iftMatrix *means = iftReadMatrix(auxFilename);

    sprintf(auxFilename, "%s/weights.npy", tmpDir);
    iftFloatArray *weights = iftReadFloatArray(auxFilename);

    iftMatrix **covariances = (iftMatrix**)iftAlloc(weights->n, sizeof(iftMatrix*));
    iftMatrix **covariancesInv = (iftMatrix**)iftAlloc(weights->n, sizeof(iftMatrix*));
    for(int g = 0; g < weights->n; g++) {
        sprintf(auxFilename, "%s/covariances_%d.npy", tmpDir, g+1);
        covariances[g] = iftReadMatrix(auxFilename);
        sprintf(auxFilename, "%s/covariances_inv_%d.npy", tmpDir, g+1);
        covariancesInv[g] = iftReadMatrix(auxFilename);
    }

    sprintf(auxFilename, "%s/covariances_det.npy", tmpDir);
    iftFloatArray *covariancesDet = iftReadFloatArray(auxFilename);

    iftBovwOPFClustModelFV *model = iftCreateBovwOPFClustModelFV(graph, means, covariances, covariancesInv, covariancesDet, weights);

    iftRemoveDir(tmpDir);

    return model;
}

/* ---------------------------------------------------------------------------------------------------------------------*/
/*                                                Learning functions                                               */
/* ---------------------------------------------------------------------------------------------------------------------*/

float iftBovwLearnImgClassif(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset, iftFileSet *dictLearnMasksFileset, char *patchDirPath, bool saveImgIntPts)
{
    bovw->dictLearnFileset = iftCopyFileSet(dictLearnFileset);
    if(dictLearnMasksFileset != NULL) bovw->dictLearnMasksFileset = iftCopyFileSet(dictLearnMasksFileset);
    int numbTrainImg = bovw->dictLearnFileset->n;

    /* create the local feats dataset */
    iftImage *img0 = iftReadImageByExt(bovw->dictLearnFileset->files[0]->path);
    iftUpdateDataSetFileInfo(bovw->dictLearnFileset->files[0]);
    iftImage *imgMask0 = NULL;
    if(dictLearnMasksFileset != NULL) imgMask0 = iftReadImageByExt(bovw->dictLearnMasksFileset->files[0]->path);
    iftRoiArray *roiArray0 = bovw->intPointDetector(img0, imgMask0, bovw->dictLearnFileset->files[0], false, bovw->funcParams);
    iftDataSet *localFeats0 = bovw->localFeatExtractor(img0, roiArray0, bovw->dictLearnFileset->files[0],patchDirPath , bovw->funcParams);
    iftDataSet *localFeats = iftCreateDataSet(numbTrainImg*localFeats0->nsamples, localFeats0->nfeats);

    /* create the label heararchy */
    int nOrderHrchies = bovw->orderHrchies ? bovw->orderHrchies->n : 0;
    iftIntMatrix *localFeatsOrderHrchy = NULL;
    if(nOrderHrchies > 0) {
        iftIntMatrix *orderHrchy0 = iftBovwComputeOrderHrchiesMatrixForImg(localFeats0, bovw->orderHrchies,
            bovw->dictLearnFileset, bovw->dictLearnMasksFileset, 0);
        localFeatsOrderHrchy = iftCreateIntMatrix(numbTrainImg*orderHrchy0->ncols, orderHrchy0->nrows);
        iftDestroyIntMatrix(&orderHrchy0);
    }

    /* execute the BoVW pipeline for every image */
    printf("- Executing Interest Point Detector and Local Feature Extractor methods ...\n");
    printf("Num. Feat vectors per image: %d (with %d feats each)\n", localFeats0->nsamples, localFeats0->nfeats);
    timer *t1 = iftTic();
    int progress = 0;
    #ifndef IFT_GPU
    #pragma omp parallel for shared(progress)
    #endif
    for(int f = 0; f < numbTrainImg; f++) {
        printf("Progress: %.2f%s\r", ((float)progress/(float)numbTrainImg*100.0), "%"); fflush(stdout);
        iftImage *img = iftReadImageByExt(bovw->dictLearnFileset->files[f]->path);
        iftUpdateDataSetFileInfo(bovw->dictLearnFileset->files[f]);
        iftImage *imgMask = NULL;
        if(dictLearnMasksFileset != NULL) imgMask = iftReadImageByExt(bovw->dictLearnMasksFileset->files[f]->path);

        /* detect and extract the local features from the image */
        iftRoiArray *roiArray = bovw->intPointDetector(img, imgMask, bovw->dictLearnFileset->files[f], saveImgIntPts, bovw->funcParams);
        iftDataSet *localFeatsAux = bovw->localFeatExtractor(img, roiArray, bovw->dictLearnFileset->files[f],patchDirPath , bovw->funcParams);
        iftCopyDataSetToDataSet(localFeatsAux, localFeats, f*localFeatsAux->nsamples);

        /* create the matrix for the label hearchy */
        if(nOrderHrchies > 0) {
            iftIntMatrix *orderHrchyAux = iftBovwComputeOrderHrchiesMatrixForImg(localFeatsAux, bovw->orderHrchies,
                bovw->dictLearnFileset, bovw->dictLearnMasksFileset, f);
            for(int col = 0; col < orderHrchyAux->ncols; col++)
                for(int row = 0; row < orderHrchyAux->nrows; row++)
                    iftMatrixElem(localFeatsOrderHrchy, f*orderHrchyAux->ncols+col, row) = iftMatrixElem(orderHrchyAux, col, row);
            iftDestroyIntMatrix(&orderHrchyAux);
        }

        progress++;
        iftDestroyImage(&img);
        if(imgMask != NULL) iftDestroyImage(&imgMask);
        iftDestroyRoiArray(&roiArray);
        iftDestroyDataSet(&localFeatsAux);
    }
    printf("Progress: 100.00%s\n", "%");

    float procTime = iftCompTime(t1, iftToc())/1000.0;
    printf("Time to extract the local feature vectors: %s\n", iftFormattedTime(procTime*1000.0));
    iftSetDistanceFunction(localFeats, bovw->distFuncId);
    iftCountNumberOfClassesDataSet(localFeats);

    if(nOrderHrchies > 0) {
        iftDict *dictAux = iftGetDictFromDict("dict_estim", bovw->funcParams);
        iftInsertIntoDict("local_feats_order_hrchy", localFeatsOrderHrchy, dictAux);
    }

    bovw->localFeats = localFeats;
    bovw->localFeatsOrderHrchy = localFeatsOrderHrchy;
    printf("Num. Image patches extracted: %d (%d per image)\n", bovw->localFeats->nsamples, bovw->localFeats->nsamples/numbTrainImg);

    /* free memory */
    iftDestroyImage(&img0);
    iftDestroyImage(&imgMask0);
    iftDestroyRoiArray(&roiArray0);
    iftDestroyDataSet(&localFeats0);

    /* dictionary estimation */
    bool isFisherVectorsEncoding = iftBovwIsFisherVectorsEncoding(bovw->codFuncId);
    if(bovw->dictEstimator != NULL) {
        printf("\n- Executing %s method ...\n", !isFisherVectorsEncoding ? "Dictionary Estimation" : "Clustering Model Estimation"); 
        timer *t2 = iftTic();
        bovw->dict = bovw->dictEstimator(localFeats, bovw->funcParams, isFisherVectorsEncoding);
        float dictTime = iftCompTime(t2, iftToc())/1000.0;
        procTime += dictTime;
        printf("Time to compute the %s: %s\n", !isFisherVectorsEncoding ? "visual dictionary" : "clustering model",
            iftFormattedTime(dictTime*1000.0));
        printf("Num. %s computed: %d\n", !isFisherVectorsEncoding ? "visual words" : "clusters", iftBovwGetNumbVisualWords(bovw));
    }

    return procTime;
}

float iftBovwLearnIntPointDetec(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset, iftFileSet *dictLearnMasksFileset, char *patchDirPath, bool saveImgIntPts)
{
    bovw->dictLearnFileset = iftCopyFileSet(dictLearnFileset);
    if(dictLearnMasksFileset != NULL) bovw->dictLearnMasksFileset = iftCopyFileSet(dictLearnMasksFileset);
    int numbTrainImg = bovw->dictLearnFileset->n;

    /* execute the BoVW interest point detection for every image */
    printf("- Executing Interest Point Detector method ...\n");
    timer *t1 = iftTic();
    int progress = 0, nIntPoints = 0;
    #ifndef IFT_GPU
    #pragma omp parallel for shared(progress)
    #endif
    for(int f = 0; f < numbTrainImg; f++) {
        printf("Progress: %.2f%s\r", ((float)progress/(float)numbTrainImg*100.0), "%"); fflush(stdout);
        iftImage *img = iftReadImageByExt(bovw->dictLearnFileset->files[f]->path);
        iftUpdateDataSetFileInfo(bovw->dictLearnFileset->files[f]);
        iftImage *imgMask = NULL;
        if(dictLearnMasksFileset != NULL) imgMask = iftReadImageByExt(bovw->dictLearnMasksFileset->files[f]->path);

        /* detect interest points from the image */
        iftRoiArray *roiArray = bovw->intPointDetector(img, imgMask, bovw->dictLearnFileset->files[f], saveImgIntPts, bovw->funcParams);
        nIntPoints += roiArray->n;

        progress++;
        iftDestroyImage(&img);
        if(imgMask != NULL) iftDestroyImage(&imgMask);
        iftDestroyRoiArray(&roiArray);
    }
    printf("Progress: 100.00%s\n", "%");

    float procTime = iftCompTime(t1, iftToc())/1000.0;
    printf("Time to detect the interest points: %s\n", iftFormattedTime(procTime*1000.0));
    printf("Num. Interest points detected: %d (%d per image)\n", nIntPoints, nIntPoints/numbTrainImg);

    return procTime;
}

iftFeatures *iftBovwComputeFeatVectFromImage(iftBagOfVisualWords *bovw, iftImage *img, iftImage *imgMask, iftFile *imgFilePtr, bool saveImgIntPts)
{
    if(bovw->dict == NULL)
        iftError("The given BoVW structure does not contain a valid Visual Dictionary.", "iftBovwComputeFeatVectFromImage");

    /* detect and extract the local features from the image */
    iftRoiArray *roiArray = bovw->intPointDetector(img, imgMask, imgFilePtr, saveImgIntPts, bovw->funcParams);
    iftDataSet *feats = bovw->localFeatExtractor(img, roiArray, imgFilePtr, NULL, bovw->funcParams);

    /* execute coding and pooling functions */
    iftFeatures *featVect = bovw->codFunc(bovw->dict, feats, bovw->funcParams);

    iftDestroyDataSet(&feats);
    iftDestroyRoiArray(&roiArray);

    return featVect;
}

iftDataSet *iftBovwComputeFeatVectsFromImgFileSet(iftBagOfVisualWords *bovw, iftFileSet *imgFileSet, iftFileSet *imgMasksFileset, bool saveImgIntPts)
{
    if(bovw->dict == NULL)
        iftError("The given BoVW structure does not contain a valid Visual Dictionary.", "iftBovwComputeFeatVectsFromImgFileSet");

    /* compute feature vectors for all the images in the imgFileSet and copy them to the dataset */
    iftDataSet *featVects = iftCreateDataSet(imgFileSet->n, iftBovwComputeNumbFeatsForCodFunc(bovw));
    featVects->nclasses = iftFileSetLabelsNumber(imgFileSet);
    featVects->function_number = bovw->distFuncId;
    iftSetDistanceFunction(featVects, bovw->distFuncId);
    featVects->ref_data_type = IFT_REF_DATA_FILESET;
    featVects->ref_data = imgFileSet;

    timer *t1 = iftTic();
    int progress = 0;
#ifndef IFT_GPU
    #pragma omp parallel for shared(progress)
#endif
    for(int s = 0; s < imgFileSet->n; s++) {
        printf("Progress: %.2f%s\r", ((float)progress/(float)imgFileSet->n*100.0), "%"); fflush(stdout);
        featVects->sample[s].id = s;
        featVects->sample[s].truelabel = imgFileSet->files[s]->label;
        iftImage *img = iftReadImageByExt(imgFileSet->files[s]->path);
        iftUpdateDataSetFileInfo(imgFileSet->files[s]);
        iftImage *imgMask = NULL;
        if(imgMasksFileset != NULL)
            imgMask = iftReadImageByExt(imgMasksFileset->files[s]->path);

        iftFeatures *feats = iftBovwComputeFeatVectFromImage(bovw, img, imgMask, imgFileSet->files[s], saveImgIntPts);
        iftCopyFloatArray(featVects->sample[s].feat, feats->val, feats->n);

        progress++;
        iftDestroyImage(&img);
        if(imgMask != NULL) iftDestroyImage(&imgMask);
        iftDestroyFeatures(&feats);
    }
    printf("Progress: 100.00%s\n", "%");
    printf("Time to compute the feature vectors: %s\n", iftFormattedTime(iftCompTime(t1, iftToc()))); 

    return featVects;
}

iftDataSet *iftBovwComputeFeatVectsFromImgFileSetBatch(iftBagOfVisualWords *bovw, iftFileSet *imgFileSet, iftFileSet *imgMasksFileset, int batchSize, bool saveImgIntPts)
{
    if(bovw->dict == NULL)
        iftError("The given BoVW structure does not contain a valid Visual Dictionary.", "iftBovwComputeFeatVectsFromImgFileSetBatch");

    /* compute feature vectors for all the images in the imgFileSet and copy them to the dataset */
    iftDataSet *featVects = iftCreateDataSet(imgFileSet->n, iftBovwComputeNumbFeatsForCodFunc(bovw));
    featVects->nclasses = iftFileSetLabelsNumber(imgFileSet);
    featVects->function_number = bovw->distFuncId;
    iftSetDistanceFunction(featVects, bovw->distFuncId);
    featVects->ref_data_type = IFT_REF_DATA_FILESET;
    featVects->ref_data = imgFileSet;

    int nImgs = imgFileSet->n;
    int nBatches = ceil((float)nImgs / (float)batchSize);

    iftImage *img0 = iftReadImageByExt(imgFileSet->files[0]->path);
    iftImage *imgMask0 = NULL;
    if(imgMasksFileset != NULL) imgMask0 = iftReadImageByExt(imgMasksFileset->files[0]->path);
    iftRoiArray *roiArray0 = bovw->intPointDetector(img0, imgMask0, imgFileSet->files[0], false, bovw->funcParams);
    iftDataSet *feats0 = bovw->localFeatExtractor(img0, roiArray0, imgFileSet->files[0], NULL, bovw->funcParams);
    int nPatches = roiArray0->n;
    int localFeatsDim = feats0->nfeats;

    timer *t1 = iftTic();
    printf("\n");
    for(int b = 0; b < nBatches; b++) {
        printf("Batch %d/%d\n", b+1, nBatches);
        int nImgsInBatch = iftMin((b+1)*batchSize, nImgs) - b*batchSize;
        iftDataSet *vectorsBatch = iftCreateDataSet(nImgsInBatch*nPatches, localFeatsDim); // it will contain the local feat vectors for all the images in the current batch
        
        /* for all the images in the current batch */
        int row = 0;
        for (int imgId = b*batchSize, i = 0; imgId < iftMin((b+1)*batchSize, nImgs); imgId++, i++) {
            printf("- computing vector distances (image %d/%d) ... \r", i+1, nImgsInBatch); fflush(stdout);
        
            featVects->sample[imgId].id = imgId;
            featVects->sample[imgId].truelabel = imgFileSet->files[imgId]->label;
            iftImage *img = iftReadImageByExt(imgFileSet->files[imgId]->path);
            iftUpdateDataSetFileInfo(imgFileSet->files[imgId]);
            iftImage *imgMask = NULL;
            if(imgMasksFileset != NULL)
                imgMask = iftReadImageByExt(imgMasksFileset->files[imgId]->path);

            /* execute interest point detection and local feat extraction methods */
            iftRoiArray *roiArray = bovw->intPointDetector(img, imgMask, imgFileSet->files[imgId], saveImgIntPts, bovw->funcParams);
            iftDataSet *feats = bovw->localFeatExtractor(img, roiArray, imgFileSet->files[imgId], NULL, bovw->funcParams);

            /* copy the local features extracted to the batch of vectors */
            for(int f = 0; f < feats->nsamples; f++) {
                iftCopyFloatArray(vectorsBatch->sample[row].feat, feats->sample[f].feat, localFeatsDim);
                row++;
            }

            iftDestroyImage(&img);
            if(imgMask != NULL) iftDestroyImage(&imgMask);
            iftDestroyRoiArray(&roiArray);
            iftDestroyDataSet(&feats);
        }
        printf("- computing vector distances (image %d/%d) ... OK\n", nImgsInBatch, nImgsInBatch);

        /* execute the batch coding function (it will return the final feature vectors for all the images in the current batch) */
        iftDict *dictAux = iftGetDictFromDict("cod_func", bovw->funcParams);
        iftInsertIntoDict("n_imgs_in_batch", nImgsInBatch, dictAux);
        iftFeatures *featVectsBatch = bovw->codFunc(bovw->dict, vectorsBatch, bovw->funcParams);
        iftRemoveValFromDict("n_imgs_in_batch", dictAux);

        /* copy the feature vectors from the batch to the final dataset */
        for (int imgId = b*batchSize, i = 0; imgId < iftMin((b+1)*batchSize, nImgs); imgId++, i++) {
            for(int w = 0; w < featVects->nfeats; w++) {
                int col = i*featVects->nfeats + w;
                featVects->sample[imgId].feat[w] = featVectsBatch->val[col];
            }
        }

        iftDestroyDataSet(&vectorsBatch);
        iftDestroyFeatures(&featVectsBatch);
    }
    printf("\nTime to compute the feature vectors: %s\n", iftFormattedTime(iftCompTime(t1, iftToc()))); 

    return featVects;
}

/* ---------------------------------------------------------------------------------------------------------------------*/
/*                                            Function pointers Implementation                                          */
/* ---------------------------------------------------------------------------------------------------------------------*/

iftRoiArray* iftBovwRandomIntPointDetector(iftImage* img, iftImage* imgMask, iftFile *imgFilePtr, bool saveImgIntPts, iftDict *funcParams)
{
    iftDict *dictAux = iftGetDictFromDict("int_point_detec", funcParams);
    iftRoiArray *roiArray = NULL;

    /* verify if the interest points were previously detected (during the learning phase) */
    if(iftDictContainKey("int_points_file", dictAux, NULL)) {
        char *intPointsFile = iftGetStrValFromDict("int_points_file", dictAux);
        roiArray = iftReadRoiArray(intPointsFile);
    }
    /* if the interest points were not previously detected, compute them */
    else {
        /* get the parameters */
        int nPoints = iftGetLongValFromDict("n_points", dictAux);
        int patchSizeInt = iftGetLongValFromDict("patch_size", dictAux);

        /* compute the interest points */
        roiArray = iftCreateRoiArray(nPoints);
        iftSize patchSize = {.x = patchSizeInt, .y = patchSizeInt, .z = iftIs3DImage(img) ? patchSizeInt : 1};

        if(imgMask) {
            iftImage *imgMaskThres = iftThreshold(imgMask,128,255,1);
            int i = 0;
            while(i < nPoints) {
                iftVoxel v;
                v.x = iftRandomInteger(0, img->xsize - iftMax(patchSize.x/2, 1));
                v.y = iftRandomInteger(0, img->ysize - iftMax(patchSize.y/2, 1));
                v.z = iftRandomInteger(0, img->zsize - iftMax(patchSize.z/2, 1));

                iftBoundingBox bb;
                if(iftImgVal(imgMaskThres, v.x, v.y, v.z) == 1) {
                    bb.begin.x = iftMax(v.x-patchSize.x/2, 0);
                    bb.begin.y = iftMax(v.y-patchSize.y/2, 0);
                    bb.begin.z = iftMax(v.z-patchSize.z/2, 0);
                    bb.end.x = iftMin(v.x+patchSize.x/2, img->xsize-1);
                    bb.end.y = iftMin(v.y+patchSize.y/2, img->ysize-1);
                    bb.end.z = iftMin(v.z+patchSize.z/2, img->zsize-1);
                    roiArray->val[i] = iftRoiFromBoundingBox(bb);
                    roiArray->val[i]->label = imgFilePtr->label;
                    i++;
                }
            }
        }
        else {
            for(int i = 0; i < nPoints; ++i) {
                iftVoxel v;
                v.x = iftRandomInteger(0, img->xsize - iftMax(patchSize.x/2, 1));
                v.y = iftRandomInteger(0, img->ysize - iftMax(patchSize.y/2, 1));
                v.z = iftRandomInteger(0, img->zsize - iftMax(patchSize.z/2, 1));

                iftBoundingBox bb;
                bb.begin.x = iftMax(v.x-patchSize.x/2, 0);
                bb.begin.y = iftMax(v.y-patchSize.y/2, 0);
                bb.begin.z = iftMax(v.z-patchSize.z/2, 0);
                bb.end.x = iftMin(v.x+patchSize.x/2, img->xsize-1);
                bb.end.y = iftMin(v.y+patchSize.y/2, img->ysize-1);
                bb.end.z = iftMin(v.z+patchSize.z/2, img->zsize-1);
                roiArray->val[i] = iftRoiFromBoundingBox(bb);
                roiArray->val[i]->label = imgFilePtr->label;
            }
        }

        /* create the interest points file and insert it into the funcParams dict */
        iftBoVWCreateIntPointsFile(roiArray, funcParams);
    }

    /* save image with interest points */
    if(saveImgIntPts)
        iftBovwSaveImgWithIntPoints(img, roiArray, imgFilePtr);

    return roiArray;
}

iftRoiArray* iftBovwGridIntPointDetector(iftImage* img, iftImage* imgMask, iftFile *imgFilePtr, bool saveImgIntPts, iftDict *funcParams)
{
    iftDict *dictAux = iftGetDictFromDict("int_point_detec", funcParams);
    iftRoiArray *roiArray = NULL;

    /* verify if the interest points were previously detected (during the learning phase) */
    if(iftDictContainKey("int_points_file", dictAux, NULL)) {
        char *intPointsFile = iftGetStrValFromDict("int_points_file", dictAux);
        roiArray = iftReadRoiArray(intPointsFile);
    }
    /* if the interest points were not previously detected, compute them */
    else {
        /* get the parameters */
        int patchSizeInt = iftGetLongValFromDict("patch_size", dictAux);
        int strideInt = iftGetLongValFromDict("stride", dictAux);

        /* compute the interest points */
        iftSize patchSize = {.x = patchSizeInt, .y = patchSizeInt, .z = iftIs3DImage(img) ? patchSizeInt : 1};
        iftSize stride = {.x = strideInt, .y = strideInt, .z = iftIs3DImage(img) ? strideInt : 0};
        iftIntArray *samplePos = NULL;
        if(imgMask) {
            samplePos = iftGridSamplingOnMask(imgMask, strideInt, -1, -1);
        }
        else {
            iftImageDomain imgDom = {.xsize = img->xsize, .ysize = img->ysize, .zsize = img->zsize};
            samplePos = iftGridSamplingForPatchExtraction(imgDom, patchSize.x, patchSize.y, patchSize.z, stride.x, stride.y, stride.z);
        }
        int nPoints = samplePos->n;
        roiArray = iftCreateRoiArray(nPoints);
        
        for(int p = 0; p < samplePos->n; p++) {
            iftVoxel v = iftGetVoxelCoord(img, samplePos->val[p]);
            iftBoundingBox bb;
            bb.begin.x = iftMax(v.x - patchSize.x/2, 0);
            bb.begin.y = iftMax(v.y - patchSize.y/2, 0);
            bb.begin.z = iftMax(v.z - patchSize.z/2, 0);
            bb.end.x = iftMin(v.x + patchSize.x/2, img->xsize-1);
            bb.end.y = iftMin(v.y + patchSize.y/2, img->ysize-1);
            bb.end.z = iftMin(v.z + patchSize.z/2, img->zsize-1);
            roiArray->val[p] = iftRoiFromBoundingBox(bb);
            roiArray->val[p]->label = imgFilePtr->label;
        }

        /* create the interest points file and insert it into the funcParams dict */
        iftBoVWCreateIntPointsFile(roiArray, funcParams);
    }

    /* save image with interest points */
    if(saveImgIntPts)
        iftBovwSaveImgWithIntPoints(img, roiArray, imgFilePtr);

    return roiArray;
}

iftRoiArray* iftBovwUnsupSuperpixelISFIntPointDetector(iftImage* img, iftImage* imgMask, iftFile *imgFilePtr, bool saveImgIntPts, iftDict *funcParams)
{
    iftDict *dictAux = iftGetDictFromDict("int_point_detec", funcParams);
    iftRoiArray *roiArray = NULL;

    /* verify if the interest points were previously detected (during the learning phase) */
    if(iftDictContainKey("int_points_file", dictAux, NULL)) {
        char *intPointsFile = iftGetStrValFromDict("int_points_file", dictAux);
        roiArray = iftReadRoiArray(intPointsFile);
    }
    /* if the interest points were not previously detected, compute them */
    else {
        /* get the parameters */
        iftDict *dictAux = iftGetDictFromDict("int_point_detec", funcParams);
        int nSuperpixels    = iftGetLongValFromDict("n_superpixels", dictAux);
        int patchSizeInt    = iftGetLongValFromDict("patch_size", dictAux);

        /* compute supervoxels by RISF (execute pre-compiled binary) */
        iftImage *imgMaskThres = iftThreshold(imgMask,128,255,1);
        char cmd1[2048], *imgFileName, *imgMaskFileName, *spixLabelsFileName, *imgId;
        imgId = iftFilename(imgFilePtr->path, iftFileExt(imgFilePtr->path));
        imgFileName = iftConcatStrings(2, imgId, "_img.png");
        imgMaskFileName = iftConcatStrings(2, imgId, "_imgMask.png");
        spixLabelsFileName = iftConcatStrings(2, imgId, "_spixLabels.png");
        iftWriteImageByExt(img, imgFileName);
        iftWriteImageByExt(imgMaskThres, imgMaskFileName);

        /* execute RISF with 3-level hierarchy (using nSuperpixels*4, nSuperpixels*2 and nSuperpixels) */
        sprintf(cmd1, "$NEWIFT_DIR/demo/BagOfVisualWords/RISF_segmentation_3-level_hrchy.sh %s %s %s %d %s", imgId, imgFileName, imgMaskFileName, nSuperpixels, spixLabelsFileName);
        if (system(cmd1) != 0)
            iftError("Error while executing $NEWIFT_DIR/demo/BagOfVisualWords/RISF_segmentation_3-level_hrchy.sh", "iftBovwUnsupSuperpixelISFIntPointDetector");

        iftImage *spixLabels = iftReadImageByExt(spixLabelsFileName);
        
        /* remove the temporary images created */
        iftRemoveFile(imgFileName);
        iftRemoveFile(imgMaskFileName);
        iftRemoveFile(spixLabelsFileName);

        /* create a mask with the resulting superpixels' boundaries */
        iftImage *spixMask = iftCreateSuperpixelBoundariesMask(spixLabels, 255);

        /* perform grid sampling on the superpixels' boundaries */
        iftIntArray *intPts = iftGridSamplingOnMask(spixMask, (float)patchSizeInt * 1.2, -1, -1);

        /* define the bounding boxes for the points extracted */
        roiArray = iftCreateRoiArray(intPts->n);
        iftSize patchSize = {.x = patchSizeInt, .y = patchSizeInt, .z = iftIs3DImage(img) ? patchSizeInt : 1};

        for(int i = 0; i < intPts->n; i++) {
            iftVoxel v = iftGetVoxelCoord(spixLabels, intPts->val[i]);
            iftBoundingBox bb;
            bb.begin.x = iftMax(v.x-patchSize.x/2, 0);
            bb.begin.y = iftMax(v.y-patchSize.y/2, 0);
            bb.begin.z = iftMax(v.z-patchSize.z/2, 0);
            bb.end.x = iftMin(v.x+patchSize.x/2, img->xsize-1);
            bb.end.y = iftMin(v.y+patchSize.y/2, img->ysize-1);
            bb.end.z = iftMin(v.z+patchSize.z/2, img->zsize-1);
            roiArray->val[i] = iftRoiFromBoundingBox(bb);
            roiArray->val[i]->label = imgFilePtr->label;
        }
        iftDestroyImage(&spixMask);

        /* create the interest points file and insert it into the funcParams dict */
        iftBoVWCreateIntPointsFile(roiArray, funcParams);
    }

    /* save image with interest points */
    if(saveImgIntPts)
        iftBovwSaveImgWithIntPoints(img, roiArray, imgFilePtr);

    return roiArray;
}

iftRoiArray* iftBovwSupSuperpixelISFIntPointDetector(iftImage* img, iftImage* imgMask, iftFile *imgFilePtr, bool saveImgIntPts, iftDict *funcParams)
{
    iftDict *dictAux = iftGetDictFromDict("int_point_detec", funcParams);
    iftRoiArray *roiArray = NULL;

    /* verify if the interest points were previously detected (during the learning phase) */
    if(iftDictContainKey("int_points_file", dictAux, NULL)) {
        char *intPointsFile = iftGetStrValFromDict("int_points_file", dictAux);
        roiArray = iftReadRoiArray(intPointsFile);
    }
    /* if the interest points were not previously detected, compute them */
    else {
        /* get the parameters */
        iftDict *dictAux = iftGetDictFromDict("int_point_detec", funcParams);
        int patchSize = iftGetLongValFromDict("patch_size", dictAux);
        int intPtsStride = iftGetLongValFromDict("int_pts_stride", dictAux);
        char *pointPosition = iftGetStrValFromDict("point_position", dictAux);
        char *patchShape = iftGetStrValFromDict("patch_shape", dictAux);

        if(iftCompareStrings(patchShape, "superpixel") && !iftCompareStrings(pointPosition, "center"))
            iftError("'Superpixel' patch shape representation can only be used with \"center\" point position",
                    "iftBovwSupSuperpixelISFIntPointDetector");

        /* extract the ROIs using the superpixel masks/labels for all classes and return the ROIs for all classes */
        iftDict *sampMasksDict = iftGetDictFromDict("sampling_masks", funcParams);
        int nClasses = iftGetLongValFromDict("n_classes", dictAux);
        
        iftRoiArray **roiArrayArray = (iftRoiArray**)iftAlloc(nClasses+1, sizeof(iftRoiArray*));

        int nIntPts = 0;
        for(int c = 1; c < nClasses+1; c++) {
            if(iftCompareStrings(pointPosition, "center")) {
                /* get superpixel labels for the given class */
                char className[2048]; sprintf(className, "class_%d", c);
                iftDict *spixLabelsPerClassDict = iftGetDictFromDict("spix_labels_per_class", sampMasksDict);
                char *spixLabelsFilename = iftGetStrValFromDict(className, spixLabelsPerClassDict);
                iftImage *spixLabels = iftReadImageByExt(spixLabelsFilename);

                /* define the ROIs using the superpixel labels */
                int nSpix = iftMaximumValue(spixLabels);
                roiArrayArray[c] = iftCreateRoiArray(nSpix);

                if(iftCompareStrings(patchShape, "superpixel")) {
                    for(int i = 0; i < nSpix; i++)
                        roiArrayArray[c]->val[i] = iftRoiFromSuperpixel(spixLabels, i+1);
                }
                else if(iftCompareStrings(patchShape, "square")) {
                    iftSize patch = {.x = patchSize, .y = patchSize, .z = iftIs3DImage(img) ? patchSize : 1};
                    for(int i = 0; i < nSpix; i++) {
                        iftRoi *roi = iftRoiFromSuperpixel(spixLabels, i+1);
                        iftVoxel v = iftRoiCenterVoxel(img, roi);
                        iftBoundingBox bb;
                        bb.begin.x = iftMax(v.x-patch.x/2, 0);
                        bb.begin.y = iftMax(v.y-patch.y/2, 0);
                        bb.begin.z = iftMax(v.z-patch.z/2, 0);
                        bb.end.x = iftMin(v.x+patch.x/2, img->xsize-1);
                        bb.end.y = iftMin(v.y+patch.y/2, img->ysize-1);
                        bb.end.z = iftMin(v.z+patch.z/2, img->zsize-1);
                        roiArrayArray[c]->val[i] = iftRoiFromBoundingBox(bb);
                    }
                }

                nIntPts += nSpix;
                iftDestroyImage(&spixLabels);
            }
            else if(iftCompareStrings(pointPosition, "boundary")) {
                /* get superpixel mask for the given class */
                char className[2048]; sprintf(className, "class_%d", c);
                iftDict *spixMasksPerClassDict = iftGetDictFromDict("spix_masks_per_class", sampMasksDict);
                char *spixMaskFilename = iftGetStrValFromDict(className, spixMasksPerClassDict);
                iftImage *spixMask = iftReadImageByExt(spixMaskFilename);
                
                /* perform grid sampling on the superpixels' boundaries */
                iftIntArray *intPts = iftGridSamplingOnMask(spixMask, intPtsStride, -1, -1); // we use all the points resulting from the given int_pts_stride

                /* define the ROIs for the points extracted */
                roiArrayArray[c] = iftCreateRoiArray(intPts->n);
                iftSize patch = {.x = patchSize, .y = patchSize, .z = iftIs3DImage(img) ? patchSize : 1};

                for(int i = 0; i < intPts->n; i++) {
                    iftVoxel v = iftGetVoxelCoord(spixMask, intPts->val[i]);
                    iftBoundingBox bb;
                    bb.begin.x = iftMax(v.x-patch.x/2, 0);
                    bb.begin.y = iftMax(v.y-patch.y/2, 0);
                    bb.begin.z = iftMax(v.z-patch.z/2, 0);
                    bb.end.x = iftMin(v.x+patch.x/2, img->xsize-1);
                    bb.end.y = iftMin(v.y+patch.y/2, img->ysize-1);
                    bb.end.z = iftMin(v.z+patch.z/2, img->zsize-1);
                    roiArrayArray[c]->val[i] = iftRoiFromBoundingBox(bb);
                }
                nIntPts += intPts->n;
                iftDestroyIntArray(&intPts);
                iftDestroyImage(&spixMask);
            }
        }
        
        /* join the ROIs in one iftRoiArray and create an array with the labels */
        roiArray = iftCreateRoiArray(nIntPts);
        int s = 0;

        for(int c = 1; c < nClasses+1; c++) {
            for(int i = 0; i < roiArrayArray[c]->n; i++) {
                roiArray->val[s] = iftCopyRoi(roiArrayArray[c]->val[i]);
                roiArray->val[s]->label = c;
                s++;
            }
        }
        for(int c = 1; c < nClasses+1; c++)
            iftDestroyRoiArray(&roiArrayArray[c]);
        iftFree(roiArrayArray);

        /* create the interest points file and insert it into the funcParams dict */
        iftBoVWCreateIntPointsFile(roiArray, funcParams);
    }

    /* save image with interest points */
    if(saveImgIntPts)
        iftBovwSaveImgWithIntPoints(img, roiArray, imgFilePtr);

    return roiArray;
}

iftDataSet* iftBovwRawFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", funcParams);
    char *colSpaceStr = iftGetStrValFromDict("col_space", dictAux);
    iftColorSpace colSpace = iftColorSpaceStrToColorSpace(colSpaceStr);

    iftColorSpace inputColSpace = YCbCr_CSPACE;
    
    /* validate the color spaces */
    bool isColorImg = iftIsColorImage(img);

    if((isColorImg && (inputColSpace == GRAY_CSPACE || inputColSpace == GRAYNorm_CSPACE)) ||
        (!isColorImg && (inputColSpace != GRAY_CSPACE && inputColSpace != GRAYNorm_CSPACE))) {
        iftError("The number of image channels and the specified col_space_in does not match",
            "iftBovwRawFeatExtractor");
    }

    if(!isColorImg && colSpace != GRAY_CSPACE && colSpace != GRAYNorm_CSPACE) {
        iftError("Color RAW image representation can only be used with color images. The given image is grayscale",
            "iftBovwRawFeatExtractor");
    }

    int nFeats = 0;
    /* the number of features depends on the size of the biggest ROI */
    if(colSpace == GRAY_CSPACE || colSpace == GRAYNorm_CSPACE)
        nFeats = roiArray->val[iftBiggestRoiInArray(roiArray)]->n;
    else
        nFeats = roiArray->val[iftBiggestRoiInArray(roiArray)]->n*3;

    int normValue = iftNormalizationValue(iftMaximumValue(img));

    /* extract features from all the ROIs */
    iftDataSet *localFeats = iftCreateDataSet((int)roiArray->n, nFeats);
    localFeats->function_number = 1;
    iftSetDistanceFunction(localFeats, localFeats->function_number);

    #pragma omp parallel for
    for(int s = 0; s < roiArray->n; s++) {
        iftRoi *roi = roiArray->val[s];
        localFeats->sample[s].id = s;
        iftAddSampleStatus(&localFeats->sample[s], IFT_SUPERVISED);
        localFeats->sample[s].truelabel = roi->label;

        int feat = 0;
        for (int v = 0; v < roi->n; v++) {
            int pixId = iftGetVoxelIndex(img, roi->val[v]);
            iftFColor colorIn = {.val[0] = 0, .val[1] = 0, .val[2] = 0};
            colorIn.val[0] = img->val[pixId];
            if(isColorImg) {
                colorIn.val[1] = img->Cb[pixId];
                colorIn.val[2] = img->Cr[pixId];
            }
            iftFColor colorOut = iftConvertPixelColorSpace(colorIn, inputColSpace, colSpace, normValue);

            /* copy the features in the dataset */
            if(colSpace == GRAY_CSPACE || colSpace == GRAYNorm_CSPACE) {
                localFeats->sample[s].feat[feat] = colorOut.val[0];
                feat++;
            }
            else {
                localFeats->sample[s].feat[feat] = colorOut.val[0];
                localFeats->sample[s].feat[feat+1] = colorOut.val[1];
                localFeats->sample[s].feat[feat+2] = colorOut.val[2];
                feat += 3;
            }
        }
    }

    /* create fileset ref data */
    if(patchDirPath != NULL) {
        localFeats->ref_data_type = IFT_REF_DATA_FILESET;
        localFeats->ref_data = iftBoVWCreateFileSetRefDataFromImgROIs(img, roiArray, imgFilePtr, patchDirPath);
    }

    return localFeats;
}

iftDataSet* iftBovwBICFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams)
{
    /* get the dynamic parameters */
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", funcParams);
    int nBins = iftGetLongValFromDict("n_bins", dictAux);
 
    /* determine the number of features */
    int nFeats;
    if(iftIsColorImage(img))
        nFeats = nBins*nBins*nBins*2;
    else
        nFeats = nBins*2;

    /* extract features from all the ROIs */
    iftDataSet *localFeats = iftCreateDataSet((int)roiArray->n, nFeats);
    localFeats->function_number = 1;
    iftSetDistanceFunction(localFeats, localFeats->function_number);

    #pragma omp parallel for
    for(int s = 0; s < roiArray->n; s++) {
        iftRoi *roi = roiArray->val[s];
        localFeats->sample[s].id = s;
        iftAddSampleStatus(&localFeats->sample[s], IFT_SUPERVISED);
        localFeats->sample[s].truelabel = roi->label;

        /* extract BIC features */
        iftImage* imgRoi = iftExtractRoiNoBkgd(img, roi);
        iftFeatures *feats = iftExtractBIC(imgRoi, NULL, nBins);

        for (int f = 0; f < feats->n; f++)
            localFeats->sample[s].feat[f] = feats->val[f];

        iftDestroyFeatures(&feats);
        iftDestroyImage(&imgRoi);
    }

    /* create fileset ref data */
    if(patchDirPath != NULL) {
        localFeats->ref_data_type = IFT_REF_DATA_FILESET;
        localFeats->ref_data = iftBoVWCreateFileSetRefDataFromImgROIs(img, roiArray, imgFilePtr, patchDirPath);
    }

    return localFeats;
}

iftDataSet* iftBovwLBPFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams)
{
    /* get the dynamic parameters */
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", funcParams);
    double circSymNeighbRad = iftGetDblValFromDict("circ_sym_neighb_rad", dictAux);

    iftAdjRel *A = iftCircular(circSymNeighbRad);
    int nFeats = roiArray->val[iftBiggestRoiInArray(roiArray)]->n; /* the number of features is equal to the size of the biggest ROI */

    /* extract features from all the ROIs */
    iftDataSet *localFeats = iftCreateDataSet((int)roiArray->n, nFeats);
    localFeats->function_number = 1;
    iftSetDistanceFunction(localFeats, localFeats->function_number);

    #pragma omp parallel for
    for(int s = 0; s < roiArray->n; s++) {
        iftRoi *roi = roiArray->val[s];
        localFeats->sample[s].id = s;
        iftAddSampleStatus(&localFeats->sample[s], IFT_SUPERVISED);
        localFeats->sample[s].truelabel = roi->label;

        /* extract LBP features */
        iftImage* imgRoi = iftExtractRoiNoBkgd(img, roi);
        iftFeatures* feats = iftExtractLBP(imgRoi, A);

        for (int f = 0; f < feats->n; f++)
            localFeats->sample[s].feat[f] = feats->val[f];

        iftDestroyFeatures(&feats);
        iftDestroyImage(&imgRoi);
    }

    /* create fileset ref data */
    if(patchDirPath != NULL) {
        localFeats->ref_data_type = IFT_REF_DATA_FILESET;
        localFeats->ref_data = iftBoVWCreateFileSetRefDataFromImgROIs(img, roiArray, imgFilePtr, patchDirPath);
    }

    iftDestroyAdjRel(&A);

    return localFeats;
}

iftDataSet* iftBovwBriefFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams)
{
    /* get the dynamic parameters */
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", funcParams);
    int patchSize = iftGetLongValFromDict("patch_size", dictAux);
    int descriptorSize = iftGetLongValFromDict("descriptor_size", dictAux);

    /* create the auxiliary files to compute BRIEF features (execute python script) */
    char *imgId = iftFilename(imgFilePtr->path, iftFileExt(imgFilePtr->path));
    char *imgFileName = iftConcatStrings(3, imgId, "_img", iftFileExt(imgFilePtr->path));
    char *intPtsFileName = iftConcatStrings(2, imgId, "_intPts.csv");
    char *localFeatsFileName = iftConcatStrings(2, imgId, "_localFeats.csv");

    iftMatrix *intPts = iftCreateMatrix(2, roiArray->n);
    for(int i = 0; i < roiArray->n; i++) {
        iftMatrixElem(intPts, 0, i) = iftRoiCenterVoxel(img, roiArray->val[i]).y;
        iftMatrixElem(intPts, 1, i) = iftRoiCenterVoxel(img, roiArray->val[i]).x;
    }

    iftWriteImageByExt(img, imgFileName);
    iftWriteMatrixCSV(intPts, intPtsFileName);

    /* compute BRIEF features (execute python script) */
    char cmd1[2048];
    sprintf(cmd1, "python3 $NEWIFT_DIR/demo/BagOfVisualWords/iftBoVWExtractBriefFeatures.py %s %s %d %d %s", imgFileName, intPtsFileName, patchSize, descriptorSize, localFeatsFileName);
    if (system(cmd1) != 0)
        iftError("Error while executing $NEWIFT_DIR/demo/BagOfVisualWords/iftBoVWExtractBriefFeatures.py", "iftBovwBriefFeatExtractor");

    iftMatrix *localFeatsM = iftReadMatrixCSV(localFeatsFileName);
    int nFeats = localFeatsM->ncols;

    if(nFeats != descriptorSize)
        iftError("The number of feats extracted does not match with the required number of feats", "iftBovwBriefFeatExtractor");

    if(localFeatsM->nrows != roiArray->n)
        iftError("The number of feats extracted does not match with the number of interest points", "iftBovwBriefFeatExtractor");
    
    /* remove the temporary images created */
    iftRemoveFile(imgFileName);
    iftRemoveFile(intPtsFileName);
    iftRemoveFile(localFeatsFileName);

    /* create the dataset with the BRIEF features */
    iftDataSet *localFeats = iftCreateDataSet((int)roiArray->n, nFeats);
    localFeats->function_number = 1;
    iftSetDistanceFunction(localFeats, localFeats->function_number);

    #pragma omp parallel for
    for(int s = 0; s < roiArray->n; s++) {
        iftRoi *roi = roiArray->val[s];
        localFeats->sample[s].id = s;
        iftAddSampleStatus(&localFeats->sample[s], IFT_SUPERVISED);
        localFeats->sample[s].truelabel = roi->label;

        for (int f = 0; f < nFeats; f++)
            localFeats->sample[s].feat[f] = iftMatrixElem(localFeatsM, f, s);
   }

    /* create fileset ref data */
    if(patchDirPath != NULL) {
        localFeats->ref_data_type = IFT_REF_DATA_FILESET;
        localFeats->ref_data = iftBoVWCreateFileSetRefDataFromImgROIs(img, roiArray, imgFilePtr, patchDirPath);
    }

    iftDestroyMatrix(&intPts);
    iftDestroyMatrix(&localFeatsM);

    return localFeats;
}

iftDataSet* iftBovwConvolutionalFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams)
{
    /* get the dynamic parameters */
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", funcParams);
    int nKernels = iftGetLongValFromDict("n_kernels", dictAux);
    int kernelSize = iftGetLongValFromDict("kernel_size", dictAux);
    int convStride = iftGetLongValFromDict("conv_stride", dictAux);
    int poolSize = iftGetLongValFromDict("pool_size", dictAux);
    int poolStride = iftGetLongValFromDict("pool_stride", dictAux);
    char *centralizedImgsDir = iftGetStrValFromDict("centralized_imgs_dir", dictAux);
    char *kernelsDir = iftGetStrValFromDict("kernels_dir", dictAux);
    char *maxPoolType = iftGetStrValFromDict("max_pool_type", dictAux);
    char *resizeImage = iftGetStrValFromDict("resize_image", dictAux);
    char *vectConstr = iftGetStrValFromDict("vect_constr", dictAux);

    if(!iftDirExists(centralizedImgsDir))
        iftError("The directory with the centralized images does not exist", "iftBovwConvolutionalFeatExtractor");

    if(!iftDirExists(kernelsDir))
        iftError("The kernels directory does not exist", "iftBovwConvolutionalFeatExtractor");

    char filename[2048];
    sprintf(filename, "%s/%s.mimg", centralizedImgsDir, iftFilename(iftBasename(imgFilePtr->path), iftFileExt(imgFilePtr->path)));
    iftMImage *mimg = iftReadMImage(filename), *aux = NULL;

    /* apply convolution between all the kernels and all the pixels in the image  */
    int nBands = iftIsColorImage(img) ? 3 : 1;
    char *kernName = iftJoinPathnames(2, kernelsDir, "kernels_layer_1.npy");
    iftMatrix *kernels = iftReadMatrix(kernName);

    if(nKernels != kernels->ncols || kernelSize != sqrt(kernels->nrows/nBands))
        iftError("The provided kernels (%d x %d) does not match the n_kernels (%d) and kernel_size (%d) parameters\n",
                 "iftBovwConvolutionalFeatExtractor", kernels->ncols, kernels->nrows, nKernels, kernelSize);

    aux = iftMConvolution(mimg, kernels, convStride);
    iftDestroyMImage(&mimg);
    mimg = iftCopyMImage(aux);
    iftDestroyMImage(&aux);

    /* relu */
    aux = iftMReLU(mimg);
    iftDestroyMImage(&mimg);
    mimg = iftCopyMImage(aux);
    iftDestroyMImage(&aux);

    /* max-pooling */
    if(iftCompareStrings(maxPoolType, "max_pool"))
        aux = iftMMaxPooling(mimg, poolSize, poolStride);
    else if(iftCompareStrings(maxPoolType, "max_pool_roi"))
        aux = iftMMaxPoolingRoi(mimg, roiArray, poolStride);

    iftDestroyMImage(&mimg);
    mimg = iftCopyMImage(aux);
    iftDestroyMImage(&aux);

    /* build the feature vectors */
    iftDataSet *localFeats = NULL;
    if(iftCompareStrings(vectConstr, "pixel_vals_in_patch") || iftCompareStrings(vectConstr, "mean_val_in_patch")) {
        localFeats = iftBovwBuildLocalFeatureVectorsFromConvolutionalMImage(mimg, img, imgFilePtr, roiArray, resizeImage, vectConstr);
    }
    else if(iftCompareStrings(vectConstr, "bic_in_patch")) {
        int nBinsPerBand = iftGetLongValFromDict("n_bins_per_band", dictAux);
        localFeats = iftBovwExtractBICFeatureVectorsFromConvolutionalMImage(mimg, img, roiArray, resizeImage, nBinsPerBand);
    }

    /* create fileset ref data */
    if(patchDirPath != NULL) {
        if(iftCompareStrings(resizeImage, "yes")) {
            localFeats->ref_data_type = IFT_REF_DATA_FILESET;
            localFeats->ref_data = iftBoVWCreateFileSetRefDataFromImgROIs(img, roiArray, imgFilePtr, patchDirPath);
        }
        else {
            iftError("The patches for each image can not be saved when resize_image=\"no\"",
                    "iftBovwConvolutionalFeatExtractor");
        }
    }

    iftDestroyMImage(&mimg);
    iftDestroyMatrix(&kernels);

    return localFeats;
}

iftDataSet* iftBovwConvolutionalMultiLayerFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams)
{
    /* get the dynamic parameters */
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", funcParams);
    iftDict *networkParams = iftGetDictFromDict("network_params", dictAux);
    char *centralizedImgsDir = iftGetStrValFromDict("centralized_imgs_dir", dictAux);
    char *kernelsDir = iftGetStrValFromDict("kernels_dir", dictAux);
    char *maxPoolType = iftGetStrValFromDict("max_pool_type", dictAux);
    char *resizeImage = iftGetStrValFromDict("resize_image", dictAux);
    char *vectConstr = iftGetStrValFromDict("vect_constr", dictAux);

    if(!iftDirExists(centralizedImgsDir))
        iftError("The directory with the centralized images does not exist", "iftBovwConvolutionalMultiLayerFeatExtractor");

    if(!iftDirExists(kernelsDir))
        iftError("The kernels directory does not exist", "iftBovwConvolutionalMultiLayerFeatExtractor");

    char filename[2048];
    sprintf(filename, "%s/%s.mimg", centralizedImgsDir, iftFilename(iftBasename(imgFilePtr->path), iftFileExt(imgFilePtr->path)));
    iftMImage *mimg = iftReadMImage(filename), *aux = NULL;

    /* apply the operations in every layer using the current batch */
    int nLayers = iftGetLongValFromDict("n_layers", networkParams);

    for(int l = 0; l < nLayers; l++) {
        char layerName[2048];
        sprintf(layerName, "layer_%d", l+1);
        iftDict *layerDict = iftGetDictFromDict(layerName, networkParams);
        
        /* convolution */
        iftDict *convDict = iftGetDictFromDict("conv", layerDict);
        int nBands = iftGetLongValFromDict("n_bands", convDict);
        int nKernels = iftGetLongValFromDict("n_kernels", convDict);
        int kernelSize = iftGetLongValFromDict("kernel_size", convDict);
        int stride = iftGetLongValFromDict("stride", convDict);

        char kernName[2048];
        sprintf(kernName, "%s/kernels_layer_%d.npy", kernelsDir, l+1);
        iftMatrix *kernels = iftReadMatrix(kernName);

        if(nKernels != kernels->ncols || kernelSize != sqrt(kernels->nrows/nBands))
            iftError("The provided kernels for layer %d (%d x %d) does not match the n_kernels (%d) and kernel_size (%d)\n",
                    "iftBovwConvolutionalMultiLayerFeatExtractor", l+1, kernels->ncols, kernels->nrows, nKernels, kernelSize);

        aux = iftMConvolution(mimg, kernels, stride);
        iftDestroyMImage(&mimg);
        mimg = iftCopyMImage(aux);
        iftDestroyMImage(&aux);

        /* relu */
        char *relu = iftGetStrValFromDict("relu", layerDict);
        if(iftCompareStrings(relu, "yes")) {
            aux = iftMReLU(mimg);
            iftDestroyMImage(&mimg);
            mimg = iftCopyMImage(aux);
            iftDestroyMImage(&aux);
        }

        /* max-pooling */
        if(iftDictContainKey("max-pool", layerDict, NULL)) {
            iftDict *maxPoolDict = iftGetDictFromDict("max-pool", layerDict);
            int poolSize = iftGetLongValFromDict("size", maxPoolDict);
            int poolStride = iftGetLongValFromDict("stride", maxPoolDict);
            
            if(iftCompareStrings(maxPoolType, "max_pool"))
                aux = iftMMaxPooling(mimg, poolSize, poolStride);
            else if(iftCompareStrings(maxPoolType, "max_pool_roi"))
                aux = iftMMaxPoolingRoi(mimg, roiArray, poolStride);

            iftDestroyMImage(&mimg);
            mimg = iftCopyMImage(aux);
            iftDestroyMImage(&aux);
        }

        iftDestroyMatrix(&kernels);
    }

    /* build the feature vectors */
    iftDataSet *localFeats = NULL;
    if(iftCompareStrings(vectConstr, "pixel_vals_in_patch") || iftCompareStrings(vectConstr, "mean_val_in_patch")) {
        localFeats = iftBovwBuildLocalFeatureVectorsFromConvolutionalMImage(mimg, img, imgFilePtr, roiArray, resizeImage, vectConstr);
    }
    else if(iftCompareStrings(vectConstr, "bic_in_patch")) {
        int nBinsPerBand = iftGetLongValFromDict("n_bins_per_band", dictAux);
        localFeats = iftBovwExtractBICFeatureVectorsFromConvolutionalMImage(mimg, img, roiArray, resizeImage, nBinsPerBand);
    }

    /* create fileset ref data */
    if(patchDirPath != NULL) {
        if(iftCompareStrings(resizeImage, "yes")) {
            localFeats->ref_data_type = IFT_REF_DATA_FILESET;
            localFeats->ref_data = iftBoVWCreateFileSetRefDataFromImgROIs(img, roiArray, imgFilePtr, patchDirPath);
        }
        else {
            iftError("The patches for each image can not be saved when resize_image=\"no\"",
                    "iftBovwConvolutionalMultiLayerFeatExtractor");
        }
    }

    iftDestroyMImage(&mimg);

    return localFeats;
}

iftDataSet* iftBovwDeepFeaturesMImageFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams)
{
    /* get the dynamic parameters */
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", funcParams);
    char *featureMimgsDir = iftGetStrValFromDict("feature_mimgs_dir", dictAux);
    char *resizeImage = iftGetStrValFromDict("resize_image", dictAux);
    char *vectConstr = iftGetStrValFromDict("vect_constr", dictAux);

    if(iftCompareStrings(iftFileExt(featureMimgsDir), ".csv")) {
        if(!iftFileExists(featureMimgsDir))
            iftError("The directory with the feature mimages does not exist", "iftBovwDeepFeaturesMImageFeatExtractor");
        iftFileSet *fs = iftLoadFileSetFromDirOrCSV(featureMimgsDir, 1, true);
        int idx = iftFindFileByLabelAndSample(fs, imgFilePtr->label, imgFilePtr->sample);
        featureMimgsDir = iftCopyString(iftDirname(fs->files[idx]->path));
        iftDestroyFileSet(&fs);
    }
    else {
        if(!iftDirExists(featureMimgsDir))
            iftError("The directory with the feature mimages does not exist", "iftBovwDeepFeaturesMImageFeatExtractor");
    }

    char filename[2048];
    sprintf(filename, "%s/%s.mimg", featureMimgsDir, iftFilename(iftBasename(imgFilePtr->path), iftFileExt(imgFilePtr->path)));
    iftMImage *mimg = iftReadMImage(filename);

    /* build the feature vectors */
    iftDataSet *localFeats = NULL;
    if(iftCompareStrings(vectConstr, "pixel_vals_in_patch") || iftCompareStrings(vectConstr, "mean_val_in_patch")) {
        localFeats = iftBovwBuildLocalFeatureVectorsFromConvolutionalMImage(mimg, img, imgFilePtr, roiArray, resizeImage, vectConstr);
    }
    else if(iftCompareStrings(vectConstr, "bic_feats_from_patch")) {
        int nBinsPerBand = iftGetLongValFromDict("n_bins_per_band", dictAux);
        localFeats = iftBovwExtractBICFeatureVectorsFromConvolutionalMImage(mimg, img, roiArray, resizeImage, nBinsPerBand);
    }

    /* create fileset ref data */
    if(patchDirPath != NULL) {
        if(iftCompareStrings(resizeImage, "yes") || (mimg->xsize == img->xsize && mimg->ysize == img->ysize)) {
            localFeats->ref_data_type = IFT_REF_DATA_FILESET;
            localFeats->ref_data = iftBoVWCreateFileSetRefDataFromImgROIs(img, roiArray, imgFilePtr, patchDirPath);
        }
        else {
            iftError("The patches for each image can not be saved when the read feature mimages has different size than the original images",
                    "iftBovwDeepFeaturesMImageFeatExtractor");
        }
    }

    iftDestroyMImage(&mimg);

    return localFeats;
}

void* iftBovwUnsupKMeansDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("dict_estim", funcParams);
    char *retRealCentStr = iftGetStrValFromDict("ret_real_cent", dictAux);
    int nGroups = iftGetLongValFromDict("n_groups", dictAux);
    int maxIter = iftGetLongValFromDict("max_iter", dictAux);
    double minImprov = iftGetDblValFromDict("min_improv", dictAux);
    bool retRealCent = iftBoolStrToBool(retRealCentStr);

    if(isFisherVectorsEncoding)
        iftError("Fisher vectors encoding is currently not supported with k-means", "iftBovwUnsupKMeansDictEstimator");

    if(nGroups <= 0)
        iftError("Number of groups should be greater than 0.", "iftBovwUnsupKMeansDictEstimator");

    /* perform k-means clustering */
    iftRandomSeed(IFT_RANDOM_SEED);
    iftSetStatus(localFeats,IFT_TRAIN);
    iftDataSet* initCentroids = iftKmeansInitCentroidsRandomNormal(localFeats, nGroups);
    iftSimpleKmeansRun(localFeats, &initCentroids, maxIter, (float)minImprov);
    iftDataSet *dict = iftExtractCentroidsFromDataSetAsDataSet(localFeats, retRealCent, false);
    return dict;
}

void* iftBovwSupOneLevelOrderHrchyKMeansDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding)
{
    /* get the dynamic parameters */
    iftDict *dictAux = iftGetDictFromDict("dict_estim", funcParams);
    int nGroups = iftGetLongValFromDict("n_groups", dictAux);
    iftIntMatrix *localFeatsOrderHrchy = iftGetIntMatrixFromDict("local_feats_order_hrchy", dictAux);

    if(isFisherVectorsEncoding)
        iftError("Fisher vectors encoding is currently not supported with k-means", "iftBovwSupOneLevelOrderHrchyKMeansDictEstimator");

    /* count the number of classes present in the one-level class-hierarchy */
    int nClasses = iftCountUniqueIntElems(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples);

    if(nGroups < nClasses || nGroups <= 0)
        iftError(
                "Number of groups should be greater than 0 and greater or equal than the number of one-level hierarchy classes.",
                "iftBovwSupOneLevelOrderHrchyKMeansDictEstimator");

    if(nClasses == 0)
        iftError("Number of pixel-level classes should be greater than 0.",
                 "iftBovwSupOneLevelOrderHrchyKMeansDictEstimator");

    /* estimate the number of groups by the one-level class-hierarchy according to their size */
    iftIntArray *groupsPerClass = iftBovwComputeNumbGroupsPerOrderHrchy(localFeatsOrderHrchy, localFeats->nsamples, nGroups);

    /* separate the dataset by the one-level class-hierarchy and run k-means */
    iftDataSet** localFeatsPerClass = (iftDataSet**) iftAlloc(nClasses, sizeof(iftDataSet*));
    iftDataSet** dictPerClass = (iftDataSet**) iftAlloc(nClasses, sizeof(iftDataSet*));
    for(int c = 1; c < nClasses+1; c++) {
        localFeatsPerClass[c-1] = NULL;
        dictPerClass[c-1] = NULL;
    }

    #pragma omp parallel for
    for(int c = 1; c < nClasses+1; c++) {
        if(iftFindIntArrayElementIndex(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples, c) != IFT_NIL) {
            int chosenOrdering[] = {c};
            localFeatsPerClass[c-1] = iftBovwSampleDataSetByOrderHrchy(localFeats, localFeatsOrderHrchy, chosenOrdering);

            // printf("Running k-Means for class %i with k=%d (nSamples = %d, nFeats = %d) ... Begin\n", c, groupsPerClass->val[c-1], localFeatsPerClass[c-1]->nsamples, localFeatsPerClass[c-1]->nfeats); 
            iftDict *dictAux = iftCopyDict(funcParams);
            iftDict *dictAux1 = iftGetDictFromDict("dict_estim", dictAux);
            iftUpdateValInDict(iftCreateGVal("n_groups"), iftCreateGVal(groupsPerClass->val[c-1]), dictAux1);
            dictPerClass[c-1] = iftBovwUnsupKMeansDictEstimator(localFeatsPerClass[c-1], dictAux, isFisherVectorsEncoding);
            // printf("Running k-Means for class %i with k=%d (nSamples = %d, nFeats = %d) ... End\n", c, groupsPerClass->val[c-1], localFeatsPerClass[c-1]->nsamples, localFeatsPerClass[c-1]->nfeats); 
        }
    }

    iftDataSet* dict = iftMergeDataSetArray(dictPerClass, nClasses, true);

    for(int c = 1; c < nClasses+1; c++) {
        iftDestroyDataSet(&localFeatsPerClass[c-1]);
        iftDestroyDataSet(&dictPerClass[c-1]);
    }
    iftFree(localFeatsPerClass);
    iftFree(dictPerClass);
    iftDestroyIntArray(&groupsPerClass);

    return dict;
}

void* iftBovwSupTwoLevelOrderHrchyKMeansDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding)
{
    return 0;
}

void* iftBovwUnsupOPFDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("dict_estim", funcParams);
    char *retRealCentStr = iftGetStrValFromDict("ret_real_cent", dictAux);
    double clustSampPerc = iftGetDblValFromDict("clust_samp_perc", dictAux);
    float knnGraphNeighbNumb = iftGetDblValFromDict("knn_graph_n_neighb", dictAux);
    char *knownNumbGroupsStr = iftGetStrValFromDict("known_n_groups", dictAux);
    int nGroups = iftGetLongValFromDict("n_groups", dictAux);
    char *useProtAsCentStr = iftGetStrValFromDict("use_prot_as_cent", dictAux);
    char *findOptmKnnForPdfStr = iftGetStrValFromDict("find_optm_knn_for_pdf", dictAux);

    bool retRealCent = iftBoolStrToBool(retRealCentStr);
    bool knownNumbGroups = iftBoolStrToBool(knownNumbGroupsStr);
    bool useProtAsCent = iftBoolStrToBool(useProtAsCentStr);
    bool findOptmKnnForPdf = iftBoolStrToBool(findOptmKnnForPdfStr);

    if(knownNumbGroups == 1 && nGroups <= 0)
        iftError("Number of groups should be greater than 0.", "iftBovwUnsupOPFDictEstimator");

    /* perform random sampling and then execute OPF clustering with the training samples */
    iftRandomSeed(IFT_RANDOM_SEED);
    int numbTrainSamp = clustSampPerc * (float)localFeats->nsamples;
    iftSampler *sampler = iftRandomSubsampling(localFeats->nsamples, 1, numbTrainSamp);
    iftDataSetSampling(localFeats, sampler, 0);

    iftKnnGraph *graph = NULL;

    if(knnGraphNeighbNumb > 1)
        graph = iftCreateKnnGraph(localFeats, iftMin(knnGraphNeighbNumb, (float)numbTrainSamp*0.3)); // knnGraphNeighbNumb is used as the number of neighbors
    else
        graph = iftCreateKnnGraph(localFeats, (float)numbTrainSamp*knnGraphNeighbNumb); // knnGraphNeighbNumb is used as a percentage of the training samples

    if(knownNumbGroups)
        iftFastUnsupTrainWithCClustersBinarySearch(graph, nGroups);
    else {
        if(findOptmKnnForPdf) {
            iftFastUnsupTrain(graph, iftNormalizedCut);
        }
        else {
            iftPDFByRange(graph);
            iftUnsupOPF(graph);
        }
    }

    /* classify the test samples */
    iftUnsupClassify(graph, localFeats);

    if(isFisherVectorsEncoding) {
        /* compute the cluster model params and create the OPF clust model*/
        iftMatrix *means = NULL, **covariances = NULL, **covariancesInv = NULL;
        iftFloatArray *covariancesDet = NULL, *weights = NULL;
        iftBovwComputeClusterModelParamsFV(graph->Z, &means, &covariances, &covariancesInv, &covariancesDet, &weights);
        iftBovwOPFClustModelFV *clustModel = iftCreateBovwOPFClustModelFV(graph, means, covariances, covariancesInv, covariancesDet, weights);
        iftDestroySampler(&sampler);

        return clustModel;
    }
    else {
        /* extract the centroids from the dataset */
        iftDataSet *dict = iftExtractCentroidsFromDataSetAsDataSet(localFeats, retRealCent, useProtAsCent);
        iftDestroyKnnGraph(&graph);
        iftDestroySampler(&sampler);

        return dict;
    }
}

void* iftBovwSupOneLevelOrderHrchyOPFDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("dict_estim", funcParams);
    // double clustSampPerc = iftGetDblValFromDict("clust_samp_perc", dictAux);
    // float knnGraphNeighbNumb = iftGetDblValFromDict("knn_graph_n_neighb", dictAux);
    char *knownNumbGroupsStr = iftGetStrValFromDict("known_n_groups", dictAux);
    int nGroups = iftGetLongValFromDict("n_groups", dictAux);
    iftIntMatrix *localFeatsOrderHrchy = iftGetIntMatrixFromDict("local_feats_order_hrchy", dictAux);
    bool knownNumbGroups = iftBoolStrToBool(knownNumbGroupsStr);

    /* count the number of classes present in the one-level class-hierarchy */
    int nClasses = iftCountUniqueIntElems(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples);

    if(knownNumbGroups && (nGroups < nClasses || nGroups <= 0))
        iftError("Number of groups should be greater than 0 and greater or equal than the number of hierarchy classes.",
                "iftBovwSupOneLevelOrderHrchyOPFDictEstimator");

    if(nClasses == 0)
        iftError("Number of one-level classes should be greater than 0.",
                 "iftBovwSupOneLevelOrderHrchyOPFDictEstimator");

    /* create the dataset for the local feats per class and the arrays for the dictionaries/clustering models per class */
    iftDataSet** dictPerClass = NULL;
    iftBovwOPFClustModelFV** clustModelPerClass = NULL;
    if(isFisherVectorsEncoding) clustModelPerClass = (iftBovwOPFClustModelFV**)iftAlloc(nClasses, sizeof(iftBovwOPFClustModelFV*));
    else dictPerClass = (iftDataSet**)iftAlloc(nClasses, sizeof(iftDataSet*));

    iftDataSet** localFeatsPerClass = (iftDataSet**) iftAlloc(nClasses, sizeof(iftDataSet*));
    for(int c = 1; c < nClasses+1; c++)
        localFeatsPerClass[c-1] = NULL;

    if(!knownNumbGroups) {
        /* separate the dataset by the one-level class-hierarchy and run OPF */
        #pragma omp parallel for
        for(int c = 1; c < nClasses+1; c++) {
            if(iftFindIntArrayElementIndex(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples, c) != IFT_NIL) {
                int chosenOrdering[] = {c};
                localFeatsPerClass[c-1] = iftBovwSampleDataSetByOrderHrchy(localFeats, localFeatsOrderHrchy, chosenOrdering);
                // printf("Running OPF for class %i (nSamples = %d, nFeats = %d, Sampling = %d, knn_graph = %d) ... Begin\n", c, localFeatsPerClass[c-1]->nsamples, localFeatsPerClass[c-1]->nfeats, (int)(clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples), (int)(knnGraphNeighbNumb*clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples)); 
                if(isFisherVectorsEncoding) clustModelPerClass[c-1] = (iftBovwOPFClustModelFV*)iftBovwUnsupOPFDictEstimator(localFeatsPerClass[c-1], funcParams, isFisherVectorsEncoding);
                else dictPerClass[c-1] = (iftDataSet*)iftBovwUnsupOPFDictEstimator(localFeatsPerClass[c-1], funcParams, isFisherVectorsEncoding);
                // printf("Running OPF for class %i (nSamples = %d, nFeats = %d, Sampling = %d, knn_graph = %d) ... End -> nGroups = %d\n", c, localFeatsPerClass[c-1]->nsamples, localFeatsPerClass[c-1]->nfeats, (int)(clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples), (int)(knnGraphNeighbNumb*clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples), dictPerClass[c-1]->nsamples); 
            }
        }
    }
    else {
        /* estimate the number of groups by the one-level class-hierarchy according to their size */
        iftIntArray *groupsPerClass = iftBovwComputeNumbGroupsPerOrderHrchy(localFeatsOrderHrchy, localFeats->nsamples, nGroups);

        /* separate the dataset by the first-level class-hierarchy and run OPF (with C clusters) */
        #pragma omp parallel for
        for(int c = 1; c < nClasses+1; c++) {
            if(iftFindIntArrayElementIndex(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples, c) != IFT_NIL) {
                int chosenOrdering[] = {c};
                localFeatsPerClass[c-1] = iftBovwSampleDataSetByOrderHrchy(localFeats, localFeatsOrderHrchy, chosenOrdering);
                // printf("Running OPF for class %i (nSamples = %d, nFeats = %d, Sampling = %d, knn_graph = %d) ... Begin -> nGroups = %d\n", c, localFeatsPerClass[c-1]->nsamples, localFeatsPerClass[c-1]->nfeats, (int)(clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples), (int)(knnGraphNeighbNumb*clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples), groupsPerClass->val[c-1]); 
                iftDict *dictAux = iftCopyDict(funcParams);
                iftDict *dictAux1 = iftGetDictFromDict("dict_estim", dictAux);
                iftUpdateValInDict(iftCreateGVal("n_groups"), iftCreateGVal(groupsPerClass->val[c-1]), dictAux1);
                if(isFisherVectorsEncoding) clustModelPerClass[c-1] = (iftBovwOPFClustModelFV*)iftBovwUnsupOPFDictEstimator(localFeatsPerClass[c-1], dictAux, isFisherVectorsEncoding);
                else dictPerClass[c-1] = (iftDataSet*)iftBovwUnsupOPFDictEstimator(localFeatsPerClass[c-1], dictAux, isFisherVectorsEncoding);
                // printf("Running OPF for class %i (nSamples = %d, nFeats = %d, Sampling = %d, knn_graph = %d) ... End\n", c, localFeatsPerClass[c-1]->nsamples, localFeatsPerClass[c-1]->nfeats, (int)(clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples), (int)(knnGraphNeighbNumb*clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples)); 
                iftDestroyDict(&dictAux);
            }
        }
        iftDestroyIntArray(&groupsPerClass);
    }

    /* return the array of cluster models or the dictionaries per class */
    if(isFisherVectorsEncoding) {
        for(int c = 1; c < nClasses+1; c++) {
            if(localFeatsPerClass[c-1] != NULL)
                iftDestroyDataSet(&localFeatsPerClass[c-1]);
        }
        iftFree(localFeatsPerClass);

        return clustModelPerClass;
    }
    else {
        iftDataSet* dict = iftMergeDataSetArray(dictPerClass, nClasses, true);

        for(int c = 1; c < nClasses+1; c++) {
            if(localFeatsPerClass[c-1] != NULL)
                iftDestroyDataSet(&localFeatsPerClass[c-1]);
            if(dictPerClass[c-1] != NULL)
                iftDestroyDataSet(&dictPerClass[c-1]);
        }
        iftFree(localFeatsPerClass);
        iftFree(dictPerClass);

        return dict;
    }
}

void* iftBovwSupTwoLevelOrderHrchyOPFDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("dict_estim", funcParams);
    // double clustSampPerc = iftGetDblValFromDict("clust_samp_perc", dictAux);
    // float knnGraphNeighbNumb = iftGetDblValFromDict("knn_graph_neighb_numb", dictAux);
    char *knownNumbGroupsStr = iftGetStrValFromDict("known_n_groups", dictAux);
    int nGroups = iftGetLongValFromDict("n_groups", dictAux);
    iftIntMatrix *localFeatsOrderHrchy = iftGetIntMatrixFromDict("local_feats_order_hrchy", dictAux);
    bool knownNumbGroups = iftBoolStrToBool(knownNumbGroupsStr);

    /* count the number of classes present in each class-hierarchy */
    int nClasses1 = iftCountUniqueIntElems(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples);
    int nClasses2 = iftCountUniqueIntElems(iftMatrixRowPointer(localFeatsOrderHrchy, 1), localFeats->nsamples);

    if(knownNumbGroups && (nGroups < nClasses1 || nGroups < nClasses2 || nGroups <= 0))
        iftError("Number of groups should be greater than 0 and greater or equal than the number of hierarchy classes.",
                "iftBovwSupTwoLevelOrderHrchyOPFDictEstimator");

    if(nClasses1 == 0)
        iftError("Number of first-level classes should be greater than 0.",
                 "iftBovwSupTwoLevelOrderHrchyOPFDictEstimator");

    if(nClasses2 == 0)
        iftError("Number of second-level classes should be greater than 0.",
                 "iftBovwSupTwoLevelOrderHrchyOPFDictEstimator");

    /* create the dataset for the local feats per class and the arrays for the dictionaries/clustering models per class */
    iftDataSet** dictPerClass = NULL;
    iftBovwOPFClustModelFV** clustModelPerClass = NULL;
    if(isFisherVectorsEncoding) clustModelPerClass = (iftBovwOPFClustModelFV**)iftAlloc(nClasses1*nClasses2, sizeof(iftBovwOPFClustModelFV*));
    else dictPerClass = (iftDataSet**)iftAlloc(nClasses1*nClasses2, sizeof(iftDataSet*));

    iftDataSet** localFeatsPerClass = (iftDataSet**) iftAlloc(nClasses1*nClasses2, sizeof(iftDataSet*));
    for(int c = 1; c < nClasses1*nClasses2+1; c++)
        localFeatsPerClass[c-1] = NULL;

    if(!knownNumbGroups) {
        /* separate the dataset by the first-level and second-level class-hierarchies and run OPF */
        #pragma omp parallel for
        for(int c1 = 0; c1 < nClasses1; c1++) {
            for(int c2 = 0; c2 < nClasses2; c2++) {
                int col = c1*nClasses2+c2;
                if(iftFindIntArrayElementIndex(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples, c1+1) != IFT_NIL ||
                    iftFindIntArrayElementIndex(iftMatrixRowPointer(localFeatsOrderHrchy, 1), localFeats->nsamples, c2+1) != IFT_NIL) {
                    int chosenOrdering[] = {c1+1, c2+1};
                    localFeatsPerClass[col] = iftBovwSampleDataSetByOrderHrchy(localFeats, localFeatsOrderHrchy, chosenOrdering);
                    // printf("Running OPF for class1 %d and class2 %d (nSamples = %d, Sampling = %d, knn_graph = %d) ... Begin\n", c1, c2, localFeatsPerClass[col]->nsamples, (int)(clustSampPerc*(float)localFeatsPerClass[col]->nsamples), (int)(knnGraphNeighbNumb*clustSampPerc*(float)localFeatsPerClass[col]->nsamples));
                    if(isFisherVectorsEncoding) clustModelPerClass[col] = (iftBovwOPFClustModelFV*)iftBovwUnsupOPFDictEstimator(localFeatsPerClass[col], funcParams, isFisherVectorsEncoding);
                    else dictPerClass[col] = (iftDataSet*)iftBovwUnsupOPFDictEstimator(localFeatsPerClass[col], funcParams, isFisherVectorsEncoding);
                    // printf("Running OPF for class1 %d and class2 %d (nSamples = %d, Sampling = %d, knn_graph = %d) ... End -> nGroups = %d\n", c1, c2, localFeatsPerClass[col]->nsamples, (int)(clustSampPerc*(float)localFeatsPerClass[col]->nsamples), (int)(knnGraphNeighbNumb*clustSampPerc*(float)localFeatsPerClass[col]->nsamples), dictPerClass[col]->nrows);
                }
                else {
                    if(!isFisherVectorsEncoding) dictPerClass[col] = iftCreateDataSet(1, localFeats->nfeats);
                }
            }
        }
    }
    else {
        /* estimate the number of groups by the one-level class-hierarchy according to their size */
        iftIntArray *groupsPerClass = iftBovwComputeNumbGroupsPerOrderHrchy(localFeatsOrderHrchy, localFeats->nsamples, nGroups);

        /* separate the dataset by the one-level class-hierarchy and run OPF (with C clusters) */
        #pragma omp parallel for
        for(int c1 = 1; c1 < nClasses1+1; c1++) {
            for(int c2 = 1; c2 < nClasses2+1; c2++) {
                if(iftFindIntArrayElementIndex(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples, c1) != IFT_NIL ||
                    iftFindIntArrayElementIndex(iftMatrixRowPointer(localFeatsOrderHrchy, 1), localFeats->nsamples, c2) != IFT_NIL) {
                    int chosenOrdering[] = {c1, c2};
                    int col = c1*nClasses2+c2-1;
                    localFeatsPerClass[col] = iftBovwSampleDataSetByOrderHrchy(localFeats, localFeatsOrderHrchy, chosenOrdering);
                    // printf("Running OPF for class %i (nSamples = %d, nFeats = %d, Sampling = %d, knn_graph = %d) ... Begin -> nGroups = %d\n", c, localFeatsPerClass[c-1]->nsamples, localFeatsPerClass[c-1]->nfeats, (int)(clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples), (int)(knnGraphNeighbNumb*clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples), groupsPerClass->val[c]); 
                    iftDict *dictAux = iftCopyDict(funcParams);
                    iftDict *dictAux1 = iftGetDictFromDict("dict_estim", dictAux);
                    iftUpdateValInDict(iftCreateGVal("n_groups"), iftCreateGVal(groupsPerClass->val[col]), dictAux1);
                    if(isFisherVectorsEncoding) clustModelPerClass[col] = (iftBovwOPFClustModelFV*)iftBovwUnsupOPFDictEstimator(localFeatsPerClass[col], dictAux, isFisherVectorsEncoding);
                    else dictPerClass[col] = (iftDataSet*)iftBovwUnsupOPFDictEstimator(localFeatsPerClass[col], dictAux, isFisherVectorsEncoding);
                    // printf("Running OPF for class %i (nSamples = %d, nFeats = %d, Sampling = %d, knn_graph = %d) ... End\n", c, localFeatsPerClass[c-1]->nsamples, localFeatsPerClass[c-1]->nfeats, (int)(clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples), (int)(knnGraphNeighbNumb*clustSampPerc*(float)localFeatsPerClass[c-1]->nsamples)); 
                    iftDestroyDict(&dictAux);
                }
            }
        }
        iftDestroyIntArray(&groupsPerClass);
    }

    /* return the array of cluster models or the dictionaries per class */
    if(isFisherVectorsEncoding) {
        for(int c = 1; c < nClasses1*nClasses2+1; c++) {
            if(localFeatsPerClass[c-1] != NULL)
                iftDestroyDataSet(&localFeatsPerClass[c-1]);
        }
        iftFree(localFeatsPerClass);

        return clustModelPerClass;
    }
    else {
        iftDataSet* dict = iftMergeDataSetArray(dictPerClass, nClasses1*nClasses2, true);

        for(int c = 1; c < nClasses1*nClasses2+1; c++) {
            if(localFeatsPerClass[c-1] != NULL)
                iftDestroyDataSet(&localFeatsPerClass[c-1]);
            if(dictPerClass[c-1] != NULL)
                iftDestroyDataSet(&dictPerClass[c-1]);
        }
        iftFree(localFeatsPerClass);
        iftFree(dictPerClass);

        return dict;
    }
}

void* iftBovwSupOneLevelOrderHrchyManualDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("dict_estim", funcParams);
    char *selCentJsonFilename = iftGetStrValFromDict("sel_cent_json_filename", dictAux);
    iftIntMatrix *localFeatsOrderHrchy  = iftGetIntMatrixFromDict("local_feats_order_hrchy", dictAux);

    if(isFisherVectorsEncoding)
        iftError("Fisher vectors encoding is currently not supported with manual dict estimation", "iftBovwSupOneLevelOrderHrchyManualDictEstimator");

    /* count the number of classes present in the one-level class-hierarchy */
    int nClasses = iftCountUniqueIntElems(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples);

    /* read the json file */
    iftDict *selCentDict = iftReadJson(selCentJsonFilename);
    char *hrchyType = iftGetStrValFromDict("hierarchy_type", selCentDict);
    iftDataSet* dict = NULL;

    /* if the hierarchy type is per position */
    if(iftCompareStrings(hrchyType, "per_position")) {
        int nPositions = iftGetLongValFromDict("n_positions", selCentDict);
        
        /* validate the number of positions in the JSON and the ones present in the local feats hierarchy */
        if(nPositions != nClasses)
            iftError("The number of positions in the JSON file and the positions extracted from the local feat vectors does not match",
                       "iftBovwSupOneLevelOrderHrchyManualDictEstimator" );
        
        /* create an array of datasets for the centroids (visual dictionary) one for each position */
        iftDataSet** dictPerPos = (iftDataSet**) iftAlloc(nPositions, sizeof(iftDataSet*));

        /* read the information (dictionary) for every position */
        for(int p = 0; p < nPositions; p++) {
            char dictName[2048];
            sprintf(dictName, "position_%d", p+1);
            iftDict *posDict = iftGetDictFromDict(dictName, selCentDict);
            // int x = iftGetLongValFromDict("x", posDict);
            // int y = iftGetLongValFromDict("y", posDict);
            char *datasetName = iftGetStrValFromDict("dataset", posDict);
            iftIntArray *selSamples = iftGetIntArrayFromDict("selected_samples", posDict);
            iftDataSet *dataset = iftReadDataSet(datasetName);

            /* get the centroids from the dataset */
            dictPerPos[p] = iftCreateDataSet(selSamples->n, localFeats->nfeats);
            for(int s = 0; s < selSamples->n; s++)
                for(int w = 0; w < localFeats->nfeats; w++)
                    dictPerPos[p]->sample[s].feat[w] = dataset->sample[selSamples->val[s]].feat[w];
            
            iftDestroyDataSet(&dataset);
        }

        dict = iftMergeDataSetArray(dictPerPos, nPositions, true);

        for(int p = 0; p < nPositions; p++)
            if(dictPerPos[p] != NULL)
                iftDestroyDataSet(&dictPerPos[p]);
        iftFree(dictPerPos);
    }
    /* if the hierarchy type is per img-class */
    else if(iftCompareStrings(hrchyType, "per_image_class")) {
        int nClassesJson = iftGetLongValFromDict("n_classes", selCentDict);
        
        /* validate the number of classes in the JSON and the ones present in the local feats hierarchy */
        if(nClassesJson != nClasses)
            iftError("The number of classes in the JSON file and the classes extracted from the local feat vectors does not match",
                       "iftBovwSupOneLevelOrderHrchyManualDictEstimator" );
        
        /* create an array of datasets for the centroids (visual dictionary) one for each class */
        iftDataSet** dictPerClass = (iftDataSet**) iftAlloc(nClasses, sizeof(iftDataSet*));

        /* read the information (dictionary) for every class */
        for(int c = 0; c < nClasses; c++) {
            char dictName[2048];
            sprintf(dictName, "class_%d", c+1);
            iftDict *classDict = iftGetDictFromDict(dictName, selCentDict);
            char *datasetName = iftGetStrValFromDict("dataset", classDict);
            iftIntArray *selSamples = iftGetIntArrayFromDict("selected_samples", classDict);
            iftDataSet *dataset = iftReadDataSet(datasetName);

            /* get the centroids from the dataset */
            dictPerClass[c] = iftCreateDataSet(selSamples->n, localFeats->nfeats);
            for(int s = 0; s < selSamples->n; s++)
                for(int w = 0; w < localFeats->nfeats; w++)
                    dictPerClass[c]->sample[s].feat[w] = dataset->sample[selSamples->val[s]].feat[w];
            
            iftDestroyDataSet(&dataset);
        }

        dict = iftMergeDataSetArray(dictPerClass, nClasses, true);

        for(int c = 0; c < nClasses; c++)
            if(dictPerClass[c] != NULL)
                iftDestroyDataSet(&dictPerClass[c]);
        iftFree(dictPerClass);
    }

    iftDestroyDict(&selCentDict);

    return dict;
}

void* iftBovwSupTwoLevelOrderHrchyManualDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("dict_estim", funcParams);
    char *selCentJsonFilename = iftGetStrValFromDict("sel_cent_json_filename", dictAux);
    iftIntMatrix *localFeatsOrderHrchy  = iftGetIntMatrixFromDict("local_feats_order_hrchy", dictAux);

    if(isFisherVectorsEncoding)
        iftError("Fisher vectors encoding is currently not supported with manual dict estimation", "iftBovwSupTwoLevelOrderHrchyManualDictEstimator");

    /* count the number of classes present in the two-level class-hierarchy */
    int nClasses1 = iftCountUniqueIntElems(iftMatrixRowPointer(localFeatsOrderHrchy, 0), localFeats->nsamples);
    int nClasses2 = iftCountUniqueIntElems(iftMatrixRowPointer(localFeatsOrderHrchy, 1), localFeats->nsamples);

    /* read the json file */
    iftDict *selCentDict = iftReadJson(selCentJsonFilename);
    char *hrchyType = iftGetStrValFromDict("hierarchy_type", selCentDict);
    iftDataSet* dict = NULL;

    /* if the hierarchy type is per img-class and position */
    if(iftCompareStrings(hrchyType, "per_image_class_and_position")) {
        int nClasses = iftGetLongValFromDict("n_classes", selCentDict);
        int nPositions = iftGetLongValFromDict("n_positions", selCentDict);
        
        /* validate the number of classes and positions in the JSON and the ones present in the local feats hierarchy */
        if((nClasses != nClasses1) || (nPositions != nClasses2))
            iftError("The number of image classes and positions in the JSON file and the ones extracted from the local feat vectors does not match",
                       "iftBovwSupTwoLevelOrderHrchyManualDictEstimator" );
        
        /* create an array of datasets for the centroids (visual dictionary) one for each position */
        iftDataSet** dictPerPos = (iftDataSet**) iftAlloc(nClasses*nPositions, sizeof(iftDataSet*));

        /* read the information (dictionary) for every position */
        for(int c = 0; c < nClasses; c++) {
            for(int p = 0; p < nPositions; p++) {
                char dictName[2048];
                sprintf(dictName, "class_%d_position_%d", c+1, p+1);
                iftDict *posDict = iftGetDictFromDict(dictName, selCentDict);
                // int x = iftGetLongValFromDict("x", posDict);
                // int y = iftGetLongValFromDict("y", posDict);
                char *datasetName = iftGetStrValFromDict("dataset", posDict);
                iftIntArray *selSamples = iftGetIntArrayFromDict("selected_samples", posDict);
                iftDataSet *dataset = iftReadDataSet(datasetName);

                /* get the centroids from the dataset */
                int col = c*nPositions+p;
                dictPerPos[col] = iftCreateDataSet(selSamples->n, localFeats->nfeats);
                for(int s = 0; s < selSamples->n; s++)
                    for(int w = 0; w < localFeats->nfeats; w++)
                        dictPerPos[col]->sample[s].feat[w] = dataset->sample[selSamples->val[s]].feat[w];
                
                iftDestroyDataSet(&dataset);
            }
        }

        dict = iftMergeDataSetArray(dictPerPos, nClasses*nPositions, true);

        for(int p = 0; p < nClasses*nPositions; p++)
            if(dictPerPos[p] != NULL)
                iftDestroyDataSet(&dictPerPos[p]);
        iftFree(dictPerPos);
    }

    iftDestroyDict(&selCentDict);

    return dict;
}

iftFeatures* iftBovwHardAsgmtCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams)
{
    iftDataSet *dict1 = (iftDataSet*)dict; // for this specific function we can safely assume that dict is a iftDataSet*
    int dictSize = dict1->nsamples;
    iftFeatures *featVect = iftCreateFeatures(dictSize);

    /* find the closest visual word to each local feat */
    for(int s = 0; s < localFeats->nsamples; s++) {
        float minDist = IFT_INFINITY_FLT;
        int minDistIdx = IFT_NIL;
        for(int r = 0; r < dictSize; r++) {
            float dist = localFeats->iftArcWeight(iftMatrixRowPointer(dict1->data, r), localFeats->sample[s].feat, localFeats->alpha, localFeats->nfeats);
            if(dist < minDist) {
                minDist = dist;
                minDistIdx = r;
            }
        }
        featVect->val[minDistIdx] += 1.0;
    }
    iftNormalizeFeatures(featVect->val, dictSize);

    return featVect;
}

iftFeatures* iftBovwSoftAsgmtCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("cod_func", funcParams);
    int numbNearNeighb = iftGetLongValFromDict("n_near_neighb", dictAux);
    char *wgtFuncStr = iftGetStrValFromDict("wgt_func", dictAux);
    int wgtFunc = iftBovwWgtFuncStrToWgtFunc(wgtFuncStr);

    iftDataSet *dict1 = (iftDataSet*)dict; // for this specific function we can safely assume that dict is a iftDataSet*
    int dictSize = dict1->nsamples;
    iftFeatures *featVect = iftCreateFeatures(dictSize);

    /* compute the distance from each local feat to each visual word using matrix multiplication */
    iftDistanceTable *distTab = iftCompEuclDistanceTable(dict1, localFeats);

    /* compute the closest words to each local feat */
    for(int s = 0; s < localFeats->nsamples; s++) {
        iftFloatArray *dist = iftCreateFloatArray(dictSize);
        iftIntArray *distIdx = iftCreateIntArray(dictSize);
        for(int r = 0; r < dictSize; r++) {
            dist->val[r] = distTab->distance_table[r][s];
            distIdx->val[r] = r;
        }
        
        iftFQuickSort(dist->val, distIdx->val, 0, dictSize-1, IFT_INCREASING);

        /* if we are using gaussian weighted distance, compute mean and variance of the distances */
        float mean = 0.0, var = 0.0;
        if(wgtFunc == BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_GAUSS) {
            mean = iftMean(dist->val, dictSize);
            var = iftVar(dist->val, dictSize);
        }

        /* add the weighted distance to the k nearest visual words */
        for(int k = 0; k < numbNearNeighb; k++) {
            switch(wgtFunc) {
                case BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_INVER:
                    //featVect->val[distIdx->val[k]] += (dist->val[distIdx->val[k]] == 0.0) ? 1.0 : 1.0 / dist->val[distIdx->val[k]];
                    featVect->val[distIdx->val[k]] += (dist->val[k] == 0.0 ? 1.0 : 1.0 / dist->val[k]);
                    break;
                case BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_GAUSS:
                    //featVect->val[distIdx->val[k]] += (var == 0.0) ? exp(-1.0) : exp(-1.0*(powf(dist->val[distIdx->val[k]]-mean, 2))/(2.0*var));
                    featVect->val[distIdx->val[k]] += (var == 0.0 ? exp(-1.0) : exp(-1.0/2.0*powf((dist->val[k]-mean)/var, 2)));
                    break;
            }
        }
        iftDestroyFloatArray(&dist);
        iftDestroyIntArray(&distIdx);
    }
    iftNormalizeFeatures(featVect->val, dictSize);
    
    iftDestroyDistanceTable(&distTab);

    return featVect;
}

iftFeatures* iftBovwHardAsgmtBatchCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams)
{
    return 0;
}

iftFeatures* iftBovwSoftAsgmtBatchCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("cod_func", funcParams);
    int numbNearNeighb = iftGetLongValFromDict("n_near_neighb", dictAux);
    char *wgtFuncStr = iftGetStrValFromDict("wgt_func", dictAux);
    int nImgsInBatch = iftGetLongValFromDict("n_imgs_in_batch", dictAux);
    int wgtFunc = iftBovwWgtFuncStrToWgtFunc(wgtFuncStr);

    iftDataSet *dict1 = (iftDataSet*)dict; // for this specific function we can safely assume that dict is a iftDataSet*
    int nPatches = localFeats->nsamples/nImgsInBatch;
    int dictSize = dict1->nsamples;
    iftFeatures *featVect = iftCreateFeatures(dictSize*nImgsInBatch);

    /* compute the distance from each local feat to each visual word using matrix multiplication */
    iftDistanceTable *distTab = iftCompEuclDistanceTable(dict1, localFeats);

    /* compute the closest words to the local feats of each image */
    for(int i = 0; i < nImgsInBatch; i++) {
        printf("- computing the closest visual words (image %d/%d) ... \r", i+1, nImgsInBatch); fflush(stdout);
        #pragma omp parallel for
        for(int p = 0; p < nPatches; p++) {
            iftFloatArray *dist = iftCreateFloatArray(dictSize);
            iftIntArray *distIdx = iftCreateIntArray(dictSize);
            int col = i*nPatches + p;
            for(int r = 0; r < dictSize; r++) {
                dist->val[r] = distTab->distance_table[r][col];
                distIdx->val[r] = r;
            }
            
            iftFQuickSort(dist->val, distIdx->val, 0, dictSize-1, IFT_INCREASING);

            /* if we are using gaussian weighted distance, compute mean and variance of the distances */
            float mean = 0.0, var = 0.0;
            if(wgtFunc == BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_GAUSS) {
                mean = iftMean(dist->val, dictSize);
                var = iftVar(dist->val, dictSize);
            }

            /* add the weighted distance to the k nearest visual words */
            for(int k = 0; k < numbNearNeighb; k++) {
                int col1 = i*dictSize + distIdx->val[k];
                switch(wgtFunc) {
                    case BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_INVER:
                        //featVect->val[distIdx->val[k]] += (dist->val[distIdx->val[k]] == 0.0) ? 1.0 : 1.0 / dist->val[distIdx->val[k]];
                        featVect->val[col1] += (dist->val[k] == 0.0 ? 1.0 : 1.0 / dist->val[k]);
                        break;
                    case BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_GAUSS:
                        //featVect->val[distIdx->val[k]] += (var == 0.0) ? exp(-1.0) : exp(-1.0*(powf(dist->val[distIdx->val[k]]-mean, 2))/(2.0*var));
                        featVect->val[col1] += (var == 0.0 ? exp(-1.0) : exp(-1.0/2.0*powf((dist->val[k]-mean)/var, 2)));
                        break;
                }
            }
            iftDestroyFloatArray(&dist);
            iftDestroyIntArray(&distIdx);
        }
        iftNormalizeFeatures(featVect->val + i*dictSize, dictSize);
    }
    printf("- computing the closest visual words (image %d/%d) ... OK\n", nImgsInBatch, nImgsInBatch);

// #ifdef IFT_GPU
//     float usedMemoryGPU =   iftMemoryUsedByArray(localFeats->data->n, sizeof(localFeats->data->val[0]), 'm') + 
//                             iftMemoryUsedByArray(dict1->n, sizeof(dict1->val[0]), 'm') +
//                             iftMemoryUsedByArray(dictDS->nsamples*localFeats->nsamples, sizeof(localFeats->data->val[0]), 'm');
//     float freeMemoryGPU = iftGetFreeMemoryGPU(0)/1024.0/1024.0;
//     printf("- GPU memory: Used: %.2f MB, Free: %.2f MB, Util: %.2f%s\n", usedMemoryGPU, freeMemoryGPU, usedMemoryGPU/freeMemoryGPU*100.0, "%");
// #endif

    iftDestroyDistanceTable(&distTab);

    return featVect;
}

iftFeatures* iftBovw2ndOrderStatsFisherVectorsCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams)
{
    /* get the parameters */
    iftDict *dictAux = iftGetDictFromDict("cod_func", funcParams);
    iftBovwDictEstimatorId dictEstimId = iftGetLongValFromDict("dict_estim_id", dictAux);

    /* compute the fisher vectors according to each type of clustering model */
    if(dictEstimId == BOVW_DICT_ESTIMATOR_UNSUP_OPF) {
        iftBovwOPFClustModelFV *clustModel = (iftBovwOPFClustModelFV*)dict;

        iftFeatures *featVect = iftBovwComputeFeatVect2ndOrderStatsFV(localFeats, clustModel->means, clustModel->covariances,
                                                                        clustModel->covariancesInv, clustModel->covariancesDet,
                                                                        clustModel->weights);
        return featVect;
    }
    else if(dictEstimId == BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL ||
        dictEstimId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS ||
        dictEstimId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE ||
        dictEstimId == BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION) {
        iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)dict;

        int nClustModels = iftGetLongValFromDict("n_clust_models", dictAux);
        iftFeatures **featVects = (iftFeatures**)iftAlloc(nClustModels, sizeof(iftFeatures*));

        for(int c = 1; c < nClustModels+1; c++)
            featVects[c-1] = iftBovwComputeFeatVect2ndOrderStatsFV(localFeats, clustModels[c]->means, clustModels[c]->covariances,
                                                                    clustModels[c]->covariancesInv, clustModels[c]->covariancesDet,
                                                                    clustModels[c]->weights);

        iftFeatures *featVect = iftJoinFeatures(featVects, nClustModels);
        return featVect;
    }
    else if(dictEstimId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS_AND_POSITION) {
        iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)dict;

        int nClustModels1 = iftGetLongValFromDict("n_clust_models_1", dictAux);
        int nClustModels2 = iftGetLongValFromDict("n_clust_models_2", dictAux);
        iftFeatures **featVects = (iftFeatures**)iftAlloc(nClustModels1*nClustModels2, sizeof(iftFeatures*));

        for(int c1 = 0; c1 < nClustModels1; c1++) {
            for(int c2 = 0; c2 < nClustModels2; c2++) {
                int col = c1*nClustModels2+c2;
                featVects[col] = iftBovwComputeFeatVect2ndOrderStatsFV(localFeats, clustModels[col]->means, clustModels[col]->covariances,
                                                                        clustModels[col]->covariancesInv, clustModels[col]->covariancesDet,
                                                                        clustModels[col]->weights);
            }
        }

        iftFeatures *featVect = iftJoinFeatures(featVects, nClustModels1*nClustModels2);
        return featVect;
    }
    else {
        iftError("Fisher vectors encoding is currently not supported with the selected dictionary estimation method", "iftBovw2ndOrderStatsFisherVectorsCodFunc");
        return 0;
    }
}

/* ---------------------------------------------------------------------------------------------------------------------*/
/*                                                Complementary functions                                               */
/* ---------------------------------------------------------------------------------------------------------------------*/

bool iftBovwIsFisherVectorsEncoding(iftBovwCodFuncId codFuncId)
{
    if(codFuncId == BOVW_COD_FUNC_2ND_ORDER_STATS_FISHER_VECTORS)
        return true;
    return false;
}

int iftBovwGetNumbVisualWords(iftBagOfVisualWords *bovw)
{
    int nVisWords = 0;
    if(iftBovwIsFisherVectorsEncoding(bovw->codFuncId)) {
        if(bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_UNSUP_OPF) {
            iftBovwOPFClustModelFV *clustModel = (iftBovwOPFClustModelFV*)(bovw->dict);
            nVisWords = clustModel->weights->n;
        }
        else if( bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION) {
            
            iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)(bovw->dict);
            int nClustModels = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 0), bovw->localFeats->nsamples);
            for(int c = 1; c < nClustModels+1; c++)
                nVisWords += clustModels[c]->weights->n;
        }
        else if( bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS_AND_POSITION) {            
            iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)(bovw->dict);
            int nClustModels1 = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 0), bovw->localFeats->nsamples);
            int nClustModels2 = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 1), bovw->localFeats->nsamples);
            for(int c1 = 0; c1 < nClustModels1; c1++) {
                for(int c2 = 0; c2 < nClustModels2; c2++) {
                    int col = c1*nClustModels2+c2;
                    nVisWords += clustModels[col]->weights->n;
                }
            }
        }
    }
    else {
        iftDataSet *dict = bovw->dict;
        nVisWords = dict->nsamples;
    }

    return nVisWords;
}

int iftBovwGetNumbFeatsPerWord(iftBagOfVisualWords *bovw)
{
    int nFeatsPerWord = 0;
    if(iftBovwIsFisherVectorsEncoding(bovw->codFuncId)) {
        if(bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_UNSUP_OPF) {
            iftBovwOPFClustModelFV *clustModel = (iftBovwOPFClustModelFV*)(bovw->dict);
            nFeatsPerWord = clustModel->nFeatsPerWord;
        }
        else if( bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION) {
            
            iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)(bovw->dict);
            int nClustModels = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 0), bovw->localFeats->nsamples);
            for(int c = 1; c < nClustModels+1; c++)
                nFeatsPerWord += clustModels[c]->nFeatsPerWord;
        }
        else if( bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS_AND_POSITION) {            
            iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)(bovw->dict);
            int nClustModels1 = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 0), bovw->localFeats->nsamples);
            int nClustModels2 = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 1), bovw->localFeats->nsamples);
            for(int c1 = 0; c1 < nClustModels1; c1++) {
                for(int c2 = 0; c2 < nClustModels2; c2++) {
                    int col = c1*nClustModels2+c2;
                    nFeatsPerWord += clustModels[col]->nFeatsPerWord;
                }
            }
        }
    }
    else {
        iftDataSet *dict = bovw->dict;
        nFeatsPerWord = dict->nfeats;
    }

    return nFeatsPerWord;
}

int iftBovwComputeNumbFeatsForCodFunc(iftBagOfVisualWords *bovw)
{
    int nGlobalFeats = 0;
    int nLocalFeats = iftBovwGetNumbFeatsPerWord(bovw);

    if(iftBovwIsFisherVectorsEncoding(bovw->codFuncId)) {
        if(bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_UNSUP_OPF) {
            iftBovwOPFClustModelFV *clustModel = (iftBovwOPFClustModelFV*)(bovw->dict);
            int nGroups = clustModel->weights->n;
            nGlobalFeats += nGroups*nLocalFeats*3;
        }
        else if( bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE ||
            bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION) {
            
            iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)(bovw->dict);
            int nClustModels = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 0), bovw->localFeats->nsamples);
            for(int c = 1; c < nClustModels+1; c++) {
                int nGroups = clustModels[c]->weights->n;
                nGlobalFeats += nGroups*nLocalFeats*3;
            }
        }
        else if( bovw->dictEstimatorId == BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS_AND_POSITION) {            
            iftBovwOPFClustModelFV **clustModels = (iftBovwOPFClustModelFV**)(bovw->dict);
            int nClustModels1 = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 0), bovw->localFeats->nsamples);
            int nClustModels2 = iftCountUniqueIntElems(iftMatrixRowPointer(bovw->localFeatsOrderHrchy, 1), bovw->localFeats->nsamples);
            for(int c1 = 0; c1 < nClustModels1; c1++) {
                for(int c2 = 0; c2 < nClustModels2; c2++) {
                    int col = c1*nClustModels2+c2;
                    int nGroups = clustModels[col]->weights->n;
                    nGlobalFeats += nGroups*nLocalFeats*3;
                }
            }
        }
    }
    else {
        iftDataSet *dict = bovw->dict;
        nGlobalFeats = dict->nsamples;
    }

    return nGlobalFeats;
}

iftFileSet *iftBoVWCreateFileSetRefDataFromImgROIs(iftImage *img, iftRoiArray *roiArray, iftFile *imgFilePtr, char *patchDirPath)
{
    iftFileSet *fileSet = iftCreateFileSet(roiArray->n);

    for(int i = 0; i < roiArray->n; i++) {
        /* create the patch image */
        iftImage *patch = iftExtractRoiNoBkgd(img, roiArray->val[i]);

        /* update the file info and write the image */
        char name[2048]; sprintf(name, "%s_%06d.png", iftFilename(imgFilePtr->path, iftFileExt(imgFilePtr->path)), i+1);
        char *filename = iftJoinPathnames(2, patchDirPath, name);
        iftWriteImageByExt(patch, filename);
        char *filenameAbs = iftAbsPathname(filename);
        fileSet->files[i] = iftCreateFile(filenameAbs);
        iftUpdateDataSetFileInfo(fileSet->files[i]);
        iftDestroyImage(&patch);
    }

    return fileSet;
}

iftIntArray *iftBovwComputeNumbGroupsPerOrderHrchy(iftIntMatrix *localFeatsOrderHrchy, int nSamples, int nGroups)
{
    int nOrderHrchies = localFeatsOrderHrchy->nrows;

    if(nOrderHrchies == 1) {
        int *classNumbers = NULL, *samplesPerClass = NULL;
        int nClasses =  iftCountElemsPerUniqueIntElem(iftMatrixRowPointer(localFeatsOrderHrchy, 0), nSamples, &classNumbers, &samplesPerClass);
        iftIntArray *groupsPerClass = iftCreateIntArray(nClasses);
        int sumGroups = 0;

        /* distribute the number of groups (nGroups) according to the size of each class */
        for (int c = 0; c < nClasses; c++) {
            groupsPerClass->val[c] = ceil((float)nGroups*(float)samplesPerClass[c]/(float)nSamples);
            sumGroups += groupsPerClass->val[c];
        }

        /* if there are remaining groups (e.g. rounding error) they are added to the first class that is not empty */
        if(sumGroups < nGroups) {
            for (int c = 0; c < nClasses; c++)
                if (groupsPerClass->val[c] > 0) {
                    groupsPerClass->val[c] += (nGroups - sumGroups);
                    break;
                }
        }
        else if(sumGroups > nGroups) { // if there are more groups, this value is substracted from the biggest group
            int maxVal = IFT_INFINITY_INT_NEG, maxValIdx = IFT_NIL;
            for (int c = 0; c < nClasses; c++)
                if (groupsPerClass->val[c] > maxVal) {
                    maxVal = groupsPerClass->val[c];
                    maxValIdx = c;
                }
            groupsPerClass->val[maxValIdx] -= (sumGroups - nGroups);
        }
        iftFree(samplesPerClass);

        return groupsPerClass;
    }

    return 0;
}

iftIntMatrix *iftBovwComputeOrderHrchiesMatrixForImg(iftDataSet *localFeats, iftIntArray *orderHrchies, iftFileSet *dictLearnFileset,
                                                        iftFileSet *dictLearnMasksFileset, int imgId)
{
    int nOrderHrchies = orderHrchies->n;
    iftIntMatrix *localFeatsOrderHrchies = iftCreateIntMatrix(localFeats->nsamples, nOrderHrchies);

    /* obtain the labels for BOVW_LOCAL_FEATS_ORDERING_PIX_LABEL (if it is the case) */
    int hrchyIdPos = iftFindIntArrayElementIndex(orderHrchies->val, nOrderHrchies, BOVW_LOCAL_FEATS_ORDERING_PIX_LABEL);
    if(hrchyIdPos != IFT_NIL) {
        if(dictLearnMasksFileset == NULL) iftError("Learn image masks need to be provided in order to use pixel-level hierarchy", "iftBovwComputeOrderHrchiesMatrixForImg");
        iftCopyIntArray(iftMatrixRowPointer(localFeatsOrderHrchies, hrchyIdPos),
                        iftBovwComputeOrderHrchyArray(localFeats, BOVW_LOCAL_FEATS_ORDERING_PIX_LABEL, 1, dictLearnMasksFileset->files[imgId]->path),
                        localFeats->nsamples);
    }

    /* obtain the labels for BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS (if it is the case) */
    hrchyIdPos = iftFindIntArrayElementIndex(orderHrchies->val, nOrderHrchies, BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS);
    if(hrchyIdPos != IFT_NIL) {
        iftCopyIntArray(iftMatrixRowPointer(localFeatsOrderHrchies, hrchyIdPos),
                        iftBovwComputeOrderHrchyArray(localFeats, BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS, 1, dictLearnFileset->files[imgId]->label),
                        localFeats->nsamples);
    }

    /* obtain the labels for BOVW_LOCAL_FEATS_ORDERING_IMAGE (if it is the case) */
    hrchyIdPos = iftFindIntArrayElementIndex(orderHrchies->val, nOrderHrchies, BOVW_LOCAL_FEATS_ORDERING_IMAGE);
    if(hrchyIdPos != IFT_NIL) {
        iftCopyIntArray(iftMatrixRowPointer(localFeatsOrderHrchies, hrchyIdPos),
                        iftBovwComputeOrderHrchyArray(localFeats, BOVW_LOCAL_FEATS_ORDERING_IMAGE, 1, imgId+1), /* to avoid imgId=0 */
                        localFeats->nsamples);
    }

    /* obtain the labels for BOVW_LOCAL_FEATS_ORDERING_POSITION (if it is the case) */
    hrchyIdPos = iftFindIntArrayElementIndex(orderHrchies->val, nOrderHrchies, BOVW_LOCAL_FEATS_ORDERING_POSITION);
    if(hrchyIdPos != IFT_NIL) {
        iftCopyIntArray(iftMatrixRowPointer(localFeatsOrderHrchies, hrchyIdPos),
                        iftBovwComputeOrderHrchyArray(localFeats, BOVW_LOCAL_FEATS_ORDERING_POSITION, 0),
                        localFeats->nsamples);
    }

    return localFeatsOrderHrchies;
}

int *iftBovwComputeOrderHrchyArray(iftDataSet *feats, int hrchyId, int nParams, ...)
{
    int *orderHrchy = iftAllocIntArray(feats->nsamples);
    va_list args;

    if(hrchyId == BOVW_LOCAL_FEATS_ORDERING_PIX_LABEL) {
        /* read dynamic parameters */
        va_start(args, nParams);
        char *imgMaskPath = va_arg(args, char * );
        va_end(args);

        iftImage *label = iftReadImageByExt(imgMaskPath);
        for (int p = 0; p < feats->nsamples; p++)
            orderHrchy[p] = label->val[feats->sample[p].id] > 127 ? 1 : 0;
    }

    if(hrchyId == BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS) {
        /* read dynamic parameters */
        va_start(args, nParams);
        int imgClass = va_arg(args, int);
        va_end(args);

        for (int p = 0; p < feats->nsamples; p++)
            orderHrchy[p] = imgClass;
    }

    if(hrchyId == BOVW_LOCAL_FEATS_ORDERING_IMAGE) {
        /* read dynamic parameters */
        va_start(args, nParams);
        int imgId = va_arg(args, int);
        va_end(args);

        for (int p = 0; p < feats->nsamples; p++)
            orderHrchy[p] = imgId;
    }

    if(hrchyId == BOVW_LOCAL_FEATS_ORDERING_POSITION) {
        /* there are no parameters for this hierarchy */

        for (int p = 0; p < feats->nsamples; p++)
            orderHrchy[p] = p+1; /* positions starts from 1 */
    }

    return orderHrchy;
}

iftDataSet *iftBovwSampleDataSetByOrderHrchy(iftDataSet *Z, iftIntMatrix *localFeatsOrderHrchy, int *chosenOrdering)
{
    /* count the number of samples for the given class combination (according to the number of hierarchies */
    int classSize=0;
    int nOrderHrchies = localFeatsOrderHrchy->nrows;
    switch(nOrderHrchies) {
        case 1:
            /* for one hierarchy, the chosenOrdering must be checked at one level */
            for (int s = 0; s < Z->nsamples; s++)
                if (iftMatrixRowPointer(localFeatsOrderHrchy,0)[s] == chosenOrdering[0]) classSize++;
            if (classSize == 0)
                iftError("There are no samples from the chosen ordering: %d",
                         "iftBovwSampleDataSetByOrderHrchy", chosenOrdering[0]);
            break;

        case 2:
            /* for two hierarchies, the chosenOrdering must be checked at two levels */
            for (int s = 0; s < Z->nsamples; s++)
                if (    iftMatrixRowPointer(localFeatsOrderHrchy,0)[s] == chosenOrdering[0] &&
                        iftMatrixRowPointer(localFeatsOrderHrchy,1)[s] == chosenOrdering[1]  ) classSize++;
            if (classSize == 0)
                iftError("There are no samples from the chosen ordering combination: %d, %d",
                         "iftBovwSampleDataSetByOrderHrchy", chosenOrdering[0], chosenOrdering[1]);
            break;
    }

    /* create the new dataset */
    iftDataSet *Z1 = iftCreateDataSet(classSize, Z->nfeats);
    Z1->function_number = Z->function_number;
    iftCopyRefData(Z1, Z->ref_data, Z->ref_data_type);
    iftSetDistanceFunction(Z1, Z->function_number);
    Z1->nclasses = Z->nclasses;
    Z1->ngroups = Z->ngroups;

    /* copy the samples that match the given class combination (according to the number of hierarquies) */
    int t = 0;
    switch(nOrderHrchies) {
        case 1:
            /* for one hierarchy, the chosenOrdering must be checked at one level */
            for (int s = 0; s < Z->nsamples; s++)
                if (iftMatrixRowPointer(localFeatsOrderHrchy,0)[s] == chosenOrdering[0]) {
                    iftCopySample(&Z->sample[s], &Z1->sample[t], Z->nfeats, true);
                    t++;
                }
            break;

        case 2:
            /* for two hierarchies, the chosenOrdering must be checked at two levels */
            for (int s = 0; s < Z->nsamples; s++)
                if (    iftMatrixRowPointer(localFeatsOrderHrchy,0)[s] == chosenOrdering[0] &&
                        iftMatrixRowPointer(localFeatsOrderHrchy,1)[s] == chosenOrdering[1]  ) {
                    iftCopySample(&Z->sample[s], &Z1->sample[t], Z->nfeats, true);
                    t++;
                }
            break;
    }

    return Z1;
}

void iftBovwPrintChosenMethods(iftBovwIntPointDetectorId intPointDetec, iftBovwLocalFeatExtractorId localFeatExtr, iftBovwDictEstimatorId dictEstim, iftBovwCodFuncId codFunc)
{
    int numbMethods = ((int)intPointDetec != IFT_NIL) + ((int)localFeatExtr != IFT_NIL) + ((int)dictEstim != IFT_NIL) + ((int)codFunc != IFT_NIL);
    printf("Chosen %s:\n", (numbMethods==1)? "method" : "methods"); 
    if(intPointDetec != IFT_NIL)
        printf("- Interest point detection: %s\n", iftBovwIntPointDetName(intPointDetec, true));

    if(localFeatExtr != IFT_NIL)
        printf("- Local feature extraction: %s\n", iftBovwLocalFeatExtrName(localFeatExtr, true));

    if(dictEstim != IFT_NIL)
        printf("- Dictionary estimation: %s\n", iftBovwDictEstimName(dictEstim, true));

    if(codFunc != IFT_NIL)
        printf("- Codification: %s\n", iftBovwCodFuncName(codFunc, true));
}

char *iftBovwIntPointDetName(iftBovwIntPointDetectorId method, bool fullName)
{
    char *name = NULL;
    switch(method) {
        case BOVW_INT_POINT_DETECTOR_RANDOM:
            name = fullName ? iftCopyString("Random") : iftCopyString("random"); break;
        case BOVW_INT_POINT_DETECTOR_GRID:
            name = fullName ? iftCopyString("Grid") : iftCopyString("grid"); break;
        case BOVW_INT_POINT_DETECTOR_UNSUP_SPIX_ISF:
            name = fullName ? iftCopyString("Unsupervised superpixel by ISF") : iftCopyString("unsup-spix-isf"); break;
        case BOVW_INT_POINT_DETECTOR_SUP_SPIX_ISF:
            name = fullName ? iftCopyString("Supervised superpixel by ISF") : iftCopyString("sup-spix-isf"); break;
        default:
            name = fullName ? iftCopyString("Random") : iftCopyString("random");
    }

    return name;
}

char *iftBovwLocalFeatExtrName(iftBovwLocalFeatExtractorId method, bool fullName)
{
    char *name = NULL;
    switch(method) {
        case BOVW_LOCAL_FEAT_EXTRACTOR_RAW:
            name = fullName ? iftCopyString("Raw") : iftCopyString("raw"); break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_BIC:
            name = fullName ? iftCopyString("BIC") : iftCopyString("bic"); break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_LBP:
            name = fullName ? iftCopyString("LBP") : iftCopyString("lbp"); break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_BRIEF:
            name = fullName ? iftCopyString("BRIEF") : iftCopyString("brief"); break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_CONV:
            name = fullName ? iftCopyString("Convolution") : iftCopyString("conv"); break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_CONV_MULTI_LAYER:
            name = fullName ? iftCopyString("Multi Layer Convolution") : iftCopyString("conv-ml"); break;
        case BOVW_LOCAL_FEAT_EXTRACTOR_DEEP_FEATS_MIMG:
            name = fullName ? iftCopyString("Deep features from Mimage") : iftCopyString("deep-feats-mimg"); break;
        default:
            name = fullName ? iftCopyString("Raw") : iftCopyString("raw");
    }

    return name;
}

char *iftBovwDictEstimName(iftBovwDictEstimatorId method, bool fullName)
{
    char *name = NULL;
    switch(method) {
        case BOVW_DICT_ESTIMATOR_UNSUP_KMEANS:
            name = fullName ? iftCopyString("Unsup k-means") : iftCopyString("unsup-kmeans"); break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_PIX_LABEL:
            name = fullName ? iftCopyString("Sup k-means with ordering by pixel label") : iftCopyString("sup-kmeans-pix-label"); break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMG_CLASS:
            name = fullName ? iftCopyString("Sup k-means with ordering by image class") : iftCopyString("sup-kmeans-img-class"); break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMAGE:
            name = fullName ? iftCopyString("Sup k-means with ordering by image") : iftCopyString("sup-kmeans-image"); break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_POSITION:
            name = fullName ? iftCopyString("Sup k-means with ordering by position") : iftCopyString("sup-kmeans-position"); break;
        case BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMG_CLASS_AND_POSITION:
            name = fullName ? iftCopyString("Sup k-means with ordering by image class and position") : iftCopyString("sup-kmeans-img-class-position"); break;
        case BOVW_DICT_ESTIMATOR_UNSUP_OPF:
            name = fullName ? iftCopyString("Unsup OPF") : iftCopyString("unsup-opf"); break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL:
            name = fullName ? iftCopyString("Sup OPF with ordering by pixel label") : iftCopyString("sup-opf-pix-label"); break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS:
            name = fullName ? iftCopyString("Sup OPF with ordering by image class") : iftCopyString("sup-opf-img-class"); break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE:
            name = fullName ? iftCopyString("Sup OPF with ordering by image") : iftCopyString("sup-opf-image"); break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION:
            name = fullName ? iftCopyString("Sup OPF with ordering by position") : iftCopyString("sup-opf-position"); break;
        case BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS_AND_POSITION:
            name = fullName ? iftCopyString("Sup OPF with ordering by image class and position") : iftCopyString("sup-opf-img-class-position"); break;
        case BOVW_DICT_ESTIMATOR_SUP_MANUAL_IMG_CLASS:
            name = fullName ? iftCopyString("Sup Manual with ordering by image class") : iftCopyString("sup-manual-img-class"); break;
        case BOVW_DICT_ESTIMATOR_SUP_MANUAL_POSITION:
            name = fullName ? iftCopyString("Sup Manual with ordering by position") : iftCopyString("sup-manual-position"); break;
        case BOVW_DICT_ESTIMATOR_SUP_MANUAL_IMG_CLASS_AND_POSITION:
            name = fullName ? iftCopyString("Sup Manual with ordering by image class and position") : iftCopyString("sup-manual-img-class-position"); break;
        default:
            name = fullName ? iftCopyString("Unsup k-Means") : iftCopyString("unsup-kmeans");
    }

    return name;
}

char *iftBovwCodFuncName(iftBovwCodFuncId method, bool fullName)
{
    char *name = NULL;
    switch(method) {
        case BOVW_COD_FUNC_HARD_ASGMT:
            name = fullName ? iftCopyString("Hard assignment") : iftCopyString("hard-asgmt"); break;
        case BOVW_COD_FUNC_SOFT_ASGMT:
            name = fullName ? iftCopyString("Soft assignment") : iftCopyString("soft-asgmt"); break;
        case BOVW_COD_FUNC_HARD_ASGMT_BATCH:
            name = fullName ? iftCopyString("Hard assignment (batch processing)") : iftCopyString("hard-asgmt-batch"); break;
        case BOVW_COD_FUNC_SOFT_ASGMT_BATCH:
            name = fullName ? iftCopyString("Soft assignment (batch processing)") : iftCopyString("soft-asgmt-batch"); break;
        case BOVW_COD_FUNC_2ND_ORDER_STATS_FISHER_VECTORS:
            name = fullName ? iftCopyString("2nd order statistics (fisher vectors)") : iftCopyString("2nd-order-stats-fv"); break;
        default:
            name = fullName ? iftCopyString("Hard assignment") : iftCopyString("hard-asgmt");
    }

    return name;
}

void iftBovwSaveImgWithIntPoints(iftImage *img, iftRoiArray *roiArray, iftFile *imgFilePtr)
{
    iftImage *aux = iftCopyImage(img);
    iftAdjRel *A = iftCircular(0.5);
    iftColor RGB, YCbCr;
    int Imax = iftNormalizationValue(iftMaximumValue(img));
    RGB.val[0] = 0;
    RGB.val[1] = Imax;
    RGB.val[2] = 0;
    YCbCr = iftRGBtoYCbCr(RGB, Imax);
    iftDrawRoiArrayBordersInPlace(aux, roiArray, YCbCr, A, true);
    char filename[2048]; sprintf(filename, "%s_intPoints.png", iftFilename(imgFilePtr->path, iftFileExt(imgFilePtr->path)));
    iftWriteImageByExt(aux, filename);
    iftDestroyImage(&aux);
    iftDestroyAdjRel(&A);
}

int iftBovwComputeBatchSize(int nVisWords, int localFeatsDim, int nPatches, float maxMemUsePerc)
{
    float nMbytesPerFloat = (float)sizeof(float)/1024.0/1024.0;
    float kernelsSize = (float)nVisWords * (float)localFeatsDim * nMbytesPerFloat;

#ifdef IFT_GPU
    float freeMemory = iftGetFreeMemoryGPU(0)/1024.0/1024.0;
#else
    float freeMemory = (float)iftGetFreePhysicalSystemMemory()/1024.0/1024.0;
#endif
    
    return (int)floor((maxMemUsePerc*freeMemory - kernelsSize) / (nPatches*nMbytesPerFloat*(localFeatsDim + (float)nVisWords)));
}

int iftBovwWgtFuncStrToWgtFunc(char *wgtFuncStr)
{
    int wgtFunc = 0;

    if(iftCompareStrings(wgtFuncStr, "inverse"))
        wgtFunc = BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_INVER;
    else if(iftCompareStrings(wgtFuncStr, "gaussian"))
        wgtFunc = BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_GAUSS;

    return wgtFunc;
}

void iftBoVWCreateIntPointsFile(iftRoiArray *roiArray, iftDict *funcParams)
{
    /* get the ouput dir and the output suffix */
    char *outputDir = iftGetStrValFromDict("output_dir", funcParams);
    char *outputSuffix = iftGetStrValFromDict("output_suffix", funcParams);
    
    /* create the file with the int points */
    char *intPointsFile = iftConcatStrings(4, outputDir, "/int_points_", outputSuffix, ".intPts");
    iftWriteRoiArray(roiArray, intPointsFile);
    iftDict *dictAux = iftGetDictFromDict("int_point_detec", funcParams);
    iftInsertIntoDict("int_points_file", intPointsFile, dictAux);
}

void iftBovwCreateRandomKernels(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset)
{
    char *outputDir = iftGetStrValFromDict("output_dir", bovw->funcParams);
    char *outSuffix = iftGetStrValFromDict("output_suffix", bovw->funcParams);
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", bovw->funcParams);
    int kernelSize = iftGetLongValFromDict("kernel_size", dictAux);
    int nKernels = iftGetLongValFromDict("n_kernels", dictAux);

    char *kernelsDir = iftJoinPathnames(3, outputDir, "kernels", outSuffix);
    if(!iftDirExists(kernelsDir))
        iftMakeDir(kernelsDir);

    char *kernName = iftJoinPathnames(2, kernelsDir, "kernels_layer_1.npy");
    iftImage *img0 = iftReadImageByExt(dictLearnFileset->files[0]->path);
    int nBands = iftIsColorImage(img0) ? 3 : 1;

    iftMatrix *kernels = iftRandomKernelBankAsMatrix(kernelSize*kernelSize, nBands, nKernels);
    iftWriteMatrix(kernels, kernName);

    kernelsDir = iftAbsPathname(kernelsDir);
    iftDestroyMatrix(&kernels);
    iftInsertIntoDict("kernels_dir", kernelsDir, dictAux);
}

void iftBovwCreateRandomKernelsMultiLayer(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset)
{
    char *outputDir = iftGetStrValFromDict("output_dir", bovw->funcParams);
    char *outSuffix = iftGetStrValFromDict("output_suffix", bovw->funcParams);
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", bovw->funcParams);
    iftDict *networkParams = iftGetDictFromDict("network_params", dictAux);

    char *kernelsDir = iftJoinPathnames(3, outputDir, "kernels", outSuffix);
    if(!iftDirExists(kernelsDir))
        iftMakeDir(kernelsDir);

    int nLayers = iftGetLongValFromDict("n_layers", networkParams);
    for(int l = 0; l < nLayers; l++) {
        char layerName[128];
        sprintf(layerName, "layer_%d", l+1);
        iftDict *layerDict = iftGetDictFromDict(layerName, networkParams);
        
        iftDict *convDict = iftGetDictFromDict("conv", layerDict);
        int nBands = iftGetLongValFromDict("n_bands", convDict);
        int nKernels = iftGetLongValFromDict("n_kernels", convDict);
        int kernelSize = iftGetLongValFromDict("kernel_size", convDict);

        char kernName[2048];
        sprintf(kernName, "%s/kernels_layer_%d.npy", kernelsDir, l+1);
        iftMatrix *kernels = iftRandomKernelBankAsMatrix(kernelSize*kernelSize, nBands, nKernels);
        iftWriteMatrix(kernels, kernName);
        iftDestroyMatrix(&kernels);
    }

    kernelsDir = iftAbsPathname(kernelsDir);
    iftInsertIntoDict("kernels_dir", kernelsDir, dictAux);
}

void iftBovwCreateJointSamplingMask(iftBagOfVisualWords *bovw, iftFileSet **dictLearnMasksFileset)
{
    printf("\n--> Creating joint sampling mask ... "); fflush(stdout);
    char *outputDir = iftGetStrValFromDict("output_dir", bovw->funcParams);
    iftDict *dictAux = iftGetDictFromDict("int_point_detec", bovw->funcParams);
    char *maskJoinOperStr = iftGetStrValFromDict("mask_join_operation", dictAux);
    iftMaskJoinOperation maskJoinOper = iftMaskJoinOperStrToMaskJoinOper(maskJoinOperStr);

    iftImage *jointMask = iftCreateJointMaskFromFileset((*dictLearnMasksFileset), maskJoinOper, 255);
    char *jointMaskDir = iftJoinPathnames(2, outputDir, "sampling_masks");
    if(!iftDirExists(jointMaskDir))
        iftMakeDir(jointMaskDir);
    char *jointMaskFilename = iftJoinPathnames(2, jointMaskDir, "joint_mask.png");
    iftWriteImageByExt(jointMask, jointMaskFilename);
    
    iftFileSet *newFileset = iftCreateFileSet((*dictLearnMasksFileset)->n);
    for(int f = 0; f < (*dictLearnMasksFileset)->n; f++)
        newFileset->files[f] = iftCreateFile(jointMaskFilename);
    
    iftDestroyFileSet(dictLearnMasksFileset);
    (*dictLearnMasksFileset) = iftCopyFileSet(newFileset);
    iftDestroyFileSet(&newFileset);
    iftDestroyImage(&jointMask);
    iftInsertIntoDict("sampling_masks", jointMaskFilename, bovw->funcParams);
    printf("OK\n");
}

void iftBovwCreateSamplingMasksPerClassISF(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset, iftFileSet *dictLearnMasksFileset)
{
    printf("\n--> Creating sampling masks per class using ISF  ...\n");
    char *outputDir = iftGetStrValFromDict("output_dir", bovw->funcParams);
    iftDict *dictAux = iftGetDictFromDict("int_point_detec", bovw->funcParams);
    char *maskJoinOperStr = iftGetStrValFromDict("mask_join_operation", dictAux);
    double adjRelRad = iftGetDblValFromDict("adj_rel_rad", dictAux);
    char *colSpaceStr = iftGetStrValFromDict("col_space", dictAux);
    int spixSeedsStride = iftGetLongValFromDict("spix_seeds_stride", dictAux);
    double alpha = iftGetDblValFromDict("alpha", dictAux);
    double beta = iftGetDblValFromDict("beta", dictAux);
    int nIter = iftGetLongValFromDict("n_iter", dictAux);

    iftMaskJoinOperation maskJoinOper = iftMaskJoinOperStrToMaskJoinOper(maskJoinOperStr);
    iftColorSpace inputColSpace = YCbCr_CSPACE;
    iftColorSpace colSpace = iftColorSpaceStrToColorSpace(colSpaceStr);

    int nClasses = iftFileSetLabelsNumber(dictLearnFileset);

    /* create dir for the sampling masks */
    char maskStr[1024];
    sprintf(maskStr, "sup-spix-isf_spix-seeds-stride%d_adjrel%.2f_alpha%.2f_beta%.2f_niter%d",
        spixSeedsStride, adjRelRad, alpha, beta, nIter);
    
    char *masksDir = iftJoinPathnames(3, outputDir, "sampling_masks", maskStr);
    if(!iftDirExists(masksDir))
        iftMakeDir(masksDir);

    /* validate the color spaces */
    iftImage *img0 = iftReadImageByExt(dictLearnFileset->files[0]->path);
    bool isColorImg = iftIsColorImage(img0);

    if((isColorImg && (inputColSpace == GRAY_CSPACE || inputColSpace == GRAYNorm_CSPACE)) ||
        (!isColorImg && (inputColSpace != GRAY_CSPACE && inputColSpace != GRAYNorm_CSPACE))) {
        iftError("The number of image channels and the specified col_space_in does not match",
            "iftBovwCreateSamplingMasksPerClassISF");
    }

    if(!isColorImg && colSpace != GRAY_CSPACE && colSpace != GRAYNorm_CSPACE) {
        iftError("Color RAW image representation can only be used with color images. The given image is grayscale",
            "iftBovwCreateSamplingMasksPerClassISF");
    }
    
    /* create a mask for each class and run the ISF algorithm with the images of each class separatedly */
    iftDict *masksPerClassDict = iftCreateDict();
    iftDict *spixMasksPerClassDict = iftCreateDict();
    iftDict *spixLabelsPerClassDict = iftCreateDict();
    int nBands = iftIsColorImage(img0) ? 3 : 1;
    int normValue = iftNormalizationValue(iftMaximumValue(img0));

    for(int c = 1; c < nClasses+1; c++) {
        printf("- class %d ... ", c); fflush(stdout);
        /* extract the samples from the given class */
        iftFileSet *fsPerClass = iftExtractClassFileset(dictLearnFileset, c);
        iftFileSet *fsMaskPerClass = iftExtractClassFileset(dictLearnMasksFileset, c);
        iftImage *maskPerClass = iftCreateJointMaskFromFileset(fsMaskPerClass, maskJoinOper, 255);

        /* create a iftMImage containing the LAB values of all the images in the given class */
        iftMImage *mimgPerClass = iftCreateMImage(img0->xsize, img0->ysize, img0->zsize, nBands*fsPerClass->n);
        for(int i = 0; i < fsPerClass->n; i++) {
            iftImage *img = iftReadImageByExt(fsPerClass->files[i]->path);

            for(int p = 0; p < img->n; p++) {
                /* convert to the chosen color space */
                iftFColor colorIn = {.val[0] = 0, .val[1] = 0, .val[2] = 0};
                colorIn.val[0] = img->val[p];
                if(isColorImg) {
                    colorIn.val[1] = img->Cb[p];
                    colorIn.val[2] = img->Cr[p];
                }
                iftFColor colorOut = iftConvertPixelColorSpace(colorIn, inputColSpace, colSpace, normValue);

                /* copy the features in the dataset */
                if(nBands == 1) {
                    mimgPerClass->val[p][i] = colorOut.val[0];
                }
                else {
                    mimgPerClass->val[p][i*nBands] = colorOut.val[0];
                    mimgPerClass->val[p][i*nBands+1] = colorOut.val[1];
                    mimgPerClass->val[p][i*nBands+2] = colorOut.val[2];
                }
            }
            iftDestroyImage(&img);
        }

        /* execute the ISF algorithm on the mimage created */
        iftAdjRel *A = iftIs3DImage(img0) ? iftSpheric(adjRelRad) : iftCircular(adjRelRad);
        iftIntArray *grid = iftGridSamplingOnMask(maskPerClass, spixSeedsStride, -1, -1); // we use all the points resulting from the given spix_seeds_stride
        iftImage *seeds = iftCreateImageFromImage(img0);
        iftIntArrayToImage(grid, seeds, 1);
        
        iftIGraph *igraph = iftExplicitIGraph(mimgPerClass, maskPerClass, NULL, A);
        iftIGraphISF_Root(igraph, seeds, alpha, beta, nIter);
        iftImage *spixLabelsPerClass = iftIGraphLabel(igraph);

        /* save the mask, the superpixel boundaries mask and the superpixel labes for the given class */
        iftImage *spixMaskPerClass = iftCreateSuperpixelBoundariesMask(spixLabelsPerClass, 255);
        char filename[2048], dictName[2048];
        
        sprintf(filename, "%s/mask_class_%d.png", masksDir, c);
        iftWriteImageByExt(maskPerClass, filename);
        sprintf(dictName, "class_%d", c);
        iftInsertIntoDict(dictName, filename, masksPerClassDict);
        iftDestroyImage(&maskPerClass);

        sprintf(filename, "%s/spix_mask_class_%d.png", masksDir, c);
        iftWriteImageByExt(spixMaskPerClass, filename);
        sprintf(dictName, "class_%d", c);
        iftInsertIntoDict(dictName, filename, spixMasksPerClassDict);
        iftDestroyImage(&spixMaskPerClass);
        
        sprintf(filename, "%s/spix_labels_class_%d.png", masksDir, c);
        iftWriteImageByExt(spixLabelsPerClass, filename);
        sprintf(dictName, "class_%d", c);
        iftInsertIntoDict(dictName, filename, spixLabelsPerClassDict);
        iftDestroyImage(&spixLabelsPerClass);

        printf("OK -> %lu superpixels extracted\n", grid->n);
    }
    
    iftDict *sampMasksDict = iftCreateDict();
    iftInsertIntoDict("masks_per_class", masksPerClassDict, sampMasksDict);
    iftInsertIntoDict("spix_masks_per_class", spixMasksPerClassDict, sampMasksDict);
    iftInsertIntoDict("spix_labels_per_class", spixLabelsPerClassDict, sampMasksDict);
    iftInsertIntoDict("sampling_masks", sampMasksDict, bovw->funcParams);
}

char *iftBovwCreateDirectoryForPatches(iftBagOfVisualWords *bovw)
{
    char patchStr[1024];
    iftDict *dictAux = iftGetDictFromDict("int_point_detec", bovw->funcParams);
    int intPointDetector = bovw->intPointDetectorId;
    char *outputDir = iftGetStrValFromDict("output_dir", bovw->funcParams);
    
    if(intPointDetector == BOVW_INT_POINT_DETECTOR_RANDOM) {
        int patchSize = iftGetLongValFromDict("patch_size_x", dictAux);
        sprintf(patchStr, "%s_patch%d", iftBovwIntPointDetName(intPointDetector, false), patchSize);
    }
    else if(intPointDetector == BOVW_INT_POINT_DETECTOR_GRID) {
        int patchSize = iftGetLongValFromDict("patch_size_x", dictAux);
        int stride = iftGetLongValFromDict("stride", dictAux);
        sprintf(patchStr, "%s_patch%d_stride%d", iftBovwIntPointDetName(intPointDetector, false), patchSize, stride);
    }
    else if(intPointDetector == BOVW_INT_POINT_DETECTOR_UNSUP_SPIX_ISF) {
        int patchSize = iftGetLongValFromDict("patch_size", dictAux);
        int nSuperpixels = iftGetLongValFromDict("n_superpixels", dictAux);
        sprintf(patchStr, "%s_patch%d_nspix%d", iftBovwIntPointDetName(intPointDetector, false), patchSize, nSuperpixels);
    }
    else if(intPointDetector == BOVW_INT_POINT_DETECTOR_SUP_SPIX_ISF) {
        int patchSize = iftGetLongValFromDict("patch_size", dictAux);
        int intPtsStride = iftGetLongValFromDict("int_pts_stride", dictAux);
        int spixSeedsStride = iftGetLongValFromDict("spix_seeds_stride", dictAux);
        float adjRelRad = iftGetDblValFromDict("adj_rel_rad", dictAux);
        float alpha = iftGetDblValFromDict("alpha", dictAux);
        float beta = iftGetDblValFromDict("beta", dictAux);
        int nIter = iftGetLongValFromDict("n_iter", dictAux);
        sprintf(patchStr, "%s_patch%d_int-pts-stride%d_spix-seeds-stride%d_adjrel%.2f_alpha%.2f_beta%.2f_niter%d",
            iftBovwIntPointDetName(intPointDetector, false), patchSize, intPtsStride, spixSeedsStride, adjRelRad,
            alpha, beta, nIter);
    }
    
    char *patchDirPath = iftJoinPathnames(3, outputDir, "patches", patchStr);
    if(!iftDirExists(patchDirPath))
        iftMakeDir(patchDirPath);

    return patchDirPath;
}

void iftBovwApplyLocalBatchCentralization(iftBagOfVisualWords *bovw, iftFileSet *imgFileset)
{
    char *outputDir = iftGetStrValFromDict("output_dir", bovw->funcParams);
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", bovw->funcParams);
    char *colSpaceStr = iftGetStrValFromDict("col_space", dictAux);
    iftColorSpace colSpace = iftColorSpaceStrToColorSpace(colSpaceStr);

    int kernelSize = IFT_NIL;
    if(iftDictContainKey("network_params", dictAux, NULL)) {
        iftDict *networkParams = iftGetDictFromDict("network_params", dictAux);
        iftDict *layerDict = iftGetDictFromDict("layer_1", networkParams);
        iftDict *convDict = iftGetDictFromDict("conv", layerDict);
        kernelSize = iftGetLongValFromDict("kernel_size", convDict);
    }
    else {
        kernelSize = iftGetLongValFromDict("kernel_size", dictAux);
    }
    
    char kernDir[2048]; sprintf(kernDir, "kernel_%dx%d", kernelSize, kernelSize);
    char *centralizedImgsDir = iftJoinPathnames(3, outputDir, "centralized_imgs", kernDir);
    
    if(!iftDirExists(centralizedImgsDir))
        iftMakeDir(centralizedImgsDir);
    centralizedImgsDir = iftAbsPathname(centralizedImgsDir);

    /* verify if the path of centralizaed images already contains the images */
    printf("\n--> Applying local batch centralization to the images ... "); fflush(stdout);
    iftFileSet *fs = iftLoadFileSetFromDirOrCSV(centralizedImgsDir, 1, true);
    
    bool allFiles = true;
    for(int i = 0; i < imgFileset->n; i++) {
        char *basename = iftFilename(imgFileset->files[i]->path, iftFileExt(imgFileset->files[i]->path));
        if(iftFindFileByBasename(fs, basename, ".mimg") == IFT_NIL) {
            allFiles = false;
            break;
        }
        iftFree(basename);
    }

    if(allFiles) {
        printf("images already centralized!\n");
    }
    else {
        iftLocalBatchCentralization(imgFileset, kernelSize, centralizedImgsDir, colSpace);
        printf("OK\n");
    }
    iftInsertIntoDict("centralized_imgs_dir", centralizedImgsDir, dictAux);
    iftDestroyFileSet(&fs);
}

iftDataSet *iftBovwBuildLocalFeatureVectorsFromConvolutionalMImage(iftMImage *mimg, iftImage *img, iftFile *imgFilePtr, iftRoiArray *roiArray, char *resizeImage, char *vectConstr)
{
    iftDataSet *localFeats = NULL;
    iftMImage *aux = NULL;
    int nBands = mimg->m;

    if(iftCompareStrings(resizeImage, "yes") || (mimg->xsize == img->xsize && mimg->ysize == img->ysize)) {
        /* resize the resulting mimage */
        if(mimg->xsize != img->xsize || mimg->ysize != img->ysize) {
            aux = iftResizeMImage(mimg, img->xsize, img->ysize, img->zsize);
            iftDestroyMImage(&mimg);
            mimg = iftCopyMImage(aux);
            iftDestroyMImage(&aux);
        }
        
        /* for each chosen position, extract the features in its neighborhood (from the resized deep features mimage) */
        int nPositions = roiArray->n;
        if(iftCompareStrings(vectConstr, "pixel_vals_in_patch"))
            localFeats = iftCreateDataSet(nPositions, nBands*roiArray->val[iftBiggestRoiInArray(roiArray)]->n);
        else if(iftCompareStrings(vectConstr, "mean_val_in_patch"))
            localFeats = iftCreateDataSet(nPositions, nBands);
        localFeats->function_number = 1;
        iftSetDistanceFunction(localFeats, localFeats->function_number);

        #pragma omp parallel for
        for(int s = 0; s < roiArray->n; s++) {
            iftRoi *roi = roiArray->val[s];
            localFeats->sample[s].id = s;
            iftAddSampleStatus(&localFeats->sample[s], IFT_SUPERVISED);
            localFeats->sample[s].truelabel = roi->label;

            if(iftCompareStrings(vectConstr, "pixel_vals_in_patch")) {
                int col = 0;
                for (int v = 0; v < roi->n; v++) {
                    int pixId = iftMGetVoxelIndex(mimg, roi->val[v]);
                    for(int b = 0; b < nBands; b++) {
                        localFeats->sample[s].feat[col] = mimg->val[pixId][b];
                        col++;
                    }
                }
            }
            else if(iftCompareStrings(vectConstr, "mean_val_in_patch")) {
                for(int b = 0; b < nBands; b++) {
                    float mean = 0;
                    for (int v = 0; v < roi->n; v++) {
                        int pixId = iftMGetVoxelIndex(mimg, roi->val[v]);
                        mean += mimg->val[pixId][b];
                    }
                    localFeats->sample[s].feat[b] = mean/=roi->n;
                }
            }
        }
    }
    else {
        /* for each pixel, extract the features in its neighborhood (from the original/reduced deep features mimage) */
        /* warning: the interest points previously detected are no longer used (all the positions from the reduced image are used instead) */
        if(iftCompareStrings(vectConstr, "pixel_vals_in_patch"))
            localFeats = iftCreateDataSet(mimg->n, nBands*roiArray->val[iftBiggestRoiInArray(roiArray)]->n);
        else if(iftCompareStrings(vectConstr, "mean_val_in_patch"))
            localFeats = iftCreateDataSet(mimg->n, nBands);
        localFeats->function_number = 1;
        iftSetDistanceFunction(localFeats, localFeats->function_number);
        int patchSize = sqrt(roiArray->val[0]->n);
        iftAdjRel *A = iftRectangular(patchSize, patchSize);

        #pragma omp parallel for
        for(int s = 0; s < mimg->n; s++) {
            localFeats->sample[s].id = s;
            iftAddSampleStatus(&localFeats->sample[s], IFT_SUPERVISED);
            localFeats->sample[s].truelabel = imgFilePtr->label;

            iftVoxel v = iftMGetVoxelCoord(mimg, s);
            if(iftCompareStrings(vectConstr, "pixel_vals_in_patch")) {
                int col = 0;
                for (int i = 0; i < A->n; i++) {
                    iftVoxel u = iftGetAdjacentVoxel(A, v, i);
                    if(iftMValidVoxel(mimg, u)) {
                        int pixId = iftMGetVoxelIndex(mimg, u);
                        for(int b = 0; b < nBands; b++) {
                            localFeats->sample[s].feat[col] = mimg->val[pixId][b];
                            col++;
                        }
                    }
                }
            }
            else if(iftCompareStrings(vectConstr, "mean_val_in_patch")) {
                for(int b = 0; b < nBands; b++) {
                    float mean = 0;
                    for (int i = 0; i < A->n; i++) {
                        iftVoxel u = iftGetAdjacentVoxel(A, v, i);
                        if(iftMValidVoxel(mimg, u)) {
                            int pixId = iftMGetVoxelIndex(mimg, u);
                            mean += mimg->val[pixId][b];
                        }
                    }
                    localFeats->sample[s].feat[b] = mean/=A->n;
                }
            }
        }
        iftDestroyAdjRel(&A);
    }

    return localFeats;
}

iftDataSet *iftBovwExtractBICFeatureVectorsFromConvolutionalMImage(iftMImage *mimg, iftImage *img, iftRoiArray *roiArray, char *resizeImage, int nBinsPerBand)
{
    /* resize the mimage */
    if(iftCompareStrings(resizeImage, "yes")) {
        if(mimg->xsize != img->xsize || mimg->ysize != img->ysize) {
            iftMImage *aux = iftResizeMImage(mimg, img->xsize, img->ysize, img->zsize);
            iftDestroyMImage(&mimg);
            mimg = iftCopyMImage(aux);
            iftDestroyMImage(&aux);
        }
    }

    /* extract features from all the ROIs */
    int nFeats = nBinsPerBand*mimg->m*2;
    iftDataSet *localFeats = iftCreateDataSet((int)roiArray->n, nFeats);
    localFeats->function_number = 1;
    iftSetDistanceFunction(localFeats, localFeats->function_number);

    //#pragma omp parallel for
    for(int s = 0; s < roiArray->n; s++) {
        iftRoi *roi = roiArray->val[s];
        localFeats->sample[s].id = s;
        iftAddSampleStatus(&localFeats->sample[s], IFT_SUPERVISED);
        localFeats->sample[s].truelabel = roi->label;

        /* extract BIC features */
        iftMImage *mimgRoi = iftMExtractRoiNoBkgd(mimg, roi);
        iftFeatures *feats = iftMExtractBIC(mimgRoi, nBinsPerBand);

        for (int f = 0; f < feats->n; f++)
            localFeats->sample[s].feat[f] = feats->val[f];

        iftDestroyFeatures(&feats);
        iftDestroyMImage(&mimgRoi);
    }

    return localFeats;
}

void iftBovwComputeMultiLayerConvolutionalFeaturesBatch(iftBagOfVisualWords *bovw, iftFileSet *imgFileset)
{
    char *outputDir = iftGetStrValFromDict("output_dir", bovw->funcParams);
    char *outSuffix = iftGetStrValFromDict("output_suffix", bovw->funcParams);
    iftDict *dictAux = iftGetDictFromDict("local_feat_extr", bovw->funcParams);
    iftDict *networkParams = iftGetDictFromDict("network_params", dictAux);
    char *centralizedImgsDir = iftGetStrValFromDict("centralized_imgs_dir", dictAux);
    char *kernelsDir = iftGetStrValFromDict("kernels_dir", dictAux);
    int batchSize = iftGetLongValFromDict("batch_size", dictAux);
    bool printDetail = true;

    if(!iftDirExists(centralizedImgsDir))
        iftError("The directory with the centralized images does not exist", "iftComputeMultiLayerConvolutionalFeatures");

    if(!iftDirExists(kernelsDir))
        iftError("The kernels directory does not exist", "iftComputeMultiLayerConvolutionalFeatures");

    /* create the directory for the filtered images */
    char *filteredImgsDir = iftJoinPathnames(3, outputDir, "filtered_imgs", outSuffix);
    if(!iftDirExists(filteredImgsDir))
        iftMakeDir(filteredImgsDir);

    /* extract the deep features */
    printf("\n--> Extracting the deep features from the images ... "); fflush(stdout);
    int nTrainImgs = imgFileset->n;
    int nBatches = ceil((float)nTrainImgs / (float)batchSize);
    char filename[2048];
    for(int b = 0; b < nBatches; b++) {
        int nImgsInBatch = iftMin((b+1)*batchSize, nTrainImgs) - b*batchSize;
        if(printDetail) { printf("\nBatch %d/%d (%d imgs): ", b+1, nBatches, nImgsInBatch); fflush(stdout);}
        else {printf("Batch %d/%d (%d imgs)\r", b+1, nBatches, nImgsInBatch); fflush(stdout);}

        /* create an array (batch) with the centralized mimages */
        iftMImageArray *mimgBatch = iftCreateMImageArray(nImgsInBatch), *aux = NULL;
        for (int imgId = b*batchSize, i = 0; imgId < iftMin((b+1)*batchSize, nTrainImgs); imgId++, i++) {
            char *basename = iftFilename(imgFileset->files[imgId]->path, iftFileExt(imgFileset->files[imgId]->path));
            sprintf(filename, "%s/%s.mimg", centralizedImgsDir, basename);
            mimgBatch->val[i] = iftReadMImage(filename);
            iftFree(basename);
        }

        /* apply the operations in every layer using the current batch */
        int nLayers = iftGetLongValFromDict("n_layers", networkParams);

        for(int l = 0; l < nLayers; l++) {
            if(printDetail) {printf("layer %d ", l+1); fflush(stdout);}
            char layerName[2048];
            sprintf(layerName, "layer_%d", l+1);
            iftDict *layerDict = iftGetDictFromDict(layerName, networkParams);
            
            /* convolution */
            iftDict *convDict = iftGetDictFromDict("conv", layerDict);
            int nBands = iftGetLongValFromDict("n_bands", convDict);
            int nKernels = iftGetLongValFromDict("n_kernels", convDict);
            int kernelSize = iftGetLongValFromDict("kernel_size", convDict);
            int stride = iftGetLongValFromDict("stride", convDict);

            char kernName[2048];
            sprintf(kernName, "%s/kernels_layer_%d.npy", kernelsDir, l+1);
            iftMatrix *kernels = iftReadMatrix(kernName);

            if(nKernels != kernels->ncols || kernelSize != sqrt(kernels->nrows/nBands))
                iftError("The provided kernels for layer %d (%d x %d) does not match the n_kernels (%d) and kernel_size (%d)\n",
                        "iftComputeMultiLayerConvolutionalFeatures", l+1, kernels->ncols, kernels->nrows, nKernels, kernelSize);

            aux = iftMConvolutionArray(mimgBatch, kernels, stride);
            iftDestroyMImageArray(&mimgBatch);
            mimgBatch = iftCopyMImageArray(aux);
            iftDestroyMImageArray(&aux);

            /* relu */
            char *relu = iftGetStrValFromDict("relu", layerDict);
            if(iftCompareStrings(relu, "yes")) {
                aux = iftMReLUArray(mimgBatch);
                iftDestroyMImageArray(&mimgBatch);
                mimgBatch = iftCopyMImageArray(aux);
                iftDestroyMImageArray(&aux);
            }

            /* max-pooling */
            if(iftDictContainKey("max-pool", layerDict, NULL)) {
                iftDict *maxPoolDict = iftGetDictFromDict("max-pool", layerDict);
                int poolSize = iftGetLongValFromDict("size", maxPoolDict);
                int poolStride = iftGetLongValFromDict("stride", maxPoolDict);

                aux = iftMMaxPoolingArray(mimgBatch, poolSize, poolStride);
                iftDestroyMImageArray(&mimgBatch);
                mimgBatch = iftCopyMImageArray(aux);
                iftDestroyMImageArray(&aux);
            }
            if(printDetail) {printf("(OK)%s", l < nLayers-1 ? ", " : " "); fflush(stdout);}

            iftDestroyMatrix(&kernels);
        }

        /* save the resulting (filtered mimages) */
        for (int imgId = b*batchSize, i = 0; imgId < iftMin((b+1)*batchSize, nTrainImgs); imgId++, i++) {
            char *basename = iftFilename(imgFileset->files[imgId]->path, iftFileExt(imgFileset->files[imgId]->path));
            sprintf(filename, "%s/%s.mimg", filteredImgsDir, basename);
            iftWriteMImage(mimgBatch->val[i], filename);
            iftFree(basename);
        }
        iftDestroyMImageArray(&mimgBatch);
    }
    
    filteredImgsDir = iftAbsPathname(filteredImgsDir);
    iftInsertIntoDict("filtered_imgs_dir", filteredImgsDir, dictAux);
}

void iftBovwComputeClusterModelParamsFV(iftDataSet *Z, iftMatrix **means, iftMatrix ***covariances, iftMatrix ***covariancesInv,
                                        iftFloatArray **covariancesDet, iftFloatArray **weights)
{
    (*means) = iftCreateMatrix(Z->nfeats, Z->ngroups);
    (*covariances) = (iftMatrix **)iftAlloc(Z->ngroups, sizeof(iftMatrix*));
    // (*covariancesRotMat) = (iftMatrix **)iftAlloc(Z->ngroups, sizeof(iftMatrix*));
    // (*covariancesFSP) = (iftFeatSpaceParam **)iftAlloc(Z->ngroups, sizeof(iftFeatSpaceParam*));
    (*covariancesInv) = (iftMatrix **)iftAlloc(Z->ngroups, sizeof(iftMatrix*));
    (*covariancesDet) = iftCreateFloatArray(Z->ngroups);
    (*weights) = iftCreateFloatArray(Z->ngroups);

    for(int g = 0; g < Z->ngroups; g++)
        (*covariances)[g] = NULL;

    for(int g = 0; g < Z->ngroups; g++) {
        iftDataSet *Zg = iftExtractGroup(Z, g+1);

        /* means */
        iftMatrix *m = iftMatrixMeanColumn(Zg->data);
        iftCopyFloatArray(iftMatrixRowPointer((*means), g), iftMatrixRowPointer(m, 0), Zg->nfeats);

        /* covariances */
        iftSetStatus(Zg, IFT_TRAIN);
        iftCentralizeDataSetInPlace(Zg);
        (*covariances)[g] = iftCovarianceMatrix(Zg->data);

        /* covariance matrices determinant */
        iftMatrix *U,*S,*Vt;
        iftSingleValueDecomp((*covariances)[g], &U, &S, &Vt);
        (*covariancesDet)->val[g] = 1.0;
        for(int i = 0; i < S->n; i++)
            if(S->val[i]*S->val[i] > 0)
                (*covariancesDet)->val[g] *= S->val[i]*S->val[i];

        /* covariance SVD rotation matrix */
        // (*covariancesRotMat)[g] = iftTransposeMatrix(Vt);

        /* covariance SVD feature space params */
        // (*covariancesFSP)[g]->mean = iftAllocFloatArray(Zg->fsp.nfeats);
        // iftCopyFloatArray((*covariancesFSP)[g]->mean, Zg->fsp.mean, Zg->fsp.nfeats);

        /* covariance inverse matrices */
        iftMatrix *V = iftTransposeMatrix(Vt);
        iftMatrix *Ut = iftTransposeMatrix(U);
        iftMatrix *Sd = iftCreateDiagonalMatrix((double*)iftMatrixRowPointer(S, 0), S->ncols);
        iftMatrix *Sinv = iftInvertDiagonalMatrix(Sd);
        (*covariancesInv)[g] = iftMultMatricesChain(3, V, Sinv, Ut);
        // Sinv -> cov matrix inverse in the Y space
        // V*Sinv*Ut -> cov matrix inverse in the X space

        /* weights */
        (*weights)->val[g] = (float)Zg->nsamples/(float)Z->nsamples;

        iftDestroyDataSet(&Zg);
        iftDestroyMatrix(&m);
        iftDestroyMatrix(&U);
        iftDestroyMatrix(&S);
        iftDestroyMatrix(&Vt);
        iftDestroyMatrix(&V);
        iftDestroyMatrix(&Ut);
        iftDestroyMatrix(&Sd);
        iftDestroyMatrix(&Sinv);
    }
}

iftFeatures *iftBovwComputeFeatVect2ndOrderStatsFV(iftDataSet *localFeats, iftMatrix *means, iftMatrix **covariances,
                                                    iftMatrix **covariancesInv, iftFloatArray *covariancesDet,
                                                    iftFloatArray *weights)
{
    int nGroups = weights->n;
    int nSamples = localFeats->nsamples;
    iftMatrix *variances = iftCreateMatrix(means->ncols, means->nrows);
    iftMatrix *Q = iftCreateMatrix(nGroups, localFeats->nsamples);

    /* compute the posterior probabilities (PDFs) */
    for(int g = 0; g < nGroups; g++) {
        float aux = 1.0/sqrtf(powf(2.0*IFT_PI, localFeats->nfeats)*covariancesDet->val[g]);

        for(int i = 0; i < localFeats->nsamples; i++) {
            iftMatrix *diff = iftCreateMatrix(localFeats->nfeats, 1);
            iftSubstractFloatArrays(localFeats->sample[i].feat,
                                    iftMatrixRowPointer(means, g),
                                    iftMatrixRowPointer(diff, 0),
                                    localFeats->nfeats);
            iftMatrix *diffT = iftTransposeMatrix(diff);
            iftMatrix *mult = iftMultMatricesChain(3, diff, covariancesInv[g], diffT);

            iftMatrixElem(Q, g, i) = aux * expf(-0.5*iftMatrixElem(mult, 0, 0));
            
            iftDestroyMatrix(&diff);
            iftDestroyMatrix(&diffT);
            iftDestroyMatrix(&mult);
        }
        iftMatrix *var = iftDiagonalFromMatrix(covariances[g], 'r');
        iftCopyFloatArray(iftMatrixRowPointer(variances,g), var->val, var->ncols);
        iftDestroyMatrix(&var);
    }

    /* Compute the sufficient statistics of descriptors */
    iftMatrix *xx = iftDataSetToFeatureMatrix(localFeats);
    iftMatrix *Qsum = iftMatrixSumColumn(Q);
    iftMultMatrixByScalar(Qsum, 1.0/(float)nSamples);
    iftMatrix *QsumT = iftTransposeMatrix(Qsum);

    iftMatrix *QT = iftTransposeMatrix(Q);
    iftMatrix *Qxx = iftMultMatrices(QT, xx);
    iftMultMatrixByScalar(Qxx, 1.0/(float)nSamples);

    iftMatrix *xx2 = iftPointWiseMatrixPow(xx, 2.0);
    iftMatrix *Qxx2 = iftMultMatrices(QT, xx2);
    iftMultMatrixByScalar(Qxx2, 1.0/(float)nSamples);

    /* Compute derivatives with respect to mixing weights, means and variances */
    /* d(Pi) */
    iftFloatArray *dPi = iftCreateFloatArray(weights->n);
    iftSubstractFloatArrays(Qsum->val, weights->val, dPi->val, weights->n);
    /* d(Mu) */
    iftMatrix *dMu = iftCreateMatrix(means->ncols, means->nrows);
    iftMatrix *QsumT_means = iftCopyMatrix(means);
    iftMultiplicationMatrixVectorByColumn(QsumT_means, QsumT);
    iftSubstractFloatArrays(Qxx->val, QsumT_means->val, dMu->val, Qxx->n);
    /* d(Sigma) */
    iftMatrix *dSigma = iftCreateMatrix(means->ncols, means->nrows);
    iftMatrix *QsumT_means2 = iftPointWiseMatrixPow(means, 2.0);
    iftMultiplicationMatrixVectorByColumn(QsumT_means2, QsumT);
    iftMatrix *QsumT_variances = iftCopyMatrix(variances);
    iftMultiplicationMatrixVectorByColumn(QsumT_variances, QsumT);
    iftMatrix *Qxx_means = iftMatricesMultiplicationPointWise(Qxx, means);
    for(int i = 0; i < dSigma->n; i++)
        dSigma->val[i] = - Qxx2->val[i] - QsumT_means2->val[i] + QsumT_variances->val[i] + 2*Qxx_means->val[i];

    /* Merge derivatives into a vector */
    iftFeatures *featVect = iftCreateFeatures(nGroups + 2*nGroups*localFeats->nfeats);
    iftCopyFloatArray(featVect->val, dPi->val, dPi->n);
    iftCopyFloatArray(featVect->val+dPi->n, dMu->val, dMu->n);
    iftCopyFloatArray(featVect->val+dPi->n+dMu->n, dSigma->val, dSigma->n);
    // iftNormalizeFeatures(featVect->val, featVect->n);

    /* free memory */
    iftDestroyMatrix(&variances);
    iftDestroyMatrix(&Q);
    iftDestroyMatrix(&xx);
    iftDestroyMatrix(&Qsum);
    iftDestroyMatrix(&QsumT);
    iftDestroyMatrix(&QT);
    iftDestroyMatrix(&Qxx);
    iftDestroyMatrix(&xx2);
    iftDestroyMatrix(&Qxx2);
    iftDestroyFloatArray(&dPi);
    iftDestroyMatrix(&dMu);
    iftDestroyMatrix(&QsumT_means);
    iftDestroyMatrix(&dSigma);
    iftDestroyMatrix(&QsumT_means2);
    iftDestroyMatrix(&QsumT_variances);
    iftDestroyMatrix(&Qxx_means);

    return featVect;
}
