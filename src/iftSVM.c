#include "iftSVM.h"

#include "ift/core/io/Dir.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"


#define MODEL_SEP "===\n"

/**
@file
@brief A description is missing here
*/
void iftWriteSVMModelTrainingParams(FILE *fout, const svmModel *model);

bool iftReadSVMModelTrainingParams(FILE *fin, svmModel *model);

void iftCopySVMModelTrainingParams(svmModel *src, svmModel *dst);

void iftReadSVMInfo(iftSVM *svm, const char *filename);

void iftWriteSVMInfo(const iftSVM *svm, const char *filename);

#if 0
void printSample(svmNode *x, double label)
{

  int index, f;
  double value;

  printf("%g ", label);
  f=0;
  while(1){
    index = x[f].index;
    value = x[f].value;

    if (index != 0 && index !=-1)
      printf("%g ", value);

    //if (index !=-1)
    //  printf("%d:%g ", index, value);

    if (index == -1)
      break;
    f++;

  }
}

void printSVMProblem(svmProblem *problem)
{

  printf("\n");
  printf("number of samples in the problem: %d\n", problem->l);

  printf("sample labels\n");
  int s;
  for(s=0; s < problem->l; s++)
    printf("%g ", problem->y[s]);
  printf("\n");

  printf("sample feature vectors\n");

  for(s=0; s < problem->l; s++){
    printSample(problem->x[s], problem->y[s]);
    printf("\n");
  }
}
#endif

/********************** PRIVATE FUNCTIONS *************************/
svmParameter *iftDefaultLibSVMParams()
{
    svmParameter *param = (svmParameter *) iftAlloc(1, sizeof(svmParameter));

    // Default libsvm parameter values
    param->svm_type = C_SVC;
    param->kernel_type = RBF;
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

    return(param);

}


iftSVM *iftCreateSVC(double C)
{

    iftSVM *svm = (iftSVM *) iftAlloc(1, sizeof(iftSVM));

    svmParameter *params = iftDefaultLibSVMParams();

    params->svm_type = C_SVC;
    params->C = C;

    svm->multiClass = IFT_OVA;
    svm->kernelization = 0;
    svm->params  = params;
    svm->nmodels = 0;
    svm->Z       = NULL;

    return(svm);
}

void iftSampleToSvmSample(iftSample ift, svmNode *svm,
                          int s, int nfeats, int precomp)
{
    if (precomp) {
        svm[0].index = 0;
        svm[0].value = (double) s + 1;
    }

    int f;
    for (f=0; f < nfeats; f++) {
        svm[precomp + f].index = f + 1;
        svm[precomp + f].value = (double) ift.feat[f];
    }

    // by convention of libsvm, the last sample node has index set to -1
    svm[precomp + nfeats].index = -1;
    svm[precomp + nfeats].value = 0.0;
}


svmProblem *iftDataSetToSVMProblem(iftSVM *svc, const iftDataSet *Z)
{
    int nsamples = Z->nsamples;
    int ntrainsamples = Z->ntrainsamples;
    int nfeats = Z->nfeats;

    int s, strain;

    //allocate problem
    svmProblem *problem = (svmProblem *) iftAlloc(1, sizeof(svmProblem));
    problem->x = (svmNode **) iftAlloc(ntrainsamples, sizeof(svmNode *));
    problem->y = (double *) iftAlloc(ntrainsamples, sizeof(double));

    int precomp = svc->params->kernel_type == IFT_PRECOMPUTED;

    if (precomp && ntrainsamples != nfeats)
        iftError("Invalid kernel matrices", "iftDataSetToSVMProblem");

    // as libsvm uses a sparse matrix representation, each sample has an
    // additional node denoting its termination. if samples represent precomputed
    // kernels, another node at the beginning of the sample vector is required to
    // identify the samples.
    int sample_svmnodes = precomp + nfeats + 1;

    for (s=0, strain=0; s < nsamples; s++){

        if (iftHasSampleStatus(Z->sample[s], IFT_TRAIN)){

            problem->x[strain] = (svmNode *) iftAlloc(sample_svmnodes,
                                                    sizeof(svmNode));

            iftSampleToSvmSample(Z->sample[s], problem->x[strain],
                                 strain, nfeats, precomp);

            // set sample label
            problem->y[strain] = (double) Z->sample[s].truelabel;

            strain++;
        }
    }

    if (strain != ntrainsamples)
        iftError("Wrong number of train samples in dataset",
                 "iftDataSetToSVMProblem");

    problem->l = ntrainsamples;

    return(problem);
}


void iftDestroyProblemAndModels(iftSVM *svm)
{
    //destroy problem
    int s;
    if (svm->problem != NULL) {
        if (svm->problem->x != NULL) {
            for (s=0; s < svm->problem->l; s++)
                iftFree(svm->problem->x[s]);

            iftFree(svm->problem->x);
        }

        if (svm->problem->y != NULL)
            iftFree(svm->problem->y);

        iftFree(svm->problem);
        svm->problem = NULL;
    }

    // destroy models
    int imodel;
    for (imodel=0; imodel < svm->nmodels; imodel++) {

        if (svm->model[imodel] != NULL) {

            // this is a bug in libsvm, which does not free this
//            iftFree(svm->model[imodel]->sv_indices);
            // bug fixed on libsvm 3.17 version
            svm_free_and_destroy_model(&(svm->model[imodel]));
        }
    }
    iftFree(svm->model);
    svm->model = NULL;
    iftFree(svm->truelabel);
    svm->truelabel = NULL;
}

void iftDestroyProblem(iftSVM *svm) {
    int s;
    if (svm->problem != NULL) {
        for (s = 0; s < svm->problem->l; s++)
            iftFree(svm->problem->x[s]);

        iftFree(svm->problem->x);
        iftFree(svm->problem->y);
        iftFree(svm->problem);
        svm->problem = NULL;
    }
}

void iftDestroyModel(iftSVM *svm, int idx_model) {
    if (svm->model[idx_model] != NULL) {

        // this is a bug in libsvm, which does not free this
        iftFree(svm->model[idx_model]->sv_indices);
        // bug fixed on libsvm 3.17 version

        svm_free_and_destroy_model(&(svm->model[idx_model]));
    }
    svm->truelabel[idx_model] = 0.0;
}
/******************************************************************/

/********************** PUBLIC FUNCTIONS *************************/

iftSVM* iftCreateSVM(iftKernelType kernel, iftMultiClass multiclass, double C, double sigma) {

    iftSVM* svm = NULL;

    switch (kernel) {
        case IFT_LINEAR:
            svm = iftCreateLinearSVC(C);
            break;
        case IFT_RBF:
            svm = iftCreateRBFSVC(C, sigma);
            break;
        case IFT_PRECOMPUTED:
            svm = iftCreatePreCompSVC(C);
            break;
        default:
            iftWarning("Invalid kernel type. Creating linear SVM.", "iftCreateSVM");
            svm = iftCreateLinearSVC(C);
    }

    svm->multiClass = multiclass;

    return svm;
}

iftSVM *iftCreateLinearSVC(double C) {
    iftSVM *svm = iftCreateSVC(C);
    svm->params->kernel_type = IFT_LINEAR;
    return svm;
}


iftSVM *iftCreateRBFSVC(double C, double sigma) {
    iftSVM *svm = iftCreateSVC(C);
    svm->params->kernel_type = RBF;
    svm->params->gamma = sigma;
    return svm;
}


iftSVM *iftCreatePreCompSVC(double C) {
    iftSVM *svm = iftCreateSVC(C);
    svm->params->kernel_type = IFT_PRECOMPUTED;
    return svm;
}


iftSVM* iftCreateSVMFromDict(const iftDict* d) {
	char* kernel = iftGetStrValFromDict("kernel", d);
	double C = iftGetDblValFromDict("C", d);

	iftSVM* svm = NULL;

	if(iftCompareStrings(kernel, "RBF")) {
		double sigma = iftGetDblValFromDict("sigma", d);
		svm = iftCreateRBFSVC(C, sigma);
	}
	else if(iftCompareStrings(kernel, "linear")) {
		svm = iftCreateLinearSVC(C);
	}
	else if(iftCompareStrings(kernel, "precomputed")){
		svm = iftCreatePreCompSVC(C);
	}

	iftFree(kernel);

	return svm;

}

void iftDestroySVM(iftSVM *svm) {
    if (svm != NULL) {
        svm_destroy_param(svm->params);
        iftFree(svm->params);
        iftDestroyProblemAndModels(svm);
        iftFree(svm);
    }
}

void iftDestroySVMClassifier(iftSVM **svm) {
    iftSVM* aux = *svm;
    svm_destroy_param(aux->params);
    iftFree(aux->params);
    iftDestroyProblemAndModels(aux);
    iftFree(aux);
    *svm = NULL;
}


iftSVM *iftReadSVM(const char *pathname) {
    char *path = NULL;
    char *tmp_dir = NULL;
    iftSVM *svm = NULL;

    if(pathname == NULL)
        iftError("Input Pathname must not be NULL", "iftReadSVM");

    if(!iftFileExists(pathname))
        iftError("Input Pathname does not refer to a valid file", "iftReadSVM", pathname);

    if(!iftCompareStrings(iftFileExt(pathname), ".zip"))
        iftError("Input Pathname for the SVM must end with a \".zip\" file extension", "iftReadSVM");

    tmp_dir = iftMakeTempDir("tmp_svm", NULL, NULL);

    if(!iftDirExists(tmp_dir))
        iftError("Error when creating temporary dir %s for unzipping the SVM file content", "iftReadSVM",
                 tmp_dir);

    iftUnzipFile(pathname, tmp_dir);

    /* Reading info */
    path = iftJoinPathnames(2, tmp_dir, "info.data");

    if (!iftFileExists(path))
        iftError("Cannot open file: \"%s\". File \"%s\" does not refer to a valid svm", "iftReadSVM",
                 path, pathname);

    /* Reading SVM info */

    svm = (iftSVM *) iftAlloc(1, sizeof(iftSVM));

    iftReadSVMInfo(svm, path);

    iftFree(path);

    path = iftJoinPathnames(2, tmp_dir, "dataset.zip");

    if (!iftFileExists(path))
        svm->Z = NULL;
    else
        svm->Z = iftReadOPFDataSet(path);

    iftFree(path);

    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);

    return svm;
}


void iftReadSVMInfo(iftSVM *svm, const char *filename) {
    FILE *fin = fopen(filename, "rb");
    char *tempfile = iftMakeTempFile(NULL, NULL, NULL);
    FILE *ft = fopen(tempfile, "wb");


    size_t nbytes = 1000000;
    char *line = iftAllocCharArray(nbytes);

    if (line == NULL) {
        iftError("Memory Error", "iftReadSVMInfo");
    }

    int nmodels = 0;
    int imodel = 0;
    int kernelization = 0;
    int ret = fscanf(fin, "%d\n", &nmodels);
    
    if (ret <= 0)
        iftError("Error Reading SVM Model 1", "iftReadSVMInfo");

    ret = fscanf(fin, "%d\n", &kernelization);
    
    if (ret <= 0)
        iftError("Error Reading Kernelization parameter", "iftReadSVMInfo");

    if (fgets(line, nbytes, fin) == NULL)
        iftError("Error Reading SVM Model 2", "iftReadSVMInfo");

    if (strcmp(line, "[]\n") == 0) {
        svm->truelabel = NULL;
    }
    else {
        svm->truelabel = iftAllocDoubleArray(nmodels);

        //line+1 to skip the [
        char *l = strtok(line + 1, " ");
        float f;
        while (l != NULL) {
            if (sscanf(l, "%f", &f) == 1) {
                svm->truelabel[imodel++] = f;
            }
            l = strtok(NULL, " ");
        }

    }
    if (fgets(line, nbytes, fin) == NULL)
        iftError("Error Reading SVM Model 3", "iftReadSVMInfo");

    svm->model = (svmModel **) iftAlloc(nmodels, sizeof(svmModel *));

    imodel = 0;

    for (; ;) {
        int *sv_indices = NULL;
        svmModel model;

        if (fgets(line, nbytes, fin) == NULL) {
            break;
        }

        if(strcmp(line, "SV_INDICES\n") == 0) {
            int nSV;

            if(fscanf(fin, "%d", &nSV) != 1)
                iftError("Error reading the number of support vectors!", "iftReadSVMInfo");

            sv_indices = iftAllocIntArray(nSV);

            for(int i = 0; i < nSV; i++) {
                if(fscanf(fin, "%d", &sv_indices[i]) != 1)
                    iftError("Error reading the support vector %03d/%3d!", "iftReadSVMInfo", i + 1, nSV);
            }

            fgetc(fin); // Eliminating \n

            if(!iftReadSVMModelTrainingParams(fin, &model))
                iftWarning("SVM model training parameters not read!", "iftReadSVMInfo");

            if(fgets(line, nbytes, fin) == NULL)
                break;

        }

        if (strcmp(line, MODEL_SEP) == 0) {
            fclose(ft);
            svm->model[imodel] = svm_load_model(tempfile);

            iftCopySVMModelTrainingParams(&model, svm->model[imodel]);

            if(sv_indices != NULL) {

                if(svm->model[imodel]->sv_indices != NULL)
                    iftFree(svm->model[imodel]->sv_indices);
                svm->model[imodel]->sv_indices = sv_indices;
            }

            imodel++;

            ft = fopen(tempfile, "wb");
            continue;
        }

        fputs(line, ft);

        if(sv_indices != NULL)
            iftFree(sv_indices);
    }

    svm->params = iftDefaultLibSVMParams();
    *(svm->params) = svm->model[0]->param;
    svm->nmodels = nmodels;
    svm->kernelization = kernelization;
    svm->multiClass = svm->truelabel==NULL? IFT_OVO : IFT_OVA;

    fclose(fin);
    fclose(ft);
    iftFree(line);

    iftRemoveFile(tempfile);
    iftFree(tempfile);
}

bool iftReadSVMModelTrainingParams(FILE *fin, svmModel *model) {/* Training parameters */
    bool success = false;
    size_t nbytes = 1000000;
    char *line = iftAllocCharArray(nbytes);

    if(fgets(line, nbytes, fin) == NULL)
        return success;

    if(iftCompareStrings(line, "TRAINING_PARAMS\n")) {

        if(fscanf(fin, "CACHE_SIZE %lf\n", &model->param.cache_size) != 1)
            iftError("Error when reading SVM 1\n", "iftReadSVMModelTrainingParams");
        if(fscanf(fin, "EPS %lf\n", &model->param.eps) != 1)
            iftError("Error when reading SVM 2\n", "iftReadSVMModelTrainingParams");
        if(fscanf(fin, "C %lf\n", &model->param.C) != 1)
            iftError("Error when reading SVM 3\n", "iftReadSVMModelTrainingParams");
        if(fscanf(fin, "NR_WEIGHT %d\n", &model->param.nr_weight) != 1)
            iftError("Error when reading SVM 4\n", "iftReadSVMModelTrainingParams");

        /* Reading weight label */

        if (fgets(line, nbytes, fin) == NULL)
            iftError("Error Reading SVM Model 5", "iftReadSVMModelTrainingParams");

        if(iftCompareStrings(line, "WEIGHT_LABEL_NULL\n")) {
            model->param.weight_label = NULL;
        } else {
            model->param.weight_label = iftAllocIntArray(model->param.nr_weight);

            iftSList *list = iftSplitString(line, " ");
            iftSNode *node = list->head->next;
            int i = 0;

            for(; node != NULL; node = node->next){
                sscanf(node->elem, "%d", &model->param.weight_label[i]);
                i++;
            }

            iftDestroySList(&list);
        }

        /* Reading weights */
        if (fgets(line, nbytes, fin) == NULL)
            iftError("Error Reading SVM Model 6", "iftReadSVMModelTrainingParams");

        if(iftCompareStrings(line, "WEIGHT_NULL\n")) {
            model->param.weight = NULL;
        } else {
            model->param.weight = iftAllocDoubleArray(model->param.nr_weight);

            iftSList *list = iftSplitString(line, " ");
            iftSNode *node = list->head->next;
            int i = 0;
            for(; node != NULL; node = node->next){
                sscanf(node->elem, "%lf", &model->param.weight[i]);
                i++;
            }
            fgetc(fin);

            iftDestroySList(&list);
        }

        if(fscanf(fin, "NU %lf\n", &model->param.nu) != 1)
            iftError("Error when reading SVM 7\n", "iftReadSVMModelTrainingParams");
        if(fscanf(fin, "P %lf\n", &model->param.p) != 1)
            iftError("Error when reading SVM 8\n", "iftReadSVMModelTrainingParams");
        if(fscanf(fin, "SHRINKING %d\n", &model->param.shrinking) != 1)
            iftError("Error when reading SVM 9\n", "iftReadSVMModelTrainingParams");
        if(fscanf(fin, "PROBABILITY %d\n", &model->param.probability) != 1)
            iftError("Error when reading SVM 10\n", "iftReadSVMModelTrainingParams");

        iftFree(line);

        success = true;
    }

    return success;
}


void iftCopySVMModelTrainingParams(svmModel *src, svmModel *dst) {
    dst->param.cache_size = src->param.cache_size;
    dst->param.eps = src->param.eps;
    dst->param.C = src->param.C;
    dst->param.nr_weight = src->param.nr_weight;
    dst->param.weight_label = src->param.weight_label;
    dst->param.weight = src->param.weight;
    dst->param.nu = src->param.nu;
    dst->param.p = src->param.p;
    dst->param.shrinking = src->param.shrinking;
    dst->param.probability = src->param.probability;
}


void iftWriteSVM(const iftSVM *svm, const char *pathname) {
    char *path = NULL;
    char *tmp_dir = NULL;

    if (svm == NULL)
        iftError("SVM is NULL", "iftWriteSVM");
    if (pathname == NULL)
        iftError("Out Pathname is NULL", "iftWriteSVM");

    //    if (svm->Z == NULL)
    //    iftError("SVM\'s reference dataset is NULL", "iftWriteSVM");

    if(!iftCompareStrings(iftFileExt(pathname), ".zip"))
        iftError("The output Pathname for the SVM model must end with extension \".zip\"", "iftWriteSVM");

    /* Creating temporary dir */
    tmp_dir = iftMakeTempDir("tmp_svm", NULL, NULL);

    /* Opening file that will contain the basic information about the graph */
    path = iftJoinPathnames(2, tmp_dir, "info.data");
    iftWriteSVMInfo(svm, path);
    iftFree(path);

    // Writing the data graph's set
    path = iftJoinPathnames(2, tmp_dir, "dataset.zip");

    if(svm->Z != NULL)
        iftWriteOPFDataSet(svm->Z, path);

    iftFree(path);

    // Zippping the folder content
    iftZipDirContent(tmp_dir, pathname);

    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
}


void iftWriteSVMInfo(const iftSVM *svm, const char *filename) {
    FILE *fout = fopen(filename, "wb");
    if (fout == NULL)
        iftError("Error when opening the Filename: %s", "iftWriteSVMInfo", filename);

    const char *tempfile = iftMakeTempFile(NULL, NULL, NULL);
    char c;

    int nmodels = svm->nmodels;
    svmModel **model = svm->model;
    int kernelization = svm->kernelization;

    fprintf(fout, "%d\n", nmodels);

    fprintf(fout, "%d\n", kernelization);
    if (svm->truelabel == NULL) {
        fprintf(fout, "[]\n");
    }
    else {
        fprintf(fout, "[");
        for (int i = 0; i < svm->nmodels; ++i) {
            fprintf(fout, " %.1lf", svm->truelabel[i]);
        }
        fprintf(fout, " ]\n");
    }

    fputs(MODEL_SEP, fout);

    for (int i = 0; i < nmodels; ++i) {

        if (model[i] == NULL)
            continue;

        int ret = svm_save_model(tempfile, model[i]);

        if (ret < 0)
            iftError("Error Saving SVM Model", "iftWriteSVMInfo");

        //copy file content
        FILE *ft = fopen(tempfile, "rb");
        for (; ;) {
            c = fgetc(ft);
            if (c == EOF)
                break;
            fputc(c, fout);
        }

        if(model[i]->sv_indices != NULL) {
            fprintf(fout, "SV_INDICES\n%d\n", model[i]->l);

            for(int j = 0; j < model[i]->l-1; j++)
                fprintf(fout, "%d ", model[i]->sv_indices[j]);
            if(model[i]->l > 0)
                fprintf(fout, "%d\n", model[i]->sv_indices[model[i]->l-1]);
        }

        iftWriteSVMModelTrainingParams(fout, model[i]);

        //append blank line
        fputs(MODEL_SEP, fout);
    }



    fclose(fout);

    iftRemoveFile(tempfile);
}

void iftWriteSVMModelTrainingParams(FILE *fout, const svmModel *model) {/* Training parameters */
    fprintf(fout, "TRAINING_PARAMS\n");
    fprintf(fout, "CACHE_SIZE %lf\n", model->param.cache_size);
    fprintf(fout, "EPS %lf\n", model->param.eps);
    fprintf(fout, "C %lf\n", model->param.C);
    fprintf(fout, "NR_WEIGHT %d\n", model->param.nr_weight);

    if(model->param.weight_label == NULL)
        fprintf(fout, "WEIGHT_LABEL_NULL\n");
    else {
        fprintf(fout, "WEIGHT_LABEL");
        for(int i = 0; i < model->param.nr_weight; i++)
            fprintf(fout, " %d", model->param.weight_label[i]);
        fprintf(fout, "\n");
    }

    if(model->param.weight == NULL)
        fprintf(fout, "WEIGHT_NULL\n");
    else {
        fprintf(fout, "WEIGHT");
        for(int i = 0; i < model->param.nr_weight; i++)
            fprintf(fout, " %lf", model->param.weight[i]);
        fprintf(fout, "\n");
    }

    fprintf(fout, "NU %lf\n", model->param.nu);
    fprintf(fout, "P %lf\n", model->param.p);
    fprintf(fout, "SHRINKING %d\n", model->param.shrinking);
    fprintf(fout, "PROBABILITY %d\n", model->param.probability);
}

iftSVMHyperplane *iftCreateSVMHyperplane(int n) {
    iftSVMHyperplane *h = (iftSVMHyperplane *) iftAlloc(1, sizeof(iftSVMHyperplane));

    h->feat = iftAllocFloatArray(n);
    h->bias = 0.0;
    h->n = n;

    return h;
}


void iftDestroySVMHyperplane(iftSVMHyperplane **h) {
    if (h != NULL) {
        iftSVMHyperplane *aux = *h;

        if (aux->feat != NULL)
            iftFree(aux->feat);
        iftFree(aux);
    }
}

iftSample iftSVMGetNormalHyperplaneAsSample(const iftSVM *svm, int model_idx, const iftDataSet *Ztrain, float *out_bias) {
    int f, sv;
    int nfeats = Ztrain->nfeats;

    if ( (svm->params->kernel_type != IFT_LINEAR) &&
         (svm->params->kernel_type != IFT_PRECOMPUTED) )
        iftError("This function is not able to compute the decision hyperplane "
                 "in the nonlinear case.",
                 "iftSVMGetNormalHyperplane");

    // hyperplane to be returned
    iftSample w;
    w.feat = (float*) iftAllocFloatArray(nfeats);

    for(sv=0; sv < svm->model[model_idx]->l; sv++) {
        //printf("sv_indices: %d\n", svm->model[idxModel]->sv_indices[sv]);

        for(f=0; f < nfeats; f++)
            w.feat[f] += svm->model[model_idx]->sv_coef[0][sv] *
                         Ztrain->sample[svm->model[model_idx]->sv_indices[sv]-1].feat[f];
    }

    for(f=0; f<nfeats; f++)
        w.feat[f] *= (float) svm->model[model_idx]->label[0];

    *out_bias = svm->model[model_idx]->label[0] * svm->model[model_idx]->rho[0];

    return w;
}

iftSVMHyperplane *iftSVMGetNormalHyperplane(iftSVM *svm, int model_idx, iftDataSet *Ztrain) {
    int f, sv;
    int nfeats = Ztrain->nfeats;

    if ( (svm->params->kernel_type != IFT_LINEAR) &&
         (svm->params->kernel_type != IFT_PRECOMPUTED) )
        iftError("This function is not able to compute the decision hyperplane "
                 "in the nonlinear case.",
                 "iftSVMGetNormalHyperplane");

    // hyperplane to be returned
    iftSVMHyperplane *hyper = iftCreateSVMHyperplane(nfeats);

    for(sv=0; sv < svm->model[model_idx]->l; sv++) {
        //printf("sv_indices: %d\n", svm->model[model_idx]->sv_indices[sv]);
        for(f=0; f < nfeats; f++)
            hyper->feat[f] += svm->model[model_idx]->sv_coef[0][sv] *
                              Ztrain->sample[svm->model[model_idx]->sv_indices[sv]-1].feat[f];
    }

    for(f=0; f<nfeats; f++)
        hyper->feat[f] *= (float) svm->model[model_idx]->label[0];

    hyper->bias = svm->model[model_idx]->label[0] * svm->model[model_idx]->rho[0];

    return hyper;
}


iftSVMHyperplane **iftSVMGetAllNormalHyperplanes(iftSVM *svm, iftDataSet *Ztrain) {
    if ( (svm->params->kernel_type != IFT_LINEAR) && (svm->params->kernel_type != IFT_PRECOMPUTED) )
        iftError("This function is not able to compute the decision hyperplane "
                 "in the nonlinear case.",
                 "iftSVMGetAllNormalHyperplanes");

    iftSVMHyperplane **hyperplanes = (iftSVMHyperplane**) iftAlloc(svm->nmodels, sizeof(iftSVMHyperplane*));

    
    if (svm->nmodels == 2) {
        hyperplanes[0] = iftSVMGetNormalHyperplane(svm, 0, Ztrain);
        hyperplanes[1] = iftSVMGetNormalHyperplane(svm, 0, Ztrain);
    }
    else {
        for (int m = 0; m < svm->nmodels; m++) {
            hyperplanes[m] = iftSVMGetNormalHyperplane(svm, m, Ztrain);
        }
    }

    return hyperplanes;
}

void iftSVMTrain(iftSVM* svm, iftDataSet* Z) {

    if(svm->multiClass==IFT_OVA) {
        iftSVMTrainOVA(svm, Z, Z);
    }
    else {
        iftSVMTrainOVO(svm, Z, Z);
    }

}

void iftSVMClassify(iftSVM* svm, iftDataSet* Z, uchar status) {

    if (svm->multiClass == IFT_OVA) {
        if (svm->params->kernel_type == IFT_LINEAR || svm->params->kernel_type == IFT_PRECOMPUTED) {
            iftSVMLinearClassifyOVA(svm, Z, status, NULL);
        }
        else {
            iftSVMClassifyOVA(svm, Z, status);
        }
    } else {
        if(svm->params->probability==1)
           iftSVMClassifyOVO_Probability(svm, Z, status);
        else
           iftSVMClassifyOVO(svm, Z, status); 
    }
}


void iftSVMTrainOVO(iftSVM *svm, const iftDataSet *Z, iftDataSet *Zref) {
    // if svm was previously trained
    if (svm->nmodels > 0)
        iftDestroyProblemAndModels(svm);

    svm->problem = iftDataSetToSVMProblem(svm, Z);

    svm->model = (svmModel **) iftAlloc(1, sizeof(svmModel *));
    svm->model[0] = svm_train(svm->problem, svm->params);

    svm->truelabel  = NULL;
    svm->nmodels    = 1;
    svm->Z          = Zref;
    svm->multiClass = IFT_OVO;
}


void iftSVMTrainOVA(iftSVM *svm, const iftDataSet *Ztrain, iftDataSet *Zref) {
    int nclasses = Ztrain->nclasses;
    int ntrainsamples = Ztrain->ntrainsamples;

    int idxtrain, idxtrained, s, istrained=1;
    double classtrain=-1.0;

    // if svm was previously trained
    if (svm->nmodels > 0)
        iftDestroyProblemAndModels(svm);
    svm->nmodels   = 0;

    svm->problem   = iftDataSetToSVMProblem(svm, Ztrain);

    svm->model     = (svmModel **) iftAlloc(nclasses, sizeof(svmModel *));
    svm->truelabel = iftAllocDoubleArray(nclasses);
    svm->Z         = Zref;
    svm->multiClass = IFT_OVA;

    //get original classes from svm problem samples
    double *origclasses = iftAllocDoubleArray(ntrainsamples);
    for (s=0; s < ntrainsamples; s++)
        origclasses[s] = svm->problem->y[s];

    //train an OVA SVM for each class
    for (idxtrain=0; idxtrain < nclasses; idxtrain++) {

        //search next class to train
        for (s=0; s < ntrainsamples; s++) {

            classtrain = origclasses[s];
            istrained = 0;
            for (idxtrained=0; idxtrained < idxtrain; idxtrained++) {

                istrained = svm->truelabel[idxtrained] == classtrain;
                if (istrained)
                    break;
            }
            if (!istrained)
                break;
        }

        // assume that every class has at least one representative in train set
        if (istrained)
            iftError("Problem on finding next class to train", "iftSVMTrainOVA");

        // set svm problem labels accordingly
        for (s=0; s < ntrainsamples; s++) {
            if (origclasses[s] == classtrain) {
                svm->problem->y[s] = +1.0;
            } else {
                svm->problem->y[s] = -1.0;
            }
        }

        // hack that avoid training two models in the originally binary case
        if (nclasses == 2 && idxtrain == 1)
            svm->model[idxtrain] = NULL;
        else
            svm->model[idxtrain] = svm_train(svm->problem, svm->params);

        svm->truelabel[idxtrain] = classtrain;
        svm->nmodels++;
    }

    if (svm->nmodels != nclasses)
        iftError("Number of trained models do not match number of classes",
                 "iftSVMTrainOVA");

    // return correct labels to svm problem just to keep it consistent
    for (s=0; s < ntrainsamples; s++)
        svm->problem->y[s] = origclasses[s];

    iftFree(origclasses);
}


int iftSVMClassifyOVO(const iftSVM *svm, iftDataSet *Z, uchar sample_status) {
    int s, nsamples = Z->nsamples;
    int nfeats = Z->nfeats;

    int nerrors = 0;
    int precomp = svm->params->kernel_type == IFT_PRECOMPUTED;

    if (svm->nmodels == 0)
        iftError("SVM not trained yet", "iftSVMClassifyOVO");

    if (svm->truelabel != NULL)
        iftError("This is not a one-versus-one SVM", "iftSVMClassifyOVO");

    // temporary buffer to get predictions from libsvm
    svmNode *x = (svmNode *) iftAlloc(precomp + nfeats + 1, sizeof(svmNode));

    for (s=0; s < nsamples; s++){

        if (iftHasSampleStatus(Z->sample[s], sample_status)) {
            iftSampleToSvmSample(Z->sample[s], x, s, nfeats, precomp);
            Z->sample[s].label = (int) svm_predict(svm->model[0], x, &Z->sample[s].weight);

            if (Z->sample[s].label != Z->sample[s].truelabel){
                nerrors++;
            }
        }
    }

    iftFree(x);

    return(nerrors);
}



/*

It performs SVM prediction using p-SVM. 
To call this function, it needs to set svm->params->probability to 1 and svm->multiClass to IFT_OVO.

*/

int iftSVMClassifyOVO_Probability(const iftSVM *svm, iftDataSet *Z, uchar sample_status) {
    int s, nsamples = Z->nsamples;
    int nfeats = Z->nfeats;
    double *prob_estimates;
    int nerrors = 0;
    int precomp = svm->params->kernel_type == IFT_PRECOMPUTED;
    prob_estimates = (double *) malloc(Z->nclasses*sizeof(double));

    if (svm->nmodels == 0)
        iftError("SVM not trained yet", "iftSVMClassifyOVO_Probability");

    if (svm->truelabel != NULL)
        iftError("This is not a p-SVM", "iftSVMClassifyOVO_Probability");

    // temporary buffer to get predictions from libsvm
    svmNode *x = (svmNode *) iftAlloc(precomp + nfeats + 1, sizeof(svmNode));

    for (s=0; s < nsamples; s++){


        if (Z->sample[s].status & sample_status) {

            iftSampleToSvmSample(Z->sample[s], x, s, nfeats, precomp);

            Z->sample[s].label = (int) svm_predict_probability(svm->model[0], x, prob_estimates);

            if (Z->sample[s].label != Z->sample[s].truelabel){
                nerrors++;
            }

            Z->sample[s].weight = prob_estimates[Z->sample[s].label-1];
        }

        
    }

    free(prob_estimates);

    iftFree(x);

    return(nerrors);
}

int iftSVMClassifyOVA(const iftSVM *svm, iftDataSet *Z, uchar sample_status) {
    int m, s, nsamples = Z->nsamples;
    int nfeats = Z->nfeats;

    int nerrors = 0;
    int precomp = svm->params->kernel_type == IFT_PRECOMPUTED;

    if (svm->nmodels == 0)
        iftError("SVM not trained yet", "iftSVMClassifyOVA");

    if (svm->truelabel == NULL)
        iftError("This is not a one-versus-all SVM", "iftSVMClassifyOVA");

    // temporary buffer to get predictions from libsvm
    svmNode *x = (svmNode *) iftAlloc(precomp + nfeats + 1, sizeof(svmNode));

    int imaxpred, nmodels = svm->nmodels;
    double maxpred, *prediction = iftAllocDoubleArray(nmodels);

    for (s = 0; s < nsamples; s++) {

        if (iftHasSampleStatus(Z->sample[s], sample_status)) {

            iftSampleToSvmSample(Z->sample[s], x, s, nfeats, precomp);

            for (m = 0; m < nmodels; m++) {

                svm_predict_values(svm->model[m], x, &(prediction[m]), NULL);
                prediction[m] *= (double) svm->model[m]->label[0];

                // hack that avoid using two models in the originally binary case
                if (nmodels == 2) {
                    prediction[m + 1] = -prediction[m];
                    break;
                }
            }

            //get index of maximum value prediction
            imaxpred = 0;
            maxpred = prediction[0];

            for (m = 1; m < nmodels; m++) {
                if (prediction[m] > maxpred) {
                    imaxpred = m;
                    maxpred = prediction[m];
                }
            }

            Z->sample[s].label = svm->truelabel[imaxpred];
            // changed by Daniel Osaku at Jul 1st 2019
            //Z->sample[s].weight = maxpred;
            Z->sample[s].weight = 1 - exp(-fabs(maxpred));

            if (Z->sample[s].label != Z->sample[s].truelabel) {
                nerrors++;

            }
        }
    }

    iftFree(prediction);
    iftFree(x);

    return (nerrors);
}


int iftSVMLinearClassifyOVA(const iftSVM *svm, iftDataSet *Z, uchar sample_status, iftMatrix **out_pred_mat) {
    int s, f, m;
    int nmodels = svm->nmodels;
    iftMatrix *pred_mat = iftCreateMatrix(svm->nmodels, Z->nsamples);
    if (svm->nmodels == 0)
        iftError("SVM not trained yet", "iftSVMLinearClassifyOVA");

    if (svm->truelabel == NULL)
        iftError("This is not a one-versus-all SVM", "iftSVMLinearClassifyOVA");

    if ((svm->params->kernel_type != IFT_LINEAR) &&
        (svm->params->kernel_type != IFT_PRECOMPUTED))
        iftError("This function is not able to compute the decision hyperplane "
                 "in the nonlinear case.",
                 "iftSVMLinearClassifyOVA");

    if(svm->Z == NULL)
        iftError("The SVM model trained for OVA classification does not contain the expected training data set",
                 "iftSVMLinearClassifyOVA");

    iftDataSet *Ztrain = svm->Z;

    iftSample *w = (iftSample *) iftAlloc(nmodels, sizeof(iftSample));
    float *rho = (float *) iftAllocFloatArray(nmodels);

    for (m = 0; m < nmodels; m++) {
        w[m] = iftSVMGetNormalHyperplaneAsSample(svm, m, Ztrain, &(rho[m]));
        if (nmodels == 2) {
            w[m + 1] = iftSVMGetNormalHyperplaneAsSample(svm, m, Ztrain, &(rho[m + 1]));
            break;
        }
    }

    float maxpred, *prediction = iftAllocFloatArray(nmodels);
    int imaxpred, nerrors = 0;

    for (s = 0; s < Z->nsamples; s++) {
        if (iftHasSampleStatus(Z->sample[s], sample_status)) {

            for (m = 0; m < nmodels; m++) {

                prediction[m] = 0.;
                for (f = 0; f < Z->nfeats; f++)
                    prediction[m] += w[m].feat[f] * Z->sample[s].feat[f];
                prediction[m] -= rho[m];

                iftMatrixElem(pred_mat, m, s) = prediction[m];

                if (nmodels == 2) {
                    prediction[m + 1] = -prediction[m];
                    iftMatrixElem(pred_mat, m+1, s) = prediction[m + 1];
                    break;
                }
            }

            //get index of maximum value prediction
            imaxpred = 0;
            maxpred = prediction[0];
            for (m = 1; m < nmodels; m++)
                if (prediction[m] > maxpred) {
                    imaxpred = m;
                    maxpred = prediction[m];
                }

            Z->sample[s].label = svm->truelabel[imaxpred];
            Z->sample[s].weight = prediction[imaxpred];

            if (Z->sample[s].label != Z->sample[s].truelabel) {
                nerrors++;
            }
        }
    }

    // returns the prediction matrix
    if (out_pred_mat != NULL)
        *out_pred_mat = pred_mat;
    else iftDestroyMatrix(&pred_mat);

    iftFree(prediction);
    for (m = 0; m < svm->nmodels; m++)
        iftFree(w[m].feat);
    iftFree(w);
    iftFree(rho);

    return nerrors;
}


iftDataSet *iftKernelizeDataSet(iftDataSet *Ztrain, iftDataSet *Ztest, int kFunction, bool traceNormalize,
                                float *ktrace) {
    int refsamples = Ztrain->nsamples;
    int insamples = Ztest->nsamples;

    if (refsamples <= 0)
        iftError("Reference dataset has no sample to kernelize from",
                 "iftKernelizeDataSet");

    if (insamples <= 0)
        iftError("input dataset has no sample to be kernelized",
                 "iftKernelizeDataSet");

    int nfeats = Ztrain->nfeats;

    if (nfeats != Ztest->nfeats)
        iftError("Number of features on both datasets do not match",
                 "iftKernelizeDataSet");

    if (kFunction != IFT_LINEAR)
        iftError("Kernel function not supporte yet", "iftKernelizeDataSet");


    // kernelization

    /* Create feature matrix of the training dataset */

    iftMatrix *X = iftDataSetToFeatureMatrix(Ztrain);
    iftMatrix *Xt = iftTransposeMatrix(X);

    iftDestroyMatrix(&X);

    /* Create feature matrix of the input dataset */

    X = iftDataSetToFeatureMatrix(Ztest);

    /* Multiply matrices: kernelization */

    iftMatrix *Xk = iftMultMatrices(X, Xt);

    iftDestroyMatrix(&X);
    iftDestroyMatrix(&Xt);

    /* Create dataset from feature matrix */

    iftDataSet *Zk = iftFeatureMatrixToDataSet(Xk);

    Zk->ngroups       = Ztest->ngroups;
    Zk->nclasses      = Ztest->nclasses;
    Zk->ntrainsamples = Ztest->ntrainsamples;
    Zk->iftArcWeight  = Ztest->iftArcWeight;
    for (int si = 0; si < insamples; si++) {
        Zk->sample[si].truelabel = Ztest->sample[si].truelabel;
        Zk->sample[si].label     = Ztest->sample[si].label;
        Zk->sample[si].id        = Ztest->sample[si].id;
        Zk->sample[si].weight    = Ztest->sample[si].weight;
        iftSetSampleStatus(&Zk->sample[si], Ztest->sample[si].status);
    }
    iftDestroyMatrix(&Xk);

    //normalize kernel matrix
    if (traceNormalize) {
        float tmpktrace = 0.;
        for (int si = 0; si < insamples; si++)
            tmpktrace += Zk->sample[si].feat[si];

        *ktrace = tmpktrace;

        if (*ktrace != 0.0)
            iftMultDataSetByScalar(Zk, 1.0 / *ktrace);
    }

    return (Zk);
}


int *iftExtractSupportVectorIndices(iftSVM *svm, int idxModel, int *n) {
    *n = svm->model[idxModel]->l;
    int *sv_indices = iftAllocIntArray(*n);

    for(int sv = 0; sv < svm->model[idxModel]->l; sv++)
        sv_indices[sv] = svm->model[idxModel]->sv_indices[sv] - 1;

    return sv_indices;
}


int **iftExtractAllSupportVectorIndices(iftSVM *svm, int **n) {
    int *nSVs = iftAllocIntArray(svm->nmodels);

    int **sv_indices = (int**) iftAlloc(svm->nmodels, sizeof(int*));

    for (int m = 0; m < svm->nmodels; m++)
        sv_indices[m] = iftExtractSupportVectorIndices(svm, m, &nSVs[m]);

    *n = nSVs;

    return sv_indices;
}


