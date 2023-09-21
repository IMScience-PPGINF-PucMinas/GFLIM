#include "iftManifold.h"


#include "ift/core/io/Stream.h"
#include "tsne.h"


iftTsne* iftCreateTsne(iftDataSet* inputDataSet) {
    iftTsne* tsne = iftAlloc(1,sizeof(iftTsne));
    tsne->inputDataSet = inputDataSet;
    tsne->lowDimensionalSpace = 2;

    tsne->max_iter = 1000;
    tsne->learningRate = 200;
    tsne->perplexity = 30;
    if(inputDataSet != NULL){
        if(tsne->perplexity*3 > inputDataSet->nsamples-1){
            tsne->perplexity = (inputDataSet->nsamples-2)/3;
            if(tsne->perplexity < 0){
                tsne->perplexity = 0;
            }
        }
    }
    tsne->startMomentum = 0.5;
    tsne->endMomentum = 0.8;
    tsne->mom_switch_iter = tsne->max_iter*0.25;
    tsne->stop_lying_iter = tsne->max_iter*0.1;
    tsne->theta = 0.5;
    tsne->exageration = 20;
    tsne->normalizeProjection = true;

    tsne->highDimensionPoints = NULL;
    tsne->lowDimensionPoints = NULL;
    tsne->destroyInputDataset = false;
    return tsne;
}

void iftDestroyTsne(iftTsne** tsnep){
    if(tsnep == NULL || *tsnep == NULL){
        return;
    }
    iftTsne* aux = *tsnep;
    if(aux->inputDataSet != NULL){
        if(aux->destroyInputDataset == true){
            iftDestroyDataSet(&(aux->inputDataSet));
            (*tsnep)->inputDataSet = NULL;
        }
    }
    if(aux->lowDimensionPoints != NULL){
        free(aux->lowDimensionPoints);
        (*tsnep)->lowDimensionPoints = NULL;
    }
    if(aux->highDimensionPoints != NULL){
        free(aux->highDimensionPoints);
        (*tsnep)->highDimensionPoints = NULL;
    }
    free((*tsnep));
    (*tsnep) = NULL;
}

void iftComputeTsneProjection(iftTsne* tsne){
    if(tsne->inputDataSet == NULL){
        iftError("Can not compute projection. High-dimensional dataset is empty",
                 "iftComputeTsneProjection");
        return;
    }
    tsne->highDimensionalSpace = tsne->inputDataSet->nfeats;
    if(tsne->highDimensionalSpace < tsne->lowDimensionalSpace){
        iftError("Can not comput projection. Low-dimensional space is greater than high-dimensional space",
                 "iftComputeTsneProjection");
        return;
    }


    if(tsne->highDimensionPoints != NULL){
        free(tsne->highDimensionPoints);
    }

    if(tsne->inputDataSet->projection != NULL){
        iftDestroyDoubleMatrix(&(tsne->inputDataSet->projection));
    }

    if(tsne->highDimensionalSpace == tsne->lowDimensionalSpace){
        iftWarning("Input dataset dimension is equal to desired dimension","iftComputeTsneProjection");
        tsne->inputDataSet->projection = iftCreateDoubleMatrix(tsne->lowDimensionalSpace,tsne->inputDataSet->nsamples);
        tsne->lowDimensionPoints = (double*)iftAlloc(tsne->inputDataSet->nsamples*tsne->lowDimensionalSpace,sizeof(double));
        int i=0;
        for (int sampleIndex = 0; sampleIndex < tsne->inputDataSet->nsamples; ++sampleIndex) {
            for (int featIndex = 0; featIndex < tsne->inputDataSet->nfeats; ++featIndex) {
                tsne->lowDimensionPoints[i] = tsne->inputDataSet->sample[sampleIndex].feat[featIndex];
                i++;
            }
        }
    }
    else if(tsne->highDimensionalSpace > tsne->lowDimensionalSpace){
        tsne->highDimensionPoints = (double*)iftAlloc(tsne->inputDataSet->nsamples*tsne->highDimensionalSpace,sizeof(double));
        tsne->lowDimensionPoints = (double*)iftAlloc(tsne->inputDataSet->nsamples*tsne->lowDimensionalSpace,sizeof(double));
        size_t totalElementsLowDimension = tsne->inputDataSet->nsamples*tsne->lowDimensionalSpace;
        size_t totalElementsHighDimension = tsne->inputDataSet->nsamples*tsne->highDimensionalSpace;

        if(tsne->inputDataSet->data != NULL){
            for (int k = 0; k < totalElementsHighDimension; ++k) {
                tsne->highDimensionPoints[k] = tsne->inputDataSet->data->val[k];
            }
        }else if(tsne->inputDataSet->sample != NULL){
            int k =0;
            for (int i = 0; i < tsne->inputDataSet->nsamples; ++i) {
                for (int j = 0; j < tsne->inputDataSet->nfeats; ++j) {
                    tsne->highDimensionPoints[k] = tsne->inputDataSet->sample[i].feat[j];
                    k++;
                }
            }
        }
        // Initialize solution (randomly)
        for(int k = 0; k < totalElementsLowDimension; k++){
            tsne->lowDimensionPoints[k] = iftRandomNormalFloat(0.0, .0001);
        }
        run_tsne(tsne->highDimensionPoints, tsne->inputDataSet->nsamples, tsne->highDimensionalSpace,
                 tsne->lowDimensionPoints, tsne->lowDimensionalSpace, tsne->perplexity, tsne->theta,
                 tsne->max_iter, tsne->stop_lying_iter, tsne->mom_switch_iter, tsne->startMomentum,
                 tsne->endMomentum, tsne->learningRate);
        tsne->inputDataSet->projection = iftCreateDoubleMatrix(tsne->lowDimensionalSpace,tsne->inputDataSet->nsamples);
    }
    int k = 0;
    if(tsne->normalizeProjection){
        double *minimus = (double*)iftAlloc(tsne->lowDimensionalSpace,sizeof(double));
        double *maximus = (double*)iftAlloc(tsne->lowDimensionalSpace,sizeof(double));
        for (int j = 0; j < tsne->lowDimensionalSpace; ++j) {
            minimus[j] = tsne->lowDimensionPoints[j];
            maximus[j] = tsne->lowDimensionPoints[j];
        }
        k=0;
        for (int i = 0; i < tsne->inputDataSet->nsamples; ++i) {
            for (int j = 0; j < tsne->lowDimensionalSpace; ++j) {
                if(minimus[j] > tsne->lowDimensionPoints[k]){
                    minimus[j] = tsne->lowDimensionPoints[k];
                }
                if(maximus[j] < tsne->lowDimensionPoints[k]){
                    maximus[j] = tsne->lowDimensionPoints[k];
                }
                k++;
            }
        }
        k=0;
        for (int i = 0; i < tsne->inputDataSet->nsamples; ++i) {
            for (int j = 0; j < tsne->lowDimensionalSpace; ++j) {
                tsne->lowDimensionPoints[k] = (tsne->lowDimensionPoints[k] - minimus[j]) / (maximus[j] - minimus[j]);
                k++;
            }
        }
        iftFree(minimus);
        iftFree(maximus);
    }
    k = 0;
    for (int i = 0; i < tsne->inputDataSet->nsamples; ++i) {
        for (int j = 0; j < tsne->lowDimensionalSpace; ++j) {
            tsne->inputDataSet->projection->val[k] = tsne->lowDimensionPoints[k];
            k++;
        }
    }

}

iftMatrix *iftTSNEInit(unsigned int numberRows, unsigned int numberColumns, double mean, double variance){
    iftMatrix *randomSolution = iftCreateMatrix(numberColumns,numberRows);

    for (size_t i = 0; i < randomSolution->n; i++) {
        randomSolution->val[i] = iftRandomNormalFloat(mean, variance);
    }

    return randomSolution;
}

iftMatrix* iftRunTSNE(iftMatrix *P, int maxIter, int ndimension) {

    int n = P->nrows;
    int N;
    double mometum = 0.5;
    double final_mometum = 0.8;
    int mom_switch_iter = 20;
    int stop_lying_iter = 100;
    double epsilon = 500;
    double min_gain = 0.01;
    int iter,k,i,j,shift;
    double sum = 0;
    //double KL_constant;
    //double cost;

    iftMatrix *y_data = NULL;
    iftMatrix *y_data_transpose = NULL;
    iftMatrix *y_data_mul_transpose = NULL;
    iftMatrix *y_squared = NULL;
    iftMatrix *y_incs = NULL;
    iftMatrix * P_Log = NULL;
    iftMatrix *gains = NULL;
    iftMatrix *sum_ydata2 = NULL;

    iftMatrix *Q = iftCreateMatrix(n,n);
    iftMatrix *L = iftCreateMatrix(n,n);
    iftMatrix * diagonalMatrix = iftCreateMatrix(n,n);
    iftMatrix *y_grads = NULL;
    iftDblArray *Lsum = NULL;
    iftDblArray *Lmean = NULL;
    //iftMatrix * KL_matrix = NULL;
    iftMatrix * P_transpose = iftTransposeMatrix(P);


    iftSetDiagonalValue(P, 0);//make sure that all diagonal elements are zero
    iftMatricesAdditionPointWiseInPlace(P,P_transpose,&P);
    iftMultMatrixByScalar(P,0.5);
    sum = iftMatrixSum(P);
    iftMultMatrixByScalar(P,1.0/sum);
    P_Log = iftCopyMatrix(P);
    iftLogarithmMatrix(P_Log);
    //KL_matrix = iftComputeOperationBetweenMatricesInPointWise(P,P_transpose,'*');
    //KL_constant = iftMatrixSum(KL_matrix);
    iftMultMatrixByScalar(P,4);
    iftMatrixMinSaturation(P, 1e-12);

    y_data = iftTSNEInit(n, ndimension, 0.0, 1.0);
    y_squared = iftCreateMatrix(ndimension,n);
    y_data_transpose = iftCreateMatrix(ndimension,n);
    double *y_data_column_mean = (double *) iftAlloc(ndimension, sizeof(double));
    double * sum_ydataSquared = (double *) iftAlloc(n, sizeof(double));
    y_incs = iftCreateSetMatrix(n, ndimension, 0);
    gains = iftCreateSetMatrix(n, ndimension, 1);
    N = n*ndimension;

    for(iter=1; iter<=maxIter; iter++){

        /**Compute joint probability that point i and j are neighbours**/
        k = 0;
        y_data_transpose->nrows = y_data->nrows;
        y_data_transpose->ncols = y_data->ncols;
        for(i=0; i<y_data->nrows; i++){
            sum_ydataSquared[i] = 0;
            for(j=0; j<y_data->ncols; j++){
                y_data_transpose->val[k] = y_data->val[k];//copying data
                y_squared->val[k] = y_data->val[k]*y_data->val[k];//pointwise square
                sum_ydataSquared[i] += y_squared->val[k];
                k++;
            }
        }
        iftTransposeMatrixInPlace(y_data_transpose);
		iftMultMatricesInPlace(y_data, y_data_transpose, 0, 0, &y_data_mul_transpose);
        k=0;
        sum = 0;
        for(i=0; i<y_data_mul_transpose->nrows; i++){
            for(j=0; j<y_data_mul_transpose->ncols; j++){
                y_data_mul_transpose->val[k] = (-2*y_data_mul_transpose->val[k]) + sum_ydataSquared[i] + sum_ydataSquared[j];//distance
                y_data_mul_transpose->val[k] = 1/(1+y_data_mul_transpose->val[k]);//student-t distribution
                k++;
            }
        }
        //y_data_transpose = num
        iftSetDiagonalValue(y_data_mul_transpose,0);
        sum = iftMatrixSum(y_data_mul_transpose);
        k=0;
        for(i=0; i<y_data_mul_transpose->nrows; i++) {
            for (j = 0; j < y_data_mul_transpose->ncols; j++) {
                Q->val[k] = (y_data_mul_transpose->val[k]/sum);
                k++;
            }
        }
        iftMatrixMinSaturation(Q, 1E-12);
        /********************************************************/

        /*********Compute the gradients****************/
        k = 0;
        for(i=0; i<Q->nrows; i++) {
            for (j = 0; j < Q->ncols; j++) {
                L->val[k] = (P->val[k]-Q->val[k])*y_data_mul_transpose->val[k];
                k++;
            }
        }
        Lsum = iftComputeSumVector(L, 'c');
        for(i=0,shift=0; i<L->nrows; i++,shift+=L->ncols){
            diagonalMatrix->val[i+shift] = Lsum->val[i];
        }
        iftSubtractMatricesInPlace(diagonalMatrix,L,L);
		iftMultMatricesInPlace(L, y_data, 0, 0, &y_grads);
        iftMultMatrixByScalar(y_grads,4);
        /********************************************/

        /*************Update solution****************/
        for(k=0; k<N ;k++){
            if((y_grads->val[k] > 0 && y_incs->val[k] > 0) || (y_grads->val[k] < 0 && y_incs->val[k] < 0) ){
                gains->val[k] = 0.8 * gains->val[k];
            }else{
                gains->val[k] = 0.2 + gains->val[k];
            }
        }
        iftMatrixMinSaturation(gains, min_gain);

        k = 0;
        for(i=0; i< y_data->nrows; i++) {
            for (j = 0; j < y_data->ncols; j++) {
                y_incs->val[k] = (mometum * y_incs->val[k]) - (epsilon * gains->val[k] * y_grads->val[k]);
                y_data->val[k] += y_incs->val[k];
                k++;
            }
        }
        for (j = 0; j < y_data->ncols; j++) {
            y_data_column_mean[j] = 0;
            for (i = 0; i < y_data->nrows; i++) {
                y_data_column_mean[j] += y_data->val[iftGetMatrixIndex(y_data,j,i)];
            }
            y_data_column_mean[j] /= y_data->nrows;
            for (i = 0; i < y_data->nrows; i++) {
                y_data->val[iftGetMatrixIndex(y_data,j,i)] -= y_data_column_mean[j];
            }
        }
        /*****************************/

        if(iter == mom_switch_iter){
            mometum = final_mometum;

        }
        if (iter == stop_lying_iter){
            iftMultMatrixByScalar(P,1.0/4);
        }


        float C = 0;

        for (int i = 0; i < P->n; ++i) {
//			printf("%f, %f\n", P->val[i], Q->val[i]);
            C += P->val[i] * log(P->val[i] / Q->val[i]);
        }

        printf("C = %f\n", C);

        iftDestroyDblArray(&Lsum);

    }



    iftDestroyMatrix(&diagonalMatrix);
    iftDestroyDblArray(&Lsum);
    iftDestroyDblArray(&Lmean);
    iftDestroyMatrix(&y_data_mul_transpose);
    iftDestroyMatrix(&y_squared);
    iftDestroyMatrix(&y_incs);
    iftDestroyMatrix(&P_Log);
    iftDestroyMatrix(&P_transpose);
    iftDestroyMatrix(&Q);
    iftDestroyMatrix(&L);
    iftDestroyMatrix(&gains);
    iftDestroyMatrix(&sum_ydata2);
    iftDestroyMatrix(&y_grads);
    iftDestroyMatrix(&y_data_transpose);
    free(sum_ydataSquared);
    free(y_data_column_mean);

    return y_data;

}

double iftHbeta(iftFloatArray *D, float beta,iftFloatArray *P){
    int i;
    double sum_P = 0;
    double sum_DP = 0;
    iftFloatArray *DP = iftCreateFloatArray(P->n);


    for (i=0; i<P->n; i++){
        P->val[i] = exp(-D->val[i]*beta);
        DP->val[i] = P->val[i] * D->val[i];
    }

    sum_P = iftSumFloatArray(P->val, P->n);
    sum_DP = iftSumFloatArray(DP->val, DP->n);

    if(sum_P < IFT_EPSILON){
        sum_P = IFT_EPSILON;
    }

    iftMultScalarFloatArray(P->val, P->n, 1.0/sum_P);
    iftDestroyFloatArray(&DP);

    return log(sum_P) + (beta*(sum_DP/sum_P));
}

void iftDistToProb(iftDistanceTable *D, double u, double tol, iftMatrix **Pp, iftFloatArray **beta){
    int n = D->nsamples1;
    iftMatrix* P = iftCreateMatrix(n, n);
    iftFloatArray *thisP = iftCreateFloatArray(n);
    iftFloatArray* distances =  iftCreateFloatArray(n);
    (*beta) = iftCreateFloatArray(n);
    double logU = log(u);
    double betaMin;
    double betaMax;
    double H,Hdiff;
    int i,j,tries;

    //doesnt exist in the original.
    //However, without normalization it can cause NAN values.
    //In Peixinho, you should trust.
    float dmax = IFT_INFINITY_FLT_NEG;

    for (int i = 0; i < D->nsamples1; ++i) {
        for (int j = 0; j < D->nsamples2; ++j) {
            if(dmax<D->distance_table[i][j]) {
                dmax = D->distance_table[i][j];
            }
        }
    }

    for (i=0; i<n; i++){
        (*beta)->val[i] = 1.0;
    }
    for (i=0; i<n; i++){

        for(j=0; j<n; j++){
            if(i == j){
                distances->val[j] = IFT_INFINITY_FLT;
            }else{
                distances->val[j] = D->distance_table[i][j]/dmax;
            }
        }

        betaMin = IFT_INFINITY_FLT_NEG;
        betaMax = IFT_INFINITY_FLT;
        H = iftHbeta(distances,(*beta)->val[i],thisP);
        Hdiff = H - logU;
        tries = 0;
        while(fabs(Hdiff) > tol && tries < 50){
            if(Hdiff > 0){
                betaMin = (*beta)->val[i];
                if(betaMax == IFT_INFINITY_FLT){
                    (*beta)->val[i] = (*beta)->val[i]*2.0;
                }else{
                    (*beta)->val[i] = ((*beta)->val[i] + betaMax)/2.0;
                }
            }
            else{
                betaMax = (*beta)->val[i];
                if(betaMin == IFT_INFINITY_FLT_NEG){
                    (*beta)->val[i] = (*beta)->val[i]/2.0;
                }else{
                    (*beta)->val[i] = ((*beta)->val[i] + betaMin)/2.0;
                }
            }

            //recompute the values
            H = iftHbeta(distances,(*beta)->val[i],thisP);
            Hdiff = H - logU;
            tries++;
        }
        for(j=0; j<n; j++){
            iftMatrixElem(P,j,i) = thisP->val[j];
        }
    }

    *Pp = P;

    iftDestroyFloatArray(&thisP);
    iftDestroyFloatArray(&distances);
}

iftMatrix *iftTSNE(iftDistanceTable *dist, unsigned int ndim, double perplexity, int maxIter) {

    iftMatrix * dataOut;
    iftFloatArray* beta;
    iftMatrix *prob;

    iftDistToProb(dist, perplexity, 0.00001, &prob, &beta);

    dataOut = iftRunTSNE(prob, maxIter, ndim);

    iftDestroyMatrix(&prob);
    iftDestroyFloatArray(&beta);

    return dataOut;
}

iftDataSet* iftDimReductionByTSNE(iftDataSet* Z, int ndim, double perplexity, size_t max_iter) {

    iftDataSet* outputDataset = NULL;
    if(Z == NULL){
        iftError("Dataset is NULL", "iftDataSetFromTSNEProjection");
        return NULL;
    }

    iftTsne* tsne             = iftCreateTsne(Z);
    tsne->lowDimensionalSpace = ndim;
    tsne->normalizeProjection = true;
    tsne->perplexity          = perplexity;
    tsne->max_iter            = max_iter;
    iftComputeTsneProjection(tsne);
    outputDataset             = iftCreateDataSet(Z->nsamples, ndim);

    outputDataset->ntrainsamples = Z->ntrainsamples;
    outputDataset->ngroups       = Z->ngroups;
    outputDataset->nclasses      = Z->nclasses;
    if(Z->ref_data != NULL)
        iftCopyRefData(outputDataset, Z->ref_data, Z->ref_data_type);

    for (int s = 0; s < Z->nsamples; s++) {
        outputDataset->sample[s].id        = Z->sample[s].id;
        outputDataset->sample[s].truelabel = Z->sample[s].truelabel;
        outputDataset->sample[s].label     = Z->sample[s].label;
        outputDataset->sample[s].group     = Z->sample[s].group;
        iftSetSampleStatus(&outputDataset->sample[s], Z->sample[s].status);
        outputDataset->sample[s].weight    = Z->sample[s].weight;
        outputDataset->sample[s].feat      = &(iftMatrixElem(outputDataset->data, 0, s));
    }
    
    for (int i = 0; i < Z->projection->n; ++i) {
      outputDataset->data->val[i] = Z->projection->val[i];
    }
    
    iftDestroyTsne(&tsne);
    return outputDataset;
}

iftImage* iftEstimateClassifierBoundariesDecision(iftDataSet* dataSet,iftGenericClassifier* classifier){
    iftWarning("The color table scheme was altered. Check if this function keeps working. Afterwards, remove this warning.", "iftEstimateClassifierBoundariesDecision");

    if(dataSet == NULL){
        return NULL;
    }
    if(dataSet->nsamples <= 0){
        return NULL;
    }
    float* minimus = iftAlloc(dataSet->nfeats,sizeof(float));
    float* maximus = iftAlloc(dataSet->nfeats,sizeof(float));
    int* sampleindexMin = iftAlloc(dataSet->nfeats,sizeof(float));
    int* sampleIndexMax = iftAlloc(dataSet->nfeats,sizeof(float));
    for (int k = 0; k < dataSet->nfeats; ++k) {
        minimus[k] = dataSet->sample[0].feat[k];
        maximus[k] = dataSet->sample[0].feat[k];
        sampleindexMin[k] = 0;
        sampleIndexMax[k] = 0;
    }
    for (int i = 0; i < dataSet->nsamples; ++i) {
        for (int j = 0; j < dataSet->nfeats; ++j) {
            if(minimus[j] > dataSet->sample[i].feat[j]){
                minimus[j] = dataSet->sample[i].feat[j];
                sampleindexMin[j] = i;
            }
            if(maximus[j] < dataSet->sample[i].feat[j]){
                maximus[j] = dataSet->sample[i].feat[j];
                sampleIndexMax[j] = i;
            }
        }
    }
    float* deltaVector = iftAlloc(dataSet->nfeats,sizeof(float));
    int nHIntervals = 256;
    //int nVIntervals = 10;
    iftColorTable* colorTable = iftCreateColorTable(dataSet->nclasses+1);
    iftImage* image = iftCreateImage((nHIntervals+1),dataSet->nfeats,1);
    iftSetCbCr(image,128);
    //int voxelIndex = 0;
    int currentFeat = 0;
    int currentX = 0;
    int currentY = 0;
    while(currentFeat < dataSet->nfeats){
        currentX = 0;
        iftDataSet* dataSetSamplesRow = iftCreateDataSet(nHIntervals+1,dataSet->nfeats);
        for (int i = 0; i < dataSet->nfeats; ++i) {
            deltaVector[i] = (dataSet->sample[sampleIndexMax[currentFeat]].feat[i] - dataSet->sample[sampleindexMin[currentFeat]].feat[i])/nHIntervals;
        }
        for (int j = 0; j <= nHIntervals; ++j) {
            for (int i = 0; i < dataSet->nfeats; ++i) {
                dataSetSamplesRow->sample[j].feat[i] = dataSet->sample[sampleindexMin[currentFeat]].feat[i] + j*deltaVector[i];
            }
        }
        iftPredictGenericClassifier(classifier,dataSetSamplesRow);
        for (int j = 0; j <= nHIntervals ; ++j) {
            int label = dataSetSamplesRow->sample[j].label;
            iftImgVal2D(image,currentX,currentY) = colorTable->color[label].val[0];
            iftImgCb2D(image,currentX,currentY) = colorTable->color[label].val[1];
            iftImgCr2D(image,currentX,currentY) = colorTable->color[label].val[2];
            currentX++;
        }
        currentY += 1;
        iftDestroyDataSet(&dataSetSamplesRow);
        currentFeat++;
    }
    iftConvertRGBImagetoYCbCrImage(image,255);
    iftWriteImageByExt(image,"a.ppm");
    iftImage *imageColorTable = iftColorTableToImage(colorTable, 20,200);
    iftWriteImageByExt(imageColorTable,"b.ppm");
    return image;
}

iftImage* iftEstimateClassifierBoundariesDecision2(iftDataSet* dataSet,iftGenericClassifier* classifier){
    float* minimus = iftAlloc(dataSet->nfeats,sizeof(float));
    float* maximus = iftAlloc(dataSet->nfeats,sizeof(float));
    int* sampleIndexMin = iftAlloc(dataSet->nfeats,sizeof(float));
    int* sampleIndexMax = iftAlloc(dataSet->nfeats,sizeof(float));
    for (int k = 0; k < dataSet->nfeats; ++k) {
        minimus[k] = dataSet->sample[0].feat[k];
        maximus[k] = dataSet->sample[0].feat[k];
        sampleIndexMin[k] = 0;
        sampleIndexMax[k] = 0;
    }
    for (int i = 0; i < dataSet->nsamples; ++i) {
        for (int j = 0; j < dataSet->nfeats; ++j) {
            if(minimus[j] > dataSet->sample[i].feat[j]){
                minimus[j] = dataSet->sample[i].feat[j];
                sampleIndexMin[j] = i;
            }
            if(maximus[j] < dataSet->sample[i].feat[j]){
                maximus[j] = dataSet->sample[i].feat[j];
                sampleIndexMax[j] = i;
            }
        }
    }
    int spaceBetweenLines = 10;
    size_t verticalSize = dataSet->nfeats*spaceBetweenLines;
    size_t horizontalSize = verticalSize/2;
    size_t actualHorizontalSize = horizontalSize+1;
    iftDataSet* gridSamples = iftCreateDataSet(actualHorizontalSize*verticalSize,dataSet->nfeats);
    gridSamples->nclasses = dataSet->nclasses;
    float* deltaVector = iftAlloc(dataSet->nfeats,sizeof(float));

    int currentFeat = 0;
    int currentRow;
    size_t sampleIndexInGrid;
    size_t shift;
    //place min/maxs on matrix
    while(currentFeat < dataSet->nfeats){
        currentRow = currentFeat*spaceBetweenLines;
        int lastColumnIndex = actualHorizontalSize-1;
        shift = currentRow*actualHorizontalSize;
        sampleIndexInGrid = shift;
        iftSample* sampleMin  = &(gridSamples->sample[sampleIndexInGrid]);
        sampleIndexInGrid = shift+lastColumnIndex;
        iftSample* sampleMax = &(gridSamples->sample[sampleIndexInGrid]);
        int indexMin = sampleIndexMin[currentFeat];
        int indexMax = sampleIndexMax[currentFeat];
        iftSample* sampleMinSrc = &(dataSet->sample[indexMin]);
        iftSample* sampleMaxSrc = &(dataSet->sample[indexMax]);
        sampleMin->truelabel = sampleMinSrc->truelabel;
        sampleMin->group = sampleMinSrc->group;
        sampleMin->isSupervised = true;

        sampleMax->truelabel = sampleMaxSrc->truelabel;
        sampleMax->group = sampleMaxSrc->group;
        sampleMax->isSupervised = true;

        for (int i = 0; i < dataSet->nfeats; ++i) {
            sampleMin->feat[i] = sampleMinSrc->feat[i];
            sampleMax->feat[i] = sampleMaxSrc->feat[i];
        }
        currentFeat++;
    }
    iftFree(sampleIndexMin);
    iftFree(sampleIndexMax);
    iftFree(minimus);
    iftFree(maximus);

    currentFeat = 0;
    //placing horizontal samples
    while(currentFeat < dataSet->nfeats){
        currentRow = currentFeat*spaceBetweenLines;
        int lastColumnIndex = actualHorizontalSize-1;
        shift = currentRow*actualHorizontalSize;
        sampleIndexInGrid = shift;
        iftSample* sampleMin  = &(gridSamples->sample[sampleIndexInGrid]);
        sampleIndexInGrid = shift+lastColumnIndex;
        iftSample* sampleMax = &(gridSamples->sample[sampleIndexInGrid]);
        for (int i = 0; i < dataSet->nfeats; ++i) {
            deltaVector[i] = (sampleMax->feat[i] - sampleMin->feat[i])/horizontalSize;
        }
        for (int j = 0; j <= horizontalSize; ++j) {
            sampleIndexInGrid = shift+j;
            iftSample* currentSample = &(gridSamples->sample[sampleIndexInGrid]);
            if(currentSample->isSupervised == true){
                continue;
            }
            for (int i = 0; i < dataSet->nfeats; ++i) {
                currentSample->feat[i] = sampleMin->feat[i] + j*deltaVector[i];
            }
            currentSample->isSupervised = true;
        }
        currentFeat++;
    }

    currentFeat = 0;
    while(currentFeat < dataSet->nfeats){
        int indexRow1 = currentFeat*spaceBetweenLines;
        int indexRow2 = (currentFeat+1)*spaceBetweenLines;
        indexRow2 = indexRow2%verticalSize;
        size_t sampleIndexRow1 = indexRow1*actualHorizontalSize;
        size_t sampleIndexRow2 = indexRow2*actualHorizontalSize;
        for (int j = 0; j <= horizontalSize; ++j) {
            sampleIndexInGrid = sampleIndexRow1+j;
            iftSample* sampleUp = &(gridSamples->sample[sampleIndexInGrid]);
            sampleIndexInGrid = sampleIndexRow2+j;
            iftSample* sampleDown = &(gridSamples->sample[sampleIndexInGrid]);
            for (int i = 0; i < dataSet->nfeats; ++i) {
                deltaVector[i] = (sampleDown->feat[i] - sampleUp->feat[i])/spaceBetweenLines;
            }
            for (int k = 0; k < spaceBetweenLines; ++k) {
                int sampleRow = indexRow1+k;
                sampleIndexInGrid = (sampleRow*actualHorizontalSize)+j;
                iftSample* currentSample = &(gridSamples->sample[sampleIndexInGrid]);
                if(currentSample->isSupervised == true){
                    continue;
                }
                for (int i = 0; i < dataSet->nfeats; ++i) {
                    currentSample->feat[i] = sampleUp->feat[i] + k*deltaVector[i];
                }
                currentSample->isSupervised = true;
            }
        }
        currentFeat++;
    }
    iftFree(deltaVector);
    iftPredictGenericClassifier(classifier,gridSamples);
    iftColorTable* colorTable = iftCreateColorTable(dataSet->nclasses+1);
    iftImage* image = iftCreateImage(actualHorizontalSize,verticalSize,1);
    iftSetCbCr(image,128);
    for (int sampleIndex = 0; sampleIndex < gridSamples->nsamples; ++sampleIndex) {
        int label = gridSamples->sample[sampleIndex].label;
        image->val[sampleIndex] = colorTable->color[label].val[0];
        image->Cb[sampleIndex] = colorTable->color[label].val[1];
        image->Cr[sampleIndex] = colorTable->color[label].val[2];
    }
    iftDestroyDataSet(&gridSamples);
    iftDestroyColorTable(&colorTable);
//    iftConvertRGBImagetoYCbCrImage(image,255);
//    iftDestroyDataSet(&dataSetSample);
//    iftDestroyGenericMatrix(&sampleMatrix);
//    iftDestroyColorTable(&colorTable);
//    iftDestroyGenericMatrix(&sampleMatrix);
    return image;
}
