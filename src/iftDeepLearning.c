#include "iftDeepLearning.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Json.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/tools/String.h"


/*------------------------------------- PRIVATE -----------------------------------*/

iftDataSet     *iftSelectPatches(iftMImage *img, iftConvNetwork *convnet, int nsamples);
void            iftCentralizePatch(const float *mean, iftMatrix *M);
iftMatrix      *iftMImageIndicesToMatrix(iftMImage *img, iftFastAdjRel *F, int stride);

iftMatrix      *iftGetVoxelValueMatrix(iftMatrix *img_ind, iftMImage *img);
void            iftActivateByBandThres(iftMImage *img);
void            iftActivateByThres(iftMImage *img, float activ_thres);
iftMImage      *iftNormalizeContrastMImage(iftMImage *img, iftAdjRel *A);
iftMImage      *iftNormalizeMImage(iftMImage *img, iftAdjRel *A);
iftMImage      *iftPooling(iftMImage *img, iftAdjRel *A, int stride, float alpha);
iftMImage      *iftMMLinearFilterConvNetwork(iftMImage *input_img, iftMatrix *input_ind, iftMMKernel *k_bank, int output_xsize, int output_ysize, int output_zsize, int output_nbands);
iftConvNetwork *iftReadConvNetworkFILE(FILE* fp);
void            iftWriteConvNetworkFILE(iftConvNetwork *convnet, FILE *fp);

iftDataSet *iftSelectPatches(iftMImage *img, iftConvNetwork *convnet, int nsamples)
{
	iftDataSet *Z;
	int dx, dy, dz, ntrials, max_ntrials;
	iftMImage *nimg;
	iftAdjRel *A;

	if (convnet->input_norm_adj_param > 0)
		nimg = iftNormalizeMImage(img,convnet->input_norm_adj);
	else
		nimg = img;

	A = convnet->k_bank[0]->A;

	iftMaxAdjShifts(A, &dx, &dy, &dz);

	int *count  = iftAllocIntArray(nimg->n);
	int *sample = iftAllocIntArray(nimg->n);

	int high=0;
	for (int p=0; p < nimg->n; p++) {
		iftVoxel u = iftMGetVoxelCoord(nimg,p);
		if ((u.x >= dx)&&(u.x < nimg->xsize-dx)&&
		    (u.y >= dy)&&(u.y < nimg->ysize-dy)&&
		    (u.z >= dz)&&(u.z < nimg->zsize-dz)){
		  count[high]  = 100;
		  sample[high] = p;
		  high++;
		}
	}
	high--;
	int last = high;

	/* Select samples based on their probability value and with minimum
	   overlapping between regions. Stop whenever the number of trials
	   is prob->maxval * prob->n to avoid infinite loop due to the
	   region overlapping constraint. */


	int s = 0; ntrials = 0; max_ntrials = 100 * high;

	while ((s < nsamples)&&(ntrials < max_ntrials)) {
		int index = iftRandomInteger(0,high);
		if (count[index] == 0) {
			iftSwap(sample[index], sample[high]);
			iftSwap(count[index], count[high]);
			high--;
		} else if (count[index]==1){
			int p           = sample[index];
			iftVoxel u      = iftMGetVoxelCoord(nimg,p);
			ntrials++;

			/* Select p only if there is no other selected voxel within the
		   same adjacency region of p */

			char select_p   = 1;

			for (int i=high+1; i <= last; i++) { /* region overlapping
						constraint */
				int q = sample[i];
				iftVoxel v = iftMGetVoxelCoord(nimg,q);
				if ( (abs(v.x-u.x)<=dx/2) &&
					 (abs(v.y-u.y)<=dy/2) &&
					 (abs(v.z-u.z)<=dz/2)    ){
					select_p = 0;
					break;
				}
			}
			if (select_p) {	/* perform selection */
				iftSwap(sample[index],sample[high]);
				iftSwap(count[index], count[high]);
				high--; s++;
			}
		}else{
			count[index]--;
		}
	}

	/* Copy selected patches to dataset */

	if (s != nsamples) {
		char msg[200];
		sprintf(msg,"Could select %d from %d required samples\n",s,nsamples);
		iftWarning(msg,"iftSelectPatches");
	}


	nsamples = s;
	Z=iftCreateDataSet(nsamples,A->n*nimg->m);
	iftCopyRefData(Z, img, IFT_REF_DATA_IMAGE);
	// Z->ref_data = img;

	s = 0;
	for (int index=high+1; index <= last; index++, s++) {
		int p           = sample[index];

		iftVoxel u      = iftMGetVoxelCoord(nimg,p);

		int j           = 0;
		Z->sample[s].id = p;
		for (int b=0; b < nimg->m; b++) {
			for (int i=0; i < A->n; i++) {
				iftVoxel v = iftGetAdjacentVoxel(A,u,i);
				if (iftMValidVoxel(nimg,v)){
					int    q = iftMGetVoxelIndex(nimg,v);
					Z->sample[s].feat[j] = nimg->val[q][b];
				}
				j++;
			}
		}
	}

	iftSetStatus(Z, IFT_TRAIN);

	if (nimg != img)
		iftDestroyMImage(&nimg);

	iftFree(count);
	iftFree(sample);

	return(Z);
}




void iftCentralizePatch(const float *mean, iftMatrix *M)
{
	int N = (M->nrows);

	for(int c2=0; c2 < M->ncols; c2++) {
		float meanL = 0.0f;
		for(int r2=0; r2 < M->nrows; r2++) {
			long index = iftGetMatrixIndex(M,c2,r2);
			meanL += M->val[index];
		}
		meanL /= N;

		float stdevL = 0.0f;
		for(int r2=0; r2 < M->nrows; r2++) {
			long index = iftGetMatrixIndex(M,c2,r2);
			M->val[index] = M->val[index] - meanL;
			stdevL = M->val[index] * M->val[index];
		}

		stdevL = sqrtf( (stdevL + 0.1f) / (N-1) ); // eps considering data in [0,1]
		for(int r2=0; r2 < M->nrows; r2++) {
			long index = iftGetMatrixIndex(M,c2,r2);
			M->val[index] = M->val[index] / stdevL;
		}
	}

	/* For every image patch, centralize it */
	for(int r2=0; r2 < M->nrows; r2++) {
		for(int c2=0; c2 < M->ncols; c2++) {
			long index = iftGetMatrixIndex(M,c2,r2);
			M->val[index] = (M->val[index]-mean[r2]);
		}
	}

}

iftMatrix *iftGetVoxelValueMatrix(iftMatrix *img_ind, iftMImage *img)
{
	iftMatrix *img_val = iftCreateMatrix(img_ind->ncols,img_ind->nrows);

#pragma omp parallel for shared(img_ind,img_val)
	for (int i=0; i < img_ind->n; i++) {
		int q = ((int)img_ind->val[i])%img->n;
		int m = ((int)img_ind->val[i])/img->n;
		img_val->val[i] = img->val[q][m];
	}
	return(img_val);
}

iftMatrix *iftMImageIndicesToMatrix(iftMImage *img, iftFastAdjRel *F, int stride)
{
	iftMatrix *M;
	int stride_z=1, zsize=img->zsize-2*F->bz;
	int xsize = ceilf(((float)(img->xsize - 2*F->bx))/stride);
	int ysize = ceilf(((float)(img->ysize - 2*F->by))/stride);

	if (img->zsize != 1){ // 3D
		zsize = ceilf(((float)zsize)/stride);
		stride_z  = stride;
	}

	M = iftCreateMatrix(xsize*ysize*zsize,F->n*img->m);

#pragma omp parallel for shared(img,M,F,stride_z,stride)
	for (int m=0; m < img->m; m++) {
		int offset_band = m * img->n;
		for (int i=0; i < F->n; i++) {
			iftVoxel u;
			int r = i+m*F->n;
			int c = 0;
			for (u.z = F->bz; u.z < img->zsize-F->bz; u.z += stride_z)
				for (u.y = F->by; u.y < img->ysize-F->by; u.y += stride)
					for (u.x = F->bx; u.x < img->xsize-F->bx; u.x += stride){
						int p = iftGetVoxelIndex(img,u);
						int q = p + F->dq[i];
						M->val[iftGetMatrixIndex(M,c,r)] = q + offset_band;
						c++;
					}
		}
	}

	return(M);
}

void iftActivateByBandThres(iftMImage *img)
{
	double *mean=iftAllocDoubleArray(img->m);

#pragma omp parallel for shared(img,mean)
	for (int b=0; b < img->m; b++) {
		for (int p=0; p < img->n; p++){
			mean[b] += img->val[p][b];
			//      mean[b] += fabs(img->val[p][b]);
		}
	}

	for (int b=0; b < img->m; b++) {
		mean[b] /= img->n;
		//    if (mean[b] < 0.0)
		mean[b] = 0.0;
	}

#pragma omp parallel for shared(img)
	for (int b=0; b < img->m; b++) {
		for (int p=0; p < img->n; p++){
			img->val[p][b] -= mean[b];
			if (img->val[p][b] < 0.)
				img->val[p][b] = 0.;
		}
	}

	iftFree(mean);
}


void iftActivateByThres(iftMImage *img, float activ_thres) {
#pragma omp parallel for shared(img, activ_thres)
	for (int b = 0; b < img->m; b++)
		for (int p = 0; p < img->n; p++) {
			img->val[p][b] -= activ_thres;
			if (img->val[p][b] < 0)
				img->val[p][b] = 0.0;
		}
}


iftMImage *iftNormalizeContrastMImage(iftMImage *img, iftAdjRel *A) {
	iftFastAdjRel *F  = iftCreateFastAdjRel(A, img->tby, img->tbz);
	iftMImage* output = iftCreateMImage(img->xsize - 2 * F->bx,
										img->ysize - 2 * F->by,
										img->zsize - 2 * F->bz,
										img->m);

	int N = img->m*F->n;

#pragma omp parallel for shared(F, img, output)
	for (int y = F->by; y < img->ysize - F->by; y++)
		for (int x = F->bx; x < img->xsize - F->bx; x++)
			for (int z = F->bz; z < img->zsize - F->bz; z++){
				int p = x + img->tby[y] + img->tbz[z];
				int output_p = (x - F->bx) + output->tby[y - F->by]
							   + output->tbz[z - F->bz];

				double mean   = 0.0;
				for (int b = 0; b < img->m; b++) {
					int i;
					for (i = 0; i < F->n; i++) {
						int q = p + F->dq[i];
						mean  += img->val[q][b];
					}
				}
				mean /= N;

				double stdev  = 0.0;
				for (int b = 0; b < img->m; b++) {
					int i;
					for (i = 0; i < F->n; i++) {
						int q = p + F->dq[i];
						stdev  += ((float)img->val[q][b]-mean)*((float)img->val[q][b]-mean);
					}
				}
				stdev = sqrtf( stdev/(N-1) );

				if (stdev > IFT_EPSILON)
					for (int b = 0; b < img->m; b++) {
						output->val[output_p][b] = ((float) img->val[p][b] - mean)/stdev;
					}

			}


	iftDestroyFastAdjRel(&F);

	return (output);
}

iftMImage *iftNormalizeMImage(iftMImage *img, iftAdjRel *A) {
	iftFastAdjRel *F = iftCreateFastAdjRel(A, img->tby, img->tbz);
	iftMImage* output = iftCreateMImage(img->xsize - 2 * F->bx,
										img->ysize - 2 * F->by,
										img->zsize - 2 * F->bz,
										img->m);


#pragma omp parallel for shared(F, output)
	for (int y = F->by; y < img->ysize - F->by; y++)
		for (int x = F->bx; x < img->xsize - F->bx; x++)
			for (int z = F->bz; z < img->zsize - F->bz; z++) {
				int p = x + img->tby[y] + img->tbz[z];
				// index of the resulting pixel with the normalized value
				int output_p = (x - F->bx) + output->tby[y - F->by]
							   + output->tbz[z - F->bz];
				double sum = 0.0;

				int b;
				for (b = 0; b < img->m; b++) {
					int i;
					for (i = 0; i < F->n; i++) {
						int q = p + F->dq[i];
						sum += (img->val[q][b] * img->val[q][b]);
					}
				}

				sum = sqrtf(sum);

				if (sum > IFT_EPSILON) {
					for (b = 0; b < img->m; b++) {
						output->val[output_p][b] = (float) img->val[p][b] / sum;
					}
				}
			}

	iftDestroyFastAdjRel(&F);

	return (output);
}

iftMImage *iftPooling(iftMImage *img, iftAdjRel *A, int stride, float alpha){
	iftFastAdjRel *F = iftCreateFastAdjRel(A, img->tby, img->tbz);

	int xsize = ceilf(((float)(img->xsize - 2*F->bx))/stride);
	int ysize = ceilf(((float)(img->ysize - 2*F->by))/stride);
	int zsize = ceilf(((float)(img->zsize - 2*F->bz))/stride);

	iftMImage* output_img = iftCreateMImage(xsize, ysize, zsize, img->m);

	int b;

	if (alpha < ALPHA_LIMIT) {
#pragma omp parallel for shared(img, output_img, A, alpha, F, xsize, ysize, zsize)
		for(b = 0; b < img->m; b++){
			iftVoxel v;
			for(v.z = F->bz; v.z < img->zsize-F->bz; v.z += stride){
				for(v.y = F->by; v.y < img->ysize-F->by; v.y += stride){
					for(v.x = F->bx; v.x < img->xsize-F->bx; v.x += stride){
						int p = v.x + img->tby[v.y] + img->tbz[v.z];
						double sum = 0.0;

						int i;
						for(i = 0; i < F->n; i++){
							int q = p + F->dq[i];
							sum += pow(img->val[q][b], alpha);
						}

						iftVoxel outputv;
						outputv.x = (v.x - F->bx)/stride;
						outputv.y = (v.y - F->by)/stride;
						outputv.z = (v.z - F->bz)/stride;
						int outputp = iftMGetVoxelIndex(output_img, outputv);
						output_img->val[outputp][b] = (float)pow(sum, 1/alpha);
					}
				}
			}
		}
	} else { /* MAX POOLING */
#pragma omp parallel for shared(img, output_img, A, alpha, F, xsize, ysize, zsize)
		for(b = 0; b < img->m; b++){
			iftVoxel v;
			for(v.z = F->bz; v.z < img->zsize-F->bz; v.z += stride){
				for(v.y = F->by; v.y < img->ysize-F->by; v.y += stride){
					for(v.x = F->bx; v.x < img->xsize-F->bx; v.x += stride){
						int p = v.x + img->tby[v.y] + img->tbz[v.z];
						float max = IFT_INFINITY_FLT_NEG;

						int i;
						for(i = 0; i < F->n; i++){
							int q = p + F->dq[i];
							if (img->val[q][b] > max) max = img->val[q][b];
						}

						iftVoxel outputv;
						outputv.x = (v.x - F->bx)/stride;
						outputv.y = (v.y - F->by)/stride;
						outputv.z = (v.z - F->bz)/stride;
						int outputp = iftMGetVoxelIndex(output_img, outputv);
						output_img->val[outputp][b] = max;
					}
				}
			}
		}
	}

	iftDestroyFastAdjRel(&F);

	return output_img;
}


iftMImage *iftMMLinearFilterConvNetwork(iftMImage *input_img, iftMatrix *input_ind, iftMMKernel *k_bank, int output_xsize, int output_ysize, int output_zsize, int output_nbands)
{
	iftMatrix   *A, *B, *M;
	iftMImage   *fimg;

	A    = iftMMKernelToMatrix(k_bank);
	B    = iftGetVoxelValueMatrix(input_ind, input_img);

	if (k_bank->W != NULL){
		iftCentralizePatch(k_bank->mean,B);
		iftMatrix *A_aux = iftMultMatrices(A,k_bank->W);
		iftDestroyMatrix(&A);
		A = A_aux;
	}

	M    = iftMultMatrices(A,B);

	// subtract bias in each band from resulting MImage
	for (int r = 0; r < M->nrows; r++) {
		for (int c = 0; c < M->ncols; c++)
			M->val[iftGetMatrixIndex(M,c,r)] -= k_bank->bias[r];
	}

	fimg = iftMatrixToMImage(M,output_xsize, output_ysize, output_zsize, output_nbands, 'r');

	iftDestroyMatrix(&A);
	iftDestroyMatrix(&B);
	iftDestroyMatrix(&M);

	return(fimg);
}



void  iftWriteConvNetworkFILE(iftConvNetwork *convnet, FILE *fp)
{

	fprintf(fp,"CONVNET\n");

	/* Write main parameters */

	iftWriteIntValue(fp,convnet->nlayers,"NLAYERS");
	iftWriteIntValue(fp,convnet->input_norm_adj_param,"INPUT_NORM_ADJ_PARAM");
	iftWriteIntValue(fp,convnet->input_xsize,"INPUT_XSIZE");
	iftWriteIntValue(fp,convnet->input_ysize,"INPUT_YSIZE");
	iftWriteIntValue(fp,convnet->input_zsize,"INPUT_ZSIZE");
	iftWriteIntValue(fp,convnet->input_nbands,"INPUT_NBANDS");
	iftWriteIntValues(fp,convnet->k_bank_adj_param,convnet->nlayers,"K_BANK_ADJ_PARAM");
	iftWriteIntValues(fp,convnet->nkernels,convnet->nlayers,"NKERNELS");
	iftWriteIntValues(fp, convnet->activ_option, convnet->nlayers, "ACTIV_OPTION");
	iftWriteFloatValues(fp, convnet->activ_thres, convnet->nlayers, "ACTIV_THRES");
	iftWriteIntValues(fp,convnet->pooling_adj_param,convnet->nlayers,"POOLING_ADJ_PARAM");
	iftWriteIntValues(fp,convnet->stride,convnet->nlayers,"STRIDE");
	iftWriteFloatValues(fp,convnet->alpha,convnet->nlayers,"ALPHA");
	iftWriteIntValues(fp,convnet->norm_adj_param,convnet->nlayers,"NORM_ADJ_PARAM");
	iftWriteIntValue(fp,convnet->rescale,"RESCALE");
	iftWriteIntValue(fp,convnet->with_weights,"WITH_WEIGHTS");

	if (convnet->with_weights) { /* write filters */
		for (int l=0; l < convnet->nlayers; l++) {
			iftWriteMMKernelFILE(convnet->k_bank[l],fp);
		}
	}

}


/*------------------------------------- PUBLIC ------------------------------------*/

float iftRelu(float f) {
	return iftMax(f, 0);
}

float iftReluDeriv(float f) {
	return (f>0)? 1.0 : 0.0;
}

float iftLeakyRelu(float f) {
	return (f>0)? f : 0.1*f;
}

float iftLeakyReluDeriv(float f) {
	return (f>0)? 1.0 : 0.1;
}

float iftSigmDeriv(float p) {
//	return p * (1.0 - p);
	float y = iftSigm(p);
	return y * (1 - y);
}

float iftSigm(float p) {
	double e = exp(-p);
	return 1.0/(1.0 + e);
}

float iftTanh(float f) {
	return 2 * iftSigm(2 * f) - 1;
}

float iftTanhDeriv(float f) {
	return 1.0 - iftFastPow(f, 2);
}

float iftNoAct(float f) {
	return f;
}

float iftNoActDeriv(float f) {
	return 1;
}

float iftForbiddenDeriv(float f) 
{
  iftError("Layer does not allow derivative computation.", "iftForbiddenDeriv");
  return 1;
}

void iftNoForward(iftNeuralNetwork* net, int l) {

}

void iftNoBackward(iftNeuralNetwork* net, int l) {}

void iftNoUpdate(iftNeuralNetwork* net, int l) {}

void iftForbiddenBackward(iftNeuralNetwork* net, int l) 
{
  iftError("Layer %d does not allow Backward operation.", "iftForbiddenDeriv", l);
}

void iftForbiddenUpdate(iftNeuralNetwork* net, int l) 
{
  iftError("Layer %d does not allow Update operation.", "iftForbiddenDeriv", l);
}

void iftPrintConvolutionalLayer(iftNeuralLayer* layer) {
	printf("- Convolutional { ");
	iftPrintTensorDim(layer->data);
	printf(" * %d x (%d, %d, %d) / %d => ", layer->noutput, layer->ninput, layer->kernelx, layer->kernely, layer->xstride);
	iftPrintTensorDim(layer->out);
	printf(" }");
}

void iftPrintFullyConnectedLayer(iftNeuralLayer* layer) {
	printf("- Fully Connected { ");
	iftPrintTensorDim(layer->data);
	printf(" * (%d, %d) => ", layer->ninput, layer->noutput);
	iftPrintTensorDim(layer->out);
	printf(" }");
}

void iftPrintCustomLayer(iftNeuralLayer* layer) {
	printf("- Custom { ");
	iftPrintTensorDim(layer->data);
	printf(" => ");
	iftPrintTensorDim(layer->out);
	printf(" }");
}

void iftPrintImageLayer(iftNeuralLayer* layer) {
    int *dim = layer->data->dimension;
	printf("- Image { %d x (%d, %d, %d) x %d band(s), %d classes => ", dim[0], dim[2], dim[3], dim[4], dim[1], layer->nclasses);
	iftPrintTensorDim(layer->out);
	printf(" }");
}

void iftPrintLossLayer(iftNeuralLayer* layer) {
	printf("- Square Loss { classes = %d }", layer->ninput);
}

void iftPrintDataSetLayer(iftNeuralLayer* layer) {
	printf("- DataSet { X = (%d samp, %d feat), %d classes => ", layer->dataset->ntrainsamples, layer->dataset->nfeats, layer->dataset->nclasses);
    iftPrintTensorDim(layer->out);
	printf(" }");
}

void iftPrintMaxPoolingLayer(iftNeuralLayer* layer) {
	printf("- Max Pooling { ");
	iftPrintTensorDim(layer->data);
	printf(" . (%d, %d) / %d => ", layer->kernelx, layer->kernely, layer->xstride);
	iftPrintTensorDim(layer->out);
	printf(" }");
}

void iftPrintNeuralLayer(iftNeuralLayer* layer) {
	switch (layer->type) {
		case IFT_FC_LAYER:
			iftPrintFullyConnectedLayer(layer);
			break;
		case IFT_CONV_LAYER:
			iftPrintConvolutionalLayer(layer);
			break;
		case IFT_LOSS_LAYER:
			iftPrintLossLayer(layer);
			break;
		case IFT_DATASET_LAYER:
			iftPrintDataSetLayer(layer);
			break;
		case IFT_MAX_POOLING_LAYER:
			iftPrintMaxPoolingLayer(layer);
			break;
        case IFT_IMAGE_LAYER:
			iftPrintImageLayer(layer);
			break;
		default:
			iftPrintCustomLayer(layer);
	}
}


void iftPrintNeuralNet(iftNeuralNetwork* net) {
    printf("General params:\n");
    printf("- Num. Layers: %d\n", net->nlayers);
    printf("- Learn rate: %.2f\n", net->learnRate);
    printf("- Minibatch size: %d\n", net->miniBatch);
    printf("Layers:\n");
	for (int l = 0; l < net->nlayers; ++l) {
		iftPrintNeuralLayer(net->layers[l]);
		printf("\n");
	}
}


void iftFullyConnectedBackward(iftNeuralNetwork* net, int l) {
	iftNeuralLayer* layer = net->layers[l];

	iftActivationFunction deriv = layer->actDeriv;

	for (int i = 0; i < layer->out->n; ++i) {
		layer->delta->val[i] *= deriv(layer->out->val[i]);
	}

	iftSetFloatArray(layer->biasUpdate->val, layer->biasUpdate->n, 0.0);

	for(int m = 0; m < net->miniBatch; ++m) {
		for(int i = 0; i < layer->noutput; ++i) {
			iftTensorElem(layer->biasUpdate, i) += iftTensorElem(layer->delta, m, i);
		}
	}

	iftMatrix* delta = iftCreateMatrixPointer(layer->delta->val, layer->noutput, net->miniBatch);
	iftMatrix* data = iftCreateMatrixPointer(layer->data->val, layer->ninput, net->miniBatch);
	iftMatrix* weightUpdate = iftCreateMatrixPointer(layer->weightUpdate->val, layer->ninput, layer->noutput);
	iftMatrix* weight = iftCreateMatrixPointer(layer->weight->val, layer->ninput, layer->noutput);
	iftMatrix* error = iftCreateMatrixPointer(layer->error->val, layer->ninput, net->miniBatch);

	iftMultMatricesInPlace(delta, data, true, false, &weightUpdate);//compute derivative to each w
	iftMultMatricesInPlace(delta, weight, false, false, &error);//compute delta for next step in backprop

//	if(net->iteration==3944)
//	{
//		printf("hue");
//
//		iftPrintMatrix(delta);
//		iftPrintMatrix(iftTransposeMatrix(delta));
//		iftPrintMatrix(data);
//	}

	iftDestroyMatrix(&delta);
	iftDestroyMatrix(&data);
	iftDestroyMatrix(&weightUpdate);
	iftDestroyMatrix(&weight);
	iftDestroyMatrix(&error);
}

iftImage* iftBatchToImage(iftTensor* data) {
	int batchSize = data->dimension[0];
	int nbands = data->dimension[1];

	int xsize = data->dimension[2];
	int ysize = data->dimension[3];


	iftImage *img = iftCreateImage(xsize * batchSize, ysize * nbands, 1);

	for (int m = 0; m < batchSize; ++m) {
		for (int b = 0; b < nbands; ++b) {
			for (int x = 0; x < xsize; ++x) {
				for (int y = 0; y < ysize; ++y) {
					iftImgVal2D(img, x + (xsize * m), y + ysize * b) = iftTensorElem(data, m, b, x, y);
				}
			}
		}
	}

	return img;
}

void iftConvolutionForward(iftNeuralNetwork* net, int l) 
{
  static const int MIN_OMP_SIZE = 5000000;
  iftNeuralLayer* layer = net->layers[l];

  iftMatrix* dataCol = iftCreateMatrixPointer(layer->buffer,
    layer->xout * layer->yout * net->miniBatch,
    layer->ninput * layer->kernelx * layer->kernely);

  // BatchImg2ColMatrix
  if (net->miniBatch > 1) {
    #pragma omp parallel for if(dataCol->n > MIN_OMP_SIZE)
    for (int m = 0; m < net->miniBatch; ++m) {
      iftImg2ColOptimizedBatch(layer->data->val + (m * layer->data->accDimension[0]), 
          layer->ninput, layer->ydata, layer->xdata, layer->kernelx, layer->xstride,
          layer->pad, net->miniBatch, m, dataCol->val);
    }
  } else {
    iftImg2ColOptimized(layer->data->val, layer->ninput, layer->ydata,
        layer->xdata, layer->kernelx, layer->xstride, layer->pad, dataCol->val);
  }

  // kernelMatrix
  iftMatrix* weight = iftCreateMatrixPointer(
      layer->weight->val,
      layer->ninput * layer->kernelx * layer->kernely,
      layer->weight->dimension[0]);

  // BatchResult
  float *outputPtr = net->miniBatch > 1 ? &(layer->buffer[dataCol->n]) : layer->out->val;
  iftMatrix *output = iftCreateMatrixPointer(
      outputPtr,
      layer->xout * layer->yout * net->miniBatch,
      weight->nrows);

  // kernelMatrix * BatchImg2ColMatrix = BatchResult
  iftMultMatricesInPlace(weight, dataCol, false, false, &output);

  // Apply bias + activation function
  iftActivationFunction activation = layer->act;
  #pragma omp parallel for if(output->n > MIN_OMP_SIZE), collapse(2)
  for (int i = 0; i < output->nrows; ++i)
    for (int j = 0; j < output->ncols; ++j)
      iftMatrixElem(output, j, i) = activation(
          iftMatrixElem(output, j, i) +
          iftTensorElem(layer->bias, i % layer->bias->dimension[0]));

  if (net->miniBatch > 1) {
    // Adjust to sample-wise memory layout
    // TODO Move to separate function
    #pragma omp parallel for if(output->n > MIN_OMP_SIZE)
    for (int m = 0; m < net->miniBatch; ++m) {
      int idx = m * (output->n / net->miniBatch);
      int nColsPerSample = output->ncols / net->miniBatch;
      for (int row = 0; row < output->nrows; ++row) {
        for (int sampleCol = 0; sampleCol < nColsPerSample; ++sampleCol) {
          int col = m * nColsPerSample + sampleCol;
          layer->out->val[idx++] = iftMatrixElem(output, col, row);
        }
      }
    }
  }

  iftDestroyMatrixPointer(&dataCol);
  iftDestroyMatrixPointer(&weight);
  iftDestroyMatrixPointer(&output);
}

void iftConvolutionBackward(iftNeuralNetwork* net, int l) {

	iftNeuralLayer* layer = net->layers[l];

	iftSetFloatArray(layer->weightUpdate->val, layer->weightUpdate->n, 0.0f);
	iftSetFloatArray(layer->biasUpdate->val, layer->biasUpdate->n, 0.0);

	iftActivationFunction deriv = layer->actDeriv;

//	iftMultScalarFloatArray(layer->delta->val, layer->delta->n, -1.0);

	for (int i = 0; i < layer->out->n; ++i) {
		layer->delta->val[i] *= deriv(layer->out->val[i]);
	}

	for(int m = 0; m < net->miniBatch; ++m) {
		#pragma omp parallel for
		for(int j = 0; j < layer->noutput; ++j) {
			iftTensorElem(layer->biasUpdate, j) += iftSumFloatArray(layer->delta->val +
																	m * layer->delta->accDimension[0]
																	+ j * layer->delta->accDimension[1], layer->xout * layer->yout);
		}
	}

	iftMatrix* dataCol = iftCreateMatrixPointer(layer->buffer, layer->xout*layer->yout, layer->ninput*layer->kernelx*layer->kernelx);
	iftMatrix* errorCol = iftCreateMatrixPointer(layer->buffer, layer->xout*layer->yout, layer->ninput*layer->kernelx*layer->kernelx);


	iftFloatArray* bufferror = iftCreateFloatArray(errorCol->n);

	for (int m = 0; m < net->miniBatch; ++m) {

		int n = iftImg2Col(layer->data->val + (m * layer->data->accDimension[0]),
						   layer->ninput, layer->ydata, layer->xdata, layer->kernelx, layer->xstride, 0, dataCol->val);

		if(n!=dataCol->n) {
			iftError("Something wrong is not right here", "iftConvolutionBackward");
		}

		iftSetFloatArray(bufferror->val, bufferror->n, 0.0);

		//#pragma omp parallel for
		for (int j = 0; j < layer->noutput; ++j) {
			iftMatrix* delta = iftCreateMatrixPointer(layer->delta->val + m*layer->delta->accDimension[0] + j*layer->delta->accDimension[1], layer->xout * layer->yout, 1);
			iftMatrix* weightUpdate = iftCreateMatrixPointer(layer->weightUpdate->val + j * layer->weightUpdate->accDimension[0], layer->ninput * layer->kernelx * layer->kernely, 1);
			iftMatrix* weight = iftCreateMatrixPointer(layer->weight->val + j * layer->weight->accDimension[0], layer->ninput * layer->kernelx * layer->kernely, 1);

//			printf("m = %d, j = %d\n", m, j);
//			iftPrintMatrix(weight);
//			iftPrintMatrix(delta);

//			printf("============\n");

			iftMultMatricesInPlace(delta, dataCol, false, true, &weightUpdate);
			iftMultMatricesInPlace(weight, delta, true, false, &errorCol);

//			iftPrintMatrix(weightUpdate);

//			iftCol2Img(errorCol->val, layer->ninput, layer->ydata, layer->xdata, layer->kernelx, layer->xstride,
//					   0, layer->error->val + (m * layer->error->accDimension[0]));

			iftLinearCombination(1.0, errorCol->val, 1.0, bufferror->val, bufferror->val, bufferror->n);

			iftDestroyMatrix(&delta);
			iftDestroyMatrix(&weightUpdate);
		}

		iftCol2Img(bufferror->val, layer->ninput, layer->ydata, layer->xdata, layer->kernelx, layer->xstride,
					   0, layer->error->val + (m * layer->error->accDimension[0]));

//		iftCopyFloatArray(layer->error->val + (m * layer->error->accDimension[0]), bufferror->val, bufferror->n);
	}

	iftDestroyFloatArray(&bufferror);

//	printf("WeightUpdate: ");
//	iftPrintFloatArray(layer->weightUpdate->val, layer->weightUpdate->n);
//	printf("Error: ");
//	iftPrintFloatArray(layer->error->val, layer->error->n);

//	iftMultScalarFloatArray(layer->weightUpdate->val, layer->weightUpdate->n, -1.0);
//	iftMultScalarFloatArray(layer->biasUpdate->val, layer->biasUpdate->n, -1.0);
//
//	iftMultScalarFloatArray(layer->delta->val, layer->delta->n, -1.0);
//	iftMultScalarFloatArray(layer->error->val, layer->error->n, -1.0);


	iftDestroyMatrix(&dataCol);
	iftDestroyMatrix(&errorCol);

}

iftNeuralLayer *iftConvolutionLayer2D(iftNeuralNetwork *net, int xsize, int ysize,
                                      int nbands, int kernelXSize, int kernelYSize,
                                      int xstride, int ystride, int nkernels, int pad,
                                      iftActivationFunction act, iftActivationFunction deriv) {

	iftNeuralLayer* layer = iftCreateNeuralLayer();

	int minibatch = net->miniBatch;

	layer->kernelx = kernelXSize;
	layer->kernely = kernelYSize;

	layer->xstride = xstride;
	layer->ystride = ystride;

	layer->ninput = nbands;
	layer->noutput = nkernels;

	layer->xdata = xsize;
	layer->ydata = ysize;

        layer->pad = pad;
	layer->yout = (ysize + 2*pad - kernelYSize) / ystride + 1;
	layer->xout = (ysize + 2*pad - kernelXSize) / xstride + 1;

	layer->act = act;
	layer->actDeriv = deriv;

	layer->data = iftCreateTensor(4, minibatch, nbands, xsize, ysize);
	layer->error = iftCreateTensor(4, minibatch, nbands, xsize, ysize);

	layer->weight = iftCreateTensor(4, nkernels, nbands, kernelXSize, kernelYSize);
	layer->weightUpdate = iftCreateTensor(4, nkernels, nbands, kernelXSize, kernelYSize);

	layer->bias = iftCreateTensor(1, nkernels);
	layer->biasUpdate = iftCreateTensor(1, nkernels);

	layer->out  = iftCreateTensor(4, minibatch, nkernels, layer->xout, layer->yout);
	layer->delta = iftCreateTensor(4, minibatch, nkernels, layer->xout, layer->yout);

//	layer->deltaW = iftCreateTensor(4, nkernels, nbands, kernelXSize, kernelYSize);

        if (net->miniBatch > 1)
          layer->buffer = iftAllocFloatArray((xsize * ysize * kernelXSize * kernelYSize * nbands * net->miniBatch) + layer->out->n);
        else
          layer->buffer = iftAllocFloatArray(xsize * ysize * kernelXSize * kernelYSize * nbands);

	float std = sqrtf(2.0f/(layer->kernelx*layer->kernely*layer->ninput + layer->kernelx*layer->kernely*layer->noutput));//Xavier initialization
	for (int i = 0; i < layer->weight->n; ++i) {
		layer->weight->val[i] = iftRandomNormalFloat(0, std*std);
//		layer->weight->val[i] = iftRandomUniform(-0.05, 0.05);
	}

	layer->forward = iftConvolutionForward;
	layer->backward = iftConvolutionBackward;
	layer->update = iftConvolutionUpdate;

	layer->type = IFT_CONV_LAYER;

	return layer;
}

iftNeuralLayer *iftConvolutionLayer2DForwardOnly(iftNeuralNetwork *net,
    int xsize, int ysize, int nbands, int kernelXSize, int kernelYSize,
    int xstride, int ystride, int nkernels, int pad, iftActivationFunction act)
{
  iftNeuralLayer* layer = iftCreateNeuralLayer();

  int minibatch = net->miniBatch;

  layer->kernelx = kernelXSize;
  layer->kernely = kernelYSize;

  layer->xstride = xstride;
  layer->ystride = ystride;

  layer->ninput = nbands;
  layer->noutput = nkernels;

  layer->xdata = xsize;
  layer->ydata = ysize;

  layer->pad = pad;
  layer->yout = (ysize + 2*pad - kernelYSize) / ystride + 1;
  layer->xout = (ysize + 2*pad - kernelXSize) / xstride + 1;

  layer->act = act;
  layer->actDeriv = iftForbiddenDeriv;

  // Actual pointer is defined later in iftNetworkAllocateSharedBuffers
  layer->data = iftCreateTensorPointer(NULL, 4, minibatch, nbands, xsize, ysize);
  layer->error = NULL;

  layer->weight = iftCreateTensor(4, nkernels, nbands, kernelXSize, kernelYSize);
  layer->weightUpdate = NULL;

  layer->bias = iftCreateTensor(1, nkernels);
  layer->biasUpdate = NULL;

  // Actual pointer is defined later in iftNetworkAllocateSharedBuffers
  layer->out  = iftCreateTensorPointer(NULL, 4, minibatch, nkernels, layer->xout, layer->yout);
  layer->delta = NULL;

  // Actual pointer is defined later in iftNetworkAllocateSharedBuffers
  layer->buffer = NULL; 
  if (net->miniBatch > 1)
    layer->bufferSize = ((xsize * ysize * kernelXSize * kernelYSize * nbands * net->miniBatch) + layer->out->n) * sizeof(float);
  else
    layer->bufferSize = (xsize * ysize * kernelXSize * kernelYSize * nbands) * sizeof(float);

  layer->forward = iftConvolutionForward;
  layer->backward = iftForbiddenBackward;
  layer->update = iftForbiddenUpdate;

  layer->type = IFT_CONV_LAYER;
  layer->isSharedBuffer = true;

  return layer;
}

void iftAppendLayer(iftNeuralNetwork* net, iftNeuralLayer* layer) {
	net->layers[net->curLayer] = layer;
	net->curLayer++;
}

iftNeuralLayer* iftFullyConnectedLayer(iftNeuralNetwork* net, int ninput, int noutput, iftActivationFunction act, iftActivationFunction deriv) {
	iftNeuralLayer* layer = iftCreateNeuralLayer();

	int minibatch = net->miniBatch;

	layer->ninput = ninput;
	layer->noutput = noutput;

	layer->data = iftCreateTensor(2, minibatch, ninput);
	layer->error= iftCreateTensor(2, minibatch, ninput);

	layer->out = iftCreateTensor(2, minibatch, noutput);
	layer->delta = iftCreateTensor(2, minibatch, noutput);

	layer->weight = iftCreateTensor(2, noutput, ninput);
	layer->weightUpdate = iftCreateTensor(2, noutput, ninput);

	layer->bias = iftCreateTensor(1, noutput);
	layer->biasUpdate = iftCreateTensor(1, noutput);

	float std = sqrtf(2.0f/(ninput + noutput));//Xavier initialization
	for (int i = 0; i < layer->weight->n; ++i) {
		layer->weight->val[i] = iftRandomNormalFloat(0, std*std);
//		layer->weight->val[i] = iftRandomUniform(-0.05, 0.05);
	}

	for (int i = 0; i < layer->bias->n; ++i) {
		layer->bias->val[i] = 0.0;
	}

	layer->act = act;
	layer->actDeriv = deriv;

	layer->forward = iftFullyConnectedForward;
	layer->backward = iftFullyConnectedBackward;
	layer->update = iftFullyConnectedUpdate;

	layer->type = IFT_FC_LAYER;

	return layer;
}

iftNeuralLayer* iftFullyConnectedLayerForwardOnly(iftNeuralNetwork* net, int ninput, int noutput, iftActivationFunction act) 
{
  iftNeuralLayer* layer = iftCreateNeuralLayer();

  int minibatch = net->miniBatch;

  layer->ninput = ninput;
  layer->noutput = noutput;

  // Actual pointer is defined later in iftNetworkAllocateSharedBuffers
  layer->data = iftCreateTensorPointer(NULL, 2, minibatch, ninput);
  layer->error = NULL;

  // Actual pointer is defined later in iftNetworkAllocateSharedBuffers
  layer->out = iftCreateTensorPointer(NULL, 2, minibatch, noutput);
  layer->delta = NULL;

  // Does not use net->buffer
  layer->bufferSize = 0; 

  layer->weight = iftCreateTensor(2, noutput, ninput);
  layer->weightUpdate = iftCreateTensor(2, noutput, ninput);

  layer->bias = iftCreateTensor(1, noutput);
  layer->biasUpdate = iftCreateTensor(1, noutput);

  layer->act = act;
  layer->actDeriv = iftForbiddenDeriv;

  layer->forward = iftFullyConnectedForward;
  layer->backward = iftForbiddenBackward;
  layer->update = iftForbiddenUpdate;

  layer->type = IFT_FC_LAYER;
  layer->isSharedBuffer = true;

  return layer;
}

iftNeuralLayer* iftEmptyImageLayer(iftNeuralNetwork* net, int xSize, int ySize, int zSize, int nBands)
{
  iftNeuralLayer* layer = iftCreateNeuralLayer();

  // Does not use fileset interface as images are loaded externally
  layer->fileset = NULL;
  layer->idx = NULL;
  layer->idxCounter = 0;

  layer->nclasses = 0;
  layer->ninput = nBands;
  layer->noutput = nBands;
  layer->xdata = xSize;
  layer->ydata = ySize;
  layer->zdata = zSize;

  // TODO check if data and delta are actually needed anywhere
  layer->data = iftCreateTensor(5, net->miniBatch, nBands, layer->xdata, layer->ydata, layer->zdata);
  layer->delta = iftCreateTensor(5, net->miniBatch, nBands, layer->xdata, layer->ydata, layer->zdata);
  layer->out = iftCreateTensor(5, net->miniBatch, nBands, layer->xdata, layer->ydata, layer->zdata);

  layer->weight = NULL;

  // TODO Create a iftForbiddenForward to avoid operation without loading any image
  layer->forward = iftNoForward;
  layer->backward = iftNoBackward;
  layer->update = iftNoUpdate;

  layer->type = IFT_IMAGE_LAYER;

  return layer;
}

iftNeuralLayer* iftEmptyImageLayerForwardOnly(iftNeuralNetwork* net, int xSize, int ySize, int zSize, int nBands)
{
  iftNeuralLayer* layer = iftCreateNeuralLayer();

  layer->fileset = NULL;
  layer->idx = NULL;
  layer->idxCounter = 0;

  layer->nclasses = 0;
  layer->ninput = nBands;
  layer->noutput = nBands;
  layer->xdata = xSize;
  layer->ydata = ySize;
  layer->zdata = zSize;

  // Not actually used, only included to simplify shared buffer interface
  layer->data = iftCreateTensorPointer(NULL, 5, net->miniBatch, nBands, layer->xdata, layer->ydata, layer->zdata);
  layer->error = NULL;

  // Actual pointer is defined later in iftNetworkAllocateSharedBuffers
  layer->out = iftCreateTensorPointer(NULL, 5, net->miniBatch, nBands, layer->xdata, layer->ydata, layer->zdata);
  layer->delta = NULL;

  layer->weight = NULL;
  layer->weightUpdate = NULL;

  // Does not use net->buffer
  layer->bufferSize = 0;

  // TODO Create a iftForbiddenForward to avoid operation without loading any image
  layer->forward = iftNoForward;
  layer->backward = iftNoBackward;
  layer->update = iftNoUpdate;

  layer->type = IFT_IMAGE_LAYER;
  layer->isSharedBuffer = true;

  return layer;
  net->layers[0]->idxCounter = 0;
}

int iftFeedEmptyImageLayerImgPath(iftNeuralNetwork* net, const char *img_path)
{
  iftNeuralLayer* layer = net->layers[0];
  assert(layer->type == IFT_IMAGE_LAYER);
  assert(layer->idxCounter < net->miniBatch);

  int i = layer->idxCounter;
  iftImage* img = iftReadImageByExt(img_path);
  assert(img->xsize == layer->xdata);
  assert(img->ysize == layer->ydata);
  assert(img->zsize == layer->zdata);
  if(iftIsColorImage(img)) {
    assert(layer->ninput == 3);
    iftImage* channel0 = iftImageRed(img);
    iftImage* channel1 = iftImageGreen(img);
    iftImage* channel2 = iftImageBlue(img);

    for (int x = 0; x < layer->xdata; ++x) {
      for (int y = 0; y < layer->ydata; ++y) {
        for (int z = 0; z < layer->zdata; ++z) {
          iftTensorElem(layer->out, i, 0, x, y, z) = iftImgVal(channel0, x, y, z);
          iftTensorElem(layer->out, i, 1, x, y, z) = iftImgVal(channel1, x, y, z);
          iftTensorElem(layer->out, i, 2, x, y, z) = iftImgVal(channel2, x, y, z);
        }
      }
    }

    iftDestroyImage(&channel0);
    iftDestroyImage(&channel1);
    iftDestroyImage(&channel2);
  } else{
    assert(layer->ninput == 1);
    for (int x = 0; x < layer->xdata; ++x) {
      for (int y = 0; y < layer->ydata; ++y) {
        for (int z = 0; z < layer->zdata; ++z) {
          iftTensorElem(layer->out, i, 0, x, y, z) = iftImgVal(img, x, y, z);
        }
      }
    }
  }

  if (net->miniBatch > 1) {
    layer->idxCounter++;
    return net->miniBatch - layer->idxCounter;
  } else {
    // Simplify operation without batches
    return 1;
  }
}

int iftFeedEmptyImageLayerColorImg(iftNeuralNetwork *net, int *R, int *G, int *B)
{
  iftNeuralLayer* layer = net->layers[0];
  assert(layer->type == IFT_IMAGE_LAYER);
  assert(layer->ninput == 3);
  assert(layer->idxCounter < net->miniBatch);

  int i = layer->idxCounter;
  int idx = 0;
  for (int z = 0; z < layer->zdata; ++z) {
    for (int y = 0; y < layer->ydata; ++y) {
      for (int x = 0; x < layer->xdata; ++x) {
        iftTensorElem(layer->out, i, 0, x, y, z) = R[idx];
        iftTensorElem(layer->out, i, 1, x, y, z) = G[idx];
        iftTensorElem(layer->out, i, 2, x, y, z) = B[idx];
        idx++;
      }
    }
  }

  if (net->miniBatch > 1) {
    layer->idxCounter++;
    return net->miniBatch - layer->idxCounter;
  } else {
    // Simplify operation without batches
    return 1;
  }
}

void iftClearEmptyImageLayerBatch(iftNeuralNetwork *net)
{
  iftNeuralLayer* layer = net->layers[0];
  assert(layer->type == IFT_IMAGE_LAYER);
  layer->idxCounter = 0;
}

iftNeuralLayer* iftImageLayer(iftNeuralNetwork* net, const iftFileSet* fileset)
{
	iftNeuralLayer* layer = iftCreateNeuralLayer();

	layer->imgs = iftCreateDict();

	int nclasses = 0;

	int* labels = iftFileSetLabels(fileset);
	nclasses = iftCountUniqueIntElems(labels, fileset->n);
	iftFree(labels);

	int ntrainsamples = 0;

	for (int i = 0; i < fileset->n; ++i) {
		if(fileset->files[i]->status & IFT_TRAIN)
			ntrainsamples++;
	}

	net->label = iftCreateTensor(2, net->miniBatch, nclasses);

	layer->fileset = iftCopyFileSet(fileset);
	layer->idx = iftCreateIntArray(ntrainsamples);
	layer->status = IFT_TRAIN;

	int c = 0;
	for (int i = 0; i < fileset->n; ++i) {
	if(fileset->files[i]->status == IFT_TRAIN)
			layer->idx->val[c++] = i;
	}

	layer->idxCounter = 0;

	iftShuffleIntArray(layer->idx->val, layer->idx->n);

	///get image size
	iftImage* img = iftReadImageByExt(fileset->files[0]->path);
	int nbands = iftIsColorImage(img)?3:1;

	layer->nclasses = nclasses;
	layer->ninput = nbands;
	layer->noutput = nbands;
	layer->xdata = img->xsize;
	layer->ydata = img->ysize;
	layer->zdata = img->zsize;
	iftDestroyImage(&img);
	///

	layer->data = iftCreateTensor(5, net->miniBatch, nbands, layer->xdata, layer->ydata, layer->zdata);
//	layer->error = iftCreateTensor(5, net->miniBatch, nbands, layer->xdata, layer->ydata, layer->zsize);
	layer->delta = iftCreateTensor(5, net->miniBatch, nbands, layer->xdata, layer->ydata, layer->zdata);
    layer->out = iftCreateTensor(5, net->miniBatch, nbands, layer->xdata, layer->ydata, layer->zdata);

	layer->weight = NULL;

	layer->forward = iftImageForward;
	layer->backward = iftNoBackward;
	layer->update = iftNoUpdate;

        layer->type = IFT_IMAGE_LAYER;

	return layer;
}

void iftRawDataForward(iftNeuralNetwork* net, int l) {

	iftNeuralLayer* layer = net->layers[l];
	iftNeuralLayer* nextLayer = net->layers[l];

	for (int m = 0; m < net->miniBatch; ++m) {

		for (int c = 0; c < layer->X->ncols; ++c) {
			int idx = layer->idx->val[layer->idxCounter];
			float x = iftMatrixElem(layer->X, c, idx);
			iftTensorElem(layer->data, m, c) = x;
			iftTensorElem(nextLayer->data, m, c) = x;
		}

		for (int c = 0; c < layer->Y->ncols; ++c) {
			int idx = layer->idx->val[layer->idxCounter];
			float y = iftMatrixElem(layer->Y, c, idx);
			iftTensorElem(net->label, m, c) = y;
		}

		layer->idxCounter++;
	}

	if(layer->idxCounter >= layer->X->nrows) {
		layer->idxCounter = 0;
		net->epoch++;
		iftShuffleIntArray(layer->idx->val, layer->idx->n);
	}

}

iftNeuralLayer* iftRawDataLayer(iftNeuralNetwork* net, const iftMatrix* X, const iftMatrix* Y) {

	if(X->nrows != Y->nrows) {
        iftError("Different number of dimensions X = %d, Y = %d.", "iftRawDataLayer", X->nrows, Y->nrows);
	}

	int nsamples = X->nrows;

	iftNeuralLayer* layer = iftCreateNeuralLayer();
	layer->idxCounter = 0;

	layer->X = iftCopyMatrix(X);
	layer->Y = iftCopyMatrix(Y);

	layer->data = iftCreateTensor(2, net->miniBatch, X->ncols);
	net->label  = iftCreateTensor(2, net->miniBatch, Y->ncols);

	layer->idx = iftIntRange(0, nsamples-1, 1);
	iftShuffleIntArray(layer->idx->val, layer->idx->n);

	layer->update   = iftNoUpdate;
	layer->forward  = iftRawDataForward;
	layer->backward = iftNoBackward;

	return layer;
}

iftNeuralLayer* iftDataSetLayer(iftNeuralNetwork* net, const iftDataSet* Z) {
	iftNeuralLayer* layer = iftCreateNeuralLayer();

	net->label = iftCreateTensor(2, net->miniBatch, Z->nclasses);

	layer->dataset = iftCopyDataSet(Z, true);
	layer->idx = iftCreateIntArray(layer->dataset->nsamples);
	layer->status = IFT_TRAIN;

	int c = 0;
	for (int i = 0; i < layer->dataset->nsamples; ++i) {
		if(iftHasSampleStatus(layer->dataset->sample[i], IFT_TRAIN))
			layer->idx->val[c++] = i;
	}

	layer->idxCounter = 0;

	iftShuffleIntArray(layer->idx->val, layer->idx->n);

	layer->data = iftCreateTensor(2, net->miniBatch, Z->nfeats);
	layer->out = iftCreateTensor(2, net->miniBatch, Z->nfeats);

	layer->mean = iftCreateTensor(1, Z->nfeats);
	layer->std = iftCreateTensor(1, Z->nfeats);

	for(int j=0; j< layer->dataset->nfeats; ++j) {
		iftTensorElem(layer->mean, j) = iftTensorElem(layer->std, j) = 0;
		for (int s = 0; s < layer->dataset->nsamples; ++s) {
			iftTensorElem(layer->mean, j) += layer->dataset->sample[s].feat[j] / (layer->dataset->nsamples);
//			iftTensorElem(layer->mean, j) = 0;
		}

		for (int s = 0; s < layer->dataset->nsamples; ++s) {
			iftTensorElem(layer->std, j) += iftFastPow(layer->dataset->sample[s].feat[j] - iftTensorElem(layer->mean, j), 2) / (layer->dataset->nsamples - 1);
		}

//		if( iftAlmostZero(iftTensorElem(layer->std, j))) {
//			iftTensorElem(layer->std, j) = 1.0;
//		}

		iftTensorElem(layer->std, j) = sqrtf(iftTensorElem(layer->std, j));
//		iftTensorElem(layer->mean, j) = 1.0;
	}


	layer->ninput = Z->nfeats;
	layer->noutput = Z->nfeats;

	layer->forward = iftDataSetForward;
	layer->backward = iftNoBackward;
	layer->update = iftNoUpdate;

	layer->type = IFT_DATASET_LAYER;

	return layer;
}

void iftNeuralNetGradientDescent(iftNeuralNetwork* net, int epochs) {

	net->epoch = 0;
	net->iteration = 0;
	int epoch = 0;
	timer *start, *end;
	double error = 0.0f;

	int i = 0;
	start = iftTic();
	while(net->epoch < epochs) {
        if(net->verbose)
            printf("Iteration: %d\n", net->iteration+1);

		iftNetForward(net);
		iftNetBackward(net);

		error += net->error/net->miniBatch;
		i++;
		iftNetUpdateWeight(net);

		net->iteration++;

		if(net->epoch > epoch){
			epoch = net->epoch;

			end = iftToc();
			
			if(net->verbose) {
				printf("Epoch %d/%d (error: %lf) ", net->epoch, epochs, error/(i));
				iftPrintFormattedTime(iftCompTime(start,end));
			}
			error = 0.0;
			i=0;

			start = iftTic();
		}
	}

}

iftNeuralLayer* iftMaxPoolingLayer2D(iftNeuralNetwork* net, int nbands, int xsize, int ysize, int poolsize, int stride) 
{
	iftNeuralLayer* layer = iftCreateNeuralLayer();

	int minibatch = net->miniBatch;

	layer->ninput = nbands;
	layer->noutput = nbands;
	layer->xdata = xsize;
	layer->ydata = ysize;
	layer->kernelx = poolsize;
	layer->kernely = poolsize;
	layer->xstride = stride;
	layer->ystride = stride;

	layer->data = iftCreateTensor(4, minibatch, nbands, xsize, ysize);
	layer->dataMask = iftCreateTensor(4, minibatch, nbands, xsize, ysize);
	layer->error = iftCreateTensor(4, minibatch, nbands, xsize, ysize);

	int pad = 0;
        layer->pad = pad;
	layer->yout = (ysize + 2*pad - poolsize) / stride + 1;
	layer->xout = (xsize + 2*pad - poolsize) / stride + 1;

	layer->out  = iftCreateTensor(4, minibatch, nbands, layer->xout, layer->yout);
	layer->delta = iftCreateTensor(4, minibatch, nbands, layer->xout, layer->yout);

	layer->ibuffer = iftAllocIntArray(xsize*ysize*poolsize*poolsize*nbands);

	layer->forward = iftMaxPoolingForward;
	layer->backward = iftMaxPoolingBackward;
	layer->update = iftNoUpdate;

	layer->type = IFT_MAX_POOLING_LAYER;

	return layer;
}

iftNeuralLayer* iftMaxPoolingLayer2DForwardOnly(iftNeuralNetwork* net, int nbands, int xsize, int ysize, int poolsize, int stride) 
{
  iftNeuralLayer* layer = iftCreateNeuralLayer();

  int minibatch = net->miniBatch;

  layer->ninput = nbands;
  layer->noutput = nbands;
  layer->xdata = xsize;
  layer->ydata = ysize;
  layer->kernelx = poolsize;
  layer->kernely = poolsize;
  layer->xstride = stride;
  layer->ystride = stride;

  // Actual pointer is defined later in iftNetworkAllocateSharedBuffers
  layer->data = iftCreateTensorPointer(NULL, 4, minibatch, nbands, xsize, ysize);
  layer->dataMask = NULL;
  layer->error = NULL;

  int pad = 0;
  layer->pad = pad;
  layer->yout = (ysize + 2*pad - poolsize) / stride + 1;
  layer->xout = (xsize + 2*pad - poolsize) / stride + 1;

  // Actual pointer is defined later in iftNetworkAllocateSharedBuffers
  layer->out = iftCreateTensorPointer(NULL, 4, minibatch, nbands, layer->xout, layer->yout);
  layer->delta = NULL;

  // Actual pointer is defined later in iftNetworkAllocateSharedBuffers
  layer->ibuffer = NULL;
  layer->bufferSize = xsize * ysize * poolsize * poolsize * nbands * sizeof(int);

  layer->forward = iftMaxPoolingForwardNoMask;
  layer->backward = iftForbiddenBackward;
  layer->update = iftNoUpdate;

  layer->type = IFT_MAX_POOLING_LAYER;
  layer->isSharedBuffer = true;

  return layer;
}

void iftMaxPoolingBackward(iftNeuralNetwork* net, int l) {
	iftNeuralLayer* layer = net->layers[l];

//	iftMatrix* dataCol = iftCreateMatrixPointer(layer->buffer, layer->xout*layer->yout, layer->ninput*layer->kernelx*layer->kernely);

	int outsize = layer->xout*layer->yout;
	int ksize   = layer->kernelx*layer->kernely;

	for(int m = 0; m < net->miniBatch; ++m) {

		iftImg2ColIdx(1, layer->ydata, layer->xdata, layer->kernelx, layer->xstride, 0, layer->ibuffer);

		#pragma omp parallel for
		for (int j = 0; j < layer->ninput; ++j) {

			float *delta = layer->delta->val + m * layer->delta->accDimension[0] + j * layer->delta->accDimension[1];
			float *mask = layer->dataMask->val + m * layer->dataMask->accDimension[0] + j * layer->dataMask->accDimension[1];
			float *error = layer->error->val + m * layer->error->accDimension[0] + j * layer->error->accDimension[1];
			int *idx = layer->ibuffer;

			for (int j = 0; j < outsize; ++j) {
				int imax = 0;
				for (int i = 0; i < ksize; ++i) {
					int idxcur = idx[j + i * outsize];
					int idxmax = idx[j + imax * outsize];
					if (idxcur != IFT_NIL && idxmax != IFT_NIL && mask[idxcur] > mask[idxmax])
						imax = i;
				}
				error[idx[j + imax * outsize]] = delta[j];
			}

		}

	}

}

void iftMaxPoolingForward(iftNeuralNetwork* net, int l) 
{
	iftNeuralLayer* layer = net->layers[l];

	int outsize = layer->xout*layer->yout;
	int ksize   = layer->kernelx*layer->kernely;

	for(int m = 0; m < net->miniBatch; ++m) {
		iftImg2ColIdx(1, layer->ydata, layer->xdata, layer->kernelx, layer->xstride, 0, layer->ibuffer);

		#pragma omp parallel for
		for (int j = 0; j < layer->ninput; ++j) {

			float *data = layer->data->val + m * layer->data->accDimension[0] + j * layer->data->accDimension[1];
			float *mask = layer->dataMask->val + m * layer->dataMask->accDimension[0] + j * layer->dataMask->accDimension[1];
			float *out = layer->out->val + m * layer->out->accDimension[0] + j * layer->out->accDimension[1];
			int *idx = layer->ibuffer;

			iftSetFloatArray(mask, outsize, 0.0);

			for (int j = 0; j < outsize; ++j) {
				int imax = 0;
				for (int i = 0; i < ksize; ++i) {
//				printf("%d %d %d %d\n", j, i, ksize, j + (i * outsize));
					int idxcur = idx[j + i * outsize];
					int idxmax = idx[j + imax * outsize];
					if (idxcur != IFT_NIL && idxmax != IFT_NIL && data[idxcur] > data[idxmax])
						imax = i;
				}
				out[j] = data[idx[j + imax * outsize]];
				mask[idx[j + imax * outsize]] = 1.0;
			}
		}
	}
}

void iftMaxPoolingForwardNoMask(iftNeuralNetwork* net, int l) 
{
  iftNeuralLayer* layer = net->layers[l];

  int outsize = layer->xout*layer->yout;
  int ksize   = layer->kernelx*layer->kernely;

  for(int m = 0; m < net->miniBatch; ++m) {
    iftImg2ColIdx(1, layer->ydata, layer->xdata, layer->kernelx, layer->xstride, 0, layer->ibuffer);

    #pragma omp parallel for
    for (int j = 0; j < layer->ninput; ++j) {
      float *data = layer->data->val + m * layer->data->accDimension[0] + j * layer->data->accDimension[1];
      float *out = layer->out->val + m * layer->out->accDimension[0] + j * layer->out->accDimension[1];
      int *idx = layer->ibuffer;

      for (int j = 0; j < outsize; ++j) {
        int imax = 0;
        for (int i = 0; i < ksize; ++i) {
          int idxcur = idx[j + i * outsize];
          int idxmax = idx[j + imax * outsize];
          if (idxcur != IFT_NIL && idxmax != IFT_NIL && data[idxcur] > data[idxmax])
            imax = i;
        }
        out[j] = data[idx[j + imax * outsize]];
      }
    }
  }
}

void iftFullyConnectedUpdate(iftNeuralNetwork* net, int l) {

	iftNeuralLayer* layer = net->layers[l];

	const double alpha = net->learnRate;
	//const double gamma = net->weightDecay;
	const double minibatch = net->miniBatch;
	const double decay = net->weightDecay;

	iftLinearCombination(1.0, layer->bias->val, -(alpha / minibatch), layer->biasUpdate->val, layer->bias->val, layer->bias->n);

	if(decay>0)
		iftLinearCombination(1.0, layer->weightUpdate->val, -decay * minibatch, layer->weight->val, layer->weightUpdate->val, layer->weightUpdate->n);

	iftLinearCombination(1.0, layer->weight->val, -(alpha/minibatch), layer->weightUpdate->val, layer->weight->val, layer->weight->n);

}

void iftNetUpdateWeight(iftNeuralNetwork *net) {

	for (int l = net->nlayers-1; l >= 0; l--) {
		net->layers[l]->update(net, l);
	}

}

void iftNetLayerBackward(iftNeuralNetwork* net, int l) {
	iftNeuralLayer* layer = net->layers[l];
	layer->backward(net, l);
	if(l>0) {
		iftNeuralLayer* prev = net->layers[l-1];
		if(prev->delta!=NULL)
			iftCopyFloatArray(prev->delta->val, layer->error->val, layer->error->n);
	}
}

void iftNetBackward(iftNeuralNetwork* net) {

	for (int l = net->nlayers - 1; l >=0 ; --l) {
//		printf("Backward #%d\n", l);
		iftNetLayerBackward(net, l);
	}

}

void iftDataSetLayerStatus(iftNeuralLayer* layer, iftSampleStatus status) {

	if(layer->idx != NULL)
		iftDestroyIntArray(&(layer->idx));

	layer->idx = iftCreateIntArray(iftCountSamples(layer->dataset, status));

	int c = 0;
	for (int i = 0; i < layer->dataset->nsamples; ++i) {
		if(iftHasSampleStatus(layer->dataset->sample[i], status))
			layer->idx->val[c++] = i;
	}

	layer->status = status;

	if(status == IFT_TRAIN)
		iftShuffleIntArray(layer->idx->val, layer->idx->n);

	layer->idxCounter = 0;

}

void iftImageLayerStatus(iftNeuralLayer* layer, iftSampleStatus status) {
	layer->status = status;

	if(layer->idx != NULL)
		iftDestroyIntArray(&(layer->idx));

	int nfiles = 0;

	for (int i = 0; i < layer->fileset->n; ++i) {
		if(layer->fileset->files[i]->status & status)
			nfiles++;
	}

	layer->idx = iftCreateIntArray(nfiles);

	int c = 0;
	for (int i = 0; i < layer->fileset->n; ++i) {
		if(layer->fileset->files[i]->status & status)
			layer->idx->val[c++] = i;
	}

	layer->status = status;

	if(status == IFT_TRAIN)
		iftShuffleIntArray(layer->idx->val, layer->idx->n);

	layer->idxCounter = 0;

}

iftTensor* iftNeuralNetClassifyImageSet(iftNeuralNetwork* net, iftFileSet* fileSet, uchar status) {
	iftNeuralLayer* lossLayer = net->layers[net->nlayers - 1];
	iftNeuralLayer* dataLayer = iftImageLayer(net, fileSet);

	iftNeuralLayer* old = net->layers[0];

	net->layers[0] = dataLayer;
	iftImageLayerStatus(dataLayer, status);

	int nsamples = dataLayer->idx->n;
	int nclasses = lossLayer->out->n / net->miniBatch;

	iftTensor* pred = iftCreateTensor(2, nsamples, nclasses);

	int offset = 0;
	int i = 0;
	int batchNum = 1;

	dataLayer->idxCounter = 0;

	while(offset < (dataLayer->idx->n)) {
		fprintf(stderr, "Starting batch %d (%d/%ld):\n", batchNum++, offset, dataLayer->idx->n);
		timer *t1 = iftTic();
		iftNetForward(net);

		for(i = 0; i < iftMin(net->miniBatch, dataLayer->idx->n - offset); ++i) {
			int idx = dataLayer->idx->val[offset + i];

			for (int c = 0; c < lossLayer->out->dimension[1]; c++) {
				pred->val[idx * pred->accDimension[0] + c] = lossLayer->out->val[i * lossLayer->out->accDimension[0] + c];
				//iftTensorElem(pred, idx, c) = iftTensorElem(lossLayer->out, i, c);
			}
		}
		offset += net->miniBatch;
		timer *t2 = iftToc();
		fprintf(stderr, "Finished in %fms\n", iftCompTime(t1,t2));
	}

	net->layers[0] = old;
	iftDestroyNeuralLayer(&dataLayer);

	return pred;
}

void iftNeuralNetClassifyDataSet(iftNeuralNetwork* net, iftDataSet* Z, uchar status) {

	iftNeuralLayer* lossLayer = net->layers[net->nlayers - 1];
//    iftDataSet* Z = dataLayer->dataset;
//    iftSampleStatus status = dataLayer->status;
//    if(status != IFT_ALL)
//        iftDataSetLayerStatus(dataLayer, IFT_ALL);
//

	iftNeuralLayer* dataLayer = iftDataSetLayer(net, Z);

	iftNeuralLayer* old = net->layers[0];
	net->layers[0] = dataLayer;
	iftDataSetLayerStatus(dataLayer, status);

//    printf("%d %d\n", Z->ntrainsamples, Z->nsamples);

//  printf("%d\n", dataLayer->idx->n);
//	iftPrintIntArray(dataLayer->idx->val, dataLayer->idx->n);

	int offset = 0;
	int i = 0;
	dataLayer->idxCounter = 0;
	while(dataLayer->idxCounter < (dataLayer->idx->n)) {
		iftNetForward(net);
		for(i = 0; i < iftMin(net->miniBatch, dataLayer->idx->n - offset); ++i) {

			int idx = dataLayer->idx->val[i + offset];
			int maxc = 0;
//            printLine(dataLayer->data, i);
//            printf("-- %d %d %d %d %d %d\n", i, dataLayer->idxCounter, idx, Z->sample[idx].id, dataLayer->idx->n, offset);
//            printf("%d\n", Z->sample[idx].id);
//            if(Z->sample[idx].id == 3176)
//			printLine(lossLayer->data, i);

			for (int c = 1; c < lossLayer->data->dimension[1]; c++) {
				if(iftTensorElem(lossLayer->data, i, maxc) < iftTensorElem(lossLayer->data, i, c))
					maxc = c;
			}

			Z->sample[idx].label = maxc+1;
//			printf("%d\\%d\n", Z->sample[idx].label, Z->sample[idx].truelabel);
		}
		offset += i;
	}

	Z->ntrainsamples = iftCountSamples(Z, IFT_TRAIN);

//    if(status != IFT_ALL)
//        iftDataSetLayerStatus(dataLayer, status);

	net->layers[0] = old;
	iftDestroyNeuralLayer(&dataLayer);

}

iftDataSet* iftNeuralNetDataSetFeatures(iftNeuralNetwork* net, int l) {

	iftNeuralLayer* dataLayer = net->layers[0];
	iftNeuralLayer* targetLayer = net->layers[l];
	iftDataSet* Z = dataLayer->dataset;
	iftDataSet* Zout = iftCreateDataSet(Z->nsamples, targetLayer->ninput);

	Zout->nclasses = Z->nclasses;
	Zout->ntrainsamples = Z->ntrainsamples;

	for (int i = 0; i < Z->nsamples; ++i) {
		Zout->sample[i].truelabel = Z->sample[i].truelabel;
		iftSetSampleStatus(&Zout->sample[i], Z->sample[i].status);
	}

	iftSampleStatus status = dataLayer->status;

	if(status != IFT_ALL)
		iftDataSetLayerStatus(dataLayer, IFT_ALL);

	dataLayer->idxCounter = 0;
	while(dataLayer->idxCounter<dataLayer->idx->n) {
		iftNetForward(net);
//    iftPrintTensor(targetLayer->data);
		int curBatch = ((dataLayer->idxCounter-1)%net->miniBatch)+1;
		for(int i = 0; i < curBatch; ++i) {
			int idx = dataLayer->idx->val[dataLayer->idxCounter - curBatch + i];
//      iftPrintIntArray(dataLayer->idx->val, dataLayer->idx->n);
//      printf("i = %d, idx = %d\n", i, idx);
			//int maxc = 0;
			for (int c = 0; c < Zout->nfeats; ++c) {
				Zout->sample[idx].feat[c] = iftTensorElem(targetLayer->data, i, c);
			}

		}
	}

//  Zout->ntrainsamples = iftCountSamples(Z, IFT_TRAIN);

	if(status != IFT_ALL)
		iftDataSetLayerStatus(dataLayer, status);

	return Zout;
}

void iftNetworkAllocateSharedBuffers(iftNeuralNetwork *net)
{
  // Find maximum sizes from layers using shared buffers
  long ioHalfSizeOrig = net->ioHalfSize;
  long long bufferSizeOrig = net->bufferSize;
  for (int i = 0; i < net->nlayers; ++i) {
    iftNeuralLayer *layer = net->layers[i];
    if (layer->isSharedBuffer) {
      net->ioHalfSize = iftMax(layer->data->n, net->ioHalfSize);
      net->ioHalfSize = iftMax(layer->out->n, net->ioHalfSize);
      net->bufferSize = iftMax(layer->bufferSize, net->bufferSize);
    }
  }

  // Actual allocation (if needed)
  if (net->ioHalfSize > ioHalfSizeOrig) {
    iftFree(net->ioBuffer);
    net->ioBuffer = iftAllocFloatArray(net->ioHalfSize * 2);
  }
  if (net->bufferSize > bufferSizeOrig) {
    iftFree(net->buffer);
    net->buffer = iftAllocFloatArray(net->bufferSize / sizeof(float));
  }

  // Make sure layers using shared buffers have the right pointers
  for (int i = 0; i < net->nlayers; ++i) {
    iftNeuralLayer *layer = net->layers[i];
    if (layer->isSharedBuffer) {
      // Alternating addresses for data/out to avoid conflict in consecutive layers
      layer->data->val = (i % 2) ? net->ioBuffer : &(net->ioBuffer[net->ioHalfSize]);
      layer->out->val = (i % 2) ? &(net->ioBuffer[net->ioHalfSize]) : net->ioBuffer;
      // buffer and ibuffer are not used simultaneously so this works for any layer
      layer->buffer = net->buffer;
      layer->ibuffer = (int *) net->buffer;
    }
  }
}

iftNeuralLayer* iftSquaredLossLayer(iftNeuralNetwork *net, int ninput) {

	int minibatch = net->miniBatch;

	iftNeuralLayer* layer = iftCreateNeuralLayer();

	layer->data = iftCreateTensor(2, minibatch, ninput);
	layer->error = iftCreateTensor(2, minibatch, ninput);

	layer->ninput = ninput;
	layer->noutput = ninput;

	layer->forward = iftNoForward;
	layer->backward = iftLossBackward;
	layer->update = iftNoUpdate;

	layer->type = IFT_LOSS_LAYER;

	return layer;
}

void iftNeuralCorrelation(iftTensor *data, iftTensor *kernel, int xstride, int ystride, iftTensor *out) {

  int batchSize = data->dimension[0];

  int nbands = data->dimension[1];
  int width = data->dimension[2];
  int height = data->dimension[3];

  int nkernels = kernel->dimension[0];
//  int kbands = kernel->dimension[1];
  int kwidth = kernel->dimension[2];
  int kheight = kernel->dimension[3];

  int kwidthHalf = (kwidth - 1)/2;
  int kheightHalf = (kheight - 1)/2;


  //image in batch
  for (int m = 0; m < batchSize; ++m) {
    //kernels
    for (int k = 0; k < nkernels; ++k) {
        //convolution
        for (int x = kwidthHalf; x < width - kwidthHalf; x += xstride) {
          for (int y = kheightHalf; y < height - kheightHalf; y += ystride) {
		  double sum = 0.0;
		  //bands
		  for (int b = 0; b < nbands; ++b) {
		          int c = 0;
		        for (int i = - kwidthHalf; i <= kwidthHalf; ++i) {
		          for (int j = - kheightHalf; j <= kheightHalf ; ++j) {
		              c++;
					  sum += iftTensorElem(data, m, b, x + i, y + j) * iftTensorElem(kernel, k, b, i+kwidthHalf, j+kheightHalf);
		          }
		        }
	//              printf("%lf %f %d\n", sum, ((float)(kwidth*kheight)), c);
		        iftTensorElem(out, m, k, (x-kwidthHalf)/xstride, (y-kheightHalf)/ystride) = sum;

		      }
		    }
		  }
    }
  }
}

void iftNeuralConvolution(iftTensor *error, iftTensor *kernel, int xstride, int ystride, iftTensor *out) {

  int batchSize = error->dimension[0];

//  int nbands = kernel->dimension[1];
  int width = error->dimension[2];
  int height = error->dimension[3];

  int nkernels = kernel->dimension[0];
  int nbands = kernel->dimension[1];
  int kwidth = kernel->dimension[2];
  int kheight = kernel->dimension[3];

  int kwidthHalf = (int)(kwidth - 1)/2;
  int kheightHalf = (int)(kheight - 1)/2;

  //image in batch
	for (int m = 0; m < batchSize; ++m) {
		//kernels
		for (int b = 0; b < nbands; ++b) {
			//convolution
			for (int x = 0; x < width * xstride; ++x) {
				for (int y = 0; y < height * ystride; ++y) {
					double sum = 0.0;
					//bands
					for (int i = 0; i < kwidth; ++i) {
						for (int j = 0; j < kheight ; ++j) {
							for (int k = 0; k < nkernels; ++k) {
								int row = x - (i - kwidthHalf)/xstride;
								int col = y - (j - kheightHalf)/ystride;
																
								if((row >= 0) && (col >= 0) && (row < width) && (col < height)){
//									printf("row: %d  col: %d \n", row, col);
									//printf("(%d %d, %d, %d)\n",  i, j, k, b);
									//printf("k = %f\n", iftTensorElem(kernel, k, b, row, col));
									sum += iftTensorElem(kernel, k, b, i, j) * iftTensorElem(error, m, b, row, col)  ;
									
								}
							}
						}
						//printf("sum: %lf\n", sum);
						iftTensorElem(out, m, b, x, y) = sum;						
          			}
        		}
			}
		}
	}
}

void iftConvolutionUpdate(iftNeuralNetwork* net, int l) {
	iftNeuralLayer* layer = net->layers[l];

	const double alpha = net->learnRate;
	//const double gamma = net->weightDecay;
	const double minibatch = net->miniBatch;
	const double decay = net->weightDecay;

	iftLinearCombination(1.0, layer->bias->val, -(alpha / minibatch), layer->biasUpdate->val, layer->bias->val, layer->bias->n);

	if(decay>0)
		iftLinearCombination(1.0, layer->weightUpdate->val, -decay * minibatch, layer->weight->val, layer->weightUpdate->val, layer->weightUpdate->n);

	iftLinearCombination(1.0, layer->weight->val, -(alpha/minibatch), layer->weightUpdate->val, layer->weight->val, layer->weight->n);

}

void iftFullyConnectedForward(iftNeuralNetwork* net, int l) {

	iftNeuralLayer* layer = net->layers[l];

	//weights sum
	iftMatrix* data = iftCreateMatrixPointer(layer->data->val, layer->ninput, net->miniBatch);
	iftMatrix* weight = iftCreateMatrixPointer(layer->weight->val, layer->ninput, layer->noutput);
	iftMatrix* out = iftCreateMatrixPointer(layer->out->val, layer->noutput, net->miniBatch);

	iftMultMatricesInPlace(data, weight, false, true, &out);

	//bias sum
	for (int i = 0; i < net->miniBatch; ++i) {
		for (int j = 0; j < layer->noutput; ++j) {
			iftTensorElem(layer->out, i, j) += iftTensorElem(layer->bias, j);
		}
	}

	//activation
	for (int i = 0; i < layer->out->n; ++i) {
//		printf("%f %f\n", layer->out->val[i], layer->act(layer->out->val[i]));
		layer->out->val[i] = layer->act(layer->out->val[i]);
	}

	iftDestroyMatrix(&data);
	iftDestroyMatrix(&weight);
	iftDestroyMatrix(&out);

}

void iftLossBackward(iftNeuralNetwork* net, int l) {

	iftNeuralLayer* layer = net->layers[l];
	net->error = 0.0f;
	for (int i = 0; i < layer->data->n; ++i) {
		layer->error->val[i] = layer->data->val[i] - net->label->val[i];
		net->error += iftFastPow(layer->error->val[i], 2)/2.0;
	}

}


int getclass(iftTensor* label, int row) {
	int i;
	for (i = 0; i < label->dimension[1]; ++i) {
		if(iftTensorElem(label, row, i)>0.0)
			break;
	}
	return i;
}

void iftDiscriminantLossBackward(iftNeuralNetwork* net, int l) {

	iftNeuralLayer* layer = net->layers[l];

	int nclasses = net->label->dimension[1];

	double* meanClasses = iftAllocDoubleArray(nclasses);
	double* varClasses = iftAllocDoubleArray(nclasses);
	int* nPerClass = iftAllocIntArray(nclasses);

	float meanAll = 0;

	net->error = 0.0;

	for (int j = 0; j < layer->ninput; ++j) {
//		iftTensorElem(layer->error, j, 0) = 0.0;

		meanAll = 0;
		for (int i = 0; i < nclasses; ++i) {
			meanClasses[i] =   0.0;
			nPerClass[i] = 0;
			varClasses[i] = 0.0;
		}

		for (int i = 0; i < net->miniBatch; ++i) {
			float x = iftTensorElem(layer->data, i, j);
			meanAll += x / net->miniBatch;
			meanClasses[getclass(net->label, i)]+= x;
			nPerClass[getclass(net->label, i)]++;
		}

		for (int c = 0; c < nclasses; ++c) {
			meanClasses[c]/=nPerClass[c];
		}

		for (int i = 0; i < net->miniBatch; ++i) {
			int c = getclass(net->label, i);
			float x = iftTensorElem(layer->data, i, j);

			varClasses[c] += iftFastPow(x-meanClasses[c], 2) / nPerClass[c];
		}

		for (int i = 0; i < net->miniBatch; ++i) {
			float x = iftTensorElem(layer->data, i, j);
			int C = getclass(net->label, i);

			float scatterWithinDeriv = 0;
			float scatterWithin = 0;
			float scatterBetweenDeriv = 0;
			float scatterBetween = 0;
			float std = 0;

			for (int c = 0; c < nclasses; ++c) {
				if(nPerClass[c]<=0) {
					printf("Caraio ...\n");
					continue;
				}
				if(c==C) {
					if(nPerClass[c]>1) {
						float newMean = iftDecrementalMean(meanClasses[c], nPerClass[c], x);
						float newVar = iftDecrementalVariance(meanClasses[c], varClasses[c], nPerClass[c], x);
						scatterWithinDeriv = x - newMean;
						scatterWithin = (iftFastPow(scatterWithinDeriv, 2))/ 2.0;
						std += newVar;

						float diff = newMean - iftDecrementalMean(meanAll, net->miniBatch, x);
						scatterBetweenDeriv += diff;
						scatterBetween += iftFastPow(diff, 2)/2.0;
					}
				}
				else {
//					float diff = meanClasses[c] - iftDecrementalMean(meanAll, net->miniBatch, x);
					float diff = x - meanClasses[c];
					scatterBetweenDeriv += diff;
					scatterBetween += (iftFastPow(diff, 4) + 1)/2.0;
//                    std += varClasses[c];
				}

			}

			//float e = (scatterWithinDeriv * scatterBetween + scatterBetweenDeriv * (scatterWithin + std) ) / iftFastPow(scatterBetween, 2);

//			float e = (scatterWithinDeriv + std) / scatterBetween;
//			printf("%f %f %f %f => %f\n", scatterWithin, scatterBetween, std, (float)net->miniBatch, e);

//			iftTensorElem(layer->error, j, 0) += e/(float)(net->miniBatch);
			net->error += (scatterWithin + std) / ((scatterBetween+1)*(float)net->miniBatch);
		}
	}

//	printf("-- %f\n", net->error);

	iftFree(meanClasses);
	iftFree(nPerClass);
}


iftNeuralLayer* iftDiscriminantLossLayer(int minibatch, int ninput) {

	iftNeuralLayer *layer = iftCreateNeuralLayer();
	layer->data = iftCreateTensor(2, minibatch, ninput);
//	layer->error = iftCreateTensor(2, ninput, 1);

	layer->ninput = ninput;
	layer->noutput = 1;

	layer->forward = iftNoForward;
	layer->backward = iftDiscriminantLossBackward;

	return layer;

}

void iftImageForward(iftNeuralNetwork* net, int l) 
{
	iftNeuralLayer* layer = net->layers[l];
	iftNeuralLayer* nextLayer = net->layers[l+1];

	iftReshapeTensor(nextLayer->data, 5, net->miniBatch, layer->ninput, layer->xdata, layer->ydata, layer->zdata);

	for (int i = 0; i < net->miniBatch; ++i) {

		int idx = layer->idx->val[layer->idxCounter];

		iftFile* file = layer->fileset->files[idx];
		iftImage* img = iftReadImageByExt(file->path);

		if(iftIsColorImage(img)) {

			iftImage* channel0 = iftImageRed(img);
			iftImage* channel1 = iftImageGreen(img);
			iftImage* channel2 = iftImageBlue(img);

			for (int x = 0; x < layer->xdata; ++x) {
				for (int y = 0; y < layer->ydata; ++y) {
					for (int z = 0; z < layer->zdata; ++z) {
						iftTensorElem(layer->data, i, 0, x, y, z) = iftImgVal(channel0, x, y, z);
						iftTensorElem(layer->data, i, 1, x, y, z) = iftImgVal(channel1, x, y, z);
						iftTensorElem(layer->data, i, 2, x, y, z) = iftImgVal(channel2, x, y, z);
					}
				}
			}

			iftDestroyImage(&channel0);
			iftDestroyImage(&channel1);
			iftDestroyImage(&channel2);
		}
		else{
			for (int x = 0; x < layer->xdata; ++x) {
				for (int y = 0; y < layer->ydata; ++y) {
					for (int z = 0; z < layer->zdata; ++z) {
						iftTensorElem(layer->data, i, 0, x, y, z) = iftImgVal(img, x, y, z);
					}
				}
			}
		}

		for (int b = 0; b < layer->ninput; ++b) {
			for (int x = 0; x < layer->xdata; ++x) {
				for (int y = 0; y < layer->ydata; ++y) {
					for (int z = 0; z < layer->zdata; ++z) {
                        iftTensorElem(layer->out, i, b, x, y, z) = iftTensorElem(layer->data, i, b, x, y, z);
					}
				}
			}
		}

		for (int c = 0; c < layer->nclasses; ++c) {
			iftTensorElem(net->label, i, c) = 0.00;
			if(file->label == (c+1))
				iftTensorElem(net->label, i, c) = 1.00;
		}

		iftDestroyImage(&img);

		layer->idxCounter++;
		if (layer->idxCounter >= layer->idx->n) { //reaches end, reshuffle and starts again
			if (layer->status == IFT_TRAIN) {
				iftShuffleIntArray(layer->idx->val, layer->idx->n);
				layer->idxCounter = 0;
				net->epoch++;
			}
			else {
				break;
			}
		}
	}
}

void iftDataSetForward(iftNeuralNetwork* net, int l) {
	iftNeuralLayer* layer = net->layers[l];

	for (int i = 0; i < net->miniBatch; ++i) {

		int idx = layer->idx->val[layer->idxCounter];

		for (int j = 0; j < layer->dataset->nfeats; ++j) {
			iftTensorElem(layer->data, i, j) = layer->dataset->sample[idx].feat[j];
//			iftTensorElem(layer->out, i, j) = iftTensorElem(layer->data, i, j);
      		iftTensorElem(layer->out, i, j) = (layer->dataset->sample[idx].feat[j] - iftTensorElem(layer->mean, j))
                                            / (iftTensorElem(layer->std, j) + 1e-15);


//        iftTensorElem(layer->out, i, j)-=127.0;
//        iftTensorElem(layer->out, i, j)/=255.0;
		}

		for (int c = 0; c < layer->dataset->nclasses; ++c) {
//			iftTensorElem(net->label, i, c) = 0.00;
//			if(layer->dataset->sample[idx].truelabel == (c+1))
//				iftTensorElem(net->label, i, c) = 1.00;
			iftTensorElem(net->label, i, c) = layer->dataset->sample[idx].truelabel;
		}

		layer->idxCounter++;
		if (layer->idxCounter >= layer->idx->n) { //reaches end, reshuffle and starts again
			if (layer->status == IFT_TRAIN) {
				iftShuffleIntArray(layer->idx->val, layer->idx->n);
				layer->idxCounter = 0;
				net->epoch++;
			}
			else {
				break;
			}
		}
	}
}

void iftPrintTensor(iftTensor* tensor) {

	if(tensor == NULL)
		return;

	if(tensor->ndimension<=1) {
		iftMatrix* m = iftCreateMatrixPointer(tensor->val, 1, tensor->dimension[0]);
		iftPrintMatrix(m);
		iftDestroyMatrix(&m);
	}
	else if (tensor->ndimension<=2) {
		iftMatrix* m = iftCreateMatrixPointer(tensor->val, tensor->dimension[1], tensor->dimension[0]);
		iftPrintMatrix(m);
		iftDestroyMatrix(&m);
	}
	else if(tensor->ndimension<=3) {
		printf("[");
		for (int i = 0; i < tensor->dimension[0]; ++i) {
			iftMatrix* m = iftCreateMatrixPointer(tensor->val + i*tensor->accDimension[0], tensor->dimension[2], tensor->dimension[1]);
			iftPrintMatrix(m);
			iftDestroyMatrix(&m);
		}
		printf("]");
	}
}

void iftNetForward(iftNeuralNetwork *net)
{
  for (int i = 0; i < net->nlayers; ++i) {
    iftNeuralLayer *layer = net->layers[i];
    layer->forward(net, i);

    // Actually exchange data with following layer
    if (i != net->nlayers - 1) {
      iftNeuralLayer *next = net->layers[i+1];

      if(next->data->n!=layer->out->n) 
        iftError("The output data dimension in layer #%d and the input data dimension in layer #%d do not match.\n", "iftNetForward", i, i+1);

      // Consecutive layers with shared buffers have matching data/out address
      if (!(layer->isSharedBuffer && next->isSharedBuffer))
        iftCopyFloatArray(next->data->val, layer->out->val, layer->out->n);
    }
  }
}

iftNeuralLayer* iftCreateNeuralLayer() {
	iftNeuralLayer* layer = iftAlloc(1, sizeof(iftNeuralLayer));
        layer->isSharedBuffer = false;
	return layer;
}

void iftDestroyNeuralLayer(iftNeuralLayer** l)
{
  assert(l != NULL);
  iftNeuralLayer* layer = *l;

  if (layer->isSharedBuffer) {
    iftDestroyTensorPointer(&(layer->data));
    iftDestroyTensorPointer(&(layer->out));
  } else {
    iftDestroyTensor(&(layer->data));
    iftDestroyTensor(&(layer->out));
    iftFree(layer->buffer);
    iftFree(layer->ibuffer);
  }

  iftDestroyDict(&(layer->imgs));
  iftDestroyTensor(&(layer->weight));
  iftDestroyTensor(&(layer->weightUpdate));
  iftDestroyTensor(&(layer->bias));
  iftDestroyTensor(&(layer->biasUpdate));
  iftDestroyTensor(&(layer->dataMask));
  iftDestroyTensor(&(layer->error));
  iftDestroyTensor(&(layer->delta));
  iftDestroyIntArray(&(layer->idx));
  iftDestroyFileSet(&(layer->fileset));
  iftDestroyDataSet(&(layer->dataset));
  iftDestroyMatrix(&(layer->X));
  iftDestroyMatrix(&(layer->Y));
  iftDestroyTensor(&(layer->mean));
  iftDestroyTensor(&(layer->std));

  iftFree(layer);
  *l = NULL;
}

iftNeuralNetwork* iftCreateNeuralNetwork(int nlayers, int minibatch) {

	iftNeuralNetwork* net = iftAlloc(1, sizeof(iftNeuralNetwork));
	net->nlayers  = nlayers;
	net->curLayer = 0;
	net->layers = (iftNeuralLayer**) iftAlloc(nlayers, sizeof(iftNeuralLayer*));

	net->learnRate = 0.1;
	net->weightDecay = 0.0;

	net->miniBatch = minibatch;

        net->error = 0.0f;
        net->epoch = 0;
        net->iteration = 0;
        net->verbose = false;
        net->label = NULL;
        net->ioBuffer = NULL;
        net->ioHalfSize = 0;
        net->buffer = NULL;
        net->bufferSize = 0;

	return net;
}

void iftDestroyNeuralNetwork(iftNeuralNetwork** n) {

	iftNeuralNetwork* net = *n;

	for (int l = 0; l < net->nlayers; ++l) {
		if(net->layers[l]!=NULL)
			iftDestroyNeuralLayer(&(net->layers[l]));
	}

	iftFree(net->layers);
	iftDestroyTensor(&(net->label));
        iftFree(net->ioBuffer);
        iftFree(net->buffer);

	iftFree(net);
}

// Redundant parts should be merged with the other function later
// Currently only difference is the assumption that fileset_path is on the json
iftNeuralNetwork* iftCreateWeightedNeuralNetworkFromJson(const char *json_path, const char *weights_path, const char *fileset_path, int miniBatchSize, bool isForwardOnly) 
{
  assert(json_path != NULL);

  /* read the general network parameters */
  iftDict *json = iftReadJson(json_path);
  FILE *fp = (weights_path != NULL) ? fopen(weights_path, "rb") : NULL;
  int nLayers = iftGetLongValFromDict("n_layers", json);
  int miniBatch = miniBatchSize > 0 ? miniBatchSize : iftGetDblValFromDict("minibatch_size", json);
  bool hasFileset = (fileset_path != NULL);
  int xSize, ySize, zSize, nBands;
  if (!hasFileset) {
    xSize = iftGetLongValFromDict("x_size", json);
    ySize = iftGetLongValFromDict("y_size", json);
    zSize = iftGetLongValFromDict("z_size", json);
    nBands = iftGetLongValFromDict("n_bands", json);
  }

  iftNeuralNetwork* net = iftCreateNeuralNetwork(nLayers, miniBatch);
  if (iftDictContainKey("learn_rate", json, NULL)) 
    net->learnRate = iftGetDblValFromDict("learn_rate", json);

  /* read the dictionaries for each layer */
  iftDict *layers = iftGetDictFromDict("layers", json);

  for(int l = 0; l < nLayers; l++) {
    char layerName[256];
    sprintf(layerName, "layer_%d", l+1);
    bool hasWeights = false;

    iftDict *layer = iftGetDictFromDict(layerName, layers);
    char *layerType = iftGetStrValFromDict("layer_type", layer);

    /* determine the activation functions (if the layer needs it) */
    iftActivationFunction activ = 0, activDeriv = 0;
    if(iftCompareStrings(layerType, "conv_layer") || iftCompareStrings(layerType, "fc_layer")) {
      char *activFunc = iftGetStrValFromDict("activ_func", layer);
      if(iftCompareStrings(activFunc, "relu")) {
        activ = iftRelu;
        activDeriv = iftReluDeriv;
      }
      else if(iftCompareStrings(activFunc, "leaky_relu")) {
        activ = iftLeakyRelu;
        activDeriv = iftLeakyReluDeriv;
      }
      else if(iftCompareStrings(activFunc, "sigmoid")) {
        activ = iftSigm;
        activDeriv = iftSigmDeriv;
      }
      else if(iftCompareStrings(activFunc, "tanh")) {
        activ = iftTanh;
        activDeriv = iftTanhDeriv;
      }
    }

    /* read the parameters specific for each type of layer */
    if(iftCompareStrings(layerType, "image_layer")) {
      if (!hasFileset) {
        if (isForwardOnly)
          net->layers[l] = iftEmptyImageLayerForwardOnly(net, xSize, ySize, zSize, nBands);
        else
          net->layers[l] = iftEmptyImageLayer(net, xSize, ySize, zSize, nBands);
      } else {
        if (isForwardOnly)
          iftError("ForwardOnly not implemented for fileset based ImageLayer",
              "iftCreatedWeightedNeuralNetworkFromJson");
        iftFileSet *fileset = iftLoadFileSetFromDirOrCSV(fileset_path, 1, true);
        net->layers[l] = iftImageLayer(net, fileset);
        iftDestroyFileSet(&fileset);
      }
    } else if(iftCompareStrings(layerType, "conv_layer")) {
      // Instantiate Convolutional Layer
      int kernelSize = iftGetLongValFromDict("kernel_size", layer);
      int stride = iftGetLongValFromDict("stride", layer);
      int nKernels = iftGetLongValFromDict("n_kernels", layer);
      int pad = (iftDictContainKey("padding", layer, NULL)) ? iftGetLongValFromDict("padding", layer) : 0;
      iftNeuralLayer *prevLayer = net->layers[l-1];
      if (isForwardOnly)
        net->layers[l] = iftConvolutionLayer2DForwardOnly(net, prevLayer->out->dimension[2], prevLayer->out->dimension[3], prevLayer->out->dimension[1], kernelSize, kernelSize, stride, stride, nKernels, pad, activ);
      else
        net->layers[l] = iftConvolutionLayer2D(net, prevLayer->out->dimension[2], prevLayer->out->dimension[3], prevLayer->out->dimension[1], kernelSize, kernelSize, stride, stride, nKernels, pad, activ, activDeriv);
      hasWeights = (weights_path != NULL);
    } else if(iftCompareStrings(layerType, "max_pooling_layer")) {
      // Instantiate Max Pooling Layer
      int poolSize = iftGetLongValFromDict("pool_size", layer);
      int stride = iftGetLongValFromDict("stride", layer);
      iftNeuralLayer *prevLayer = net->layers[l-1];
      if (isForwardOnly)
        net->layers[l] = iftMaxPoolingLayer2DForwardOnly(net, prevLayer->noutput, prevLayer->out->dimension[2], prevLayer->out->dimension[3], poolSize, stride);
      else
        net->layers[l] = iftMaxPoolingLayer2D(net, prevLayer->noutput, prevLayer->out->dimension[2], prevLayer->out->dimension[3], poolSize, stride);
    } else if(iftCompareStrings(layerType, "fc_layer")) {
      // Note that previous implementation assumes only a single FC layer exists
      // We are adding support for multiple ones with the "units" field
      iftNeuralLayer *prevLayer = net->layers[l-1];
      int units = iftGetLongValFromDict("units", layer);
      int nInput = prevLayer->noutput;

      // Implicitly flatten the input for the first FC layer 
      if (prevLayer->xout > 0)
        nInput *= prevLayer->xout;
      if (prevLayer->yout > 0)
        nInput *= prevLayer->yout;
      if (prevLayer->zout > 0)
        nInput *= prevLayer->zout;
      if (isForwardOnly)
        net->layers[l] = iftFullyConnectedLayerForwardOnly(net, nInput, units, activ);
      else
        net->layers[l] = iftFullyConnectedLayer(net, nInput, units, activ, activDeriv);
      hasWeights = (weights_path != NULL);
    }

    if (hasWeights) {
      // Debug, check if we will read the correct amount of data
      int nBytes = 0;
      nBytes += sizeof(*(net->layers[l]->bias->val)) * net->layers[l]->bias->n;
      nBytes += sizeof(*(net->layers[l]->weight->val)) * net->layers[l]->weight->n;
      int expectedNBytes = iftGetLongValFromDict("weights_nbytes", layer);
      if (nBytes != expectedNBytes) {
        printf("Allocated %d bytes but expected %d on layer %d\n", nBytes, expectedNBytes, l);
        exit(-1);
      }

      // Load bias
      int nRead;
      nRead = fread(
          net->layers[l]->bias->val,
          sizeof(*(net->layers[l]->bias->val)),
          net->layers[l]->bias->n,
          fp);
      // Load kernels
      nRead = fread(
          net->layers[l]->weight->val,
          sizeof(*(net->layers[l]->weight->val)),
          net->layers[l]->weight->n,
          fp);
      (void) nRead; // Avoid warnings
    }

    iftDestroyDict(&layer);
  }

  if (net->label == NULL)
    net->label = iftCreateTensor(2, net->miniBatch, net->layers[nLayers-1]->noutput);

  iftDestroyDict(&layers);
  iftDestroyDict(&json);
  fclose(fp);

  if (isForwardOnly)
    iftNetworkAllocateSharedBuffers(net);

  return net;
}

iftNeuralNetwork* iftCreateNeuralNetworkFromJson(const char *filename) {

    /* read the general network parameters */
    iftDict *json = iftReadJson(filename);
    int nLayers = iftGetLongValFromDict("n_layers", json);
    int miniBatch = iftGetLongValFromDict("minibatch_size", json);
    float learnRate = iftGetDblValFromDict("learn_rate", json);
    iftFileSet *fileset = iftLoadFileSetFromDirOrCSV(iftGetStrValFromDict("img_fileset", json), 1, true);
    iftImage *img0 = iftReadImageByExt(fileset->files[0]->path);
    int nLabels = iftFileSetLabelsNumber(fileset);
    
    iftNeuralNetwork* net = iftCreateNeuralNetwork(nLayers, miniBatch);
	net->learnRate = learnRate;

    /* read the dictionaries for each layer */
    iftDict *layers = iftGetDictFromDict("layers", json);
    for(int l = 0; l < nLayers; l++) {
        char layerName[256]; sprintf(layerName, "layer_%d", l+1);
        iftDict *layer = iftGetDictFromDict(layerName, layers);
        char *layerType = iftGetStrValFromDict("layer_type", layer);

        /* determine the activation functions (if the layer needs it) */
        iftActivationFunction activ = 0, activDeriv = 0;
        if(iftCompareStrings(layerType, "conv_layer") || iftCompareStrings(layerType, "fc_layer")) {
            char *activFunc = iftGetStrValFromDict("activ_func", layer);
            if(iftCompareStrings(activFunc, "relu")) {
                activ = iftRelu;
                activDeriv = iftReluDeriv;
            }
            else if(iftCompareStrings(activFunc, "leaky_relu")) {
                activ = iftLeakyRelu;
                activDeriv = iftLeakyReluDeriv;
            }
            else if(iftCompareStrings(activFunc, "sigmoid")) {
                activ = iftSigm;
                activDeriv = iftSigmDeriv;
            }
            else if(iftCompareStrings(activFunc, "tanh")) {
                activ = iftTanh;
                activDeriv = iftTanhDeriv;
            }
        }

        /* read the parameters specific for each type of layer */
        if(iftCompareStrings(layerType, "image_layer")) {
            net->layers[l] = iftImageLayer(net, fileset);
        }
        else if(iftCompareStrings(layerType, "conv_layer")) {
            int kernelSize = iftGetLongValFromDict("kernel_size", layer);
            int stride = iftGetLongValFromDict("stride", layer);
            int nKernels = iftGetLongValFromDict("n_kernels", layer);
            iftNeuralLayer *prevLayer = net->layers[l-1];
            net->layers[l] = iftConvolutionLayer2D(net, prevLayer->out->dimension[2], prevLayer->out->dimension[3], prevLayer->out->dimension[1], kernelSize, kernelSize, stride, stride, nKernels, 0, activ, activDeriv);
        }
        else if(iftCompareStrings(layerType, "max_pooling_layer")) {
            int poolSize = iftGetLongValFromDict("pool_size", layer);
            int stride = iftGetLongValFromDict("stride", layer);
            iftNeuralLayer *prevLayer = net->layers[l-1];
            net->layers[l] = iftMaxPoolingLayer2D(net, prevLayer->noutput, prevLayer->out->dimension[2], prevLayer->out->dimension[3], poolSize, stride);
        }
        else if(iftCompareStrings(layerType, "fc_layer")) {
            iftNeuralLayer *prevLayer = net->layers[l-1];
            net->layers[l] = iftFullyConnectedLayer(net, prevLayer->out->dimension[2] * prevLayer->out->dimension[3] * prevLayer->noutput, nLabels, activ, activDeriv);
        }
        else if(iftCompareStrings(layerType, "loss_layer")) {
            net->layers[l] = iftSquaredLossLayer(net, nLabels);
        }

        iftDestroyDict(&layer);
    }

    iftDestroyFileSet(&fileset);
    iftDestroyImage(&img0);
    iftDestroyDict(&layers);
    iftDestroyDict(&json);

	return net;
}

iftConvNetwork*   iftCreateConvNetwork(int nlayers)
{
	iftConvNetwork *convnet = (iftConvNetwork *) iftAlloc(1,sizeof(iftConvNetwork));

	convnet->nlayers        = nlayers;
	convnet->nstages        = 3*nlayers + 2;
	convnet->input_norm_adj = NULL;
	convnet->input_norm_adj_param = 0;
	convnet->input_xsize    = 0;
	convnet->input_ysize    = 0;
	convnet->input_zsize    = 0;
	convnet->input_nbands   = 0;
	convnet->k_bank         = (iftMMKernel **) iftAlloc(nlayers,sizeof(iftMMKernel *));
	convnet->k_bank_adj_param  = iftAllocIntArray(nlayers);
	convnet->nkernels          = iftAllocIntArray(nlayers);
	convnet->activ_option      = iftAllocIntArray(nlayers);
	convnet->activ_thres		 = iftAllocFloatArray(nlayers);
	convnet->pooling_adj       = (iftAdjRel **) iftAlloc(nlayers, sizeof(iftAdjRel *));
	convnet->pooling_adj_param = iftAllocIntArray(nlayers);
	convnet->stride            = iftAllocIntArray(nlayers);
	convnet->alpha             = iftAllocFloatArray(nlayers);
	convnet->norm_adj          = (iftAdjRel **) iftAlloc(nlayers, sizeof(iftAdjRel *));
	convnet->norm_adj_param    = iftAllocIntArray(nlayers);
	convnet->img_ind           = (iftMatrix **) iftAlloc(convnet->nlayers, sizeof(iftMatrix *));
	convnet->rescale           = 0;
	convnet->with_weights      = 0;

	if ((convnet->k_bank == NULL) ||
		(convnet->pooling_adj == NULL)||
		(convnet->norm_adj == NULL)||
		(convnet->img_ind == NULL)||
		(convnet->activ_option == NULL)||
		(convnet->activ_thres == NULL))
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateConvNetwork");

	for (int l=0; l < nlayers; l++) {
		convnet->k_bank[l]           = NULL;
		convnet->pooling_adj[l]      = NULL;
		convnet->norm_adj[l]         = NULL;
		convnet->img_ind[l]          = NULL;
		convnet->k_bank_adj_param[l]  = 0;
		convnet->norm_adj_param[l]    = 0;
		convnet->pooling_adj_param[l] = 0;
		convnet->activ_thres[l] = -ACTIV_THRES_LIMIT; // nullify threshold - it will be the threshold by band mean
	}

	return(convnet);
}


void iftPrintConvNetArch(iftConvNetwork* convNet) {
	printf("ConvNet { \n\tnlayers = %d\n\tinput = (%d, %d, %d)", convNet->nlayers, convNet->input_xsize, convNet->input_ysize, convNet->input_zsize);
	for (int l = 0; l < convNet->nlayers; ++l) {
		printf("\n\tLayer%d {\n\t\tnkernels = %d\n\t\tkernelSize = %d\n\t\tpoolingSize = %d\n\t\tstrideSize = %d\n\t }", l+1, convNet->nkernels[l], convNet->k_bank_adj_param[l], convNet->pooling_adj_param[l], convNet->stride[l]);
	}
	printf("\n}\n");
}

void iftDestroyConvNetwork(iftConvNetwork **convnet)
{
	iftConvNetwork *aux = *convnet;

	if (aux->input_norm_adj != NULL) iftDestroyAdjRel(&aux->input_norm_adj);

	for (int l=0; l < aux->nlayers; l++) {
		if (aux->k_bank[l] != NULL)      iftDestroyMMKernel(&aux->k_bank[l]);
		if (aux->pooling_adj[l] != NULL) iftDestroyAdjRel(&aux->pooling_adj[l]);
		if (aux->norm_adj[l] != NULL)    iftDestroyAdjRel(&aux->norm_adj[l]);
	}

	for (int i=0; i < aux->nlayers; i++) {
		if (aux->img_ind[i] != NULL) iftDestroyMatrix(&aux->img_ind[i]);
	}

	iftFree(aux->activ_option);
	iftFree(aux->activ_thres);
	iftFree(aux->k_bank);
	iftFree(aux->pooling_adj);
	iftFree(aux->norm_adj);
	iftFree(aux->img_ind);

	iftFree(aux->stride);
	iftFree(aux->alpha);
	iftFree(aux->nkernels);
	iftFree(aux->k_bank_adj_param);
	iftFree(aux->pooling_adj_param);
	iftFree(aux->norm_adj_param);

	iftFree(aux);
	*convnet = NULL;
}


iftMSConvNetwork* iftCreateMSConvNetwork(int nscales)
{
	iftMSConvNetwork *msconvnet = (iftMSConvNetwork *)iftAlloc(1,sizeof(iftMSConvNetwork));

	msconvnet->nscales = nscales;
	msconvnet->convnet = (iftConvNetwork **) iftAlloc(nscales, sizeof(iftConvNetwork *));
	msconvnet->reduction_factor = iftAllocFloatArray(nscales);
	msconvnet->output_norm_adj  = NULL;
	msconvnet->output_norm_adj_param  = 0;

	if (msconvnet->convnet == NULL)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateMSConvNetwork");

	for (int n=0; n < nscales; n++)
		msconvnet->convnet[n] = NULL;

	return(msconvnet);
}

void  iftDestroyMSConvNetwork(iftMSConvNetwork **msconvnet)
{
	iftMSConvNetwork *aux = *msconvnet;

	for (int n=0; n < aux->nscales; n++)
		if (aux->convnet[n] != NULL) iftDestroyConvNetwork(&(aux->convnet[n]));
	iftFree(aux->convnet);
	iftFree(aux->reduction_factor);
	if (aux->output_norm_adj != NULL) iftDestroyAdjRel(&(aux->output_norm_adj));
	iftFree(aux);
	*msconvnet = NULL;
}

iftConvNetwork* iftCreateConvNetworkFromDict(const iftDict* params) {
	int nlayers = iftGetDblValFromDict("n_layers", params);

	iftConvNetwork* conv = NULL;

	conv = iftCreateConvNetwork(nlayers);
	conv->input_xsize = iftGetDblValFromDict("xsize", params);
	conv->input_ysize = iftGetDblValFromDict("ysize", params);
	conv->input_zsize = iftGetDblValFromDict("zsize", params);
	conv->input_nbands = iftGetDblValFromDict("n_bands", params);
	conv->input_norm_adj_param = iftGetDblValFromDict("normalization_size", params);

	char layerKey[IFT_STR_DEFAULT_SIZE];
	for (int l = 0; l < conv->nlayers; ++l) {
		sprintf(layerKey, "layer%d", l+1);
		iftDict* layer = iftGetDictFromDict(layerKey, params);
		conv->nkernels[l] = iftGetDblValFromDict("n_filter", layer);
		conv->k_bank_adj_param[l] = iftGetDblValFromDict("filter_size", layer);
		conv->activ_option[l] = iftGetDblValFromDict("activation", layer);
		conv->pooling_adj_param[l] = iftGetDblValFromDict("pool_size", layer);
		if(conv->activ_option[l])
			conv->activ_thres[l] = iftGetDblValFromDict("activation_thresh", layer);
		conv->alpha[l] = iftGetDblValFromDict("pool_alpha", layer);
		conv->norm_adj_param[l] = iftGetDblValFromDict("normalization_size", layer);
		conv->stride[l] = iftGetDblValFromDict("stride_size", layer);
	}

	conv->with_weights = iftDictContainKey("weight", params, NULL);

	if (!conv->with_weights){ /* create filters with random coefficients */
		iftCreateRandomKernelBanks(conv);
		iftCreateAdjRelAlongNetwork(conv);
		iftMImageIndexMatrices(conv);
	}

	return(conv);
}


iftConvNetwork* iftReadConvNetworkJson(const char *filename) {

	iftConvNetwork* conv = NULL;

	iftDict* params = iftReadJson(filename);

	conv = iftCreateConvNetworkFromDict(params);

	if (conv->with_weights){ /* read kernel weights */
		char* weight = iftGetStrValFromDict("weight", params);
		const char* dir = iftParentDir(filename);
		const char* weightpath = iftJoinPathnames(2, dir, weight);
		FILE* fp = fopen(weightpath, "rb");
		iftLoadKernelBanks(conv, fp);
		fclose(fp);

		iftCreateAdjRelAlongNetwork(conv);
		iftMImageIndexMatrices(conv);
	}

	iftDestroyDict(&params);

	return conv;
}

iftConvNetwork *iftReadConvNetwork(const char *filename)
{
	iftConvNetwork *convnet = NULL;

	const char *ext = iftFileExt(filename);

	if(iftCompareStrings(ext, ".json")) {
		convnet = iftReadConvNetworkJson(filename);
	}
	else if(iftCompareStrings(ext, ".convnet")) {

		FILE *fp = fopen(filename, "rb");
		if (fp == NULL)
            iftError(MSG_FILE_OPEN_ERROR, "iftReadConvNetwork", filename);
		convnet = iftReadConvNetworkFILE(fp);
		fclose(fp);
		return convnet;
	}
	else {
        iftError("Invalid extension: %s", "iftReadConvNetwork", ext);
	}

	return convnet;
}

void  iftWriteConvNetwork(iftConvNetwork *convnet, const char *filename)
{
	FILE *fp = fopen(filename, "wb");
	if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteConvNetwork", filename);
	iftWriteConvNetworkFILE(convnet, fp);
	fclose(fp);
}

bool iftValidArchitecture(iftConvNetwork* convnet) {
	int xsize[convnet->nstages], ysize[convnet->nstages], zsize[convnet->nstages], nbands[convnet->nstages];

	return iftImageDimensionsAlongNetwork(convnet,xsize,ysize,zsize,nbands);
}


// To see the documentation and usage of this method, please check the file: ift_root/demo/iftMergeConvNetworks.c
iftConvNetwork *iftMergeConvNetwork(iftConvNetwork *convnet1, iftConvNetwork *convnet2) {
	int nlayers = convnet1->nlayers + convnet2->nlayers;

	if (convnet1->rescale != convnet2->rescale)
        iftError("Error: Parameter \"RESCALE\" is not equal between convnet_1.convnet and convnet_2.convnet",
                 "iftMergeConvNetwork");

	if (convnet1->with_weights != convnet2->with_weights)
        iftError("Error: Parameter \"WITH_WEIGHTS\" is not equal between convnet_1.convnet and convnet_2.convnet",
                 "iftMergeConvNetwork");


	iftConvNetwork *mconvnet = iftCreateConvNetwork(nlayers);

	puts("#### Copying the convnet_1.convnet Parameters ####");
	mconvnet->input_norm_adj_param 	= convnet1->input_norm_adj_param;
	mconvnet->input_xsize 			= convnet1->input_xsize;
	mconvnet->input_ysize 			= convnet1->input_ysize;
	mconvnet->input_zsize 			= convnet1->input_zsize;
	mconvnet->input_nbands 			= convnet1->input_nbands;
	mconvnet->rescale				    = convnet1->rescale;
	mconvnet->with_weights			= convnet1->with_weights;
	for (int l = 0; l < convnet1->nlayers; l++) {
		mconvnet->k_bank_adj_param[l] 	= convnet1->k_bank_adj_param[l];
		mconvnet->nkernels[l] 		  	= convnet1->nkernels[l];
		mconvnet->activ_option[l] 	  	= convnet1->activ_option[l];
		mconvnet->activ_thres[l] 	  	= convnet1->activ_thres[l];
		mconvnet->pooling_adj_param[l]  = convnet1->pooling_adj_param[l];
		mconvnet->stride[l] 			= convnet1->stride[l];
		mconvnet->alpha[l] 				= convnet1->alpha[l];
		mconvnet->norm_adj_param[l] 	= convnet1->norm_adj_param[l];
	}

	puts("#### Copying the convnet_2.convnet Parameters ####");
	for (int l = convnet1->nlayers; l < nlayers; l++) {
		int idx = l - convnet1->nlayers;
		mconvnet->k_bank_adj_param[l] 	= convnet2->k_bank_adj_param[idx];
		mconvnet->nkernels[l] 			= convnet2->nkernels[idx];
		mconvnet->activ_option[l] 		= convnet2->activ_option[idx];
		mconvnet->activ_thres[l] 		= convnet2->activ_thres[idx];
		mconvnet->pooling_adj_param[l]  = convnet2->pooling_adj_param[idx];
		mconvnet->stride[l] 			= convnet2->stride[idx];
		mconvnet->alpha[l] 				= convnet2->alpha[idx];
		mconvnet->norm_adj_param[l] 	= convnet2->norm_adj_param[idx];
	}

	puts("#### Copying the Kernel weights of the convnet_1.convnet and convnet_2.convnet ####");
	if (mconvnet->with_weights == 1) {
		for (int l = 0; l < convnet1->nlayers; l++)
			mconvnet->k_bank[l] = iftCopyMMKernel(convnet1->k_bank[l]);

		for (int l = convnet1->nlayers; l < nlayers; l++) {
			int idx = l - convnet1->nlayers;
			mconvnet->k_bank[l] = iftCopyMMKernel(convnet2->k_bank[idx]);
		}
	}

	return mconvnet;
}


iftMSConvNetwork* iftReadMSConvNetwork(char *filename)
{
	iftMSConvNetwork *msconvnet=NULL;
	int   s, nscales;
	FILE *fp = fopen(filename,"rb");

	if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadMSConvNetwork", filename);

	iftVerifyToken(fp,"MSCONVNET","iftReadMSConvNetwork");

	/* Read main parameters */

	iftReadIntValue(fp,&nscales,"NSCALES","iftReadMSConvNetwork");
	if (nscales <= 0)
        iftError("Invalid number of scales", "iftReadMSConvNetwork");

	msconvnet = iftCreateMSConvNetwork(nscales);

	iftReadFloatValues(fp,&msconvnet->reduction_factor,nscales,"REDUCTION_FACTOR","iftReadMSConvNetwork");
	iftReadFloatValue(fp,&msconvnet->output_norm_adj_param,"OUTPUT_NORM_ADJ_PARAM","iftReadMSConvNetwork");
	for (int s=0; s < msconvnet->nscales; s++) {
		msconvnet->convnet[s] = iftReadConvNetworkFILE(fp);
	}

	if (msconvnet->convnet[0]->input_zsize != 1 ){ // 3D
		for (s=1; s < nscales; s++) {
			if ((iftRound(msconvnet->convnet[s]->input_xsize / msconvnet->reduction_factor[s]) != msconvnet->convnet[0]->input_xsize / msconvnet->reduction_factor[0]) ||
				(iftRound(msconvnet->convnet[s]->input_ysize / msconvnet->reduction_factor[s]) != msconvnet->convnet[0]->input_ysize / msconvnet->reduction_factor[0]) ||
				(iftRound(msconvnet->convnet[s]->input_zsize / msconvnet->reduction_factor[s]) != msconvnet->convnet[0]->input_zsize / msconvnet->reduction_factor[0])    )
			{
				char msg[200];
				sprintf(msg,"Reduction factor %f at scale %d will create incompatible input image for the convolution network of the next scale",msconvnet->reduction_factor[s],s);
                iftError(msg, "iftReadMSConvNetwork");
			}
		}
	}else { // 2D
		for (s=1; s < nscales; s++) {
			if ((iftRound(msconvnet->convnet[s]->input_xsize / msconvnet->reduction_factor[s]) != msconvnet->convnet[0]->input_xsize / msconvnet->reduction_factor[0]) ||
				(iftRound(msconvnet->convnet[s]->input_ysize / msconvnet->reduction_factor[s]) != msconvnet->convnet[0]->input_ysize / msconvnet->reduction_factor[0])   )
			{
				char msg[200];
				sprintf(msg,"Reduction factor %f at scale %d will create incompatible input image for the convolution network of the next scale",msconvnet->reduction_factor[s],s);
                iftError(msg, "iftReadMSConvNetwork");
			}
		}
	}

	if (msconvnet->output_norm_adj_param > 0.0) {
		if (msconvnet->convnet[0]->input_zsize != 1 )
			msconvnet->output_norm_adj = iftCuboid(msconvnet->output_norm_adj_param,msconvnet->output_norm_adj_param,msconvnet->output_norm_adj_param);
		else
			msconvnet->output_norm_adj = iftRectangular(msconvnet->output_norm_adj_param,msconvnet->output_norm_adj_param);
	}

	fclose(fp);

	return(msconvnet);
}

void  iftWriteMSConvNetwork(iftMSConvNetwork *msconvnet, const char *filename)
{
	FILE *fp = fopen(filename,"wb");

	if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteMSConvNetwork", filename);

	/* writing the token */
	fprintf(fp,"MSCONVNET\n");

	/* Read main parameters */
	iftWriteIntValue(fp,msconvnet->nscales,"NSCALES");

	iftWriteFloatValues(fp,msconvnet->reduction_factor,msconvnet->nscales,"REDUCTION_FACTOR");
	iftWriteFloatValue(fp,msconvnet->output_norm_adj_param,"OUTPUT_NORM_ADJ_PARAM");
	for (int s=0; s < msconvnet->nscales; s++) {
		iftWriteConvNetworkFILE(msconvnet->convnet[s],fp);
	}

	fclose(fp);
}

void iftVerifyNan(iftMImage *img, char *msg, int layer)
{

#pragma omp parallel for shared(img, msg, layer)
	for (int b=0; b < img->m; b++)
		for (int p=0; p < img->n; p++)
			if (isnan(img->val[p][b])){
				char MSG[200];
				sprintf(MSG,"%s in layer %d",msg,layer);
                iftError(MSG, "iftVerifyNan");
			}
}

iftMImage *iftApplyConvNetwork(iftMImage *img, iftConvNetwork *convnet)
{
	iftMImage  *norm_input, *filt_img, *pool_img, *norm_output=NULL;
	int xsize[convnet->nstages], ysize[convnet->nstages], zsize[convnet->nstages], nbands[convnet->nstages];
	int l, i = 0;

	if(!iftImageDimensionsAlongNetwork(convnet,xsize,ysize,zsize,nbands)) {
		char msg[200];
		sprintf(msg,"High network parameters require higher input image resolution. Image at stage %d will have insufficient dimensions: xsize=%d, ysize=%d, zsize=%d, nbands=%d",i,xsize[i],ysize[i],zsize[i],nbands[i]);
        iftError(msg, "iftImageDimensionsAlongNetwork");
	}

	if ( (convnet->input_xsize != img->xsize) ||
		 (convnet->input_ysize != img->ysize) ||
		 (convnet->input_zsize != img->zsize) ||
		 (convnet->input_nbands != img->m) )
        iftError(
                "Image dimensions (%d %d %d, %d bands) are not compatible with the convolution network dimensions (%d %d %d, %d bands)",
                "iftApplyConvNetwork", img->xsize, img->ysize, img->zsize, img->m, convnet->input_xsize,
                convnet->input_ysize,
                convnet->input_zsize, convnet->input_nbands);


	if (convnet->input_norm_adj_param > 0){
		norm_input = iftNormalizeMImage(img, convnet->input_norm_adj);
	}else
		norm_input = img;

	//  iftVerifyNan(norm_input,"after input normalization",0);

	for (l=0, i=2; l < convnet->nlayers; l++, i += 3) {

		filt_img = iftMMLinearFilterConvNetwork(norm_input, convnet->img_ind[l], convnet->k_bank[l], xsize[i], ysize[i], zsize[i], nbands[i]);

		if (convnet->activ_option[l] == 0)
			iftActivateByThres(filt_img, convnet->activ_thres[l]);
		else if (convnet->activ_option[l] == 1)
			iftActivateByBandThres(filt_img);

		if (convnet->pooling_adj_param[l] > 0) {
			pool_img = iftPooling(filt_img, convnet->pooling_adj[l],convnet->stride[l], convnet->alpha[l]);
			iftDestroyMImage(&filt_img);
		}else{
			pool_img = filt_img;
		}

		// iftVerifyNan(pool_img,"after pooling",l);


		if (convnet->norm_adj_param[l] > 0){
			norm_output = iftNormalizeMImage(pool_img, convnet->norm_adj[l]);
			iftDestroyMImage(&pool_img);
		}else{
			norm_output = pool_img;
		}

		//    iftVerifyNan(norm_output,"after output normalization",l);

		if (norm_input != img)
			iftDestroyMImage(&norm_input);

		norm_input = norm_output;

	}

	if (convnet->rescale) {
		iftFastAdjRel *F;
		iftMImage *aux;
		int bx, by, bz;


		for (l=convnet->nlayers-1, i=3*l+2; l >= 0; l--, i -= 3) {

			if (convnet->norm_adj_param[l] > 0) {
				iftMaxAdjShifts(convnet->norm_adj[l],&bx,&by,&bz);
				aux = iftCopyMImage(norm_output);
				iftDestroyMImage(&norm_output);
				norm_output = iftMAddFrame(aux,bx,by,bz,0);
				iftDestroyMImage(&aux);
			}

			if (convnet->stride[l] > 1){ // increase scale
				aux = iftCopyMImage(norm_output);
				iftDestroyMImage(&norm_output);
				if (convnet->input_zsize != 1){ // 3D
					norm_output = iftMInterp(aux,convnet->stride[l],convnet->stride[l],convnet->stride[l]);
				}else{ // 2D
					norm_output = iftMInterp2D(aux,convnet->stride[l],convnet->stride[l]);
				}
				iftDestroyMImage(&aux);
			}

			if (convnet->pooling_adj_param[l] > 0) {
				iftMaxAdjShifts(convnet->pooling_adj[l],&bx,&by,&bz);
				aux = iftCopyMImage(norm_output);
				iftDestroyMImage(&norm_output);
				norm_output = iftMAddFrame(aux,bx,by,bz,0);
				iftDestroyMImage(&aux);
			}

			aux = iftCreateMImage(xsize[i-1],ysize[i-1],zsize[i-1],nbands[i-1]);
			F   = iftCreateFastAdjRel(convnet->k_bank[l]->A,img->tby,img->tbz);
			iftDestroyMImage(&aux);
			aux = iftCopyMImage(norm_output);
			iftDestroyMImage(&norm_output);
			norm_output = iftMAddFrame(aux,F->bx,F->by,F->bz,0);
			iftDestroyFastAdjRel(&F);
			iftDestroyMImage(&aux);
		}

		if (convnet->input_norm_adj_param > 0) {
			iftMaxAdjShifts(convnet->input_norm_adj,&bx,&by,&bz);
			aux = iftCopyMImage(norm_output);
			iftDestroyMImage(&norm_output);
			norm_output = iftMAddFrame(aux,bx,by,bz,0);
			iftDestroyMImage(&aux);
		}

	}

	return(norm_output);

}

iftMImage *iftApplyMSConvNetwork(iftMImage *mimg, iftMSConvNetwork *msconvnet)
{
	int scale;

	iftMImage ** arrmimgs = (iftMImage**) calloc (msconvnet->nscales, sizeof(iftMImage*));
	if (!arrmimgs)
        iftError("Could not allocate array of M Images", "iftExtractDeepFeaturesMultiScale");

	timer *tic, *toc;
	tic = iftTic();
	for (scale = 0; scale < msconvnet->nscales; scale++)
	{
		if (msconvnet->convnet[scale]->input_zsize != 1){ // 3D
			arrmimgs[scale] = iftMInterp   (mimg    , msconvnet->reduction_factor[scale], msconvnet->reduction_factor[scale], msconvnet->reduction_factor[scale]);
		} else {
			arrmimgs[scale] = iftMInterp2D (mimg, msconvnet->reduction_factor[scale], msconvnet->reduction_factor[scale]);
		}
		printf("(%d,%d) -> (%d,%d)\n",iftMXSize(mimg),iftMYSize(mimg),iftMXSize(arrmimgs[scale]),iftMYSize(arrmimgs[scale]));

#ifdef IFT_CONVNET_DEBUG_WRITE_ALL
		iftImage* rimg = iftMImageToImage(arrmimgs[scale],4095,0);

		sprintf(buffer,"output/rescaled%d.pgm",scale);
		iftWriteImageP2(rimg,buffer);
		iftDestroyImage(&rimg);
#endif
		iftMImage* simg = iftApplyConvNetwork(arrmimgs[scale],msconvnet->convnet[scale]);
		iftDestroyMImage(&(arrmimgs[scale]));
		arrmimgs[scale] = simg;
	}
	toc = iftToc();
	printf("multi-scale convnetwork: %f\n", iftCompTime(tic,toc));


	/* rescaling */

	tic = iftTic();
	int m,b,totBands=0;
	for (scale = 0; scale < msconvnet->nscales ; scale++)
		totBands += arrmimgs[scale]->m;
	iftMImage* output = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, totBands);

	for (b=0,scale = 0; scale < msconvnet->nscales ; scale++)
	{
//		fprintf(stderr,"(%3d,%3d,%3d) ->",arrmimgs[scale]->xsize,arrmimgs[scale]->ysize,arrmimgs[scale]->m);
		iftMImage *scaled_image;
		if (msconvnet->convnet[scale]->input_zsize != 1){ // 3D
			scaled_image = iftMInterp  (arrmimgs[scale], (float) mimg->xsize/arrmimgs[scale]->xsize, (float) mimg->ysize/arrmimgs[scale]->ysize, (float) mimg->zsize/arrmimgs[scale]->zsize);
		} else {
			scaled_image = iftMInterp2D(arrmimgs[scale], (float) mimg->xsize/arrmimgs[scale]->xsize, (float) mimg->ysize/arrmimgs[scale]->ysize);
		}
		iftDestroyMImage(&(arrmimgs[scale]));
		arrmimgs[scale] = scaled_image;
        for (m = 0; m < arrmimgs[scale]->m; b++, m++) {
            for (int i = 0; i < output->n; i++) {
                output->val[i][b] = arrmimgs[scale]->val[i][m];
            }
		}
//		fprintf(stderr," (%3d,%3d,%3d)\n",arrmimgs[scale]->xsize,arrmimgs[scale]->ysize,arrmimgs[scale]->m);
	}
	toc = iftToc();
	printf("rescaling: %f\n",iftCompTime(tic,toc));

#ifdef IFT_CONVNET_DEBUG_WRITE_ALL
	sprintf(buffer, "output/final");
	iftWriteMImage(output,buffer);
	buffer[0] = 0;
#endif

	if (msconvnet->output_norm_adj_param > 0)
	{
		tic = iftTic();
		iftMImage *output_normalized = iftNormalizeMImage(output, msconvnet->output_norm_adj);
//		iftMImage *output_normalized = iftNormalizeContrastMImage(output, msconvnet->output_norm_adj);		
//		printf("resolution: (%dx%d)\n",iftMXSize(output_normalized),iftMYSize(output_normalized));
#ifdef IFT_CONVNET_DEBUG_WRITE_ALL
		iftWriteMImage(output_normalized,"output/output_normalized");
#endif
		iftDestroyMImage(&output);
		output = output_normalized;
		toc = iftToc();
		printf("output normalization: %f\n",iftCompTime(tic,toc));
	}

	for(scale=0; scale < msconvnet->nscales; scale++)
		iftDestroyMImage(&(arrmimgs[scale]));
	iftFree(arrmimgs);

	return output;
}


bool iftImageDimensionsAlongNetwork(iftConvNetwork *convnet, int *xsize, int *ysize, int *zsize, int *nbands)
{
	int bx, by, bz, l, i, j, k;
	iftMImage     *img;
	iftFastAdjRel *F;

	/* Compute the number of bands at each stage of each layer along network */

	nbands[0] = convnet->input_nbands;
	nbands[1] = convnet->input_nbands;
	for (l=0, i=2, j=i+1, k=i+2; l < convnet->nlayers; l++,i=i+3,j=i+1,k=i+2){
		nbands[i] = nbands[j] = nbands[k] = convnet->nkernels[l];
	}

	/* Compute the image dimensions at each stage of each layer along network */

	if (convnet->input_norm_adj_param > 0)
		iftMaxAdjShifts(convnet->input_norm_adj,&bx,&by,&bz);
	else
		bx = by = bz = 0;

	/* Before and after input normalization stage */

	xsize[0] = convnet->input_xsize;
	ysize[0] = convnet->input_ysize;
	zsize[0] = convnet->input_zsize;

	xsize[1] = convnet->input_xsize - 2*bx;
	ysize[1] = convnet->input_ysize - 2*by;
	zsize[1] = convnet->input_zsize - 2*bz;

	for (l=0, i=2, j=i+1, k=i+2; l < convnet->nlayers; l++,i=i+3,j=i+1,k=i+2){

		if(xsize[i-1]<=0 || ysize[i-1]<=0 || zsize[i-1]<=0) {
			break;
		}

		/* After filtering */

		img = iftCreateMImage(xsize[i-1],ysize[i-1],zsize[i-1],nbands[i-1]);
		F   = iftCreateFastAdjRel(convnet->k_bank[l]->A,img->tby,img->tbz);

		xsize[i] = xsize[i-1] - 2*F->bx;
		ysize[i] = ysize[i-1] - 2*F->by;
		zsize[i] = zsize[i-1] - 2*F->bz;

		iftDestroyMImage(&img);
		iftDestroyFastAdjRel(&F);

		/* After pooling */

		if (convnet->pooling_adj_param[l] > 0)
			iftMaxAdjShifts(convnet->pooling_adj[l],&bx,&by,&bz);
		else
			bx = by = bz = 0;

		xsize[j] = xsize[i] - 2*bx;
		ysize[j] = ysize[i] - 2*by;
		zsize[j] = zsize[i] - 2*bz;

		if (convnet->stride[l] > 1){ // reduce scale
			xsize[j] = ceilf(((float)xsize[j])/convnet->stride[l]);
			ysize[j] = ceilf(((float)ysize[j])/convnet->stride[l]);
			if (convnet->input_zsize != 1) // 3D
				zsize[j] = ceilf(((float)zsize[j])/convnet->stride[l]);
		}

		/* After normalization */

		if (convnet->norm_adj_param[l] > 0)
			iftMaxAdjShifts(convnet->norm_adj[l],&bx,&by,&bz);
		else
			bx = by = bz = 0;

		xsize[k] = xsize[j] - 2*bx;
		ysize[k] = ysize[j] - 2*by;
		zsize[k] = zsize[j] - 2*bz;

	}

	for (i=0; i < convnet->nstages; i++){
		if ((nbands[i] <= 0) || (xsize[i] <= 0) || (ysize[i] <= 0)|| (zsize[i] <= 0))
		{
			return 0;
		}
	}

	for (l = 0; l < convnet->nlayers; ++l) {
		size_t x = xsize[2+l*3];
		size_t y = ysize[2+l*3];
		size_t z = zsize[2+l*3];
		size_t nbands = (size_t)convnet->k_bank[l]->nbands;
		size_t im_memory = x*y*z;
		size_t kernel_memory= nbands*(size_t)(convnet->k_bank[l]->A->n);
		if(iftLog(im_memory, 2)+iftLog(kernel_memory, 2)>31)//exceeds the maximum contiguous memory limit
			return 0;
	}

	for (l = 0; l < convnet->nlayers; ++l) {
		size_t x = xsize[2+l*3];
		size_t y = ysize[2+l*3];
		size_t z = zsize[2+l*3];
		size_t nbands = (size_t)convnet->k_bank[l]->nbands;
		size_t im_memory = x*y*z;
		size_t kernel_memory= nbands*(size_t)(convnet->k_bank[l]->A->n);
		if(iftLog(im_memory, 2)+iftLog(kernel_memory, 2)+iftLog(sizeof(double), 2)>=32)//exceeds the maximum contiguous memory limit
			return 0;
	}

	return 1;

}

iftMMKernel *iftUnsupLearnKernelsFromDataSet(iftDataSet* Z,iftAdjRel* A, float kmax_perc, bool whitening)
{
	iftDataSet  *Zp,*Zk=NULL;
	iftKnnGraph *graph=NULL;
	iftMMKernel *K;
	int          nkernels,nbands;

	if ((kmax_perc >  1.0)||(kmax_perc < 0))
        iftError("Invalid percentage of adjacent samples", "iftUnsupLearnKernelsFromDataSet");


	if (whitening)
		Zp = iftWhiteningTransform(Z);
	else
		Zp = Z;

	iftSetDistanceFunction(Zp,1);

	nbands = Zp->nfeats/A->n;

	/* Compute clustering using normalized graph-cut for best k estimation */

	if (Zp->nsamples <= 500) {
		graph = iftCreateKnnGraph(Zp,(int) iftMax((kmax_perc * Zp->nsamples), 5));
		iftUnsupTrain(graph, iftNormalizedCut);
	} else {
		iftSelectUnsupTrainSamples(Zp,500.0/(float)Zp->nsamples,100);
		graph = iftUnsupLearn(Zp, kmax_perc, iftNormalizedCut, 100, false);
	}

	/* Select kernels from root samples */
	nkernels = 0;
	for (int u=0; u < graph->nnodes; u++)  {
		if (graph->node[u].root==u) {
			nkernels ++;
		}
	}

	//  printf("Creating bank with %d kernels \n",nkernels);

	K = iftCreateMMKernel(A,nbands,nkernels); //img->m
	int k=0;
	for (int u=0; u < graph->nnodes; u++)  {
		if (graph->node[u].root == u) {
			int s = graph->node[u].sample;
			int j = 0;
			for (int b=0; b < nbands; b++) { // img->m
				for (int i=0; i < A->n; i++) {
					K->weight[k][b].val[i] = Zp->sample[s].feat[j];
					j++;
				}
			}
			k++;
		}
	}


	if (Zp->fsp.W != NULL){
		K->W     = iftCopyMatrix(Zp->fsp.W);
		K->mean  = iftAllocFloatArray(Zp->nfeats);
		for (int i=0; i < Zp->nfeats; i++) {
			K->mean[i]  = Zp->fsp.mean[i];
		}
	}

	if (graph != NULL)
		iftDestroyKnnGraph(&graph);
	if (whitening)
		iftDestroyDataSet(&Zp);

	if (Zk != NULL)
		iftDestroyDataSet(&Zk);

	return(K);
}

iftMMKernel *iftUnsupLearnKernelsByKmeansFromDataSet(iftDataSet* Z,iftAdjRel* A, int k, char whitening)
{
	iftDataSet  *Zp,*Zk=NULL;
	iftMMKernel *K;
	int          nkernels,nbands;

	if ((k <= 1))
        iftError("Invalid number of kernels", "iftUnsupLearnKernelsByKmensFromDataSet");


	if (whitening)
		Zp = iftWhiteningTransform(Z);
	else
		Zp = Z;

	iftSetDistanceFunction(Zp,1);

	nbands = Zp->nfeats/A->n;

	/* Initialization of Centroids */

	if (whitening)
		Zk = iftKmeansInitCentroidsRandomNormal(Zp,k);
	else
		Zk = iftKmeansInitCentroidsFromSamples(Zp,k);

	/* Kmeans */

	iftKmeansRun(0,Zp, &Zk, 50,1E-5);

	nkernels = Zk->nsamples;

	//  printf("Creating bank with %d kernels \n",nkernels);

	K = iftCreateMMKernel(A,nbands,nkernels); //img->m

	for (int s=0; s < Zk->nsamples; s++)  {
		int j = 0;
		for (int b=0; b < nbands; b++) {
			for (int i=0; i < A->n; i++) {
				K->weight[s][b].val[i] = Zk->sample[s].feat[j];
				j++;
			}
		}

	}

	if (Zp->fsp.W != NULL){
		K->W     = iftCopyMatrix(Zp->fsp.W);
		K->mean  = iftAllocFloatArray(Zp->nfeats);
		for (int i=0; i < Zp->nfeats; i++) {
			K->mean[i]  = Zp->fsp.mean[i];
		}
	}

	if (whitening)
		iftDestroyDataSet(&Zp);

	if (Zk != NULL)
		iftDestroyDataSet(&Zk);

	return(K);
}

iftMMKernel *iftUnsupLearnKernelsByKmedoidsFromDataSet(iftDataSet* Z,iftAdjRel* A, int k, char whitening)
{
	iftDataSet  *Zp,*Zk=NULL;
	iftMMKernel *K;
	int          nkernels,nbands;

	if ((k <= 1))
        iftError("Invalid number of kernels", "iftUnsupLearnKernelsByKmedoidsFromDataSet");


	if (whitening)
		Zp = iftWhiteningTransform(Z);
	else
		Zp = Z;

	iftSetDistanceFunction(Zp,1);

	nbands = Zp->nfeats/A->n;

	/* Initialization of Centroids */

	if (whitening)
		Zk = iftKmeansInitCentroidsRandomNormal(Zp,k);
	else
		Zk = iftKmeansInitCentroidsFromSamples(Zp,k);

	/* Kmeans */

	iftKmeansRun(1,Zp, &Zk, 50,1E-5);

	nkernels = Zk->nsamples;

	//  printf("Creating bank with %d kernels \n",nkernels);

	K = iftCreateMMKernel(A,nbands,nkernels); //img->m

	for (int s=0; s < Zk->nsamples; s++)  {
		int j = 0;
		for (int b=0; b < nbands; b++) {
			for (int i=0; i < A->n; i++) {
				K->weight[s][b].val[i] = Zk->sample[s].feat[j];
				j++;
			}
		}
	}

	if (Zp->fsp.W != NULL){
		K->W     = iftCopyMatrix(Zp->fsp.W);
		K->mean  = iftAllocFloatArray(Zp->nfeats);
		for (int i=0; i < Zp->nfeats; i++) {
			K->mean[i]  = Zp->fsp.mean[i];
		}
	}

	if (whitening)
		iftDestroyDataSet(&Zp);

	if (Zk != NULL)
		iftDestroyDataSet(&Zk);

	return(K);
}

iftMMKernel *iftUnsupLearnKernelsBySphericalKmeansFromDataSet(iftDataSet* Z,iftAdjRel* A, int k, char whitening)
{
	iftDataSet  *Zp,*Zk=NULL;
	iftMMKernel *K;
	int          nkernels,nbands;

	if ((k <= 1))
        iftError("Invalid number of kernels", "iftUnsupLearnKernelsByKmedoidsFromDataSet");


	if (whitening)
		Zp = iftWhiteningTransform(Z);
	else
		Zp = Z;

	iftSetDistanceFunction(Zp,1);

	nbands = Zp->nfeats/A->n;

	/* Initialization of Centroids */

	if (whitening)
		Zk = iftKmeansInitCentroidsRandomNormal(Zp,k);
	else
		Zk = iftKmeansInitCentroidsFromSamples(Zp,k);

	/* Kmeans */

	iftSphericalKmeansRun(Zp, &Zk, 50);

	nkernels = Zk->nsamples;

	//  printf("Creating bank with %d kernels \n",nkernels);

	K = iftCreateMMKernel(A,nbands,nkernels); //img->m

	for (int s=0; s < Zk->nsamples; s++)  {
		int j = 0;
		for (int b=0; b < nbands; b++) {
			for (int i=0; i < A->n; i++) {
				K->weight[s][b].val[i] = Zk->sample[s].feat[j];
				j++;
			}
		}
	}

	if (Zp->fsp.W != NULL){
		K->W     = iftCopyMatrix(Zp->fsp.W);
		K->mean  = iftAllocFloatArray(Zp->nfeats);
		for (int i=0; i < Zp->nfeats; i++) {
			K->mean[i]  = Zp->fsp.mean[i];
		}
	}

	if (whitening)
		iftDestroyDataSet(&Zp);

	if (Zk != NULL)
		iftDestroyDataSet(&Zk);

	return(K);
}


void iftUnsupLearnKernels(iftMImage *img, iftConvNetwork *convnet, int nsamples, float kmax_perc, bool whitening)
{
	iftDataSet  *Z;
	iftMMKernel *K;

	if (convnet->nlayers != 1)
        iftError("It requires a single-layer convolution network", "iftUnsupLearnKernels");


	Z  = iftSelectPatches(img, convnet,nsamples);


	K = iftUnsupLearnKernelsFromDataSet(Z,convnet->k_bank[0]->A,kmax_perc,whitening);


	iftDestroyMMKernel(&convnet->k_bank[0]);
	convnet->k_bank[0] = K;
	convnet->nkernels[0] = K->nkernels;

	iftDestroyDataSet(&Z);
}

void iftUnsupLearnKernelsByKmeans(iftMImage *img, iftConvNetwork *convnet, int nsamples, int k, char whitening)
{
	iftDataSet  *Z;
	iftMMKernel *K;

	if (convnet->nlayers != 1)
        iftError("It requires a single-layer convolution network", "iftUnsupLearnKernelsByKmedoids");


	Z  = iftSelectPatches(img, convnet,nsamples);


	K = iftUnsupLearnKernelsByKmeansFromDataSet(Z,convnet->k_bank[0]->A,k,whitening);


	iftDestroyMMKernel(&convnet->k_bank[0]);
	convnet->k_bank[0] = K;
	convnet->nkernels[0] = K->nkernels;

	iftDestroyDataSet(&Z);
}

void iftUnsupLearnKernelsByKmedoids(iftMImage *img, iftConvNetwork *convnet, int nsamples, int k, char whitening)
{
	iftDataSet  *Z;
	iftMMKernel *K;

	if (convnet->nlayers != 1)
        iftError("It requires a single-layer convolution network", "iftUnsupLearnKernelsByKmedoids");


	Z  = iftSelectPatches(img, convnet,nsamples);


	K = iftUnsupLearnKernelsByKmedoidsFromDataSet(Z,convnet->k_bank[0]->A,k,whitening);


	iftDestroyMMKernel(&convnet->k_bank[0]);
	convnet->k_bank[0] = K;
	convnet->nkernels[0] = K->nkernels;

	iftDestroyDataSet(&Z);
}

void iftSelectRandomKernels(iftMImage *img, iftConvNetwork *convnet, char whitening)
{
	iftDataSet  *Z, *Zp;
	iftMMKernel *K;

	if (convnet->nlayers != 1)
        iftError("It requires a single-layer convolution network", "iftSelectRandomKernels");

	Z  = iftSelectPatches(img, convnet, convnet->nkernels[0]);

	if (whitening)
		Zp = iftWhiteningTransform(Z);
	else
		Zp = Z;

	convnet->nkernels[0] = Zp->nsamples;

	K = iftCreateMMKernel(convnet->k_bank[0]->A,convnet->k_bank[0]->nbands,convnet->nkernels[0]);

	int s = 0;
	for (int k=0; k < K->nkernels; k++) {
		int j = 0;
		for (int b=0; b < K->nbands; b++) {
			for (int i=0; i < K->A->n; i++) {
				K->weight[k][b].val[i] = Zp->sample[s].feat[j];
				j++;
			}
		}
		s++;
	}

	if (Zp->fsp.W != NULL){
		K->W     = iftCopyMatrix(Zp->fsp.W);
		K->mean  = iftAllocFloatArray(Zp->nfeats);
		for (int i=0; i < Zp->nfeats; i++) {
			K->mean[i]  = Zp->fsp.mean[i];
		}
	}


	iftDestroyMMKernel(&convnet->k_bank[0]);
	convnet->k_bank[0]   = K;
	if (Zp != Z)
		iftDestroyDataSet(&Zp);
	iftDestroyDataSet(&Z);
}




void iftConvertHyperplaneInKernelBand(iftSVMHyperplane *h, iftBand *kernel, int band_size) {

	for (int i = 0; i < h->n; i++) {
		int b = i/band_size;
		int idx = i % band_size;

		kernel[b].val[idx] = h->feat[i];
	}
}



// Return the n learned hyperplanes, where n is the number of classes from dataset.
// Compute the sum of all predicted values from samples in the hyperplanes
iftSVMHyperplane **iftSupLearnKernelsFromDataSet(iftDataSet* Z, float *predict_sum) {
	iftSVM *svm = NULL;
	float C = 1e5;
	int nW = Z->nfeats; //size of hyperplane
	int nclasses = Z->nclasses; // num_hyperplanes = nclasses

	iftSVMHyperplane **hps = (iftSVMHyperplane**) iftAlloc(nclasses, sizeof(iftSVMHyperplane*));
	for (int i = 0; i < nclasses; i++)
		hps[i] = iftCreateSVMHyperplane(nW);

	iftSetDistanceFunction(Z, 1);
	iftSetStatus(Z, IFT_TRAIN);

	svm = iftCreateLinearSVC(C);
	iftSVMTrainOVA(svm, Z, Z);

	// Extract the optimal hyperplane to each class
//	#pragma omp paralell for
	for (int classIdx = 0; classIdx < nclasses; classIdx++) {
		float bias;
		iftSample s = iftSVMGetNormalHyperplaneAsSample(svm, classIdx, Z, &bias);
		iftCopyFloatArray(hps[classIdx]->feat, s.feat, nW);
		hps[classIdx]->bias = bias;
	}

	iftSetStatus(Z, IFT_TEST);
	iftMatrix *predict_matrix = NULL;
	iftSVMLinearClassifyOVA(svm, Z, IFT_TEST, &predict_matrix);
	float sum = 0.0;
//	#pragma omp paralell for reduction(+:sum)
	for (int s = 0; s < Z->nsamples; s++)
		for (int c = 0; c < Z->nclasses; c++) {
			sum += predict_matrix->val[iftGetMatrixIndex(predict_matrix, c, s)];
		}
	*predict_sum = sum;

	// Deallocators
	iftDestroyMatrix(&predict_matrix);

	iftFree(predict_matrix);

	iftDestroySVM(svm);

	return hps;
}

void iftCreateAdjRelAlongNetwork(iftConvNetwork *convnet)
{
	if (convnet->input_zsize != 1) { // 3D
		if (convnet->input_norm_adj_param > 0)
			convnet->input_norm_adj = iftCuboid(convnet->input_norm_adj_param,convnet->input_norm_adj_param,convnet->input_norm_adj_param);
		for (int l=0; l < convnet->nlayers; l++) {
			if (convnet->pooling_adj_param[l] > 0)
				convnet->pooling_adj[l] = iftCuboid(convnet->pooling_adj_param[l],convnet->pooling_adj_param[l],convnet->pooling_adj_param[l]);
			if (convnet->norm_adj_param[l] > 0)
				convnet->norm_adj[l]    = iftCuboid(convnet->norm_adj_param[l],convnet->norm_adj_param[l],convnet->norm_adj_param[l]);
		}
	}else{ // 2D
		if (convnet->input_norm_adj_param > 0)
			convnet->input_norm_adj = iftRectangular(convnet->input_norm_adj_param,convnet->input_norm_adj_param);
		for (int l=0; l < convnet->nlayers; l++) {
			if (convnet->pooling_adj_param[l] > 0)
				convnet->pooling_adj[l] = iftRectangular(convnet->pooling_adj_param[l],convnet->pooling_adj_param[l]);
			if (convnet->norm_adj_param[l] > 0)
				convnet->norm_adj[l]    = iftRectangular(convnet->norm_adj_param[l],convnet->norm_adj_param[l]);
		}
	}
}


void iftMImageIndexMatrices(iftConvNetwork *convnet)
{
	int xsize[convnet->nstages], ysize[convnet->nstages], zsize[convnet->nstages], nbands[convnet->nstages];
	iftMImage     *img;
	iftFastAdjRel *F;
	int            l,i=0;

	/* Compute image dimensions along the network */

	if(!iftImageDimensionsAlongNetwork(convnet,xsize,ysize,zsize,nbands)) {
		char msg[200];
		sprintf(msg,"High network parameters require higher input image resolution. Image at stage %d will have insufficient dimensions: xsize=%d, ysize=%d, zsize=%d, nbands=%d",i,xsize[i],ysize[i],zsize[i],nbands[i]);
        iftError(msg, "iftImageDimensionsAlongNetwork");
	}

	for (l=0, i=2; l < convnet->nlayers; l++,i=i+3)
	{

		/* For filtering */

		img = iftCreateMImage(xsize[i-1],ysize[i-1],zsize[i-1],nbands[i-1]);
		F   = iftCreateFastAdjRel(convnet->k_bank[l]->A,img->tby,img->tbz);
		convnet->img_ind[l] = iftMImageIndicesToMatrix(img,F,1);
		iftDestroyMImage(&img);
		iftDestroyFastAdjRel(&F);

	}

}

void iftCreateRandomKernelBanks(iftConvNetwork *convnet)
{
	iftAdjRel *A;
	int nbands[convnet->nlayers];

	nbands[0] = convnet->input_nbands;
	for (int l=1; l < convnet->nlayers; l++)
		nbands[l] = convnet->nkernels[l-1];

	for (int l=0; l < convnet->nlayers; l++) {
		if (convnet->input_zsize != 1)
			A = iftCuboid(convnet->k_bank_adj_param[l],convnet->k_bank_adj_param[l],convnet->k_bank_adj_param[l]);
		else
			A = iftRectangular(convnet->k_bank_adj_param[l],convnet->k_bank_adj_param[l]);

		convnet->k_bank[l] = iftRandomMMKernel(A, nbands[l], convnet->nkernels[l]);
		iftDestroyAdjRel(&A);
	}

}



void iftLoadKernelBanks(iftConvNetwork *convnet, FILE *fp)
{
	for (int l=0; l < convnet->nlayers; l++) {
		convnet->k_bank[l] = iftReadMMKernelFILE(fp);
	}
}

iftConvNetwork *iftReadConvNetworkFILE(FILE* fp)
{
	iftConvNetwork *convnet=NULL;
	int   nlayers;

	iftVerifyToken(fp,"CONVNET","iftReadConvNetwork");

	/* Read main parameters */

	iftReadIntValue(fp,&nlayers,"NLAYERS","iftReadConvNetwork");
	if (nlayers <= 0)
        iftError("Invalid number of layers", "iftReadConvNetwork");


	convnet = iftCreateConvNetwork(nlayers);

	iftReadIntValue(fp,&convnet->input_norm_adj_param,"INPUT_NORM_ADJ_PARAM","iftReadConvNetwork");
	iftReadIntValue(fp,&convnet->input_xsize,"INPUT_XSIZE","iftReadConvNetwork");
	iftReadIntValue(fp,&convnet->input_ysize,"INPUT_YSIZE","iftReadConvNetwork");
	iftReadIntValue(fp,&convnet->input_zsize,"INPUT_ZSIZE","iftReadConvNetwork");
	iftReadIntValue(fp,&convnet->input_nbands,"INPUT_NBANDS","iftReadConvNetwork");
	iftReadIntValues(fp,&convnet->k_bank_adj_param,nlayers,"K_BANK_ADJ_PARAM","iftReadConvNetwork");
	iftReadIntValues(fp,&convnet->nkernels,nlayers,"NKERNELS","iftReadConvNetwork");
	iftReadIntValues(fp, &convnet->activ_option, nlayers, "ACTIV_OPTION", "iftReadConvNetwork");
	iftReadFloatValues(fp, &convnet->activ_thres, nlayers, "ACTIV_THRES", "iftReadConvNetwork");
	iftReadIntValues(fp,&convnet->pooling_adj_param,nlayers,"POOLING_ADJ_PARAM","iftReadConvNetwork");
	iftReadIntValues(fp,&convnet->stride,nlayers,"STRIDE","iftReadConvNetwork");
	iftReadFloatValues(fp,&convnet->alpha,nlayers,"ALPHA","iftReadConvNetwork");
	iftReadIntValues(fp,&convnet->norm_adj_param,nlayers,"NORM_ADJ_PARAM","iftReadConvNetwork");
	iftReadIntValue(fp,&convnet->rescale,"RESCALE","iftReadConvNetwork");
	iftReadIntValue(fp,&convnet->with_weights,"WITH_WEIGHTS","iftReadConvNetwork");

	if ((convnet->input_norm_adj_param < 0) ||
		(convnet->input_nbands <= 0) ||
		(convnet->input_xsize < 3) ||
		(convnet->input_ysize < 3) ||
		(convnet->input_zsize < 1) ||
		(convnet->rescale < 0) ||
		(convnet->with_weights < 0))
        iftError("Invalid input parameters", "iftReadConvNetwork");

	for (int l=0; l < convnet->nlayers; l++) {
		if ((convnet->activ_option[l] != 0) && (convnet->activ_option[l] != 1)) {
			char msg[200];
			sprintf(msg, "Invalid option for Activation in the layer %d (ACTIV_OPTION)", l);
            iftError(msg, "iftReadConvNetwork");
		}

		if ((convnet->k_bank_adj_param[l] <= 0) ||
			(convnet->nkernels[l] < 1) ||
			(convnet->pooling_adj_param[l] < 0) ||
			(convnet->stride[l] < 1) ||
			(convnet->alpha[l] < 1) ||
			(convnet->norm_adj_param[l] < 0)){
			char msg[100];
			sprintf(msg,"Invalid input parameters for layer %d",l);
            iftError(msg, "iftReadConvNetwork");
		}
	}

	/* generate random filters or read kernel weights */

	if (convnet->with_weights == 0){ /* create filters with random coefficients */
		iftCreateRandomKernelBanks(convnet);
	}else if (convnet->with_weights == 1){ /* read kernel weights */
		iftLoadKernelBanks(convnet, fp);
	}


	/* Create adjacency relations */

	iftCreateAdjRelAlongNetwork(convnet);


	/* Create image index matrices */

	iftMImageIndexMatrices(convnet);

	return(convnet);
}

iftConvNetwork *iftCopyFirstLayersConvNetwork(iftConvNetwork* convnet, int nlayers)
{
	iftConvNetwork *out_convnet;

	if (nlayers <= 0)
        iftError("Invalid number of layers", "iftReadConvNetwork");

	out_convnet = iftCreateConvNetwork(nlayers);

	out_convnet->input_norm_adj_param  = convnet->input_norm_adj_param;
	out_convnet->input_xsize           = convnet->input_xsize;
	out_convnet->input_ysize           = convnet->input_ysize;
	out_convnet->input_zsize           = convnet->input_zsize;
	out_convnet->input_nbands          = convnet->input_nbands;

	for(int l=0;l<nlayers;l++) {
		out_convnet->k_bank_adj_param [l] = convnet->k_bank_adj_param [l];
		out_convnet->nkernels         [l] = convnet->nkernels         [l];
		out_convnet->pooling_adj_param[l] = convnet->pooling_adj_param[l];
		out_convnet->stride           [l] = convnet->stride           [l];
		out_convnet->alpha            [l] = convnet->alpha            [l];
		out_convnet->norm_adj_param   [l] = convnet->norm_adj_param   [l];

		out_convnet->k_bank           [l] = iftCopyMMKernel(convnet->k_bank[l]);
	}
	out_convnet->rescale         = convnet->rescale;
	out_convnet->with_weights    = convnet->with_weights;
	out_convnet->activ_option    = convnet->activ_option;
	out_convnet->activ_thres     = convnet->activ_thres;

	/* Create adjacency relations */
	iftCreateAdjRelAlongNetwork(out_convnet);

	/* Create image index matrices */
	iftMImageIndexMatrices(out_convnet);

	return(out_convnet);
}


float iftImg2ColGetPixel(float *im, int height, int width, int channels,
					   int row, int col, int channel, int pad)
{
	row -= pad;
	col -= pad;

	if (row < 0 || col < 0 ||
		row >= height || col >= width) return 0;
	return im[col + width*(row + height*channel)];
}

int iftImg2ColGetIdx(int height, int width, int channels,
						 int row, int col, int channel, int pad)
{
	row -= pad;
	col -= pad;

	if (row < 0 || col < 0 ||
		row >= height || col >= width) return IFT_NIL;
	return col + width*(row + height*channel);
}

//This piece of code is from Caffe
//https://github.com/BVLC/caffe/blob/master/LICENSE
int iftImg2Col(float* data_im,
				int channels,  int height,  int width,
				int ksize,  int stride, int pad, float* data_col)
{
	int c,h,w;
	int height_col = (height + 2*pad - ksize) / stride + 1;
	int width_col = (width + 2*pad - ksize) / stride + 1;

	int n = 0;
	int channels_col = channels * ksize * ksize;
	for (c = 0; c < channels_col; ++c) {
		int w_offset = c % ksize;
		int h_offset = (c / ksize) % ksize;
		int c_im = c / ksize / ksize;
		for (h = 0; h < height_col; ++h) {
			for (w = 0; w < width_col; ++w) {
				int im_row = h_offset + h * stride;
				int im_col = w_offset + w * stride;
				int col_index = (c * height_col + h) * width_col + w;
//				printf("(%d, %d) col_index: %d\n", im_row, im_col, col_index);
				float value = iftImg2ColGetPixel(data_im, height, width, channels,
												  im_row, im_col, c_im, pad);
				data_col[col_index] = value;
				n++;
			}
		}
	}

	return n;
}

void iftImg2ColOptimized(float* inputImg, int nBands,  int ySize,  int xSize,
    int kSize, int stride, int pad, float* outMxBuffer)
{
  // Effective amount of kernel placements along each axis
  int ySizeMod = (ySize + 2 * pad - kSize) / stride + 1;
  int xSizeMod = (xSize + 2 * pad - kSize) / stride + 1;

  int nRowsOut = nBands * kSize * kSize;
  // int nColsOut = xSizeMod * ySizeMod * batchSize;

  for (int outRow = 0; outRow < nRowsOut; ++outRow) {
    // Each row corresponds to position within the kernel
    int xOffsetFromKernel = outRow % kSize;
    int yOffsetFromKernel = (outRow / kSize) % kSize;
    int imgBand = outRow / (kSize * kSize);
    
    // Separate intervals with padded image values
    // This avoids a conditional in the innermost loop
    int yNoPadStart = iftMax(pad - yOffsetFromKernel, 0);
    int yNoPadEnd = iftMin(ySizeMod + pad - yOffsetFromKernel, ySizeMod);
    int xNoPadStart = iftMax(pad - xOffsetFromKernel, 0);
    int xNoPadEnd = iftMin(xSizeMod + pad - xOffsetFromKernel, xSizeMod);

    // Fill actual image values
    for (int yOut = yNoPadStart; yOut < yNoPadEnd; ++yOut) {
      for (int xOut = xNoPadStart; xOut < xNoPadEnd; ++xOut) {
        int xImg = (xOffsetFromKernel + xOut * stride) - pad;
        int yImg = (yOffsetFromKernel + yOut * stride) - pad;
        // (xImg, yImg, imgBand) collapsed to linear index
        int imgIdx = xImg + xSize * (yImg + imgBand * ySize);

        // Each col corresponds to a kernel placement over the input image
        // (xOut, yOut, outRow) collapsed to linear index
        int outIdx = (outRow * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = inputImg[imgIdx];
      }
    }

    // Fill pad values
    for (int yOut = 0; yOut < yNoPadStart; ++yOut) {
      for (int xOut = 0; xOut < xSizeMod; ++xOut) {
        int outIdx = (outRow * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = 0.0f;
      }
    }
    for (int yOut = yNoPadEnd; yOut < ySizeMod; ++yOut) {
      for (int xOut = 0; xOut < xSizeMod; ++xOut) {
        int outIdx = (outRow * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = 0.0f;
      }
    }
    for (int yOut = 0; yOut < ySizeMod; ++yOut) {
      for (int xOut = 0; xOut < xNoPadStart; ++xOut) {
        int outIdx = (outRow * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = 0.0f;
      }
      for (int xOut = xNoPadEnd; xOut < xSizeMod; ++xOut) {
        int outIdx = (outRow * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = 0.0f;
      }
    }
  }
}

void iftImg2ColOptimizedBatch(float* inputImg, int nBands,  int ySize,  int xSize,
    int kSize, int stride, int pad, int batchSize, int batchPos, float* outMxBuffer)
{
  // Effective amount of kernel placements along each axis
  int ySizeMod = (ySize + 2 * pad - kSize) / stride + 1;
  int xSizeMod = (xSize + 2 * pad - kSize) / stride + 1;

  int nRowsOut = nBands * kSize * kSize;
  // int nColsOut = xSizeMod * ySizeMod * batchSize;

  for (int outRow = 0; outRow < nRowsOut; ++outRow) {
    // Each row corresponds to position within the kernel
    int xOffsetFromKernel = outRow % kSize;
    int yOffsetFromKernel = (outRow / kSize) % kSize;
    int imgBand = outRow / (kSize * kSize);
    
    // Separate intervals with padded image values
    // This avoids a conditional in the innermost loop
    int yNoPadStart = iftMax(pad - yOffsetFromKernel, 0);
    int yNoPadEnd = iftMin(ySizeMod + pad - yOffsetFromKernel, ySizeMod);
    int xNoPadStart = iftMax(pad - xOffsetFromKernel, 0);
    int xNoPadEnd = iftMin(xSizeMod + pad - xOffsetFromKernel, xSizeMod);

    // Fill actual image values
    for (int yOut = yNoPadStart; yOut < yNoPadEnd; ++yOut) {
      for (int xOut = xNoPadStart; xOut < xNoPadEnd; ++xOut) {
        int xImg = (xOffsetFromKernel + xOut * stride) - pad;
        int yImg = (yOffsetFromKernel + yOut * stride) - pad;
        // (xImg, yImg, imgBand) collapsed to linear index
        int imgIdx = xImg + xSize * (yImg + imgBand * ySize);

        // Each col corresponds to a kernel placement over the input image
        // (xOut, yOut, batchPos, outRow) collapsed to linear index
        int outIdx = ((outRow * batchSize + batchPos) * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = inputImg[imgIdx];
      }
    }

    // Fill pad values
    for (int yOut = 0; yOut < yNoPadStart; ++yOut) {
      for (int xOut = 0; xOut < xSizeMod; ++xOut) {
        int outIdx = ((outRow * batchSize + batchPos) * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = 0.0f;
      }
    }
    for (int yOut = yNoPadEnd; yOut < ySizeMod; ++yOut) {
      for (int xOut = 0; xOut < xSizeMod; ++xOut) {
        int outIdx = ((outRow * batchSize + batchPos) * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = 0.0f;
      }
    }
    for (int yOut = 0; yOut < ySizeMod; ++yOut) {
      for (int xOut = 0; xOut < xNoPadStart; ++xOut) {
        int outIdx = ((outRow * batchSize + batchPos) * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = 0.0f;
      }
      for (int xOut = xNoPadEnd; xOut < xSizeMod; ++xOut) {
        int outIdx = ((outRow * batchSize + batchPos) * ySizeMod + yOut) * xSizeMod + xOut;
        outMxBuffer[outIdx] = 0.0f;
      }
    }
  }
}


int iftImg2ColIdx(int channels,  int height,  int width,
			   int ksize,  int stride, int pad, int* idx_col)
{
	int c,h,w;
	int height_col = (height + 2*pad - ksize) / stride + 1;
	int width_col = (width + 2*pad - ksize) / stride + 1;

	int n = 0;
	int channels_col = channels * ksize * ksize;
	for (c = 0; c < channels_col; ++c) {
		int w_offset = c % ksize;
		int h_offset = (c / ksize) % ksize;
		int c_im = c / ksize / ksize;
		for (h = 0; h < height_col; ++h) {
			for (w = 0; w < width_col; ++w) {
				int im_row = h_offset + h * stride;
				int im_col = w_offset + w * stride;
				int col_index = (c * height_col + h) * width_col + w;
				int value = iftImg2ColGetIdx(height, width, channels,
												 im_row, im_col, c_im, pad);
//				printf("(%d, %d) col_index: %d, %d\n", im_row, im_col, col_index, value);

				idx_col[col_index] = value;
				n++;
			}
		}
	}

	return n;
}


void iftCol2ImgSetPixel(float *im, int height, int width, int channels,
					  int row, int col, int channel, int pad, float val)
{
	row -= pad;
	col -= pad;

	if (row < 0 || col < 0 ||
		row >= height || col >= width) return;
	im[col + width*(row + height*channel)] += val;
}


void iftCol2Img(float* data_col,
				int channels,  int height,  int width,
				int ksize,  int stride, int pad, float* data_im)
{
	int c,h,w;
	int height_col = (height + 2*pad - ksize) / stride + 1;
	int width_col = (width + 2*pad - ksize) / stride + 1;

	int channels_col = channels * ksize * ksize;
	for (c = 0; c < channels_col; ++c) {
		int w_offset = c % ksize;
		int h_offset = (c / ksize) % ksize;
		int c_im = c / ksize / ksize;
		for (h = 0; h < height_col; ++h) {
			for (w = 0; w < width_col; ++w) {
				int im_row = h_offset + h * stride;
				int im_col = w_offset + w * stride;
				int col_index = (c * height_col + h) * width_col + w;
				double val = data_col[col_index];
				iftCol2ImgSetPixel(data_im, height, width, channels,
								 im_row, im_col, c_im, pad, val);
			}
		}
	}
}
