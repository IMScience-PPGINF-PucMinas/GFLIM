//
// Created by Alan Peixinho on 4/12/15.
//

#include "iftBagOfFeatures.h"

#include "ift/core/io/Stream.h"


#define MAX_FEATS 1000000

bool iftBowSetup(iftBagOfFeatures *bow) {

    iftImage* testImage = NULL;
    float* testFeats = NULL;

    int *patchesTest = NULL;
    int nSamplesTest = 0;

    int nPatchesFeatsTest = 0;
    iftMImage* mImage = NULL;

    /**
     * Check if the patch fits inside image
     */
    if(bow->patchSize.x>bow->imgSize.x || bow->patchSize.y>bow->imgSize.y || bow->patchSize.z>bow->imgSize.z) {
        iftError("Patches size must be less or equal image size.", "iftValidateParamsBOW");
        return 1;
    }

    /**
     * Check if the number of features returned by FeatureExtractor matches the indicated in nPatchFeats
     */
    if(bow->isColored)
      testImage = iftCreateColorImage(bow->imgSize.x, bow->imgSize.y, bow->imgSize.z,8);
    else
        testImage = iftCreateImage(bow->imgSize.x, bow->imgSize.y, bow->imgSize.z);

    testFeats = iftAllocFloatArray(MAX_FEATS);

    iftVoxel zero = (iftVoxel) {.x=0, .y=0, .z=0};

    nPatchesFeatsTest = bow->featsExtractor(testImage, zero, bow->patchSize, testFeats);

    bow->nPatchfeats = nPatchesFeatsTest;

    /**
     * Check if number of patches extracted by Sampler matches the indicated in samplesPerImage
     */
    patchesTest = iftAllocIntArray(bow->samplesPerImage);
    nSamplesTest = bow->sampler(testImage, bow->patchSize, bow->samplesPerImage, patchesTest);

    if(nSamplesTest !=bow->samplesPerImage)
    {
        iftError("Wrong number of samples per image. Indicated %d, but found %d", "iftValidateParamsBOW",
                 bow->samplesPerImage,
                 nSamplesTest);
        return 1;
    }

    /**
     * Check if the feats number returned by the Coding matches the indicated in nFeats
     */
    mImage = iftCreateMImage((bow->imgSize.x - bow->patchSize.x)/bow->stride.x,
                             (bow->imgSize.y - bow->patchSize.y)/bow->stride.y,
                             iftMax((bow->imgSize.z - bow->patchSize.z) / bow->stride.z, 1),
                             bow->dictionarySize);

    iftBowCodingPooling coding = bow->codingPooling;
    int nFeatsTest = coding(mImage, testFeats, bow->nfeats);

    bow->nfeats = nFeatsTest;

    iftDestroyMImage(&mImage);
    iftDestroyImage(&testImage);
    iftFree(patchesTest);
    iftFree(testFeats);

    return 0;
}

int iftBowPatchesRandomSampler(iftImage* img, iftVoxel size, int nsamples, int *patchesOut)
{
    iftVoxel v;

    for (int i = 0; i < nsamples; ++i) {
        v.x = iftRandomInteger(0, img->xsize - size.x);
        v.y = iftRandomInteger(0, img->ysize - size.y);
        v.z = iftRandomInteger(0, img->zsize - size.z);

        patchesOut[i] = iftGetVoxelIndex(img, v);
    }

    return nsamples;
}

int iftBowPatchesDenseSampler(iftImage* img, iftVoxel size, int nsamples, int *patchesOut)
{
    int count = 0;
    iftVoxel v;

    for (v.z = 0; v.z <= img->zsize - size.z; v.z+=size.z) {
        for (v.y = 0; v.y <= img->ysize - size.y; v.y+=size.y) {
            for (v.x = 0; v.x <= img->xsize - size.x; v.x+=size.x) {

                int idx = iftGetVoxelIndex(img, v);

                if(count>=nsamples)
                {
                    count++;
                    continue;
                }

                patchesOut[count++] = idx;
            }
        }
    }

    return count;
}

iftDataSet* iftBowKMeansKernels(iftDataSet* z, int nkernels){

    iftDataSet* kernels = iftKmeansInitCentroidsRandomNormal(z, nkernels);
    iftSimpleKmeansRun(z, &kernels, 20, .001);

    return kernels;
}

iftDataSet* iftBowSupKMeansKernels(iftDataSet* z, int nkernels) {

    int nclasses = z->nclasses;

    if(nkernels < nclasses) {
        iftError("Number of kernels should greater or equal than the number of classes.", "iftBowSupKMeansKernels");
    }

    int *kernelsPerClass = iftAllocIntArray(nclasses);
    int reminder = nkernels - nclasses*(nkernels/nclasses);

    for (int i = 0; i < nclasses; ++i) {
        kernelsPerClass[i] = nkernels/nclasses;
    }

    for (int i = 0; i < reminder; ++i) {
        kernelsPerClass[i]++;
    }

    int sum = 0;
    for (int i = 0; i < nclasses; ++i) {
        sum+=kernelsPerClass[i];
    }
    assert(sum==nkernels);

    iftDataSet** zc = (iftDataSet**) iftAlloc(nclasses, sizeof(iftDataSet*));
    iftDataSet** classKernels = (iftDataSet**) iftAlloc(nclasses, sizeof(iftDataSet*));

    //#pragma omp parallel for
    for (int i = 0; i < nclasses; ++i) {
        zc[i] = iftExtractClass(z, i+1);
        classKernels[i] = iftKmeansInitCentroidsRandomNormal(z, kernelsPerClass[i]);
        printf("Running %dth Kmeans k=%d ... (nsamples = %d, nfeats = %d)\n", i+1, kernelsPerClass[i], zc[i]->nsamples, zc[i]->nfeats);
        iftKmeansRun(kernelsPerClass[i], zc[i], &classKernels[i], 1000, .001);
    }

    iftDataSet* kernels = iftMergeDataSetArray(classKernels, nclasses, true);

    for (int i = 0; i < nclasses; ++i) {
        iftDestroyDataSet(&classKernels[i]);
        iftDestroyDataSet(&zc[i]);
    }

    iftFree(classKernels);
    iftFree(zc);

    return kernels;
}

int iftBowPatchRawFeatsExtractor(iftImage* img, iftVoxel begin, iftVoxel patchSize, float* featsOut)
{
    int count = 0;
    iftVoxel v;

    for(v.z = begin.z; v.z < begin.z+patchSize.z; ++v.z) {
        for (v.y = begin.y; v.y < begin.y + patchSize.y; ++v.y) {
            for (v.x = begin.x; v.x < begin.x + patchSize.x; ++v.x) {
                featsOut[count++] = iftImgVoxelVal(img, v) / 255.0;
            }
        }
    }

    return count;
}


int iftBowPatchColorRawFeatsExtractor(iftImage* img, iftVoxel begin, iftVoxel patchSize, float* featsOut) {
    int count = 0, idx;
    iftVoxel v;

    for(v.z = begin.z; v.z < begin.z+patchSize.z; ++v.z) {
        for (v.y = begin.y; v.y < begin.y + patchSize.y; ++v.y) {
            for (v.x = begin.x; v.x < begin.x + patchSize.x; ++v.x) {
                idx = iftGetVoxelIndex(img, v);
                featsOut[count++] = img->val[idx] / 255.0;
                featsOut[count++] = img->Cb[idx] / 255.0;
                featsOut[count++] = img->Cr[idx] / 255.0;
            }
        }
    }

    return count;
}

int iftBowPatchBICFeatsExtractor(iftImage *img, iftVoxel begin, iftVoxel patchSize, float *featsOut) {

    iftVoxel end = (iftVoxel) iftVectorSum(begin, patchSize);

    int nbins = 6;
    int nfeats;

    iftBoundingBox bb = {.begin = begin, .end = end};
    iftImage* roi = iftExtractROI(img, bb);
    iftFeatures* feats = iftExtractBIC(roi, NULL, nbins);
    nfeats = feats->n;

    for (int i = 0; i < feats->n; ++i) {
        featsOut[i] = feats->val[i];
    }

    iftDestroyFeatures(&feats);
    iftDestroyImage(&roi);

    return nfeats;
}


int iftNoCoding(iftMImage *mImg, float *featsOut, int n) {

    int nbands = mImg->m;
    int count = 0;

    for (int band = 0; band < nbands; ++band) {
        for (int p = 0; p < mImg->n; ++p) {
            featsOut[count++] = mImg->val[p][band];
        }
    }

    return count;
}

int iftBowHardCoding(iftMImage *mImg, float *featsOut, int n) {

    int nbands = mImg->m;

    for (int band = 0; band < nbands; ++band) {
        featsOut[band] = 0.0;
    }

    for (int p = 0; p < mImg->n; ++p) {
        int maxIdx = 0;
        for (int band = 0; band < nbands; ++band) {
            if(mImg->val[p][band]>mImg->val[p][maxIdx])
                maxIdx = band;
        }
        featsOut[maxIdx]++;
    }

    iftUnitNorm(featsOut, nbands);

    return nbands;
}

int iftBowSoftCoding(iftMImage *mImg, float *featsOut, int n) {

    for (int b = 0; b < mImg->m; b++)
        featsOut = 0;

    for (int i = 0; i < mImg->n; i++) {
        for (int b = 0; b < mImg->m; b++) {
            featsOut[b] += mImg->val[i][b];
        }
    }

    iftUnitNorm(featsOut, mImg->m);

    return mImg->m;
}



int iftBowHVSumPooling(iftMImage *mImg, float *featsOut, int n) {
    int nbands = mImg->m;
    // int count = 0;

    memset(featsOut, 0, n*sizeof(float));

    for (int band = 0; band < nbands; ++band) {

        int hstart = band*(mImg->xsize+mImg->ysize);
        int vstart = band*(mImg->xsize+mImg->ysize) + mImg->xsize;
        iftVoxel v;
        v.z = 0;
        for (v.x = 0; v.x < mImg->xsize; ++v.x) {
            for (v.y = 0; v.y < mImg->ysize; ++v.y) {
                float value = mImg->val[iftMGetVoxelIndex(mImg, v)][band];
                featsOut[hstart + v.x] += value;
                featsOut[vstart + v.y] += value;
            }
        }
    }
    
    return nbands*(mImg->xsize+mImg->ysize);
}

iftBagOfFeatures *iftCreateBow(bool color, int xImgSize, int yImgSize, int zImgSize,
                     int xPatchSize, int yPatchSize, int zPatchSize,
                     int xStrideSize, int yStrideSize, int zStrideSize, iftPatchesSampler sampler,
                     iftPatchFeatsExtractor featsExtractor, iftBowKernelsEstimator dictionaryEstimator,
                     iftBowCodingPooling coding, int samplesPerImage, int dictSize) {

    iftBagOfFeatures* bow = (iftBagOfFeatures *) iftAlloc(1, sizeof(iftBagOfFeatures));

    bow->isColored = color;

    bow->imgSize.x = xImgSize; bow->imgSize.y = yImgSize; bow->imgSize.z = zImgSize;
    bow->patchSize.x = xPatchSize; bow->patchSize.y = yPatchSize; bow->patchSize.z = zPatchSize;
    bow->stride.x = xStrideSize; bow->stride.y = yStrideSize; bow->stride.z = zStrideSize;

    bow->dictionarySize = dictSize;
    bow->featsExtractor = featsExtractor;
    bow->samplesPerImage = samplesPerImage;
    bow->sampler = sampler;


    bow->dictionary = NULL;
    bow->dictionaryEstimator = dictionaryEstimator;
    bow->codingPooling = coding;

    return bow;
}

void iftDestroyBow(iftBagOfFeatures **bow) {

    if((*bow)->dictionary!=NULL)
    {
        iftDestroyDataSet(&((*bow)->dictionary));
    }

    iftFree(*bow);
}

void iftBowWrite(iftBagOfFeatures *bow) {
    iftError("Not Implemented yet. Shame on you lazy programmer.", "iftBowWrite");
}

iftBagOfFeatures *iftBowRead(const char *filename) {
    
    iftError("Not Implemented yet. Shame on you lazy programmer.", "iftBowRead");

    return NULL;
}

void iftBowLearn(iftBagOfFeatures *bow, iftFileSet *images) {

    int* imgSamplesBuffer = iftAllocIntArray(bow->samplesPerImage);

    int nImages = 0;

    for (int i = 0; i < images->n; ++i)
        if(images->files[i]->status == IFT_TRAIN) nImages++;

    if(nImages<=0) {
        iftError("There are no training images in the dataset.", "iftBOWLearn");
    }

    iftDataSet* samples = iftCreateDataSet(bow->samplesPerImage *nImages, bow->nPatchfeats);

    int idx, label;
    int samplesNumber = 0;

//    printf("Create samples dataset.\n");

    /**
     * Creates samples dataset
     */
    #pragma omp parallel for
    for (int imgIdx = 0; imgIdx < images->n; ++imgIdx) {

        if(images->files[imgIdx]->status!=IFT_TRAIN)
            continue;

        label = images->files[imgIdx]->label;
        idx = images->files[imgIdx]->sample;
        iftImage* img = iftReadImageByExt(images->files[imgIdx]->path);

        if(img->xsize != bow->imgSize.x || img->ysize!=bow->imgSize.y) {
            iftError("Wrong input size. Expected (%d,%d,%d) image.", "iftBowLearn");//TODO: lancar mensagem de erro
        }

        bow->sampler(img, bow->patchSize, bow->samplesPerImage, imgSamplesBuffer);

        for (int i = 0; i < bow->samplesPerImage; ++i) {
            iftVoxel v = iftGetVoxelCoord(img, imgSamplesBuffer[i]);
            bow->featsExtractor(img, v, bow->patchSize, samples->sample[samplesNumber].feat);

            samples->sample[samplesNumber].id = idx;
            samples->sample[samplesNumber].truelabel = label;

            samplesNumber++;
        }

        iftDestroyImage(&img);
    }

    iftCountNumberOfClassesDataSet(samples);

//    printf("Run Dictionary Learn.\n");
    iftBowKernelsEstimator kernelsEstimator = bow->dictionaryEstimator;
    //create dictionary
    iftDataSet* dictionary = kernelsEstimator(samples, bow->dictionarySize);

    bow->dictionary = dictionary;
}

//TODO: speed up with GPU
iftMImage *iftBowTransform(iftBagOfFeatures *bow, iftImage *img) {

    iftMImage* mImage = iftCreateMImage((bow->imgSize.x - bow->patchSize.x)/bow->stride.x,
                                        (bow->imgSize.y - bow->patchSize.y)/bow->stride.y,
                                        iftMax((bow->imgSize.z - bow->patchSize.z) / bow->stride.z, 1), bow->dictionarySize);

    float *bufferFeats = iftAllocFloatArray(bow->nPatchfeats);
    iftVoxel vmimg, vimg;
    int p;

    for (int kernelIdx = 0; kernelIdx < bow->dictionarySize; ++kernelIdx) {

        float* kernel = bow->dictionary->sample[kernelIdx].feat;

        for(vmimg.z = 0; vmimg.z < mImage->zsize; vmimg.z+=1) {
            vimg.z = vmimg.z*bow->stride.z;
            for (vmimg.y = 0; vmimg.y < mImage->ysize; vmimg.y+=1) {
                vimg.y = vmimg.y*bow->stride.y;
                for (vmimg.x = 0; vmimg.x < mImage->xsize; vmimg.x+=1) {
                    vimg.x = vmimg.x*bow->stride.x;
                    bow->featsExtractor(img, vimg, bow->patchSize, bufferFeats);
                    float dist = iftSquaredFeatDistance(kernel, bufferFeats, bow->nPatchfeats);
                    p = iftGetVoxelIndex(mImage, vmimg);
                    mImage->val[p][kernelIdx] = 1.0 / (1.0 + dist);/** Inverse of distance to the kernels */
                }
            }
        }
    }

    iftFree(bufferFeats);

    return mImage;
}

iftDataSet *iftBowExtractFeaturesBatch(iftBagOfFeatures *bow, iftFileSet *imgArray) {
    iftBowCodingPooling coding = bow->codingPooling;
    iftDataSet* Z = iftCreateDataSet(imgArray->n, bow->nfeats);

    #pragma omp parallel for
    for (int i = 0; i < imgArray->n; ++i) {
        iftImage* img = iftReadImageByExt(imgArray->files[i]->path);
        iftMImage* mImg = iftBowTransform(bow, img);
        coding(mImg, Z->sample[i].feat, Z->nfeats);
        Z->sample[i].truelabel = imgArray->files[i]->label;
        iftSetSampleStatus(&Z->sample[i], imgArray->files[i]->status);
        iftDestroyImage(&img);
        iftDestroyMImage(&mImg);
    }

    iftCountNumberOfClassesDataSet(Z);

    int ntrain = 0;
    for (int i = 0; i < Z->nsamples; ++i) {
        if(iftHasSampleStatus(Z->sample[i], IFT_TRAIN))
            ntrain++;
    }
    Z->ntrainsamples = ntrain;

    return Z;
}

iftFeatures *iftBowExtractFeatures(iftBagOfFeatures *bow, iftImage *img) {
    iftFeatures* features = iftCreateFeatures(bow->nfeats);

    if(img->xsize != bow->imgSize.x || img->ysize!=bow->imgSize.y) {
        iftError("Wrong input size. Expected (%d,%d,%d) images.", "iftBowExtractFeatures", bow->imgSize.x,
                 bow->imgSize.y, bow->imgSize.z);
        return NULL;
    }

    iftMImage* mImage = iftBowTransform(bow, img);

    iftBowCodingPooling coding = bow->codingPooling;
    coding(mImage, features->val, features->n);
    iftDestroyMImage(&mImage);

    return features;
}
