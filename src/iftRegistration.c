#include "iftRegistration.h"

#include "ift/core/io/Dir.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include <math.h>


//MSPS parameters
#define R_X 0
#define R_Y 1
#define R_Z 2
#define T_X 3
#define T_Y 4
#define T_Z 5

#define MAX_SCORE 1000000

/********************** PRIVATE FUNCTIONS *************************/




/********************** PUBLIC FUNCTIONS *************************/

iftFileSet *iftCentralizeLabelImages(const iftFileSet *label_imgs_files, const char *output_dir) {
    iftVoxel *gcs        = NULL;
    iftBoundingBox *mbbs = iftAllMinBoundingBox(label_imgs_files, &gcs);

    iftImageDomain dom = {0, 0, 0};

    // finds the resulting image domain, which is the maximum bounding box that fits all min. bounding boxes
    // centralized by their geometric centers
    for (size_t i = 0; i < label_imgs_files->n; i++) {
        int xsize = iftRound(2 * iftMax(mbbs[i].end.x - gcs[i].x + 1, gcs[i].x - mbbs[i].begin.x + 1));
        int ysize = iftRound(2 * iftMax(mbbs[i].end.y - gcs[i].y + 1, gcs[i].y - mbbs[i].begin.y + 1));
        int zsize = iftRound(2 * iftMax(mbbs[i].end.z - gcs[i].z + 1, gcs[i].z - mbbs[i].begin.z + 1));

        if (dom.xsize < xsize)
            dom.xsize = xsize;
        if (dom.ysize < ysize)
            dom.ysize = ysize;
        if (dom.zsize < zsize)
            dom.zsize = zsize;
    }

    iftFileSet *out_label_img_files = iftCreateFileSet(label_imgs_files->n);
    char *tmp_dir = NULL;
    if (output_dir == NULL)
        tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);
    else tmp_dir = iftCopyString(output_dir);


    // centralize the labeled images
    #pragma omp parallel for
    for (size_t i = 0; i < label_imgs_files->n; i++) {
        char *filename                = iftFilename(label_imgs_files->files[i]->path, NULL);
        char *out_path                = iftJoinPathnames(2, tmp_dir, filename);
        out_label_img_files->files[i] = iftCreateFile(out_path);

        iftImage *label_img         = iftReadImageByExt(label_imgs_files->files[i]->path);
        iftImage *central_label_img = iftCreateImage(dom.xsize, dom.ysize, dom.zsize);
        iftCopyVoxelSize(label_img, central_label_img);
        iftVoxel center             = iftImageCenter(central_label_img);

        iftVoxel u;
        for (u.z = mbbs[i].begin.z; u.z <= mbbs[i].end.z; u.z++) {
            for (u.y = mbbs[i].begin.y; u.y <= mbbs[i].end.y; u.y++) {
                for (u.x = mbbs[i].begin.x; u.x <= mbbs[i].end.x; u.x++) {
                    iftVoxel disp = iftVectorSub(center, gcs[i]);
                    iftVoxel v    = iftVectorSum(u, disp);

                    iftImgVoxelVal(central_label_img, v) = iftImgVoxelVal(label_img, u);
                }
            }
        }
        iftWriteImageByExt(central_label_img, out_path);

        iftDestroyImage(&label_img);
        iftDestroyImage(&central_label_img);
        iftFree(filename);
        iftFree(out_path);
    }


    iftFree(mbbs);
    iftFree(gcs);
    iftFree(tmp_dir);


    return out_label_img_files;
}
/*************************************************************************************/



typedef struct ift_ShapeRegisterProb
{
    int npoints;
    iftPoint *points;
    iftPoint center;

    iftImage *score;
    iftMatrix*_bufferT;/** Speed up allocating transform matrix only once **/

} iftShapeRegisterProb;
//
//iftImage* euclideanDistanceTransform(iftSet* img){
//
//   iftImage* gradient = iftImageGradientMagnitude(img, iftCircular(25.0));
//   printf("%d, %d\n", iftMinimumValue(gradient), iftMaximumValue(gradient));
//
//   float mean = 0.0f;
//   for (int i = 0; i < gradient->n; ++i) {
//      mean += (float)gradient->val[i]/gradient->n;
//   }
//
//   for (int i = 0; i < gradient->n; ++i) {
//      gradient->val[i] = 255*(gradient->val[i] > mean + 10);
//   }
//
//   iftImage* edt = iftEuclDistTrans(gradient, iftCircular(1.5), IFT_BOTH, NULL, NULL, NULL);
//   iftImage* normedt = iftNormalize(edt, 0, 255);
//
//   normedt = iftCopyImage(gradient);
//
//   iftDestroyImage(&gradient);
//   iftDestroyImage(&edt);
//
//   return normedt;
//}

iftPoint iftGetPointsCenter(const iftPoint* p, int npoints)
{
    int i;
    iftPoint center;
    double x, y, z;
    x = y = z = 0.0f;

    for (i = 0; i<npoints; ++i) {
        x += p[i].x;
        y += p[i].y;
        z += p[i].z;
    }

    center.x=x/npoints; center.y=y/npoints; center.z=z/npoints;

    return center;
}


iftImage* iftEuclideanScoreImage(iftImage* img, float decay){

    iftAdjRel *A, *B;

    if(iftIs3DImage(img)) {
        A = iftSpheric(1.0);
        B = iftSpheric(1.8);
    }
    else {
        A = iftCircular(1.0);
        B = iftCircular(1.5);
    }

    iftImage* shape = iftObjectBorders(img, A, false, true);
    iftImage* eucl = iftEuclDistTrans(img, B, IFT_BOTH, NULL, NULL, NULL);

    float max = iftMaximumValue(eucl);
    for (int i = 0; i < eucl->n; ++i) {
        eucl->val[i] = max - eucl->val[i];

        if(decay>1.0)
            eucl->val[i] = MAX_SCORE*pow(eucl->val[i]/max, decay);
    }

    iftDestroyImage(&shape);

    return eucl;
}

char iftImageValidPoint(iftImage* img, iftPoint p)
{
    return ((p.x >= 0)&&(p.x < img->xsize)&&
            (p.y >= 0)&&(p.y < img->ysize)&&
            (p.z >= 0)&&(p.z < img->zsize));
}

int iftScoreValueAtPoint(iftImage* img, iftPoint p)
{
    if(iftImageValidPoint(img, p))
    if(img->dz>0)
        return iftImageValueAtPoint(img, p);
    else
        return iftImageValueAtPoint2D(img, p);
    else
        return 0;
}

void iftGetTransformationMatrix(void *prob, float *theta){
    iftShapeRegisterProb* shapeRegister = (iftShapeRegisterProb*) prob;
    iftMatrix* T = shapeRegister->_bufferT;

    float tx = shapeRegister->center.x;
    float ty = shapeRegister->center.y;
    float tz = shapeRegister->center.z;

    float rx = theta[R_X];
    float ry = theta[R_Y];
    float rz = theta[R_Z];

    rx = rx * IFT_PI / 180.0;
    ry = ry * IFT_PI / 180.0;
    rz = rz * IFT_PI / 180.0;

    float cosrx = cos(rx);
    float cosry = cos(ry);
    float cosrz = cos(rz);

    float sinrx = sin(rx);
    float sinry = sin(ry);
    float sinrz = sin(rz);

    T->val[iftGetMatrixIndex(T, 0, 0)] = cosry*cosrz;
    T->val[iftGetMatrixIndex(T, 1, 0)] = -sinrz;
    T->val[iftGetMatrixIndex(T, 2, 0)] = sinry*cosrz;
    T->val[iftGetMatrixIndex(T, 3, 0)] = -tx*cosry*cosrz + tx + ty*sinrz - tz*sinry*cosrz;
    T->val[iftGetMatrixIndex(T, 0, 1)] = sinrx*sinry + sinrz*cosrx*cosry;
    T->val[iftGetMatrixIndex(T, 1, 1)] = cosrx*cosrz;
    T->val[iftGetMatrixIndex(T, 2, 1)] = -sinrx*cosry + sinry*sinrz*cosrx;
    T->val[iftGetMatrixIndex(T, 3, 1)] = -tx*(sinrx*sinry + sinrz*cosrx*cosry) - ty*cosrx*cosrz + ty - tz*(-sinrx*cosry + sinry*sinrz*cosrx);
    T->val[iftGetMatrixIndex(T, 0, 2)] = sinrx*sinrz*cosry - sinry*cosrx;
    T->val[iftGetMatrixIndex(T, 1, 2)] = sinrx*cosrz;
    T->val[iftGetMatrixIndex(T, 2, 2)] = sinrx*sinry*sinrz + cosrx*cosry;
    T->val[iftGetMatrixIndex(T, 3, 2)] = -tx*(sinrx*sinrz*cosry - sinry*cosrx) - ty*sinrx*cosrz - tz*(sinrx*sinry*sinrz + cosrx*cosry) + tz;
    T->val[iftGetMatrixIndex(T, 0, 3)] = 0;
    T->val[iftGetMatrixIndex(T, 1, 3)] = 0;
    T->val[iftGetMatrixIndex(T, 2, 3)] = 0;
    T->val[iftGetMatrixIndex(T, 3, 3)] = 1;
}

//
//void _getTransformationMatrix(void *prob, float *theta)
//{
//
//    iftShapeRegisterProb* shapeRegister = (iftShapeRegisterProb*) prob;
//
//   //TODO: optimize this
//   iftMatrix* r1 = iftRotationMatrix(IFT_AXIS_X, theta[R_X]);
//   iftMatrix* r2 = iftRotationMatrix(IFT_AXIS_Z, theta[R_Z]);
//   iftMatrix* r3 = iftRotationMatrix(IFT_AXIS_Y, theta[R_Y]);
//
//   iftVector tv;
//
//   tv.x = theta[T_X];
//   tv.y = theta[T_Y];
//   tv.z = theta[T_Z];
//
////   tv.x = tv.y = tv.z = 0;
//
//   iftMatrix* tc = iftTranslationMatrix(shapeRegister->center);
//   iftMatrix* mtc = iftTranslationMatrix(iftVectorScalarProd(shapeRegister->center, -1.0f));
//   iftMatrix* t = iftTranslationMatrix(tv);
//   iftMatrix* Tf = iftMultMatricesChain(6, t, tc, r1, r2, r3, mtc);
//
//   for (int i = 0; i < shapeRegister->_bufferT->n; ++i) {
//      shapeRegister->_bufferT->val[i] = Tf->val[i];
//   }
//
//   iftDestroyMatrix(&r1);
//   iftDestroyMatrix(&r2);
//   iftDestroyMatrix(&t);
//   iftDestroyMatrix(&Tf);
//   iftDestroyMatrix(&tc);
//   iftDestroyMatrix(&mtc);
//}

//this function computes the transformation matrix for each comparison (debug purposes)
void iftComputeTransformationMatrix(void *prob, float *theta)
{

    iftShapeRegisterProb* shapeRegister = (iftShapeRegisterProb*) prob;

    iftMatrix* r1 = iftRotationMatrix(IFT_AXIS_X, theta[R_X]);
    iftMatrix* r2 = iftRotationMatrix(IFT_AXIS_Z, theta[R_Z]);
    iftMatrix* r3 = iftRotationMatrix(IFT_AXIS_Y, theta[R_Y]);

    iftVector tv;

    tv.x = theta[T_X];
    tv.y = theta[T_Y];
    tv.z = theta[T_Z];

    iftMatrix* tc = iftTranslationMatrix(shapeRegister->center);
    iftMatrix* mtc = iftTranslationMatrix((iftVector)iftVectorScalarProd(shapeRegister->center, -1.0f));
    iftMatrix* t = iftTranslationMatrix(tv);
    iftMatrix* Tf = iftMultMatricesChain(6, t, tc, r1, r2, r3, mtc);

    for (int i = 0; i < shapeRegister->_bufferT->n; ++i) {
        shapeRegister->_bufferT->val[i] = Tf->val[i];
    }

    iftDestroyMatrix(&r1);
    iftDestroyMatrix(&r2);
    iftDestroyMatrix(&t);
    iftDestroyMatrix(&Tf);
    iftDestroyMatrix(&tc);
    iftDestroyMatrix(&mtc);
}

float iftAffineTransformFitness(void *prob, float* x)
{
    int i;
    double scoreSum = 0.0f;
    iftShapeRegisterProb* shapeRegister = (iftShapeRegisterProb*)prob;
    iftPoint* points = shapeRegister->points;
    iftImage* score = shapeRegister->score;

    iftGetTransformationMatrix(prob, x);

    for (i = 0; i < shapeRegister->npoints; ++i) {
        iftPoint pt = iftTransformPoint(shapeRegister->_bufferT, points[i]);
        float f = iftScoreValueAtPoint(score, pt);
        scoreSum+=f/shapeRegister->npoints;
    }

    return (float)scoreSum;
}

void iftInitDeltaShapeRegister(iftMSPS* msps)
{
    //TODO: get variation as a parameter
    float angles[] = {1, 5, 7};//rotation
    //float steps[] = {1, 5, 17, 49, 121};//translation
    float angle_dev[] = {0.1, 0.5, 1.0};//angle deviation

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < msps->m; ++j) {
            int idx = iftGetMatrixIndex(msps->delta, i, j);
            msps->delta->val[idx] = angles[j];
            msps->sigma->val[idx] = angles[j]* angle_dev[j];
        }
    }
}

iftMatrix *iftRotationMatrixToAlignByPrincipalAxes(iftImage *bin)
{ 
  iftDataSet *Z  = iftObjectToDataSet(bin);
  iftSetStatus(Z,IFT_TRAIN);

  iftDataSet *Zc = iftCentralizeDataSet(Z);

  iftMatrix *A = iftDatasetCovarianceMatrix(Zc);
  iftMatrix *U,*Vt,*S;
  
  iftSingleValueDecomp(A,&U,&S,&Vt);

  iftDestroyMatrix(&U);
  iftDestroyMatrix(&S);
  iftDestroyMatrix(&A);
  iftDestroyDataSet(&Z);
  iftDestroyDataSet(&Zc);
  
  return(Vt);
}

iftMatrix *iftShapeRegister(iftPoint *orig, int norig, iftImage *score) {

    iftShapeRegisterProb* prob = (iftShapeRegisterProb*) iftAlloc(1, sizeof(iftShapeRegisterProb));

    prob->npoints = norig;
    prob->points = orig;
    prob->score = score;
    prob->center = iftGetPointsCenter(orig, norig);
    prob->_bufferT = iftCreateMatrix(4,4);

    iftMSPS* msps = iftCreateMSPS(3, 3, iftAffineTransformFitness, prob);

    //random MSPS
    //msps->iftPerturbation = iftMSPSLinearRandomPerturbation;

    iftInitDeltaShapeRegister(msps);
    iftMSPSMax(msps);
    iftGetTransformationMatrix(prob, msps->theta);

    iftFree(prob);

    return prob->_bufferT;
}

float iftRegistrationRMSE(iftImage *fixed_img, iftImage *moved_img)
{
    float mse = 0.;

    if (!iftIsDomainEqual(fixed_img,moved_img))
        iftError("Images do not have the same domain","iftMSE");


    for (int p = 0; p < fixed_img->n; p++){
        mse += (fixed_img->val[p] - moved_img->val[p])*(fixed_img->val[p] - moved_img->val[p]);
    }

    mse /= fixed_img->n;
    mse = sqrt(mse);

    return mse;
}

/******************************************************************/
