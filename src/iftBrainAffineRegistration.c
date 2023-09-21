
#include "iftBrainAffineRegistration.h" 

#include "ift/core/io/Stream.h"



//---------------------------------------------------------------------
//
// PRIVATE - Old-ift legacy code 
//
//---------------------------------------------------------------------

typedef struct _voxel {
  int x,y,z;
} Voxel;

typedef struct _scene {
  int *data;
  int xsize,ysize,zsize;
  float dx,dy,dz;
  int *tby, *tbz;
  int maxval, n;
  //nifti_image *nii_hdr;
} Scene;

#define GetVoxel(s,n) ((s)->data[(n)])
#define SceneLen(s) ((s)->xsize * (s)->ysize * (s)->zsize)
#define VoxelAddress(s,x,y,z) ((x)+(s)->tby[(y)]+(s)->tbz[(z)])

#define MSG1  "Cannot allocate memory space"
#define MSG2  "Cannot open file"
#define MSG3  "Invalid option"
#define MSG4  "Could not locate nifti header"
#define MSG5  "Nifti-1 data type not supported"


static void Error(const char *msg, const char *func){ /* It prints error message and exits
                                    the program. */
  fprintf(stderr,"Error:%s in %s\n",msg,func);
  exit(-1);
}


typedef unsigned char uchar;
typedef float real;
//typedef enum msp_boolean {msp_false,msp_true} msp_bool;

static uchar *AllocUCharArray(int n)
{
  uchar *v=NULL;
  v = (uchar *) calloc(n,sizeof(uchar));
  if (v == NULL)
    Error(MSG1,"AllocUCharArray");
  return(v);
}
   
static int *AllocIntArray(int n)
{
  int *v=NULL;
  v = (int *) calloc(n,sizeof(int));
  if (v == NULL)
    Error(MSG1,"AllocIntArray");
  return(v);
}

static ushort *AllocUShortArray(int n)
{
  ushort *v=NULL;
  v = (ushort *) calloc(n,sizeof(ushort));
  if (v == NULL)
    Error(MSG1,"AllocUShortArray");
  return(v);
}

static double *AllocDoubleArray(int n)
{
  double *v=NULL;
  v = (double *) calloc(n,sizeof(double));
  if (v == NULL)
    Error(MSG1,"AllocDoubleArray");
  return(v);
}

static float *AllocFloatArray(int n)
{
  float *v=NULL;
  v = (float *) calloc(n,sizeof(float));
  if (v == NULL)
    Error(MSG1,"AllocFloatArray");
  return(v);
}


static Scene  *CreateScene(int xsize,int ysize,int zsize)
{
  Scene *scn=NULL;
  int i,xysize;

  scn = (Scene *) calloc(1,sizeof(Scene));
  if (scn == NULL){
    Error(MSG1,"CreateScene");
  }

  scn->data    = AllocIntArray(xsize*ysize*zsize);
  scn->xsize   = xsize;
  scn->ysize   = ysize;
  scn->zsize   = zsize;
  scn->dx      = 1.0;
  scn->dy      = 1.0;
  scn->dz      = 1.0;
  scn->tby     = AllocIntArray(ysize);
  scn->tbz     = AllocIntArray(zsize);
  scn->maxval  = 0;
  scn->n       = xsize*ysize*zsize;
  //scn->nii_hdr = NULL;

  if (scn->data==NULL || scn->tbz==NULL || scn->tby==NULL) {
    Error(MSG1,"CreateScene");
  }

  scn->tby[0]=0;
  for (i=1; i < ysize; i++)
    scn->tby[i]=scn->tby[i-1] + xsize;

  scn->tbz[0]=0; xysize = xsize*ysize;
  for (i=1; i < zsize; i++)
    scn->tbz[i]=scn->tbz[i-1] + xysize;

  return(scn);
}


/* static Scene *CopyScene(Scene *scn){ */
/*   Scene *aux = CreateScene(scn->xsize,scn->ysize,scn->zsize); */
/*   aux->dx     = scn->dx; */
/*   aux->dy     = scn->dy; */
/*   aux->dz     = scn->dz; */
/*   aux->maxval = scn->maxval; */
/*   //aux->nii_hdr = NULL; */

/*   //if( scn->nii_hdr != NULL ) { */
/*   //  aux->nii_hdr = nifti_copy_nim_info( scn->nii_hdr ); */
/*   //} */
/*   aux->n=scn->xsize*scn->ysize*scn->zsize; */
/*   memcpy(aux->data,scn->data,sizeof(int) * aux->n); */
/*   return aux; */
/* } */

static Voxel   Transform_Voxel(float M[4][4],Voxel v){
  Voxel nv;

  nv.x = (int) (M[0][0]*v.x + M[0][1]*v.y + M[0][2]*v.z + M[0][3]);
  nv.y = (int) (M[1][0]*v.x + M[1][1]*v.y + M[1][2]*v.z + M[1][3]);
  nv.z = (int) (M[2][0]*v.x + M[2][1]*v.y + M[2][2]*v.z + M[2][3]);

  return(nv);

}

static void     DestroyScene(Scene **scn)
{
  Scene *aux;

  aux = *scn;
  if(aux != NULL){
    if (aux->data != NULL)  free(aux->data);
    if (aux->tby  != NULL)  free(aux->tby);
    if (aux->tbz  != NULL)  free(aux->tbz);
    //if (aux->nii_hdr != NULL) nifti_image_free( aux->nii_hdr );
    free(aux);
    *scn = NULL;
  }
}

static int ValidVoxel(Scene *scn, int vx, int vy, int vz)
{
  if ((vx >= 0)&&(vx < scn->xsize)&&
      (vy >= 0)&&(vy < scn->ysize)&&
      (vz >= 0)&&(vz < scn->zsize))
    return(1);
  else
    return(0);
}

static int MaximumValue3(Scene *scn)
{
  unsigned int i, n, r;
  int Imax;

  Imax = 0;
  n = scn->xsize*scn->ysize*scn->zsize - 1;
  r = n%4;
  n -= r;
  for (i=0; i < n; i+=4) {
    if (scn->data[i] > Imax)
      Imax = scn->data[i];
    if (scn->data[i+1] > Imax)
      Imax = scn->data[i+1];
    if (scn->data[i+2] > Imax)
      Imax = scn->data[i+2];
    if (scn->data[i+3] > Imax)
      Imax = scn->data[i+3];
  }
  while (r != 0) {
    if (scn->data[i+r-1] > Imax)
      Imax = scn->data[i+r-1];
    --r;
  }

  scn->maxval = Imax;

  return(Imax);
}

static void WriteScene(Scene *scn, char *filename)
{
  FILE *fp=NULL;
  int Imax;
  int i,n;
  uchar  *data8 =NULL;
  ushort *data16=NULL;

  // Checking file type
  int len = strlen(filename);
  if ( (len>=4) && ((strcasecmp(filename + len - 4, ".hdr")==0) || (strcasecmp(filename + len - 4, ".img")==0))) {
    //WriteScene_Analyze(scn, filename);
    //WriteScene_Nifti1(scn, filename);
    Error(MSG2,"WriteScene: Invalid file name or extension.");
    return;
  }
  if ( (len>=8) && (strcasecmp(filename + len - 8, ".scn.bz2")==0)) {
    //WriteCompressedScene(scn, filename);
    Error(MSG2,"WriteScene: Invalid file name or extension.");
    return;
  }
  if ( (len>=7) && (strcasecmp(filename + len - 7, ".nii.gz")==0)) {
    //WriteScene_Nifti1( scn, filename );
    Error(MSG2,"WriteScene: Invalid file name or extension.");
    return;
  }
  if ( (len>=4) && (strcasecmp(filename + len - 4, ".nii")==0)) {
    //WriteScene_Nifti1( scn, filename );
    return;
  }
  if ( (len<=4) || (strcasecmp(filename + len - 4, ".scn")!=0)) {
    Error(MSG2,"WriteScene: Invalid file name or extension.");
  }

  // Writing the scn file
  fp = fopen(filename,"wb");
  if (fp == NULL)
    Error(MSG2,"WriteScene");

  fprintf(fp,"SCN\n");
  fprintf(fp,"%d %d %d\n",scn->xsize,scn->ysize,scn->zsize);
  fprintf(fp,"%f %f %f\n",scn->dx,scn->dy,scn->dz);

  Imax = MaximumValue3(scn);

  n = scn->xsize*scn->ysize*scn->zsize;
  if (Imax < 256) {
    fprintf(fp,"%d\n",8);
    data8 = AllocUCharArray(n);
    for (i=0; i < n; i++)
      data8[i] = (uchar) scn->data[i];
    fwrite(data8,sizeof(uchar),n,fp);
    free(data8);
  } else if (Imax < 65536) {
    fprintf(fp,"%d\n",16);
    data16 = AllocUShortArray(n);
    for (i=0; i < n; i++)
      data16[i] = (ushort) scn->data[i];
    fwrite(data16,sizeof(ushort),n,fp);
    free(data16);
  } else {
    fprintf(fp,"%d\n",32);
    fwrite(scn->data,sizeof(int),n,fp);
  }
  fclose(fp);
}



typedef struct _curve { /* Curve */
  double *X;
  double *Y;
  int n;
} Curve;

static Curve *CreateCurve(int n)
{
  Curve *curve=NULL;

  curve = (Curve *) calloc(1,sizeof(Curve));
  if (curve != NULL) {
    curve->X = AllocDoubleArray(n);
    curve->Y = AllocDoubleArray(n);
    curve->n = n;
  } else {
    Error(MSG1,"CreateCurve");
  }
  return(curve);
}

static void DestroyCurve(Curve **curve)
{
  Curve *aux;

  aux = *curve;
  if (aux != NULL){
    if (aux->X != NULL) free(aux->X);
    if (aux->Y != NULL) free(aux->Y);
    free(aux);
    *curve = NULL;
  }
}

static Curve *Histogram3(Scene *scn)
{
  int i, n, nbins;
  Curve *hist = NULL;

  nbins = MaximumValue3(scn)+1;
  hist  = CreateCurve(nbins);
  n = scn->xsize * scn->ysize * scn->zsize;
  for (i = 0; i < n; i++)
    hist->Y[scn->data[i]]++;
  for (i = 0; i < nbins; i++)
    hist->X[i] = i;

  return (hist);
}


static Curve *NormHistogram3(Scene *scn)
{
  int i, sum;
  Curve *nhist;

  nhist = Histogram3(scn);
  sum = scn->xsize * scn->ysize * scn->zsize;
  for (i = 0; i < nhist->n; i++){
    nhist->Y[i] /= sum;
    nhist->X[i]=i;
  }

  return (nhist);
}

static Curve *NormAccHistogram3(Scene *scn)
{
  int i;
  Curve *ahist;

  ahist = NormHistogram3(scn);
  for (i = 1; i < ahist->n; i++){
    ahist->Y[i] = ahist->Y[i-1] + ahist->Y[i];
    ahist->X[i] = i;
  }
  return (ahist);
}

static Scene *Add3(Scene *scn1, int value)
{
  Scene *sum;
  int p, n;

  sum = CreateScene(scn1->xsize, scn1->ysize, scn1->zsize);
  n = scn1->xsize * scn1->ysize * scn1->zsize;
  for (p = 0; p < n; p++)
    sum->data[p] = scn1->data[p] + value;

  return (sum);
}

typedef struct _vector{
  float x;
  float y;
  float z;
} Vector, Point, Vertex;

/* static void VectorNormalize(Vector *v) { */

/*   float norm; */

/*   norm= sqrt(v->x*v->x + v->y*v->y + v->z*v->z); */

/*   if (norm != 0.0) { */

/*     v->x = v->x/norm; */
/*     v->y = v->y/norm; */
/*     v->z = v->z/norm; */
/*   } */
/* } */


typedef struct _adjrel3 {
  int *dx;
  int *dy;
  int *dz;
  int n;
} AdjRel3;
static AdjRel3 *CreateAdjRel3(int n)
{
  AdjRel3 *A=NULL;

  A = (AdjRel3 *) calloc(1,sizeof(AdjRel3));
  if (A != NULL){
    A->dx = AllocIntArray(n);
    A->dy = AllocIntArray(n);
    A->dz = AllocIntArray(n);
    A->n  = n;
  } else {
    Error(MSG1,"CreateAdjRel3");
  }

  return(A);
}

static void DestroyAdjRel3(AdjRel3 **A)
{
  AdjRel3 *aux;

  aux = *A;
  if (aux != NULL){
    if (aux->dx != NULL) free(aux->dx);
    if (aux->dy != NULL) free(aux->dy);
    if (aux->dz != NULL) free(aux->dz);
    free(aux);
    *A = NULL;
  }
}


static AdjRel3 *Spheric(float r)
{
  AdjRel3 *A=NULL;
  int i,n,r0,r2,dx,dy,dz,i0=0;

  n=0;
  r0 = (int)r;
  r2  = (int)(r*r + 0.5);
  for(dz=-r0;dz<=r0;dz++)
    for(dy=-r0;dy<=r0;dy++)
      for(dx=-r0;dx<=r0;dx++)
      if(((dx*dx)+(dy*dy)+(dz*dz)) <= r2)
        n++;

  A = CreateAdjRel3(n);
  i=0;
  for(dz=-r0;dz<=r0;dz++)
    for(dy=-r0;dy<=r0;dy++)
      for(dx=-r0;dx<=r0;dx++)
        if(((dx*dx)+(dy*dy)+(dz*dz)) <= r2){
          A->dx[i]=dx;
          A->dy[i]=dy;
          A->dz[i]=dz;
          if ((dx==0)&&(dy==0)&&(dz==0))
            i0 = i;
          i++;
        }

  /* shift to right and place central voxel at first */

  for (i=i0; i > 0; i--) {
    dx = A->dx[i];
    dy = A->dy[i];
    dz = A->dz[i];
    A->dx[i] = A->dx[i-1];
    A->dy[i] = A->dy[i-1];
    A->dz[i] = A->dz[i-1];
    A->dx[i-1] = dx;
    A->dy[i-1] = dy;
    A->dz[i-1] = dz;
  }

  return(A);
}

typedef struct _realmatrix {
  real **val;
  int ncols,nrows;
} RealMatrix;

static RealMatrix *CreateRealMatrix(int ncols,int nrows){
  RealMatrix *mat=NULL;
  real *aux;
  int i;

  mat = (RealMatrix *) calloc(1,sizeof(RealMatrix));
  if(mat == NULL)
    Error(MSG1,"CreateRealMatrix");

  aux = (real *)calloc(nrows*ncols, sizeof(real));
  mat->val = (real **) calloc(nrows, sizeof(real *));
  if(mat->val == NULL || aux == NULL)
    Error(MSG1,"CreateRealMatrix");

  mat->val[0] = aux;
  for(i=1; i<nrows; i++)
    mat->val[i] = mat->val[i-1] + ncols;

  mat->ncols = ncols;
  mat->nrows = nrows;

  return(mat);
}


static void        DestroyRealMatrix(RealMatrix **mat){
  RealMatrix *aux;

  aux = *mat;
  if(aux != NULL){
    if(aux->val != NULL){
      if(*(aux->val) != NULL)
        free(*(aux->val));
      free(aux->val);
    }
    free(aux);
    *mat = NULL;
  }
}

static RealMatrix *CloneRealMatrix(RealMatrix *mat){
  RealMatrix *matc;

  matc = CreateRealMatrix(mat->ncols, mat->nrows);
  memcpy(matc->val[0], mat->val[0],
         mat->ncols*mat->nrows*sizeof(real));

  return(matc);
}

RealMatrix *InvertRealMatrix(RealMatrix *A){
  RealMatrix *B=NULL;
  RealMatrix *X=NULL;
  int i,j,k;
  double m;

  if(A->ncols!=A->nrows)
    Error("Matrix dimension error","InvertRealMatrix");

  X = CreateRealMatrix(A->ncols, A->nrows);
  B = CloneRealMatrix(A);

  for(i=0; i<A->nrows; i++)
    X->val[i][i] = 1.0;

  for(k=0; k<A->nrows; k++){
    m = B->val[k][k];
    //if(m < 0.0000000000001)
      // Error("Singular matrix","InvertRealMatrix");   THIS IS WRONG

    B->val[k][k] = 1.0;
    for(j=0; j<A->ncols; j++){
      if(j!=k)
        B->val[k][j] = B->val[k][j]/m;
      X->val[k][j] = X->val[k][j]/m;
    }

    for(i=0; i<A->nrows; i++){
      if(i!=k){
        m = B->val[i][k]/B->val[k][k];

        B->val[i][k] = 0.0;
        for(j=0; j<A->ncols; j++){
          if(j!=k)
            B->val[i][j] = B->val[i][j] - m*B->val[k][j];
          X->val[i][j] = X->val[i][j] - m*X->val[k][j];
        }
      }
    }
  }

  DestroyRealMatrix(&B);
  return X;
}

static RealMatrix *MultRealMatrix(RealMatrix *A,
                           RealMatrix *B){
  RealMatrix *M = NULL;
  int i,j,k;

  if(A->ncols!=B->nrows)
    Error("Matrix dimension error","MultRealMatrix");

  M = CreateRealMatrix(B->ncols, A->nrows);
  for(i=0; i<M->nrows; i++){
    for(j=0; j<M->ncols; j++){
      M->val[i][j] = 0.0;
      for (k=0; k<A->ncols; k++)
        M->val[i][j] += A->val[i][k]*B->val[k][j];
    }
  }
  return(M);
}

static RealMatrix* RotationMatrix3(int axis, // options: 0 (x) / 1 (y) / 2 (z)
                            double th)
{
  RealMatrix *m;
  m = CreateRealMatrix(4,4);
  if (axis==0) {
    m->val[0][0] = 1.0;    m->val[0][1] = 0.0;    m->val[0][2] = 0.0;    m->val[0][3] = 0.0;
    m->val[1][0] = 0.0;    m->val[1][1] = cos(th);    m->val[1][2] = -sin(th);    m->val[1][3] = 0.0;
    m->val[2][0] = 0.0;    m->val[2][1] = sin(th);    m->val[2][2] = cos(th);    m->val[2][3] = 0.0;
    m->val[3][0] = 0.0;    m->val[3][1] = 0.0;    m->val[3][2] = 0.0;    m->val[3][3] = 1.0;
  }
  if (axis==1) {
    m->val[0][0] = cos(th);    m->val[0][1] = 0.0;    m->val[0][2] = sin(th);    m->val[0][3] = 0.0;
    m->val[1][0] = 0.0;    m->val[1][1] = 1;    m->val[1][2] = 0.0;    m->val[1][3] = 0.0;
    m->val[2][0] = -sin(th);    m->val[2][1] = 0;    m->val[2][2] = cos(th);    m->val[2][3] = 0.0;
    m->val[3][0] = 0.0;    m->val[3][1] = 0.0;    m->val[3][2] = 0.0;    m->val[3][3] = 1.0;

  }
  if (axis==2) {
    m->val[0][0] = cos(th);    m->val[0][1] = -sin(th);    m->val[0][2] = 0.0;    m->val[0][3] = 0.0;
    m->val[1][0] = sin(th);    m->val[1][1] = cos(th);    m->val[1][2] = 0.0;    m->val[1][3] = 0.0;
    m->val[2][0] = 0.0;    m->val[2][1] = 0.0;    m->val[2][2] = 1.0;    m->val[2][3] = 0.0;
    m->val[3][0] = 0.0;    m->val[3][1] = 0.0;    m->val[3][2] = 0.0;    m->val[3][3] = 1.0;

  }
  return m;
}

static RealMatrix* TranslationMatrix3(float dx, float dy, float dz)
{
  RealMatrix *m;
  m = CreateRealMatrix(4,4);
  m->val[0][0] = 1.0;    m->val[0][1] = 0.0;    m->val[0][2] = 0.0;    m->val[0][3] = dx;
  m->val[1][0] = 0.0;    m->val[1][1] = 1.0;    m->val[1][2] = 0.0;    m->val[1][3] = dy;
  m->val[2][0] = 0.0;    m->val[2][1] = 0.0;    m->val[2][2] = 1.0;    m->val[2][3] = dz;
  m->val[3][0] = 0.0;    m->val[3][1] = 0.0;    m->val[3][2] = 0.0;    m->val[3][3] = 1.0;
  return m;
}

static RealMatrix* ScaleMatrix3(float Sx, float Sy, float Sz)
{
  RealMatrix *m;
  m = CreateRealMatrix(4,4);
  m->val[0][0] = Sx;    m->val[0][1] = 0.0;    m->val[0][2] = 0.0;    m->val[0][3] = 0;
  m->val[1][0] = 0.0;    m->val[1][1] = Sy;    m->val[1][2] = 0.0;    m->val[1][3] = 0;
  m->val[2][0] = 0.0;    m->val[2][1] = 0.0;    m->val[2][2] = Sz;    m->val[2][3] = 0;
  m->val[3][0] = 0.0;    m->val[3][1] = 0.0;    m->val[3][2] = 0.0;    m->val[3][3] = 1.0;
  return m;
}

static RealMatrix* ShearMatrix3(float SHxy,float SHxz,float SHyx,float SHyz,float SHzx,float SHzy)
{
  RealMatrix *m;
  m = CreateRealMatrix(4,4);
  m->val[0][0] = 1.0;    m->val[0][1] = SHxy;    m->val[0][2] = SHxz;    m->val[0][3] = 0;
  m->val[1][0] = SHyx;    m->val[1][1] = 1.0;    m->val[1][2] = SHyz;    m->val[1][3] = 0;
  m->val[2][0] = SHzx;    m->val[2][1] = SHzy;    m->val[2][2] = 1.0;    m->val[2][3] = 0;
  m->val[3][0] = 0.0;    m->val[3][1] = 0.0;    m->val[3][2] = 0.0;    m->val[3][3] = 1.0;
  return m;
}


/* static RealMatrix* TransformVoxel(RealMatrix *m, Voxel v) */
/* { */
/*   RealMatrix *vm,*res; */
/*   vm = CreateRealMatrix(1,4); */
/*   vm->val[0][0]=v.x; */
/*   vm->val[1][0]=v.y; */
/*   vm->val[2][0]=v.z; */
/*   vm->val[3][0]=1.0; */
/*   res=MultRealMatrix(m,vm); */
/*   DestroyRealMatrix(&vm); */
/*   return res; */
/* } */




#define INTERIOR    0
#define EXTERIOR    1
#define BOTH        2
#define WHITE       0
#define GRAY        1
#define BLACK       2
#define NIL        -1
#define INCREASING  1
#define DECREASING  0
#ifndef MAX
#define MAX(x,y) (((x) > (y))?(x):(y))
#endif

#ifndef MIN
#define MIN(x,y) (((x) < (y))?(x):(y))
#endif

#define ROUND(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))

#define SIGN(x) ((x >= 0)?1:-1)

#define MINVALUE   0 /* define queue to remove node with minimum value */
#define MAXVALUE   1 /* define queue to remove node with maximum value */
#define FIFOBREAK 0  /* define queue to solve ambiguity by FIFO */
#define LIFOBREAK 1  /* define queue to solve ambiguity by LIFO */
#define QSIZE     32768

#define SetTieBreak(a,b) a->C.tiebreak=b 
#define SetRemovalPolicy(a,b) a->C.removal_policy=b 

typedef struct _gqnode { 
  int  next;  /* next node */
  int  prev;  /* prev node */
  char color; /* WHITE=0, GRAY=1, BLACK=2 */ 
} GQNode;

typedef struct _gdoublylinkedlists {
  GQNode *elem;  /* all possible doubly-linked lists of the circular queue */
  int nelems;  /* total number of elements */
  int *value;   /* the value of the nodes in the graph */
} GDoublyLinkedLists; 

typedef struct _gcircularqueue { 
  int  *first;   /* list of the first elements of each doubly-linked list */
  int  *last;    /* list of the last  elements of each doubly-linked list  */
  int  nbuckets; /* number of buckets in the circular queue */
  int  minvalue;  /* minimum value of a node in queue */
  int  maxvalue;  /* maximum value of a node in queue */
  char tiebreak; /* 1 is LIFO, 0 is FIFO (default) */
  char removal_policy; /* 0 is MINVALUE and 1 is MAXVALUE */
} GCircularQueue;

typedef struct _gqueue { /* Priority queue by Dial implemented as
                           proposed by A. Falcao */
  GCircularQueue C;
  GDoublyLinkedLists L;
} GQueue;


static void ResetGQueue(GQueue *Q)
{
    int i;

    Q->C.minvalue = INT_MAX;
    Q->C.maxvalue = INT_MIN;
    SetTieBreak(Q,FIFOBREAK);
    SetRemovalPolicy(Q,MINVALUE);
    for (i=0; i < Q->C.nbuckets+1; i++)
        Q->C.first[i]=Q->C.last[i]=NIL;

    for (i=0; i < Q->L.nelems; i++)
    {
        Q->L.elem[i].next =  Q->L.elem[i].prev = NIL;
        Q->L.elem[i].color = WHITE;
    }

}

static GQueue *CreateGQueue(int nbuckets, int nelems, int *value)
{
    GQueue *Q=NULL;

    Q = (GQueue *) malloc(1*sizeof(GQueue));

    if (Q != NULL)
    {
        Q->C.first = (int *)malloc((nbuckets+1) * sizeof(int));
        Q->C.last  = (int *)malloc((nbuckets+1) * sizeof(int));
        Q->C.nbuckets = nbuckets;
        if ( (Q->C.first != NULL) && (Q->C.last != NULL) )
        {
            Q->L.elem = (GQNode *)malloc(nelems*sizeof(GQNode));
            Q->L.nelems = nelems;
            Q->L.value   = value;
            if (Q->L.elem != NULL)
            {
                ResetGQueue(Q);
            }
            else
                Error(MSG1,"CreateGQueue");
        }
        else
            Error(MSG1,"CreateGQueue");
    }
    else
        Error(MSG1,"CreateGQueue");

    return(Q);
}


static void DestroyGQueue(GQueue **Q)
{
    GQueue *aux;

    aux = *Q;
    if (aux != NULL)
    {
        if (aux->C.first != NULL) free(aux->C.first);
        if (aux->C.last  != NULL) free(aux->C.last);
        if (aux->L.elem  != NULL) free(aux->L.elem);
        free(aux);
        *Q = NULL;
    }
}

static GQueue *GrowGQueue(GQueue **Q, int nbuckets)
{
    GQueue *Q1=CreateGQueue(nbuckets,(*Q)->L.nelems,(*Q)->L.value);
    int i,bucket;

    Q1->C.minvalue  = (*Q)->C.minvalue;
    Q1->C.maxvalue  = (*Q)->C.maxvalue;
    Q1->C.tiebreak = (*Q)->C.tiebreak;
    Q1->C.removal_policy = (*Q)->C.removal_policy;
    for (i=0; i<(*Q)->C.nbuckets; i++)
        if ((*Q)->C.first[i]!=NIL)
        {
            bucket = (*Q)->L.value[(*Q)->C.first[i]]%Q1->C.nbuckets;
            Q1->C.first[bucket] = (*Q)->C.first[i];
            Q1->C.last[bucket]  = (*Q)->C.last[i];
        }
    if ((*Q)->C.first[(*Q)->C.nbuckets]!=NIL)
    {
        bucket = Q1->C.nbuckets;
        Q1->C.first[bucket] = (*Q)->C.first[(*Q)->C.nbuckets];
        Q1->C.last[bucket]  = (*Q)->C.last[(*Q)->C.nbuckets];
    }

    for (i=0; i < (*Q)->L.nelems; i++)
        Q1->L.elem[i]  = (*Q)->L.elem[i];

    DestroyGQueue(Q);
    return(Q1);
}


static void InsertGQueue(GQueue **Q, int elem)
{
    int bucket,minvalue=(*Q)->C.minvalue,maxvalue=(*Q)->C.maxvalue;

    if (((*Q)->L.value[elem]==INT_MAX)||((*Q)->L.value[elem]==INT_MIN))
        bucket=(*Q)->C.nbuckets;
    else
    {
        if ((*Q)->L.value[elem] < minvalue)
            minvalue = (*Q)->L.value[elem];
        if ((*Q)->L.value[elem] > maxvalue)
            maxvalue = (*Q)->L.value[elem];
        if ((maxvalue-minvalue) > ((*Q)->C.nbuckets-1))
        {
            (*Q) = GrowGQueue(Q,2*(maxvalue-minvalue)+1);
            fprintf(stdout,"Warning: Doubling queue size\n");
        }
        if ((*Q)->C.removal_policy==MINVALUE)
        {
            bucket=(*Q)->L.value[elem]%(*Q)->C.nbuckets;
        }
        else
        {
            bucket=(*Q)->C.nbuckets-1-((*Q)->L.value[elem]%(*Q)->C.nbuckets);
        }
        (*Q)->C.minvalue = minvalue;
        (*Q)->C.maxvalue = maxvalue;
    }
    if ((*Q)->C.first[bucket] == NIL)
    {
        (*Q)->C.first[bucket]   = elem;
        (*Q)->L.elem[elem].prev = NIL;
    }
    else
    {
        (*Q)->L.elem[(*Q)->C.last[bucket]].next = elem;
        (*Q)->L.elem[elem].prev = (*Q)->C.last[bucket];
    }

    (*Q)->C.last[bucket]     = elem;
    (*Q)->L.elem[elem].next  = NIL;
    (*Q)->L.elem[elem].color = GRAY;
}

static int RemoveGQueue(GQueue *Q)
{
    int elem=NIL, next, prev;
    int last, current;

    if (Q->C.removal_policy==MINVALUE)
        current=Q->C.minvalue%Q->C.nbuckets;
    else
        current=Q->C.nbuckets-1-(Q->C.maxvalue%Q->C.nbuckets);

    /** moves to next element **/

    if (Q->C.first[current] == NIL)
    {
        last = current;

        current = (current + 1) % (Q->C.nbuckets);

        while ((Q->C.first[current] == NIL) && (current != last))
        {
            current = (current + 1) % (Q->C.nbuckets);
        }

        if (Q->C.first[current] != NIL)
        {
            if (Q->C.removal_policy==MINVALUE)
                Q->C.minvalue = Q->L.value[Q->C.first[current]];
            else
                Q->C.maxvalue = Q->L.value[Q->C.first[current]];
        }
        else
        {
            if (Q->C.first[Q->C.nbuckets] != NIL)
            {
                current = Q->C.nbuckets;
                if (Q->C.removal_policy==MINVALUE)
                    Q->C.minvalue = Q->L.value[Q->C.first[current]];
                else
                    Q->C.maxvalue = Q->L.value[Q->C.first[current]];
            }
            else
            {
                Error("GQueue is empty\n","RemoveGQueue");
            }
        }
    }

    if (Q->C.tiebreak == LIFOBREAK)
    {
        elem = Q->C.last[current];
        prev = Q->L.elem[elem].prev;
        if (prev == NIL)           /* there was a single element in the list */
        {
            Q->C.last[current] = Q->C.first[current]  = NIL;
        }
        else
        {
            Q->C.last[current]   = prev;
            Q->L.elem[prev].next = NIL;
        }
    }
    else   /* Assume FIFO policy for breaking ties */
    {
        elem = Q->C.first[current];
        next = Q->L.elem[elem].next;
        if (next == NIL)           /* there was a single element in the list */
        {
            Q->C.first[current] = Q->C.last[current]  = NIL;
        }
        else
        {
            Q->C.first[current] = next;
            Q->L.elem[next].prev = NIL;
        }
    }

    Q->L.elem[elem].color = BLACK;

    return elem;
}

static void RemoveGQueueElem(GQueue *Q, int elem)
{
    int prev,next,bucket;

    if ((Q->L.value[elem] == INT_MAX)||(Q->L.value[elem] == INT_MIN))
        bucket = Q->C.nbuckets;
    else
    {
        if (Q->C.removal_policy == MINVALUE)
            bucket = Q->L.value[elem]%Q->C.nbuckets;
        else
            bucket = Q->C.nbuckets-1-(Q->L.value[elem]%Q->C.nbuckets);
    }

    prev = Q->L.elem[elem].prev;
    next = Q->L.elem[elem].next;

    /* if elem is the first element */
    if (Q->C.first[bucket] == elem)
    {
        Q->C.first[bucket] = next;
        if (next == NIL) /* elem is also the last one */
            Q->C.last[bucket] = NIL;
        else
            Q->L.elem[next].prev = NIL;
    }
    else    /* elem is in the middle or it is the last */
    {
        Q->L.elem[prev].next = next;
        if (next == NIL) /* if it is the last */
            Q->C.last[bucket] = prev;
        else
            Q->L.elem[next].prev = prev;
    }

    Q->L.elem[elem].color = BLACK;

}

static void UpdateGQueue(GQueue **Q, int elem, int newvalue)
{
    RemoveGQueueElem(*Q, elem);
    (*Q)->L.value[elem] = newvalue;
    InsertGQueue(Q, elem);
}

static int EmptyGQueue(GQueue *Q)
{
    int last,current;

    if (Q->C.removal_policy == MINVALUE)
        current=Q->C.minvalue%Q->C.nbuckets;
    else
        current=Q->C.nbuckets - 1 - (Q->C.maxvalue%Q->C.nbuckets);

    if (Q->C.first[current] != NIL)
        return 0;

    last = current;

    current = (current + 1) % (Q->C.nbuckets);

    while ((Q->C.first[current] == NIL) && (current != last))
    {
        current = (current + 1) % (Q->C.nbuckets);
    }

    if (Q->C.first[current] == NIL)
    {
        if (Q->C.first[Q->C.nbuckets] == NIL)
        {
	    ResetGQueue(Q);
            return(1);
        }
    }

    return (0);
}




static Scene *WaterGray3(Scene *scn, Scene *marker, AdjRel3 *A)
{
  Scene  *cost=NULL,*label=NULL;
  GQueue *Q=NULL;
  int i,p,q,tmp,n,r=1,sz;
  Voxel u,v;
  
  n     = scn->xsize*scn->ysize*scn->zsize;
  sz    = (scn->xsize * scn->ysize);
  cost  = CreateScene(scn->xsize,scn->ysize,scn->zsize);
  label = CreateScene(scn->xsize,scn->ysize,scn->zsize);
  Q     = CreateGQueue(MaximumValue3(marker)+2,n,cost->data);
  for (p=0; p < n; p++) {
    cost->data[p]=marker->data[p]+1;
    InsertGQueue(&Q,p);
  }

  while(!EmptyGQueue(Q)) {
    p=RemoveGQueue(Q);
    if (label->data[p]==0) {
      cost->data[p]=scn->data[p];
      label->data[p]=r;
      r++;
    }
    u.z =  p/sz;
    u.x = (p%sz)%scn->xsize;
    u.y = (p%sz)/scn->xsize;
    for (i=1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      v.z = u.z + A->dz[i];
      if (ValidVoxel(scn,v.x,v.y,v.z)){	
	q = v.x + scn->tby[v.y] + scn->tbz[v.z];
	if (cost->data[q] > cost->data[p]){
	  tmp = MAX(cost->data[p],scn->data[q]);
	  if (tmp < cost->data[q]){
	    UpdateGQueue(&Q,q,tmp);
	    label->data[q] = label->data[p];
	  }
	}
      }
    }
  }
  
  DestroyGQueue(&Q);
  DestroyScene(&cost);

  return(label);
}



typedef struct _set {
  int elem;
  struct _set *next;
} Set;



static void InsertSet(Set **S, int elem)
{
  Set *p=NULL;

  p = (Set *) calloc(1,sizeof(Set));
  if (p == NULL) Error(MSG1,"InsertSet");
  if (*S == NULL){
    p->elem  = elem;
    p->next  = NULL;
  }else{
    p->elem  = elem;
    p->next  = *S;
  }
  *S = p;
}

static void DestroySet(Set **S)
{
  Set *p;
  while(*S != NULL){
    p = *S;
    *S = p->next;
    free(p);
  }
}




static Scene *GetBorder3(Scene *scn, AdjRel3 *A)
{
  Scene *hscn=NULL;
  int p,q,i;
  Voxel u,v;

  hscn = CreateScene(scn->xsize,scn->ysize,scn->zsize);
  for (u.z=0; u.z < hscn->zsize; u.z++){
    for (u.y=0; u.y < hscn->ysize; u.y++){
      for (u.x=0; u.x < hscn->xsize; u.x++){
	p = u.x + hscn->tby[u.y] + hscn->tbz[u.z];
	if (scn->data[p] != 0) {
	  for (i=1; i < A->n; i++){
	    v.x = u.x + A->dx[i];
	    v.y = u.y + A->dy[i];
	    v.z = u.z + A->dz[i];
	    if (ValidVoxel(hscn,v.x,v.y,v.z)){
	      q = v.x + hscn->tby[v.y] + hscn->tbz[v.z];
	      if (scn->data[p] != scn->data[q]){
		hscn->data[p] = scn->data[p];
	        break;
	      }
	    } else {
	      hscn->data[p] = scn->data[p];
	      break;
	    }
	  }
	}
      }
    }
  }
  return(hscn);
}




/* Type for a fitness function. That is, an optimization problem. */
typedef double (*problem) ( double *x, void *context);

struct optContext
{
  problem f;                 /* The fitness function to be optimized. */
  void *fContext;            /* The problem's context. */
  size_t fDim;               /* Dimensionality of problem. */
  double *theta;             /* Theta to be initially considered and to be return
                                as the best one in the end of the optimization. */
  double * lowerInit;   /* Lower initialization boundary. */
  double * upperInit;   /* Upper initialization boundary. */
  double * lowerBound;  /* Lower search-space boundary. */
  double * upperBound;  /* Upper search-space boundary. */
  size_t numIterations;      /* Number of fitness-evaluations. */
  char verbose;              /* Whether it prints information from the optmiation or not. */
};


/***************************************************************
 * Multiscale Parameter Search (MSPS)
*****************/

/*Whether the method is for minimization or maximization. */
/* > is maximizatio, < is minimization */
#define MinMaxOp >

static void InitVector(double *v,  double value, size_t n);
static void CopyVector(double *dest,  double *src, size_t n);
static void SumVector(double *dest,  double *a,  double *b, size_t n);
static char BoundTest(double *x, size_t i,  double *lower,  double *upper);
static void PrintFunctionEval( char *label,  double *theta, size_t n, double v );

static double MSPS( double *param, void* context)
{
  /* Cast void-ptr context to correct struct-type. */
  struct optContext *c = (struct optContext*) context;

  /* Clone context to local variables for easier reference. */
  problem f = c->f;                         /* Fitness function. */
  void *fContext = c->fContext;             /* Context for fitness function. */
  size_t n = c->fDim;                       /* Dimensionality of problem. */
  double * lowerInit  = c->lowerInit;  /* Lower initialization boundary. */
  double * upperInit  = c->upperInit;  /* Upper initialization boundary. */
  double * lowerBound = c->lowerBound; /* Lower search-space boundary. */
  double * upperBound = c->upperBound; /* Upper search-space boundary. */
  size_t numIterations = c->numIterations;  /* Number of iterations to perform. */
  char verbose = c->verbose;

  /* Retrieve the parameters of the method. */
  double scales = param[0];
  double degree = param[1];
  double gamma  = param[2];

  /* Initialize search-range factor and decrease-factor. */
  double r = 1;                    /* Search-range. */
  double q = 1 / pow(2.0, gamma);   //Decrease-factor. 

  /* Iteration variables. */
  size_t i, j, k=0, bigIter=1, bestP=0;

  /* Rounded scales*/
  size_t rScales = (size_t) (scales + 0.5);

  double scaleNorm;

  /* Allocate and set up the delta matrix. */
  double **delta = (double**) malloc(sizeof(double*)*n);
  for ( i=0; i < n; i++ )
  {
    delta[i] = AllocDoubleArray(rScales);
    scaleNorm = pow( rScales, degree ) / (( upperInit[i] - lowerInit[i] ) / 2 );
    for ( j=1; j <= rScales; j++ )
    {
      delta[i][j-1] = pow( j, degree ) / scaleNorm;
    }
  }

  /* Allocate thetas and deltas. */
  double *theta      = AllocDoubleArray(n);
  double *thetaTmp   = AllocDoubleArray(n);
  double *thetaPrime = AllocDoubleArray(n);

  double *deltaPrimeS = AllocDoubleArray(n);
  double *deltaPrimeP = AllocDoubleArray(n);

  double *VP          = AllocDoubleArray(n);

  /* Temporary store. */
  double pivot;
  char bounded;

  /* Auxiliary flags to avoid unnecessary function evaluations */
  char testIntra, testInter;

  /* Function value variables. */
  double V, V0, Vneg, Vpos, Vprime;

  /* Initial position. */
  CopyVector(theta, c->theta, n);
  CopyVector(thetaPrime, theta, n);

  /* Compute the function value at the initial position. */
  Vprime = f( theta, fContext );

  do
  {
    if ( verbose )
    {
      printf( "iter #%lu\n", bigIter++ );
      PrintFunctionEval( "  init:", thetaPrime, n, Vprime );
    }

    V0 = Vprime;
    CopyVector(theta, thetaPrime, n);
    InitVector(deltaPrimeP, 0.0, n);
    InitVector(VP, Vprime, n);
    testInter = 0;

    for (j=0; j < rScales && k < numIterations; j++) //Scales
    {
      if ( verbose ) printf( "\n  scale #%lu\n", j+1 );

      testIntra = 0;

      for (i=0; i < n && k < numIterations; i++) //Parameters
      {
        if ( verbose ) printf( "    par #%lu\n", i+1 );

        deltaPrimeS[i] = 0;
        V = V0;
        pivot = theta[i];
        theta[i] =  pivot + ( delta[i][j] * r );

        bounded = BoundTest(theta, i, lowerBound, upperBound);
        if (!bounded)
        {
          Vpos = f( theta, fContext ); k++;
          if ( Vpos MinMaxOp V )      { V = Vpos;     deltaPrimeS[i] = delta[i][j] * r; testIntra = 1; }
          if ( Vpos MinMaxOp VP[i] )  { VP[i] = Vpos; deltaPrimeP[i] = delta[i][j] * r; testInter = 1; }
          if ( Vpos MinMaxOp Vprime ) { Vprime = Vpos; }

          if (verbose ) PrintFunctionEval( "      pos:", theta, n, Vpos );
        }

        theta[i] =  pivot - ( delta[i][j] * r );
        bounded = BoundTest(theta, i, lowerBound, upperBound);
        if (!bounded)
        {
          Vneg = f( theta, fContext ); k++;
          if ( Vneg MinMaxOp V )      { V = Vneg;      deltaPrimeS[i] = -delta[i][j] * r; testIntra = 1; }
          if ( Vneg MinMaxOp VP[i] )  { VP[i] = Vneg;  deltaPrimeP[i] = -delta[i][j] * r; testInter = 1; }
          if ( Vneg MinMaxOp Vprime ) { Vprime = Vneg; }

          if (verbose ) PrintFunctionEval( "      neg:", theta, n, Vneg );
        }
        theta[i] = pivot;
      }
      //Test the scale's composed displacement
      if ( testIntra )
      {
        SumVector(thetaTmp, theta, deltaPrimeS, n);
        V = f( thetaTmp, fContext ); k++;
        if ( V MinMaxOp Vprime )  { Vprime = V;  CopyVector(thetaPrime, thetaTmp, n); }

        if (verbose ) PrintFunctionEval( "\n    comp:", thetaTmp, n, V );
      }
    }

    //Test the interscale composed displacement
    if ( rScales > 1 && testInter )
    {
      SumVector(thetaTmp, theta, deltaPrimeP, n);
      V = f(thetaTmp, fContext ); k++;
      if ( V MinMaxOp Vprime )  { Vprime = V; CopyVector(thetaPrime, thetaTmp, n); }

      if (verbose ) PrintFunctionEval( "\n  inter-comp.:", thetaTmp, n, V );
    }

    /***************************************************************
     * Test the best scale or interscale composed displacement
     * against the individual parameter best fitness.
     *
     * Here there is a difference between the pseudo-code and the
     * implementation due to efficiency.
     *
     * The initial sum is just to ensure that if V==Vprime after the
     * loop, then the best displacement is certainly based on only
     * one of the parameters.
     **************************************************************/
    if ( Vprime MinMaxOp Vprime + 1.0 )
      V = Vprime + 1.0; //minimization
    else
      V = Vprime - 1.0; //maximization

    for (i=0; i< n; i++) //parameters (dimensions)
    {
      if ( VP[i] MinMaxOp V ) { bestP = i; V = VP[i]; }
    }
    if ( V == Vprime ) //implies that the best fitness comes from the axis
    {
      CopyVector(thetaTmp, theta, n);
      thetaTmp[bestP] = thetaTmp[bestP] + deltaPrimeP[bestP];
      CopyVector(thetaPrime, thetaTmp, n);
    }

    if (verbose )
    {
      PrintFunctionEval( "best:", thetaPrime, n, Vprime ); printf("\n");
    }

    //If no improvement was obtained, decrease the sampling-range.
    if ( Vprime == V0)
      r *= q;

  } while ( k < numIterations );

  CopyVector( c->theta, thetaPrime, n );

  /* Delete thetas and deltas. */
  for (i=0; i < n; i++ )
    free(delta[i]);

  free(delta);

  free(theta);
  free(thetaTmp);
  free(thetaPrime);

  free(deltaPrimeS);
  free(deltaPrimeP);
  free(VP);

  /* Return best-found function value. */
  return Vprime;
}

static void InitVector(double *v,  double value, size_t n)
{
  size_t i;

  assert(v);

  for (i=0; i<n; i++)
  {
    v[i] = value;
  }
}
static void CopyVector(double *dest,  double *src, size_t n)
{
  size_t i;

  assert(dest);
  assert(src);

  for (i=0; i<n; i++)
  {
    dest[i] = src[i];
  }
}

static void SumVector(double *dest,  double *a,  double *b, size_t n)
{
   size_t i;

   assert(dest);
   assert(a);
   assert(b);

   for (i=0; i<n; i++)
   {
      dest[i] = a[i]+b[i];
   }
}
static char BoundTest(double *x, size_t i,  double *lower,  double *upper)
{
   assert(upper[i] >= lower[i]);

   if (x[i] < lower[i])
   {
      x[i] = lower[i];
      return 1;
   }
   else if (x[i] > upper[i])
   {
      x[i] = upper[i];
      return 1;
   }
   return 0;
}

static void PrintFunctionEval( char *label,  double *theta, size_t n, double v )
{
  size_t i;

  printf("%s ", label);
  for ( i=0; i<n; i++ )
  {
    printf( "%+02.2f ",theta[i] );
    fflush(stdout);
  }

  printf("- v: %f \n", v );
}



static float Cofator(float M[4][4], int l, int c){

	float aux[3][3];
	int linhas[3]={0,0,0},colunas[3]={0,0,0};
	int i,j;
	float cof=0;
	int aux_sinal;

	if(l==0){
		linhas[0]=1; linhas[1]=2; linhas[2]=3;
	} 
	else if(l==1){
		linhas[0]=0; linhas[1]=2; linhas[2]=3;
	}
	else if(l==2){
		linhas[0]=0; linhas[1]=1;linhas[2]=3;
	}
	else if(l==3){
		linhas[0]=0; linhas[1]=1;linhas[2]=2;
	}
	if(c==0){
		colunas[0]=1; colunas[1]=2; colunas[2]=3;
	} 
	else if(c==1){
		colunas[0]=0; colunas[1]=2;colunas[2]=3;
	}
	else if(c==2){
		colunas[0]=0; colunas[1]=1;colunas[2]=3;
	}
	else if(c==3){
		colunas[0]=0; colunas[1]=1;colunas[2]=2;
	}
	
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			aux[i][j]=M[linhas[i]][colunas[j]];
		}
	}

	cof=aux[0][0]*aux[1][1]*aux[2][2] + aux[0][1]*aux[1][2]*aux[2][0] + aux[0][2]*aux[1][0]*aux[2][1]
		 -aux[0][0]*aux[1][2]*aux[2][1] - aux[0][1]*aux[1][0]*aux[2][2] - aux[0][2]*aux[1][1]*aux[2][0];
	
	aux_sinal=l+c;
	if ((aux_sinal%2) ==1)
		cof=-1*cof;
	return(cof);
}

static void TransMatrix(float M1[4][4],float M2[4][4])
{
  int i,j;

  for (j=0; j < 4; j++)
    for (i=0; i < 4; i++)
      M2[i][j]=M1[j][i];
}



static void inversa(float M[4][4], float IM[4][4]){

	float detM;
	float aux[4][4],M_cofatores[4][4],Trans_cofat[4][4];
	int i,j;
	
	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			aux[i][j]=M[i][j];
			//printf("[%f] ",M[i][j]);
		}
		//printf("\n");
	}
			//calcula o determinante de M (a ultima linha de M é 0 0 0 1 - por laplace faco 1*Cofator M3x3)
	detM = aux[0][0]*aux[1][1]*aux[2][2] + aux[0][1]*aux[1][2]*aux[2][0] + aux[0][2]*aux[1][0]*aux[2][1]
			  -aux[0][0]*aux[1][2]*aux[2][1] - aux[0][1]*aux[1][0]*aux[2][2] - aux[0][2]*aux[1][1]*aux[2][0];
	
	
	if(detM==0)
		Error(MSG1,"Matriz Inversa (determinante==0)");
 
	// calculo da matriz dos cofatores
	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			M_cofatores[i][j]=Cofator(aux,i,j);
		}
	}

	//Transposta da Matriz dos cofatores
	TransMatrix(M_cofatores,Trans_cofat);

	//dividir a Matriz Transposta pelo determinante de M
	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			IM[i][j]=Trans_cofat[i][j]/detM;
		}

	}

}

static Point TransformPoint(float M[4][4], Point p)
{
  Point np;

  np.x = M[0][0]*p.x + M[0][1]*p.y + M[0][2]*p.z + M[0][3];
  np.y = M[1][0]*p.x + M[1][1]*p.y + M[1][2]*p.z + M[1][3];
  np.z = M[2][0]*p.x + M[2][1]*p.y + M[2][2]*p.z + M[2][3];

  return(np);
}

static int valueInterpol3(Point p, Scene *scn){

  float d_q12_q2, d_q12_q1, d_q34_q3, d_q34_q4, d_q1_q2, d_q3_q4,d_q12_q34;
  float d_q56_q6, d_q56_q5, d_q78_q7, d_q78_q8, d_q5_q6, d_q7_q8,d_q56_q78;
  float d_p1_q12,d_p1_q34,d_p2_q56,d_p2_q78,d_p1_p2;
  float d_p_p1,d_p_p2,d_p_q12,d_p_q34;
  float I_q12,I_q34,I_q56,I_q78;
  float I_p1,I_p2;
  int I_p,I_q1,I_q2,I_q3, I_q4,I_q5,I_q6,I_q7,I_q8;
  int x,X,y,Y,z,Z; //X,Y maior inteiro proximo, x,y:menor inteiro proximo

  x=(int)p.x;
  X=(int)p.x+1;
  y=(int)p.y;
  Y=(int)p.y+1;
  z=(int)p.z;
  Z=(int)p.z+1;

  //if point is in the center of the slice
  if(z==p.z){
     if((x==p.x) && (y==p.y))
        return scn->data[scn->tbz[z]+scn->tby[y]+x];
    
    // point has 2 neighboors in x direction
    if(y==p.y){
      I_q1=scn->data[scn->tbz[z]+scn->tby[y]+x];
      I_q2=scn->data[scn->tbz[z]+scn->tby[y]+X];
      d_q12_q1=p.x-x; 
      d_q12_q2=X-p.x;
      d_q1_q2=X-x;
      I_q12=(d_q12_q2*I_q1 + d_q12_q1*I_q2)/d_q1_q2;
      return ((int)I_q12);
    }
    
    // point has 2 neighboors in y direction
    if (x==p.x){
      I_q1=scn->data[scn->tbz[z]+scn->tby[y]+x];
      I_q2=scn->data[scn->tbz[z]+scn->tby[Y]+x];
      d_q12_q1=p.y-y; 
      d_q12_q2=Y-p.y;
      d_q1_q2=Y-y;
      I_q12=(d_q12_q2*I_q1 + d_q12_q1*I_q2)/d_q1_q2;
      return ((int)I_q12);;
    }
    
    //if point has 4 neighboors
    I_q1=scn->data[scn->tbz[z]+scn->tby[y]+x]; I_q2=scn->data[scn->tbz[z]+scn->tby[y]+X];
    I_q3=scn->data[scn->tbz[z]+scn->tby[Y]+x]; I_q4=scn->data[scn->tbz[z]+scn->tby[Y]+X];

    d_q12_q1=d_q34_q3=p.x-x; 
    d_q34_q4=d_q12_q2=X-p.x;

    d_p_q12=p.y-y;
    d_p_q34=Y-p.y;

    d_q1_q2=d_q3_q4=X-x;
    d_q12_q34=Y-y;

    I_q12=(d_q12_q2*I_q1 + d_q12_q1*I_q2)/d_q1_q2;
    I_q34=(d_q34_q4*I_q3 + d_q34_q3*I_q4)/d_q3_q4;
    I_p= (int)((d_p_q34*I_q12 + d_p_q12*I_q34)/d_q12_q34);

    if(I_p<0){

        printf("Negativo: %d  \nPonto[%f][%f]  Vizinhos[%d][%d]  [%d][%d] [%d][%d] [%d][%d] \n",I_p,p.x,p.y,x,y,X,y,x,Y,X,Y);
    }
    return I_p;
  }
  
  //if point has 2 neighbors in z direction
  if(y==p.y && x==p.x && z!=p.z){
      I_q1=scn->data[scn->tbz[(int)z]+ scn->tby[(int)y]+(int)x];
      I_q2=scn->data[scn->tbz[(int)Z]+ scn->tby[(int)y]+(int)x];
      d_q12_q1 =p.z-z;
      d_q12_q2 =Z-p.z;
      d_q1_q2=Z-z;
      I_q12=(d_q12_q2*I_q1 + d_q12_q1*I_q2)/d_q1_q2;
      return((int)I_q12);
  }
  
 
  
  //if point has the 8 neighbors 
  I_q1=scn->data[scn->tbz[(int)z]+ scn->tby[(int)y]+(int)x]; I_q2=scn->data[scn->tbz[(int)z]+ scn->tby[(int)y]+(int)X];
  I_q3=scn->data[scn->tbz[(int)z]+ scn->tby[(int)Y]+(int)x]; I_q4=scn->data[scn->tbz[(int)z]+ scn->tby[(int)Y]+(int)X];
  
  I_q5=scn->data[scn->tbz[(int)Z]+ scn->tby[(int)y]+(int)x]; I_q6=scn->data[scn->tbz[(int)Z]+ scn->tby[(int)y]+(int)X];
  I_q7=scn->data[scn->tbz[(int)Z]+ scn->tby[(int)Y]+(int)x]; I_q8=scn->data[scn->tbz[(int)Z]+ scn->tby[(int)Y]+(int)X];
  
  
  d_q12_q1 = d_q34_q3 = d_q56_q5 = d_q78_q7 = p.x-x; 
  d_q34_q4 = d_q12_q2 = d_q56_q6 = d_q78_q8 = X-p.x;
  
  d_p1_q12 = d_p2_q56 = p.y-y;
  d_p1_q34 = d_p2_q78 = Y-p.y;
  
  d_q1_q2 = d_q3_q4 = d_q5_q6 = d_q7_q8 = X-x;
  d_q12_q34 = d_q56_q78= Y-y;
  
  I_q12=(d_q12_q2*I_q1 + d_q12_q1*I_q2)/d_q1_q2;
  I_q34=(d_q34_q4*I_q3 + d_q34_q3*I_q4)/d_q3_q4;
  I_p1=(d_p1_q34*I_q12 + d_p1_q12*I_q34)/d_q12_q34;
  
  I_q56=(d_q56_q6*I_q5 + d_q56_q5*I_q6)/d_q5_q6;
  I_q78=(d_q78_q8*I_q7 + d_q78_q7*I_q8)/d_q7_q8;
  I_p2=(d_p2_q78*I_q56 + d_p2_q56*I_q78)/d_q56_q78;
  
  d_p_p1=p.z-z;
  d_p_p2=Z-p.z;
  d_p1_p2=Z-z;
  
  I_p= (int)((d_p_p2*I_p1 + d_p_p1*I_p2)/d_p1_p2);
  if(I_p<0){
    printf("Negativo: %d  \nPonto[%f][%f]  Vizinhos[%d][%d]  [%d][%d] [%d][%d] [%d][%d]\n",I_p,p.x,p.y,x,y,X,y,x,Y,X,Y);
  }
  
  return I_p; 
}



static void transformScene(Scene *scn, float T[4][4], Scene *scn_out){
	int i, j, k;
      //  Scene *scn_out=CreateScene(scn->xsize,scn->ysize,scn->zsize);
        Point orig, dest;
        float IM[4][4];
        inversa(T, IM);
        for(k=0;k<scn_out->zsize;k++){
            for(j=0;j<scn_out->ysize;j++){
                for(i=0;i<scn_out->xsize;i++){
                    orig.x=i;
                    orig.y=j;
                    orig.z=k;
                    dest=TransformPoint(IM, orig);
                    if((dest.x>=0) && (dest.y>=0) && (dest.z>=0) && (dest.x<=scn->xsize-1) && (dest.y<=scn->ysize-1) && (dest.z<=scn->zsize-1)
                            && ValidVoxel(scn, (int)dest.x, (int)dest.y, (int)dest.z)){
                                       scn_out->data[scn_out->tbz[(int)orig.z]+scn_out->tby[(int)orig.y]+(int)orig.x]=(int)valueInterpol3(dest, scn);
                    }
                }
            }
        }
        scn_out->maxval=MaximumValue3(scn_out);
	scn_out->dx=scn->dx;
	scn_out->dy=scn->dy;
	scn_out->dz=scn->dz;
}
















//--------------------------------------------------------------------------------
// PRIVATE - Registration-specific funcions - Using the old ift library





static Scene *NormalizeAccHist3(Scene *scn) {
    Scene *out = CreateScene(scn->xsize, scn->ysize, scn->zsize);
    Curve *acc = NormAccHistogram3(scn);
    int max, i, n, min;
    for (i = acc->n - 1; acc->Y[i] > 0.991; i--);
    max = i;
    for (i = 1; acc->Y[i] < 0.1; i++);
    min = i;
  

    n = scn->xsize * scn->ysize * scn->zsize;
    for (i = 0; i < n; i++) {
        if (scn->data[i] < min)
            out->data[i] = 0;

        else if (scn->data[i] <= max) {
            out->data[i] = ((scn->data[i] - min)*4095) / (max - min);
        } else
            out->data[i] = 4095;
    }
    DestroyCurve(&acc);
    return (out);
}

static float Distance(float *f1, float *f2, int n){
  int i;
  float dist;

  dist = 0;
  for (i=0; i < n; i++)
    dist += (f2[i]-f1[i]);//(f2[i]-f1[i])/2.0;
  //dist /= n;

  return(dist);//exp(-(dist-0.5)*(dist-0.5)/2.0));
}

static Scene *TextGradient(Scene *scn){
  float   dist,gx,gy,gz;
  int     i,p,q,n=scn->xsize*scn->ysize*scn->zsize;//,Imax=MaximumValue3(scn);
  Voxel   u,v;
  AdjRel3 *A=Spheric(1.0),*A6=Spheric(1.0);
  float   *mg=AllocFloatArray(A6->n);
  Scene   *grad=CreateScene(scn->xsize,scn->ysize,scn->zsize);


  typedef struct _features {
    float *f;
  } Features;

  Features *feat=(Features *)calloc(n,sizeof(Features));
  for (p=0; p < n; p++)
    feat[p].f = AllocFloatArray(A->n);

  for (u.z=0; u.z < scn->zsize; u.z++)
    for (u.y=0; u.y < scn->ysize; u.y++)
      for (u.x=0; u.x < scn->xsize; u.x++) {
        p = u.x + scn->tby[u.y] + scn->tbz[u.z];
        for (i=0; i < A->n; i++) {
          v.x = u.x + A->dx[i];
          v.y = u.y + A->dy[i];
          v.z = u.z + A->dz[i];
          if (ValidVoxel(scn,v.x,v.y,v.z)){
            q = v.x + scn->tby[v.y] + scn->tbz[v.z];
            feat[p].f[i]=(float)scn->data[q];///(float)Imax;
          }
        }
      }

  for (i=0; i < A6->n; i++)
    mg[i]=sqrt(A6->dx[i]*A6->dx[i]+A6->dy[i]*A6->dy[i]+A6->dz[i]*A6->dz[i]);

  for (u.z=0; u.z < scn->zsize; u.z++)
    for (u.y=0; u.y < scn->ysize; u.y++)
      for (u.x=0; u.x < scn->xsize; u.x++) {
        p = u.x + scn->tby[u.y] + scn->tbz[u.z];
        gx = gy = gz = 0.0;
        for (i=1; i < A6->n; i++) {
          v.x = u.x + A6->dx[i];
          v.y = u.y + A6->dy[i];
          v.z = u.z + A6->dz[i];
          if (ValidVoxel(scn,v.x,v.y,v.z)){
            q = v.x + scn->tby[v.y] + scn->tbz[v.z];
            dist = Distance(feat[p].f,feat[q].f,A->n);
            gx  += dist*A6->dx[i]/mg[i];
            gy  += dist*A6->dy[i]/mg[i];
            gz  += dist*A6->dz[i]/mg[i];
          }
        }
        grad->data[p]=(int)sqrt(gx*gx + gy*gy + gz*gz);//(100000.0*sqrt(gx*gx + gy*gy + gz*gz));
      }


  for (p=0; p < n; p++)
    free(feat[p].f);
  free(feat);

  free(mg);
  DestroyAdjRel3(&A);
  DestroyAdjRel3(&A6);
  return(grad);
}


static Scene * getWaterGray3(Scene *scn, float di, float factor){
  Scene *grad=NULL,*handicap=NULL;
  Scene    *label=NULL;

  AdjRel3   *A=NULL;

  grad  = TextGradient(scn);

  A = Spheric(di); 
  handicap = Add3(grad,(int)(factor*MaximumValue3(grad)));
  label = WaterGray3(grad,handicap,A);
  
  DestroyAdjRel3(&A);
  DestroyScene(&grad);  
  DestroyScene(&handicap);  
  
  return(label);
}


static Scene *getWGBorder(Scene *label){

  //Scene *label=getWaterGray3(scn,1.8,0.02);
	
  Scene *hscn=CreateScene(label->xsize,label->ysize,label->zsize);
  int p,q,i;
  AdjRel3 *A=NULL;
  Voxel u,v;

  A    = Spheric(1.0);
  for (u.z=0; u.z < hscn->zsize; u.z++){
    for (u.y=0; u.y < hscn->ysize; u.y++){
      for (u.x=0; u.x < hscn->xsize; u.x++){
	p = u.x + hscn->tby[u.y] + hscn->tbz[u.z];
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  v.z = u.z + A->dz[i];
	  if (ValidVoxel(hscn,v.x,v.y,v.z)){
	    q = v.x + hscn->tby[v.y] + hscn->tbz[v.z];
	    if (label->data[p] < label->data[q]){
	      hscn->data[p] = 1;
	      break;
	    }
	  }
	}
      }
    }
  }
  DestroyAdjRel3(&A);
  return(hscn);
}	

static Set *getSetBorder(Scene *scn){
    
    AdjRel3 *adj=Spheric(1.0);
    Scene *borda=GetBorder3(scn, adj);
    Set *S=NULL;
    int i, n;
    n=scn->xsize*scn->ysize*scn->zsize;
    for(i=0;i<n;i++){
        if(borda->data[i]==1){
            InsertSet(&S, i);
        }
    }
    DestroyScene(&borda);
    DestroyAdjRel3(&adj);
    return S;
}





// Parameters: 
// theta - array with the  parameters of the transform
// n - # of elements in theta
// center - the center of the rotation
// RT - is the resulting transformation matrix
static void calculateTransformation(float *theta, int n, Point center, float RT[4][4]) {
  // Compute transformation in this order: centered rotation (center is cx,cy,cz), translation, scaling and shear. Parameters theta are Rx, Ry, Rz, Tx, Ty, Tz, Sx, Sy, Sz, SHxy, SHxz, SHyz, SHyz, SHzx, SHzy.
  // The matrix multiplication is done backwards:
  // M = T SH S -Tc Rz Ry Rx Tc
  float Rx=0,Ry=0,Rz=0,Tx=0,Ty=0,Tz=0;
  float Sx=1,Sy=1,Sz=1,SHxy=0,SHxz=0,SHyx=0,SHyz=0,SHzx=0,SHzy=0;
  RealMatrix *trans1,*rot,*rot1,*rot2,*rot3,*trans2;
  RealMatrix *aux1,*aux2;
  RealMatrix *trans, *scale, *shear, *final;

  if (n>=3) {
    Rx = theta[0] * PI / 180.0;
    Ry = theta[1] * PI / 180.0;
    Rz = theta[2] * PI / 180.0;
  }
  if (n>=6) {
    Tx = theta[3];
    Ty = theta[4];
    Tz = theta[5];
  }
  if (n>=9) {
    Sx = theta[6];
    Sy = theta[7];
    Sz = theta[8];
  }
  if (n==15) {
    SHxy = theta[9];
    SHxz = theta[10];
    SHyx = theta[11];
    SHyz = theta[12];
    SHzx = theta[13];
    SHzy = theta[14];
  }

  // Rotation
  trans1 = TranslationMatrix3(-center.x,-center.y,-center.z);
  rot1 = RotationMatrix3(0,Rx);
  rot2 = RotationMatrix3(1,Ry);
  rot3 = RotationMatrix3(2,Rz);
  trans2 = TranslationMatrix3(center.x,center.y,center.z);
  // Compose transform rot=trans2 x rot3 x rot2 x rot1 x trans1
  aux1 = MultRealMatrix(trans2,rot3);
  aux2 = MultRealMatrix(aux1,rot2);
  DestroyRealMatrix(&aux1);
  aux1 = MultRealMatrix(aux2,rot1);
  DestroyRealMatrix(&aux2);
  rot = MultRealMatrix(aux1,trans1);
  DestroyRealMatrix(&aux1);
  DestroyRealMatrix(&rot1);
  DestroyRealMatrix(&rot2);
  DestroyRealMatrix(&rot3);
  DestroyRealMatrix(&trans1);
  DestroyRealMatrix(&trans2);

  // Translation
  trans = TranslationMatrix3(Tx,Ty,Tz);

  // Scale
  scale = ScaleMatrix3(Sx,Sy,Sz);

  // Shear
  shear = ShearMatrix3(SHxy,SHxz,SHyx,SHyz,SHzx,SHzy);

  // compose final final = trans x shear x scale x rot
  aux1 = MultRealMatrix(trans,shear);
  aux2 = MultRealMatrix(aux1,scale);
  DestroyRealMatrix(&aux1);
  final = MultRealMatrix(aux2,rot);
  DestroyRealMatrix(&aux2);
  DestroyRealMatrix(&rot);
  DestroyRealMatrix(&trans);
  DestroyRealMatrix(&scale);
  DestroyRealMatrix(&shear);

  RT[0][0]=final->val[0][0];
  RT[0][1]=final->val[0][1];
  RT[0][2]=final->val[0][2];
  RT[0][3]=final->val[0][3];
  RT[1][0]=final->val[1][0];
  RT[1][1]=final->val[1][1];
  RT[1][2]=final->val[1][2];
  RT[1][3]=final->val[1][3];
  RT[2][0]=final->val[2][0];
  RT[2][1]=final->val[2][1];
  RT[2][2]=final->val[2][2];
  RT[2][3]=final->val[2][3];
  RT[3][0]=final->val[3][0];
  RT[3][1]=final->val[3][1];
  RT[3][2]=final->val[3][2];
  RT[3][3]=final->val[3][3];
  DestroyRealMatrix(&final);  
}




/*

void calcTransformation(float *theta, float RT[4][4], float xc, float yc, float zc) {
    //calcula a matriz de tranformacao (rotacao e translacao) RT a partir dos parametros theta[i]
    // a transformacao é centrada no ponto de centro de gravidade (xc,yc,zc)
    
    //M= T -Tc R Tc
    float thx, thy, thz, X, Y, Z;
    thx = theta[0];
    thy = theta[1];
    thz = theta[2];
    if (thx < 0)
        thx = 360 + thx;
    if (thy < 0)
        thy = 360 + thy;
    if (thz < 0)
        thz = 360 + thz;
    //radianos
    X = thx * PI / 180;
    Y = thy * PI / 180;
    Z = thz * PI / 180;

    float T[4][4];
    T[0][0] = 1.0;
    T[0][1] = 0.0;
    T[0][2] = 0.0;
    T[0][3] = theta[3];
    T[1][0] = 0.0;
    T[1][1] = 1.0;
    T[1][2] = 0.0;
    T[1][3] = theta[4];
    T[2][0] = 0.0;
    T[2][1] = 0.0;
    T[2][2] = 1.0;
    T[2][3] = theta[5];
    T[3][0] = 0.0;
    T[3][1] = 0.0;
    T[3][2] = 0.0;
    T[3][3] = 1.0;
    float S=theta[6];
    
    float R[4][4];

    R[0][0] = S*(cos(Z) * cos(Y));
    R[1][0] = S*(sin(Z) * cos(Y) * cos(X) - sin(Y) * sin(X));
    R[2][0] = S*(-1. * sin(Z) * cos(Y) * sin(X) - sin(Y) * cos(X));
    R[3][0] = 0.0;

    R[0][1] = S*(-1. * sin(Z));
    R[1][1] = S*(cos(Z) * cos(X));
    R[2][1] = S*(-1. * cos(Z) * sin(X));
    R[3][1] = 0.0;

    R[0][2] = S*(cos(Z) * sin(Y));
    R[1][2] = S*(sin(Z) * sin(Y) * cos(X) + cos(Y) * sin(X));
    R[2][2] = S*(-1. * sin(Z) * sin(Y) * sin(X) + cos(Y) * cos(X));
    R[3][2] = 0.0;

    R[0][3] = 0.0;
    R[1][3] = 0.0;
    R[2][3] = 0.0;
    R[3][3] = 1.0;

    float TC[4][4];
    TC[0][0] = 1.0;
    TC[0][1] = 0.0;
    TC[0][2] = 0.0;
    TC[0][3] = (float) - 1 * xc;
    TC[1][0] = 0.0;
    TC[1][1] = 1.0;
    TC[1][2] = 0.0;
    TC[1][3] = (float) - 1 * yc;
    TC[2][0] = 0.0;
    TC[2][1] = 0.0;
    TC[2][2] = 1.0;
    TC[2][3] = (float) - 1 * zc;
    TC[3][0] = 0.0;
    TC[3][1] = 0.0;
    TC[3][2] = 0.0;
    TC[3][3] = 1.0;

    float T_C[4][4]; // -TC
    T_C[0][0] = 1.0;
    T_C[0][1] = 0.0;
    T_C[0][2] = 0.0;
    T_C[0][3] = (float) xc;
    T_C[1][0] = 0.0;
    T_C[1][1] = 1.0;
    T_C[1][2] = 0.0;
    T_C[1][3] = (float) yc;
    T_C[2][0] = 0.0;
    T_C[2][1] = 0.0;
    T_C[2][2] = 1.0;
    T_C[2][3] = (float) zc;
    T_C[3][0] = 0.0;
    T_C[3][1] = 0.0;
    T_C[3][2] = 0.0;
    T_C[3][3] = 1.0;

    float m1[4][4], m2[4][4];
 
    //M= T -Tc R Tc
    MultMatrices(T, T_C, m1);
    MultMatrices(m1, R, m2);
    MultMatrices(m2, TC, RT);

   
}

*/


/* The context-struct for benchmark problems. */
  struct RegContext
  {
    size_t   n;          /* Dimensionality. */
    Scene *fixed_grad; // gradient of fixed image
    Scene *moving; // moving image
    Set *moving_set; // set of watershed line points moving image)
    Point center;
  };

/*************** Function to be optimized *********************/
static double CriteryFunction(const double *x, void *context)
{
  struct RegContext const* c = (struct RegContext const*) context;
  float T[4][4];

  float *theta_aux=AllocFloatArray(9); // 7 parameters (3 rot, 3 trans, 1 scale)
  theta_aux[0]=(float)x[0];
  theta_aux[1]=(float)x[1];
  theta_aux[2]=(float)x[2];
  theta_aux[3]=(float)x[3];
  theta_aux[4]=(float)x[4];
  theta_aux[5]=(float)x[5];
  theta_aux[6]=(float)x[6]; // Sx=Sy=Sz
  theta_aux[7]=(float)x[6];
  theta_aux[8]=(float)x[6];
  calculateTransformation(theta_aux, 9, c->center, T);

  //calcula T(Sj) e a distancia D ao mesmo tempo
  Set *S = c->moving_set;
  Scene *fixed_grad = c->fixed_grad;
  Scene *moving = c->moving;
  int p;
  double D=0.0; 
  Voxel v, v_transf;
  int counter=0;
  while (S != NULL) {
    counter++;
    p = S->elem;
    v.z = p / (moving->xsize * moving->ysize);
    v.y = (p - moving->tbz[v.z]) / (moving->xsize);
    v.x = (p - moving->tbz[v.z]) % (moving->xsize);
    v_transf = Transform_Voxel(T, v);
    if (ValidVoxel(fixed_grad, v_transf.x, v_transf.y, v_transf.z))
      D += fixed_grad->data[fixed_grad->tbz[v_transf.z] + fixed_grad->tby[v_transf.y] + v_transf.x];
    S = S->next;
  }
  free(theta_aux);
  if (counter!=0) 
    D = (D/counter)/4096;
  else D=0;
  //printf("- Score: %f\n",D);

  return D;
}






// input: fixa, movel
// output: T,best_theta
static Scene *Register3(Scene *fixa, Scene *movel, float T[4][4], float **best_theta)
{
  int i;
  int debug=0;
  AdjRel3 *adj=NULL;
  Scene *B=NULL;
  Set *S=NULL;
  Point centro;
  Scene *fixa_norm=NULL, *movel_norm=NULL;
  Scene *movel_wg=NULL, *movel_wglines=NULL;
  Scene *out=NULL;
  
  printf("[gradient] "); fflush(stdout);
  //imagem de borda realcada 
  adj = Spheric(1.0);
  fixa_norm = NormalizeAccHist3(fixa);
  //B = MorphGrad3(fixa_norm,adj);
  Scene *tmp = TextGradient(fixa_norm);
  B  = NormalizeAccHist3(tmp);
  DestroyScene(&tmp);
  //WriteScene(B,"B.scn");
  DestroyAdjRel3(&adj);
  DestroyScene(&fixa_norm);
 
  printf("[watergray] ");fflush(stdout);
  //Conjunto de pontos (S) da imagem movel
  movel_norm = NormalizeAccHist3(movel); 
  movel_wg = getWaterGray3(movel_norm, 1.8, 0.07); 
  DestroyScene(&movel_norm);
  movel_wglines = getWGBorder(movel_wg);
  DestroyScene(&movel_wg);
  S = getSetBorder(movel_wglines);
  DestroyScene(&movel_wglines);
  
  //centro = CalcCenterOfGravity(S, movel);
  centro.x = movel->xsize/2;
  centro.y = movel->ysize/2;
  centro.z = movel->zsize/2;

  printf("[search]\n "); fflush(stdout);

  /********************************************
   * Problem layer
   ********************************************/
  int dimensions = 7; //7 parameters
  struct RegContext regContext;
  regContext.n = dimensions;
  regContext.fixed_grad = B;
  regContext.moving = movel;
  regContext.moving_set = S;
  regContext.center = centro;

  //theta from which the optimization starts. it can be a random point as well.
  double theta0[dimensions];

  double lowerBound[dimensions];
  double upperBound[dimensions];
  double lowerInit[dimensions];
  double upperInit[dimensions];

  for ( i=0; i<dimensions; i++)
    theta0[i] = 0.0;
  theta0[6] = 1.0; //scale

  lowerBound[0]=-15;  upperBound[0]=15;
  lowerBound[1]=-15;  upperBound[1]=15;
  lowerBound[2]=-15;  upperBound[2]=15;
  lowerBound[3]=-(B->xsize/12);  upperBound[3]=B->xsize/12;
  lowerBound[4]=-(B->ysize/12);  upperBound[4]=B->ysize/12;
  lowerBound[5]=-(B->zsize/12);  upperBound[5]=B->zsize/12;
  lowerBound[6]=0.95;  upperBound[6]=1.05;

  lowerInit[0]=-15;  upperInit[0]=15;
  lowerInit[1]=-15;  upperInit[1]=15;
  lowerInit[2]=-15;  upperInit[2]=15;
  lowerInit[3]=-(B->xsize/12);  upperInit[3]=B->xsize/12;
  lowerInit[4]=-(B->ysize/12);  upperInit[4]=B->ysize/12;
  lowerInit[5]=-(B->zsize/12);  upperInit[5]=B->zsize/12;
  lowerInit[6]=0.95;  upperInit[6]=1.05;
  

  /********************************************
   * Optimization layer
   ********************************************/
  double settings[] = {6,3,1}; /* number of scales and its increasing degree */
  size_t  numIterations = 1000;    /* Number of iterations. */

  /* Define the context for the MSPS optimization */
  struct optContext mspsContext;

  mspsContext.f          =  (problem) CriteryFunction;
  mspsContext.fContext   = (void*) &regContext;
  mspsContext.fDim       = regContext.n;
  mspsContext.theta      = theta0;
  mspsContext.lowerInit = lowerInit;
  mspsContext.upperInit = upperInit;
  mspsContext.lowerBound = lowerBound;
  mspsContext.upperBound = upperBound;
  mspsContext.numIterations = numIterations;
  mspsContext.verbose = 0;

  //double Vprime;
  MSPS( settings, (void*) &mspsContext);

  DestroyScene(&B);
  DestroySet(&S);

  *best_theta = AllocFloatArray(9);
  (*best_theta)[0]=mspsContext.theta[0];
  (*best_theta)[1]=mspsContext.theta[1];
  (*best_theta)[2]=mspsContext.theta[2];
  (*best_theta)[3]=mspsContext.theta[3];
  (*best_theta)[4]=mspsContext.theta[4];
  (*best_theta)[5]=mspsContext.theta[5];
  (*best_theta)[6]=mspsContext.theta[6];
  (*best_theta)[7]=mspsContext.theta[6];
  (*best_theta)[8]=mspsContext.theta[6];
  printf("Parametros finais: RX=%.2f RY=%.2f RZ=%.2f TX=%.2f TY=%.2f TZ=%.2f S=%.2f\n",(*best_theta)[0],(*best_theta)[1],(*best_theta)[2],(*best_theta)[3],(*best_theta)[4],(*best_theta)[5],(*best_theta)[6]);fflush(stdout);

  calculateTransformation(*best_theta,9,centro,T);
  out = CreateScene(fixa->xsize, fixa->ysize, fixa->zsize );    
  transformScene( movel, T, out );
  if (debug==1) WriteScene(out,"out.scn");
  return out;
}







//---------------------------------------------------------------------
// PRIVATE -  Old-ift to New-iftiftImage *image, iftImage *mask,  convertion functions



static Scene *Img2Scn(iftImage *img) {
  Scene *scn=NULL;
  if (img!=NULL) {
    scn=CreateScene(img->xsize,img->ysize,img->zsize);
    if (scn!=NULL) {
      memcpy(scn->data,img->val,sizeof(int) * scn->n);
      scn->dx=img->dx;
      scn->dy=img->dy;
      scn->dz=img->dz;
    }
  }
  return scn;
}


static iftImage *Scn2Img(Scene *scn) {
  iftImage *img=NULL;
  if (scn!=NULL) {
    img=iftCreateImageFromBuffer(scn->xsize,scn->ysize,scn->zsize, scn->data);
    if (img!=NULL) {
      img->dx=scn->dx;
      img->dy=scn->dy;
      img->dz=scn->dz;
    }
  }
  return img;
}


   







//--------------------------------------------------------------------------
//
//  PUBLIC - MSP function for the new-ift code

iftImage *iftBrainAffineRegistration(iftImage *fixed, iftImage *moving){

  float T[4][4];
  float *best_theta;
  Scene *out_scn, *fixed_scn, *moving_scn;
  iftImage *out;

  fixed_scn = Img2Scn(fixed);
  moving_scn = Img2Scn(moving);
  out_scn = Register3(fixed_scn,moving_scn,T, &best_theta);
  out = Scn2Img(out_scn);
  //DestroyScene(&out_scn); // do not destroy it
  DestroyScene(&fixed_scn);
  DestroyScene(&moving_scn);
  return out;
  
}





