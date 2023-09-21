#include "iftGraphics.h"

#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftMathMorph.h"


/**
@file
@brief A description is missing here
*/
/* ------------------------------- Private Methods -------------------------*/

/* Methods used for surface and volume rendering */

float     iftPhongShading(iftGraphicalContext *gc, int p, iftVector N, float dist);
float     iftDepthShading(iftGraphicalContext *gc, int p, iftVector N, float dist);
iftVector iftObjectNormalOnTheFly(iftGraphicalContext *gc, iftImage *image, int po, iftAdjRel *A);
void      iftRenderFromSRBuffers(iftGraphicalContext *gc, iftImage *image);
void      iftResetSRBuffers(iftGraphicalContext *gc);
float     iftOpacityValueAtPoint(iftGraphicalContext *gc, iftPoint P);
iftVector iftNormalVectorAtPoint(iftGraphicalContext *gc, iftPoint P);
void      iftSurfaceShadingAlongRay(iftGraphicalContext *gc, iftPoint P0, iftPoint P1, iftPoint Pn, float *red, float *green, float *blue, iftAdjRel *A);
void      iftVoxelInfoByRayCasting(iftGraphicalContext *gc, iftPoint P0, iftPoint P1, iftPoint Pn, int *voxel, float *depth, int *object, iftAdjRel *A);
void      iftVolumeShadingAlongRay(iftGraphicalContext *gc, iftPoint P0, iftPoint P1, iftPoint Pn, float *value);
void      iftVoxelSplatting(iftGraphicalContext *gc, iftImage *image, iftPoint P, int p, float phong_val, float opac);
iftImage *iftSurfaceRenderingByRayCasting(iftGraphicalContext *gc);
iftImage *iftSurfaceRenderingBySplatting(iftGraphicalContext *gc);
iftImage *iftVolumeRenderingByRayCasting(iftGraphicalContext *gc);

/* Methods to create the graphical context */

void                 iftSetFTBaccess(iftGraphicalContext *gc);
iftPhongModel       *iftCreatePhongModel(iftFImage *scene);
iftSRBuffers        *iftCreateSRBuffers(iftImage *label);
iftObjectAttributes *iftCreateObjectAttributes(iftImage *label, int *nobjs);
int                  iftGetNormalIndex(iftVector normal);
iftVector           *iftCreateNormalTable();
void                 iftSetSceneFaces(iftGraphicalContext *gc);


float iftPhongShading(iftGraphicalContext *gc, int p, iftVector N, float dist)
{
    float cos_theta;
    float cos_2theta, pow, phong_val=0.0;

    cos_theta   = iftVectorInnerProduct(gc->viewdir->V, N);

    /* |angle| <= 90° */

    if (cos_theta >= IFT_EPSILON){

        cos_2theta = 2*cos_theta*cos_theta - 1;

        if (cos_2theta <= IFT_EPSILON){  /* |angle| >= 45° */
            pow = 0.;
        }else {
            pow = 1.;
            for (int k=0; k < gc->phong->ns; k++)
                pow = pow*cos_2theta;
        }

        phong_val = gc->phong->ka + gc->phong->Idist[(int)dist]*(gc->phong->kd*cos_theta + gc->phong->ks*pow);
    }

    return(phong_val);
}

float iftDepthShading(iftGraphicalContext *gc, int p, iftVector N, float dist)
{
    float cos_theta;
    float depth_val=0.0;

    cos_theta  = iftVectorInnerProduct(gc->viewdir->V, N);
    if (cos_theta > 1.0)
        cos_theta=1.0;

    if (cos_theta > IFT_EPSILON){  /* |angle| <= 90° */

        depth_val = gc->phong->Idist[(int)dist];

    }

    return(depth_val);
}


iftVector iftObjectNormalOnTheFly(iftGraphicalContext *gc, iftImage *image, int po, iftAdjRel *A)
{
    iftVector N, V, valid_neighbor[A->n];
    int i, p, q, n, qo;
    iftVoxel u,v,uo,vo;

    p  = gc->surf_render[po].voxel;
    u  = iftGetVoxelCoord(gc->label,p);

    uo = iftGetVoxelCoord(image,po);

    n = 0;
    for (i=1; i < A->n; i++) {
        vo = iftGetAdjacentVoxel(A,uo,i);
        if (iftValidVoxel(image,vo)){
            qo = iftGetVoxelIndex(image,vo);
            q  = gc->surf_render[qo].voxel;
            if (q != IFT_NIL) {
                v  = iftGetVoxelCoord(gc->label,q);
                int dist = iftSquaredVoxelDistance(u,v);
                if ((dist <= 100)&&(gc->label->val[q]==gc->label->val[p])){
                    valid_neighbor[n].x = v.x - u.x;
                    valid_neighbor[n].y = v.y - u.y;
                    valid_neighbor[n].z = v.z - u.z;
                    n++;
                }
            }
        }
    }

    N.x = N.y = N.z = 0.0;
    if ( n > 2 ) {
        for (i=0; i < n-1; i++) {
            V = (iftVector)iftVectorCrossProd(valid_neighbor[i],valid_neighbor[i+1]);
            N.x += V.x;
            N.y += V.y;
            N.z += V.z;
        }
        N.x /= -n; N.y /= -n; N.z /= -n;
        N = iftNormalizeVector(N);
    }

    return(N);
}

void iftRenderFromSRBuffers(iftGraphicalContext *gc, iftImage *image)
{
    iftAdjRel *A=iftClockCircular(3.0);

    for (int po=0; po < image->n; po++) {

        iftVector N;
        float     opac, dist, phong_val;
        int       p;
        iftColor  RGB, YCbCr;

        if (gc->surf_render[po].voxel != IFT_NIL) {

            N          = iftObjectNormalOnTheFly(gc,image,po,A);
            opac       = gc->object[gc->surf_render[po].object].opacity;
            p          = gc->surf_render[po].voxel;
            dist       = gc->surf_render[po].depth;

            phong_val  = opac * iftPhongShading(gc,p,N,dist);
            RGB.val[0] = (int)(255.0 * phong_val * gc->object[gc->surf_render[po].object].red);
            RGB.val[1] = (int)(255.0 * phong_val * gc->object[gc->surf_render[po].object].green);
            RGB.val[2] = (int)(255.0 * phong_val * gc->object[gc->surf_render[po].object].blue);
            YCbCr      = iftRGBtoYCbCr(RGB,255);

            image->val[po] = YCbCr.val[0];
            image->Cb[po]  = YCbCr.val[1];
            image->Cr[po]  = YCbCr.val[2];

        }
    }

    iftDestroyAdjRel(&A);
}

void iftResetSRBuffers(iftGraphicalContext *gc)
{
    int n = iftDiagonalSize(gc->label);

    n = n*n;
    for (int p=0; p < n; p++) {
        gc->surf_render[p].depth    = IFT_INFINITY_FLT;
        gc->surf_render[p].opacity  = 1.0;
        gc->surf_render[p].voxel    = IFT_NIL;
        gc->surf_render[p].object   = IFT_NIL;
    }
}

float iftOpacityValueAtPoint(iftGraphicalContext *gc, iftPoint P)
{
    iftVoxel u[8];
    int      p[8], i;
    float    dx=1.0,dy=1.0,dz=1.0, val[6];
    float    value;
    float    o[8];

    if ((int)(P.x+1.0)==gc->scene->xsize) dx = 0.0;
    if ((int)(P.y+1.0)==gc->scene->ysize) dy = 0.0;
    if ((int)(P.z+1.0)==gc->scene->zsize) dz = 0.0;

    u[0].x = (int)P.x;      u[0].y = (int)P.y;       u[0].z = (int)P.z;
    u[1].x = (int)(P.x+dx); u[1].y = (int)P.y;       u[1].z = (int)P.z;
    u[2].x = (int)P.x;      u[2].y = (int)(P.y+dy);  u[2].z = (int)P.z;
    u[3].x = (int)(P.x+dx); u[3].y = (int)(P.y+dy);  u[3].z = (int)P.z;
    u[4].x = (int)P.x;      u[4].y = (int)P.y;       u[4].z = (int)(P.z+dz);
    u[5].x = (int)(P.x+dx); u[5].y = (int)P.y;       u[5].z = (int)(P.z+dz);
    u[6].x = (int)P.x;      u[6].y = (int)(P.y+dy);  u[6].z = (int)(P.z+dz);
    u[7].x = (int)(P.x+dx); u[7].y = (int)(P.y+dy);  u[7].z = (int)(P.z+dz);

    if (gc->opacity != NULL){ /* volume rendering */
        for (i=0; i < 8; i++) {
            if (iftFValidVoxel(gc->scene,u[i])){
                p[i] = iftFGetVoxelIndex(gc->scene,u[i]);
                o[i]=gc->opacity->val[p[i]];
            }else{
                u[0].x = iftRound(P.x);
                u[0].y = iftRound(P.y);
                u[0].z = iftRound(P.z);
                p[0]   = iftFGetVoxelIndex(gc->scene,u[0]);
                return(gc->opacity->val[p[0]]);
            }
        }
    }else{ /* surface rendering */
        if (gc->nobjs > 0){
            for (i=0; i < 8; i++) {
                if (iftFValidVoxel(gc->scene,u[i])){
                    p[i] = iftFGetVoxelIndex(gc->scene,u[i]);
                    o[i] = gc->object[gc->label->val[p[i]]].opacity;
                }else{
                    u[0].x = iftRound(P.x);
                    u[0].y = iftRound(P.y);
                    u[0].z = iftRound(P.z);
                    p[0]   = iftFGetVoxelIndex(gc->scene,u[0]);
                    return(gc->object[gc->label->val[p[0]]].opacity);
                }
            }
        } else {
            iftError("Cannot interpolate opacities", "iftOpacityValueAtPoint");
        }
    }

    val[0] =(float)o[1]*(P.x-u[0].x)+(float)o[0]*(u[1].x-P.x);
    val[1] =(float)o[3]*(P.x-u[2].x)+(float)o[2]*(u[3].x-P.x);
    val[2] =(float)o[5]*(P.x-u[4].x)+(float)o[4]*(u[5].x-P.x);
    val[3] =(float)o[7]*(P.x-u[6].x)+(float)o[6]*(u[7].x-P.x);
    val[4] = val[1]*(P.y-u[0].y) + val[0]*(u[2].y-P.y);
    val[5] = val[3]*(P.y-u[0].y) + val[2]*(u[2].y-P.y);
    value  = (val[5]*(P.z-u[0].z) + val[4]*(u[4].z-P.z));

    return(value);
}

iftVector iftNormalVectorAtPoint(iftGraphicalContext *gc, iftPoint P)
{
    iftVoxel  u[8];
    int       p[8], i;
    float     dx=1.0,dy=1.0,dz=1.0,val[8];
    iftVector V[8];

    if ((int)(P.x+1.0)==gc->scene->xsize) dx = 0.0;
    if ((int)(P.y+1.0)==gc->scene->ysize) dy = 0.0;
    if ((int)(P.z+1.0)==gc->scene->zsize) dz = 0.0;

    u[0].x = (int)P.x;      u[0].y = (int)P.y;       u[0].z = (int)P.z;
    u[1].x = (int)(P.x+dx); u[1].y = (int)P.y;       u[1].z = (int)P.z;
    u[2].x = (int)P.x;      u[2].y = (int)(P.y+dy);  u[2].z = (int)P.z;
    u[3].x = (int)(P.x+dx); u[3].y = (int)(P.y+dy);  u[3].z = (int)P.z;
    u[4].x = (int)P.x;      u[4].y = (int)P.y;       u[4].z = (int)(P.z+dz);
    u[5].x = (int)(P.x+dx); u[5].y = (int)P.y;       u[5].z = (int)(P.z+dz);
    u[6].x = (int)P.x;      u[6].y = (int)(P.y+dy);  u[6].z = (int)(P.z+dz);
    u[7].x = (int)(P.x+dx); u[7].y = (int)(P.y+dy);  u[7].z = (int)(P.z+dz);

    for (i=0; i < 8; i++) {
        if (iftFValidVoxel(gc->scene,u[i])){
            p[i] = iftFGetVoxelIndex(gc->scene,u[i]);
            V[i].x = gc->phong->normal[gc->normal->val[p[i]]].x;
            V[i].y = gc->phong->normal[gc->normal->val[p[i]]].y;
            V[i].z = gc->phong->normal[gc->normal->val[p[i]]].z;
        }else{
            u[0].x = iftRound(P.x);
            u[0].y = iftRound(P.y);
            u[0].z = iftRound(P.z);
            p[0]   = iftFGetVoxelIndex(gc->scene,u[0]);
            V[0].x = gc->phong->normal[gc->normal->val[p[i]]].x;
            V[0].y = gc->phong->normal[gc->normal->val[p[i]]].y;
            V[0].z = gc->phong->normal[gc->normal->val[p[i]]].z;
            return(V[0]);
        }
    }

    val[0] =(float)V[1].x*(P.x-u[0].x)+(float)V[0].x*(u[1].x-P.x);
    val[1] =(float)V[3].x*(P.x-u[2].x)+(float)V[2].x*(u[3].x-P.x);
    val[2] =(float)V[5].x*(P.x-u[4].x)+(float)V[4].x*(u[5].x-P.x);
    val[3] =(float)V[7].x*(P.x-u[6].x)+(float)V[6].x*(u[7].x-P.x);
    val[4] = val[1]*(P.y-u[0].y) + val[0]*(u[2].y-P.y);
    val[5] = val[3]*(P.y-u[0].y) + val[2]*(u[2].y-P.y);
    V[0].x = (val[5]*(P.z-u[0].z) + val[4]*(u[4].z-P.z));

    val[0] =(float)V[1].y*(P.x-u[0].x)+(float)V[0].y*(u[1].x-P.x);
    val[1] =(float)V[3].y*(P.x-u[2].x)+(float)V[2].y*(u[3].x-P.x);
    val[2] =(float)V[5].y*(P.x-u[4].x)+(float)V[4].y*(u[5].x-P.x);
    val[3] =(float)V[7].y*(P.x-u[6].x)+(float)V[6].y*(u[7].x-P.x);
    val[4] = val[1]*(P.y-u[0].y) + val[0]*(u[2].y-P.y);
    val[5] = val[3]*(P.y-u[0].y) + val[2]*(u[2].y-P.y);
    V[0].y = (val[5]*(P.z-u[0].z) + val[4]*(u[4].z-P.z));

    val[0] =(float)V[1].z*(P.x-u[0].x)+(float)V[0].z*(u[1].x-P.x);
    val[1] =(float)V[3].z*(P.x-u[2].x)+(float)V[2].z*(u[3].x-P.x);
    val[2] =(float)V[5].z*(P.x-u[4].x)+(float)V[4].z*(u[5].x-P.x);
    val[3] =(float)V[7].z*(P.x-u[6].x)+(float)V[6].z*(u[7].x-P.x);
    val[4] = val[1]*(P.y-u[0].y) + val[0]*(u[2].y-P.y);
    val[5] = val[3]*(P.y-u[0].y) + val[2]*(u[2].y-P.y);
    V[0].z = (val[5]*(P.z-u[0].z) + val[4]*(u[4].z-P.z));

    return(V[0]);
}

void iftSurfaceShadingAlongRay(iftGraphicalContext *gc, iftPoint P0, iftPoint P1, iftPoint Pn, float *red, float *green, float *blue, iftAdjRel *A)
{
    iftPoint     u;
    int          i,k, n, p, obj_flag[gc->nobjs+1];
    float        Dx=(Pn.x-P1.x), Dy=(Pn.y-P1.y), Dz=(Pn.z-P1.z), dist;
    float        dx=0, dy=0, dz=0, phong_val, acum_opacity, opac;
    iftVoxel     v,w;
    iftVector    N;
    char         flag;
    
    *red = *green = *blue = 0.0;

    /* DDA - Digital Differential Analyzer */

    if (iftVectorsAreEqual(P1, Pn)) {
        n = 1;
    }else{ /* process points from P1 to Pn */
        if ((fabs(Dx) >= fabs(Dy))&&(fabs(Dx) >= fabs(Dz))) { /* Dx is the maximum projection of
  							     vector P1Pn */
            n  = (int)(fabs(Dx)+1);
            dx = iftSign(Dx);
            dy = dx*Dy/Dx;
            dz = dx*Dz/Dx;
        }else{
            if ((fabs(Dy) >= fabs(Dx))&&(fabs(Dy) >= fabs(Dz))) { /* Dy is the maximum projection of
  							       vector P1Pn */
                n  = (int)(fabs(Dy)+1);
                dy = iftSign(Dy);
                dx = dy*Dx/Dy;
                dz = dy*Dz/Dy;
            } else { /* Dz is the maximum projection of vector P1Pn */
                n  = (int)(fabs(Dz)+1);
                dz = iftSign(Dz);
                dx = dz*Dx/Dz;
                dy = dz*Dy/Dz;
            }
        }
    }

    /* Execute shading model along the viewing ray */

    u.x  = P1.x;  u.y = P1.y;  u.z = P1.z;

    for (k=1; k <= gc->nobjs; k++)
        obj_flag[k]=0;

    for (k=0, acum_opacity=1.0; (k < n) && (acum_opacity > IFT_EPSILON); k++) {

        v.x = iftRound(u.x);
        v.y = iftRound(u.y);
        v.z = iftRound(u.z);

        flag=0;
	for (i=0; (i < A->n)&&(flag==0); i++) { /* it avoids to skip the
						   surface of the
						   objects */
            w = iftGetAdjacentVoxel(A,v,i);
            if(iftFValidVoxel(gc->scene, w)){
                p = iftFGetVoxelIndex(gc->scene,w);
                if ((gc->border->val[p]!=0)&&(obj_flag[gc->label->val[p]]==0)){
                    flag=1;
                    if (gc->object[gc->label->val[p]].visibility != 0){
                        if (gc->object[gc->label->val[p]].opacity > IFT_EPSILON){
                            dist = iftPointDistance(u,P0);
                            N.x  = gc->phong->normal[gc->normal->val[p]].x;
                            N.y  = gc->phong->normal[gc->normal->val[p]].y;
                            N.z  = gc->phong->normal[gc->normal->val[p]].z;
                            opac = gc->object[gc->label->val[p]].opacity;
                            phong_val = opac * iftPhongShading(gc,p,N,dist) * acum_opacity;
                            *red   += phong_val * gc->object[gc->label->val[p]].red;
                            *green += phong_val * gc->object[gc->label->val[p]].green;
                            *blue  += phong_val * gc->object[gc->label->val[p]].blue;
                            acum_opacity = acum_opacity * (1.0 -  opac);
                            obj_flag[gc->label->val[p]]=1; /* it avoids multiple
					      renditions of a same
					      object in a same ray */
                        }
                    }
                }
            }
	}

        u.x = u.x + dx;
        u.y = u.y + dy;
        u.z = u.z + dz;
    }

    if (*red<0) *red=0;
    if (*green<0) *green=0;
    if (*blue<0) *blue=0;
    if (*red>1) {*red=1;}
    if (*green>1) {*green=1;}
    if (*blue>1) {*blue=1;}

}

void iftVoxelInfoByRayCasting(iftGraphicalContext *gc, iftPoint P0, iftPoint P1, iftPoint Pn, int *voxel, float *depth, int *object, iftAdjRel *A) /* only for scenes with opaque objects */
{
    iftPoint     u;
    int          i,k, n, p;
    float        Dx=(Pn.x-P1.x), Dy=(Pn.y-P1.y), Dz=(Pn.z-P1.z);
    float        dx=0, dy=0, dz=0;
    iftVoxel     v, w;
    char         flag1, flag2;

    *voxel = IFT_NIL; *depth = IFT_INFINITY_FLT; *object = IFT_NIL;

    /* DDA - Digital Differential Analyzer */

    if (iftVectorsAreEqual(P1, Pn)) {
        n = 1;
    }else{ /* process points from P1 to Pn */

        if ((fabs(Dx) >= fabs(Dy))&&(fabs(Dx) >= fabs(Dz))) { /* Dx is the maximum projection of
  							     vector P1Pn */
            n  = (int)(fabs(Dx)+1);
            dx = iftSign(Dx);
            dy = dx*Dy/Dx;
            dz = dx*Dz/Dx;
        }else{
            if ((fabs(Dy) >= fabs(Dx))&&(fabs(Dy) >= fabs(Dz))) { /* Dy is the maximum projection of
  							       vector P1Pn */
                n  = (int)(fabs(Dy)+1);
                dy = iftSign(Dy);
                dx = dy*Dx/Dy;
                dz = dy*Dz/Dy;
            } else { /* Dz is the maximum projection of vector P1Pn */
                n  = (int)(fabs(Dz)+1);
                dz = iftSign(Dz);
                dx = dz*Dx/Dz;
                dy = dz*Dy/Dz;
            }
        }
    }


    /* Find the closest voxel and its depth to the viewing plane */

    u.x  = P1.x;  u.y = P1.y;  u.z = P1.z;

    flag1=0;
    for (k=0; (k < n)&&(flag1==0); k++) {

        v.x = iftRound(u.x); v.y = iftRound(u.y); v.z = iftRound(u.z);

        flag2=0;
        for (i=0; (i < A->n)&&(flag2==0); i++) { /* it avoids to skip the
					       surface of the
					       objects */
            w = iftGetAdjacentVoxel(A,v,i);
            p = iftFGetVoxelIndex(gc->scene,w);
            if (gc->border->val[p]!=0){
                flag2=1;
                if (gc->object[gc->label->val[p]].visibility != 0){
                    if (gc->object[gc->label->val[p]].opacity > 0.0){
                        *depth  = iftPointDistance(u,P0);
                        *voxel  = p;
                        *object = gc->label->val[p];
                        flag1   = 1;
                    }
                }
            }
        }

        u.x = u.x + dx;
        u.y = u.y + dy;
        u.z = u.z + dz;
    }

}


void iftVolumeShadingAlongRay(iftGraphicalContext *gc, iftPoint P0, iftPoint P1, iftPoint Pn, float *value)
{
    iftPoint     u;
    int          k, n, p;
    float        Dx=(Pn.x-P1.x), Dy=(Pn.y-P1.y), Dz=(Pn.z-P1.z);
    float        dx=0, dy=0, dz=0, acum_opacity=1.0, dist, opac;
    iftVector    N;
    iftVoxel     v;

    *value = 0.0;

    /* DDA - Digital Differential Analyzer */

    if (iftVectorsAreEqual(P1, Pn)) {
        n = 1;
    }else{ /* process points from P1 to Pn */
        if ((fabs(Dx) >= fabs(Dy))&&(fabs(Dx) >= fabs(Dz))) { /* Dx is the maximum projection of
							     vector P1Pn */
            n  = (int)(fabs(Dx)+1);
            dx = iftSign(Dx);
            dy = dx*Dy/Dx;
            dz = dx*Dz/Dx;
        }else{
            if ((fabs(Dy) >= fabs(Dx))&&(fabs(Dy) >= fabs(Dz))) { /* Dy is the maximum projection of
							       vector P1Pn */
                n  = (int)(fabs(Dy)+1);
                dy = iftSign(Dy);
                dx = dy*Dx/Dy;
                dz = dy*Dz/Dy;
            } else { /* Dz is the maximum projection of vector P1Pn */
                n  = (int)(fabs(Dz)+1);
                dz = iftSign(Dz);
                dx = dz*Dx/Dz;
                dy = dz*Dy/Dz;
            }
        }
    }

    u.x  = P1.x;  u.y = P1.y;  u.z = P1.z;

    /* Execute Phong model along the viewing ray */

    for (k=0; (k < n)&&(acum_opacity > IFT_EPSILON); k++) {
        v.x = iftRound(u.x); v.y = iftRound(u.y); v.z = iftRound(u.z);
        p = iftFGetVoxelIndex(gc->scene,v);
        if (gc->opacity->val[p] > 0.0){
            dist = iftPointDistance(u,P0);
            N.x  = gc->phong->normal[gc->normal->val[p]].x;
            N.y  = gc->phong->normal[gc->normal->val[p]].y;
            N.z  = gc->phong->normal[gc->normal->val[p]].z;
            N    = iftNormalVectorAtPoint(gc,u);
            opac = iftOpacityValueAtPoint(gc,u);
            *value += opac * iftPhongShading(gc,p,N,dist) * acum_opacity;
            acum_opacity = acum_opacity * (1.0 -  opac);
        }
        u.x = u.x + dx;
        u.y = u.y + dy;
        u.z = u.z + dz;
    }

}


void iftVoxelSplatting(iftGraphicalContext *gc, iftImage *image, iftPoint P, int p, float phong_val, float opac)
{
    iftColor   RGB, YCbCr;
    float      acum_phong_val;
    iftVoxel   v;
    int        po;

    v.x = (int)P.x; v.y = (int)P.y; v.z = 0;
    if (iftValidVoxel(image,v)){
        po  = iftGetVoxelIndex(image,v);
        if (gc->surf_render[po].opacity > IFT_EPSILON){
            acum_phong_val = 255.0 * phong_val * gc->surf_render[po].opacity;
            RGB.val[0] = (int)( acum_phong_val * gc->object[gc->label->val[p]].red);
            RGB.val[1] = (int)( acum_phong_val * gc->object[gc->label->val[p]].green);
            RGB.val[2] = (int)( acum_phong_val * gc->object[gc->label->val[p]].blue);
            YCbCr      = iftRGBtoYCbCr(RGB,255);
            image->val[po] = YCbCr.val[0]; image->Cb[po] = YCbCr.val[1]; image->Cr[po] = YCbCr.val[2];
            gc->surf_render[po].opacity = gc->surf_render[po].opacity * (1.0 - opac);
        }
    }

    v.x = (int)(P.x+1); v.y = (int)P.y; v.z = 0;
    if (iftValidVoxel(image,v)){
        po  = iftGetVoxelIndex(image,v);
        if (gc->surf_render[po].opacity > IFT_EPSILON){
            acum_phong_val = 255.0 * phong_val * gc->surf_render[po].opacity;
            RGB.val[0] = (int)( acum_phong_val * gc->object[gc->label->val[p]].red);
            RGB.val[1] = (int)( acum_phong_val * gc->object[gc->label->val[p]].green);
            RGB.val[2] = (int)( acum_phong_val * gc->object[gc->label->val[p]].blue);
            YCbCr      = iftRGBtoYCbCr(RGB,255);
            image->val[po] = YCbCr.val[0]; image->Cb[po] = YCbCr.val[1]; image->Cr[po] = YCbCr.val[2];
            gc->surf_render[po].opacity = gc->surf_render[po].opacity * (1.0 - opac);
        }
    }

    v.x = (int)P.x; v.y = (int)(P.y+1); v.z = 0;
    if (iftValidVoxel(image,v)){
        po  = iftGetVoxelIndex(image,v);
        if (gc->surf_render[po].opacity > IFT_EPSILON){
            acum_phong_val = 255.0 * phong_val * gc->surf_render[po].opacity;
            RGB.val[0] = (int)( acum_phong_val * gc->object[gc->label->val[p]].red);
            RGB.val[1] = (int)( acum_phong_val * gc->object[gc->label->val[p]].green);
            RGB.val[2] = (int)( acum_phong_val * gc->object[gc->label->val[p]].blue);
            YCbCr      = iftRGBtoYCbCr(RGB,255);
            image->val[po] = YCbCr.val[0]; image->Cb[po] = YCbCr.val[1]; image->Cr[po] = YCbCr.val[2];
            gc->surf_render[po].opacity = gc->surf_render[po].opacity * (1.0 - opac);
        }
    }

    v.x = (int)(P.x+1); v.y = (int)(P.y+1); v.z = 0;
    if (iftValidVoxel(image,v)){
        po  = iftGetVoxelIndex(image,v);
        if (gc->surf_render[po].opacity > IFT_EPSILON){
            acum_phong_val = 255.0 * phong_val * gc->surf_render[po].opacity;
            RGB.val[0] = (int)( acum_phong_val * gc->object[gc->label->val[p]].red);
            RGB.val[1] = (int)( acum_phong_val * gc->object[gc->label->val[p]].green);
            RGB.val[2] = (int)( acum_phong_val * gc->object[gc->label->val[p]].blue);
            YCbCr      = iftRGBtoYCbCr(RGB,255);
            image->val[po] = YCbCr.val[0]; image->Cb[po] = YCbCr.val[1]; image->Cr[po] = YCbCr.val[2];
            gc->surf_render[po].opacity = gc->surf_render[po].opacity * (1.0 - opac);
        }
    }

}

iftImage *iftSurfaceRenderingByRayCasting(iftGraphicalContext *gc)
{
    iftImage  *image;
    iftVector  n;
    int        diag  = iftFDiagonalSize(gc->scene);
    iftAdjRel *A=iftSpheric(1.0);

    n.x  = 0; n.y = 0; n.z = 1;
    n    = iftTransformVector(gc->viewdir->Rinv,n);

    image =  iftCreateImage(diag,diag,1);
    iftSetCbCr(image,128);

    if (gc->normal != NULL) {

#pragma omp parallel for shared(image,A,gc,n)
        for (int po=0; po < image->n; po++) {

            iftPoint  P0,P1,Pn;
            iftColor  RGB, YCbCr;
            float     red=0.0, green=0.0, blue=0.0;

            P0.x   =  iftGetXCoord(image,po);
            P0.y   =  iftGetYCoord(image,po);
            P0.z   = -diag/2.0;
            P0     =  iftTransformPoint(gc->viewdir->Tinv,P0);

            if (iftIntersectionPoints(gc,P0,n,&P1,&Pn)){

                iftSurfaceShadingAlongRay(gc,P0,P1,Pn,&red,&green,&blue,A);

                RGB.val[0]     = (int)(255.0*red);
                RGB.val[1]     = (int)(255.0*green);
                RGB.val[2]     = (int)(255.0*blue);
                YCbCr          = iftRGBtoYCbCr(RGB,255);
                image->val[po] = YCbCr.val[0];
                image->Cb[po]  = YCbCr.val[1];
                image->Cr[po]  = YCbCr.val[2];
            }
        }
    } else { /* Estimate normals on-the-fly */

        if (iftAlmostZero(gc->overall_opac - 1.0)) {

#pragma omp parallel for shared(image,A,gc,n)
            for (int po=0; po < image->n; po++) {

                iftPoint  P0,P1,Pn;

                P0.x   =  iftGetXCoord(image,po);
                P0.y   =  iftGetYCoord(image,po);
                P0.z   = -diag/2.0;
                P0     =  iftTransformPoint(gc->viewdir->Tinv,P0);

                if (iftIntersectionPoints(gc,P0,n,&P1,&Pn)){

                    iftVoxelInfoByRayCasting(gc,P0,P1,Pn,&gc->surf_render[po].voxel,&gc->surf_render[po].depth,&gc->surf_render[po].object,A);


                }
            }

            iftRenderFromSRBuffers(gc,image);

        } else {
            iftError("Semi-transparent objects require to set object normals", "iftSurfaceRenderingByRayCasting");
        }

    }

    iftDestroyAdjRel(&A);

    return(image);
}

iftImage *iftSurfaceRenderingBySplatting(iftGraphicalContext *gc)
{
    int        xo,yo,zo,xf,yf,zf,dx,dy,dz,p,diag=iftFDiagonalSize(gc->scene);
    iftVoxel   u;
    iftPoint   P;
    iftVector  N;
    float      opac, phong_val;
    iftImage  *image = iftCreateImage(diag,diag,1);

    iftSetCbCr(image,128);

    if (gc->overall_opac != 1.0)
        iftError(
                "This function needs to be extended for semi-transparent voxel splatting by creating a buffer for each object in order to avoid multiple voxels of a same object to be splatted onto a same pixel",
                "iftSurfaceRenderingBySplatting");

    xo = gc->viewdir->ftb[gc->viewdir->octant].xo; yo = gc->viewdir->ftb[gc->viewdir->octant].yo; zo = gc->viewdir->ftb[gc->viewdir->octant].zo;
    xf = gc->viewdir->ftb[gc->viewdir->octant].xf; yf = gc->viewdir->ftb[gc->viewdir->octant].yf; zf = gc->viewdir->ftb[gc->viewdir->octant].zf;
    dx = gc->viewdir->ftb[gc->viewdir->octant].dx; dy = gc->viewdir->ftb[gc->viewdir->octant].dy; dz = gc->viewdir->ftb[gc->viewdir->octant].dz;

    switch (gc->viewdir->paxis) {

        case IFT_AXIS_X:

            for (u.x=xo; (u.x != xf); u.x = u.x + dx)
                for (u.y=yo; (u.y != yf); u.y = u.y + dy)
                    for (u.z=zo; (u.z != zf); u.z = u.z + dz){
                        p = iftFGetVoxelIndex(gc->scene,u);
                        if ((gc->border->val[p]!=0)&&
                            (gc->opacity->val[p] > 0.0)&&
                            (gc->object[gc->label->val[p]].opacity > 0.0)&&
                            (gc->object[gc->label->val[p]].visibility != 0)){
                            P.x       = u.x; P.y = u.y; P.z = u.z;
                            P         = iftTransformPoint(gc->viewdir->T,P);
                            N.x       = gc->phong->normal[gc->normal->val[p]].x;
                            N.y       = gc->phong->normal[gc->normal->val[p]].y;
                            N.z       = gc->phong->normal[gc->normal->val[p]].z;
                            N         = iftTransformVector(gc->viewdir->R,N);
                            opac      = gc->object[gc->label->val[p]].opacity;
                            phong_val = opac * iftPhongShading(gc,p,N,P.z);
                            iftVoxelSplatting(gc, image, P, p, phong_val, opac);
                        }
                    }

            break;

        case IFT_AXIS_Y:

            for (u.y=yo; (u.y != yf); u.y = u.y + dy)
                for (u.x=xo; (u.x != xf); u.x = u.x + dx)
                    for (u.z=zo; (u.z != zf); u.z = u.z + dz){
                        p = iftFGetVoxelIndex(gc->scene,u);
                        if ((gc->border->val[p]!=0)&&
                            (gc->opacity->val[p]!=0)&&
                            (gc->object[gc->label->val[p]].opacity > 0.0)&&
                            (gc->object[gc->label->val[p]].visibility != 0)){
                            P.x       = u.x; P.y = u.y; P.z = u.z;
                            P         = iftTransformPoint(gc->viewdir->T,P);
                            N.x       = gc->phong->normal[gc->normal->val[p]].x;
                            N.y       = gc->phong->normal[gc->normal->val[p]].y;
                            N.z       = gc->phong->normal[gc->normal->val[p]].z;
                            N         = iftTransformVector(gc->viewdir->R,N);
                            opac      = gc->object[gc->label->val[p]].opacity;
                            phong_val = opac * iftPhongShading(gc,p,N,P.z);
                            iftVoxelSplatting(gc, image, P, p, phong_val, opac);
                        }
                    }

            break;

        case IFT_AXIS_Z:

            for (u.z=zo; (u.z != zf); u.z = u.z + dz)
                for (u.x=xo; (u.x != xf); u.x = u.x + dx)
                    for (u.y=yo; (u.y != yf); u.y = u.y + dy){
                        p = iftFGetVoxelIndex(gc->scene,u);
                        if ((gc->border->val[p]!=0)&&
                            (gc->opacity->val[p]!=0)&&
                            (gc->object[gc->label->val[p]].opacity > 0.0)&&
                            (gc->object[gc->label->val[p]].visibility != 0)){
                            P.x       = u.x; P.y = u.y; P.z = u.z;
                            P         = iftTransformPoint(gc->viewdir->T,P);
                            N.x       = gc->phong->normal[gc->normal->val[p]].x;
                            N.y       = gc->phong->normal[gc->normal->val[p]].y;
                            N.z       = gc->phong->normal[gc->normal->val[p]].z;
                            N         = iftTransformVector(gc->viewdir->R,N);
                            opac      = gc->object[gc->label->val[p]].opacity;
                            phong_val = opac * iftPhongShading(gc,p,N,P.z);
                            iftVoxelSplatting(gc, image, P, p, phong_val, opac);
                        }
                    }
            break;

        default:
            iftError("Invalid principal axis", "iftSurfaceRenderingBySplatting");

    }


    return(image);
}

iftImage *iftVolumeRenderingByRayCasting(iftGraphicalContext *gc)
{
    iftImage  *image;
    iftVector  n;
    int        diag  = iftFDiagonalSize(gc->scene);

    image =  iftCreateImage(diag,diag,1);

    n.x  = 0; n.y = 0; n.z = 1;
    n    = iftTransformVector(gc->viewdir->Rinv,n);

    for (int po=0; po < image->n; po++) {
        iftPoint  P0,P1,Pn;
        float     value;

        P0.x   =  iftGetXCoord(image,po);
        P0.y   =  iftGetYCoord(image,po);
        P0.z   = -diag/2.0;
        P0     =  iftTransformPoint(gc->viewdir->Tinv,P0);

        if (iftIntersectionPoints(gc,P0,n,&P1,&Pn)){
            iftVolumeShadingAlongRay(gc,P0,P1,Pn,&value);
            image->val[po] = (int)(255.0*value);
        }
    }

    return(image);
}

void iftSetSceneFaces(iftGraphicalContext *gc)
{
    iftFImage *scene=gc->scene;

    gc->face[0].pos.x = scene->xsize/2.0;
    gc->face[0].pos.y = scene->ysize/2.0;
    gc->face[0].pos.z = 0;
    gc->face[1].pos.x = scene->xsize/2.0;
    gc->face[1].pos.y = scene->ysize/2.0;
    gc->face[1].pos.z = scene->zsize-1;
    gc->face[2].pos.x = scene->xsize/2.0;
    gc->face[2].pos.y = 0;
    gc->face[2].pos.z = scene->zsize/2.0;
    gc->face[3].pos.x = scene->xsize/2.0;
    gc->face[3].pos.y = scene->ysize-1;
    gc->face[3].pos.z = scene->zsize/2.0;
    gc->face[4].pos.x = 0;
    gc->face[4].pos.y = scene->ysize/2.0;
    gc->face[4].pos.z = scene->zsize/2.0;
    gc->face[5].pos.x = scene->xsize-1;
    gc->face[5].pos.y = scene->ysize/2.0;
    gc->face[5].pos.z = scene->zsize/2.0;

    gc->face[0].normal.x =  0;
    gc->face[0].normal.y =  0;
    gc->face[0].normal.z = -1;
    gc->face[1].normal.x =  0;
    gc->face[1].normal.y =  0;
    gc->face[1].normal.z =  1;
    gc->face[2].normal.x =  0;
    gc->face[2].normal.y = -1;
    gc->face[2].normal.z =  0;
    gc->face[3].normal.x =  0;
    gc->face[3].normal.y =  1;
    gc->face[3].normal.z =  0;
    gc->face[4].normal.x = -1;
    gc->face[4].normal.y =  0;
    gc->face[4].normal.z =  0;
    gc->face[5].normal.x =  1;
    gc->face[5].normal.y =  0;
    gc->face[5].normal.z =  0;

}

iftVector *iftCreateNormalTable()
{
    int        i,gamma,alpha;
    float      gamma_rad,alpha_rad;
    iftVector *normaltable;

    /* creates normal look-up table */

    normaltable = (iftVector*)iftAlloc(NUM_OF_NORMALS,sizeof(iftVector));

    if (normaltable == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateNormalTable");
    }

    normaltable[0].x = 0.0;
    normaltable[0].y = 0.0;
    normaltable[0].z = 0.0;

    i=1;
    for (gamma=-90; gamma <= 90; gamma++){
        gamma_rad = (IFT_PI * gamma) / 180.0;
        for (alpha=0; alpha < 360; alpha++){
            alpha_rad = (IFT_PI * alpha) / 180.0;
            normaltable[i].x = cos(gamma_rad)*cos(alpha_rad);
            normaltable[i].y = cos(gamma_rad)*sin(alpha_rad);
            normaltable[i].z = sin(gamma_rad);
            i++;
        }
    }

    return(normaltable);
}

int iftGetNormalIndex(iftVector normal)
{
    int gamma,alpha,index;

    if ((normal.x == 0.0) && (normal.y == 0.0) && (normal.z == 0.0)) {
        return(0);
    }

    gamma = (int)(asin(normal.z) * 180.0 / IFT_PI); /* [-90,90] */
    alpha = (int)(atan2(normal.y,normal.x) * 180.0 / IFT_PI); /* [-180,180] */
    if (alpha < 0)
        alpha += 360;
    index = ((gamma+90)*360) + alpha + 1;

    return(index);
}

iftObjectAttributes *iftCreateObjectAttributes(iftImage *label, int *nobjs)
{
    iftObjectAttributes *object;
    *nobjs = iftMaximumValue(label);

    object = (iftObjectAttributes *)iftAlloc(*nobjs+1,sizeof(iftObjectAttributes));


    /* background */

    object[0].opacity    = 0;
    object[0].red        = 0;
    object[0].green      = 0;
    object[0].blue       = 0;
    object[0].visibility = 0;

    /* default for objects */

    for (int i=1; i <= *nobjs; i++){
        object[i].opacity    = 1;
        object[i].red        = 1;
        object[i].green      = 1;
        object[i].blue       = 1;
        object[i].visibility = 1;
    }

    return(object);
}

iftSRBuffers *iftCreateSRBuffers(iftImage *label)
{
    iftSRBuffers *surf_render;
    int n = iftDiagonalSize(label);

    n = n*n;
    surf_render            = (iftSRBuffers *) iftAlloc(n,sizeof(iftSRBuffers));

    for (int p=0; p < n; p++) {
        surf_render[p].depth    = IFT_INFINITY_FLT;
        surf_render[p].opacity  = 1.0;
        surf_render[p].voxel    = IFT_NIL;
        surf_render[p].object   = IFT_NIL;
    }

    return(surf_render);
}

iftPhongModel *iftCreatePhongModel(iftFImage *scene)
{
    iftPhongModel *phong=(iftPhongModel *)iftAlloc(1,sizeof(iftPhongModel));

    phong->ka     = 0.1;
    phong->kd     = 0.7;
    phong->ks     = 0.2;
    phong->ns     = 5.0;
    phong->normal = iftCreateNormalTable();
    phong->ndists = (int)(2.0 * iftFDiagonalSize(scene) + 1);
    phong->Idist  = iftAllocFloatArray(phong->ndists);

    for(int d=0; d < phong->ndists; d++){
        phong->Idist[d] = (float)0.8*(phong->ndists - d)/(float)phong->ndists + 0.2;
    }

    return(phong);
}

void iftSetFTBaccess(iftGraphicalContext *gc)
{
    gc->viewdir->ftb    = (iftFTBaccess *)iftAlloc(8,sizeof(iftFTBaccess));

    /* Set the front-to-back voxel access from the closest to the farthest octant */

    int dx[8], dy[8], dz[8];

    dx[0] = 1; dx[1] = -1; dx[2] =  1; dx[3] = -1; dx[4] =  1; dx[5] = -1; dx[6] =  1; dx[7] = -1;
    dy[0] = 1; dy[1] =  1; dy[2] = -1; dy[3] = -1; dy[4] =  1; dy[5] =  1; dy[6] = -1; dy[7] = -1;
    dz[0] = 1; dz[1] =  1; dz[2] =  1; dz[3] =  1; dz[4] = -1; dz[5] = -1; dz[6] = -1; dz[7] = -1;

    int xo[8], yo[8], zo[8];

    xo[0] = 0; xo[1] = gc->scene->xsize-1; xo[2] = 0;                  xo[3] = gc->scene->xsize-1; xo[4] = 0;                  xo[5] = gc->scene->xsize-1; xo[6] = 0;                  xo[7] = gc->scene->xsize-1;
    yo[0] = 0; yo[1] = 0;                  yo[2] = gc->scene->ysize-1; yo[3] = gc->scene->ysize-1; yo[4] = 0;                  yo[5] = 0;                  yo[6] = gc->scene->ysize-1; yo[7] = gc->scene->ysize-1;
    zo[0] = 0; zo[1] = 0;                  zo[2] = 0;                  zo[3] = 0;                  zo[4] = gc->scene->zsize-1; zo[5] = gc->scene->zsize-1; zo[6] = gc->scene->zsize-1; zo[7] = gc->scene->zsize-1;

    int xf[8], yf[8], zf[8];

    for (int i=0; i < 8; i++) {
        if (dx[i]==1)
            xf[i]=gc->scene->xsize;
        else
            xf[i]=-1;
        if (dy[i]==1)
            yf[i]=gc->scene->ysize;
        else
            yf[i]=-1;
        if (dz[i]==1)
            zf[i]=gc->scene->zsize;
        else
            zf[i]=-1;
    }

    for (int i=0; i < 8; i++) {
        gc->viewdir->ftb[i].dx = dx[i];    gc->viewdir->ftb[i].dy = dy[i];    gc->viewdir->ftb[i].dz = dz[i];
        gc->viewdir->ftb[i].xo = xo[i];    gc->viewdir->ftb[i].yo = yo[i];    gc->viewdir->ftb[i].zo = zo[i];
        gc->viewdir->ftb[i].xf = xf[i];    gc->viewdir->ftb[i].yf = yf[i];    gc->viewdir->ftb[i].zf = zf[i];
    }

}

void iftDrawBoundingBoxBorders2DInPlace(iftImage *img, iftBoundingBox bb, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint) {
    int Imax = iftNormalizationValue(iftMaximumValue(img));
    // draws x-lines
    //#pragma omp parallel for
    for (int y = bb.begin.y; y <= bb.end.y; y += (bb.end.y-bb.begin.y))
        for (int x = bb.begin.x; x <= bb.end.x; x++) {
            iftVoxel u = {.x = x, .y = y, .z = 0};
            iftDrawPoint(img, u, YCbCr, A, Imax);
        }

    // draws y-lines
    //#pragma omp parallel for
    for (int x = bb.begin.x; x <= bb.end.x; x += (bb.end.x-bb.begin.x))
        for (int y = bb.begin.y; y <= bb.end.y; y++) {
            iftVoxel u = {.x = x, .y = y, .z = 0};
            iftDrawPoint(img, u, YCbCr, A, Imax);
        }
    
    // draws the central point
    if(drawCentralPoint)
        iftDrawPoint(img, iftBoundingBoxCenterVoxel(bb), YCbCr, A, Imax);
}


void iftDrawBoundingBoxBorders3DInPlace(iftImage *img, iftBoundingBox bb, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint) {
    int Imax = iftNormalizationValue(iftMaximumValue(img));
    // draws x-lines
#pragma omp parallel for
    for (int z = bb.begin.z; z <= bb.end.z; z += (bb.end.z-bb.begin.z))
        for (int y = bb.begin.y; y <= bb.end.y; y += (bb.end.y-bb.begin.y))
            for (int x = bb.begin.x; x <= bb.end.x; x++) {
                iftVoxel u = {.x = x, .y = y, .z = z};
                iftDrawPoint(img, u, YCbCr, A, Imax);
            }

    // draws y-lines
#pragma omp parallel for
    for (int z = bb.begin.z; z <= bb.end.z; z += (bb.end.z-bb.begin.z))
        for (int x = bb.begin.x; x <= bb.end.x; x += (bb.end.x-bb.begin.x))
            for (int y = bb.begin.y; y <= bb.end.y; y++) {
                iftVoxel u = {.x = x, .y = y, .z = z};
                iftDrawPoint(img, u, YCbCr, A, Imax);
            }

    // draws z-lines
#pragma omp parallel for
    for (int y = bb.begin.y; y <= bb.end.y; y += (bb.end.y-bb.begin.y))
        for (int x = bb.begin.x; x <= bb.end.x; x += (bb.end.x-bb.begin.x))
            for (int z = bb.begin.z; z <= bb.end.z; z++) {
                iftVoxel u = {.x = x, .y = y, .z = z};
                iftDrawPoint(img, u, YCbCr, A, Imax);
            }
    
    // draws the central point
    if(drawCentralPoint)
        iftDrawPoint(img, iftBoundingBoxCenterVoxel(bb), YCbCr, A, Imax);
}


/* ------------------------------- Public Methods --------------------------*/

void iftSetSceneNormal(iftGraphicalContext *gc)
{
    iftAdjRel *A   = iftSpheric(3.0);
    float     *mag = iftAllocFloatArray(A->n);
    float      Delta;
    int        i,p,q;
    iftVoxel   u,v;
    iftVector  N;


    if (gc->opacity == NULL)
        iftError("Set scene opacity first", "iftSetSceneNormal");

    if (gc->normal != NULL)
        iftDestroyImage(&gc->normal);

    gc->normal = iftCreateImage(gc->scene->xsize,gc->scene->ysize,gc->scene->zsize);

    for (i=0; i < A->n; i++)
        mag[i]=sqrtf(A->dx[i]*A->dx[i]+A->dy[i]*A->dy[i]+A->dz[i]*A->dz[i]);

    for (u.z=0; u.z < gc->scene->zsize; u.z++)
        for (u.y=0; u.y < gc->scene->ysize; u.y++)
            for (u.x=0; u.x < gc->scene->xsize; u.x++) {
                p = iftFGetVoxelIndex(gc->scene,u);
                if (gc->opacity->val[p]>0.0) {
                    N.x = N.y = N.z = 0.0;
                    for (i=1; i < A->n; i++) {
                        v = iftGetAdjacentVoxel(A,u,i);
                        if (iftFValidVoxel(gc->scene,v)){
                            q = iftFGetVoxelIndex(gc->scene,v);
                            Delta = gc->scene->val[q]-gc->scene->val[p];
                            N.x  += Delta*A->dx[i]/mag[i];
                            N.y  += Delta*A->dy[i]/mag[i];
                            N.z  += Delta*A->dz[i]/mag[i];
                        }
                    }
                    /* it assumes scenes with objects brighter than the background */
                    N = iftNormalizeVector(N);
                    N.x = -N.x; N.y = -N.y; N.z = -N.z;
                    gc->normal->val[p] = iftGetNormalIndex(N);
                }
            }


    iftFree(mag);
    iftDestroyAdjRel(&A);

}

void iftSetProjectionMode(iftGraphicalContext *gc, char proj_mode)
{
    gc->proj_mode = proj_mode;
    gc->viewdir->V.x = 0; gc->viewdir->V.y = 0; gc->viewdir->V.z = -1;
    if (gc->proj_mode == RAYCASTING){
        gc->viewdir->V   = iftTransformVector(gc->viewdir->Rinv,gc->viewdir->V);
    }
}

void iftSetObjectNormal(iftGraphicalContext *gc, char normal_type)
{
    iftAdjRel *A = NULL;
    iftImage  *dist = NULL;
    int        i,p,q;
    iftVoxel   u,v;
    iftVector  N;
    iftSet    *S=NULL;

    if (gc->label == NULL)
        iftError("Object labels are required", "iftSetObjectNormal");

    if (gc->normal != NULL)
        iftDestroyImage(&gc->normal);

    if (gc->opacity != NULL)
        iftDestroyFImage(&gc->opacity);

    A             = iftSpheric(1.0);
    S             = iftObjectBorderSet(gc->label,A);
    iftDestroyAdjRel(&A);

    /* estimate object-based normal vectors and set opacity scene for the shell */

    gc->normal  = iftCreateImage(gc->label->xsize,gc->label->ysize,gc->label->zsize);
    gc->opacity = iftCreateFImage(gc->label->xsize,gc->label->ysize,gc->label->zsize);

    if (normal_type == OBJECT_NORMAL){
        /* extract shell inside the object and estimate normals for the
         border voxels */

        A             = iftSpheric(sqrtf(3.0));
        dist          = iftShellDistTrans(gc->label,A,IFT_INTERIOR,5);
        iftDestroyAdjRel(&A);

        A   = iftSpheric(5.0);

        while (S != NULL) {

            p = iftRemoveSet(&S);

            gc->opacity->val[p] = 1.0;

            u = iftGetVoxelCoord(dist,p);
            N.x = N.y = N.z = 0.0;

            for (i=1; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A,u,i);
                if (iftValidVoxel(dist,v)){
                    q = iftGetVoxelIndex(dist,v);
                    if ((gc->label->val[p]==gc->label->val[q])&&
                        (dist->val[q] != IFT_INFINITY_INT)){
                        N.x  += dist->val[q]*A->dx[i];
                        N.y  += dist->val[q]*A->dy[i];
                        N.z  += dist->val[q]*A->dz[i];
                    }
                }
            }

            /* force normal to point outward the object */

            N   = iftNormalizeVector(N);
            N.x = -N.x; N.y = -N.y; N.z = -N.z;
            gc->normal->val[p] = iftGetNormalIndex(N);
        }

    } else { /* SCENE_NORMAL */
        A   = iftSpheric(5.0);

        while (S != NULL) {

            p = iftRemoveSet(&S);

            gc->opacity->val[p] = 1.0;

            u = iftFGetVoxelCoord(gc->scene,p);

            N.x = N.y = N.z = 0.0;
            for (i=1; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A,u,i);
                if (iftFValidVoxel(gc->scene,v)){
                    q = iftFGetVoxelIndex(gc->scene,v);
                    if (gc->label->val[q]==gc->label->val[p]){
                        N.x  += gc->scene->val[q]*A->dx[i];
                        N.y  += gc->scene->val[q]*A->dy[i];
                        N.z  += gc->scene->val[q]*A->dz[i];
                    }
                }
            }

            /* it assumes scenes with objects brighter than the
           background, THEN forces normal to point outward the
           object */

            N = iftNormalizeVector(N);
            N.x = -N.x; N.y = -N.y; N.z = -N.z;

            gc->normal->val[p] = iftGetNormalIndex(N);
        }
    }


    iftDestroyAdjRel(&A);
    iftDestroyImage(&dist);
}

void iftSetObjectColor(iftGraphicalContext *gc, int object, float red, float green, float blue)
{

    if ((object > gc->nobjs)||(object <= 0))
        iftError("Invalid object", "iftSetObjectColor");

    if ((red < 0)||(red > 1)||(green < 0)||(green > 1)||(blue < 0)||(blue > 1))
        iftError("Invalid color value", "iftSetObjectColor");

    gc->object[object].red   = red;
    gc->object[object].green = green;
    gc->object[object].blue  = blue;

}

void iftSetObjectOpacity(iftGraphicalContext *gc, int object, float opacity)
{

    if ((object > gc->nobjs)||(object <= 0))
        iftError("Invalid object", "iftSetObjectOpacity");

    if ((opacity < 0)||(opacity > 1))
        iftError("Invalid opacity value", "iftSetObjectOpacity");

    gc->object[object].opacity   = opacity;

    gc->overall_opac=1.0;

    for (int i=1; i <= gc->nobjs; i++)
        gc->overall_opac *=  gc->object[i].opacity;

}

void iftSetObjectVisibility(iftGraphicalContext *gc, int object, char visibility)
{

    if ((object > gc->nobjs)||(object <= 0))
        iftError("Invalid object", "iftSetObjectVisibility");

    if ((visibility != 0)&&(visibility != 1))
        iftError("Invalid visibility value", "iftSetObjectVisibility");

    gc->object[object].visibility = visibility;

}

void iftSetViewDir(iftGraphicalContext *gc, float tilt, float spin)
{

    iftVector  t;
    iftMatrix *Rx,*Ry,*Txyz,*Tuv,*aux;
    int        diag=iftFDiagonalSize(gc->scene);

    if (gc->viewdir == NULL) {
        gc->viewdir = (iftViewDir *)iftAlloc(1,sizeof(iftViewDir));
        iftSetFTBaccess(gc);
    }else{
        iftDestroyMatrix(&gc->viewdir->T);
        iftDestroyMatrix(&gc->viewdir->Tinv);
        iftDestroyMatrix(&gc->viewdir->R);
        iftDestroyMatrix(&gc->viewdir->Rinv);
    }

    /* Set scene transformation and rotation matrices */

    t.x = -gc->scene->xsize/2.0; t.y = -gc->scene->ysize/2.0; t.z = -gc->scene->zsize/2.0;
    Txyz            =  iftTranslationMatrix(t);
    Rx              =  iftRotationMatrix(IFT_AXIS_X, tilt);
    Ry              =  iftRotationMatrix(IFT_AXIS_Y, spin);
    gc->viewdir->R  =  iftMultMatrices(Ry,Rx);
    t.x =  diag/2.0; t.y = diag/2.0; t.z = diag/2.0;
    Tuv             =  iftTranslationMatrix(t);

    aux             =  iftMultMatrices(gc->viewdir->R,Txyz);
    gc->viewdir->T  =  iftMultMatrices(Tuv,aux);

    iftDestroyMatrix(&aux);
    iftDestroyMatrix(&Txyz);
    iftDestroyMatrix(&Tuv);
    iftDestroyMatrix(&Rx);
    iftDestroyMatrix(&Ry);

    /* Set scene transformation and rotation inverse matrices */

    t.x      = -diag/2.0; t.y = -diag/2.0; t.z = -diag/2.0;
    Tuv      =  iftTranslationMatrix(t);
    t.x      =  gc->scene->xsize/2.0; t.y = gc->scene->ysize/2.0; t.z = gc->scene->zsize/2.0;
    Txyz     =  iftTranslationMatrix(t);
    Rx       =  iftRotationMatrix(IFT_AXIS_X, -tilt);
    Ry       =  iftRotationMatrix(IFT_AXIS_Y, -spin);
    gc->viewdir->Rinv  =  iftMultMatrices(Rx,Ry);

    aux               =  iftMultMatrices(gc->viewdir->Rinv,Tuv);
    gc->viewdir->Tinv =  iftMultMatrices(Txyz,aux);

    iftDestroyMatrix(&aux);
    iftDestroyMatrix(&Txyz);
    iftDestroyMatrix(&Tuv);
    iftDestroyMatrix(&Rx);
    iftDestroyMatrix(&Ry);

    /* Set viewing direction */

    gc->viewdir->V.x = 0; gc->viewdir->V.y = 0; gc->viewdir->V.z = -1;
    if (gc->proj_mode == RAYCASTING){
        gc->viewdir->V   = iftTransformVector(gc->viewdir->Rinv,gc->viewdir->V);
    }

    /* Set principal axis */

    float max_proj = IFT_INFINITY_FLT_NEG;
    iftVector N;

    N.x = 0; N.y = 0; N.z = 1;
    N = iftTransformVector(gc->viewdir->Rinv,N);
    if (fabs(N.x) > max_proj){
        gc->viewdir->paxis  = IFT_AXIS_X;
        max_proj            = fabs(N.x);
    }
    if (fabs(N.y) > max_proj){
        gc->viewdir->paxis  = IFT_AXIS_Y;
        max_proj            = fabs(N.y);
    }
    if (fabs(N.z) > max_proj){
        gc->viewdir->paxis  = IFT_AXIS_Z;
        max_proj            = fabs(N.z);
    }

    /* Set closest octant */

    float depth= IFT_INFINITY_FLT;

    for (int i=0; i < 8; i++) {
        iftPoint P;

        P.x = gc->viewdir->ftb[i].xo; P.y = gc->viewdir->ftb[i].yo; P.z = gc->viewdir->ftb[i].zo;
        P   = iftTransformPoint(gc->viewdir->T,P);
        if (P.z < depth) {
            gc->viewdir->octant = i;
            depth = P.z;
        }
    }

}

void iftSetSceneOpacity(iftGraphicalContext *gc, float min_val, float max_val, iftImage *grad, int grad_thres, float max_opac)
{
    if (gc->opacity == NULL)
        gc->opacity = iftCreateFImage(gc->scene->xsize,gc->scene->ysize,gc->scene->zsize);

    int Gt = grad_thres; 
    if (Gt < 0){
      Gt = iftOtsu(grad);
      Gt = Gt - 0.1*Gt;
    }
    int Gmax = iftMaximumValue(grad);
    float K = Gmax-Gt;
    
    for (int p = 0; p < grad->n; p++) {
      if ((gc->scene->val[p] >= min_val)&&(gc->scene->val[p] <= max_val)){
	float x = (float)grad->val[p]-(float)Gt;
	  if (x > 0)
	    gc->opacity->val[p] = max_opac * (1.0 - x/K);
      }
    }
}

void iftSetAmbientReflection(iftGraphicalContext *gc, float ka)
{
    if ((ka <= 1.0) && (ka >= 0.0))
        gc->phong->ka = ka;
}

void iftSetDiffuseReflection(iftGraphicalContext *gc, float kd)
{
    if ((kd <= 1.0) && (kd >= 0.0))
        gc->phong->ka = kd;
}

void iftSetSpecularReflection(iftGraphicalContext *gc, float ks)
{
    if ((ks <= 1.0) && (ks >= 0.0))
        gc->phong->ka = ks;
}

char iftIntersectionPoints(iftGraphicalContext *gc, iftPoint P0, iftVector n, iftPoint *P1, iftPoint *Pn)
{
    iftVector V;
    iftPoint  P;
    int i;
    float lambda, a, b, lambda_1= IFT_INFINITY_FLT, lambda_n= IFT_INFINITY_FLT_NEG;

    P1->x=Pn->x=0;
    P1->y=Pn->y=0;
    P1->z=Pn->z=0;

    for (i=0; i < 6; i++) {
        b = iftVectorInnerProduct(gc->face[i].normal, n);
        if (fabs(b) > IFT_EPSILON){
            V.x    = gc->face[i].pos.x - P0.x;
            V.y    = gc->face[i].pos.y - P0.y;
            V.z    = gc->face[i].pos.z - P0.z;
            a      = iftVectorInnerProduct(gc->face[i].normal, V);
            lambda = a/b;

            P.x = iftRound(P0.x + lambda * n.x);
            P.y = iftRound(P0.y + lambda * n.y);
            P.z = iftRound(P0.z + lambda * n.z);

            if (iftFValidPoint(gc->scene,P)){
                if (lambda < lambda_1) {
                    lambda_1 = lambda;
                }
                if (lambda > lambda_n) {
                    lambda_n = lambda;
                }
            }
        }
    }

    if (lambda_1 < lambda_n){

        P1->x = iftRound(P0.x + lambda_1 * n.x);
        P1->y = iftRound(P0.y + lambda_1 * n.y);
        P1->z = iftRound(P0.z + lambda_1 * n.z);

        Pn->x = iftRound(P0.x + lambda_n * n.x);
        Pn->y = iftRound(P0.y + lambda_n * n.y);
        Pn->z = iftRound(P0.z + lambda_n * n.z);

        return(1);
    }else{
        return(0);
    }
}

iftGraphicalContext *iftCreateGraphicalContext(iftFImage *scene, iftImage *label)
{
    iftGraphicalContext *gc;

    gc = (iftGraphicalContext *) iftAlloc(1, sizeof(iftGraphicalContext));

    gc->scene          = iftFCopyImage(scene);
    gc->phong          = iftCreatePhongModel(scene);
    gc->viewdir        = NULL;
    gc->proj_mode      = RAYCASTING;
    gc->nobjs          = 0;
    gc->overall_opac   = 1.0;
    iftSetViewDir(gc,0,0);
    iftSetSceneFaces(gc);

    if (label != NULL){ /* for surface rendering */
        gc->label       = iftCopyImage(label);
        iftAdjRel *A    = iftSpheric(1.0);
        gc->border      = iftObjectBorders(label,A, true, true);
        iftDestroyAdjRel(&A);
        gc->object      = iftCreateObjectAttributes(label,&gc->nobjs);
        gc->surf_render = iftCreateSRBuffers(label);
    }else{ /* for volume rendering */
        gc->label       = NULL;
        gc->border      = NULL;
        gc->object      = NULL;
        gc->surf_render = NULL;
    }

    gc->opacity       = NULL;
    gc->normal        = NULL;

    return(gc);
}

void iftDestroyGraphicalContext(iftGraphicalContext *gc)
{
    if (gc != NULL) {
        if (gc->nobjs != 0){
            iftFree(gc->object);
            iftFree(gc->surf_render);
            iftDestroyImage(&gc->label);
            iftDestroyImage(&gc->border);
        }
        if (gc->opacity != NULL)
            iftDestroyFImage(&gc->opacity);
        if (gc->normal != NULL)
            iftDestroyImage(&gc->normal);
        iftFree(gc->phong->normal);
        iftFree(gc->phong->Idist);
        iftFree(gc->phong);
        iftDestroyMatrix(&gc->viewdir->T);
        iftDestroyMatrix(&gc->viewdir->Tinv);
        iftDestroyMatrix(&gc->viewdir->R);
        iftDestroyMatrix(&gc->viewdir->Rinv);
        iftFree(gc->viewdir->ftb);
        iftFree(gc->viewdir);
        iftDestroyFImage(&gc->scene);
        iftFree(gc);
    }else{
        iftWarning("Graphical context is already NULL","iftDestroyGraphicalContext");
    }
}

iftImage *iftSurfaceRender(iftGraphicalContext *gc)
{
    iftImage  *image=NULL;

    if (gc->nobjs==0)
        iftError("There are no objects for visualization", "iftSurfaceRender");

    if (gc->normal == NULL){ /* normals estimated on the fly */
        iftResetSRBuffers(gc);
        if (gc->proj_mode == SPLATTING) {
            iftWarning("Changing projection mode to ray casting","iftSurfaceRender");
            iftSetProjectionMode(gc,RAYCASTING);
        }
        image = iftSurfaceRenderingByRayCasting(gc);
    }else{ /* pre-computed normals */
        if (gc->proj_mode == SPLATTING)
            image = iftSurfaceRenderingBySplatting(gc);
        else
            image = iftSurfaceRenderingByRayCasting(gc);
    }

    return(image);
}

iftImage *iftVolumeRender(iftGraphicalContext *gc)
{
    iftImage  *image=NULL;

    if (gc->opacity == NULL) {
        iftError("Set opacity scene", "iftVolumeRenderingByRayCasting");
    }

    if (gc->normal == NULL) {
        iftSetSceneNormal(gc);
    }

    if (gc->proj_mode == SPLATTING) {
        iftWarning("Changing projection mode to ray casting","iftVolumeRender");
        iftSetProjectionMode(gc,RAYCASTING);
    }

    image = iftVolumeRenderingByRayCasting(gc);

    return(image);
}

iftImage *iftGraphicImage(iftImage *img)
{
    iftImage *gimg;
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));


    gimg = iftNormalize(img,0,normalization_value);

    if (!iftIsColorImage(gimg)){
        gimg->Cb = iftAllocUShortArray(gimg->n);
        gimg->Cr = iftAllocUShortArray(gimg->n);
        for (int q=0; q < gimg->n; q++) {
            gimg->Cb[q] = normalization_value/2;
            gimg->Cr[q] = normalization_value/2;
        }
    }

    return(gimg);
}

iftImage *iftDrawPath(iftImage *img, iftImage *pred, int last, iftFColor normRGB, iftAdjRel *B)
{
    int q, Vmax;
    int maxRangevalue = iftNormalizationValue(iftMaximumValue(img));
    iftVoxel u;
    iftImage *gimg = iftGraphicImage(img);
    iftColor  YCbCr, RGB;

    Vmax       = iftMaximumValue(gimg);

    RGB.val[0] = Vmax*normRGB.val[0];
    RGB.val[1] = Vmax*normRGB.val[1];
    RGB.val[2] = Vmax*normRGB.val[2];

    YCbCr      = iftRGBtoYCbCr(RGB,Vmax);

    q = last;
    while (q != IFT_NIL){
        u = iftGetVoxelCoord(gimg,q);
        iftDrawPoint(gimg, u, YCbCr, B, maxRangevalue);
        q = pred->val[q];
    }

    return(gimg);
}

void iftDrawPoint(iftImage *img, iftVoxel u, iftColor YCbCr, iftAdjRel *B, int normValue)
{
    int q,i;
    iftVoxel v;

    if (!iftValidVoxel(img,u))
        iftError("Point is outside the image domain", "iftDrawPoint");

    if (!iftIsColorImage(img)){
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
        for (q=0; q < img->n; q++) {
            img->Cb[q] = normValue/2;
            img->Cr[q] = normValue/2;
        }
    }

    for (i=0; i < B->n; i++) {
        v.x = u.x + B->dx[i];
        v.y = u.y + B->dy[i];
        v.z = u.z + B->dz[i];
        if (iftValidVoxel(img,v)){
            q = iftGetVoxelIndex(img,v);
            img->val[q]=YCbCr.val[0];
            img->Cb[q]=(ushort) YCbCr.val[1];
            img->Cr[q]=(ushort) YCbCr.val[2];
        }
    }

}

void iftDrawPointAlpha(iftImage *img, iftVoxel u, iftColor YCbCr, iftAdjRel *B, int normValue, float alpha)
{
    int q,i;
    iftVoxel v;

    if (!iftValidVoxel(img,u))
        iftError("Point is outside the image domain", "iftDrawPointAlpha");

    if (!iftIsColorImage(img)){
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
        for (q=0; q < img->n; q++) {
            img->Cb[q] = normValue/2;
            img->Cr[q] = normValue/2;
        }
    }

    iftColor RGB = iftYCbCrtoRGB(YCbCr, normValue);

    for (i=0; i < B->n; i++) {
        v.x = u.x + B->dx[i];
        v.y = u.y + B->dy[i];
        v.z = u.z + B->dz[i];
        if (iftValidVoxel(img,v)){
            q = iftGetVoxelIndex(img,v);

            iftColor color = iftGetRGB(img, q, normValue);
            color.val[0] = color.val[0]*(1.0-alpha) + RGB.val[0]*alpha;
            color.val[1] = color.val[1]*(1.0-alpha) + RGB.val[1]*alpha;
            color.val[2] = color.val[2]*(1.0-alpha) + RGB.val[2]*alpha;

            iftSetRGB(img, q, color.val[0], color.val[1], color.val[2], normValue);
        }
    }
}

void iftDrawLine(iftImage *img, iftVoxel u1, iftVoxel u2, iftColor YCbCr, iftAdjRel *B)
{
    iftPoint u;
    iftVoxel v;
    int p, k, n;
    float dx=0, dy=0, dz=0, Dx=u2.x-u1.x, Dy=u2.y-u1.y, Dz=u2.z-u1.z;
    int maxRangeValue = iftNormalizationValue(iftMaximumValue(img));

    if ( (!iftValidVoxel(img,u1)) || (!iftValidVoxel(img,u2)) )
        iftError("Line has end point(s) outside the image domain", "iftDrawLine");

    if (img->Cb == NULL){
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
        for (p=0; p < img->n; p++) {
            img->Cb[p] = maxRangeValue/2;
            img->Cr[p] = maxRangeValue/2;
        }
    }

    /* DDA - Digital Differential Analyzer */

    if (iftVoxelsAreEqual(u1,u2)) {
        n = 1;
    }else{ /* draw line from u1 to u2 */
        if ((fabs(Dx) >= fabs(Dy))&&(fabs(Dx) >= fabs(Dz))) { /* Dx is the maximum projection of
				     vector u1u2 */
            n  = (int)(fabs(Dx)+1);
            dx = iftSign(Dx);
            dy = (float)dx*Dy/(float)Dx;
            dz = (float)dx*Dz/(float)Dx;
        }else{
            if ((fabs(Dy) >= fabs(Dx))&&(fabs(Dy) >= fabs(Dz))) { /* Dy is the maximum projection of
				       vector u1u2 */
                n  = (int)(fabs(Dy)+1);
                dy = iftSign(Dy);
                dx = (float)dy*Dx/(float)Dy;
                dz = (float)dy*Dz/(float)Dy;
            } else { /* Dz is the maximum projection of vector u1u2 */
                n  = (int)(fabs(Dz)+1);
                dz = iftSign(Dz);
                dx = (float)dz*Dx/(float)Dz;
                dy = (float)dz*Dy/(float)Dz;
            }
        }
    }

    u.x = u1.x;  u.y = u1.y;  u.z = u1.z;
    for (k=0; k < n; k++) {
        v.x = iftRound(u.x); v.y = iftRound(u.y); v.z = iftRound(u.z);
        iftDrawPoint(img, v, YCbCr, B, maxRangeValue);
        u.x = (u.x + dx);
        u.y = (u.y + dy);
        u.z = (u.z + dz);
    }
}

void iftDrawPoints(iftImage *img, iftSet *S, iftColor YCbCr, iftAdjRel *B)
{
    iftSet *Saux=NULL;
    iftVoxel u;
    int p;
    int maxRangeValue = iftNormalizationValue(iftMaximumValue(img));

    if (img->Cb == NULL){
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
        for (p=0; p < img->n; p++) {
            img->Cb[p] = maxRangeValue/2;
            img->Cr[p] = maxRangeValue/2;
        }
    }

    Saux = S;
    while (Saux!=NULL) {
        p   = Saux->elem;
        u.x = iftGetXCoord(img,p);
        u.y = iftGetYCoord(img,p);
        u.z = iftGetZCoord(img,p);
        iftDrawPoint(img, u, YCbCr, B, maxRangeValue);
        Saux = Saux->next;
    }

}

void iftDrawCurve(iftImage *img, iftCurve *curve, iftColor YCbCr, iftAdjRel *B)
{
    iftVoxel u;
    int p;
    int maxRangeValue = iftNormalizationValue(iftMaximumValue(img));

    if (img->Cb == NULL){
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
        for (p=0; p < img->n; p++) {
            img->Cb[p] = maxRangeValue/2;
            img->Cr[p] = maxRangeValue/2;
        }
    }

    for (int i=0; i < curve->npts; i++) {
        u.x = (int) curve->pt[i].x;
        u.y = (int) curve->pt[i].y;
        u.z = (int) curve->pt[i].z;
        iftDrawPoint(img, u, YCbCr, B, maxRangeValue);
    }
}

void iftDrawBorders(iftImage *img, iftImage *label, iftAdjRel *A, iftColor YCbCr, iftAdjRel *B)
{
    iftVoxel u,v;
    int i,p,q;
    int maxRangeValue = iftNormalizationValue(iftMaximumValue(img));

    if ((img->xsize != label->xsize)||
        (img->ysize != label->ysize)||
        (img->zsize != label->zsize))
        iftError("Images must have the same domain", "iftDrawBorders");


    if (!iftIsColorImage(img))
        iftSetCbCr(img,maxRangeValue/2);

    if (A->n > 1){
        for (p=0; p < img->n; p++) {
            u = iftGetVoxelCoord(label,p);
            for (i=0; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A,u,i);
                if (iftValidVoxel(label,v)){
                    q = iftGetVoxelIndex(label,v);
                    if (label->val[p] != label->val[q]){
                        iftDrawPoint(img, u, YCbCr, B, maxRangeValue);
                        break;
                    }
                }
            }
        }
    } else {
        for (p=0; p < img->n; p++)
            if (label->val[p] != 0){
                u = iftGetVoxelCoord(label,p);
                iftDrawPoint(img, u, YCbCr, B, maxRangeValue);
            }
    }
}

void iftDrawObject(iftImage *img, iftImage *bin, iftColor YCbCr, iftAdjRel *B)
{
    int p;
    iftVoxel u;
    int maxRangeValue = iftNormalizationValue(iftMaximumValue(img));

    if ((img->xsize != bin->xsize)||
        (img->ysize != bin->ysize)||
        (img->zsize != bin->zsize))
        iftError("Images must have the same domain", "iftDrawObject");

    if(!iftIsColorImage(img))
        iftSetCbCr(img,maxRangeValue/2);

    for (p=0; p < bin->n; p++) {
        if (bin->val[p]!=0){
            u.x = iftGetXCoord(bin,p);
            u.y = iftGetYCoord(bin,p);
            u.z = iftGetZCoord(bin,p);
            iftDrawPoint(img, u, YCbCr, B, maxRangeValue);
        }
    }

}

void iftDrawRoots(iftImage *img, iftImage *root, iftColor YCbCr, iftAdjRel *B)
{
    int p;
    iftVoxel u;
    int maxRangeValue = iftNormalizationValue(iftMaximumValue(img));

    iftVerifyImageDomains(img,root,"iftDrawObject");

    iftSetCbCr(img,maxRangeValue/2);

    for (p=0; p < root->n; p++) {
        if (root->val[p]==p){
            u.x = iftGetXCoord(root,p);
            u.y = iftGetYCoord(root,p);
            u.z = iftGetZCoord(root,p);
            iftDrawPoint(img, u, YCbCr, B, maxRangeValue);
        }
    }
}

void iftDrawLabels(iftImage *img, const iftImage *label_img, const iftColorTable *cmap,
                   const iftAdjRel *A, bool fill, float alpha) {
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));
    iftVerifyImageDomains(img, label_img, "iftDrawLabels");
    
    if (cmap == NULL)
        iftError("Color Table is NULL", "iftDrawLabels");
    
    int max_label = iftMaximumValue(label_img);
    if (max_label >= cmap->ncolors)
        iftError("There is no color for the label %d\nNumber of color = %d\n", "iftDrawLabels", max_label,
                 cmap->ncolors);
    
    if (img->Cb == NULL)
        iftSetCbCr(img,(iftMaxImageRange(iftImageDepth(img))+1)/2);

    for (int p = 0; p < img->n; p++) {
        iftVoxel v = iftGetVoxelCoord(img, p);

        if (label_img->val[p] > 0) {
            iftColor color_lb = cmap->color[label_img->val[p]];

            color_lb = iftYCbCrtoRGB(color_lb, normalization_value);

            if (fill) {
                iftColor color = iftGetRGB(img, p, normalization_value);
                color.val[0] = (color.val[0] * (1.0 - alpha)) + (color_lb.val[0] * alpha);
                color.val[1] = (color.val[1] * (1.0 - alpha)) + (color_lb.val[1] * alpha);
                color.val[2] = (color.val[2] * (1.0 - alpha)) + (color_lb.val[2] * alpha);

                iftSetRGB(img, p, color.val[0], color.val[1], color.val[2], normalization_value);
            }
            
            if (A) {
                for (int i = 1; i < A->n; i++) {
                    iftVoxel u = iftGetAdjacentVoxel(A, v, i);
    
                    if (iftValidVoxel(img, u) && label_img->val[p] != label_img->val[iftGetVoxelIndex(img, u)]) {
                        iftSetRGB(img, p, color_lb.val[0], color_lb.val[1], color_lb.val[2], normalization_value);
                        break;
                    }
                }
            }
        }
    }
}


iftImage *iftDrawVoxelSamples(iftDataSet *Z, iftFeatureLabel opt, char *ref_data_type) {
    iftWarning("The color table scheme was altered. Check if this function keeps working. Afterwards, remove this warning.", "iftDrawVoxelSamples");

    iftImage      *img = NULL;
    int           s, p;
    iftVoxel      u;
    iftAdjRel     *B;
    iftColor      RGB, YCbCr;
    iftColorTable *ctb = NULL;
    int           max_val;

    if ((Z->ref_data == NULL) ||
        ((strcmp(ref_data_type, "iftImage") != 0) &&
         (strcmp(ref_data_type, "iftMImage") != 0) &&
         (strcmp(ref_data_type, "iftFImage") != 0))) {
        iftError("Reference data must be an image", "iftDrawVoxelSamples");
    }

    if (strcmp(ref_data_type, "iftImage") == 0) {
        iftImage *aux = (iftImage *) Z->ref_data;
        img = iftCreateImage(aux->xsize, aux->ysize, aux->zsize);
        iftSetCbCr(img, 128);
        if (iftIs3DImage(aux)) {
            max_val = iftMaximumValue(aux);
            for (p = 0; p < aux->n; p++)
                img->val[p] = (int) (255. * log(aux->val[p] + 2) / log(max_val + 2));
        } else {
            for (p = 0; p < aux->n; p++)
                img->val[p] = aux->val[p];
        }
    } else {
        if (strcmp(ref_data_type, "iftFImage") == 0) {
            iftFImage *aux = (iftFImage *) Z->ref_data;
            img = iftFImageToImage(aux, 255);
            iftSetCbCr(img, 128);
        } else {
            if (strcmp(ref_data_type, "iftMImage") == 0) {
                iftMImage *aux = (iftMImage *) Z->ref_data;
                img = iftMImageToImage(aux, 255, 0);
                iftSetCbCr(img, 128);
            }
        }
    }

    if (iftIs3DImage(img))
        B = iftSpheric(3.0);
    else
        B = iftCircular(3.0);

    switch (opt) {

        case IFT_LABEL:

            ctb = iftCreateColorTable(iftMax(Z->ngroups + 2, 1000));

            for (s = 0; s < Z->nsamples; s++) {
                p = Z->sample[s].id;
                u = iftGetVoxelCoord(img, p);
                iftDrawPoint(img, u, ctb->color[Z->sample[s].label], B, 255);
            }
            iftDestroyColorTable(&ctb);

            break;

        case IFT_CLASS:

            ctb = iftCreateColorTable(iftMax(Z->nclasses + 2, 1000));

            for (s = 0; s < Z->nsamples; s++) {
                p = Z->sample[s].id;
                u = iftGetVoxelCoord(img, p);
                iftDrawPoint(img, u, ctb->color[Z->sample[s].truelabel], B, 255);
            }
            iftDestroyColorTable(&ctb);

            break;

        case IFT_POINT:

            RGB.val[0] = 255;
            RGB.val[1] = 255;
            RGB.val[2] = 255;
            YCbCr = iftRGBtoYCbCr(RGB, 255);
            for (s = 0; s < Z->nsamples; s++) {
                p = Z->sample[s].id;
                u = iftGetVoxelCoord(img, p);
                iftDrawPoint(img, u, YCbCr, B, 255);
            }

            break;

        default:
    
            iftError("Invalid option", "iftDrawVoxelSamples");

    }

    iftDestroyAdjRel(&B);

    return (img);
}


iftImage *iftDraw2DFeatureSpace(iftDataSet *Z, iftFeatureLabel opt, iftSampleStatus status)
{

    iftImage *img=NULL;
    int       s,p;
    iftVoxel  u;
    float     xmin,xmax,ymin,ymax,minw,maxw;
    iftAdjRel *B=iftCircular(sqrtf(7.0));
    iftColor   RGB, YCbCr;
    iftColorTable *ctb=NULL;
    float *x=NULL, *y=NULL;

    if (Z->nfeats != 2)
        iftError("Feature space must be 2D", "iftDraw2DFeatureSpace");

    minw=xmin=ymin= IFT_INFINITY_FLT; maxw= xmax= ymax= IFT_INFINITY_FLT_NEG;
    for (s=0; s < Z->nsamples; s++){
        if (Z->sample[s].feat[0] < xmin)
            xmin = Z->sample[s].feat[0];
        if (Z->sample[s].feat[1] < ymin)
            ymin = Z->sample[s].feat[1];
        if (Z->sample[s].feat[0] > xmax)
            xmax = Z->sample[s].feat[0];
        if (Z->sample[s].feat[1] > ymax)
            ymax = Z->sample[s].feat[1];
        if (Z->sample[s].weight < minw)
            minw = Z->sample[s].weight;
        if (Z->sample[s].weight > maxw)
            maxw = Z->sample[s].weight;
    }

    x = iftAllocFloatArray(Z->nsamples);
    y = iftAllocFloatArray(Z->nsamples);

    for (s=0; s < Z->nsamples; s++){
        x[s] = 520*(Z->sample[s].feat[0]-xmin + 1e-10)/(xmax - xmin + 1e-10);
        y[s] = 520*(Z->sample[s].feat[1]-ymin + 1e-10)/(ymax - ymin + 1e-10);
    }

    img=iftCreateImage(540,540,1);
    iftSetCbCr(img,128);

    for (p=0; p < img->n; p++) {
        RGB.val[0] = 255;
        RGB.val[1] = 255;
        RGB.val[2] = 255;
        YCbCr = iftRGBtoYCbCr(RGB,255);
        img->val[p] = YCbCr.val[0];
        img->Cb[p]  = YCbCr.val[1];
        img->Cr[p]  = YCbCr.val[2];
    }

    int ncolors = 1000;
	
    if(opt == IFT_GROUP) {
      ncolors = Z->ngroups;
    }
    
    if ((opt == IFT_CLASS)||(opt == IFT_LABEL))
      ncolors = Z->nclasses;

    ctb = iftCategoricalColorTable(ncolors);

    float value, Red, Green, Blue;

    for (s=0; s < Z->nsamples; s++) {

      if (iftHasSampleStatus(Z->sample[s], status)) {
 
	u.x = (int)x[s] + 5;
        u.y = (int)y[s] + 5;
        u.z = 0;
        p = iftGetVoxelIndex(img,u);

        switch (opt) {

            case IFT_WEIGHT:
                value = (Z->sample[s].weight - minw) / (maxw - minw);
                iftHeatColorMapping(value,&Red,&Green,&Blue);
                RGB.val[0] = (int) 255*Red;
                RGB.val[1] = (int) 255*Green;
                RGB.val[2] = (int) 255*Blue;
                YCbCr = iftRGBtoYCbCr(RGB,255);
                iftDrawPoint(img, u, YCbCr, B, 255);
                break;

            case IFT_GROUP:
	      if (Z->sample[s].group > 0){
                iftDrawPoint(img, u, ctb->color[Z->sample[s].group - 1], B, 255);
	      }else{
		RGB.val[0] = 0;
                RGB.val[1] = 0;
                RGB.val[2] = 0;
                YCbCr = iftRGBtoYCbCr(RGB,255);
                iftDrawPoint(img, u, YCbCr, B, 255);
	      }    
	      break;

            case IFT_LABEL:
	      if (Z->sample[s].label > 0){
                iftDrawPoint(img, u, ctb->color[Z->sample[s].label - 1], B, 255);
	      }else{
		RGB.val[0] = 0;
                RGB.val[1] = 0;
                RGB.val[2] = 0;
                YCbCr = iftRGBtoYCbCr(RGB,255);
                iftDrawPoint(img, u, YCbCr, B, 255);
	      }    
	      break;

            case IFT_CLASS:
	      if (Z->sample[s].truelabel > 0){
		iftDrawPoint(img, u, ctb->color[Z->sample[s].truelabel-1], B, 255);
	      }else{
		RGB.val[0] = 0;
                RGB.val[1] = 0;
                RGB.val[2] = 0;
                YCbCr = iftRGBtoYCbCr(RGB,255);
                iftDrawPoint(img, u, YCbCr, B, 255);
	      }    
	      break;

            case IFT_POINT:
                RGB.val[0] = 0;
                RGB.val[1] = 0;
                RGB.val[2] = 0;
                YCbCr = iftRGBtoYCbCr(RGB,255);
                iftDrawPoint(img, u, YCbCr, B, 255);
                break;

            default:
                RGB.val[0] = 0;
                RGB.val[1] = 0;
                RGB.val[2] = 0;
                YCbCr = iftRGBtoYCbCr(RGB,255);
                iftDrawPoint(img, u, YCbCr, B, 255);
        }
      }
    }

    iftFree(x);
    iftFree(y);
    iftDestroyAdjRel(&B);
    iftDestroyColorTable(&ctb);

    return(img);
}



iftImage *iftProjectMeanValue(iftImage *img, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize)
{
    iftImage *mip,*aux;
    iftPoint  P1,P0;
    iftVoxel  u;
    int       i,xyviewsize,offset,p;

    /* Initialize output image */

    mip   = iftCreateImage(xviewsize,yviewsize,1);

    /* Starting at slabsize slices before, reslice the scene until
       slabsize slices after the center of the cut plane and along
       the direction of v. */

    xyviewsize = xviewsize*yviewsize;
    P0   = cutplane->pos;
    for (i=-slabsize,offset=0; i <= slabsize; i++,offset+=xyviewsize) {
        P1.x = P0.x + i*cutplane->normal.x;
        P1.y = P0.y + i*cutplane->normal.y;
        P1.z = P0.z + i*cutplane->normal.z;
        u.x  = iftRound(P1.x);
        u.y  = iftRound(P1.y);
        u.z  = iftRound(P1.z);

        if (iftValidVoxel(img,u)){
            iftSetPlanePos(cutplane,P1.x,P1.y,P1.z);
            aux  = iftGetSlice(img,cutplane,xviewsize,yviewsize);
            for (p=0; p < mip->n; p++) {
                mip->val[p] += aux->val[p];
            }
            iftDestroyImage(&aux);
        }
    }

    for (p=0; p < mip->n; p++) {
        mip->val[p] /= (2*slabsize+1);
    }

    return(mip);
}

iftImage *iftProjectMinValue(iftImage *img, iftPlane *cutplane, int xviewsize, int yviewsize,
                             int slabsize) {
    iftImage *mip, *aux;
    iftPoint P1, P0;
    iftVoxel u;
    int      i, xyviewsize, offset, p;
    int      img_max_val;


    /* Initialize output image */

    mip         = iftCreateImage(xviewsize, yviewsize, 1);
    img_max_val = iftMaximumValue(img);
    iftSetImage(mip, img_max_val + 1);

    /* Starting at slabsize slices before, reslice the scene until
       slabsize slices after the center of the cut plane and along the
       direction of v. */

    xyviewsize = xviewsize * yviewsize;
    P0         = cutplane->pos;
    for (i = -slabsize, offset = 0; i <= slabsize; i++, offset += xyviewsize) {
        P1.x = P0.x + i * cutplane->normal.x;
        P1.y = P0.y + i * cutplane->normal.y;
        P1.z = P0.z + i * cutplane->normal.z;
        u.x  = iftRound(P1.x);
        u.y  = iftRound(P1.y);
        u.z  = iftRound(P1.z);

        if (iftValidVoxel(img, u)) {
            iftSetPlanePos(cutplane, P1.x, P1.y, P1.z);
            aux = iftGetSlice(img, cutplane, xviewsize, yviewsize);
            for (p = 0; p < mip->n; p++) {
                if (mip->val[p] > aux->val[p])
                    mip->val[p] = aux->val[p];
            }
            iftDestroyImage(&aux);
        }
    }

    for (p = 0; p < mip->n; p++) {
        if (mip->val[p] == (img_max_val + 1))
            mip->val[p] = 0;
    }

    return (mip);
}

iftImage *iftProjectMaxValue(iftImage *img, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize)
{
    iftImage *mip,*aux;
    iftPoint  P1,P0;
    iftVoxel  u;
    int       i,xyviewsize,offset,p;

    /* Initialize output image */

    mip   = iftCreateImage(xviewsize,yviewsize,1);

    /* Starting at slabsize slices before, reslice the scene until
       slabsize slices after the center of the cut plane and along
       the direction of v. */

    xyviewsize = xviewsize*yviewsize;
    P0   = cutplane->pos;
    for (i=-slabsize,offset=0; i < slabsize; i++,offset+=xyviewsize) {
        P1.x = P0.x + i*cutplane->normal.x;
        P1.y = P0.y + i*cutplane->normal.y;
        P1.z = P0.z + i*cutplane->normal.z;
        u.x  = iftRound(P1.x);
        u.y  = iftRound(P1.y);
        u.z  = iftRound(P1.z);

        if (iftValidVoxel(img,u)){
            iftSetPlanePos(cutplane,P1.x,P1.y,P1.z);
            aux  = iftGetSlice(img,cutplane,xviewsize,yviewsize);
            for (p=0; p < mip->n; p++) {
                if (mip->val[p] < aux->val[p])
                    mip->val[p] = aux->val[p];
            }
            iftDestroyImage(&aux);
        }
    }

    return(mip);
}

iftImage *iftProjectMeanObjectValue(iftImage *img, iftImage *obj, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize)
{
    iftImage  *mip,*aux1,*aux2;
    iftPoint  P1,P0;
    iftVoxel  u;
    int       i,xyviewsize,offset,p;

    /* Initialize output image */

    mip   = iftCreateImage(xviewsize,yviewsize,1);

    /* Starting at slabsize slices before, reslice the scene until
       slabsize slices after the center of the cut plane and along
       the direction of v. */

    xyviewsize = xviewsize*yviewsize;
    P0 = cutplane->pos;
    for (i=-slabsize,offset=0; i <= slabsize; i++,offset+=xyviewsize) {
        P1.x = P0.x + i*cutplane->normal.x;
        P1.y = P0.y + i*cutplane->normal.y;
        P1.z = P0.z + i*cutplane->normal.z;
        u.x  = iftRound(P1.x);
        u.y  = iftRound(P1.y);
        u.z  = iftRound(P1.z);

        if (iftValidVoxel(img,u)){
            iftSetPlanePos(cutplane,P1.x,P1.y,P1.z);
            aux1  = iftGetSlice(img,cutplane,xviewsize,yviewsize);
            aux2  = iftGetSlice(obj,cutplane,xviewsize,yviewsize);
            for (p=0; p < mip->n; p++) {
                if (aux2->val[p] != 0){
                    mip->val[p] += aux1->val[p];
                }
            }
            iftDestroyImage(&aux1);
            iftDestroyImage(&aux2);
        }
    }

    for (p=0; p < mip->n; p++) {
        mip->val[p]/=(2*slabsize+1);
    }

    return(mip);
}


iftImage *iftProjectMinObjectValue(iftImage *img, iftImage *obj, iftPlane *cutplane, int xviewsize,
                                   int yviewsize, int slabsize) {
    iftImage *mip, *aux1, *aux2;
    iftPoint P1, P0;
    iftVoxel u;
    int      i, xyviewsize, offset, p;
    int      img_max_val;


    /* Initialize output image */

    mip         = iftCreateImage(xviewsize, yviewsize, 1);
    img_max_val = iftMaximumValue(img);

    iftSetImage(mip, img_max_val + 1);

    /* Starting at slabsize slices before, reslice the scene until
       slabsize slices after the center of the cut plane and along
       the direction of v. */

    xyviewsize = xviewsize * yviewsize;
    P0         = cutplane->pos;
    for (i = -slabsize, offset = 0; i <= slabsize; i++, offset += xyviewsize) {
        P1.x = P0.x + i * cutplane->normal.x;
        P1.y = P0.y + i * cutplane->normal.y;
        P1.z = P0.z + i * cutplane->normal.z;
        u.x  = iftRound(P1.x);
        u.y  = iftRound(P1.y);
        u.z  = iftRound(P1.z);

        if (iftValidVoxel(img, u)) {
            iftSetPlanePos(cutplane, P1.x, P1.y, P1.z);
            aux1 = iftGetSlice(img, cutplane, xviewsize, yviewsize);
            aux2 = iftGetSlice(obj, cutplane, xviewsize, yviewsize);
            for (p = 0; p < mip->n; p++) {
                if (aux2->val[p] != 0) {
                    if (mip->val[p] > aux1->val[p])
                        mip->val[p] = aux1->val[p];
                }
            }
            iftDestroyImage(&aux1);
            iftDestroyImage(&aux2);
        }
    }

    for (p = 0; p < mip->n; p++) {
        if (mip->val[p] == (img_max_val + 1))
            mip->val[p] = 0;
    }

    return (mip);
}


iftImage *iftProjectMaxObjectValue(iftImage *img, iftImage *obj, iftPlane *cutplane, int xviewsize,
                                   int yviewsize, int slabsize) {
    iftImage *mip, *aux1, *aux2;
    iftPoint P1, P0;
    iftVoxel u;
    int      i, xyviewsize, offset, p;

    /* Initialize output image */

    mip = iftCreateImage(xviewsize, yviewsize, 1);

    /* Starting at slabsize slices before, reslice the scene until
       slabsize slices after the center of the cut plane and along
       the direction of v. */

    xyviewsize = xviewsize * yviewsize;
    P0         = cutplane->pos;
    for (i = -slabsize, offset = 0; i <= slabsize; i++, offset += xyviewsize) {
        P1.x = P0.x + i * cutplane->normal.x;
        P1.y = P0.y + i * cutplane->normal.y;
        P1.z = P0.z + i * cutplane->normal.z;
        u.x  = iftRound(P1.x);
        u.y  = iftRound(P1.y);
        u.z  = iftRound(P1.z);

        if (iftValidVoxel(img, u)) {
            iftSetPlanePos(cutplane, P1.x, P1.y, P1.z);
            aux1 = iftGetSlice(img, cutplane, xviewsize, yviewsize);
            aux2 = iftGetSlice(obj, cutplane, xviewsize, yviewsize);
            for (p = 0; p < mip->n; p++) {
                if (aux2->val[p] != 0) {
                    if (mip->val[p] < aux1->val[p])
                        mip->val[p] = aux1->val[p];
                }
            }
            iftDestroyImage(&aux1);
            iftDestroyImage(&aux2);
        }
    }

    return (mip);
}

iftImage *iftCurvilinearSplatting(iftImage *img, iftFImage *dist, iftPlane *cutplane, int viewsize, float depth)
{
    iftMatrix  *R = iftRotationMatrixToAlignVectorWithZ(cutplane->normal);
    iftImage *proj;
    iftFImage *zbuff;
    iftVoxel u,v;
    iftPoint P1,P2;
    int p,q;
    float dmin= iftMax(depth - 0.5, 0.5),dmax= depth + 1.0;


    /* Initialize output image */

    proj   = iftCreateImage(viewsize,viewsize,1);

    /* Perform voxel splatting, using the direction of the viewing
       plane's normal and depth information*/

    zbuff   = iftCreateFImage(viewsize,viewsize,1);
    iftFSetImage(zbuff, IFT_INFINITY_FLT);

    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++){

                p   = iftGetVoxelIndex(img,u);

                if ((dist->val[p] >= dmin)&&(dist->val[p] <= dmax)) {

                    /* Set voxel coordinate to the rotated system with respect
                       to the center of the object */

                    P1.x = u.x - cutplane->pos.x;
                    P1.y = u.y - cutplane->pos.y;
                    P1.z = u.z - cutplane->pos.z;
                    P2   = iftTransformPoint(R,P1);
                    P1.x = P2.x + viewsize/2.0;
                    P1.y = P2.y + viewsize/2.0;
                    P1.z = P2.z + viewsize/2.0;

                    /* Perform voxel splatting */

                    v.z = 0;

                    v.y = (int) (P1.y);
                    v.x = (int) (P1.x);

                    if (iftValidVoxel(proj,v)){
                        q   = iftGetVoxelIndex(proj,v);
                        if (zbuff->val[q]>P1.z){
                            proj->val[q]  = img->val[p];
                            zbuff->val[q] = P1.z;
                        }
                    }

                    v.y = (int) (P1.y);
                    v.x = (int) (P1.x+1.0);

                    if (iftValidVoxel(proj,v)){
                        q   = iftGetVoxelIndex(proj,v);
                        if (zbuff->val[q]>P1.z){
                            proj->val[q]  = img->val[p];
                            zbuff->val[q] = P1.z;
                        }
                    }

                    v.y = (int) (P1.y+1.0);
                    v.x = (int) (P1.x);

                    if (iftValidVoxel(proj,v)){
                        q   = iftGetVoxelIndex(proj,v);
                        if (zbuff->val[q]>P1.z){
                            proj->val[q]  = img->val[p];
                            zbuff->val[q] = P1.z;
                        }
                    }

                    v.y = (int) (P1.y+1.0);
                    v.x = (int) (P1.x+1.0);

                    if (iftValidVoxel(proj,v)){
                        q   = iftGetVoxelIndex(proj,v);
                        if (zbuff->val[q]>P1.z){
                            proj->val[q]  = img->val[p];
                            zbuff->val[q] = P1.z;
                        }
                    }
                }
            }

    iftDestroyMatrix(&R);
    iftDestroyFImage(&zbuff);

    return(proj);
}


iftImage *iftColorizeImageByLabels(iftImage *orig, iftImage *label, iftColorTable *ctb)
{
  iftImage *cimg = NULL;

  iftVerifyImageDomains(orig,label,"iftColorizeImageByLabels");
 
  cimg = iftCreateImage(label->xsize,label->ysize,label->zsize);
  for (int p=0; p < cimg->n; p++)
    cimg->val[p]=orig->val[p];

  iftSetCbCr(cimg,iftNormalizationValue(iftMaximumValue(cimg))/2);

  for (int p=0; p < cimg->n; p++) {
    if (label->val[p] != 0){
      cimg->Cb[p] = ctb->color[label->val[p]%(ctb->ncolors-1)].val[1];
      cimg->Cr[p] = ctb->color[label->val[p]%(ctb->ncolors-1)].val[2];
    }
  }

  return(cimg);
}

iftImage *iftColorizeComp(iftImage *label)
{

    iftColorTable *ctb;
    iftImage *img;
    int p;

    ctb = iftCreateRandomColorTable(iftMaximumValue(label)+1);
    img = iftCreateImage(label->xsize,label->ysize,label->zsize);
    iftSetImage(img,255);
    iftSetCbCr(img,128);

    for (p=0; p < label->n; p++)
        if (label->val[p]>0){
            img->val[p] = ctb->color[label->val[p]].val[0];
            img->Cb[p]  = ctb->color[label->val[p]].val[1];
            img->Cr[p]  = ctb->color[label->val[p]].val[2];
        }
    iftDestroyColorTable(&ctb);
    return(img);
}

iftImage *iftColorizeCompOverImage(iftImage *orig, iftImage *label)
{
    iftColorTable *ctb;
    iftImage *img;
    int p;

    iftVerifyImageDomains(orig,label,"iftColorizeCompOverImage");

    ctb = iftCreateRandomColorTable(iftMaximumValue(label)+1);
    img = iftCreateImage(label->xsize,label->ysize,label->zsize);
    for (p=0; p < img->n; p++)
        img->val[p]=orig->val[p];

    iftSetCbCr(img,iftNormalizationValue(iftMaximumValue(img))/2);

    for (p=0; p < label->n; p++)
        if (label->val[p]>0){
	  img->Cb[p] = ctb->color[label->val[p]].val[1];
	  img->Cr[p] = ctb->color[label->val[p]].val[2];
        }
    
    iftDestroyColorTable(&ctb);
    return(img);
}

void iftTiltSpinFromPrincipalAxis(iftImage *bin, float *tilt, float *spin)
{
    iftVector  paxis;

    paxis   = iftPrincipalAxis(bin);
    if (fabs(paxis.z) > IFT_EPSILON){
        *tilt   = atan(paxis.y/paxis.z);
        paxis.z = paxis.z/cos(*tilt);
        *tilt   = (*tilt) * 180 / IFT_PI;
        *spin   = atan(paxis.x/paxis.z) * 180 / IFT_PI;
    }else{ // paxis.z = 0
        if (fabs(paxis.x) < IFT_EPSILON) // paxis.x = 0
            *spin = 0;
        else
            *spin = 90;
        if (fabs(paxis.y) < IFT_EPSILON) // paxis.y = 0
            *tilt = 0;
        else
            *tilt = 90;
    }
}



void iftDrawLineRecursiveAux(iftImage *img, iftVoxel px1, iftVoxel pxn, iftColor YCbCr, iftAdjRel *B,int maxRangeValue) {
    iftVoxel px;

    px.x = (px1.x + pxn.x)/2;
    px.y = (px1.y + pxn.y)/2;
    px.z = (px1.z + pxn.z)/2;

    if (iftValidVoxel(img, px)) {
        iftDrawPoint(img, px, YCbCr, B, maxRangeValue);
    }

    if(iftVoxelDistance(px1, pxn) <= 1.5)
        return;

    iftDrawLineRecursiveAux(img, px1, px, YCbCr, B, maxRangeValue);
    iftDrawLineRecursiveAux(img, px, pxn, YCbCr, B, maxRangeValue);
}

void iftDrawLineRecursive(iftImage *img, iftVoxel px1, iftVoxel pxn, iftColor YCbCr, iftAdjRel *B) {
    int maxRangeValue = iftNormalizationValue(iftMaximumValue(img));
    if (iftValidVoxel(img, px1)) {
        iftDrawPoint(img, px1, YCbCr, B, maxRangeValue);
    }

    if (iftValidVoxel(img, pxn)) {
        iftDrawPoint(img, pxn, YCbCr, B, maxRangeValue);
    }

    if((iftValidVoxel(img, px1)) && (iftValidVoxel(img, pxn)))
        iftDrawLineRecursiveAux(img, px1, pxn, YCbCr, B,maxRangeValue);
}

iftImage *iftMaskImageFromPolygon2D(iftVoxel *points, int n, int xsize, int ysize)
{
    int i;
    iftImage *closed = NULL;
    iftImage *mask = iftCreateImage(xsize, ysize, 1);
    iftAdjRel *B = iftCircular(0.0);
    iftColor color;

    color.val[0] = 1;
    color.val[1] = 0;
    color.val[2] = 0;

    for(i = 0; i < n; i++)
    {
        iftVoxel p, q;

        p = points[i];
        q = points[(i+1) % n];

        iftDrawLineRecursive(mask, p, q, color, B);
    }

    closed = iftCloseBasins(mask,NULL,NULL);
    iftDestroyImage(&mask);

    mask = closed;

    iftDestroyAdjRel(&B);

    return mask;
}

void iftDrawBoundingBoxBordersInPlace(iftImage *img, iftBoundingBox bb, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint)
{
    if (img == NULL)
        iftError("Image is NULL", "iftDrawBoundingBoxBordersInPlace");
    if (!iftValidVoxel(img, bb.begin))
        iftError("Invalid begin voxel (%d, %d, %d) in Image Domain: %d, %d, %d\n",
                 "iftDrawBoundingBoxBordersInPlace", bb.begin.x, bb.begin.y, bb.begin.z,
                 img->xsize, img->ysize, img->zsize);
    if (!iftValidVoxel(img, bb.end))
        iftError("Invalid end voxel (%d, %d, %d) in Image Domain: %d, %d, %d\n",
                 "iftDrawBoundingBoxBordersInPlace", bb.end.x, bb.end.y, bb.end.z,
                 img->xsize, img->ysize, img->zsize);

    if (iftIs3DImage(img))
        iftDrawBoundingBoxBorders3DInPlace(img, bb, YCbCr, A, drawCentralPoint);
    else
        iftDrawBoundingBoxBorders2DInPlace(img, bb, YCbCr, A, drawCentralPoint);
}

void iftDrawBoundingBoxArrayBordersInPlace(iftImage *img, iftBoundingBoxArray *bbs, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint)
{
    if (img == NULL)
        iftError("Image is NULL", "iftDrawBoundingBoxArrayBordersInPlace");

    if (bbs->n == 0)
        iftError("Bounding box array is empty", "iftDrawBoundingBoxArrayBordersInPlace");

    for(int i = 0; i < bbs->n; i++)
        iftDrawBoundingBoxBordersInPlace(img, bbs->val[i], YCbCr, A, drawCentralPoint);
}

void iftDrawRoiBordersInPlace(iftImage *img, iftRoi *roi, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint)
{
    if (img == NULL)
        iftError("Image is NULL", "iftDrawRoiBordersInPlace");

    if (roi->n == 0)
        iftError("ROI is empty", "iftDrawRoiBordersInPlace");

    for(int i = 0; i < roi->n; i++)
        if (!iftValidVoxel(img, roi->val[i]))
            iftError("Invalid ROI voxel (%d, %d, %d) in Image Domain: %d, %d, %d\n",
                 "iftDrawRoiBordersInPlace", roi->val[i].x, roi->val[i].y, roi->val[i].z,
                 img->xsize, img->ysize, img->zsize);

    int Imax = iftNormalizationValue(iftMaximumValue(img));
    iftAdjRel *adj = iftIs3DImage(img) ? iftSpheric(1.0) : iftCircular(1.0);
    iftSize maskSize = {.x = img->xsize, .y = img->ysize, .z = img->zsize};
    iftImage *mask = iftMaskFromRoi(roi, maskSize);

    /* draw all the border pixels (i.e. the ones who have at least one neighbor which is not part of the ROI) */
    #pragma omp parallel for
    for(int i = 0; i < roi->n; i++) {
        iftVoxel v = roi->val[i];
        bool drawPx = false;
        for(int j = 0; j < adj->n; j++) {
            iftVoxel u = iftGetAdjacentVoxel(adj, v, j);
            if(iftValidVoxel(img, u) && iftImgVoxelVal(mask, v) != iftImgVoxelVal(mask, u))
                drawPx = true;
        }
        if(drawPx)
            iftDrawPoint(img, v, YCbCr, A, Imax);
    }

    /* draw the central point */
    if(drawCentralPoint)
        iftDrawPoint(img, iftRoiCenterVoxel(img, roi), YCbCr, A, Imax);
    
}

void iftDrawRoiArrayBordersInPlace(iftImage *img, iftRoiArray *roiArray, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint)
{
    if (img == NULL)
        iftError("Image is NULL", "iftDrawRoiArrayBordersInPlace");

    if (roiArray->n == 0)
        iftError("ROI array is empty", "iftDrawRoiArrayBordersInPlace");

    for(int i = 0; i < roiArray->n; i++)
        iftDrawRoiBordersInPlace(img, roiArray->val[i], YCbCr, A, drawCentralPoint);
}

void iftDrawBinaryLabeledSeeds(iftImage *img,iftLabeledSet *seeds,iftColor YCbCr, iftAdjRel *A)
{
    iftColor YCbCr_compl;
    int Imax = iftNormalizationValue(iftMaximumValue(img));

    YCbCr_compl.val[0] = YCbCr.val[0];
    YCbCr_compl.val[1] = Imax-YCbCr.val[1];
    YCbCr_compl.val[2] = Imax-YCbCr.val[2];

    iftLabeledSet *S = seeds;
    while (S != NULL) {
        int p = S->elem;
        iftVoxel u = iftGetVoxelCoord(img,p);
        if (S->label == 1){
            iftDrawPoint(img,u,YCbCr,A,Imax);
        } else {
            iftDrawPoint(img,u,YCbCr_compl,A,Imax);
        }
        S = S->next;
    }
}

void iftDrawColorLegend(iftImage *img, iftColorTable *ctb, int firstColor, int lastColor, char *location, int colorBoxSize)
{
    int nColors = lastColor - firstColor + 1;
    int Imax = iftNormalizationValue(iftMaximumValue(img));

    /* dimensions of the legend */
    int legendHeight = (colorBoxSize-1)*nColors;
    int legendWidth = colorBoxSize;

    /*determine the top-left and bottom-right corners */
    iftVoxel legendTopLeftCorner = {.x = 0, .y = 0, .z = 0}, legendBottomRightCorner = {.x = 0, .y = 0, .z = 0};

    if(iftCompareStrings(location, "top-left"))
    {
        legendTopLeftCorner.x = 0;
        legendTopLeftCorner.y = 0;
        legendBottomRightCorner.x = legendWidth;
        legendBottomRightCorner.y = legendHeight;
    }

    if(iftCompareStrings(location, "top-right"))
    {
        legendTopLeftCorner.x = img->xsize - legendWidth;
        legendTopLeftCorner.y = 0;
        legendBottomRightCorner.x = img->xsize;
        legendBottomRightCorner.y = legendHeight;
    }

    if(iftCompareStrings(location, "bottom-left"))
    {
        legendTopLeftCorner.x = 0;
        legendTopLeftCorner.y = img->ysize - legendHeight;
        legendBottomRightCorner.x = legendWidth;
        legendBottomRightCorner.y = img->ysize;
    }

    if(iftCompareStrings(location, "bottom-right"))
    {
        legendTopLeftCorner.x = img->xsize - legendWidth;
        legendTopLeftCorner.y = img->ysize - legendHeight;
        legendBottomRightCorner.x = img->xsize;
        legendBottomRightCorner.y = img->ysize;
    }

    /* colors for background and edges */
    iftColor whiteColor;
    whiteColor.val[0] = Imax;
    whiteColor.val[1] = Imax;
    whiteColor.val[2] = Imax;
    whiteColor = iftRGBtoYCbCr(whiteColor, Imax);

    iftColor blackColor;
    blackColor.val[0] = 0;
    blackColor.val[1] = 0;
    blackColor.val[2] = 0;
    blackColor = iftRGBtoYCbCr(blackColor, Imax);

    iftVoxel pos = {.x = 0, .y = 0, .z = 0};

    /* draw the box for the legend */
    for(pos.y = legendTopLeftCorner.y; pos.y < legendBottomRightCorner.y; pos.y++)
    {
        for(pos.x = legendTopLeftCorner.x; pos.x < legendBottomRightCorner.x; pos.x++)
        {
            int p = iftGetVoxelIndex(img, pos);
            if(pos.y == legendTopLeftCorner.y || pos.y == legendBottomRightCorner.y-1 || pos.x == legendTopLeftCorner.x || pos.x == legendBottomRightCorner.x-1)
                iftSetYCbCr(img, p, blackColor);
            else
                iftSetYCbCr(img, p, whiteColor);
        }
    }

    /* draw the boxes for every color */
    int i, j;
    iftVoxel colorLeftCorner  = {.x = legendTopLeftCorner.x, .y = legendTopLeftCorner.y, .z = 0};
    for(int c = firstColor; c <= lastColor; c++)
    {
        iftColor color = ctb->color[c];
        for(pos.y = colorLeftCorner.y, i = 0; i < colorBoxSize; pos.y++, i++)
        {
            for(pos.x = colorLeftCorner.x, j = 0; j < colorBoxSize; pos.x++, j++)
            {
                int p = iftGetVoxelIndex(img, pos);
                if(i == 0 || i == colorBoxSize-1 || j == 0 || j == colorBoxSize-1)
                    iftSetYCbCr(img, p, blackColor);
                else
                    iftSetYCbCr(img, p, color);
            }
        }
        colorLeftCorner.y += colorBoxSize-1;
    }
}

void iftDrawFancyColorLegend(iftImage *img, iftColorTable *ctb, int firstColor, int lastColor, char *location, int colorBoxSize)
{
    int nColors = lastColor - firstColor + 1;
    int Imax = iftNormalizationValue(iftMaximumValue(img));

    /* dimensions of the color boxes */
    int spaceSize = ceilf((float)colorBoxSize/4.0);

    /* dimensions of the legend */
    int legendHeight = (colorBoxSize*nColors) + (spaceSize*(nColors+1));
    int legendWidth = colorBoxSize + spaceSize*2;

    /*determine the top-left and bottom-right corners */
    iftVoxel legendTopLeftCorner = {.x = 0, .y = 0, .z = 0}, legendBottomRightCorner = {.x = 0, .y = 0, .z = 0};

    if(iftCompareStrings(location, "top-left"))
    {
        legendTopLeftCorner.x = spaceSize;
        legendTopLeftCorner.y = spaceSize;
        legendBottomRightCorner.x = spaceSize + legendWidth;
        legendBottomRightCorner.y = spaceSize + legendHeight;
    }

    if(iftCompareStrings(location, "top-right"))
    {
        legendTopLeftCorner.x = img->xsize-1 - legendWidth - spaceSize;
        legendTopLeftCorner.y = spaceSize;
        legendBottomRightCorner.x = img->xsize-1 - spaceSize;
        legendBottomRightCorner.y = spaceSize + legendHeight;
    }

    if(iftCompareStrings(location, "bottom-left"))
    {
        legendTopLeftCorner.x = spaceSize;
        legendTopLeftCorner.y = img->ysize-1 - legendHeight - spaceSize;
        legendBottomRightCorner.x = spaceSize + legendWidth;
        legendBottomRightCorner.y = img->ysize-1 - spaceSize;
    }

    if(iftCompareStrings(location, "bottom-right"))
    {
        legendTopLeftCorner.x = img->xsize-1 - legendWidth - spaceSize;
        legendTopLeftCorner.y = img->ysize-1 - legendHeight - spaceSize;
        legendBottomRightCorner.x = img->xsize-1 - spaceSize;
        legendBottomRightCorner.y = img->ysize-1 - spaceSize;
    }

    /* colors for background and edges */
    iftColor whiteColor;
    whiteColor.val[0] = Imax;
    whiteColor.val[1] = Imax;
    whiteColor.val[2] = Imax;
    whiteColor = iftRGBtoYCbCr(whiteColor, Imax);

    iftColor blackColor;
    blackColor.val[0] = 0;
    blackColor.val[1] = 0;
    blackColor.val[2] = 0;
    blackColor = iftRGBtoYCbCr(blackColor, Imax);

    iftVoxel pos = {.x = 0, .y = 0, .z = 0};

    /* draw the box for the legend */
    for(pos.y = legendTopLeftCorner.y; pos.y < legendBottomRightCorner.y; pos.y++)
    {
        for(pos.x = legendTopLeftCorner.x; pos.x < legendBottomRightCorner.x; pos.x++)
        {
            int p = iftGetVoxelIndex(img, pos);
            if(pos.y == legendTopLeftCorner.y || pos.y == legendBottomRightCorner.y-1 || pos.x == legendTopLeftCorner.x || pos.x == legendBottomRightCorner.x-1)
                iftSetYCbCr(img, p, blackColor);
            else
                iftSetYCbCr(img, p, whiteColor);
        }
    }

    /* draw the boxes for every color */
    int i, j;
    iftVoxel colorLeftCorner  = {.x = legendTopLeftCorner.x + spaceSize, .y = legendTopLeftCorner.y + spaceSize, .z = 0};
    for(int c = firstColor; c <= lastColor; c++)
    {
        iftColor color = ctb->color[c];
        for(pos.y = colorLeftCorner.y, i = 0; i < colorBoxSize; pos.y++, i++)
        {
            for(pos.x = colorLeftCorner.x, j = 0; j < colorBoxSize; pos.x++, j++)
            {
                int p = iftGetVoxelIndex(img, pos);
                if(i == 0 || i == colorBoxSize-1 || j == 0 || j == colorBoxSize-1)
                    iftSetYCbCr(img, p, blackColor);
                else
                    iftSetYCbCr(img, p, color);
            }
        }
        colorLeftCorner.y += colorBoxSize + spaceSize;
    }
}
