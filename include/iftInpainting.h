#ifndef IFT_INPAINTING_H_
#define IFT_INPAINTING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftAdjacency.h"
#include "iftImage.h"
#include "iftFImage.h"
#include "iftImageMath.h"
#include "iftMathMorph.h"


typedef struct ift_texture {
    int   *feat;
    char   valid;
} iftTexture;

typedef struct ift_texture_samples {
    iftTexture *sample;
    int nfeats, nsamples;
} iftTextureSamples;

/* Method for inpainting implemented as proposed in

A. Criminisi, P. Perez and K. Toyama, Region Filling and Object
Removal by Exemplar-Based Image Inpainting, \emph{IEEE TRANSACTIONS ON
IMAGE PROCESSING}, VOL. 13, NO. 9, SEP 2004.

*/
iftImage *iftInpainting(iftImage *img, iftImage *mask, iftTextureSamples *texture,
                        iftAdjRel *A);

/* This function samples texture from a shell around the mask's boundary, ranging
	from 2*adj_radius to 4*adj_radius outside it, and uses iftInpainting for
	foreground completion.
 */
iftImage *iftInpaintingFromSamplesAroundMaskBoundaryOutside(iftImage *img, iftImage *mask, float adj_radius);
/* This function samples texture inside the mask's boundary and inpaints it
	to decontamine	background colors from the object's border.
 */
iftImage *iftInpaintBoundary(iftImage *img, iftImage *mask);

/* Texture sampling */

iftTextureSamples *iftCreateTextureSamples(int nsamples, int nfeats);
void               iftDestroyTextureSamples(iftTextureSamples **texture);

/* This function samples texture on the foreground portion of the mask. A sample is
	only valid if it is *entirely* contained within the mask, hence, those near the
	mask's boundary are disconsidered because they may overlap with the filling area.
	This function can be used to determine that the entire background should be used
	for sampling, for example.
 */
iftTextureSamples *iftTextureSamplesOnMask(iftImage *img, iftImage *mask, iftAdjRel *A);
iftTextureSamples *iftTextureSamplesAroundMaskBoundaryOutside(iftImage *img, iftImage *mask, float adj_radius);
iftTextureSamples *iftTextureSamplesAroundMaskBoundaryInside(iftImage *img, iftImage *mask, float adj_radius);


#ifdef __cplusplus
}
#endif

#endif
