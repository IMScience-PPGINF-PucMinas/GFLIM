#include "ift/core/dtypes/Color.h"

#include "ift/core/dtypes/Dict.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/IntQueue.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/tools/Regex.h"


/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Definition for a function point that converts a color from a color space to another one,
 * @author Samuel Martins
 * @date Jun 12, 2018
 */
typedef iftColor (*iftConvertColorFunc)(iftColor cin, int normalization_value);

/**
 * @brief Create a color table with <n_colors> YCbCr colors of 3*8 bits resulting from the gradient between
 * the beginning RGB color <begin_RGB> and ending RGB color <end_RGB>.
 * @author Samuel Martins
 * @date Jun 12, 2018
 */
iftColorTable *_iftCreateRGBGradientColorTable(iftColor begin_RGB, iftColor end_RGB, int n_colors) {
    if (n_colors < 2)
        iftWarning("Number of colors: %d < 2", "iftCreateRGBGradientColorTable", n_colors);

    iftColorTable *ctb = iftCreateColorTable(n_colors);

    for (int c = 0; c < n_colors; c++) {
        float alpha = ((float) c) / (n_colors-1);

        iftColor RGB;

        RGB.val[0] = iftRound(((1-alpha) * begin_RGB.val[0]) + (alpha * end_RGB.val[0]));
        RGB.val[1] = iftRound(((1-alpha) * begin_RGB.val[1]) + (alpha * end_RGB.val[1]));
        RGB.val[2] = iftRound(((1-alpha) * begin_RGB.val[2]) + (alpha * end_RGB.val[2]));

        ctb->color[c] = iftRGBtoYCbCr(RGB, 255); 
    }

    return ctb;
}

/**
 * @brief Applies a breadth-first search on the Intermediate Color Ticks in order to put the color
 * ticks in a good order to generate different colors.
 * @param inter_ticks Intermediate Color Ticks.
 * @param n           Number of Intermediate Color Ticks.
 * @param ticks       Ticks according a BFS.
 */
void iftBFSColorTicks(int *inter_ticks, int n, iftIntArray *ticks) {
    if (n <= 0)
        return;

    int i = 1; // index from where the ticks will be added into ticks array

    iftIntQueue *bfs = iftCreateIntQueue((n+1)*4); // stores the beginning and ending position of the sub-array

    iftInsertIntQueue(bfs, 0);
    iftInsertIntQueue(bfs, n-1);

    while (!iftIsIntQueueEmpty(bfs)) {
        int i0, i1;
        iftRemoveIntQueue(bfs, &i0);
        iftRemoveIntQueue(bfs, &i1);

        if (i0 == i1) {
            ticks->val[i++] = inter_ticks[i0];
        }
        else if (i0 < i1) {
            int c = (i0 + i1) / 2;
            ticks->val[i++] = inter_ticks[c];

            iftInsertIntQueue(bfs, i0);
            iftInsertIntQueue(bfs, c-1);
            iftInsertIntQueue(bfs, c+1);
            iftInsertIntQueue(bfs, i1);
        }
    }

    iftDestroyIntQueue(&bfs);
}


iftColor *iftGetAllColors(int n_max_colors, iftIntArray *ticks) {
    // base to generate colors
    iftColor mask[7] = {{.val={1,1,0}}, {.val={0,0,1}}, {.val={1,0,0}},
                        {.val={0,1,0}}, {.val={1,0,1}}, {.val={0,1,1}},
                        {.val={1,1,1}}};

    iftColor *all_colors = (iftColor*) iftAlloc(n_max_colors, sizeof(iftColor));
    iftDict *visited     = iftCreateDict(); // checks if the colors were already visited

    int c = 0;
    for (int r = 0; r < ticks->n; r++) {
        int R = ticks->val[r];

        for (int g = 0; g < ticks->n; g++) {
            int G = ticks->val[g];

            for (int b = 0; b < ticks->n; b++) {
                int B = ticks->val[b];

                for (int i = 0; i < 7; i++) {
                    iftColor color;
                    color.val[0] = R * mask[i].val[0];
                    color.val[1] = G * mask[i].val[1];
                    color.val[2] = B * mask[i].val[2];

                    char color_key[32];
                    sprintf(color_key, "%d%d%d", color.val[0], color.val[1], color.val[2]);

                    // if the color was not still visited
                    if (!iftDictContainKey(color_key, visited, NULL)) {
                        iftInsertIntoDict(color_key, NULL, visited); // visits the key

                        all_colors[c].val[0] = R * mask[i].val[0];
                        all_colors[c].val[1] = G * mask[i].val[1];
                        all_colors[c].val[2] = B * mask[i].val[2];
                        c++;
                    }
                }
            }
        }
    }

    iftDestroyDict(&visited);


    return all_colors;
}
/******************************************************************/

/********************** PUBLIC FUNCTIONS *************************/
iftColorSpace iftColorSpaceStrToColorSpace(char *colSpaceStr) {
    iftColorSpace colSpace = 0;

    if(!strcmp(colSpaceStr, "ycbcr"))
        colSpace = YCbCr_CSPACE;
    else if(!strcmp(colSpaceStr, "ycbcr_norm"))
        colSpace = YCbCrNorm_CSPACE;
    else if(!strcmp(colSpaceStr, "rgb"))
        colSpace = RGB_CSPACE;
    else if(!strcmp(colSpaceStr, "rgb_norm"))
        colSpace = RGBNorm_CSPACE;
    else if(!strcmp(colSpaceStr, "gray"))
        colSpace = GRAY_CSPACE;
    else if(!strcmp(colSpaceStr, "gray_norm"))
        colSpace = GRAYNorm_CSPACE;
    else if(!strcmp(colSpaceStr, "weighted_ycbcr"))
        colSpace = WEIGHTED_YCbCr_CSPACE;
    else if(!strcmp(colSpaceStr, "lab"))
        colSpace = LAB_CSPACE;
    else if(!strcmp(colSpaceStr, "lab_norm"))
        colSpace = LABNorm_CSPACE;
    else if(!strcmp(colSpaceStr, "lab_norm_2"))
        colSpace = LABNorm2_CSPACE;
    else if(!strcmp(colSpaceStr, "hsv"))
        colSpace = HSV_CSPACE;
    
    return colSpace;
}

iftColorTable *iftCreateColorTable(int n_colors) {
  
  if (n_colors <= 0)
      iftError("Invalid num of colors: %d <= 0", "iftCreateColorTable", n_colors);
  
    iftColorTable *ctb = (iftColorTable*) iftAlloc(1, sizeof(iftColorTable));
    ctb->ncolors       = n_colors;
    ctb->color         = (iftColor*) iftAlloc(n_colors, sizeof(iftColor));
    if (ctb->color == NULL)
        iftError("Cannot allocate Color Table", "iftCreateColorTable");

    int S = 255;

    iftRandomSeed(time(NULL));

    for (int c = 0, h = iftRandomInteger(0,359); c < n_colors; c++, h = iftRandomInteger(0,359)) {
      ctb->color[c].val[0] = h;
      ctb->color[c].val[1] = S;
      ctb->color[c].val[2] = 255;
      if (c%10 == 0) {
    	S = S - 10;
    	if (S < 20)
    	  S = 255;
      }
      ctb->color[c] = iftRGBtoYCbCr(iftHSVtoRGB(ctb->color[c], 255), 255);
    }

    return ctb;
}

iftColorTable *iftCreateRandomColorTable(int n_colors) {
  
  if (n_colors <= 0)
      iftError("Invalid num of colors: %d <= 0", "iftCreateRandomColorTable", n_colors);
  
    iftColorTable *ctb = (iftColorTable*) iftAlloc(1, sizeof(iftColorTable));
    ctb->ncolors       = n_colors;
    ctb->color         = (iftColor*) iftAlloc(n_colors, sizeof(iftColor));
    if (ctb->color == NULL)
        iftError("Cannot allocate Color Table", "iftCreateColorTable");

    iftRandomSeed(time(NULL));

    for (int c = 0; c < n_colors; c++) {
      ctb->color[c].val[0] = iftRandomInteger(0,255);
      ctb->color[c].val[1] = iftRandomInteger(0,255);
      ctb->color[c].val[2] = iftRandomInteger(0,255);
      ctb->color[c]        = iftRGBtoYCbCr(ctb->color[c], 255);
    }

    return ctb;
}

iftFColorMatrix *iftCreateFColorMatrix(int n_colors) {

    if (n_colors <= 0)
        iftError("Invalid num of colors: %d <= 0", "iftCreateColorTable", n_colors);

    iftFColorMatrix *cm = (iftFColorMatrix*) iftAlloc(1, sizeof(iftFColorMatrix));
    cm->ncolors       = n_colors;
    cm->color         = (iftFColor*) iftAlloc(n_colors, sizeof(iftFColor));
    if (cm->color == NULL)
        iftError("Cannot allocate Color Table", "iftCreateColorTable");

    return cm;
}

iftColorTable *iftCreateHeatMapColorTable(float delta, int scalingFactor){
    if(delta <= 0 || delta > 1){
        iftWarning("delta must be between (0.0 - 1.0]","iftCreateHeatMapColorTable");
    }
    float begin = 0.0;
    float end = 1.0;
    float currentValue = begin;
    float R,G,B;
    int ncolors = ((end-begin)/(delta)) + 1;
    iftColorTable *colorTable = (iftColorTable*) iftAlloc(1, sizeof(iftColorTable));
    colorTable->ncolors = ncolors;
    colorTable->color = (iftColor*) iftAlloc(ncolors, sizeof(iftColor));
    int currentColorTableRow = 0;
    while(currentValue <= end){
        iftHeatColorMapping(currentValue, &R, &G, &B);
        colorTable->color[currentColorTableRow].val[0] = R*scalingFactor;
        colorTable->color[currentColorTableRow].val[1] = G*scalingFactor;
        colorTable->color[currentColorTableRow].val[2] = B*scalingFactor;
        currentValue += delta;
        currentColorTableRow++;
    }
    return colorTable;
}

void iftConvertRGBColorTableToYCbCrColorTable(iftColorTable *colorTable_rgb,int normalization_value){
    iftColor yCbCr;
    iftColor rgb;
    for (int i = 0; i < colorTable_rgb->ncolors; ++i) {
        rgb.val[0] = colorTable_rgb->color[i].val[0];
        rgb.val[1] = colorTable_rgb->color[i].val[1];
        rgb.val[2] = colorTable_rgb->color[i].val[2];
        yCbCr = iftRGBtoYCbCr(rgb, normalization_value);
        colorTable_rgb->color[i].val[0] = yCbCr.val[0];
        colorTable_rgb->color[i].val[1] = yCbCr.val[1];
        colorTable_rgb->color[i].val[2] = yCbCr.val[2];
    }
}

iftColorTable *iftConvertYCbCrColorTableToRGBColorTable(const iftColorTable *ctb, int normalization_value) {
    iftColorTable *rgb_ctb = iftCreateColorTable(ctb->ncolors);

    #pragma omp parallel for
    for (int c = 0; c < ctb->ncolors; c++)
        rgb_ctb->color[c] = iftYCbCrtoRGB(ctb->color[c], normalization_value);


    return rgb_ctb;
}

void iftDestroyColorTable(iftColorTable **ctb)
{
    iftColorTable *aux=*ctb;

    if (aux != NULL) {
        iftFree(aux->color);
        iftFree(aux);
        *ctb = NULL;
    }

}

void iftDestroyFColorMatrix(iftFColorMatrix **cm)
{
    iftFColorMatrix *aux=*cm;

    if (aux != NULL) {
        iftFree(aux->color);
        iftFree(aux);
        *cm = NULL;
    }

}

iftColor iftRGBColor(int R, int G, int B) {
    iftColor color;
    color.val[0] = R;
    color.val[1] = G;
    color.val[2] = B;

    return color;
}

iftFColor iftConvertPixelColorSpace(iftFColor colorIn, char colSpaceIn, char colSpaceOut, int normValue)
{
    iftColor auxColor = {.val[0] = 0, .val[1] = 0, .val[2] = 0};
    iftFColor auxColorF = {.val[0] = 0, .val[1] = 0, .val[2] = 0}, colorOut = {.val[0] = 0, .val[1] = 0, .val[2] = 0};

    if(colSpaceIn == colSpaceOut) {
        colorOut.val[0] = colorIn.val[0];
        colorOut.val[1] = colorIn.val[1];
        colorOut.val[2] = colorIn.val[2];
    } else {
        /* convert from any color space to YCbCr */
        switch(colSpaceIn) {
            case YCbCr_CSPACE:
                auxColor.val[0] = colorIn.val[0];
                auxColor.val[1] = colorIn.val[1];
                auxColor.val[2] = colorIn.val[2];
                break;

            case YCbCrNorm_CSPACE:
                auxColor.val[0] = colorIn.val[0]*(float)normValue;
                auxColor.val[1] = colorIn.val[1]*(float)normValue;
                auxColor.val[2] = colorIn.val[2]*(float)normValue;
                break;

            case RGB_CSPACE:
                auxColor.val[0] = colorIn.val[0];
                auxColor.val[1] = colorIn.val[1];
                auxColor.val[2] = colorIn.val[2];
                auxColor = iftRGBtoYCbCr(auxColor, normValue);
                break;

            case RGBNorm_CSPACE:
                auxColor.val[0] = colorIn.val[0]*(float)normValue;
                auxColor.val[1] = colorIn.val[1]*(float)normValue;
                auxColor.val[2] = colorIn.val[2]*(float)normValue;
                auxColor = iftRGBtoYCbCr(auxColor, normValue);
                break;

            case GRAY_CSPACE:
                auxColor.val[0] = colorIn.val[0];
                break;

            case GRAYNorm_CSPACE:
                auxColor.val[0] = colorIn.val[0]*(float)normValue;
                break;

            case HSV_CSPACE:
                auxColor.val[0] = colorIn.val[0];
                auxColor.val[1] = colorIn.val[1];
                auxColor.val[2] = colorIn.val[2];
                auxColor = iftHSVtoRGB(auxColor, normValue);
                auxColor = iftRGBtoYCbCr(auxColor, normValue);
                break;

            case LAB_CSPACE:
                auxColorF.val[0] = colorIn.val[0];
                auxColorF.val[1] = colorIn.val[1];
                auxColorF.val[2] = colorIn.val[2];
                auxColor = iftLabtoRGB(auxColorF, normValue);
                auxColor = iftRGBtoYCbCr(auxColor, normValue);
                break;

            case LABNorm_CSPACE:
                auxColorF.val[0] = colorIn.val[0]*(float)normValue;
                auxColorF.val[1] = colorIn.val[1]*(float)normValue;
                auxColorF.val[2] = colorIn.val[2]*(float)normValue;
                auxColor = iftLabtoRGB(auxColorF, normValue);
                auxColor = iftRGBtoYCbCr(auxColor, normValue);
                break;

            case LABNorm2_CSPACE:
                auxColorF.val[0] = colorIn.val[0]*(float)normValue;
                auxColorF.val[1] = colorIn.val[1]*(float)normValue;
                auxColorF.val[2] = colorIn.val[2]*(float)normValue;
                auxColor = iftLabtoRGB(auxColorF, normValue);
                auxColor = iftRGBtoYCbCr(auxColor, normValue);
                break;

            default:
                iftError("Unexpected color space", "iftConvertPixelColorSpace");
                break;
        }

        /* convert from YCbCr to the desired color space */
        switch(colSpaceOut) {
            case YCbCr_CSPACE:
                colorOut.val[0] = auxColor.val[0];
                colorOut.val[1] = auxColor.val[1];
                colorOut.val[2] = auxColor.val[2];
                break;

            case YCbCrNorm_CSPACE:
                colorOut.val[0] = (float)auxColor.val[0]/(float)normValue;
                colorOut.val[1] = (float)auxColor.val[1]/(float)normValue;
                colorOut.val[2] = (float)auxColor.val[2]/(float)normValue;
                break;

            case GRAY_CSPACE:
                colorOut.val[0] = auxColor.val[0];
                break;

            case GRAYNorm_CSPACE:
                colorOut.val[0] = (float)auxColor.val[0]/(float)normValue;
                break;

            case RGB_CSPACE:
                auxColor = iftYCbCrtoRGB(auxColor, normValue);
                colorOut.val[0] = auxColor.val[0];
                colorOut.val[1] = auxColor.val[1];
                colorOut.val[2] = auxColor.val[2];
                break;

            case RGBNorm_CSPACE:
                auxColor = iftYCbCrtoRGB(auxColor, normValue);
                colorOut.val[0] = (float)auxColor.val[0]/(float)normValue;
                colorOut.val[1] = (float)auxColor.val[1]/(float)normValue;
                colorOut.val[2] = (float)auxColor.val[2]/(float)normValue;
                break;

            case HSV_CSPACE:
                auxColor = iftYCbCrtoRGB(auxColor, normValue);
                auxColor = iftRGBtoHSV(auxColor, normValue);
                colorOut.val[0] = auxColor.val[0];
                colorOut.val[1] = auxColor.val[1];
                colorOut.val[2] = auxColor.val[2];
                break;

            case LAB_CSPACE:
                auxColorF = iftYCbCrtoLab(auxColor, normValue);
                colorOut.val[0] = auxColorF.val[0];
                colorOut.val[1] = auxColorF.val[1];
                colorOut.val[2] = auxColorF.val[2];
                break;

            case LABNorm_CSPACE:
                auxColorF = iftYCbCrtoLabNorm(auxColor, normValue);
                colorOut.val[0] = auxColorF.val[0];
                colorOut.val[1] = auxColorF.val[1];
                colorOut.val[2] = auxColorF.val[2];
                break;

            case LABNorm2_CSPACE:
                auxColorF = iftYCbCrtoLabNorm2(auxColor, normValue);
                colorOut.val[0] = auxColorF.val[0];
                colorOut.val[1] = auxColorF.val[1];
                colorOut.val[2] = auxColorF.val[2];
                break;
                
            default:
                iftError("Unexpected color space", "iftConvertPixelColorSpace");
                break;
        }
    }

    return colorOut;
}

iftColor iftRGBtoYCbCr(iftColor cin, int normalization_value)
{
    iftColor cout;
    float a = (16.0/255.0)*(float)normalization_value;
    float b = (128.0/255.0)*(float)normalization_value;

    cout.val[0]=(int)(0.256789062*(float)cin.val[0]+
                      0.504128906*(float)cin.val[1]+
                      0.09790625*(float)cin.val[2]+a);
    cout.val[1]=(int)(-0.148222656*(float)cin.val[0]+
                      -0.290992187*(float)cin.val[1]+
                      0.439214844*(float)cin.val[2]+b);
    cout.val[2]=(int)(0.439214844*(float)cin.val[0]+
                      -0.367789063*(float)cin.val[1]+
                      -0.071425781*(float)cin.val[2]+b);

    for(int i=0; i < 3; i++) {
        if (cout.val[i] < 0) cout.val[i] = 0;
        if (cout.val[i] > normalization_value) cout.val[i] = normalization_value;
    }

    return(cout);
}

iftColor iftYCbCrtoRGB(iftColor cin, int normalization_value)
{
    iftColor cout;
    float a = (16.0/255.0)*(float)normalization_value;
    float b = (128.0/255.0)*(float)normalization_value;

    cout.val[0]=(int)(1.164383562*((float)cin.val[0]-a)+
                      1.596026786*((float)cin.val[2]-b));

    cout.val[1]=(int)(1.164383562*((float)cin.val[0]-a)+
                      -0.39176229*((float)cin.val[1]-b)+
                      -0.812967647*((float)cin.val[2]-b));

    cout.val[2]=(int)(1.164383562*((float)cin.val[0]-a)+
                      2.017232143*((float)cin.val[1]-b));

    for(int i=0; i < 3; i++) {
        if (cout.val[i] < 0) cout.val[i] = 0;
        if (cout.val[i] > normalization_value) cout.val[i] = normalization_value;
    }

    return(cout);
}

/**
 * Convert a RGB color value to YCbCr BT.2020 color space.
 *
 * http://http://en.wikipedia.org/wiki/Rec._2020
 *
 * Book Digital Video Concepts, Methods, and Metrics: Quality, Compression, Performance,
 * and Power Trade-off Analysis, ISBN 9781430267133, page 27
 * https://books.google.com.br/books?id=fr08BQAAQBAJ
 *
 * @param cin Input RGB color
 * @param rgbBitDepth Input bit depth
 * @param yCbCrBitDepth Output bit depth (10 or 12)
 * @return Color value in YCbCr space
 **/
iftColor iftRGBtoYCbCrBT2020(iftColor cin, const int rgbBitDepth, const int yCbCrBitDepth)
{
    int minLum, minChr, quantLum, quantChr;
    iftColor cout;

    switch (yCbCrBitDepth) {
        case 8:
            minLum = 16; // 16 * 2^(bitDepth-8)
            minChr = 128; // 128 * 2^(bitDepth-8)
            quantLum = 219.0; // 219 * 2^(bitDepth-8)
            quantChr = 224.0;  // 224 * 2^(bitDepth-8)
            break;
        case 10:
            minLum = 64; // 16 * 2^(bitDepth-8)
            minChr = 512; // 128 * 2^(bitDepth-8)
            quantLum = 876; // 219 * 2^(bitDepth-8)
            quantChr = 896;  // 224 * 2^(bitDepth-8)
            break;
        case 12:
            minLum = 256; // 16 * 2^(bitDepth-8)
            minChr = 2048; // 128 * 2^(bitDepth-8)
            quantLum = 3504; // 219 * 2^(bitDepth-8)
            quantChr = 3584;  // 224 * 2^(bitDepth-8)
            break;
        case 16:
            minLum = 4096; // 16 * 2^(bitDepth-8)
            minChr = 32768; // 128 * 2^(bitDepth-8)
            quantLum = 56064.0; // 219 * 2^(bitDepth-8)
            quantChr = 57344.0;  // 224 * 2^(bitDepth-8)
            break;
        default:
            iftError("Bit depth not specified in BT.2020", "iftRGBtoYCbCrBT2020");
            cout.val[0] = cout.val[1] = cout.val[2] = 0;
            return cout;
    }

    double maxRgbValue = (double) ((1 << rgbBitDepth) - 1);
    double r = cin.val[0] / maxRgbValue;
    double g = cin.val[1] / maxRgbValue;
    double b = cin.val[2] / maxRgbValue;

    double y = 0.2627 * r + 0.6780 * g + 0.0593 * b;
    double cb = (b - y) / 1.8814;
    double cr = (r - y) / 1.4746;

    // clip luminance to [0..1] and chrominance to [-0.5..0.5]
    if (y < 0.0) y = 0.0;
    else if (y > 1.0) y = 1.0;
    if (cb < -0.5) cb = -0.5;
    else if (cb > 0.5) cb = 0.5;
    if (cr < -0.5) cr = -0.5;
    else if (cr > 0.5) cr = 0.5;

    // perform quantization
    cout.val[0] = (int) (y * quantLum) + minLum;
    cout.val[1] = (int) (cb * quantChr) + minChr;
    cout.val[2] = (int) (cr * quantChr) + minChr;

    return cout;
}

/**
 * Convert a YCbCr BT.2020 color value to RGB color space.
 *
 * http://http://en.wikipedia.org/wiki/Rec._2020
 *
 * Book Digital Video Concepts, Methods, and Metrics: Quality, Compression, Performance,
 * and Power Trade-off Analysis, ISBN 9781430267133, page 27
 * https://books.google.com.br/books?id=fr08BQAAQBAJ
 *
 * @param cin Input YCbCr color
 * @param yCbCrBitDepth Input bit depth (10 or 12)
 * @param rgbBitDepth Output bit depth
 * @return Color value in RGB space
 **/
iftColor iftYCbCrBT2020toRGB(iftColor cin, const int yCbCrBitDepth, const int rgbBitDepth)
{
    int minLum, minChr;
    double quantLum, quantChr;
    iftColor cout;

    switch (yCbCrBitDepth) {
        case 8:
            minLum = 16; // 16 * 2^(bitDepth-8)
            minChr = 128; // 128 * 2^(bitDepth-8)
            quantLum = 219.0; // 219 * 2^(bitDepth-8)
            quantChr = 224.0;  // 224 * 2^(bitDepth-8)
            break;
        case 10:
            minLum = 64; // 16 * 2^(bitDepth-8)
            minChr = 512; // 128 * 2^(bitDepth-8)
            quantLum = 876.0; // 219 * 2^(bitDepth-8)
            quantChr = 896.0;  // 224 * 2^(bitDepth-8)
            break;
        case 12:
            minLum = 256; // 16 * 2^(bitDepth-8)
            minChr = 2048; // 128 * 2^(bitDepth-8)
            quantLum = 3504.0; // 219 * 2^(bitDepth-8)
            quantChr = 3584.0;  // 224 * 2^(bitDepth-8)
            break;
        case 16:
            minLum = 4096; // 16 * 2^(bitDepth-8)
            minChr = 32768; // 128 * 2^(bitDepth-8)
            quantLum = 56064.0; // 219 * 2^(bitDepth-8)
            quantChr = 57344.0;  // 224 * 2^(bitDepth-8)
            break;
        default:
            iftError("Bit depth not specified in BT.2020", "iftYCbCrBT2020toRGB");
            cout.val[0] = cout.val[1] = cout.val[2] = 0;
            return cout;
    }

    double y = (cin.val[0] - minLum) / quantLum;
    double cb = (cin.val[1] - minChr) / quantChr;
    double cr = (cin.val[2] - minChr) / quantChr;

    double r = cr * 1.4746 + y;
    double b = cb * 1.8814 + y;
    double g = (y - 0.2627 * r - 0.0593 * b) / 0.6780;

    // clip rgb values to [0..1]
    if (r < 0.0) r = 0.0;
    else if (r > 1.0) r = 1.0;
    if (g < 0.0) g = 0.0;
    else if (g > 1.0) g = 1.0;
    if (b < 0.0) b = 0.0;
    else if (b > 1.0) b = 1.0;

    // perform quantization
    double maxRgbValue = (double) ((1 << rgbBitDepth) - 1);
    cout.val[0] = (int) (r * maxRgbValue);
    cout.val[1] = (int) (g * maxRgbValue);
    cout.val[2] = (int) (b * maxRgbValue);

    return cout;
}

iftFColor iftYCbCrtoLab(iftColor ycbcr, int normalization_value) {
    iftColor rgb = iftYCbCrtoRGB(ycbcr, normalization_value);

    return iftRGBtoLab(rgb, normalization_value);
}

iftFColor iftYCbCrtoLabNorm2(iftColor ycbcr, int normalization_value) {
    iftColor rgb = iftYCbCrtoRGB(ycbcr, normalization_value);

    return iftRGBtoLabNorm2(rgb, normalization_value);
}

iftFColor iftRGBtoLabNorm(iftColor rgb, int normalization_value)
{
    //RGB to XYZ

    float R = rgb.val[0]/(float)normalization_value;
    float G = rgb.val[1]/(float)normalization_value;
    float B = rgb.val[2]/(float)normalization_value;

    if(R <= 0.04045)	R = R/12.92;
    else	        R = pow((R+0.055)/1.055,2.4);

    if(G <= 0.04045)	G = G/12.92;
    else		G = pow((G+0.055)/1.055,2.4);

    if(B <= 0.04045)	B = B/12.92;
    else		B = pow((B+0.055)/1.055,2.4);

    float X = (0.4123955889674142161*R + 0.3575834307637148171*G + 0.1804926473817015735*B);
    float Y = (0.2125862307855955516*R + 0.7151703037034108499*G + 0.07220049864333622685*B);
    float Z = (0.01929721549174694484*R + 0.1191838645808485318*G + 0.9504971251315797660*B);

    //XYZ to lab
    X /= WHITEPOINT_X;
    Y /= WHITEPOINT_Y;
    Z /= WHITEPOINT_Z;
    X = LABF(X);
    Y = LABF(Y);
    Z = LABF(Z);
    float L = 116*Y - 16;
    float a = 500*(X - Y);
    float b = 200*(Y - Z);

    iftFColor lab;
    lab.val[0] = L;
    lab.val[1] = a;
    lab.val[2] = b;

    return lab;
}

iftFColor iftYCbCrtoLabNorm(iftColor ycbcr, int normalization_value) {
    iftColor rgb = iftYCbCrtoRGB(ycbcr, normalization_value);

    return iftRGBtoLabNorm(rgb, normalization_value);
}

iftColor iftLabtoRGB(iftFColor lab, int normalization_value)
{
    //Lab to XYZ
    float L = lab.val[0];
    float a = lab.val[1];
    float b = lab.val[2];

    L = (L + 16)/116;
    a = L + a/500;
    b = L - b/200;

    float X = WHITEPOINT_X*LABINVF(a);
    float Y = WHITEPOINT_Y*LABINVF(L);
    float Z = WHITEPOINT_Z*LABINVF(b);

    //XYZ to RGB
    float R = ( 3.2406*X - 1.5372*Y - 0.4986*Z);
    float G = (-0.9689*X + 1.8758*Y + 0.0415*Z);
    float B = ( 0.0557*X - 0.2040*Y + 1.0570*Z);

    iftColor rgb;
    rgb.val[0] = R*(float)normalization_value;
    rgb.val[1] = G*(float)normalization_value;
    rgb.val[2] = B*(float)normalization_value;

    return rgb;
}

iftColor iftLabtoQLab(iftFColor lab, int normalization_value){
    iftColor cout;

    cout.val[0] = (lab.val[0] / 99.998337)*normalization_value;
    cout.val[1] = ((lab.val[1] + 86.182236)/(86.182236 + 98.258614) )*normalization_value;
    cout.val[2] = ((lab.val[2] + 107.867744)/(107.867744 + 94.481682))*normalization_value;

    return cout;
}

iftFColor iftQLabToLab(iftColor qlab, int normalization_value){
    iftFColor cout;

    cout.val[0] = ((float)qlab.val[0]/normalization_value) * 99.998337;
    cout.val[1] = (((float)qlab.val[1]/normalization_value) * (86.182236 + 98.258614)) -86.182236;
    cout.val[2] = (((float)qlab.val[2]/normalization_value) * (107.867744 + 94.481682)) -107.867744;

    return cout;
}

iftColor iftHexaColorToRGB(const char *hexa_color) {
    if (!iftRegexMatch(hexa_color, "^#[a-fA-F0-9]{6}$"))
        iftError("Invalid hexa color: %s\nTry: #XXXXXX", "iftHexaColorToRGB", hexa_color);
    
    char hexa[3];
    hexa[2] = '\0';

    iftColor RGB;
    
    hexa[0] = hexa_color[1];
    hexa[1] = hexa_color[2];
    RGB.val[0] = (int) strtol(hexa, NULL, 16);

    hexa[0] = hexa_color[3];
    hexa[1] = hexa_color[4];
    RGB.val[1] = (int) strtol(hexa, NULL, 16);
    
    hexa[0] = hexa_color[5];
    hexa[1] = hexa_color[6];
    RGB.val[2] = (int) strtol(hexa, NULL, 16);

    return RGB;
}

char *iftConvertColorToHexaColor(iftColor RGB) {
    char *hexa_color = iftAlloc(7, sizeof(char));

    char channel[3];

    sprintf(channel, "%X", RGB.val[0]);
    if (strlen(channel) == 1) {
        hexa_color[0] = '0';
        hexa_color[1] = channel[0];
    }
    else {
        hexa_color[0] = channel[0];
        hexa_color[1] = channel[1];
    }

    sprintf(channel, "%X", RGB.val[1]);
    if (strlen(channel) == 1) {
        hexa_color[2] = '0';
        hexa_color[3] = channel[0];
    }
    else {
        hexa_color[2] = channel[0];
        hexa_color[3] = channel[1];
    }

    sprintf(channel, "%X", RGB.val[2]);
    if (strlen(channel) == 1) {
        hexa_color[4] = '0';
        hexa_color[5] = channel[0];
    }
    else {
        hexa_color[4] = channel[0];
        hexa_color[5] = channel[1];
    }
    
    hexa_color[6] = '\0';

    return hexa_color;
}

iftColor iftGrayScaleToBlueToRedColor(float intensity, float norm_value)
{
    float value = 4*(intensity/norm_value)+1;

    iftColor rgb_color;
    rgb_color.val[0] = norm_value * iftMax(0,(3-(float)fabs(value-4)-(float)fabs(value-5))/2);
    rgb_color.val[1] = norm_value * iftMax(0,(4-(float)fabs(value-2)-(float)fabs(value-4))/2);
    rgb_color.val[2] = norm_value * iftMax(0,(3-(float)fabs(value-1)-(float)fabs(value-2))/2);

    iftColor ycbcr = iftRGBtoYCbCr(rgb_color, norm_value);

    return ycbcr;
}

iftColorTable *iftBlueToRedColorTable(int ncolors)
{
    iftColorTable *ctb = iftCreateColorTable(ncolors);

    if (ncolors <= 50){

      float dR, dG, dB, val[3];
      int   point[4];

      if (ncolors < 7)
	iftError("It requires at least seven colors", "iftBlueToRedColorTable");

      point[0] =   0;
      point[1] =   ncolors/3;
      point[2] = 2*ncolors/3;
      point[3] =   ncolors-1;
      
      /* Build a blue-to-red lookup table */
      
      /* blue */
      val[0] = ctb->color[point[0]].val[0] = 0;
      val[1] = ctb->color[point[0]].val[1] = 183;
      val[2] = ctb->color[point[0]].val[2] = 255;
      
      /* green */
      ctb->color[point[1]].val[0] = 0;
      ctb->color[point[1]].val[1] = 169;
      ctb->color[point[1]].val[2] = 6;
      
      /* yellow */
      ctb->color[point[2]].val[0] = 255;
      ctb->color[point[2]].val[1] = 255;
      ctb->color[point[2]].val[2] = 0;
      
      /* red */
      ctb->color[point[3]].val[0] = 173;
      ctb->color[point[3]].val[1] = 0;
      ctb->color[point[3]].val[2] = 20;
      
      for (int j=0; j < 3; j++) {
        dR = (float)(ctb->color[point[j+1]].val[0]-ctb->color[point[j]].val[0])/
	  (float)(point[j+1]-point[j]);
        dG = (float)(ctb->color[point[j+1]].val[1]-ctb->color[point[j]].val[1])/
	  (float)(point[j+1]-point[j]);
        dB = (float)(ctb->color[point[j+1]].val[2]-ctb->color[point[j]].val[2])/
	  (float)(point[j+1]-point[j]);
        for (int i = point[j]+1; i < point[j+1]; i++){
	  val[0] += dR;
	  val[1] += dG;
	  val[2] += dB;
	  ctb->color[i].val[0] = (int)(val[0]);
	  ctb->color[i].val[1] = (int)(val[1]);
	  ctb->color[i].val[2] = (int)(val[2]);
        }
      }
    } else {
      for (int i=0; i < ncolors; i++){
	ctb->color[i] = iftGrayScaleToBlueToRedColor(i, ncolors);
      }
    }

    return(ctb);
}

iftColor iftRGBtoHSV(iftColor cin, int normalization_value) {
    float r = ((float)cin.val[0]/normalization_value),
            g = ((float)cin.val[1]/normalization_value),
            b = ((float)cin.val[2]/normalization_value), v, x, f;
    float a[3];
    int   i;
    iftColor cout;

    // RGB are each on [0, 1]. S and V are returned on [0, 1] and H is
    // returned on [0, 6].

    x = iftMin(iftMin(r, g), b);
    v = iftMax(iftMax(r, g), b);
    if (v == x) {
        a[0]=0.0;
        a[1]=0.0;
        a[2]=v;
    } else {
        f = (r == x) ? g - b : ((g == x) ? b - r : r - g);
        i = (r == x) ? 3 : ((g == x) ? 5 : 1);
        a[0]=((float)i)-f/(v-x);
        a[1]=(v-x)/v;
        a[2]=0.299*r+0.587*g+0.114*b;
    }

    // (un)normalize

    cout.val[0] = (int)(a[0]*60.0);
    cout.val[1] = (int)(a[1]*normalization_value);
    cout.val[2] = (int)(a[2]*normalization_value);

    return(cout);
}

iftColor iftHSVtoRGB(iftColor cin, int normalization_value) {
    // H is given on [0, 6]. S and V are given on [0, 1].
    // RGB are each returned on [0, 1].
    float h = ((float)cin.val[0]/60.0),
            s = ((float)cin.val[1]/normalization_value),
            v = ((float)cin.val[2]/normalization_value), m, n, f;
    float a[3]={0,0,0};
    int i;
    iftColor cout;

    if (s==0.0) {
        a[0]=a[1]=a[2]=v;
    } else {
        i = (int) floor(h);
        f = h - (float)i;
        if(!(i & 1)) f = 1 - f; // if i is even
        m = v * (1 - s);
        n = v * (1 - s * f);
        switch (i) {
            case 6:
            case 0: a[0]=v; a[1]=n; a[2]=m; break;
            case 1: a[0]=n; a[1]=v; a[2]=m; break;
            case 2: a[0]=m; a[1]=v; a[2]=n; break;
            case 3: a[0]=m; a[1]=n; a[2]=v; break;
            case 4: a[0]=n; a[1]=m; a[2]=v; break;
            case 5: a[0]=v; a[1]=m; a[2]=n; break;
        }
    }

    // (un)normalize
    for(i=0;i<3;i++)
        cout.val[i]=a[i]*normalization_value;

    return(cout);
}

iftFColor iftRGBTo(iftColor RGB, int normalization_value, iftColorSpace cspace) {
    iftConvertColorFunc convert_func = NULL;

    if (cspace == RGB_CSPACE)
        return iftColorToFColor(RGB);
    else {
        if (cspace == YCbCr_CSPACE)
            convert_func = iftRGBtoYCbCr;
        else if (cspace == HSV_CSPACE)
            convert_func = iftRGBtoHSV;
        else iftError("Required Color Space not supported yet!", "iftRGBTo");
    
        return iftColorToFColor(convert_func(RGB, normalization_value));
    }
}

/**
 * @brief Converts values from [0,1] to R, G, B colors according to the heat color map (from blue to red). 
 * @param value: input value 
 * @param R: output red value
 * @param G: output green value 
 * @param B: output blue value 
 * @author Alexandre Falcao
 */

void iftHeatColorMapping(float value, float *R, float *G, float *B)
{
    value = (6-2)*value+1;
    *R = iftMax(0.0f,(3-(float)fabs(value-4)-(float)fabs(value-5))/2);
    *G = iftMax(0.0f,(4-(float)fabs(value-2)-(float)fabs(value-4))/2);
    *B = iftMax(0.0f,(3-(float)fabs(value-1)-(float)fabs(value-2))/2);
}

iftColorTable *iftCreateGradientColorTable(const iftStrArray *hexa_colors, int n_colors) {    
    int n_intervals = hexa_colors->n - 1;
    int interval_size = ceil(((float) n_colors) / n_intervals);

    iftColorTable *ctb = iftCreateColorTable(n_colors);
    int c = 0; // index for populating the resulting color table

    for (int i = 0; i < n_intervals; i++) {
        char *begin_hexa_color = hexa_colors->val[i];
        char *end_hexa_color = hexa_colors->val[i+1];

        iftColor begin_RGB = iftHexaColorToRGB(begin_hexa_color);
        iftColor end_RGB = iftHexaColorToRGB(end_hexa_color);

        int n_colors_interval = interval_size;

        // if it is in the last iteration
        if (i == (n_intervals-1)) {
            int last_interval_size = (n_colors - c); // == (n_colors-1) - c +1

            if (interval_size != last_interval_size)
                n_colors_interval = last_interval_size;
        }

        iftColorTable *ctb_interval = _iftCreateRGBGradientColorTable(begin_RGB, end_RGB, n_colors_interval);
        for (int j = 0; j < ctb_interval->ncolors; j++)
            ctb->color[c++] = ctb_interval->color[j];
        iftDestroyColorTable(&ctb_interval);

        c--;
    }

    // #pragma omp parallel for
    // for (int c = 0; c < ctb->ncolors; c++)
    //     ctb->color[c] = iftRGBtoYCbCr(ctb->color[c], 255);

    return ctb;
}

iftColorTable *iftGrayColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(2);

    strcpy(hexa_colors->val[0], "#000000");
    strcpy(hexa_colors->val[1], "#FFFFFF");

    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);

    iftDestroyStrArray(&hexa_colors);

    return ctb;
}

iftColorTable *iftIronColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(2);

    strcpy(hexa_colors->val[0], "#FF0000");
    strcpy(hexa_colors->val[1], "#FFFF00");

    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);

    iftDestroyStrArray(&hexa_colors);

    return ctb;
}

iftColorTable *iftHotColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(5);

    strcpy(hexa_colors->val[0], "#000000");
    strcpy(hexa_colors->val[1], "#8b0000");
    strcpy(hexa_colors->val[2], "#ffa500");
    strcpy(hexa_colors->val[3], "#ffff00");
    strcpy(hexa_colors->val[4], "#ffffff");

    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);

    iftDestroyStrArray(&hexa_colors);

    return ctb;    
}

iftColorTable *iftRainbowColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(6);

    strcpy(hexa_colors->val[0], "#FF0000");
    strcpy(hexa_colors->val[1], "#ffed00");
    strcpy(hexa_colors->val[2], "#00FF00");
    strcpy(hexa_colors->val[3], "#00ffea");
    strcpy(hexa_colors->val[4], "#0020ff");
    strcpy(hexa_colors->val[5], "#c800ff");

    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, iftMax(n_colors,6));

    iftDestroyStrArray(&hexa_colors);

    return ctb;
}

iftColorTable *iftCategoricalColorTable(int n_colors) {
    if (n_colors <= 21)
        return iftCategory21ColorTable(n_colors);
    else {
        iftColorTable *ctb_categ = iftCategory21ColorTable(21);

        iftStrArray *hexa_colors = iftCreateStrArray(7);

        strcpy(hexa_colors->val[0], "#00FFFF");
        strcpy(hexa_colors->val[1], "#FF00FF");
        strcpy(hexa_colors->val[2], "#FFFF00");
        strcpy(hexa_colors->val[3], "#00FF00");
        strcpy(hexa_colors->val[4], "#FF8000");
        strcpy(hexa_colors->val[5], "#008080");
        strcpy(hexa_colors->val[6], "#800080");

        int n_remaining_colors = n_colors - 21;
        iftColorTable *ctb_grad = iftCreateGradientColorTable(hexa_colors, n_remaining_colors);

        iftIntArray *indices = iftIntRange(0, n_remaining_colors-1, 1);
        iftShuffleIntArray(indices->val, indices->n);


        iftColorTable *ctb = iftCreateColorTable(n_colors);
        for (int c = 0; c < 21; c++)
            ctb->color[c] = ctb_categ->color[c];
        
        for (int c = 0; c < n_remaining_colors; c++)
            ctb->color[c + 21] = ctb_grad->color[indices->val[c]];

        iftDestroyStrArray(&hexa_colors);
        iftDestroyIntArray(&indices);
        iftDestroyColorTable(&ctb_grad);
        iftDestroyColorTable(&ctb_categ);

        return ctb;
    }
}

iftColorTable *iftReverseRainbowColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(2);

    strcpy(hexa_colors->val[0], "#cd00ff");
    strcpy(hexa_colors->val[1], "#FF0000");

    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);

    iftDestroyStrArray(&hexa_colors);

    return ctb;    
}

iftColorTable *iftHeatMapColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(5);

    strcpy(hexa_colors->val[0], "#0000FF");
    strcpy(hexa_colors->val[1], "#00FFFF");
    strcpy(hexa_colors->val[2], "#00FF00");
    strcpy(hexa_colors->val[3], "#FFFF00");
    strcpy(hexa_colors->val[4], "#FF0000");

    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);

    iftDestroyStrArray(&hexa_colors);

    return ctb;    
}

iftColorTable *iftRedHotColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(3);
    
    strcpy(hexa_colors->val[0], "#000000");
    strcpy(hexa_colors->val[1], "#FF0000");
    strcpy(hexa_colors->val[2], "#FFFFFF");
    
    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);
    
    iftDestroyStrArray(&hexa_colors);
    
    return ctb;    
}

iftColorTable *iftGreenHotColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(3);
    
    strcpy(hexa_colors->val[0], "#000000");
    strcpy(hexa_colors->val[1], "#00FF00");
    strcpy(hexa_colors->val[2], "#FFFFFF");
    
    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);
    
    iftDestroyStrArray(&hexa_colors);

    return ctb;    
}

iftColorTable *iftBlueHotColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(3);
    
    strcpy(hexa_colors->val[0], "#000000");
    strcpy(hexa_colors->val[1], "#0000FF");
    strcpy(hexa_colors->val[2], "#FFFFFF");
    
    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);
    
    iftDestroyStrArray(&hexa_colors);

    return ctb;    
}

iftColorTable *iftRedYellowHotColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(4);
    
    strcpy(hexa_colors->val[0], "#000000");
    strcpy(hexa_colors->val[1], "#FF0000");
    strcpy(hexa_colors->val[2], "#FFFF00");
    strcpy(hexa_colors->val[3], "#FFFFFF");
    
    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);
    
    iftDestroyStrArray(&hexa_colors);

    return ctb;    
}

iftColorTable *iftRedToBlueColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(5);

    strcpy(hexa_colors->val[0], "#a60126");
    strcpy(hexa_colors->val[1], "#fb9c59");
    strcpy(hexa_colors->val[2], "#ecebd0");
    strcpy(hexa_colors->val[3], "#78b1d3");
    strcpy(hexa_colors->val[4], "#43489e");

    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);
    
    iftDestroyStrArray(&hexa_colors);

    return ctb;    
}


iftColorTable *iftViridisColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(8);

    strcpy(hexa_colors->val[0], "#480054");
    strcpy(hexa_colors->val[1], "#4f307f");
    strcpy(hexa_colors->val[2], "#435b8d");
    strcpy(hexa_colors->val[3], "#347f8e");
    strcpy(hexa_colors->val[4], "#22a287");
    strcpy(hexa_colors->val[5], "#3dc36c");
    strcpy(hexa_colors->val[6], "#93db35");
    strcpy(hexa_colors->val[7], "#f3e91c");

    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);
    
    iftDestroyStrArray(&hexa_colors);

    return ctb;    
}


iftColorTable *iftPlasmaColorTable(int n_colors) {
    iftStrArray *hexa_colors = iftCreateStrArray(8);

    strcpy(hexa_colors->val[0], "#2f0087");
    strcpy(hexa_colors->val[1], "#6200a4");
    strcpy(hexa_colors->val[2], "#9200a6");
    strcpy(hexa_colors->val[3], "#ba2f8a");
    strcpy(hexa_colors->val[4], "#d85b69");
    strcpy(hexa_colors->val[5], "#ee8949");
    strcpy(hexa_colors->val[6], "#f6bd27");
    strcpy(hexa_colors->val[7], "#e4fa15");

    iftColorTable *ctb = iftCreateGradientColorTable(hexa_colors, n_colors);
    
    iftDestroyStrArray(&hexa_colors);

    return ctb;    
}


iftColorTable *iftCategory21ColorTable(int n_colors) {
    if (n_colors > 21)
        iftError("Number of colors for this color table must be <= 21. Input ncolors was %d",
                  "iftCategory21ColorTable", n_colors);

    iftStrArray *hexa_colors = iftCreateStrArray(21);

    strcpy(hexa_colors->val[0], "#e41a1c");
    strcpy(hexa_colors->val[1], "#4daf4a");
    strcpy(hexa_colors->val[2], "#377eb8");
    strcpy(hexa_colors->val[3], "#984ea3");
    strcpy(hexa_colors->val[4], "#ff7f00");
    strcpy(hexa_colors->val[5], "#ffff33");
    strcpy(hexa_colors->val[6], "#a65628");
    strcpy(hexa_colors->val[7], "#f781bf");
    strcpy(hexa_colors->val[8], "#999999");
    strcpy(hexa_colors->val[9], "#bcbc35");
    strcpy(hexa_colors->val[10], "#29bece");
    strcpy(hexa_colors->val[11], "#fd9898");
    strcpy(hexa_colors->val[12], "#afc8e7");
    strcpy(hexa_colors->val[13], "#9ade8d");
    strcpy(hexa_colors->val[14], "#c5b1d4");
    strcpy(hexa_colors->val[15], "#fdbb7d");
    strcpy(hexa_colors->val[16], "#dbda91");
    strcpy(hexa_colors->val[17], "#c39c95");
    strcpy(hexa_colors->val[18], "#f6b7d2");
    strcpy(hexa_colors->val[19], "#c7c7c7");
    strcpy(hexa_colors->val[20], "#a0dae4");

    iftColorTable *ctb = iftCreateColorTable(n_colors);
    for (int c = 0; c < n_colors; c++) {
        ctb->color[c] = iftRGBtoYCbCr(iftHexaColorToRGB(hexa_colors->val[c]), 255);
    }
    iftDestroyStrArray(&hexa_colors);

    return ctb;
}

iftColorTable *iftCategory10ColorTable(int n_colors) {
    if (n_colors > 10)
        iftError("Number of colors for this color table must be <= 10. Input ncolors was %d",
                  "iftCategory10ColorTable", n_colors);

    iftStrArray *hexa_colors = iftCreateStrArray(10);

    strcpy(hexa_colors->val[0], "#1f77b4");
    strcpy(hexa_colors->val[1], "#ff7f0e");
    strcpy(hexa_colors->val[2], "#2ca02c");
    strcpy(hexa_colors->val[3], "#d62728");
    strcpy(hexa_colors->val[4], "#9467bd");
    strcpy(hexa_colors->val[5], "#8c564b");
    strcpy(hexa_colors->val[6], "#e377c2");
    strcpy(hexa_colors->val[7], "#7f7f7f");
    strcpy(hexa_colors->val[8], "#bcbd22");
    strcpy(hexa_colors->val[9], "#17becf");

    iftColorTable *ctb = iftCreateColorTable(n_colors);
    for (int c = 0; c < n_colors; c++)
        ctb->color[c] = iftRGBtoYCbCr(iftHexaColorToRGB(hexa_colors->val[c]), 255);
    iftDestroyStrArray(&hexa_colors);

    return ctb;
}
