/**
 * @file iftColor.h
 * @brief Structs and function prototypes for color.
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 * @note Programs:
 * @ref iftGenerateColorTable.c (demo/Miscellaneous/iftGenerateColorTable.c) = Generates a Color Table. 
 */

#ifndef IFT_COLOR_H_
#define IFT_COLOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/StrArray.h"
#include "iftCommon.h"


/**
 * @brief Color space enumeration
 * @note added by Adan Echemendia. The three color components from LABNorm2_CSPACE are normalized from 0 to 1. Not so
 * the components from LABNorm_CSPACE.
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
//! swig()
typedef enum ift_color_space {
    YCbCr_CSPACE,
    YCbCrNorm_CSPACE,
    RGB_CSPACE,
    RGBNorm_CSPACE,
    GRAY_CSPACE,
    GRAYNorm_CSPACE,
    WEIGHTED_YCbCr_CSPACE,
    LAB_CSPACE,
    LABNorm_CSPACE,
	LABNorm2_CSPACE,
    HSV_CSPACE,
    UNDEFINED_CSPACE
} iftColorSpace;

#define WHITEPOINT_X	0.950456
#define WHITEPOINT_Y	1.0
#define WHITEPOINT_Z	1.088754

#define LABF(t)	\
	((t >= 8.85645167903563082e-3) ? \
	pow(t,0.333333333333333) : (841.0/108.0)*(t) + (4.0/29.0))

#define LABINVF(t)	\
	((t >= 0.206896551724137931) ? \
	((t)*(t)*(t)) : (108.0/841.0)*((t) - (4.0/29.0)))


/**
 * @brief Float color space structure
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
typedef struct ift_fcolor {
	float val[3];
} iftFColor;


/**
 * @brief Integer color space structure
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 */
//! swig(extend = iftColorExt.i)
typedef struct ift_color {
  int val[3];
	float alpha;
} iftColor;


/**
 * @brief Array of colors
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 */
 //! swig(destroyer = iftDestroyColorTable, extend = iftColorTableExt.i)
typedef struct ift_colortable {
  iftColor *color;
  int ncolors;
} iftColorTable;

typedef struct ift_fcolormatrix {
    iftFColor *color;
    int ncolors;
} iftFColorMatrix;

static inline iftFColor iftColorToFColor(iftColor color) {
    iftFColor fcolor;
    fcolor.val[0] = color.val[0];
    fcolor.val[1] = color.val[1];
    fcolor.val[2] = color.val[2];
    return fcolor;
}


/**
 * @brief Converts a iftColorSpace name (str) into its corresponding iftColorSpace value
 * @author Cesar Castelo
 * @date Mar 19, 2019
 * @ingroup Color
 */
iftColorSpace iftColorSpaceStrToColorSpace(char *colSpaceStr);


/**
 * @brief Creates an array by generating colors in the color space YCbCr by random.
 * @author Samuel Martins
 * @date Jun 15, 2016
 * @ingroup Color
 */
 //! swig(newobject, stable)
iftColorTable *iftCreateColorTable(int ncolors);


/**
 * @brief Creates a float color matrix allocating 0 to its value.
 * @author Leonardo de Melo
 * @date Jun 15, 2016
 * @ingroup Color
 */
//! swig(newobject, stable)
iftFColorMatrix *iftCreateFColorMatrix(int ncolors);


iftColorTable *iftCreateHeatMapColorTable(float delta, int scalingFactor);
  
//! swig(newobject)
void iftConvertRGBColorTableToYCbCrColorTable(iftColorTable *colorTable_rgb, int normalization_value);


/**
 * @brief Convert a YCbCr Color Table to a RGB Color Table.
 * @author Samuel Martins
 * @date Aug 10, 2018.
 */
//! swig(newobject)
iftColorTable *iftConvertYCbCrColorTableToRGBColorTable(const iftColorTable *ctb, int normalization_value);

/**
 * @brief Destroy color array from memory
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
void iftDestroyColorTable(iftColorTable **ctb);

/**
 * @brief Destroy float color matrix from memory
 *
 * @author Leonardo de Melo
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
void iftDestroyFColorMatrix(iftFColorMatrix **cm);

/**
 * @brief Create a color array from blue to red
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
//! swig(newobject)
iftColorTable *iftBlueToRedColorTable(int ncolors);



/**
 * @brief Create a color array with random YCbCr colors
 * 
 * @author  Alexandre Falcao
 * @param   number of colors in the array
 * @date    Nov, 2020
 * @ingroup Color
 *
 */
//! swig(newobject)
iftColorTable *iftCreateRandomColorTable(int n_colors);

/**
 * @brief Converts an RGB color to the iftColor data structure (does not convert it to YCbCr).
 *
 * @author Thiago Vallin Spina
 * @date Feb 14, 2016
 *
 * @param R The Red component.
 * @param G The Green component.
 * @param B The Blue component.
 *
 * @return The iftColor structure with the RGB color.
 */
  iftColor  iftRGBColor(int R, int G, int B);


/**
 * @brief Converts a grayscale value into a RGB color in a BlueToRed scale.
 *
 * @author Alexandre Falcão
 * @date Mar 10th, 2021
 *
 * @param grayscale (intensity) value.
 * @param normalization value (maximum in the scale).
 *
 * @return The iftColor structure with the RGB color.
 */
  iftColor iftGrayScaleToBlueToRedColor(float intensity, float norm_value);

/**
 * @brief Converts color values between color spaces
 * @warning It receives and returns an iftFColor struct since some color spaces use float values
 * @param colorIn Input color
 * @param colSpaceIn Input color space
 * @param colSpaceOut Output color space
 * @param normValue Normalization value
 *
 * @author Cesar Castelo
 * @date Jan 14, 2019
 * @ingroup Color
 */
iftFColor iftConvertPixelColorSpace(iftFColor colorIn, char colSpaceIn, char colSpaceOut, int normValue);

/**
 * @brief Convert color from RGB to YCbCr
 * https://en.wikipedia.org/wiki/YCbCr (ITU-R BT.601 conversion)
 *
 * @author Alexandre Falcao, modified by Cesar Castelo
 * @date Jul 30, 2016 (Nov 05, 2018)
 * @ingroup Color
 *
 */
//! swig()
iftColor iftRGBtoYCbCr(iftColor cin, int normalization_value);

/**
 * @brief Convert color from YCbCr to RGB
 * https://en.wikipedia.org/wiki/YCbCr (ITU-R BT.601 conversion)
 *
 * @author Alexandre Falcao, modified by Cesar Castelo
 * @date Jun 10, 2016 (Nov 05, 2018)
 * @ingroup Color
 *
 */
//! swig()
iftColor iftYCbCrtoRGB(iftColor cin, int normalization_value);

/**
 * @brief Convert color from RGB to HSV
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftRGBtoHSV(iftColor cin, int normalization_value);

/**
 * @brief Convert color from HSV to RGB
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftHSVtoRGB(iftColor cin, int normalization_value);


/**
 * @brief Convert a RGB color <RGB> to the color space <cspace> and return a iftFColor since some
 *        color spaces can return float values.
 * @author Samuka Martins
 * @date Jun 12, 2018
 */
iftFColor iftRGBTo(iftColor RGB, int normalization_value, iftColorSpace cspace);


// iftFColor iftYCbCrTo(iftColor )



/**
 * @brief Convert color from RGB to Lab
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
static inline iftFColor iftRGBtoLab(iftColor rgb, int normalization_value)
{
        /* with no normalization. */
    // L in [0, 99.998337]
    // a in [-86.182236, 98.258614]
    // b in [-107.867744, 94.481682]

    //RGB to XYZ
    float R = rgb.val[0]/(float)normalization_value;
    float G = rgb.val[1]/(float)normalization_value;
    float B = rgb.val[2]/(float)normalization_value;

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

/**
 * @brief Convert color from RGB to Normalized Lab. The function calls iftRGBtoLab function and normalizes each value to get in the interval [0..1]
 * @author Adan Echemendia
 * @date Jul 26, 2017
 * @ingroup Color
 *
 */
static inline iftFColor iftRGBtoLabNorm2(iftColor rgb, int normalization_value)
{
    /* get lab values*/
    iftFColor lab=iftRGBtoLab(rgb,normalization_value);
    /*normalize each value*/
    lab.val[0]=lab.val[0]/99.998337f;
    lab.val[1]=(lab.val[1]+86.182236f)/(86.182236f+98.258614f);
    lab.val[2]=(lab.val[2]+107.867744f)/(107.867744f+94.481682f);

    return lab;
}

/**
 * @brief Convert color from YcbCr to Lab
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */  
iftFColor iftYCbCrtoLab(iftColor ycbcr, int normalization_value);

/**
 * @brief Convert color from YCbCr to LabNorm
 *
 * @author Cesar Castelo
 * @date Jan 10, 2018
 * @ingroup Color
 *
 */
iftFColor iftYCbCrtoLabNorm(iftColor ycbcr, int normalization_value);

/**
 * @brief Convert color from YcbCr to LabNorm (all components are in the interval [0..1])
 * @author Adan Echemendia
 * @date Jul, 2017
 * @ingroup Color
 *
 */
iftFColor iftYCbCrtoLabNorm2(iftColor ycbcr, int normalization_value);
/**
 * @brief Convert color from RGB to LabNorm
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftFColor iftRGBtoLabNorm(iftColor rgb, int normalization_value);

/**
 * @brief Convert color from Lab to RGB
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftLabtoRGB(iftFColor lab, int normalization_value);

/**
 * @brief Convert color from RGB to YCbCr BT.2020
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
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftRGBtoYCbCrBT2020(iftColor cin, const int rgbBitDepth, const int yCbCrBitDepth);

/**
 * @brief Convert color from YCbCr BT.2020 to RGB
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
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftYCbCrBT2020toRGB(iftColor cin, const int yCbCrBitDepth, const int rgbBitDepth);

/**
 * @brief Convert color from Lab to QLab
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftLabtoQLab(iftFColor lab,int normalization_value);

/**
 * @brief Convert color from QLab to Lab
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftFColor iftQLabToLab(iftColor qlab, int normalization_value);


/**
 * @brief Convert a hexadecimal color (#XXXXXX) to a RGB iftColor.
 * @author Samuka Martins
 * @date Jun 12, 2018
 */
//! swig()
iftColor iftHexaColorToRGB(const char *hexa_color);

/**
 * @brief Convert a RGB color to hexadecimal color (#XXXXXX).
 * @author Samuka Martins
 * @date Jun 12, 2018
 */
char *iftConvertColorToHexaColor(iftColor RGB);


/**
 * @brief Converts values from [0,1] to R, G, B colors according to the heat color map (from blue to red). 
 * @param value: input value 
 * @param R: output red value
 * @param G: output green value 
 * @param B: output blue value 
 * @author Alexandre Falcao
 * @date Aug, 5th, 2016
 * @ingroup Color
 */
void iftHeatColorMapping(float value, float *R, float *G, float *B);


/**
 * @brief Create a color table with <n_colors> colors of 3*8 bits in the color space YCbCr
 * resulting from the gradient of the hexa colors.
 *
 * The parameter <hexa_colors> stores a list of hexa colors which will be our "control colors" for the
 * gradient computation.
 * The color table will have n_colors-1 intervals of gradient colors, so that each interval has
 * the same number of colors.
 * 
 * Eg: hexa_colors = ["#000000", "#FF0000", "#FFFFFF"], n_colors = 15
 * 
 * n_control_colors = 3
 * number of intervals = 3 - 1 = 2
 * number of colors per interval (gradient) = ceil(15 / 2) = 8
 *
 * The color table will have two intervals of gradient colors.
 * The first one with colors from the gradient from [0] = #000000 (black) to [8] = #FF0000 (red)
 * The second one with colors from the gradient from [8] = #FF0000 (red) to [14] = #FFFFFF (white).
 *
 * The last interval can have a different number of colors compared to the others, according to the number
 * of control colors and required n_colors.
 *
 * From this function, it is possible to create several color tables.
 * Good site for examples are:
 * https://jpgraph.net/download/manuals/chunkhtml/ch22s08.html
 * https://www.slicer.org/wiki/Documentation/4.0/Modules/Colors
 * 
 * @param  hexa_colors Array of strings with the control colors in hexa decimal: #XXXXXX
 * @param  n_colors    Required number of colors for the resulting color table.
 * @return             Resulting color table of gradient colors.
 *
 * @author Samuka Martins
 * @date Jun 12, 2018.
 */
//! swig(newobject)
iftColorTable *iftCreateGradientColorTable(const iftStrArray *hexa_colors, int n_colors);


/**
 * @brief Create a Gray Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 12, 2018
 */
//! swig(newobject)
iftColorTable *iftGrayColorTable(int n_colors);

/**
 * @brief Create an Iron Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 12, 2018
 */
//! swig(newobject)
iftColorTable *iftIronColorTable(int n_colors);

/**
 * @brief Create a Hot Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 12, 2018
 */
//! swig(newobject)
iftColorTable *iftHotColorTable(int n_colors);

/**
 * @brief Create a Rainbow Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 12, 2018
 */
//! swig(newobject)
iftColorTable *iftRainbowColorTable(int n_colors);

/**
 * @brief Create a Categorical Color Table in the color space YCbCr..
 *
 * If the number of colors is <= 21, the Category21 color table is returned (@see iftCategory21ColorTable).
 * Otherwise, the colors from the 22nd are shuffled based on the gradient of some colors.
 * 
 * @author Alexandre Falcao, Samuel Martins
 * @date Jul 20th, 2018
 */
//! swig(newobject)
iftColorTable *iftCategoricalColorTable(int n_colors);

  
/**
 * @brief Create a Reverse Rainbow Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 12, 2018
 */
//! swig(newobject)
iftColorTable *iftReverseRainbowColorTable(int n_colors);

/**
 * @brief Create a Heat Map Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 12, 2018
 */
//! swig(newobject)
iftColorTable *iftHeatMapColorTable(int n_colors);


/**
 * @brief Create a Red Hot Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 2, 2018
 */
//! swig(newobject)
iftColorTable *iftRedHotColorTable(int n_colors);


/**
 * @brief Create a Green Hot Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 2, 2018
 */
//! swig(newobject)
iftColorTable *iftGreenHotColorTable(int n_colors);


/**
 * @brief Create a Blue Hot Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 2, 2018
 */
//! swig(newobject)
iftColorTable *iftBlueHotColorTable(int n_colors);


/**
 * @brief Create a Red-Yellow Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Jun 5, 2018
 */
//! swig(newobject)
iftColorTable *iftRedYellowHotColorTable(int n_colors);


/**
 * @brief Create a Red-to-Blue Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Oct 23, 2019
 */
//! swig(newobject)
iftColorTable *iftRedToBlueColorTable(int n_colors);


/**
 * @brief Create a Viridis Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Feb 4, 2020
 * 
 * https://www.kennethmoreland.com/color-advice/
 */
//! swig(newobject)
iftColorTable *iftViridisColorTable(int n_colors);


/**
 * @brief Create a Plasma Color Table in the color space YCbCr.
 * @author Samuka Martins
 * @date Feb 4, 2020
 * 
 * https://www.kennethmoreland.com/color-advice/
 */
//! swig(newobject)
iftColorTable *iftPlasmaColorTable(int n_colors);


/**
 * @brief Create a Categorical Color Table with up to 21 colors in the color space YCbCr..
 *
 * This color table is inspired at the color tables Set1 and Category20 from
 * https://vega.github.io/vega/docs/schemes/
 * 
 * @author Samuka Martins
 * @date Jun 22, 20180
 */
//! swig(newobject)
iftColorTable *iftCategory21ColorTable(int n_colors);

/**
 * @brief Create a Categorical Color Table with up to 10 colors in the color space YCbCr..
 *
 * This color table is based on the color table Category10 from
 * https://vega.github.io/vega/docs/schemes/
 * 
 * @author Cesar Castelo
 * @date Set 15, 2018
 */
//! swig(newobject)
iftColorTable *iftCategory10ColorTable(int n_colors);


#ifdef __cplusplus
}
#endif

#endif
