/* libtiff/tiffconf.h.  Generated from tiffconf.h.in by configure.  */
/*
  Configuration defines for installed libtiff.
  This file maintained for backward compatibility. Do not use definitions
  from this file in your programs.
*/
#include <stdint.h>

#ifndef _TIFFCONF_
#define _TIFFCONF_

/* Signed 16-bit type */
#ifndef TIFF_INT16_T
#define TIFF_INT16_T int16_t
#endif

/* Signed 32-bit type */
#ifndef TIFF_INT32_T
#define TIFF_INT32_T int32_t
#endif

/* Signed 64-bit type */
#ifndef TIFF_INT64_T
#define TIFF_INT64_T int64_t
#endif

/* Signed 8-bit type */
#ifndef TIFF_INT8_T
#define TIFF_INT8_T int8_t
#endif

/* Unsigned 16-bit type */
#ifndef TIFF_UINT16_T
#define TIFF_UINT16_T uint16_t
#endif

/* Unsigned 32-bit type */
#ifndef TIFF_UINT32_T
#define TIFF_UINT32_T uint32_t
#endif

/* Unsigned 64-bit type */
#ifndef TIFF_UINT64_T
#define TIFF_UINT64_T uint64_t
#endif

/* Unsigned 8-bit type */
#ifndef TIFF_UINT8_T
#define TIFF_UINT8_T uint8_t
#endif

/* Signed size type */
#define TIFF_SSIZE_T signed long

/* Pointer difference type */
#define TIFF_PTRDIFF_T ptrdiff_t

/* Define to 1 if the system has the type `int16'. */
/* #undef HAVE_INT16 */

/* Define to 1 if the system has the type `int32'. */
/* #undef HAVE_INT32 */

/* Define to 1 if the system has the type `int8'. */
/* #undef HAVE_INT8 */

/* Compatibility stuff. */

/* Define as 0 or 1 according to the floating point format suported by the
   machine */
#define HAVE_IEEEFP 1

/* Set the native cpu bit order (FILLORDER_LSB2MSB or FILLORDER_MSB2LSB) */
#define HOST_FILLORDER FILLORDER_LSB2MSB

/* Native cpu byte order: 1 if big-endian (Motorola) or 0 if little-endian
   (Intel) */
#define HOST_BIGENDIAN 0

/* Support CCITT Group 3 & 4 algorithms */
#define CCITT_SUPPORT 1

/* Support JPEG compression (requires IJG JPEG library) */
/* #undef JPEG_SUPPORT */

/* Support JBIG compression (requires JBIG-KIT library) */
/* #undef JBIG_SUPPORT */

/* Support LogLuv high dynamic range encoding */
#define LOGLUV_SUPPORT 1

/* Support LZW algorithm */
#define LZW_SUPPORT 1

/* Support NeXT 2-bit RLE algorithm */
#define NEXT_SUPPORT 1

/* Support Old JPEG compresson (read contrib/ojpeg/README first! Compilation
   fails with unpatched IJG JPEG library) */
/* #undef OJPEG_SUPPORT */

/* Support Macintosh PackBits algorithm */
#define PACKBITS_SUPPORT 1

/* Support Pixar log-format algorithm (requires Zlib) */
/* #undef PIXARLOG_SUPPORT */

/* Support ThunderScan 4-bit RLE algorithm */
#define THUNDER_SUPPORT 1

/* Support Deflate compression */
/* #undef ZIP_SUPPORT */

/* Support strip chopping (whether or not to convert single-strip uncompressed
   images to mutiple strips of ~8Kb to reduce memory usage) */
#define STRIPCHOP_DEFAULT TIFF_STRIPCHOP

/* Enable SubIFD tag (330) support */
#define SUBIFD_SUPPORT 1

/* Treat extra sample as alpha (default enabled). The RGBA interface will
   treat a fourth sample with no EXTRASAMPLE_ value as being ASSOCALPHA. Many
   packages produce RGBA files but don't mark the alpha properly. */
#define DEFAULT_EXTRASAMPLE_AS_ALPHA 1

/* Pick up YCbCr subsampling info from the JPEG data stream to support files
   lacking the tag (default enabled). */
#define CHECK_JPEG_YCBCR_SUBSAMPLING 1

/* Support MS MDI magic number files as TIFF */
#define MDI_SUPPORT 1

/*
 * Feature support definitions.
 * XXX: These macros are obsoleted. Don't use them in your apps!
 * Macros stays here for backward compatibility and should be always defined.
 */
#define COLORIMETRY_SUPPORT
#define YCBCR_SUPPORT
#define CMYK_SUPPORT
#define ICC_SUPPORT
#define PHOTOSHOP_SUPPORT
#define IPTC_SUPPORT

#endif /* _TIFFCONF_ */