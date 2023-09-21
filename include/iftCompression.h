

/**
 * @file iftCompression.h
 * @brief Low level i/o interface to compressed and noncompressed files.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 22, 2016
 * @ingroup Compression
 *
 * @note Programs:
 * * @ref iftTestCompression.c = It loads a 3D image file and compresses it into a second file with .zscn or .scn.gz file extension,
 * * and does the reverse process.
 */

#ifndef IFT_COMPRESSION_H
#define IFT_COMPRESSION_H

/*
 * iftCompression.h  (zipped or non-zipped library)
 * *****            This code is released to the public domain.            *****

 * *****  Author: Mark Jenkinson, FMRIB Centre, University of Oxford       *****
 * *****  Date:   September 2004                                           *****

 * *****  Neither the FMRIB Centre, the University of Oxford, nor any of   *****
 * *****  its employees imply any warranty of usefulness of this software  *****
 * *****  for any purpose, and do not assume any liability for damages,    *****
 * *****  incidental or otherwise, caused by any use of this document.     *****
 */


/**
 *
 * This library provides an interface to both compressed (gzip/zlib) and
 * uncompressed (normal) file IO.  The functions are written to have the
 * same interface as the standard file IO functions.
 *
 * To use this library instead of normal file IO, the following changes
 * are required:
 *  - replace all instances of FILE* with iftGZipFile
 *  - change the name of all function calls, replacing the initial character
 *    f with the znz  (e.g. fseek becomes iftGZipSeek)
 *    one exception is rewind() -> iftGZipRewind()
 *  - add a third parameter to all calls to znzopen (previously fopen)
 *    that specifies whether to use compression (1) or not (0)
 *  - use znz_isnull rather than any (pointer == NULL) comparisons in the code
 *    for znzfile types (normally done after a return from iftGZipOpen)
 *
 * NB: seeks for writable files with compression are quite restricted
 *
 */


/*=================*/
#ifdef  __cplusplus
extern "C" {
#endif
/*=================*/

#include "iftCommon.h"

/* include optional check for HAVE_FDOPEN here, from deleted config.h:

   uncomment the following line if fdopen() exists for your compiler and
   compiler options
*/
/* #define HAVE_FDOPEN */

#ifndef HAVE_ZLIB
#define HAVE_ZLIB
#endif

#ifdef HAVE_ZLIB
#if defined(ITKZLIB)
#include "itk_zlib.h"
#else
#include "zlib.h"
#endif
#endif


/**
 * @brief Structure that may be used to store a file with or without compression in seamless mode.
 */
struct iftGZipPtr {
    bool withz;
    FILE* nzfptr;
#ifdef HAVE_ZLIB
    gzFile zfptr;
#endif
};

/* the type for all file pointers */
typedef struct iftGZipPtr *iftGZipFile;


/* int znz_isnull(iftGZipFile f); */
/* int znzclose(iftGZipFile f); */
#define znz_isnull(f) ((f) == NULL)
#define znzclose(f)   Xznzclose(&(f))

/**
 * @brief Opens a file for storing compressed or uncompressed data. Replaces function fopen.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 * @param path The path to the file
 * @param mode The reading mode ("r", "rb", "w", "wb", and all other available standard modes)
 * @param use_compression If true, then the file will be compressed. Otherwise, a standard file is written/read.
 * @ingroup Compression
 *
 * @return The file pointer
 */
iftGZipFile iftGZipOpen(const char *path, const char *mode, bool use_compression);

/**
 * @brief Opens a file for storing compressed or uncompressed data. Replaces function open.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 * @param fd The file descriptor number.
 * @param mode The reading mode ("r", "rb", "w", "wb", and all other available standard modes)
 * @param use_compression If true, then the file will be compressed. Otherwise, a standard file is written/read.
 *
 * @return The file pointer
 */
iftGZipFile iftGZipDOpen(int fd, const char *mode, bool use_compression);

/**
 * @brief Closes a zipped or not zipped file. Replaces function fclose.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 * @note The <file> variable is freed and set to NULL afterwards
 */
int iftGZipClose(iftGZipFile *file);


/**
 * @brief Replaces function fread. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
size_t iftGZipRead(void *buf, size_t size, size_t nmemb, iftGZipFile file);

/**
 * @brief Replaces function fwrite. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
size_t iftGZipWrite(const void *buf, size_t size, size_t nmemb, iftGZipFile file);

/**
 * @brief Replaces function fseek. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
long iftGZipSeek(iftGZipFile file, long offset, int whence);

/**
 * @brief Replaces function frewind. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
int iftGZipRewind(iftGZipFile stream);

/**
 * @brief Replaces function ftell. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
long iftGZipTell(iftGZipFile file);

/**
 * @brief Replaces function fputs. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
int iftGZipPuts(const char *str, iftGZipFile file);

/**
 * @brief Replaces function fgets. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
char *iftGZipGets(char *str, int size, iftGZipFile file);

/**
 * @brief Replaces function fputc. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
int iftGZipPutc(int c, iftGZipFile file);

/**
 * @brief Replaces function fgetc. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
int iftGZipGetc(iftGZipFile file);

#if !defined(WIN32)

/**
 * @brief Replaces function fprintf. Same parameters except for <file>.
 *
 * @author Mark Jenkinson
 * @note Adapted by Thiago Vallin Spina
 * @date Mar 31, 2016
 * @ingroup Compression
 *
 */
int iftGZipPrintf(iftGZipFile stream, const char *format, ...);
#endif

#endif //IFT_IFTCOMPRESS_H

#ifdef __cplusplus
}

#endif


