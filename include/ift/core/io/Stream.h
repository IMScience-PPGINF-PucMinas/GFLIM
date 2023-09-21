//
// Created by Samuel Martins on 2018-12-16.
//

#ifndef IFT_STREAM_H
#define IFT_STREAM_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/tools/Dialog.h"


/**
 * @brief Allocates a contigous memory block of <n> chunks of size <size>. See also iftFree().
 * @param n Number of memory chunks.
 * @param size Size of the memory chunk.
 * @return Pointer to the first chunk of the allocated stream data.
 * @author Peixinho
 */
static inline void *iftAlloc(size_t n, size_t size) {
    return calloc(n, size);
}


/**
 * @brief Releases a data stream (contigous memory block) previously allocated. See also iftAlloc().
 * @param data Pointer to the first chunk of the stream data.
 * @author Peixinho
 */
static inline void iftFree(void *data) {
    //so you can just call free for NULL pointer, and be happy :)
    if (data != NULL)
      free(data);
}


/**
 * @brief Reallocates a stream data (contigous memory block) to a new one with n chunks of size <size>. See also iftFree().
 * @param data Pointer to the first chunk of the stream data.
 * @param size Size of the memory chunk.
 * @return Pointer to the first chunk of the reallocated stream data.
 * @author Peixinho
 */
static inline void *iftRealloc(void *data, size_t size) {
    return realloc(data, size);
}


/**
 * @brief Copies a source stream data (contigous memory block) <src> with <n> memory chunks of size <chunk_size>
 * to a pre-allocated destination stream data <dst>.
 *
 * @warning It expects that <dst> has enough memory space to copy the source data stream.
 *
 * @param dst Destination stream data.
 * @param src Source destination stream data.
 * @param n Number of memory chunks from src.
 * @param chunk_size Memory chuck size.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
static inline void iftCopyData(void *dst, const void *src, size_t n, size_t chunk_size) {
    size_t nbytes = n * chunk_size;
    memmove(dst, src, nbytes);
}


/**
 * @brief Returns the content of a FILE with pathname <code><b>pathname</code></b> in a string.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname Pathname from the FILE to be read.
 * @return The FILE content as a String.
 *
 * @note All the content of the file is stored in an only string.
 */
char *iftReadFileAsString(const char *pathname, ...);


/**
 * @brief Writes a data stream <data> with <n> memory chunks of size <chunk_size> bytes into a file <file>.
 *
 * @param file File where the stream will be written.
 * @param data Data stream to be written.
 * @param n Number of memory chunks from the data stream.
 * @param chuck_size Size (in bytes) of the memory chunks from the data stream.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteDataInFile(FILE *file, const void *data, size_t n, size_t chuck_size);


/**
 * @brief Writes a string to a a file.
 *
 * @param str String to be written to a File.
 * @param pathname Pathname from the file (it can be built on time like printf)
 *
 * @author Samuka Martins
 * @date Nov 4, 2018
 * @ingroup File
 */
void iftWriteStringToFile(const char *str, const char *pathname, ...);


/**
 * @brief Zip the Content of a Dir.
 * @author Samuel Martins
 * @date Jun 9, 2016
 * @ingroup File
 *
 * @param dir_path     Directory's pathname.
 * @param out_zip_path Pathname from the Output Zip file.
 */
void iftZipDirContent(const char *dir_path, const char *out_zip_path);


/**
 * @brief Unzip a zip file to a directory.
 * @author Samuel Martins
 * @date Jun 9, 2016
 * @ingroup File
 *
 * @param zip_path Zip file to be unzipped.
 * @param out_dir  Output directory where the zip file will be unzipped
 */
void iftUnzipFile(const char *zip_path, const char *out_dir);


/**
 * @brief Finds a string in a File. If the string exists, it returns true. Otherwise, returns false.
 * @author Samuel Martins
 * @date Dec 8th, 2015
 * @ingroup File
 *
 * @param pathname Pathname from the FILE to be read.
 * @param str String to be searched in the FILE.
 * @return True if the string exists in the File. False, otherwise.
 *
 * @exception pathname does not exist or is not a File.
 */
bool iftFindInFile(const char *pathname, const char *str);


/**
 * @brief Check if a file is empty.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param fp File pointer to the file.
 * @return True if it is a directory, false otherwise.
 */
bool iftIsFileContentEmpty(FILE *fp);


/**
 * @brief Get a line (without \n) from a File. If there is no line to get (cursor is in the end of the file),
 * it returns NULL.
 * @param  stream Open file whose line will be read.
 * @return        The desired line.
 *
 * @author Samuka Martins
 * @date Oct 31, 2017
 */
char *iftGetLine(FILE *stream);


/**
 * @brief Sed command: Replace the string old_str with new_str into file file_path.
 * @param  file_path Pathname from the target file.
 * @param  old       Old string to be replaced.
 * @param  new       New string.
 *
 * @attention The input file must is read as ASCII
 *
 * @author Samuka Martins
 * @date Oct 31, 2017
 */
void iftSed(const char *file_path, const char *old_str, const char *new_str);


/**
 * @brief Returns the size (in bytes) from a given datatype.
 * @author Samuka
 * @date Dec 15, 2018
 */
size_t iftCDataTypeSize(iftCDataType cdtype);

#ifdef __cplusplus
}
#endif

#endif //IFT_STREAM_H
