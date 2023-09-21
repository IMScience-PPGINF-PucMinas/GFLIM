//
// Created by Samuel Martins on 2018-12-19.
//

#ifndef IFT_FILE_H
#define IFT_FILE_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#if defined(__WIN32) || defined(__WIN64)
#include <Windows.h>
#endif


// DIRECTORY SEPARATOR STRING
/**
 * @brief Directory Separator String.
 * @author Samuel Martins
 * @ingroup File
 */
#if defined(__WIN32) || defined(__WIN64)
#define IFT_SEP_C "\\"
#else
#define IFT_SEP_C "/"
#endif

/**
 * @brief An abstract representation of file.
 *
 * \ref iftFile contains the file info for creating, using and managing disk files.
 *
 * @author Samuel Martins
 * @ingroup File
 */
typedef struct ift_file {
    /** File path */
    char *path;
    /**Label for machine learning dataset files*/
    int label;
    /** Sample id following the new convention for IFT datasets (Thiago Vallin Spina, Mar 09, 2016) */
    int sample;
    /** Suffix following the new convention for IFT datasets (Thiago Vallin Spina, Mar 09, 2016) */
    char *suffix;
    /**Status for machine learning dataset files
     * iftFile can be marked as IFT_TRAIN, TEST, IFT_OUTLIER as the iftSample*/
    char status;
} iftFile;



/**
 * @brief Creates an iftFile.
 * @author Samuel Martins
 * @date Nov 10, 2015
 * @ingroup File
 *
 * Creates an iftFile from a <code><b>pathname</b></code> on memory, BUT IT DOES NOT CREATE IT ON DISK.
 * If the <code><b>pathname</b><code> exists, the field <exists> is set to 1.
 * Otherwise, it is set to 0.
 * If the file exists and is a directory, an error will be raised.
 *
 *
 * @param pathname The file to be read.
 * @return An iftFile with the <code><b>pathname</b><code>.
 * @sa iftLoadDir()
 */
iftFile *iftCreateFile(const char *pathname, ...);


/**
 * @brief Destroys an iftFile from memory, <b>BUT IT DOES NOT DELETE IT FROM THE DISK</b>.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param file The iftFile to be destroyed.
 */
void iftDestroyFile(iftFile **file);


/**
 * @brief Prints the information of an iftFile.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param f The iftFile to be printed.
 *
 * @warning It does not print the information of sub-directories from an iftFile dir.
 */
void iftPrintFileInfo(const iftFile *f);


/**
 * @brief Creates a valid temporary file on <b>DISK</b>.
 * @author Samuel Martins
 * @date Dec 11, 2015
 * @ingroup File
 *
 * The temporary file is created in a given directory <dir_name>. If <dir_name> is NULL, the directory /tmp is considered.
 *
 * @param prefix Prefix from the tmp file. If it is NULL, no prefix is added.
 * @param suffix Suffix from the tmp file. If it is NULL, no prefix is added.
 * @param dir_name Directory where the temp file is created. If it is NULL, nothe /tmp is considered.
 *
 * @return The temporary file pathname.
 * @sa iftTempPathname(), iftMakeTempDir()
 */
char *iftMakeTempFile(const char *prefix, const char *suffix, const char *dir_name);


/**
 * @brief Creates a valid temporary pathname. If a prefix and/or a suffix are given (!= NULL), they are added into the temp pathname.
 *
 * The temporary pathname is created in a given directory <dir_name>. If <dir_name> is NULL, the directory /tmp is considered.
 *
 * @param prefix Prefix from the temporary pathname. If it is NULL, no prefix is added.
 * @param suffix Suffix from the temporary pathname. If it is NULL, no suffix is added.
 * @param dir_name Directory where the temp pathname is created. If it is NULL, nothe /tmp is considered.
 * @return The temporary pathname.
 *
 * @attention No file or directory is created on disk.
 * @sa iftMakeTempFile(), iftMakeTempDir()

 * @author Samuel Martins
 * @date Dec 11th, 2015
 * @ingroup File
 */
char *iftMakeTempPathname(const char *prefix, const char *suffix, const char *dir_name);


/**
 * @brief Copies an iftFile.
 * @author Peixinho
 * @ingroup File
 *
 * @param file The iftFile to be copied.
 * @return A copy of the input iftFile.
 */
iftFile *iftCopyFile(const iftFile* file);


/**
 * @brief Copies any <b>File</b> (text, binary, pdf, doc, png, ...) or <b>Directory</b> from DISK.
 * @author Samuel Martins
 * @date Dec 1, 2015
 * @ingroup File
 *
 * @param src Filename from the source file/directory to be copied.
 * @param dst Destination Filename from the copy. It can be a directory filename where the src will be copied.
 *
 * @attention If the <code><b>dst</b></code> filename is an existing directory, the <b>src</b> (file or directory) will be copied into the <b>dst</b> directory with the same src basename.
 *
 * @exception src filename does not exist.
 * @exception src filename is a <i>Directory</i> and the <code><b>dst</b></code> is a <i>File</i>.
 */
void iftCopyFromDisk(const char *src, const char *dst);


/**
 * @brief Removes the file @p <b>pathname</b> from the <b>DISK</b>.
 * @author Samuel Martins
 * @date Oct 7th, 2015
 * @ingroup File
 *
 * @param pathname The file to be removed from the DISK.
 *
 * @warning If the <code><b>pathname</code></b> is a directory or it does not exist, nothing will be done.
 */
void iftRemoveFile(const char *pathname);


/**
 * @brief Rename the file @p <b>oldfilename</b> to <b>newfilename</b>.
 * If newfilename points to a different path, the file will be moved
 * @author Cesar Castelo
 * @date Jan 17, 2019
 * @ingroup File
 *
 * @param oldfilename The file to be renamed
 * @param newfilename The new filename for the file
 *
 */
void iftRenameOrMoveFile(const char *oldfilename, const char *newfilename);


/**
 * @brief Checks if exists a file/directory with the <pathname>.
 *
 * This function checks if exists a file/directory with the <pathname>.
 *
 * @param pathname The pathname to be checked.
 * @return 1 if the file/directory exists, 0 otherwise.
 *
 * @author Samuel Martins
 * @date Dec 1, 2015
 */
static inline bool iftPathnameExists(const char *pathname) {
    struct stat buffer;
    return (stat (pathname, &buffer) == 0);
}


/**
 * @brief Checks if the <pathname> is a file (not dir) and exists on the disk.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname to be checked.
 * @return True if it is a file, false otherwise.
 */
bool iftFileExists(const char *pathname);


/**
 * @brief Gets the filename from a pathname. If <code><b>suffix</code></b> is != NULL and != "", the suffix will be removed from the filename.
 * @author Samuel Martins
 * @date Dec 11, 2015
 * @ingroup File
 *
 * E.g: \n
 * ("/home/batman/Pictures/avatar.png", NULL) = "avatar.png" \n
 * ("/home/batman/Pictures/avatar.png", .png) = "avatar"
 * ("/home/batman/Pictures", .png) = "/home/batman/Pictures"
 *
 * @param pathname The Input Pathname.
 * @param suffix The suffix to be removed from the Oathaname. If it is NULL or "", no suffix is removed.
 * @return The Filename from <b>pathname</b>.
 *
 * @note If the given <code><b>suffix</code></b> is not a suffix from the pathname, nothing will be removed.
 * @attention If <code><b>suffix</code></b> is != NULL and != "", the suffix will be removed from the filename.
 */
//! swig(newobject)
char *iftFilename(const char *pathname, const char *suffix);

/**
 * @brief Gets the dirname from a pathname.
 * @author Cesar Castelo
 * @date Jun 3, 2019
 * @ingroup File
 *
 * E.g: \n
 * ("/home/batman/Pictures/avatar.png") = "/home/batman/Pictures/"
 * ("/home/batman/Pictures/avatar") = "/home/batman/Pictures/"
 * ("/home/batman/Pictures/avatar/") = "/home/batman/Pictures/avatar/"
 *
 * @param pathname The Input Pathname.
 * @return The dirname from <b>pathname</b>.
 *
 * @note If the last character from pathname is "/", then the same string is returned, since it is already a directory
 */
//! swig(newobject)
char *iftDirname(const char *pathname);

/**
 * @brief Gets the Extension of file <code><b>pathname</code></b>. It there is no extension, a blank string "" will be returned.
 * @author Samuel Martins
 * @date Dec 11, 2015
 * @ingroup String
 *
 * @param str Pathname from the file.
 * @return The file extension or a blank string "", if there is no extension.
 *
 * @note The extension is returned getting the strings after the last '.' (including it)
 * @warning Only files with extensions "*.tar.gz", "*.tar.bz", and "*.tar.bz2" has the correct extension returned.
 */
//! swig()
const char *iftFileExt(const char *pathname);


/**
 * @brief Checks if the File Extensions are the same.
 * @author Azael Sousa
 * @date Nov 29th, 2017
 * @ingroup File
 *
 * @param file1 is the pathname of the first file
 * @param file2 is the pathname of the second file
 * @return True if the extensions are the same, false otherwise.
 *
 * @exception not case-sensitive.
 */
bool iftAreFileExtensionsEqual(const char *file1,const char *file2);


/**
 * @brief Joins two pathnames. It automatically treats the '/' (separation char) of the end of <code><b>pathname1</code></b>
 * and the beginning of the <code><b>pathname2</code></b>.
 * @author Samuel Martins
 * @date Dec 13rd, 2015
 * @ingroup File
 *
 * @param The number of paths to be joined.
 * @return The joined pathname.
 */
char *iftJoinPathnames(long n, ...);


/**
 * @brief Gets the absolute path from a pathname.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname.
 * @return The absolute pathname from <code><b>pathname</code></b>.
 */
char *iftAbsPathname(const char *pathname);


/**
 * @brief Gets the Basename from a Pathname.
 * @author Samuel Martins
 * @date Mar 21, 2016
 * @ingroup File
 *
 * E.g: \n
 * ("/home/batman/Pictures/avatar.png") = "/home/batman/Pictures/avatar" \n
 * ("my_image.zip", .png) = "my_image"
 *
 * @param  pathname Input Pathname
 * @return          The Basename from the Input Pathname
 */
//! swig(newobject)
char *iftBasename(const char *pathname);


/**
 * @brief Checks if the Image Pathname has a valid suported Extension (scn, pgm, ppm, png).
 * @author Samuel Martins
 * @date Mar 3, 2016
 * @ingroup File
 *
 * @param  img_pathname Image Pathname to be checked.
 * @return              True if it's valid, false otherwise;
 *
 * @warning It DOESN'T check if it exists a file with this pathname.
 */
bool iftIsImagePathnameValid(const char *img_pathname);


/**
 * @brief Checks if the Pathname is from a valid Image Pathname (*.scn, *.pgm, *.ppm, *.png).
 * @author Samuel Martins
 * @date Dec 9th, 2015
 * @ingroup File
 *
 * @param img_pathname The Image Pathname to be checked.
 * @return True if the image pathname is valid, false otherwise.
 *
 * @exception img_pathname does not exist or is a directory.
 */
bool iftIsImageFile(const char *img_pathname);


/**
 * @brief Adds Escape Characters to '(', ')', and '/' (linux, mac) or '\' (windows) to a Pathname
 *        (it should work on Linux, Mac and Windows).
 * @author Samuel Martins
 * @date Dec 2nd, 2015
 * @ingroup File
 *
 * E.g (Linux): pathname = "/usr/local/(abc)"
 *              escaped  = "\/usr\/local\/\(abc\)"
 *
 * @param pathname Pathname to be used.
 * @return The escaped pathname.
 *
 * @warning Function limitation: The function just works to pathnames with at most 1000 characteres to escape.
 */
char *iftAddEscapeChars(const char *pathname);


/**
 * @brief Replace the "~/" from an input path with the absolute home path
 * @param path Input path
 * @author Azael Sousa, Samuel Botter
 * @date May 18, 2018
 */
char *iftExpandUser(const char *path);


/**
 * @brief Get the file index, according to the path.
 * @author Peixinho
 *
 * The image index is stored as part of the image file name.
 * The IFT default nomenclature is <label>_<sample>_<suffix>.<ext>
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param path Image path.
 * @return The image index.
 * @sa iftValidDataSetFilenameConvention
 */
int iftImageSampleId(char *path);


/**
 * @brief Get the file label, according to the path.
 * @author Peixinho
 *
 * The image label is stored as part of the image file name.
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param path Image path.
 * @return The image label.
 * @sa iftValidDataSetFilenameConvention
 */
int iftImageSampleLabel(char *path);


/**
 * @brief Get the file label and id, according to the path.
 * @author Deangeli
 *
 * The image label is stored as part of the image file name.
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param path Image path.
 * @return A int pointer which the first position (index 0) stores the image label and the second one (index 1) stores
 * the sample id.
 * @sa iftValidDataSetFilenameConvention
 */
int* iftGetImageSampleInfo(char *path);


/**
 * @brief Get the file suffix, if the path follows IFT's convention.
 * @author Thiago Vallin Spina
 *
 * The image suffix is stored as part of the image file name.
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param path Image path.
 * @return The image suffix.
 * @sa iftValidDataSetFilenameConvention
 */
char* iftImageSuffix(char* path);


/**
 * @brief Checks if the path's filename conforms to IFT's default nomenclature.
 *
 * @author Thiago Vallin Spina
 * @date Mar 09, 2016
 *
 * The IFT default nomenclature is <label>_<sample>_<suffix>.<ext>
 *
 * Ex1.: 000002_000015_LeftLung.pgm
 *
 * The suffix may also be empty to conform with the old convention, or be separated by multiple '_' characters
 *
 * Ex2.: 000002_000015.pgm
 * Ex3.: 000002_000015_Right_Lung.pgm
 *
 * @param path The path name.
 *
 * @note The old convention 000002_0000015.pgm is also respected.
 */
bool iftValidDataSetFilenameConvention(const char *path);


/**
 * @brief Updstes the dataset information of a file.
 * @param file Input file.
 * @date Dec 1, 2015
 */
void iftUpdateDataSetFileInfo(iftFile* file);



#ifdef __cplusplus
}
#endif

#endif //IFT_FILE_H
