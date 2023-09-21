//
// Created by Samuel Martins on 2018-12-19.
//

#ifndef IFT_DIR_H
#define IFT_DIR_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/io/File.h"


/**
 * @brief An abstract representation of directories. Contains a list of files and subdirectories under that directory.
 * @author Samuel Martins
 * @ingroup File
 */
typedef struct ift_dir {
    char *path;
    /** Number of files from the directory */
    long nfiles;
    /** Vector of iftFiles with the children */
    iftFile **files;
    /** Number of files from the directory */
    long nsubdirs;
    /** Vector of iftDir with the subdirs */
    struct ift_dir **subdirs;
} iftDir;



/**
 * @brief Destroys an iftDir.
 * @author Samuel Martins
 * @ingroup File
 *
 * Destroys the iftDir <b><code>dir</code></b> from memory, BUT IT DOES NOT DELETE IT FROM THE DISK.
 * If its sub-directories also were listed/loaded, the function also destroys them.
 *
 * @param dir The iftDir to be destroyed.
 */
void iftDestroyDir(iftDir **dir);


/**
 * @brief Prints the information of an iftDir.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param dir The iftDir to be printed.
 *
 * @attention It does not print the information of its sub-directories.
 */
void iftPrintDirInfo(const iftDir *dir);


/**
 * @brief Prints all files and sub-directories from an iftFile directory as a tree.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param dir The iftDir to be printed.
 */
void iftPrintDirAsTree(const iftDir *dir);


/**
 * @brief Creates a temporary directory on <b>DISK</b>. If a prefix and/or a suffix are given (!= NULL), they are added into the temp dir.
 * @author Samuel Martins
 * @date Dec 11th, 2015
 * @ingroup File
 *
 * The temporary directory is created in a given directory <dir_name>. If <dir_name> is NULL, the directory /tmp is considered.
 *
 * @param prefix Prefix from the temporary directory. If it is NULL, no prefix is added.
 * @param suffix Suffix from the temporary directory. If it is NULL, no suffix is added.
 * @param dir_name Directory where the temp directory is created. If it is NULL, nothe /tmp is considered.
 *
 * @return The temporary directory pathname.
 * @sa iftMakeTempFile(), iftMakeTempPathname()
 */
char *iftMakeTempDir(const char *prefix, const char *suffix, const char *dir_name);


/**
 * @brief Copies a <b>Directory</b> and all its content (files and subdirectories) from DISK.
 *
 * @author Samuel Martins
 * @date December 01, 2015
 * @ingroup File
 *
 * @attention All exception checks were done on the function iftCopyFromDisk().
 */
void iftCopyDirFromDisk(const char *src, const char *dst);


/**
 * @brief Removes the directory <code><b>dir_pathname</code></b> from the <b>DISK</b>.
 * @author Samuel Martins
 * @date Oct 7th, 2015
 * @ingroup File
 *
 * @param dir_pathname The directory to be removed from the DISK.
 *
 * @warning If the <code><b>dir_pathname</code></b> is a file or it does not exist, nothing will be done.
 */
void iftRemoveDir(const char *dir_pathname);


/**
 * @brief Loads a Directory <dir_pathname>.
 * @author Samuel Martins
 * @date Aug 14, 2015
 * @ingroup File
 *
 * Loads a Directory from a <dir_pathname> on memory.
 * THE DIRECTORY MUST EXIST.
 * The parameter @p hier_levels indicates until which hierarchical level will be considered in the
 * loading.
 *
 * If <code>hier_levels=2</code>, all files and subdirs until the 2º level from the root directory are loaded, but
 * the content of the subdirs from the 2º level are not loaded.
 * If  <code>hier_levels=0</code>, all files and subdirs from all levels are loaded.
 * If <code>hier_levels is</code> greater than the real number of levels from the root directory, all files and subdirs
 * will be loaded.
 *
 * A separator char '/' (linux) or '\' (windows) is put at the end of the pathname.
 * If the file exists and is not a directory (it is a file), an error will be raised.
 *
 * @param dir_pathname The directory to be read.
 * @param hier_levels Number of levels that will be loaded. The vale 0 loads all files from the whole hierarchy.
 * @return An iftDir with the directory.
 */
iftDir *iftLoadDir(const char *dir_pathname, long hier_levels);


/**
 * @brief Loads all files (not subdirs) in the directory <dir_pathname> from a given <extension>.
 * @author Samuel Martins
 * @ingroup File
 *
 * Loads all files (not subdirs) in the directory <b><code>dir_pathname</code></b> from a given <b><code>extension</code></b>.
 * The filenames are sorted in ascending order.
 * It only gets the files from the 1st file hierarchical level.
 * If <b><code>extension</code></b> = "", it gets all files (not subdirs).
 *
 * @param dir_pathname The directory to be read.
 * @param suffix Suffix from the files.
 * @return An iftDir with all files with <extension> inside <dir_pathname>.
 */
iftDir *iftLoadFilesFromDirBySuffix(const char *dir_pathname, const char *suffix);


/**
 * @brief Loads all files (not subdirs) in the directory <dir_pathname> from a given <regex>.
 * @author Samuel Martins
 * @author Jul 7, 2016
 * @ingroup File
 *
 * Loads all files (not subdirs) in the directory <b><code>dir_pathname</code></b> from a given <b><code>regex</code></b>. \n
 * The filenames are sorted in ascending order. \n
 * It only gets the files from the 1st file hierarchical level. \n
 * Example of Regex: \n
 * regex = ".*Pos.*\\.scn", equivalent to ls *Pos*.scn in bash-linux
 *
 * @param dir_pathname The directory to be read.
 * @param regex        Target Regex.
 * @return             An iftDir with all files with <extension> inside <dir_pathname> filtered by <regex>.
 */
iftDir *iftLoadFilesFromDirByRegex(const char *dir_pathname, const char *regex);


/**
 * @brief Sorts the file list in the directory by pathname
 * @author Thiago Vallin Spina
 * @ingroup File
 *
 * @param files The file list.
 */
void iftSortDir(iftDir *files);


/**
 * @brief Get the file labels, according to the filed path under that directory.
 * @author Peixinho
 *
 * The image label is stored as part of the image file name.
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param dir Directory containing the files.
 * @return The image labels.
 * @sa iftValidDataSetFilenameConvention
 */
int* iftImageLabels(iftDir *dir);


/**
 * @brief Checks if the <b><code>pathname</code></b> is a directory on the disk.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname to be checked.
 * @return True if it is a directory, false otherwise.
 */
bool iftDirExists(const char *pathname, ...);


/**
 * @brief Gets the parent dir from a file or directory.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname of the file/directory.
 * @return The parent directory from <pathname>.
 *
 * @warning It does not check if the parent dir exists.
 * @warning Return the parent dir WITHOUT THE SLASH AT THE END: e.g: /home/samuel --> parent_dir = /home
 */
char *iftParentDir(const char *pathname);


/**
 * @brief Creates all directories in a hierarchy on the disk if they do not exist.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname of the directory to be created on disk.
 *
 * @warning If the directory already exists, do nothing.
 * @note http://www.thegeekstuff.com/2012/06/c-directory/
 */
void iftMakeDir(const char *dir_path);



#ifdef __cplusplus
}
#endif

#endif //IFT_DIR_H
