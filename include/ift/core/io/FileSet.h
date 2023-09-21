//
// Created by Samuel Martins on 2018-12-19.
//

#ifndef IFT_FILESET_H
#define IFT_FILESET_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/Dict.h"
#include "ift/core/io/CSV.h"
#include "ift/core/io/File.h"



/**
 * @brief A struct that contains an array to iftFile*.
 * @author Samuel Martins
 * @date Oct 6th, 2015
 * @ingroup File
 */
//! swig(destroyer = iftDestroyFileSet, extend = iftFileSetExt.i)
typedef struct ift_fileset {
    /** Number of Files */
    long n;
    /** Array of Files */
    iftFile **files;
} iftFileSet;



/**
 * @brief Allocates an iftFileSet structure but does NOT allocate the elements of the file array.
 * @author Thiago Vallin Spina
 *
 * @param n The number of files.
 * @return The allocated structure
 */
iftFileSet *iftCreateFileSet(long nfiles);


/**
 * @brief Creates a file set and initializes it with given number of strings.
 *
 * @param nfiles Number of strings to be copied into the file set.
 * @param ... Sequence of <nfiles> strings that will be copied (in order) into the file set.
 * @return Initialized file set.
 *
 * @author Samuel Martins
 * @date Jan 24, 2019
 */
iftFileSet *iftInitFileSet(long nfiles, ...);


/**
 * @brief Destroys a File Array (iftFileSet) from memory, <b>BUT IT DOES NOT DELETE IT FROM THE DISK</b>.
 * @author Samuel Martins
 * @date Oct 6th, 2015
 * @ingroup File
 *
 * @param file The iftFile to be destroyed.
 */
void iftDestroyFileSet(iftFileSet **farr);


/**
 * @brief Copies a fileset
 * @author Peixinho
 * @param Original fileset
 * @return Fileset copy
 * @date Jun, 2016
 */
iftFileSet *iftCopyFileSet(const iftFileSet* fileset);


/**
 * @brief Loads an Array of iftFiles from a Directory (in the entire directory hierarchy).
 *
 * The parameter @p hier_levels indicates until which hierarchical level will be considered in the
 * loading.
 *
 * If <code>hier_levels=2</code>, all files and subdirs until the 2º level from the root directory are loaded, but
 * the content of the subdirs from the 2º level are not loaded.
 * If  <code>hier_levels=0</code>, all files and subdirs from all levels are loaded.
 * If <code>hier_levels is</code> greater than the real number of levels from the root directory, all files and subdirs
 * will be loaded.
 *
 * @author Samuel Martins
 * @date Oct 6th, 2015
 * @ingroup File
 *
 * @param dir_pathname The pathname from the directory to be read.
 * @param hier_levels Number of folder levels that will be loaded.
 * @return The array of iftFiles.
 *
 * @exception dir_pathname does not exist or it is a File.
 */
iftFileSet *iftLoadFileSetFromDir(const char *dir_pathname, long hier_level);


/**
 * @brief Loads an Array of iftFiles with extension <b><code>ext</code></b> from a Directory.
 * @author Samuel Martins
 * @date Oct 7th, 2015
 * @ingroup File
 *
 * Loads an Array of iftFiles with extension <ext> from a Directory.
 * The function gets the files in the entire directory hierarchy.
 *
 * @param dir_pathname The pathname from the directory to be read.
 * @return The array of iftFiles.
 */
//! swig(newobject)
iftFileSet *iftLoadFileSetFromDirBySuffix(const char *dir_pathname, const char *ext, int depth);


/**
 * @brief Loads an Array of iftFiles with Regex <b><code>regex</code></b> from a Directory.
 * @author Samuel Martins
 * @date Jul 7, 2015
 * @ingroup File
 *
 * @note The function gets the files in the entire directory hierarchy.
 *
 * @param dir_pathname The pathname from the directory to be read.
 * @param regex        Target Regex.
 * @param sort_pathnames Sort the pathnames.
 * @return The array of iftFiles filtered by the regex.
 */
iftFileSet *iftLoadFileSetFromDirByRegex(const char *dir_pathname, const char *regex, bool sort_pathnames);


/**
* @brief Reads an Array of iftFiles from a CSV.
* @author Samuel Martins
* @date Sep 25th, 2015
* @ingroup File
*
* Reads an iftCSV from file.
* The CSV must have exactly the following separating pattern:
*
* Brazil,100,20.50
* France,256,35.00
*
* @param csv_pathname The pathname from the CSV to be read.
* @param sort_pathnames It indicates if the pathnames should be sorted.
* @return The array of iftFiles.
*/
iftFileSet *iftLoadFileSetFromCSV(const char *csv_pathname, bool sort_pathnames);


/**
 * @brief Reads an Array of iftFiles from a Directory or a CSV file.
 *
 * The parameter @p hier_levels indicates until which hierarchical level will be considered in the
 * loading.
 *
 * If <code>hier_levels=2</code>, all files and subdirs until the 2º level from the root directory are loaded, but
 * the content of the subdirs from the 2º level are not loaded.
 * If  <code>hier_levels=0</code>, all files and subdirs from all levels are loaded.
 * If <code>hier_levels is</code> greater than the real number of levels from the root directory, all files and subdirs
 * will be loaded.
 *
 * @author Samuel Martins
 * @date Sep 25th, 2015
 * @ingroup File
 *
 * @param file_entry The pathname from the directory of the CSV file to be read.
 * @param hier_levels Number of levels that will be loaded.
 * @param sort_pathnames It indicates if the pathnames should be sorted.
 *
 * @return The array of iftFiles.
 *
 * @exception The file entry <code><b>file_entry</b></code> is neither a directory nor a CSV file.
 */
//! swig(newobject)
iftFileSet *iftLoadFileSetFromDirOrCSV(const char *file_entry, long hier_levels, bool sort_pathnames);


/**
 * @brief      Writes a File Set as a CSV, one pathname per row.
 *
 * @param[in]  fset  File set to be written.
 * @param[in]  path  Output CSV file.
 *
 * @author Samuka
 * @date   Jan, 27, 2017
 */
void iftWriteFileSetAsCSV(const iftFileSet *fset, const char *path, ...);


/**
 * @brief Removes all files in the file set @p <b>fset</b> from the <b>DISK</b>.
 * @author Samuel Martins
 * @date Dec 10, 2018
 * @ingroup File
 *
 * @param fset File set whose files will be removed from the DISK.
 */
void iftRemoveFileSet(const iftFileSet *fset);


/**
 * @brief Merge two file sets (fset1 comes first than fset2)
 * @author Samuel Martins
 * @date Feb 20, 2016
 */
iftFileSet *iftMergeFileSet(const iftFileSet *fset1, const iftFileSet *fset2);


/**
 * @brief      Filters the files from the file set <dst> by matching the filenames of each pathname from <src> with <dst>.
 *
 * For each pathname from <src>, this functions tries to match its filename with the filename of one pathname in <dst>. \n
 * If matched, the pathname from <dst> is copied. Otherwise an exception will be raised.
 *
 * @param[in]  src   Source file set whose filenames will be filter the destination file set.
 * @param[in]  dst   Destination file set which will be filtered
 * @return           Filtered file set.
 *
 * @author Samuka
 * @date Jan 27, 2017
 * @ingroup File
 */
iftFileSet *iftFilterFileSetByFilename(const iftFileSet *src, const iftFileSet *dst);


/**
 * @brief      Filters the files by indexes.
 *
 * @param[in]  fset  File set to be filtered.
 * @param      idxs  Indexes from the files to be chosen.
 *
 * @return     Filtered file set.
 *
 * @author Samuka
 * @date Jan 30, 2017
 * @ingroup File
 */
iftFileSet *iftFilterFileSetByIndexes(const iftFileSet *fset, const iftIntArray *idxs);


/**
 * @brief Finds the image file by label and returns its index in the list.
 * @author Thiago Vallin Spina
 *
 * @param files File list.
 * @param label Image label.
 *
 * @return The index of the file in the file list or NIL if it is not present.
 */
int iftFindFileByLabel(const iftFileSet *files, int label);


/**
 * @brief Finds the image file by label and sample id and returns its index in the list.
 *
 * @author Thiago Vallin Spina
 * @date Mar 09, 2016
 *
 * @param files File list.
 * @param label Image label.
 *
 * @return The index of the file in the file list or NIL if it is not present.
 */
int iftFindFileByLabelAndSample(const iftFileSet *files, int label, int sample);


/**
 * @brief Finds the file by the path's basename'sand returns its index in the list.
 * @author Thiago Vallin Spina
 *
 * @param files File list.
 * @param prefix The path's basename.
 *
 * @return The index of the file in the file list or NIL if it is not present.
 */
int iftFindFileByBasename(const iftFileSet *files, const char *bname, const char *suffix);


/**
 * Count number of unique label on a file set.
 * @author Peixinho
 * @date Jun, 2017
 * @param fileset
 * @return Number of unique labels
 */
int iftFileSetLabelsNumber(const iftFileSet* fileset);


/**
 * Counts the number of samples per class in a fileset
 * @author Cesar Castelo
 * @date Ago 08, 2018
 * @param fileset
 * @return Array with the number of samples per class
 */
int *iftCountSamplesPerClassFileSet(iftFileSet *fileset);


/**
 * @brief Get the file labels, according to the IFT file format.
 * @author Peixinho
 * @ingroup File
 *
 * The image label is stored as part of the image file name.
 * The IFT default nomenclature is <label>_<id>.<ext>
 *
 * Ex.: 000002_000015.pgm
 *
 * @param dir FileSet containing the files.
 * @return The file labels.
 */
int* iftFileSetLabels(const iftFileSet* fileSet);


/**
 * @brief Converts a FileSet to a CSV file structure.
 *
 * @author Thiago Vallin Spina
 * @date Mar 9, 2016
 *
 * @param files The input file set.
 * @return A CSV file where each row has a path from the file set.
 */
iftCSV *iftFileSetToCSV(const iftFileSet *files);


/**
 * @brief Creates a dictionary from a File Set, where for each, the pair (key, value) corresponds to the
 * File Key and Pathname.
 * The File Key is the filename without extension.
 * E.g: For the pathname: "/user/vader/images/000001.png"
 * (key, value) = ("000001", "/user/vader/images/000001.png")
 *
 *
 * @param  fset The File Set to be used.
 * @return      The resulting Dictionary of Files.
 */
iftDict *iftFileSetToDict(const iftFileSet *fset);


/**
 * @brief Sorts the file list by pathname
 * @author Thiago Vallin Spina
 * @ingroup File
 *
 * @param files The file list.
 */
void iftSortFileSet(iftFileSet *files);


/**
 * @brief Checks if the File Set only have supported images (*.scn, *.pgm, *.ppm, *.png), raising an error if false.
 * @author Samuel Martins
 * @date Aug 12, 2016
 * @ingroup File
 *
 * @param img_paths File Set with Image's pathnames.
 * @exception At least a pathname is not from an Image.
 */
void iftValidateFileSet(const iftFileSet *img_files);


/**
 * @brief Extracts all the samples that belong to a given class in the fileset
 *
 * @author Cesar Castelo
 * @date Mar 14, 2019
 *
 * @param fileset The input fileset.
 * @param label The class label to be extracted.
 * @return A fileset containing the samples from the given class
 */
iftFileSet *iftExtractClassFileset(const iftFileSet *fileset, int label);


#ifdef __cplusplus
}
#endif

#endif //IFT_FILESET_H
