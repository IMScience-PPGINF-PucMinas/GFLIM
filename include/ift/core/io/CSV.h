//
// Created by Samuel Martins on 2018-12-20.
//

#ifndef IFT_CSV_H
#define IFT_CSV_H


#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Comma-Separated Values struct.
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 */
//! swig(destroyer = iftDestroyCSV)
typedef struct ift_csv {
    /** Number of rows of the CSV matrix. */
    long nrows;
    /** Number of columns of the CSV matrix. */
    long ncols;
    /** Header of the CSV. */
    char **header;
    /** CSV matrix of strings. Each string has 512 characters. */
    char ***data;
} iftCSV;


/**
 * @brief Creates an iftCSV (and allocates their string matrix).
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 *
 * @param nrows Number of rows of the CSV matrix.
 * @param ncols Number of columns of the CSV matrix.
 * @return An iftCSV.
 */
//! swig(newobject)
iftCSV *iftCreateCSV(long nrows, long ncols);


/**
 * @brief Destroys an iftCSV.
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 *
 * @param csv Reference to the iftCSV to be destroyed.
 */
void iftDestroyCSV(iftCSV **csv);


/**
 * @brief Reads an iftCSV from file.
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 *
 * Reads an iftCSV from file.
 * The CSV can be a header in the first row. However, it is only indeed considered a header, if
 * all its columns start with letters, and at least one column from the second row is a number.
 *
 * @param csv_pathname The pathname from the CSV to be read.
 * @return The iftCSV.
 */
//! swig(newobject)
iftCSV *iftReadCSV(const char *csv_pathname, const char separator);


/**
 * @brief Writes a CSV structure into a file.
 * @author Samuel Martins
 * @date Feb 23, 2016
 * @ingroup CSV
 *
 * @param csv The original CSV structure.
 * @param filename The output filename.
 * @param separator The selected separating character.
 *
 * @warning This function issues an error if separator != ',' and separator != ';'
 */
//! swig()
void iftWriteCSV(const iftCSV *csv, const char *filename, const char separator);


/**
 * @brief Prints an iftCSV.
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 *
 * @param csv iftCSV to be printed.
 */
//! swig()
void iftPrintCSV(const iftCSV *csv);


/**
 * @brief Merge two CSV files.
 * @author Samuel Martins
 * @date May 18, 2018
 */
//! swig()
iftCSV *iftMergeCSVs(const iftCSV *csv1, const iftCSV *csv2);


/**
 * @brief Copy a CSV.
 * @author Samuel Martins
 * @date Jun 21, 2018
 */
//! swig()
iftCSV *iftCopyCSV(const iftCSV *csv);


/**
 * @brief Insert Header <header> into csv <csv>. The header must have all column names separated by comma.
 * Eg: name,city,address
 * @author Samuel Martins
 * @date Aug 13, 2018
 */
//! swig()
void iftSetCSVHeader(iftCSV *csv, const char *header);


#ifdef __cplusplus
}
#endif

#endif //IFT_CSV_H
