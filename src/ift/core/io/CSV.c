#include "ift/core/io/CSV.h"

#include "ift/core/dtypes/SList.h"
#include "ift/core/io/File.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/Regex.h"
#include "ift/core/tools/String.h"
#include "iftMemory.h"



/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Check if the csv file <csv_pathname> has a header.
 *
 * The first row of the csv is only  considered a header if
 * all its columns start with letters and at least one column from the second row of CSV is a number.
 *
 * @author Samuka Martins
 * @date Jun 4, 2018.
 */
bool _iftHasCSVHeader(const char *csv_pathname, char separator) {
    bool has_header = false;
    
    FILE *fp = fopen(csv_pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "_iftHasCSVHeader", csv_pathname);
    
    // reads the first line for checking if it is a CSV header
    char *line = iftGetLine(fp);
    if (line != NULL) {
        bool is_first_line_ok = true;
        
        // if all columns of first line have letters in the beginning of the string, the row can be the header
        char strSeparator[2] = {separator, '\0'};
        iftSList *SL = iftSplitString(line, strSeparator);
        while (!iftIsSListEmpty(SL)) {
            char *column = iftRemoveSListHead(SL);
            
            if (!iftRegexMatch(column, "^[a-zA-Z]+.*$", separator)) {
                is_first_line_ok = false;
                iftFree(column);
                break;
            }
            iftFree(column);
        }
        iftDestroySList(&SL);
        
        if (is_first_line_ok) {
            iftFree(line);
            
            line = iftGetLine(fp);
            if (line != NULL) {
                iftSList *SL = iftSplitString(line, strSeparator);
                iftFree(line);
                
                while (SL->n != 0) {
                    char *column = iftRemoveSListHead(SL);
                    
                    // if at least one column of the second row is a number (integer or real)
                    // the first row is a header
                    if (iftRegexMatch(column, "^[0-9]+(.[0-9]+)?$", separator)) {
                        iftFree(column);
                        has_header = true;
                        break;
                    }
                    iftFree(column);
                }
                iftDestroySList(&SL);
            }
        }
    }
    fclose(fp);
    
    return has_header;
}



/**
 * @brief Counts the number of Rows and Columns of the CSV File.
 *
 * @author Samuel Martins
 * @date September 16, 2015
 *
 * @param csv_pathname The pathname from the CSV to be read.
 * @param nrows Number of rows which will be figured out (passed by reference).
 * @param ncols Number of columns which will be figured out (passed by reference).
 */
void _iftCountNumOfRowsAndColsFromCSVFile(const char *csv_pathname, long *nrows, long *ncols, char separator) {
    char strSeparator[2] = {separator, '\0'};
    
    FILE *fp = fopen(csv_pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "_iftCountNumOfRowsAndColsFromCSVFile", csv_pathname);
    
    *nrows = 0;
    *ncols = 0;
    
    // reads the first line from the file to get the number of cols
    iftSList *SL = NULL;
    char *line = iftGetLine(fp);
    
    // gets the number of columns from the first line, because such number must be the same for
    // the entire csv
    if (line != NULL) {
        SL = iftSplitString(line, strSeparator);
        (*nrows)++;
        *ncols = SL->n;
        iftDestroySList(&SL);
    }
    
    iftFree(line);
    line = iftGetLine(fp);
    
    // gets each line of the file
    while (line != NULL) {
        SL = iftSplitString(line, strSeparator);
        
        if (*ncols != SL->n)
            iftError("Number of Columns is different in the lines: %d - %d",
                     "_iftCountNumOfRowsAndColsFromCSVFile", *ncols, SL->n);
        
        iftDestroySList(&SL);
        (*nrows)++;
        
        iftFree(line);
        line = iftGetLine(fp);
    }
    
    fclose(fp);
}


/**
 * @brief Creates an iftCSV WITHOUT allocating its String Matrix.
 *
 * @author Samuel Martins
 * @date September 15, 2015
 *
 * Creates an iftCSV. The CSV matrix IS NOT allocated on memory.
 *
 * @param nrows Number of rows of the CSV matrix.
 * @param ncols Number of columns of the CSV matrix.
 * @return An iftCSV.
 */
iftCSV *_iftCreateCSVWithoutStringAllocation(long nrows, long ncols) {
    iftCSV *csv = (iftCSV*) iftAlloc(1, sizeof(iftCSV));
    
    csv->nrows = nrows;
    csv->ncols = ncols;
    
    // allocates the CSV string matrix
    csv->data = (char***) iftAlloc(nrows, sizeof(char**));
    for (long i = 0; i < nrows; i++) {
        csv->data[i] = (char**) iftAlloc(ncols, sizeof(char*));
    }
    
    return csv;
}



/********************** PUBLIC FUNCTIONS *************************/
iftCSV *iftCreateCSV(long nrows, long ncols) {
    iftCSV *csv = (iftCSV*) iftAlloc(1, sizeof(iftCSV));
    csv->header = NULL;
    
    csv->nrows = nrows;
    csv->ncols = ncols;
    
    // allocates the CSV string matrix
    csv->data = (char***) iftAlloc(nrows, sizeof(char**));
    for (long i = 0; i < nrows; i++) {
        csv->data[i] = (char**) iftAlloc(ncols, sizeof(char*));
        
        for (long j = 0; j < ncols; j++) {
            csv->data[i][j] = iftAllocCharArray(IFT_STR_DEFAULT_SIZE);
        }
    }
    
    return csv;
}


void iftDestroyCSV(iftCSV **csv) {
    iftCSV *csv_aux = *csv;
    
    if (csv_aux != NULL) {
        if (csv_aux->data != NULL) {
            // deallocates the CSV string matrix
            for (long i = 0; i < csv_aux->nrows; i++) {
                if (csv_aux->data[i] != NULL)
                    for (long j = 0; j < csv_aux->ncols; j++) {
                        iftFree(csv_aux->data[i][j]);
                    }
                iftFree(csv_aux->data[i]);
            }
        }
        iftFree(csv_aux->data);
        
        if (csv_aux->header != NULL) {
            for (int c = 0; c < csv_aux->ncols; c++)
                iftFree(csv_aux->header[c]);
            iftFree(csv_aux->header);
        }
        iftFree(csv_aux);
        *csv = NULL;
    }
}


iftCSV *iftReadCSV(const char *csv_pathname, const char separator) {
    if (!iftFileExists(csv_pathname))
        iftError("The CSV file pathname \"%s\" does not exists!", "iftReadCSV", csv_pathname);
    
    char strSeparator[2] = {separator, '\0'};
    
    bool has_header = _iftHasCSVHeader(csv_pathname, separator);
    
    long nrows, ncols;
    _iftCountNumOfRowsAndColsFromCSVFile(csv_pathname, &nrows, &ncols, separator);
    if (has_header)
        nrows--;
    
    iftCSV *csv = _iftCreateCSVWithoutStringAllocation(nrows, ncols);
    
    FILE *fp = fopen(csv_pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "_iftCountNumOfRowsAndColsFromCSVFile", csv_pathname);
    
    // copies the values from the CSV file
    iftSList *SL = NULL;
    char *line = iftGetLine(fp);
    
    if (has_header) {
        csv->header = iftAlloc(csv->ncols, sizeof(char*));
        
        SL = iftSplitString(line, strSeparator);
        
        for (long j = 0; j < csv->ncols; j++) {
            csv->header[j] = iftRemoveSListHead(SL); // just points to string
            // removes the '\n' and '\r' from the paths
            iftRightTrim(csv->header[j], '\n');
            iftRightTrim(csv->header[j], '\r');
        }
        iftDestroySList(&SL);
        
        iftFree(line);
        line = iftGetLine(fp);
    }
    
    long i = 0;
    while (line != NULL) {
        SL = iftSplitString(line, strSeparator);
        
        for (long j = 0; j < csv->ncols; j++) {
            csv->data[i][j] = iftRemoveSListHead(SL); // just points to string
            // removes the '\n' and '\r' from the paths
            iftRightTrim(csv->data[i][j], '\n');
            iftRightTrim(csv->data[i][j], '\r');
        }
        i++;
        iftDestroySList(&SL);
        
        iftFree(line);
        line = iftGetLine(fp);
    }
    
    fclose(fp);
    
    return csv;
}


void iftWriteCSV(const iftCSV *csv, const char *filename, const char separator) {
    if (csv != NULL) {
        FILE *f = fopen(filename, "wb");
        
        if(f == NULL)
            iftError("Unable to open file %s for writing!", "iftWriteCSV", filename);
        
        if(separator != ',' && separator != ';')
            iftError("The CSV items may only be separated with ',' or ';'. Please prefer the former whenever possible.",
                     "iftWriteCSV");
        
        if (csv->header != NULL) {
            for (long j = 0; j < (csv->ncols-1); j++) {
                fprintf(f, "%s%c", csv->header[j], separator);
            }
            fprintf(f, "%s\n", csv->header[csv->ncols-1]);
        }
        
        for (long i = 0; i < csv->nrows; i++) {
            for (long j = 0; j < (csv->ncols-1); j++) {
                fprintf(f, "%s%c", csv->data[i][j], separator);
            }
            fprintf(f, "%s\n", csv->data[i][csv->ncols-1]);
        }
        
        fclose(f);
    }
}


void iftPrintCSV(const iftCSV *csv) {
    if (csv != NULL) {
        printf("nrows: %lu, ncols: %lu\n", csv->nrows, csv->ncols);
        
        if (csv->header != NULL) {
            for (long j = 0; j < (csv->ncols-1); j++)
                printf("%s,", csv->header[j]);
            printf("%s\n", csv->header[csv->ncols-1]);
        }
        
        for (long i = 0; i < csv->nrows; i++) {
            for (long j = 0; j < (csv->ncols-1); j++) {
                printf("%s,", csv->data[i][j]);
            }
            printf("%s\n", csv->data[i][csv->ncols-1]);
        }
        puts("");
    }
}


iftCSV *iftMergeCSVs(const iftCSV *csv1, const iftCSV *csv2) {
    if (csv1->ncols != csv2->ncols)
        iftError("CSV files with different number of columns: %d != %d", "iftMergeCSVs", csv1->ncols, csv2->ncols);
    
    iftCSV *csv_merge = iftCreateCSV(csv1->nrows + csv2->nrows, csv1->ncols);
    
    if (csv1->header != NULL) {
        csv_merge->header = iftAlloc(csv_merge->ncols, sizeof(char *));
        #pragma omp parallel for
        for (int c = 0; c < csv1->ncols; c++)
            csv_merge->header[c] = iftCopyString(csv1->header[c]);
    }

    #pragma omp parallel for
    for (int r = 0; r < csv1->nrows; r++)
        for (int c = 0; c < csv1->ncols; c++)
            strcpy(csv_merge->data[r][c], csv1->data[r][c]);
    
    long offset = csv1->nrows;
    #pragma omp parallel for
    for (int r = 0; r < csv2->nrows; r++)
        for (int c = 0; c < csv2->ncols; c++)
            strcpy(csv_merge->data[r+offset][c], csv2->data[r][c]);
    
    return csv_merge;
}



iftCSV *iftCopyCSV(const iftCSV *csv) {
    iftCSV *copy = iftCreateCSV(csv->nrows, csv->ncols);
    
    if (csv->header != NULL) {
        copy->header = iftAlloc(copy->ncols, sizeof(char *));
        #pragma omp parallel for
        for (int c = 0; c < csv->ncols; c++)
            copy->header[c] = iftCopyString(csv->header[c]);
    }

    #pragma omp parallel for
    for (int r = 0; r < csv->nrows; r++)
        for (int c = 0; c < csv->ncols; c++)
            strcpy(copy->data[r][c], csv->data[r][c]);
    
    return copy;
}


void iftSetCSVHeader(iftCSV *csv, const char *header) {
    iftSList *SL = iftSplitString(header, ",");
    if (SL->n != csv->ncols)
        iftError("Number of column names is != the number of column in CSV: %d != %d", "iftSetCSVHeader",
                 SL->n, csv->ncols);
    
    if (csv->header == NULL)
        csv->header = iftAlloc(csv->ncols, sizeof(char *));
    
    for (int c = 0; c < csv->ncols; c++)
        csv->header[c] = iftRemoveSListHead(SL);
    iftDestroySList(&SL);
}
