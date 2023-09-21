#include "ift/core/io/FileSet.h"

#include "ift/core/dtypes/Dict.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/SList.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftCommon.h"


/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Function that is used to compare 2 iftFile in qsort.
 */
int _iftCmpFilesSortFileSet(const void *a, const void *b) {
    iftFile **f1 = (iftFile**) a;
    iftFile **f2 = (iftFile**) b;
    
    return strcmp((*f1)->path, (*f2)->path);
}


/**
 * @brief Gets the Pathnames from a hierarchical level in the dir.
 *
 * @author Samuel Martins
 * @date October 06, 2015
 *
 * Gets the Pathnames from a hierarchical level in the dir.
 * After that, gets the files from its subdirectories recursevely.
 *
 * @param dir The directory to be read.
 * @return The array of iftFiles.
 */
void _iftGetFilesFromDirRec(iftDir *dir, iftSList *SL) {
    for (long i = 0; i < dir->nfiles; i++)
        iftInsertSListIntoTail(SL, dir->files[i]->path);
    
    for (long i = 0; i < dir->nsubdirs; i++)
        _iftGetFilesFromDirRec(dir->subdirs[i], SL);
}


/**
 * @brief Gets the Pathnames with the extension <ext> from a hierarchical level in the dir.
 *
 * @author Samuel Martins
 * @date October 06, 2015
 *
 * Gets the Pathnames with the extension <ext> from a hierarchical level in the dir.
 * After that, gets the files from its subdirectories recursevely.
 *
 * @param dir The directory to be read.
 * @return The array of iftFiles.
 */
void _iftGetFilesFromDirByExtRec(iftDir *dir, iftSList *SL, const char *ext) {
    for (long i = 0; i < dir->nfiles; i++) {
        if (iftEndsWith(dir->files[i]->path, ext))
            iftInsertSListIntoTail(SL, dir->files[i]->path);
    }
    
    for (long i = 0; i < dir->nsubdirs; i++)
        _iftGetFilesFromDirByExtRec(dir->subdirs[i], SL, ext);
}





/********************** PUBLIC FUNCTIONS *************************/
iftFileSet *iftCreateFileSet(long nfiles) {
    iftFileSet *farr = NULL;
    farr        = (iftFileSet *) iftAlloc(1, sizeof(iftFileSet));
    farr->files = (iftFile**) iftAlloc(nfiles, sizeof(iftFile*));
    farr->n = nfiles;
    
    for(long i = 0; i < farr->n; i++)
        farr->files[i] = NULL;
    
    return farr;
}


iftFileSet *iftInitFileSet(long nfiles, ...) {
    iftFileSet *fset = iftCreateFileSet(nfiles);
    
    va_list strings;
    va_start(strings, nfiles);
    for (long i = 0; i < fset->n; i++)
        fset->files[i] = iftCreateFile(va_arg(strings, char*));
    va_end(strings);
    
    return fset;
}


void iftDestroyFileSet(iftFileSet **farr) {
    if (farr != NULL) {
        iftFileSet *faux = *farr;
        
        if (faux != NULL) {
            if (faux->files != NULL) {
                for (long i = 0; i < faux->n; i++)
                    iftDestroyFile(&(faux->files[i]));
                iftFree(faux->files);
            }
            iftFree(faux);
            *farr = NULL;
        }
    }
}


iftFileSet *iftCopyFileSet(const iftFileSet* fileset) {
    if (fileset == NULL)
        return NULL;
    
    iftFileSet* copy = iftCreateFileSet(fileset->n);
    
    for (int i = 0; i < fileset->n; ++i)
        copy->files[i] = iftCopyFile(fileset->files[i]);
    
    return copy;
}


iftFileSet *iftLoadFileSetFromDir(const char *dir_pathname, long hier_level) {
    if (dir_pathname == NULL)
        iftError("Directory \"%s\" is NULL", "iftLoadFileSetFromDir");
    if (!iftDirExists(dir_pathname))
        iftError("Directory \"%s\" does not exist", "iftLoadFileSetFromDir", dir_pathname);
    
    iftDir *dir        = NULL; // loads all files from all dir levels
    iftFileSet *farr   = NULL;
    iftSList *SL       = NULL;
    iftSNode *snode    = NULL, *tnode = NULL;
    
    dir = iftLoadDir(dir_pathname, hier_level);
    SL  = iftCreateSList();
    
    // Gets all files in the entire from the directory and its subdirs
    _iftGetFilesFromDirRec(dir, SL);
    
    /**** Converts the String List into a File Array ****/
    farr        = iftCreateFileSet(SL->n);
    
    // copies each node and destroys it
    snode = SL->head;
    long i = 0;
    while (snode!= NULL) {
        farr->files[i] = iftCreateFile(snode->elem);
        i++;
        
        tnode = snode;
        snode = snode->next;
        iftFree(tnode->elem);
        iftFree(tnode);
    }
    SL->head = SL->tail = NULL;
    iftDestroySList(&SL);
    /*******************/
    
    iftDestroyDir(&dir);
    
    return farr;
}


iftFileSet *iftLoadFileSetFromDirBySuffix(const char *dir_pathname, const char *ext, int depth) {
    if (!iftDirExists(dir_pathname))
        iftError("Directory \"%s\" does not exist", "iftReadFilesFromdDirectory", dir_pathname);
    
    iftDir *dir      = NULL; // loads all files from all dir levels
    iftFileSet *farr = NULL;
    iftSList *SL     = NULL;
    iftSNode *snode  = NULL, *tnode = NULL;
    
    dir = iftLoadDir(dir_pathname, depth);
    SL  = iftCreateSList();
    
    // Gets all files with extension <ext> in the entire from the directory and its subdirs
    _iftGetFilesFromDirByExtRec(dir, SL, ext);
    
    /**** Converts the String List into a File Array ****/
    farr        = iftCreateFileSet(SL->n);
    
    // copies each node and destroys it
    snode = SL->head;
    long i = 0;
    while (snode!= NULL) {
        farr->files[i] = iftCreateFile(snode->elem);
        i++;
        
        // destroys node
        tnode = snode;
        snode = snode->next;
        iftFree(tnode->elem);
        iftFree(tnode);
    }
    SL->head = SL->tail = NULL;
    iftDestroySList(&SL);
    /*******************/
    
    iftDestroyDir(&dir);
    
    return farr;
}


iftFileSet *iftLoadFileSetFromDirByRegex(const char *dir_pathname, const char *regex, bool sort_pathnames) {
    if (dir_pathname == NULL)
        iftError("Dir's pathname is NULL", "iftLoadFilesFromDirByRegex");
    if (regex == NULL)
        iftError("Regex is NULL", "iftLoadFilesFromDirByRegex");
    if (!iftDirExists(dir_pathname))
        iftError("Directory \"%s\" does not exist", "iftReadFilesFromdDirectory", dir_pathname);
    
    iftDir *dir      = iftLoadFilesFromDirByRegex(dir_pathname, regex);
    iftFileSet *farr = iftCreateFileSet(dir->nfiles);
    
    for (long i = 0; i < dir->nfiles; i++)
        farr->files[i] = iftCopyFile(dir->files[i]);
    
    if (sort_pathnames)
        iftSortFileSet(farr);
    
    iftDestroyDir(&dir);
    
    return farr;
}


iftFileSet *iftLoadFileSetFromCSV(const char *csv_pathname, bool sort_pathnames) {
    iftCSV *csv   = iftReadCSV(csv_pathname, ',');
    long nfiles = csv->nrows*csv->ncols;
    
    iftFileSet *farr = iftCreateFileSet(nfiles);

    #pragma omp parallel for
    for (long i = 0; i < csv->nrows; i++)
        for (long j = 0; j < csv->ncols; j++) {
            long p = j + i*csv->ncols;
            char *aux = iftExpandUser(csv->data[i][j]);
            farr->files[p] = iftCreateFile(aux);
            iftFree(aux);
        }
    iftDestroyCSV(&csv);
    
    if (sort_pathnames)
        iftSortFileSet(farr);
    
    return farr;
}


iftFileSet *iftLoadFileSetFromDirOrCSV(const char *file_entry, long hier_levels, bool sort_pathnames) {
    iftFileSet *fset = NULL;
    if (iftDirExists(file_entry))
        fset = iftLoadFileSetFromDir(file_entry, hier_levels); // it also returns a sorted list
    else {
        char *lower_file_entry = iftLowerString(file_entry);
        
        if (iftFileExists(file_entry) && iftEndsWith(lower_file_entry, ".csv")) {
            fset = iftLoadFileSetFromCSV(file_entry, sort_pathnames);
            iftFree(lower_file_entry);
            
            // if (sort_pathnames)
            //     iftSortFileSet(fset);
        }
        else
            iftError("Invalid File Entry: %s\nIt is neither a directory nor a CSV file",
                     "iftLoadFileSetFromDirOrCSV", file_entry);
    }
    
    return fset;
}


void iftWriteFileSetAsCSV(const iftFileSet *fset, const char *format, ...) {
    va_list args;
    char path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(path, format, args);
    va_end(args);
    
    iftCSV *csv = iftFileSetToCSV(fset);
    iftWriteCSV(csv, path, ',');
    iftDestroyCSV(&csv);
}


void iftRemoveFileSet(const iftFileSet *fset) {
    if (fset != NULL) {
        for (long i = 0; i < fset->n; i++)
            iftRemoveFile(fset->files[i]->path);
    }
}


iftFileSet *iftMergeFileSet(const iftFileSet *fset1, const iftFileSet *fset2) {
    if (fset1 == NULL)
        iftError("File Set 1 is NULL", "iftMergeFileSet");
    if (fset2 == NULL)
        iftError("File Set 2 is NULL", "iftMergeFileSet");
    
    iftFileSet *fset3 = iftCreateFileSet(fset1->n + fset2->n);
    
    for (long i = 0; i < fset1->n; i++) {
        iftDestroyFile(&fset3->files[i]);
        fset3->files[i] = iftCopyFile(fset1->files[i]);
    }
    
    for (long j = 0; j < fset2->n; j++) {
        long i = j + fset1->n;
        
        iftDestroyFile(&fset3->files[i]);
        fset3->files[i] = iftCopyFile(fset2->files[j]);
    }
    
    return fset3;
}


iftFileSet *iftFilterFileSetByFilename(const iftFileSet *src, const iftFileSet *dst) {
    iftFileSet *filt_set = iftCreateFileSet(src->n);

    #pragma omp parallel for
    for (long i = 0; i < src->n; i++) {
        const char *img_path = src->files[i]->path;
        char *img_filename   = iftFilename(img_path, iftFileExt(img_path));
        
        for (long j = 0; j < dst->n; j++) {
            const char *label_path = dst->files[j]->path;
            char *label_filename   = iftFilename(label_path, iftFileExt(label_path));
            
            if (iftCompareStrings(img_filename, label_filename)) {
                filt_set->files[i] = iftCreateFile(label_path);
            }
            
            iftFree(label_filename);
        }
        iftFree(img_filename);
        
        // not found
        if (filt_set->files[i] == NULL) {
            iftError("Not found the corresponding Label Image for the Image \"%s\"",
                     "iftBuildLabelSetFromImageSet", img_path);
        }
    }
    
    return filt_set;
}


iftFileSet *iftFilterFileSetByIndexes(const iftFileSet *fset, const iftIntArray *idxs) {
    iftFileSet *filt_set = iftCreateFileSet(idxs->n);

    #pragma omp parallel for
    for (int i = 0; i < idxs->n; i++)
        filt_set->files[i] = iftCreateFile(fset->files[idxs->val[i]]->path);
    
    return filt_set;
}


int iftFindFileByLabel(const iftFileSet *files, int label) {
    int i, file = IFT_NIL;
    
    for(i = 0; i < files->n && file == IFT_NIL; i++) {
        if(files->files[i]->label == label)
            file = i;
    }
    
    return file;
}


int iftFindFileByLabelAndSample(const iftFileSet *files, int label, int sample) {
    int i, file = IFT_NIL;
    
    for(i = 0; i < files->n && file == IFT_NIL; i++) {
        if(files->files[i]->label == label && files->files[i]->sample == sample)
            file = i;
    }
    
    return file;
}


int iftFindFileByBasename(const iftFileSet *files, const char *bname, const char *suffix) {
    int i, file = IFT_NIL;
    
    for(i = 0; i < files->n && file == IFT_NIL; i++) {
        char *cur_bname = iftFilename(files->files[i]->path, suffix);
        if(iftCompareStrings(bname, cur_bname))
            file = i;
        
        iftFree(cur_bname);
    }
    
    return file;
}


int iftFileSetLabelsNumber(const iftFileSet* fileset) {
    int* labels = iftFileSetLabels(fileset);
    iftIntArray *unique = iftIntArrayUnique(labels, fileset->n);
    int nlabels = iftMaxValIntArray(unique->val, unique->n, NULL);
    iftFree(labels);
    
    return nlabels;
}


int* iftCountSamplesPerClassFileSet(iftFileSet *fileset)
{
    int nClasses = iftFileSetLabelsNumber(fileset);
    int *sampPerClass = iftAllocIntArray(nClasses+1);
    
    for(int i = 0; i < fileset->n; i++) {
        iftUpdateDataSetFileInfo(fileset->files[i]);
        if (fileset->files[i]->label > 0)
            sampPerClass[fileset->files[i]->label]++;
    }
    
    return sampPerClass;
}


int* iftFileSetLabels(const iftFileSet* fileSet) {
    int n = fileSet->n;
    int i;
    int *labelArray = (int*) iftAlloc(n, sizeof(int));
    
    for (i = 0; i < n; ++i) {
        labelArray[i] = fileSet->files[i]->label;
    }
    
    return labelArray;
}


iftCSV *iftFileSetToCSV(const iftFileSet *files) {
    iftCSV *csv = iftCreateCSV(files->n, 1);
    
    for(int i = 0; i < files->n; i++) {
        /** Freeing the standard-size string allocated in iftCreateCSV just to ensure that the string will have the right size **/
        iftFree(csv->data[i][0]);
        csv->data[i][0] = iftCopyString(files->files[i]->path);
    }
    
    return csv;
}


iftDict *iftFileSetToDict(const iftFileSet *fset) {
    if (fset == NULL)
        iftError("File Set is NULL", "iftFileSetToDict");
    
    iftDict *dict = iftCreateDictWithApproxSize(fset->n);
    
    for (long i = 0; i < fset->n; i++) {
        char *path = fset->files[i]->path;
        char *fkey = iftFilename(path, iftFileExt(path));
        iftInsertIntoDict(fkey, path, dict);
        iftFree(fkey);
    }
    
    return dict;
}


void iftSortFileSet(iftFileSet *files) {
    qsort(files->files, files->n, sizeof(iftFile*), _iftCmpFilesSortFileSet);
}


void iftValidateFileSet(const iftFileSet *img_files) {
    if (img_files->n == 0)
        iftError("There are no Files in the File Set", "iftValidateFileSet");
    else {
        #pragma omp parallel for
        for (long i = 0; i < img_files->n; i++)
            if (!iftIsImageFile(img_files->files[i]->path))
                iftError("File \"%s\" is not an Image.\nCheck your File Set", "iftValidateFileSet",
                         img_files->files[i]->path);
    }
}

iftFileSet *iftExtractClassFileset(const iftFileSet *fileset, int label)
{
    int nSamp = 0;
    for(int i = 0; i < fileset->n; i++) {
        iftUpdateDataSetFileInfo(fileset->files[i]);
        if(fileset->files[i]->label == label)
            nSamp++;
    }

    if(nSamp == 0)
        iftError("There are no samples from the selected class in the fileset", "iftExtractClassFileset");

    iftFileSet *newFileset = iftCreateFileSet(nSamp);

    int s = 0;
    for(int i = 0; i < fileset->n; i++) {
        if(fileset->files[i]->label == label) {
            newFileset->files[s] = iftCopyFile(fileset->files[i]);
            s++;
        }
    }

    return newFileset;
}
