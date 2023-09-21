#include "ift/core/io/Dir.h"

#include "ift/core/tools/Regex.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftMemory.h"
#include <dirent.h>


/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Function that is used to compare 2 iftFile in qsort.
 */
int _iftCmpFiles(const void *a, const void *b) {
    iftFile **f1 = (iftFile**) a;
    iftFile **f2 = (iftFile**) b;
    
    return strcmp((*f1)->path, (*f2)->path);
}


/**
 * @brief Function that is used to compare 2 iftFile in qsort.
 */
int _iftCmpDirs(const void *a, const void *b) {
    iftDir **dir1 = (iftDir**) a;
    iftDir **dir2 = (iftDir**) b;
    
    return strcmp((*dir1)->path, (*dir2)->path);
}


/**
 * @brief Counts the Number of Files and SubDirs in the directory <pathname>.
 *
 * This function counts the number of files and subdirs in the directory <pathname>.
 *
 * @param pathname The directory to be read.
 * @param n Reference that will stores the number of files from <pathname>.
 * @param nsubdirs Reference that will stores the number of files from <pathname>.
 */
void _iftCountFilesInDirectory(const char *dir_pathname, long *nfiles, long *nsubdirs) {
    //http://pubs.opengroup.org/onlinepubs/007908799/xsh/dirent.h.html
    //http://www.delorie.com/gnu/docs/glibc/libc_270.html
    DIR *system_dir;
    struct dirent *entry;
    char msg[512];
    char *pathname = NULL;
    *nfiles = 0;
    *nsubdirs = 0;
    
    system_dir = opendir(dir_pathname);
    if (system_dir) {
        while ((entry = readdir(system_dir)) != NULL)
            // it excludes the system_dir . and ..
            if ((strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0)) {
                pathname = iftJoinPathnames(2, dir_pathname, entry->d_name);
                
                if (iftDirExists(pathname)) {
                    (*nsubdirs)++;
                }
                else {
                    (*nfiles)++;
                }
                
                iftFree(pathname);
                pathname = NULL;
            }
        closedir(system_dir);
    }
    else {
        sprintf(msg, "Error opening directory path: \"%s\"", dir_pathname);
        iftError(msg, "_iftCountFilesInDirectory");
    }
}


/**
 * @brief Counts ONLY the Number of Files with (not subdirs) <extension) in the directory <dir_pathname>.
 *
 * Counts ONLY the Number of Files with (not subdirs) <extension) in the directory <dir_pathname>.
 * If <extension> = "", the function gets all the files in the directory.
 *
 * @param dir_pathname The directory to be read.
 * @param extension Extension from the files.
 * @return Number of files (not subdirs) in the directory.
 */
int _iftCountOnlyFilesInDirectory(const char *dir_pathname, const char *extension) {
    DIR *dir;
    struct dirent *entry;
    
    long file_count = 0;
    
    dir = opendir(dir_pathname);
    if (dir) {
        while ((entry = readdir(dir)) != NULL)
            if (!iftDirExists(entry->d_name) && (iftEndsWith(entry->d_name, extension)))
                file_count++;
        closedir(dir);
    }
    else {
        char msg[512];
        sprintf(msg, "Error opening directory path: \"%s\"", dir_pathname);
        iftError(msg, "_iftCountOnlyFilesInDirectory");
    }
    
    return file_count;
}


/**
 * @brief Lists all childs (files/subdirs) from the curr_level.
 *
 * Lists all childs (files/subdirs) from the <curr_level> and calls itself for the subdirs recursevely.
 * It is used by the method iftListDirectory.
 *
 * @param dir The iftDir of the root directory.
 * @param hier_levels Number of levels which will be loaded.
 * @param curr_level Level to be loaded.
 */
void _iftListDirectoryRec(iftDir *dir, long hier_levels, long curr_level) {
    DIR *system_dir;
    struct dirent *entry;
    char *pathname = NULL;
    
    dir->files   = NULL;
    dir->subdirs = NULL;
    
    if (curr_level <= hier_levels) {
        system_dir = opendir(dir->path);
        if (system_dir) {
            _iftCountFilesInDirectory(dir->path, &dir->nfiles, &dir->nsubdirs);
            if (dir->nfiles != 0)
                dir->files   = (iftFile**) iftAlloc(dir->nfiles, sizeof(iftFile*));
            if (dir->nsubdirs != 0)
                dir->subdirs = (iftDir**) iftAlloc(dir->nsubdirs, sizeof(iftDir*));
            
            long i = 0, j = 0;
            while ((entry = readdir(system_dir)) != NULL) {
                // it excludes the dir . and ..
                if ((strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0)) {
                    pathname = iftJoinPathnames(2, dir->path, entry->d_name);
                    
                    if (iftDirExists(pathname)) { // it is a directory
                        iftDir *subdir = (iftDir*) iftAlloc(1, sizeof(iftDir));
                        subdir->path = pathname;
                        
                        subdir->nfiles   = 0;
                        subdir->nsubdirs = 0;
                        subdir->files    = NULL;
                        subdir->subdirs  = NULL;
                        
                        _iftListDirectoryRec(subdir, hier_levels, curr_level+1);
                        dir->subdirs[j++] = subdir;
                    }
                    else { // it is a File
                        iftFile *f = (iftFile*) iftAlloc(1, sizeof(iftFile));
                        f->path = pathname;
                        dir->files[i++] = f;
                        f->suffix = NULL;
                    }
                }
            }
            closedir(system_dir);
            
            /* sorts the pathnames using qsort functions */
            qsort(dir->files, dir->nfiles, sizeof(iftFile*), _iftCmpFiles);
            qsort(dir->subdirs, dir->nsubdirs, sizeof(iftDir*), _iftCmpDirs);
        }
        else {
            char msg[512];
            sprintf(msg, "Error opening directory path: \"%s\"", dir->path);
            iftError(msg, "_iftListDirectoryRec");
        }
    }
}


/**
 * @brief Lists all childs (files/subdirs) from the root directory <dir> until the level <hier_level>.
 *
 * This function lists all childs (files/subdirs) from the root directory <dir> until the level <hier_level>.
 * The parameter <hier_levels> indicates until which hierarchical level will be considered in the
 * loading.
 *
 * The loaded files and subdirs are stored in the iftDir dir->files, dir->subdirs.
 *
 * @param root_dir The iftDir of the root directory.
 * @param hier_levels Number of levels which will be loaded.
 */
void _iftListDirectory(iftDir *root_dir, long hier_levels) {
    if (root_dir == NULL) {
        iftError("Directory is NULL", "_iftListDirectory");
    }
    else {
        if (hier_levels == 0)
            hier_levels = IFT_INFINITY_INT; // trick to set the hier_levels as the possible maximum
        
        root_dir->nfiles   = 0;
        root_dir->nsubdirs = 0;
        root_dir->files    = NULL;
        root_dir->subdirs  = NULL;
        
        long curr_level = 1;
        _iftListDirectoryRec(root_dir, hier_levels, curr_level);
    }
}


/**
 * @brief Lists all files (not subdirs) in the directory <dir_pathname> with a given <suffix>.
 *
 * This function lists all files (not subdirs) in the directory <dir_pathname> of a given <suffix>.
 * The filenames are sorted in ascending order.
 * It only gets the files from the 1st file hierarchical level.
 * If <suffix> = "", it gets all files (not subdirs).
 *
 * @param dir_pathname The directory to be read.
 * @param nchilds The number of files (not subdirs) from the <pathname> is returned in this reference.
 * @param suffix suffix from the files (with extension).
 * @return A vector of iftFiles* with all files/subdirs inside <dir_pathname>.
 */
void _iftListFilesFromDirectory(iftDir *dir, const char *suffix) {
    DIR *system_dir;
    struct dirent *entry;
    char *pathname = NULL;
    
    if (dir != NULL) {
        dir->files    = NULL;
        dir->subdirs  = NULL;
        dir->nfiles   = 0;
        dir->nsubdirs = 0;
        
        system_dir = opendir(dir->path);
        if (system_dir) {
            dir->nfiles = _iftCountOnlyFilesInDirectory(dir->path, suffix);
            
            if (dir->nfiles != 0) {
                dir->files = (iftFile**) iftAlloc(dir->nfiles, sizeof(iftFile*));
                
                long i = 0;
                while ((entry = readdir(system_dir)) != NULL) {
                    pathname = iftJoinPathnames(2, dir->path, entry->d_name);
                    
                    if (!iftDirExists(pathname) && (iftEndsWith(pathname, suffix))) {
                        iftFile *f = iftCreateFile(pathname);
                        dir->files[i++] = f;
                    }
                }
                closedir(system_dir);
                
                /* sorts the pathnames using qsort functions */
                iftSortDir(dir);
            }
        }
        else {
            char msg[512];
            sprintf(msg, "Error opening directory path: \"%s\"", dir->path);
            iftError(msg, "_iftListFilesFromDirectory");
        }
    }
}


/**
 * @brief Prints all childs (files/sub-directories) as a Tree.
 *
 * This function prints all childs (files/sub-directories) from a directory.
 * It also computes the total number of files and subdirs from this vector.
 *
 * @param dir The parent directory to be printed.
 * @param node_level The hierarchic level from the files. Initially, the nodel_level is expected to be 1.
 * @param tn_files Reference to the variable that will store the total number of files
 * from this file hierarchy. Initially, n is expected to be 0
 * @param tn_subdirs Reference to the variable that will store the total number of sub-directories
 * from this file hierarchy. Initially, nsubdirs is expected to be 0.
 */
void _iftPrintPathnameAsNode(const iftDir *dir, long node_level, long *tn_files,
                            long *tn_subdirs) {
    if (dir != NULL) {
        // prints the sub-directories
        for (long i = 0; i < dir->nsubdirs; i++) {
            for (long j = 1; j < node_level; j++)
                printf("│   ");
            
            if (i == (dir->nsubdirs-1) && (dir->nfiles == 0))
                printf("└── ");
            else
                printf("├── ");
            
            char *pathname = iftSplitStringAt(dir->subdirs[i]->path, IFT_SEP_C, -1);
            printf(IFT_COLOR_BLUE "%s\n" IFT_COLOR_RESET, pathname);
            iftFree(pathname);
            (*tn_subdirs)++;
            
            _iftPrintPathnameAsNode(dir->subdirs[i], node_level+1, tn_files, tn_subdirs);
        }
        
        // prints the files
        for (long i = 0; i < dir->nfiles; i++) {
            for (long j = 1; j < node_level; j++)
                printf("│   ");
            
            if (i == (dir->nfiles-1))
                printf("└── ");
            else
                printf("├── ");
            
            char *pathname = iftSplitStringAt(dir->files[i]->path, IFT_SEP_C, -1);
            printf("%s\n", pathname);
            iftFree(pathname);
            (*tn_files)++;
        }
    }
}


/**
 * @brief Removes a directory represented by the iftDir <code><b>dir</code></b> from the <b>DISK</b>.
 * @author Samuel Martins
 * @date Oct 8th, 2015
 * @ingroup File
 *
 * @param dir The iftDir that represents the directory to be removed from the DISK.
 *
 * @warning If the pathname is a file or it does not exist, nothing will be done.
 * @attention The iftDir object/struct IS NOT DESTROYED, then it must be done out this function.
 */
void _iftRemoveDirFromDirObject(iftDir *dir) {
    if (iftPathnameExists(dir->path)) {
        if (iftFileExists(dir->path))
            iftWarning("The directory \"%s\" is a FILE... Nothing to do", "_iftRemoveDirFromDirObject", dir->path);
        else {
            // Removes all files inside the directory
            for (long i = 0; i < dir->nfiles; i++)
                iftRemoveFile(dir->files[i]->path);
            
            // Removes the subdirectories
            for (long i = 0; i < dir->nsubdirs; i++)
                _iftRemoveDirFromDirObject(dir->subdirs[i]);
            
            
            /**
             * If the directory <dir->path> is into an external HD, shared folder by sshfs, etc,
             * the OS can create hidden files inside it to indicate this (Eg: the files .fuse_hiddenXXXXXX).
             * This creation can happen after the loading of the directory into iftDir.
             * I (Samuel) tried a lot of things (rmdir, rm -rf, etc) and nothing worked.
             * It seems that temporary directories created by your program have these hidden files
             * only while your program runs.
             */
            if (rmdir(dir->path) == -1)
                iftWarning("Problem to remove the directory: \"%s\"... It is likely to have files and/or subdirs inside it.",
                           "_iftRemoveDirFromDirObject", dir->path);
        }
    }
}




/********************** PUBLIC FUNCTIONS *************************/
void iftDestroyDir(iftDir **dir) {
    if (dir != NULL) {
        iftDir *dir_aux = *dir;
        
        if (dir_aux != NULL) {
            if (dir_aux->path != NULL)
                iftFree(dir_aux->path);
            
            // deallocates the files
            if (dir_aux->files != NULL) {
                for (long i = 0; i < dir_aux->nfiles; i++)
                    iftDestroyFile(&dir_aux->files[i]);
                
                iftFree(dir_aux->files);
                dir_aux->files = NULL;
            }
            
            if (dir_aux->subdirs != NULL) {
                // deallocates the subdirs
                for (long j = 0; j < dir_aux->nsubdirs; j++)
                    iftDestroyDir(&dir_aux->subdirs[j]);
                
                iftFree(dir_aux->subdirs);
                dir_aux->subdirs = NULL;
            }
            iftFree(dir_aux);
            *dir = NULL;
        }
    }
}


void iftPrintDirInfo(const iftDir *dir) {
    if (dir != NULL) {
        printf("Pathname: %s\n", dir->path);
        printf("Num of Files: %lu\n", dir->nfiles);
        printf("Num of Subdirs: %lu\n", dir->nsubdirs);
    }
}


void iftPrintDirAsTree(const iftDir *dir) {
    printf(IFT_COLOR_BLUE "\n%s\n" IFT_COLOR_RESET, dir->path);
    long n_total_files = 0, n_total_subdirs = 0;
    _iftPrintPathnameAsNode(dir, 1, &n_total_files, &n_total_subdirs);
    printf("\n%lu sub-directorie(s), %lu file(s)\n\n", n_total_subdirs, n_total_files);
}


char *iftMakeTempDir(const char *prefix, const char *suffix, const char *dir_name) {
    char *dir_pathname = iftMakeTempPathname(prefix, suffix, dir_name);
    iftMakeDir(dir_pathname);
    return dir_pathname;
}


void iftCopyDirFromDisk(const char *src, const char *dst) {
    iftDir *src_dir = iftLoadDir(src, 0); // loads the dir and all its file hierarchy
    char *dir_pathname = NULL;
    
    if (iftDirExists(dst)) {
        char *src_basename = iftFilename(src, NULL);
        dir_pathname       = iftJoinPathnames(2, dst, src_basename);
        iftFree(src_basename);
    }
        // creates the destination directory
    else {
        dir_pathname = iftCopyString(dst);
    }
    
    iftMakeDir(dir_pathname);
    
    // copies the files inside the src directory
    for (long i = 0; i < src_dir->nfiles; i++)
        iftCopyFromDisk(src_dir->files[i]->path, dir_pathname);
    
    // copies the subdirs inside the src directory
    for (long i = 0; i < src_dir->nsubdirs; i++)
        iftCopyDirFromDisk(src_dir->subdirs[i]->path, dir_pathname);
    
    iftFree(dir_pathname);
    iftDestroyDir(&src_dir);
    puts("");
}


void iftRemoveDir(const char *dir_pathname) {
    if (dir_pathname != NULL) {
        if (iftPathnameExists(dir_pathname)) {
            if (iftFileExists(dir_pathname))
                iftWarning("The directory \"%s\" is a FILE... Nothing to do", "iftRemoveDir", dir_pathname);
            else {
                iftDir *dir = iftLoadDir(dir_pathname, 0);
                _iftRemoveDirFromDirObject(dir);
                iftDestroyDir(&dir);
            }
        }
    }
}


iftDir *iftLoadDir(const char *dir_pathname, long hier_levels) {
    char msg[512];
    iftDir *dir = NULL;
    
    
    if (iftPathnameExists(dir_pathname)) {
        // it is really a directory and it exists
        if (iftDirExists(dir_pathname)) {
            dir = (iftDir*) iftAlloc(1, sizeof(iftDir));
            dir->path = iftAllocCharArray(strlen(dir_pathname) + 2); // one more char to put the separation '/'
            strcpy(dir->path, dir_pathname);
            
            // puts the '/' at the end of the pathname
            if (dir->path[strlen(dir->path) - 1] != IFT_SEP_C[0])
                strcat(dir->path, IFT_SEP_C);
            
            _iftListDirectory(dir, hier_levels);
        }
            // it is a File instead of a Directory
        else {
            sprintf(msg, "Pathname \"%s\" is a File", dir_pathname);
            iftError(msg, "iftLoadDir");
        }
    }
    else {
        sprintf(msg, "Pathname \"%s\" does not exist!", dir_pathname);
        iftError(msg, "iftLoadDir");
    }
    
    return dir;
}


iftDir *iftLoadFilesFromDirBySuffix(const char *dir_pathname, const char *suffix) {
    char msg[512];
    
    iftDir *dir = NULL;
    
    if (iftPathnameExists(dir_pathname)) {
        // it is really a directory and it exists
        if (iftDirExists(dir_pathname)) {
            dir = (iftDir*) iftAlloc(1, sizeof(iftDir));
            dir->path = iftAllocCharArray(strlen(dir_pathname) + 2); // one more char to put the separation '/'
            strcpy(dir->path, dir_pathname);
            
            // puts the '/' at the end of the pathname
            if (dir->path[strlen(dir->path) - 1] != IFT_SEP_C[0])
                strcat(dir->path, IFT_SEP_C);
            
            _iftListFilesFromDirectory(dir, suffix);
        }
        else { // it is really a directory and it exists
            sprintf(msg, "\"%s\" is a file!", dir->path);
            iftError(msg, "iftLoadFilesFromDirBySuffix");
        }
    }
    else {
        sprintf(msg, "Pathname \"%s\" does not exist!", dir_pathname);
        iftError(msg, "iftLoadFilesFromDirBySuffix");
    }
    
    return dir;
}


iftDir *iftLoadFilesFromDirByRegex(const char *dir_pathname, const char *regex) {
    if (dir_pathname == NULL)
        iftError("Dir's Pathname is NULL", "iftLoadFilesFromDirByRegex");
    if (regex == NULL)
        iftError("Regex is NULL", "iftLoadFilesFromDirByRegex");
    
    iftDir *dir          = iftLoadDir(dir_pathname, 1);
    long n_all_files   = dir->nfiles;
    long n_final_files = 0;
    
    for (long f = 0; f < n_all_files; f++) {
        char *filename = iftFilename(dir->files[f]->path, NULL);
        
        if (iftRegexMatch(filename, regex))
            n_final_files++;
        
        iftFree(filename);
    }
    
    iftFile **file_array = dir->files;
    dir->files           = NULL;
    
    dir->files  = (iftFile**) iftAlloc(n_final_files, sizeof(iftFile*));
    dir->nfiles = n_final_files;
    
    long i = 0;
    for (long f = 0; f < n_all_files; f++) {
        char *filename = iftFilename(file_array[f]->path, NULL);
        
        if (iftRegexMatch(filename, regex)) {
            dir->files[i++] = file_array[f];
            file_array[f]   = NULL;
        }
        else iftDestroyFile(&file_array[f]);
        
        iftFree(filename);
    }
    iftFree(file_array);
    
    return dir;
}


void iftSortDir(iftDir *files) {
    qsort(files->files, files->nfiles, sizeof(iftFile*), _iftCmpFiles);
}


int *iftImageLabels(iftDir *dir) {
    
    int n = dir->nfiles;
    int i;
    int *labelArray = (int*) iftAlloc(n, sizeof(int));
    
    for (i = 0; i < n; ++i) {
        labelArray[i] = iftImageSampleLabel(dir->files[i]->path);
    }
    
    return labelArray;
}


bool iftDirExists(const char *format, ...) {
    va_list args;
    char pathname[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(pathname, format, args);
    va_end(args);
    
    struct stat st;
    
    if (stat(pathname, &st) == 0) {
        if (S_ISDIR(st.st_mode))
            return true; //it's a directory
    }
    return false ;
}


char *iftParentDir(const char *pathname) {
    char *filename   = NULL;
    char *parent_dir = NULL;
    
    filename   = iftSplitStringAt(pathname, IFT_SEP_C, -1);
    parent_dir = iftRemoveSuffix(pathname, filename);
    iftFree(filename);
    
    // if the parent_dir is empty
    if (strcmp(parent_dir, "") == 0) {
        strcpy(parent_dir, ".");
    }
    else { // eliminates the slash (dir. separator char) at the end
        parent_dir[strlen(parent_dir)-1] = '\0';
    }
    
    return (parent_dir);
}


void iftMakeDir(const char *dir_path) {
    if (!iftDirExists(dir_path)) {
        char *parent_dir = iftAllocCharArray(IFT_STR_DEFAULT_SIZE);
        strcpy(parent_dir, "");
        
        iftSList *SL    = iftSplitString(dir_path, IFT_SEP_C);
        char *inter_dir = iftRemoveSListHead(SL);
        
        while (inter_dir != NULL) {
            strcat(parent_dir, inter_dir);
            strcat(parent_dir, IFT_SEP_C);
            
            if (!iftDirExists(parent_dir)) {
            #if defined(__linux) || defined(__APPLE__)
                if (mkdir(parent_dir, 0777) == -1) // Create the directory
                    iftError("Problem to create the directory: %s", "iftMakeDir", dir_path);
            #else
                if (!CreateDirectory(parent_dir, NULL))
                        iftError("Problem to create the directory", "iftMakeDir");
            #endif
            }
            
            iftFree(inter_dir);
            inter_dir = iftRemoveSListHead(SL);
        }
        iftDestroySList(&SL);
        iftFree(parent_dir);
    }
}


