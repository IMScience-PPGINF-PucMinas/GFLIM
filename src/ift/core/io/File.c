#include "ift/core/io/File.h"

#include "ift/core/dtypes/SList.h"
#include "ift/core/io/Dir.h"
#include "ift/core/tools/OS.h"
#include "ift/core/tools/Regex.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftCommon.h"
#include "iftDataSet.h"
#include "iftMemory.h"



#if defined(__WIN32) || defined(__WIN64)
/* Implementation of mkstemp for Windows from http://stackoverflow.com/questions/6036227/mkstemp-implementation-for-win32 */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

/* mkstemp extracted from libc/sysdeps/posix/tempname.c.  Copyright
   (C) 1991-1999, 2000, 2001, 2006 Free Software Foundation, Inc.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.  */

static const char letters[] =
"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

/* Generate a temporary file name based on TMPL.  TMPL must match the
   rules for mk[s]temp (i.e. end in "XXXXXX").  The name constructed
   does not exist at the time of the call to mkstemp.  TMPL is
   overwritten with the result.  */
int
mkstemp_custom (char *tmpl)
{
  int len;
  char *XXXXXX;
  static unsigned long long value;
  unsigned long long random_time_bits;
  unsigned int count;
  int fd = -1;
  int save_errno = errno;

  /* A lower bound on the number of temporary files to attempt to
     generate.  The maximum total number of temporary file names that
     can exist for a given template is 62**6.  It should never be
     necessary to try all these combinations.  Instead if a reasonable
     number of names is tried (we define reasonable as 62**3) fail to
     give the system administrator the chance to remove the problems.  */
#define ATTEMPTS_MIN (62 * 62 * 62)

  /* The number of times to attempt to generate a temporary file.  To
     conform to POSIX, this must be no smaller than TMP_MAX.  */
#if ATTEMPTS_MIN < TMP_MAX
  unsigned int attempts = TMP_MAX;
#else
  unsigned int attempts = ATTEMPTS_MIN;
#endif

  len = strlen (tmpl);
  if (len < 6 || strcmp (&tmpl[len - 6], "XXXXXX"))
    {
      errno = EINVAL;
      return -1;
    }

/* This is where the Xs start.  */
  XXXXXX = &tmpl[len - 6];

  /* Get some more or less random data.  */
  {
    SYSTEMTIME      stNow;
    FILETIME ftNow;

    // get system time
    GetSystemTime(&stNow);
    stNow.wMilliseconds = 500;
    if (!SystemTimeToFileTime(&stNow, &ftNow))
    {
        errno = -1;
        return -1;
    }

    random_time_bits = (((unsigned long long)ftNow.dwHighDateTime << 32)
                        | (unsigned long long)ftNow.dwLowDateTime);
  }
  value += random_time_bits ^ (unsigned long long)GetCurrentThreadId ();

  for (count = 0; count < attempts; value += 7777, ++count)
    {
      unsigned long long v = value;

      /* Fill in the random bits.  */
      XXXXXX[0] = letters[v % 62];
      v /= 62;
      XXXXXX[1] = letters[v % 62];
      v /= 62;
      XXXXXX[2] = letters[v % 62];
      v /= 62;
      XXXXXX[3] = letters[v % 62];
      v /= 62;
      XXXXXX[4] = letters[v % 62];
      v /= 62;
      XXXXXX[5] = letters[v % 62];

      fd = open (tmpl, O_RDWR | O_CREAT | O_EXCL, _S_IREAD | _S_IWRITE);
      if (fd >= 0)
    {
      errno = save_errno;
      return fd;
    }
      else if (errno != EEXIST)
    return -1;
    }

  /* We got out of the loop because we ran out of combinations to try.  */
  errno = EEXIST;
  return -1;
}
#endif




/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Copies any <b>File</b> (text, binary, pdf, doc, png, ...) from DISK.
 *
 * @author Samuel Martins
 * @date December 01, 2015
 * @ingroup File
 *
 * @attention All exception checks were done on the function iftCopyFromDisk().
 */
void _iftCopyFileFromDisk(const char *src, const char *dst) {
    long BUFFER_SIZE = 1024; // for speed, let's copy a block of 1024 chars at a time, not just a char
    long n;
    char buffer[BUFFER_SIZE];
    char *dst_pathname = NULL;
    
    // copies src file inside the dst directory with the same src basename
    if (iftDirExists(dst)) {
        char *src_basename = iftFilename(src, NULL);
        dst_pathname       = iftJoinPathnames(2, dst, src_basename);
        iftFree(src_basename);
    }
    else {
        dst_pathname = iftCopyString(dst);
    }
    
    FILE *src_fp = fopen(src, "rb");
    FILE *dst_fp = fopen(dst_pathname, "wb");
    
    while ((n = fread(buffer, 1, BUFFER_SIZE, src_fp)) != 0) {
        fwrite(buffer, 1, n, dst_fp);
    }
    
    iftFree(dst_pathname);
    fclose(src_fp);
    fclose(dst_fp);
}







/********************** PUBLIC FUNCTIONS *************************/
/**
 * If the file follows the ift file format, updates the File properties. (id and label)
 */
iftFile *iftCreateFile(const char *format, ...) {
    va_list args;
    char pathname[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(pathname, format, args);
    va_end(args);
    
    // it is a Directory instead of a File
    if (iftDirExists(pathname))
        iftError("Pathname \"%s\" is a directory", "iftCreateFile", pathname);
    
    iftFile *f = (iftFile*) iftAlloc(1, sizeof(iftFile));
    f->path = iftCopyString(pathname);
    f->suffix = NULL;
    iftUpdateDataSetFileInfo(f);
    f->status = IFT_TRAIN;
    
    return f;
}


void iftDestroyFile(iftFile **f) {
    if (f != NULL) {
        iftFile *f_aux = *f;
        
        if (f_aux != NULL) {
            if (f_aux->path != NULL) {
                iftFree(f_aux->path);
                f_aux->path = NULL;
            }
            if(f_aux->suffix != NULL) {
                iftFree(f_aux->suffix);
                f_aux->suffix = NULL;
            }
            iftFree(f_aux);
            *f = NULL;
        }
    }
}


void iftPrintFileInfo(const iftFile *f) {
    if (f != NULL) {
        printf("Pathname: %s\n", f->path);
        printf("Stored on disk: ");
        if (iftFileExists(f->path))
            printf("TRUE\n");
        else
            printf("FALSE\n");
    }
}


/**
 * @brief Makes a Temp File on disk.
 *
 * Makes a Temp File on disk and returns your filename.
 *
 * @return The pathname from the temp file.
 */
char *iftMakeTempFile(const char *prefix, const char *suffix, const char *dir_name) {
    int id;
    char temp[12];
    
    strcpy(temp, "XXXXXX");

#if defined(__linux) || defined(__APPLE__)
    // Create Temporary Names
    if ((id = mkstemp(temp)) == -1) {
        iftError("MKSTEMP Error", "iftMakeTempFile");
    }
#endif

#if defined(__WIN32) || defined(__WIN64)
    if((id = mkstemp_custom(temp)) == -1) {
        iftError("MKSTEMP Error", "iftMakeTempFile");
    }
#endif
    close(id);
    
    char *final_temp = iftCopyString(temp);
    char *aux        = NULL;
    if (prefix != NULL) {
        aux = iftConcatStrings(2, prefix, final_temp);
        iftFree(final_temp);
        final_temp = aux;
    }
    
    if (suffix != NULL) {
        aux = iftConcatStrings(2, final_temp, suffix);
        iftFree(final_temp);
        final_temp = aux;
    }
    
    if (dir_name == NULL)
        aux = iftJoinPathnames(2, "/tmp", final_temp);
    else
        aux = iftJoinPathnames(2, dir_name, final_temp);
    
    iftFree(final_temp);
    final_temp = aux;
    
    // temporary file already exists... generate another one
    if (iftFileExists(final_temp)) {
        iftFree(final_temp);
        return iftMakeTempFile(prefix, suffix, dir_name);
    }
    else {
        if (!iftCompareStrings(temp, final_temp)) {
            iftRunProgram("mv", "%s %s", temp, final_temp);
        }
    }
    
    return final_temp;
}


char *iftMakeTempPathname(const char *prefix, const char *suffix, const char *dir_name) {
    char *tmp_pathname = iftMakeTempFile(prefix, suffix, dir_name);
    iftRemoveFile(tmp_pathname);
    
    return tmp_pathname;
}


iftFile *iftCopyFile(const iftFile* file) {
    iftFile* f = NULL;
    
    if (file == NULL)
        iftError("The file to be copied is NULL", "iftCopyFile");
    else {
        f         = (iftFile*) iftAlloc(1, sizeof(iftFile));
        f->path   = iftCopyString(file->path);
        f->sample = file->sample;
        f->label  = file->label;
        f->status = file->status;
        f->suffix = NULL;
        if(file->suffix != NULL)
            f->suffix   = iftCopyString(file->suffix);
    }
    
    return f;
}


void iftCopyFromDisk(const char *src, const char *dst) {
    if (src == NULL)
        iftError("Source Filename is NULL", "iftCopyFromDisk");
    if (dst == NULL)
        iftError("Target Filename is NULL", "iftCopyFromDisk");
    if (!iftFileExists(src) && !iftDirExists(src))
        iftError("Source filename \"%s\" does not exist", "iftCopyFromDisk", src);
    if (iftDirExists(src) && iftFileExists(dst))
        iftError("Cannot overwrite a non-directory \"%s\" with directory \"%s\"", "iftCopyFromDisk", dst, src);
    
    
    if (iftFileExists(src))
        _iftCopyFileFromDisk(src, dst);
    else
        iftCopyDirFromDisk(src, dst);
}


void iftRemoveFile(const char *pathname) {
    if (pathname != NULL) {
        if (iftPathnameExists(pathname)) {
            if (iftDirExists(pathname))
                iftWarning("The pathname \"%s\" is a DIRECTORY... Nothing to do", "iftRemoveFile", pathname);
            else if (remove(pathname) != 0)
                iftError("Problem to remove the file: \"%s\"", "iftRemoveFile", pathname);
        }
    }
}

void iftRenameOrMoveFile(const char *oldfilename, const char *newfilename)
{
    if (oldfilename != NULL && newfilename != NULL) {
        if (iftPathnameExists(oldfilename)) {
            if (iftDirExists(oldfilename))
                iftWarning("The pathname \"%s\" is a DIRECTORY... Nothing to do", "iftRenameOrMoveFile", oldfilename);
            else if (rename(oldfilename, newfilename) != 0)
                iftError("Problem to rename/move the file: \"%s\"", "iftRenameOrMoveFile", oldfilename);
        }
    }
}

bool iftFileExists(const char *pathname) {
    return (iftPathnameExists(pathname) && !iftDirExists(pathname));
}

char *iftFilename(const char *pathname, const char *suffix) {
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftFilename");
    
    char *base = iftSplitStringAt(pathname, IFT_SEP_C, -1);
    
    if ((suffix != NULL) && (!iftCompareStrings(suffix, ""))) {
        char *out_base = iftRemoveSuffix(base, suffix);
        iftFree(base);
        base = out_base;
    }
    
    return base;
}

char *iftDirname(const char *pathname) {
    char *dirname = NULL;
    
    if(pathname[strlen(pathname)-1] == '/') {
        dirname = iftCopyString(pathname);
    }
    else {
        iftSList *SL = iftSplitString(pathname, "/");
        dirname = iftCopyString("");
        char *aux = NULL;
        iftSNode *ptr = SL->head;
        for(int i = 0; i < SL->n-1; i++) {
            aux = iftConcatStrings(3, dirname, ptr->elem, "/");
            iftFree(dirname);
            dirname = iftCopyString(aux);
            iftFree(aux);
            ptr = ptr->next;
        }
        iftDestroySList(&SL);
    }
    iftRightTrim(dirname, '/');
    
    return dirname;
}

const char *iftFileExt(const char *pathname) {
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftFileExt");
    
    const char *dot = strrchr(pathname, '.'); // returns a pointer to the last occurrence of '.'
    
    if ((!dot) || (dot == pathname)) {
        return ("");
    }
    else {
        if (iftRegexMatch(pathname, "^.*\\.tar\\.(gz|bz|bz2)$") || iftRegexMatch(pathname, "^.*\\.(scn|nii)\\.gz$")) {
            dot -= 4; // points to the penultimate dot '.'
        }
        
        return dot; // returns the extension with '.'
    }
}


bool iftAreFileExtensionsEqual(const char *file1,const char *file2) {
    if ((file1 == NULL) || (file2 == NULL))
        iftError("One or both files are NULL", "iftAreFileExtensionsEqual");
    
    if (iftCompareStrings(iftFileExt(file1), iftFileExt(file2)))
        return true;
    else
        return false;
}


char *iftJoinPathnames(long n, ...) {
    if (n <= 0)
        iftError("Number of pathnames to be concatenated is <= 0", "iftJoinPathnames");
    
    long out_str_size = 1; // '\0'
    
    // Counts the size of the concatenated string
    va_list path_list;
    va_start(path_list, n);
    for (int i = 0; i < n; i++)
        out_str_size += strlen(va_arg(path_list, char*)) + 1; // one char for '/' (separation char)
    va_end(path_list);
    
    char *joined_path = iftAllocCharArray(out_str_size);
    char *aux = iftAllocCharArray(out_str_size);
    
    va_start(path_list, n);
    strcpy(joined_path, va_arg(path_list, char*));
    
    for (int i = 1; i < n; i++) {
        char *path = va_arg(path_list, char*);
        if (iftStartsWith(path, IFT_SEP_C))
            path++; // skip the first char, which is the directory separator
        
        if (iftEndsWith(joined_path, IFT_SEP_C))
            sprintf(aux, "%s%s", joined_path, path);
        else
            sprintf(aux, "%s%s%s", joined_path, IFT_SEP_C, path);

        iftFree(joined_path);
        joined_path = iftCopyString(aux);
    }
    iftFree(aux);
    
    return joined_path;
}


char *iftAbsPathname(const char *pathname) {
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftAbsPathname");
    
    char *abs_pathname = iftAllocCharArray(IFT_STR_DEFAULT_SIZE);
    char *ptr = NULL;
    
    char *final_path = NULL;
    
    if (strlen(pathname) == 1) {
        if (pathname[0] == '~') {
            final_path = getenv("HOME");
        }
    }
    else {
        if ((pathname[0] == '~') && (pathname[1] == IFT_SEP_C[0])) {
            const char *pt = pathname + 2;
            
            char *home = getenv("HOME");
            final_path = iftJoinPathnames(2, home, pt);
        }
    }
    
    if (final_path == NULL)
        final_path = iftCopyString(pathname);


#if defined(__linux) || (__APPLE__)
    ptr = realpath(final_path, abs_pathname);
    
    if (ptr != abs_pathname)
        iftError("Problems with the function realpath\nPath: %s", "iftAbsPathname", final_path);
#else
    if (GetFullPathName(final_path, IFT_STR_DEFAULT_SIZE, abs_pathname, NULL) == 0)
        iftError("Problems with the function realpath\nPath: %s", "iftAbsPathname", final_path);
#endif
    
    iftFree(final_path);
    
    return abs_pathname;
}


char *iftBasename(const char *pathname) {
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftBasename");
    
    char *base = NULL;
    
    // If the pathname does not have a file extension, we use it entirely as the basename
    if(iftCompareStrings(iftFileExt(pathname), "")) {
        base = iftCopyString(pathname);
    } else {
        base = iftSplitStringAt(pathname, iftFileExt(pathname), 0);
    }
    return base;
}


bool iftIsImagePathnameValid(const char *img_pathname) {
    if (img_pathname == NULL)
        iftError("Image Pathname is NULL", "iftIsImagePathnameValid");
    
    char *lower_pathname = iftLowerString(img_pathname);
    bool is_valid        = iftRegexMatch(lower_pathname, "^.+\\.(jpg|jpeg|pgm|ppm|scn|png|scn\\.gz|zscn|hdr|nii|nii\\.gz)$");
    iftFree(lower_pathname);
    
    return is_valid;
}


bool iftIsImageFile(const char *img_pathname) {
    if (img_pathname == NULL)
        iftError("Image Pathname is NULL", "iftIsImageFile");
    
    return (iftFileExists(img_pathname) && iftIsImagePathnameValid(img_pathname));
}


char *iftAddEscapeChars(const char *pathname) {
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftAddEscapeChars");
    
    char *out = iftAllocCharArray(strlen(pathname)+1000); // max num of escapes: 1000
    long j = 0; // index of the <out>
    
    for (long i = 0; i < strlen(pathname); i++) {
        if ((pathname[i] == IFT_SEP_C[0]) || (pathname[i] == '(') || (pathname[i] == ')'))
            out[j++] = '\\';
        
        out[j++] = pathname[i];
    }
    out[j] = '\0';
    
    return out;
}


char *iftExpandUser(const char *path) {
    return iftReplaceString(path,"~", getenv("HOME"));
}


int iftImageSampleId(char *path) {
    int sample = IFT_NIL;
    
    if(iftValidDataSetFilenameConvention(path))
    {
        const char *FILE_SEP = "_";
        char *bname = iftFilename(path, iftFileExt(path));
        char *str = iftSplitStringAt(bname, FILE_SEP, 1);
        
        sample = atoi(str);
        
        iftFree(str);
        iftFree(bname);
    } else {
        iftError("String %s not following the IFT file format.", "iftImageSampleLabel", path);
    }
    
    return sample;
}


int iftImageSampleLabel(char *path)
{
    int label = IFT_NIL;
    
    if(iftValidDataSetFilenameConvention(path))
    {
        const char *FILE_SEP = "_";
        char *bname = iftFilename(path, iftFileExt(path));
        char *str = iftSplitStringAt(bname, FILE_SEP, 0);
        
        label = atoi(str);
        
        iftFree(str);
        iftFree(bname);
    } else {
        iftError("String %s not following the IFT file format.", "iftImageSampleLabel", path);
    }
    
    return label;
}


int* iftGetImageSampleInfo(char *path){
    int* sampleInfo = NULL;
    if(iftValidDataSetFilenameConvention(path)){
        sampleInfo = calloc(2,sizeof(int));
        const char *FILE_SEP = "_";
        char *bname = iftFilename(path, iftFileExt(path));
        char *label_str = iftSplitStringAt(bname, FILE_SEP, 0);
        char *id_str = iftSplitStringAt(bname, FILE_SEP, 1);
        
        sampleInfo[0] = atoi(label_str);
        sampleInfo[1] = atoi(id_str);
        
        iftFree(label_str);
        iftFree(id_str);
        iftFree(bname);
    }else{
        iftError("String %s not following the IFT file format.", "iftGetImageSampleInfo", path);
    }
    return sampleInfo;
}


char* iftImageSuffix(char* path)
{
    char *suffix       = NULL;
    
    if(iftValidDataSetFilenameConvention(path)) {
        const char *FILE_SEP = "_";
        
        char *bname = iftFilename(path, iftFileExt(path));
        iftSList *SL = iftSplitString(bname, FILE_SEP);
        
        if (SL->n >= 3) {
            long i, suffix_len = 0;
            iftSNode *aux = NULL;
            
            for (aux = SL->head->next->next; aux != NULL; aux = aux->next) {
                suffix_len += strlen(aux->elem);
            }
            
            suffix_len += SL->n - 3; // Considering all of the separators '_' that will be added
            
            suffix = iftAllocCharArray(suffix_len + 1); // Considering '\0'
            
            i = 0;
            for (aux = SL->head->next->next; aux != NULL; aux = aux->next) {
                strcpy(suffix + i, aux->elem);
                i += strlen(aux->elem);
                
                // Appending the separator '_'
                if (aux->next != NULL) {
                    strcpy(suffix + i, FILE_SEP);
                    i++;
                }
            }
        } else {
            suffix = iftCopyString("");
        }
        
        iftFree(bname);
        iftDestroySList(&SL);
    } else {
        iftError("String %s not following the IFT file format.", "iftImageSuffix", path);
    }
    
    return suffix;
}


bool iftValidDataSetFilenameConvention(const char *path) {
    if (path == NULL)
        iftError("Pathname is NULL", "iftValidDataSetFilenameConvention");
    
    char *filename = iftFilename(path, NULL);
    bool valid_dataset_file = iftRegexMatch(filename, "^[0-9]+_[0-9]+.+$");
    iftFree(filename);
    
    return valid_dataset_file;
}


void iftUpdateDataSetFileInfo(iftFile* file) {
    if (iftValidDataSetFilenameConvention(file->path)) {
        file->sample = iftImageSampleId(file->path);
        file->label  = iftImageSampleLabel(file->path);
        file->suffix = iftImageSuffix(file->path);
    }
}

