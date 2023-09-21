#include "ift/core/io/Stream.h"

#include "ift/core/io/File.h"
#include "ift/core/io/Dir.h"
#include "ift/core/tools/Regex.h"
#include "ift/core/tools/String.h"
#include "iftMemory.h"


char *iftReadFileAsString(const char *format, ...) {
    va_list args;
    char pathname[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(pathname, format, args);
    va_end(args);
    
    FILE *fp = fopen(pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadFileAsString", pathname);
    
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    
    char *data = (char*) iftAlloc((len+1), sizeof(char));
    if (fread(data, 1, len, fp) != len)
        iftError("Number of read chars is different from the allocated string size", "iftReadFileAsString");
    
    data[len] = '\0';
    
    fclose(fp);
    
    return data;
}


void iftWriteDataInFile(FILE *file, const void *data, size_t n, size_t chuck_size) {
    if (file == NULL)
        iftError("File is NULL", "iftWriteDataInFile");
 
    if (fwrite(data, chuck_size, n, file) != n)
        iftError("Cannot save all data elements", "iftWriteDataInFile");
}


void iftWriteStringToFile(const char *str, const char *format, ...) {
    va_list args;
    char pathname[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(pathname, format, args);
    va_end(args);
    
    FILE *fp = fopen(pathname, "w");
    fprintf(fp, "%s", str);
    fclose(fp);
}


void iftZipDirContent(const char *dir_path, const char *out_zip_path) {
    if (!iftDirExists(dir_path))
        iftError("Dir \"%s\" does not exist", "iftZipDirContent", dir_path);
    /*
    if (!iftEndsWith(out_zip_path, ".zip"))
        iftError("Invalid Output Zip file \"%s\"", "iftZipDirContent", out_zip_path);
    */
    
    // ZIP ALL FILES
    char cmd[2048];
    char *tmp_path = iftMakeTempPathname("tmp_", ".zip", NULL);
    
    sprintf(cmd, ">/dev/null cd %s; zip -rq %s *;>/dev/null cd -; cp %s %s; rm %s", dir_path, tmp_path, tmp_path, out_zip_path, tmp_path);
    if (system(cmd) == -1) {
        iftError("System Cmd Error when zipping the content of the directory \"%s\"", "iftZipDirContent", dir_path);
    }
    iftFree(tmp_path);
}


void iftUnzipFile(const char *zip_path, const char *out_dir) {
    // checkers
    if (zip_path == NULL)
        iftError("Zip File is NULL", "iftUnzipFile");
    /*
    if (!iftRegexMatch(zip_path, "^.+\\.zip$"))
        iftError("Invalid Zip File \"%s\"", "iftUnzipFile", zip_path);
    */

    if (out_dir == NULL)
        iftError("Output Directory is NULL", "iftUnzipFile");
    if (!iftDirExists(out_dir)) {
        iftMakeDir(out_dir);
    }
    
    char cmd[2048];
    sprintf(cmd, "unzip -q %s -d %s", zip_path, out_dir);
    if (system(cmd) == -1) {
        iftError("System Cmd Error when unzipping the Zip File  \"%s\"", "iftReadBlockClassifiers", zip_path);
    }
}



bool iftFindInFile(const char *pathname, const char *str) {
    if (pathname == NULL)
        iftError("Pathname from the File is NULL", "iftSearchInFile");
    if (!iftFileExists(pathname))
        iftError("Pathname is not a File or does not exists", "iftSearchInFile");
    if (str == NULL)
        iftError("String is NULL", "iftSearchInFile");
    
    FILE *fp        = fopen(pathname, "rb");
    char *tmp       = iftAllocCharArray(IFT_STR_DEFAULT_SIZE);
    bool find_result = false;
    
    while (fgets(tmp, IFT_STR_DEFAULT_SIZE, fp) != NULL)
        if ((strstr(tmp, str)) != NULL) {
            find_result = true;
            break;
        }
    
    fclose(fp);
    iftFree(tmp);
    
    return find_result;
}


bool iftIsFileContentEmpty(FILE *fp) {
    long saved_off_set = ftell(fp);
    fseek(fp, 0, SEEK_END);
    
    if (ftell(fp) == 0)
        return true;
    
    // it is not empty - move the pointer to the original position
    fseek(fp, saved_off_set, SEEK_SET);
    
    return false;
}


char *iftGetLine(FILE *stream) {
    if (stream == NULL || feof(stream))
        return NULL;
    
    size_t buffer_size = 0;
    char *line = NULL;
    
    // getline reads an entire line from stream, storing the text (including the newline and a terminating null character)
    // in a buffer and storing the buffer address in *line
    // If you set *line to a null pointer, and *buffer_size to zero, before the call, then getline
    // allocates the initial buffer for you by calling malloc. This buffer should be freed by the
    // user program even if getline() failed.
    // If an error occurs or end of file is reached without any bytes read, getline returns -1.
    if (getline(&line, &buffer_size, stream) == -1) {
        if (feof(stream)) {
        	if (line != NULL)
        		free(line);
        	return NULL;
        } else
            iftError("Error with getline command", "iftGetLine");
    }
    
    size_t str_len = strlen(line);
    if (line[str_len-1] == '\n')
        line[str_len-1] = '\0'; // without this, it reads the \n and does not put the \0 at the end
    
    return line;
}


void iftSed(const char *file_path, const char *old_str, const char *new_str) {
    char *tmpfile = iftMakeTempFile("tmpfile_", NULL, NULL);
    
    FILE *fp = fopen(file_path, "r");
    FILE *fw = fopen(tmpfile, "w");
    
    
    char *line = iftGetLine(fp);
    while (line != NULL) {
        char *replaced_line = iftReplaceString(line, old_str, new_str);
        fprintf(fw, "%s\n", replaced_line);
        
        iftFree(line);
        iftFree(replaced_line);
        
        line = iftGetLine(fp);
    }
    
    fclose(fp);
    fclose(fw);
    iftCopyFromDisk(tmpfile, file_path);
    iftRemoveFile(tmpfile);
    iftFree(tmpfile);
}


size_t iftCDataTypeSize(iftCDataType cdtype) {
    switch (cdtype) {
        case IFT_CHAR_TYPE:
            return sizeof(char);
        case IFT_SHORT_TYPE:
            return sizeof(short);
        case IFT_INT_TYPE:
            return sizeof(int);
        case IFT_LONG_TYPE:
            return sizeof(long);
        case IFT_UCHAR_TYPE:
            return sizeof(uchar);
        case IFT_USHORT_TYPE:
            return sizeof(ushort);
        case IFT_UINT_TYPE:
            return sizeof(uint);
        case IFT_ULONG_TYPE:
            return sizeof(ulong);
        case IFT_FLT_TYPE:
            return sizeof(float);
        case IFT_DBL_TYPE:
            return sizeof(double);
        case IFT_PTR_TYPE:
            return sizeof(void *);
        default:
            iftError("Unsupported CDataType", "iftCDataTypeSize");
            return 0;
    }
}


