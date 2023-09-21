#include "iftCompression.h"

#include "ift/core/io/Stream.h"

/**
 * @file iftCompression.c
 * @brief Low level i/o interface to compressed and noncompressed files.
 *        Written by Mark Jenkinson, FMRIB
 *
 * This library provides an interface to both compressed (gzip/zlib) and
 * uncompressed (normal) file IO.  The functions are written to have the
 * same interface as the standard file IO functions.
 *
 * To use this library instead of normal file IO, the following changes
 * are required:
 *  - replace all instances of FILE* with iftGZipFile
 *  - change the name of all function calls, replacing the initial character
 *    f with the znz  (e.g. fseek becomes iftGZipSeek)
 *    one exception is rewind() -> iftGZipRewind()
 *  - add a third parameter to all calls to znzopen (previously fopen)
 *    that specifies whether to use compression (1) or not (0)
 *  - use znz_isnull rather than any (pointer == NULL) comparisons in the code
 *    for znzfile types (normally done after a return from iftGZipOpen)
 *
 * NB: seeks for writable files with compression are quite restricted
 *
 */

/*
 * iftCompression.c  (zipped or non-zipped library)
 * *****            This code is released to the public domain.            *****

 * *****  Author: Mark Jenkinson, FMRIB Centre, University of Oxford       *****
 * *****  Date:   September 2004                                           *****

 * *****  Neither the FMRIB Centre, the University of Oxford, nor any of   *****
 * *****  its employees imply any warranty of usefulness of this software  *****
 * *****  for any purpose, and do not assume any liability for damages,    *****
 * *****  incidental or otherwise, caused by any use of this document.     *****
 */



iftGZipFile iftGZipOpen(const char *path, const char *mode, bool use_compression)
{
    iftGZipFile file;
    file = (iftGZipFile) iftAlloc(1, sizeof(struct iftGZipPtr));
    if( file == NULL ){
        fprintf(stderr,"** IFT_ERROR: iftGZipOpen failed to alloc znzptr\n");
        return NULL;
    }

    file->nzfptr = NULL;

#ifdef HAVE_ZLIB
    file->zfptr = NULL;

    if (use_compression) {
        file->withz = 1;
        if((file->zfptr = gzopen(path,mode)) == NULL) {
            iftFree(file);
            file = NULL;
        }
    } else {
#endif

        file->withz = 0;
        if((file->nzfptr = fopen(path,mode)) == NULL) {
            iftFree(file);
            file = NULL;
        }

#ifdef HAVE_ZLIB
    }
#endif

    return file;
}


iftGZipFile iftGZipDOpen(int fd, const char *mode, bool use_compression)
{
    iftGZipFile file;
    file = (iftGZipFile) iftAlloc(1, sizeof(struct iftGZipPtr));
    if( file == NULL ){
        fprintf(stderr,"** IFT_ERROR: iftGZipDOpen failed to alloc znzptr\n");
        return NULL;
    }
#ifdef HAVE_ZLIB
    if (use_compression) {
        file->withz = 1;
        file->zfptr = gzdopen(fd,mode);
        file->nzfptr = NULL;
    } else {
#endif
        file->withz = 0;
#ifdef HAVE_FDOPEN
        file->nzfptr = fdopen(fd,mode);
#endif
#ifdef HAVE_ZLIB
        file->zfptr = NULL;
    };
#endif
    return file;
}


int iftGZipClose(iftGZipFile *file)
{
    int retval = 0;
    if (*file!=NULL) {
#ifdef HAVE_ZLIB
        if ((*file)->zfptr!=NULL)  { retval = gzclose((*file)->zfptr); }
#endif
        if ((*file)->nzfptr!=NULL) { retval = fclose((*file)->nzfptr); }

        iftFree(*file);
        *file = NULL;
    }
    return retval;
}


/* we already assume ints are 4 bytes */
#undef GZip_MAX_BLOCK_SIZE
#define GZip_MAX_BLOCK_SIZE (1<<30)

size_t iftGZipRead(void *buf, size_t size, size_t nmemb, iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    size_t     remain = size*nmemb;
    char     * cbuf = (char *)buf;
    unsigned   n2read;
    int        nread;

    if (file->zfptr!=NULL) {
        /* gzread/write take unsigned int length, so maybe read in int pieces
           (noted by M Hanke, example given by M Adler)   6 July 2010 [rickr] */
        while( remain > 0 ) {
            n2read = (remain < GZip_MAX_BLOCK_SIZE) ? remain : GZip_MAX_BLOCK_SIZE;
            nread = gzread(file->zfptr, (void *)cbuf, n2read);
            if( nread < 0 ) return nread; /* returns -1 on error */

            remain -= nread;
            cbuf += nread;

            /* require reading n2read bytes, so we don't get stuck */
            if( nread < (int)n2read ) break;  /* return will be short */
        }

        /* warn of a short read that will seem complete */
        if( remain > 0 && remain < size )
            fprintf(stderr,"** iftGZipRead: read short by %u bytes\n",(unsigned)remain);

        return nmemb - remain/size;   /* return number of members processed */
    }
#endif
    return fread(buf,size,nmemb,file->nzfptr);
}

size_t iftGZipWrite(const void *buf, size_t size, size_t nmemb, iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    size_t     remain = size*nmemb;
    char     * cbuf = (char *)buf;
    unsigned   n2write;
    int        nwritten;

    if (file->zfptr!=NULL) {
        while( remain > 0 ) {
            n2write = (remain < GZip_MAX_BLOCK_SIZE) ? remain : GZip_MAX_BLOCK_SIZE;
            nwritten = gzwrite(file->zfptr, (void *)cbuf, n2write);

            /* gzread returns 0 on error, but in case that ever changes... */
            if( nwritten < 0 ) return nwritten;

            remain -= nwritten;
            cbuf += nwritten;

            /* require writing n2write bytes, so we don't get stuck */
            if( nwritten < (int)n2write ) break;
        }

        /* warn of a short write that will seem complete */
        if( remain > 0 && remain < size )
            fprintf(stderr,"** iftGZipWrite: write short by %u bytes\n",(unsigned)remain);

        return nmemb - remain/size;   /* return number of members processed */
    }
#endif
    return fwrite(buf,size,nmemb,file->nzfptr);
}

long iftGZipSeek(iftGZipFile file, long offset, int whence)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return (long) gzseek(file->zfptr,offset,whence);
#endif
    return fseek(file->nzfptr,offset,whence);
}

int iftGZipRewind(iftGZipFile stream)
{
    if (stream==NULL) { return 0; }
#ifdef HAVE_ZLIB
    /* On some systems, gzrewind() fails for uncompressed files.
     Use gzseek(), instead.               10, May 2005 [rickr]

     if (stream->zfptr!=NULL) return gzrewind(stream->zfptr);
  */

    if (stream->zfptr!=NULL) return (int)gzseek(stream->zfptr, 0L, SEEK_SET);
#endif
    rewind(stream->nzfptr);
    return 0;
}

long iftGZipTell(iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return (long) gztell(file->zfptr);
#endif
    return ftell(file->nzfptr);
}

int iftGZipPuts(const char *str, iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return gzputs(file->zfptr,str);
#endif
    return fputs(str,file->nzfptr);
}


char *iftGZipGets(char *str, int size, iftGZipFile file)
{
    if (file==NULL) { return NULL; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return gzgets(file->zfptr,str,size);
#endif
    return fgets(str,size,file->nzfptr);
}


int iftGZipFlush(iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return gzflush(file->zfptr,Z_SYNC_FLUSH);
#endif
    return fflush(file->nzfptr);
}


int iftGZipEOF(iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return gzeof(file->zfptr);
#endif
    return feof(file->nzfptr);
}


int iftGZipPutc(int c, iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return gzputc(file->zfptr,c);
#endif
    return fputc(c,file->nzfptr);
}


int iftGZipGetc(iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return gzgetc(file->zfptr);
#endif
    return fgetc(file->nzfptr);
}

#if !defined (WIN32)
int iftGZipPrintf(iftGZipFile stream, const char *format, ...)
{
    int retval=0;
    va_list va;
    if (stream==NULL) { return 0; }
    va_start(va, format);
#ifdef HAVE_ZLIB
    char *tmpstr;
    if (stream->zfptr!=NULL) {
        int size;  /* local to HAVE_ZLIB block */
        size = strlen(format) + 1000000;  /* overkill I hope */
        tmpstr = (char *)iftAlloc(1, size);
        if( tmpstr == NULL ){
            fprintf(stderr,"** IFT_ERROR: iftGZipPrintf failed to alloc %d bytes\n", size);
            return retval;
        }
        vsprintf(tmpstr,format,va);
        retval=gzprintf(stream->zfptr,"%s",tmpstr);
        iftFree(tmpstr);
    } else
#endif
    {
        retval=vfprintf(stream->nzfptr,format,va);
    }
    va_end(va);
    return retval;
}

#endif

//void iftCompressDirContents(const char *input, const char *output) {
//    iftFileSet *files = NULL;
//
//    if(!iftCompareStrings(iftFileExt(output), ".gz") && !iftCompareStrings(iftFileExt(output), ".GZ")){
//        iftError("Output compressed file extension must be .gz or .GZ", "iftCompressPath");
//    }
//
//    if(iftDirExists(input)) {
//        files = iftLoadFileSetFromDir(input, 0);
//    } else if(iftFileExists(input)) {
//        files = iftCreateFileSet(1);
//        files->files[0] = iftCreateFile(input);
//    } else {
//        iftError("A valid path to a directory or path must be passed as input. Path \"%s\" does not exist!",
//                 "iftCompressPath", input);
//    }
//
//    for(size_t i = 0; i < files->n; i++) {
//
//    }
//
//}
