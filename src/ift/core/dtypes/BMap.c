#include "ift/core/dtypes/BMap.h"

#include "ift/core/io/Stream.h"


iftBMap *iftCreateBMap(int n) {
    iftBMap *b;
    b= (iftBMap *) iftAlloc(1,sizeof(iftBMap));
    b->n        = n;
    b->nbytes   = n/8;
    if (n%8) b->nbytes++;
    b->val = (char *) iftAlloc(b->nbytes,sizeof(char));
    if (b->val==NULL){
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateBMap");
    }
    return b;
}


void iftDestroyBMap(iftBMap **bmap) {
    iftBMap *aux=*bmap;
    
    if (aux != NULL) {
        iftFree(aux->val);
        iftFree(aux);
        *bmap=NULL;
    }
    
}


void iftFillBMap(iftBMap *bmap, int value) {
    int p;
    for (p=0; p < bmap->nbytes; p++)
        bmap->val[p] = (value?0xff:0);
}


iftBMap *iftCopyBMap(const iftBMap *src) {
    iftBMap *dst=iftCreateBMap(src->n);
    int p;
    for (p=0; p < src->nbytes; p++)
        dst->val[p] = src->val[p];
    return(dst);
}


iftBMap *iftReadBMap(const char *path) {
    if (path == NULL)
        iftError("Pathname is NULL", "iftReadBMap");
    
    FILE *fp = fopen(path, "rb");
    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftReadBMap", path);
    
    int n;
    if (fscanf(fp, "%d\n", &n) != 1)
        iftError("Error when reading the Bit Map Size", "iftReadBMap");
    
    iftBMap *bmap = iftCreateBMap(n);
    if (fread(bmap->val, sizeof(char), bmap->nbytes, fp) != bmap->nbytes)
        iftError("Error when reading the Bit Map Values", "iftReadBMap");
    
    fclose(fp);
    
    return bmap;
}


void iftWriteBMap(const iftBMap *bmap, const char *path) {
    if (bmap == NULL)
        iftError("Bit Map is NULL", "iftWriteBMap");
    if (path == NULL)
        iftError("Pathname is NULL", "iftWriteBMap");
    
    FILE *fp = fopen(path, "wb");
    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteBMap", path);
    
    fprintf(fp, "%d\n", bmap->n);
    fwrite(bmap->val, sizeof(char), bmap->nbytes, fp);
    
    fclose(fp);
}
