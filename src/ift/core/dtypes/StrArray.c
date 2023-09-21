#include "ift/core/dtypes/StrArray.h"

#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"

#include "iftCommon.h"
#include "iftMemory.h"



iftStrArray *iftCreateStrArray(long n) {
    iftStrArray *sarr = (iftStrArray *) iftAlloc(1, sizeof(iftStrArray));
    
    sarr->n   = n;
    sarr->val = (char**) iftAlloc(n, sizeof(char*));
    for (int i = 0; i < n; i++) {
        sarr->val[i] = iftAllocCharArray(2048);
    }
    
    return sarr;
}


void iftDestroyStrArray(iftStrArray **sarr) {
    if (sarr != NULL && *sarr != NULL) {
        iftStrArray *sarr_aux = *sarr;
        
        if (sarr_aux->val != NULL) {
            for (int i = 0; i < sarr_aux->n; i++) {
                iftFree(sarr_aux->val[i]);
            }
            iftFree(sarr_aux->val);
        }
        iftFree(sarr_aux);
        *sarr = NULL;
    }
}

void iftResizeStrArray(iftStrArray **iarr, long n) {
    if (n == 0)
        iftError("The new size must be greater than 0", "iftResizeStrArray");
    
    if ((*iarr) == NULL) {
        (*iarr) = iftCreateStrArray(n);
    } else {
        if(n != (*iarr)->n) {
            iftStrArray *iarr1 = iftCopyStrArray((*iarr)->val, (*iarr)->n);
            iftDestroyStrArray(iarr);
            (*iarr) = iftCreateStrArray(n);
            for (int i = 0; i < iftMin(iarr1->n, (*iarr)->n); i++) {
                iftFree((*iarr)->val[i]);
                (*iarr)->val[i] = iftCopyString(iarr1->val[i]);
            }
            iftDestroyStrArray(&iarr1);
        }
    }
}

iftStrArray *iftCopyStrArray(char **src_arr, long n) {
    iftStrArray *copy = iftCreateStrArray(n);
    
    for (int i = 0; i < copy->n; i++) {
        iftFree(copy->val[i]);
        copy->val[i] = iftCopyString(src_arr[i]);
    }
    
    return copy;
}

char *iftStrArrayAsString(const iftStrArray *sarr) {
    // 2 chars for square brackets + 2050 chars per number (at most) + n-1 commas
    char *str = iftAlloc(2 + (2050 * sarr->n) + (sarr->n-1), sizeof(char));
    
    sprintf(str, "[");
    for (long i = 0; i <= (sarr->n-2); i++)
        sprintf(str, "%s\"%s\", ", str, sarr->val[i]);
    sprintf(str, "%s\"%s\"]", str, sarr->val[sarr->n-1]);
    
    return str;
}


iftStrArray *iftStrValues(int n, ...) {
    int i;
    const char* str;
    iftStrArray * array = iftCreateStrArray(n);
    
    va_list arguments;
    va_start ( arguments, n );
    
    for (i = 0; i < n; ++i) {
        str = va_arg ( arguments, const char* );
        strcpy(array->val[i], str);
    }
    
    va_end ( arguments );                  // Cleans up the list
    
    return array;
}
