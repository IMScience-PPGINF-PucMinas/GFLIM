#include "ift/core/dtypes/ULongArray.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"


iftULongArray *iftCreateULongArray(long n) {
    iftULongArray *iarr = (iftULongArray*) iftAlloc(1, sizeof(iftULongArray));
    
    iarr->n = n;
    iarr->val  = iftAllocULongArray(n);
    
    return iarr;
}


void iftDestroyULongArray(iftULongArray **iarr) {
    if (iarr != NULL && *iarr != NULL) {
        iftULongArray *iarr_aux = *iarr;
        
        if (iarr_aux->val != NULL)
            iftFree(iarr_aux->val);
        iftFree(iarr_aux);
        *iarr = NULL;
    }
}

void iftResizeULongArray(iftULongArray **iarr, long n) {
    if (n == 0)
        iftError("The new size must be greater than 0", "iftResizeULongArray");
    
    if ((*iarr) == NULL) {
        (*iarr) = iftCreateULongArray(n);
    } else {
        if(n != (*iarr)->n) {
            iftULongArray *iarr1 = iftCreateULongArray((*iarr)->n);
            iftCopyData(iarr1->val, (*iarr)->val, (*iarr)->n, sizeof(ulong));
            iftDestroyULongArray(iarr);
            (*iarr) = iftCreateULongArray(n);
            iftCopyData((*iarr)->val, iarr1->val, iftMin(iarr1->n, (*iarr)->n), sizeof(ulong));
            iftDestroyULongArray(&iarr1);
        }
    }
}
