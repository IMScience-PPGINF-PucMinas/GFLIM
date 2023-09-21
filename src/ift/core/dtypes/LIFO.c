#include "ift/core/dtypes/LIFO.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"
#include "iftMemory.h"


iftLIFO *iftCreateLIFO(int n)
{
    iftLIFO *L=(iftLIFO *)iftAlloc(1,sizeof(iftLIFO));
    
    L->LIFO   = iftAllocIntArray(n);
    L->color  = iftAllocCharArray(n);
    L->n      = n;
    L->last   = IFT_NIL;
    
    return(L);
}

void     iftDestroyLIFO(iftLIFO **L)
{
    iftLIFO *aux=*L;
    
    if (aux != NULL) {
        iftFree(aux->LIFO);
        iftFree(aux->color);
        iftFree(aux);
        *L = NULL;
    }
}

char iftInsertLIFO(iftLIFO *L, int node)
{
    if (iftFullLIFO(L)){
        iftWarning("LIFO is full","iftInsertLIFO");
        return 0;
    }
    L->last++; L->LIFO[L->last]=node;
    L->color[node]= IFT_GRAY;
    
    return 1;
}

int      iftRemoveLIFO(iftLIFO *L)
{
    int node= IFT_NIL;
    
    if (!iftEmptyLIFO(L)){
        node = L->LIFO[L->last];  L->last--;
        L->color[node]= IFT_BLACK;
    }else{
        iftWarning("LIFO is empty","iftRemoveLIFO");
    }
    
    return node;
}

char     iftEmptyLIFO(iftLIFO *L)
{
    if (L->last == IFT_NIL) {
        // Changed by Falcao and Nikolas
        // iftResetLIFO(L);
        return(1);
    }else
        return(0);
}

char     iftFullLIFO(iftLIFO *L)
{
    if (L->last==L->n-1) {
        return(1);
    }else
        return(0);
}

void iftResetLIFO(iftLIFO *L)
{
    int p;
    
    L->last= IFT_NIL;
    for (p=0; p < L->n; p++)
        L->color[p]= IFT_WHITE;
}

