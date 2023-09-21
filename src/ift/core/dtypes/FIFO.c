#include "ift/core/dtypes/FIFO.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"


iftFIFO *iftCreateFIFO(int n)
{
    iftFIFO *F=(iftFIFO *)iftAlloc(1,sizeof(iftFIFO));
    
    F->FIFO  = iftAllocIntArray(n);
    F->color = iftAllocCharArray(n);
    F->n     = n;
    F->first=F->last=0;
    
    return(F);
}

void iftDestroyFIFO(iftFIFO **F) {
    if (F != NULL) {
        iftFIFO *aux=*F;
        
        if (aux != NULL) {
            iftFree(aux->FIFO);
            iftFree(aux->color);
            iftFree(aux);
            *F = NULL;
        }
    }
}

char iftInsertFIFO(iftFIFO *F, int elem)
{
    if (iftFullFIFO(F)){
        iftWarning("FIFO is full","iftInsertFIFO");
        return 0;
    }
    F->color[elem]= IFT_GRAY;
    F->FIFO[F->last]=elem;  F->last++;
    
    return 1;
}

int      iftRemoveFIFO(iftFIFO *F)
{
    int node= IFT_NIL;
    
    if (!iftEmptyFIFO(F)){
        node = F->FIFO[F->first];  F->first++;
        F->color[node]= IFT_BLACK;
    }else{
        iftWarning("FIFO is empty","iftRemoveFIFO");
    }
    
    return node;
}

bool iftEmptyFIFO(iftFIFO *F)
{
    if (F->first == F->last) {
        // Changed by Falcao and Nikolas
        // iftResetFIFO(F);
        return(1);
    }else
        return(0);
}

bool iftFullFIFO(iftFIFO *F)
{
    if (F->last==F->n) {
        return(1);
    }else
        return(0);
}

void     iftResetFIFO(iftFIFO *F)
{
    int p;
    for (p=0; p < F->n; p++)
        F->color[p] = IFT_WHITE;
    F->first=F->last=0;
}
int    iftColorFIFO(iftFIFO *F, int pos)
{
    return F->color[pos];
}

