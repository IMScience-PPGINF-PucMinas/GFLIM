#include "ift/core/dtypes/FHeap.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"


iftFHeap *iftCreateFHeap(int n, float *value) {
    iftFHeap *H = NULL;
    int i;
    
    if (value == NULL) {
        iftError("Cannot create heap without priority value map", "iftCreateFHeap");
    }
    
    H = (iftFHeap *) iftAlloc(1, sizeof(iftFHeap));
    if (H != NULL) {
        H->n       = n;
        H->value    = value;
        H->color   = (char *) iftAlloc(n, sizeof(char));
        H->node   = (int *) iftAlloc(n, sizeof(int));
        H->pos     = (int *) iftAlloc(n, sizeof(int));
        H->last    = -1;
        H->removal_policy = IFT_MINVALUE;
        if (H->color == NULL || H->pos == NULL || H->node == NULL)
            iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateFHeap");
        for (i = 0; i < H->n; i++) {
            H->color[i] = IFT_WHITE;
            H->pos[i]   = -1;
            H->node[i] = -1;
        }
    }
    else
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateFHeap");
    
    return H;
}

void iftDestroyFHeap(iftFHeap **H) {
    iftFHeap *aux = *H;
    if (aux != NULL) {
        if (aux->node != NULL) iftFree(aux->node);
        if (aux->color != NULL) iftFree(aux->color);
        if (aux->pos != NULL)   iftFree(aux->pos);
        iftFree(aux);
        *H = NULL;
    }
}

void  iftGoUpFHeap(iftFHeap *H, int i) {
    int j = iftDad(i);
    
    if(H->removal_policy == IFT_MINVALUE){
        
        while ((j >= 0) && (H->value[H->node[j]] > H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
    else{ /* removal_policy == IFT_MAXVALUE */
        
        while ((j >= 0) && (H->value[H->node[j]] < H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
}

void iftGoDownFHeap(iftFHeap *H, int i) {
    int j, left = iftLeftSon(i), right = iftRightSon(i);
    
    j = i;
    if(H->removal_policy == IFT_MINVALUE){
        
        if ((left <= H->last) &&
            (H->value[H->node[left]] < H->value[H->node[i]]))
            j = left;
        if ((right <= H->last) &&
            (H->value[H->node[right]] < H->value[H->node[j]]))
            j = right;
    }
    else{ /* removal_policy == IFT_MAXVALUE */
        
        if ((left <= H->last) &&
            (H->value[H->node[left]] > H->value[H->node[i]]))
            j = left;
        if ((right <= H->last) &&
            (H->value[H->node[right]] > H->value[H->node[j]]))
            j = right;
    }
    
    if(j != i) {
        iftSwap(H->node[j], H->node[i]);
        H->pos[H->node[i]] = i;
        H->pos[H->node[j]] = j;
        iftGoDownFHeap(H, j);
    }
}

char iftFullFHeap(iftFHeap *H) {
    if (H->last == (H->n - 1))
        return 1;
    else
        return 0;
}

char iftEmptyFHeap(iftFHeap *H) {
    if (H->last == -1){
        return 1;
    }else{
        return 0;
    }
}


char iftInsertFHeap(iftFHeap *H, int node) {
    
    if (!iftFullFHeap(H)) {
        H->last++;
        H->node[H->last] = node;
        H->color[node]   = IFT_GRAY;
        H->pos[node]     = H->last;
        iftGoUpFHeap(H, H->last);
        return 1;
    } else {
        iftWarning("FHeap is full","iftInsertFHeap");
        return 0;
    }
    
}

int iftRemoveFHeap(iftFHeap *H) {
    int node= IFT_NIL;
    
    if (!iftEmptyFHeap(H)) {
        node = H->node[0];
        H->pos[node]   = -1;
        H->color[node] = IFT_BLACK;
        H->node[0]     = H->node[H->last];
        H->pos[H->node[0]] = 0;
        H->node[H->last] = -1;
        H->last--;
        iftGoDownFHeap(H, 0);
    }else{
        iftWarning("FHeap is empty","iftRemoveFHeap");
    }
    
    return node;
    
}

void iftRemoveFHeapElem(iftFHeap *H, int pixel){
    
    if(H->pos[pixel] == -1)
        iftError("Element is not in the Heap", "iftRemoveFHeapElem");
    
    float aux = H->value[pixel];
    
    if(H->removal_policy == IFT_MINVALUE)
        H->value[pixel] = IFT_INFINITY_FLT_NEG;
    else
        H->value[pixel] = IFT_INFINITY_FLT;
    
    iftGoUpFHeap(H, H->pos[pixel]);
    iftRemoveFHeap(H);
    
    H->value[pixel] = aux;
    H->color[pixel] = IFT_WHITE;
    
}

void iftResetFHeap(iftFHeap *H)
{
    int i;
    
    for (i=0; i < H->n; i++) {
        H->color[i] = IFT_WHITE;
        H->pos[i]   = -1;
        H->node[i] = -1;
    }
    H->last = -1;
}


