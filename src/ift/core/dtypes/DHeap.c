#include "ift/core/dtypes/DHeap.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"


iftDHeap *iftCreateDHeap(int n, double *value) {
    iftDHeap *H = NULL;
    int i;
    
    if (value == NULL) {
        iftError("Cannot create heap without priority value map", "iftCreateDHeap");
    }
    
    H = (iftDHeap *) iftAlloc(1, sizeof(iftDHeap));
    if (H != NULL) {
        H->n       = n;
        H->value   = value;
        H->color   = (char *) iftAlloc(sizeof(char), n);
        H->node    = (int *) iftAlloc(sizeof(int), n);
        H->pos     = (int *) iftAlloc(sizeof(int), n);
        H->last    = -1;
        H->removal_policy = MINVALUE;
        if (H->color == NULL || H->pos == NULL || H->node == NULL)
            iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDHeap");
        for (i = 0; i < H->n; i++) {
            H->color[i] = IFT_WHITE;
            H->pos[i]   = -1;
            H->node[i] = -1;
        }
    }
    else
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDHeap");
    
    return H;
}

void iftDestroyDHeap(iftDHeap **H) {
    iftDHeap *aux = *H;
    if (aux != NULL) {
        if (aux->node != NULL) iftFree(aux->node);
        if (aux->color != NULL) iftFree(aux->color);
        if (aux->pos != NULL)   iftFree(aux->pos);
        iftFree(aux);
        *H = NULL;
    }
}


void  iftGoUpDHeap_(iftDHeap *H, int i, bool show) {
    int j = iftDad(i);
    
    if(show) printf("i:%d j:%d, H->node[i]:%d, H->node[j]:%d \n",i,j, H->node[i],H->node[j]);
    if(show) printf("H->value[H->node[i]]:%lf H->value[H->node[j]]:%lf \n",H->value[H->node[i]], H->value[H->node[j]]);

    if(H->removal_policy == MINVALUE){
        while ((j >= 0) && (H->value[H->node[j]] > H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
    else{ /* removal_policy == MAXVALUE */
        
        while ((j >= 0) && (H->value[H->node[j]] < H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
}

void  iftGoUpDHeap(iftDHeap *H, int i) {
    int j = iftDad(i);
    
    if(H->removal_policy == MINVALUE){
        while ((j >= 0) && (H->value[H->node[j]] > H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
    else{ /* removal_policy == MAXVALUE */
        
        while ((j >= 0) && (H->value[H->node[j]] < H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
}

void iftGoDownDHeap(iftDHeap *H, int i) {
    int j, left = iftLeftSon(i), right = iftRightSon(i);
    
    j = i;
    if(H->removal_policy == MINVALUE){
        
        if ((left <= H->last) &&
            (H->value[H->node[left]] < H->value[H->node[i]]))
            j = left;
        if ((right <= H->last) &&
            (H->value[H->node[right]] < H->value[H->node[j]]))
            j = right;
    }
    else{ /* removal_policy == MAXVALUE */
        
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
        iftGoDownDHeap(H, j);
    }
}

char iftFullDHeap(iftDHeap *H) {
    if (H->last == (H->n - 1))
        return 1;
    else
        return 0;
}

char iftEmptyDHeap(iftDHeap *H) {
    if (H->last == -1){
        return 1;
    }else{
        return 0;
    }
}


char iftInsertDHeap(iftDHeap *H, int node) {
    
    if (!iftFullDHeap(H)) {
        H->last++;
        H->node[H->last] = node;
        H->color[node]   = IFT_GRAY;
        H->pos[node]     = H->last;
        iftGoUpDHeap(H, H->last);
        return 1;
    } else {
        iftWarning("DHeap is full","iftInsertDHeap");
        return 0;
    }
    
}

int iftRemoveDHeap(iftDHeap *H) {
    int node= IFT_NIL;
    
    if (!iftEmptyDHeap(H)) {
        node = H->node[0];
        H->pos[node]   = -1;
        H->color[node] = IFT_BLACK;
        H->node[0]     = H->node[H->last];
        H->pos[H->node[0]] = 0;
        H->node[H->last] = -1;
        H->last--;
        iftGoDownDHeap(H, 0);
    }else{
        iftWarning("DHeap is empty","iftRemoveDHeap");
    }
    
    return node;
    
}

void    iftRemoveDHeapElem(iftDHeap *H, int pixel){
    
    if(H->pos[pixel] == -1)
        iftError("Element is not in the Heap", "iftRemoveDHeapElem");
    
    double aux = H->value[pixel];
    
    if(H->removal_policy == MINVALUE)
        H->value[pixel] = IFT_INFINITY_DBL_NEG;
    else
        H->value[pixel] = IFT_INFINITY_DBL;
    
    iftGoUpDHeap(H, H->pos[pixel]);
    iftRemoveDHeap(H);
    
    H->value[pixel] = aux;
    H->color[pixel] = IFT_WHITE;
    
}

void iftResetDHeap(iftDHeap *H)
{
    int i;
    
    for (i=0; i < H->n; i++) {
        H->color[i] = IFT_WHITE;
        H->pos[i]   = -1;
        H->node[i] = -1;
    }
    H->last = -1;
}


