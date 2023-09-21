//
// Created by Samuel Martins on 2019-01-04.
//

#ifndef IFT_FSET_H
#define IFT_FSET_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


typedef struct ift_fset {
    float elem;
    struct ift_fset *next;
} iftFSet;

void    iftInsertFSet(iftFSet **S, float elem);
float   iftRemoveFSet(iftFSet **S);
void    iftRemoveFSetElem(iftFSet **S, float elem);
void    iftDestroyFSet(iftFSet **S);

#ifdef __cplusplus
}
#endif

#endif //IFT_FSET_H
