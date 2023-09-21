//
// Created by deangeli on 4/30/17.
//

#ifndef IFT_IFTDATASETCONTEXT_H
#define IFT_IFTDATASETCONTEXT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftGenericLinkedList.h"
#include "iftDataSet.h"

typedef struct _iftDataSetContext {
    iftDataSet* dataSet;
    bool* eraseSample;
    char* contextName;
} iftDataSetContext;

typedef struct _iftDataSetContextManager {
    iftLinkedList *contextList;
    iftDataSetContext* contextReference;
    iftDataSetContext* currentContext;
} iftDataSetContextManager;


iftDataSetContextManager* createDataSetContextManager();

void iftDestroyDataSetContext(void *DataSetContext);
void copyDataSetInfo(iftDataSetContext* context1,iftDataSetContext* context2);


#ifdef __cplusplus
}
#endif

#endif //IFT_IFTDATASETCONTEXT_H
