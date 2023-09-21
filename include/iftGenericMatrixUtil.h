

#ifndef _IFTGENERICMATRIXUTIL_H
#define _IFTGENERICMATRIXUTIL_H

#include "iftGenericMatrix.h"
//#include "iftDistanceFunctions.h"


//inline double computeDistanceBetweenRows(iftGenericMatrix* iftGenericMatrix1, iftGenericMatrix* iftGenericMatrix2,
//                                         size_t indexRow_Source ,size_t indexRow_Target,
//                                         iftDistanceFunction distanceFunction,
//                                         iftArgumentList* argumentList){
//
//    size_t nCols = iftGenericMatrix2->numberColumns;
//    float *vec_source = ((float*)iftGenericMatrix1->matrixData->data) + (indexRow_Source*iftGenericMatrix1->numberColumns);
//    float *vec_target = ((float*)iftGenericMatrix2->matrixData->data) + (indexRow_Target*iftGenericMatrix2->numberColumns);
//
//    return distanceFunction(vec_source,vec_target,nCols,argumentList);
//}
//
//
//inline size_t findNearestRow(iftGenericMatrix* source, iftGenericMatrix* target,
//                             size_t indexRow_Target,
//                             iftDistanceFunction distanceFunction,
//                             iftArgumentList* argumentList){
//
//    double minDistance = DBL_MAX;
//    double currentDistance = 0;
//    size_t indexMin = 0;
//    for (size_t i = 0; i < source->numberRows; ++i) {
//        //printf("eiei2\n");
//        currentDistance = computeDistanceBetweenRows(source,target,i,indexRow_Target,distanceFunction,argumentList);
//        //printf("eiei1\n");
//        if(currentDistance < minDistance){
//            minDistance = currentDistance;
//            indexMin = i;
//        }
//    }
//    return indexMin;
//}
//
//inline double* computeAllDistancesBetweenRowAndiftGenericMatrix(iftGenericMatrix* source, iftGenericMatrix* target,
//                                                     size_t indexRow_Target,
//                                                     iftDistanceFunction distanceFunction,
//                                                     iftArgumentList* argumentList){
//    double* distances = (double*)calloc(source->numberRows,sizeof(double));
//    for (size_t i = 0; i < source->numberRows; ++i) {
//        distances[i] = computeDistanceBetweenRows(source,target,i,indexRow_Target,distanceFunction,argumentList);
//    }
//    return distances;
//}
//
////inline double* myInsertionSort(double* vector, size_t n){
////    double* auxVector = (double*)calloc(n,sizeof(double));
////    double aux;
////
////    for (size_t i = 0; i < n; ++i) {
////        auxVector[i] = vector[i];
////        for (size_t j = i; j > 0; --j) {
////
////            if(auxVector[j-1] > auxVector[j]){
////                aux = auxVector[j];
////                auxVector[j] = auxVector[j-1];
////                auxVector[j-1] =  aux;
////            }else{
////                break;
////            }
////
////        }
////    }
////    return auxVector;
////}
//
////option = 0 : vertical
////option = 1 : horizontal;
//inline iftGenericMatrix* computeiftGenericMatrixMean(iftGenericMatrix* matrix, int option){
//
//    iftGenericMatrix* output = NULL;
//
//    if(option == 0){
//        output = iftCreateGenericMatrix(1,matrix->numberColumns,matrix->numberElements);
//        for (size_t i = 0; i < matrix->numberRows; ++i) {
//            for (size_t j = 0; j < matrix->numberColumns; ++j) {
//                MATRIX_GET_ELEMENT_PO_AS(float,output,0,j) += MATRIX_GET_ELEMENT_PO_AS(float,matrix,i,j);
//            }
//        }
//
//        for (size_t j = 0; j < matrix->numberColumns; ++j) {
//            MATRIX_GET_ELEMENT_PO_AS(float,output,0,j) /= matrix->numberColumns;
//        }
//
//    }else if(option == 1){
//        output = iftCreateGenericMatrix(matrix->numberRows,1,matrix->numberElements);
//        for (size_t i = 0; i < matrix->numberRows; ++i) {
//            for (size_t j = 0; j < matrix->numberColumns; ++j) {
//                MATRIX_GET_ELEMENT_PO_AS(float,output,i,0) += MATRIX_GET_ELEMENT_PO_AS(float,matrix,i,j);
//            }
//        }
//
//        for (size_t j = 0; j < matrix->numberColumns; ++j) {
//            MATRIX_GET_ELEMENT_PO_AS(float,output,0,j) /= matrix->numberColumns;
//        }
//
//
//    }
//
//    return output;
//}


#endif //LIBFL_iftGenericMatrixUTIL_H
