//
// Created by deangeli on 5/21/17.
//

#ifndef LIBFL_VECTORUTIL_H
#define LIBFL_VECTORUTIL_H

#include "iftGenericVector.h"

#define mergeVectorsSameType_macro(type, vector1, size1, vector2, size2, vector3) \
    if(vector1 == NULL && vector2 == NULL){ \
         vector3 = NULL; \
    }\
    else if(vector1 == NULL){ \
        vector3 = (type*)calloc(size2,sizeof(type)); \
        memmove(vector3,vector2, size2*(sizeof(type))); \
    }\
    else if(vector2 == NULL){\
        vector3 = (type*)calloc(size1,sizeof(type));\
        memmove(vector3,vector1,size1* (sizeof(type))); \
    } \
    else{\
        vector3 = (type *) calloc(size1+size2,sizeof(vector3)); \
        memmove(vector3,vector1, size1*(sizeof(type)) ); \
        type* shiftedVector = vector3 + size1; \
        memmove(shiftedVector,vector2, size2*(sizeof(type)) );\
    }

#define appendVectorsSameType_macro(type, vector_destination, vector_source, destinationShiftIndex, sizeSource) \
    if(vector_destination != NULL && vector_source != NULL){ \
        type* shiftedVector = vector_destination + destinationShiftIndex;\
        memmove(shiftedVector,vector_source, sizeSource*(sizeof(type)) );\
    }\



#define mergeVectorsDifferentTypes_macro(type1, vector1, size1, type2, vector2, size2, type3, vector3) \
    if(vector1 == NULL && vector2 == NULL){ \
         vector3 = NULL; \
    }\
    else if(vector1 == NULL){ \
        vector3 = (type3*)calloc(size2,sizeof(type3)); \
        for(int i=0; i<size2; i++){\
            vector3[i] = vector2[i];\
        }\
    }\
    else if(vector2 == NULL){\
        vector3 = (type3*)calloc(size1,sizeof(type3));\
        for(int i=0; i<size1; i++){\
            vector3[i] = vector1[i];\
        }\
    } \
    else{\
        vector3 = (type3 *) calloc(size1+size2,sizeof(type3)); \
        for(int i=0; i < size1; i++){\
            vector3[i] = vector1[i];\
        }\
        for(int i=0; i<size2; i++){\
            vector3[size1+i] = vector2[i];\
        }\
    }


#define writeGenericVectorPointerAsText(type, symbol_print_string, symbol_separator_char, symbol_end_char,vector, vectorSize, filePointer)\
    {\
        if (filePointer == NULL){iftError(MSG_FILE_OPEN_ERROR, "writeGenericVectorPointerAsText macro");}\
        if(vectorSize > 0){\
            bool isEmptySeparator = (symbol_separator_char == '\0') ;\
            if(isEmptySeparator == true){\
                for(int i=0; i< vectorSize; i++){\
                    fprintf(filePointer,symbol_print_string,vector[i]);\
                }\
            }else{\
                for(int i=0; i< vectorSize-1; i++){\
                    fprintf(filePointer,symbol_print_string,vector[i]);\
                    fprintf(filePointer,"%c",symbol_separator_char);\
                }\
                fprintf(filePointer,symbol_print_string, vector[vectorSize-1]);\
            }\
            bool isEmptyEnd = (symbol_end_char == '\0') ;\
            if(isEmptyEnd != true){\
                fprintf(filePointer,"%c",symbol_end_char);\
            }\
        }\
    }


#define readGenericVectorPointerAsText(type, symbol_print_string, symbol_separator_char, symbol_end_char,vector, vectorSize, filePointer)\
    {\
        if (filePointer == NULL){iftError(MSG_FILE_OPEN_ERROR, "readGenericVectorPointerAsText macro");}\
        long int startPosition = ftell(filePointer);\
        unsigned long long iter=0;unsigned long long iterMax = 10000000000;\
        type aux;char aux_char;\
        if(vectorSize > 0){\
            do{\
                if(fscanf (fp,symbol_print_string, &aux)!= 1){iftError("Error","readGenericVectorPointerAsText macro 1");}\
                if(fscanf (fp,"%c",&aux_char)!= 1){iftError("Error","readGenericVectorPointerAsText macro 2");}\
                vector[iter] = aux;\
                iter++;\
            }while(aux_char != symbol_end_char);\
        }else{\
            if(vector != NULL){iftWarning("vector is not NULL\n","readGenericVectorPointerAsText macro");}\
            do{\
                if(fscanf (fp,symbol_print_string, &aux)!= 1){iftWarning("Error","readGenericVectorPointerAsText macro 3");}\
                if(fscanf (fp,"%c",&aux_char)!= 1){iftWarning("Error","readGenericVectorPointerAsText macro 4");}\
                iter++;\
            }while(aux_char != symbol_end_char && iter < iterMax);\
            if(iter >= iterMax){\
                iftWarning("Maximum number of iterations was reached","readGenericVectorPointerAsText macro 5");\
            }else{\
                vectorSize = iter;\
                vector = iftAlloc(iter,sizeof(type));\
                fseek ( filePointer , startPosition , SEEK_SET );\
                iter=0;\
                do{\
                    if(fscanf (fp,symbol_print_string, &aux)!= 1){iftWarning("Error","readGenericVectorPointerAsText macro 6");}\
                    if(fscanf (fp,"%c",&aux_char)!= 1){iftWarning("Error","readGenericVectorPointerAsText macro 7");}\
                    vector[iter] = aux;\
                    iter++;\
                }while(aux_char != symbol_end_char);\
            }\
        }\
    }

#define writeGenericVectorPointerAsBinary(type, vector, vetorSize, filePointer)\
    {\
        if (filePointer == NULL){iftError(MSG_FILE_OPEN_ERROR, "writeGenericVectorPointerAsBinary macro");}\
        if(vetorSize > 0){\
            if(fwrite ((void*)vector , sizeof(type), vetorSize, filePointer) != vetorSize){iftError("Error","writeGenericVectorPointerAsBinary macro 1");};\
        }else{\
            while(fwrite ((void*)vector , sizeof(type), 1, filePointer) != 1);\
        }\
    }

#define readGenericVectorPointerAsBinary(type, vector, vetorSize, filePointer)\
    {\
        if (filePointer == NULL){iftError(MSG_FILE_OPEN_ERROR, "writeGenericVectorPointerAsBinary macro");}\
        unsigned long long iter = 0;\
        if(vetorSize > 0){\
            if(fread((void*)vector,sizeof(type),vetorSize,filePointer) != vetorSize){\
                iftWarning("Error while reading data","readGenericVectorPointerAsBinary macro");\
            }\
        }else{\
            while(fread( (&(vector[iter])),sizeof(type),1,filePointer) != 1){iter++;}\
            vetorSize = iter;\
        }\
    }


inline double * myInsertionSort(double* vector, size_t n,size_t * indecesOrdered){
    double* auxVector = (double*)calloc(n,sizeof(double));
    double aux;
    size_t auxIndex;
    for (size_t i = 0; i < n; ++i) {
        auxVector[i] = vector[i];
        indecesOrdered[i] = i;
        for (size_t j = i; j > 0; --j) {

            if(auxVector[j-1] > auxVector[j]){
                aux = auxVector[j];
                auxVector[j] = auxVector[j-1];
                auxVector[j-1] =  aux;

                auxIndex = indecesOrdered[j];
                indecesOrdered[j] = indecesOrdered[j-1];
                indecesOrdered[j-1] = auxIndex;
            }else{
                break;
            }

        }
    }
    return auxVector;
}

inline double* mergeVectors(double* vector1, size_t size1, double* vector2, size_t size2){
    double *mergedVector = NULL;
    mergeVectorsSameType_macro( double, vector1, size1, vector2, size2, mergedVector);
    return mergedVector;
}


inline void myInsertionSortInplace(double* vector, size_t n,size_t * indecesOrdered){
    double aux;
    size_t auxIndex;
    for (size_t i = 0; i < n; ++i) {
        indecesOrdered[i] = i;
        for (size_t j = i; j > 0; --j) {

            if(vector[j-1] > vector[j]){
                aux = vector[j];
                vector[j] = vector[j-1];
                vector[j-1] =  aux;

                auxIndex = indecesOrdered[j];
                indecesOrdered[j] = indecesOrdered[j-1];
                indecesOrdered[j-1] = auxIndex;
            }else{
                break;
            }

        }
    }
}

/*do{\
//            fscanf (fp,"%d%c", &aux_int,&aux_char);\
//            printf("%d%c\n",aux_int,aux_char);\
//            iter++;\
//        }while(aux_char != symbol_end_char && iter < iterMax)\
//        if(iter >= iterMax){\
//            iftWarning("Maximum number of iterations was reached","readGenericVectorPointerAsText macro");\
//        }\ */



//inline void myInsertionSortInplace(int* vector, size_t n){
//    int aux;
//    for (size_t i = 0; i < n; ++i) {
//        for (size_t j = i; j > 0; --j) {
//
//            if(vector[j-1] > vector[j]){
//                aux = vector[j];
//                vector[j] = vector[j-1];
//                vector[j-1] =  aux;
//            }else{
//                break;
//            }
//
//        }
//    }
//}



//inline iftGenericVector* findUniquesInIntegerVector(iftGenericVector* vector,bool ordered){
//    iftGenericVector* uniques = iftCreateNullVector(vector->size,sizeof(int));
//
//    VECTOR_GET_ELEMENT_AS(int,uniques,0) = VECTOR_GET_ELEMENT_AS(int,vector,0);
//    size_t sizeUniquesVector = 1;
//    size_t countEquals;
//    for (size_t vectorIndex = 1; vectorIndex < vector->size; ++vectorIndex) {
//        countEquals = 0;
//        for (size_t uniquesVectorIndex = 0; uniquesVectorIndex < sizeUniquesVector; ++uniquesVectorIndex) {
//            if(VECTOR_GET_ELEMENT_AS(int,vector,vectorIndex) ==
//               VECTOR_GET_ELEMENT_AS(int,uniques,uniquesVectorIndex)){
//                countEquals++;
//            }
//        }
//        if(countEquals == 0){
//            VECTOR_GET_ELEMENT_AS(int,uniques,sizeUniquesVector) = VECTOR_GET_ELEMENT_AS(int,vector,vectorIndex);
//            sizeUniquesVector++;
//        }
//    }
//    uniques->size = sizeUniquesVector;
//    iftShrinkToFit(uniques);
//    if(ordered){
//        myInsertionSortInplace((int*)uniques->data,uniques->size);
//    }
//    return uniques;
//}



#endif //LIBFL_VECTORUTIL_H
