//
// Created by deangeli on 5/5/17.
//

#ifndef _DATATRANSFERENCE_H_
#define _DATATRANSFERENCE_H_

//for some reasons memcpy is too slow
#define TRANSFER_DATA_COPY(buffDest,buffSrc,nBytes) \
        unsigned char* buffDestination = (unsigned char*)buffDest; \
        unsigned char* buffSource = (unsigned char*)buffSrc; \
        switch (nBytes) { \
            case 1:\
                buffDestination[0] = buffSource[0]; \
                break;\
            case 2:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                break; \
            case 3:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                break; \
            case 4:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                break; \
            case 5:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                break; \
            case 6:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                break; \
            case 7:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                break; \
            case 8:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                buffDestination[7] = buffSource[7]; \
                break; \
            case 9:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                buffDestination[7] = buffSource[7]; \
                buffDestination[8] = buffSource[8]; \
                break; \
            case 10:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                buffDestination[7] = buffSource[7]; \
                buffDestination[8] = buffSource[8]; \
                buffDestination[9] = buffSource[9]; \
                break; \
            case 11:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                buffDestination[7] = buffSource[7]; \
                buffDestination[8] = buffSource[8]; \
                buffDestination[9] = buffSource[9]; \
                buffDestination[10] = buffSource[10]; \
                break; \
            case 12:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                buffDestination[7] = buffSource[7]; \
                buffDestination[8] = buffSource[8]; \
                buffDestination[9] = buffSource[9]; \
                buffDestination[10] = buffSource[10]; \
                buffDestination[11] = buffSource[11]; \
                break; \
            case 13:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                buffDestination[7] = buffSource[7]; \
                buffDestination[8] = buffSource[8]; \
                buffDestination[9] = buffSource[9]; \
                buffDestination[10] = buffSource[10]; \
                buffDestination[11] = buffSource[11]; \
                buffDestination[12] = buffSource[12]; \
                 break; \
            case 14:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                buffDestination[7] = buffSource[7]; \
                buffDestination[8] = buffSource[8]; \
                buffDestination[9] = buffSource[9]; \
                buffDestination[10] = buffSource[10]; \
                buffDestination[11] = buffSource[11]; \
                buffDestination[12] = buffSource[12]; \
                buffDestination[13] = buffSource[13]; \
                break; \
            case 15:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                buffDestination[7] = buffSource[7]; \
                buffDestination[8] = buffSource[8]; \
                buffDestination[9] = buffSource[9]; \
                buffDestination[10] = buffSource[10]; \
                buffDestination[11] = buffSource[11]; \
                buffDestination[12] = buffSource[12]; \
                buffDestination[13] = buffSource[13]; \
                buffDestination[14] = buffSource[14]; \
                break; \
            case 16:\
                buffDestination[0] = buffSource[0]; \
                buffDestination[1] = buffSource[1]; \
                buffDestination[2] = buffSource[2]; \
                buffDestination[3] = buffSource[3]; \
                buffDestination[4] = buffSource[4]; \
                buffDestination[5] = buffSource[5]; \
                buffDestination[6] = buffSource[6]; \
                buffDestination[7] = buffSource[7]; \
                buffDestination[8] = buffSource[8]; \
                buffDestination[9] = buffSource[9]; \
                buffDestination[10] = buffSource[10]; \
                buffDestination[11] = buffSource[11]; \
                buffDestination[12] = buffSource[12]; \
                buffDestination[13] = buffSource[13]; \
                buffDestination[14] = buffSource[14]; \
                buffDestination[15] = buffSource[15]; \
                break; \
            default:\
                for (size_t i = 0; i < nBytes; ++i) { \
                    buffDestination[i] = buffSource[i]; \
                }\
                break; \
        }


inline void transferDataCopy(unsigned  char* buffDestination,const unsigned  char* buffSource, size_t nBytes){
    switch (nBytes){
        case 1:
            buffDestination[0] = buffSource[0];
            break;
        case 2:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            break;
        case 3:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            break;
        case 4:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            break;
        case 5:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            break;
        case 6:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            break;
        case 7:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            break;
        case 8:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            buffDestination[7] = buffSource[7];
            break;
        case 9:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            buffDestination[7] = buffSource[7];
            buffDestination[8] = buffSource[8];
            break;
        case 10:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            buffDestination[7] = buffSource[7];
            buffDestination[8] = buffSource[8];
            buffDestination[9] = buffSource[9];
            break;
        case 11:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            buffDestination[7] = buffSource[7];
            buffDestination[8] = buffSource[8];
            buffDestination[9] = buffSource[9];
            buffDestination[10] = buffSource[10];
            break;
        case 12:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            buffDestination[7] = buffSource[7];
            buffDestination[8] = buffSource[8];
            buffDestination[9] = buffSource[9];
            buffDestination[10] = buffSource[10];
            buffDestination[11] = buffSource[11];
            break;
        case 13:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            buffDestination[7] = buffSource[7];
            buffDestination[8] = buffSource[8];
            buffDestination[9] = buffSource[9];
            buffDestination[10] = buffSource[10];
            buffDestination[11] = buffSource[11];
            buffDestination[12] = buffSource[12];
            break;
        case 14:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            buffDestination[7] = buffSource[7];
            buffDestination[8] = buffSource[8];
            buffDestination[9] = buffSource[9];
            buffDestination[10] = buffSource[10];
            buffDestination[11] = buffSource[11];
            buffDestination[12] = buffSource[12];
            buffDestination[13] = buffSource[13];
            break;
        case 15:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            buffDestination[7] = buffSource[7];
            buffDestination[8] = buffSource[8];
            buffDestination[9] = buffSource[9];
            buffDestination[10] = buffSource[10];
            buffDestination[11] = buffSource[11];
            buffDestination[12] = buffSource[12];
            buffDestination[13] = buffSource[13];
            buffDestination[14] = buffSource[14];
            break;
        case 16:
            buffDestination[0] = buffSource[0];
            buffDestination[1] = buffSource[1];
            buffDestination[2] = buffSource[2];
            buffDestination[3] = buffSource[3];
            buffDestination[4] = buffSource[4];
            buffDestination[5] = buffSource[5];
            buffDestination[6] = buffSource[6];
            buffDestination[7] = buffSource[7];
            buffDestination[8] = buffSource[8];
            buffDestination[9] = buffSource[9];
            buffDestination[10] = buffSource[10];
            buffDestination[11] = buffSource[11];
            buffDestination[12] = buffSource[12];
            buffDestination[13] = buffSource[13];
            buffDestination[14] = buffSource[14];
            buffDestination[15] = buffSource[15];
            break;
        default:
            for (size_t i = 0; i < nBytes; ++i) {
                buffDestination[i] = buffSource[i];
            }
            break;
    }
}

//do not use it...memmove is more effcient
#define TRANSFER_DATA_RIGHT_SHIFT(_buffDest,_buffSrc,nBytes,elementSize) \
    for (int k = nBytes-elementSize; k >= 0 ; k -= elementSize) {\
        TRANSFER_DATA_COPY(_buffDest+k, _buffSrc+k,elementSize); \
    }

inline void transferDataRightShift(void* _buffDest,void* _buffSrc, size_t nBytes, size_t elementSize){
    for (int k = nBytes-elementSize; k >= 0 ; k -= elementSize) {
        TRANSFER_DATA_COPY(_buffDest+k, _buffSrc+k,elementSize);
    }
}

#endif //LIBFL_DATATRANSFERENCE_H
