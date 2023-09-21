#include "ift/core/dtypes/BasicDataTypes.h"

#include "ift/core/dtypes/DblArray.h"
#include "ift/core/dtypes/Dict.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/IntMatrix.h"
#include "ift/core/dtypes/StrArray.h"
#include "iftMatrix.h"


////////////// GENERIC TYPE MANIPULATION //////////////


iftGVal iftInitBoolGVal(bool val) {
    iftGVal gv;
    gv.type = IFT_BOOL_TYPE;
    gv.bool_val = val;
    return gv;
}


iftGVal iftInitCharGVal(char val) {
    iftGVal gv;
    gv.type = IFT_CHAR_TYPE;
    gv.char_val = val;
    return gv;
}


iftGVal iftInitUCharGVal(uchar val) {
    iftGVal gv;
    gv.type = IFT_UCHAR_TYPE;
    gv.uchar_val = val;
    return gv;
}


iftGVal iftInitStrGVal(const char *val) {
    char *str = (char *) calloc(strlen(val)+1, sizeof(char));
    strcpy(str, val);
    iftGVal gv;
    gv.type = IFT_STR_TYPE;
    gv.str_val = str;
    return gv;
}


iftGVal iftInitLongGVal(long val) {
    iftGVal gv;
    gv.type = IFT_LONG_TYPE;
    gv.long_val = val;
    return gv;
}


iftGVal iftInitULongGVal(ulong val) {
    iftGVal gv;
    gv.type = IFT_ULONG_TYPE;
    gv.ulong_val = val;
    return gv;
}


iftGVal iftInitDblGVal(double val) {
    iftGVal gv;
    gv.type = IFT_DBL_TYPE;
    gv.dbl_val = val;
    return gv;
}


iftGVal iftInitIntArrayGVal(iftIntArray *array) {
    iftGVal gv;
    gv.type = IFT_INT_ARRAY_TYPE;
    gv.int_array_val = array;
    return gv;
}


iftGVal iftInitDblArrayGVal(iftDblArray *array) {
    iftGVal gv;
    gv.type = IFT_DBL_ARRAY_TYPE;
    gv.dbl_array_val = array;
    return gv;
}


iftGVal iftInitStrArrayGVal(iftStrArray *array) {
    iftGVal gv;
    gv.type = IFT_STR_ARRAY_TYPE;
    gv.str_array_val = array;
    return gv;
}


iftGVal iftInitIntMatrixGVal(iftIntMatrix *mat) {
    iftGVal gv;
    gv.type = IFT_INT_MATRIX_TYPE;
    gv.int_matrix_val = mat;
    return gv;
}


iftGVal iftInitDblMatrixGVal(iftMatrix *mat) {
    iftGVal gv;
    gv.type = IFT_DBL_MATRIX_TYPE;
    gv.dbl_matrix_val = mat;
    return gv;
}


iftGVal iftInitStrMatrixGVal(iftStrMatrix *mat) {
    iftGVal gv;
    gv.type = IFT_STR_MATRIX_TYPE;
    gv.str_matrix_val = mat;
    return gv;
}


iftGVal iftInitDictGVal(iftDict *dict) {
    iftGVal gv;
    gv.type = IFT_DICT_TYPE;
    gv.dict_val = dict;
    return gv;
}


iftGVal iftInitPtrGVal(void *val) {
    iftGVal gv;
    gv.type = IFT_PTR_TYPE;
    gv.ptr_val = val;
    return gv;
}


void iftFreeGVal(iftGVal gv) {
    if (gv.type == IFT_STR_TYPE)
        free(gv.str_val);
    else if (gv.type == IFT_INT_ARRAY_TYPE)
        iftDestroyIntArray(&gv.int_array_val);
    else if (gv.type == IFT_DBL_ARRAY_TYPE)
        iftDestroyDblArray(&gv.dbl_array_val);
    else if (gv.type == IFT_STR_ARRAY_TYPE)
        iftDestroyStrArray(&gv.str_array_val);
    else if (gv.type == IFT_INT_MATRIX_TYPE)
        iftDestroyIntMatrix(&gv.int_matrix_val);
    else if (gv.type == IFT_DBL_MATRIX_TYPE)
        iftDestroyMatrix(&gv.dbl_matrix_val);
    else if (gv.type == IFT_STR_MATRIX_TYPE)
        iftDestroyStrMatrix(&gv.str_matrix_val);

    gv.ptr_val = NULL;
}


const char *iftCDataTypeToString(iftCDataType datatype) {
    switch(datatype) {
        case IFT_UNTYPED:
            return "untyped";
        case IFT_BOOL_TYPE:
            return "boolean";
        case IFT_CHAR_TYPE:
            return "char";
        case IFT_UCHAR_TYPE:
            return "unsigned char";
        case IFT_STR_TYPE:
            return "string (char*)";
        case IFT_INT_TYPE:
            return "int";
        case IFT_UINT_TYPE:
            return "unsigned int";
        case IFT_LONG_TYPE:
            return "long";
        case IFT_ULONG_TYPE:
            return "unsigned long";
        case IFT_FLT_TYPE:
            return "float";
        case IFT_DBL_TYPE:
            return "double";
        case IFT_INT_ARRAY_TYPE:
            return "iftIntArray";
        case IFT_DBL_ARRAY_TYPE:
            return "iftDblArray";
        case IFT_STR_ARRAY_TYPE:
            return "iftStrArray";
        case IFT_INT_MATRIX_TYPE:
            return "iftIntMatrix";
        case IFT_DBL_MATRIX_TYPE:
            return "iftMatrix";
        case IFT_STR_MATRIX_TYPE:
            return "iftStrMatrix";
        case IFT_DICT_TYPE:
            return "iftDict*";
        case IFT_PTR_TYPE:
            return "void*";
        default:
            return "unknown type"; // just to avoid compilation warnings
    }
}


char *iftGValToString(iftGVal gval) {
    char *str = (char *) calloc(IFT_STR_DEFAULT_SIZE, sizeof(char));

    switch(gval.type) {
        case IFT_BOOL_TYPE:
            sprintf(str, "%s", gval.bool_val ? "true" : "false");
            break;
        case IFT_CHAR_TYPE:
            sprintf(str, "\'%c\'", gval.char_val);
            break;
        case IFT_UCHAR_TYPE:
            sprintf(str, "\'%c\'", gval.char_val);
            break;
        case IFT_STR_TYPE:
            sprintf(str, "\"%s\"", gval.str_val);
            break;
        case IFT_LONG_TYPE:
            sprintf(str, "%ld", gval.long_val);
            break;
        case IFT_ULONG_TYPE:
            sprintf(str, "%lu", gval.ulong_val);
            break;
        case IFT_DBL_TYPE:
            sprintf(str, "%lf", gval.dbl_val);
            break;
        case IFT_PTR_TYPE:
            if (gval.ptr_val == 0x0) {
                sprintf(str, "null");
            }
            else sprintf(str, "%p", gval.ptr_val);
            break;
        default:
            strcpy(str, "");
            break;
    }

    return str;
}



bool iftCompareGVal(iftGVal val1, iftGVal val2) {
    if (val1.type == val2.type) {
        switch(val1.type) {
            case IFT_BOOL_TYPE:
                return (val1.bool_val == val2.bool_val);
            case IFT_CHAR_TYPE:
                return (val1.char_val == val2.char_val);
            case IFT_UCHAR_TYPE:
                return (val1.uchar_val == val2.uchar_val);
            case IFT_STR_TYPE:
                return (strcmp(val1.str_val, val2.str_val) == 0);
            case IFT_LONG_TYPE:
                return (val1.long_val == val2.long_val);
            case IFT_ULONG_TYPE:
                return (val1.ulong_val == val2.ulong_val);
            case IFT_DBL_TYPE:
                return (val1.dbl_val == val2.dbl_val);
            case IFT_PTR_TYPE:
                return (val1.ptr_val == val2.ptr_val);
            default:
                return false; // IFT_UNTYPED
        }
    }
    else return false;
}

iftGVal iftCopyGVal(iftGVal gval) {
    iftGVal out;
    out.type = gval.type;

    switch (gval.type) {
        case IFT_BOOL_TYPE:
            out = iftInitBoolGVal(gval.bool_val);
            break;
        case IFT_CHAR_TYPE:
            out = iftInitCharGVal(gval.char_val);
            break;
        case IFT_UCHAR_TYPE:
            out = iftInitUCharGVal(gval.uchar_val);
            break;
        case IFT_STR_TYPE:
            out = iftInitStrGVal(gval.str_val);
            break;
        case IFT_INT_TYPE:
        case IFT_LONG_TYPE:
            out = iftInitLongGVal(gval.long_val);
            break;
        case IFT_UINT_TYPE:
        case IFT_ULONG_TYPE:
            out = iftInitULongGVal(gval.ulong_val);
            break;
        case IFT_FLT_TYPE:
        case IFT_DBL_TYPE:
            out = iftInitDblGVal(gval.dbl_val);
            break;
        case IFT_INT_ARRAY_TYPE:
            out.int_array_val = iftCreateIntArray(gval.int_array_val->n);
            iftCopyIntArray(out.int_array_val->val, gval.int_array_val->val, gval.int_array_val->n);
            break;
        case IFT_DBL_ARRAY_TYPE:
            out.dbl_array_val = iftCopyDblArray(gval.dbl_array_val->val, gval.dbl_array_val->n);
            break;
        case IFT_STR_ARRAY_TYPE:
            out.str_array_val = iftCopyStrArray(gval.str_array_val->val, gval.str_array_val->n);
            break;
        case IFT_INT_MATRIX_TYPE:
            out.int_matrix_val = iftCopyIntMatrix(gval.int_matrix_val->val, gval.int_matrix_val->nrows, gval.int_matrix_val->ncols);
            break;
        case IFT_DBL_MATRIX_TYPE:
            out.dbl_matrix_val = iftCopyMatrix(gval.dbl_matrix_val);
            break;
        case IFT_STR_MATRIX_TYPE:
            out.str_matrix_val = iftCopyStrMatrix(gval.str_matrix_val->val, gval.str_matrix_val->nrows, gval.str_matrix_val->ncols);
            break;
        case IFT_DICT_TYPE:
            out.dict_val = iftCopyDict(gval.dict_val);
            break;
        case IFT_UNTYPED:
        case IFT_PTR_TYPE:
        default:
            out.ptr_val = gval.ptr_val;
            break;
    }

    return out;
}




bool iftGetBoolVal(iftGVal gval) {
    return gval.bool_val;
}


char iftGetCharVal(iftGVal gval) {
    return gval.char_val;
}


uchar iftGetUCharVal(iftGVal gval) {
    return gval.uchar_val;
}


char *iftGetStrVal(iftGVal gval) {
    return gval.str_val;
}


const char *iftGetConstStrVal(iftGVal gval) {
    return gval.str_val;
}


long iftGetLongVal(iftGVal gval) {
    return gval.long_val;
}


ulong iftGetULongVal(iftGVal gval) {
    return gval.ulong_val;
}


double iftGetDblVal(iftGVal gval) {
    return gval.dbl_val;
}


iftIntArray *iftGetIntArrayVal(iftGVal gval) {
    return gval.int_array_val;
}


iftDblArray *iftGetDblArrayVal(iftGVal gval) {
    return gval.dbl_array_val;
}


iftStrArray *iftGetStrArrayVal(iftGVal gval) {
    return gval.str_array_val;
}


iftIntMatrix *iftGetIntMatrixVal(iftGVal gval){
    return gval.int_matrix_val;
}


iftMatrix *iftGetDblMatrixVal(iftGVal gval) {
    return gval.dbl_matrix_val;
}


iftStrMatrix *iftGetStrMatrixVal(iftGVal gval) {
    return gval.str_matrix_val;
}


iftDict *iftGetDictVal(iftGVal gval) {
    return gval.dict_val;
}


void *iftGetPtrVal(iftGVal gval) {
    return gval.ptr_val;
}



void iftCopyVoxel(iftVoxel *src, iftVoxel *dst) {
    (*dst).x = (*src).x;
    (*dst).y = (*src).y;
    (*dst).z = (*src).z;
    (*dst).t = (*src).t;    
}


