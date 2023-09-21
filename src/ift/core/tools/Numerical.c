#include "ift/core/tools/Numerical.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"

iftIntArray *iftIntRange(int begin, int end, int inc) {
    int n = ((end - begin) / inc) + 1;
    
    iftIntArray *space = iftCreateIntArray(n);
    
    for (int i = 0; i < n; ++i) {
        space->val[i] = begin + (inc*i);
    }
    
    return space;
}


iftIntArray *iftIntRepeat(int x, int n_repetitions) {
    if (n_repetitions <= 0)
        iftError("Number of repetitions %d is <= 0", "iftIntRepeat", n_repetitions);
    
    iftIntArray *rep_arr = iftCreateIntArray(n_repetitions);

    #pragma omp parallel for
    for (long i = 0; i < rep_arr->n; i++)
        rep_arr->val[i] = x;
    
    return rep_arr;
}


iftDblArray *iftRange(double begin, double end, double inc) {
    int n = ((end-begin)/inc) + 1;
    iftDblArray * space = iftCreateDblArray(n);
    for (int i = 0; i < n; ++i) {
        space->val[i] = begin + inc*i;
    }
    return space;
}


iftDblArray *iftGeomRange(double begin, double end, double mul) {
    int n = ceil(iftLog(end/begin, mul)) + 1;
    iftDblArray * space = iftCreateDblArray(n);
    
    for (int i = 0; i < n; ++i) {
        space->val[i] = begin * pow(mul,i);
    }
    
    return space;
}
