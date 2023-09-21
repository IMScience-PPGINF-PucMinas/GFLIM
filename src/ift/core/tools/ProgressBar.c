#include "ift/core/tools/ProgressBar.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"
#include "iftMemory.h"


iftProgressBar* iftCreateProgressBar(float start, float end) {
    return iftCreateProgressBarSize(start, end, 50);
}


iftProgressBar* iftCreateProgressBarSize(float start, float end, int size) {
    iftProgressBar* bar = (iftProgressBar*) iftAlloc(1, sizeof(iftProgressBar));
    
    if(end<=start) {
        iftError("Invalid interval", "iftCreateProgressBarSize");
    }
    
    bar->start = bar->cur = start;
    bar->end = end;
    bar->size = size;
    
    return bar;
}


void iftDestroyProgressBar(iftProgressBar** pbar) {
    iftFree(*pbar);
    *pbar = NULL;
}


void iftProgressBarUpdate(iftProgressBar* bar, int cur) {
    bar->cur = iftMin(iftMax(cur, bar->start), bar->end);
    fprintf(stdout, "\r|");
    float step = (bar->end-bar->start)/bar->size;
    for (int i = 0; i < bar->size; ++i) {
        if ((step*i)>=bar->cur) {
            fprintf(stdout, " ");
        }
        else {
            fprintf(stdout, "=");
        }
    }
    fprintf(stdout, "|\t%4.2f%%\t(%4.2f/%4.2f)", 100.0*(bar->cur - bar->start)/(bar->end - bar->start), bar->cur, bar->end);
    fflush(stdout);
}

