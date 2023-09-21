#include "ift/core/tools/Dialog.h"

void iftAssertAux(bool condition, const char* msg, const char* func, ...) {
    if(!condition) {
        va_list args;
        char final_msg[4096];
        va_start(args, func);
        vsprintf(final_msg, msg, args);
        va_end(args);
        iftError(final_msg, func);
    }
}


void iftError(const char *msg, const char *func, ...) {
    va_list args;
    char    final_msg[4096];
    
    va_start(args, func);
    vsprintf(final_msg, msg, args);
    va_end(args);
    
    fprintf(stderr, "\nError in %s: \n%s\n", func, final_msg);
    fflush(stdout);
    exit(-1);
}


void iftWarning(const char *msg, const char *func, ...) {
    va_list args;
    char    final_msg[4096];
    
    va_start(args, func);
    vsprintf(final_msg, msg, args);
    va_end(args);
    
    fprintf(stdout, "\nWarning in %s: \n%s\n", func, final_msg);
}


void iftDeprecated(const char *old_function, const char *new_function, const char *message) {
    if(message != NULL) {
        fprintf(stderr,
                IFT_COLOR_YELLOW "\n-----\nThe function \"%s\" is deprecated.\nUse the function \"%s\".\n%s\n\n-----\n\n" IFT_COLOR_RESET,
                old_function, new_function, message);
    } else {
        fprintf(stderr,
                IFT_COLOR_YELLOW "\n-----\nThe function \"%s\" is deprecated.\nUse the function \"%s\".\n-----\n\n" IFT_COLOR_RESET,
                old_function, new_function);
    }
}

