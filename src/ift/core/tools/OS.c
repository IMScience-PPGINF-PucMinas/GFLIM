#include "ift/core/tools/OS.h"

#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftMemory.h"



iftOS iftGetOS(void) {
    iftOS os;
    
    #if defined(__linux) || defined (linux) || defined (__linux__)
        os = IFT_LINUX;
    #elif defined(__MACH__) || defined (__APPLE__)
        os = IFT_MACOS;
    #elif defined(_WIN32)
        os = IFT_WIN32;
    #elif defined(_WIN64)
        os = IFT_WIN64;
    #else
        iftError("Operational System not supported!", "iftGetOS");
    #endif
    
    return os;
}



void iftRunProgram(const char *prog, const char *format, ...) {
    va_list vaargs;
    char prog_args[IFT_STR_DEFAULT_SIZE];

    va_start(vaargs, format);
    vsprintf(prog_args, format, vaargs);
    va_end(vaargs);

    char *cmd = iftConcatStrings(3, prog, " ", prog_args);
    if (system(cmd) != 0)
        iftError("System command error executing program %s. Try to compile it", "iftRunProgram", prog);
    iftFree(cmd);
}

void iftGetTerminalWindowSize(int *n_rows, int *n_cols)
{
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

    *n_rows = w.ws_row;
    *n_cols = w.ws_col;
}

void iftPrintSeparatorLineInTerminal(char character)
{
    int nTermimalRows, nTerminalCols;
    iftGetTerminalWindowSize(&nTermimalRows, &nTerminalCols);

    for(int i = 0; i < nTerminalCols; i++)
        printf("%c", character);
    printf("\n");
}
