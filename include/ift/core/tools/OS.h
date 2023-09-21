//
// Created by Samuel Martins on 08/12/18.
//

#ifndef IFT_OS_H
#define IFT_OS_H

#include <sys/ioctl.h>
#include <unistd.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef enum {
    IFT_LINUX, IFT_MACOS, IFT_WIN32, IFT_WIN64
} iftOS;



/**
* @brief Get the Operational System. Operation Systems supported: IFT_LINUX, IFT_MACOS, IFT_WIN32, IFT_WIN64
* @author Samuka Martins
* @date Jul 18, 2017
* @return The enumeration of the Operation System
*/
iftOS iftGetOS(void);


/**
 * @brief Runs a program.
 * @param prog Program name.
 * @param args Arguments of the program.
 *
 * @note args can be built on-the-fly.
 */
void iftRunProgram(const char *prog, const char *args, ...);

/**
 * @brief Get the size of the terminal in characters (n_rows and n_cols).
 * @author Cesar Castelo
 * @date Mar 04, 2019
 * @return n_rows and n_rcols.
 */
void iftGetTerminalWindowSize(int *n_rows, int *n_cols);

/**
 * @brief Prints a line of characters which occupies the entire terminal width.
 * @author Cesar Castelo
 * @date Mar 04, 2019
 * @param character to be printed.
 */
void iftPrintSeparatorLineInTerminal(char character);


#ifdef __cplusplus
}
#endif

#endif //IFT_OS_H
