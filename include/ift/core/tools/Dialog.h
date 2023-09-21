//
// Created by samuel on 13/12/18.
//

#ifndef IFT_DIALOG_H
#define IFT_DIALOG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


#define MSG_MEMORY_ALLOC_ERROR  "Cannot allocate memory space"
#define MSG_FILE_OPEN_ERROR  "Cannot open file %s"
#define MSG_FILE_FORMAT_ERROR  "Format error in line %d of file %s"

/**
 * @brief Color tags used to print messages with a different color.
 * @author Samuel Martins
 * @date Feb 29, 2016
 * @ingroup Dialog
 * @{
 */
#define IFT_COLOR_YELLOW "\x1B[33m"
#define IFT_COLOR_BLUE "\x1B[34m"
#define IFT_COLOR_RESET "\033[0m"
/** }@ */

/**
 * @brief Tries to assert if a condition is satisfied.
 * @warning This method is ignored in non debug mode.
 * @author Peixinho
 * @date Nov, 2016
 */
#ifdef IFT_DEBUG
#define iftAssert(cond, msg, func, ...) iftAssertAux(cond, msg, func, __VA_ARGS__)
#else
#define iftAssert(cond, msg, func, ...) (void)0
#endif


/**
 * @brief An exception is raised from function <b>func</b> with the error message <b>msg</b>. The program exits abnormally.
 * @date May, 2017
 * @ingroup Dialog
 * @author Peixinho
 *
 * @note The function accepts formatted string like printf. E.g: ("Invalid value %d", "function_name", value)
 */
void iftError(const char *msg, const char *func, ...);


/**
 * @brief A warning message is printed from function <b>func</b>.
 * @date Feb 29, 2016
 * @ingroup Dialog
 *
 * @note The function accepts formatted string like printf. E.g: ("Invalid value %d", "function_name", value)
 */
void iftWarning(const char *msg, const char *func, ...);


/**
 * @brief A warning message is printed indicating the depreciated function and the newest function.
 * @author Samuel Martins
 * @date Sep 14, 2015
 *
 * @param old_function Depreciated function.
 * @param new_function The newest function.
 */
void iftDeprecated(const char *old_function, const char *new_function, const char *message);

#ifdef __cplusplus
}
#endif

#endif //IFT_DIALOG_H
