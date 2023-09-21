//
// Created by Samuel Martins on 08/12/18.
//

#ifndef IFT_REGEX_H
#define IFT_REGEX_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include <stdbool.h>

/**
 * @brief Checks if a string matches with a regular expression.
 * @author Samuel Martins
 * @date Dec 11, 2015
 * @ingroup Regex
 * @note See a complete demo in @ref iftRegexExamples.c and @ref iftTestRegex.c
 *
 * @param str The string to be checked.
 * @param regex_pattern The regular expression.
 * @return True/False.
 *
 *
 * @exception str is NULL.
 * @exception regex_pattern is NULL.
 * @exception regex_pattern does not compile.
 */
bool iftRegexMatch(const char *str, const char *regex_pattern, ...);

#ifdef __cplusplus
}
#endif

#endif //IFT_REGEX_H
