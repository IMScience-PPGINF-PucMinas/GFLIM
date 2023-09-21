//
// Created by Samuel Martins on 10/12/18.
//
#ifndef IFT_STRING_H
#define IFT_STRING_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/SList.h"



/**
 * @brief Removes the occurences of a character c in the end of the string.
 * @author Peixinho
 * @param s
 * @param c
 */
void iftRightTrim(char* s, char c);


/**
 * @brief Removes the occurences of a character c in the beginning of the string.
 * @author Peixinho
 * @param s
 * @param c
 */
void iftLeftTrim(char* s, char c);


/**
 * @brief Removes the occurences of a character c in both ends of the string.
 * @author Peixinho
 * @param s
 * @param c
 */
void iftTrim(char* s, char c);


/**
 * @brief Creates a copy of the string @p <b>str</b> returning it. If @p <b>str</b> is NULL, returns NULL.
 * @author Samuel Martins
 * @date Oct 23, 2015
 * @ingroup String
 *
 * @param str The string to be copied.
 * @return The copied string.
 */
char *iftCopyString(const char *str, ...);


/**
 * @brief Concatenates strings.
 * @author Samuel Martins
 * @date Nov 11st, 2015
 * @ingroup String
 *
 * @param n Number of strings to be concatenated.
 * @param ... The @p <b>n</b> strings to be concatenated.
 * @return The concatenated string.

 * @warning If the number of strings is Zero, a NULL pointer will be returned.
 */
char *iftConcatStrings(int n, ...);



/**
 * @brief Compares two strings.
 * @author Samuel Martins
 * @date Dec 10, 2015
 * @ingroup String
 *
 * @param str1 First String.
 * @param str2 Second String.
 * @return True if they are equal. False, otherwise.
 */
bool iftCompareStrings(const char *str1, const char *str2);


/**
 * @brief Checks if the strings ends with a suffix.
 * @author Samuel Martins
 * @date Dec 10, 2015
 * @ingroup String
 *
 * @param string The string to be compared
 * @param string The suffix to be compared
 * @return True/False.
 */
bool iftEndsWith(const char *str, const char *suffix);


/**
 * @brief Checks if the strings starts with a prefix.
 * @author Samuel Martins
 * @date Dec 10, 2015
 * @ingroup String
 *
 * @param string The string to be compared
 * @param string The prefix to be compared
 * @return True/False.
 */
bool iftStartsWith(const char *str, const char *prefix);


/**
 * @brief Returns a String without the given @p <b>prefix</b>. If the prefix does not match or is "", returns the original string.\n
 * @author Samuel Martins
 * @date Dec 10, 2015
 * @ingroup String
 *
 * @param string The string to be compared
 * @param prefix The prefix to be compared
 * @return The string without prefix.
 */
char *iftRemovePrefix(const char *str, const char *prefix);


/**
 * @brief Returns a String (copy) without the given <code><b>suffix</b></code>. If the suffix does not match or is "", returns the original string.\n
 * @author Samuel Martins
 * @date Dec 10, 2015
 * @ingroup String
 *
 * @param string The string to be compared
 * @param suffix The suffix to be compared
 * @return The string (copy) without suffix.
 */
char *iftRemoveSuffix(const char *str, const char *suffix);


/**
 * @brief Splits a String and returns all substrings into a String Doubly Linked List.
 * @author Samuel Martins
 * @date Sep 16, 2015
 * @ingroup String
 *
 * Splits the string <phrase> for the <delimiter> and returns all substrings into a
 * String Doubly Linked List.
 * LIMITATION: the function supports substrings of length 512 chars
 *
 * @param phrase String to be splitted.
 * @param delimiter Substring where the <phrase> will be splitted.
 * @return The Linked List with all substrings.
 */
iftSList *iftSplitString(const char *phrase, const char *delimiter);


/**
 * @brief Splits a String and returns the sub-string of a given position.
 * @author Samuel Martins
 * @ingroup String
 *
 * Splits the key <phrase> for the <delimiter> and return the sub-string of
 * the position <position>.
 * YOU MUST DEALLOCATE THE RETURNED STRING AFTER.
 * LIMITATION: the function supports at most 4096 substrings of length 512 chars
 *
 * Position means the segment position in the sub-strings. The first segment is 0.
 * The method accepts negative position, since abs(position) <= number of splitted sub-strings
 * Negative positions access from back to front
 *    ---> e.g. abc_def, '_' --> [0] - abc, [1] = def, [-1] = def, [-2] = abc.
 *
 * @param phrase String to be splitted.
 * @param delimiter Sub-string where the <phrase> will be splitted.
 * @param position Position of the required sub-string.
 * @return Required sub-string (char*).
 *
 * @attention See a demo in <i>ift_dir/demo/Miscellaneous/iftTestSplitString</i>
 */
char *iftSplitStringAt(const char *phrase, const char *delimiter, long position);


/**
 * @brief Converts a string to Uppercase.
 * @author Samuel Martins
 * @date Dec 7, 2015
 * @ingroup String
 *
 * @param str The string to be converted.
 * @return The uppercase string.
 *
 * @exception str is NULL.
 */
char *iftUpperString(const char *str);


/**
 * @brief Converts a string to Lowercase.
 * @author Samuel Martins
 * @date Dec 7, 2015
 * @ingroup String
 *
 * @param str The string to be converted.
 * @return The lowercase string.
 */
char *iftLowerString(const char *str);


/**
 * @brief Replace all occurrences of the substring <b>old_sub</b> with <b>new_sub</b> in the string <b>str</b>.
 * @author Samuel Martins
 * @date Mar 9, 2016
 * @ingroup String
 * @note See a demo in @ref iftTestReplaceString.c
 * 
 * @param str The string to be replaced.
 * @param old_sub Old substring.
 * @param new_sub New substring.
 * @return The replaced string.
 */
char *iftReplaceString(const char *str, const char *old_sub, const char *new_sub);


#ifdef __cplusplus
}
#endif

#endif

