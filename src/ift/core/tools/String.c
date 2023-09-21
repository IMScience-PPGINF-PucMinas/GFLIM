#include "ift/core/tools/String.h"

#include "ift/core/io/Stream.h"
#include "iftMemory.h"



void iftRightTrim(char* s, char c) {

    int idx = strlen(s) - 1;

    while(s[idx] == c) {
        idx--;
    }

    s[idx+1] = '\0';

}


void iftLeftTrim(char* s, char c) {
    int idx = 0;

    while (s[idx] == c) {
        idx++;
    }

    if(idx>0)
        memmove(s, s+idx, strlen(s) - idx);
}


void iftTrim(char* s, char c) {
    iftRightTrim(s, c);
    iftLeftTrim(s, c);
}


char *iftCopyString(const char *format, ...) {
    va_list args;
    char str[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(str, format, args);
    va_end(args);
    
    char *copy = iftAllocCharArray(strlen(str) + 1);
    strcpy(copy, str);

    return copy;
}


char *iftConcatStrings(int n, ...) {
    if (n <= 0)
        iftError("Number of Strings to be concatenated is <= 0", "iftConcatStrings");

    size_t out_str_size = 1; // '\0'

    // Counts the size of the concatenated string
    va_list strings;
    va_start(strings, n);
    for (int i = 0; i < n; i++)
        out_str_size += strlen(va_arg(strings, char*));
    va_end(strings);

    char *concat_str = iftAllocCharArray(out_str_size);

    va_start(strings, n);
    for (int i = 0; i < n; i++)
        strcat(concat_str, va_arg(strings, char*));
    va_end(strings);

    return concat_str;
}


bool iftCompareStrings(const char *str1, const char *str2) {
    if (str1 == NULL)
        iftError("First String is NULL", "iftCompareStrings");
    if (str2 == NULL)
        iftError("Second String is NULL", "iftCompareStrings");

    return (strcmp(str1, str2) == 0);
}


bool iftEndsWith(const char *str, const char *suffix) {
    if (str == NULL)
        iftError("String is NULL", "iftEndsWith");
    if (suffix == NULL)
        iftError("Suffix is NULL", "iftEndsWith");

    size_t len_suffix = strlen(suffix);
    size_t len_str    = strlen(str);

    if (len_suffix <= len_str) {
        size_t shift = len_str - len_suffix;
        return (strncmp(str+shift, suffix, len_suffix) == 0);
    }
    else
        return false;
}


bool iftStartsWith(const char *str, const char *prefix) {
    if (str == NULL)
        iftError("String is NULL", "iftStartsWith");
    if (prefix == NULL)
        iftError("Prefix is NULL", "iftStartsWith");

    size_t len_prefix = strlen(prefix);
    size_t len_str    = strlen(str);

    if (len_prefix <= len_str)
        return (strncmp(str, prefix, len_prefix) == 0);        
    else
        return false;
}


char *iftRemovePrefix(const char *str, const char *prefix) {
    if (str == NULL)
        iftError("String is NULL", "iftRemovePrefix");
    if (prefix == NULL)
        iftError("Prefix is NULL", "iftRemovePrefix");

    if (!iftCompareStrings(prefix, "") && iftStartsWith(str, prefix)) {
        size_t shift = strlen(prefix);
        return (iftCopyString(str + shift));   
    }
    else 
        return (iftCopyString(str));
}


char *iftRemoveSuffix(const char *str, const char *suffix) {
    if (str == NULL)
        iftError("String is NULL", "iftRemoveSuffix");
    if (suffix == NULL)
        iftError("Suffix is NULL", "iftRemoveSuffix");

    if (!iftCompareStrings(suffix, "") && iftEndsWith(str, suffix)) {
        size_t shift = strlen(str) - strlen(suffix);
        char *out_str = iftCopyString(str);
        out_str[shift] = '\0';

        return (out_str);   
    }
    else {
        return iftCopyString(str);
    }
}


iftSList *iftSplitString(const char *phrase, const char *delimiter) {
    if (phrase == NULL)
        iftError("String to be splitted is NULL", "iftSplitString");
    if (delimiter == NULL)
        iftError("Delimiter is NULL", "iftSplitString");

    char *buf       = iftAllocString(strlen(phrase)+1);
    const char *pr  = phrase;
    const char *loc = strstr(pr, delimiter); // pointer to first delimiter occurrence in the string

    size_t length = strlen(delimiter);
    size_t bytes;

    iftSList *SL = iftCreateSList();

    // build a list of sub-strings
    while (loc != NULL) {
        bytes = loc - pr;
        strncpy(buf, pr, bytes);
        buf[bytes] = '\0'; // \0 character must be added manually because strncpy does not do that, as opposed to other functions such as strcpy and sprintf

        iftInsertSListIntoTail(SL, buf);

        pr = loc + length;
        loc = strstr(pr, delimiter);
    }

    // Copies the last substring to the left of the last delimiter found OR
    // Copies the whole string if it doesn't have the delimiter 
    strcpy(buf, pr);
    iftInsertSListIntoTail(SL, buf);

    iftFree(buf);

    return SL;
}


char *iftSplitStringAt(const char *phrase, const char *delimiter, long position) {
    if (phrase == NULL)
        iftError("String to be splitted is NULL", "iftSplitStringAt");
    if (delimiter == NULL)
        iftError("Delimiter is NULL", "iftSplitStringAt");

    iftSList *SL = iftSplitString(phrase, delimiter);
    iftSNode *snode = NULL;

    // Copies the split sub-string of the position
    if (position >= 0) {
        if ((position+1) <= SL->n) {

            snode = SL->head;
            for (size_t i = 0; i < position; i++)
                snode = snode->next;
        }
        else {
            iftError("Invalid Position %ld\n-> Position Index must be < %ld\n",
                     "iftSplitStringAt", position, SL->n);
        }
    } else {
        if (labs(position) <= SL->n) {
            long real_pos = SL->n + position;
            snode = SL->tail;
            for (size_t i = SL->n-1; i > real_pos; i--) {
                snode = snode->prev;
            }
        }
        else {
            iftError("Invalid Negative Position %ld\n-> Negative Position Index must be >= %ld\n",
                     "iftSplitStringAt", position, -1 * SL->n);
        }
    }

    char *str = iftCopyString(snode->elem);

    iftDestroySList(&SL);

    return str;
}


char *iftUpperString(const char *str) {
    if (str == NULL)
        iftError("Input string is NULL", "iftUpperString");

    char *out_str = iftAllocCharArray(strlen(str)+1);

    for (size_t c = 0; c < strlen(str); c++)
        out_str[c] = toupper(str[c]);

    return out_str;
}


char *iftLowerString(const char *str) {
    if (str == NULL)
        iftError("Input string is NULL", "iftLowerString");

    char *out_str = iftAllocCharArray(strlen(str)+1);

    for (size_t c = 0; c < strlen(str); c++)
        out_str[c] = tolower(str[c]);

    return out_str;
}


char *iftReplaceString(const char *str, const char *old_sub, const char *new_sub) {
    if (str == NULL)
        iftError("String is NULL", "iftReplaceString");
    if (old_sub == NULL)
        iftError("Old Substring is NULL", "iftReplaceString");
    if (new_sub == NULL)
        iftError("New Substring is NULL", "iftReplaceString");

    iftSList *SL = iftSplitString(str, old_sub);

    long n_sub   = SL->n-1; // number of old subtrings found in str
    size_t str_len = strlen(str) + (n_sub * strlen(new_sub)) + 1; // size of the replaced string (with '\0')

    // adds the number of chars of the new substring + '\0'
    str_len += (n_sub * strlen(new_sub)) + 1;
    char *rep_str = iftAllocCharArray(str_len);

    // builds the replaced String - the list has always at least one element
    char *elem = iftRemoveSListHead(SL);
    while ((elem != NULL) && (SL->n >= 1)) {
        strcat(rep_str, elem);
        strcat(rep_str, new_sub);
        elem = iftRemoveSListHead(SL);
    }
    strcat(rep_str, elem); // copies the last element from the List

    iftDestroySList(&SL);

    return rep_str;
}
