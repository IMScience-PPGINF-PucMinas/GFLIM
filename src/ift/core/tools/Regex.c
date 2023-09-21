#include "ift/core/tools/Regex.h"

#include <regex.h>


bool iftRegexMatch(const char *str, const char *regex_pattern, ...) {
    if (str == NULL)
        iftError("String is NULL", "iftRegexMatch");
    if (regex_pattern == NULL)
        iftError("Regular Expression is NULL", "iftRegexMatch");
    
    char error[IFT_STR_DEFAULT_SIZE];
    regex_t regex;
    int reti;
    
    va_list args;
    char final_regex_pattern[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, regex_pattern);
    vsprintf(final_regex_pattern, regex_pattern, args);
    va_end(args);
    
    // Compile Regex
    if ((reti = regcomp(&regex, final_regex_pattern, REG_EXTENDED|REG_NOSUB)) != 0) {
        regerror(reti, &regex, error, sizeof(error));
        iftError("Regex Compilation Failed: \"%s\"\n" \
                 "IFT_ERROR: %s", "iftRegexMatch", final_regex_pattern, error);
    }
    
    // Execute Regex
    reti = regexec(&regex, str, (size_t) 0, NULL, 0);
    regfree(&regex);
    
    return (reti == 0);
}
