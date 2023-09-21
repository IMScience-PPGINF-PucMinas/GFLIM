//
// Created by Samuel Martins on 08/12/18.
//

#ifndef IFT_COMMANDLINEPARSER_H
#define IFT_COMMANDLINEPARSER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/Dict.h"



/**
* @brief An Option from the Command Line.
* @author Samuel Martins
* @date Feb 15, 2016
* @ingroup CommandLineParser
* @warning At least one name must be defined. Otherwise, an error will be raised in iftCreateCmdLineParser().
*/
typedef struct ift_cmd_line_opt {
    /** Short option name. IT MUST HAVE THE PREFIX '-' */
    char short_name[4];
    /** Long option name. IT MUST HAVE THE PREFIX '--' */
    char long_name[128];
    /** Tells if the option has an argument. */
    bool has_arg;
    /** Datatype of the Argument (if the option have it). */
    iftCDataType arg_type;
    /** Tells if the option is required. */
    bool required;
    /** Description of the option (optional). It is used in the iftPrintUsage(). */
    char help[2048];
} iftCmdLineOpt;


/**
 * @brief The Command Line Parser.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 */
typedef struct ift_cmd_line_parser {
    /** Program name */
    char program_name[128];
    /** Number of options. */
    int n_opts;
    /** An array with the command line options. */
    iftCmdLineOpt *opts;
    /** Description of the program usage. It is used in the iftPrintUsage(). */
    char *description;
} iftCmdLineParser;



/**
 * @brief Creates a Command Line Parser from a set of command line options.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 *
 * @param description Description of the program usage (optional). It is used in iftPrintUsage().
 * @param Number of options.
 * @param Array of Command Line Options.
 *
 * @exception Number of options is negative.
 * @exception Command Line Option with no name defined.
 * @exception Option Short name without the prefix '-'.
 * @exception Option Short name is "-h".
 * @exception Option Long name without the prefix '--'.
 * @exception Option Long name is "--help".
 */
iftCmdLineParser *iftCreateCmdLineParser(const char *description, int n_opts, iftCmdLineOpt cmd_line_opts[]);


/**
 * @brief Destroys a Command Line Parser.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 */
void iftDestroyCmdLineParser(iftCmdLineParser **parser);


/**
 * @brief Prints the Usage Message from the command line options.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 *
 * @param parser Command Line Parser that contains all command line options.
 *
 * @note This function is called by the iftParseCmdLine().
 */
void iftPrintUsage(const iftCmdLineParser *parser);


/**
 * @brief Parses an Input Command Line.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 *
 * @param argc Size of the array <code><b>argv</b></code>.
 * @param argv String array.
 * @param parser The command line parser.
 * @return The dictionary with the passed arguments.
 *
 * @note If no options is passed or the options "-h", "--h", or "--help" are passed, the function is going to call the
 * usage printing function iftPrintUsage().
 * @note See a complete demo in demo/Miscellaneous/iftTestCmdLineParser*.c
 */
iftDict *iftParseCmdLine(int argc, const char *argv[], iftCmdLineParser *parser);



#ifdef __cplusplus
}
#endif

#endif //IFT_COMMANDLINEPARSER_H
