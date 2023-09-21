#include "ift/core/tools/CommandLineParser.h"

#include "ift/core/io/Stream.h"
#include "ift/core/tools/Regex.h"
#include "ift/core/tools/String.h"


/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Validates the command line options are valid.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 *
 * @exception Number of options is negative.
 * @exception Command Line Option with no name defined.
 * @exception Option Short name without the prefix '-'.
 * @exception Option Short name is "-h".
 * @exception Option Long name without the prefix '--'.
 * @exception Option Long name is "--help".
 * @exception Number of arguments of an option is negative.
 */
void _iftValidateCmdLineOpts(int n_opts, const iftCmdLineOpt cmd_line_opts[]) {
    // EXCEPTION CHECKERS
    if (n_opts < 0)
        iftError("Number of command line options %d is < 0", "_iftValidateCmdLineOpts", n_opts);
    // checks each option
    for (int i = 0; i < n_opts; i++) {
        // short option name
        if (iftCompareStrings(cmd_line_opts[i].short_name, "")) {
            if (iftCompareStrings(cmd_line_opts[i].long_name, ""))
                iftError("Short and Long names from the option %d are empty \"\"",
                         "_iftValidateCmdLineOpts", i + 1);
        }
        else {
            if (iftRegexMatch(cmd_line_opts[i].short_name, "^-[a-zA-Z]([a-zA-Z]|[0-9])?$")) {
                if (iftCompareStrings(cmd_line_opts[i].short_name, "-h"))
                    iftError("The short option name \"-h\" is not allowed.\n" \
                             "It is only used to call the Usage Printing function",
                             "_iftValidateCmdLineOpts");
            }
            else
                iftError("Invalid command line short name: %s\n" \
                          "Try: -[a-zA-Z]([a-zA-Z]|[0-9])+\n\n" \
                          "Exs: -f, -a, -D, -Y, -x1", "_iftValidateCmdLineOpts", cmd_line_opts[i].short_name);
        }
        // short option name
        if (!iftCompareStrings(cmd_line_opts[i].long_name, "")) {
            if (iftRegexMatch(cmd_line_opts[i].long_name, "^--[a-zA-Z]([a-zA-Z]|[0-9]|-|_)+$")) {
                if (iftCompareStrings(cmd_line_opts[i].long_name, "--h") ||
                    iftCompareStrings(cmd_line_opts[i].long_name, "--help"))
                    iftError("The long option names \"--h\" and \"--help\" are not allowed.\n" \
                             "The areonly used to call the Usage Printing function.",
                             "iftCreateCmdLineParser");
            }
            else
                iftError("Invalid command line long name: %s\n" \
                          "Try: --[a-zA-Z]([a-zA-Z]|-|_)+\n\n" \
                          "Exs: --input-image, --outImage, --normalize, --score_file, --BOW",
                         "iftCreateCmdLineParser", cmd_line_opts[i].long_name);
        }
        // Check the Argument Type
        if (cmd_line_opts[i].has_arg) {
            if ((cmd_line_opts[i].arg_type != IFT_LONG_TYPE) &&
                (cmd_line_opts[i].arg_type != IFT_DBL_TYPE) &&
                (cmd_line_opts[i].arg_type != IFT_STR_TYPE)) {
                iftError("Invalid Argument Type. Try:\nIFT_LONG_TYPE or IFT_DBL_TYPE or IFT_STR_TYPE",
                         "iftCreateCmdLineParser");
            }
        }
    }
}


/**
 * @brief Checks if there are required options in the command line parser.
 * @author Samuel Martins
 * @date Feb 16, 2016
 * @ingroup CmdLineParse
 */
bool _iftHasRequiredOpts(const iftCmdLineParser *parser) {
    for (int i = 0; i < parser->n_opts; i++)
        if (parser->opts[i].required)
            return true;
    
    return false;
}

/**
 * @brief Checks if there are optional options in the command line parser.
 * @author Samuel Martins
 * @date Feb 16, 2016
 * @ingroup CmdLineParse
 */
bool _iftHasOptionalOpts(const iftCmdLineParser *parser) {
    for (int i = 0; i < parser->n_opts; i++)
        if (!parser->opts[i].required)
            return true;
    
    return false;
}


/**
 * @brief Checks if all required/mandatory options were passed.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 *
 * @param parser The used parser.
 * @param args   The resulting option dictionary.
 */
void _iftCheckRequiredCmdLineOpt(const iftCmdLineParser *parser, const iftDict *args) {
    for (int i = 0; i < parser->n_opts; i++) {
        if (parser->opts[i].required) {
            
            char *key = NULL;
            
            // it sets the key which will be used to find out the option in the dict
            if (!iftCompareStrings(parser->opts[i].short_name, ""))
                key = parser->opts[i].short_name;
            else
                key = parser->opts[i].long_name;
            
            
            // the required option is missing
            if (!iftDictContainKey(key, args, NULL)) {
                char missing_opt_names[1024]; // just used to print the missing option names
                
                if (!iftCompareStrings(parser->opts[i].short_name, "") &&
                    !iftCompareStrings(parser->opts[i].long_name, "")) {
                    sprintf(missing_opt_names, "\'%s\' or \'%s\'",
                            parser->opts[i].short_name, parser->opts[i].long_name);
                }
                else if (!iftCompareStrings(parser->opts[i].short_name, ""))
                    strcpy(missing_opt_names, parser->opts[i].short_name);
                else
                    strcpy(missing_opt_names, parser->opts[i].long_name);
                
                fprintf(stderr, "Required option %s is missing.\nRun \'%s -h\' or \'%s --help\' to see a full list of " \
                "available command line options", missing_opt_names, parser->program_name, parser->program_name);
                exit(EXIT_FAILURE);
            }
        }
    }
}


/**
 * @brief Get the command line option from the parser with a give (short or long) name.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 *
 * @param  cmd_line_name The (short or long) name to be found out.
 * @param  parser The command line parser.
 * @param  out_opt Reference to the cmd line option to be returned.
 * @return True if the option was found. False, otherwise.
 */
bool _iftGetCmdLineOpt(const char *cmd_line_name, iftCmdLineParser *parser, iftCmdLineOpt *out_opt) {
    for (int i = 0; i < parser->n_opts; i++)
        if (iftCompareStrings(cmd_line_name, parser->opts[i].short_name) ||
            iftCompareStrings(cmd_line_name, parser->opts[i].long_name)) {
            *out_opt = parser->opts[i];
            return true;
        }
    
    return false;
}



/********************** PUBLIC FUNCTIONS *************************/
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
iftCmdLineParser *iftCreateCmdLineParser(const char *description, int n_opts, iftCmdLineOpt cmd_line_opts[]) {
    _iftValidateCmdLineOpts(n_opts, cmd_line_opts);
    
    iftCmdLineParser *parser = (iftCmdLineParser*) iftAlloc(1, sizeof(iftCmdLineParser));
    parser->n_opts           = n_opts;
    parser->opts             = cmd_line_opts;
    parser->description      = iftCopyString(description);
    
    return parser;
}


void iftDestroyCmdLineParser(iftCmdLineParser **parser) {
    iftCmdLineParser *aux = *parser;
    
    if (aux != NULL) {
        iftFree(aux->description);
        iftFree(aux);
        *parser = NULL;
    }
}


void iftPrintUsage(const iftCmdLineParser *parser) {
    if (parser == NULL)
        iftError("Parser is NULL", "iftPrintUsage");
    
    int length_for_printing_opts = 0; // used to align the option descriptions (helps)
    for (int i = 0; i < parser->n_opts; i++) {
        int length = strlen(parser->opts[i].short_name) + strlen(parser->opts[i].long_name);
        
        if (length_for_printing_opts < length)
            length_for_printing_opts = length;
    }
    length_for_printing_opts += 8; // first 2 spaces, 1 comma between short and long names, 4 spaces after long name
    
    ////////////////////// PRINTING //////////////////////
    char opt_names[1024];
    int n_spaces = 0;
    fprintf(stdout, "--------------\n%s\n\n", parser->description);
    fprintf(stdout, "Usage\n  %s [OPTIONS]\n\n", parser->program_name);
    
    // Prints the Required Options
    if (_iftHasRequiredOpts(parser)) {
        fprintf(stdout, "Required Options\n");
        
        for (int i = 0; i < parser->n_opts; i++) {
            if (parser->opts[i].required) {
                if (!iftCompareStrings(parser->opts[i].short_name, "")) {
                    sprintf(opt_names, "  %s", parser->opts[i].short_name);
                    
                    if (!iftCompareStrings(parser->opts[i].long_name, ""))
                        sprintf(opt_names, "%s, %s", opt_names, parser->opts[i].long_name);
                }
                else sprintf(opt_names, "  %s", parser->opts[i].long_name);
                
                n_spaces = length_for_printing_opts-strlen(opt_names);
                sprintf(opt_names, "%s%*s", opt_names, n_spaces, ""); // aligning the option descriptions (helps)
                fprintf(stdout, "%s", opt_names);
                if (parser->opts[i].has_arg)
                    fprintf(stdout, "[HAS ARG] %s\n", parser->opts[i].help);
                else
                    fprintf(stdout, "%s\n", parser->opts[i].help);
            }
        }
        puts("");
    }
    
    // Prints the Optional Options
    if (_iftHasOptionalOpts(parser)) {
        fprintf(stdout, "Optional\n");
        
        for (int i = 0; i < parser->n_opts; i++) {
            if (!parser->opts[i].required) {
                if (!iftCompareStrings(parser->opts[i].short_name, "")) {
                    sprintf(opt_names, "  %s", parser->opts[i].short_name);
                    
                    if (!iftCompareStrings(parser->opts[i].long_name, ""))
                        sprintf(opt_names, "%s, %s", opt_names, parser->opts[i].long_name);
                }
                else sprintf(opt_names, "  %s", parser->opts[i].long_name);
                
                n_spaces = length_for_printing_opts-strlen(opt_names);
                sprintf(opt_names, "%s%*s", opt_names, n_spaces, ""); // aligning the option descriptions (helps)
                fprintf(stdout, "%s", opt_names);
                if (parser->opts[i].has_arg)
                    fprintf(stdout, "[HAS ARG] %s\n", parser->opts[i].help);
                else
                    fprintf(stdout, "%s\n", parser->opts[i].help);
            }
        }
        puts("");
    }
    
    // help options
    fprintf(stdout, "Help Options\n");
    sprintf(opt_names, "  -h, --help");
    n_spaces = length_for_printing_opts-strlen(opt_names);
    sprintf(opt_names, "%s%*s", opt_names, n_spaces, ""); // aligning the option descriptions (helps)
    fprintf(stdout, "%s", opt_names);
    fprintf(stdout, "%s", "Show help options\n\n");
    
    exit(EXIT_FAILURE);
    ///////////////////////////////////////////////////////
}


iftDict *iftParseCmdLine(int argc, const char *argv[], iftCmdLineParser *parser) {
    // EXCEPTION CHECKERS
    if (argv == NULL)
        iftError("argv is NULL", "iftParseCmdLine");
    if (parser == NULL)
        iftError("Command Line Parser is NULL", "iftParseCmdLine");
    
    iftDict *arg_dict = NULL;
    strcpy(parser->program_name, argv[0]); // copies the program name to the parser
    
    char default_error_msg[IFT_STR_DEFAULT_SIZE];
    sprintf(default_error_msg, "- Run \'%s -h\' or \'%s --help\' to see a full list of available " \
                               "command line options\n", parser->program_name, parser->program_name);
    
    
    // checks if some help option was passed
    for (int i = 1; i < argc; i++)
        if (iftCompareStrings(argv[i], "-h") || iftCompareStrings(argv[i], "--h") ||
            iftCompareStrings(argv[i], "--help")) {
            iftPrintUsage(parser);
        }
    
    
    // gets the number of required options
    // this will enables to use a program that contains only optional options
    int n_req_opts = 0;
    for (int i = 0; i < parser->n_opts; i++)
        n_req_opts += (parser->opts[i].required);
    
    
    // Parse the Command Line
    if ((argc == 1) && (n_req_opts != 0))
        iftPrintUsage(parser);
    else {
        // the num of expected elements to be inserted into dict is this due to the possibility of two option names
        arg_dict = iftCreateDict();
        
        // for each cmd line option and argument (it ignores the program name in argv[0])
        for (int i = 1; i < argc; i++) {
            iftCmdLineOpt opt;
            if (_iftGetCmdLineOpt(argv[i], parser, &opt)) {
                iftGVal arg = {.type=IFT_UNTYPED, .ptr_val=NULL};
                
                if (opt.has_arg) {
                    i++;
                    if (i >= argc) {
                        fprintf(stderr, "- Argument for the option \'%s\' is missing\n%s", argv[i-1], default_error_msg);
                        fflush(stdout);
                        fflush(stderr);
                        exit(EXIT_FAILURE);
                    }
                    else if (iftRegexMatch(argv[i], "^-[a-zA-Z]([a-zA-Z]|[0-9])?$") || iftRegexMatch(argv[i], "^--[a-zA-Z]([a-zA-Z]|[0-9]|-|_)+$")) {
                        fprintf(stderr, "- Invalid argument \'%s\' for the option \'%s\'\n%s", argv[i], argv[i-1], default_error_msg);
                        fflush(stdout);
                        fflush(stderr);
                        exit(EXIT_FAILURE);
                    }
                    
                    if (opt.arg_type == IFT_STR_TYPE) {
                        arg = iftCreateGVal(argv[i]);
                    }
                        // checks if it's a valid number
                    else if (iftRegexMatch(argv[i], "^-?[0-9]+(\\.[0-9]*)?$")) {
                        if (opt.arg_type == IFT_LONG_TYPE) {
                            arg = iftCreateGVal(atol(argv[i]));
                        }
                        else if (opt.arg_type == IFT_DBL_TYPE)
                            arg = iftCreateGVal(atof(argv[i]));
                        else
                            iftError("Invalid Datatype for the Argument of the option: %s",
                                     "iftParseCmdLine", argv[i - 1]);
                    }
                    else
                        iftError("Invalid Argument for the of the option: %s\nArgument: %s",
                                 "iftParseCmdLine", argv[i - 1], argv[i]);
                }
                
                if (!iftCompareStrings(opt.short_name, ""))
                    iftInsertKeyValIntoDict(iftCreateGVal(opt.short_name), iftCopyGVal(arg), arg_dict);
                if (!iftCompareStrings(opt.long_name, ""))
                    iftInsertKeyValIntoDict(iftCreateGVal(opt.long_name), iftCopyGVal(arg), arg_dict);
                
                if (arg.type == IFT_STR_TYPE)
                    iftFree(arg.str_val);
            }
            else {
                fprintf(stderr, "- Unknown option: %s\n%s", argv[i], default_error_msg);
                fflush(stdout);
                fflush(stderr);
                exit(EXIT_FAILURE);
            }
        }
        
        // Checks if all required/mandatory option were passed
        _iftCheckRequiredCmdLineOpt(parser, arg_dict);
    }
    
    return arg_dict;
}


