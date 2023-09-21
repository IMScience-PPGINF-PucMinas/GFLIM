#include "ift/core/io/NumPy.h"

#include "ift/core/tools/Regex.h"
#include "ift/core/tools/String.h"


/**
 * @brief Allocates a Numpy Array.
 * @author Samuel Martins
 */
iftNumPyHeader *_iftAllocNumPyHeader() {
    return iftAlloc(1, sizeof(iftNumPyHeader));
}


/**
 * @brief Skips a given char from an input numpy file buffer.
 * @author Samuel Martins
 */
const char *_iftSkipChar(const char *buff, const char c) {
    size_t buff_size = strlen(buff);
    size_t i = 0;
    
    for (; (buff[0] != c) && (i < buff_size); buff++, i++) {}
    
    if (i == buff_size)
        iftError("Not found char '%c' to skip", "_iftSkipChar", c);
    
    buff++; // skip the char
    
    return buff;
}


/**
 * @brief Skips qutoes from an input numpy file buffer.
 * @author Samuel Martins
 */
const char *_iftSkipQuote(const char *buff) {
    size_t buff_size = strlen(buff);
    size_t i = 0;
    
    for (size_t i = 0; (buff[0] != '\'') && (buff[0] != '\"') && (i < buff_size); buff++, i++) {}
    
    if (i == buff_size)
        iftError("Not found quote ' or \" to skip", "_iftSkipQuote");
    
    buff++; // skip the char
    
    return buff;
}


/**
 * @brief Gets a numpy file buffer after skipping some characters.
 * @author Samuel Martins
 */
char *_iftGetStrNumPyFile(const char **buff_in) {
    const char *buff = *buff_in;
    
    buff              = _iftSkipQuote(buff); // pointer in the first char after opening quote
    const char *buff2 = _iftSkipQuote(buff); // pointer in the first char after closing quote
    
    size_t str_size = strlen(buff) - strlen(buff2) - 1; // remove the closing quote
    char *str = iftAlloc(str_size+1, sizeof(char));
    memcpy(str, buff, str_size);
    
    buff += (str_size+1); // move the buffer to the next char after the closing ' or "
    *buff_in = buff2;
    
    return str;
}


/**
 * @brief Skips the key NUMPY from an input numpy file buffer.
 * @author Samuel Martins
 */
char *_iftGetKeyNumPyFile(const char **buff_in) {
    char *key = _iftGetStrNumPyFile(buff_in);
    *buff_in  = _iftSkipChar(*buff_in, ':');
    
    return key;
}


/**
 * @brief Gets the numpy dtype string from an input numpy file buffer.
 * @author Samuel Martins
 */
char *_iftGetNumpyDType(const char **buff) {
    char *key = _iftGetKeyNumPyFile(buff);
    if (!iftCompareStrings(key, "descr"))
        iftError("Invalid key: \"%s\"\nIt should be \"descr\"", "_iftReadNumPyHeader", key);
    iftFree(key);
    
    char *datatype = _iftGetStrNumPyFile(buff);
    
    if (datatype[0] == '>')
        iftError("Big-Endian not support yet", "_iftGetNumpyDType");
    else if (!iftRegexMatch(datatype, "^(<|\\|)?((i|u)(1|2|4|8)|f(4|8))$"))
        iftError("Unsupported datatype: %s\nSupported dtypes: i1, i2, i4, i8, u1, u2, u4, u8, f4, f8 (without prefix or " \
                 "with prefix < or |)", "_iftGetNumpyDType", datatype);
    
    char *dtype = NULL;
    
    // ignores the byte order symbol
    if (datatype[0] == '<' || datatype[0] == '|')
        dtype = iftCopyString(&datatype[1]);
    else dtype = iftCopyString(datatype);
    
    iftFree(datatype);
    
    return dtype;
}


/**
 * @brief Gets the Fortran Order from an input numpy file buffer.
 * @author Samuel Martins
 */
bool _iftGetFortranOrder(const char **buff_in) {
    const char *buff1 = *buff_in;
    
    char *key = _iftGetKeyNumPyFile(&buff1);
    if (!iftCompareStrings(key, "fortran_order"))
        iftError("Invalid key: \"%s\"\nIt should be \"fortran_order\"", "_iftReadNumPyHeader", key);
    iftFree(key);
    
    // skip all consecutive spaces - buff1 will point to the first char after the spaces
    while (buff1[0] == ' ')
        buff1 = _iftSkipChar(buff1, ' ');
    
    const char *buff2 = _iftSkipChar(buff1, ','); // pointer in the first char after comma
    
    size_t str_size = strlen(buff1) - strlen(buff2) - 1; // remove the comma
    char *bool_str = iftAlloc(str_size+1, sizeof(char));
    memcpy(bool_str, buff1, str_size);
    
    bool bool_val = true;
    if (iftCompareStrings(bool_str, "False"))
        bool_val = false;
    else if (iftCompareStrings(bool_str, "True"))
        bool_val = true;
    else iftError("Invalid Boolean Value: %s\n", "_iftGetFortranOrder", bool_str);
    
    iftFree(bool_str);
    
    *buff_in = buff2;
    
    return bool_val;
}


/**
 * @brief Gets the numpy shape from an input numpy file buffer.
 * @author Samuel Martins
 */
long *_iftGetShapeNumPyFile(const char **buff_in, size_t *n_dims_out) {
    const char *buff1 = *buff_in;
    
    char *key = _iftGetKeyNumPyFile(&buff1);
    if (!iftCompareStrings(key, "shape"))
        iftError("Invalid key: \"%s\"\nIt should be \"shape\"", "_iftReadNumPyHeader", key);
    iftFree(key);
    
    buff1             = _iftSkipChar(buff1, '(');
    const char *buff2 = _iftSkipChar(buff1, ')');
    
    size_t str_size = strlen(buff1) - strlen(buff2) - 1; // remove the )
    char *shape_str = iftAlloc(str_size+1, sizeof(char));
    memcpy(shape_str, buff1, str_size);
    *buff_in = buff2;
    
    // numpy allows shapes with size until 32
    // we consider that the number can be at most 36 digits
    char shape_str_numbers[32][37];
    size_t n_dims = 0;
    int begin = 0, end = 0;
    
    while (begin < strlen(shape_str)) {
        for (; shape_str[begin] == ' ' && begin < strlen(shape_str); begin++) {}
        end = begin;
        for (; shape_str[end] != ',' && end < strlen(shape_str); end++) {}
        
        int str_len = end - begin;
        // [end] is a comma
        if (str_len != 0) {
            memcpy(shape_str_numbers[n_dims], &shape_str[begin], str_len);
            shape_str_numbers[n_dims][str_len] = '\0';
            n_dims++;
            begin = end + 1; // begin is the first char after the comma comma
        }
    }
    iftFree(shape_str);
    
    long *shape = iftAlloc(n_dims, sizeof(long));
    for (int i = 0; i < n_dims; i++)
        shape[i] = atol(shape_str_numbers[i]);
    
    *n_dims_out = n_dims;
    
    // return shape;
    return shape;
}


/**
 * @brief Gets the beginning header bytes from an input numpy file.
 * @author Samuel Martins
 */
unsigned short _iftReadNumPyHeaderBeginningBytes(FILE *npy_file) {
    // The first 6 bytes are a magic string: exactly "\\x93NUMPY".
    fgetc(npy_file); // ignoring the first byte
    
    char magic[6];
    if (fread(magic, sizeof(char), 5, npy_file) != 5)
        iftError("It could not read the magic string", "_iftReadNumPyHeaderBeginningBytes");
    magic[5] = '\0';
    
    if (!iftCompareStrings(magic, "NUMPY"))
        iftError("Magic string %s is not NUMPY", "_iftReadNumPyHeaderBeginningBytes", magic);
    
    // The next 1 byte is an unsigned byte: the major version number of the file
    // format, e.g. \\x01.
    fgetc(npy_file);
    
    // The next 1 byte is an unsigned byte: the minor version number of the file
    // format, e.g. \\x00.
    fgetc(npy_file);
    
    // The next 2 bytes form a little-endian unsigned short int: the length of the
    // header data HEADER_LEN.
    unsigned short hdr_len;
    if (fread(&hdr_len, sizeof(unsigned short), 1, npy_file) != 1)
        iftError("It could not read the header length", "_iftReadNumPyHeaderBeginningBytes");
    
    return hdr_len;
}


/**
 * @brief Reads the numpy header from a numpy file.
 * @author Samuel Martins
 */
iftNumPyHeader *_iftReadNumPyHeader(FILE *npy_file) {
    unsigned short hdr_len = _iftReadNumPyHeaderBeginningBytes(npy_file);
    
    char *format = iftAlloc(hdr_len+1, sizeof(char));
    if (fread(format, sizeof(char), hdr_len, npy_file) != hdr_len)
        iftError("It could not read the array format", "_iftReadNumPyHeader");
    format[hdr_len] = '\0';
    
    const char *buff = format;
    buff             = _iftSkipChar(buff, '{');
    
    iftNumPyHeader *header = _iftAllocNumPyHeader();
    
    header->dtype         = _iftGetNumpyDType(&buff);
    header->fortran_order = _iftGetFortranOrder(&buff);
    header->shape         = _iftGetShapeNumPyFile(&buff, &header->n_dims);
    header->size = 1;
    for (int i = 0; i < header->n_dims; i++)
        header->size *= header->shape[i];
    iftFree(format);
    
    return header;
}


/**
 * @brief Writes the numpy shape into a numpy file.
 * @author Samuel Martins
 */
void _iftWriteNumPyHeader(FILE *file, const iftNumPyHeader *header) {
    char magic[7];
    magic[0] = -109; // -109 is the number stored in the first byte... why? I don't know
    strcpy(&magic[1], "NUMPY");
    fprintf(file, "%s", magic);
    
    char major_version = 1;
    char minor_version = 0;
    fwrite(&major_version, sizeof(char), 1, file);
    fwrite(&minor_version, sizeof(char), 1, file);
    
    char format[512];
    
    sprintf(format, "{'descr': '<%s', ", header->dtype);
    sprintf(format, "%s'fortran_order': %s, ", format, header->fortran_order ? "True" : "False");
    sprintf(format, "%s'shape': (%ld,", format, header->shape[0]);
    if (header->n_dims > 1) {
        for (int i = 1; i < (header->n_dims - 1); i++)
            sprintf(format, "%s %ld,", format, header->shape[i]);
        sprintf(format, "%s %ld", format, header->shape[header->n_dims - 1]);
    }
    sprintf(format, "%s)}", format);
    
    int total_hdr_len = strlen(format) + 10 + 1; // 10 = 6 (magic word) + 2 (major/minor) + 2 (header len)
    // 1 because the last header byte will be '\n'
    
    // it must be divisible by 16 for alignment purposes
    // finds the next greater number divisible by 16
    int next_divisible = 16;
    while (next_divisible < total_hdr_len)
        next_divisible += 16;
    
    // fill the header until reaching the next_divisible bytes, padding it with space,
    // and with the last header byte with \n
    int padding = next_divisible - total_hdr_len + 1;
    for (int i = 0; i < (padding-1); i++)
        strcat(format, " ");
    strcat(format, "\n");
    
    unsigned short format_len = strlen(format);
    fwrite(&format_len, sizeof(unsigned short), 1, file);
    fwrite(format, sizeof(char), format_len, file);
}


/************ PUBLIC FUNCTIONS ***************/
iftNumPyHeader *iftCreateNumPyHeader(iftCDataType cdtype, const long *shape, size_t n_dims) {
    if (shape == NULL)
        iftError("Shape is NULL", "iftCreateNumPyHeader");
    if (n_dims >= 5)
        iftError("Unsupported number of array dimensions %ld >= 5.", "iftCreateNumPyHeader", n_dims);
    
    iftNumPyHeader *header = iftAlloc(1, sizeof(iftNumPyHeader));
    header->n_dims = n_dims;
    header->shape = iftAlloc(n_dims, sizeof(long));
    header->size = 1;
    for (size_t i = 0; i < n_dims; i++) {
        header->shape[i] = shape[i];
        header->size *= header->shape[i];
    }
    header->fortran_order = false;
    header->dtype = iftCDataTypeAsNumPyDType(cdtype);

    return header;
}


void iftDestroyNumPyHeader(iftNumPyHeader **header) {
    if (header) {
        iftNumPyHeader *aux = *header;
        
        if (aux != NULL) {
            iftFree(aux->dtype);
            iftFree(aux->shape);
            iftFree(aux);
        }
    }
}


void *iftReadNumPy(const char *format, iftNumPyHeader **header_out, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, header_out);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    if (!iftEndsWith(npy_path, ".npy"))
        iftError("Invalid numpy extension for path \"%s\". Try *.npy", "iftReadNumPy", npy_path);
    
    FILE *npy_file         = fopen(npy_path, "rb");
    iftNumPyHeader *header = _iftReadNumPyHeader(npy_file);
    
    void *data = NULL;
    if (iftCompareStrings(header->dtype, "i1")) {
        data = iftAlloc(header->size, sizeof(char));
        if (fread(data, sizeof(char), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else if (iftCompareStrings(header->dtype, "i2")) {
        data = iftAlloc(header->size, sizeof(short));
        if (fread(data, sizeof(short), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else if (iftCompareStrings(header->dtype, "i4")) {
        data = iftAlloc(header->size, sizeof(int));
        if (fread(data, sizeof(int), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else if (iftCompareStrings(header->dtype, "i8")) {
        data = iftAlloc(header->size, sizeof(long));
        if (fread(data, sizeof(long), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else if (iftCompareStrings(header->dtype, "u1")) {
        data = iftAlloc(header->size, sizeof(uchar));
        if (fread(data, sizeof(uchar), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else if (iftCompareStrings(header->dtype, "u2")) {
        data = iftAlloc(header->size, sizeof(ushort));
        if (fread(data, sizeof(ushort), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else if (iftCompareStrings(header->dtype, "u4")) {
        data = iftAlloc(header->size, sizeof(uint));
        if (fread(data, sizeof(uint), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else if (iftCompareStrings(header->dtype, "u8")) {
        data = iftAlloc(header->size, sizeof(ulong));
        if (fread(data, sizeof(ulong), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else if (iftCompareStrings(header->dtype, "f4")) {
        data = iftAlloc(header->size, sizeof(float));
        if (fread(data, sizeof(float), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else if (iftCompareStrings(header->dtype, "f8")) {
        data = iftAlloc(header->size, sizeof(double));
        if (fread(data, sizeof(double), header->size, npy_file) != header->size)
            iftError("Error when reading data", "iftReadNumPy");
    }
    else iftError("Dataype \"%s\" not supported yet", "iftReadNumPy", header->dtype);
    
    fclose(npy_file);
    
    *header_out = header;
    
    return data;
}


void iftWriteNumPy(const iftNumPyHeader *header, const void *data, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);

    iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);
    size_t cdtype_size = iftCDataTypeSize(cdtype);
    
    FILE *file = fopen(npy_path, "wb");
    _iftWriteNumPyHeader(file, header);
    iftWriteDataInFile(file, data, header->size, cdtype_size);
    fclose(file);
}


iftCDataType iftNumPyDTypeAsCDataType(const char *dtype) {
    if (iftCompareStrings(dtype, "i1"))
        return(IFT_CHAR_TYPE);
    else if (iftCompareStrings(dtype, "i2"))
        return(IFT_SHORT_TYPE);
    else if (iftCompareStrings(dtype, "i4"))
        return(IFT_INT_TYPE);
    else if (iftCompareStrings(dtype, "i8"))
        return(IFT_LONG_TYPE);
    else if (iftCompareStrings(dtype, "u1"))
        return(IFT_UCHAR_TYPE);
    else if (iftCompareStrings(dtype, "u2"))
        return(IFT_USHORT_TYPE);
    else if (iftCompareStrings(dtype, "u4"))
        return(IFT_UINT_TYPE);
    else if (iftCompareStrings(dtype, "u8"))
        return(IFT_ULONG_TYPE);
    else if (iftCompareStrings(dtype, "f4"))
        return(IFT_FLT_TYPE);
    else if (iftCompareStrings(dtype, "f8"))
        return(IFT_DBL_TYPE);
    else {
        iftWarning("Unsupported numpy datatype %s", "iftNumPyDTypeAsCDataType", dtype);
        return(IFT_UNTYPED);
    }
}


char *iftCDataTypeAsNumPyDType(iftCDataType cdtype) {
    switch (cdtype) {
        case IFT_CHAR_TYPE:
            return iftCopyString("i1");
        case IFT_SHORT_TYPE:
            return iftCopyString("i2");
        case IFT_INT_TYPE:
            return iftCopyString("i4");
        case IFT_LONG_TYPE:
            return iftCopyString("i8");
        case IFT_UCHAR_TYPE:
            return iftCopyString("u1");
        case IFT_USHORT_TYPE:
            return iftCopyString("u2");
        case IFT_UINT_TYPE:
            return iftCopyString("u4");
        case IFT_ULONG_TYPE:
            return iftCopyString("u8");
        case IFT_FLT_TYPE:
            return iftCopyString("f4");
        case IFT_DBL_TYPE:
            return iftCopyString("f8");
        default:
            iftError("Unsupported CDataType to NumPy dtype", "iftCDataTypeAsNumPyDType");
            return NULL;
    }
}

