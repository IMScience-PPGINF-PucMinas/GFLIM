#include "ift/core/io/Json.h"

#include "ift/core/dtypes/FloatArray.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/IntMatrix.h"
#include "ift/core/dtypes/StrArray.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftMatrix.h"


/********************** cJSON Sources *************************/
/**
 * @brief Json Data Types.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 */
typedef enum {
    IFT_JSON_FALSE,
    IFT_JSON_TRUE,
    IFT_JSON_NULL,
    IFT_JSON_INT,
    IFT_JSON_DBL,
    IFT_JSON_STRING,
    IFT_JSON_ARRAY,
    IFT_JSON_DICT
} iftJsonType;


typedef struct cJSON_Hooks {
    void *(*malloc_fn)(size_t sz);
    
    void (*free_fn)(void *ptr);
} cJSON_Hooks;

#define cJSON_IsReference 256
#define cJSON_StringIsConst 512

/* Macros for creating things quickly. */
#define cJSON_AddNullToObject(object, name)        cJSON_AddItemToObject(object, name, cJSON_CreateNull())
#define cJSON_AddTrueToObject(object, name)        cJSON_AddItemToObject(object, name, cJSON_CreateTrue())
#define cJSON_AddFalseToObject(object, name)        cJSON_AddItemToObject(object, name, cJSON_CreateFalse())
#define cJSON_AddBoolToObject(object, name, b)    cJSON_AddItemToObject(object, name, cJSON_CreateBool(b))
#define cJSON_AddNumberToObject(object, name, n)    cJSON_AddItemToObject(object, name, cJSON_CreateNumber(n))
#define cJSON_AddStringToObject(object, name, s)    cJSON_AddItemToObject(object, name, cJSON_CreateString(s))

/* When assigning an integer value, it needs to be propagated to dbl_val too. */
#define cJSON_SetIntValue(object, val)            ((object)?(object)->valueint=(object)->valuedouble=(val):(val))
#define cJSON_SetNumberValue(object, val)        ((object)?(object)->valueint=(object)->valuedouble=(val):(val))


static const char *ep;

const char *cJSON_GetErrorPtr(void) { return ep; }

static int cJSON_strcasecmp(const char *s1, const char *s2) {
    if (!s1) return (s1 == s2) ? 0 : 1;
    if (!s2) return 1;
    for (; tolower(*s1) == tolower(*s2); ++s1, ++s2) if (*s1 == 0) return 0;
    return tolower(*(const unsigned char *) s1) - tolower(*(const unsigned char *) s2);
}

static void *(*cJSON_malloc)(size_t sz) = malloc;

static void (*cJSON_free)(void *ptr) = free;

static char *cJSON_strdup(const char *str) {
    size_t len;
    char *copy;
    
    len = strlen(str) + 1;
    if (!(copy = (char *) cJSON_malloc(len))) return 0;
    memcpy(copy, str, len);
    return copy;
}

void cJSON_InitHooks(cJSON_Hooks *hooks) {
    if (!hooks) { /* Reset hooks */
        cJSON_malloc = malloc;
        cJSON_free = free;
        return;
    }
    
    cJSON_malloc = (hooks->malloc_fn) ? hooks->malloc_fn : malloc;
    cJSON_free = (hooks->free_fn) ? hooks->free_fn : free;
}

/* Internal constructor. */
static iftJson *cJSON_New_Item(void) {
    iftJson *node = (iftJson *) cJSON_malloc(sizeof(iftJson));
    if (node) memset(node, 0, sizeof(iftJson));
    return node;
}

/* Delete a cJSON structure. */
void cJSON_Delete(iftJson *c) {
    iftJson *next;
    while (c) {
        next = c->next;
        if (!(c->type & cJSON_IsReference) && c->first_child) cJSON_Delete(c->first_child);
        if (!(c->type & cJSON_IsReference) && c->str_val) cJSON_free(c->str_val);
        if (!(c->type & cJSON_StringIsConst) && c->key) cJSON_free(c->key);
        cJSON_free(c);
        c = next;
    }
}

/* Parse the input text to generate a number, and populate the result into item. */
static const char *parse_number(iftJson *item, const char *num) {
    iftJsonType jtype = IFT_JSON_INT;
    
    double n = 0, sign = 1, scale = 0;
    int subscale = 0, signsubscale = 1;
    
    if (*num == '-') sign = -1, num++;    /* Has sign? */
    if (*num == '0') num++;            /* is zero */
    if (*num >= '1' && *num <= '9')
        do n = (n * 10.0) + (*num++ - '0'); while (*num >= '0' && *num <= '9');    /* Number? */
    if (*num == '.' && num[1] >= '0' && num[1] <= '9') {
        jtype = IFT_JSON_DBL;
        num++;
        do n = (n * 10.0) + (*num++ - '0'), scale--; while (*num >= '0' && *num <= '9');
    }    /* Fractional part? */
    if (*num == 'e' || *num == 'E')        /* Exponent? */
    {
        num++;
        if (*num == '+') num++; else if (*num == '-') signsubscale = -1, num++;        /* With sign? */
        while (*num >= '0' && *num <= '9') subscale = (subscale * 10) + (*num++ - '0');    /* Number? */
    }
    
    n = sign * n * pow(10.0, (scale + subscale * signsubscale));    /* number = +/- number.fraction * 10^+/- exponent */
    
    item->dbl_val = n;
    item->int_val = n;
    item->type = jtype;
    return num;
}

static int pow2gt(int x) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x + 1;
}

typedef struct {
    char *buffer;
    int length;
    int offset;
} printbuffer;

static char *ensure(printbuffer *p, int needed) {
    char *newbuffer;
    int newsize;
    if (!p || !p->buffer) return 0;
    needed += p->offset;
    if (needed <= p->length) return p->buffer + p->offset;
    
    newsize = pow2gt(needed);
    newbuffer = (char *) cJSON_malloc(newsize);
    if (!newbuffer) {
        cJSON_free(p->buffer);
        p->length = 0, p->buffer = 0;
        return 0;
    }
    if (newbuffer) memcpy(newbuffer, p->buffer, p->length);
    cJSON_free(p->buffer);
    p->length = newsize;
    p->buffer = newbuffer;
    return newbuffer + p->offset;
}

static int update(printbuffer *p) {
    char *str;
    if (!p || !p->buffer) return 0;
    str = p->buffer + p->offset;
    return p->offset + strlen(str);
}

/* Render the number nicely from the given item into a key. */
static char *print_number(iftJson *item, printbuffer *p) {
    char *str = 0;
    double d = item->dbl_val;
    if (d == 0) {
        if (p) str = ensure(p, 2);
        else str = (char *) cJSON_malloc(2);    /* special case for 0. */
        if (str) strcpy(str, "0");
    }
    else if (fabs(((double) item->int_val) - d) <= DBL_EPSILON && d <= INT_MAX && d >= INT_MIN) {
        if (p) str = ensure(p, 21);
        else str = (char *) cJSON_malloc(21);    /* 2^64+1 can be represented in 21 chars. */
        if (str) sprintf(str, "%d", item->int_val);
    }
    else {
        if (p) str = ensure(p, 64);
        else str = (char *) cJSON_malloc(64);    /* This is a nice tradeoff. */
        if (str) {
            if (fabs(floor(d) - d) <= DBL_EPSILON && fabs(d) < 1.0e60)sprintf(str, "%.0f", d);
            else if (fabs(d) < 1.0e-6 || fabs(d) > 1.0e9) sprintf(str, "%e", d);
            else sprintf(str, "%f", d);
        }
    }
    return str;
}

static unsigned parse_hex4(const char *str) {
    unsigned h = 0;
    if (*str >= '0' && *str <= '9') h += (*str) - '0'; else if (*str >= 'A' && *str <= 'F')
        h += 10 + (*str) - 'A';
    else if (*str >= 'a' && *str <= 'f') h += 10 + (*str) - 'a'; else return 0;
    h = h << 4;
    str++;
    if (*str >= '0' && *str <= '9') h += (*str) - '0'; else if (*str >= 'A' && *str <= 'F')
        h += 10 + (*str) - 'A';
    else if (*str >= 'a' && *str <= 'f') h += 10 + (*str) - 'a'; else return 0;
    h = h << 4;
    str++;
    if (*str >= '0' && *str <= '9') h += (*str) - '0'; else if (*str >= 'A' && *str <= 'F')
        h += 10 + (*str) - 'A';
    else if (*str >= 'a' && *str <= 'f') h += 10 + (*str) - 'a'; else return 0;
    h = h << 4;
    str++;
    if (*str >= '0' && *str <= '9') h += (*str) - '0'; else if (*str >= 'A' && *str <= 'F')
        h += 10 + (*str) - 'A';
    else if (*str >= 'a' && *str <= 'f') h += 10 + (*str) - 'a'; else return 0;
    return h;
}

/* Parse the input text into an unescaped cstring, and populate item. */
static const unsigned char firstByteMark[7] = {0x00, 0x00, 0xC0, 0xE0, 0xF0, 0xF8, 0xFC};

static const char *parse_string(iftJson *item, const char *str) {
    const char *ptr = str + 1;
    char *ptr2;
    char *out;
    int len = 0;
    unsigned uc, uc2;
    if (*str != '\"') {
        ep = str;
        return 0;
    }    /* not a key! */
    
    while (*ptr != '\"' && *ptr && ++len) if (*ptr++ == '\\') ptr++;    /* Skip escaped quotes. */
    
    out = (char *) cJSON_malloc(len + 1);    /* This is how long we need for the key, roughly. */
    if (!out) return 0;
    
    ptr = str + 1;
    ptr2 = out;
    while (*ptr != '\"' && *ptr) {
        if (*ptr != '\\') *ptr2++ = *ptr++;
        else {
            ptr++;
            switch (*ptr) {
                case 'b':
                    *ptr2++ = '\b';
                    break;
                case 'f':
                    *ptr2++ = '\f';
                    break;
                case 'n':
                    *ptr2++ = '\n';
                    break;
                case 'r':
                    *ptr2++ = '\r';
                    break;
                case 't':
                    *ptr2++ = '\t';
                    break;
                case 'u':     /* transcode utf16 to utf8. */
                    uc = parse_hex4(ptr + 1);
                    ptr += 4;    /* get the unicode char. */
                    
                    if ((uc >= 0xDC00 && uc <= 0xDFFF) || uc == 0) break;    /* check for invalid.  */
                    
                    if (uc >= 0xD800 && uc <= 0xDBFF)    /* UTF16 surrogate pairs.  */
                    {
                        if (ptr[1] != '\\' || ptr[2] != 'u') break;    /* missing second-half of surrogate. */
                        uc2 = parse_hex4(ptr + 3);
                        ptr += 6;
                        if (uc2 < 0xDC00 || uc2 > 0xDFFF) break;    /* invalid second-half of surrogate.    */
                        uc = 0x10000 + (((uc & 0x3FF) << 10) | (uc2 & 0x3FF));
                    }
                    
                    len = 4;
                    if (uc < 0x80) len = 1; else if (uc < 0x800) len = 2; else if (uc < 0x10000) len = 3;
                    ptr2 += len;
                    
                    switch (len) {
                        case 4:
                            *--ptr2 = ((uc | 0x80) & 0xBF);
                            uc >>= 6;
                        case 3:
                            *--ptr2 = ((uc | 0x80) & 0xBF);
                            uc >>= 6;
                        case 2:
                            *--ptr2 = ((uc | 0x80) & 0xBF);
                            uc >>= 6;
                        case 1:
                            *--ptr2 = (uc | firstByteMark[len]);
                    }
                    ptr2 += len;
                    break;
                default:
                    *ptr2++ = *ptr;
                    break;
            }
            ptr++;
        }
    }
    *ptr2 = 0;
    if (*ptr == '\"') ptr++;
    item->str_val = out;
    item->type = IFT_JSON_STRING;
    return ptr;
}

/* Render the cstring provided to an escaped version that can be printed. */
static char *print_string_ptr(const char *str, printbuffer *p) {
    const char *ptr;
    char *ptr2, *out;
    int len = 0, flag = 0;
    unsigned char token;
    
    for (ptr = str; *ptr; ptr++) flag |= ((*ptr > 0 && *ptr < 32) || (*ptr == '\"') || (*ptr == '\\')) ? 1 : 0;
    if (!flag) {
        len = ptr - str;
        if (p) out = ensure(p, len + 3);
        else out = (char *) cJSON_malloc(len + 3);
        if (!out) return 0;
        ptr2 = out;
        *ptr2++ = '\"';
        strcpy(ptr2, str);
        ptr2[len] = '\"';
        ptr2[len + 1] = 0;
        return out;
    }
    
    if (!str) {
        if (p) out = ensure(p, 3);
        else out = (char *) cJSON_malloc(3);
        if (!out) return 0;
        strcpy(out, "\"\"");
        return out;
    }
    ptr = str;
    while ((token = *ptr) && ++len) {
        if (strchr("\"\\\b\f\n\r\t", token)) len++; else if (token < 32) len += 5;
        ptr++;
    }
    
    if (p) out = ensure(p, len + 3);
    else out = (char *) cJSON_malloc(len + 3);
    if (!out) return 0;
    
    ptr2 = out;
    ptr = str;
    *ptr2++ = '\"';
    while (*ptr) {
        if ((unsigned char) *ptr > 31 && *ptr != '\"' && *ptr != '\\') *ptr2++ = *ptr++;
        else {
            *ptr2++ = '\\';
            switch (token = *ptr++) {
                case '\\':
                    *ptr2++ = '\\';
                    break;
                case '\"':
                    *ptr2++ = '\"';
                    break;
                case '\b':
                    *ptr2++ = 'b';
                    break;
                case '\f':
                    *ptr2++ = 'f';
                    break;
                case '\n':
                    *ptr2++ = 'n';
                    break;
                case '\r':
                    *ptr2++ = 'r';
                    break;
                case '\t':
                    *ptr2++ = 't';
                    break;
                default:
                    sprintf(ptr2, "u%04x", token);
                    ptr2 += 5;
                    break;    /* escape and print */
            }
        }
    }
    *ptr2++ = '\"';
    *ptr2++ = 0;
    return out;
}

/* Invote print_string_ptr (which is useful) on an item. */
static char *print_string(iftJson *item, printbuffer *p) { return print_string_ptr(item->str_val, p); }

/* Predeclare these prototypes. */
static const char *parse_value(iftJson *item, const char *value);

static char *print_value(iftJson *item, int depth, int fmt, printbuffer *p);

static const char *parse_array(iftJson *item, const char *value);

static char *print_array(iftJson *item, int depth, int fmt, printbuffer *p);

static const char *parse_object(iftJson *item, const char *value);

static char *print_object(iftJson *item, int depth, int fmt, printbuffer *p);

/* Utility to jump whitespace and cr/lf */
static const char *skip(const char *in) {
    while (in && *in && (unsigned char) *in <= 32) in++;
    return in;
}

/* Parse an object - create a new root, and populate. */
iftJson *cJSON_ParseWithOpts(const char *value, const char **return_parse_end, int require_null_terminated) {
    const char *end = 0;
    iftJson *c = cJSON_New_Item();
    ep = 0;
    if (!c) return 0;       /* memory fail */
    
    end = parse_value(c, skip(value));
    if (!end) {
        cJSON_Delete(c);
        return 0;
    }    /* parse failure. ep is set. */
    
    /* if we require null-terminated JSON without appended garbage, skip and then check for a null terminator */
    if (require_null_terminated) {
        end = skip(end);
        if (*end) {
            cJSON_Delete(c);
            ep = end;
            return 0;
        }
    }
    if (return_parse_end) *return_parse_end = end;
    return c;
}


/* Render a cJSON item/entity/structure to text. */
char *cJSON_Print(iftJson *item) { return print_value(item, 0, 1, 0); }

char *cJSON_PrintUnformatted(iftJson *item) { return print_value(item, 0, 0, 0); }

char *cJSON_PrintBuffered(iftJson *item, int prebuffer, int fmt) {
    printbuffer p;
    p.buffer = (char *) cJSON_malloc(prebuffer);
    p.length = prebuffer;
    p.offset = 0;
    return print_value(item, 0, fmt, &p);
    return p.buffer;
}


/* Parser core - when encountering text, process appropriately. */
static const char *parse_value(iftJson *item, const char *value) {
    if (!value) return 0;    /* Fail on null. */
    if (!strncmp(value, "null", 4)) {
        item->type = IFT_JSON_NULL;
        return value + 4;
    }
    if (!strncmp(value, "false", 5)) {
        item->type = IFT_JSON_FALSE;
        return value + 5;
    }
    if (!strncmp(value, "true", 4)) {
        item->type = IFT_JSON_TRUE;
        item->int_val = 1;
        return value + 4;
    }
    if (*value == '\"') { return parse_string(item, value); }
    if (*value == '-' || (*value >= '0' && *value <= '9')) { return parse_number(item, value); }
    if (*value == '[') { return parse_array(item, value); }
    if (*value == '{') { return parse_object(item, value); }
    
    ep = value;
    return 0;    /* failure. */
}

/* Render a value to text. */
static char *print_value(iftJson *item, int depth, int fmt, printbuffer *p) {
    char *out = 0;
    if (!item) return 0;
    if (p) {
        switch ((item->type) & 255) {
            case IFT_JSON_NULL: {
                out = ensure(p, 5);
                if (out) strcpy(out, "null");
                break;
            }
            case IFT_JSON_FALSE: {
                out = ensure(p, 6);
                if (out) strcpy(out, "false");
                break;
            }
            case IFT_JSON_TRUE: {
                out = ensure(p, 5);
                if (out) strcpy(out, "true");
                break;
            }
            case IFT_JSON_INT:
            case IFT_JSON_DBL:
                out = print_number(item, p);
                break;
            case IFT_JSON_STRING:
                out = print_string(item, p);
                break;
            case IFT_JSON_ARRAY:
                out = print_array(item, depth, fmt, p);
                break;
            case IFT_JSON_DICT:
                out = print_object(item, depth, fmt, p);
                break;
        }
    }
    else {
        switch ((item->type) & 255) {
            case IFT_JSON_NULL:
                out = cJSON_strdup("null");
                break;
            case IFT_JSON_FALSE:
                out = cJSON_strdup("false");
                break;
            case IFT_JSON_TRUE:
                out = cJSON_strdup("true");
                break;
            case IFT_JSON_INT:
            case IFT_JSON_DBL:
                out = print_number(item, 0);
                break;
            case IFT_JSON_STRING:
                out = print_string(item, 0);
                break;
            case IFT_JSON_ARRAY:
                out = print_array(item, depth, fmt, 0);
                break;
            case IFT_JSON_DICT:
                out = print_object(item, depth, fmt, 0);
                break;
        }
    }
    return out;
}

/* Build an array from input text. */
static const char *parse_array(iftJson *item, const char *value) {
    iftJson *aux;
    if (*value != '[') {
        ep = value;
        return 0;
    }    /* not an array! */
    
    item->type = IFT_JSON_ARRAY;
    value = skip(value + 1);
    if (*value == ']') return value + 1;    /* empty array. */
    
    item->first_child = aux = cJSON_New_Item();
    if (!item->first_child) return 0;         /* memory fail */
    value = skip(parse_value(aux, skip(value)));    /* skip any spacing, get the value. */
    if (!value) return 0;
    
    while (*value == ',') {
        iftJson *new_item;
        if (!(new_item = cJSON_New_Item())) return 0;    /* memory fail */
        aux->next = new_item;
        new_item->prev = aux;
        aux = new_item;
        value = skip(parse_value(aux, skip(value + 1)));
        if (!value) return 0;    /* memory fail */
    }
    
    /* end of array */
    if (*value == ']') {
        item->last_child = aux; // the variable aux points to the last element of the array
        return value + 1;
    }
    ep = value;
    return 0;    /* malformed. */
}

/* Render an array to text */
static char *print_array(iftJson *item, int depth, int fmt, printbuffer *p) {
    char **entries;
    char *out = 0, *ptr, *ret;
    int len = 5;
    iftJson *child = item->first_child;
    int numentries = 0, i = 0, fail = 0;
    size_t tmplen = 0;
    
    /* How many entries in the array? */
    while (child) numentries++, child = child->next;
    /* Explicitly handle numentries==0 */
    if (!numentries) {
        if (p) out = ensure(p, 3);
        else out = (char *) cJSON_malloc(3);
        if (out) strcpy(out, "[]");
        return out;
    }
    
    if (p) {
        /* Compose the output array. */
        i = p->offset;
        ptr = ensure(p, 1);
        if (!ptr) return 0;
        *ptr = '[';
        p->offset++;
        child = item->first_child;
        while (child && !fail) {
            print_value(child, depth + 1, fmt, p);
            p->offset = update(p);
            if (child->next) {
                len = fmt ? 2 : 1;
                ptr = ensure(p, len + 1);
                if (!ptr) return 0;
                *ptr++ = ',';
                if (fmt)*ptr++ = ' ';
                *ptr = 0;
                p->offset += len;
            }
            child = child->next;
        }
        ptr = ensure(p, 2);
        if (!ptr) return 0;
        *ptr++ = ']';
        *ptr = 0;
        out = (p->buffer) + i;
    }
    else {
        /* Allocate an array to hold the values for each */
        entries = (char **) cJSON_malloc(numentries * sizeof(char *));
        if (!entries) return 0;
        memset(entries, 0, numentries * sizeof(char *));
        /* Retrieve all the results: */
        child = item->first_child;
        while (child && !fail) {
            ret = print_value(child, depth + 1, fmt, 0);
            entries[i++] = ret;
            if (ret) len += strlen(ret) + 2 + (fmt ? 1 : 0); else fail = 1;
            child = child->next;
        }
        
        /* If we didn't fail, try to malloc the output key */
        if (!fail) out = (char *) cJSON_malloc(len);
        /* If that fails, we fail. */
        if (!out) fail = 1;
        
        /* Handle failure. */
        if (fail) {
            for (i = 0; i < numentries; i++) if (entries[i]) cJSON_free(entries[i]);
            cJSON_free(entries);
            return 0;
        }
        
        /* Compose the output array. */
        *out = '[';
        ptr = out + 1;
        *ptr = 0;
        for (i = 0; i < numentries; i++) {
            tmplen = strlen(entries[i]);
            memcpy(ptr, entries[i], tmplen);
            ptr += tmplen;
            if (i != numentries - 1) {
                *ptr++ = ',';
                if (fmt)*ptr++ = ' ';
                *ptr = 0;
            }
            cJSON_free(entries[i]);
        }
        cJSON_free(entries);
        *ptr++ = ']';
        *ptr++ = 0;
    }
    return out;
}

/* Build an object from the text. */
static const char *parse_object(iftJson *item, const char *value) {
    iftJson *aux;
    if (*value != '{') {
        ep = value;
        return 0;
    }    /* not an object! */
    
    item->type = IFT_JSON_DICT;
    value = skip(value + 1);
    if (*value == '}') return value + 1;    /* empty array. */
    
    item->first_child = aux = cJSON_New_Item();
    if (!item->first_child) return 0;
    value = skip(parse_string(aux, skip(value)));
    if (!value) return 0;
    aux->key = aux->str_val;
    aux->str_val = 0;
    if (*value != ':') {
        ep = value;
        return 0;
    }    /* fail! */
    value = skip(parse_value(aux, skip(value + 1)));    /* skip any spacing, get the value. */
    if (!value) return 0;
    
    while (*value == ',') {
        iftJson *new_item;
        if (!(new_item = cJSON_New_Item())) return 0; /* memory fail */
        aux->next = new_item;
        new_item->prev = aux;
        aux = new_item;
        value = skip(parse_string(aux, skip(value + 1)));
        if (!value) return 0;
        aux->key = aux->str_val;
        aux->str_val = 0;
        if (*value != ':') {
            ep = value;
            return 0;
        }    /* fail! */
        value = skip(parse_value(aux, skip(value + 1)));    /* skip any spacing, get the value. */
        if (!value) return 0;
    }
    
    /* end of array */
    if (*value == '}') {
        item->last_child = aux; // the variable aux points to the last element of the dictionary
        return value + 1;
    }
    ep = value;
    return 0;    /* malformed. */
}

/* Render an object to text. */
static char *print_object(iftJson *item, int depth, int fmt, printbuffer *p) {
    char **entries = 0, **names = 0;
    char *out = 0, *ptr, *ret, *str;
    int len = 7, i = 0, j;
    iftJson *child = item->first_child;
    int numentries = 0, fail = 0;
    size_t tmplen = 0;
    /* Count the number of entries. */
    while (child) numentries++, child = child->next;
    /* Explicitly handle empty object case */
    if (!numentries) {
        if (p) out = ensure(p, fmt ? depth + 4 : 3);
        else out = (char *) cJSON_malloc(fmt ? depth + 4 : 3);
        if (!out) return 0;
        ptr = out;
        *ptr++ = '{';
        if (fmt) {
            *ptr++ = '\n';
            for (i = 0; i < depth - 1; i++) *ptr++ = '\t';
        }
        *ptr++ = '}';
        *ptr++ = 0;
        return out;
    }
    if (p) {
        /* Compose the output: */
        i = p->offset;
        len = fmt ? 2 : 1;
        ptr = ensure(p, len + 1);
        if (!ptr) return 0;
        *ptr++ = '{';
        if (fmt) *ptr++ = '\n';
        *ptr = 0;
        p->offset += len;
        child = item->first_child;
        depth++;
        while (child) {
            if (fmt) {
                ptr = ensure(p, depth);
                if (!ptr) return 0;
                for (j = 0; j < depth; j++) *ptr++ = '\t';
                p->offset += depth;
            }
            print_string_ptr(child->key, p);
            p->offset = update(p);
            
            len = fmt ? 2 : 1;
            ptr = ensure(p, len);
            if (!ptr) return 0;
            *ptr++ = ':';
            if (fmt) *ptr++ = '\t';
            p->offset += len;
            
            print_value(child, depth, fmt, p);
            p->offset = update(p);
            
            len = (fmt ? 1 : 0) + (child->next ? 1 : 0);
            ptr = ensure(p, len + 1);
            if (!ptr) return 0;
            if (child->next) *ptr++ = ',';
            if (fmt) *ptr++ = '\n';
            *ptr = 0;
            p->offset += len;
            child = child->next;
        }
        ptr = ensure(p, fmt ? (depth + 1) : 2);
        if (!ptr) return 0;
        if (fmt) for (i = 0; i < depth - 1; i++) *ptr++ = '\t';
        *ptr++ = '}';
        *ptr = 0;
        out = (p->buffer) + i;
    }
    else {
        /* Allocate space for the names and the objects */
        entries = (char **) cJSON_malloc(numentries * sizeof(char *));
        if (!entries) return 0;
        names = (char **) cJSON_malloc(numentries * sizeof(char *));
        if (!names) {
            cJSON_free(entries);
            return 0;
        }
        memset(entries, 0, sizeof(char *) * numentries);
        memset(names, 0, sizeof(char *) * numentries);
        
        /* Collect all the results into our arrays: */
        child = item->first_child;
        depth++;
        if (fmt) len += depth;
        while (child) {
            names[i] = str = print_string_ptr(child->key, 0);
            entries[i++] = ret = print_value(child, depth, fmt, 0);
            if (str && ret) len += strlen(ret) + strlen(str) + 2 + (fmt ? 2 + depth : 0); else fail = 1;
            child = child->next;
        }
        
        /* Try to allocate the output key */
        if (!fail) out = (char *) cJSON_malloc(len);
        if (!out) fail = 1;
        
        /* Handle failure */
        if (fail) {
            for (i = 0; i < numentries; i++) {
                if (names[i]) cJSON_free(names[i]);
                if (entries[i]) cJSON_free(entries[i]);
            }
            cJSON_free(names);
            cJSON_free(entries);
            return 0;
        }
        
        /* Compose the output: */
        *out = '{';
        ptr = out + 1;
        if (fmt)*ptr++ = '\n';
        *ptr = 0;
        for (i = 0; i < numentries; i++) {
            if (fmt) for (j = 0; j < depth; j++) *ptr++ = '\t';
            tmplen = strlen(names[i]);
            memcpy(ptr, names[i], tmplen);
            ptr += tmplen;
            *ptr++ = ':';
            if (fmt) *ptr++ = ' ';
            strcpy(ptr, entries[i]);
            ptr += strlen(entries[i]);
            if (i != numentries - 1) *ptr++ = ',';
            if (fmt) *ptr++ = '\n';
            *ptr = 0;
            cJSON_free(names[i]);
            cJSON_free(entries[i]);
        }
        
        cJSON_free(names);
        cJSON_free(entries);
        if (fmt) for (i = 0; i < depth - 1; i++) *ptr++ = '\t';
        *ptr++ = '}';
        *ptr++ = 0;
    }
    return out;
}

/* Get Array size/item / object item. */
int cJSON_GetArraySize(iftJson *array) {
    iftJson *c = array->first_child;
    int i = 0;
    while (c)i++, c = c->next;
    return i;
}

iftJson *cJSON_GetArrayItem(iftJson *array, int item) {
    iftJson *c = array->first_child;
    while (c && item > 0) item--, c = c->next;
    return c;
}


iftJson *cJSON_GetObjectItem(iftJson *object, const char *string) {
    iftJson *c = object->first_child;
    while (c && cJSON_strcasecmp(c->key, string))
        c = c->next;
    
    return c;
}

/* Utility for array list handling. */
static void suffix_object(iftJson *prev, iftJson *item) {
    prev->next = item;
    item->prev = prev;
}

/* Utility for handling references. */
static iftJson *create_reference(iftJson *item) {
    iftJson *ref = cJSON_New_Item();
    if (!ref) return 0;
    memcpy(ref, item, sizeof(iftJson));
    ref->key = 0;
    ref->type |= cJSON_IsReference;
    ref->next = ref->prev = 0;
    return ref;
}

/* Add item to array/object. */
void cJSON_AddItemToArray(iftJson *array, iftJson *item) {
    iftJson *c = array->first_child;
    if (!item) return;
    if (!c) { array->first_child = item; } else {
        while (c && c->next) c = c->next;
        suffix_object(c, item);
    }
}

void cJSON_AddItemToObject(iftJson *object, const char *string, iftJson *item) {
    if (!item) return;
    if (item->key) cJSON_free(item->key);
    item->key = cJSON_strdup(string);
    cJSON_AddItemToArray(object, item);
}

void cJSON_AddItemToObjectCS(iftJson *object, const char *string, iftJson *item) {
    if (!item) return;
    if (!(item->type & cJSON_StringIsConst) && item->key) cJSON_free(item->key);
    item->key = (char *) string;
    item->type |= cJSON_StringIsConst;
    cJSON_AddItemToArray(object, item);
}

void cJSON_AddItemReferenceToArray(iftJson *array, iftJson *item) {
    cJSON_AddItemToArray(array, create_reference(item));
}

void cJSON_AddItemReferenceToObject(iftJson *object, const char *string, iftJson *item) {
    cJSON_AddItemToObject(object, string, create_reference(item));
}

iftJson *cJSON_DetachItemFromArray(iftJson *array, int which) {
    iftJson *c = array->first_child;
    while (c && which > 0) c = c->next, which--;
    if (!c) return 0;
    if (c->prev) c->prev->next = c->next;
    if (c->next) c->next->prev = c->prev;
    if (c == array->first_child) array->first_child = c->next;
    c->prev = c->next = 0;
    return c;
}

void cJSON_DeleteItemFromArray(iftJson *array, int which) { cJSON_Delete(cJSON_DetachItemFromArray(array, which)); }

iftJson *cJSON_DetachItemFromObject(iftJson *object, const char *string) {
    int i = 0;
    iftJson *c = object->first_child;
    while (c && cJSON_strcasecmp(c->key, string)) i++, c = c->next;
    if (c) return cJSON_DetachItemFromArray(object, i);
    return 0;
}

void cJSON_DeleteItemFromObject(iftJson *object, const char *string) {
    cJSON_Delete(cJSON_DetachItemFromObject(object, string));
}

/* Replace array/object items with new ones. */
void cJSON_InsertItemInArray(iftJson *array, int which, iftJson *newitem) {
    iftJson *c = array->first_child;
    while (c && which > 0) c = c->next, which--;
    if (!c) {
        cJSON_AddItemToArray(array, newitem);
        return;
    }
    newitem->next = c;
    newitem->prev = c->prev;
    c->prev = newitem;
    if (c == array->first_child) array->first_child = newitem; else newitem->prev->next = newitem;
}

void cJSON_ReplaceItemInArray(iftJson *array, int which, iftJson *newitem) {
    iftJson *c = array->first_child;
    while (c && which > 0) c = c->next, which--;
    if (!c) return;
    newitem->next = c->next;
    newitem->prev = c->prev;
    if (newitem->next) newitem->next->prev = newitem;
    if (c == array->first_child) array->first_child = newitem; else newitem->prev->next = newitem;
    c->next = c->prev = 0;
    cJSON_Delete(c);
}

void cJSON_ReplaceItemInObject(iftJson *object, const char *string, iftJson *newitem) {
    int i = 0;
    iftJson *c = object->first_child;
    while (c && cJSON_strcasecmp(c->key, string))i++, c = c->next;
    if (c) {
        newitem->key = cJSON_strdup(string);
        cJSON_ReplaceItemInArray(object, i, newitem);
    }
}


/* Duplication */
iftJson *cJSON_Duplicate(iftJson *item, int recurse) {
    iftJson *newitem, *cptr, *nptr = 0, *newchild;
    /* Bail on bad ptr */
    if (!item) return 0;
    /* Create new item */
    newitem = cJSON_New_Item();
    if (!newitem) return 0;
    /* Copy over all vars */
    newitem->type =
            item->type & (~cJSON_IsReference), newitem->int_val = item->int_val, newitem->dbl_val = item->dbl_val;
    if (item->str_val) {
        newitem->str_val = cJSON_strdup(item->str_val);
        if (!newitem->str_val) {
            cJSON_Delete(newitem);
            return 0;
        }
    }
    if (item->key) {
        newitem->key = cJSON_strdup(item->key);
        if (!newitem->key) {
            cJSON_Delete(newitem);
            return 0;
        }
    }
    /* If non-recursive, then we're done! */
    if (!recurse) return newitem;
    /* Walk the ->next chain for the first_child. */
    cptr = item->first_child;
    while (cptr) {
        newchild = cJSON_Duplicate(cptr, 1);        /* Duplicate (with recurse) each item in the ->next chain */
        if (!newchild) {
            cJSON_Delete(newitem);
            return 0;
        }
        if (nptr) {
            nptr->next = newchild, newchild->prev = nptr;
            nptr = newchild;
            newitem->last_child = newchild;
        }    /* If newitem->first_child already set, then crosswire ->prev and ->next and move on */
        else {
            newitem->first_child = newitem->last_child = newchild;
            nptr = newchild;
        }                    /* Set newitem->first_child and move to it */
        cptr = cptr->next;
    }
    return newitem;
}

void cJSON_Minify(char *json) {
    char *into = json;
    while (*json) {
        if (*json == ' ') json++;
        else if (*json == '\t') json++;    /* Whitespace characters. */
        else if (*json == '\r') json++;
        else if (*json == '\n') json++;
        else if (*json == '/' && json[1] == '/')
            while (*json && *json != '\n')
                json++;    /* double-slash comments, to end of line. */
        else if (*json == '/' && json[1] == '*') {
            while (*json && !(*json == '*' && json[1] == '/')) json++;
            json += 2;
        }    /* multiline comments. */
        else if (*json == '\"') {
            *into++ = *json++;
            while (*json && *json != '\"') {
                if (*json == '\\') *into++ = *json++;
                *into++ = *json++;
            }
            *into++ = *json++;
        } /* key literals, which are \" sensitive. */
        else *into++ = *json++;            /* All other characters. */
    }
    *into = 0;    /* and null-terminate. */
}
/******************************************************************/


/* Default options for cJSON_Parse */
iftJson *cJSON_Parse(const char *value) {
    return cJSON_ParseWithOpts(value, 0, 0);
}


/* Create basic types: */
iftJson *cJSON_CreateNull(void) {
    iftJson *item = cJSON_New_Item();
    if (item)item->type = IFT_JSON_NULL;
    return item;
}

iftJson *cJSON_CreateTrue(void) {
    iftJson *item = cJSON_New_Item();
    if (item)item->type = IFT_JSON_TRUE;
    return item;
}

iftJson *cJSON_CreateFalse(void) {
    iftJson *item = cJSON_New_Item();
    if (item)item->type = IFT_JSON_FALSE;
    return item;
}

iftJson *cJSON_CreateBool(int b) {
    iftJson *item = cJSON_New_Item();
    if (item)item->type = b ? IFT_JSON_TRUE : IFT_JSON_FALSE;
    return item;
}

iftJson *cJSON_CreateNumber(double num) {
    iftJson *item = cJSON_New_Item();
    if (item) {
        item->type = IFT_JSON_DBL;
        item->dbl_val = num;
        item->int_val = (int) num;
    }
    return item;
}

iftJson *cJSON_CreateString(const char *string) {
    iftJson *item = cJSON_New_Item();
    if (item) {
        item->type = IFT_JSON_STRING;
        item->str_val = cJSON_strdup(string);
    }
    return item;
}

iftJson *cJSON_CreateArray(void) {
    iftJson *item = cJSON_New_Item();
    if (item)item->type = IFT_JSON_ARRAY;
    return item;
}

iftJson *cJSON_CreateObject(void) {
    iftJson *item = cJSON_New_Item();
    if (item)item->type = IFT_JSON_DICT;
    return item;
}

/* Create Arrays: */
iftJson *cJSON_CreateIntArray(const int *numbers, int count) {
    int i;
    iftJson *n = 0, *p = 0, *a = cJSON_CreateArray();
    for (i = 0; a && i < count; i++) {
        n = cJSON_CreateNumber(numbers[i]);
        if (!i)a->first_child = n; else suffix_object(p, n);
        p = n;
    }
    return a;
}

iftJson *cJSON_CreateFloatArray(const float *numbers, int count) {
    int i;
    iftJson *n = 0, *p = 0, *a = cJSON_CreateArray();
    for (i = 0; a && i < count; i++) {
        n = cJSON_CreateNumber(numbers[i]);
        if (!i)a->first_child = n; else suffix_object(p, n);
        p = n;
    }
    return a;
}

iftJson *cJSON_CreateDoubleArray(const double *numbers, int count) {
    int i;
    iftJson *n = 0, *p = 0, *a = cJSON_CreateArray();
    for (i = 0; a && i < count; i++) {
        n = cJSON_CreateNumber(numbers[i]);
        if (!i)a->first_child = n; else suffix_object(p, n);
        p = n;
    }
    return a;
}

iftJson *cJSON_CreateStringArray(const char **strings, int count) {
    int i;
    iftJson *n = 0, *p = 0, *a = cJSON_CreateArray();
    for (i = 0; a && i < count; i++) {
        n = cJSON_CreateString(strings[i]);
        if (!i)a->first_child = n; else suffix_object(p, n);
        p = n;
    }
    return a;
}

/************************** My Custom Wrappers/Functions *********************************/
/*********** MY PRIVATE FUNCTIONS AND DEFINITIONS *************/

/**
 * @brief Prints an iftJson following a printing format.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param json iftJson to be printed.
 */
void iftPrintJson(iftJson *json);


/**
 * @brief Prints an iftJson following a minified format.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param json iftJson to be printed.
 */
void iftPrintJsonMinified(iftJson *json);



/**
 * @brief Copies a Json (or a JDict or a JNode).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * @param json The iftJson to be copied.
 * @return The copied iftJson.
 */
iftJson *iftCopyJson(iftJson *json);

/**
 * @brief Copies a Json as a string that can be used for display and/or printing to a file.
 * @author Thiago V. Spina
 * @date Oct 28, 2016
 *
 * @param json The iftJson to be copied.
 * @return A string representing the json file.
 */
char *iftCopyJsonAsString(iftJson *json);

/**
 * @brief Returns an iftJNode REFERENCE from the key <key>. BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 * @sa iftGetDictReference
 *
 * Returns an iftJNode reference from the key <key>.
 * BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY. Then, if you destroy the <jnode> outside, the Json will be inconsistent.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the an intermediate key is not a dict, an Error message is going to be raised.
 *
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The reference of the iftJNode required.
 */
iftJNode *iftGetJNodeReference(iftJson *json, const char *key);


/**
 * @brief Returns an iftJNode (COPYINT IT) from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 * @sa iftGetDict
 *
 * Returns an iftJNode (COPYINT IT) from the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the an intermediate key is not a dict, an Error message is going to be raised.
 *
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The iftJNode required.
 */
iftJNode *iftGetJNode(iftJson *json, const char *key);


/**
 * @brief Returns a boolean value from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a boolean value from the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a Boolean, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The boolean value required.
 */
bool iftGetJBool(iftJson *json, const char *key);


/**
 * @brief Returns a NULL pointer from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a NULL pointer from the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not "null", an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return NULL, if the value for the key is "null".
 */
void *iftGetJNull(iftJson *json, const char *key);


/**
 * @brief Returns an Integer value from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns an Integer value with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key do not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an Integer, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Integer value required.
 */
int iftGetJInt(iftJson *json, const char *key);


/**
 * @brief Returns a double value from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a double value with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a Double, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The double value required.
 */
double iftGetJDouble(iftJson *json, const char *key);


/**
 * @brief Returns a string value (COPYING IT) from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a string value (COPYING IT) with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a String, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The string value required.
 */
char *iftGetJString(iftJson *json, const char *key);


/**
 * @brief Returns the reference of the iftJDict with the key <key>. BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns the reference of the iftJDict with the key <key>.
 * BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY. Then, if you destroy the <jnode> outside, the Json will be inconsistent.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an iftJDict, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The reference of the required iftJDict.
 */
iftJDict *iftGetJDictReference(iftJson *json, const char *key);


/**
 * @brief Returns the iftJDict (COPYING IT) with the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns the the iftJDict (COPYING IT) with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an iftJDict, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The required iftJDict.
 */
iftJDict *iftGetJDict(iftJson *json, const char *key);


/**
 * @brief Check if a Json Node is an Integer
 * @param  jnode JNode must be a number
 * @return       true/false
 */
static inline bool iftIsJInt(const iftJNode *jnode) {
    return (jnode->int_val == jnode->dbl_val);
}


/**
 * @brief Check if all elements from an JArray is integer.
 * @param  jarray JNode must be an array with NUMBER elements.
 * @return      true/false
 *
 * @author Samuka Martins
 * @date Sep 17, 2017
 */
bool iftIsJIntArray(const iftJNode *jarray);



/**
 * @brief Returns an Integer Array from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns an Integer Array with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an Integer Array, an Error message is also going to be raised.
 * If some value of the array is not an Integer, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Integer Array required.
 */
iftIntArray *iftGetJIntArray(iftJson *json, const char *key);


/**
 * @brief Returns a Double Array from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a Double Array with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a Double Array, an Error message is also going to be raised.
 * If some value of the array is not a Double, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Double Array required.
 */
iftDblArray *iftGetJDblArray(iftJson *json, const char *key);


/**
 * @brief Returns a String Array from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a String Array with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a String Array, an Error message is also going to be raised.
 * If some value of the array is not a String, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The String Array required.
 */
iftStrArray *iftGetJStrArray(iftJson *json, const char *key);


/**
 * @brief Returns an Integer Matrix from the key <key>.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * Returns an Integer Matrix with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an Integer Matrix, an Error message is also going to be raised.
 * If some value of the matrix is not an Integer, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Integer Matrix required.
 */
iftIntMatrix *iftGetJIntMatrix(iftJson *json, const char *key);


/**
 * @brief Returns a Double Matrix from the key <key>.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * Returns a Double Matrix with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a Double Matrix, an Error message is also going to be raised.
 * If some value of the matrix is not a Double, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Double Matrix required.
 */
iftMatrix *iftGetJDblMatrix(iftJson *json, const char *key);


/**
 * @brief Returns a String Matrix (COPYING each string) from the key <key>.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * Returns a String Matrix (COPYING each string) with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a String Matrix, an Error message is also going to be raised.
 * If some value of the matrix is not a String, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The String Matrix required.
 */
iftStrMatrix *iftGetJStrMatrix(iftJson *json, const char *key);


// TODO: iftGetJDictArray


/**
 * @brief Creates a Json root which consists of an only Json node with its elements pointing to NULL.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * Creates a Json root which consists of an only Json node with its elements pointing to NULL.
 * It is necessary to add new elements to Json.
 *
 * @return The created iftJson root.
 */
iftJson *iftCreateJsonRoot(void);


/**
 * @brief Creates an EMPTY Json Node with its elements pointing to NULL (Similar to iftCreateJsonRoot()).
 *
 * @author Samuel Martins
 * @date October 20, 2015
 * @sa iftCreateJsonRoot()
 *
 * @return The created iftJNode.
 */
iftJNode *iftCreateJNode();


/**
 * @brief Creates a Json Node with type IFT_JSON_ARRAY and its elements pointing to NULL.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @return The created Json Array Node.
 */
iftJNode *iftCreateJArray();


/**
 * @brief Creates a JDict with type IFT_JSON_DICT and its elements pointing to NULL (Similar to iftCreateJsonRoot()).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 * @sa iftCreateJsonRoot()
 *
 * @return The created iftJDict.
 */
iftJDict *iftCreateJDict();


/**
 * @brief Adds a JNode REFERENCE to a Json <json> in the key <key>. BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY.
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds a JNode <jnode> REFERENCE to Json in the key <key>.
 * BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY. Then, if you destroy the <jnode> outside, the Json will be inconsistent.
 * If the key already exists already exists, it is going to be destroyed and the new JNode will be assigned instead.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @param jnode The reference of the JNode to added.
 */
void iftAddJNodeReferenceToJson(iftJson *json, const char *key, iftJNode *jnode);


/**
 * @brief Adds a JNode (COPYING IT) to a Json <json> in the key <key>.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds a JNode <jnode> (COPYING IT) to Json in the key <key>.
 * If the key already exists already exists, it is going to be destroyed and the new JNode will be assigned instead.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @param jnode The JNode to be copied and added to Json.
 */
void iftAddJNodeToJson(iftJson *json, const char *key, iftJNode *jnode);


/**
 * @brief Adds or sets a new Boolean value to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets a new Boolean value <bool_val> to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param bool_val The Boolean value to be added into.
 * @param key The key to be used.
 */
void iftAddBoolToJson(iftJson *json, const char *key, unsigned int bool_val);


/**
 * @brief Adds or sets a new Null value to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets a new Null value to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param key The key to be used.
 */
void iftAddNullToJson(iftJson *json, const char *key);


/**
 * @brief Adds or sets a new integer value to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets a new integer value <int_val> to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param int_val The integer to be added into.
 * @param key The key to be used.
 */
void iftAddIntToJson(iftJson *json, const char *key, int int_val);


/**
 * @brief Adds or sets a new double value to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets a new double value <dbl_val> to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param dbl_val The double to be added into.
 * @param key The key to be used.
 */
void iftAddDoubleToJson(iftJson *json, const char *key, double dbl_val);


/**
 * @brief Adds or sets a new string value (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets a new string value <str_val> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param str_val The string to be added into.
 * @param key The key to be used.
 */
void iftAddStringToJson(iftJson *json, const char *key, const char *str_val);


/**
 * @brief Adds or sets the REFERENCE of a JDict to a Json (or a JDict). BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY.
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets the REFERENCE of the JDict <jdict> to the iftJson or iftJDict <json> in the key <key>.
 * BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY. Then, if you destroy the <jdict> outside, the Json will be inconsistent.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param jdict The reference of the iftJDict to be added.
 * @param key The key to be used.
 */
void iftAddJDictReferenceToJson(iftJson *json, const char *key, iftJDict *jdict_ref);


/**
 * @brief Adds or sets a JDict (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets the JDict <jdict> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param jdict The iftJDict to be added.
 * @param key The key to be used.
 */
void iftAddJDictToJson(iftJson *json, const char *key, iftJDict *jdict);


/**
 * @brief Adds or sets an Integer Array (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the Integer Array <iarr> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param iarr The iftIntArray to be added.
 * @param key The key to be used.
 */
void iftAddIntArrayToJson(iftJson *json, const char *key, const iftIntArray *iarr);


/**
 * @brief Adds or sets a Double Array (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the Double Array <darr> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param darr The iftDblArray to be added.
 * @param key The key to be used.
 */
void iftAddFloatArrayToJson(iftJson *json, const char *key, iftFloatArray *darr);


/**
 * @brief Adds or sets a String Array (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the String Array <sarr> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param sarr The iftStrArray to be added.
 * @param key The key to be used.
 */
void iftAddStringArrayToJson(iftJson *json, const char *key, iftStrArray *sarr);


/**
 * @brief Adds or sets an Integer Matrix (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the Integer Matrix <iM> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param iM The iftIntArray to be added.
 * @param key The key to be used.
 */
void iftAddIntMatrixToJson(iftJson *json, const char *key, iftIntMatrix *iM);


/**
 * @brief Adds or sets a Double Matrix (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the Double Matrix <dM> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param dM The iftDblArray to be added.
 * @param key The key to be used.
 */
void iftAddDoubleMatrixToJson(iftJson *json, const char *key, iftMatrix *dM);


/**
 * @brief Adds or sets a String Matrix (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the String Matrix <sM> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param sM The iftStrArray to be added.
 * @param key The key to be used.
 */
void iftAddStrMatrixToJson(iftJson *json, const char *key, iftStrMatrix *sM);


/**
 * @brief Adds (APPENDING IT) a new Integer into Json Array Node (node with type = IFT_JSON_ARRAY). BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new Integer into Json Array Node (node with type = IFT_JSON_ARRAY).
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * No key is necessary.
 * The value is inserted at the end of the array.
 * If the JNode used is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 *
 * @param jarray The JArray node where the value will be added.
 * @param int_val The Integer to be added.
 */
void iftAddIntToJArrayNode(iftJNode *jarray, int int_val);


/**
 * @brief Adds (APPENDING IT) a new Double into Json Array Node (node with type = IFT_JSON_ARRAY). BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new Double into Json Array Node (node with type = IFT_JSON_ARRAY).
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * No key is necessary.
 * The value is inserted at the end of the array.
 * If the JNode used is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 *
 * @param jarray The JArray node where the value will be added.
 * @param dbl_val The Double to be added.
 */
void iftAddDoubleToJArrayNode(iftJNode *jarray, double dbl_val);


/**
 * @brief Adds (APPENDING IT) a new String into Json Array Node (node with type = IFT_JSON_ARRAY). BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new String into Json Array Node (node with type = IFT_JSON_ARRAY).
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * No key is necessary.
 * The value is inserted at the end of the array.
 * If the JNode used is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 *
 * @param jarray The JArray node where the value will be added.
 * @param str_val The String to be added.
 */
void iftAddStringToJArrayNode(iftJNode *jarray, const char *str_val);


/**
 * @brief Adds (APPENDING IT) a new Integer into Json Array with key <key>. BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new Integer into Json Array with key <key>.
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * The value is inserted at the end of the array.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the JNode with key <key> is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The Json (or JDict) to be used.
 * @param key The key of the JArray to be used.
 * @param int_val The Integer to be added.
 */
void iftAddIntToJArray(iftJson *json, const char *key, int int_val);


/**
 * @brief Adds (APPENDING IT) a new Double into Json Array with key <key>. BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new Double into Json Array with key <key>.
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * The value is inserted at the end of the array.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the JNode with key <key> is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The Json (or JDict) to be used.
 * @param key The key of the JArray to be used.
 * @param dbl_val The Double to be added.
 */
void iftAddDoubleToJArray(iftJson *json, const char *key, double dbl_val);


/**
 * @brief Adds (APPENDING IT) a new String into Json Array with key <key>. BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new String into Json Array with key <key>.
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * The value is inserted at the end of the array.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the JNode with key <key> is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The Json (or JDict) to be used.
 * @param key The key of the JArray to be used.
 * @param str_val The String to be added.
 */
void iftAddStringToJArray(iftJson *json, const char *key, const char *str_val);


/**
 * @brief Deletes (deallocates) the JNode with the key <key> and its descendant from a Json or a JDict.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Deletes (deallocates) the JNode with the key <key> and its descendant from a Json or a JDict.
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param key The key from the node to be deleted.
 */
void iftDeleteJNode(iftJson *json, const char *key);


bool iftJsonContainKey(const char *key, iftJson *json);







/****************************** FUNCTION IMPLEMENTATION ******************************/



/**
 * @brief Counts the number of elements of a Json Array.
 *
 * @author Samuel Martins
 * @date October 16, 2015
 *
 * @param jarray Json array to be analyzed.
 * @param key Json key just used in the error printing.
 * @return The number of elements.
 */
int _iftCountJNumElems(iftJson *jarray, iftJsonType jtype, const char *key) {
    if (jarray->type != IFT_JSON_ARRAY)
        iftError("The key \"%s\" is not an Array", "_iftCountJNumElems", key);
    
    int n = 0;
    iftJNode *jaux = jarray->first_child; // gets the first element of the array
    
    // counts the number of elements and checks if all of them are integers
    while (jaux != NULL) {
        if (jaux->type != jtype)
            iftError("Elements with different types in the array with key: \"%s\"", "_iftCountJNumElems", key);
        
        n++;
        jaux = jaux->next;
    }
    
    return n;
}


/**
 * @brief Counts the number of Rows and Columns of a Json Matrix.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * @param jmatrix Json matrix to be analyzed.
 * @param key Json key just used in the error printing.
 * @param nrows Reference of the Number of Rows to be counted.
 * @param ncols Reference of the Number of Cols to be counted.
 */
void _iftCountJNumRowsAndCols(iftJson *jmatrix, iftJsonType jtype, const char *key, int *nrows, int *ncols) {
    if (jmatrix->type != IFT_JSON_ARRAY)
        iftError("The key \"%s\" is not an Array", "_iftCountJNumRowsAndCols", key);
    
    *nrows = 0;
    *ncols = 0;
    iftJNode *jrow = jmatrix->first_child; // gets the first sub-array
    iftJNode *jcol = NULL;
    
    // gets the first row just to compute the number of cols
    if (jrow != NULL) {
        if (jrow->type != IFT_JSON_ARRAY)
            iftError("Invalid elements (not sub-array) in the matrix with key: \"%s\"",
                     "iftGet_iftCountJNumRowsAndColsJIntMatrix", key);
        (*nrows)++;
        
        jcol = jrow->first_child; // gets the row array
        
        while (jcol != NULL) {
            if (jcol->type != jtype)
                iftError("Elements with different types in the matrix with key: \"%s\"", "_iftCountJNumRowsAndCols",
                         key);
            (*ncols)++;
            jcol = jcol->next;
        }
    }
    
    jrow = jrow->next; // gets the second row
    
    // counts the number of rows, checks if the matrix is regular, and checks if all of them are integers
    while (jrow != NULL) {
        if (jrow->type != IFT_JSON_ARRAY)
            iftError("Invalid elements (not sub-array) in the matrix with the key: \"%s\"",
                     "iftGet_iftCountJNumRowsAndColsJIntMatrix", key);
        (*nrows)++;
        
        jcol = jrow->first_child; // gets the row array
        long int _ncols = 0;
        
        while (jcol != NULL) {
            if (jcol->type != jtype)
                iftError("Matrix with elements non-integers", "_iftCountJNumRowsAndCols");
            _ncols++;
            jcol = jcol->next;
        }
        
        jrow = jrow->next;
        
        // checks if the matrix is regular
        if (*ncols != _ncols)
            iftError("The integer matrix is not regular: Rows with different number of columns",
                     "_iftCountJNumRowsAndCols");
    }
}


/**
 * @brief Creates a Json Boolean node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param bool_val The boolean value used.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJBool(unsigned int bool_val) {
    if ((bool_val != 0) && (bool_val != 1))
        iftError("Invalid Boolean value: %d... Try 0 or 1", "iftCreateJBool", bool_val);
    
    iftJNode *jbool = iftCreateJNode();
    jbool->type     = bool_val;
    jbool->int_val  = bool_val;
    jbool->dbl_val  = (double) bool_val;
    
    return jbool;
}


/**
 * @brief Creates a Json Null node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @return The created iftJNode.
 */
iftJNode *iftCreateJNull() {
    iftJNode *jnull = iftCreateJNode();
    jnull->type     = IFT_JSON_NULL;
    
    return jnull;
}


/**
 * @brief Creates a Json Integer node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param int_val The integer value used.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJInt(int int_val) {
    iftJNode *jint = iftCreateJNode();
    jint->type     = IFT_JSON_INT;
    jint->int_val  = int_val;
    jint->dbl_val  = (double) int_val;
    
    return jint;
}


/**
 * @brief Creates a Json Double node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param dbl_val The double value used.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJDouble(double dbl_val) {
    iftJNode *jdbl = iftCreateJNode();
    jdbl->type     = IFT_JSON_DBL;
    jdbl->dbl_val  = dbl_val;
    jdbl->int_val  = (int) dbl_val;
    
    return jdbl;
}


/**
 * @brief Creates a Json String node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param str_val The string value used.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJString(const char *str_val) {
    iftJNode *jstr = iftCreateJNode();
    jstr->type     = IFT_JSON_STRING;
    if (str_val != NULL) {
        jstr->str_val = iftAllocCharArray(strlen(str_val)+1);
        strcpy(jstr->str_val, str_val);
    }
    
    return jstr;
}


/**
 * @brief Creates a Json Integer Array node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param iarr The Integer Array to be copied.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJIntArray(const iftIntArray *iarr) {
    iftJNode *jarray = iftCreateJNode();
    jarray->type     = IFT_JSON_ARRAY;
    
    // inserts the first child element
    // code duplication in order to avoid intermediate ifs
    iftJNode *jchild = NULL;
    if (iarr->n != 0) {
        jchild              = iftCreateJInt(iarr->val[0]);
        jarray->first_child = jchild;
        
        for (int i = 1; i < iarr->n; i++) {
            iftJNode *jnew_child = iftCreateJInt(iarr->val[i]);
            
            jchild->next     = jnew_child;
            jnew_child->prev = jchild;
            jchild           = jnew_child;
        }
    }
    
    jarray->last_child = jchild;
    
    return jarray;
}


/**
 * @brief Creates a Json Double Array node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param darr The Double Array to be copied.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJFloatArray(iftFloatArray *darr) {
    iftJNode *jarray = iftCreateJNode();
    jarray->type     = IFT_JSON_ARRAY;
    
    // inserts the first child element
    // code duplication in order to avoid intermediate ifs
    iftJNode *jchild = NULL;
    if (darr->n != 0) {
        jchild              = iftCreateJDouble(darr->val[0]);
        jarray->first_child = jchild;
        
        for (int i = 1; i < darr->n; i++) {
            iftJNode *jnew_child = iftCreateJDouble(darr->val[i]);
            
            jchild->next     = jnew_child;
            jnew_child->prev = jchild;
            jchild           = jnew_child;
        }
    }
    
    jarray->last_child = jchild;
    
    return jarray;
}


/**
 * @brief Creates a Json String Array node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param sarr The String Array to be copied.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJStringArray(iftStrArray *sarr) {
    iftJNode *jarray = iftCreateJNode();
    jarray->type     = IFT_JSON_ARRAY;
    
    // inserts the first child element
    // code duplication in order to avoid intermediate ifs
    iftJNode *jchild = NULL;
    if (sarr->n != 0) {
        jchild              = iftCreateJString(sarr->val[0]);
        jarray->first_child = jchild;
        
        for (int i = 1; i < sarr->n; i++) {
            iftJNode *jnew_child = iftCreateJString(sarr->val[i]);
            
            jchild->next     = jnew_child;
            jnew_child->prev = jchild;
            jchild           = jnew_child;
        }
    }
    
    jarray->last_child = jchild;
    
    return jarray;
}


/**
 * @brief Converts a linear JArray into a JMatrix IN PLACE.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param jmatrix The linear JArray to be converted Into a JMatrix IN PLACE.
 * @param nrows The number of Rows from the Matrix.
 * @param nrows The number of Columns from the Matrix.
 */
void iftConvertJArrayToJMatrixInPlace(iftJNode **jmatrix, int nrows, int ncols) {
    iftJNode *jrow    = NULL;
    iftJNode *jcol    = NULL;
    
    iftJNode *jmatrix_aux   = *jmatrix;
    *jmatrix                = iftCreateJNode();
    (*jmatrix)->type        = IFT_JSON_ARRAY;
    (*jmatrix)->first_child = jmatrix_aux;
    
    /***** Transform the linear JArray into a JMatrix ****/
    // inserts the first child element
    // code duplication in order to avoid intermediate ifs
    jrow       = jmatrix_aux;
    jcol       = jrow->first_child;
    jcol->prev = NULL;
    // scans until finding the linear array last node of the row from the JMatrix
    for (int j = 1; j < ncols; j++) {
        jcol = jcol->next;
    }
    jrow->last_child = jcol;
    
    
    for (int i = 1; i < nrows; i++) {
        jrow->next       = iftCreateJNode();
        jrow->next->type = IFT_JSON_ARRAY;
        jrow->next->prev = jrow;
        jrow             = jrow->next;
        
        jcol              = jcol->next;
        jcol->prev->next  = NULL;
        jcol->prev        = NULL;
        jrow->first_child = jcol;
        
        for (int j = 1; j < ncols; j++) {
            jcol = jcol->next;
        }
        jrow->last_child = jcol;
    }
    jcol->next = NULL;
    jrow->next = NULL;
    (*jmatrix)->last_child = jrow;
}


/**
 * @brief Creates a Json Integer Matrix node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param sarr The Integer Matrix to be copied.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJIntMatrix(iftIntMatrix *iM) {
    int nrows = iM->nrows;
    int ncols = iM->ncols;
    
    // "converts" the Matrix into a Linear Array
    iftIntArray *iarr = (iftIntArray*) iftAlloc(1, sizeof(iftIntArray));
    iarr->val         = iM->val;
    iarr->n = iM->n;
    
    // JArray of the Matrix as a linear array
    iftJNode *jmatrix = iftCreateJIntArray(iarr);
    iftConvertJArrayToJMatrixInPlace(&jmatrix, nrows, ncols);
    
    iftFree(iarr); // just deallocates the point to iftIntArray
    
    return jmatrix;
}


/**
 * @brief Creates a Json Double Matrix node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param sarr The Double Matrix to be copied.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJFloatMatrix(iftMatrix *dM) {
    int nrows = dM->nrows;
    int ncols = dM->ncols;
    
    // "converts" the Matrix into a Linear Array
    iftFloatArray *darr = (iftFloatArray *) iftAlloc(1, sizeof(iftFloatArray));
    darr->val            = dM->val;
    darr->n              = dM->n;
    
    // JArray of the Matrix as a linear array
    iftJNode *jmatrix = iftCreateJFloatArray(darr);
    iftConvertJArrayToJMatrixInPlace(&jmatrix, nrows, ncols);
    
    iftFree(darr); // just deallocates the point to iftIntArray
    
    return jmatrix;
}


/**
 * @brief Creates a Json Double Matrix node (without a key).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @param sarr The Double Matrix to be copied.
 * @return The created iftJNode.
 */
iftJNode *iftCreateJStrMatrix(iftStrMatrix *sM) {
    int nrows = sM->nrows;
    int ncols = sM->ncols;
    
    // "converts" the Matrix into a Linear Array
    iftStrArray *sarr = (iftStrArray *) iftAlloc(1, sizeof(iftStrArray));
    sarr->val            = sM->val;
    sarr->n              = sM->n;
    
    // JArray of the Matrix as a linear array
    iftJNode *jmatrix = iftCreateJStringArray(sarr);
    iftConvertJArrayToJMatrixInPlace(&jmatrix, nrows, ncols);
    
    iftFree(sarr); // just deallocates the point to iftIntArray
    
    return jmatrix;
}
/****************************************************/


iftJson *iftReadJsonExternalLibStruct(const char *json_path) {
    iftJson *json = NULL;
    
    char *text = iftReadFileAsString(json_path);
    json = cJSON_Parse(text);
    
    if (!json)
        iftError("JSON Error before: [%s]\n", "iftReadJsonExternalLibStruct", cJSON_GetErrorPtr());
    iftFree(text);
    
    return json;
}


void iftDestroyJson(iftJson **json) {
    iftJson *jaux = *json;
    
    if (jaux != NULL) {
        cJSON_Delete(jaux);
        *json = NULL;
    }
}


void iftPrintJson(iftJson *json) {
    char *text = cJSON_Print(json);
    puts(text);
    iftFree(text);
}

void iftPrintJsonMinified(iftJson *json) {
    char *text = cJSON_PrintUnformatted(json);
    puts(text);
    iftFree(text);
}


iftJson *iftCopyJson(iftJson *json) {
    return (cJSON_Duplicate(json, 1));
}


char *iftCopyJsonAsString(iftJson *json) {
    return cJSON_Print(json);
}

iftJNode *iftGetJNodeReference(iftJson *json, const char *key) {
    char keys_error_msg[1024] = ""; // used to show the key hierarchy if a key is not found
    char str[128];
    
    if (json == NULL)
        iftError("Json is NULL", "iftGetJNode");
    
    iftJNode *jnode    = json;
    iftSList *key_list = iftSplitString(key, ":");
    char *elem_key     = iftRemoveSListTail(key_list); // gets the key of the element
    
    // checking if all intermediate keys are dicts and gets the penultimate dict
    while (key_list->n != 0) {
        char *intermediate_key = iftRemoveSListHead(key_list);
        
        sprintf(str, "[%s]", intermediate_key);
        strcat(keys_error_msg, str);
        
        jnode = cJSON_GetObjectItem(jnode, intermediate_key);
        
        if (jnode == NULL)
            iftError("Invalid key: \"%s\".\nFull keys: \"%s\"", "iftGetJNode", intermediate_key, keys_error_msg);
        if (jnode->type != IFT_JSON_DICT)
            iftError("Intermediate key: \"%s\" is not a Dictionary.\nFull keys: \"%s\"", "iftGetJNode",
                     intermediate_key,
                     keys_error_msg);
        iftFree(intermediate_key);
    }
    iftDestroySList(&key_list);
    
    // gets the last element
    jnode = cJSON_GetObjectItem(jnode, elem_key);
    if (jnode == NULL) {
        sprintf(str, "[%s]", elem_key);
        strcat(keys_error_msg, str);
        iftError("Invalid key: \"%s\".\nFull keys: \"%s\"", "iftGetJNode", elem_key, keys_error_msg);
    }
    iftFree(elem_key);
    
    return jnode;
}


iftJNode *iftGetJNode(iftJson *json, const char *key) {
    iftJNode *jnode = iftGetJNodeReference(json, key);
    return (iftCopyJson(jnode));
}


bool iftGetJBool(iftJson *json, const char *key) {
    iftJNode *jnode = iftGetJNodeReference(json, key);
    
    if ((jnode->type != IFT_JSON_FALSE) && (jnode->type != IFT_JSON_TRUE))
        iftError("The key \"%s\" is not a Boolean value", "iftGetJBool", key);
    
    return ((bool) jnode->type);
}


void *iftGetJNull(iftJson *json, const char *key) {
    iftJNode *jnode = iftGetJNodeReference(json, key);
    
    if (jnode->type != IFT_JSON_NULL)
        iftError("The key \"%s\" is not a Null value", "iftGetJNull", key);
    
    return NULL;
}


int iftGetJInt(iftJson *json, const char *key) {
    iftJNode *jnode = iftGetJNodeReference(json, key);
    
    if (jnode->type != IFT_JSON_INT)
        iftError("The key \"%s\" is not an Integer value", "iftGetJInt", key);
    
    return jnode->int_val;
}


double iftGetJDouble(iftJson *json, const char *key) {
    iftJNode *jnode = iftGetJNodeReference(json, key);
    
    if (jnode->type != IFT_JSON_DBL)
        iftError("The key \"%s\" is not a Double value", "iftGetJDouble", key);
    
    return jnode->dbl_val;
}


char *iftGetJString(iftJson *json, const char *key) {
    iftJNode *jnode = iftGetJNodeReference(json, key);
    
    if (jnode->type != IFT_JSON_STRING)
        iftError("The key \"%s\" is not a String value", "iftGetJString", key);
    
    char *str = iftAllocCharArray(strlen(jnode->str_val)+1);
    strcpy(str, jnode->str_val);
    
    return str;
}


iftJDict *iftGetJDictReference(iftJson *json, const char *key) {
    iftJNode *jdict = iftGetJNodeReference(json, key);
    
    if (jdict->type != IFT_JSON_DICT)
        iftError("The key \"%s\" is not a Double value", "iftGetJDouble", key);
    
    return jdict;
}


bool iftIsJIntArray(const iftJNode *jarray) {
    for (iftJNode *node = jarray->first_child; node != NULL ; node = node->next) {
        printf("(%d != %lf) = %s\n", node->int_val, node->dbl_val, iftBoolAsString(node->int_val != node->dbl_val));
        
        if (node->int_val != node->dbl_val) {
            printf("%d != %lf\n", node->int_val, node->dbl_val);
            return false;
        }
    }
    
    return true;
}


iftJDict *iftGetJDict(iftJson *json, const char *key) {
    iftJNode *jdict = iftGetJDictReference(json, key);
    
    return (iftCopyJson(jdict));
}


iftIntArray *iftGetJIntArray(iftJson *json, const char *key) {
    iftJNode *jarray = iftGetJNodeReference(json, key);
    int n            = _iftCountJNumElems(jarray, IFT_JSON_INT, key);
    
    iftIntArray *iarr = iftCreateIntArray(n);
    iftJNode *jaux    = jarray->first_child; // gets the first element of the array
    int i             = 0;
    
    while (jaux != NULL) {
        if ((((double) jaux->int_val) - jaux->dbl_val) > DBL_EPSILON)
            iftError("Invalid elements (non-integers) in the array with key: \"%s\"", "iftGetJIntArray", key);
        
        iarr->val[i++] = jaux->int_val;
        jaux           = jaux->next;
    }
    
    return iarr;
}


iftDblArray *iftGetJDblArray(iftJson *json, const char *key) {
    iftJNode *jarray = iftGetJNodeReference(json, key);
    int n       = _iftCountJNumElems(jarray, IFT_JSON_DBL, key);
    
    iftDblArray *darr = iftCreateDblArray(n);
    iftJNode *jaux       = jarray->first_child; // gets the first element of the array
    int i                = 0;
    
    while (jaux != NULL) {
        // printf("(%d == %lf) ? = %s\n", jaux->int_val, jaux->dbl_val, iftBoolAsString(jaux->int_val == jaux->dbl_val));
        darr->val[i++] = jaux->dbl_val;
        jaux           = jaux->next;
    }
    
    return darr;
}


iftStrArray *iftGetJStrArray(iftJson *json, const char *key) {
    iftJNode *jarray = iftGetJNodeReference(json, key);
    int n            = _iftCountJNumElems(jarray, IFT_JSON_STRING, key);
    
    iftStrArray *sarr = (iftStrArray *) iftAlloc(1, sizeof(iftStrArray));
    sarr->n              = n;
    sarr->val            = (char**) iftAlloc(n, sizeof(char*));
    
    iftJNode *jaux = jarray->first_child; // gets the first element of the array
    int i          = 0;
    
    while (jaux != NULL) {
        sarr->val[i] = iftAllocCharArray(strlen(jaux->str_val)+1);
        strcpy(sarr->val[i++], jaux->str_val);
        jaux         = jaux->next;
    }
    
    return sarr;
}


iftIntMatrix *iftGetJIntMatrix(iftJson *json, const char *key) {
    int nrows, ncols;
    
    iftJNode *jmatrix = iftGetJNodeReference(json, key);
    
    _iftCountJNumRowsAndCols(jmatrix, IFT_JSON_INT, key, &nrows, &ncols);
    
    iftIntMatrix *iM = iftCreateIntMatrix(ncols, nrows);
    
    iftJNode *jrow = jmatrix->first_child; // gets the first sub-array
    iftJNode *jcol = NULL;
    int i = 0;
    
    while (jrow != NULL) {
        jcol = jrow->first_child;
        
        while (jcol != NULL) {
            if ((((double) jcol->int_val) - jcol->dbl_val) > DBL_EPSILON)
                iftError("Invalid elements (non-integers) in the matrix with key: \"%s\"", "iftGetJIntMatrix", key);
            
            iM->val[i] = jcol->int_val;
            jcol = jcol->next;
            i++;
        }
        jrow = jrow->next;
    }
    
    return iM;
}


iftMatrix *iftGetJDblMatrix(iftJson *json, const char *key) {
    int nrows, ncols;
    
    iftJNode *jmatrix = iftGetJNodeReference(json, key);
    
    _iftCountJNumRowsAndCols(jmatrix, IFT_JSON_DBL, key, &nrows, &ncols);
    
    iftMatrix *dM = iftCreateMatrix(ncols, nrows);
    
    iftJNode *jrow = jmatrix->first_child; // gets the first sub-array
    iftJNode *jcol = NULL;
    int i = 0;
    
    while (jrow != NULL) {
        jcol = jrow->first_child;
        
        while (jcol != NULL) {
            dM->val[i] = jcol->dbl_val;
            jcol = jcol->next;
            i++;
        }
        jrow = jrow->next;
    }
    
    return dM;
}


iftStrMatrix *iftGetJStrMatrix(iftJson *json, const char *key) {
    int nrows, ncols;
    
    iftJNode *smatrix = iftGetJNodeReference(json, key);
    
    _iftCountJNumRowsAndCols(smatrix, IFT_JSON_STRING, key, &nrows, &ncols);
    iftStrMatrix *sM = (iftStrMatrix*) iftAlloc(1, sizeof(iftStrMatrix));
    sM->nrows = nrows;
    sM->ncols = ncols;
    sM->n = nrows*ncols;
    sM->tbrow = iftAllocIntArray(nrows);
    for (int r = 0; r < nrows; r++) {
        sM->tbrow[r] = r*ncols;
    }
    
    sM->val = (char**) iftAlloc(sM->n, sizeof(char*));
    
    iftJNode *jrow = smatrix->first_child; // gets the first sub-array
    iftJNode *jcol = NULL;
    int i = 0;
    
    while (jrow != NULL) {
        jcol = jrow->first_child;
        
        while (jcol != NULL) {
            sM->val[i] = iftAllocCharArray(strlen(jcol->str_val)+1);
            strcpy(sM->val[i], jcol->str_val);
            jcol = jcol->next;
            i++;
        }
        jrow = jrow->next;
    }
    
    return sM;
}


iftJson *iftCreateJsonRoot(void) {
    return (cJSON_CreateObject());
}


iftJNode *iftCreateJNode() {
    return (cJSON_New_Item());
}


iftJNode *iftCreateJArray() {
    return (cJSON_CreateArray());
}


iftJDict *iftCreateJDict() {
    return (cJSON_CreateObject());
}


void iftAddJNodeReferenceToJson(iftJson *json, const char *key, iftJNode *jnode) {
    char keys_error_msg[1024] = ""; // used to show the key hierarchy if a key is not found
    char str[128];
    
    if (json == NULL)
        iftError("Json is NULL", "iftAddJNodeReferenceToJson");
    if (jnode == NULL) {
        iftWarning("Input JNode is NULL... Nothing to Do", "iftAddJNodeReferenceToJson");
        return;
    }
    
    iftSList *key_list = iftSplitString(key, ":");
    char *elem_key     = iftRemoveSListTail(key_list); // gets the key of the element to added (it can not exist).
    jnode->key         = elem_key;
    
    iftJson *jdict = json;
    
    while (key_list->n != 0) {
        char *intermediate_key = iftRemoveSListHead(key_list);
        
        sprintf(str, "[%s]", intermediate_key);
        strcat(keys_error_msg, str);
        
        jdict = cJSON_GetObjectItem(jdict, intermediate_key);
        
        if (jdict == NULL)
            iftError("Invalid key: \"%s\".\nFull keys: \"%s\"", "iftAddJNodeReferenceToJson", intermediate_key,
                     keys_error_msg);
        if (jdict->type != IFT_JSON_DICT)
            iftError("Intermediate key: \"%s\" is not a Dictionary.\nFull keys: \"%s\"", "iftAddJNodeReferenceToJson",
                     intermediate_key, keys_error_msg);
        iftFree(intermediate_key);
    }
    
    
    iftJNode *jnode_curr = cJSON_GetObjectItem(jdict, elem_key);
    
    // the key does not exist
    if (jnode_curr == NULL) {
        // empty dict
        if (jdict->first_child == NULL) {
            jdict->first_child = jdict->last_child = jnode;
        }
        else {
            jdict->last_child->next = jnode;
            jnode->prev             = jdict->last_child;
            jdict->last_child       = jnode;
        }
    }
        // the key exists - the input JNode is positioned in place of the current JNode
        // then, the current JNode is destroyed
    else {
        if (jdict->first_child == jnode_curr)
            jdict->first_child = jnode;
        if (jdict->last_child == jnode_curr)
            jdict->last_child = jnode;
        if (jnode_curr->prev != NULL)
            jnode_curr->prev->next = jnode;
        if (jnode_curr->next != NULL)
            jnode_curr->next->prev = jnode;
        
        jnode->prev      = jnode_curr->prev;
        jnode->next      = jnode_curr->next;
        jnode_curr->prev = NULL;
        jnode_curr->next = NULL;
        
        iftDestroyJson(&jnode_curr);
    }
    
    iftDestroySList(&key_list);
}


void iftAddJNodeToJson(iftJson *json, const char *key, iftJNode *jnode) {
    iftJNode *jcopy = iftCopyJson(jnode);
    iftAddJNodeReferenceToJson(json, key, jcopy);
}


void iftAddBoolToJson(iftJson *json, const char *key, unsigned int bool_val) {
    if ((bool_val != 0) && (bool_val != 1))
        iftError("Invalid Boolean value: %d... Try 0 or 1", "iftAddBoolToJson", bool_val);
    
    iftJNode *jbool = iftCreateJBool(bool_val);
    iftAddJNodeReferenceToJson(json, key, jbool);
}


void iftAddNullToJson(iftJson *json, const char *key) {
    iftJNode *jnull = iftCreateJNull();
    iftAddJNodeReferenceToJson(json, key, jnull);
}


void iftAddIntToJson(iftJson *json, const char *key, int int_val) {
    iftJNode *jint = iftCreateJInt(int_val);
    iftAddJNodeReferenceToJson(json, key, jint);
}


void iftAddDoubleToJson(iftJson *json, const char *key, double dbl_val) {
    iftJNode *jdbl = iftCreateJDouble(dbl_val);
    iftAddJNodeReferenceToJson(json, key, jdbl);
}


void iftAddStringToJson(iftJson *json, const char *key, const char *str_val) {
    iftJNode *jstr = iftCreateJString(str_val);
    iftAddJNodeReferenceToJson(json, key, jstr);
}


void iftAddJDictReferenceToJson(iftJson *json, const char *key, iftJDict *jdict_ref) {
    iftAddJNodeReferenceToJson(json, key, jdict_ref);
}


void iftAddJDictToJson(iftJson *json, const char *key, iftJDict *jdict) {
    iftJDict *jcopy    = iftCopyJson(jdict);
    iftAddJNodeReferenceToJson(json, key, jcopy);
}


void iftAddIntArrayToJson(iftJson *json, const char *key, const iftIntArray *iarr) {
    iftJNode *jiarr = iftCreateJIntArray(iarr);
    iftAddJNodeReferenceToJson(json, key, jiarr);
}


void iftAddFloatArrayToJson(iftJson *json, const char *key, iftFloatArray *darr) {
    iftJNode *jdarr = iftCreateJFloatArray(darr);
    iftAddJNodeReferenceToJson(json, key, jdarr);
}


void iftAddStringArrayToJson(iftJson *json, const char *key, iftStrArray *sarr) {
    iftJNode *jsarr = iftCreateJStringArray(sarr);
    iftAddJNodeReferenceToJson(json, key, jsarr);
}


void iftAddIntMatrixToJson(iftJson *json, const char *key, iftIntMatrix *iM) {
    iftJNode *jint_matrix = iftCreateJIntMatrix(iM);
    iftAddJNodeReferenceToJson(json, key, jint_matrix);
}


void iftAddDoubleMatrixToJson(iftJson *json, const char *key, iftMatrix *dM) {
    iftJNode *jdbl_matrix = iftCreateJFloatMatrix(dM);
    iftAddJNodeReferenceToJson(json, key, jdbl_matrix);
}


void iftAddStrMatrixToJson(iftJson *json, const char *key, iftStrMatrix *sM) {
    iftJNode *jstr_matrix = iftCreateJStrMatrix(sM);
    iftAddJNodeReferenceToJson(json, key, jstr_matrix);
}


void iftAddIntToJArrayNode(iftJNode *jarray, int int_val) {
    if (jarray == NULL)
        iftError("Json Array Node is NULL", "iftAddIntToJArrayNode");
    if (jarray->type != IFT_JSON_ARRAY)
        iftError("Json Node considered is not an Array (type != IFT_JSON_ARRAY)", "iftAddIntToJArrayNode");
    
    iftJNode *jint = iftCreateJInt(int_val);
    
    // empty array
    if (jarray->first_child == NULL) {
        jarray->first_child = jarray->last_child = jint;
    }
    else {
        jarray->last_child->next = jint;
        jint->prev               = jarray->last_child;
        jarray->last_child       = jint;
    }
}


void iftAddDoubleToJArrayNode(iftJNode *jarray, double dbl_val) {
    if (jarray == NULL)
        iftError("Json Array Node is NULL", "iftAddDoubleToJArrayNode");
    if (jarray->type != IFT_JSON_ARRAY)
        iftError("Json Node considered is not an Array (type != IFT_JSON_ARRAY)", "iftAddDoubleToJArrayNode");
    
    iftJNode *jdbl = iftCreateJDouble(dbl_val);
    
    // empty array
    if (jarray->first_child == NULL) {
        jarray->first_child = jarray->last_child = jdbl;
    }
    else {
        jarray->last_child->next = jdbl;
        jdbl->prev               = jarray->last_child;
        jarray->last_child       = jdbl;
    }
}


void iftAddStringToJArrayNode(iftJNode *jarray, const char *str_val) {
    if (jarray == NULL)
        iftError("Json Array Node is NULL", "iftAddStringToJArrayNode");
    if (jarray->type != IFT_JSON_ARRAY)
        iftError("Json Node considered is not an Array (type != IFT_JSON_ARRAY)", "iftAddStringToJArrayNode");
    
    iftJNode *jstr = iftCreateJString(str_val);
    
    // empty array
    if (jarray->first_child == NULL) {
        jarray->first_child = jarray->last_child = jstr;
    }
    else {
        jarray->last_child->next = jstr;
        jstr->prev               = jarray->last_child;
        jarray->last_child       = jstr;
    }
}


void iftAddIntToJArray(iftJson *json, const char *key, int int_val) {
    iftJNode *jarray = iftGetJNodeReference(json, key);
    iftAddIntToJArrayNode(jarray, int_val);
}


void iftAddDoubleToJArray(iftJson *json, const char *key, double dbl_val) {
    iftJNode *jarray = iftGetJNodeReference(json, key);
    iftAddDoubleToJArrayNode(jarray, dbl_val);
}


void iftAddStringToJArray(iftJson *json, const char *key, const char *str_val) {
    iftJNode *jarray = iftGetJNodeReference(json, key);
    iftAddStringToJArrayNode(jarray, str_val);
}


void iftDeleteJNode(iftJson *json, const char *key) {
    char keys_error_msg[1024] = ""; // used to show the key hierarchy if a key is not found
    char str[128];
    
    if (json == NULL)
        iftError("Json is NULL", "iftDeleteJNode");
    
    iftSList *key_list = iftSplitString(key, ":");
    char *elem_key     = iftRemoveSListTail(key_list); // gets the key of the element to added (it can not exist).
    
    iftJson *jdict = json;
    
    // gets the parent dictionary from the node to be deleted
    while (key_list->n != 0) {
        char *intermediate_key = iftRemoveSListHead(key_list);
        
        sprintf(str, "[%s]", intermediate_key);
        strcat(keys_error_msg, str);
        
        jdict = cJSON_GetObjectItem(jdict, intermediate_key);
        
        if (jdict == NULL)
            iftError("Invalid key: \"%s\".\nFull keys: \"%s\"", "iftDeleteJNode", intermediate_key, keys_error_msg);
        if (jdict->type != IFT_JSON_DICT)
            iftError("Intermediate key: \"%s\" is not a Dictionary.\nFull keys: \"%s\"", "iftDeleteJNode",
                     intermediate_key, keys_error_msg);
        iftFree(intermediate_key);
    }
    
    // gets the desired JNode to be deleted
    iftJNode *jnode = cJSON_GetObjectItem(jdict, elem_key);
    if (jnode == NULL) {
        sprintf(str, "[%s]", elem_key);
        strcat(keys_error_msg, str);
        iftError("The key does not exists: \"%s\".\nFull keys: \"%s\"", "iftDeleteJNode", elem_key, keys_error_msg);
    }
    else {
        if (jdict->first_child == jnode)
            jdict->first_child = jnode->next;
        if (jdict->last_child == jnode)
            jdict->last_child = jnode->prev;
        if (jnode->prev != NULL)
            jnode->prev->next = jnode->next;
        if (jnode->next != NULL)
            jnode->next->prev = jnode->prev;
        
        jnode->prev = jnode->next = NULL;
        iftDestroyJson(&jnode);
    }
    
    iftFree(elem_key);
    iftDestroySList(&key_list);
}


bool iftJsonContainKey(const char *key, iftJson *json) {
    return (cJSON_GetObjectItem(json, key) != NULL);
}
/******************************************************************/



iftDict *iftJsonToDict(iftJson *json) {
    iftDict *dict = iftCreateDict();
    
    for (iftJNode *node = json->first_child; node != NULL ; node = node->next) {
        if (node->type == IFT_JSON_DICT) {
            iftInsertIntoDict(node->key, iftJsonToDict(node), dict);
        }
        else if (node->type == IFT_JSON_ARRAY) {
            iftJNode *first_elem = node->first_child;
            
            if (first_elem->type == IFT_JSON_INT)
                iftInsertIntoDict(node->key, iftGetJIntArray(json, node->key), dict);
            else if (first_elem->type == IFT_JSON_DBL)
                iftInsertIntoDict(node->key, iftGetJDblArray(json, node->key), dict);
            else if (first_elem->type == IFT_JSON_STRING) {
                iftInsertIntoDict(node->key, iftGetJStrArray(json, node->key), dict);
            }
                // matrix
            else if (first_elem->type == IFT_JSON_ARRAY) {
                iftJNode *first_matrix_elem = first_elem->first_child;
                
                if (first_matrix_elem->type == IFT_JSON_INT)
                    iftInsertIntoDict(node->key, iftGetJIntMatrix(json, node->key), dict);
                else if (first_matrix_elem->type == IFT_JSON_DBL)
                    iftInsertIntoDict(node->key, iftGetJDblMatrix(json, node->key), dict);
                else if (first_matrix_elem->type == IFT_JSON_STRING) {
                    iftInsertIntoDict(node->key, iftGetJStrMatrix(json, node->key), dict);
                }
            }
        }
        else if (node->type == IFT_JSON_STRING) {
            iftInsertIntoDict(node->key, iftGetJString(json, node->key), dict);
        }
        else if (node->type == IFT_JSON_INT)
            iftInsertIntoDict(node->key, iftGetJInt(json, node->key), dict);
        else if (node->type == IFT_JSON_DBL)
            iftInsertIntoDict(node->key, iftGetJDouble(json, node->key), dict);
        else if ((node->type == IFT_JSON_FALSE) || (node->type == IFT_JSON_TRUE)) {
            iftInsertIntoDict(node->key, (bool) iftGetJBool(json, node->key), dict);
        }
        else if (node->type == IFT_JSON_NULL) {
            iftInsertIntoDict(node->key, iftGetJNull(json, node->key), dict);
        }
    }
    
    return dict;
}




iftDict *iftReadJson(const char *format, ...) {
    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);
    
    iftJson* json = iftReadJsonExternalLibStruct(filename);
    iftDict* dict = iftJsonToDict(json);
    iftDestroyJson(&json);
    
    return dict;
}



/**
 * @brief Transforms the pair Key and Val in strings to print after. Keys will be surrounded by quotes ""
 * @author Samuel Martins
 * @date Sep 17, 2017
 * @ingroup Dictionary
 */
void iftKeyAsQuoteAndValToString(iftKeyVal kv, char **out_key, char **out_val) {
    char *key = iftAlloc(1024, sizeof(char));
    char *val = iftAlloc(1024, sizeof(char));
    
    // prints the key
    switch(kv.key.type) {
        case IFT_UNTYPED:
            strcpy(key, "");
        case IFT_BOOL_TYPE:
            sprintf(key, "%s", kv.key.bool_val ? "true" : "false");
            break;
        case IFT_CHAR_TYPE:
        case IFT_UCHAR_TYPE:
            // since char and uchar have the same size and the fields shares the same memory location,
            // we can just use the field char_val for printing
            sprintf(key, "\'%c\'", kv.key.char_val);
            break;
        case IFT_STR_TYPE:
            iftFree(key);
            key = iftConcatStrings(3, "\"", kv.key.str_val, "\"");
            break;
        case IFT_LONG_TYPE:
            sprintf(key, "\"%ld\"", kv.key.long_val);
            break;
        case IFT_ULONG_TYPE:
            sprintf(key, "\"%lu\"", kv.key.ulong_val);
            break;
        case IFT_DBL_TYPE:
            sprintf(key, "\"%lf\"", kv.key.dbl_val);
            break;
        default:
            iftError("Invalid Type for Key", "iftKeyValToString");
    }
    
    // prints the val
    switch(kv.val.type) {
        case IFT_UNTYPED:
            strcpy(val, "");
        case IFT_BOOL_TYPE:
            sprintf(val, "%s", kv.val.bool_val ? "true" : "false");
            break;
        case IFT_CHAR_TYPE:
        case IFT_UCHAR_TYPE:
            // since char and uchar have the same size and the fields shares the same memory location,
            // we can just use the field char_val for printing
            sprintf(val, "\'%c\'", kv.val.char_val);
            break;
        case IFT_STR_TYPE:
            iftFree(val);
            val = iftConcatStrings(3, "\"", kv.val.str_val, "\"");
            break;
        case IFT_LONG_TYPE:
            sprintf(val, "%ld", kv.val.long_val);
            break;
        case IFT_ULONG_TYPE:
            sprintf(val, "%lu", kv.val.ulong_val);
            break;
        case IFT_DBL_TYPE:
            sprintf(val, "%lf", kv.val.dbl_val);
            break;
        case IFT_INT_ARRAY_TYPE:
            iftFree(val);
            val = iftIntArrayAsString(kv.val.int_array_val);
            break;
        case IFT_DBL_ARRAY_TYPE:
            iftFree(val);
            val = iftDblArrayAsString(kv.val.dbl_array_val);
            break;
        case IFT_STR_ARRAY_TYPE:
            iftFree(val);
            val = iftStrArrayAsString(kv.val.str_array_val);
            break;
        case IFT_INT_MATRIX_TYPE:
            iftFree(val);
            val = iftIntMatrixAsString(kv.val.int_matrix_val);
            break;
        case IFT_DBL_MATRIX_TYPE:
            iftFree(val);
            val = iftMatrixAsString(kv.val.dbl_matrix_val);
            break;
        case IFT_STR_MATRIX_TYPE:
            iftFree(val);
            val = iftStrMatrixAsString(kv.val.str_matrix_val);
            break;
        case IFT_PTR_TYPE:
            if (kv.val.ptr_val == 0x0) {
                sprintf(val, "null");
            }
            else sprintf(val, "%p", kv.val.ptr_val);
            break;
        default:
            iftError("Invalid Type for Val", "iftKeyValToString");
    }
    
    *out_key = key;
    *out_val = val;
}



void iftWriteJsonAux(const iftDict *dict, FILE *fp, int tab, int newline) {
    fprintf(fp, "{");
    unsigned long elem = 0;
    if (!iftIsDictEmpty(dict)) {
        for (iftKeyVal keyval = iftDictFirst(dict); iftDictValidKey(keyval); keyval = iftDictNext(dict, keyval), ++elem) {
            if (elem > 0)//do not add comma to the first element
                fprintf(fp, ",");
            for (int i = 0; i < newline; i++) fprintf(fp, "\n");
            for (int i = 0; i < tab; i++) fprintf(fp, "\t");
                if (keyval.val.type == IFT_DICT_TYPE) {
                char *dictname = iftGValToString(keyval.key);
                fprintf(fp, "%s: ", dictname);
                iftWriteJsonAux(iftGetDictVal(keyval.val), fp, tab+1, newline);
                iftFree(dictname);
            }
            else {
                char *key = NULL;
                char *val = NULL;
                iftKeyAsQuoteAndValToString(keyval, &key, &val);
                fprintf(fp, "%s: %s", key, val);
                iftFree(key);
                iftFree(val);
            }
        }
    }
    for (int i=0; i<newline; i++) fprintf(fp, "\n");
    for (int i=0; i<tab-1; i++) fprintf(fp, "\t");
    fprintf(fp, "}");
    fflush(stdout);
}


void iftWriteJson(const iftDict *dict, const char *format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);
    FILE *fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteJson", filename);
    iftWriteJsonAux(dict, fp, 1, 1);
    fprintf(fp, "\n");
    
    fclose(fp);
}


void iftWriteJsonMinified(const iftDict *dict, const char *format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);
    
    FILE *fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteJson", filename);
    
    iftWriteJsonAux(dict, fp, -1000, -1000);
    fprintf(fp, "\n");
    
    fclose(fp);
}


iftDict *iftDictFromJsonString(const char *value) {
    iftJson *json = cJSON_Parse(value);
    
    if (!json)
        iftError("JSON Error before: [%s]\n", "iftLearnBoVW.c", cJSON_GetErrorPtr());
    
    iftDict *dict = iftJsonToDict(json);
    iftDestroyJson(&json);

    return dict;
}
