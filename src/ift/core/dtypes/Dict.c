#include "ift/core/dtypes/Dict.h"

#include "ift/core/dtypes/DblArray.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/IntMatrix.h"
#include "ift/core/dtypes/StrArray.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftMatrix.h"


/*********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Finds the Power of 2 greater than a given number @p <b>n</b>.
 * @ingroup Dictionary
 * @param n Number to be checked.
 * @param expo Reference that stores the exponent of the found power of two.
 * @return The found power of 2.
 *
 * @author Samuel Martins
 * @date Jan 17th, 2016
 */
ulong iftPowerOf2GtN(size_t n, ulong *expo) {
    ulong pow2 = 2;
    
    (*expo) = 1;
    while (pow2 < n) {
        pow2 = pow2 << 1; // (pow2 << 1) = pow2 * 2^(1)
        (*expo)++;
    }
    
    return pow2;
}


/**
 * @brief Finds the farthest Prime Number in a Range of Powers of 2.
 *
 * It assumes that lower_pow_of_2 and higher_pow_of_2 are indeed power of 2 and that the former is
 * less then the latter.
 *
 * @param lower_pow_of_2 The lower value of the range, which is a power of 2.
 * @param higher_pow_of_2 The higher value of the range, which is a power of 2.
 * @return The farthest Prime Number in the range.
 *
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 */
ulong iftFastestPrimeFromPowerOf2(ulong lower_pow_of_2, ulong higher_power_of_2) {
    // Exceptional case
    if ((lower_pow_of_2 == 2) && (higher_power_of_2 == 4))
        return 4;
    
    ulong mean = (lower_pow_of_2 + higher_power_of_2) /2;;
    ulong candidates[2];
    
    candidates[0] = mean - 1;
    candidates[1] = mean + 1;
    
    while((candidates[0] > lower_pow_of_2) && (candidates[1] < higher_power_of_2)) {
        if (iftIsPrime(candidates[0]))
            return(candidates[0]);
        if (iftIsPrime(candidates[1]))
            return(candidates[1]);
        candidates[0]--;
        candidates[1]++;
    }
    
    iftError("There is no prime number between %ld and %ld", "iftFastestPrimeFromPowerOf2",
             lower_pow_of_2, higher_power_of_2);
    
    return -1; // just to avoid warnings... Never will be reached
}


/**
 * @brief Transforms a string to a long value by using the algorithm Jenkins One At A Time Hash.
 *
 * @note http://en.wikipedia.org/wiki/Jenkins_hash_function.
 *
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 */
long iftStrToLong(const char *str) {
    if (str == NULL)
        iftError("String is NULL", "iftStrToLong");
    
    long lval = 0;
    
    for(size_t i = 0; i < strlen(str); i++) {
        lval += str[i];
        lval += (lval << 10);
        lval ^= (lval >> 6);
    }
    lval += (lval << 3);
    lval ^= (lval >> 11);
    lval += (lval << 15);
    
    return lval;
}


/**
 * @brief Transforms a key into a hash key.
 * @author Samuel Martins
 * @date Jan 6th, 2016
 * @ingroup Dictionary
 */
size_t iftGetHashKey(iftGVal key) {
    long hash_key = 0;
    
    switch(key.type) {
        case IFT_CHAR_TYPE:
            hash_key = (size_t) key.char_val;
            break;
        case IFT_UCHAR_TYPE:
            hash_key = (size_t) key.uchar_val;
            break;
        case IFT_STR_TYPE:
            hash_key = iftStrToLong(key.str_val);
            break;
        case IFT_LONG_TYPE:
            hash_key = key.long_val;
            break;
        case IFT_ULONG_TYPE:
            hash_key = (long) key.ulong_val; // it could have loss of data
            break;
        case IFT_DBL_TYPE:
            hash_key = (key.dbl_val * 10000); // truncates the abs double value to int by increasing four decimal numbers
            break;
        default:
            iftError("Invalid key datatype %s for hashing", "iftGetHashKey",
                     iftCDataTypeToString(key.type));
    }
    
    return labs(hash_key);
}


/**
 * @brief Applies a Double Hashing over a key with the datatypes: char, unsigned char, string (char*), long, and double.
 * @author Samuel Martins
 * @date Jan 6th, 2016
 * @ingroup Dictionary
 */
size_t iftHashing(iftGVal key, iftDict *dict) {
    if (dict == NULL)
        iftError("Dict is NULL", "iftHashing");
    
    size_t hash_key = iftGetHashKey(key);
    
    size_t delta = 1;
    size_t h0, h;
    h0 = h = (hash_key % dict->capacity);
    
    while ((dict->table[h].key.type != IFT_UNTYPED) && (delta <= dict->capacity)) {
        // found a table position with the same hash_key... the new value will replace this bucket
        if (iftCompareGVal(dict->table[h].key, key)) {
            return h;
        }
        else {
            delta++;
            h = (h0 + (delta - 1) * (1 + (hash_key % (dict->capacity - 2)))) % dict->capacity;
        }
    }
    
    if (delta > dict->capacity)
        return IFT_HASH_FAIL;
    
    return h; // successful
}


/**
 * @brief Transforms the pair Key and Val in strings to print after. Returns "" if there is no key val, that is, if the bucket is available.
 * @warning It does not check if the key and val parameters has sufficient space to store the strings.
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 */
void iftKeyValToString(iftKeyVal kv, char **out_key, char **out_val) {
    char *key = iftAlloc(1024, sizeof(char));
    char *val = iftAlloc(1024, sizeof(char));
    
    // prints the key
    switch(kv.key.type) {
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
            sprintf(key, "%ld", kv.key.long_val);
            break;
        case IFT_ULONG_TYPE:
            sprintf(key, "%lu", kv.key.ulong_val);
            break;
        case IFT_DBL_TYPE:
            sprintf(key, "%lf", kv.key.dbl_val);
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
            if (kv.val.ptr_val == 0x0)
                sprintf(val, "null");
            else sprintf(val, "\"%p\"", kv.val.ptr_val);
            break;
        case IFT_DICT_TYPE:
            iftFree(val);
            val = iftConcatStrings(3, "\"", iftDictToStringMinified(kv.val.dict_val), "\"");
            break;
        default:
            iftError("Invalid Type for Val", "iftKeyValToString");
    }
    
    *out_key = key;
    *out_val = val;
}


/**
 * @brief Resize a dict and rehash all elements. It's used in the function iftInsertKeyValIntoDict().
 * @author Samuel Martins
 * @date Feb 8th, 2016
 * @ingroup Dictionary
 *
 * @warning It considers that the dict is full. Then, it rehashes all elements from the dict without checking if the element really exists.
 * @warning The new number of expected elements to be inserted is (2*dict->size + 1).
 */
void iftDictResizing(iftDict **dict) {
    if (dict == NULL)
        iftError("Dict Reference is NULL", "iftDictResizing");
    if (*dict == NULL)
        iftError("Dict is NULL", "iftDictResizing");
    
    iftDict *aux = *dict;
    size_t n = 2*aux->capacity + 1; // new number of elements to be inserted
    
    // Rehashing the elements
    iftDict *new_dict = iftCreateDictWithApproxSize(n); // results in a dict with a better size to avoid colisions
    for (iftKeyVal k = iftDictFirst(aux); iftDictValidKey(k); k = iftDictNext(aux, k))
        iftInsertKeyValIntoDict(k.key, k.val, new_dict);
    
    // just deallocates the table of the original dict in order to keep the reference to it
    // ps: the strings (char*) does not deallocated, because the new_dict points to them
    iftFree(aux->table);
    
    aux->capacity = new_dict->capacity;
    aux->table    = new_dict->table;
    aux->firstKey = new_dict->firstKey;
    aux->lastKey = new_dict->lastKey;
    new_dict->table = NULL; // avoid that the destroyer deallocates the table now pointed by the original dict
    iftDestroyDict(&new_dict);
}
/*******************************************************************/


/*********************** PUBLIC FUNCTIONS *************************/
iftDict *iftCreateDict() {
    return iftCreateDictWithInitialSize(IFT_INITIAL_DICT_SIZE);
}


iftDict *iftCreateDictWithApproxSize(size_t n) {
    if ((long) n <= 0)
        iftError("The Approximate Dictionary Size %ld is <= 0", "iftCreateDict", (long) n);
    
    ulong expo;
    ulong lower_pow_of_2  = iftPowerOf2GtN(n, &expo);
    ulong higher_pow_of_2 = 1 << (expo+1); // 1 * 2^(expo+1)
    size_t size           = iftFastestPrimeFromPowerOf2(lower_pow_of_2, higher_pow_of_2);
    
    return iftCreateDictWithInitialSize(size);
}


iftDict *iftCreateDictWithInitialSize(size_t n) {
    iftDict *dict = (iftDict*) iftAlloc(1, sizeof(iftDict));
    
    dict->capacity             = n;
    dict->n_busy_buckets   = 0;
    dict->table            = (iftKeyVal*) iftAlloc(dict->capacity, sizeof(iftKeyVal));
    
    for (size_t i = 0; i < dict->capacity; i++) {
        dict->table[i].key.type = IFT_UNTYPED;
        dict->table[i].val.type = IFT_UNTYPED;
    }
    
    dict->erase_elements = false;
    dict->firstKey       = dict->lastKey = IFT_NIL;
    
    return dict;
}


void iftDestroyDict(iftDict **dict) {
    if (dict != NULL) {
        iftDict *aux = *dict;
        
        if (aux != NULL) {
            if (aux->table != NULL) {
                if (!iftIsDictEmpty(aux) && (aux->erase_elements)) {
                    for (iftKeyVal keyval = iftDictFirst(aux); iftDictValidKey(keyval); keyval = iftDictNext(aux, keyval)) {
                        iftFreeGVal(keyval.key);
                        if (keyval.val.type == IFT_DICT_TYPE)
                            iftDestroyDict(&keyval.val.dict_val);
                        else iftFreeGVal(keyval.val);
                    }
                }
                iftFree(aux->table);
            }
            iftFree(aux);
            *dict = NULL;
        }
    }
}


iftDictArray *iftCreateDictArray(size_t n) {
    iftDictArray *dict_arr = iftAlloc(1, sizeof(iftDictArray));
    dict_arr->n            = n;
    
    for (size_t i = 0; i < dict_arr->n; i++)
        dict_arr->val[i] = iftCreateDict();
    
    
    return dict_arr;
}


void iftDestroyDictArray(iftDictArray **dict_arr) {
    iftDictArray *aux = *dict_arr;
    
    if (aux != NULL) {
        for (size_t i = 0; i < aux->n; i++)
            iftDestroyDict(&aux->val[i]);
        iftFree(aux);
        *dict_arr = NULL;
    }
}


iftDict *iftCopyDict(const iftDict *dict) {
    iftDict *copy = NULL;
    
    if (dict != NULL) {
        copy                 = (iftDict*) iftAlloc(1, sizeof(iftDict));
        copy->capacity       = dict->capacity;
        copy->n_busy_buckets = dict->n_busy_buckets;
        copy->table          = (iftKeyVal*) iftAlloc(copy->capacity, sizeof(iftKeyVal));
        copy->firstKey       = dict->firstKey;
        copy->lastKey        = dict->lastKey;
        copy->erase_elements = dict->erase_elements;
        
        int i = dict->firstKey;
        iftKeyVal keyval = iftDictFirst(dict);
        while (iftDictValidKey(keyval)) {
            copy->table[i].key  = iftCopyGVal(keyval.key);
            copy->table[i].val  = iftCopyGVal(keyval.val);
            copy->table[i].prev = keyval.prev;
            copy->table[i].next = keyval.next;
            
            i = keyval.next;
            keyval = iftDictNext(dict, keyval);
        }
    }
    
    return copy;
}


bool iftIsDictEmpty(const iftDict *dict) {
    if (dict == NULL)
        iftError("Dict is NULL", "iftIsDictEmpty");
    
    return (dict->n_busy_buckets == 0);
}


bool iftIsDictFull(const iftDict *dict) {
    if (dict == NULL)
        iftError("Dict is NULL", "iftIsDictFull");
    
    return (dict->n_busy_buckets == dict->capacity);
}


bool iftDictContainGKey(iftGVal key, const iftDict *dict, size_t *key_index) {
    if (dict == NULL)
        iftError("Dict is NULL", "iftDictContainGKey");
    
    bool contain_gkey = false;
    
    if (!iftIsDictEmpty(dict)) {
        // It is basically the same hashing function of iftHashing(), however,
        // it ignores empty buckets and continue seeking the key (by double hashing) until
        // all buckets have been visited.
        // This enables that a key B rehashed due to a colision with a key A is able to reach
        // its real position even if the key A have been removed.
        size_t hash_key = iftGetHashKey(key);
        size_t delta    = 1;
        size_t h0, h;
        h0 = h = (hash_key % dict->capacity);
        
        while ((delta <= dict->capacity) && !iftCompareGVal(dict->table[h].key, key)) {
            delta++;
            h = (h0 + (delta - 1) * (1 + (hash_key % (dict->capacity - 2)))) % dict->capacity;
        }
        
        if ((delta > dict->capacity) || !iftCompareGVal(dict->table[h].key, key))
            contain_gkey = false;
        else {
            if (key_index != NULL)
                *key_index = h;
            contain_gkey = true;
        }
    }
    
    // deallocates if the input key is a string
    if (key.type == IFT_STR_TYPE) {
        iftFree(key.str_val);
        key.str_val = NULL;
    }
    
    return contain_gkey;
}


void iftPrintDictAux(const iftDict *dict, int tab, int newline) {
    if (dict != NULL) {
        printf("{");
        unsigned long elem = 0;
        if (!iftIsDictEmpty(dict)) {
            for (iftKeyVal keyval = iftDictFirst(dict); iftDictValidKey(keyval); keyval = iftDictNext(dict, keyval), ++elem) {
                
                if(elem>0)//do not add comma to the first element
                    printf(",");
                for (int i=0; i<newline; i++) printf("\n");
                for (int i=0; i<tab; i++) printf("\t");
                
                if(keyval.val.type == IFT_DICT_TYPE) {
                    char* dictname = iftGValToString(keyval.key);
                    printf("%s: ", dictname);
                    iftPrintDictAux(iftGetDictVal(keyval.val), tab+1, newline);
                    iftFree(dictname);
                }
                else {
                    char *key = NULL;
                    char *val = NULL;
                    iftKeyValToString(keyval, &key, &val);
                    printf("%s: %s", key, val);
                    iftFree(key);
                    iftFree(val);
                }
            }
        }
        for (int i=0; i<newline; i++) printf("\n");
        for (int i=0; i<tab-1; i++) printf("\t");
        printf("}\n");
        fflush(stdout);
    }
}


void iftPrintDict(const iftDict *dict) {
    iftPrintDictAux(dict, 1, 1);
}


void iftPrintDictMinified(const iftDict *dict) {
    iftPrintDictAux(dict, -1000, -1000);
}


void iftPrintDictAsArray(const iftDict *dict) {
    if (dict != NULL) {
        if (!iftIsDictEmpty(dict)) {
            //printf("%lu/%lu\n", (dict->n_busy_buckets), dict->capacity);
            
            long unsigned i=0;
            for (iftKeyVal keyval = iftDictFirst(dict) ; iftDictValidKey(keyval); keyval = iftDictNext(dict, keyval)) {
                char *key = NULL;
                char *val = NULL;
                iftKeyValToString(keyval, &key, &val);
                if(dict->n_busy_buckets < 10)
                    printf("[%1lu] %s: %s\n", i, key, val);
                else if(dict->n_busy_buckets < 100)
                    printf("[%2lu] %s: %s\n", i, key, val);
                else if(dict->n_busy_buckets < 1000)
                    printf("[%3lu] %s: %s\n", i, key, val);
                else if(dict->n_busy_buckets < 10000)
                    printf("[%4lu] %s: %s\n", i, key, val);
                else if(dict->n_busy_buckets < 100000)
                    printf("[%5lu] %s: %s\n", i, key, val);
                else
                    printf("[%lu] %s: %s\n", i, key, val);
                i++;
                iftFree(key);
                iftFree(val);
            }
        }
    }
}


bool iftInsertKeyValIntoDict(iftGVal key, iftGVal val, iftDict *dict) {
    // INSERTIONS
    if (iftIsDictFull(dict)) {
        iftDictResizing(&dict);
    }
    
    // if the key is a string, it could have intermediate dicts.sub-keys (e.g: abc:def - abc is a dict)
    if (key.type == IFT_STR_TYPE) {
        iftSList *SL = iftSplitString(key.str_val, ":");
        
        
        // string key has a sub-key, which must be an intermediate dict
        if (SL->n > 1) {
            char *sub_key = iftRemoveSListHead(SL);
            iftDestroySList(&SL);
            
            iftDict *sub_dict = NULL;
            
            // gets the sub-dict
            size_t h;
            if (iftDictContainKey(sub_key, dict, &h)) {
                if (dict->table[h].val.type == IFT_DICT_TYPE) {
                    sub_dict = dict->table[h].val.dict_val;
                }
                else {
                    sub_dict = iftCreateDict();
                    iftInsertIntoDict(sub_key, sub_dict, dict);
                }
            }
            else {
                sub_dict = iftCreateDict();
                iftInsertIntoDict(sub_key, sub_dict, dict);
            }
            
            
            char *prefix = iftConcatStrings(2, sub_key, ":");
            iftFree(sub_key);
            
            sub_key = iftRemovePrefix(key.str_val, prefix);
            iftFree(prefix);
            
            // insert the remaining sub keys
            iftFree(key.str_val);
            key = iftCreateGVal(sub_key);
            iftFree(sub_key);
            return iftInsertKeyValIntoDict(key, val, sub_dict);
        }
        else iftDestroySList(&SL);
    }
    
    size_t h = iftHashing(key, dict);
    
    if (h == IFT_HASH_FAIL)
        return false;
    
    // if the chosen bucket is empty, it increments the n_busy_buckets
    dict->n_busy_buckets += (dict->table[h].key.type == IFT_UNTYPED);
    
    // bucket is empty, update the iterator
    if (dict->table[h].key.type == IFT_UNTYPED) {
        if (dict->lastKey!=IFT_NIL)
            dict->table[dict->lastKey].next = h;
        
        dict->table[h].next = IFT_NIL;
        dict->table[h].prev = dict->lastKey;
        
        dict->lastKey = h;
        
        if (dict->firstKey == IFT_NIL)
            dict->firstKey = h;
    }
        // busy bucket - deallocate its element
    else {
        // just to avoid memory leak when overwriting existent key elements
        if (dict->table[h].key.type == IFT_STR_TYPE) {
            iftFree(dict->table[h].key.str_val);
            dict->table[h].key.str_val = NULL;
        }
        if (dict->table[h].val.type == IFT_STR_TYPE)
            iftFree(dict->table[h].val.str_val);
        else if (dict->erase_elements)
            iftFreeGVal(dict->table[h].val);
    }
    
    dict->table[h].key = key;
    dict->table[h].val = val;
    
    return true;
}


iftGVal iftGetGValFromDict(iftGVal key, const iftDict *dict, iftCDataType val_type) {
    size_t h;
    iftGVal val = {.type=IFT_UNTYPED};
    
    // since we can have compound keys (eg: "sub-dict:sub-sub-dict:test", we need to access the final/last dict to the element
    const iftDict *final_dict = dict;
    
    // if the key is a string, it could have intermediate dicts/sub-keys (e.g: abc:def - abc is a dict)
    if (key.type == IFT_STR_TYPE) {
        iftSList *SL = iftSplitString(key.str_val, ":");
        
        int n = SL->n;
        
        // finds the final/last dict to insert the element
        
        if (n > 1) {
            // scan up to the penultimate node
            for (int i = 0; i <= (n-2); i++) {
                char *sub_key = iftRemoveSListHead(SL);
                
                if (iftDictContainKey(sub_key, final_dict, &h)) {
                    if (final_dict->table[h].val.type != IFT_DICT_TYPE)
                        iftError("Sub-key \"%s\" is not a dict", "iftGetGValFromDict", sub_key);
                    final_dict = final_dict->table[h].val.dict_val;
                }
                else
                    iftError("Sub-key \"%s\" does not have a dict element", "iftGetGValFromDict", sub_key);
                
                iftFree(sub_key);
            }
        }
        char *final_key = iftRemoveSListHead(SL);
        iftDestroySList(&SL);
        
        iftFree(key.str_val);
        key = iftCreateGVal(final_key);
        iftFree(final_key);
    }
    
    if (iftDictContainKey(key, final_dict, &h))
        val = final_dict->table[h].val;
    else
        iftError("KeyError: %s (%s) does not exist", "iftGetGValFromDict", iftGValToString(key),
                 iftCDataTypeToString(key.type));
    
    if (val.type != val_type)
        iftError("Value %s from key %s is actually %s and not %s",
                 "iftGetValFromDict", iftGValToString(val), iftGValToString(key),
                 iftCDataTypeToString(val.type), iftCDataTypeToString(val_type));
    
    // deallocates if the input key is a string
    if (key.type == IFT_STR_TYPE) {
        iftFree(key.str_val);
        key.str_val = NULL;
    }
    
    return val;
}


bool iftRemoveGValFromDict(iftGVal key, iftDict *dict) {
    if (dict == NULL)
        iftError("Dict is NULL", "iftRemoveGValFromDict");
    
    size_t h = 0;
    bool removed;
    
    if (iftDictContainKey(key, dict, &h)) {
        if (dict->table[h].key.type == IFT_STR_TYPE) {
            iftFree(dict->table[h].key.str_val);
            dict->table[h].key.str_val = NULL;
        }
        if (dict->table[h].val.type == IFT_STR_TYPE) {
            iftFree(dict->table[h].val.str_val);
            dict->table[h].val.str_val = NULL;
        }
        else if (dict->erase_elements)
            iftFreeGVal(dict->table[h].val);
        
        dict->table[h].key.type = IFT_UNTYPED; // REMOVAL - the value keeps there, but it is not accessed via dict due to this flag
        dict->table[h].val.type = IFT_UNTYPED; // REMOVAL - the value keeps there, but it is not accessed via dict due to this flag
        dict->n_busy_buckets--;
        removed = true;
    } else removed = false;
    
    // deallocates if the input key is a string
    if (key.type == IFT_STR_TYPE) {
        iftFree(key.str_val);
        key.str_val = NULL;
    }
    
    if(dict->lastKey==h)
        dict->lastKey = dict->table[h].prev;
    else
        dict->table[dict->table[h].next].prev = dict->table[h].prev;
    
    if(dict->firstKey==h)
        dict->firstKey = dict->table[h].next;
    else
        dict->table[dict->table[h].prev].next = dict->table[h].next;
    
    return removed;
}


void iftUpdateValInDict(iftGVal key, iftGVal val, iftDict *dict)
{
    iftRemoveValFromDict(key, dict);
    iftInsertIntoDict(key, val, dict);
}


iftKeyVal iftDictFirst(const iftDict* dict) {
    return dict->table[dict->firstKey];
}


iftKeyVal iftDictNext(const iftDict* dict, const iftKeyVal keyval) {
    if(keyval.next!=IFT_NIL)
        return dict->table[keyval.next];
    else {
        iftKeyVal keyval;
        keyval.key.type = IFT_UNTYPED;
        keyval.val.type = IFT_UNTYPED;
        keyval.prev = IFT_NIL;
        keyval.next = IFT_NIL;
        return keyval;
    }
}


bool iftDictValidKey(const iftKeyVal keyval) {
    return keyval.key.type != IFT_UNTYPED && keyval.val.type != IFT_UNTYPED;
}


iftDict *iftIntArrayToDict(const int *arr, size_t n) {
    iftDict *dict = iftCreateDict();
    
    for (int idx = 0; idx < n; idx++)
        iftInsertIntoDict(arr[idx], idx, dict);
    
    return dict;
}


iftDict *iftMergeDicts(iftDict *dict1, iftDict *dict2)
{
    iftDict *dict3 = iftCreateDict();
    
    for(iftKeyVal k = iftDictFirst(dict1); iftDictValidKey(k); k = iftDictNext(dict1, k))
        iftInsertKeyValIntoDict(k.key, k.val, dict3);
    
    for(iftKeyVal k = iftDictFirst(dict2); iftDictValidKey(k); k = iftDictNext(dict2, k))
        iftInsertKeyValIntoDict(k.key, k.val, dict3);
    
    return dict3;
}


char *iftDictToStringMinified(const iftDict *dict) {
    char *strVal = NULL, *strAux = NULL;
    
    if (dict != NULL) {
        strVal = iftCopyString("{");
        unsigned long elem = 0;
        if (!iftIsDictEmpty(dict)) {
            for (iftKeyVal keyval = iftDictFirst(dict); iftDictValidKey(keyval); keyval = iftDictNext(dict, keyval), ++elem) {
                
                if(elem>0) {//do not add comma to the first element
                    strAux = iftConcatStrings(2, strVal, ", "); iftFree(strVal);
                    strVal = iftCopyString(strAux); iftFree(strAux);
                }
                
                if(keyval.val.type == IFT_DICT_TYPE) {
                    char* dictname = iftGValToString(keyval.key);
                    strAux = iftConcatStrings(4, strVal, dictname, ":", iftDictToStringMinified(iftGetDictVal(keyval.val))); iftFree(strVal);
                    strVal = iftCopyString(strAux); iftFree(strAux);
                    iftFree(dictname);
                }
                else {
                    char *key = NULL;
                    char *val = NULL;
                    iftKeyValToString(keyval, &key, &val);
                    strAux = iftConcatStrings(4, strVal, key, ": ", val); iftFree(strVal);
                    strVal = iftCopyString(strAux); iftFree(strAux);
                    iftFree(key);
                    iftFree(val);
                }
            }
        }
        strAux = iftConcatStrings(2, strVal, "}"); iftFree(strVal);
        strVal = iftCopyString(strAux); iftFree(strAux);
    }
    
    return strVal;
}
