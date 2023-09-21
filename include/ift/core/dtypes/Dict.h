/**
 * @file iftDict.h
 * @brief Definitions and functions of a Generic Dict based on Hash Tables.
 * @author Samuel Martins
 * @date Jan 17, 2016
 * @ingroup Dictionary
 *
 * - An example about how to create, insert, and print elements can be found in demo/Miscellaneous/iftTestDict1.c
 * - An example about how to get elements can be found in demo/Miscellaneous/iftTestDict2.c
 * - An example about how to remove elements can be found in demo/Miscellaneous/iftTestDict3.c
 */

#ifndef IFT_DICT_H
#define IFT_DICT_H


#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"



/**
 * @brief When creating a new dict, it has this initial size.
 * @ingroup Dictionary
 */
#define IFT_INITIAL_DICT_SIZE 103

/**
 * @brief When creating a new dict, it has this initial size.
 * @ingroup Dictionary
 */
#define IFT_HASH_FAIL -1 // "flag" to indicate that there is no available buckets in the dict


/**
 * @brief A pair of Key Value of generic datatypes used as iteration for iftDict.
 * @author Samuel Martins
 * @date Jan 6th, 2016
 * @ingroup Dictionary
 */
typedef struct ift_key_val {
    /** Key as Generic Value */
    iftGVal key;
    /** Value as Generic Value */
    iftGVal val;
    /** Previous KeyVal of the iterator */
    int prev;
    /** Next KeyVal of the iterator */
    int next;
} iftKeyVal;


/**
 * @brief A generic dictionary based on Hash Table that stores pairs of key values. All keys will have a common datatype, the same holds to all values.
 * @author Samuel Martins
 * @date Jan 6th, 2016
 * @ingroup Dictionary
 *
 * @note Its type definition in iftBasicDataType.h
 */
struct ift_dict {
    /** Size of the dict. When inserting more elements than this size, the dict is expanded with a new size. */
    size_t capacity;
    /** Number of busy buckets of the dict. */
    size_t n_busy_buckets;
    /** Dictionary. */
    iftKeyVal *table;
    /** It indicates if the Destroy Function will deallocate all elements from the dict */
    bool erase_elements;
    /** First key, for iteration (and world domination) purposes.*/
    int firstKey;
    /** Last key, for iteration (and world domination) purposes.*/
    int lastKey;
};


/**
 * @brief Array of Dictionaries. Useful to insert/get array of dicts in Json files
 * @author Samuel Martins
 * @date Sep 15, 2017
 */
typedef struct ift_dict_array {
    /** Size of Array */
    size_t n;
    /** Array of Dicts */
    iftDict **val;
} iftDictArray;


/**
 * @brief Creates a Dictionary structure, based on Hash Tables, that stores pairs of keys and values of any type.
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 *
 * @note The initial size of the dict is IFT_INITIAL_DICT_SIZE. When the dict is FULL, it will be resized automatically.
 * @note Initially, empty buckets of the dictionary has key type and val type equals to IFT_UNTYPED.
 *
 * @return The created dictionary.
 *
 * @sa iftCreateDictWithApproxSize(), iftCreateDictWithInitialSize()
 */
iftDict *iftCreateDict();


/**
 * @brief Creates a Dictionary with an approximate size.
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 *
 * Initially, empty buckets of the dictionary has key type and val type equals to IFT_UNTYPED.
 *
 * Actually, the <b>real dictionary size</b> is computed from the passed approximate size <b>n</b>,
 * and it is a bigger and better size used to minimize collisions due to the hashing.\n
 * The real dictionary size is found according to the following rules:
 * -# size >= n;
 * -# size is the farthest prime number of the powers of 2 that surround <b>n</b>;
 * -# size and (1 + (k % m)), where m = size - 2, do not have GCD (greatest common divisor);
 * -# (1 + (k % m)) > size, where m = size - 2.
 *
 * @param n The approximate size of the Dict.
 * @return The created dictionary.
 *
 * @note When the dict is FULL, it will be resized automatically.
 * @note Initially, empty buckets of the dictionary has key type and val type equals to IFT_UNTYPED.
 *
 * @sa iftCreateDict(), iftCreateDictWithInitialSize()
 */
iftDict *iftCreateDictWithApproxSize(size_t n);


/**
 * @brief Creates a Dictionary with the initial size <b>n</n>. USE iftCreateDictWithApproxSize() INSTEAD OF.
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 *
 * @param n The initial size of the Dict.
 * @return The created dictionary.
 *
 * @note When the dict is FULL, it will be resized automatically.
 * @note Initially, empty buckets of the dictionary has key type and val type equals to IFT_UNTYPED.
 *
 * @sa iftCreateDict(), iftCreateDictWithApproxSize()
 */
iftDict *iftCreateDictWithInitialSize(size_t n);


/**
 * @brief Destroys a dictionary. If the filed erase_elements is true, it will deallocate EVERYTHING dynamically allocated.
 * @author Samuel Martins
 * @date Jan 17, 2016
 * @ingroup Dictionary
 */
void iftDestroyDict(iftDict **dict);


/**
 * @brief Create an array of n dictionaries.
 * @author Samuka Martins
 * @date Sep 15, 2017
 */
iftDictArray *iftCreateDictArray(size_t n);


/**
 * @brief Destroy an array of n dictionaries.
 * @author Samuka Martins
 * @date Sep 15, 2017
 */
void iftDestroyDictArray(iftDictArray **dict_arr);


/**
 * @brief Copies a dictionary.
 * @author Samuel Martins
 * @date Feb 29, 2016
 * @ingroup Dictionary
 */
iftDict *iftCopyDict(const iftDict *dict);


/**
 * @brief Checks if the dict is empty.
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 */
bool iftIsDictEmpty(const iftDict *dict);


/**
 * @brief Checks if the dict is full.
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 */
bool iftIsDictFull(const iftDict *dict);


/**
 * @brief Checks if the dict contains the key. SUPER SLOW FUNCTION IF USED INSIDE LOOPS
 * @author Samuel Martins
 * @date Feb 16, 2016
 * @ingroup Dictionary
 *
 * When seeking the key, it ignores empty buckets and continue the search until all buckets
 * have been visited.\n
 * This enables that a key B rehashed (i.e., assigned to a new position), due to a colision with
 * a key A, is able to reach its real position even if the key A have been removed.\n
 *
 * @param KEY The key to be checked.
 * @param DICT The considered dict.
 * @param KEY_INDEX Returns the dict table index of the key, if it exists. Otherwise, the IFT_HASH_FAIL is returned.
 * @return True if the dict contains the key and false otherwise.
 *
 * @attention If the <code>key_index</code> is NULL, nothing is returned to it.
 * @warning If the GVal key is indeed a string, it will be deallocated at the end of the function.
 */
#define iftDictContainKey(KEY, DICT, KEY_INDEX) iftDictContainGKey(iftCreateGVal((KEY)), (DICT), (KEY_INDEX))


/**
 * @brief Checks if the dict contains the key (of the type GVal). USE THE MACRO FUNCTION iftDictContainKey() INSTEAD OF.
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 *
 * When seeking the key, it ignores empty buckets and continue the search until all buckets
 * have been visited.\n
 * This enables that a key B rehashed (i.e., assigned to a new position), due to a colision with
 * a key A, is able to reach its real position even if the key A have been removed.\n
 *
 * @param key The GVal key to be checked.
 * @param dict The considered dict.
 * @param key_index Returns the dict table index of the key, if it exists. Otherwise, the IFT_HASH_FAIL is returned.
 * @return True if the dict contains the key and false otherwise.
 *
 * @note USE THE MACRO FUNCTION iftDictContainKey() INSTEAD OF.
 * @warning If the <code>key_index</code> is NULL, nothing is returned to it.
 * @warning If the key is a string, it is deallocated at the end of the function.
 */
bool iftDictContainGKey(iftGVal key, const iftDict *dict, size_t *key_index);


/**
 * @brief Prints the dictionary in a fancy style.
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @ingroup Dictionary
 * @sa iftPrintDictMinified(), iftPrintDictAsArray()
 */
void iftPrintDict(const iftDict *dict);


/**
 * @brief Prints the dictionary in a minified view.
 * @author Samuel Martins
 * @date Feb 4th, 2016
 * @ingroup Dictionary
 * @sa iftPrintDict(), iftPrintDictAsArray()
 */
void iftPrintDictMinified(const iftDict *dict);


/**
 * @brief Prints the dictionary similarly to an array, by starting from the beginning index to the end.
 * @author Samuel Martins
 * @date Feb 4th, 2016
 * @ingroup Dictionary
 * @sa iftPrintDict(), iftPrintDictMinified()
 */
void iftPrintDictAsArray(const iftDict *dict);


/**
 * @brief Generic function to insertion of elements into a dictionary. All elements (KEY and VAL) are passed as REFERENCE (not a copy) to dict.
 * @author Samuel Martins
 * @date Feb 4th, 2016
 * @ingroup Dictionary
 *
 * @param KEY The key of any datatype to be inserted.
 * @param VAL The value of any datatype to be inserted.
 * @param DICT The target dictionary where the pair key value will be inserted.
 *
 * @note A datatype #GVal is created for the key and val and then both are inserted in the dictionary.
 */
#define iftInsertIntoDict(KEY,VAL,DICT) iftInsertKeyValIntoDict(iftCreateGVal(KEY), iftCreateGVal(VAL), DICT)


/**
 * @brief Insert a pair of generic key val into the dictionary. USE THE MACRO FUNCTION iftInsertIntoDict() INSTEAD OF.
 * @author Samuel Martins
 * @date Feb 4th, 2016
 * @ingroup Dictionary
 *
 * @param key The generic key to be inserted.
 * @param val The generic val to be inserted.
 * @param dict The target dictionary where the generic pair key value will be inserted.
 * @return True if the pair were inserted, False otherwise.
 *
 * @attention USE THE MACRO FUNCTION iftInsertIntoDict() INSTEAD OF.
 * @warning The allowed Datatypes for the keys are:\n
 * <b>char, unsigned char, string (char, unsigned char, string (char*), short, int, long, unsigned short,
 * unsigned int, "unsigned long, float, double.</b>
 */
bool iftInsertKeyValIntoDict(iftGVal key, iftGVal val, iftDict *dict);


/**
 * @brief Gets the REFERENCE (except the string which is passed as a COPY) of the value from the dict based on its key.
 *
 * It is possible to pass multiple keys, according to dict hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn"
 *
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup Dictionary
 *
 * @param KEY The key used to get the value.
 * @param DICT The used dictionary.
 * @return The desired value.
 *
 * @warning The value is not removed from the dict.
 * @warning If the key is a string, it is deallocated at the end.
 * @exception If the input key does not exist.
 * @{
 */
#define iftGetBoolValFromDict(KEY,DICT)     iftGetBoolVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_BOOL_TYPE))
#define iftGetCharValFromDict(KEY,DICT)     iftGetCharVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_CHAR_TYPE))
#define iftGetUCharValFromDict(KEY,DICT)    iftGetUCharVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_UCHAR_TYPE))
#define iftGetStrValFromDict(KEY,DICT)      iftGetStrVal(iftCopyGVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_STR_TYPE)))
#define iftGetConstStrValFromDict(KEY,DICT) iftGetConstStrVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_STR_TYPE))
#define iftGetLongValFromDict(KEY,DICT)     iftGetLongVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_LONG_TYPE))
#define iftGetULongValFromDict(KEY,DICT)    iftGetULongVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_ULONG_TYPE))
#define iftGetDblValFromDict(KEY,DICT)      iftGetDblVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_DBL_TYPE))
#define iftGetPtrValFromDict(KEY,DICT)      iftGetPtrVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_PTR_TYPE))
#define iftGetIntArrayFromDict(KEY, DICT)   iftGetIntArrayVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_INT_ARRAY_TYPE))
#define iftGetDblArrayFromDict(KEY, DICT)   iftGetDblArrayVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_DBL_ARRAY_TYPE))
#define iftGetStrArrayFromDict(KEY, DICT)   iftGetStrArrayVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_STR_ARRAY_TYPE))
#define iftGetIntMatrixFromDict(KEY, DICT)  iftGetIntMatrixVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_INT_MATRIX_TYPE))
#define iftGetDblMatrixFromDict(KEY, DICT)  iftGetDblMatrixVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_DBL_MATRIX_TYPE))
#define iftGetStrMatrixFromDict(KEY, DICT)  iftGetStrMatrixVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_STR_MATRIX_TYPE))
#define iftGetDictFromDict(KEY, DICT)       iftGetDictVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_DICT_TYPE))
/** @} */


/**
 * @brief Gets (without removing) the generic value (#iftGVal) from the dict based on its key, which is also an #iftGVal.\n
 * Use the Getter Macro functions (iftGet*ValFromDict()) instead of.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup Dictionary
 *
 * @param KEY The key used to get the value.
 * @param DICT The used dictionary.
 * @return The desired value.
 *
 * @attention Use the Getter Macro functions (iftGet*ValFromDict()) instead of.
 * @warning The value is not removed from the dict.
 * @warning If the key is a string, it is deallocated at the end of the function.
 * @exception If the input key does not exist.
 */
iftGVal iftGetGValFromDict(iftGVal key, const iftDict *dict, iftCDataType true_val_type);


/**
 * @brief Removes the value from the dict based on its key.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup Dictionary
 *
 * @param KEY The key used to get the value.
 * @param DICT The used dictionary.
 * @return The desired value.
 *
 * @warning If the key is a string, it is deallocated at the end of the function.
 * @exception If the input key does not exist.
 * @{
 */
#define iftRemoveValFromDict(KEY,DICT) iftRemoveGValFromDict(iftCreateGVal(KEY), DICT)
/** @} */


/**
 * @brief Removes the generic value (#iftGVal) from the dict based on its key, which is also an #iftGVal.
 * Use the Getter Macro functions (iftRemove*ValFromDict()) instead of.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup Dictionary
 *
 * @param KEY The key used to get the value.
 * @param DICT The used dictionary.
 * @return The desired value.
 *
 * @attention Use the Getter Macro functions (iftRemove*ValFromDict()) instead of.
 * @warning If the key is a string, it is deallocated in the end of the function.
 * @exception If the input key does not exist.
 */
bool iftRemoveGValFromDict(iftGVal key, iftDict *dict);


/**
 * @brief Updates the generic value (#iftGVal) from the dict based on its key, which is also an #iftGVal.
 *
 * @author Cesar Castelo
 * @date Feb 22, 2018
 * @ingroup Dictionary
 *
 * @param KEY The key used to get the value.
 * @param VAL The value to be inserted.
 * @param DICT The used dictionary.
 * @return The desired value.
 *
 * @warning This is NOT an efficient implementation (suitable only for small dictionaries).
 * It just calls the iftRemoveValFromDict and iftInsertIntoDict functions
 */
void iftUpdateValInDict(iftGVal key, iftGVal val, iftDict *dict);


/**
 * @brief Returns the first element in the dictionary.
 * @author Peixinho
 * @April, 2016
 */
iftKeyVal iftDictFirst(const iftDict* dict);


/**
 * @brief Gets the next element in the dictionary.
 * @author Peixinho
 * @date April, 2016
 */
iftKeyVal iftDictNext(const iftDict* dict, const iftKeyVal keyval);


/**
 * @brief Check if this is a valid element in the dictionary, for iteration.
 * @author Peixinho
 * @date April, 2016
 */
bool iftDictValidKey(const iftKeyVal keyval);


/**
 * @brief Convert an Integer Array to a Dict. The elements from the array will be the keys from the dict,
 * whereas theirs indices will be the corresponding values from the dict.
 *
 * @author Samuka Martins
 * @date Sep 13, 2017
 */
iftDict *iftIntArrayToDict(const int *arr, size_t n);


/**
 * @brief Merges two dictionaries and return a new one
 *
 * @author Cesar Castelo
 * @date Feb 27, 2018
 */
iftDict *iftMergeDicts(iftDict *dict1, iftDict *dict2);


/**
 * @brief Creates a string from a dictionary in a minifield view.
 * @author Cesar Castelo
 * @date March 01, 2018
 * @ingroup Dictionary
 */
char *iftDictToStringMinified(const iftDict *dict);


#ifdef __cplusplus
}
#endif


#endif //IFT_DICT_H
