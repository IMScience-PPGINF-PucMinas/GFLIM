//
// Created by Samuel Martins on 2018-12-20.
//
/**
 * @file iftJson.h
 * @brief Json parser and manipulation.
 * @author Samuel Martins
 * @date Oct 8, 2015
 * @ingroup Json
 *
 * @note An example of Json reading and writing and how to get elements from it can be found in
 * ift/demo/Miscellaneous/DataStrutuctures/iftTestJson1.c
 * @note An example of how to create a Json and add elements to it can be found in
 * ift/demo/Miscellaneous/DataStrutuctures/iftTestJson2.c
 *
 * This is an adaption and extension from the cJSON parser:
 * https://github.com/kbranigan/cJSON
 */
/*
  Copyright (c) 2009 Dave Gamble

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/


#ifndef IFT_JSON_H
#define IFT_JSON_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/Dict.h"


/**
 * @brief Adapted cJson Structure.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * This Json structure is basically a Linked List. For sake of clarity, we have defined three types to reference it:
 * --> iftJson references a complete Json, that is, the root of the Json.
 * --> iftJDict references a dictionary of the Json. Note that the root of the Json is also a dictionary.
 * --> iftJNode references a node/element of the Json. Each node has a iftJsonType that indicates the type of the element.
 */
typedef struct _ift_json {
    /** Type of Json node, see above - they also store the BOOL values. */
    int type;
    /** Key of the json node. The root of the json has key = "(null)" */
    char *key;
    /** Integer value, if type == IFT_JSON_NUMBER */
    int int_val;
    /** Double value, if type == IFT_JSON_NUMBER */
    double dbl_val;
    /** String value, if type == IFT_JSON_STRING */
    char *str_val;
    
    /** An Array or Dictionary will have a first child and a last child pointer pointing to a chain of the items in the array/dict. */
    struct _ift_json *first_child, *last_child;
    /* Next/Prev allow to walk array/dict chains. */
    struct _ift_json *next, *prev;
} iftJson, iftJDict, iftJNode;


/**
 * @brief Parses a JSON string and creates a iftDict to be returned
 *
 * @author Cesar Castelo
 * @date Mar 14, 2019
 *
 * @param value JSON string to be parsed
 */
iftDict *iftDictFromJsonString(const char *value);

/**
 * @brief Converts a JSON structure into a iftDict structure
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param json JSON structure to be converted
 */
iftDict *iftJsonToDict(iftJson *json);

/**
 * @brief Destroys an iftJson.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param json Reference to the iftJson to be destroyed.
 */
void iftDestroyJson(iftJson **json);


/**
 * @brief Reads a Json file from the disk and returns a dictionary
 * @author Samuel Martins
 * @date Sep 17, 2017
 *
 * @param json_pathname The pathname from the Json to be read.
 * @return Resulting dictionary.
 */
iftDict *iftReadJson(const char *json_path, ...);


/**
 * @brief Writes a dict to a JSON file.
 * @author Samuka Martins
 * @date Sep 17, 2017
 */
void iftWriteJson(const iftDict *dict, const char *json_path, ...);


/**
 * @brief Writes a dict to a JSON file in a minified way.
 * @author Samuka Martins
 * @date Sep 17, 2017
 */
void iftWriteJsonMinified(const iftDict *dict, const char *json_path, ...);

#ifdef __cplusplus
}
#endif

#endif //IFT_JSON_H
