//
// Created by Samuel Martins on 2018-12-20.
//

#ifndef IFT_BMAP_H
#define IFT_BMAP_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Definition of a Bitmap data struct.
 *
 * Since we the less datasize is 1 byte (char), we need to store the bits in this way. Thus, we have to allocate
 * an array of bytes, so that a single (byte) position in this bitmap can store 8 bitmap elements (8 bits).
 * Note that the less significant bits are the right-most bits from a number.
 *
 * Ex: A bitmap with a byte/char array with 2 values can store 16 bitmap elements (bits).
 *
 *             7 6 5 4 3 2 1 0 (storing order of the bitmap elements)
 * [0] =  25 = 0 0 0 1 1 0 0 1
 *
 *
 *             15 14 13 12 11 10 9  8 (storing order of the bitmap elements)
 * [1] = 100 = 0  1  1  0  0  1  0  0
 */
//! swig(destroyer = iftDestroyBMap)
typedef struct ift_bitmap {
    /** Array of chars (int8), the less datatype size (1 byte). */
    char *val;
    /** Array size: Number char elements (bytes) from bitmap (char array). */
    int nbytes;
    /** Number of available bits (bitmap elements) to be set: nbytes * 8 */
    int n;
} iftBMap;


/**
 * @brief Creates a Bitmap of n bits (bitmap elements).
 *
 * @warning that if the number of bitmap elements is not multiple of 8, it cannot allocate fraction of bytes.
 * Then, the next multiple of 8 is used what results in more available bitmap elements than thos required.
 *
 * @param n Number of bitmap elements (bits) from the array.
 * @return Created Bitmap.
 */
//! swig(newobject)
iftBMap *iftCreateBMap(int n);


/**
 * @brief Destroys a Bitmap.
 */
void iftDestroyBMap(iftBMap **bmap);


/**
 * @brief Read a Bit Map from File.
 * @author Samuel Martins
 * @date Apr 12, 2016
 * @ingroup BitMap
 *
 * A string is read from the file and is assigned at the char* in the bit map structure. \n
 * Thus, we have one byte more allocated due to the '\0', but the number of bytes from the structure
 * ignores it.
 */
//! swig(newobject)
iftBMap *iftReadBMap(const char *path);


/**
 * @brief Write a Bit Map.
 * @author Samuel Martins
 * @date Apr 12, 2016
 * @ingroup BitMap
 *
 * In order to avoid to write multiple chars in file, we build a string and write it. \n
 * For that, we need to allocate one byte more for the '\0', but nothing changes because of this.
 */
//! swig()
void iftWriteBMap(const iftBMap *bmap, const char *path);


/**
 * @brief Fills the 8 bitmap elements (bits) stored into each char/int8 elements in the bitmap with the
 * bits from the value <value>.
 *
 * @note For each char element (byte), set its 8 bits (bitmap elements) with the bits from the value <value>.
 */
//! swig()
void iftFillBMap(iftBMap *bmap, int value);


/**
 * @brief Copies a bitmap.
 */
//! swig(newobject)
iftBMap *iftCopyBMap(const iftBMap *src);


/**
 * @brief Sets value 0 to the b-th bitmap element.
 *
 * This function has a similar logic to function iftBMapSet1. See it for a better explanation.
 * The only different is we want a 8-bit mask where only the desired bit has value 0, and the remaining bits have value 1.
 * Then, we perform an AND operation between the 8-bits from the char/byte position [i] with this mask.
 *
 * Note that: ~0 = 11111111
 *
 * @param bmap Bitmap.
 * @param b Required b-th bitmap element.
 *
 * @author Alexandre Falcão, Samuel Martins
 * @date Dec 20, 2018
 */
//! swig()
static inline void iftBMapSet0(iftBMap *bmap, int b) {
    bmap->val[b >> 3] &= ((~0) ^ (1 << (b & 0x07)));
}


/**
 * @brief Sets value 1 to the b-th bitmap element.
 *
 * Given that the bitmap has bmap->n elements (bits), we need to find the position [i] in the char/int8 array where
 * the bitmap element is. For that, we have two options:
 * i = floor(b / 8), or
 * i = b >> 3 // shifting 3 bits from a number to the right equals the first option.
 *
 * Once we have found the position in the byte/char array that stores the b-th bitmap elements, we need to set 1
 * in its bit and keep the same bit in the others 7 bits.
 * For that, we want to perform a OR operation between the current bits from [i] and a bit mask with only the required bit as 1.
 *
 * We need to find the bitmap element position in the 8 bits from [i] (remember: the bitmap elements are stored from
 * right to left: see the Struct definition for a better explanation).
 * The bitmap element bis a int32 value. Its 3 less significant bits defines such position in the bits of [i].
 *
 * Ex: b = 35 = 00100011 (bits)
 * i = floor(35 / 8) = 4
 *
 * Suppose we have the following bits stored in this byte position:
 *
 *       39 38 37 36 35 34 33 32 (bitmap elements stored into this byte position)
 * [4] = 0  0  1  0  0  0  1  0
 *       0  0  0  0  1  0  0  0 (desired mask to only set 1 in the 35th bitmap element (keeping the remaining bits) after OR operation)
 *
 * Looking for the 3 less significant bits from 35, we have 011 = 3 (operation: b & 0x07).
 * The (3 + 1)th bit in the this byte position is exactly the desired mask.
 * Therefore, we create the desired mask by left-shifting x bits from the number 1 (00000001), where x is the decimal
 * number resulting from the 3 less significant bits from the required bitmap element <b>.
 * x = (b & 0x07)
 * mask = 1 << x
 *
 * Finally, we perform a OR operation between the 8bits from [i] with this mask.
 *
 * @param bmap Bitmap.
 * @param b Required b-th bitmap element.
 *
 * @author Alexandre Falcão, Samuel Martins
 * @date Dec 20, 2018
 */
//! swig()
static inline void iftBMapSet1(iftBMap *bmap, int b) {
    bmap->val[b >> 3] |= (1 << (b & 0x07));
}


/**
 * @brief Check if the b-th bitmap element has value 1.
 *
 * It has a similar logic to the functions iftBMapSet0 and iftBMapSet1.
 *
 * @param bmap Bitmap.
 * @param b Required b-th bitmap element.
 *
 * @author Alexandre Falcão, Samuel Martins
 * @date Dec 20, 2018
 */
//! swig()
static inline bool iftBMapValue(const iftBMap *bmap, int b) {
    return ((bmap->val[b >> 3] & (1 << (b & 0x07))) != 0);
}


/**
 * @brief Toogles/inverts the value (0 to 1 and viceversa) of the b-th bitmap element.
 *
 * It has a similar logic to the functions iftBMapSet0 and iftBMapSet1.
 *
 * @param bmap Bitmap.
 * @param b Required b-th bitmap element.
 *
 * @author Alexandre Falcão, Samuel Martins
 * @date Dec 20, 2018
 */
//! swig()
static inline void iftBMapToggle(iftBMap *bmap, int b) {
    bmap->val[b >> 3] ^= (1 << (b & 0x07));
}


#ifdef __cplusplus
}
#endif


#endif //IFT_BMAP_H
