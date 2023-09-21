#ifndef IFT_LABELEDSET_H_
#define IFT_LABELEDSET_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/Set.h"
#include "iftCommon.h"
#include "iftImage.h"


//! swig(extend = iftLabeledSetExt.i, destroyer = iftDestroyLabeledSet)
typedef struct ift_labeledset {
  int elem;
  int label;
  int marker;
  int handicap;
  struct ift_labeledset *next;
} iftLabeledSet;

void iftInsertLabeledSet(iftLabeledSet **S, int elem, int label);
void iftInsertLabeledSetMarkerAndHandicap(iftLabeledSet **S, int elem, int label, int marker, int handicap);
int  iftRemoveLabeledSet(iftLabeledSet **S, int *label);
void iftRemoveLabeledSetElem(iftLabeledSet **S, int elem);
void iftDestroyLabeledSet(iftLabeledSet **S);

/**
 * @brief Returns c+1 where c is the number of objects from 1 to c in the set and 0 is the background. 
 * 
 * 
 * @param S     Labeled set by subsequent integer numbers from 0 to c
 *
 * @author Falcao
 * @date Set 25th, 2017
 */

int  iftNumberOfLabels(const iftLabeledSet *S);

/**
 * @brief Returns c+1 where c is the number of markers from 1 to c in the set;
 *
 *
 * @param S     Marked set by subsequent integer numbers from 0 to c
 *
 * @author Jordão
 * @date June 30th, 2018
 */

int  iftNumberOfMarkers(const iftLabeledSet *S);

/**
 * @brief Inserts a Set with label <label> into a Labeled Set. All elements are removed from Set and it
 * is assigned to NULL.
 * 
 * @param S     Set to be inserted.
 * @param label Label of the elements from set.
 * @param T     Reference of the Labeled Set
 *
 * @author samuka
 * @date Jan 2, 2017
 */
void iftInsertSetIntoLabeledSet(iftSet **S, int label, iftLabeledSet **T);


/**
 * @brief Translates a Labeled Set by a displacement vector.
 * @param  S        Labeled Set to be translated.
 * @param  img      Image where the labeled set is.
 * @param  disp_vec Displacement vector to translate the labeled set.
 * @return          Translated Labeled Set.
 * 
 * @author Samuka Martins
 */
iftLabeledSet *iftTranslateLabeledSet(const iftLabeledSet *S, const iftImage *img, iftVector disp_vec);


/**
 * @brief  Inserts a new element in a labeled set if it is not there. 
 * @param  S        Labeled Set.
 * @param  elem     New element for insertion.
 * @param  label    Its object/region label.
 * @return True/False to indicate the success or failure of the operation.
 * @date   Sep 10, 2021
 * @author Alexandre Falcao
 */
  char iftUnionLabeledSetElem(iftLabeledSet **S, int elem, int label);
  
/**
 * @brief  Inserts a new element in a labeled set if it is not there. 
 * @param  S        Labeled Set.
 * @param  elem     New element for insertion.
 * @param  label    Its object/region label.
 * @param  marker   Its marker label.
 * @param  handicap Its handicap value.
 * @return True/False to indicate the success or failure of the operation.
 * @date   Sep 10, 2021
 * @author Alexandre Falcao
 */
  char iftUnionLabeledSetElemMarkerAndHandicap(iftLabeledSet **S, int elem, int label, int marker, int handicap);

/**
 * @brief Translates a Set by a displacement vector.
 * @param  S        Set to be translated.
 * @param  img      Image where the labeled set is.
 * @param  disp_vec Displacement vector to translate the labeled set.
 * @return          Translated Labeled Set.
 * 
 * @author Samuka Martins
 * @date Nov 16, 2017
 */
iftSet *iftTranslateSet(const iftSet *S, const iftImage *img, iftVector disp_vec);



iftLabeledSet* iftCopyLabeledSet(iftLabeledSet *s);
iftLabeledSet* iftCopyOrderedLabeledSet(iftLabeledSet *s);
  void iftPrintLabeledSet(iftLabeledSet *S);

/**
 * @brief Copies the labeled set in reverse order
 *
 * @author Thiago Vallin Spina
 *
 * @param s The input seed set.
 * @return The reversed seed set.
 */
iftLabeledSet* iftReverseLabeledSet(iftLabeledSet *s);

void iftConcatLabeledSet(iftLabeledSet **S1,iftLabeledSet **S2);
void iftRemoveSubsetLabeledSet(iftLabeledSet **S1,iftLabeledSet **S2);

int iftLabeledSetSize(const iftLabeledSet *s);

iftSet* iftLabeledSetToSet(iftLabeledSet *S, int lb);


/**
 * @brief Copies the elements in the labeled set ot and iftSet* regardless of their labels.
 *
 * @author Thiago Vallin Spina
 *
 * @param S labeled set
 * @return An iftSet* with all elements in <S>
 */
iftSet* iftLabeledSetElemsToSet(iftLabeledSet *S);
iftLabeledSet* iftCopyLabels(iftLabeledSet *S, int lb);

iftSet* iftLabeledSetMarkersToSet(iftLabeledSet *S, int marker);

/**
 * @brief Copies the set of markers into another labeled set.
 *
 * @author Thiago Vallin Spina
 * @date Jan 25, 2016
 *
 * @param S The original set of seeds.
 * @param marker The marker ID that should be copied.
 *
 * @return The filtered labeled set.
 */
iftLabeledSet* iftCopyLabeledSetMarkers(iftLabeledSet *S, int marker);

/**
 * @brief Verifies if an element belongs to the labeled seed set.
 *
 * @param S The labeled seed set.
 * @param elem The desired element.
 * @return 1 if found, 0 otherwise
 * @sa iftSetHasElement
 */
int iftLabeledSetHasElement(iftLabeledSet *S, int elem);


/**
 * @brief Insert a Labeled Set into an Image.
 * @author Samuka Martins
 * @date Nov 16, 2017
 */
void iftInsertLabeledSetIntoImage(const iftLabeledSet *S, iftImage *img);

#ifdef __cplusplus
}
#endif

#endif

