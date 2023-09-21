#ifndef IFT_BRAINAFFINEREGISTRATION_H_
#define IFT_BRAINAFFINEREGISTRATION_H_

#ifdef __cplusplus
extern "C" {
#endif
  
#include "iftCommon.h"
#include "iftImage.h"

  
/**
 * @brief Register the moving image to the fixed image
 *
 * @author Guilherme C.S. Ruppert
 * @date Sep 17, 2017
 *
 * @param image The fixed (reference) image
 * @param image The input image (the one to be registered)
 *
 * @return The resulting moving image registered to the fixed image domain.
*/


iftImage *iftBrainAffineRegistration(iftImage *fixed, iftImage *moving);

   
 


// Private
  


  

  
#ifdef __cplusplus
}
#endif

#endif
