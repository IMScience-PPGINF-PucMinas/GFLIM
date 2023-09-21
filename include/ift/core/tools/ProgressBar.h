//
// Created by Samuel Martins on 08/12/18.
//

#ifndef IFT_PROGRESSBAR_H
#define IFT_PROGRESSBAR_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"



/**
 * @brief Keep information of progress of a task
 * @author Peixinho
 * @date Apr, 2017
 */
typedef struct ift_progress_bar {
    float start;
    float end;
    float cur;
    int size;
} iftProgressBar;


/**
 * @brief Creates a progress bar
 * @param start starting point
 * @param end ending point
 * @return iftProgressBar object
 * @author Peixinho
 * @date Apr, 2017
 */
iftProgressBar* iftCreateProgressBar(float start, float end);


/**
 * @brief Creates a progress bar with specified size
 * @param start starting point
 * @param end ending point
 * @param size progress bar print size
 * @return iftProgressBar object
 * @author Peixinho
 * @date Apr, 2017
 */
iftProgressBar* iftCreateProgressBarSize(float start, float end, int size);


/**
 * @brief Destroys a progress bar
 * @param pbar Pointer to progress bar object
 * @return iftProgressBar object
 * @author Peixinho
 * @date Apr, 2017
 */
void iftDestroyProgressBar(iftProgressBar** pbar);


/**
 * Updates the current position of a progressbar object.
 * @param bar iftProgressBar object
 * @param step current position
 * @date Apr, 2017
 * @author Peixinho
 */
void iftProgressBarUpdate(iftProgressBar* bar, int cur);



#ifdef __cplusplus
}
#endif

#endif //IFT_PROGRESSBAR_H
