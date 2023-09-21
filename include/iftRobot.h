#ifndef IFT_IFTROBOT_H
#define IFT_IFTROBOT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftImage.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "iftRepresentation.h"
#include "iftSimilarity.h"


/**
 * @brief Robot for segmentation experiments automatic evaluation.
 * @note Default functions should not be changed to receive a bunch of parameters, all extra features
 * added to the struct should be changed outside to avoid functions with unnecessary parameters
 * @author Jordão Bragantini
 * @ingroup Robot
 * @date Aug. 31, 2018
 */
//! swig(extend = iftRobotExt.i, destroyer = iftDestroyRobot)
typedef struct ift_robot {

    int stop_iter; /* default */
    int cur_iter; /* default */

    float stop_acur; /* default */
    float cur_acur; /* default */

    int max_seed_per_iter; /* default */
    bool converged; /* default */

    iftImage *orig; /* default */
    iftImage *gt_img; /* default */
    iftImage *segm; /* default */
    iftImage *error; /* default */

    iftLabeledSet* seeds; /* default */

    iftAdjRel *mrk_radius; /* default */
    iftAdjRel *limit_mrk_rad; /* default */

} iftRobot;


/**
 * @note This function should not be change to support extra parameters to iftRobot, changes should be made outside
 * or in other function to avoid, unnecessary parameters of specific tasks
 * @param orig_path path to original image
 * @param gt_path   path to original image ground truth
 * @param stop_iter maximum number of iterations
 * @param stop_acur minimum accuracy
 * @author Jordão Bragantini
 * @ingroup Robot
 * @date Aug. 31, 2018
 * @return new segmentation robot
 */
//! swig(newobject)
iftRobot *iftCreateRobot(const char *orig_path, const char *gt_path, int stop_iter, float stop_acur);

/**
 * @brief Reset segmentation information, not every parameter set previously
 * @param bot
 * @author Jordão Bragantini
 * @ingroup Robot
 * @date Sep. 14, 2018
 */
//! swig()
void iftResetRobotSegmentation(iftRobot *bot);

/**
 * @brief Destroy robot
 * @param bot
 * @author Jordão Bragantini
 * @ingroup Robot
 * @date Aug. 31, 2018
 */
void iftDestroyRobot(iftRobot **bot);

/**
 * @brief Updates robot error, accuracy and iterations
 * @param bot
 * @author Jordão Bragantini
 * @ingroup Robot
 * @date Aug. 31, 2018
 */
//! swig()
void iftRobotUpdateError(iftRobot *bot);

/**
 * @brief Append additional markers to the robot, the location is selected according to biggest error component area
 * and the additional markers are set on the Multi Scale Skeleton on the position with longest distance from the borders
 * @author Jordão Bragantini
 * @ingroup Robot
 * @date Aug. 31, 2018
 */
//! swig()
void iftRobotFindSeedsMSSkel(iftRobot *bot);

/**
 * @brief Append additional markers to the robot, the location is select according to biggest error component area
 * and the additional markers are set on the component center of mass.
 * @param bot
 * @author Jordão Bragantini
 * @ingroup Robot
 * @date Aug. 31, 2018
 */
//! swig()
void iftRobotFindSeedsCenterOfMass(iftRobot *bot);

/**
 * @brief Print robot data
 * @param bot
 * @author Jordão Bragantini
 * @ingroup Robot
 * @date Aug. 31, 2018
 */
//! swig()
void iftRobotPrintInfo(iftRobot *bot);

/**
 * @brief Returns if the robot should continue or not operating
 * @param bot
 * @author Jordão Bragantini
 * @ingroup Robot
 * @date Aug. 31, 2018
 * @return true if the robot should continue, otherwise false
 */
//! swig()
bool iftRobotContinue(iftRobot *bot);

#ifdef __cplusplus
}
#endif

#endif //IFT_IFTROBOT_H
