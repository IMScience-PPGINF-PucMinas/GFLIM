//
// Created by Peixinho on 12/12/16.
//

#ifdef __cplusplus
extern "C" {
#endif

#include "iftDataSet.h"
#include "iftGenericMatrix.h"
#include "iftClassification.h"
//#include "tsne.h"


#ifndef IFT_IFTMANIFOLD_H
#define IFT_IFTMANIFOLD_H


typedef struct ift_tsne {
    iftDataSet* inputDataSet;
    size_t lowDimensionalSpace;
    size_t highDimensionalSpace;

    /*The perplexity is related to the number of nearest neighbors that is
     * used in other manifold learning algorithms. Larger datasets usually
     * require a larger perplexity. Consider selecting a value between 5 and 50.
     * The choice is not extremely critical since t-SNE is quite insensitive
     * to this parameter.*/
    double perplexity;
    /*This is the trade-off between speed and accuracy for Barnes-Hut T-SNE.
     * ‘theta’ is the angular size (referred to as theta in [3]) of a distant
     * node as measured from a point. If this size is below ‘angle’ then it is
     * used as a summary node of all points contained within it. This method is
     * not very sensitive to changes in this parameter in the range of 0.2 - 0.8.
     * Angle less than 0.2 has quickly increasing computation time and angle greater
     * 0.8 has quickly increasing error.*/
    double theta;
    /*Maximum number of iterations for the optimization*/
    size_t max_iter;

    // iteration at Stop lying about the P-values
    size_t stop_lying_iter;

    // iteration at momemtum changes value
    size_t mom_switch_iter;
    double startMomentum;
    double endMomentum;

    //tipical range 100-600
    double learningRate;
    /*Controls how tight natural clusters in the original space are in the embedded space and how much space
     * will be between them. For larger values, the space between natural clusters will be larger in the embedded space.
     * Again, the choice of this parameter is not very critical. If the cost function increases during initial
     * optimization, the early exaggeration factor or the learning rate might be too high. default 20
     * */
    double exageration;

    double* highDimensionPoints;
    double* lowDimensionPoints;

    bool normalizeProjection;
    bool  destroyInputDataset;

} iftTsne;


iftTsne* iftCreateTsne(iftDataSet* inputDataSet);
void iftDestroyTsne(iftTsne** tsnep);
void iftComputeTsneProjection(iftTsne* tsne);


iftMatrix *iftTSNE(iftDistanceTable *dist, unsigned int ndim, double perplexity, int maxIter);


/**
 * @brief Projects a dataset into a smaller dimension using t-SNE and returns a new dataset containing the projected data.
 * It copies the projected data in the projection matrix of the new dataset.
 * Furthermore, if copyHighDimFeats is true it copies the original (high dim) data to the new dataset's feat and data pointers.
 * If false it copies the reduced data to the new dataset's feat and data pointers.
 * @param Z Input dataset
 * @param ndim Number of dimension
 * @param perplexity Perplexity to execute t-SNE
 * @param max_iter Maximum number of iterations for t-SNE
 * @return Projected dataset
 * @author Cesar Castelo
 * @date Apr 26, 2018
 */
//! swig(newobject)
iftDataSet* iftDimReductionByTSNE(iftDataSet* Z, int ndim, double perplexity, size_t max_iter);

iftImage* iftEstimateClassifierBoundariesDecision(iftDataSet* dataSet,iftGenericClassifier* classifier);

iftImage* iftEstimateClassifierBoundariesDecision2(iftDataSet* dataSet,iftGenericClassifier* classifier);


#endif //IFT_IFTMANIFOLD_H_H


#ifdef __cplusplus
}
#endif
