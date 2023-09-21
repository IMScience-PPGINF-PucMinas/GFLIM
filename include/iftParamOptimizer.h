/**
 * @file
 * @brief Parameter Optimization utilities
 * @author Peixinho
 */

#ifndef IFT_PARAMOPTIMIZER_H
#define IFT_PARAMOPTIMIZER_H

#include "ift/core/dtypes/Dict.h"
#include "iftDataSet.h"
#include "iftFitnessFunc.h"
#include "iftMetrics.h"

typedef double (*iftDictFitnessFunc) (iftDict*, iftDict *params);
typedef double (*iftClassifierFitnessFunc)   (iftDict*, iftDataSet*, iftDict *params);
typedef double (*iftDescriptorFitnessFunc)   (iftDict*, iftFileSet*, iftDict *params);

//TODO: change all the void* problem in a params dictionary

/*
 * @brief Searchs for the best score parameters among all the possible values.
 * @param opt Space search params
 * @param f Fitness function, evaluates the classifier given the specified parameters.
 * @param problem Extra parameters to an specific problem
 * @author Peixinho
 * */
iftDict* iftGridSearch(iftDict *opt, iftDictFitnessFunc f, iftDict *problem);

/*
 * @brief Searchs for the best classifier parameters in a given dataset using the iftGridSearch() method.
 * @param opt Space search params
 * @param f Fitness function, evaluates the classifier given the specified parameters.
 * @param dataset The train and test dataset, with thre proper IFT_TRAIN/IFT_TEST split.
 * @param problem Extra parameters to an specific problem
 * @author Peixinho
 * */
iftDict * iftGridSearchClassifier(iftDict *opt, iftClassifierFitnessFunc f, iftDataSet *dataset, iftSampler *sampler,
                                  iftDict *problem);

/*
 * @brief Searchs for the best descriptor parameters in a given fileset using the iftGridSearch() method.
 * @param opt Space search params
 * @param f Fitness function, evaluates the descriptor given the specified parameters.
 * @param dataset The fileset to be used in descriptor evaluation.
 * @param problem Extra parameters to an specific problem
 * @author Peixinho
 * */
iftDict * iftGridSearchDescriptor(iftDict *, iftDescriptorFitnessFunc f, iftFileSet *dataset, iftDict *problem);

/*
 * @brief Searchs for the best score parameters using random parameters trials.
 * @param opt Space search params
 * @param f Fitness function, evaluates the classifier given the specified parameters.
 * @ntrials The number of random trials.
 * @param problem Extra parameters to an specific problem
 * @author Peixinho
 * @return The best params
 * */
iftDict* iftRandomSearch(iftDict *opt, iftDictFitnessFunc f, int ntrials, iftDict *problem);

iftDict * iftRandomSearchClassifier(iftDict *, iftClassifierFitnessFunc, iftDataSet *dataset, iftSampler *sampler,
                                    int ntrials, iftDict *);
iftDict * iftRandomSearchDescriptor(iftDict *, iftDescriptorFitnessFunc fitnessFunc, iftFileSet *, int ntrials,
                                    iftDict *);

/*
 * @brief Creates a uniform (integer) search space with the given begin, end and step parameters
 * @param begin First value
 * @param end Last value
 * @param step Increment for the sequence
 * @author Cesar Castelo
 * @return iftIntArray containing the search space (parameter values)
 * */
iftIntArray *iftUniformIntSearchSpace(int begin, int end, int step);

/*
 * @brief Creates a uniform (double) search space with the given begin, end and step parameters
 * @param begin First value
 * @param end Last value
 * @param step Increment for the sequence
 * @author Cesar Castelo
 * @return iftDblArray containing the search space (parameter values)
 * */
iftDblArray *iftUniformDblSearchSpace(double begin, double end, double step);

#endif //IFT_PARAMOPTIMIZER_H
