/**
 * @file iftSVM.h
 * @brief Definitions and functions about Support Vector Machines.
 * @author Giovani Chiachia
 *
 * @note This is an encapsulation from the LibSVM, which is found in ift_root/libsvm.
 */

#ifndef _IFT_SVM_H_
#define _IFT_SVM_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef enum ift_multi_class {IFT_OVO = 0, IFT_OVA = 1} iftMultiClass;
typedef enum ift_kernel_type {IFT_LINEAR = 0, IFT_RBF = 1, IFT_PRECOMPUTED = 2} iftKernelType;

#include "iftDataSet.h"
#include "svm.h"

/**
 * @brief Redefinition of some variables from LibSVM for reasons of ease.
 * @ingroup SVM
 * @{
 */
typedef struct svm_node      svmNode;
typedef struct svm_problem   svmProblem;
typedef struct svm_parameter svmParameter;
typedef struct svm_model     svmModel;
/** @} */


/**
 * @brief SVM Hyperplane.
 * @author Samuel Martins
 * @ingroup SVM
 */
typedef struct ift_svm_hyperplane {
    /** Number of Features/Dimensions from the Hyperplane */
    int n;
    /** Feature Vector */
    float *feat;
    /** Bias from the Hyperplane */
    float bias;
} iftSVMHyperplane;


/**
 * @brief SVM structure.
 * @author Giovani Chiachia
 * @ingroup SVM
 */
typedef struct ift_svm {
    iftMultiClass multiClass;
    /** LibSVM Structure for SVM Parameters */
    svmParameter *params;
    /** The target Problem to be solved by SVM */
    svmProblem *problem;
    /** // LibSVM structure to represent SVM Models */
    svmModel **model;
    /** True labels from each SVM model */
    double *truelabel;
    /** Number of SVM Models */
    int nmodels;
    /** Data set for classifying with OVA */
    int kernelization;
    iftDataSet *Z;
} iftSVM;

/**
 * @brief Utility function to encapsulate the multiclass classification strategy
 * @param svm
 * @param Z
 * @author Peixinho
 * @date May, 2017
 */
void iftSVMClassify(iftSVM* svm, iftDataSet* Z, uchar status);

/**
 * @brief Utility function to encapsulate the multiclass training strategy
 * @param svm
 * @param Z
 * @author Peixinho
 * @date May, 2017
 */
void iftSVMTrain(iftSVM* svm, iftDataSet* Z);

/**
 * @brief Utility function to create a SVM for any parameter
 * @param kernel
 * @param multiclass
 * @param C
 * @param sigma
 * @author Peixinho
 * @return
 */
iftSVM* iftCreateSVM(iftKernelType kernel, iftMultiClass multiclass, double C, double sigma);

/////////////////////////////////////
/**
 * @brief Create a Linear Support Vector Classification (SVC).
 * @author Giovani Chiachia
 * @ingroup SVM
 *
 * @param  C Penalty parameter C of the error term.
 * @return   The created Linear SVC.
 */
iftSVM *iftCreateLinearSVC(double C);

/**
 * @brief Create Support Vector Classification (SVC) with Radial Basis Function (RBF) kernel.
 * @author Giovani Chiachia
 * @ingroup SVM
 *
 * @param  C     Penalty parameter C of the error term.
 * @param  sigma Kernel coefficient for RBF.
 * @return       The created RBF SVC.
 */
iftSVM *iftCreateRBFSVC(double C, double sigma);

/**
 * @brief Create a Linear Support Vector Classification (SVC) pre-computing the kernels.
 * @author Giovani Chiachia
 * @ingroup SVM
 *
 * @note It has a faster training step than the normal Linear SVC.
 *
 * @param  C Penalty parameter C of the error term.
 * @return   The created Linear SVC with pre-computed the kernels.
 */
iftSVM *iftCreatePreCompSVC(double C);

/**
 * @brief Create a SVM classifier according to the params described in the dictionary.
 * @param d Params dictionary. Keys => {"C":double, "kernel":string, "sigma":double}
 * @author Peixinho
 * @return
 */
iftSVM* iftCreateSVMFromDict(const iftDict* d);

/**
 * @brief Destroy a SVM.
 * @author Giovani Chiachia
 * @ingroup SVM
 *
 * @param svm SVM to be destroyed.
 */
void iftDestroySVM(iftSVM *svm);

/**
 * @brief Loads a trained SVM model from disk.
 * @author Thiago Vallin Spina
 * @ingroup SVM
 * @date May 17, 2016
 *
 * The support vector indices are also read, even though that is not done by LibSVM.
 *
 * @param filename File path to load.
 * @return         The SVM model.
 */
iftSVM *iftReadSVM(const char *pathname);

/**
 * @brief Saves a SVM model on disk for later usage.
 * @author Thiago Vallin Spina
 * @ingroup SVM
 * @date May 17, 2016
 *
 * The support vector indices are also saved, even though that is not done by LibSVM.
 *
 * @param svm      SVM model/classifier to be saved.
 * @param filename File path to save.
 */
void iftWriteSVM(const iftSVM *svm, const char *pathname);

/**
 * @brief Creates a SVM Hyperplane structure.
 * @author Samuel Martins
 * @ingroup SVM
 *
 * @param  n Number of Features/Dimensions from the Hyperplane
 * @return   The allocated/created SVM Hyperplane.
 */
iftSVMHyperplane *iftCreateSVMHyperplane(int n);

/**
 * @brief Destroy a SVM Hyperplane.
 * @author Samuel Martins
 * @ingroup SVM
 */
void iftDestroySVMHyperplane(iftSVMHyperplane **h);

/**
 * @brief Get the Normal Hyperplane (as a Sample) of the Model with index <b>model_idx</b> from a Trained Linear SVM as a Sample.
 * @author Samuel Martins
 * @ingroup SVM
 *
 * @param svm       Linear SVM (previously trained) used to get the required hyperplane.
 * @param model_idx Index from the required model whose hyperplane will be returned.
 * @param Ztrain    Dataset used to train the input linear SVM that is used to get the required hyperplane
 *                  via inner product.
 * @param out_bias  Returns by Reference (if != NULL) the bias from the required hyperplane.
 * @return          The hyperplane of the Model with index <b>model_idx</b> as a sample.
 */
iftSample iftSVMGetNormalHyperplaneAsSample(const iftSVM *svm, int model_idx, const iftDataSet *Ztrain, float *out_bias);


/**
 * @brief Get the Normal Hyperplane of the Model with index <b>model_idx</b> from a Trained Linear SVM.
 * @author Samuel Martins
 * @ingroup SVM
 *
 * @param svm       Linear SVM (previously trained) used to get the required hyperplane.
 * @param model_idx Index from the required model whose hyperplane will be returned.
 * @param Ztrain    Dataset used to train the input linear SVM that is used to get the required hyperplane
 *                  via inner product.
 * @return          The hyperplane of the Model with index <b>model_idx</b>.
 */
iftSVMHyperplane *iftSVMGetNormalHyperplane(iftSVM* svm, int model_idx, iftDataSet* Ztrain);


/**
 * @brief Get the Normal Hyperplanes of all Models from a Trained Linear SVM.
 * @author Samuel Martins
 * @ingroup SVM
 *
 * @param svm       Linear SVM (previously trained) used to get the hyperplanes.
 * @param Ztrain    Dataset used to train the input linear SVM that is used to get all hyperplanes
 *                  via inner product.
 * @return          Array with all hyperplanes.
 */
iftSVMHyperplane **iftSVMGetAllNormalHyperplanes(iftSVM *svm, iftDataSet *Ztrain);
/////////////////////////////////////


/**
 * @brief Train an One-versus-One SVM classifier from a Dataset with training samples.
 * @author Giovani Chiachia
 * @ingroup SVM
 *
 * @param svm SVM struct allocated where the classifier is going to be stored.
 * @param Z   Dataset (with training samples) used to train the SVM.
 */
void iftSVMTrainOVO(iftSVM *svm, const iftDataSet *Z, iftDataSet *Zref);

/**
 * @brief Train an One-versus-All SVM classifier from a Dataset with training samples.
 * @author Thiago Vallin Spina
 * @ingroup SVM
 * @date July 07, 2016
 *
 * @param svm SVM struct allocated where the classifier is going to be stored.
 * @param Z   Data set (with training samples) used to train the SVM.
 * @param Ztrain Data set
 * @param Zref The data set used by function iftSVMLinearClassifyOVA for classification. It may be NULL for iftSVMClassifyOVA
 */
void iftSVMTrainOVA(iftSVM *svm, const iftDataSet *Ztrain, iftDataSet *Zref);

/**
 * @brief Classify the samples with status <b>sample_status</b> from the dataset <b>Z</b> with a
 * Trained One-versus-One SVM classifier.
 * @author Giovani Chiachia
 * @ingroup SVM
 *
 * @note The predicted label of each sample is compared with its true label.
 *
 * @param svm           One-versus-One SVM used to classify the samples.
 * @param Z             Dataset (with training samples) used to train the SVM.
 * @param sample_status Status of the samples to be classified.
 * @return              Number of samples incorrectly predicted.
 */
int iftSVMClassifyOVO(const iftSVM *svm, iftDataSet *Z, uchar sample_status);

/**
 * @brief Classify the samples with status <b>sample_status</b> from the dataset <b>Z</b> with a
 * Trained One-versus-All SVM classifier.
 * @author Giovani Chiachia
 * @ingroup SVM
 *
 * @note The predicted label of each sample is compared with its true label.
 *
 * @param svm           One-versus-All SVM used to classify the samples.
 * @param Z             Dataset used to train the SVM - it uses the samples with status IFT_TRAIN.
 * @param sample_status Status of the samples to be classified.
 * @return              Number of samples incorrectly predicted.
 */



int iftSVMClassifyOVO_Probability(const iftSVM *svm, iftDataSet *Z, uchar sample_status);

/**
 * @brief Classify the samples with status <b>sample_status</b> from the dataset <b>Z</b> with a
 * Trained p-SVM classifier.
 * @author Daniel Osaku 
 * @ingroup SVM
 *
 * @note The predicted label of each sample is compared with its true label.
 *
 * @param svm           p-SVM used to classify the samples.
 * @param Z             Dataset used to train the SVM - it uses the samples with status IFT_TRAIN.
 * @param sample_status Status of the samples to be classified.
 * @return              Number of samples incorrectly predicted.
 */

int iftSVMClassifyOVA(const iftSVM *svm, iftDataSet *Z, uchar sample_status);

/**
 * @brief Optimized function to classify the samples with status <b>predStatus</b> from the dataset <b>Z</b>
 * with a trained One-versus-All LINEAR SVM classifier.
 * @author Giovani Chiachia
 * @ingroup SVM
 *
 * @warning This functions only works with Linear/Pre-Computed SVM classifiers.
 * @note <b>Prediction matrix</b> - one sample per row, one svm model per column.
 * @note The predicted label of each sample is compared with its true label.
 *
 * @param svm           One-versus-All Linear SVM used to classify the samples.
 * @param Z             Dataset to be classified.
 * @param sample_status Status of the samples to be classified.
 * @param out_pred_mat  Return by Reference (it != NULL) the matrix with the predicted values resulting
 *                      of the SVM classification.
 * @return              Number of samples incorrectly predicted.
 */
int iftSVMLinearClassifyOVA(const iftSVM *svm, iftDataSet *Z, uchar sample_status, iftMatrix **out_pred_mat);


/**
 * @brief Computes the inner product between each sample from the training set <b>Zref</b> and
 * the samples from testing set <b>Zin</b>
 * @author Giovani Chiachia
 * @ingroup SVM
 *
 * @param  Ztrain         Training set.
 * @param  Ztest          Testing set.
 * @param  kFunction      Kernel Function: [LINEAR].
 * @param  traceNormalize If true, normalizes the kernel matrix by its trace.
 * @param  ktrace         Return by reference the kernel matrix's trace.
 * @return                Kernel matrix as a Data set.
 */
iftDataSet *iftKernelizeDataSet(iftDataSet *Ztrain, iftDataSet *Ztest, int kFunction, bool traceNormalize,
                                float *ktrace);


/**
 * @todo: comment
 */
svmProblem *iftDataSetToSVMProblem(iftSVM *svc, const iftDataSet *Z);


void iftDestroyProblem(iftSVM *svm);

void iftDestroyModel(iftSVM *svm, int idx_model);

void iftDestroyProblemAndModels(iftSVM *svm);
void iftSampleToSvmSample(iftSample ift, svmNode *svm,
                          int s, int nfeats, int precomp);
void iftDestroySVMClassifier(iftSVM **svm);

#ifdef __cplusplus
}
#endif

#endif
