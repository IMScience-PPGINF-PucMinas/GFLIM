/**
 * @file iftDeepLearning.h
 * @brief Deep Learning structures and algorithms
 * @author Peixinho
 * @date Jun, 2016
 */


#ifndef IFT_DEEPLEARNING_H_
#define IFT_DEEPLEARNING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftMImage.h"
#include "iftAdjacency.h"
#include "iftGraphics.h"
#include "ift/core/dtypes/Color.h"
#include "iftFiltering.h"
#include "iftInterpolation.h"
#include "iftClustering.h"
#include "iftSVM.h"
#include "iftKmeans.h"
#include "iftTensor.h"
#ifdef IFT_GPU
#include "iftDeepLearning.cuh"
#endif

#define ACTIV_THRES_LIMIT  1000000
#define ALPHA_LIMIT        100.

typedef struct ift_neural_network iftNeuralNetwork;
typedef struct ift_neural_layer iftNeuralLayer;

typedef void (*iftNeuralForward) (iftNeuralNetwork*, int);
typedef void (*iftNeuralBackward) (iftNeuralNetwork*, int);
typedef void (*iftNeuralUpdate) (iftNeuralNetwork*, int);

typedef float (*iftActivationFunction) (float);

typedef enum ift_neural_layer_type { IFT_UNKNOWN_LAYER,
	IFT_MAX_POOLING_LAYER,IFT_FC_LAYER, IFT_CONV_LAYER,
	IFT_LOSS_LAYER, IFT_DATASET_LAYER, IFT_IMAGE_LAYER } iftNeuralLayerType;

float iftRelu(float f);

float iftReluDeriv(float f);

float iftLeakyRelu(float f);

float iftLeakyReluDeriv(float f);

float iftSigm(float p);

float iftSigmDeriv(float p);

float iftTanh(float f);

float iftTanhDeriv(float f);

float iftForbiddenDeriv(float f);

/**
 * @brief Neural net layer parameters. The same structure can represent all different layers. The layer currently instantied is defined in iftNeuralLayer.type parameter.
 * @author Peixinho
 * @date Nov, 2017
 */
struct ift_neural_layer {

  iftDict* imgs;

  iftTensor* weight;
  iftTensor* weightUpdate;

  iftTensor* bias;
  iftTensor* biasUpdate;

  iftTensor* data;
  iftTensor* dataMask;
  
  iftTensor* error;
  
  iftTensor* out;
  iftTensor* delta;
  
  float* buffer;
  int *ibuffer;
  long long bufferSize; // Size of either buffer in number of bytes

  int idxCounter;
  iftIntArray* idx;
  iftFileSet* fileset;
  iftDataSet* dataset;
  iftMatrix * X;
  iftMatrix * Y;
  iftSampleStatus status;
  int nsamples;
  
  iftTensor* mean;
  iftTensor* std;
  
  iftActivationFunction act;
  iftActivationFunction actDeriv;

  iftNeuralForward forward;
  iftNeuralBackward backward;
  iftNeuralUpdate update;
  
  iftNeuralLayerType type;
  bool isSharedBuffer;

  int nclasses;
  int ninput, noutput;//nbands when in convolution

  //convolutional only
  int kernelx, kernely, kernelz;
  int xout, yout, zout;
  int xdata, ydata, zdata;
  int xstride, ystride, zstride;
  int pad;
};

/**
 * @brief Neural network global parameters.
 * @date Nov, 2017
 * @author Peixinho
 */
struct ift_neural_network {
  // Number of layers
  int nlayers;

  int curLayer;

  // Ordered layer data from input to output (0 to nlayers-1)
  iftNeuralLayer** layers;

  // Learning rate for gradient descent method
  float learnRate;
  float weightDecay;

  // Minibatch size
  //   minibatch=1 is equivalent to stochastic gradient descent,
  //   minibatch=ntrain is equivalent to the batch gradient descent.
  int miniBatch;

  // Network loss error
  float error;

  // Current epoch
  //   Each epoch refers to one forward/backward pass on the entire train dataset.
  int epoch;

  // Current iteration
  //   Each iteration refers to one forward/backward pass.
  int iteration;

  // Print information about convergence on stdout
  bool verbose;

  // Expected output for the network.
  //   We aim to minimize the difference between the network output and label.
  iftTensor* label;

  // Pre-allocated memory for layers using shared buffers.
  float *ioBuffer; // for data and out
  long ioHalfSize; // Size in num of elements
  float *buffer; // cast to (int *) when using layer ibuffer
  long long bufferSize; // Size in bytes
};

/**
 * @brief Prints the neural network architecture. Usefull to check consistence of defined schema.
 * @param net
 * @date Nov, 2017
 * @author Peixinho
 */
void iftPrintNeuralNet(iftNeuralNetwork* net);

/**
 * @brief Forward all layers (from input to output) of a neural network.
 * @param net Neural net
 * @date Nov, 2017
 * @author Peixinho
 */
void iftNetForward(iftNeuralNetwork* net);

/**
 * @brief Backward all layers (from output to input) of a neural network.
 * @param net Neural net
 * @date Nov, 2017
 * @author Peixinho
 */
void iftNetBackward(iftNeuralNetwork* net);

/**
 * @brief Updates the weights of all layers, with the derivative computed in iftNetBackward().
 * @warning Should be called after iftNetBackward()
 * @param net
 */
void iftNetUpdateWeight(iftNeuralNetwork *net);

iftNeuralLayer* iftMaxPoolingLayer2D(iftNeuralNetwork* net, int nbands,
    int xsize, int ysize, int poolsize, int stride);

/**
 * @brief Creates a 2D max pooling layer optimized for testing.
 * @author Felipe Galvao
 * @date Nov, 2019
 *
 * @details Analogous to iftMaxPoolingLayer2D but removing everything
 *   not used in the forward function and using shared buffers.
 *   The function iftNetworkAllocateSharedBuffers must be called over
 *   the network before running iftNetForward.
 */
iftNeuralLayer* iftMaxPoolingLayer2DForwardOnly(iftNeuralNetwork* net,
    int nbands, int xsize, int ysize, int poolsize, int stride);

/**
 * @brief Forward function of max-pooling layers.
 * @author Peixinho
 * @param net Neural network
 * @date Nov, 2017
 * @param l
 */
void iftMaxPoolingForward(iftNeuralNetwork* net, int l);

/**
 * @brief Forward function of max-pooling layers optimized for testing.
 * @author Felipe Galvao
 * @date Nov, 2019
 *
 * @details Meant for ForwardOnly operation, see iftMaxPoolingLayer2DForwardOnly.
 */
void iftMaxPoolingForwardNoMask(iftNeuralNetwork* net, int l);

/**
 * @brief Backward (derivative) function of max-pooling layers.
 * @author Peixinho
 * @param net Neural network
 * @date Nov, 2017
 * @param l
 */
void iftMaxPoolingBackward(iftNeuralNetwork* net, int l);

/**
 * @brief Learns the optimum neural net weights,
 * with a stochastic gradient descent approach.
 * @param net Neural network.
 * @param epochs Maximum number of epochs.
 * @Peixinho
 * @date Nov, 2017
 */
void iftNeuralNetGradientDescent(iftNeuralNetwork* net, int epochs);

/**
 * @brief Creates an empty neural net layer. This is the base structure for any neural net layer.
 * @date Nov, 2017
 * @author Peixinho
 */
iftNeuralLayer* iftCreateNeuralLayer();

/**
 * @brief Destroy any neural net layer.
 * @warning This function should be able to destroy all kind of layers.
 * @param l Neural net layer.
 */
void iftDestroyNeuralLayer(iftNeuralLayer** l);

/**
 * @brief Create a neural network with the specified number of layers.
 * @param nlayers Number of layers
 * @param minibatch Number of samples considered in minibatch optimization approach.
 * @return The neural net.
 * @date Nov, 2017
 * @author Peixinho
 */
iftNeuralNetwork* iftCreateNeuralNetwork(int nlayers, int minibatch);

/**
 * @brief Destroys the neural network, and associated layers.
 * @param net Neural net to destroy.
 * @date Nov, 2017
 * @author Peixinho
 */
void iftDestroyNeuralNetwork(iftNeuralNetwork** net);

/**
 * @brief Creates a neural network using a given architecture and weights.
 * @date Sep 2019
 * @author Felipe Galvao
 *
 * TODO Add specification of how binary file is structured.
 *
 * @param json_path Json file containing the Neural Network architecture.
 * @param weights_path (Optional) Binary file containing the Neural Network weights.
 * @param fileset_path (Optional) Pre-defined set of images for the image layer.
 * @param miniBatchSize If included in json, this is optional and overrides that value.
 * @param isForwardOnly Optimize network if only NetForward will be used.
 * @return Neural network
 */
iftNeuralNetwork* iftCreateWeightedNeuralNetworkFromJson(const char *json_path, const char *weights_path, const char *fileset_path, int miniBatchSize, bool isForwardOnly); 

/**
 * @brief Creates a Neural Network using a given architecture (json file).
 * @param filename Json file containing the Neural Network architecture.
 * @return Neural network
 * @date Oct 15, 2018
 * @author Cesar Castelo
 */
iftNeuralNetwork *iftCreateNeuralNetworkFromJson(const char *filename);

/**
 * @brief A Loss layer similar to the linear discriminant analysis objective function.
 * @warning Prototype only
 * @param minibatch number of batches
 * @param ninput number of input neurons
 * @return A driscriminant loss layer
 */
iftNeuralLayer* iftDiscriminantLossLayer(int minibatch, int ninput);

/**
 * @brief Backward (derivative) function of discriminant loss layers.
 * @author Peixinho
 * @date Nov, 2017
 * @param net Neural network
 * @param l layer index
 */
void iftDiscriminantLossBackward(iftNeuralNetwork* net, int l);

/**
 * @brief Forward function for image input layer.
 * @date Nov, 2017
 * @author Peixinho
 * @param net Neural net
 * @param l layer index
 */
void iftImageForward(iftNeuralNetwork* net, int l);

/**
 * @brief Forward function for \ref iftDataSet input layer.
 * @date Nov, 2017
 * @author Peixinho
 * @param net Neural net
 * @param l layer index
 */
void iftDataSetForward(iftNeuralNetwork* net, int l);

/**
 * @brief An empty function, to define layers without a forward function. Null pattern is awesome, you should try.
 * @author Peixinho
 * @date Nov, 2017
 * @param net Neural net
 * @param l layer index
 */
void iftNoForward(iftNeuralNetwork* net, int l);

/**
 * @brief An empty function, to define layers without a backward function. Null pattern is awesome, you should try.
 * @author Peixinho
 * @date Nov, 2017
 * @param net Neural net
 * @param l layer index
 */
void iftNoBackward(iftNeuralNetwork* net, int l);

/**
 * @brief An empty function, to define layers without an update function. Null pattern is awesome, you should try.
 * @author Peixinho
 * @date Nov, 2017
 * @param net Neural net
 * @param l layer index
 */
void iftNoUpdate(iftNeuralNetwork* net, int l);

/**
 * @brief Throws error when backward layer function is used.
 * @author Felipe Galvao
 * @date Nov, 2019
 *
 * @details Ensures that ForwardOnly layers are not used for network
 *   backward. Unlike iftNoBackward usage, this is for layers that
 *   conceptually _should_ perform the operation, but were created
 *   without the appropriate data for optimization reasons.
 */
void iftForbiddenBackward(iftNeuralNetwork *net, int l);

/**
 * @brief Throws error when update layer function is used.
 * @author Felipe Galvao
 * @date Nov, 2019
 *
 * @details Ensures that ForwardOnly layers are not used for network
 *   update. Unlike iftNoUpdate usage, this is for layers that
 *   conceptually _should_ perform the operation, but were created
 *   without the appropriate data for optimization reasons.
 */
void iftForbiddenUpdate(iftNeuralNetwork *net, int l);

/**
 * @brief The backward (derivative) of a squared loss layer.
 * @param net
 * @param l
 */
void iftSquaredLossBackward(iftNeuralNetwork *net, int l);

/**
 * @brief Update the weights of a fully connected layer.
 * @param net Neural net
 * @param l layer index.
 * @author Peixinho
 * @date Nov, 2017
 */
void iftFullyConnectedUpdate(iftNeuralNetwork* net, int l);

/**
 * @brief Update the weights of a convolution layer.
 * @param net Neural net
 * @param l layer index.
 * @author Peixinho
 * @date Nov, 2017
 */
void iftConvolutionUpdate(iftNeuralNetwork* net, int l);

/**
 * @brief Forward function for image input layer.
 * @date Nov, 2017
 * @author Peixinho
 * @param net Neural net
 * @param l layer index
 */
void iftFullyConnectedForward(iftNeuralNetwork* net, int l);

/**
 * @brief Backward (derivative) function for image input layer.
 * @date Nov, 2017
 * @author Peixinho
 * @param net Neural net
 * @param l layer index
 */
void iftFullyConnectedBackward(iftNeuralNetwork* net, int l);

void iftLossBackward(iftNeuralNetwork* net, int l);

/**
 * @brief Convert windowed iteration of image to an extended matrix, where each column refers to one pixel and its neighbors
 * @sa iftCol2Img()
 * @param data_im Data array
 * @param channels Number of channels
 * @param height Data height
 * @param width Data width
 * @param ksize Number of kernels
 * @param stride Stride size
 * @param pad Padding size
 * @param data_col data array in extended column format
 * @return Output size array
 */
int iftImg2Col(float* data_im,
			   int channels,  int height,  int width,
			   int ksize,  int stride, int pad, float* data_col);

// TODO Add doc -Felipe Galvao
void iftImg2ColOptimized(float* inputImg, int nBands,  int ySize,  int xSize,
    int kSize, int stride, int pad, float* outMxBuffer);
void iftImg2ColOptimizedBatch(float* inputImg, int nBands,  int ySize,  int xSize,
    int kSize, int stride, int pad, int batchSize, int batchPos, float* outMxBuffer);

/**
 * @brief Convert windowed iteration of image to an extended matrix, where each column refers to one pixel and its neighbors
 * @sa iftImg2Col()
 * @param data_col data array in extended column format
 * @param channels Number of channels
 * @param height Data height
 * @param width Data width
 * @param ksize Number of kernels
 * @param stride Stride size
 * @param pad Padding size
 * @param data_im Ouput data array
 * @return Output size array
 */
void iftCol2Img(float* data_col,
				int channels,  int height,  int width,
				int ksize,  int stride, int pad, float* data_im);

/**
 * @brief Retrieve the index of the column extended matrix in the original image space.
 * @sa iftImg2Col()
 * @param channels Number of channels
 * @param height Data height
 * @param width Data width
 * @param ksize Number of kernels
 * @param stride Stride size
 * @param pad Padding size
 * @param idx_col Index
 * @return Output size array
 */
int iftImg2ColIdx(int channels,  int height,  int width,
				  int ksize,  int stride, int pad, int* idx_col);

/**
 * @brief Creates a squared loss layer.
 * @warning Please it should go as the last layer of a neural net.
 * @param net Neural net
 * @param ninput Number of input layers
 * @return Return the squared loss layer
 * @author Peixinho
 * @date Nov, 2017
 */
iftNeuralLayer* iftSquaredLossLayer(iftNeuralNetwork *net, int ninput);

/**
 * @brief Creates a 2D convolution layer
 * @author Peixinho
 * @date Nov, 2017
 * @param net Neural net
 * @param xsize input x size
 * @param ysize input y size
 * @param nbands number of input channels
 * @param kernelXSize kernel x size
 * @param kernelYSize kernel y size
 * @param xstride stride x size
 * @param ystride stride y size
 * @param deriv derivative of activation function
 * @param nkernels Number of convolution kernels
 * @param pad Size of 2D padding added around input
 * @param act activation function
 * @return The convolution layer
 */
iftNeuralLayer *iftConvolutionLayer2D(iftNeuralNetwork *net, int xsize,
    int ysize, int nbands, int kernelXSize, int kernelYSize, int xstride,
    int ystride, int nkernels, int pad, iftActivationFunction act,
    iftActivationFunction deriv);

/**
 * @brief Creates a 2D convolution layer optimized for testing.
 * @author Felipe Galvao
 * @date Nov, 2019
 *
 * @details Analogous to iftConvolutionLayer2D but removing everything
 *   not used in the forward function and using shared buffers.
 *   The function iftNetworkAllocateSharedBuffers must be called over
 *   the network before running iftNetForward.
 */
iftNeuralLayer *iftConvolutionLayer2DForwardOnly(iftNeuralNetwork *net,
    int xsize, int ysize, int nbands, int kernelXSize, int kernelYSize,
    int xstride, int ystride, int nkernels, int pad, iftActivationFunction act);

/**
 * @brief Creates an incomplete input image layer.
 * @author Felipe Galvao
 * @date Sep, 2019
 *
 * @details This is for cases where you can not define a fileset
 *   in advance (e.g., online processing). Note that the function
 *   iftFeedEmptyImageLayer must be called to define input image
 *   before forwarding the net for classification. Also useful
 *   to test the network on individual images.
 *
 * @param net Neural network
 * @param xSize Dimension for future images
 * @param ySize Dimension for future images
 * @param zSize Dimension for future images
 * @param nBands Dimension for future images
 * @return The input image layer
 */
iftNeuralLayer* iftEmptyImageLayer(iftNeuralNetwork* net, int xSize, int ySize, int zSize, int nBands);

/**
 * @brief Creates an empty image layer optimized for testing.
 * @author Felipe Galvao
 * @date Nov, 2019
 *
 * @details Analogous to iftEmptyImageLayer but using shared buffers.
 *   The function iftNetworkAllocateSharedBuffers must be called over
 *   the network before running iftNetForward.
 */
iftNeuralLayer* iftEmptyImageLayerForwardOnly(iftNeuralNetwork* net, int xSize, int ySize, int zSize, int nBands);

/**
 *  @brief Sets image from path when using iftEmptyImageLayer operation.
 *  @author Felipe Galvao
 *  @date Sep, 2019
 *
 *  @details Any additional normalization must be performed afterwards
 *    over net->layer[0]->out data.
 *
 *  @param net Neural Network starting with an EmptyImageLayer
 *  @param img_path Image to be fed into the network
 *  @return Remaining slots in current batch
 */
int iftFeedEmptyImageLayerImgPath(iftNeuralNetwork* net, const char *img_path);

/**
 *  @brief Sets color image from buffer when using iftEmptyImageLayer operation.
 *  @author Felipe Galvao
 *  @date Nov, 2019
 *
 *  TODO create float version to allow proper normalization beforehand. 
 *
 *  @param net Neural Network starting with an EmptyImageLayer
 *  @param R Red channel
 *  @param G Green channel
 *  @param B Blue channel
 *  @return Remaining slots in current batch
 **/
int iftFeedEmptyImageLayerColorImg(iftNeuralNetwork *net, int *R, int *G, int *B);

/**
 *  @brief Signals EmptyImageLayer functions to clear existing batch.
 *  @author Felipe Galvao
 *  @date Nov, 2019
 **/
void iftClearEmptyImageLayerBatch(iftNeuralNetwork *net);

/**
 * @brief Creates an input image layer
 * @author Peixinho
 * @date Nov, 2017
 * @param net Neural net
 * @param fileset The list of files to input
 * @return The input image layer
 */
iftNeuralLayer* iftImageLayer(iftNeuralNetwork* net, const iftFileSet* fileset);

/**
 * @brief Creates a raw input layer
 * @author Peixinho
 * @date Nov, 2017
 * @param net Neural net
 * @param X Data input matrix (nsamples x nfeats)
 * @param Y Label input matrix (nsamples x noutput)
 * @return The raw input layer
 */
iftNeuralLayer* iftRawDataLayer(iftNeuralNetwork* net, const iftMatrix* X, const iftMatrix* Y);

/**
 * @brief Creates a dataset input layer
 * @author Peixinho
 * @date Nov, 2017
 * @param net Neural net
 * @param Z Input dataset
 * @return The input dataset layer
 */
iftNeuralLayer* iftDataSetLayer(iftNeuralNetwork*, const iftDataSet* Z);

/**
 * @brief Creates a fully connected layer.
 * @param net Neural net
 * @param ninput number of input neurons
 * @param noutput number of output neurons
 * @param act Activation function
 * @param deriv Derivative of activation function
 * @return Fully connected layer
 */
iftNeuralLayer* iftFullyConnectedLayer(iftNeuralNetwork* net, int ninput,
    int noutput, iftActivationFunction act, iftActivationFunction deriv);

/**
 * @brief Creates a fully connected layer optimized for testing.
 * @author Felipe Galvao
 * @date Nov, 2019
 *
 * @details Analogous to iftFullyConnectedLayer but removing everything
 *   not used in the forward function and using shared buffers.
 *   The function iftNetworkAllocateSharedBuffers must be called over
 *   the network before running iftNetForward.
 */
iftNeuralLayer* iftFullyConnectedLayerForwardOnly(iftNeuralNetwork* net,
    int ninput, int noutput, iftActivationFunction act);

/**
 * @brief Classify images on a trained neural net
 * @param net Neural net
 * @param fileSet Images input
 * @param status samples to be classified
 */
iftTensor* iftNeuralNetClassifyImageSet(iftNeuralNetwork* net, iftFileSet* fileSet, uchar status);

/**
 * @brief Classify dataset on a trained neural net
 * @param net Neural net
 * @param Z Images input
 * @param status samples to be classified
 */
void iftNeuralNetClassifyDataSet(iftNeuralNetwork* net, iftDataSet* Z, uchar status);

/**
 * @brief Extracts as a dataset, the output of any layer.
 * @param net Neural net
 * @param l layer index
 * @return Hidden layer dataset output.
 */
iftDataSet* iftNeuralNetDataSetFeatures(iftNeuralNetwork* net, int l);

/**
 * @brief Allocation and setup for layers using shared buffers.
 * @author Felipe Galvao
 * @date Nov, 2019
 *
 * @details This function first guarantees that enough memory is
 *   allocated for the shared buffers, and then goes through each
 *   layer using shared buffers to guarantee they point to the right
 *   memory addresses.
 **/
void iftNetworkAllocateSharedBuffers(iftNeuralNetwork *net);

/**
 * @brief Convolution Network Descriptor
 * Contains the ConvNet architecture information along with the kernel weights.
 * @author Peixinho
 * @date Jun, 2016
 **/
typedef struct ift_convnetwork {
    int nlayers; // number of layers
    int nstages; // number of stages is 3*layers + 2
    iftAdjRel *input_norm_adj; // normalization adjacency of the input image
    int input_norm_adj_param; // parameter to create the normalization adjacency of the input image
    int input_xsize, input_ysize, input_zsize, input_nbands; // input image dimensions
    iftMMKernel **k_bank; // one kernel bank per layer
    int *k_bank_adj_param; // parameter to create adjacency of each kernel bank
    int *nkernels; // number of kernels per layer
    iftAdjRel **pooling_adj; // one pooling adjacency per layer
    int *pooling_adj_param; // parameters to create one pooling adjacency per layer
    int *stride; // one pooling stride per layer
    float *alpha; // one parameter of the pooling metric per layer
    iftAdjRel **norm_adj; // one normalization adjacency for the end of each layer
    int *norm_adj_param; // parameter to create one normalization adjacency for the end of each layer
    iftMatrix **img_ind;  // one image index matrix per layer for fast filtering
    int rescale;  // Output the image of the last layer with its 0- actual resolution or 1- with the resolution of the input image
    int with_weights; // Write 1-with or 0-without kernel weights */
    float *activ_thres; // one activation threshold per layer
    int *activ_option; // one option per layer
} iftConvNetwork;

/**
 * @brief Multi Scale Convolution Network.
 * Contains many ConvNets, analysing the image in different input scales.
 * @date Jun, 2016
 * @author Peixinho
 */
typedef struct ift_msconvnetwork {
    int nscales; // number of scales
    float *reduction_factor; // one reduction factor per scale
    iftAdjRel *output_norm_adj;  // normalization adjacency of the output image
    float output_norm_adj_param; // parameter to create the normalization adjacency of the output image
    iftConvNetwork **convnet; // one convolution network per scale
} iftMSConvNetwork;

/**
 * @brief Creates a Convolution Network with nlayers.
 * @param nlayers Number of layers in the ConvNet.
 * @return The Convolution Network structure.
 * @author Peixinho
 */
iftConvNetwork *iftCreateConvNetwork(int nlayers);

/**
 * @brief Destroys the Convolution Network.
 * @param convnet
 * @author Peixinho
 * @date Jun, 2016
 */
void iftDestroyConvNetwork(iftConvNetwork **convnet);

/**
 * @brief Create Multi Scale Convolution Network
 * @param nscales Number of input scales
 * @return Multi Scale Convolution Network
 * @author Peixinho
 * @date Junt, 2016
 */
iftMSConvNetwork *iftCreateMSConvNetwork(int nscales);

/**
 * @brief Destroy Multi Scale Convolution Network
 * @param msconvnet
 * @author Peixinho
 * @date Jun, 2016
 */
void iftDestroyMSConvNetwork(iftMSConvNetwork **msconvnet);

/**
 * @brief Create Convolution Network from a dictionary of parameters
 * @param params Dictionary of parameters
 * @return Convolution Network
 * @author Peixinho
 * @date Nov, 2016
 */
iftConvNetwork* iftCreateConvNetworkFromDict(const iftDict* params);

/**
 * @brief Read Convolution Network from disk.
 * The Convolution Network can be defined as a .convnet or .json (Describes the )
 * @param filename Convolution Network definition.
 * @return
 */
iftConvNetwork *iftReadConvNetwork(const char *filename);

/**
 * @brief Print the Convolution Network architecture.
 * @param convNet Convolution Network
 * @author Peixinho
 * @date Jun, 2016
 */
void iftPrintConvNetArch(iftConvNetwork *convnet);

/**
 * @brief Check if the architecture is feasible.
 * The pooling and convolution layers may cause the image to have the size reduced.
 * This funcion checks if all the hidden layers will have a valid size (greater than zero), and also checks if the ammount of memory needed to represent this network can be allocated.
 * @param convnet
 * @return True if the architecture is feasible, False otherwise
 */
bool iftValidArchitecture(iftConvNetwork *convnet);

/**
 * @brief Saves the Convolution Network in disk.
 * The ConvNet can be saved in two formats (.convnet/.json)
 * @warning .json for now is only able to save the ConvNet architecture, to save the ConvNet weights you must use .convnet format
 * @param convnet Convolution Network
 * @param filename Filepath
 */
void iftWriteConvNetwork(iftConvNetwork *convnet, const char *filename);

/**
 * @brief Appends two Convolution Networks, as one single serialized Convolution Network.
 * The second ConvNet is placed at the end of the first, taking as input the first ConvNet output.
 * @param convnet1 First ConvNet
 * @param convnet2 Second ConvNet
 * @return Joined ConvNet.
 * @author Samuel
 * @date Jun, 2016
 */
iftConvNetwork *iftMergeConvNetwork(iftConvNetwork *convnet1, iftConvNetwork *convnet2);

/**
 * @brief Read a Multi Scale Convolution Network from disk.
 * @param filename Filepath
 * @return The ConvNet.
 * @author Peixinho
 * @date Jun, 2016
 */
iftMSConvNetwork *iftReadMSConvNetwork(char *filename);

/**
 * @brief Stores a Multi Scale Convolution Network on disk.
 * @param msconvnet Multi Scale ConvNet
 * @param filename Filepath
 * @author Peixinho
 * @date Jun, 2016
 */
void iftWriteMSConvNetwork(iftMSConvNetwork *msconvnet, const char *filename);

/**
 * @brief Apply ConvNet in image.
 * The ConvNet output is a multiband image
 * @warning The high dimensional output is better suited for linear classifier. See also ::iftSVM and ::iftLogReg.
 * @param img Input image.
 * @param convnet Convolution Network
 * @return The multibandimage.
 * @author Peixinho
 * @date Jun, 2016
 */
iftMImage *iftApplyConvNetwork(iftMImage *img, iftConvNetwork *convnet);


/**
 * @brief Apply Multi Scale ConvNet in image.
 * The ConvNet output is a multiband image
 * @warning The high dimensional output is better suited for linear classifier. See also ::iftSVM and ::iftLogReg.
 * @param img Input image.
 * @param convnet Multi Scale Convolution Network
 * @return The multibandimage.
 * @date Jun, 2016
 * @author Peixinho
 */
iftMImage *iftApplyMSConvNetwork(iftMImage *img, iftMSConvNetwork *convnet);

/**
 * @brief Computes the multiband image dimensions along the hidden layers.
 * @param convnet Convolution Network
 * @param xsize
 * @param ysize
 * @param zsize
 * @param nbands
 * @return True if the architecture is valid, False otherwise. See also iftValidArchitecture()
 */
bool iftImageDimensionsAlongNetwork(iftConvNetwork *convnet, int *xsize, int *ysize, int *zsize, int *nbands);

/**
 * @brief Learns the kernel weights of a one layer Convolution Network.
 * @warning Only works with 1 layer ConvNets.
 * @param img Input image
 * @param convnet Convolution Network
 * @param nsamples Number of sampled patches for clustering
 * @param kmax_perc OPF unsup maximum k
 * @param whitening Apply, or not, whitening in the dataset
 * @author Peixinho
 * @date Jun, 2016
 */
void iftUnsupLearnKernels(iftMImage *img, iftConvNetwork *convnet, int nsamples, float kmax_perc, bool whitening);

/**
 * Learns the kernel weights of a one layer Convolution Network with KMeans clustering.
 * @warning Only works with 1 layer ConvNets.
 * @param img Input image
 * @param convnet Convolution Network
 * @param nsamples Number of sampled patches for clustering
 * @param k Number of clusters in KMeans
 * @param whitening Apply, or not, whitening in the dataset
 * @date Jun, 2016
 * @author Peixinho
 */
void iftUnsupLearnKernelsByKmeans(iftMImage *img, iftConvNetwork *convnet, int nsamples, int k, char whitening);

/**
 * Learns the kernel weights of a one layer Convolution Network with KMedoids clustering.
 * @warning Only works with 1 layer ConvNets.
 * @param img Input image
 * @param convnet Convolution Network
 * @param nsamples Number of sampled patches for clustering
 * @param k Number of clusters in KMedoids
 * @param whitening Apply, or not, whitening in the dataset
 * @author Peixinho
 * @date Jun, 2016
 */
void iftUnsupLearnKernelsByKmedoids(iftMImage *img, iftConvNetwork *convnet, int nsamples, int k, char whitening);


/**
 * Learns a multiband kernel from a dataset of image patches with unsupervised OPF clustering.
 * @param Z Patches dataset
 * @param A Patches adjacency
 * @param kmax_perc Maximum k percentage for OPF clustering
 * @param whitening Apply, or not, whitening in the dataset
 * @return The multiband kernel.
 * @author Peixinho
 * @date Jun, 2016
 */
iftMMKernel *iftUnsupLearnKernelsFromDataSet(iftDataSet *Z, iftAdjRel *A, float kmax_perc, bool whitening);

/**
 * Learns a multiband kernel from dataset of image patches with unsupervised KMedoids clustering.
 * @param Z Patches dataset
 * @param A Patches adjacency
 * @param k Number of clusters for KMedoids clustering
 * @param whitening Apply, or not, whitening in the dataset
 * @return The multiband kernel.
 * @author Peixinho
 * @date Jun, 2016
 */
iftMMKernel *iftUnsupLearnKernelsByKmedoidsFromDataSet(iftDataSet *Z, iftAdjRel *A, int k, char whitening);


/**
 * @brief Supervised Learn dataset with SVM
 * @param Z Dataset
 * @param[out] predict_sum TODO: discover what this variable returns
 * @return The hyperplanes learned
 * @author Peixinho
 * @date Jun, 2016
 */
iftSVMHyperplane **iftSupLearnKernelsFromDataSet(iftDataSet *Z, float *predict_sum);


/**
 * @brief Learn weights for a one layered Convolution Network, with randomly selected patchs from image
 * @param img Input image
 * @param convnet Convolution Network
 * @param whitening Apply, or not, whitening in dataset
 */
void iftSelectRandomKernels(iftMImage *img, iftConvNetwork *convnet, char whitening);

/**
 * @brief Creates the adjacency objects, according to kernel size in convolution layers
 * @param convnet Convolution Network
 */
void iftCreateAdjRelAlongNetwork(iftConvNetwork *convnet);

/**
 * @brief Creates the fast convolution matrices for each layer
 * @param convnet Convolution Network
 */
void iftMImageIndexMatrices(iftConvNetwork *convnet);

/**
 * @brief Creates randomly initialized kernels
 * @param convnet Convolution Network
 * @author Peixinho
 * @date Jun, 2016
 */
void iftCreateRandomKernelBanks(iftConvNetwork *convnet);

/**
 * @brief Creates a ConvNet from the first layers of a ConvNet
 * @param convnet Convolution Network
 * @param nlayers Number of layers to be copied
 * @author Peixinho
 * @date Jun, 2016
 * @return
 */
iftConvNetwork *iftCopyFirstLayersConvNetwork(iftConvNetwork *convnet, int nlayers);

/**
 * @brief Loads kernel weights from file
 * @param convnet Convolution Network
 * @param fp Pointer for file with kernel weights
 * TODO: I believe this shouldbe a private function
 * @author Peixinho
 * @date Jun, 2016
 */
void iftLoadKernelBanks(iftConvNetwork *convnet, FILE *fp);

/**
 * @brief Converts SVM hyperplanes to kernel
 * @param h SVM hyperplanes
 * @param kernel Output kernel
 * @param band_size Input image bands
 * @author Peixinho
 * @date Jun, 2016
 */
void iftConvertHyperplaneInKernelBand(iftSVMHyperplane *h, iftBand *kernel, int band_size);

/**
 * @brief Prints Convolution Network architecture
 * @param convNet Convolution Network
 * @author Peixinho
 * @date Jun, 2016
 */
void iftPrintConvNetArch(iftConvNetwork *convNet);

#ifdef __cplusplus
}
#endif

#endif
