#ifndef IFT_H_
#define IFT_H_

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/BMap.h"
#include "ift/core/dtypes/CharArray.h"
#include "ift/core/dtypes/Color.h"
#include "ift/core/dtypes/DblArray.h"
#include "ift/core/dtypes/DHeap.h"
#include "ift/core/dtypes/Dict.h"
#include "ift/core/dtypes/DList.h"
#include "ift/core/dtypes/FHeap.h"
#include "ift/core/dtypes/FIFO.h"
#include "ift/core/dtypes/FloatArray.h"
#include "ift/core/dtypes/FSet.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/IntQueue.h"
#include "ift/core/dtypes/IntMatrix.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "ift/core/dtypes/LIFO.h"
#include "ift/core/dtypes/List.h"
#include "ift/core/dtypes/LongArray.h"
#include "ift/core/dtypes/Matrix.h"
#include "ift/core/dtypes/Set.h"
#include "ift/core/dtypes/SList.h"
#include "ift/core/dtypes/StrArray.h"
#include "ift/core/dtypes/UCharArray.h"
#include "ift/core/dtypes/ULongArray.h"
#include "ift/core/io/CSV.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/File.h"
#include "ift/core/io/FileSet.h"
#include "ift/core/io/Json.h"
#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/CommandLineParser.h"
#include "ift/core/tools/Dialog.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/tools/OS.h"
#include "ift/core/tools/ProgressBar.h"
#include "ift/core/tools/Regex.h"
#include "ift/core/tools/String.h"
#include "ift/imgproc/basic/Histogram.h"
#include "ift/imgproc/dtypes/BoundingBox.h"
#include "ift/imgproc/dtypes/VoxelArray.h"
#include "ift/imgproc/dtypes/Roi.h"
#include "ift/medical/brain/AnomalyDetection.h"
#include "ift/medical/brain/AsymmetryAnalysis.h"
#include "ift/medical/brain/MSP.h"
#include "ift/medical/brain/PreProcessing.h"
#include "ift/medical/registration/Elastix.h"
#include "ift/medical/segm/AdaPro.h"
#include "ift/medical/segm/MALF.h"
#include "ift/medical/segm/SOSM-S.h"
#include "ift/metriclearn/MetricLearnCore.h"
#include "ift/metriclearn/LargeMargin.h"
#include "ift/metriclearn/NeighCompAnalysis.h"
#include "ift/segm/DynamicTrees.h"
#include "ift/segm/iftDynamicForest.h"
#include "ift/segm/GrabCut.h"
#include "iftIterativeOPF.h"
#include "iftFLIM.h"
#include "iftSpectrum.h"
#include "iftImageSequence.h"
#include "iftActiveLearning.h"
#include "iftAdjacency.h"
#include "iftAdjSet.h"
#include "iftBagOfFeatures.h"
#include "iftBagOfVisualWords.h"
#include "iftClassification.h"
#include "iftClustering.h"
#include "iftCommon.h"
#include "iftCompression.h"
#include "iftCompTree.h"
#include "iftCurve.h"
#include "iftDataSet.h"
#include "iftSaliency.h"
#include "iftSaliencyPriors.h"
#include "iftDeepLearning.h"
#include "iftDescriptors.h"
#include "iftDicom.h"
#include "iftDisjointSet.h"
#include "iftFiltering.h"
#include "iftFImage.h"
#include "iftFImageForest.h"
#include "iftFunctions.h"
#include "iftGenericLinkedList.h"
#include "iftGeometric.h"
#include "iftGraphics.h"
#include "iftHierarchicCluster.h"
#include "iftIGraph.h"
#include "iftImage.h"
#include "iftImageForest.h"
#include "iftImageMath.h"
#include "iftInpainting.h"
#include "iftInterpolation.h"
#include "iftKernel.h"
#include "iftKmeans.h"
#include "iftKmeans_v2.h"
#include "iftManifold.h"
#include "iftMathMorph.h"
#include "iftMatrix.h"
#include "iftMaxflow.h"
#include "iftMemory.h"
#include "iftMetrics.h"
#include "iftMImage.h"
#include "iftMSPS.h"
#include "iftMST.h"
#include "iftParamOptimizationProblems.h"
#include "iftParamOptimizer.h"
#include "iftPlane.h"
#include "iftRadiometric.h"
#include "iftReconstruction.h"
#include "iftRegion.h"
#include "iftRegistration.h"
#include "iftRepresentation.h"
#include "iftRobot.h"
#include "iftSeeds.h"
#include "iftSegmentation.h"
#include "iftSegmentationResuming.h"
#include "iftSimilarity.h"
#include "iftSlic.h"
#include "iftSort.h"
#include "iftSVM.h"
#include "iftTensor.h"
#include "iftVideo.h"
#include "iftBrainAffineRegistration.h"
#include "iftRISF.h"
#include "iftSupervoxelAtlas.h"

#ifdef IFT_GPU
    #include "iftDeepLearning.cuh"
    #include "iftMemory.cuh"
    #include "iftMatrix.cuh"
    #include "iftFLIM.cuh"
#endif

//generic data structures. Private use only in order to keep Type safety.
//#include "iftGenericVector.h"
//#include "iftGenericVectorUtil.h"
//#include "iftGenericMatrix.h"
//#include "iftGenericMatrixUtil.h"
//#include "iftArgumentList.h"
//#include "iftGenericLinkedList.h"
//#include "iftDataTransference.h"

// OS can have a library that defines a global variable as I, then we remove this definition
// for the C library compilation.
// therefore, we can declare a variable I
#undef I


/**
 * @defgroup Basic
 * @brief Structures and functions used for basic operations.
 * @{
 *
 * @defgroup DataTypes
 * @brief Definition of Basic Data Types.
 * @note Programs:
 * * @note iftTestGVal.c = It shows how to use the Generic Value (iftGVal)
 * @{ @}
 * 
 * @defgroup Dialog
 * @brief Functions used to show Dialog Messages.
 * @{ @}
 *
 * @defgroup Geometry
 * @brief Analytic Geometry operations.
 * @{ @}
 *
 * @defgroup Memory
 * @brief Functions and definitions for memory management.
 * @{ @}
 *
 * @defgroup Numeric
 * @brief Functions about Numeric operations.
 * @{ @}
 *
 * @defgroup Regex
 * @brief Regular Expression functions.
 * @note Programs:
 * * @ref iftRegexExamples.c = Examples and a quick explanation about Regex
 * * @ref iftTestRegex.c = Checks if an expression/string matches with a Regex
 * @{ @}
 *
 * @defgroup String
 * @brief String manipulation functions.
 * @note Programs:
 * * @ref iftTestSplitString.c = This program splits a string at a given position
 * * @ref iftTestReplaceString.c = It returns a copy of the string in which the occurrences of old substring have been replaced
 * with new substring.
 * @{ @}
 * 
 * @defgroup Time
 * @brief A set of functions to compute and show the runtime of programs and/or code blocks.
 * @{ @}
 * 
 * @} 
 */
/**
 * @defgroup DataStructures
 * @brief Module with Data Structure definitions and functions.
 * @{
 *
 * @defgroup BitMap
 * @brief Bit Map structure.
 * @{ @}
 * 
 * @defgroup Dictionary
 * @brief Dictionary structure, based on Hash Tables, that stores pairs of generic keys and values.
 *
 * @note Examples can be found in demo/Miscellaneous/DataStructures/iftTestDict*.c
 *
 * The allowed Datatypes for the keys are: <b>IFT_CHAR_TYPE, IFT_UCHAR_TYPE, IFT_STR_TYPE, SHORT_TYPE, INT_TYPE, IFT_LONG_TYPE, USHORT_TYPE, INT_TYPE, IFT_LONG_TYPE,
 * FLOAT_TYPE, IFT_DBL_TYPE, IFT_PTR_TYPE.</b>\n
 * The allowed Datatypes for the value are the same of the key plus <b>IFT_PTR_TYPE</b>\n.
 * Initially, empty buckets of the dictionary has key type and val type equals to IFT_UNTYPED.
 * @{ @}
 *
 * @defgroup Image
 * @brief Image definition and basic processing functions.
 * @note <b>Programs:</b>
 * * @ref iftExtractROI.c = Extract ROIs of an Image
 * * @ref iftInsertROI.c = Inserts/Copies an Image (Region Of Interest (ROI)) inside a Target Image.
 * * @ref iftExtractObject.c = Extract an Object from an Image
 * * @ref iftExtractAllObjects.c = Extract All Objects (Labels) from an Image
 * * @ref iftExtractImageROI.c = Extract the Min. Bounding Box (ROI) from an Image with all objects/voxels
 *                               (non-zero) inside it. \n
 * *
 * * @ref iftLabelImageAreaVolume.c = Computes the area/volume from each Object from an Label Image. 
 * * @ref iftSetVoxelSize.c = Overwrites the Pixel/Voxel Sizes of an Image or a Set of Images to some given Values.
 * @{ @}
 *
 * @defgroup StringDoublyLinkedList
 * @brief String Doubly Linked List definitions and functions.
 * @note Programs:
 * * @ref iftTestSList.c = It shows how to use the String Doubly Linked List.
 * @{ @}
 * 
 * @}
 */

/**
 * @defgroup MachineLearning
 * @brief Algorithms and data structures for machine learning and pattern recognition.
 * @{
 * @defgroup Supervised
 * @brief Labeled based algorithms.
 * @{ @}
 * @defgroup Algorithms not relying on labeled information.
 * @brief .
 * @{ @}
 * @}
 */

/**
 * @defgroup Dataset
 * @brief All Management function and definitions of Dataset.
 * @{
 *
 * @defgroup Metrics
 * @brief Several functions used to compute Metrics between Samples, Feature Vectors, etc, from a Dataset.
 * @{ @}
 *
 * @}
 */
/**
 * @defgroup ImageProcessing
 * @brief Image processing functions.
 * @{ 
 * 
 * @defgroup Filtering
 * @brief Several Image Filters.
 * @{ @}
 * 
 * @defgroup ImageConverters
 * @brief Several Image Converter functions.
 * @{ @}
 *
 * @defgroup ImageMetrics
 * @brief Functions that computes some metric (distance or similarity) between images.
 * @{ @}
 *
 * @defgroup Interpolation
 * @brief Image reslicing and other interpolation functions.
 * @{ @}
 *
 * 
 * @defgroup Registration
 * @brief Image Registration functions and approaches.
 * @{
 *
 * @defgroup Elastix
 * @brief Functions and definitions of the Elastix Registration approach.
 * @{ @}
 * 
 * @}
 *
 * @}
 */
/**
 * @defgroup ObjectModels
 * @brief Functions and definitions of any type of Object Model, such as Prob. Atlas, Fuzzy models, etc.
 * @{
 *
 * @defgroup Atlas
 * @brief Functions and Definitions about Probabilistic Atlases.
 * @note Programs
 * * @ref iftBuildProbAtlas.c = Builds a Probabilistic Atlas from a set of Images and their Labels.
 * @{ @}
 *
 * @}
 */
/**
 * @defgroup Utilities
 * @{
 *
 * @defgroup Allocation
 * @brief Common functions to allocate memory.
 * @{ @}
 *
 * 
 * @defgroup CommandLineParser
 * @brief Common functions to parser and manipulate Command Line Arguments from the Input.
 * @note Programs
 * * @ref iftTestCmdLineParser1.c, @ref iftTestCmdLineParser2.c = Show how to parse input command line
 *   arguments from the keyboard
 * @{ @}
 *
 * 
 * @defgroup CSV
 * @brief Generic CSV manipulations.
 * @note Programs:
 * * @ref iftTestCSV.c = It shows how to open a CSV file.
 * @{ @}
 *
 * 
 * @defgroup File
 * @brief File and path management functions.
 * @note Programs:
 * * @ref iftTestLoadDir.c = It loads a directory.
 * * @ref iftMakeDirsExamples.c = It shows examples of how to create Pathnames, Files, and Directories.
 * @{ @}
 *
 * 
 * @defgroup Messages
 * @brief Console user messages.
 * @{ @}
 *
 * @defgroup Compression
 * @brief Allows the reading and writing of compressed data into/from files.
 * @{ @}
 *
 * @}
 */

#endif
