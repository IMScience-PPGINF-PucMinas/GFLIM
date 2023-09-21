#ifndef _IFT_VIDEO_H_
#define _IFT_VIDEO_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"

#define IFT_VIDEO_FOLDER_FRAME_NZEROES 5

// This function counts the number of frame folders
int iftCountNumberOfFrameFolders(const char *path);

// This function returns the entire frame's path as <path>/<frame_id>/<name>
// using 'IFT_VIDEO_FOLDER_FRAME_NZEROES' zeroes to the left of frame_id
char* iftFramePath(const char *path, const  char *name, int frame_id);

// This function reads a frame (or other image stored in a frame folder)
iftImage* iftReadFrame(const char *path, char *name, int frame_id);

// This function creates a frame folder and writes the image into it
void iftWriteFrame(iftImage *frame, const char *path, char *name,
				   int frame_id);

// This function reads a video folder as a volumetric (color) image.
// Set start_frame = NIL and/or end_frame = NIL to read from the beginning/end
// of the video automatically. Example: iftReadVideoFolderAsVolume("GaTech/yunakim_long2/frames", NIL, NIL, "frame.ppm");
// By passing another frame name such as "label.pgm" a grayscale volume is read instead with
// the segmentation result.
iftImage* iftReadVideoFolderAsVolume(const char* folder_name, int start_frame,
									int end_frame, char *frame_name);

/**
 * @brief Loads the images in a folder as a volume.
 *
 * @author Thiago Vallin Spina
 * @date Mar 23, 2016
 *
 * @param folder_name The input folder.
 * @return The volume.
 *
 * @note This only works for 2D slices, which may be colored or not.
 */
iftImage* iftReadImageFolderAsVolume(const char* folder_name);

// This function writes a (color) volume as a video folder
void iftWriteVolumeAsVideoFolder(iftImage *video, const char *folder_name,
								 char *frame_name);

// This function first converts the video file frames into images and then
// stores them into the video folder format we use (<folder>/00001/frame.ppm)
void iftConvertVideoFileToVideoFolder(const char *video_path, const char *output_folder,
										int rotate);

// This function converts a video file into images
void iftConvertVideoFramesToImages(const char *video_path, const char *output_folder,
										const char *frame_name, const char *extension,
										int rotate);

// This function stores all frames stores as image files in a video folder
// with format <ouput_folder>/00001/frame.ppm, starting from frame id 1
void iftStoreFramesInVideoFolder(const char *input_folder, const char *input_extension,
									const char *output_folder, const char *output_filename);


#ifdef __cplusplus
}
#endif

#endif
