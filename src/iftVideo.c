#include "iftVideo.h"

#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"

int iftCountNumberOfFrameFolders(const char *path)
{
	int file_count = 0;

#ifdef __linux
	DIR * dirp;
	struct dirent * entry;

	dirp = opendir(path); /* There should be error handling after this */
	if(dirp != NULL)
	{
		while ((entry = readdir(dirp)) != NULL) {
			if (entry->d_type == DT_DIR) { /* If the entry is a regular file */
				file_count++;
			}
		}

		closedir(dirp);
	}

	return file_count-2; // eliminating direcotires '.' and '..'
#else
	char command[IFT_STR_DEFAULT_SIZE], frame[IFT_STR_DEFAULT_SIZE];
	FILE *f = NULL;

	sprintf(command, "ls -v %s > _tmp_.txt", path);
	system(command);
	
	f = fopen("_tmp_.txt", "r");
	while(!feof(f)){
	  fscanf(f, "%s", frame);
	  file_count++;
	}
	fclose(f);
	file_count--; // Removing '\n'
	
	system("rm -f _tmp_.txt");

	return file_count;
#endif
}

char* iftFramePath(const char *folder, const  char *name, int frame_id)
{
	char path[IFT_STR_DEFAULT_SIZE], *output_path = NULL;

	sprintf(path, "%s/%0*d/%s", folder, IFT_VIDEO_FOLDER_FRAME_NZEROES,
			frame_id, name);

	output_path = (char*)iftAlloc(strlen(path), sizeof(char));

	strcpy(output_path, path);

	return output_path;
}

iftImage* iftReadFrame(const char *path, char *name, int frame_id)
{
	char *filename = iftFramePath(path, name, frame_id);
	iftImage *frame = NULL;

	frame = iftReadImageByExt(filename);

	if(filename == NULL)
        iftError("Invalid frame file format", "iftReadFrame");

	iftFree(filename);

	return frame;
}

void iftWriteFrame(iftImage *frame, const char *path, char *name,
				   int frame_id)
{
	int success = 0;
	char *pos, ext[10];
	char *filename = iftFramePath(path, name, frame_id);
	char *folder_name = iftFramePath(path,"", frame_id);
	char cmd[16+strlen(folder_name)];

	sprintf(cmd,"mkdir -p %s", folder_name);

	if(system(cmd) == -1) iftError("Command error", "iftWriteFrame");

	pos = strrchr(name,'.');
	pos++;
	sscanf(pos,"%s",ext);

	if(!iftIs3DImage(frame))
	{
		if(strcmp(ext,"ppm") == 0)
		{
			if(iftIsColorImage(frame))
			{
				iftWriteImageP6(frame, filename);
				success = 1;
			}
		}
		else if(strcmp(ext,"pgm") == 0)
		{
			if(iftMaximumValue(frame) > 255)
				iftWriteImageP2(frame, filename);
			else
				iftWriteImageP5(frame, filename);
			success = 1;
		}
	}

	if(success == 0)
        iftError("Invalid frame file format", "iftWriteFrame");

	iftFree(filename);
	iftFree(folder_name);
}

iftImage* iftReadVideoFolderAsVolume(const char* folder_name, int start_frame,
									 int end_frame, char *frame_name)
{
	int i;
	char *pos, ext[10];
	iftImage *volume = NULL, *frame = NULL;

	if(start_frame < 0)
		start_frame = 1;

	if(end_frame < 0)
		end_frame = iftCountNumberOfFrameFolders(folder_name);

	if(end_frame < start_frame)
	{
        iftError("Start frame must be smaller than end frame", "iftReadVideoFolderAsVolume");
	}


	frame = iftReadFrame(folder_name, frame_name, start_frame);

	pos = strrchr(frame_name,'.') + 1;
	sscanf(pos,"%s",ext);

	if(strcmp(ext,"ppm") == 0)
	{
		volume = iftCreateColorImage(frame->xsize, frame->ysize, end_frame - start_frame+1, 8);
	}
	else if(strcmp(ext,"pgm") == 0)
	{
		volume = iftCreateImage(frame->xsize, frame->ysize, end_frame - start_frame+1);
	}
	else
	{
        iftError("Invalid frame file format", "iftReadVideoFolderAsVolume");
	}

	iftPutXYSlice(volume, frame, 0);
	iftDestroyImage(&frame);

	for(i = start_frame+1; i <= end_frame; i++)
	{
		frame = iftReadFrame(folder_name, frame_name, i);

		iftPutXYSlice(volume, frame, i-start_frame);

		iftDestroyImage(&frame);
	}

	return volume;
}

iftImage* iftReadImageFolderAsVolume(const char* folder_name)
{
	int z;
	iftImage *volume = NULL;
	iftFileSet *files = NULL;

	files = iftLoadFileSetFromDirOrCSV(folder_name, 0, true);

	// Verifying if all
	for(z = 0; z < files->n; z++)
	{
		const char *path = files->files[z]->path;
		if(!iftIsImageFile(path) || iftCompareStrings(iftFileExt(path),".scn")) {
            iftError("File %s cannot be read as an image for constructing a 3D volume!", "iftReadImageFolderAsVolume",
                     path);
		}
	}

	for(z = 0; z < files->n; z++)
	{
		const char *path = files->files[z]->path;
		iftImage *slice = iftReadImageByExt(path);

		if(volume == NULL) {
			if(iftIsColorImage(slice))
				volume = iftCreateColorImage(slice->xsize, slice->ysize, files->n, iftImageDepth(slice));
			else
				volume = iftCreateImage(slice->xsize, slice->ysize, files->n);
		}

		if(slice->xsize != volume->xsize || slice->ysize != volume->ysize)
            iftError("The XY dimensions of the current image being loaded (%s) differs from the first image selected" \
            "for initializing the volume! First image: %d %d %d, current image: %d %d %d\n",
                     "iftReadImageFolderAsVolume",
                     path, volume->xsize, volume->ysize, volume->zsize,
                     slice->xsize, slice->ysize, slice->zsize);

		if(iftIsColorImage(slice) != iftIsColorImage(volume)) {
            iftError(
                    "Both the current slice (%s) image being loaded and the first one must either be colored or grayscaled!"
                    "iftReadImageFolderAsVolume", path);
		}

		iftPutXYSlice(volume, slice, z);

		iftDestroyImage(&slice);
	}

	return volume;
}


void iftWriteVolumeAsVideoFolder(iftImage *video, const char *folder_name,
								 char *frame_name)
{
	int i;
	char *pos, ext[10];
	iftImage *frame = NULL;

	pos = strrchr(frame_name,'.');
	pos++;
	sscanf(pos,"%s",ext);

	if((strcmp(ext,"ppm") == 0 && !iftIsColorImage(video))
	   || (strcmp(ext,"pgm") == 0 && iftIsColorImage(video)))
	{
        iftError("Invalid frame file format", "iftWriteVolumeAsVideoFolder");
	}

	for(i = 0; i < video->zsize; i++)
	{
		frame = iftGetXYSlice(video, i);
		iftWriteFrame(frame, folder_name, frame_name, i+1);

		iftDestroyImage(&frame);
	}
}


void iftConvertVideoFileToVideoFolder(const char *video_path, const char *output_folder,
									  int rotate)
{
	iftConvertVideoFramesToImages(video_path, "_tmp", "frame", "ppm", rotate);

	iftStoreFramesInVideoFolder("_tmp","ppm", output_folder, "frame");

	if(system("rm -rf _tmp") == -1) iftError("Command error", "iftConvertVideoFileToVideoFolder");
}

void iftConvertVideoFramesToImages(const char *video_path, const char *output_folder,
								   const char *frame_name, const char *extension,
								   int rotate)
{
	char cmd[65536];
	char output_fname[65536];

	sprintf(output_fname,"%s%%0%dd.%s",frame_name, IFT_VIDEO_FOLDER_FRAME_NZEROES, extension);

	sprintf(cmd,"mkdir -p %s", output_folder);
	if(system(cmd) == -1) iftError("Command error", "iftConvertVideoFramesToImages");

	if(rotate < 0)
		sprintf(cmd, "ffmpeg -i %s -f image2 %s", video_path, output_fname);
	else
		sprintf(cmd, "ffmpeg -i %s -f image2 -vf \"transpose=%d\" %s", video_path,
				rotate, output_fname);

	if(system(cmd) == -1) iftError("Command error", "iftConvertVideoFramesToImages");

	sprintf(cmd, "mv %s*.%s %s",frame_name, extension, output_folder);
	if(system(cmd) == -1) iftError("Command error", "iftConvertVideoFramesToImages");

}

void iftStoreFramesInVideoFolder(const char *input_folder, const char *input_extension,
								 const char *output_folder, const char *output_filename)
{
	char cmd[65536];
	char expr[] = "j=1;" \
				"for i in $(ls %s/*.%s); do " \
					"number=`printf \"%%0%dd\" $j`;"
			"mkdir -p %s/$number; " \
					"mv $i %s/$number/%s.%s;" \
					"j=`expr $j + 1`; " \
				"done";
	sprintf(cmd, expr, input_folder, input_extension, IFT_VIDEO_FOLDER_FRAME_NZEROES,
			output_folder, output_folder, output_filename,	input_extension);

	if(system(cmd) == -1) iftError("Command error", "iftStoreFramesInVideoFolder");
}

