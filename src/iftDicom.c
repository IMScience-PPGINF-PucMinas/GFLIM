#include "iftDicom.h"

#include "ift/core/io/Dir.h"
#include "ift/core/io/Stream.h"


/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Counts the number of Dicom Images (slices) inside the directory <dir>.
 *
 * Counts the number of Dicom Images (slices) inside the directory <dir>.
 *
 * @param dir The directory to be used.
 * @return The number of Dicom Images (slices).
 */
int iftCountDicomSlices(iftDir *dir) {
    FILE *fp = NULL;
    char cmd[512];

    int nslices = 0;

    char *temp_file = iftMakeTempFile(NULL, NULL, NULL);
    
    for (int i = 0; i < dir->nfiles; i++) {
        if (iftFileExists(dir->files[i]->path)) { // if file exists
            // Get the informations from gdcminfo, storing them into a temp file
            sprintf(cmd, "gdcminfo -i %s > %s", dir->files[i]->path, temp_file);
            if (system(cmd) == -1)
                iftError("System command error: Check if \"libgdcm-tools\" is installed!",
                         "iftCountDicomSlices");

            fp = fopen(temp_file, "r");

            if (!iftIsFileContentEmpty(fp)) // if the file is a dicom slice (2D image)
                nslices++;

            fclose(fp);
        }
    }
    unlink(temp_file); // remove the temp_file from disk
    iftFree(temp_file);

    return nslices;
}


/**
 * @brief Reads the dicom slice informations from a Text File.
 *
 * Reads the dicom slice informations from a Text File.
 * All informations are "returned" into the references of the variables passed as parameters.
 *
 * @param fp File pointer from the Dicom informations. It must be opened.
 * @param xsize Reference to the variable that will receive the xsize value.
 * @param ysize Reference to the variable that will receive the ysize value.
 * @param bits_per_voxel Reference to the variable that will receive the bits_per_voxel value.
 * @param z_position Reference to the variable that will receive the z_position value.
 * @param dx Reference to the variable that will receive the dx value.
 * @param dy Reference to the variable that will receive the dy value.
 * @param orientation Reference to the variable that will receive the orientation value.
 */
void iftReadDicomSliceInfo(FILE *fp, int *xsize, int *ysize, int *bits_per_voxel, float *z_position, float *dx, float *dy, char **orientation, float orientation_vectors[2][3]) {

    char line[512], aux[512];
    char *pos = NULL, *orientation_aux = NULL;
    int trash1;
    float trash2, trash3;

    /* Read DICOM data */
    while (fgets(line, 512, fp) != NULL) { // read line by line with at most 512 characters per line
        if ((strstr(line, "Dimensions") != NULL) && (line[0] == 'D')) {
            pos = strstr(line, "("); // shift the '('
            sscanf(pos, "(%d,%d,%d)", xsize, ysize, &trash1);
        }
        else if (strstr(line, "BitsAllocated") != NULL) {
            pos = strstr(line, ":");
            sscanf((pos+1), "%d", bits_per_voxel);
        }
        else if (strstr(line, "Origin:") != NULL) {
            pos = strstr(line, "(");
            sscanf(pos, "(%f,%f,%f)", &trash2, &trash3, z_position);
        }
        else if (strstr(line, "Spacing") != NULL) {
            pos = strstr(line, "(");
            sscanf(pos, "(%f,%f,%f)", dx, dy, &trash2);
        }
        else if (strstr(line, "Orientation Label") != NULL) {
            pos = strstr(line, ":");
            sscanf((pos+2), "%s", aux);
            orientation_aux              = iftAllocCharArray(strlen(aux)+1);
            strcpy(orientation_aux, aux);
            orientation_aux[strlen(aux)] = '\0';
            *orientation                 = orientation_aux;
        }
	//Gets the orientation vectors as made available by gdcminfo.
	else if (strstr(line, "DirectionCosines") != NULL) {
            pos = strstr(line, "(");
            sscanf(pos, "(%f,%f,%f,%f,%f,%f)", &orientation_vectors[0][0], &orientation_vectors[0][1], &orientation_vectors[0][2],
						   &orientation_vectors[1][0], &orientation_vectors[1][1], &orientation_vectors[1][2]);
	        /*printf("Vector was read: (%f,%f,%f,%f,%f,%f)", orientation_vectors[0][0], orientation_vectors[0][1], orientation_vectors[0][2],
						   orientation_vectors[1][0], orientation_vectors[1][1], orientation_vectors[1][2]);*/
        }
    }
}


/**
 * @brief Function that is used to compare 2 dicom slices in qsort.
 */
int iftCompareDicomSlices(const void *a, const void *b) {
    iftDicomSlice **slice1 = (iftDicomSlice**) a;    
    iftDicomSlice **slice2 = (iftDicomSlice**) b;

    if ((*slice1)->z_position < (*slice2)->z_position)
        return -1;
    else if ((*slice1)->z_position > (*slice2)->z_position)
        return 1;
    else
        return 0;
}


/**
 * @brief Sorts the dicom slices from a Dicom image.
 *
 * Sorts the dicom slices from a Dicom image in ascending order, according to their z_position.
 * qsort is used to sort the vector.
 * Obs: dicom->zsize == nslices
 *
 * @param dicom Dicom image.
 */
void iftSortDicomSlices(iftDicomImage *dicom) {
    if (dicom != NULL) {
        if (dicom->zsize > 0) {
            qsort(dicom->slices, dicom->zsize, sizeof(iftDicomSlice*), iftCompareDicomSlices);
        }
        else {
            iftWarning("Number of dicom slices (zsize) is <= 0. Nothing to do!", "iftSortDicomSlices");
            puts("");
        }
    }
}
/******************************************************************/



/********************** PUBLIC FUNCTIONS *************************/
/**
 * @brief Creates an iftDicomSlice.
 *
 * Allocates an iftDicomSlice on memory.
 *
 * @return An iftDicomSlice.
 */
iftDicomSlice *iftCreateDicomSlice() {
    iftDicomSlice *dicom_slice = (iftDicomSlice*) iftAlloc(1, sizeof(iftDicomSlice));
    return dicom_slice;
}


/**
 * @brief Destroys an iftDicomSlice.
 *
 * Destroys the iftDicomSlice <dicom_slice> from the memory. 
 *
 * @param dicom_slice The reference from the iftDicomSlice.
 * @param dicom_slice The reference from the iftDicomSlice.
 */
void iftDestroyDicomSlice(iftDicomSlice **dicom_slice) {
    iftDicomSlice *ds_aux = *dicom_slice;

    if (ds_aux != NULL) {
        iftFree(ds_aux);
        *dicom_slice = NULL;
    }
}



/**
 * @brief Creates an iftDicomImage.
 *
 * Allocates an iftDicomImage on memory.
 *
 * @param nslices The number of dicom slices xy from the dicom image. zsize = nslices.
 * @return An iftDicomImage.
 */
iftDicomImage *iftCreateDicomImage(int nslices) {
    if (nslices < 0)
        iftError("Number of dicom slices < 0", "iftCreateDicomImage");

    iftDicomImage *dicom  = (iftDicomImage*) iftAlloc(1, sizeof(iftDicomImage));
    dicom->xsize          = 0;
    dicom->ysize          = 0;
    dicom->zsize          = nslices;
    dicom->orientation    = NULL;
    dicom->dx             = 0;
    dicom->dy             = 0;
    dicom->dz             = 0;
    dicom->bits_per_voxel = 0;

    dicom->slices = (iftDicomSlice**) iftAlloc(nslices, sizeof(iftDicomSlice*));
    for (int z = 0; z < nslices; z++)
        dicom->slices[z] = iftCreateDicomSlice();

    return dicom;
}

/**
 * @brief Determines if two sets of orientation vectors are the same (have the same elements).
 *
 * @param vecs1, vecs2 The two vector sets to be compared.
 *
 * @return 0 if different, 1 if equal.
 */
int iftCompareOrientationVectors(float vec1[2][3], float vec2[2][3]) {
	int i, j;

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 3; j++) {
			if (vec1[i][j] != vec2[i][j]) return 0; 
		}
	}

	return 1;
}

/**
 * @brief Reads all dicom images (slices) from a directory.
 *
 * Reads and sorts all dicom images (slices) from a directory and returns an iftDicomImage with all of them.
 * The dicom slices are sorted by their z_position values in ascending order.
 *
 * @param dir_pathname The directory pathname to be read.
 * @return The iftDicomImage with the informations from the dicom image.
 */
iftDicomImage *iftReadDicomImage(char *dir_pathname) {
    FILE *fp = NULL;
    char msg[512], cmd[512];
    int xsize, ysize, bits_per_voxel;
    float dx, dy;
    char *orientation;
    float orientation_vectors[2][3];

    iftDicomImage *dicom = NULL;

    // load all files (not subdirs) from the 1st hierarchical level from the directory <dir_pathname>
    // this function already checks if the directory exists
    iftDir *dicom_dir = iftLoadFilesFromDirBySuffix(dir_pathname, "");

    int nslices = iftCountDicomSlices(dicom_dir);
    
    if (nslices == 0) {
        sprintf(msg, "There is no dicom image in the directory \"%s\". It will return NULL", dicom_dir->path);
        iftWarning(msg, "iftReadDicomImage");
        puts("");
        iftDestroyDir(&dicom_dir);
        return NULL;
    }
    else {
        dicom = iftCreateDicomImage(nslices);
        char *temp_file = iftMakeTempFile(NULL, NULL, NULL);

        /* reads the first dicom slice out of the main loop in order to compare if all dicom slices
           have the same properties, such as xsize, ysize, ...
        */
        int first_slice_idx = -1;
        // for each file inside the directory
        for (int i = 0; i < dicom_dir->nfiles; i++) {

            // Get the informations from gdcminfo, storing them into a temp file
            sprintf(cmd, "gdcminfo -i %s > %s", dicom_dir->files[i]->path, temp_file);
            if (system(cmd) == -1)
                iftError("System command error: Check if \"libgdcm-tools\" is installed!",
                         "iftReadDicomImage");

            fp = fopen(temp_file, "r");

            if (!iftIsFileContentEmpty(fp)) { // it is a dicom slice
                first_slice_idx = i;
                
                strcpy(dicom->slices[0]->filename, dicom_dir->files[i]->path);

                iftReadDicomSliceInfo(fp, &dicom->xsize, &dicom->ysize, &dicom->bits_per_voxel,
                    &dicom->slices[0]->z_position, &dicom->dx, &dicom->dy, &dicom->orientation, dicom->orientation_vectors);

                /*
                 * An error will raise if dicom->xsize or dicom->ysize is 0, every source code using dicom images has
                 * to check if they are 0 and disregard that slice
                 */

                fclose(fp);
                break;
            }
            fclose(fp);
        }

        int z = 1;
        /* finds the remaining dicom slices and compare their informations with the first */
        for (int i = first_slice_idx+1; (i < dicom_dir->nfiles) && (z < nslices); i++) {
            // Get the informations from gdcminfo, storing them into a temp file
            sprintf(cmd, "gdcminfo -i %s > %s", dicom_dir->files[i]->path, temp_file);
            if (system(cmd) == -1)
                iftError("System command error: Check if \"libgdcm-tools\" is installed!",
                         "iftReadDicomImage");
            
            fp = fopen(temp_file, "r");

            if (!iftIsFileContentEmpty(fp)) { // it is a dicom slice
                strcpy(dicom->slices[z]->filename, dicom_dir->files[i]->path);

                iftReadDicomSliceInfo(fp, &xsize, &ysize, &bits_per_voxel,
                    &dicom->slices[z]->z_position, &dx, &dy, &orientation, orientation_vectors);

                // Checking the slice information with the first dicom slice

                if (dicom->xsize != xsize) {
                    printf("\nError in Dicom file: %s \nSlice X Size: %d\nX Size from reference slice: %d\n",dicom_dir->files[i]->path, dicom->xsize, xsize);
                    iftError("Multiple values of xsize in the dicom images", "iftReadDicomImage");
                }
                if (dicom->ysize != ysize) {
                    printf("\nError in Dicom file: %s \nSlice Y Size: %d\nY Size from reference slice: %d\n",dicom_dir->files[i]->path, dicom->ysize, ysize);
                    iftError("Multiple values of ysize in the dicom images", "iftReadDicomImage");
                }
                if (dicom->bits_per_voxel != bits_per_voxel) {
                    printf("\nError in Dicom file: %s \nBits per voxel in Slice: %d\nBits per voxel in reference slice: %d\n",dicom_dir->files[i]->path, dicom->bits_per_voxel, bits_per_voxel);
                    iftError("Multiple values of bits_per_voxel in the dicom images", "iftReadDicomImage");
                }
                if (fabs(dicom->dx-dx) > IFT_EPSILON) {
                    printf("\nError in Dicom file: %s \nSlice Dx: %f\nDx in reference slice: %f\n",dicom_dir->files[i]->path, dicom->dx, dx);
                    iftError("Multiple values of dx in the dicom images", "iftReadDicomImage");
                }
                if (fabs(dicom->dy-dy) > IFT_EPSILON) {
                    printf("\nError in Dicom file: %s \nSlice Dy: %f\nDy in reference slice: %f\n",dicom_dir->files[i]->path, dicom->dy, dy);
                    iftError("Multiple values of dy in the dicom images", "iftReadDicomImage");
                }
                if (strcmp(dicom->orientation, orientation) != 0) {
                    printf("\nError in Dicom file: %s \nSlice Orientation: %s\nOrientation in reference slice: %s\n",dicom_dir->files[i]->path, dicom->orientation, orientation);
                    iftError("Multiple values of orientation in the dicom images", "iftReadDicomImage");
                }
                if (iftCompareOrientationVectors(dicom->orientation_vectors, orientation_vectors) == 0) {
                    printf("\nError in Dicom file: %s \n",dicom_dir->files[i]->path);
                    iftError("Multiple values of orientation vectors in the dicom images", "iftReadDicomImage");
                }

                
                iftFree(orientation);
                z++;
            }
            fclose(fp);
        }

        iftDestroyDir(&dicom_dir);
        unlink(temp_file); // remove the temp_file from the disk
        iftFree(temp_file);

        int nvalid_slices = 0;
        if (nslices > 1) {
            // Sorts the dicom slices according to their z_positions in increasing order
            iftSortDicomSlices(dicom);
            dicom->dz = 0.0;
            int i;
            float dz[nslices];
            for (i = 0; ((i < (dicom->zsize-1)) && (iftAlmostZero(dicom->dz))); i++){
                dicom->dz = fabs(dicom->slices[i+1]->z_position - dicom->slices[i]->z_position);
            }
            printf("===========\ndicom->dz = %f\n===========\n",dicom->dz);
            // checks if the "dz" between each subsequent pair of slices are the same
            // remember: dicom->zsize == nslices
            for (int z = 0; z < (dicom->zsize-1); z++) {
                dz[z] = fabs(dicom->slices[z+1]->z_position - dicom->slices[z]->z_position);
                //printf("%s - %s = dz = %f\n",dicom->slices[z+1]->filename,dicom->slices[z]->filename,dz);
                if (dz[z] > IFT_EPSILON)
                    nvalid_slices++;
            }

            if (nvalid_slices != nslices) {

                iftDicomSlice **new_slices = NULL;
                new_slices = (iftDicomSlice **) iftAlloc(nvalid_slices, sizeof(iftDicomSlice *));

                for (i = 0; i < nvalid_slices; i++)
                    new_slices[i] = iftCreateDicomSlice();
                int pos = 0;
                for (int z = 0; z < (dicom->zsize - 1); z++) {
                    if (dz[z] > IFT_EPSILON) {
                        strncpy(new_slices[pos]->filename,dicom->slices[z]->filename,512);
                        new_slices[pos++]->z_position = dicom->slices[z]->z_position;
                    }
                }
                for (int z = 0; z < dicom->zsize; z++)
                    iftDestroyDicomSlice(&dicom->slices[z]);
                dicom->slices = new_slices;
                dicom->zsize = nvalid_slices;
            }
        }
    }

    return dicom;
}

/*
 * @brief Destroys an iftDicomImage.
 *
 * Destroys the iftDicomImage <dicom> from the memory. 
 *
 * @param dicom The reference from the iftDicomImage.
 */
void iftDestroyDicomImage(iftDicomImage **dicom) {
    iftDicomImage *dicom_aux = *dicom;

    if (dicom_aux != NULL) {
        iftFree(dicom_aux->orientation);

        for (int z = 0; z < dicom_aux->zsize; z++)
            iftDestroyDicomSlice(&dicom_aux->slices[z]);
        iftFree(dicom_aux);
        *dicom = NULL;
    }
}

/**
 * @brief Prints the informations from an iftDicomImage.
 */
void iftPrintDicomInfo(iftDicomImage *dicom) {
    if (dicom != NULL) {
        puts("------------------------------------------");
        printf("- xsize = %d, ysize = %d, zsize = %d\n", dicom->xsize, dicom->ysize, dicom->zsize);
        printf("- dx = %f, dy = %f, dz = %f\n", dicom->dx, dicom->dy, dicom->dz);
        printf("- orientation = %s\n", dicom->orientation);
	    printf("- direction vectors = (%f,%f,%f), (%f,%f,%f)\n", 
					dicom->orientation_vectors[0][0],dicom->orientation_vectors[0][1],dicom->orientation_vectors[0][2],
					dicom->orientation_vectors[1][0],dicom->orientation_vectors[1][1],dicom->orientation_vectors[1][2]);
        for (int z = 0; z < dicom->zsize; z++) {
            printf("--> %s - z_position = %f\n", dicom->slices[z]->filename, dicom->slices[z]->z_position);
        }
        puts("------------------------------------------\n");
    }
}

/*
 * @brief Converts an ift dicom image to the scene format (ift image)
 *
 * @param dicom The dicom image to be converted.
 * @return The converted scene image.
 */
iftImage *iftConvertDicom2Scene(iftDicomImage * dicom) {
  iftImage *img; 
  uchar *data8  = NULL; /* 8-bit data are always from [0-255] */ 
  short *data16 = NULL; /* 16-bit data might be negative */
  int delta_z, n_read_pixels;
  char cmd[2048], msg[2048];
  FILE *fp2;
  
  int  xysize = dicom->xsize * dicom->ysize;

  // Initial setup for the conversion.
  if (dicom->bits_per_voxel == 8)
    data8  = iftAllocUCharArray(xysize);
  else if (dicom->bits_per_voxel == 16)
    data16 = iftAllocShortArray(xysize);
  else {
    sprintf(msg, "Number of bits must be 8 or 16: found %d", dicom->bits_per_voxel);
      iftError(msg, "iftConvertDicom2Scene");
  }
  
  img = iftCreateImage(dicom->xsize, dicom->ysize, dicom->zsize);
  img->dx = dicom->dx;
  img->dy = dicom->dy;
  img->dz = dicom->dz;
  
    char *temp_file1 = iftMakeTempFile(NULL, NULL, NULL);
    char *temp_file2 = iftMakeTempFile(NULL, NULL, NULL);
    
    // Conversion process per se: goes through different conversion phases,
    // considers the appropriate type of integer.
    // #pragma omp parallel for shared(dicom) private(fp2)
    for (int z = 0; z < dicom->zsize; z++) {        
      printf("Converting %s\n",dicom->slices[z]->filename);
      sprintf(cmd, "gdcmconv --raw %s %s", dicom->slices[z]->filename, temp_file1);
      
      if (system(cmd) == -1)
          iftError("System command error: gdcmconv", "iftConvertDicom2Scene");
      
      sprintf(cmd, "gdcmraw  -i %s -o %s", temp_file1, temp_file2);
      if (system(cmd) == -1)
          iftError("System command error: gdcmraw", "iftConvertDicom2Scene");
      
      fp2     = fopen(temp_file2, "rb");
      delta_z = z*xysize;
      if (dicom->bits_per_voxel == 8) {
          n_read_pixels = fread(data8,sizeof(uchar),xysize,fp2);
          if (n_read_pixels != xysize) {
              sprintf(msg, "Reading error (dicom image of 8 bits): Number of Pixels from slice %d is different from xysize: %d, %d\n" \
		  "Slice: %s", z, n_read_pixels, xysize, dicom->slices[z]->filename);
              iftError(msg, "iftConvertDicom2Scene");
          }
          for (int p = 0; p < xysize; p++)
              img->val[p+delta_z] = data8[p];
      } else {
          n_read_pixels = fread(data16,sizeof(short),xysize,fp2);
          if (n_read_pixels != xysize) {
              sprintf(msg, "Reading error (dicom image of 16 bits): Number of Pixels from slice %d is different from xysize: %d, %d\n" \
		  "Slice: %s", z, n_read_pixels, xysize, dicom->slices[z]->filename);
              iftError(msg, "iftConvertDicom2Scene");
          }
          for (int p = 0; p < xysize; p++)
              img->val[p+delta_z] = data16[p];
      }
      fclose(fp2);
    }
    
    unlink(temp_file1);
    unlink(temp_file2);

    iftFree(temp_file1);
    iftFree(temp_file2);

    if (dicom->bits_per_voxel==8)
        iftFree(data8);
    else
        iftFree(data16);

    /* 16-bit images should be normalized within [0-4095] */
    
    if (dicom->bits_per_voxel==16){
      iftImage *nimg;
      int img_min_val = iftMinimumValue(img);
      for (int p=0; p < img->n; p++) /* make sure all values are
					non-negatives */
	img->val[p] = img->val[p]-img_min_val;

      nimg = iftNormalizeWithNoOutliers(img,0,4095,0.98);
      iftDestroyImage(&img);    
      return(nimg);
    } else {
      return(img);
    }

}

/*
 * @brief Performs the cross product between two 3D vectors.
 *
 * Uses the formula V X W = (Vy.Wz - Vz Wy) i + (Vz Wx - Vx Wz) j + 
 * (Vx Wy - Vy Wx) k, where i, j and k are the unity vectors.
 * Source: http://math.oregonstate.edu/bridge/papers/dot+cross.pdf
 *
 * @param vec1, vec2 The vectors between which the function will get the 
 * cross product.
 *
 * @return The vector resulting from the cross product.
 */
iftMatrix * ift3DVectorCrossProduct(iftMatrix * vec1, iftMatrix * vec2) {

    iftMatrix * result = iftCreateMatrix(1,3);

    result->val[0] = (vec1->val[1] * vec2->val[2]) - 
                     (vec1->val[2] * vec2->val[1]);
    result->val[1] = (vec1->val[2] * vec2->val[0]) - 
                     (vec1->val[0] * vec2->val[2]);
    result->val[2] = (vec1->val[0] * vec2->val[1]) - 
                     (vec1->val[1] * vec2->val[0]);

    return result;
}

/*
 * @brief Adds the values in two matrices and returns the resultant matrix.
 *
 * @param mat1 mat2 The matrices to be added.
 *
 * @return The resultant matrix.
 */
iftMatrix * iftAddMatrices (iftMatrix * mat1, iftMatrix * mat2) {
    int i;
    iftMatrix * sum;

    if (mat1->ncols != mat2->ncols || mat1->nrows != mat2->nrows) {
        iftWarning("Matrices cannot be added because they have different dimensions.",
                   "iftAddMatrices");
    }

    sum = iftCreateMatrix(mat1->ncols, mat1->nrows);

    for (i = 0; i < mat1->n; i++) {
        sum->val[i] = mat1->val[i] + mat2->val[i];
    }

    return sum;
}

/*
 * @brief Converts an ift image ("img") to the reference coordinate system using 
 * the orientation vectors of the dicom image ("dicom") which originated "img".
 *
 * The reference coordinate system is used so that all images can be visualized 
 * the same way.
 *
 * @param img The image to be converted.
 * @param dicom The dicom image equivalent to img.
 *
 * @return The image converted to the new coordinate system.
 */
iftImage * iftConvertDicomCoordinates(iftImage * img, iftDicomImage * dicom) {
    int i, j, p, q;
    iftMatrix * coordVec, * new_coordVec; //Coordinate vectors.
    iftVoxel v, u;
    iftImage * new_img = NULL;
    int new_xsize = 0, new_ysize = 0, new_zsize = 0;

    /* Approximate the values in the orientation vectors to either -1.0, 0.0,or 1.0 */ 

    for (i = 0; i < 2; i++) {
      for (j = 0; j < 3; j++) {
	if (dicom->orientation_vectors[i][j]<0.0){
	  if (fabs(dicom->orientation_vectors[i][j]+1.0)<fabs(dicom->orientation_vectors[i][j])){
	    dicom->orientation_vectors[i][j]=-1.0;
	  }else{
	    dicom->orientation_vectors[i][j]=+0.0;
	  }
	} else {
	  if (fabs(dicom->orientation_vectors[i][j]-1.0)<fabs(dicom->orientation_vectors[i][j])){
	    dicom->orientation_vectors[i][j]=+1.0;
	  }else{
	    dicom->orientation_vectors[i][j]=+0.0;
	  }
	}
      }
    }

    /* Compute transformation matrix */

    iftMatrix * transformation_matrix = iftCreateMatrix(3, 3);

    for (i = 0; i < 2; i++) {
      for (j = 0; j < 3; j++) {
	transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, j, i)] = dicom->orientation_vectors[i][j];
      }
    }

    if (strcmp(dicom->orientation,"AXIAL")==0){
      transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, 0, 2)] = 0;
      transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, 1, 2)] = 0;
      transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, 2, 2)] = 1;
    } else {
      if (strcmp(dicom->orientation,"CORONAL")==0){
    	transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, 0, 2)] = 0;
    	transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, 1, 2)] = 1;
    	transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, 2, 2)] = 0;
      } else { /* SAGITAL */
    	transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, 0, 2)] = 1;
    	transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, 1, 2)] = 0;
    	transformation_matrix->val[iftGetMatrixIndex(transformation_matrix, 2, 2)] = 0;
      }
    }

    iftPrintMatrix(transformation_matrix);
    printf("%s\n",dicom->orientation);

    //Allocates memory for the image in the new coordinates.

    //Find which axes correspond to x, y and z in the reference coordinate 
    //system, and assign the sizes to be allocated accordingly.

    if (fabs(transformation_matrix->val[iftGetMatrixIndex
        (transformation_matrix, 0, 0)]) == 1) {
        new_xsize = img->xsize;
    } else if (fabs(transformation_matrix->val[iftGetMatrixIndex
        (transformation_matrix, 0, 1)]) == 1) {
        new_xsize = img->ysize;
    } else if (fabs(transformation_matrix->val[iftGetMatrixIndex
        (transformation_matrix, 0, 2)]) == 1) {
        new_xsize = img->zsize;
    }

    if (fabs(transformation_matrix->val[iftGetMatrixIndex
        (transformation_matrix, 1, 0)]) == 1) {
        new_ysize = img->xsize;
    } else if (fabs(transformation_matrix->val[iftGetMatrixIndex
        (transformation_matrix, 1, 1)]) == 1) {
        new_ysize = img->ysize;
    } else if (fabs(transformation_matrix->val[iftGetMatrixIndex
        (transformation_matrix, 1, 2)]) == 1) {
        new_ysize = img->zsize;
    }

    if (fabs(transformation_matrix->val[iftGetMatrixIndex
        (transformation_matrix, 2, 0)]) == 1) {
        new_zsize = img->xsize;
    } else if (fabs(transformation_matrix->val[iftGetMatrixIndex
        (transformation_matrix, 2, 1)]) == 1) {
        new_zsize = img->ysize;
    } else if (fabs(transformation_matrix->val[iftGetMatrixIndex
        (transformation_matrix, 2, 2)]) == 1) {
        new_zsize = img->zsize;
    }

    new_img = iftCreateImage(new_xsize, new_ysize, new_zsize);
    iftCopyVoxelSize(img,new_img);

    /*Transforms the original dicom image by multiplying each pixel coordinate 
     *by the transformation matrix that was found and assigning that pixel's 
     *value to the resulting coordinate.*/

    coordVec = iftCreateMatrix(3,1);
    for (u.z = 0; u.z < img->zsize; u.z++) {
        for (u.y = 0; u.y < img->ysize; u.y++) {
            for (u.x = 0; u.x < img->xsize; u.x++) {
                //The vector which indicates the current coordinates.
                coordVec->val[iftGetMatrixIndex(coordVec, 0, 0)] = u.x;
                coordVec->val[iftGetMatrixIndex(coordVec, 1, 0)] = u.y;
                coordVec->val[iftGetMatrixIndex(coordVec, 2, 0)] = u.z; 
                //Vector after multiplication and addition: new coordinate system.
                new_coordVec = iftMultMatrices(coordVec, transformation_matrix);
                v.x = new_coordVec->val[iftGetMatrixIndex(new_coordVec, 0, 0)];
                v.y = new_coordVec->val[iftGetMatrixIndex(new_coordVec, 1, 0)];
                v.z = new_coordVec->val[iftGetMatrixIndex(new_coordVec, 2, 0)];
		if (v.x < 0) v.x = new_xsize - 1 + v.x;
		if (v.y < 0) v.y = new_ysize - 1 + v.y;
		if (v.z < 0) v.z = new_zsize - 1 + v.z;
                //Assigns the value of the current voxel to its 
                //correspondent coordinate.
                p = iftGetVoxelIndex(img,u);
                q = iftGetVoxelIndex(new_img,v);
                new_img->val[q] = img->val[p];
		iftDestroyMatrix(&new_coordVec);
            }
        }
    }

    iftDestroyMatrix(&transformation_matrix);		
    iftDestroyMatrix(&coordVec);		

    return new_img;
}


/*
 * @brief Build the Dicom Transformation Matrix.
 *
 * The coordinates from the 3rd vector always will have positive, because the dicom slices
 * are read correctly. 
 *
 * @param v1 Unit Vector.
 * @param v2 Unit Vector.
 */
iftMatrix *iftDicomTransMatrix(iftVector v1, iftVector v2) {
    // iftVector v3;

    // // Cross Product: returns an orthogonal unit vector v3 with respect to plan v1v2
    // // We force the unit vector v3 to be positive, because the dicom slices are read of the correct way
    // v3.x = fabs((v1.y * v2.z) - (v1.z * v2.y));
    // v3.y = fabs((v1.z * v2.x) - (v1.x * v2.z));
    // v3.z = fabs((v1.x * v2.y) - (v1.y * v2.x));

    // iftMatrix *T = iftCreateMatrix(3,3);

    return NULL;
}


/*
 * @brief Validate the coordinates from the Orientation Vectors of the Dicom Image. 
 *
 * The coordinate must have values -1 or 0 or 1.
 *
 * @param dicom The dicom image.
 */
void iftValidateDicomOrientationVectors(iftDicomImage *dicom) {
    float coord;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++) {
            coord = dicom->orientation_vectors[0][0];
            if ((coord != -1) && (coord != 0) && (coord != 1)) {
                char msg[512];
                sprintf(msg, "Invalid Orientation Vector Coord: [%d][%d] = %f\n" \
                             "Value must be: -1 or 0 or 1", i, j, coord);
                iftError(msg, "iftValidateDicomOrientationVectors");
            }
        }
}







