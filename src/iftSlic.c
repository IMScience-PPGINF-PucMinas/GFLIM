#include "iftSlic.h"

#include "ift/core/io/Stream.h"



/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Create the Slic Parameter Structure
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 */
iftSlic *iftCreateSlic(double compactness, int step, int img_max_range, bool is_3D_img, bool is_color_img) {
    iftSlic *slic = (iftSlic*) iftAlloc(1, sizeof(iftSlic));

    slic->compactness   = compactness;
    slic->step          = step;
    slic->img_max_range = img_max_range;
    slic->is_3D_img     = is_3D_img;
    slic->is_color_img  = is_color_img;

    return slic;
}


/**
 * @brief Allocates the Feature Vectors for each cell's centers
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 */
void iftAllocFeatVectorsForClusterCenters(iftSlic *slic) {
    slic->color_feats   = (iftFColor*) iftAlloc(slic->n_cells, sizeof(iftFColor));
    slic->spatial_feats = (iftPoint*) iftAlloc(slic->n_cells, sizeof(iftPoint));
}


/**
 * @brief Gets the Color LAB (or Gray) Feature Vector from the Voxel <b>voxel_idx</b> that is stored
 * in the Multi-Band Image <b>mimg</b>. 
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 *
 * @warning If the Image is Colored, the multi-band image has 3 bands. Otherwise, it has only 1 band.
 * 
 * @param mimg          Multi-Band Image that contains the converted (3 bands) LAB Image
 *                      (if the YCbCr Image is colored) or the original Grayscale Image (1 band) 
 * @param voxel_idx     Index from the requested voxel
 * @param img_max_range Maximum range of the image to be oversegmented: e.g 255 (for 8 bits images)
 * @return              The Color (or Gray) Feature Vector.
 */
iftFColor iftGetColorFeats(const iftMImage *mimg, int voxel_idx, int img_max_range) {
    iftFColor cfeats;

    // LAB Color Image
    if (mimg->m == 3) {
        cfeats.val[0] = mimg->val[voxel_idx][0]; // L
        cfeats.val[1] = mimg->val[voxel_idx][1]; // A
        cfeats.val[2] = mimg->val[voxel_idx][2]; // B
    }
    // GRAYSCALE IMAGE
    else {
        cfeats.val[0] = mimg->val[voxel_idx][0];
    }

    return cfeats;
}



/**
 * @brief Finds the voxels with the Minimum Gradient around each cell's centers, within a region of
 * 3x3 (2D image) or 3x3x3 (3D image), and moves the centers to there.
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 *
 * @note The cell's center coordinates are in the spatial feature vector in slic. 
 * @warning It is expected that the input image has at least 3x3 (2D) or 3x3x3 (3D)
 * 
 * @param slic Slic Parameters that contains the cell's centers and feature vectors.
 * @param mimg Multi-Band Image that contains the converted (3 bands) LAB Image
 *             (if the YCbCr Image is colored) or the original Grayscale Image (1 band) 
 * @param grad Image Gradient used.
 */
void iftFindLocalMinimumForSlicSeeds(iftSlic *slic, const iftMImage *mimg, const iftImage *grad) {
    iftAdjRel *A = NULL;
    
    // Adjacency of 3x3 (2D) or 3x3x3
    if (slic->is_3D_img)
        A = iftSpheric(1.74); // sqrt(3.0)
    else A = iftCircular(1.42); // sqrt(2.0)

    for (int cell = 0; cell < slic->n_cells; cell++) {
        iftVoxel center = iftPointToVoxel(slic->spatial_feats[cell]);
        if (!slic->is_3D_img)
            center.z = 0;

        int c        = iftGetVoxelIndex(grad, center);
        int min_grad = grad->val[c]; // gradient from the cluster center voxel

        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, center, i);

            if (iftValidVoxel(grad, v)) {
                int q = iftGetVoxelIndex(grad, v);
                
                if (grad->val[q] < min_grad) {
                    min_grad = grad->val[q];
                    c = q;
                }       
            }
        }

        // assigns the final center voxel from the cell <cell>
        center = iftGetVoxelCoord(grad, c);
        slic->spatial_feats[cell] = iftVoxelToPoint(center);

        // assigns the color space features from the final center voxel of the cluster
        slic->color_feats[cell] = iftGetColorFeats(mimg, c, slic->img_max_range);
    }

    // DESTROYERS
    iftDestroyAdjRel(&A);
}

/**
 * @brief Performs Kmeans segmentation to find the superpixels.
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 *
 * @param  slic Slic Parameter with cell's center feature vectors
 * @param  mimg Multi-Band Image that contains the converted (3 bands) LAB Image
 *              (if the YCbCr Image is colored) or the original Grayscale Image (1 band)
 * @return      The Superpixel Image.
 */
iftImage *iftPerformSuperpixelsSlic(iftSlic *slic, const iftMImage *mimg) {
    iftImage *superpixels      = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    double spatial_dist_factor = (slic->compactness*slic->compactness) / (slic->step*slic->step);

    // minimum distance from each pixel to best neighbor cell center
    double *min_dists = iftAllocDoubleArray(mimg->n);

    // runs Slic for 10 iteration at most
    for (int it = 1; it <= 10; it++) { //10

      for (int p = 0; p < mimg->n; p++)
	min_dists[p] = IFT_INFINITY_DBL;

      for (int cell = 0; cell < slic->n_cells; cell++) {
	iftPoint center         = slic->spatial_feats[cell];
	iftFColor cfeats_center = slic->color_feats[cell];
	
	// finds the region of size (step x step) or (step x step x step) around the cell's center
	iftBoundingBox bb;
	bb.begin.x = iftMax(0, center.x - slic->step);
	bb.end.x   = iftMin(mimg->xsize-1, center.x + slic->step);
	
	bb.begin.y = iftMax(0, center.y - slic->step);
	bb.end.y   = iftMin(mimg->ysize-1, center.y + slic->step);
	
	if (slic->is_3D_img) {
	  bb.begin.z = iftMax(0, center.z - slic->step);
	  bb.end.z   = iftMin(mimg->zsize-1, center.z + slic->step);
	}
	else {
	  bb.begin.z = 0;
	  bb.end.z   = 0;
	}

	// for each voxel around the cluster center within a region of step x step
#pragma omp parallel for
	for (int z = bb.begin.z; z <= bb.end.z; z++) {
	  for (int y = bb.begin.y; y <= bb.end.y; y++) {
	    for (int x = bb.begin.x; x <= bb.end.x; x++) {
	      iftPoint u = {.x = x, .y = y, .z = z};
	      int p      = iftGetVoxelIndex(mimg, iftPointToVoxel(u));
	      
	      iftFColor cfeats_u = iftGetColorFeats(mimg, p, slic->img_max_range);
	      
	      double dist = iftSlicDistance(center, u, cfeats_center, cfeats_u, spatial_dist_factor,
					    slic->is_3D_img, slic->is_color_img);

	      // new minimum distance for the voxel p
	      // assigns the superpixel label for the voxel p
	      if (dist < min_dists[p]) {
		min_dists[p] = dist;
		superpixels->val[p] = cell+1; // the label from the cell (it starts from 1 to n_cells)
	      }
	    }
	  }
	}
      }
      iftUpdateCellCenters(slic, superpixels, mimg);
    }

    printf("slic->ncells: %d\n",slic->n_cells);

    // DESTROYERS
    iftFree(min_dists);

    return superpixels;
}


/**
 * @brief Initializes the initial cell's center and their feature vectors for an Input 2D Image
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 *
 * @note The initial cells are regular grids
 */
void iftInit2DSeeds(iftSlic *slic, const iftMImage *mimg, const iftImage *grad) {
    int n_cells_x = iftRound(mimg->xsize / (1.0 * slic->step));
    int n_cells_y = iftRound(mimg->ysize / (1.0 * slic->step));

    /*
    int n_regionsX;
    int n_regionsY;
    int rs = 2*slic->step;
    if ((mimg->xsize % rs) == 0){
        n_regionsX = mimg->xsize / rs;
    } else {
        n_regionsX = mimg->xsize / rs + 1;
    }
    if ((mimg->ysize % rs) == 0){
        n_regionsY = mimg->ysize / rs;
    } else {
        n_regionsY = mimg->ysize / rs + 1;
    }
    slic->region_map = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    slic->mini_map = iftCreateImage(n_regionsX,n_regionsX,1);
    slic->index_table = iftCreateImage(n_regionsX*n_regionsY,250,1);
    slic->table_size = iftAllocIntArray(n_regionsX*n_regionsY);

    for (int p=0; p<slic->region_map->n; p++) {
        iftVoxel u = iftGetVoxelCoord(slic->region_map,p);
        iftVoxel ur;
        ur.x = u.x/rs;
        ur.y = u.y/rs;
        ur.z = 0;
        int reg_num = iftGetVoxelIndex(slic->region_map,ur);
        slic->region_map->val[p] = reg_num;
    }
    */

    // num of remainder voxels in the axes - the regular grids could not cover all image domain
    int n_rem_voxels_x = mimg->xsize - slic->step * n_cells_x;
    int n_rem_voxels_y = mimg->ysize - slic->step * n_cells_y;

    // fixes the numbers of cells in the axes if they overflow their domains
    if (n_rem_voxels_x < 0) {
        n_cells_x--;
        n_rem_voxels_x = mimg->xsize - slic->step * n_cells_x;
    }
    if (n_rem_voxels_y < 0) {
        n_cells_y--;
        n_rem_voxels_y = mimg->ysize - slic->step * n_cells_y;
    }

    // the remainder voxels is distributed into the cells in each axis
    double n_add_voxels_per_cell_x = n_rem_voxels_x / (1.0 * n_cells_x);
    double n_add_voxels_per_cell_y = n_rem_voxels_y / (1.0 * n_cells_y);

    // offset = displacement from the first voxel in a cell (initially, it's a regular grid) to its center
    iftVoxel off;
    off.x = (slic->step + n_add_voxels_per_cell_x) / 2;
    off.y = (slic->step + n_add_voxels_per_cell_y) / 2;
    off.z = 0;

    slic->n_cells = n_cells_x * n_cells_y;
    iftAllocFeatVectorsForClusterCenters(slic);

    #pragma omp parallel for
    for (int cell_y = 0; cell_y < n_cells_y; cell_y++) {
        int acc_rem_y = cell_y * n_add_voxels_per_cell_y; // accumulated remainder voxel in y-axis

        for (int cell_x = 0; cell_x < n_cells_x; cell_x++) {
            int acc_rem_x = cell_x * n_add_voxels_per_cell_x; // accumulated remainder voxel in x-axis

            int cell = (cell_y * n_cells_x) + cell_x;

            // assigns center voxels for initial regular cells
            slic->spatial_feats[cell].x = (cell_x * slic->step) + (off.x + acc_rem_x);
            slic->spatial_feats[cell].y = (cell_y * slic->step) + (off.y + acc_rem_y);

            iftVoxel center = iftPointToVoxel(slic->spatial_feats[cell]);
            center.z = 0;

            int c = iftGetVoxelIndex(grad, center);

            // assigns the color space features for initial regular cells
            slic->color_feats[cell] = iftGetColorFeats(mimg, c, slic->img_max_range);
        }
    }

    /*
    for (int cell=0; cell< slic->n_cells; cell++){
        iftVoxel u, t;
        u.x = slic->spatial_feats[cell].x;
        u.y = slic->spatial_feats[cell].y;
        u.z = 0;
        int p = iftGetVoxelIndex(slic->region_map,u);
        int reg_num = slic->region_map->val[p];
        t.x = reg_num;
        t.y = slic->table_size[reg_num];
        t.z = 0;
        int pt = iftGetVoxelCoord(slic->region_map,t);
        slic->index_table->val[pt] = cell;
        slic->table_size[reg_num]++;
    }
     */

    iftFindLocalMinimumForSlicSeeds(slic, mimg, grad);
}


/**
 * @brief Initializes the initial cell's center and their feature vectors for an Input 3D Image 
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 *
 * @note The initial cells are regular grids
 */
void iftInit3DSeeds(iftSlic *slic, const iftMImage *mimg, const iftImage *grad) {
    int n_cells_x = iftRound(mimg->xsize / (1.0 * slic->step));
    int n_cells_y = iftRound(mimg->ysize / (1.0 * slic->step));
    int n_cells_z = iftRound(mimg->zsize / (1.0 * slic->step));

    // num of remainder voxels in the axes - the regular grids could not cover all image domain
    int n_rem_voxels_x = mimg->xsize - slic->step * n_cells_x;
    int n_rem_voxels_y = mimg->ysize - slic->step * n_cells_y;
    int n_rem_voxels_z = mimg->zsize - slic->step * n_cells_z;

    // fixes the numbers of cells in the axes if they overflow their domains
    if (n_rem_voxels_x < 0) {
        n_cells_x--;
        n_rem_voxels_x = mimg->xsize - slic->step * n_cells_x;
    }
    if (n_rem_voxels_y < 0) {
        n_cells_y--;
        n_rem_voxels_y = mimg->ysize - slic->step * n_cells_y;
    }
    if (n_rem_voxels_z < 0) {
        n_cells_z--;
        n_rem_voxels_z = mimg->zsize - slic->step * n_cells_z;
    }

    // the remainder voxels is distributed into the cells in each axis
    double n_add_voxels_per_cell_x = n_rem_voxels_x / (1.0 * n_cells_x);
    double n_add_voxels_per_cell_y = n_rem_voxels_y / (1.0 * n_cells_y);
    double n_add_voxels_per_cell_z = n_rem_voxels_z / (1.0 * n_cells_z);

    // offset = displacement from the first voxel in a cell (initially, it's a regular grid) to its center
    iftVoxel off;
    off.x = (slic->step + n_add_voxels_per_cell_x) / 2;
    off.y = (slic->step + n_add_voxels_per_cell_y) / 2;
    off.z = (slic->step + n_add_voxels_per_cell_z) / 2;

    slic->n_cells = n_cells_x * n_cells_y * n_cells_z;
    iftAllocFeatVectorsForClusterCenters(slic);

    // image with the center points
    // iftImage *cimg = iftCreateImage(grad->xsize, grad->ysize, grad->zsize);

    for (int cell_z = 0; cell_z < n_cells_z; cell_z++) {
        int acc_rem_z = cell_z * n_add_voxels_per_cell_z; // accumulated remainder voxel in z-axis
     
        for (int cell_y = 0; cell_y < n_cells_y; cell_y++) {
            int acc_rem_y = cell_y * n_add_voxels_per_cell_y; // accumulated remainder voxel in y-axis

            for (int cell_x = 0; cell_x < n_cells_x; cell_x++) {
                int acc_rem_x = cell_x * n_add_voxels_per_cell_x; // accumulated remainder voxel in x-axis

                int cell = cell_z * (n_cells_x * n_cells_y) + (cell_y * n_cells_x) + cell_x;

                // assigns center voxels for initial regular cells  
                slic->spatial_feats[cell].x = (cell_x * slic->step) + (off.x + acc_rem_x);
                slic->spatial_feats[cell].y = (cell_y * slic->step) + (off.y + acc_rem_y);
                slic->spatial_feats[cell].z = (cell_z * slic->step) + (off.z + acc_rem_z);

                iftVoxel center = iftPointToVoxel(slic->spatial_feats[cell]);
                int c           = iftGetVoxelIndex(grad, center); 
                // cimg->val[c] = 255;

                // assigns the color space features for initial regular cells
                slic->color_feats[cell] = iftGetColorFeats(mimg, c, slic->img_max_range);
            }
        }
    }
    iftFindLocalMinimumForSlicSeeds(slic, mimg, grad);
    // iftWriteImage(cimg, "centers.scn");
    // exit(-1);
}


/**
 * @brief Returns the Label from the First Adjacent Cell to a voxel u that was already Relabeled, or
 * or -1 otherwise.
 * @authors Samuel Martins, John Vargas
 * @date May 25, 2016
 *
 * @note Voxels from non-relabeled cells have value 0.
 * 
 * @param  u               Target voxel whose neighborhood is analyzed.
 * @param  new_superpixels Relabeled Superpixel Image.
 * @param  A               Adjacency Relation to the Neighborhood.
 * @return                 The label from the First Adjacent Relabeled Cell.
 */
int iftGetLabelFromAdjRelabeledCell(iftVoxel u, const iftImage *new_superpixels,
                                    const iftAdjRel *A) {
    int adj_label = -1;

    // gets the label of the first adjacent cell that was already relabeled
    for (int i = 1; i < A->n; i++) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, i);

        if (iftValidVoxel(new_superpixels, v)) {
            int q = iftGetVoxelIndex(new_superpixels, v);

            // if the adjacent voxel belongs to a relabeled cell 
            if (new_superpixels->val[q] != 0) {
                adj_label = new_superpixels->val[q]; // gets the label from the adjacent cell
                return adj_label;
            }
        }
    }

    return adj_label;
}


/**
 * @brief Relabels a Cell from <superpixels> to <new_superpixels>.
 *
 * Relabels all voxels from a given cell with the initial voxel <initial_voxel_idx> with the new label <new_label> \n
 * The parameter <idxs_in_cell> is a pre-allocated array with size superpixels->n where will store
 * all indices from the cell's voxels. \n
 * This arrays simulates a queue during the relabeling.
 * 
 * @authors Samuel Martins, John Vargas
 * @date May 25, 2016
 * 
 * @param  superpixels       Image with the old cell's labels. 
 * @param  new_superpixels   Image where the new cell's labels are stored.
 * @param  initial_voxel_idx Initial voxel index from the target cell.
 * @param  new_label         New label for the Cell.
 * @param  A                 Adjacency Relation that corresponds to the cell's voxel connectivity
 * @param  idxs_in_cell      Pre-allocated array of indices from the cell's voxels used as queue during cell relabeling.
 * @return                   The number of voxels from the cell.
 */
int iftCellRelabeling(const iftImage *superpixels, iftImage *new_superpixels, int initial_voxel_idx,
                      int new_label, const iftAdjRel *A, int *idxs_in_cell) {
    int old_label = superpixels->val[initial_voxel_idx];

    // inserts into Queue the initial cell's voxel
    int cell_size   = 1;
    idxs_in_cell[0] = initial_voxel_idx;

    for (int c = 0; c < cell_size; c++) {
        int p      = idxs_in_cell[c]; // removes from the Queue
        iftVoxel u = iftGetVoxelCoord(superpixels, p);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(superpixels, v)) {
                int q = iftGetVoxelIndex(superpixels, v);

                // if q is adjacent (with same old label) to p and it was not relabeled
                if ((superpixels->val[q] == old_label) && (new_superpixels->val[q] == 0)) {
                    new_superpixels->val[q] = new_label;
                    idxs_in_cell[cell_size] = q; // inserts into queue the adjacent voxel
                    cell_size++;
                }
            }
        }
    }

    return cell_size;
}
/******************************************************************/




/********************** PUBLIC FUNCTIONS *************************/
iftSlic *iftInitSlic(const iftImage *grad, const iftMImage *mimg, int input_n_cells, double comp,
                     int img_max_range) {
    // checkers
    if (grad == NULL)
        iftError("Gradient Image is NULL", "iftInitSlic");
    if (iftIs3DImage(grad)) {
      if ((grad->xsize < 3) || (grad->ysize < 3) || (grad->zsize < 3)) {
          iftError("Invalid Gradient Image Domain %dx%dx%d\n"        \
         "It must be at least 3x3x3", "iftInitSlic", grad->xsize, grad->ysize, grad->zsize);
      }
    }
    else {
      if ((grad->xsize < 3) || (grad->ysize < 3))
          iftError("Invalid Gradient Image Domain %dx%d\n"        \
         "It must be at least 3x3", "iftInitSlic", grad->xsize, grad->ysize);
    }
    if (input_n_cells <= 0)
        iftError("Invalid Input Number of Superpixels: %d... Try > 0",
                 "iftGenerateSuperpixelsBySlic", input_n_cells);
    if (comp <= 0.0)
        iftError("Invalid Compactness: %lf... Try > 0 (0.01, 0.1, 1, 10)",
                 "iftGenerateSuperpixelsBySlic", comp);
    if (grad->n < input_n_cells) {
      iftWarning("Input Number of Superpixels %d is greater than Image size: %d",
		 "iftGenerateSuperpixelsBySlic", input_n_cells, grad->n);
      input_n_cells = grad->n;
    }

    // initial regular grid (cluster) side = distance between two initial adjacent grid centers
    int step;
    if (iftIs3DImage(grad))
      step = pow(grad->n / (1.0 * input_n_cells), 1.0 / 3.0);
    else 
      step = iftRound(sqrtf(grad->n / (1.0 * input_n_cells)));
    
    bool is_color_img = (mimg->m == 3);
    iftSlic *slic     = iftCreateSlic(comp, step, img_max_range, iftIs3DImage(grad), is_color_img);

    // initializes the Slic Seeds
    if (iftIs3DImage(grad)) {
      iftInit3DSeeds(slic, mimg, grad);
    }
    else {
      iftInit2DSeeds(slic, mimg, grad);
    }

    return slic;
}


void iftDestroySlic(iftSlic **slic) {
    if (slic != NULL) {
        iftSlic *aux = *slic;
        if (aux->color_feats != NULL)
            iftFree(aux->color_feats);
        if (aux->spatial_feats != NULL)
            iftFree(aux->spatial_feats);
        if (aux->region_map != NULL)
            iftDestroyImage(&(aux->region_map));
        if (aux->mini_map != NULL)
            iftDestroyImage(&(aux->mini_map));
        if (aux->index_table != NULL)
            iftDestroyImage(&(aux->index_table));
        if (aux->table_size != NULL)
            iftFree(aux->table_size);


        iftFree(aux);
        *slic = NULL;
    }
}


// distance from voxel center to u (without the square roots for optimization)
// spatial dist factor = m^2 / S^2 in the formula of the Distance in the SLIC paper
double iftSlicDistance(iftPoint center, iftPoint u, iftFColor cfeats_center, iftFColor cfeats_u,
                       double spatial_dist_factor, bool is_3D_img, bool is_color_img) {
    // spatial distance from u to center to the power of 2
    double ds = ((center.x - u.x) * (center.x - u.x)) + ((center.y - u.y) * (center.y - u.y));
    if (is_3D_img)
        ds += ((center.z - u.z) * (center.z - u.z));

    // color dist to the power of 2
    double dc = (cfeats_center.val[0] - cfeats_u.val[0]) * (cfeats_center.val[0] - cfeats_u.val[0]);
    if (is_color_img)
        dc += (cfeats_center.val[1] - cfeats_u.val[1]) * (cfeats_center.val[1] - cfeats_u.val[1]) +
              (cfeats_center.val[2] - cfeats_u.val[2]) * (cfeats_center.val[2] - cfeats_u.val[2]);

    double dist = dc + (ds * spatial_dist_factor);

    return dist;
}


void iftEnforceLabelConnectivity(iftImage *superpixels, iftSlic *slic) {
    if (superpixels == NULL)
        iftError("Superpixel Image is NULL", "iftEnforceLabelConnectivity");

    // 4-neighborhood
    iftAdjRel *A = NULL;
    if (iftIs3DImage(superpixels))
      A = iftSpheric(1.0);
    else 
      A = iftCircular(1.0);

    // stores the index of each voxel from a given cell (at most n voxels) - it simulates a Queue
    int *idxs_in_cell         = iftAllocIntArray(superpixels->n);
    int min_cell_size         = superpixels->n / (4 * slic->n_cells);
    iftImage *new_superpixels = iftCreateImage(superpixels->xsize, superpixels->ysize, superpixels->zsize);

    int new_label = 1; // initial new label

    // relabels all voxels from the superpixels
    for (int p = 0; p < superpixels->n; p++) {
        // if the voxel belongs to a non-relabeled cell
        if (new_superpixels->val[p] == 0) {
            iftVoxel u    = iftGetVoxelCoord(superpixels, p);
            int adj_label = iftGetLabelFromAdjRelabeledCell(u, new_superpixels, A);

            // no adjacent relabeled cell was found
            if (adj_label == -1)
	      adj_label = new_label;

            int cell_size = iftCellRelabeling(superpixels, new_superpixels, p, new_label, A, idxs_in_cell);
	    
            // incorporates the cell into the chosen adjacent relabeled cell
            if (cell_size < min_cell_size) {
	      for (int c = 0; c < cell_size; c++) {
		int q = idxs_in_cell[c];
		
		new_superpixels->val[q] = adj_label;
	      }
	      new_label--; // it will keep the current new label for the next iteration
            }
            new_label++; // the entire cell has the voxel p was relabeled
        }
    }
    
    slic->n_cells = new_label;
    
    // assigns (in place) the new superpixels into the superpixel image
    iftFree(superpixels->val);
    superpixels->val     = new_superpixels->val;
    new_superpixels->val = NULL;
    

    // DESTROYERS
    iftDestroyAdjRel(&A);
    iftFree(idxs_in_cell);
    iftDestroyImage(&new_superpixels);
}


// TODO: change the accumulator to ulong instead of iftPoint and iftFColor
void iftUpdateCellCenters(iftSlic *slic, const iftImage *superpixels, const iftMImage *mimg) {
    // CHECKERS
    if (slic == NULL)
        iftError("Slic is NULL", "iftUpdateCellCenters");
    if (slic->n_cells <= 0) {
        printf("n->cells is equal to 0. Failed\n");
        return;
    }
    if (slic->spatial_feats == NULL)
        iftError("Slic Spatial Feats are NULL", "iftUpdateCellCenters");
    if (slic->color_feats == NULL)
        iftError("Slic Color Feats are NULL", "iftUpdateCellCenters");
    if (superpixels == NULL)
        iftError("Superpixel Image is NULL", "iftUpdateCellCenters");
    if (mimg == NULL)
        iftError("Multi-band Color Image (LAB or Grayscale) is NULL", "iftUpdateCellCenters");

    // creates the accumulator arrays including the LAB and z-axis
    double *acc_color_feats[3]   = {NULL, NULL, NULL};
    double *acc_spatial_feats[3] = {NULL, NULL, NULL};

    for (int i = 0; i < 3; i++) {
        acc_color_feats[i]   = iftAllocDoubleArray(slic->n_cells);
        acc_spatial_feats[i] = iftAllocDoubleArray(slic->n_cells);
    }
    int *n_voxels_per_cell = iftAllocIntArray(slic->n_cells);

    if (slic->is_color_img) {
        for (int p = 0; p < superpixels->n; p++) {
            int cell = superpixels->val[p] - 1; // e.g: the label 2 corresponds to the cell 1

            acc_color_feats[0][cell] += mimg->val[p][0]; // L
            acc_color_feats[1][cell] += mimg->val[p][1]; // A
            acc_color_feats[2][cell] += mimg->val[p][2]; // B
            
            // it is summing up the Z-axis, but it will be ignored if the oversegmented image is 2D
            iftVoxel u = iftGetVoxelCoord(superpixels, p);
            acc_spatial_feats[0][cell] += u.x;
            acc_spatial_feats[1][cell] += u.y;
            acc_spatial_feats[2][cell] += u.z;

            n_voxels_per_cell[cell]++;
        }
    }
    else {
        for (int p = 0; p < superpixels->n; p++) {
            int cell = superpixels->val[p] - 1; // e.g: the label 2 corresponds to the cell 1

            acc_color_feats[0][cell] += mimg->val[p][0]; // gray
            
            // it is summing up the Z-axis, but it will be ignored if the oversegmented image is 2D
            iftVoxel u = iftGetVoxelCoord(superpixels, p);
            acc_spatial_feats[0][cell] += u.x;
            acc_spatial_feats[1][cell] += u.y;
            acc_spatial_feats[2][cell] += u.z;

            n_voxels_per_cell[cell]++;
        }   
    }
    
    // updates the cell's center in the slic structure
    for (int cell = 0; cell < slic->n_cells; cell++) {
        // if the image is gray, the positions [1] and [2] will be ignored by the Slic Distance Function
        slic->color_feats[cell].val[0] = acc_color_feats[0][cell] / (1.0 * n_voxels_per_cell[cell]);
        slic->color_feats[cell].val[1] = acc_color_feats[1][cell] / (1.0 * n_voxels_per_cell[cell]);
        slic->color_feats[cell].val[2] = acc_color_feats[2][cell] / (1.0 * n_voxels_per_cell[cell]);

        // if the image is 3D, the z-axis will be ignored by the Slic Distance Function
        slic->spatial_feats[cell].x = acc_spatial_feats[0][cell] / (1.0 * n_voxels_per_cell[cell]);
        slic->spatial_feats[cell].y = acc_spatial_feats[1][cell] / (1.0 * n_voxels_per_cell[cell]);
        slic->spatial_feats[cell].z = acc_spatial_feats[2][cell] / (1.0 * n_voxels_per_cell[cell]);
    }


    // DESTROYERS
    for (int i = 0; i < 3; i++) {
        iftFree(acc_color_feats[i]);
        iftFree(acc_spatial_feats[i]);
    }
    iftFree(n_voxels_per_cell);
}


iftImage *iftGenerateSuperpixelsBySlic(const iftImage *img, const iftImage *input_grad, int input_n_cells,
                                       double comp, int img_max_range, int *out_n_cells) {
    // CHECKERS
  if (img == NULL)
      iftError("Image is NULL", "iftGenerateSuperpixelsBySlic");
  if (!iftIsIsotropic(img)) {
    char msg[512];
    if (iftIs3DImage(img))
      sprintf(msg, "Image should be isotropic (same pixel/voxel sizes)\n" \
	      "(dx, dy, dz) = (%.2f, %.2f, %.2f)", img->dx, img->dy, img->dz);
    else sprintf(msg, "Image should be isotropic (same pixel/voxel sizes)\n" \
		 "(dx, dy) = (%.2f, %.2f)", img->dx, img->dy);
      iftError(msg, "iftGenerateSuperpixelsBySlic");
  }

  // multi-band image that stores the LAB image or Grayscale image
  iftMImage *mimg = NULL;
  if (iftIsColorImage(img))
    mimg = iftImageToMImage(img, LABNorm_CSPACE); // 3 Bands - L, A, B
  else 
    mimg = iftImageToMImage(img, GRAY_CSPACE); 

  // uses the input gradient or computes another one if it is NULL 
  iftImage *grad = NULL;
  if (input_grad == NULL) {
      // Compute the Image Gradient
      iftAdjRel *A = NULL;
      
      // Adjacency of 3x3 (2D) or 3x3x3
      if (iftIs3DImage(img))
	A = iftSpheric(1.74); // sqrt(3.0)
      else 
	A = iftCircular(1.42); // sqrt(2.0)

      // we don't convert the color image space from YCbCr to LAB for reasons of simplicity
      grad = iftImageBasins(img, A);

      iftDestroyAdjRel(&A);
    }
    else 
      grad = iftCopyImage(input_grad);

  iftSlic *slic = iftInitSlic(grad, mimg, input_n_cells, comp, img_max_range);
  iftDestroyImage(&grad);
  iftImage *superpixels = iftPerformSuperpixelsSlic(slic, mimg);
  iftEnforceLabelConnectivity(superpixels, slic);
  iftCopyVoxelSize(img, superpixels);

  if (out_n_cells != NULL)
    *out_n_cells = slic->n_cells;
  
  // DESTROYERS
  iftDestroySlic(&slic);
  iftDestroyMImage(&mimg);
  
  return superpixels;
}


iftImage *iftGenerateSuperpixelsBySlicInBlocks(const iftImage *img, const iftImage *input_grad, int n_blocks,
                                               int input_n_cells, double comp, int img_max_range,
                                               iftImageTiles **out_blocks, int *out_n_cells) {
    if (n_blocks <= 0)
        iftError("Number of Blocks %d <= 0", "iftGenerateSuperpixelsBySlicInBlocks", n_blocks);
    if (img == NULL)
        iftError("Image is NULL", "iftGenerateSuperpixelsBySlicInBlocks");

    // tries to find regular blocks for the given number of blocks
    int block_vol = img->n / (1.0 * n_blocks);
    
    int n_blocks_x, n_blocks_y, n_blocks_z;
    iftNumberOfEquallyDimensionedTilesWithGivenMaximumSize(img->xsize, img->ysize, img->zsize,
                                                           block_vol, &n_blocks_x, &n_blocks_y, &n_blocks_z);

    iftImageTiles *blocks = iftImageTilesByEquallyDividingAxes(img, n_blocks_x, n_blocks_y, n_blocks_z);
    printf("--> Number of Blocks generated: %d\n", blocks->ntiles);
    iftImage *super_img   = iftCreateImage(img->xsize, img->ysize, img->zsize);
    
    for (int b = 0; b < blocks->ntiles; b++) {
      iftBoundingBox bb     = blocks->tile_coords[b];
      iftImage *block_img   = iftExtractROI(img, bb);
      iftImage *block_super = iftGenerateSuperpixelsBySlic(block_img, input_grad, input_n_cells,
                                                             comp, img_max_range, NULL);
      iftInsertROI(block_super, super_img, bb.begin); 

      iftDestroyImage(&block_img);
      iftDestroyImage(&block_super);
    }
    
    int n_cells;
    iftRenumerateBlockSuperpixels(super_img, blocks, &n_cells);
    
    if (out_blocks != NULL)
      *out_blocks = blocks;
    else iftDestroyImageTiles(&blocks);
    
    if (out_n_cells != NULL)
      *out_n_cells = n_cells;
    
    return super_img;
}
/******************************************************************/
