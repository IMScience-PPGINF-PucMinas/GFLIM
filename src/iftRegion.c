#include "iftRegion.h"

#include "ift/core/dtypes/FIFO.h"
#include "ift/core/io/Stream.h"
#include "iftSegmentation.h"
#include "iftGraphics.h"
#include "iftCommon.h"



/*----------------PRIVATE FUNCTIONS--------------- */
float iftGetStdevFromTwoSamples(float mean1, float mean2, float stdev1, float stdev2, float n1, float n2)
{
	float var1, var2, var_merge, mean_merge;
	var1 = stdev1*stdev1;
	var2 = stdev2*stdev2;
	mean_merge = (n1*mean1+n2*mean2)/(n1+n2);
	var_merge = (n1*(var1+mean1*mean1) + n2*(var2+mean2*mean2))/(n1+n2) - (mean_merge*mean_merge);
	return sqrt(var_merge);
}

/*------------PUBLIC FUNCTIONS------------------------*/
iftRegionGraph* iftRegionGraphFromLabelImage(iftImage* label_image, iftDataSet* dataset, iftAdjRel* adjacency){

	int n_regions = dataset->nsamples;
	iftRegionGraph* region_graph = iftCreateRegionGraph(dataset,n_regions);

#pragma omp parallel for
	for(int p = 0; p < label_image->n; p++){
		iftVoxel pixel = iftGetVoxelCoord(label_image,p);
		int pixel_label = label_image->val[p];

		int i;
		for(i = 1; i < adjacency->n; i++){
			iftVoxel neighbor = iftGetAdjacentVoxel(adjacency,pixel,i);

			if(iftValidVoxel(label_image, neighbor)){
				int neighbor_label = label_image->val[iftGetVoxelIndex(label_image,neighbor)];

				if(pixel_label != neighbor_label){
#pragma omp critical
					iftUnionSetElem(&(region_graph->node[pixel_label - 1].adjacent),neighbor_label - 1);
				}
			}
		}
	}

	return region_graph;
}

iftRegionGraph* iftCreateRegionGraph(iftDataSet* dataset, int n_regions){
	iftRegionGraph *region_graph = (iftRegionGraph*)iftAlloc(1,sizeof(iftRegionGraph));

	region_graph->dataset = dataset;
	region_graph->nnodes = n_regions;

	region_graph->node = (iftRegionGraphNode*)iftAlloc(n_regions, sizeof(iftRegionGraphNode));

	region_graph->root = iftAllocIntArray(n_regions);
	region_graph->pred = iftAllocIntArray(n_regions);
	region_graph->pathval = iftAllocFloatArray(n_regions);

	for(int r = 0; r < n_regions; r++){
		region_graph->root[r] = r;
		region_graph->pred[r] = IFT_NIL;
		region_graph->pathval[r] = IFT_INFINITY_FLT;
		region_graph->node[r].node_type=0;
	}

	region_graph->heap = iftCreateFHeap(n_regions, region_graph->pathval);

	return region_graph;
}

void iftDestroyRegionGraph(iftRegionGraph** rg){
	iftRegionGraph *region_graph = *rg;
	if(region_graph == NULL)
		return;

	int i;
	for(i = 0; i < region_graph->nnodes; i++){
		iftDestroySet(&(region_graph->node[i].adjacent));
	}

	iftFree(region_graph->node);
	iftFree(region_graph->root);
	iftFree(region_graph->pred);
	iftFree(region_graph->pathval);
	iftDestroyFHeap(&region_graph->heap);

	iftFree(region_graph);

	*rg = NULL;
}

void iftSuperpixelClassification(iftRegionGraph *graph, iftImage *label_image, iftLabeledSet *pixel_seeds){
	iftDataSet* dataset = graph->dataset;

	iftResetFHeap(graph->heap);

	int r;
	for(r = 0; r < graph->nnodes; r++){
		graph->pathval[r] = IFT_INFINITY_FLT;
		graph->root[r] = r;
		graph->pred[r] = IFT_NIL;
		dataset->sample[r].label = -1;
	}

	iftLabeledSet* seed = pixel_seeds;
	while(seed != NULL){
		r = label_image->val[seed->elem] - 1;
		if(dataset->sample[r].label == -1){
			graph->pathval[r] = 0;
			dataset->sample[r].label = seed->label;
			iftInsertFHeap(graph->heap,r);
		}
		seed = seed->next;
	}

	while( !iftEmptyFHeap(graph->heap)){
		r = iftRemoveFHeap(graph->heap);
		dataset->sample[r].weight = IFT_MAXWEIGHT * graph->pathval[r];

		iftSet* adj;
		for(adj = graph->node[r].adjacent; adj != NULL; adj = adj->next){
			int v = adj->elem;

			if(graph->heap->color[v] != IFT_BLACK){
				float weight = dataset->iftArcWeight(dataset->sample[r].feat,dataset->sample[v].feat, dataset->alpha,dataset->nfeats);

				float tmp = iftMax(graph->pathval[r], weight);
				if(tmp < graph->pathval[v]){
					dataset->sample[v].label = dataset->sample[r].label;
					graph->pathval[v] = tmp;
					graph->root[v] = graph->root[r];
					graph->pred[v] = r;

					if(graph->heap->color[v] == IFT_WHITE)
						iftInsertFHeap(graph->heap, v);
					else
						iftGoUpFHeap(graph->heap,graph->heap->pos[v]);
				}
			}
		}
	}

	iftResetFHeap(graph->heap);
}

//This function does not implement seed removal
void iftDiffSuperpixelClassification(iftRegionGraph *graph, iftImage *label_image, iftLabeledSet *new_seeds){
	iftDataSet *dataset = graph->dataset;
	int r;

	iftResetFHeap(graph->heap);

	iftLabeledSet *S = new_seeds;
	while(S != NULL){
		if(S->label == IFT_NIL)
            iftError("This function does not implement seed removal", "iftDiffSuperpixelClassification");

		r = label_image->val[S->elem] - 1;
		if(graph->heap->color[r] == IFT_WHITE){
			graph->pathval[r] = 0;
			graph->root[r] = r;
			graph->pred[r] = IFT_NIL;
			dataset->sample[r].label = S->label;

			iftInsertFHeap(graph->heap,r);
		}


		S = S->next;
	}

	while( !iftEmptyFHeap(graph->heap)){
		r = iftRemoveFHeap(graph->heap);
		dataset->sample[r].weight = IFT_MAXWEIGHT * graph->pathval[r];

		iftSet* adj;
		for(adj = graph->node[r].adjacent; adj != NULL; adj = adj->next){
			int v = adj->elem;

			if(graph->heap->color[v] != IFT_BLACK){
				float weight = dataset->iftArcWeight(dataset->sample[r].feat,dataset->sample[v].feat, dataset->alpha,dataset->nfeats);

				float tmp = iftMax(graph->pathval[r], weight);
				if(tmp < graph->pathval[v] || graph->pred[v] == r){
					dataset->sample[v].label = dataset->sample[r].label;
					graph->pathval[v] = tmp;
					graph->root[v] = graph->root[r];
					graph->pred[v] = r;

					if(graph->heap->color[v] == IFT_WHITE)
						iftInsertFHeap(graph->heap, v);
					else
						iftGoUpFHeap(graph->heap,graph->heap->pos[v]);
				}
			}
		}
	}

	iftResetFHeap(graph->heap);
}

//The range of the distance function NEEDS to be [0,1]
//This can be achieved by choosing the appropriate arc weight function for the dataset.
void iftFastSuperpixelClassification(iftRegionGraph *graph, iftImage *label_image, iftLabeledSet *pixel_seeds){
	iftDataSet *dataset = graph->dataset;

	int *pathval = iftAllocIntArray(graph->nnodes);

	iftGQueue *Q = iftCreateGQueue(IFT_MAXWEIGHT + 1, graph->nnodes, pathval);
	int r, p;
	for(r = 0; r < graph->nnodes; r++){
		pathval[r] = IFT_INFINITY_INT;
		graph->root[r] = r;
		graph->pred[r] = IFT_NIL;
		dataset->sample[r].label = -1;
	}

	iftLabeledSet *s = pixel_seeds;
	while(s != NULL){
		p = s->elem;

		int r = label_image->val[p] - 1;
		if(dataset->sample[r].label == -1){
			dataset->sample[r].label = s->label;
			pathval[r] = 0;

			iftInsertGQueue(&Q, r);
		}

		s = s->next;
	}

	while(!iftEmptyGQueue(Q)){
		r = iftRemoveGQueue(Q);
		dataset->sample[r].weight = pathval[r];

		iftSet *adjs = graph->node[r].adjacent;

		while(adjs != NULL){
			int s = adjs->elem;

			if(pathval[s] > pathval[r]){
				float weight = dataset->iftArcWeight(dataset->sample[r].feat,dataset->sample[s].feat, dataset->alpha,dataset->nfeats);

				int tmp = iftMax((int)(IFT_MAXWEIGHT * weight), pathval[r]);
				if(tmp < pathval[s]){
					dataset->sample[s].label = dataset->sample[r].label;
					pathval[s] = tmp;
					graph->root[s] = graph->root[r];
					graph->pred[s] = r;

					iftInsertGQueue(&Q,s);
				}
			}
			adjs = adjs->next;
		}
	}
}


void iftSuperpixelClassificationGeodesic(iftRegionGraph *graph, iftImage *label_image, iftLabeledSet *pixel_seeds, float beta){
	iftDataSet* dataset = graph->dataset;

	iftResetFHeap(graph->heap);

	int r;
	for(r = 0; r < graph->nnodes; r++){
		graph->pathval[r] = IFT_INFINITY_FLT;
		dataset->sample[r].label = -1;
	}

	iftLabeledSet* seed = pixel_seeds;
	while(seed != NULL){
		r = label_image->val[seed->elem] - 1;
		if(dataset->sample[r].label == -1){
			graph->pathval[r] = 0;
			dataset->sample[r].label = seed->label;
			iftInsertFHeap(graph->heap,r);
		}
		seed = seed->next;
	}

	float max_dist = IFT_INFINITY_FLT_NEG;
	float min_dist = IFT_INFINITY_FLT;
	for(r = 0; r < graph->nnodes; r++){
		iftSet *adj = graph->node[r].adjacent;
		while(adj){
			int v = adj->elem;
			float dist = dataset->iftArcWeight(dataset->sample[r].feat,dataset->sample[v].feat, dataset->alpha,dataset->nfeats);
			if(dist > max_dist)
				max_dist = dist;
			if(dist < min_dist)
				min_dist = dist;
			adj = adj->next;
		}
	}

	while( !iftEmptyFHeap(graph->heap)){
		r = iftRemoveFHeap(graph->heap);
		dataset->sample[r].weight = IFT_MAXWEIGHT * graph->pathval[r];

		iftSet* adj;
		for(adj = graph->node[r].adjacent; adj != NULL; adj = adj->next){
			int v = adj->elem;

			if(graph->heap->color[v] != IFT_BLACK){
				float dist = dataset->iftArcWeight(dataset->sample[r].feat,dataset->sample[v].feat, dataset->alpha,dataset->nfeats);
				dist = (dist - min_dist)/(max_dist - min_dist);

				float similarity = exp(-dist/0.5);
				dist = 1 - pow(similarity, beta);

				float tmp = graph->pathval[r] + dist;
				if(tmp < graph->pathval[v]){
					dataset->sample[v].label = dataset->sample[r].label;
					graph->pathval[v] = tmp;

					if(graph->heap->color[v] == IFT_WHITE)
						iftInsertFHeap(graph->heap, v);
					else
						iftGoUpFHeap(graph->heap,graph->heap->pos[v]);
				}
			}
		}
	}

	iftResetFHeap(graph->heap);

}

/* Given a label image with n labels different from 0, returns an image relabeled from 1 to n and leaves region 0 untouched */
iftImage* iftRelabelRegions(iftImage* labelled, iftAdjRel* adj_rel){

	iftImage *relabelled = iftCreateImage(labelled->xsize,labelled->ysize,labelled->zsize);
	iftCopyVoxelSize(labelled,relabelled);
	iftFIFO *F = iftCreateFIFO(labelled->n);

	int nlabels = 1;

	int i;
	for(i = 0; i < labelled->n; i++){
		if( (labelled->val[i] != 0) && (relabelled->val[i] == 0)){
			relabelled->val[i] = nlabels;
			iftInsertFIFO(F,i);

			while(!iftEmptyFIFO(F)){
				int p = iftRemoveFIFO(F);
				iftVoxel u = iftGetVoxelCoord(labelled,p);

				int j;
				for(j = 1; j < adj_rel->n; j++){
					iftVoxel v = iftGetAdjacentVoxel(adj_rel,u,j);

					if(iftValidVoxel(labelled,v)){
						int q = iftGetVoxelIndex(labelled,v);

//                        if((relabelled->val[q] == 0) && (labelled->val[p] == labelled->val[q]) ){
//                            relabelled->val[q] = nlabels;
//                            iftInsertFIFO(F,q);
//                        }

						if(relabelled->val[q] == 0){
                            if(labelled->val[p] == labelled->val[q])
                                iftInsertFIFO(F,q);
                            relabelled->val[q] = nlabels;
						}
					}
				}
			}

			nlabels++;
		}
	}

	iftDestroyFIFO(&F);

	return relabelled;
}

/* Receives an image labelled from 0 to n and returns a labeled set with the n geodesic centers, ordered by decreasing distance to their respective borders.*/
iftLabeledSet* iftGeodesicCenters(const iftImage* label_image){
	iftLabeledSet* geo_centers = NULL;

	int nregions = iftMaximumValue(label_image);

	iftImage *dist = iftCreateImage(label_image->xsize,label_image->ysize,label_image->zsize);
	iftImage *root = iftCreateImage(label_image->xsize,label_image->ysize,label_image->zsize);
	iftGQueue *Q = iftCreateGQueue(IFT_QSIZE,label_image->n,dist->val);

	iftAdjRel *adj_rel = NULL;
	if(iftIs3DImage(label_image))
		adj_rel = iftSpheric(sqrt(1.0));
	else
		adj_rel = iftCircular(sqrt(1.0));

	//Creates a distance map that is 0 in the (internal) borders and +infinity otherwise
	int p;
	for(p = 0; p < label_image->n; p++){
		dist->val[p] = IFT_INFINITY_INT;

		iftVoxel u = iftGetVoxelCoord(label_image, p);

		int i;
		for(i = 1; i < adj_rel->n; i++){
			iftVoxel v = iftGetAdjacentVoxel(adj_rel,u,i);

			if(iftValidVoxel(label_image,v)){
				int q = iftGetVoxelIndex(label_image,v);
				if((label_image->val[p] != label_image->val[q]) && (dist->val[p] > 0)){
					dist->val[p] = 0;
					root->val[p] = p;
					iftInsertGQueue(&Q,p);
				}
			}else if(dist->val[p] > 0){
				dist->val[p] = 0;
				root->val[p] = p;
				iftInsertGQueue(&Q,p);
			}
		}
	}

	iftDestroyAdjRel(&adj_rel);

	if(iftIs3DImage(label_image))
		adj_rel = iftSpheric(sqrtf(3.0));
	else
		adj_rel = iftCircular(sqrtf(2.0));

	while(!iftEmptyGQueue(Q)){
		p = iftRemoveGQueue(Q);

		iftVoxel u = iftGetVoxelCoord(label_image,p);
		iftVoxel r = iftGetVoxelCoord(label_image, root->val[p]);

		int i;
		for(i = 1; i < adj_rel->n; i++){
			iftVoxel v = iftGetAdjacentVoxel(adj_rel,u,i);

			if(iftValidVoxel(label_image,v)){
				int q = iftGetVoxelIndex(label_image,v);

				if ((dist->val[q] > dist->val[p]) && (label_image->val[p] == label_image->val[q])){
					int tmp = (v.x-r.x)*(v.x-r.x) + (v.y-r.y)*(v.y-r.y) + (v.z-r.z)*(v.z-r.z);
					if (tmp < dist->val[q]){
						if(dist->val[q] != IFT_INFINITY_INT)
						  iftRemoveGQueueElem(Q,q);

						dist->val[q] = tmp;
						root->val[q] = root->val[p];
						iftInsertGQueue(&Q,q);
					}
				}
			}
		}
	}

    /* printf("nregions %d\n", nregions); */
    /* 	iftWriteImageByExt(dist, "tmp/what.pgm"); */
	iftDestroyImage(&root);
	iftDestroyGQueue(&Q);
	iftDestroyAdjRel(&adj_rel);

	//Finds, for each region, the most distant pixel from its border
	int centers[nregions];
	int centers_distance[nregions];
	int i;
	for(i = 0; i < nregions; i++)
		centers_distance[i] = -1;

	for(p = 0; p < label_image->n; p++){
		int region = label_image->val[p] - 1;
		if(centers_distance[region] < dist->val[p]){
			centers[region] = p;
			centers_distance[region] = dist->val[p];
		}
	}

	iftDestroyImage(&dist);

	//This array will be scrambled by quicksort
	int centers_index[nregions];
	for(i = 0; i < nregions; i++)
		centers_index[i] = i;

	//Sorts the centers in increasing order, so inserting them in order in the iftLabeledSet geo_centers gives us a decreasing order
	iftQuickSort(centers_distance, centers_index, 0,nregions - 1, IFT_INCREASING);

	for(i = 0; i < nregions; i++){
		int index = centers_index[i];

		iftInsertLabeledSet(&geo_centers,centers[index],index + 1);
	}

	return geo_centers;
}


iftVoxel *iftGeodesicCenterVoxelCoords(const iftImage *super_img, int *out_n_clusters) {
	int n_clusters = iftMaximumValue(super_img);

	iftVoxel *roots            = (iftVoxel*) iftAlloc(n_clusters, sizeof(iftVoxel)); // one root for each supervoxel
    iftLabeledSet *geo_centers = iftGeodesicCenters(super_img);
    iftLabeledSet *center      = geo_centers;
    while (center != NULL) {
        int c         = center->label - 1; // ex: supervoxel with label 2 has idx [1]
        int voxel_idx = center->elem;
        roots[c]      = iftGetVoxelCoord(super_img, voxel_idx);
        center        = center->next;
    }
    iftDestroyLabeledSet(&geo_centers);

    if (out_n_clusters != NULL)
    	*out_n_clusters = n_clusters;

    return roots;
}


iftLabeledSet* iftGeometricCenters(iftImage* label_image){
	iftLabeledSet* centers = NULL;

	int nregions = iftMaximumValue(label_image);

	ulong *xc = iftAllocULongArray(nregions);
	ulong *yc = iftAllocULongArray(nregions);
	ulong *zc = iftAllocULongArray(nregions);

	int *np = iftAllocIntArray(nregions);

	int p;
	for(p = 0; p < label_image->n; p++){
		int region = label_image->val[p] - 1;

		if(region >= 0) {
			iftVoxel v = iftGetVoxelCoord(label_image, p);

			xc[region] += v.x;
			yc[region] += v.y;
			zc[region] += v.z;

			np[region]++;
		}
	}

	int i;
	for(i = 0; i < nregions; i++){
	  if(np[i] > 0) {
		iftVoxel v;
		v.x = xc[i]/np[i];
		v.y = yc[i]/np[i];
		v.z = zc[i]/np[i];

		iftInsertLabeledSet(&centers,iftGetVoxelIndex(label_image,v),i+1);
	  }
	}

	iftFree(xc);
	iftFree(yc);
	iftFree(zc);

	iftFree(np);

	return centers;
}

//Costly
iftImage* iftCreateRefinedLabelImage(iftImage* image, iftLabeledSet* seed, float spatial_radius, int volume_threshold, int steps, float vol_ratio){
		if(steps < 0 || vol_ratio >= 1)
            iftError("Invalid parameters.", "iftCreateRefinedLabelImage");

        iftAdjRel* adj_relation = iftCircular(spatial_radius);
        iftImage* basins_original = iftImageBasins(image,adj_relation);
        iftImage* basins = iftCopyImage(basins_original);

        iftImage* marker = iftVolumeClose(basins,volume_threshold,NULL);
		iftImage* label_image = iftWaterGray(basins,marker,adj_relation);

		int last_volume = volume_threshold;

		iftAdjRel *adj_relabeling = iftCircular(1.0);

        int i;
        for(i = 0; (i < steps) && (last_volume > 1) ; i++){
        	int nregions = iftMaximumValue(label_image);
        	int *region_class = iftAllocIntArray(nregions);

        	iftSet *multiple_labeled = 0;

        	iftLabeledSet* s = seed;
        	while(s != NULL){
        		int p = s->elem;
        		int c = s->label;
        		int region = label_image->val[p];

        		if(region_class[region - 1] != 0){
        			if(region_class[region - 1] != c)
        				iftInsertSet(&multiple_labeled,region);
        		}else{
        			region_class[region - 1] = c;
        		}

        		s = s->next;
        	}

        	if(multiple_labeled == 0){
        		iftFree(region_class);
        		break;
        	}

        	iftDestroyImage(&basins);
        	iftDestroyImage(&marker);

        	last_volume = iftMax(1, last_volume * vol_ratio);

            basins = iftCopyImage(basins_original);
	    marker = iftVolumeClose(basins, last_volume,NULL);
			iftImage *refined = iftWaterGray(basins,marker,adj_relation);

			int p;
			for(p = 0; p < image->n; p++){
				iftSet *r = multiple_labeled;
				while(r != NULL){
					if(r->elem == label_image->val[p]){
						label_image->val[p] = nregions + refined->val[p];
						break;
					}

					r = r->next;
				}
			}

			iftImage* old = label_image;
			label_image = iftRelabelRegions(label_image, adj_relabeling);
			iftDestroyImage(&old);

			iftDestroyImage(&refined);
			iftDestroySet(&multiple_labeled);
        }

        iftDestroyAdjRel(&adj_relabeling);
        iftDestroyAdjRel(&adj_relation);
        iftDestroyImage(&basins);
        iftDestroyImage(&basins_original);
        iftDestroyImage(&marker);

        return label_image;
}

void iftSetSuperpixelClassesFromGroundTruth(iftImage* label_image, iftImage* groud_truth, iftDataSet* dataset){
	int *n_obj_pixels = iftAllocIntArray(dataset->nsamples);
	int *n_bkg_pixels = iftAllocIntArray(dataset->nsamples);

	int p;
	for(p = 0; p < groud_truth->n; p++){
		int superpixel = label_image->val[p] - 1;

		if(groud_truth->val[p] != 0)
			n_obj_pixels[superpixel]++;
		else
			n_bkg_pixels[superpixel]++;
	}

	int r;
	for(r = 0; r < dataset->nsamples; r++){
		if(n_obj_pixels[r] > n_bkg_pixels[r])
			dataset->sample[r].truelabel = 2;
		else
			dataset->sample[r].truelabel = 1;
	}
}

iftImage* iftSuperpixelLabelImageFromDataset(iftImage *label_image, iftDataSet *dataset){
	iftImage *output = iftCreateImage(label_image->xsize,label_image->ysize,label_image->zsize);

	int p;
	for(p = 0; p < output->n; p++){
		int superpixel = label_image->val[p] - 1;
		output->val[p] = dataset->sample[superpixel].label;
	}

	return output;
}

//Finds the regions that are adjacent to regions with different label in the /graph/'s dataset
iftBMap		*iftAdjacentToDifferentLabels(iftRegionGraph *graph){
	iftSample *samples = graph->dataset->sample;
	iftBMap *selected = iftCreateBMap(graph->nnodes);

	int s;
	for(s = 0; s < graph->nnodes; s++){
		iftSet *adj = graph->node[s].adjacent;

		while(adj){
			int t = adj->elem;

			if(samples[s].label != samples[t].label){
				iftBMapSet1(selected,s);
				iftBMapSet1(selected,t);
			}

			adj = adj->next;
		}
	}

	return selected;
}

iftRegionHierarchyNode* iftCreateRegionHierarchyNode(int merging_time, iftRegionHierarchyNode* left, iftRegionHierarchyNode *right, int node_id){
	iftRegionHierarchyNode *node = (iftRegionHierarchyNode*)iftAlloc(1, sizeof(iftRegionHierarchyNode));

    node->id = node_id;
	node->merging_time = merging_time;
	node->left = left;
	node->right = right;
	node->region = IFT_NIL;

	node->xmax = iftMax(left->xmax, right->xmax);
	node->ymax = iftMax(left->ymax, right->ymax);
	node->zmax = iftMax(left->zmax, right->zmax);
	node->xmin = iftMin(left->xmin, right->xmin);
	node->ymin = iftMin(left->ymin, right->ymin);
	node->zmin = iftMin(left->zmin, right->zmin);

	return node;
}

iftRegionHierarchyNode* iftCreateRegionHierarchyLeaf(iftImage* label_image, iftSet* pixels, int region,
                                                     int node_id){
	iftRegionHierarchyNode *node = (iftRegionHierarchyNode*)iftAlloc(1, sizeof(iftRegionHierarchyNode));

    node->id = node_id;
	node->merging_time = 0;
	node->left = 0;
	node->right = 0;
	node->region = region;

	node->xmax = 0;
	node->ymax = 0;
	node->zmax = 0;
	node->xmin = IFT_INFINITY_INT;
	node->ymin = IFT_INFINITY_INT;
	node->zmax = IFT_INFINITY_INT;

	while(pixels){
		iftVoxel v = iftGetVoxelCoord(label_image, pixels->elem);

		node->xmax = iftMax(node->xmax, v.x);
		node->ymax = iftMax(node->ymax, v.y);
		node->zmax = iftMax(node->zmax, v.z);
		node->xmin = iftMin(node->xmin, v.x);
		node->ymin = iftMin(node->ymin, v.y);
		node->zmin = iftMin(node->zmin, v.z);

		pixels = pixels->next;
	}

	return node;
}


void iftDestroyRegionHierarchyNode(iftRegionHierarchyNode **region_hnode){
	iftRegionHierarchyNode *rh_node = *region_hnode;
	if(rh_node){
		iftDestroyRegionHierarchyNode(&rh_node->left);
		iftDestroyRegionHierarchyNode(&rh_node->right);

		iftFree(rh_node);
		*region_hnode = 0;
	}
}

void iftDestroyRegionHierarchy(iftRegionHierarchy** region_h){
	iftRegionHierarchy *rh = *region_h;
	if(rh){
		iftDestroyImage(&rh->label_image);

		int i;
		for(i = 0; i < rh->nleaves; i++){
			iftDestroySet(&rh->pixels_by_node[i]);
		}

		iftFree(rh->pixels_by_node);
		iftDestroyRegionHierarchyNode(&rh->root);
		iftFree(rh);

		*region_h = 0;
	}
}

iftRegionHierarchy* iftCreateRegionHierarchy(iftDataSet *dataset, iftAdjRel* adjacency, iftMergingFun merge_func){

	iftImage *label_image = (iftImage*)dataset->ref_data;
	iftRegionHierarchy* rh = (iftRegionHierarchy*)iftAlloc(1, sizeof(iftRegionHierarchy));

	//Creates a region graph
	iftRegionGraph *region_graph = iftRegionGraphFromLabelImage(label_image, dataset, adjacency);
	rh->nleaves = region_graph->nnodes;

	//Stores the pixels belonging to each region
	rh->pixels_by_node = (iftSet**)iftAlloc(region_graph->nnodes,sizeof(iftSet*));
	int p;
	for(p = 0; p < label_image->n; p++){
		iftInsertSet(&rh->pixels_by_node[label_image->val[p] - 1],p);
	}
	rh->label_image = iftCopyImage(label_image);

	iftRegionHierarchyNode **hierarchy_nodes = (iftRegionHierarchyNode**)iftAlloc(region_graph->nnodes,sizeof(iftRegionHierarchyNode*));
	int r;
	int node_id = 1;
	for(r = 0; r < region_graph->nnodes; r++){
		hierarchy_nodes[r] = iftCreateRegionHierarchyLeaf(label_image, rh->pixels_by_node[r], r, node_id++);
	}

	iftResetFHeap(region_graph->heap);

	//Finds the distance to the nearest neighbor, stores it in "pathval" and the nearest neighbor in "destination"
	int *dest = iftAllocIntArray(region_graph->nnodes);
	for(r = 0; r < region_graph->nnodes; r++){
		iftSet *adj = region_graph->node[r].adjacent;

		region_graph->pathval[r] = IFT_INFINITY_FLT;
		dest[r] = -1;

		while(adj){
			int neighbor = adj->elem;

			float dist = dataset->iftArcWeight(region_graph->dataset->sample[r].feat,region_graph->dataset->sample[neighbor].feat,dataset->alpha,region_graph->dataset->nfeats);
			if(region_graph->pathval[r] > dist){
				region_graph->pathval[r] = dist;
				dest[r] = neighbor;
			}

			adj = adj->next;
		}

		if(dest[r] == -1)
            iftError("The graph is disconnected", "iftRegionMerging");

		iftInsertFHeap(region_graph->heap,r);
	}

	//Repeat until we achieve the desired number of iterations
	int i;
	int u;
	for(i = 0 ;; i++){
		//The vertices u and v will be merged
		u = iftRemoveFHeap(region_graph->heap);
		int v = dest[u];

		//Creates a new node in the hierarchy
		iftRegionHierarchyNode *new_node = iftCreateRegionHierarchyNode(i+1, hierarchy_nodes[u], hierarchy_nodes[v], node_id++);
		hierarchy_nodes[u] = new_node;

		//Remove the node v from consideration
		dest[v] = -1;
		iftRemoveFHeapElem(region_graph->heap, v);

		//Compute the features for the new node
		float *old_feats_u = region_graph->dataset->sample[u].feat;
		region_graph->dataset->sample[u].feat = merge_func(region_graph->dataset->sample[u].feat,region_graph->dataset->sample[v].feat, dataset->alpha, region_graph->dataset->nfeats);
		iftFree(old_feats_u);

		//Nobody should be adjacent to v
		iftSet *v_adj = region_graph->node[v].adjacent;
		while(v_adj){
			iftRemoveSetElem(&region_graph->node[v_adj->elem].adjacent,v);
			v_adj = v_adj->next;
		}

		//Removes u from v's adjacents
		iftRemoveSetElem(&region_graph->node[v].adjacent,u);
		//Joins u's and v's adjacents
		iftSet *old_adj_u = region_graph->node[u].adjacent;
		region_graph->node[u].adjacent = iftSetUnion(region_graph->node[u].adjacent,region_graph->node[v].adjacent);
		iftDestroySet(&old_adj_u);

		//Removes v's adjacents
		iftDestroySet(&region_graph->node[v].adjacent);

		//If u is adjacent to w, w is adjacent to u (it would be faster to consider only v's adjacents)
		iftSet* u_adj = region_graph->node[u].adjacent;
		while(u_adj){
			if(!iftSetHasElement(region_graph->node[u_adj->elem].adjacent, u))
				iftInsertSet(&region_graph->node[u_adj->elem].adjacent, u);

			u_adj = u_adj->next;
		}

		//The nearest node to u will be recomputed
		region_graph->pathval[u] = IFT_INFINITY_FLT;

		u_adj = region_graph->node[u].adjacent;
		while(u_adj){
			//w is adjacent to u
			int w = u_adj->elem;

			//Finds the nearest adjacent to u
			float dist = dataset->iftArcWeight(region_graph->dataset->sample[w].feat,region_graph->dataset->sample[u].feat,dataset->alpha,region_graph->dataset->nfeats);
			if(region_graph->pathval[u] > dist){
				region_graph->pathval[u] = dist;
				dest[u] = w;
			}

			//Finds the nearest to w
			iftSet *w_adj = region_graph->node[w].adjacent;
			iftRemoveFHeapElem(region_graph->heap, w);
			region_graph->pathval[w] = IFT_INFINITY_FLT;
			while(w_adj){
				dist = dataset->iftArcWeight(region_graph->dataset->sample[w].feat, region_graph->dataset->sample[w_adj->elem].feat,dataset->alpha,region_graph->dataset->nfeats);
				if(region_graph->pathval[w] > dist){
					region_graph->pathval[w] = dist;
					dest[w] = w_adj->elem;
				}

				w_adj = w_adj->next;
			}
			iftInsertFHeap(region_graph->heap,w);

			u_adj = u_adj->next;
		}

		//Tests whether there is only one region left
		if(iftEmptyFHeap(region_graph->heap)){
			break;
		}

		//Re-inserts u in the queue
		iftInsertFHeap(region_graph->heap, u);
	}

	iftFree(dest);

	//Sets the root of the hierarchy
	rh->root = hierarchy_nodes[u];


	//Frees the array, not its contents (which now belong to rh)
	iftFree(hierarchy_nodes);
	iftDestroyRegionGraph(&region_graph);

	return rh;
}

//Private function: draw label on pixels recursively
void iftRHDrawSubregions(iftRegionHierarchy *rh, iftRegionHierarchyNode *node, iftImage* label_image, int label){
	if(node->region != IFT_NIL){
		iftSet *pixels = rh->pixels_by_node[node->region];

		while(pixels){
			int p = pixels->elem;
			label_image->val[p] = label;
			pixels = pixels->next;
		}
	}else{
		iftRHDrawSubregions(rh, node->left, label_image, label);
		iftRHDrawSubregions(rh, node->right, label_image, label);
	}
}

//Private function: flattens the hierarchy recursively
int iftRHFlatten(iftRegionHierarchy *rh, iftRegionHierarchyNode* node, iftImage *label_image, int cutoff, int current_label){
	if(node->merging_time <= cutoff){
		iftRHDrawSubregions(rh, node, label_image, current_label);
		current_label++;
	}else{
		current_label = iftRHFlatten(rh, node->left, label_image, cutoff, current_label);
		current_label = iftRHFlatten(rh, node->right, label_image, cutoff, current_label);
	}

	return current_label;
}

iftImage* iftFlattenRegionHierarchy(iftRegionHierarchy* rh, int nregions){
	if(nregions >= rh->nleaves)
		return iftCopyImage(rh->label_image);

	iftImage *label_image = iftCreateImage(rh->label_image->xsize,rh->label_image->ysize,rh->label_image->zsize);
	iftCopyVoxelSize(rh->label_image,label_image);
	if(nregions <= 1){
		iftSetImage(label_image, 1);
		return label_image;
	}

	//Nodes with merging_time <= cutoff will create a new region
	int cutoff = rh->nleaves - nregions;

	iftRegionHierarchyNode *root = rh->root;

	int current_label = 1;
	current_label = iftRHFlatten(rh, root->left, label_image, cutoff, current_label);
	current_label = iftRHFlatten(rh, root->right, label_image, cutoff, current_label);

	//These tests are somewhat expensive
	//if(iftMinimumValue(label_image) == 0)
	//	iftError("Inconsistent label image","iftFlattenRegionHierarchy");
//	if(iftMaximumValue(label_image) != nregions || current_label != nregions + 1)
//		iftError("Inconsistent number of regions", "iftFlattenRegionHierarchy");
	if(current_label != nregions + 1)
        iftError("Inconsistent number of regions", "iftFlattenRegionHierarchy");

	return label_image;
}

iftImage* iftEliminateRegionsByArea(iftDataSet* dataset, iftAdjRel *adj, int threshold){
	iftImage *label_image = (iftImage*)dataset->ref_data;

	int nregions = iftMaximumValue(label_image);
	int *npixels = iftAllocIntArray(nregions);
	int *representatives = iftAllocIntArray(nregions);

	int p;
	for(p = 0; p < label_image->n; p++){
		int region = label_image->val[p] - 1;

		if(npixels[region] == 0)
			representatives[region] = p;

		npixels[region]++;
	}

	iftLabeledSet *pixel_seeds = 0;
	int n_output_regions = 0;
	int r;
	for(r = 0; r < nregions; r++){
		if(npixels[r] >= threshold){
			n_output_regions++;
			iftInsertLabeledSet(&pixel_seeds, representatives[r], n_output_regions);
		}
	}

	if (n_output_regions == 0){
		iftImage *single_region = iftCreateImage(label_image->xsize, label_image->ysize, label_image->zsize);
		iftSetImage(single_region, 1);
		return single_region;
	}

	iftRegionGraph *rg = iftRegionGraphFromLabelImage(label_image, dataset, adj);
	iftSuperpixelClassification(rg, label_image, pixel_seeds);

	iftFree(npixels);
	iftFree(representatives);
	iftDestroyLabeledSet(&pixel_seeds);
	iftDestroyRegionGraph(&rg);

	return iftSuperpixelLabelImageFromDataset(label_image, dataset);
}

iftRegionHierarchyNode** iftGetRegionHierarchyNodes(iftRegionHierarchy *rh){
	int nnodes = rh->root->merging_time + rh->nleaves;
	iftRegionHierarchyNode** rha = (iftRegionHierarchyNode**)iftAlloc(nnodes, sizeof(iftRegionHierarchyNode*));

	rha[0] = rh->root;
	int next_i = 1;

	int i;
	for(i = 0; i < nnodes; i++){
        if( i >= next_i )
            iftError("Inconsistent region hierarchy", "iftGetRegionHierarchyNodes");

		iftRegionHierarchyNode* node = rha[i];

		if(node->region == IFT_NIL){
			rha[next_i] = node->left;
			next_i++;
			rha[next_i] = node->right;
			next_i++;
		}
	}

	return rha;
}

void iftSetRegionGraphNodeTypeForBorderSuperpixelsInTiles(iftRegionGraph *region_graph, iftImage *label, iftImage *tiles_img, iftAdjRel *A) {

 #pragma omp parallel for
 	for (int p = 0; p < label->n; ++p) {
 		iftVoxel u = iftGetVoxelCoord(label, p);
 		iftVoxel v;
		for (int i = 1; i < A->n; ++i) {
			v = iftGetAdjacentVoxel(A, u, i);
			if (iftValidVoxel(label,v)) {
				int q = iftGetVoxelIndex(label, v);
				if (tiles_img->val[p] != tiles_img->val[q] && label->val[p] != label->val[q]) {
#pragma omp critical
					region_graph->node[label->val[p] - 1].node_type = 1;
					break;
				}
			}
		}
 	}

}

iftRegionHierarchy* iftMergeBorderRegions(iftDataSet *dataset, iftRegionGraph *region_graph, iftMergingFun merge_func, float threshold,iftRegionHierarchyNode ***output_hierarchy_nodes){

	iftImage *label_image = (iftImage*)dataset->ref_data;
	iftRegionHierarchy* rh = (iftRegionHierarchy*)iftAlloc(1, sizeof(iftRegionHierarchy));

	//Creates a region graph
	rh->nleaves = region_graph->nnodes;

	//Stores the pixels belonging to each region
	rh->pixels_by_node = (iftSet**)iftAlloc(region_graph->nnodes,sizeof(iftSet*));
	int p;

	for(p = 0; p < label_image->n; p++){
		iftInsertSet(&rh->pixels_by_node[label_image->val[p] - 1],p);
	}
	rh->label_image = iftCopyImage(label_image);

	iftRegionHierarchyNode **hierarchy_nodes = (iftRegionHierarchyNode**)iftAlloc(region_graph->nnodes,sizeof(iftRegionHierarchyNode*));
	int r;
	int node_id = 1;
	for(r = 0; r < region_graph->nnodes; r++){
		hierarchy_nodes[r] = iftCreateRegionHierarchyLeaf(label_image, rh->pixels_by_node[r], r, node_id++);
	}

	iftResetFHeap(region_graph->heap);

	//Finds the distance to the nearest neighbor, stores it in "pathval" and the nearest neighbor in "destination"
	int *dest = iftAllocIntArray(region_graph->nnodes);
	for(r = 0; r < region_graph->nnodes; r++) {

		region_graph->pathval[r] = IFT_INFINITY_FLT;
		dest[r] = -1;
		if (region_graph->node[r].node_type == 1) {
			iftSet *adj = region_graph->node[r].adjacent;

			while (adj) {
				int neighbor = adj->elem;

				if (region_graph->node[neighbor].node_type == 1) {
					float dist = dataset->iftArcWeight(region_graph->dataset->sample[r].feat,
																						 region_graph->dataset->sample[neighbor].feat, dataset->alpha,
																						 region_graph->dataset->nfeats);
					if (region_graph->pathval[r] > dist) {
						region_graph->pathval[r] = dist;
						dest[r] = neighbor;
					}
				}

				adj = adj->next;
			}

//			if (dest[r] == -1)
//				iftError("The graph is disconnected", "iftRegionMerging");

			iftInsertFHeap(region_graph->heap, r);
		}
	}

	//Repeat until we achieve the desired number of iterations
	int i;
	int u;
	for(i = 0 ;; i++){

//		printf("iteration: %d\n",i);

		//The vertices u and v will be merged
		u = iftRemoveFHeap(region_graph->heap);
		int v = dest[u];

		if (region_graph->pathval[u] > threshold)
			break;

		//Creates a new node in the hierarchy
		iftRegionHierarchyNode *new_node = iftCreateRegionHierarchyNode(i+1, hierarchy_nodes[u], hierarchy_nodes[v], node_id++);
		hierarchy_nodes[u] = new_node;

		//Remove the node v from consideration
		dest[v] = -1;
		iftRemoveFHeapElem(region_graph->heap, v);

		//Compute the features for the new node
//		float *old_feats_u = region_graph->dataset->sample[u].feat;
//		iftFree(old_feats_u);
		region_graph->dataset->sample[u].feat = merge_func(region_graph->dataset->sample[u].feat,region_graph->dataset->sample[v].feat, dataset->alpha, region_graph->dataset->nfeats);


		//Nobody should be adjacent to v
		iftSet *v_adj = region_graph->node[v].adjacent;
		while(v_adj){
			iftRemoveSetElem(&region_graph->node[v_adj->elem].adjacent,v);
			v_adj = v_adj->next;
		}

		//Removes u from v's adjacents
		iftRemoveSetElem(&region_graph->node[v].adjacent,u);
		//Joins u's and v's adjacents
		iftSet *old_adj_u = region_graph->node[u].adjacent;
		region_graph->node[u].adjacent = iftSetUnion(region_graph->node[u].adjacent,region_graph->node[v].adjacent);
		iftDestroySet(&old_adj_u);

		//Removes v's adjacents
		iftDestroySet(&region_graph->node[v].adjacent);

		//If u is adjacent to w, w is adjacent to u (it would be faster to consider only v's adjacents)
		iftSet* u_adj = region_graph->node[u].adjacent;
		while(u_adj){
			if(!iftSetHasElement(region_graph->node[u_adj->elem].adjacent, u))
				iftInsertSet(&region_graph->node[u_adj->elem].adjacent, u);

			u_adj = u_adj->next;
		}

		//The nearest node to u will be recomputed
		region_graph->pathval[u] = IFT_INFINITY_FLT;

		u_adj = region_graph->node[u].adjacent;
		while(u_adj){
			//w is adjacent to u
			int w = u_adj->elem;
			float dist;
			if (region_graph->node[w].node_type == 1) {
				//Finds the nearest adjacent to u
				dist = dataset->iftArcWeight(region_graph->dataset->sample[w].feat, region_graph->dataset->sample[u].feat,
																		 dataset->alpha, region_graph->dataset->nfeats);
				if (region_graph->pathval[u] > dist) {
					region_graph->pathval[u] = dist;
					dest[u] = w;
				}

				//Finds the nearest to w
				iftSet *w_adj = region_graph->node[w].adjacent;
				iftRemoveFHeapElem(region_graph->heap, w);
				region_graph->pathval[w] = IFT_INFINITY_FLT;
				while (w_adj) {
					if (region_graph->node[w_adj->elem].node_type == 1) {
						dist = dataset->iftArcWeight(region_graph->dataset->sample[w].feat,
																				 region_graph->dataset->sample[w_adj->elem].feat, dataset->alpha,
																				 region_graph->dataset->nfeats);
						if (region_graph->pathval[w] > dist) {
							region_graph->pathval[w] = dist;
							dest[w] = w_adj->elem;
						}

					}
					w_adj = w_adj->next;
				}
				iftInsertFHeap(region_graph->heap, w);
			}

			u_adj = u_adj->next;
		}

		//Tests whether there is only one region left
		if(iftEmptyFHeap(region_graph->heap)){
			break;
		}

		//Re-inserts u in the queue
		iftInsertFHeap(region_graph->heap, u);
	}

	iftFree(dest);

	//Sets the root of the hierarchy
	rh->root = hierarchy_nodes[u];
	// set output array of hierarchy nodes
	*output_hierarchy_nodes = hierarchy_nodes;

	return rh;
}

float* iftMergeLabColorMeanStdAndSizeFeatures(float *f1, float *f2, float *alpha, int n)
{
	float *merged_feat;
	int nbands = 3;
	assert(n == (2*nbands+1));
	merged_feat = iftAllocFloatArray(n);
	// Compute merge features
	for (int i = 0; i < nbands; ++i)
	{
		merged_feat[i] = (f1[2*nbands]*f1[i]+f2[2*nbands]*f2[i])/(f1[2*nbands]+f2[2*nbands]);
		merged_feat[i+nbands] = iftGetStdevFromTwoSamples(f1[i], f2[i], f1[i+nbands], f2[i+nbands], f1[2*nbands], f2[2*nbands]);
	}
	merged_feat[2*nbands] = f1[2*nbands] + f2[2*nbands];
	return merged_feat;
}

iftImage *iftJoinProbabilityRegions(iftMImage *prob, iftImage *label, int norm_val, bool decrement)
{
	iftImage *out = iftCreateImage(label->xsize, label->ysize, label->zsize);

	int n_labels = iftMaximumValue(label);

	if (decrement == false)
		n_labels++;


	if (prob->m != n_labels) {
        iftError("Number of labels must equal to the probability map number of bands (classes)",
                 "iftJoinProbabilityRegions");
	}

	for (int p = 0; p < label->n; p++)
	{
		int idx = label->val[p] - decrement;
		if (idx >= 0) /* this add the possibility to penalize a region */
		{
			out->val[p] = (1.0f - prob->val[p][idx]) * norm_val;
		} else {
			out->val[p] = norm_val;
		}
	}

	return out;
}

