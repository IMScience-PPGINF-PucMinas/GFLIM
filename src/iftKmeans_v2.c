#include "iftKmeans_v2.h"
#include "ift/core/dtypes/IntArray.h"

/**
 * @brief Check if centroids have converged according to given tolerance.
 * @param iftDataSet *Z --- input dataset
 * @param curr_centroids --- current centroids
 * @param prev_centroids --- previous centroids
 * @param tol --- tolerance of the difference of cluster centers of two consecutive iterations
 * @param distance function in the feature space. It assumes the
 * function was already chosen.
 * @author David Aparco Cardenas 
 * @date Jun 2nd, 2020
 */

bool iftCentroidsConverged(iftDataSet *Z, float **curr_centroids, float **prev_centroids, float tol, iftKmeansDistance distance) {
    int n_clusters = Z->ngroups;
    int n_features = Z->nfeats;
    float dist;

    for (int i = 0; i < n_clusters; i++) {    
        dist = distance(curr_centroids[i], prev_centroids[i], n_features);
        if (dist > tol) {
            return false;
        }
    }

    return true;
}

/* -------------PUBLIC FUNCTIONS--------------------------------------- */

void iftKMeans(iftDataSet *Z, int n_clusters, int max_iter, float tol, int *init_centroids, iftCentroidsHist *centroids_hist, iftKmeansDistance distance) {
    
    if (distance == NULL)
        distance = iftEuclideanDistance;

    if (Z == NULL) {
        iftError("Dataset is NULL", "iftKMeans");
    }

    if (Z->nsamples == 0) {
        iftError("Dataset is empty", "iftKMeans");
    }

    if (n_clusters > Z->nsamples) {
        iftError("Number of clusters is greater than number of samples", "iftKMeans");
    }

    if (tol <= 0) {
        iftError("Tolerance must be greater than zero", "iftKMeans");
    }

    Z->ngroups = n_clusters;

    // Write validation for init_centroids
    if (init_centroids == NULL) {
        init_centroids = iftRandomIntegers(0, Z->nsamples - 1, n_clusters);
    }

    if (centroids_hist != NULL) {
        centroids_hist->n_clusters = n_clusters;
        centroids_hist->centroid   = iftAlloc(n_clusters, sizeof(*centroids_hist->centroid));
        centroids_hist->n_feats    = Z->nfeats;
    }

    float d, d_min;
    int s, c_min;
    int iter               = 0;
    int n_feats            = Z->nfeats;
    float *curr_distances  = iftAlloc(Z->nsamples, sizeof(*curr_distances));
    float **curr_centroids = iftAlloc(n_clusters, sizeof(*curr_centroids));
    float **prev_centroids = iftAlloc(n_clusters, sizeof(*prev_centroids));

    for (int i = 0; i < n_clusters; i++) {
        s                 = init_centroids[i];
        curr_centroids[i] = iftAlloc(n_feats, sizeof(*curr_centroids[i]));
        prev_centroids[i] = iftAlloc(n_feats, sizeof(*prev_centroids[i]));

        if (centroids_hist != NULL) {
            centroids_hist->centroid[i]            = iftAlloc(1, sizeof(*centroids_hist->centroid[i]));
            centroids_hist->centroid[i]->id        = i + 1;
            centroids_hist->centroid[i]->n_hist    = 0;
            centroids_hist->centroid[i]->feat      = NULL;
            centroids_hist->centroid[i]->hist_feat = NULL;
        }
    }
    
    iftDynamicSet **S = iftAlloc(n_clusters, sizeof(*S));
    
    // Initializing centroids
    for (int i = 0; i < n_clusters; i++) {
        s = init_centroids[i];
        iftCopyFloatArray(curr_centroids[i], Z->sample[s].feat, n_feats);
        if (centroids_hist != NULL) {
            iftInsertCentroidHistory(centroids_hist->centroid[i], Z->sample[s].feat, n_feats);
        }
    }

    while (iter < max_iter && !iftCentroidsConverged(Z, curr_centroids, prev_centroids, tol, distance)) {
        // Initializing dynamic sets (groups)
        for (int i = 0; i < n_clusters; i++) {
            S[i] = iftCreateDynamicSet(n_feats);
            if (iter == 0) {
                s = init_centroids[i];
                Z->sample[s].group = i + 1;
                iftInsertSampleDynamicSet(Z, S[i], s);
            }
        }
        // Assigning samples to dynamic sets
        for (s = 0; s < Z->nsamples; s++) {
            d_min = IFT_INFINITY_FLT;
            c_min = IFT_NIL;
            
            for (int i = 0; i < n_clusters; i++) {
                d = distance(curr_centroids[i], Z->sample[s].feat, n_feats);
                if (d < d_min) {
                    d_min = d;
                    c_min = i;
                }
            }
            curr_distances[s]  = d_min;
            Z->sample[s].group = c_min + 1;
            iftInsertSampleDynamicSet(Z, S[c_min], s);
        }

        // Updating centroids
        for (int i = 0; i < n_clusters; i++) {
            // Updating previous centroids
            iftCopyFloatArray(prev_centroids[i], curr_centroids[i], n_feats);
            // Updating current centroids
            iftCopyDblArrayToFloatArray(curr_centroids[i], S[i]->mean, n_feats);
            if (centroids_hist != NULL) {
                iftInsertCentroidHistory(centroids_hist->centroid[i], curr_centroids[i], n_feats);
            }
        }
        
        iter += 1;
    }

    // Set the status IFT_PROTOTYPE to samples closest to the actual centroids
    if (S[0] != NULL) {
        iftSetDataSetPrototypeStatus(Z, S);
    } else {
        for (int s=0; s < Z->nsamples; s++)
	        iftRemoveSampleStatus(&Z->sample[s], IFT_PROTOTYPE);

        for (int i = 0; i < n_clusters; i++) {
	        s = init_centroids[i];
	        iftAddSampleStatus(&Z->sample[s], IFT_PROTOTYPE);
        }
        
        for (int s = 0; s < Z->nsamples; s++) {
	        Z->sample[s].group = 1;
        }
    }

    // Destroying dynamic sets
    for (int i = 0; i < n_clusters; i++) {
        iftDestroyDynamicSet(&S[i]);
    }

    // Saving centroids to iftCentroids
    if (centroids_hist != NULL) {
        for (int i = 0; i < centroids_hist->n_clusters; i++) {
            iftCentroid *aux = centroids_hist->centroid[i];
            aux->feat = iftAlloc(centroids_hist->n_feats, sizeof(*aux->feat));
            iftCopyFloatArray(aux->feat, curr_centroids[i], centroids_hist->n_feats);
        }
        centroids_hist->n_iter = iter;
    }
}


void iftSetDataSetPrototypeStatus(iftDataSet *Z, iftDynamicSet **S) {
    float min_distance, distance;
    int min_sample, s;
    int n_clusters = Z->ngroups;
    int n_features = Z->nfeats;
    iftSet *aux = NULL;
    float *mean = iftAlloc(n_features, sizeof(*mean));


    for (int s=0; s < Z->nsamples; s++)
      iftRemoveSampleStatus(&Z->sample[s], IFT_PROTOTYPE);
      
    
    for (int i = 0; i < n_clusters; i++) {
        min_distance = IFT_INFINITY_FLT;
        min_sample   = IFT_NIL;
        aux          = S[i]->begin;
        while (aux != NULL) {
            s = aux->elem;
            iftCopyDblArrayToFloatArray(mean, S[i]->mean, n_features);
            distance = iftEuclideanDistance(mean, Z->sample[s].feat, n_features);
            if (distance < min_distance) {
                min_distance = distance;
                min_sample = s;
            }
            aux = aux->next;
        }
        
        if (min_sample != IFT_NIL) {
            iftAddSampleStatus(&Z->sample[min_sample], IFT_PROTOTYPE);
        }
    }  
}

void iftInsertCentroidHistory(iftCentroid *centroid, float *centroid_feat, int n_feats) {
    centroid->hist_feat = iftRealloc(centroid->hist_feat, (centroid->n_hist + 1) * sizeof(*centroid->hist_feat));
    centroid->hist_feat[centroid->n_hist] = iftAlloc(n_feats, sizeof(*centroid->hist_feat[centroid->n_hist]));
    iftCopyFloatArray(centroid->hist_feat[centroid->n_hist], centroid_feat, n_feats);
    centroid->n_hist += 1;
}

int *iftGetDataSetPrototypes(iftDataSet *Z) {
    if (Z == NULL) {
        iftError("Dataset is NULL", "iftKMeans");
    }

    if (Z->ngroups <= 0) {
        iftError("Number of groups equal to 0", "iftKMeans");
    }

    int *prototypes = iftAllocIntArray(Z->ngroups);

    for (int s = 0; s < Z->nsamples; s++) {
        if (iftHasStatus(IFT_PROTOTYPE, Z->sample[s].status)) {
            prototypes[Z->sample[s].group - 1] = s;
        }
    }

    return prototypes; 
}

void iftInsertSampleDynamicSet(iftDataSet *Z, iftDynamicSet *S, int s) {
    if (S->size) {
        iftSet *a = (iftSet*) iftAlloc(1, sizeof *a);
        if (!a)
            iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertDynamicSet");
        a->elem = s;
        S->end->next = a;
        S->end = a;
    } else {
        S->begin = (iftSet*) iftAlloc(1, sizeof *S->begin);
        S->begin->elem = s;
        S->end = S->begin;
    }
    
    S->size += 1;

    for (int i = 0; i < S->dim; i++) {
        S->mean[i] += (Z->sample[s].feat[i] - S->mean[i]) / S->size;
    }   
}

void iftKMedoids(iftDataSet *Z, int n_clusters, int max_iter, iftKmeansDistance distance) {

    if (distance == NULL)
        distance = iftEuclideanDistance;

    if (Z == NULL) {
        iftError("Dataset is NULL", "iftKMedoids");
    }

    if (Z->nsamples == 0) {
        iftError("Dataset is empty", "iftKMedoids");
    }

    if (n_clusters > Z->nsamples) {
        iftError("Number of clusters is greater than number of samples", "iftKMedoids");
    }

    Z->ngroups = n_clusters;

    iftMatrix *distance_matrix = iftCreateMatrix(Z->nsamples, Z->nsamples);

    int i, j, s, t;

    // calculate the distance between every pair of all objects
    for (int row = 0; row < Z->nsamples; row++) {
        for (int col = 0; col < Z->nsamples; col++) {
            i = iftGetMatrixIndex(distance_matrix, col, row);
            if (row == col) {
                distance_matrix->val[i] = .0;
            } else if (col < row) {
                j = iftGetMatrixIndex(distance_matrix, row, col);
                distance_matrix->val[i] = distance_matrix->val[j];
            } else {
                distance_matrix->val[i] = distance(Z->sample[row].feat, Z->sample[col].feat, Z->nfeats);
            }
        }
    }

    iftMatrix *sumrow = iftMatrixSumRow(distance_matrix);

    float *v = iftAllocFloatArray(Z->nsamples);
    int *index = iftAllocIntArray(Z->nsamples);

    for (int obj = 0; obj < Z->nsamples; obj++) {
        index[obj] = obj;
        v[obj] = .0;
        for (int row = 0; row < Z->nsamples; row++) {
            i = iftGetMatrixIndex(distance_matrix, obj, row);
            j = iftGetMatrixIndex(sumrow, 0, row);
            v[obj] += distance_matrix->val[i] / sumrow->val[j];
        }
    }

    iftFQuickSort(v, index, 0, Z->nsamples - 1, IFT_INCREASING);

    // select initial medoids
    int *medoids = iftAllocIntArray(n_clusters);
    
    for (i = 0; i < n_clusters; i++) {
        medoids[i] = index[i];
    }

    int min_index, nearest_medoid, iter = 0;
    float min_distance;
    float sum_distances = IFT_INFINITY_FLT;
    float prev_sum_distances = IFT_INFINITY_FLT_NEG;
    float *cluster_sum_distances = iftAllocFloatArray(n_clusters);

    iftSet **cluster_set = iftAlloc(n_clusters, sizeof(*cluster_set));

    while (iter < max_iter && sum_distances != prev_sum_distances) {

        // initializing sum of distances from all objects to their medoids
        for (i = 0; i < n_clusters; i++) {
            cluster_sum_distances[i] = .0;
            cluster_set[i] = NULL;
        }

        // assigning each object to the nearest medoid
        for (s = 0; s < Z->nsamples; s++) {
            min_distance = IFT_INFINITY_FLT;
            nearest_medoid = s;
            for (int med = 0; med < n_clusters; med++) {
                i = iftGetMatrixIndex(distance_matrix, medoids[med], s);
                if (distance_matrix->val[i] < min_distance) {
                    min_distance = distance_matrix->val[i];
                    nearest_medoid = med;        
                }
            }
            Z->sample[s].group = nearest_medoid + 1;
            cluster_sum_distances[nearest_medoid] += min_distance;
            iftInsertSet(&cluster_set[nearest_medoid], s);
        }

        // update medoids
        prev_sum_distances = sum_distances;
        sum_distances = .0;
        for (int c = 0; c < n_clusters; c++) {
            iftIntArray *cluster_array = iftSetToArray(cluster_set[c]);
            float *cluster_sum = iftAllocFloatArray(cluster_array->n);
            for (s = 0; s < cluster_array->n; s++) {
                cluster_sum[s] = .0;
                for (t = 0; t < cluster_array->n; t++) {
                    i = iftGetMatrixIndex(distance_matrix, cluster_array->val[t], cluster_array->val[s]);
                    cluster_sum[s] += distance_matrix->val[i];
                }
            }
            min_index = iftMinIndexFloatArray(cluster_sum, cluster_array->n);
            sum_distances += cluster_sum[min_index];
            medoids[c] = cluster_array->val[min_index];
        }

        iter += 1;
    }

    for (s = 0; s < Z->nsamples; s++)
	    iftRemoveSampleStatus(&Z->sample[s], IFT_PROTOTYPE);

    for (int c = 0; c < n_clusters; c++) {
        s = medoids[c];
        iftAddSampleStatus(&Z->sample[s], IFT_PROTOTYPE);
    }

    iftDestroyMatrix(&distance_matrix);
    iftDestroyMatrix(&sumrow);
   
}
