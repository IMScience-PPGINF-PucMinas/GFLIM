#include "iftKmeans.h"

#include "ift/core/io/Stream.h"


// Adapted from - http://www.medphysics.wisc.edu/~ethan/kmeans/
// kmeans.c
// Ethan Brodsky - October 2011

#define _SILENCE

/*--------------------------- Private methods ---------------------------*/
float iftDistanceKmeans(float *f1, float *f2, int dim){
  float dist=0.0f;

  for (int i=0; i < dim; i++)
    dist += (f1[i]-f2[i])*(f1[i]-f2[i]);

  return dist;
}

float iftDistanceSphericalKmeans (float *f1, float *f2, int dim) {
	float inner_prod = 0.0;

	for (int i = 0; i < dim; i++)
		inner_prod += f1[i] * f2[i];

	return inner_prod;
}

void iftCalcAllDistances(iftDataSet* Z,iftDataSet* Zk, float* distance_output) {

  int dim = Z->nfeats;
  int n = Z->nsamples;
  int k = Zk->nsamples;

#pragma omp parallel for
  for (int ii = 0; ii < n; ii++) // for each point
    for (int kk = 0; kk < k; kk++) // for each cluster
      // calculate distance between point and cluster centroid
      distance_output[ii*k + kk] = iftDistanceKmeans( (Z->sample[ii].feat), (Zk->sample[kk].feat), dim);
}
  
float iftCalcTotalDistance(iftDataSet* Z,iftDataSet* Zk,int* label) {

  int dim = Z->nfeats;
  int n   = Z->nsamples;

  // NOTE: a point with cluster assignment -1 is ignored
  float tot_D = 0.0;
    
  // for every point,
  for (int ii = 0; ii < n; ii++) {
    // which cluster is it in?
    int active_cluster = label[ii];
        
    // sum distance
    if (active_cluster != -1)
      tot_D += iftDistanceKmeans( (Z->sample[ii].feat), (Zk->sample[active_cluster].feat), dim);
  }
   
  return tot_D;
}

void iftChooseAllClustersFromDistances(iftDataSet* Z,iftDataSet* Zk,float* distance_array,int* label) {

  int n = Z->nsamples;
  int k = Zk->nsamples;

  // for each point
#pragma omp parallel for
  for (int ii = 0; ii < n; ii++) {
    int best_index = 0;
    float closest_distance = distance_array[ii*k + 0];
  
    // for each cluster - initiallized on "0"
    for (int kk = 1; kk < k; kk++) {
      // distance between point and cluster centroid
      if (distance_array[ii*k + kk] <= closest_distance) {
        best_index = kk;
        closest_distance = distance_array[ii*k + kk];
      }
    }

    // record in array
    label[ii] = best_index;
  }
}

int iftCalcClusterCentroids(int bKmedoids,iftDataSet* Z,iftDataSet* Zk,int* ClAss) {
  int dim = Z->nfeats;
  int n = Z->nsamples;
  int k = Zk->nsamples;

  int* cluster_member_count = iftAllocIntArray(k);
  
  // initialize cluster centroid coordinate sums to zero
  for (int kk = 0; kk < k; kk++) {
    cluster_member_count[kk] = 0;
        
    for (int ff = 0; ff < dim; ff++)
      Zk->sample[kk].feat[ff] = 0.;
  }

  for (int ii = 0; ii < n; ii++) {
    // which cluster is it in?
    int active_cluster = ClAss[ii];

    // update count of members in that cluster
    cluster_member_count[active_cluster]++;
  }

  // verifying for empty clusters
  for (int kk = 0; kk < k; kk++) {
    if (cluster_member_count[kk] == 0) {
#ifndef _SILENCE
      fprintf(stderr,"empty cluster!\n");
#endif
      iftFree(cluster_member_count);
      return 0;
    }
  }

  // sum all points
  // for every point
  for (int ii = 0; ii < n; ii++) {
    // which cluster is it in?
    int active_cluster = ClAss[ii];

    // sum point coordinates for finding centroid
    for (int ff = 0; ff < dim; ff++)
      Zk->sample[active_cluster].feat[ff] += Z->sample[ii].feat[ff];
  }

  // now divide each coordinate sum by number of members to find mean/centroid
  // for each cluster
  for (int kk = 0; kk < k; kk++) {
    // for each dimension
    for (int ff = 0; ff < dim; ff++)
      Zk->sample[kk].feat[ff] /= cluster_member_count[kk]; 
    /// XXXX will divide by zero here for any empty clusters!
  }

  iftFree(cluster_member_count);

  if (bKmedoids) {
    int*   closest_member   = iftAllocIntArray(k);
    float* closest_distance = iftAllocFloatArray(k);

    for (int kk = 0; kk < k; kk++) {
      closest_member[kk]   =-1;
      closest_distance[kk] = IFT_INFINITY_FLT;
    }

    // Kmedoids - seraching for the closest point to centroids
    for (int ii = 0; ii < n; ii++) {
      // which cluster is it in?
      int active_cluster = ClAss[ii];
    
      float cur_distance =  iftDistanceKmeans( (Z->sample[ii].feat), (Zk->sample[active_cluster].feat), dim);    
      if ( cur_distance < closest_distance[active_cluster] ) {
	closest_distance[active_cluster] = cur_distance;
	closest_member[active_cluster]   = ii;
      }
    }

    // defining new centroid based on the closest point
    for (int kk = 0; kk < k; kk++) {
      // for each dimension
      for (int ff = 0; ff < dim; ff++)
	Zk->sample[kk].feat[ff] = Z->sample[closest_member[kk]].feat[ff]; 
      /// XXXX will divide by zero here for any empty clusters!
    }
  
    iftFree(closest_member);
    iftFree(closest_distance);
  }

  return 1; // OK!
}


/*
 * Strategy to avoid empty cluster - founded in http://www.academypublisher.com/ijrte/vol01/no01/ijrte0101220226.pdf
 */
int iftCalcClusterCentroids2(int bKmedoids, iftDataSet* Z, iftDataSet* Zk, int* ClAss) {
	printf("------- iftCalcClusterCentroids2\n");
	int dim = Z->nfeats;
	int n = Z->nsamples;
	int k = Zk->nsamples;

	int* cluster_member_count = iftAllocIntArray(k);

	for (int ii = 0; ii < n; ii++) {
		// which cluster is it in?
		int active_cluster = ClAss[ii];

		// update count of members in that cluster
		cluster_member_count[active_cluster]++;
	}

	// verifying for empty clusters
	for (int kk = 0; kk < k; kk++) {
		if (cluster_member_count[kk] == 0) {
#ifndef _SILENCE
			fprintf(stderr,"empty cluster!\n");
#endif
			iftFree(cluster_member_count);
			return 0;
		}
	}

	// sum all points
	// for every point
	for (int ii = 0; ii < n; ii++) {
		// which cluster is it in?
		int active_cluster = ClAss[ii];

		// sum point coordinates for finding centroid
		for (int ff = 0; ff < dim; ff++)
			Zk->sample[active_cluster].feat[ff] += Z->sample[ii].feat[ff];
	}

	// now divide each coordinate sum by number of members to find mean/centroid
	// for each cluster
	for (int kk = 0; kk < k; kk++) {
		// for each dimension
		for (int ff = 0; ff < dim; ff++)
			Zk->sample[kk].feat[ff] /= (cluster_member_count[kk]+1);
		/// XXXX will divide by zero here for any empty clusters!
	}

	iftFree(cluster_member_count);

	if (bKmedoids) {
		int* closest_member = iftAllocIntArray(k);
		float* closest_distance = iftAllocFloatArray(k);

		for (int kk = 0; kk < k; kk++) {
			closest_member[kk] = -1;
			closest_distance[kk] = IFT_INFINITY_FLT;
		}

		// Kmedoids - seraching for the closest point to centroids
		for (int ii = 0; ii < n; ii++) {
			// which cluster is it in?
			int active_cluster = ClAss[ii];

			float cur_distance = iftDistanceKmeans((Z->sample[ii].feat),
					(Zk->sample[active_cluster].feat), dim);
			if (cur_distance < closest_distance[active_cluster]) {
				closest_distance[active_cluster] = cur_distance;
				closest_member[active_cluster] = ii;
			}
		}

		// defining new centroid based on the closest point
		for (int kk = 0; kk < k; kk++) {
			// for each dimension
			for (int ff = 0; ff < dim; ff++)
				Zk->sample[kk].feat[ff] =
						Z->sample[closest_member[kk]].feat[ff];
			/// XXXX will divide by zero here for any empty clusters!
		}

		iftFree(closest_member);
		iftFree(closest_distance);
	}

	return 1; // OK!
}

void iftGetClusterMemberCount(iftDataSet* Z,iftDataSet* Zk,int* ClAss,int *cluster_member_count) {
  int n = Z->nsamples;
  int k = Zk->nsamples;

  // initialize cluster member counts
  for (int kk = 0; kk < k; kk++) 
    cluster_member_count[kk] = 0;
  
  // count members of each cluster    
  for (int ii = 0; ii < n; ii++)
    cluster_member_count[ClAss[ii]]++;
}

void iftClusterDiag(iftDataSet* Z,iftDataSet* Zk,int* ClAss)
{
  int dim = Z->nfeats;
  int k = Zk->nsamples;

  int* cluster_member_count = iftAllocIntArray(k);

  iftGetClusterMemberCount(Z,Zk,ClAss,cluster_member_count);
     
  printf("  Final clusters \n");
  for (int kk = 0; kk < k; kk++) { 

    float norm = 0.; 
    for (int ff=0;ff < dim ; ff++)
      norm += Zk->sample[kk].feat[ff] * Zk->sample[kk].feat[ff];
    norm = sqrt(norm);

    printf("    cluster %d:     members: %8d, centroid (%.1f %.1f), norm: %f \n", kk, cluster_member_count[kk], Zk->sample[kk].feat[0], Zk->sample[kk].feat[1],norm);
  }

  iftFree(cluster_member_count);
}

void copy_assignment_array(int n, int *src, int *tgt) {
  for (int ii = 0; ii < n; ii++)
    tgt[ii] = src[ii];
}
  
int assignment_change_count(int n, int a[], int b[]) {
  int change_count = 0;

  for (int ii = 0; ii < n; ii++)
    if (a[ii] != b[ii])
      change_count++;
        
  return change_count;
}

/*--------------------------- Public methods ---------------------------*/
void iftClusterDataSetByKMeans(iftDataSet *Z, int nb_cluster, int maxIterations, float minImprovement, char kmedoids,
                               char grid_sampling, char debug) {

  if (nb_cluster > Z->nsamples)
      iftError("The desired number of clusters can't be greather than the number of samples of the dataset",
               "iftClusterDataSetByKMeans");

  /* check if the nb of desired clusters is the same as the nb of samples*/
  if (nb_cluster == Z->nsamples){
    int l=0;
    for (int i=0;i<Z->nsamples;i++)
      Z->sample[i].group=++l;
    Z->ngroups=nb_cluster;
    return;
  }

  int nsamples = Z->nsamples;
  int nfeats = Z->nfeats;
  /* initialize the centers randomly from the dataset*/
  iftSample **centers;

  if (grid_sampling && Z->ref_data_type == IFT_REF_DATA_MIMAGE){
    iftMImage *mimg=(iftMImage *)Z->ref_data;
    iftImage *mask1 = iftSelectImageDomain(mimg->xsize,mimg->ysize,mimg->zsize);
    iftImage *mask_sampling = iftGridSampling(mimg, mask1,nb_cluster);

    nb_cluster=0;
    for (int i=0;i<mask_sampling->n;i++){
      if (mask_sampling->val[i]){
        nb_cluster++;
      }
    }

    centers=(iftSample **)iftAlloc(nb_cluster, sizeof(iftSample *));

    for (int kk=0;kk<nb_cluster;kk++){
      centers[kk]=(iftSample*) iftAlloc(1, sizeof(iftSample));
    }

    int t=0;
    for (int s=0;s<Z->nsamples;s++){
      int voxel=Z->sample[s].id;
      if (mask_sampling->val[voxel]){
        iftCopySample(&Z->sample[s],centers[t],nfeats,false);
        centers[t]->id = s;
        t++;
      }
    }
    iftDestroyImage(&mask1);
    iftDestroyImage(&mask_sampling);
  }
  else{
    centers=(iftSample **)iftAlloc(nb_cluster, sizeof(iftSample *));

    for (int kk=0;kk<nb_cluster;kk++){
      centers[kk]=(iftSample*) iftAlloc(1, sizeof(iftSample));
    }

    // randomly select K samples
    int *sample = iftAllocIntArray(nsamples);
    int *count  = iftAllocIntArray(nsamples);
    for (int s=0; s < nsamples; s++) {
      sample[s]=s;
      count[s]=5;
    }
    int t = 0;int high = nsamples-1;
    int i,s;
    while (t < nb_cluster) {
      i = iftRandomInteger(0,high);
      if (count[i]==0){
        s=sample[i];
        iftCopySample(&Z->sample[s],centers[t],nfeats,false);
        iftSwap(sample[i], sample[high]);
        iftSwap(count [i], count [high]);
        centers[t]->id = s;
        t++; high--;
      }
      count[i]--;
    }
    iftFree(count);
    iftFree(sample);
  }

  /*begin to execute the Algorithm*/
  /*---------------------------------------------------------------*/

  int *clusterCount = iftAllocIntArray(nb_cluster);

  float *sample_near_dist_to_cluster  = iftAllocFloatArray(nsamples);
  int   *current_labels = iftAllocIntArray(nsamples);
  int   *temp_labels = iftAllocIntArray(nsamples);
  float **temp_feats=(float **)iftAlloc(nb_cluster, sizeof(float *));
  for (int i=0; i < nb_cluster; i++)
    temp_feats[i]=iftAllocFloatArray(nfeats);

  float curr_error, prevError = FLT_MAX;
  /* vectors to compute the new centroids*/
  int*   closest_member=NULL;
  float* closest_distance=NULL;

  closest_member   = iftAllocIntArray(nb_cluster);
  closest_distance = iftAllocFloatArray(nb_cluster);

  int s;
  for (int it = 0; it < maxIterations; ++it) {

    if (debug){
      printf("iter nb: %d\n",it);
    }

    for (int s=0;s<nsamples;s++){
      sample_near_dist_to_cluster[s]=IFT_INFINITY_FLT;
    }

    /* compute the distances from each sample to centers*/
#pragma omp parallel for
    for (int ii = 0; ii < nsamples; ii++) // for each point
      for (int kk = 0; kk < nb_cluster; kk++)
      {
        float dist = iftDistanceKmeans( (Z->sample[ii].feat), (centers[kk]->feat), nfeats);
        if (dist < sample_near_dist_to_cluster[ii]){
          sample_near_dist_to_cluster[ii]=dist;
          temp_labels[ii]=kk;

        }
      }

    /* compute the total error*/
    curr_error = 0.0;
    for (int ii = 0; ii < nsamples; ii++) {
      curr_error += sample_near_dist_to_cluster[ii];
    }

    if((fabs(curr_error - prevError) < minImprovement))
      break;

    iftCopyIntArray(current_labels, temp_labels, nsamples);
    prevError = curr_error;

    //reset cluster values
    for (int kk = 0; kk < nb_cluster; ++kk) {
      clusterCount[kk] = 0;
      for (int j=0;j<nfeats;j++)
        temp_feats[kk][j]=0.0;
    }

    /* update the new cluster center values*/
    for (int ii = 0; ii < nsamples; ++ii) {
      int c = current_labels[ii];
      clusterCount[c]++;

      for (int j = 0; j < nfeats; ++j) {
        temp_feats[c][j] += Z->sample[ii].feat[j];
      }
    }

    /*if there is an empty cluster we remove it. We can think of another strategy here*/
    bool flag=true;
    bool decrement_cluster;
    while (flag){
      decrement_cluster=false;
      for (int c=0;c<nb_cluster;c++){
        if (clusterCount[c]==0){
          clusterCount[c]=clusterCount[nb_cluster-1];
          clusterCount[nb_cluster-1]=0;

          for (int j=0;j<nfeats;j++){
            temp_feats[c][j]=temp_feats[nb_cluster-1][j];
          }

          for (int ii=0;ii<nsamples;ii++){
            if (current_labels[ii] == nb_cluster-1)
              current_labels[ii]=c;
          }
          decrement_cluster=true;
          break;
        }
      }
      if (decrement_cluster)
        nb_cluster--;
      else
        flag=false;
    }

    if (debug){
      for (int c=0;c<nb_cluster;c++)
        if (clusterCount[c] == 0)
          iftError("Cluster %d emtpy","iftClusterDataSetByKMeans",c);
    }



    /* divide the feature values by the cluster count*/
#pragma omp parallel for
    for (int kk = 0; kk < nb_cluster; ++kk) {
      for (int j=0;j<nfeats;j++){
        temp_feats[kk][j]=temp_feats[kk][j]/clusterCount[kk];
      }
    }

    if (kmedoids){
      /* for each cluster find the new centroid, it would be the nearest sample to the fiction features computed for the center*/
      for (int kk = 0; kk < nb_cluster; kk++) {
        closest_member[kk]   =-1;
        closest_distance[kk] = IFT_INFINITY_FLT;
      }
      for (int ii = 0; ii < nsamples; ii++) {
        // which cluster is it in?
        int active_cluster = current_labels[ii];
        float cur_distance = iftDistanceKmeans( (Z->sample[ii].feat), temp_feats[active_cluster], nfeats);
        if ( cur_distance < closest_distance[active_cluster] ) {
          closest_distance[active_cluster] = cur_distance;
          closest_member[active_cluster]   = ii;
        }
      }

      for (int kk = 0; kk < nb_cluster; kk++) {
        s=closest_member[kk];
        iftCopySample(&Z->sample[s],centers[kk],nfeats,false);
        centers[kk]->id = s;
      }
    }
    else{
      /*the cluster centers aren't real samples, so we only update the feature vectors*/
      for (int kk=0;kk<nb_cluster;kk++)
        centers[kk]->feat=temp_feats[kk];
    }
  }

  /* mark the corresponding prototypes in the dataset for a possible posterior propagation of labels, if the method was kmedoid they are already selected, otherwise find the sample nearby each cluster center*/
  if (!kmedoids){
    /* we need to find the real representative of each cluster */
    for (int kk = 0; kk < nb_cluster; kk++) {
      closest_member[kk]   =-1;
      closest_distance[kk] = IFT_INFINITY_FLT;
    }
    for (int ii = 0; ii < nsamples; ii++) {
      // which cluster is it in?
      int active_cluster = current_labels[ii];
      float cur_distance = iftDistanceKmeans( (Z->sample[ii].feat), temp_feats[active_cluster], nfeats);
      if ( cur_distance < closest_distance[active_cluster] ) {
        closest_distance[active_cluster] = cur_distance;
        closest_member[active_cluster]   = ii;
      }
    }

    for (int kk = 0; kk < nb_cluster; kk++) {
      s=closest_member[kk];
      iftCopySample(&Z->sample[s],centers[kk],nfeats,false);
      centers[kk]->id = s;
    }
  }

  for (int kk=0;kk<nb_cluster;kk++) {
    iftAddSampleStatus(&Z->sample[centers[kk]->id], IFT_PROTOTYPE);
    iftAddSampleStatus(&Z->sample[centers[kk]->id], IFT_CENTROID);
  }

  /* just assign the labels*/
  for (int ii = 0; ii < Z->nsamples; ++ii) {
    Z->sample[ii].group = current_labels[ii]+1;
  }
  Z->ngroups=nb_cluster;

  iftFree(closest_member);
  iftFree(closest_distance);
  for (int kk=0;kk<nb_cluster;kk++)
    iftFree(centers[kk]);
  iftFree(centers);
  iftFree(clusterCount);
  iftFree(temp_labels);
  iftFree(current_labels);
  iftFree(sample_near_dist_to_cluster);
  for (int i=0; i < nb_cluster; i++)
    iftFree(temp_feats[i]);
  iftFree(temp_feats);
}

void iftSimpleKmeansRun(iftDataSet* Z, iftDataSet** pZk, int maxIterations,float minImprovement) {

  int n = Z->nsamples;
  int k = (*pZk)->nsamples;
  int nfeats = Z->nfeats;

  int *clusterCount = iftAllocIntArray(k);
  float *dist      = iftAllocFloatArray(n*k);
  int   *classCur = iftAllocIntArray(n);
  int   *classPrev = iftAllocIntArray(n);

  float error, prevError = FLT_MAX;

  int it;
  for (it = 0; ; ++it) {

    iftCalcAllDistances(Z, *pZk, dist);
    iftChooseAllClustersFromDistances(Z, *pZk, dist, classCur);
    error = iftCalcTotalDistance(Z, *pZk, classCur);

    error/=n;

    if((fabs(error - prevError) < minImprovement) || (it >= maxIterations))
      break;

    iftCopyIntArray(classPrev, classCur, n);
    prevError = error;

    //update clusters ignoring empty clusters
    for (int i = 0; i < k; ++i) {
      clusterCount[i] = 0;

      for (int j = 0; j < nfeats; ++j) {
        (*pZk)->sample[i].feat[j] = 0.0;
      }
    }

    for (int i = 0; i < n; ++i) {
      clusterCount[classCur[i]]++;
    }

    for (int i = 0; i < n; ++i) {
      int c = classCur[i];
      for (int j = 0; j < nfeats; ++j) {
        (*pZk)->sample[c].feat[j] += Z->sample[i].feat[j] / clusterCount[c];
      }
    }
  }

    for (int i = 0; i < Z->nsamples; ++i) {
      Z->sample[i].label = classCur[i];
      Z->sample[i].group = classCur[i];
    }
  Z->ngroups = k;

  iftFree(clusterCount);
  iftFree(classPrev);
  iftFree(classCur);
  iftFree(dist);
}

void iftKmeansRun(int bKmedoids,iftDataSet* Z,iftDataSet** pZk,int maxIterations,float minImprovement) {

  int n = Z->nsamples;
  int k = (*pZk)->nsamples;

  if (n == k) {
    for (int ii=0; ii < n ; ii++)
      for(int ff=0; ff < Z->nfeats; ff++ )
	(*pZk)->sample[ii].feat[ff] = Z->sample[ii].feat[ff];
    return;

  } else if (n < k) {
    char msg[300];
    sprintf(msg,"Number of samples (%d) smaller than the number of required centroids (%d)",n,k);
      iftError(msg, "iftKmeansRun");
  }

  float *dist      = iftAllocFloatArray(n*k);
  int   *ClAssCur  = iftAllocIntArray(n);
  int   *ClAssPrev = iftAllocIntArray(n);

  if (!dist || !ClAssCur || !ClAssPrev )
      iftError("Error allocating dist arrays", "iftKmeansRun");
    
  // initial setup  
  iftCalcAllDistances(Z,*pZk,dist);
  iftChooseAllClustersFromDistances(Z,*pZk,dist,ClAssCur);
  copy_assignment_array(n, ClAssCur,ClAssPrev);

  // BATCH UPDATE
  float prev_totD = IFT_INFINITY_FLT;
  int   batch_iteration = 0;
  while (batch_iteration < maxIterations) {
#ifndef _SILENCE
    fprintf(stdout,"batch iteration %d \n", batch_iteration);fflush(stdout);
    //iftClusterDiag(Z,*pZk,ClAssCur);
#endif

    // update cluster centroids - no empty cluster
    int it=0;
    while (!iftCalcClusterCentroids(bKmedoids,Z,*pZk,ClAssCur))  {
      // deal with empty clusters
      if (it++ == maxIterations) {
        // restarting
        prev_totD = IFT_INFINITY_FLT;
        batch_iteration = 0;
        it = 0;
        //
        k = (int)(0.9*k+0.5);
        fprintf(stdout,"reducing k: %d\n",k);fflush(stdout);
      }        
      iftDestroyDataSet(pZk);
      *pZk = iftKmeansInitCentroidsFromSamples(Z,k);
      copy_assignment_array(n, ClAssCur,ClAssPrev);
      iftCalcAllDistances(Z,*pZk,dist);
      iftChooseAllClustersFromDistances(Z,*pZk,dist,ClAssCur);
    }

    // see if we've failed to improve
    float totD = iftCalcTotalDistance(Z,*pZk,ClAssCur);
    if (prev_totD - totD < minImprovement) {
      // failed to improve - currently solution worse than previous
      // restore old assignments
      copy_assignment_array(n, ClAssPrev, ClAssCur);

      // recalc centroids
      iftCalcClusterCentroids(bKmedoids,Z,*pZk,ClAssCur);
      // printf("  negative progress made on this step - iteration completed (%.2f) \n", totD - prev_totD);
      // done with this phase
      break;
    }
           
    // save previous step
    copy_assignment_array(n, ClAssCur, ClAssPrev);
         
    // move all points to nearest cluster
    iftCalcAllDistances(Z,*pZk,dist);
    iftChooseAllClustersFromDistances(Z,*pZk,dist,ClAssCur);
         
    int change_count = assignment_change_count(n, ClAssCur, ClAssPrev);
       
    //printf("%3d   %u   %9d  %16.2f %17.2f\n", batch_iteration, 1, change_count, totD, totD - prev_totD);
    //fflush(stdout);
         
    // done with this phase if nothing has changed
    if (change_count == 0) {
      // printf("  no change made on this step - iteration completed \n");
      break;
    }

    printf("Iteration: %d with error: %f\n",batch_iteration,totD);
    prev_totD = totD;
    batch_iteration++;
  }

  //iftClusterDiag(Z,*pZk,ClAssCur);

  // write to output array
  for(int ii=0; ii < Z->nsamples ; ii++)
    Z->sample[ii].label = ClAssCur[ii]+1;
  Z->ngroups = k;

  iftFree(dist);
  iftFree(ClAssCur);
  iftFree(ClAssPrev);
}          


void iftSphericalKmeansRun(iftDataSet* Z, iftDataSet** pZk, int maxIterations) {
	int n = Z->nsamples;
	int k = (*pZk)->nsamples;
	int dim = (*pZk)->nfeats;

	if (n == k) {
		for (int ii = 0; ii < n; ii++)
			for (int ff = 0; ff < dim; ff++)
				(*pZk)->sample[ii].feat[ff] = Z->sample[ii].feat[ff];
		return;
	} else if (n < k) {
		char msg[300];
		sprintf(msg,
				"Number of samples (%d) smaller than the number of required centroids (%d)",
				n, k);
        iftError(msg, "iftKmeansRun");
	}

	iftDataSet* Zknew = iftCreateDataSet(k, dim);
	iftDataSet* Ztmp = NULL;

	float *innerprod = iftAllocFloatArray(k);
	int *cluster_member_count = iftAllocIntArray(k);
	int *cluster_assigned = iftAllocIntArray(n);

	if (!innerprod || !cluster_member_count || !cluster_assigned)
        iftError("Error allocating dist arrays", "iftKmeansRun");

	// normalizing to unit norm the centroids
	Ztmp = iftUnitNormDataSet(*pZk);
	iftDestroyDataSet(pZk);
	*pZk = Ztmp;

	int batch_iteration = 0;
	while (batch_iteration < maxIterations) {
#ifndef _SILENCE
		fprintf(stdout,"batch iteration %d \n", batch_iteration);fflush(stdout);
#endif

		for (int kk = 0; kk < k; kk++) {
			cluster_member_count[kk] = 0;
			for (int ff = 0; ff < dim; ff++)
				Zknew->sample[kk].feat[ff] = 0.;
		}

		// For each sample
		for (int s = 0; s < n; s++) {
			// Determining the kk centroid that achieves the maximum inner product with sample s
			int maxkk = 0;
			for (int kk = 0; kk < k; kk++) {
				innerprod[kk] = 0.;
				for (int ff = 0; ff < dim; ff++)
					innerprod[kk] += Z->sample[s].feat[ff]
							* (*pZk)->sample[kk].feat[ff];
				if ((kk != 0) && (fabsf(innerprod[maxkk]) < fabsf(innerprod[kk])))
					maxkk = kk;
			}
			cluster_assigned[s] = maxkk;
			cluster_member_count[maxkk]++;

			// updating new centroid
			for (int ff = 0; ff < dim; ff++)
				Zknew->sample[maxkk].feat[ff] += innerprod[maxkk]
						* Z->sample[s].feat[ff];
		}

		// Zk <- Zknew & reinitializing empty clusters
		for (int kk = 0; kk < k; kk++) {
			if (cluster_member_count[kk] == 0) {
#ifndef _SILENCE
				fprintf(stdout,"batch iteration %d/kernel %d\n", batch_iteration,kk);
#endif
				for (int ff = 0; ff < dim; ff++)
					(*pZk)->sample[kk].feat[ff] = iftRandomNormalFloat(0.0, 1.0);

			} else {
				for (int ff = 0; ff < dim; ff++)
					(*pZk)->sample[kk].feat[ff] = Zknew->sample[kk].feat[ff];
			}
		}

		// normalizing to unit norm the centroids
		Ztmp = iftUnitNormDataSet(*pZk);
		iftDestroyDataSet(pZk);
		*pZk = Ztmp;

		batch_iteration++;
	}

	// write to output array
	for (int ii = 0; ii < Z->nsamples; ii++)
		Z->sample[ii].label = cluster_assigned[ii] + 1;
	Z->ngroups = k;

	iftFree(innerprod);
	iftFree(cluster_member_count);
	iftFree(cluster_assigned);

	iftDestroyDataSet(&Zknew);
}          

iftDataSet* iftKmeansInitCentroidsFromSamples(iftDataSet* Z, int k) {

  int s, t, i, high, *sample=NULL, *count=NULL;

  int n = Z->nsamples;
  int dim = Z->nfeats;

  iftDataSet* Zk = iftCreateDataSet(k,dim);

  // randomly select K samples 
  sample = iftAllocIntArray(n); 
  count  = iftAllocIntArray(n); 
  for (s=0; s < n; s++) {
    sample[s]=s;
    count[s]=100;
  }
  t = 0;high = n-1;
  while ( t < k ) {
    i = iftRandomInteger(0,high);
    if (count[i]==0){
      s=sample[i];
      iftCopySample(&Z->sample[s],&Zk->sample[t],Z->nfeats,true);
      iftSwap(sample[i], sample[high]);
      iftSwap(count [i], count [high]);
      t++; high--;
    }else{
      count[i]--;
    }
  }

  iftFree(count);
  iftFree(sample);

  return Zk;
}

iftDataSet* iftKmeansInitCentroidsRandomNormal(iftDataSet* Z, int k) {
  int dim = Z->nfeats;

  iftDataSet* Zk = iftCreateDataSet(k,dim);

  // randomly selecting centroids from a normal distribution - one feature at a time
  for(int ff=0; ff < dim; ff++) { 
    for(int kk=0; kk < k; kk++) {
      Zk->sample[kk].feat[ff] = iftRandomNormalFloat(0.0, 1.0);
    }
  }

  // normalizing each centroid to unit length
  for(int kk=0; kk < k; kk++) {
    float norm=0;
    for(int ff=0; ff < dim; ff++) { 
      norm+= Zk->sample[kk].feat[ff] * Zk->sample[kk].feat[ff];
    }
    norm=sqrt(norm);
    if (norm > 1.) {
      for(int ff=0; ff < dim; ff++) { 
	Zk->sample[kk].feat[ff] /= norm;
      }
    }
  }

  return Zk;
}



/* Build a iftDataSet with the most representative samples (nearest samples from centroids)
 *
 * WARNING: Zorig must have just one class
 *
 * *** Inputs ***
 * Zorig = Dataset with the original samples
 * Zc    = Dataset with the clustered samples (ps: the samples can be reduced spatially
 * Zk	 = Dataset with the k centroids
 * m	 = number of representative samples to be chosen in each cluster
 *
 * if include_centroids, check if the sample feat space were Reduced
 */
iftDataSet *iftNearestSamplesFromCentroids(iftDataSet* Zorig, iftDataSet *Zc, iftDataSet* Zk, int m, int include_centroids) {
	char msg[256];

	if (Zorig->nclasses != 1)
        iftError("Invalid Number of Classes from the Original Dataset... Zorig->nclasses must be 1",
                 "iftNearestSamplesFromCentroids");

	if (Zorig->nsamples != Zc->nsamples) {
		sprintf(msg, "Incompatible Original and Clustered Dataset: num of samples is different (%d - %d)", Zorig->nsamples, Zc->nsamples);
        iftError(msg, "iftNearestSamplesFromCentroids");
	}

	if (Zc->nfeats != Zk->nfeats) {
		sprintf(msg, "Incompatible Clustered and Centroid Dataset: num of feats is different (%d - %d)", Zc->nfeats, Zk->nfeats);
        iftError(msg, "iftNearestSamplesFromCentroids");
	}

	int k 		 = Zk->nsamples; // number of clusters
	int nsamples = Zorig->nsamples;

	float **dist  = (float**) iftAlloc(k, sizeof(float*)); // matrix KxM that stores the distances of the M nearest samples from the clusteres
	int **indices = (int**) iftAlloc(k, sizeof(int*)); // matrix KxM that stores the index of the M nearest samples from the clusteres

	// Initialize the distance and index matrices
	for (int i = 0; i < k; i++) {
		dist[i]    = iftAllocFloatArray(m);
		indices[i] = iftAllocIntArray(m);

		for (int j = 0; j < m; j++) {
			dist[i][j] 	  = IFT_INFINITY_FLT;
			indices[i][j] = -1;
		}
	}

	// It obtains the nearest samples from each cluster
	for (int s = 0; s < nsamples; s++) {
		int cluster = Zc->sample[s].label - 1;
		float d = iftDistanceKmeans(Zc->sample[s].feat, Zk->sample[cluster].feat, Zc->nfeats);
		int idx = -1;


		for (int i = 0; i < m; i++)
			if (d < dist[cluster][i]) {
				idx = i;
				break;
			}

		if (idx >= 0) {
			for (int i = m-1; i > idx; i--) {
				dist[cluster][i] 	= dist[cluster][i-1];
				indices[cluster][i] = indices[cluster][i-1];
			}
			dist[cluster][idx] 	  = d;
			indices[cluster][idx] = s;
		}

//		printf("---- sample: %d, cluster: %d, d: %f\n", s, cluster, d);
//		printf("Dist - [%d]\n", s);
//		for (int ii = 0; ii < k; ii++) {
//			for (int jj = 0; jj < m; jj++)
//				printf("%f ", dist[ii][jj]);
//			puts("");
//		}
//		printf("\nMat - [%d]\n", s);
//		for (int ii = 0; ii < k; ii++) {
//			for (int jj = 0; jj < m; jj++)
//				printf("%d ", indices[ii][jj]);
//			puts("");
//		}
//		puts("");
	}

	// Check if there are m samples in each cluster
	int ns; // num of samples per cluster
	for (int i = 0; i < k; i++) {
		ns = 0;
		for (int j = 0; j < m; j++)
			if (indices[i][j] == -1) {
				char msg[256];
				sprintf(msg, "Cluster %d doesn't have m = %d samples --- %d/%d\n", i+1, m, ns, m);
                iftError(msg, "iftNearestSamplesFromCentroids");
			} else
				ns++;
	}

	iftDataSet *Zout = iftCreateDataSet(k*m, Zorig->nfeats);
	Zout->nclasses = 1;
	for (int i = 0; i < k; i++)
		for (int j = 0; j < m; j++) {
			int s = i*m + j;
			int idx = indices[i][j];

			Zout->sample[s].id 		= Zorig->sample[idx].id;
			Zout->sample[s].truelabel 	= Zorig->sample[idx].truelabel;
			Zout->sample[s].label 	= Zorig->sample[idx].label;
			iftSetSampleStatus(&Zout->sample[s], Zorig->sample[idx].status);
			Zout->sample[s].weight 	= Zorig->sample[idx].weight;

			for (int f = 0; f < Zorig->nfeats; f++) {
				Zout->sample[s].feat[f] = Zorig->sample[idx].feat[f];
			}
		}

	if (include_centroids) {
		puts("--> Including Centroids");
		iftDataSet *Zaux;
		if (iftIsTransformedDataSetByPCA(Zc)) { // REDUCED
			puts("--> Applying the PCA Inverse Transformation in Centroids\n");
			iftSetStatus(Zc, IFT_TRAIN);
			Zaux = iftInverseTransformTestDataSetByPCA(Zc, Zk);
		} else // NOT REDUCED
			Zaux = iftCopyDataSet(Zk, true);

		for (int i = 0; i < k; i++) {
			int s = i*m + (m-1);

			Zout->sample[s].id 		= Zaux->sample[i].id;
			Zout->sample[s].truelabel 	= Zaux->sample[i].truelabel;
			Zout->sample[s].label 	= Zaux->sample[i].label;
			iftSetSampleStatus(&Zout->sample[s], Zaux->sample[i].status);
			Zout->sample[s].weight 	= Zaux->sample[i].weight;

			for (int f = 0; f < Zorig->nfeats; f++)
				Zout->sample[s].feat[f] = Zaux->sample[i].feat[f];
		}

		iftDestroyDataSet(&Zaux);
	}


	return Zout;
}


/* Build a iftDataSet by following an alternating order of the nearest samples from centroids.
 * The sequence of the dataset building is: one nearest sample of each cluster per time
 *
 * WARNING: Zorig must have just one class
 *
 * *** Inputs ***
 * Zorig = Dataset with the original samples
 * Zc    = Dataset with the clustered samples (ps: the samples can be reduced spatially
 * Zk	 = Dataset with the k centroids
 * method = "kmeans" or "spherical"
 *
 * if include_centroids, check if the sample feat space were Reduced
 */
iftDataSet *iftBuildDataSetAsNearestSamplesFromCentroids(iftDataSet* Zorig, iftDataSet *Zc,
			iftDataSet* Zk, int include_centroids, char *method) {
	char msg[256];

	if (Zorig->nsamples != Zc->nsamples) {
		sprintf(msg, "Incompatible Original and Clustered Dataset: num of samples is different (%d - %d)", Zorig->nsamples, Zc->nsamples);
        iftError(msg, "iftBuildDataSetAsNearestSamplesFromCentroids");
	}

	if (Zc->nfeats != Zk->nfeats) {
		sprintf(msg, "Incompatible Clustered and Centroid Dataset: num of feats is different (%d - %d)", Zc->nfeats, Zk->nfeats);
        iftError(msg, "iftBuildDataSetAsNearestSamplesFromCentroids");
	}

	int k = Zk->nsamples; // number of clusters

	int *samples_per_cluster = iftAllocIntArray(k); // number of samples per cluster
	int *idx_cluster = iftAllocIntArray(k); // it stores the current idx of the distance and index matrices of each cluster

	for (int s = 0; s < Zc->nsamples; s++)
		samples_per_cluster[Zc->sample[s].label-1]++;

	for (int i = 0; i < k; i++)
		printf("Cluster: %d - %d samples\n", i+1, samples_per_cluster[i]);
	puts("");

	// matrix that stores the distance of each sample with its centroid
	float **dist  = (float**) iftAlloc(k, sizeof(float*));
	// matrix that stores the sample indices of each cluster
	int **indices = (int**) iftAlloc(k, sizeof(int*));
	for (int i = 0; i < k; i++) {
		dist[i]    = iftAllocFloatArray(samples_per_cluster[i]);
		indices[i] = iftAllocIntArray(samples_per_cluster[i]);
	}

	if (!strcmp(method, "kmeans")) {
		// Compute the distance to the centroid for every sample
		for (int s = 0; s < Zc->nsamples; s++) {
			int cluster = Zc->sample[s].label - 1;
			dist[cluster][idx_cluster[cluster]] = iftDistanceKmeans(Zc->sample[s].feat, Zk->sample[cluster].feat, Zc->nfeats);
			indices[cluster][idx_cluster[cluster]] = s;
			idx_cluster[cluster]++;
		}

		// Sort the distances of each cluster in increasing order
		for (int i = 0; i < k; i++) {
			idx_cluster[i] = 0;
			iftFQuickSort(dist[i], indices[i], 0, samples_per_cluster[i]-1, IFT_INCREASING);
		}
	}
	else if (!strcmp(method, "spherical")) {
		for (int s = 0; s < Zc->nsamples; s++) {
			int cluster = Zc->sample[s].label - 1;
			dist[cluster][idx_cluster[cluster]] = iftDistanceSphericalKmeans(Zc->sample[s].feat, Zk->sample[cluster].feat, Zc->nfeats);
			indices[cluster][idx_cluster[cluster]] = s;
			idx_cluster[cluster]++;
		}

		// Sort the distances of each cluster in decreasing order
		for (int i = 0; i < k; i++) {
			idx_cluster[i] = 0;
			iftFQuickSort(dist[i], indices[i], 0, samples_per_cluster[i]-1, IFT_DECREASING);
		}
	}

	// just to check if the distance matrix is in IFT_INCREASING ORDER
//	for (int i = 0; i < k; i++) {
//		printf("Cluster: %d\n", i);
//		for (int j = 0; j < samples_per_cluster[i]; j++)
//			printf("[%d] = %f\n", indices[i][j], dist[i][j]);
//		puts("\n");
//	}

	iftDataSet *Zout = NULL;
	int s = 0;
	if (include_centroids) {
		Zout = iftCreateDataSet(k + Zorig->nsamples, Zorig->nfeats);

		puts("--> Including Centroids");
		iftDataSet *Zaux;
		if (iftIsTransformedDataSetByPCA(Zc)) { // REDUCED
			puts("--> Applying the PCA Inverse Transformation in Centroids\n");
			iftSetStatus(Zc, IFT_TRAIN);
			Zaux = iftInverseTransformTestDataSetByPCA(Zc, Zk);
		} else // NOT REDUCED
			Zaux = iftCopyDataSet(Zk, true);

		for (int i = 0; i < k; i++) {
			Zout->sample[s].id 		= Zaux->sample[i].id;
			Zout->sample[s].truelabel 	= Zaux->sample[i].truelabel;
			Zout->sample[s].label 	= Zaux->sample[i].label;
			iftSetSampleStatus(&Zout->sample[s], Zaux->sample[i].status);
			Zout->sample[s].weight 	= Zaux->sample[i].weight;

			for (int f = 0; f < Zorig->nfeats; f++)
				Zout->sample[s].feat[f] = Zaux->sample[i].feat[f];
			s++;
		}

		iftDestroyDataSet(&Zaux);
	} else {
		Zout = iftCreateDataSet(Zorig->nsamples, Zorig->nfeats);
	}

	while (s < Zout->nsamples) {
		for (int cluster = 0; cluster < k; cluster++) {
			if (idx_cluster[cluster] < samples_per_cluster[cluster]) {
				int idx = indices[cluster][idx_cluster[cluster]++];

//				printf("[%d] = %d\n", s, idx);
				Zout->sample[s].id     = Zorig->sample[idx].id;
				Zout->sample[s].truelabel  = Zorig->sample[idx].truelabel;
				Zout->sample[s].label  = Zorig->sample[idx].label;
				iftSetSampleStatus(&Zout->sample[s], Zorig->sample[idx].status);
				Zout->sample[s].weight = Zorig->sample[idx].weight;

				for (int f = 0; f < Zorig->nfeats; f++)
					Zout->sample[s].feat[f] = Zorig->sample[idx].feat[f];
				s++;
			}
		}
	}

	for (int i = 0; i < k; i++) {
		iftFree(dist[i]);
		iftFree(indices[i]);
	}
	iftFree(dist);
	iftFree(indices);
	iftFree(idx_cluster);
	iftFree(samples_per_cluster);

	return Zout;
}
