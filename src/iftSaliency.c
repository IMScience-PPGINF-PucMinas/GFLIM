#include "iftSaliency.h"
#include "iftSaliencyPriors.h"
#include "ift/core/dtypes/SaliencyArray.h"
#include "ift/core/dtypes/IntQueue.h"

iftITSELFParameters *iftInitializeITSELFParametersByDefault(){
    iftITSELFParameters *params = iftCreateITSELFParameters();
    params->itself_iterations = 12;
    params->number_superpixels = 2500;
    params->obj_seeds_number = 30;
    params->oisf_gamma = (float)10.0;
    params->oisf_iterations = 3;
    params->oisf_beta = 12;
    params->superpixel_increase_ratio = (float)0.8;
    params->oisf_alpha = (float)0.8;
    params->query_importance = (float)1;
    params->normalization_value = (float)0.01;
    params->integration_lambda = (float)0.0001;
    params->integration_iteration = 1;
    params->enhancement_type = ITSELF_BY_SALIENCY;

    params->prior_params.center_variance = (float)0.0;
    params->prior_params.color_variance = (float)0.0;
    params->prior_params.ellipse_variance = (float)0.02;
    params->prior_params.focus_variance = (float)0.2;
    params->prior_params.red_yellow_variance = (float)0.0;
    params->prior_params.white_variance = (float)0.0;
    params->prior_params.foreground_distance_variance = (float)0.02;
    params->prior_params.foreground_feature_variance = (float)0.05;

    return params;
}

iftITSELFParameters *iftCreateITSELFParameters(){
    iftITSELFParameters *darr = (iftITSELFParameters *) iftAlloc(1, sizeof(iftITSELFParameters));
    return darr;
}


void iftDestroyITSELFParameters(iftITSELFParameters **darr){
    if (darr != NULL && *darr != NULL) {
        iftITSELFParameters *darr_aux = *darr;

        iftFree(darr_aux);
        *darr = NULL;
    }
}

iftITSELFParameters *iftInitializeITSELFParametersByDefaultU2Net(){
    iftITSELFParameters *params = iftCreateITSELFParameters();
    params->itself_iterations = 12;
    params->number_superpixels = 2500;
    params->obj_seeds_number = 10;
    params->oisf_gamma = (float)10.0;
    params->oisf_iterations = 5;
    params->oisf_beta = 12;
    params->superpixel_increase_ratio = (float)0.8;
    params->oisf_alpha = (float)0.8;
    params->query_importance = (float)1;
    params->normalization_value = (float)0.01;
    params->integration_lambda = (float)0.0001;
    params->integration_iteration = 1;
    params->enhancement_type = ITSELF_BY_SALIENCY;

    params->prior_params.center_variance = (float)0.0;
    params->prior_params.color_variance = (float)0.0;
    params->prior_params.ellipse_variance = (float)0.02;
    params->prior_params.focus_variance = (float)0.2;
    params->prior_params.red_yellow_variance = (float)0.0;
    params->prior_params.white_variance = (float)0.0;
    params->prior_params.foreground_distance_variance = (float)0.02;
    params->prior_params.foreground_feature_variance = (float)0.05;

    return params;
}

iftITSELFParameters *iftInitializeITSELFParametersByDefaultScribbles(){
    iftITSELFParameters *params = iftCreateITSELFParameters();
    params->itself_iterations = 12;
    params->number_superpixels = 2400;
    params->obj_seeds_number = 10;
    params->oisf_gamma = (float)1.0;
    params->oisf_iterations = 5;
    params->oisf_beta = 12;
    params->superpixel_increase_ratio = (float)0.8;
    params->oisf_alpha = (float)0.8;
    params->query_importance = (float)1;
    params->normalization_value = (float)0.01;
    params->integration_lambda = (float)0.0001;
    params->integration_iteration = 1;
    params->enhancement_type = ITSELF_BY_SCRIBBLES;

    params->prior_params.center_variance = (float)0.0;
    params->prior_params.color_variance = (float)0.0;
    params->prior_params.ellipse_variance = (float)0.0;
    params->prior_params.focus_variance = (float)0.0;
    params->prior_params.red_yellow_variance = (float)0.0;
    params->prior_params.white_variance = (float)0.0;
    params->prior_params.foreground_distance_variance = (float)0.00;
    params->prior_params.foreground_feature_variance = (float)0.00;

    return params;
}

void iftITSELFParametersSetParam(iftITSELFParameters *params, char *parameter_name, float value){
    if(strcmp(parameter_name, "itself_iterations") == 0)
        params->itself_iterations = (int)value;
    if(strcmp(parameter_name, "number_superpixels") == 0)
        params->number_superpixels = (int)value;
    if(strcmp(parameter_name, "obj_seeds_number") == 0)
        params->obj_seeds_number = (int)value;
    if(strcmp(parameter_name, "oisf_gamma") == 0)
        params->oisf_gamma = value;
    if(strcmp(parameter_name, "oisf_iterations") == 0)
        params->oisf_iterations = (int)value;
    if(strcmp(parameter_name, "oisf_beta") == 0)
        params->oisf_beta = value;
    if(strcmp(parameter_name, "oisf_alpha") == 0)
        params->oisf_alpha = value;
    if(strcmp(parameter_name, "oisf_dist_penalty") == 0)
        params->oisf_dist_penalty = value;
    if(strcmp(parameter_name, "superpixel_increase_ratio") == 0)
        params->superpixel_increase_ratio = value;
    if(strcmp(parameter_name, "query_importance") == 0)
        params->query_importance = value;
    if(strcmp(parameter_name, "normalization_value") == 0)
        params->normalization_value = value;
    if(strcmp(parameter_name, "integration_lambda") == 0)
        params->integration_lambda = value;
    if(strcmp(parameter_name, "integration_iteration") == 0)
        params->integration_iteration = (int)value;
    if(strcmp(parameter_name, "enhancement_type") == 0)
        params->enhancement_type = (int)value;

}

void iftCopyITSELFParameters(iftITSELFParameters *copied_to, iftITSELFParameters *to_copy){
    copied_to->itself_iterations = to_copy->itself_iterations;
    copied_to->number_superpixels = to_copy->number_superpixels;
    copied_to->obj_seeds_number = to_copy->obj_seeds_number;
    copied_to->oisf_gamma = to_copy->oisf_gamma;
    copied_to->oisf_iterations = to_copy->oisf_iterations;
    copied_to->oisf_beta = to_copy->oisf_beta;
    copied_to->oisf_alpha = to_copy->oisf_alpha;
    copied_to->superpixel_increase_ratio = to_copy->superpixel_increase_ratio;
    copied_to->query_importance = to_copy->query_importance;
    copied_to->normalization_value = to_copy->normalization_value;
    copied_to->integration_lambda = to_copy->integration_lambda;
    copied_to->integration_iteration = to_copy->integration_iteration;

    copied_to->prior_params.center_variance = to_copy->prior_params.center_variance ;
    copied_to->prior_params.color_variance = to_copy->prior_params.color_variance;
    copied_to->prior_params.ellipse_variance = to_copy->prior_params.ellipse_variance;
    copied_to->prior_params.focus_variance = to_copy->prior_params.focus_variance;
    copied_to->prior_params.red_yellow_variance = to_copy->prior_params.red_yellow_variance;
    copied_to->prior_params.white_variance = to_copy->prior_params.white_variance;
    copied_to->prior_params.foreground_distance_variance = to_copy->prior_params.foreground_distance_variance;
    copied_to->prior_params.foreground_feature_variance = to_copy->prior_params.foreground_feature_variance;
}

void iftWriteITSELFParamsToFile(iftITSELFParameters *params, char *out_file){
    FILE *fp = fopen(out_file,"wb");

    fprintf(fp, "itself_iterations,%d\n", params->itself_iterations);
    fprintf(fp, "number_superpixels,%d\n", params->number_superpixels);
    fprintf(fp, "obj_seeds_number,%d\n", params->obj_seeds_number);
    fprintf(fp, "oisf_gamma,%.2f\n", params->oisf_gamma);
    fprintf(fp, "oisf_iterations,%d\n", params->oisf_iterations);
    fprintf(fp, "oisf_beta,%.2f\n", params->oisf_beta);
    fprintf(fp, "oisf_alpha,%.2f\n", params->oisf_alpha);
    fprintf(fp, "superpixel_increase_ratio,%.2f\n", params->superpixel_increase_ratio);
    fprintf(fp, "query_importance,%.2f\n", params->query_importance);
    fprintf(fp, "normalization_value,%.2f\n", params->normalization_value);
    fprintf(fp, "integration_lambda,%.6f\n", params->integration_lambda);
    fprintf(fp, "integration_iteration,%d\n", params->integration_iteration);

    fprintf(fp, "prior_params.center_variance,%.2f\n", params->prior_params.center_variance);
    fprintf(fp, "prior_params.color_variance,%.2f\n", params->prior_params.color_variance);
    fprintf(fp, "prior_params.ellipse_variance,%.2f\n", params->prior_params.ellipse_variance);
    fprintf(fp, "prior_params.focus_variance,%.2f\n", params->prior_params.focus_variance);
    fprintf(fp, "prior_params.red_yellow_variance,%.2f\n", params->prior_params.red_yellow_variance);
    fprintf(fp, "prior_params.white_variance,%.2f\n", params->prior_params.white_variance);
    fprintf(fp, "prior_params.foreground_distance_variance,%.2f\n", params->prior_params.foreground_distance_variance);
    fprintf(fp, "prior_params.foreground_feature_variance,%.2f\n", params->prior_params.foreground_feature_variance);

    fclose(fp);
}

iftSaliencyArray *iftGBSSingleForegroundMap(iftSaliencyGraph **saliency_graph_ref, iftMImage *lab, iftImage *superpixel_img, float query_importance, iftImage *query_image, float normalization_value){
    float max_saliency, min_saliency;
    iftSaliencyGraph *saliency_graph = *saliency_graph_ref;
    int query_type = IFT_SALIENCY_MAP_QUERY;

    /*----------------------------------- Create maps for foreground seeds ------------------------------------------- */

    iftSetSaliencyQuery(superpixel_img, saliency_graph, query_type, 2, query_image);

    iftSaliencyArray *saliency;

    iftSaliencyArray *old_saliency = iftSaliencyPriorFromImage(superpixel_img, query_image);
    old_saliency->superpixel_image = superpixel_img;
    saliency = iftComputeSuperpixelGraphDissimilarityNew(saliency_graph, query_importance, old_saliency, normalization_value);

    iftDestroySaliencyArray(&old_saliency);

    saliency_graph->saliency = saliency;

    max_saliency = 0;
    min_saliency = IFT_INFINITY_FLT;

    for(int rj = 0; rj < saliency_graph->n; rj++){ //Get Max
        if(saliency->val[rj] > max_saliency)
            max_saliency = saliency->val[rj];
        if(saliency->val[rj] < min_saliency)
            min_saliency = saliency->val[rj];
    }
    for(int rj = 0; rj < saliency_graph->n; rj++){
        saliency->val[rj] = (saliency->val[rj] - min_saliency) /(max_saliency-min_saliency);
    }
    *saliency_graph_ref = saliency_graph;

    return saliency;
}

iftSaliencyArray *iftGBSSingleBackgroundMap(iftSaliencyGraph **saliency_graph_ref, iftMImage *lab, iftImage *superpixel_img, float query_importance, iftImage *query_image, float normalization_value){
    float max_saliency, min_saliency;
    iftSaliencyGraph *saliency_graph = *saliency_graph_ref;
    int query_type = IFT_SALIENCY_MAP_BACKGROUND_QUERY;

    /*-------- Create maps for background seeds ------- */

    iftSetSaliencyQuery(superpixel_img, saliency_graph, query_type, 2, query_image);
    normalization_value+=0.03;

    iftSaliencyArray *saliency;
    iftSaliencyArray *old_saliency = iftSaliencyPriorFromImage(superpixel_img, query_image);
    saliency = iftComputeSuperpixelGraphDissimilarityNew(saliency_graph, query_importance, old_saliency, normalization_value);

    iftDestroySaliencyArray(&old_saliency);

    saliency_graph->saliency = saliency;
    max_saliency = 0;
    min_saliency = IFT_INFINITY_FLT;

    for(int rj = 0; rj < saliency_graph->n; rj++){ //Get Max
        if(saliency->val[rj] > max_saliency)
            max_saliency = saliency->val[rj];
        if(saliency->val[rj] < min_saliency)
            min_saliency = saliency->val[rj];
    }
    for(int rj = 0; rj < saliency_graph->n; rj++){
        saliency->val[rj] = 1 * (saliency->val[rj] - min_saliency) /(max_saliency-min_saliency);
    }

    *saliency_graph_ref = saliency_graph;

    return saliency;
}

void iftCreateOrUpdateSaliencyGraph(iftSaliencyGraph **prev_saliency_graph, iftMImage *lab, iftImage *label_img, float adjacency_radius, iftMImage *features){
    iftIntArray *region_sizes = NULL;
    iftFloatArray *region_sizes_prob = NULL;
    iftSaliencyGraph *saliency_graph = NULL;
    int nfeats = features->m;

    /* For new graphs, create colors and color distances */
    if(*prev_saliency_graph == NULL){
        saliency_graph = (iftSaliencyGraph*) iftAlloc(1, sizeof(iftSaliencyGraph));
        saliency_graph->combined_priors = NULL;
        saliency_graph->prior_list = NULL;
        saliency_graph->improved = 1;
    } else {
        saliency_graph = *prev_saliency_graph;
        iftDestroySaliencyGraphSuperpixels(&saliency_graph);
        iftDestroyFloatArray(&(saliency_graph->region_sizes_prob));
        iftDestroyIntArray(&(saliency_graph->region_sizes));
    }

    /* Because the superpixels change, recompute every superpixel related feature */
    region_sizes = iftCountLabelSpels(label_img);
    region_sizes_prob = iftCreateFloatArray(region_sizes->n-1);
    saliency_graph->region_sizes = region_sizes;
    saliency_graph->region_sizes_prob = region_sizes_prob;
    saliency_graph->n = (int) region_sizes->n - 1;
    saliency_graph->superpixel_mean_color = (iftFloatArray **) iftAlloc((size_t) saliency_graph->n, sizeof(iftFloatArray *));
    saliency_graph->superpixel_mean_feat = (iftFloatArray **) iftAlloc((size_t) saliency_graph->n, sizeof(iftFloatArray *));
    saliency_graph->feature_number = features->m;
    for (int r = 0; r < saliency_graph->n; r++){
        saliency_graph->superpixel_mean_color[r] = iftCreateFloatArray(lab->m);
        saliency_graph->superpixel_mean_feat[r] = iftCreateFloatArray(features->m);
    }

    saliency_graph->nfeats = nfeats;

    iftFloatArray **adjacency = (iftFloatArray**) iftAlloc((size_t)saliency_graph->n, sizeof(iftFloatArray*));
    for(int r = 0; r < saliency_graph->n; r++) {
        adjacency[r] = iftCreateFloatArray(saliency_graph->n);
    }

    /*------ Create superpixel color relation and adjacency ------*/
    iftAdjRel *A = iftCircular(adjacency_radius);

    for(int p = 0; p < label_img->n; p++) {
        int superpixelP, superpixelQ;
        superpixelP = label_img->val[p] - 1;
        iftVoxel u = iftGetVoxelCoord(label_img,p);

        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftValidVoxel(label_img, v)) {
                int q = iftGetVoxelIndex(label_img, v);
                superpixelQ = label_img->val[q] - 1;
                if(superpixelP != superpixelQ) {
                    adjacency[superpixelP]->val[superpixelQ]=1;
                }
            }
        }
    }

    for(int p = 0; p < lab->n; p++){
        for(int b = 0; b < features->m; b++)
            saliency_graph->superpixel_mean_feat[label_img->val[p] -1]->val[b] += features->val[p][b];
        for(int b = 0; b < lab->m; b++)
            saliency_graph->superpixel_mean_color[label_img->val[p] -1]->val[b] += lab->val[p][b];

    }

    for(int r = 0; r < saliency_graph->n; r++){
        for(int b = 0; b < features->m; b++)
            saliency_graph->superpixel_mean_feat[r]->val[b] /= region_sizes->val[r + 1];
        for(int b = 0; b < lab->m; b++)
            saliency_graph->superpixel_mean_color[r]->val[b] /= region_sizes->val[r + 1];
    }

    iftFloatArray **extended_adjacency = (iftFloatArray**) iftAlloc((size_t) saliency_graph->n, sizeof(iftFloatArray*));
    for(int r = 0; r < saliency_graph->n; r++)
        extended_adjacency[r] = iftCreateFloatArray(saliency_graph->n);

    for(int r = 0; r < saliency_graph->n; r++) {
        for(int rj = 0; rj < saliency_graph->n; rj++) {
            extended_adjacency[r]->val[rj]=adjacency[r]->val[rj];
        }
    }
    for(int r =0; r < saliency_graph->n; r++){
        for(int q=0; q < saliency_graph->n; q++){
            if(adjacency[r]->val[q] > 0.0)
                for (int a = 0; a < saliency_graph->n; a++)
                    if (r != a && adjacency[q]->val[a] > 0.0){
                        extended_adjacency[r]->val[a]=1;
                    }
        }
    }

    saliency_graph->adjacency = extended_adjacency;
    iftSetSaliencyQuery(label_img, saliency_graph, IFT_NO_QUERY, 0, NULL);

    saliency_graph->saliency = NULL;

    iftDestroyAdjRel(&A);
    for(int r = 0; r < saliency_graph->n; r++) {
        iftFloatArray *aux = adjacency[r];
        iftDestroyFloatArray(&aux);
    }

    iftFree(adjacency);
    *prev_saliency_graph = saliency_graph;
}


iftSaliencyArray *iftComputeSuperpixelGraphDissimilarityNew(iftSaliencyGraph *saliency_graph, float query_importance, iftSaliencyArray *old_saliency, float variance) {
    iftSaliencyArray *saliency;
    iftIntArray *query;
    if (saliency_graph->query != NULL)
        query = saliency_graph->query;
    else {
        query = iftCreateIntArray(saliency_graph->n);
        for (int s = 0; s < saliency_graph->n; s++)
            query->val[s] = 0;
    }
    int query_number = 0;
    for(int q = 0; q < query->n; q++)
        if(query->val[q] > 0)
            query_number++;
    float epsilon = 0.000000001;
    iftSaliencyArray *saliency_by_adjacency = iftCreateSaliencyArray(saliency_graph->n);
    iftSaliencyArray *saliency_by_query = iftCreateSaliencyArray(saliency_graph->n);
    saliency = NULL;

    iftFloatArray **w = (iftFloatArray **) iftAlloc((size_t) saliency_graph->n,
                                                    sizeof(iftFloatArray *));
    for (int r = 0; r < saliency_graph->n; r++)
        w[r] = iftCreateFloatArray(saliency_graph->n);
    iftFloatArray **w_query = (iftFloatArray **) iftAlloc((size_t) saliency_graph->n,
                                                          sizeof(iftFloatArray *));
    for (int r = 0; r < saliency_graph->n; r++)
        w_query[r] = iftCreateFloatArray(saliency_graph->n);
    for (int s = 0; s < saliency_graph->n; s++) {
        for (int r = 0; r < saliency_graph->n; r++) {
            if (saliency_graph->adjacency[s]->val[r] || s == r) {
                /* ---------- COMPUTE SIMILARITY BETWEEN COLORS ------------ */
                w[s]->val[r] = expf(-((float) computeFeatureDistance(saliency_graph->superpixel_mean_feat[s]->val,saliency_graph->superpixel_mean_feat[r]->val,saliency_graph->feature_number)/variance));
            }
            if(query->val[r] || query->val[s]){
                w_query[s]->val[r] = expf(-((float) computeFeatureDistance(saliency_graph->superpixel_mean_feat[s]->val,saliency_graph->superpixel_mean_feat[r]->val,saliency_graph->feature_number)/variance));
                w_query[r]->val[s] = w_query[s]->val[r];
                //printf("Feat %f\n", saliency_graph->superpixel_mean_color[s]->val[0]);
            }
        }
    }

    saliency = iftCreateSaliencyArray(saliency_graph->n);

    for (int r = 0; r < saliency_graph->n; r++) {
        saliency_by_query->val[r] = 0;
        float adjacency_size_dist = 0;
        for (int q = 0; q < saliency_graph->n; q++) {
            if (query->val[q])
                if (q != r) {
                    if(saliency_graph->query_weight == 0)
                        saliency_by_query->val[r] += w_query[r]->val[q];
                    else
                        saliency_by_query->val[r] = iftMax(w_query[r]->val[q], saliency_by_query->val[r]);
                }
            if(saliency_graph->adjacency[r]->val[q]) {
                adjacency_size_dist += w[r]->val[q];
                saliency_by_adjacency->val[r] += w[r]->val[q]* ((old_saliency->val[q]+epsilon));
            }
        }

        if(saliency_graph->query_weight == 0)
            saliency_by_query->val[r] = expf(-saliency_by_query->val[r]*variance);
        if(adjacency_size_dist == 0)
            adjacency_size_dist = epsilon;
        saliency_by_adjacency->val[r] = saliency_by_adjacency->val[r] / adjacency_size_dist;
    }

    float max_adj_sal = 0, min_adj_sal = IFT_INFINITY_FLT;
    float max_query_sal = 0, min_query_sal = IFT_INFINITY_FLT;
    for (int r = 0; r < saliency_graph->n; r++) {
        if (saliency_by_adjacency->val[r] > max_adj_sal)
            max_adj_sal = saliency_by_adjacency->val[r];
        if (saliency_by_adjacency->val[r] < min_adj_sal)
            min_adj_sal = saliency_by_adjacency->val[r];

        if (saliency_by_query->val[r] > max_query_sal)
            max_query_sal = saliency_by_query->val[r];
        if (saliency_by_query->val[r] < min_query_sal)
            min_query_sal = saliency_by_query->val[r];
    }
    for (int r = 0; r < saliency_graph->n; r++){
        saliency_by_adjacency->val[r] = (saliency_by_adjacency->val[r] - min_adj_sal) / (max_adj_sal - min_adj_sal);
        if(min_query_sal == max_query_sal){
            saliency_by_query->val[r] = saliency_by_adjacency->val[r];
        }
        else{
            if(query->val[r] && saliency_graph->query_weight == 1)
                saliency_by_query->val[r] = max_query_sal;
            saliency_by_query->val[r] = (saliency_by_query->val[r] - min_query_sal) / (max_query_sal - min_query_sal);
        }
        saliency->val[r] = query_importance*saliency_by_query->val[r] + (1-query_importance)*saliency_by_adjacency->val[r];
    }
    float stdev_in_image = 0;
    for(int s = 0; s < saliency_graph->n; s++){
        stdev_in_image += square(saliency->val[s] - 0.5) * saliency_graph->region_sizes->val[s+1];
    }
    stdev_in_image = (stdev_in_image/old_saliency->superpixel_image->n);

    if(stdev_in_image < 0.15){
        iftImage *sal_image = iftSaliencyImageFromArray(saliency_graph->prior_list[0], saliency_graph->prior_list[0]->superpixel_image);
        iftSaliencyArray *sal_array = iftSaliencyPriorFromImage(old_saliency->superpixel_image, sal_image);
        for(int s = 0; s < saliency_graph->n; s++)
            saliency->val[s] = sal_array->val[s];
        iftDestroySaliencyArray(&sal_array);
        iftDestroyImage(&sal_image);
    }


    for (int r = 0; r < saliency_graph->n; r++) {
        iftFloatArray *aux = w_query[r];
        iftDestroyFloatArray(&aux);
    }
    iftFree(w_query);
    for (int r = 0; r < saliency_graph->n; r++) {
        iftFloatArray *aux = w[r];
        iftDestroyFloatArray(&aux);
    }

    iftFree(w);
    iftDestroySaliencyArray(&saliency_by_adjacency);
    iftDestroySaliencyArray(&saliency_by_query);
    return saliency;
}

void iftSetSaliencyQuery(iftImage *label_img, iftSaliencyGraph *saliency_graph, int query_type, int border_radius, iftImage *query_image){
    if(saliency_graph->query != NULL)
        iftDestroyIntArray(&(saliency_graph->query));
    saliency_graph->query = iftCreateIntArray(saliency_graph->n);
    int threshold;
    saliency_graph->query_weight = 0;
    iftImage *thresholded;
    saliency_graph->query_size=0;

    iftIntArray *region_saliency_sizes;

    switch(query_type){
        case IFT_TOP_BORDER_QUERY:
            for(int p=0; p < label_img->n; p++){
                iftVoxel u;
                u = iftGetVoxelCoord(label_img, p);
                if(u.y <= border_radius)
                    saliency_graph->query->val[label_img->val[p] -1] = 1;
            }
            break;

        case IFT_RIGHT_BORDER_QUERY:
            for(int p=0; p < label_img->n; p++){
                iftVoxel u;
                u = iftGetVoxelCoord(label_img, p);
                if(u.x >= label_img->xsize - border_radius)
                    saliency_graph->query->val[label_img->val[p] -1] = 1;
            }
            break;

        case IFT_BOTTOM_BORDER_QUERY:
            for(int p=0; p < label_img->n; p++){
                iftVoxel u;
                u = iftGetVoxelCoord(label_img, p);
                if(u.y >= label_img->ysize - border_radius)
                    saliency_graph->query->val[label_img->val[p] -1] = 1;
            }
            break;

        case IFT_LEFT_BORDER_QUERY:
            for(int p=0; p < label_img->n; p++){
                iftVoxel u;
                u = iftGetVoxelCoord(label_img, p);
                if(u.x <= border_radius)
                    saliency_graph->query->val[label_img->val[p] -1] = 1;
            }
            break;

        case IFT_ALL_BORDERS_QUERY:
            for(int p=0; p < label_img->n; p++){
                iftVoxel u;
                u = iftGetVoxelCoord(label_img, p);
                if(u.y <= border_radius || u.y >= label_img->ysize - border_radius || u.x <= border_radius || u.x >= label_img->xsize - border_radius)
                    saliency_graph->query->val[label_img->val[p] -1] = 1;
            }
            break;
        case IFT_NO_QUERY:
            for(int r=0; r < saliency_graph->n; r++){
                saliency_graph->query->val[r] = 0;
            }
            break;
        case IFT_SALIENCY_MAP_FOREGROUND_QUERY:

            threshold = iftOtsu(query_image);
            int max_image_sal = iftMaximumValue(query_image);
            if(threshold > max_image_sal)
                threshold = max_image_sal-1;
            region_saliency_sizes = iftCreateIntArray(saliency_graph->region_sizes->n-1);
            saliency_graph->query_threshold = threshold;
            for(int p=0; p < query_image->n; p++){
                if(query_image->val[p] >= threshold)
                    region_saliency_sizes->val[label_img->val[p]-1]++;
            }
            for(int r=0; r < region_saliency_sizes->n; r++) {
                if ((float) region_saliency_sizes->val[r] / (float)saliency_graph->region_sizes->val[r + 1] >= 0.5){
                    saliency_graph->query->val[r] = 1;
                }

            }

            thresholded = iftCreateImage(label_img->xsize, label_img->ysize, label_img->zsize);

            for(int p = 0; p < label_img->n; p++)
                if(saliency_graph->query->val[label_img->val[p]-1])
                    thresholded->val[p] = 255;
                else
                    thresholded->val[p] = 0;

            iftDestroyImage(&thresholded);
            iftDestroyIntArray(&region_saliency_sizes);
            saliency_graph->query_weight = 1;

            break;
        case IFT_SALIENCY_MAP_BACKGROUND_QUERY:\

            threshold = iftOtsu(query_image);
            saliency_graph->query_threshold = threshold;

            region_saliency_sizes = iftCreateIntArray(saliency_graph->region_sizes->n-1);
            for(int p=0; p < query_image->n; p++){
                if(query_image->val[p] < threshold) {
                    region_saliency_sizes->val[label_img->val[p] - 1]++;
                }
            }
            for(int r=0; r < region_saliency_sizes->n; r++) {
                if (((float) region_saliency_sizes->val[r] / (float)saliency_graph->region_sizes->val[r + 1]) >= 0.5) {
                    saliency_graph->query->val[r] = 1;
                }
            }

            iftDestroyIntArray(&region_saliency_sizes);

            break;
        case IFT_SEED_FOREGROUND_QUERY:
            for(int p=0; p < query_image->n; p++){
                if(query_image->val[p] == 2)
                    saliency_graph->query->val[label_img->val[p] - 1] = 1;
            }
            saliency_graph->query_weight = 1;
            break;
        case IFT_SEED_BACKGROUND_QUERY:
            for(int p=0; p < query_image->n; p++){
                if(query_image->val[p] == 1)
                    saliency_graph->query->val[label_img->val[p] - 1] = 1;
            }
            break;

        case IFT_SALIENCY_MAP_JOINT_QUERY:

            threshold = iftOtsu(query_image);
            max_image_sal = iftMaximumValue(query_image);
            if(threshold > max_image_sal)
                threshold = max_image_sal-1;
            saliency_graph->query_threshold = threshold;
            region_saliency_sizes = iftCreateIntArray(saliency_graph->region_sizes->n-1);
            for(int p=0; p < query_image->n; p++){
                if(query_image->val[p] >= threshold)
                    region_saliency_sizes->val[label_img->val[p]-1]++;
            }
            for(int r=0; r < region_saliency_sizes->n; r++) {
                if ((float) region_saliency_sizes->val[r] / (float)saliency_graph->region_sizes->val[r + 1] >= 0.5){
                    saliency_graph->query->val[r] = 2;
                }

            }

            for(int s=0; s < region_saliency_sizes->n; s++)
                region_saliency_sizes->val[s] = 0;

            for(int p=0; p < query_image->n; p++){
                if(query_image->val[p] < threshold)
                    region_saliency_sizes->val[label_img->val[p]-1]++;
            }
            for(int r=0; r < region_saliency_sizes->n; r++) {
                if ((float) region_saliency_sizes->val[r] / (float)saliency_graph->region_sizes->val[r + 1] >= 0.5){
                    saliency_graph->query->val[r] = 1;
                }

            }

            iftDestroyIntArray(&region_saliency_sizes);
            saliency_graph->query_weight = 1;

            break;
        default:
            for(int r=0; r < saliency_graph->n; r++){
                saliency_graph->query->val[r] = 0;
            }
    }
    for(int p = 0; p < label_img->n; p++)
        if(saliency_graph->query->val[label_img->val[p]-1])
            saliency_graph->query_size+=1;
}

void iftDestroySaliencyGraphSuperpixels(iftSaliencyGraph **saliency_graph){
    iftSaliencyGraph *aux = *saliency_graph;
    iftFloatArray *aux2;
    if(aux != NULL){
        if(aux->adjacency != NULL){
            for(int r = 0; r < aux->n; r++) {
                aux2 = aux->adjacency[r];
                iftDestroyFloatArray(&aux2);
            }
            iftFree(aux->adjacency);
        }

        if(aux->superpixel_mean_color != NULL){
            for(int r = 0; r < aux->n; r++) {
                aux2 = aux->superpixel_mean_color[r];
                iftDestroyFloatArray(&aux2);
            }
            iftFree(aux->superpixel_mean_color);
        }

        if(aux->superpixel_mean_feat != NULL){
            for(int r = 0; r < aux->n; r++) {
                aux2 = aux->superpixel_mean_feat[r];
                iftDestroyFloatArray(&aux2);
            }
            iftFree(aux->superpixel_mean_feat);
        }
    }

}

void iftDestroySaliencyGraph(iftSaliencyGraph **saliency_graph){
    iftSaliencyGraph *aux = *saliency_graph;
    iftFloatArray *aux2;
    if(aux != NULL){
        for(int r = 0; r < aux->n; r++) {
            aux2 = aux->adjacency[r];
            iftDestroyFloatArray(&aux2);
        }
        if(aux->adjacency != NULL)
            iftFree(aux->adjacency);

        for(int r = 0; r < aux->n; r++) {
            aux2 = aux->superpixel_mean_color[r];
            iftDestroyFloatArray(&aux2);
        }

        iftFree((aux->superpixel_mean_color));

        for(int r = 0; r < aux->n; r++) {
            aux2 = aux->superpixel_mean_feat[r];
            iftDestroyFloatArray(&aux2);
        }

        iftFree((aux->superpixel_mean_feat));

        iftDestroyIntArray(&(aux->query));

        if(aux->saliency != NULL)
            iftDestroySaliencyArray(&(aux->saliency));

        iftDestroyPriorList(saliency_graph);

        if(aux->combined_priors!=NULL)
            iftDestroySaliencyArray(&(aux->combined_priors));
        aux->combined_priors = NULL;

        iftDestroyIntArray(&aux->region_sizes);
        iftDestroyFloatArray(&aux->region_sizes_prob);

        iftFree(aux);
        *saliency_graph = NULL;
    }

}

double computeColorDistance(iftFColor color1, iftFColor color2){
    double distance=0;

    for (int i =0; i<3; i++) {
        distance += square((double) (color1.val[i] - color2.val[i]));
    }
    distance = sqrt(distance);
    return distance;
}

double computeFeatureDistance(float *feat1, float *feat2, int feat_number){
    double distance=0;

    for (int i =0; i<feat_number; i++)
        distance+= square((double)(feat1[i] - feat2[i]));
    distance = sqrt(distance);
    return distance;
}


void iftCuboidPriorIntegrationImplicit(iftSaliencyGraph **ref_saliency_graph, iftAdjRel *A, iftImage *label_img, iftITSELFParameters *params){
    iftSaliencyGraph *saliencyGraph = *ref_saliency_graph;
    int number_maps = saliencyGraph->npriors;
    iftImage **prior_images = (iftImage**) iftAlloc((size_t) number_maps, sizeof(iftImage*));
    for(int i = 0; i < number_maps; i++)
        prior_images[i] = iftCreateImageFromImage(label_img);

    for(int i = 0; i < number_maps; i++)
        for(int p = 0; p < label_img->n; p++){
            int superpixel = label_img->val[p]-1;
            prior_images[i]->val[p] = (int)(255 * saliencyGraph->prior_list[i]->val[superpixel]);
        }

    iftFImage **log_saliency = (iftFImage**) iftAlloc((size_t) number_maps, sizeof(iftImage*));
    for(int i = 0; i < number_maps; i++)
        log_saliency[i] = iftCreateFImage(label_img->xsize, label_img->ysize, label_img->zsize);

    iftIntArray *thresholds = iftCreateIntArray(number_maps);

    /* --------- Compute Automata States ---------- */
    for(int t = 0; t < params->integration_iteration; t++) {
        for(int i = 0; i < number_maps; i++) {
            thresholds->val[i] = iftOtsu(prior_images[i]);
        }

        iftImage *summed_image = iftCreateImageFromImage(label_img);
        for (int p = 0; p < label_img->n; p++) {
            iftVoxel u = iftGetVoxelCoord(label_img, p);
            for (int i = 0; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(label_img, v)) {
                    int q = iftGetVoxelIndex(label_img, v);
                    for (int sm = 0; sm < number_maps; sm++) {
                        if (prior_images[sm]->val[q] > thresholds->val[sm])
                            summed_image->val[p] += 1;
                        else
                            summed_image->val[p] -= 1;
                    }
                }
            }
        }

        float epsilon = (float)0.00000001;
        for (int p = 0; p < label_img->n; p++)
            for (int i = 0; i < number_maps; i++){
                float sal = (float)prior_images[i]->val[p]/255;
                float log_sal = logf(sal+ epsilon)/(1-sal + epsilon);
                log_saliency[i]->val[p] = (1 * (log_sal + params->integration_lambda * (float)summed_image->val[p]));
            }

        for (int p = 0; p < label_img->n; p++)
            for (int i = 0; i < number_maps; i++){
                float ex = expf(log_saliency[i]->val[p]);
                prior_images[i]->val[p] = (int)(255 * ex/(1+ex));
            }

        iftDestroyImage(&summed_image);
    }

    iftImage *combined_prior_image = iftCreateImageFromImage(label_img);

    for(int p = 0; p < label_img->n; p++)
        for(int i = 0; i < number_maps; i++)
            combined_prior_image->val[p] += prior_images[i]->val[p];


    int max = IFT_INFINITY_INT_NEG, min = IFT_INFINITY_INT;

    for(int p = 0; p < label_img->n; p++){
        if(combined_prior_image->val[p] > max)
            max = combined_prior_image->val[p];
        if(combined_prior_image->val[p] < min)
            min = combined_prior_image->val[p];
    }

    for(int p = 0; p < label_img->n; p++)
        combined_prior_image->val[p] = 255*(combined_prior_image->val[p] - min) / (max - min);

    iftSaliencyArray *combined_priors = iftSaliencyPriorFromImage(label_img, combined_prior_image);

    for(int r = 0; r < number_maps; r++) {
        iftImage *aux2 = prior_images[r];
        iftDestroyImage(&aux2);
        iftFImage *aux = log_saliency[r];
        iftDestroyFImage(&aux);
    }
    iftFree(prior_images);
    iftFree(log_saliency);
    iftDestroyIntArray(&thresholds);
    iftDestroyImage(&combined_prior_image);

    iftDestroySaliencyArray(&(saliencyGraph->combined_priors));

    saliencyGraph->combined_priors = combined_priors;

}

iftImage *iftCuboidSaliencyIntegration(iftImage **maps, int number_maps, iftAdjRel *A, int number_iterations, float lambda){
    iftImage **binary_prior_images = (iftImage**) iftAlloc((size_t) number_maps, sizeof(iftImage*));
    for(int i = 0; i < number_maps; i++)
        binary_prior_images[i] = iftCreateImageFromImage(maps[0]);

    iftImage **new_maps = (iftImage**) iftAlloc((size_t) number_maps, sizeof(iftImage*));
    for(int i = 0; i < number_maps; i++)
        new_maps[i] = iftCreateImage(maps[0]->xsize, maps[0]->ysize, maps[0]->zsize);

    iftFImage **log_images = (iftFImage**) iftAlloc((size_t) number_maps, sizeof(iftFImage*));
    for(int i = 0; i < number_maps; i++)
        log_images[i] = iftCreateFImage(maps[0]->xsize, maps[0]->ysize, maps[0]->zsize);

    iftFImage **log_saliency = (iftFImage**) iftAlloc((size_t) number_maps, sizeof(iftFImage*));
    for(int i = 0; i < number_maps; i++)
        log_saliency[i] = iftCreateFImage(maps[0]->xsize, maps[0]->ysize, maps[0]->zsize);

    iftIntArray *thresholds = iftCreateIntArray(number_maps);

    for(int p = 0; p < maps[0]->n; p++){
        for(int i = 0; i < number_maps; i++) {
            new_maps[i]->val[p] = maps[i]->val[p];
        }
    }

    /* ---------- Compute Automata States ------------ */
    float mean_threshold = 0;
    for(int t = 0; t < number_iterations; t++) {
        for(int i = 0; i < number_maps; i++) {
            thresholds->val[i] = iftMin(iftMaximumValue(new_maps[i]) - 1, iftOtsu(new_maps[i]));
            mean_threshold+=(float)thresholds->val[i];
        }
        mean_threshold/=(float)number_maps;
        for (int i = 0; i < number_maps; i++) {
            for (int p = 0; p < maps[0]->n; p++)
                if ((float)maps[i]->val[p] > mean_threshold)
                    binary_prior_images[i]->val[p] = 1;
                else
                    binary_prior_images[i]->val[p] = -1;
        }

        iftImage *summed_image = iftCreateImageFromImage(maps[0]);


        for (int p = 0; p < maps[0]->n; p++) {
            iftVoxel u = iftGetVoxelCoord(maps[0], p);
            for (int i = 0; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(maps[0], v)) {
                    int q = iftGetVoxelIndex(maps[0], v);
                    for (int sm = 0; sm < number_maps; sm++) {
                        summed_image->val[p] += binary_prior_images[sm]->val[q];
                    }
                }
            }
        }

        float epsilon = (float)0.0000000000001;
        for (int p = 0; p < maps[0]->n; p++)
            for (int i = 0; i < number_maps; i++) {
                float sal = (float)maps[i]->val[p]/255.0;
                double log_sal = logf(sal+epsilon)/(1-sal+epsilon);
                log_images[i]->val[p] = (float)log_sal;
            }

        for (int p = 0; p < maps[0]->n; p++)
            for (int i = 0; i < number_maps; i++)
                log_saliency[i]->val[p] = (1 * (log_images[i]->val[p] + lambda * (float)summed_image->val[p]));

        for (int p = 0; p < maps[0]->n; p++)
            for (int i = 0; i < number_maps; i++){
                float ex = expf(log_saliency[i]->val[p]);
                new_maps[i]->val[p] = (int)(255 * ex/(1+ex));
            }

        for(int i = 0; i < number_maps; i++)
            for(int p=0; p < maps[0]->n; p++)
                log_images[i]->val[p] = 0;
        iftDestroyImage(&summed_image);
    }

    iftImage *combined_prior_image = iftCreateImageFromImage(maps[0]);

    for(int p = 0; p < maps[0]->n; p++)
        for(int i = 0; i < number_maps; i++)
            combined_prior_image->val[p] += new_maps[i]->val[p];


    int max = IFT_INFINITY_INT_NEG, min = IFT_INFINITY_INT;

    for(int p = 0; p < maps[0]->n; p++){
        if(combined_prior_image->val[p] > max)
            max = combined_prior_image->val[p];
        if(combined_prior_image->val[p] < min)
            min = combined_prior_image->val[p];
    }

    for(int p = 0; p < maps[0]->n; p++)
        combined_prior_image->val[p] = (int)(255*(float)(combined_prior_image->val[p] - min) / (float)(max - min));

    for(int r = 0; r < number_maps; r++) {
        iftImage *aux2 = binary_prior_images[r];
        iftDestroyImage(&aux2);
        aux2 = new_maps[r];
        iftDestroyImage(&aux2);
        iftFImage *aux = log_images[r];
        iftDestroyFImage(&aux);
        aux = log_saliency[r];
        iftDestroyFImage(&aux);
    }
    iftFree(log_images);
    iftFree(log_saliency);
    iftFree(new_maps);
    iftFree(binary_prior_images);
    iftDestroyIntArray(&thresholds);

    return combined_prior_image;
}

void iftWriteSaliencyImage(iftImage *label_img, iftSaliencyArray *saliency, char *filename, bool heatmap){
    iftImage *prior_image;
    if(heatmap) {
        prior_image = iftCreateColorImage(label_img->xsize, label_img->ysize, label_img->zsize, 3);

        /*---------- Convert to heatmap ------------*/
        for (int p = 0; p < label_img->n; p++) {
            float R, G, B;
            iftHeatColorMapping(saliency->val[label_img->val[p] - 1], &R, &G, &B);

            iftColor rgb_color = iftRGBColor(iftRound(R * 255), iftRound(G * 255), iftRound(B * 255));
            iftColor ycbcr_color = iftRGBtoYCbCr(rgb_color, 255);

            prior_image->val[p] = ycbcr_color.val[0];
            prior_image->Cb[p] = (ushort) ycbcr_color.val[1];
            prior_image->Cr[p] = (ushort) ycbcr_color.val[2];
        }
    }else{
        prior_image = iftCreateImageFromImage(label_img);
        for(int p = 0; p < label_img->n; p++)
            prior_image->val[p] = (int) (255 * saliency->val[label_img->val[p] - 1]);
    }

    iftWriteImageByExt(prior_image, filename);

    iftDestroyImage(&prior_image);
}

iftImage *updateSuperpixels(iftImage *original_image, iftImage *saliency_img, iftITSELFParameters *params){
    params->number_superpixels = (int)iftMax((params->number_superpixels)*params->superpixel_increase_ratio, 30);

    iftImage *segmented = iftCreateImageFromImage(saliency_img);
    int threshold = iftOtsu(saliency_img);
    for(int p = 0; p < saliency_img->n; p++)
        segmented->val[p] = saliency_img->val[p] >= threshold;
    int n_components = iftConnCompCounting(segmented);
    iftDestroyImage(&segmented);
    float obj_seeds_number = n_components*params->obj_seeds_number;

    float obj_perc=iftMin(obj_seeds_number / (float)(params->number_superpixels), 0.8);
    float alpha = params->oisf_alpha, beta = params->oisf_beta, gamma = params->oisf_gamma;
    int iter = params->oisf_iterations;

    /* ------------ Compute Superpixels ------------- */
    iftImage *mask = iftSelectImageDomain(original_image->xsize, original_image->ysize, original_image->zsize);
    iftImage *seeds = iftSamplingByOSMOX(saliency_img, mask, (params->number_superpixels), obj_perc, (float)4.0);
    iftIGraph *igraph = iftInitOISFIGraph(original_image, mask, saliency_img);
    iftIGraph *out_igraph = iftIGraphOISF(igraph, seeds, alpha, beta, gamma, iter);
    iftImage *superpixel_img = iftIGraphLabel(out_igraph);
    iftDestroyIGraph(&out_igraph);
    iftDestroyIGraph(&igraph);
    iftDestroyImage(&seeds);
    iftDestroyImage(&mask);

    return superpixel_img;
}

iftImage *computeAndUpdateSaliency(iftImage ***saliency_list_ref, int i, iftSaliencyGraph *saliency_graph, iftMImage *lab, iftMImage *mimg, iftImage *superpixel_img, iftImage *query_image, iftITSELFParameters *params) {
    iftImage **saliency_list = *saliency_list_ref;
    iftSaliencyArray *saliency_foreground = NULL;
    iftSaliencyArray *saliency_background = NULL;
    iftImage *saliency_img;
    iftAdjRel *A = iftCircular(sqrtf((float)2.0));
    saliency_foreground = iftGBSSingleForegroundMap(&saliency_graph, lab, superpixel_img, params->query_importance, query_image, params->normalization_value);
    saliency_background = iftGBSSingleBackgroundMap(&saliency_graph, lab, superpixel_img, params->query_importance, query_image, params->normalization_value);
    saliency_foreground->superpixel_image = superpixel_img;
    saliency_background->superpixel_image = superpixel_img;
    for(int s = 0; s < saliency_foreground->n; s++)
        saliency_foreground->val[s]= saliency_background->val[s] * saliency_foreground->val[s];

    saliency_img = iftSaliencyImageFromArray(saliency_foreground, superpixel_img);
    saliency_list[i+1] = saliency_img;
    iftImage *updatedSaliency = iftCuboidSaliencyIntegration(saliency_list, i+2, A, params->integration_iteration, params->integration_lambda);
    saliency_graph->saliency = iftSaliencyPriorFromImage(superpixel_img, updatedSaliency);

    iftDestroyImage(&updatedSaliency);
    saliency_img = iftSaliencyImageFromArray(saliency_graph->saliency, superpixel_img);
    iftDestroySaliencyArray(&saliency_graph->saliency);
    iftDestroyAdjRel(&A);

    return saliency_img;
}


iftImage *computeAndUpdateSaliencyByPriors(iftImage ***saliency_list_ref, int i, iftSaliencyGraph *saliency_graph, iftMImage *lab, iftMImage *mimg, iftImage *superpixel_img, iftImage *query_image, iftITSELFParameters *params) {
    iftImage **saliency_list = *saliency_list_ref;
    iftSaliencyArray *saliency_foreground = NULL;
    iftSaliencyArray *saliency_background = NULL;
    iftImage *saliency_img;
    iftAdjRel *A = iftCircular(sqrtf((float)2.0));
    saliency_foreground = iftGBSSingleForegroundMap(&saliency_graph, lab, superpixel_img, params->query_importance, query_image, params->normalization_value);
    saliency_background = iftGBSSingleBackgroundMap(&saliency_graph, lab, superpixel_img, params->query_importance, query_image, params->normalization_value);
    saliency_foreground->superpixel_image = superpixel_img;
    saliency_background->superpixel_image = superpixel_img;
    iftAddPriorToGraph(&saliency_graph, saliency_foreground);
    iftAddPriorToGraph(&saliency_graph, saliency_background);
    for(int s = 0; s < saliency_foreground->n; s++)
        saliency_foreground->val[s]= saliency_background->val[s] * saliency_foreground->val[s];

    saliency_img = iftSaliencyImageFromArray(saliency_foreground, superpixel_img);
    saliency_list[i+1] = saliency_img;
    iftImage *updatedSaliency = iftCuboidSaliencyIntegration(saliency_list, i+2, A, params->integration_iteration, params->integration_lambda);
    saliency_graph->saliency = iftSaliencyPriorFromImage(superpixel_img, updatedSaliency);

    iftDestroyImage(&updatedSaliency);
    saliency_img = iftSaliencyImageFromArray(saliency_graph->saliency, superpixel_img);
    iftDestroySaliencyArray(&saliency_graph->saliency);
    iftDestroyAdjRel(&A);

    return saliency_img;
}

void implicitSalientRegionSmoothing(iftSaliencyGraph **saliency_graph_ref){
    iftSaliencyGraph *saliency_graph = *saliency_graph_ref;
    iftSaliencyArray *smoothed_sal = iftCreateSaliencyArray(saliency_graph->n);
    iftIntArray *closest_query = iftCreateIntArray(smoothed_sal->n);

    for(int s = 0; s < smoothed_sal->n; s++)
        smoothed_sal->val[s] = saliency_graph->saliency->val[s];

    for(int i = 0; i < 5; i ++) {
        for (int s = 0; s < smoothed_sal->n; s++) {
            float min_distance = IFT_INFINITY_FLT;
            if (saliency_graph->saliency->val[s] > 0.05) {
                for (int q = 0; q < smoothed_sal->n; q++)
                    if (saliency_graph->query->val[q] && q != s) {
                        float distance = computeFeatureDistance(saliency_graph->superpixel_mean_color[s]->val,
                                                                saliency_graph->superpixel_mean_color[q]->val, 3);
                        if (distance < min_distance) {
                            closest_query->val[s] = q;
                            min_distance = distance;
                        }
                    }
                float s_sal_certainty = fabs(0.5 - saliency_graph->saliency->val[s]);
                float q_sal_certainty = fabs(0.5 - saliency_graph->saliency->val[closest_query->val[s]]);
                //            smoothed_sal->val[s] = iftMax(saliency_graph->saliency->val[closest_query->val[s]],saliency_graph->saliency->val[s]);
                if (q_sal_certainty > s_sal_certainty-0.4)
                    smoothed_sal->val[s] = smoothed_sal->val[closest_query->val[s]];
            }
        }
    }
    iftDestroySaliencyArray(&(saliency_graph->saliency));
    saliency_graph->saliency = smoothed_sal;
}

iftImage *getFinalSaliency(iftImage **saliency_list, iftSaliencyGraph *saliency_graph, iftImage *superpixel_img, iftITSELFParameters *params){
    iftImage *final_saliency;
    iftAdjRel *A = iftCircular(sqrtf((float)2.0));
    if(params->itself_iterations > 1 && saliency_graph->improved == 1){
        final_saliency = iftCuboidSaliencyIntegration(saliency_list, (params->itself_iterations+1), A, params->integration_iteration, params->integration_lambda);

        iftImage *foreground = iftSaliencyImageFromArray(saliency_graph->prior_list[0], saliency_graph->prior_list[0]->superpixel_image);
        for(int p = 0; p < final_saliency->n; p++)
            final_saliency->val[p] = iftMax(final_saliency->val[p], foreground->val[p]);
        iftSetSaliencyQuery(superpixel_img, saliency_graph, IFT_SALIENCY_MAP_BACKGROUND_QUERY, 2.0, final_saliency);
        iftSaliencyArray *old_sal = iftSaliencyPriorFromImage(superpixel_img, final_saliency);
        iftSaliencyArray *new_sal = iftComputeSuperpixelGraphDissimilarityNew(saliency_graph, 1, old_sal, params->normalization_value+0.03);

        for(int s = 0; s < saliency_graph->n; s++)
            new_sal->val[s] *= old_sal->val[s];

        iftDestroySaliencyArray(&old_sal);

        iftImage *new_sal_image = iftSaliencyImageFromArray(new_sal, superpixel_img);
        for(int p = 0; p < final_saliency->n; p++)
            final_saliency->val[p] = iftMax(new_sal_image->val[p], foreground->val[p]);
        iftDestroySaliencyArray(&new_sal);
        iftDestroyImage(&new_sal_image);
        iftDestroyImage(&foreground);
    }else
        final_saliency = iftSaliencyImageFromArray(saliency_graph->prior_list[0], saliency_graph->prior_list[0]->superpixel_image);
    iftDestroyAdjRel(&A);

    return final_saliency;
}

iftImage *iftITSELF_Priors(iftImage *orig, iftImage *initial_saliency, iftMImage *features, iftITSELFParameters *params){
    iftImage *superpixel_img = NULL;
    iftImage *saliency_img = NULL;
    float adjacency_radius = sqrtf((float)2.0);
    iftSaliencyGraph *saliency_graph = NULL;

    iftMImage *lab = convertToLab(orig);

    /* Initializing variables for first iteration */
    saliency_img = iftCopyImage(initial_saliency);
    iftImage **saliency_list = (iftImage**) iftAlloc(params->itself_iterations+1, sizeof(iftImage*));
    saliency_list[0] = initial_saliency;
    iftAdjRel *A = iftCircular(sqrtf((float)2.0));

    /* --------------------------------------- ITERATIONS -----------------------------------------------------*/
    for (int i = 0; i < params->itself_iterations; i++) {
        fflush(stdout);

        superpixel_img = updateSuperpixels(orig, saliency_img, params);

        /* ------ (Re)Compute Saliency Graph ----- */
        if(features == NULL)
            iftCreateOrUpdateSaliencyGraph(&saliency_graph, lab, superpixel_img, adjacency_radius, lab);
        else
            iftCreateOrUpdateSaliencyGraph(&saliency_graph, lab, superpixel_img, adjacency_radius, features);

        if(i==0){
            iftSaliencyArray *is = iftSaliencyPriorFromImage(superpixel_img, initial_saliency);
            is->superpixel_image = superpixel_img;
            iftAddPriorToGraph(&saliency_graph, is);
        }
        if(saliency_graph->improved == 0)
            break;
        //iftSetSaliencyQuery(superpixel_img, saliency_graph, IFT_ALL_BORDERS_QUERY, sqrt(2.0), NULL, 1);

        /* --------------- Compute Priors -------------- */

        iftDestroySaliencyArray(&(saliency_graph->saliency));

        /* ----------------- Compute Saliency --------------------*/
        /* Destroy past saliency map */
        iftImage *query_image = iftCopyImage(saliency_img);
        iftDestroyImage(&saliency_img);
        if(features == NULL)
            saliency_img = computeAndUpdateSaliency(&saliency_list, i, saliency_graph, lab, lab, superpixel_img, query_image, params);
        else
            saliency_img = computeAndUpdateSaliency(&saliency_list, i, saliency_graph, lab, features, superpixel_img, query_image, params);

        iftDestroyImage(&query_image);

//        if(i >= params.itself_iterations/2)
//            params.query_importance=0.2;
    }

    /*----------------------------------------------------------End iterations -------------------------------------------------------*/
    iftDestroyAdjRel(&A);
    params->number_superpixels = 2600;
    superpixel_img = updateSuperpixels(orig, saliency_img, params);

    iftCreateOrUpdateSaliencyGraph(&saliency_graph, lab, superpixel_img, adjacency_radius, lab);
    iftImage *final_saliency = getFinalSaliency(saliency_list, saliency_graph, superpixel_img, params);

    iftDestroyMImage(&lab);
    for(int r = 1; r < (params->itself_iterations+1); r++) {
        iftImage *aux2 = saliency_list[r];
        iftDestroyImage(&aux2);
    }
    iftFree(saliency_list);
    iftDestroyImage(&saliency_img);
    iftDestroyImage(&superpixel_img);
    iftDestroySaliencyGraph(&saliency_graph);


    return final_saliency;
}

iftImage *iftITSELF(iftImage *orig, iftImage *initial_saliency, iftMImage *features, iftITSELFParameters *params){
    iftImage *superpixel_img = NULL;
    iftImage *saliency_img = NULL;
    float adjacency_radius = sqrtf((float)2.0);
    iftSaliencyGraph *saliency_graph = NULL;

    iftMImage *lab = convertToLab(orig);

    /* Initializing variables for first iteration */
    saliency_img = iftCopyImage(initial_saliency);
    iftImage **saliency_list = (iftImage**) iftAlloc(params->itself_iterations+1, sizeof(iftImage*));
    saliency_list[0] = initial_saliency;
    iftAdjRel *A = iftCircular(sqrtf((float)2.0));

    /* --------------------------------------- ITERATIONS -----------------------------------------------------*/
    for (int i = 0; i < params->itself_iterations; i++) {
        fflush(stdout);

        superpixel_img = updateSuperpixels(orig, saliency_img, params);

        /* ------ (Re)Compute Saliency Graph ----- */
        if(features == NULL)
            iftCreateOrUpdateSaliencyGraph(&saliency_graph, lab, superpixel_img, adjacency_radius, lab);
        else
            iftCreateOrUpdateSaliencyGraph(&saliency_graph, lab, superpixel_img, adjacency_radius, features);

        if(i==0){
            iftSaliencyArray *is = iftSaliencyPriorFromImage(superpixel_img, initial_saliency);
            is->superpixel_image = superpixel_img;
            iftAddPriorToGraph(&saliency_graph, is);
        }
        if(saliency_graph->improved == 0)
            break;

        iftDestroySaliencyArray(&(saliency_graph->saliency));

        /* ----------------- Compute Saliency --------------------*/
        /* Destroy past saliency map */
        iftImage *query_image = iftCopyImage(saliency_img);
        iftDestroyImage(&saliency_img);
        if(features == NULL)
            saliency_img = computeAndUpdateSaliency(&saliency_list, i, saliency_graph, lab, lab, superpixel_img, query_image, params);
        else
            saliency_img = computeAndUpdateSaliency(&saliency_list, i, saliency_graph, lab, features, superpixel_img, query_image, params);

        iftDestroyImage(&query_image);

    }

    /*----------------------------------------------------------End iterations -------------------------------------------------------*/
    iftDestroyAdjRel(&A);
    params->number_superpixels = 2600;
    superpixel_img = updateSuperpixels(orig, saliency_img, params);

    iftCreateOrUpdateSaliencyGraph(&saliency_graph, lab, superpixel_img, adjacency_radius, lab);
    iftImage *final_saliency = getFinalSaliency(saliency_list, saliency_graph, superpixel_img, params);


    iftDestroyMImage(&lab);
    for(int r = 1; r < (params->itself_iterations+1); r++) {
        iftImage *aux2 = saliency_list[r];
        iftDestroyImage(&aux2);
    }
    iftFree(saliency_list);
    iftDestroyImage(&saliency_img);
    iftDestroyImage(&superpixel_img);
    iftDestroySaliencyGraph(&saliency_graph);


    return final_saliency;
}

iftImage *iftSESS(iftImage *original_image, iftImage *initial_saliency, iftITSELFParameters *params, iftMImage *features){
    iftImage *output = iftCreateColorImage(initial_saliency->xsize, initial_saliency->ysize, initial_saliency->zsize,
                                           iftImageDepth(initial_saliency));
    int normalization_value = iftNormalizationValue(iftMaximumValue(initial_saliency));
    iftImage *saliency_img = NULL;

    if(normalization_value < 7){
        return initial_saliency;
    }

    saliency_img = iftITSELF(original_image, initial_saliency, features, params);
    for(int p = 0; p < saliency_img->n; p++){
        output->val[p] = saliency_img->val[p];
        output->Cb[p] = 128;
        output->Cr[p] = 128;
    }

    iftDestroyImage(&saliency_img);

    return(output);
}

// ADDITIONAL STUFF THAT SHOULDN'T BE HERE
iftMImage *iftNormalizeByBand(iftMImage *feats){
    float max_feat, min_feat;
    iftMImage *normalized_feats = iftCreateMImage(feats->xsize, feats->ysize, feats->zsize, (int)feats->m);
    for(int b = 0; b < feats->m; b++){
        max_feat = 0;
        min_feat = IFT_INFINITY_FLT;
        for(int p = 0; p < feats->n; p++){
            if(feats->val[p][b] > max_feat)
                max_feat = feats->val[p][b];
            if(feats->val[p][b] < min_feat)
                min_feat = feats->val[p][b];
        }
        for(int p = 0; p < feats->n; p++){
            if(max_feat-min_feat != 0)
                normalized_feats->val[p][b] = (feats->val[p][b] - min_feat) / (max_feat - min_feat);
            else
                normalized_feats->val[p][b] = 0;
        }
    }

    return normalized_feats;
}

iftMImage *iftExtendMImageByGrayObjSalMap(iftMImage *mimg, iftImage* objsm){
    int max, norm;
    iftMImage *emimg;

    // Init
    max = iftMaximumValue(objsm);
    norm = (int) iftNormalizationValue(max);
    emimg = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, (int)mimg->m+1);

    // For all pixels
    for (int p = 0; p < mimg->n; p++)  {
        // In all bands
        for(int b = 0; b < mimg->m; b++ )  emimg->val[p][b] = mimg->val[p][b];
        // Put the value of the object saliency map in the last band
        emimg->val[p][mimg->m] = (float)objsm->val[p]/(float)norm;
    }

    return emimg;
}

iftIGraph *iftInitOISFIGraph(iftImage *img, iftImage *mask, iftImage *objsm){
    iftMImage *mimg, *obj_mimg;
    iftAdjRel *A;
    iftIGraph *igraph;
    if (iftIs3DImage(img)) A = iftSpheric((float)1.0);
    else A = iftCircular((float)1.0);
    if (iftIsColorImage(img)) mimg   = iftImageToMImage(img,LAB_CSPACE);
    else mimg   = iftImageToMImage(img,GRAY_CSPACE);
    obj_mimg = iftExtendMImageByGrayObjSalMap(mimg, objsm);
    igraph = iftImplicitIGraph( obj_mimg, mask, A);
    iftDestroyMImage(&mimg);
    iftDestroyMImage(&obj_mimg);
    iftDestroyAdjRel(&A);

    return igraph;
}

iftMImage *convertToLab(iftImage *image){
    iftMImage *lab_image;
    lab_image = iftCreateMImage(image->xsize, image->ysize, image->zsize, 3);
    float normalization_value = (float)iftNormalizationValue(iftMaximumValue(image));

    iftColor origin_color;
    iftFColor lab_color;
    int p;
    for(p = 0; p < image->n; p++){
        if(iftIsColorImage(image)){
            origin_color.val[0] = image->val[p];
            origin_color.val[1] = image->Cb[p];
            origin_color.val[2] = image->Cr[p];

            lab_color = iftYCbCrtoLabNorm(origin_color,(int) normalization_value);

            lab_image->val[p][0] = lab_color.val[0];
            lab_image->val[p][1] = lab_color.val[1];
            lab_image->val[p][2] = lab_color.val[2];
        } else{
            lab_image->val[p][0] = image->val[p];
            lab_image->val[p][1] = image->val[p];
            lab_image->val[p][2] = image->val[p];
        }
    }
    return lab_image;
}

/*
 * Return a quantized dataset and change the value of cluster_number to the number of clusters created
 */
iftDataSet *computeQuantizedDatasetKMeans(iftMImage *mimg, int max_number_of_clusters, int *cluster_number,  int maxIterations, float minImprovement){
    iftDataSet *Z;
    Z = iftMImageToDataSet(mimg, NULL, 0);
//    int old_stdout = dup(1);
//    freopen ("/dev/null", "w", stdout); // or "nul" instead of "/dev/null"
/* FOR THE PUBLISHED VERSION USE THIS CLUSTERING */
    iftClusterDataSetByKMeans(Z, max_number_of_clusters, maxIterations, minImprovement, 0, '0', 0);
//    iftKMeans(Z, max_number_of_clusters, maxIterations, minImprovement, NULL, NULL);

    *cluster_number = Z->ngroups;

    int min = IFT_INFINITY_INT;
    for(int p = 0; p < Z->nsamples; p++)
        if(Z->sample[p].group < min)
            min = Z->sample[p].group;
//    fclose(stdout);
//    stdout = fdopen(old_stdout, "w");

    return Z;
}

void iftWriteOverlay (iftImage* orig, iftImage *label, const char *filename) {
    int normvalue;
    iftImage *overlay;
    iftAdjRel *A;
    iftColor RGB, YCbCr;

    normvalue = (int)iftNormalizationValue(iftMaximumValue(orig));

    A = iftCircular((float)0.0);

    overlay = iftBorderImage(label,0);

    RGB = iftRGBColor(normvalue, normvalue, 0);
    YCbCr = iftRGBtoYCbCr(RGB, normvalue);

    iftImage *printable = iftCopyImage(orig);

    iftDrawBorders(printable,overlay,A,YCbCr,A);

    iftWriteImageByExt(printable, filename);

    iftDestroyImage(&overlay);
    iftDestroyImage(&printable);
    iftDestroyAdjRel(&A);
}

void iftWriteLabAsRGB(iftMImage *lab, char *filename){
    iftColor ycbcr_color;
    iftColor rgb_color;
    iftFColor lab_color;
    int p;
    iftImage *writable = iftCreateColorImage(lab->xsize, lab->ysize, lab->zsize, 3);
    for(p = 0; p < writable->n; p++){
        lab_color.val[0] = lab->val[p][0] * 150;
        lab_color.val[1] = lab->val[p][1] * 256 - 128;
        lab_color.val[2] = lab->val[p][2] * 256 - 128;

        rgb_color = iftLabtoRGB(lab_color,(int) 255);
        ycbcr_color = iftRGBtoYCbCr(rgb_color, 255);

        writable->val[p] = ycbcr_color.val[0];
        writable->Cb[p] = ycbcr_color.val[1];
        writable->Cr[p] = ycbcr_color.val[2];
    }
    iftWriteImageByExt(writable, filename);
    iftDestroyImage(&writable);

}


int iftConnCompCounting(const iftImage *label_img) {
    iftImage *out_label_img = iftCreateImageFromImage(label_img);

    iftIntQueue *Q = iftCreateIntQueue(label_img->n);
    iftAdjRel *A = iftSpheric((float)sqrt(3.0));

    int label = 1;

    for (int t = 0; t < label_img->n; t++) {
        if ((label_img->val[t] != 0) && (out_label_img->val[t] == 0)) {
            out_label_img->val[t] = label;
            iftInsertIntQueue(Q, t);

            while (!iftIsIntQueueEmpty(Q)) {
                int p;
                iftRemoveIntQueue(Q, &p);
                iftVoxel u = iftGetVoxelCoord(label_img, p);

                for (int i = 1; i < A->n; i++) {
                    iftVoxel v = {.x = u.x + A->dx[i], .y = u.y + A->dy[i], .z = u.z + A->dz[i]};

                    if (iftValidVoxel(label_img, v)) {
                        int q = iftGetVoxelIndex(label_img, v);

                        if ((out_label_img->val[q] == 0) && (label_img->val[p] == label_img->val[q])) {
                            out_label_img->val[q] = label;
                            iftInsertIntQueue(Q, q);
                        }
                    }
                }
            }
            label++;
        }
    }

    iftDestroyIntQueue(&Q);
    iftDestroyAdjRel(&A);

    iftDestroyImage(&out_label_img);
    return label - 1;
}

iftImage *iftSamplingByOSMOX_new(iftImage* objsm, iftImage *mask, int num_seeds, float obj_perc, float dist_penalty)
{
    if(objsm->n < num_seeds || num_seeds < 0) {
        iftError("Invalid number of seeds!", "iftSamplingByOSMOX");
    }
    else if(obj_perc < 0.0 || obj_perc > 1.0) {
        iftError("Invalid object percentage!", "iftSamplingByOSMOX");
    }
    if(mask != NULL) iftVerifyImageDomains(objsm, mask, "iftSamplingByOSMOX");

    int obj_seeds, bkg_seeds, max_val, min_val;
    iftSet *obj_set, *bkg_set, *s;
    iftImage *seed_img, *mask_copy, *invsm;

    obj_set = NULL;
    bkg_set = NULL;
    s = NULL;

    iftMinMaxValues(objsm, &min_val, &max_val);

    if( mask == NULL ) mask_copy = iftSelectImageDomain(objsm->xsize, objsm->ysize, objsm->zsize);
    else mask_copy = iftCopyImage(mask);

    seed_img = iftCreateImage(objsm->xsize, objsm->ysize, objsm->zsize);

    // Establish the number of seeds
    if(max_val == min_val)
    {
        obj_seeds = num_seeds;
        bkg_seeds = 0;
    }
    else
    {
        obj_seeds = iftRound(num_seeds * obj_perc);
        bkg_seeds = num_seeds - obj_seeds;
    }

    obj_set = iftObjSalMapSamplByHighestValue(objsm, mask_copy, obj_seeds, dist_penalty);

    invsm = iftComplement(objsm);

    s = obj_set;
    while( s != NULL ) {
        mask_copy->val[s->elem] = 0;
        seed_img->val[s->elem] = 1;
        s = s->next;
    }

    bkg_set = iftObjSalMapSamplByHighestValue(invsm, mask_copy, bkg_seeds, dist_penalty);

    s = bkg_set;
    while( s != NULL ) {
        seed_img->val[s->elem] = 1;
        s = s->next;
    }

    // Free
    iftDestroyImage(&invsm);
    iftDestroySet(&obj_set);
    iftDestroySet(&bkg_set);
    iftDestroySet(&s);
    iftDestroyImage(&mask_copy);

    return (seed_img);
}

int iftIsEmptyMap(iftImage *saliency){
    int sal = saliency->val[0];
    int is_empty = 1;
    for(int p = 0; p < saliency->n; p++)
        if(saliency->val[p] != sal){
            is_empty = 0;
            break;
        }
    return is_empty;
}
