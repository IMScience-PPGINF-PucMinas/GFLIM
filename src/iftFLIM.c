#include "iftFLIM.h"
#include "iftFLIM.cuh"

iftMImage *ReadInputMImage(char *filename);

iftLabeledSet *LabelMarkers(iftMImage *img, iftLabeledSet *S);

iftMatrix *ComputeKernelBank(iftDataSet *Z, int *ngroups, uchar grouping_method);

iftDataSet *ComputeSeedDataSet(iftMImage *img, iftLabeledSet *S, iftAdjRel *A, int nsamples);

void NormalizeImageByZScore(iftMImage *img, float *mean, float *stdev);

void NormalizeImagePerBand(iftMImage *img, iftImage *label);

iftMatrix *LearnKernelBank(iftMImage *input, iftLabeledSet *M, iftAdjRel *A, int nsamples, int nkernels_per_image, int nkernels_per_marker);

iftMatrix *ConsensusKernelbank(iftFileSet *fs_seeds, char *inputdata_dir, int noutput_channels);

void StatisticsFromAllSeeds(iftFileSet *fs_seeds, char *inputdata_dir, float *mean, float *stdev, float stdev_factor);

iftMatrix *SelectRelevantKernelsByPCA(iftMatrix *kernels, int number_of_kernels);

iftMatrix *SelectRelevantKernelsByKmeans(iftMatrix *K, int ngroups);

void ReadMeanStdev(char *basepath, float *mean, float *stdev, int ninput_channels);

void WriteMeanStdev(char *basepath, float *mean, float *stdev, int ninput_channels);

float *ReadBias(char *basepath);

void WriteBias(char *basepath, float *bias, int number_of_kernels);

iftImage *SubsampleImage(iftImage *img, int stride);

void MImageToFeatureMatrixAtRow(iftMImage *mimg, iftAdjRel *A, iftMatrix *matrix, int initRow);

void FLIMExtractFeaturesPerImage(iftFileSet *fs_images, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                 char *param_dir, char *feat_dir, char *object_dir);

void FLIMExtractFeaturesFromLayerPerImage(iftFileSet *fs_images, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                          char *param_dir, int layer_index, char *feat_dir, char *object_dir);

void FLIMExtractFeaturesPerBatch(iftFileSet *fs_images, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                 char *param_dir, char *feat_dir, char *object_dir);

void FLIMExtractFeaturesFromLayerPerBatch(iftFileSet *fs_images, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                          char *param_dir, int layer_index, char *feat_dir, char *object_dir);

iftMImage **
FLIMConvolutionalLayer(iftMImage **input, int nimages, iftImage **mask, iftFLIMArch *arch, int layer, int layer_index, int atrous_factor, char *param_dir);

iftMatrix *KernelBankByPCA(iftDataSet *Z, int *nkernels);

/* ------------ private functions --------------------------------- */

iftMatrix *KernelBankByPCA(iftDataSet *Z, int *nkernels)
{
    /* Estimate the initial kernel bank by PCA */

    iftMatrix *kern_pca = iftRotationMatrixByPCA(Z);
    iftMatrix *kern_rot = iftTransposeMatrix(kern_pca);
    iftDestroyMatrix(&kern_pca);
    *nkernels = iftMin(*nkernels, kern_rot->ncols);
    iftMatrix *kern_bank = iftCreateMatrix(*nkernels, kern_rot->nrows);
    for (int col = 0; col < *nkernels; col++)
    {
        for (int row = 0; row < kern_rot->nrows; row++)
        {
            iftMatrixElem(kern_bank, col, row) = iftMatrixElem(kern_rot, col, row);
        }
    }
    iftDestroyMatrix(&kern_rot);

    return (kern_bank);
}

void MImageToFeatureMatrixAtRow(iftMImage *mimg, iftAdjRel *A, iftMatrix *matrix, int initRow)
{
#pragma omp parallel for
    for (int p = 0; p < mimg->n; p++)
    {
        iftVoxel u = iftMGetVoxelCoord(mimg, p);
        for (int i = 0; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);
                for (int b = 0; b < mimg->m; b++)
                {
                    iftMatrixElem(matrix, b + i * mimg->m, initRow + p) = mimg->val[q][b];
                }
            }
        }
    }
}

void FLIMExtractFeaturesPerBatch(iftFileSet *fs_images, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                 char *param_dir, char *feat_dir, char *object_dir)
{
    iftMImage **input = NULL, **output = NULL;
    iftImage **mask = NULL;
    char filename[200];
    int nimages = last - first + 1;
    char ext[20];
    sprintf(ext, "%s", iftFileExt(fs_images->files[0]->path));
    char **basename = (char **)calloc(nimages, sizeof(char *));
    input = (iftMImage **)calloc(nimages, sizeof(iftMImage *));
    output = NULL;
    if (object_dir != NULL)
        mask = (iftImage **)calloc(nimages, sizeof(iftImage *));

    /* Load batch */

    for (int i = first, k = 0; i <= last; i++, k++)
    {
        sprintf(filename, "%s/%s", orig_dir, fs_images->files[i]->path);
        if (strcmp(ext, ".mimg") == 0)
        {
            input[k] = iftReadMImage(filename);
        }
        else
        {
            input[k] = ReadInputMImage(filename);
        }
        basename[k] = iftFilename(fs_images->files[i]->path, ext);
        if (mask != NULL)
        {
            sprintf(filename, "%s/%s%s", object_dir, basename[k], ext);
            if (iftFileExists(filename))
            {
                mask[k] = iftReadImageByExt(filename);
                int stride = mask[k]->xsize / input[k]->xsize;
                if (stride > 1)
                {
                    iftImage *aux = SubsampleImage(mask[k], stride);
                    iftDestroyImage(&mask[k]);
                    mask[k] = aux;
                }
            }
        }
    }

    /* For each layer do */

    for (int l = 0; l < arch->nlayers; l++)
    {
        output = FLIMConvolutionalLayer(input, nimages, mask, arch, l, l + 1, 1, param_dir);
        for (int k = 0; k < nimages; k++)
        {
            iftDestroyMImage(&input[k]);
            input[k] = output[k];
            output[k] = NULL;
        }
        iftFree(output);
    }

    for (int k = 0; k < nimages; k++)
    {
        sprintf(filename, "%s/%s.mimg", feat_dir, basename[k]);
        iftFree(basename[k]);
        iftWriteMImage(input[k], filename);
        iftDestroyMImage(&input[k]);
        if (mask != NULL)
            iftDestroyImage(&mask[k]);
    }
    iftFree(mask);
    iftFree(input);
    iftFree(basename);
}

void FLIMExtractFeaturesFromLayerPerBatch(iftFileSet *fs_images, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                          char *param_dir, int layer_index, char *feat_dir, char *object_dir)
{
    iftMImage **input = NULL, **output = NULL;
    iftImage **mask = NULL;
    char filename[200];
    int nimages = last - first + 1;
    char ext[20];
    sprintf(ext, "%s", iftFileExt(fs_images->files[0]->path));
    char **basename = (char **)calloc(nimages, sizeof(char *));
    input = (iftMImage **)calloc(nimages, sizeof(iftMImage *));
    output = NULL;
    if (object_dir != NULL)
        mask = (iftImage **)calloc(nimages, sizeof(iftImage *));

    /* Load batch */

    for (int i = first, k = 0; i <= last; i++, k++)
    {
        sprintf(filename, "%s/%s", orig_dir, fs_images->files[i]->path);
        if (strcmp(ext, ".mimg") == 0)
        {
            input[k] = iftReadMImage(filename);
        }
        else
        {
            input[k] = ReadInputMImage(filename);
        }
        basename[k] = iftFilename(fs_images->files[i]->path, ext);
        if (mask != NULL)
        {
            if (iftIs3DMImage(input[k]))
            {
                sprintf(filename, "%s/%s.nii.gz", object_dir, basename[k]);
            }
            else
            {
                sprintf(filename, "%s/%s.png", object_dir, basename[k]);
            }
            if (iftFileExists(filename))
            {
                mask[k] = iftReadImageByExt(filename);
                int stride = mask[k]->xsize / input[k]->xsize;
                if (stride > 1)
                {
                    iftImage *aux = SubsampleImage(mask[k], stride);
                    iftDestroyImage(&mask[k]);
                    mask[k] = aux;
                }
            }
        }
    }

    /* For each layer do */

    output = FLIMConvolutionalLayer(input, nimages, mask, arch, 0, layer_index, 1, param_dir);
    for (int k = 0; k < nimages; k++)
    {
        iftDestroyMImage(&input[k]);
        input[k] = output[k];
        output[k] = NULL;
    }
    iftFree(output);

    for (int k = 0; k < nimages; k++)
    {
        sprintf(filename, "%s/%s.mimg", feat_dir, basename[k]);
        iftFree(basename[k]);
        iftWriteMImage(input[k], filename);
        iftDestroyMImage(&input[k]);
        if (mask != NULL)
            iftDestroyImage(&mask[k]);
    }
    iftFree(mask);
    iftFree(input);
    iftFree(basename);
}

void FLIMGraphExtractFeaturesFromLayerPerBatch(iftFileSet *fs_graphs, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                               char *param_dir, int layer_index, char *feat_dir, char *object_dir)
{
    iftFLIMGraph **input = NULL, **output = NULL;
    iftImage **mask = NULL;
    char filename[200];
    int ngraphs = last - first + 1;
    char ext[20];
    sprintf(ext, "%s", iftFileExt(fs_graphs->files[0]->path));
    char **basename = (char **)calloc(ngraphs, sizeof(char *));
    input = (iftFLIMGraph **)calloc(ngraphs, sizeof(iftFLIMGraph *));
    output = NULL;
    if (object_dir != NULL)
        iftWarning("Mask not implemented for graphs", "FLIMGraphExtractFeaturesFromLayerPerBatch");

    /* Load batch */

    for (int i = first, k = 0; i <= last; i++, k++)
    {
        sprintf(filename, "%s/%s", orig_dir, fs_graphs->files[i]->path);
        input[k] = iftReadFLIMGraph(filename);
        basename[k] = iftFilename(fs_graphs->files[i]->path, ext);
    }

    /* For each layer do */

    output = FLIMGraphConvolutionalLayer(input, ngraphs, mask, arch, 0, layer_index, 1, param_dir);
    for (int k = 0; k < ngraphs; k++)
    {
        iftDestroyFLIMGraph(&input[k]);
        input[k] = output[k];
        output[k] = NULL;
    }
    iftFree(output);

    for (int k = 0; k < ngraphs; k++)
    {
        sprintf(filename, "%s/%s.json", feat_dir, basename[k]);
        iftFree(basename[k]);
        iftWriteFLIMGraph(input[k], filename);
        iftDestroyFLIMGraph(&input[k]);
    }

    iftFree(input);
    iftFree(basename);
}

void FLIMExtractFeaturesPerImage(iftFileSet *fs_images, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                 char *param_dir, char *feat_dir, char *object_dir)
{
    iftMImage **input = NULL, **output = NULL;
    iftImage **mask = NULL;
    char filename[200];

    // printf("\n");

    for (int i = first; i <= last; i++)
    {
        input = (iftMImage **)calloc(1, sizeof(iftMImage *));
        char ext[20];
        sprintf(ext, "%s", iftFileExt(fs_images->files[i]->path));
        sprintf(filename, "%s/%s", orig_dir, fs_images->files[i]->path);
        if (strcmp(ext, ".mimg") == 0)
        {
            input[0] = iftReadMImage(filename);
        }
        else
        {
            input[0] = ReadInputMImage(filename);
        }
        char *basename = iftFilename(fs_images->files[i]->path, ext);
        if (object_dir != NULL)
        {
            mask = (iftImage **)calloc(1, sizeof(iftImage *));
            sprintf(filename, "%s/%s%s", object_dir, basename, ext);
            if (iftFileExists(filename))
            {
                mask[0] = iftReadImageByExt(filename);
                int stride = mask[0]->xsize / input[0]->xsize;
                if (stride > 1)
                {
                    iftImage *aux = SubsampleImage(mask[0], stride);
                    iftDestroyImage(&mask[0]);
                    mask[0] = aux;
                }
            }
        }

        // printf("Processing file %s: %d of %d files\r", basename, i + 1, last - first + 1);
        // fflush(stdout);

        for (int l = 0; l < arch->nlayers; l++)
        {
            output = FLIMConvolutionalLayer(input, 1, mask, arch, l, l + 1, 1, param_dir);
            iftDestroyMImage(&input[0]);
            input[0] = output[0];
            output[0] = NULL;
            iftFree(output);
        }

        sprintf(filename, "%s/%s.mimg", feat_dir, basename);
        iftFree(basename);
        iftWriteMImage(input[0], filename);
        iftDestroyMImage(&input[0]);
        iftFree(input);
        if (mask != NULL)
        {
            iftDestroyImage(&mask[0]);
            iftFree(mask);
        }
    }
}

void FLIMExtractFeaturesFromLayerPerImage(iftFileSet *fs_images, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                          char *param_dir, int layer_index, char *feat_dir, char *object_dir)
{
    iftMImage **input = NULL, **output = NULL;
    iftImage **mask = NULL;
    char filename[200];

    // printf("\n");

    for (int i = first; i <= last; i++)
    {

        char ext[20];
        sprintf(ext, "%s", iftFileExt(fs_images->files[i]->path));
        sprintf(filename, "%s/%s", orig_dir, fs_images->files[i]->path);
        input = (iftMImage **)calloc(1, sizeof(iftMImage *));
        if (strcmp(ext, ".mimg") == 0)
        {
            input[0] = iftReadMImage(filename);
        }
        else
        {
            input[0] = ReadInputMImage(filename);
        }
        char *basename = iftFilename(fs_images->files[i]->path, ext);
        if (object_dir != NULL)
        {
            mask = (iftImage **)calloc(1, sizeof(iftImage *));
            if (iftIs3DMImage(input[0]))
            {
                sprintf(filename, "%s/%s.nii.gz", object_dir, basename);
            }
            else
            {
                sprintf(filename, "%s/%s.png", object_dir, basename);
            }
            if (iftFileExists(filename))
            {
                mask[0] = iftReadImageByExt(filename);
                int stride = mask[0]->xsize / input[0]->xsize;
                if (stride > 1)
                {
                    iftImage *aux = SubsampleImage(mask[0], stride);
                    iftDestroyImage(&mask[0]);
                    mask[0] = aux;
                }
            }
        }

        // printf("Processing file %s: %d of %d files\r", basename, i + 1, last - first + 1);
        // fflush(stdout);

        output = FLIMConvolutionalLayer(input, 1, mask, arch, 0, layer_index, 1, param_dir);
        iftDestroyMImage(&input[0]);
        input[0] = output[0];
        output[0] = NULL;
        iftFree(output);

        sprintf(filename, "%s/%s.mimg", feat_dir, basename);
        iftFree(basename);
        iftWriteMImage(input[0], filename);
        iftDestroyMImage(&input[0]);
        iftFree(input);
        if (mask != NULL)
        {
            iftDestroyImage(&mask[0]);
            iftFree(mask);
        }
    }
    // printf("end \n");
}

void FLIMGraphExtractFeaturesFromLayerPerImage(iftFileSet *fs_graphs, int first, int last, char *orig_dir, iftFLIMArch *arch,
                                               char *param_dir, int layer_index, char *feat_dir, char *object_dir)
{
    iftFLIMGraph **input = NULL, **output = NULL;
    iftImage **mask = NULL;
    char filename[200];

    // printf("\n");

    for (int i = first; i <= last; i++)
    {

        char ext[20];
        sprintf(ext, "%s", iftFileExt(fs_graphs->files[i]->path));
        sprintf(filename, "%s/%s", orig_dir, fs_graphs->files[i]->path);
        input = (iftFLIMGraph **)calloc(1, sizeof(iftFLIMGraph *));

        input[0] = iftReadFLIMGraph(filename);

        char *basename = iftFilename(fs_graphs->files[i]->path, ext);
        if (object_dir != NULL)
        {
            iftWarning("Object mask not supported for graphs", "FLIMGraphExtractFeaturesFromLayerPerImage");
        }

        // printf("Processing file %s: %d of %d files\r", basename, i + 1, last - first + 1);
        // fflush(stdout);

        output = FLIMGraphConvolutionalLayer(input, 1, mask, arch, 0, layer_index, 1, param_dir);
        iftDestroyFLIMGraph(&input[0]);
        input[0] = output[0];
        output[0] = NULL;
        iftFree(output);

        sprintf(filename, "%s/%s.json", feat_dir, basename);
        iftFree(basename);
        iftWriteFLIMGraph(input[0], filename);
        iftDestroyFLIMGraph(&input[0]);
        iftFree(input);
    }
}


iftMImage *ReadInputMImage(char *filename)
{
    iftImage *img = iftReadImageByExt(filename);
    iftMImage *input = iftImageToMImage(img, LABNorm2_CSPACE);
    iftDestroyImage(&img);
    return (input);
}

/*
Dado um conjunto S de sementes em uma imagem img,
retorna o mesmo conjunto com rótulos iguais para sementes vizinhas.
*/
iftLabeledSet *LabelMarkers(iftMImage *img, iftLabeledSet *S)
{
    iftLabeledSet *M = NULL, *seed = S;
    iftImage *markers = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftAdjRel *A = NULL;

    if (iftIs3DMImage(img))
        A = iftSpheric(sqrtf(3.0));
    else
        A = iftCircular(sqrtf(2.0));

    while (seed != NULL)
    {
        int p = seed->elem;
        markers->val[p] = 255;
        seed = seed->next;
    }

    /*
      Atribui o mesmo rótulo a sementes conectadas
      no espaço da imagem conforme a adjacência A.
      Os rótulos são valores inteiros sequenciais a partir de 1.
      */
    iftImage *lbmarkers = iftLabelComp(markers, A);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&markers);

    /*
    Constroi um conjunto de sementes rotuladas
    */
    seed = S;
    while (seed != NULL)
    {
        int p = seed->elem;
        iftInsertLabeledSet(&M, p, lbmarkers->val[p]);
        seed = seed->next;
    }

    iftDestroyImage(&lbmarkers);

    return (M);
}

void iftInsertFLIMGraphNeighborNode(unsigned long node_index, unsigned long neighbor_index, iftFLIMGraph *graph)
{
    iftFLIMGraphNodeList *tmp = (iftFLIMGraphNodeList *)malloc(sizeof(iftFLIMGraphNodeList));
    tmp->next = NULL;
    tmp->node = &(graph->nodes[neighbor_index]);
    if (graph->nodes[node_index].numNeighbors == 0)
    {
        graph->nodes[node_index].neighbors_list_head = tmp;
        graph->nodes[node_index].neighbors_list_tail = tmp;
        graph->nodes[node_index].numNeighbors++;
        return;
    }

    graph->nodes[node_index].neighbors_list_tail->next = tmp;
    graph->nodes[node_index].neighbors_list_tail = tmp;
    graph->nodes[node_index].numNeighbors++;
}

/*
Dado um conjunto S de sementes em um grafo graph,
retorna o mesmo conjunto com rótulos iguais para sementes vizinhas.
*/
iftLabeledSet *LabelMarkersGraph(iftFLIMGraph *graph, iftLabeledSet *S)
{
    iftLabeledSet *M = NULL, *seed = S;
    iftDHeap *queue;

    double *markers = (double *)malloc((int)(graph->num_nodes) * sizeof(double));
    double max_label = (double)(graph->num_nodes) + 1.0;
    queue = iftCreateDHeap((int)(graph->num_nodes), markers);

    if (markers == NULL)
    {
        iftError("Error to allocate memory for markers", "LabelMarkersGraph");
    }

    for (int i = 0; i < graph->num_nodes; i++)
        markers[i] = max_label + 1;

    while (seed != NULL)
    {
        int p = seed->elem;
        markers[p] = max_label;
        iftInsertDHeap(queue, p);
        seed = seed->next;
    }

    double curr_label = 1.0;
    while (!iftEmptyDHeap(queue))
    {
        int node_index = iftRemoveDHeap(queue);

        if (markers[node_index] == max_label)
        { // unconnected seed
            markers[node_index] = curr_label;
            curr_label++;
        }
        double node_label = markers[node_index];

        iftFLIMGraphNodeList *node_list_ptr = graph->nodes[node_index].neighbors_list_head;

        while (node_list_ptr != NULL)
        {
            int adj_index = (int)(node_list_ptr->node->index);

            if (markers[adj_index] == max_label && queue->color[adj_index] == IFT_GRAY)
            {
                markers[adj_index] = node_label;
                iftGoUpDHeap(queue, queue->pos[adj_index]);
            }
            node_list_ptr = node_list_ptr->next;
        }
    }

    /*
    Constroi um conjunto de sementes rotuladas
    */
    seed = S;
    while (seed != NULL)
    {
        int p = seed->elem;
        iftInsertLabeledSet(&M, p, (int)(markers[p]));
        seed = seed->next;
    }

    iftFree(markers);
    iftDestroyDHeap(&queue);

    return (M);
}

iftDataSet *ComputeSeedDataSet(iftMImage *img, iftLabeledSet *S, iftAdjRel *A, int nsamples)
{
    int tensor_size = img->m * A->n;
    iftDataSet *Z = iftCreateDataSet(nsamples, tensor_size);
    int ninput_channels = img->m;

    Z->nclasses = 0;
    int s = 0;
    iftLabeledSet *seed = S;
    while (seed != NULL)
    {
        int p = seed->elem;
        Z->sample[s].id = s;
        Z->sample[s].truelabel = seed->label;
        if (Z->sample[s].truelabel > Z->nclasses)
            Z->nclasses = Z->sample[s].truelabel;
        iftVoxel u = iftMGetVoxelCoord(img, p);
        int j = 0;
        for (int k = 0; k < A->n; k++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, k);
            if (iftMValidVoxel(img, v))
            {
                int q = iftMGetVoxelIndex(img, v);
                for (int b = 0; b < ninput_channels; b++)
                {
                    Z->sample[s].feat[j] = img->val[q][b];
                    j++;
                }
            }
            else
            {
                for (int b = 0; b < img->m; b++)
                {
                    Z->sample[s].feat[j] = 0;
                    j++;
                }
            }
        }
        s++;
        seed = seed->next;
    }

    iftSetStatus(Z, IFT_TRAIN);
    iftAddStatus(Z, IFT_SUPERVISED);

    return (Z);
}


/**
 * @brief Insert in "maxHeap" at most k immediate neighbors to the nodes in "roots" that are
   most similar (Euclidean distance) with the features in "feature_reference".
   These neighbors must have a white state in "maxHeap".
 * @author Isabela
 * @date abr, 2023
 * @param graph: graph with node relations and features
 * @param roots: initial nodes indexes to get the neighbors
 * @param k: number of neighbors to get
 * @param dissimilarity: array to store the dissimilarity between the reference node features and the found neighbors.
 * @param feature_reference: reference node features
 * @param maxHeap: heap of nodes to store the immediate neighbors found
 */
void getNeighborsFLIMGraph(iftFLIMGraph *graph, iftSet *roots, int k, double *dissimilarity, double *feature_reference, iftDHeap *maxHeap)
{
    if (roots == NULL)
        printf("ALGUMA COISA ESTÁ ERRADA AQUI... \n");
    iftSet *tmp_roots = roots;
    while (tmp_roots != NULL)
    {

        iftFLIMGraphNodeList *neigbors_list = graph->nodes[tmp_roots->elem].neighbors_list_head;
        while (neigbors_list != NULL)
        {
            if (maxHeap->color[neigbors_list->node->index] == IFT_WHITE)
            {
                double dist = 0.0;
                for (int i = 0; i < graph->num_feats; i++)
                {
                    dist += (feature_reference[i] - neigbors_list->node->feats[i]) * (feature_reference[i] - neigbors_list->node->feats[i]);
                }

                if (maxHeap->last + 1 == k)
                {
                    int last = iftRemoveDHeap(maxHeap);
                    if (dist < dissimilarity[last])
                    {
                        dissimilarity[neigbors_list->node->index] = dist;
                        iftInsertDHeap(maxHeap, neigbors_list->node->index);
                    }
                    else
                        iftInsertDHeap(maxHeap, last);
                }
                else
                {
                    dissimilarity[neigbors_list->node->index] = dist;
                    iftInsertDHeap(maxHeap, neigbors_list->node->index);
                }
            }
            neigbors_list = neigbors_list->next;
        }
        tmp_roots = tmp_roots->next;
    }
}

/**
 * @brief Insert in minHeap the k neighbors of node_index according to the graph.
 * roots, dissimilarity, and maxHeap may have NULL value.
 * @author Isabela
 * @date abr, 2023
 * @param node_index: reference node to get the neighbors
 * @param graph: graph with node relations and features
 * @param k: number of neighbors
 * @param dissimilarity: (optional) auxiliary structure to compute the dissimilarity between the node and the neighbors
 * @param maxHeap: (optional) heap of nodes to store the immediate neighbors
 * @param minHeap: heap of nodes to store the k neighbors
 * @param dilation: number of skips to get neighbors (>= 1) (similar to dilation in convolutional network)
 */
void getKNeigborsFLIMGraph(int node_index, iftFLIMGraph *graph, int k,
                           double *dissimilarity, iftSet *roots,
                           iftDHeap *maxHeap, iftDHeap *minHeap,
                           int dilation)
{
    bool has_maxHeap = false;
    bool has_minHeap = false;
    bool has_dissimilarity = false;

    if (maxHeap != NULL)
        has_maxHeap = true;
    if (minHeap != NULL)
        has_minHeap = true;
    if (dissimilarity != NULL)
        has_dissimilarity = true;

    // MaxHeap para armazenar até k vizinhos mais similares.
    // Se a quantidade de vizinhos for maior que k, remove o menos similar
    // do heap (que está no topo) usando "pop" e verifica quem é o mais similar
    // para nova inserção.
    if (!has_maxHeap)
        maxHeap = iftCreateDHeap(graph->num_nodes, dissimilarity);
    maxHeap->removal_policy = MAXVALUE;

    // MinHeap armazena todos os vizinhos que serão incluídos na matriz de saída.
    // O minHeap fornece a ordem de inclusão na matriz de saída
    // (do mais similar para o menos similar).
    // A raíz do mineap deve sempre ser o vértice inicial,
    // pois possui maior similaridade consigo mesmo.
    if (!has_minHeap)
        minHeap = iftCreateDHeap((int)graph->num_nodes, dissimilarity);

    if (!has_dissimilarity)
        dissimilarity = (double *)malloc(graph->num_nodes * sizeof(double));

    dissimilarity[node_index] = 0.0;
    iftInsertDHeap(maxHeap, node_index);
    iftRemoveDHeap(maxHeap);
    iftInsertDHeap(minHeap, node_index);
    iftInsertSet(&roots, node_index); // initial nodes

    int found_neighbors = 1;
    int skips_performed = 1;

    while (found_neighbors < k)
    {
        int miss_neighbors = k - found_neighbors;
        getNeighborsFLIMGraph(graph, roots, miss_neighbors, dissimilarity, graph->nodes[node_index].feats, maxHeap);

        iftDestroySet(&roots);
        if (iftEmptyDHeap(maxHeap))
            break; // there are no neighbors to include

        bool new_iter = (maxHeap->last + 1 < miss_neighbors || skips_performed != dilation);

        // adiciona os k vizinhos selecionados (no maxHeap) em um minHeap
        while (!iftEmptyDHeap(maxHeap))
        {
            int node_index = iftRemoveDHeap(maxHeap);
            if (new_iter)
                iftInsertSet(&roots, node_index);
            if (skips_performed == dilation)
            {
                iftInsertDHeap(minHeap, node_index);
                found_neighbors++;
            }
        }
        if (skips_performed != dilation)
            skips_performed++;
    }
    if (maxHeap != NULL && has_maxHeap)
        iftResetDHeap(maxHeap);

    if (!has_dissimilarity)
        free(dissimilarity);
    if (!has_maxHeap)
        iftDestroyDHeap(&maxHeap);
}


iftDataSet *ComputeSeedDataSetGraph(iftFLIMGraph *graph, iftLabeledSet *S, int *kernel_size, int nsamples)
{

    int kernel_neighbors = iftMin(graph->num_nodes, kernel_size[0] * kernel_size[1]);
    int tensor_size = graph->num_feats * kernel_neighbors;

    iftDataSet *Z = iftCreateDataSet(nsamples, tensor_size);

    double *dissimilarity = (double *)malloc(graph->num_nodes * sizeof(double));

    // MaxHeap para armazenar até k vizinhos mais similares.
    // Se a quantidade de vizinhos for maior que k, remove o menos similar
    // do heap (que está no topo) usando "pop" e verifica quem é o mais similar
    // para nova inserção.
    iftDHeap *maxHeap = NULL;
    maxHeap = iftCreateDHeap((int)(graph->num_nodes), dissimilarity);
    maxHeap->removal_policy = MAXVALUE;

    // MinHeap armazena todos os vizinhos que serão incluídos na matriz de saída.
    // O minHeap fornece a ordem de inclusão na matriz de saída
    // (do mais similar para o menos similar).
    // A raíz do mineap deve sempre ser o vértice inicial,
    // pois possui maior similaridade consigo mesmo.
    iftDHeap *minHeap;
    minHeap = iftCreateDHeap((int)(graph->num_nodes), dissimilarity);

    iftSet *roots = NULL;

    Z->nclasses = 0;
    int s = 0;
    iftLabeledSet *seed = S;

    while (seed != NULL)
    {
        int p = seed->elem;

        Z->sample[s].id = s;
        Z->sample[s].truelabel = seed->label;
        if (Z->sample[s].truelabel > Z->nclasses)
            Z->nclasses = Z->sample[s].truelabel;

        getKNeigborsFLIMGraph(p, graph, kernel_neighbors,
                              dissimilarity, roots,
                              maxHeap, minHeap, 1);

        int j = 0;

        while (!iftEmptyDHeap(minHeap))
        {
            int node_index = iftRemoveDHeap(minHeap);
            for (int i = 0; i < graph->num_feats; i++)
            {
                Z->sample[s].feat[j] = graph->nodes[node_index].feats[i];
                j++;
            }
        }
        iftResetDHeap(minHeap);

        // caso não tenha vizinhos suficientes, completa o kernel com zeros
        while (j / graph->num_feats < kernel_neighbors)
        {
            for (int i = 0; i < graph->num_feats; i++)
            {
                Z->sample[s].feat[j] = 0.0;
                j++;
            }
        }
        s++;
        seed = seed->next;
    }

    iftDestroyDHeap(&maxHeap);
    iftDestroyDHeap(&minHeap);
    iftFree(dissimilarity);
    iftDestroySet(&roots);

    iftSetStatus(Z, IFT_TRAIN);
    iftAddStatus(Z, IFT_SUPERVISED);

    return (Z);
}


/* Para amostras de mesma label em Z
    (cada amostra contém as features de uma semente e k vizinhos mais próximos),
    constrói kernels (um por linha em uma matriz) no formato [Z->nfeats][ngroups].
    Caso o número de amostras (Z->nsamples) seja menor ou igual ao número de kernels
    desejado (ngroups), os kernels retornados conterão todas as amostras de Z.
    Caso contrário, faz uma aproximação de valores com kmeans, iterated watershed ou OPF.
*/
iftMatrix *ComputeKernelBank(iftDataSet *Z, int *ngroups, uchar grouping_method)
{
    iftMatrix *kernels = NULL;
    iftGraph *graph = NULL; /* For IOPF and IW */
    iftKnnGraph *knngraph = NULL;
    int i = 0;
    int *groupsize = NULL, mingroupsize = iftMax(0.05 * Z->nsamples, 1);

    iftRandomSeed(42);

    if (*ngroups >= Z->nsamples)
    { /* use all markers instead */
        *ngroups = Z->nsamples;

        // ngroups é a quantidade de kernels por marcador
        // Z->nfeats = graph->num_feats * kernel_neighbors
        // kernels[Z->nfeats][ngroups]

        kernels = iftCreateMatrix(*ngroups, Z->nfeats);

        for (int s = 0; s < Z->nsamples; s++)
        {
            iftUnitNorm(Z->sample[s].feat, Z->nfeats);
            for (int j = 0; j < Z->nfeats; j++) // Z->nfeats = graph->num_feats * kernel_neighbors
            {
                // kernels contém em cada coluna o kernel de uma semente
                // e seu numero de linhas corresponde ao número de vizinhos*features do grafo
                iftMatrixElem(kernels, i, j) = Z->sample[s].feat[j];
            }
            i++;
        }
        return (kernels);
    }

    switch (grouping_method)
    {

    case 0: /* KMeans */

        /* avoid small clusters (outliers) with less than 5% of samples */
        kernels = iftCreateMatrix(*ngroups, Z->nfeats);
        iftKMeans(Z, *ngroups, 100, 0.001, NULL, NULL, iftEuclideanDistance);
        groupsize = iftAllocIntArray(iftCountNumberOfGroupsDataSet(Z) + 1);
        for (int s = 0; s < Z->nsamples; s++)
        {
            groupsize[Z->sample[s].group]++;
        }
        for (int s = 0; s < Z->nsamples; s++)
        {
            if (iftHasSampleStatus(Z->sample[s], IFT_PROTOTYPE))
            {
                if (groupsize[Z->sample[s].group] >= mingroupsize)
                {
                    iftUnitNorm(Z->sample[s].feat, Z->nfeats);
                    for (int j = 0; j < Z->nfeats; j++)
                    {
                        iftMatrixElem(kernels, i, j) = Z->sample[s].feat[j];
                    }
                    i++;
                }
            }
        }

        break;

    case 1: /* IW in a complete graph */

        iftAddStatus(Z, IFT_TRAIN);
        graph = iftCreateGraph(Z);
        graph->is_complete = true;
        iftIteratedWatersheds(graph, *ngroups, 100, IFT_MAX, IFT_DATA_FEAT, true);
        kernels = iftCreateMatrix(*ngroups, Z->nfeats);
        for (int i = 0; i < *ngroups; i++)
        {
            int u = graph->centroids[i];
            int s = graph->node[u].sample;
            iftUnitNorm(Z->sample[s].feat, Z->nfeats);
            for (int j = 0; j < Z->nfeats; j++)
                iftMatrixElem(kernels, i, j) = Z->sample[s].feat[j];
        }
        iftDestroyGraph(&graph);

        break;

    case 2: /* OPF C Clusters */

        iftAddStatus(Z, IFT_TRAIN);
        knngraph = iftCreateKnnGraph(Z, iftMax(Z->nsamples - 1, 1));
        *ngroups = iftUnsupTrainWithCClusters(knngraph, *ngroups);

        kernels = iftCreateMatrix(*ngroups, Z->nfeats);
        for (int s = 0; s < Z->nsamples; s++)
        {
            if (iftHasSampleStatus(Z->sample[s], IFT_PROTOTYPE))
            {
                iftUnitNorm(Z->sample[s].feat, Z->nfeats);
                for (int j = 0; j < Z->nfeats; j++)
                {
                    iftMatrixElem(kernels, i, j) = Z->sample[s].feat[j];
                }
                i++;
            }
        }
        iftDestroyKnnGraph(&knngraph);

        break;

    default:
        iftError("Grouping method must be 0 for KMeans, 1 for IW", "ComputeKernelBank");
    }
    iftFree(groupsize);

    return (kernels);
}


void NormalizeImageByZScore(iftMImage *img, float *mean, float *stdev)
{
    int ninput_channels = img->m;

#pragma omp parallel for
    for (int p = 0; p < img->n; p++)
    {
        for (int b = 0; b < ninput_channels; b++)
        {
            img->val[p][b] = (img->val[p][b] - mean[b]) / stdev[b];
        }
    }
}

void NormalizeGraphByZScore(iftFLIMGraph *graph, float *mean, float *stdev)
{
    int num_feats = graph->num_feats;

    for (int p = 0; p < graph->num_nodes; p++)
    {
        for (int b = 0; b < num_feats; b++)
        {
            graph->nodes[p].feats[b] = ((float)graph->nodes[p].feats[b] - mean[b]) / stdev[b];
        }
    }
}


void NormalizeImagePerBand(iftMImage *img, iftImage *label)
{
    if (label == NULL)
    {
        for (int b = 0; b < img->m; b++)
        {
            float Imin = IFT_INFINITY_FLT, Imax = IFT_INFINITY_FLT_NEG;
            for (int p = 0; p < img->n; p++)
            {
                if (img->val[p][b] > Imax)
                    Imax = img->val[p][b];
                if ((img->val[p][b] < Imin) && (img->val[p][b] > 0))
                    Imin = img->val[p][b];
            }
            if (Imin < Imax)
                for (int p = 0; p < img->n; p++)
                {
                    img->val[p][b] = (img->val[p][b] - Imin) / (Imax - Imin);
                    if (img->val[p][b] < 0)
                        img->val[p][b] = 0;
                }
        }
    }
    else
    {
        for (int b = 0; b < img->m; b++)
        {
            float Imin = IFT_INFINITY_FLT, Imax = IFT_INFINITY_FLT_NEG;
            for (int p = 0; p < img->n; p++)
            {
                if (label->val[p] != 0)
                {
                    if (img->val[p][b] > Imax)
                        Imax = img->val[p][b];
                    if ((img->val[p][b] < Imin) && (img->val[p][b] > 0))
                        Imin = img->val[p][b];
                }
            }
            if (Imin < Imax)
                for (int p = 0; p < img->n; p++)
                {
                    if (label->val[p] != 0)
                    {
                        img->val[p][b] = (img->val[p][b] - Imin) / (Imax - Imin);
                        if (img->val[p][b] < 0)
                            img->val[p][b] = 0;
                    }
                    else
                    {
                        img->val[p][b] = 0;
                    }
                }
        }
    }
}

iftMatrix *SelectRelevantKernelsByPCA(iftMatrix *kernels, int number_of_kernels)
{
    iftMatrix *kern_t = iftTransposeMatrix(kernels);
    iftDataSet *Z = iftFeatureMatrixToDataSet(kern_t);
    iftDestroyMatrix(&kern_t);

    iftSetStatus(Z, IFT_TRAIN);
    iftNormalizeDataSetByZScoreInPlace(Z, NULL, IFT_EPSILON);

    iftMatrix *kern_pca = iftRotationMatrixByPCA(Z);
    iftMatrix *kern_rot = iftTransposeMatrix(kern_pca);
    iftDestroyMatrix(&kern_pca);

    number_of_kernels = iftMin(number_of_kernels, kern_rot->ncols);
    iftMatrix *M = iftCreateMatrix(number_of_kernels, kern_rot->nrows);

    for (int k = 0; k < number_of_kernels; k++)
    {
        for (int j = 0; j < kern_rot->nrows; j++)
        {
            iftMatrixElem(M, k, j) = iftMatrixElem(kern_rot, k, j);
        }
    }

    iftDestroyMatrix(&kern_rot);
    iftDestroyDataSet(&Z);

    return (M);
}

iftMatrix *SelectRelevantKernelsByKmeans(iftMatrix *K, int ngroups)
{
    iftMatrix *Kt = iftTransposeMatrix(K);
    iftDataSet *Z = iftFeatureMatrixToDataSet(Kt);
    iftDestroyMatrix(&Kt);

    iftRandomSeed(42);

    iftMatrix *kernels = iftCreateMatrix(ngroups, Z->nfeats);

    iftKMeans(Z, ngroups, 100, 0.001, NULL, NULL, iftCosineDistance);
    int i = 0;
    for (int s = 0; s < Z->nsamples; s++)
    {
        if (iftHasSampleStatus(Z->sample[s], IFT_PROTOTYPE))
        {
            iftUnitNorm(Z->sample[s].feat, Z->nfeats);
            for (int j = 0; j < Z->nfeats; j++)
                iftMatrixElem(kernels, i, j) = Z->sample[s].feat[j];
            i++;
        }
    }
    iftDestroyDataSet(&Z);
    return (kernels);
}


iftMatrix *LearnKernelBank(iftMImage *input, iftLabeledSet *M, iftAdjRel *A, int nsamples, int nkernels_per_image, int nkernels_per_marker)
{
    int tensor_size = input->m * A->n;

    // Z contém as informações dos kernels.
    // Z->sample contém um vetor com os kernels.
    // Cada posição é um kernel e contém o rótulo da semente que o gerou.
    // As features de cada kernel correspondem às cores da semente e de suas posições vizinhas
    iftDataSet *Z = ComputeSeedDataSet(input, M, A, nsamples);

    iftMatrix **kernels = (iftMatrix **)calloc(Z->nclasses + 1, sizeof(iftMatrix *));
    int total_nkernels = 0;

    for (int c = 1; c <= Z->nclasses; c++)
    {
        iftDataSet *Z1 = iftExtractSamplesFromClass(Z, c);
        int nkernels = nkernels_per_marker;
        /* 0: kmeans, 1: iterated watershed and 2: OPF c clusters */
        // retorna uma matriz em que
        // cada coluna contém o kernel de uma semente
        // e seu numero de linhas corresponde ao número de vizinhos*features do grafo
        kernels[c] = ComputeKernelBank(Z1, &nkernels, 0);
        total_nkernels += nkernels;
        iftDestroyDataSet(&Z1);
    }

    iftMatrix *kernelbank = iftCreateMatrix(total_nkernels, tensor_size);

    int k = 0;
    for (int c = 1; c <= Z->nclasses; c++)
    {
        for (int col = 0; col < kernels[c]->ncols; col++, k++)
        {
            for (int row = 0; row < kernels[c]->nrows; row++)
            {
                iftMatrixElem(kernelbank, k, row) = iftMatrixElem(kernels[c], col, row);
            }
        }
    }

    for (int c = 0; c <= Z->nclasses; c++)
    {
        iftDestroyMatrix(&kernels[c]);
    }
    iftFree(kernels);
    iftDestroyDataSet(&Z);

    if (kernelbank->ncols > nkernels_per_image)
    { /* force a number of kernels per image */
        iftMatrix *Rkernels = SelectRelevantKernelsByKmeans(kernelbank, nkernels_per_image);
        iftDestroyMatrix(&kernelbank);
        kernelbank = Rkernels;
    }

    return (kernelbank);
}

iftMatrix *LearnKernelBankGraph(iftFLIMGraph *graph, iftLabeledSet *M, int *kernel_size, int nsamples, int nkernels_per_image, int nkernels_per_marker)
{
    int tensor_size = graph->num_feats * kernel_size[0] * kernel_size[1];

    // Z contém as informações dos kernels.
    // Z->sample contém um vetor com os kernels.
    // Z->nclasses = numero de marcadores
    // Cada posição é um kernel e contém o rótulo da semente que o gerou. Ou seja, pode conter kernels de mesmo rótulo
    // As features de cada kernel correspondem às cores da semente e de suas posições vizinhas
    iftDataSet *Z = ComputeSeedDataSetGraph(graph, M, kernel_size, nsamples);

    // vetor de matrizes de kernels, Cada matriz é banco de kernels
    // cuja quantidade de matrizes (kernels) é igual ao número de marcadores
    iftMatrix **kernels = (iftMatrix **)calloc(Z->nclasses + 1, sizeof(iftMatrix *));
    int total_nkernels = 0;

    // para cada rótulo c, ou seja, para cada marcador
    // calcula seus kernels (uma matriz), armazenando-os em kernels[c]
    for (int c = 1; c <= Z->nclasses; c++)
    {
        iftDataSet *Z1 = iftExtractSamplesFromClass(Z, c);
        int nkernels = nkernels_per_marker;
        /* 0: kmeans, 1: iterated watershed and 2: OPF c clusters */
        kernels[c] = ComputeKernelBank(Z1, &nkernels, 0);
        total_nkernels += nkernels;
        iftDestroyDataSet(&Z1);
    }

    // kernelbank une os kernels de todos os marcadores
    // Seu tamanho é [nkernels][graph->num_feats * kernel_size[0] * kernel_size[1]],
    // em que nkernels é
    iftMatrix *kernelbank = iftCreateMatrix(total_nkernels, tensor_size);

    int k = 0;
    for (int c = 1; c <= Z->nclasses; c++)
    {
        for (int col = 0; col < kernels[c]->ncols; col++, k++)
        {
            for (int row = 0; row < kernels[c]->nrows; row++)
            {
                iftMatrixElem(kernelbank, k, row) = iftMatrixElem(kernels[c], col, row);
            }
        }
    }

    for (int c = 0; c <= Z->nclasses; c++)
    {
        iftDestroyMatrix(&kernels[c]);
    }
    iftFree(kernels);
    iftDestroyDataSet(&Z);

    if (kernelbank->ncols > nkernels_per_image)
    { /* force a number of kernels per image */
        iftMatrix *Rkernels = SelectRelevantKernelsByKmeans(kernelbank, nkernels_per_image);
        iftDestroyMatrix(&kernelbank);
        kernelbank = Rkernels;
    }

    return (kernelbank);
}

/* Create from the graph's nodes a matrix with size [num nodes][num_adj * num features],
   in which each matrix row contains the features of a node and its [num_adj] most similar neigbors.
   For nodes with degree < num_adj, the search for neighbors include neighbors of neighbors (search degree > 2).
   While the number of neighbors is lower than num_adj, search degree increases.
   The search for neighbors also stops when there aren't new nodes to increase the search.
   The similarity measure is euclidean distance in feature (e.g, color) space.
*/
iftMatrix *iftGraphToFeatureMatrix(iftFLIMGraph *graph, int num_adj, float *fill_band_with)
{
    // cria a matriz de saída
    iftMatrix *matrix = iftCreateMatrix(num_adj * graph->num_feats, graph->num_nodes);

    double *dissimilarity = (double *)malloc(graph->num_nodes * sizeof(double));

    // MaxHeap para armazenar até k vizinhos mais similares.
    // Se a quantidade de vizinhos for maior que k, remove o menos similar
    // do heap (que está no topo) usando "pop" e verifica quem é o mais similar
    // para nova inserção.
    iftDHeap *maxHeap = NULL;
    maxHeap = iftCreateDHeap(graph->num_nodes, dissimilarity);
    maxHeap->removal_policy = MAXVALUE;

    // MinHeap armazena todos os vizinhos que serão incluídos na matriz de saída.
    // O minHeap fornece a ordem de inclusão na matriz de saída
    // (do mais similar para o menos similar).
    // A raíz do mineap deve sempre ser o vértice inicial,
    // pois possui maior similaridade consigo mesmo.
    iftDHeap *minHeap;
    minHeap = iftCreateDHeap((int)graph->num_nodes, dissimilarity);
    iftSet *roots = NULL;

    for (int p = 0; p < graph->num_nodes; p++)
    {

        getKNeigborsFLIMGraph(p, graph, num_adj,
                              dissimilarity, roots,
                              maxHeap, minHeap, 1);

        // move as features dos vértices incluídos no minHeap para a matriz de saída
        int i = 0;

        while (!iftEmptyDHeap(minHeap))
        {
            int node_index = iftRemoveDHeap(minHeap);
            for (int b = 0; b < graph->num_feats; b++)
            {
                iftMatrixElem(matrix, b + i * graph->num_feats, p) = (float)(graph->nodes[node_index].feats[b]);
            }
            i++;
        }
        iftResetDHeap(minHeap);

        // caso não tenha vizinhos suficientes,
        // completa o kernel com zeros ou com o valor de fill_band_with
        float tmp[graph->num_feats];
        if (fill_band_with != NULL)
            memcpy(tmp, fill_band_with, graph->num_feats * sizeof(float));
        else
            memset(tmp, 0, graph->num_feats * sizeof(float));

        while (i < num_adj)
        {
            for (int b = 0; b < graph->num_feats; b++)
            {
                iftMatrixElem(matrix, b + i * graph->num_feats, p) = (float)(tmp[b]);
            }
            i++;
        }
    }
    iftDestroyDHeap(&maxHeap);
    iftDestroyDHeap(&minHeap);
    iftFree(dissimilarity);
    iftDestroySet(&roots);

    return matrix;
}



iftMImage **
FLIMConvolutionalLayer(iftMImage **input, int nimages, iftImage **mask, iftFLIMArch *arch, int layer, int layer_index, int atrous_factor, char *param_dir)
{
    /* Set input parameters */

    iftAdjRel *A = iftFLIMAdaptiveAdjRelFromKernel(arch->layer[layer], atrous_factor, iftIs3DMImage(input[0]));
    int ninput_channels = input[0]->m;
    char *basename = NULL;
    iftMImage **output = (iftMImage **)calloc(nimages, sizeof(iftMImage *));

    /* Read consensus parameters of the current layer */

    char filename[200];
    sprintf(filename, "%s/conv%d-kernels.npy", param_dir, layer_index);
    iftMatrix *kernelbank = iftReadMatrix(filename);
    float *mean = iftAllocFloatArray(ninput_channels);
    float *stdev = iftAllocFloatArray(ninput_channels);
    sprintf(filename, "%s/conv%d", param_dir, layer_index);
    ReadMeanStdev(filename, mean, stdev, ninput_channels);

    /* BIAS: use environment variable to use the bias idea */
    float *bias = NULL;
    char *use_bias = getenv("USE_BIAS");

    if (use_bias != NULL)
    {
        bias = ReadBias(filename);
    }

    /* Apply convolutional layer */
    for (int i = 0; i < nimages; i++)
    {
        /* marker-based normalization */

        /* BIAS: bias do not need normalization, but fill the matrix with mean */
        iftMatrix *imgM;
        if (use_bias == NULL)
        {
            NormalizeImageByZScore(input[i], mean, stdev);
            imgM = iftMImageToFeatureMatrix(input[i], A, NULL);
        }
        else
        {
            imgM = iftMImageToFeatureMatrix(input[i], A, mean);
        }

        /* convolution */
        // printf("iftMultMatrices: input[i][%d][%d][%ld] imgM[%d][%d] kernelbank[%d][%d] \n", input[i]->xsize, input[i]->ysize, input[i]->m, imgM->nrows, imgM->ncols, kernelbank->nrows, kernelbank->ncols);
        iftMatrix *conv = iftMultMatrices(imgM, kernelbank);
        iftDestroyMatrix(&imgM);
        iftMImage *activ = iftMatrixToMImage(conv, input[i]->xsize, input[i]->ysize, input[i]->zsize, kernelbank->ncols,
                                             'c');
        iftDestroyMatrix(&conv);

        if (mask != NULL)
        { /* apply mask whenever it is the case */
            int stride = mask[i]->xsize / activ->xsize;
            if (stride > 1)
            {
                iftImage *aux = SubsampleImage(mask[i], stride);
                iftDestroyImage(&mask[i]);
                mask[i] = aux;
            }

#pragma omp parallel for
            for (int p = 0; p < activ->n; p++)
            {
                if (mask[i]->val[p] == 0)
                    for (int b = 0; b < activ->m; b++)
                    {
                        activ->val[p][b] = 0;
                    }
            }
        }

        /* activation in place */

        if (arch->layer[layer].relu)
        { /* ReLU in place */
#pragma omp parallel for
            for (int p = 0; p < activ->n; p++)
            {
                for (int b = 0; b < activ->m; b++)
                {
                    /* BIAS: check for the environment variable */
                    if (use_bias != NULL)
                    {
                        activ->val[p][b] += bias[b];
                    }
                    if (activ->val[p][b] < 0)
                        activ->val[p][b] = 0;
                }
            }
        }

        /* pooling */

        if (strcmp(arch->layer[layer].pool_type, "no_pool") != 0)
        {
            iftMImage *pool = NULL;
            if (strcmp(arch->layer[layer].pool_type, "avg_pool") == 0)
            { /* ignore the stride to learn the model */
                pool = iftFLIMAtrousAveragePooling(activ, arch->layer[layer].pool_size[0], arch->layer[layer].pool_size[1], arch->layer[layer].pool_size[2], atrous_factor, arch->layer[layer].pool_stride);
                iftDestroyMImage(&activ);
                activ = pool;
            }
            else
            {
                if (strcmp(arch->layer[layer].pool_type, "max_pool") == 0)
                { /* ignore the stride to learn the model */
                    pool = iftFLIMAtrousMaxPooling(activ, arch->layer[layer].pool_size[0], arch->layer[layer].pool_size[1], arch->layer[layer].pool_size[2], atrous_factor, arch->layer[layer].pool_stride);
                    iftDestroyMImage(&activ);
                    activ = pool;
                }
                else
                {
                    iftWarning("Invalid pooling type has been ignored", "FLIMConvolutionalLayer");
                }
            }
        }

        output[i] = iftCopyMImage(activ);
        iftDestroyMImage(&activ);
    }

    iftDestroyAdjRel(&A);
    iftDestroyMatrix(&kernelbank);
    iftFree(mean);
    iftFree(stdev);
    /* BIAS: free bias array */
    if (use_bias != NULL)
    {
        iftFree(bias);
    }
    return output;
}

iftFLIMGraph **
FLIMGraphConvolutionalLayer(iftFLIMGraph **graph, int ngraphs, iftImage **mask, iftFLIMArch *arch, int layer, int layer_index, int atrous_factor, char *param_dir)
{
    // Set input parameters

    int ninput_channels = graph[0]->num_feats;
    char *basename = NULL;
    iftFLIMGraph **output_graph = (iftFLIMGraph **)calloc(ngraphs, sizeof(iftFLIMGraph *));

    // Read consensus parameters of the current layer

    char filename[200];
    sprintf(filename, "%s/conv%d-kernels.npy", param_dir, layer_index);
    iftMatrix *kernelbank = iftReadMatrix(filename);
    float *mean = iftAllocFloatArray(ninput_channels);
    float *stdev = iftAllocFloatArray(ninput_channels);
    sprintf(filename, "%s/conv%d", param_dir, layer_index);
    ReadMeanStdev(filename, mean, stdev, ninput_channels);

    // BIAS: use environment variable to use the bias idea
    float *bias = NULL;
    char *use_bias = getenv("USE_BIAS");

    if (use_bias != NULL)
    {
        bias = ReadBias(filename);
    }

    // Apply convolutional layer
    for (int i = 0; i < ngraphs; i++)
    {
        // marker-based normalization

        // BIAS: bias do not need normalization, but fill the matrix with mean
        // cria uma matriz de features de tamanho [n][k*m], em que n, k e m são
        // a quantidade de vértices, vizinhos e features, respectivamente
        iftMatrix *imgM;

        if (use_bias == NULL)
        {
            NormalizeGraphByZScore(graph[i], mean, stdev);
            imgM = iftGraphToFeatureMatrix(graph[i], arch->layer[layer].kernel_size[0] * arch->layer[layer].kernel_size[1], NULL);
        }
        else
        {
            imgM = iftGraphToFeatureMatrix(graph[i], arch->layer[layer].kernel_size[0] * arch->layer[layer].kernel_size[1], mean);
        }

        // convolution
        //printf("iftMultMatrices: graph[i]->num_nodes=%ld,->num_feats=%ld imgM[%d][%d] kernelbank[%d][%d] \n", graph[i]->num_nodes, graph[i]->num_feats, imgM->nrows, imgM->ncols, kernelbank->nrows, kernelbank->ncols);
        iftMatrix *conv = iftMultMatrices(imgM, kernelbank); // (n, m*k) x (m*k, k*k) = (n, k*k)
        iftDestroyMatrix(&imgM);

        /*if (mask != NULL)
        { // apply mask whenever it is the case
            int stride = mask[i]->xsize / activ->xsize;
            if (stride > 1)
            {
                iftImage *aux = SubsampleImage(mask[i], stride);
                iftDestroyImage(&mask[i]);
                mask[i] = aux;
            }

#pragma omp parallel for
            for (int p = 0; p < activ->n; p++)
            {
                if (mask[i]->val[p] == 0)
                    for (int b = 0; b < activ->m; b++)
                    {
                        activ->val[p][b] = 0;
                    }
            }
        }*/

        // activation in place
        if (arch->layer[layer].relu)
        { // ReLU in place
            for (int p = 0; p < conv->nrows; p++)
            {
                for (int b = 0; b < conv->ncols; b++)
                {
                    // BIAS: check for the environment variable
                    if (use_bias != NULL)
                        iftMatrixElem(conv, b, p) += bias[b];
                    if (iftMatrixElem(conv, b, p) < 0)
                        iftMatrixElem(conv, b, p) = 0;
                }
            }
        }
        // pooling

        if (strcmp(arch->layer[layer].pool_type, "no_pool") != 0)
        {
            output_graph[i] = iftMatrixToFLIMGraph(conv, graph[i]);
            isValidFLIMGraph(output_graph[i]);

            if (strcmp(arch->layer[layer].pool_type, "avg_pool") == 0)
            { // ignore the stride to learn the model
                iftFLIMGraphAtrousAveragePooling(conv, output_graph[i], arch->layer[layer].pool_size[0], arch->layer[layer].pool_size[1], arch->layer[layer].pool_size[2], atrous_factor, arch->layer[layer].pool_stride);
            }
            else
            {
                if (strcmp(arch->layer[layer].pool_type, "max_pool") == 0)
                { // ignore the stride to learn the model
                    iftFLIMGraphAtrousMaxPooling(conv, output_graph[i], arch->layer[layer].pool_size[0], arch->layer[layer].pool_size[1], arch->layer[layer].pool_size[2], atrous_factor, arch->layer[layer].pool_stride);
                }
                else
                {
                    iftWarning("Invalid pooling type has been ignore", "FLIMConvolutionalLayer");
                }
            }

            for (int j = 0; j < output_graph[i]->num_nodes; j++)
            {
                for (int b = 0; b < output_graph[i]->num_feats; b++)
                {
                    output_graph[i]->nodes[j].feats[b] = (double)(iftMatrixElem(conv, b, j));
                }
            }
        }
        else
        {
            output_graph[i] = iftMatrixToFLIMGraph(conv, graph[i]);
        }
        iftDestroyMatrix(&conv);
    }
    iftDestroyMatrix(&kernelbank);
    iftFree(mean);
    iftFree(stdev);

    // BIAS: free bias array
    if (use_bias != NULL)
    {
        iftFree(bias);
    }

    return output_graph;
}


iftMatrix *ConsensusKernelbank(iftFileSet *fs_seeds, char *inputdata_dir, int noutput_channels)
{
    int nkernels = fs_seeds->n;
    iftMatrix **kernels = (iftMatrix **)calloc(nkernels, sizeof(iftMatrix *));
    int ncols = 0;
    int nrows = 0;
    char filename[200];

    /* Load kernels from the training images */
    for (int i = 0; i < nkernels; i++)
    {
        char *basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
        sprintf(filename, "%s/%s-kernels.npy", inputdata_dir, basename);
        kernels[i] = iftReadMatrix(filename);
        ncols += kernels[i]->ncols;
        iftFree(basename);
    }
    nrows = kernels[0]->nrows;

    /* Copy all kernels into a single matrix */

    iftMatrix *Mkernels = iftCreateMatrix(ncols, nrows);
    int k = 0;
    for (int i = 0; i < nkernels; i++)
    {
        for (int col = 0; col < kernels[i]->ncols; col++, k++)
        {
            for (int row = 0; row < kernels[i]->nrows; row++)
            {
                iftMatrixElem(Mkernels, k, row) = iftMatrixElem(kernels[i], col, row);
            }
        }
    }

    for (int i = 0; i < nkernels; i++)
    {
        iftDestroyMatrix(&kernels[i]);
    }
    iftFree(kernels);

    /* Reduce the number of kernels into the desired number of output channels */

    iftMatrix *Rkernels = NULL;

    if (Mkernels->ncols <= noutput_channels)
    {
        Rkernels = iftCopyMatrix(Mkernels);
    }
    else
    {
        if (Mkernels->ncols > Mkernels->nrows)
        {
            Rkernels = SelectRelevantKernelsByPCA(Mkernels, noutput_channels);
        }
        else
        {
            Rkernels = SelectRelevantKernelsByKmeans(Mkernels, noutput_channels);
        }
    }

    iftDestroyMatrix(&Mkernels);

    return (Rkernels);
}

iftMatrix *ConsensusKernelbankGraph(iftFileSet *fs_seeds, char *inputdata_dir, int noutput_channels)
{
    int nkernels = fs_seeds->n;
    iftMatrix **kernels = (iftMatrix **)calloc(nkernels, sizeof(iftMatrix *));
    int ncols = 0;
    int nrows = 0;
    char filename[200];

    /* Load kernels from the training images */
    for (int i = 0; i < nkernels; i++)
    {
        char *basename = iftFilename(fs_seeds->files[i]->path, "-seeds_graph.txt");
        sprintf(filename, "%s/%s-kernels.npy", inputdata_dir, basename);
        kernels[i] = iftReadMatrix(filename);
        ncols += kernels[i]->ncols;
        iftFree(basename);
    }
    nrows = kernels[0]->nrows;

    /* Copy all kernels into a single matrix */

    iftMatrix *Mkernels = iftCreateMatrix(ncols, nrows);
    int k = 0;
    for (int i = 0; i < nkernels; i++)
    {
        for (int col = 0; col < kernels[i]->ncols; col++, k++)
        {
            for (int row = 0; row < kernels[i]->nrows; row++)
            {
                iftMatrixElem(Mkernels, k, row) = iftMatrixElem(kernels[i], col, row);
            }
        }
    }

    for (int i = 0; i < nkernels; i++)
    {
        iftDestroyMatrix(&kernels[i]);
    }
    iftFree(kernels);

    /* Reduce the number of kernels into the desired number of output channels */

    iftMatrix *Rkernels = NULL;

    if (Mkernels->ncols <= noutput_channels)
    {
        Rkernels = iftCopyMatrix(Mkernels);
    }
    else
    {
        if (Mkernels->ncols > Mkernels->nrows)
        {
            Rkernels = SelectRelevantKernelsByPCA(Mkernels, noutput_channels);
        }
        else
        {
            Rkernels = SelectRelevantKernelsByKmeans(Mkernels, noutput_channels);
        }
    }

    iftDestroyMatrix(&Mkernels);

    return (Rkernels);
}


void StatisticsFromAllSeeds(iftFileSet *fs_seeds, char *inputdata_dir, float *mean, float *stdev, float stdev_factor)
{
    int nseeds = 0, ninput_channels = 0;
    char *basename = NULL;
    char filename[200];
    iftMImage *input = NULL;

    /*
    Para cada arquivo (imagem):
        - lê o arquivo de sementes como iftMImage
        - transforma de iftMImage para iftLabeledSet
        - itera sobre as sementes, contando seu número e o somatório de suas features (cor)
    */
    for (int i = 0; i < fs_seeds->n; i++)
    {
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
        sprintf(filename, "tmp/%s.mimg", basename);
        input = iftReadMImage(filename);
        iftFree(basename);

        ninput_channels = input->m;
        iftLabeledSet *S = iftMReadSeeds(input, fs_seeds->files[i]->path);
        while (S != NULL)
        {
            int l;
            int p = iftRemoveLabeledSet(&S, &l);
            nseeds += 1;
            for (int b = 0; b < ninput_channels; b++)
            {
                mean[b] += input->val[p][b];
            }
        }
        iftDestroyMImage(&input);
    }

    // termina de calcular a média
    for (int b = 0; b < ninput_channels; b++)
    {
        mean[b] = mean[b] / nseeds;
    }

    /*
    Para cada arquivo (imagem):
        - lê o arquivo de sementes como iftMImage
        - transforma de iftMImage para iftLabeledSet
        - itera sobre as sementes, contando seu número e o desvio padrão de suas features (cor)
    */
    for (int i = 0; i < fs_seeds->n; i++)
    {
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
        sprintf(filename, "tmp/%s.mimg", basename);
        input = iftReadMImage(filename);
        iftFree(basename);

        ninput_channels = input->m;
        iftLabeledSet *S = iftMReadSeeds(input, fs_seeds->files[i]->path);
        while (S != NULL)
        {
            int l;
            int p = iftRemoveLabeledSet(&S, &l);
            for (int b = 0; b < ninput_channels; b++)
            {
                stdev[b] += (input->val[p][b] - mean[b]) * (input->val[p][b] - mean[b]);
            }
        }
        iftDestroyMImage(&input);
    }

    // termina de calcualr o desvio padrão
    for (int b = 0; b < ninput_channels; b++)
    {
        stdev[b] = sqrtf(stdev[b] / nseeds) + stdev_factor;
    }
}

void ReadMeanStdev(char *basepath, float *mean, float *stdev, int ninput_channels)
{
    char filename[2][200];
    FILE *fp[2];

    sprintf(filename[0], "%s-mean.txt", basepath);
    sprintf(filename[1], "%s-stdev.txt", basepath);
    fp[0] = fopen(filename[0], "r");
    fp[1] = fopen(filename[1], "r");
    for (int b = 0; b < ninput_channels; b++)
    {
        if (fscanf(fp[0], "%f", &mean[b]) != 1)
            ;
        if (fscanf(fp[1], "%f", &stdev[b]) != 1)
            ;
    }
    fclose(fp[0]);
    fclose(fp[1]);
}

void WriteMeanStdev(char *basepath, float *mean, float *stdev, int ninput_channels)
{
    char filename[2][200];
    FILE *fp[2];

    sprintf(filename[0], "%s-mean.txt", basepath);
    sprintf(filename[1], "%s-stdev.txt", basepath);
    fp[0] = fopen(filename[0], "w");
    fp[1] = fopen(filename[1], "w");
    for (int b = 0; b < ninput_channels; b++)
    {
        fprintf(fp[0], "%f ", mean[b]);
        fprintf(fp[1], "%f ", stdev[b]);
    }
    fclose(fp[0]);
    fclose(fp[1]);
}

float *ReadBias(char *basepath)
{
    char filename[200];
    FILE *fp;
    int number_of_kernels;

    sprintf(filename, "%s-bias.txt", basepath);
    fp = fopen(filename, "r");
    if (fscanf(fp, "%d", &number_of_kernels) != 1)
        ;
    float *bias = iftAllocFloatArray(number_of_kernels);
    for (int k = 0; k < number_of_kernels; k++)
    {
        if (fscanf(fp, "%f", &bias[k]) != 1)
            ;
    }
    fclose(fp);

    return (bias);
}

void WriteBias(char *basepath, float *bias, int number_of_kernels)
{
    char filename[200];
    FILE *fp;

    sprintf(filename, "%s-bias.txt", basepath);
    fp = fopen(filename, "w");
    fprintf(fp, "%d\n", number_of_kernels);
    for (int k = 0; k < number_of_kernels; k++)
    {
        fprintf(fp, "%f ", bias[k]);
    }
    fclose(fp);
}

iftImage *SubsampleImage(iftImage *img, int stride)
{
    iftImage *simg = iftCreateImage(ceilf(img->xsize / (float)stride), ceilf(img->ysize / (float)stride),
                                    iftMax(ceilf(img->zsize / (float)stride), 1));
    int q = 0;
    iftVoxel u;
    for (u.z = 0; u.z < img->zsize; u.z = u.z + stride)
    {
        for (u.y = 0; u.y < img->ysize; u.y = u.y + stride)
        {
            for (u.x = 0; u.x < img->xsize; u.x = u.x + stride)
            {
                int p = iftGetVoxelIndex(img, u);
                simg->val[q] = img->val[p];
                q++;
            }
        }
    }

    return (simg);
}

/* --------------------- public functions ----------------------------- */
/* Create a new empty graph (graph_orig == NULL) or a copy of graph_orig */
iftFLIMGraph *iftCreateFLIMGraph(int num_feats, iftFLIMGraph *graph_orig)
{
    iftFLIMGraph *graph;

    graph = (iftFLIMGraph *)malloc(sizeof(iftFLIMGraph));
    graph->num_feats = num_feats;

    if (graph_orig == NULL)
    {
        graph->num_nodes = 0;
        return graph;
    }

    graph->num_nodes = graph_orig->num_nodes;
    graph->nodes = (iftFLIMGraphNode *)calloc(graph->num_nodes, sizeof(iftFLIMGraphNode));

    for (unsigned long node_index = 0; node_index < graph->num_nodes; node_index++)
    {
        graph->nodes[node_index].numNeighbors = graph_orig->nodes[node_index].numNeighbors;
        graph->nodes[node_index].numPixels = graph_orig->nodes[node_index].numPixels;
        graph->nodes[node_index].index = node_index;

        graph->nodes[node_index].neighbors_list_head = NULL;
        graph->nodes[node_index].neighbors_list_tail = NULL;

        graph->nodes[node_index].feats = (double *)malloc(num_feats * sizeof(double));
        for (unsigned long b = 0; b < graph_orig->num_feats; b++)
            graph->nodes[node_index].feats[b] = graph_orig->nodes[node_index].feats[b];
        for (unsigned long b = graph_orig->num_feats; b < graph->num_feats; b++)
            graph->nodes[node_index].feats[b] = 0;

        iftFLIMGraphNodeList *neighbors = graph_orig->nodes[node_index].neighbors_list_head;
        while (neighbors != NULL)
        {
            // allocate a new neigbor for node graph->nodes[node_index]
            iftFLIMGraphNodeList *tmp = (iftFLIMGraphNodeList *)malloc(sizeof(iftFLIMGraphNodeList));

            tmp->next = NULL;
            tmp->node = &(graph->nodes[neighbors->node->index]);
            if (graph->nodes[node_index].numNeighbors == 0)
            {
                graph->nodes[node_index].neighbors_list_head = tmp;
                graph->nodes[node_index].neighbors_list_tail = tmp;
                graph->nodes[node_index].numNeighbors++;
            }
            else
            {
                graph->nodes[node_index].neighbors_list_tail->next = tmp;
                graph->nodes[node_index].neighbors_list_tail = tmp;
                graph->nodes[node_index].numNeighbors++;
            }
            neighbors = neighbors->next;
        }
    }
    return graph;
}

void isValidFLIMGraph(iftFLIMGraph *graph)
{
    if (graph == NULL)
    {
        iftError("Graph is NULL", "isValidFLIMGraph");
    }
    if (graph->num_nodes == 0)
    {
        iftError("Graph has no nodes", "isValidFLIMGraph");
    }
    if (graph->nodes == NULL)
    {
        iftError("Graph has no nodes", "isValidFLIMGraph");
    }
}

/*
    Create a graph with n vertices, each one with m features.
    The graph edges come from a reference graph
    and the node features come from a matrix with size [n][m].
*/
iftFLIMGraph *iftMatrixToFLIMGraph(iftMatrix *matrix, iftFLIMGraph *graph_ref)
{
    iftFLIMGraph *graph;

    assert(matrix->nrows == graph_ref->num_nodes);

    graph = (iftFLIMGraph *)malloc(sizeof(iftFLIMGraph));
    graph->num_feats = matrix->ncols;

    graph->num_nodes = graph_ref->num_nodes;
    graph->nodes = (iftFLIMGraphNode *)calloc(graph->num_nodes, sizeof(iftFLIMGraphNode));

    for (int node_index = 0; node_index < graph->num_nodes; node_index++)
    {
        graph->nodes[node_index].numNeighbors = 0;
        graph->nodes[node_index].numPixels = graph_ref->nodes[node_index].numPixels;
        graph->nodes[node_index].index = node_index;

        graph->nodes[node_index].neighbors_list_head = NULL;
        graph->nodes[node_index].neighbors_list_tail = NULL;

        graph->nodes[node_index].feats = (double *)malloc(matrix->ncols * sizeof(double));
        for (int b = 0; b < matrix->ncols; b++)
            graph->nodes[node_index].feats[b] = iftMatrixElem(matrix, b, node_index);

        iftFLIMGraphNodeList *neighbors = graph_ref->nodes[node_index].neighbors_list_head;
        while (neighbors != NULL)
        {
            // allocate a new neigbor for node graph->nodes[node_index]
            iftFLIMGraphNodeList *tmp = (iftFLIMGraphNodeList *)malloc(sizeof(iftFLIMGraphNodeList));

            tmp->next = NULL;
            tmp->node = &(graph->nodes[neighbors->node->index]);
            if (graph->nodes[node_index].numNeighbors == 0)
            {
                graph->nodes[node_index].neighbors_list_head = tmp;
                graph->nodes[node_index].neighbors_list_tail = tmp;
                graph->nodes[node_index].numNeighbors++;
            }
            else
            {
                graph->nodes[node_index].neighbors_list_tail->next = tmp;
                graph->nodes[node_index].neighbors_list_tail = tmp;
                graph->nodes[node_index].numNeighbors++;
            }
            neighbors = neighbors->next;
        }
    }
    isValidFLIMGraph(graph);
    return graph;
}


iftFLIMGraph *iftReadFLIMGraph(char *filename)
{
    iftDict *graph_dict = iftReadJson(filename);

    unsigned long num_nodes = strtoul(iftGetStrValFromDict("num_nodes", graph_dict), NULL, 0);
    unsigned long num_feats = strtoul(iftGetStrValFromDict("num_feats", graph_dict), NULL, 0);

    iftFLIMGraph *graph = (iftFLIMGraph *)calloc(1, sizeof(iftFLIMGraph));

    graph->num_nodes = num_nodes;
    graph->num_feats = num_feats;

    graph->nodes = (iftFLIMGraphNode *)calloc(graph->num_nodes, sizeof(iftFLIMGraphNode));

    for (unsigned long spx_index = 0; spx_index < graph->num_nodes; spx_index++)
    {
        graph->nodes[spx_index].numNeighbors = 0;
        graph->nodes[spx_index].neighbors_list_head = NULL;
        graph->nodes[spx_index].neighbors_list_tail = NULL;
        graph->nodes[spx_index].numPixels = 0;
        graph->nodes[spx_index].index = spx_index;
        graph->nodes[spx_index].feats = (double *)malloc(num_feats * sizeof(double));
        for (unsigned long feat_index = 0; feat_index < num_feats; feat_index++)
            graph->nodes[spx_index].feats[feat_index] = 0.0;
    }

    for (unsigned long l = 0; l < graph->num_nodes; l++)
    {
        char name[50];
        sprintf(name, "node%lu", l);
        iftDict *node_dict = iftGetDictFromDict(name, graph_dict);

        unsigned long node_index = strtoul(iftGetStrValFromDict("node_index", node_dict), NULL, 0);
        unsigned long num_pixels = strtoul(iftGetStrValFromDict("num_pixels", node_dict), NULL, 0);
        unsigned long num_neighbors = strtoul(iftGetStrValFromDict("num_neighbors", node_dict), NULL, 0);

        graph->nodes[node_index].numPixels = num_pixels;

        iftDblArray *feats = iftGetDblArrayFromDict("feats", node_dict);
        for (unsigned long feat_index = 0; feat_index < graph->num_feats; feat_index++)
        {
            graph->nodes[node_index].feats[feat_index] = feats->val[feat_index];
        }

        iftStrArray *input = iftGetStrArrayFromDict("neighbors", node_dict);

        unsigned long added_neighb = 0;
        for (unsigned long neighb = 0; neighb < num_neighbors; neighb++)
        {
            unsigned long neighb_val = strtoul(input->val[neighb], NULL, 0);
            if (node_index > neighb_val)
            {
                added_neighb++;
                iftInsertFLIMGraphNeighborNode(node_index, neighb_val, graph);
                iftInsertFLIMGraphNeighborNode(neighb_val, node_index, graph);
            }
        }

        // iftDestroyStrArray(&input);
        // iftDestroyDict(&node_dict);
    }
    iftDestroyDict(&graph_dict);
    isValidFLIMGraph(graph);
    return (graph);
}


iftImage *iftGraphToImage(iftFLIMGraph *graph, iftImage *labels, int Imax, int band)
{
    int p, b = band;

    iftImage *img = iftCreateImage(labels->xsize, labels->ysize, labels->zsize);
    double min = IFT_INFINITY_FLT, max = IFT_INFINITY_FLT_NEG;

    if ((band < 0) || (band >= (int)graph->num_feats))
        iftError("Invalid band", "iftGraphToImage");

    for (p = 0; p < graph->num_nodes; p++)
    {
        if (graph->nodes[p].feats[b] < min)
            min = graph->nodes[p].feats[b];

        if (graph->nodes[p].feats[b] > max)
            max = graph->nodes[p].feats[b];
    }

    double final_max = 0;
    double final_min = 255;

    if (max > min)
    {
        for (p = 0; p < img->n; p++)
        {
            img->val[p] = (int)((double)(Imax) * ((graph->nodes[labels->val[p] - 1].feats[b] - min) / (max - min)));

            if (img->val[p] < final_min)
                final_min = (double)(img->val[p]);
            if (img->val[p] > final_max)
                final_max = (double)(img->val[p]);
        }
    }

    img->dx = 1;
    img->dy = 1;
    img->dz = 1;

    return img;
}


iftMImage *iftGraphToMImage(iftFLIMGraph *graph, iftImage *labels)
{
    int p, b;

    iftMImage *img = iftCreateMImage(labels->xsize, labels->ysize, labels->zsize, graph->num_feats);

    for (b = 0; b < img->m; b++)
    {
        for (p = 0; p < img->n; p++)
        {
            img->val[p][b] = graph->nodes[labels->val[p] - 1].feats[b];
        }
    }

    img->dx = 1;
    img->dy = 1;
    img->dz = 1;

    return img;
}

/**
void iftgetArchVal3DIntVector(iftFLIMArch *arch, char *key, FILE *fp, int *vec){
    char tmp[255];

    char k;
    int id = 0;
    do{
        k = fgetc(fp);
        if(k == ',' || k == ']') {
            vec[id] = atoi(tmp);
            tmp[0] = '\0';
            id++;
        }
        else{
            if(k >= 48 && k <= 57){ // is a number in [0,9]
                tmp[strlen(tmp)] = k;
                tmp[strlen(tmp)+1] = '\0';
            }
        }
        id++;
    }while(id < 3);

    return;
}

void iftgetArchVal1Level(iftFLIMArch *arch, char *key, FILE *fp){
    char tmp[255];

    if(strcmp(key, "\"apply_intrinsic_atrous\"") == 0){
        fscanf(fp, "%s[^\n]", tmp); // [true,]  getValue
        while(tmp[strlen(tmp)] == ',' || tmp[strlen(tmp)] == '\n' || tmp[strlen(tmp)] == ' ')
            tmp[strlen(tmp)-1] = '\0';
        arch->apply_intrinsic_atrous = (strcmp(tmp, "true") == 0);
        return;
    }

    if(strcmp(key, "\"nlayers\"") == 0){
        fscanf(fp, "%s[^\n]", tmp);
        while(tmp[strlen(tmp)] == ',' || tmp[strlen(tmp)] == '\n' || tmp[strlen(tmp)] == ' ')
            tmp[strlen(tmp)-1] = '\0';
        arch->nlayers = atoi(tmp);
        return;
    }

    if(strcmp(key, "\"stdev_factor\"") == 0){
        fscanf(fp, "%s[^\n]", tmp);
        while(tmp[strlen(tmp)] == ',' || tmp[strlen(tmp)] == '\n' || tmp[strlen(tmp)] == ' ')
            tmp[strlen(tmp)-1] = '\0';
        arch->stdev_factor = atof(tmp);
        return;
    }
}

iftFLIMArch * iftReadFLIMArch2(char *filename){

    iftFLIMArch *arch = (iftFLIMArch *)calloc(1, sizeof(iftFLIMArch));
    char tmp[255];

    FILE *fp = fopen(filename, "r");
    int space_tabs = 8, ntabs = 0;

    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftReadFLIMArch", filename);

    fscanf(fp, "%s[^\n]", tmp); // [{\n ]

    fscanf(fp, "%s[^ ]", tmp); // ["apply_intrinsic_atrous":] getKey
    tmp[strlen(tmp)-1] = '\0';
    printf("key: %s \n", tmp);

    fscanf(fp, "%s[^ ]", tmp); // [true,]  getValue
    if(tmp[strlen(tmp)] == ',') tmp[strlen(tmp)-1] = '\0';
    printf("value: %s \n", tmp, strlen(tmp));

    exit(1);


    fscanf(fp, "        \"stdev_factor\": %f,\n", &arch->stdev_factor);

    //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
    fscanf(fp, "        \"nlayers\": %d,\n", &arch->nlayers);

    //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
    fscanf(fp, "        \"apply_intrinsic_at_each_layer\": %s,\n", tmp);
    printf("intrinsic_at_each_layer: %s \n", tmp);
    arch->apply_intrinsic_atrous = strcmp(tmp, "true") == 0 ? true : false;

    for (int l = 1; l < arch->nlayers+1; l++){
        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "%s[^\n]", tmp);
        printf("tmp ( \"layer_n\": { ): %s \n", tmp);

        //ntabs++;
        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "%s[^\n]", tmp);
        printf("tmp ( \"conv\": { ): %s \n", tmp);

        //ntabs++;
        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "                        \"kernel_size\": [%d, %d, %d],\n",
                    &(arch->layer[l-1].kernel_size[0]),
                    &(arch->layer[l-1].kernel_size[1]),
                    &(arch->layer[l-1].kernel_size[2]));

        printf("                        \"kernel_size\": [%d, %d, %d],\n",
                    (arch->layer[l-1].kernel_size[0]),
                    (arch->layer[l-1].kernel_size[1]),
                    (arch->layer[l-1].kernel_size[2]));

        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "                        \"nkernels_per_marker\": %d,\n", &(arch->layer[l-1].nkernels_per_marker));
        printf("                        \"nkernels_per_marker\": %d,\n", (arch->layer[l-1].nkernels_per_marker));

        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "                        \"dilation_rate\": [%d, %d, %d],\n",
                    &(arch->layer[l-1].dilation_rate[0]),
                    &(arch->layer[l-1].dilation_rate[1]),
                    &(arch->layer[l-1].dilation_rate[2]));
        printf("                        \"dilation_rate\": [%d, %d, %d],\n",
                    (arch->layer[l-1].dilation_rate[0]),
                    (arch->layer[l-1].dilation_rate[1]),
                    (arch->layer[l-1].dilation_rate[2]));

        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "                        \"nkernels_per_image\": %d,\n", &(arch->layer[l-1].nkernels_per_image));
        printf("                        \"nkernels_per_image\": %d,\n", (arch->layer[l-1].nkernels_per_image));

        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "                        \"noutput_channels\": %d\n", &(arch->layer[l-1].noutput_channels));
        printf("                        \"noutput_channels\": %d\n", (arch->layer[l-1].noutput_channels));

        //ntabs--;
        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "%s[^\n]", tmp);
        printf("%s \n", tmp);

        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "                \"relu\": %s,\n", tmp);
        printf("                \"relu\": %s,\n", tmp);
        arch->layer[l-1].relu = strcmp(tmp, "true") == 0 ? true : false;

        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "%s[^\n]", tmp);
        printf("%s \n", tmp);

        //ntabs++;
        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "                        \"type\": \"%s\",\n", arch->layer[l-1].pool_type);
        printf("                        \"type\": \"%s\",\n", arch->layer[l-1].pool_type);

        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "                        \"size\": [%d, %d, %d],\n",
                    &(arch->layer[l-1].pool_size[0]),
                    &(arch->layer[l-1].pool_size[1]),
                    &(arch->layer[l-1].pool_size[2]));
        printf("                        \"size\": [%d, %d, %d],\n",
                    (arch->layer[l-1].pool_size[0]),
                    (arch->layer[l-1].pool_size[1]),
                    (arch->layer[l-1].pool_size[2]));

        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "                        \"stride\": %d\n", &(arch->layer[l-1].pool_stride));
        printf("                        \"stride\": %d\n", (arch->layer[l-1].pool_stride));

        //ntabs--;
        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        fscanf(fp, "%s[^\n]", tmp);
        printf("%s \n", tmp);

        //ntabs--;
        //for(int i=0; i < space_tabs*ntabs; i++) fgetc(fp);
        if(l != arch->nlayers) {
            fscanf(fp, "%s[^\n]", tmp);
            printf("%s \n", tmp);
        }
    }
    fclose(fp);
    return arch;
}

void iftWriteFLIMArch2(iftFLIMArch *arch, char *filename)
{
    FILE *fp = fopen(filename, "wb");
    int space_tabs = 8, ntabs = 0;

    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteFLIMArch2", filename);

    fprintf(fp, "{\n");

    ntabs++;
    for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
    fprintf(fp, "\"stdev_factor\": %.6f,\n", arch->stdev_factor);

    for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
    fprintf(fp, "\"nlayers\": %d,\n", arch->nlayers);

    for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
    bool x = arch->apply_intrinsic_atrous;
    fprintf(fp, "\"apply_intrinsic_atrous\": %s,\n", x ? "true":"false");

    for (int l = 1; l < arch->nlayers+1; l++){
        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"layer%d\": {\n", l);

        ntabs++;
        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"conv\": {\n");

        ntabs++;
        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"kernel_size\": [");
        for(int i=0; i < 2; i++){
            fprintf(fp, "%d, ", arch->layer[l-1].kernel_size[i]);
        }
        fprintf(fp, "%d],\n", arch->layer[l-1].kernel_size[2]);

        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"nkernels_per_marker\": %d,\n", arch->layer[l-1].nkernels_per_marker);

        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"dilation_rate\": [");
        for(int i=0; i < 2; i++){
            fprintf(fp, "%d, ", arch->layer[l-1].dilation_rate[i]);
        }
        fprintf(fp, "%d],\n", arch->layer[l-1].dilation_rate[2]);

        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"nkernels_per_image\": %d,\n", arch->layer[l-1].nkernels_per_image);
        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"noutput_channels\": %d\n", arch->layer[l-1].noutput_channels);

        ntabs--;
        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "},\n");

        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        x = arch->layer[l-1].relu;
        fprintf(fp, "\"relu\": %s,\n", x ? "true":"false");

        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"pooling\": {\n");

        ntabs++;
        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"type\": \"%s\",\n", arch->layer[l-1].pool_type);

        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"size\": [");
        for(int i=0; i < 2; i++){
            fprintf(fp, "%d, ", arch->layer[l-1].pool_size[i]);
        }
        fprintf(fp, "%d],\n", arch->layer[l-1].pool_size[2]);

        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "\"stride\": %d\n", arch->layer[l-1].pool_stride);

        ntabs--;
        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        fprintf(fp, "}\n");

        ntabs--;
        for(int i=0; i < space_tabs*ntabs; i++) fprintf(fp, " ");
        if(l == arch->nlayers) fprintf(fp, "}\n");
        else fprintf(fp, "},\n");
    }
    fprintf(fp, "}");

    fclose(fp);
}
*/

void iftWriteFLIMGraph(iftFLIMGraph *graph, char *filename)
{

    FILE *fp = fopen(filename, "wb");

    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteFLIMGraph", filename);

    fprintf(fp, "{\n");

    for (int l = 0; l < graph->num_nodes; l++)
    {
        for (int i = 0; i < 4; i++)
            fprintf(fp, " ");
        fprintf(fp, "\"node%d\": {\n", l);

        for (int i = 0; i < 8; i++)
            fprintf(fp, " ");
        fprintf(fp, "\"feats\": [\n");
        for (int i = 0; i < graph->num_feats; i++)
        {
            for (int j = 0; j < 12; j++)
                fprintf(fp, " ");
            fprintf(fp, "%lf", graph->nodes[l].feats[i]);
            if (i != graph->num_feats - 1)
                fprintf(fp, ",");
            fprintf(fp, "\n");
        }
        for (int i = 0; i < 8; i++)
            fprintf(fp, " ");
        fprintf(fp, "],\n");

        for (int i = 0; i < 8; i++)
            fprintf(fp, " ");
        fprintf(fp, "\"neighbors\": [\n");
        iftFLIMGraphNodeList *node_list_ptr = graph->nodes[l].neighbors_list_head;
        for (unsigned long neighb = 0; neighb < graph->nodes[l].numNeighbors; neighb++)
        {
            for (int i = 0; i < 12; i++)
                fprintf(fp, " ");
            fprintf(fp, "\"%lu\"", node_list_ptr->node->index);
            node_list_ptr = node_list_ptr->next;
            if (neighb != graph->nodes[l].numNeighbors - 1)
                fprintf(fp, ",");
            fprintf(fp, "\n");
        }
        for (int i = 0; i < 8; i++)
            fprintf(fp, " ");
        fprintf(fp, "],\n");

        for (int i = 0; i < 8; i++)
            fprintf(fp, " ");
        fprintf(fp, "\"node_index\": \"%lu\",\n", graph->nodes[l].index);
        for (int i = 0; i < 8; i++)
            fprintf(fp, " ");
        fprintf(fp, "\"num_neighbors\": \"%lu\",\n", graph->nodes[l].numNeighbors);
        for (int i = 0; i < 8; i++)
            fprintf(fp, " ");
        fprintf(fp, "\"num_pixels\": \"%lu\"\n", graph->nodes[l].numPixels);

        for (int i = 0; i < 4; i++)
            fprintf(fp, " ");
        fprintf(fp, "},\n");
    }
    for (int i = 0; i < 4; i++)
        fprintf(fp, " ");
    fprintf(fp, "\"num_feats\": \"%lu\",\n", graph->num_feats);
    for (int i = 0; i < 4; i++)
        fprintf(fp, " ");
    fprintf(fp, "\"num_nodes\": \"%lu\"\n", graph->num_nodes);
    fprintf(fp, "}");

    fclose(fp);
}

void iftWriteFLIMGraph_2(iftFLIMGraph *graph, char *filename)
{
    iftDict *new_graph = iftCreateDict();
    iftDict *nodes[graph->num_nodes];

    /* Dictionary does not copy nested dicts, and therefore all nested dicts must be allocated first */

    char *num_nodes_char = (char *)malloc(50 * sizeof(char));
    sprintf(num_nodes_char, "%lu", graph->num_nodes);
    iftInsertIntoDict("num_nodes", num_nodes_char, new_graph);

    char *num_feats_char = (char *)malloc(50 * sizeof(char));
    sprintf(num_feats_char, "%lu", graph->num_feats);
    iftInsertIntoDict("num_feats", num_feats_char, new_graph);

    iftStrArray *node_index_char = iftCreateStrArray(graph->num_nodes);
    iftStrArray *numPixels_char = iftCreateStrArray(graph->num_nodes);
    iftStrArray *numNeighbors_char = iftCreateStrArray(graph->num_nodes);
    iftDblArray **feats = (iftDblArray **)malloc(graph->num_nodes * sizeof(iftDblArray *));
    iftStrArray **neighbors = (iftStrArray **)malloc(graph->num_nodes * sizeof(iftStrArray *));

    for (int l = 0; l < graph->num_nodes; l++)
    {
        nodes[l] = iftCreateDict();
        feats[l] = iftCreateDblArray(graph->num_feats);
    }
    for (int l = 0; l < graph->num_nodes; l++)
    {
        char name[100];
        sprintf(name, "node%d", l);

        for (int i = 0; i < graph->num_feats; i++)
            (feats[l])->val[i] = graph->nodes[l].feats[i];
        iftInsertIntoDict("feats", feats[l], nodes[l]);

        sprintf(node_index_char->val[l], "%lu", graph->nodes[l].index);
        iftInsertIntoDict("node_index", node_index_char->val[l], nodes[l]);

        sprintf(numPixels_char->val[l], "%lu", graph->nodes[l].numPixels);
        iftInsertIntoDict("num_pixels", numPixels_char->val[l], nodes[l]);

        sprintf(numNeighbors_char->val[l], "%lu", graph->nodes[l].numNeighbors);
        iftInsertIntoDict("num_neighbors", numNeighbors_char->val[l], nodes[l]);

        neighbors[l] = iftCreateStrArray(graph->nodes[l].numNeighbors);
        iftFLIMGraphNodeList *node_list_ptr = graph->nodes[l].neighbors_list_head;
        for (unsigned long neighb = 0; neighb < graph->nodes[l].numNeighbors; neighb++)
        {
            sprintf(neighbors[l]->val[neighb], "%lu", node_list_ptr->node->index);
            node_list_ptr = node_list_ptr->next;
        }
        iftInsertIntoDict("neighbors", neighbors[l], nodes[l]);

        char buf[255];
        sprintf(buf, "node%d", l);
        iftInsertIntoDict(name, nodes[l], new_graph);
    }

    iftWriteJson(new_graph, filename);
    // for (int l = 0; l < graph->num_nodes; l++) iftDestroyDict(&nodes[l]);

    iftDestroyDict(&new_graph);

    iftFree(num_nodes_char);
    iftFree(num_feats_char);
    iftDestroyStrArray(&node_index_char);
    iftDestroyStrArray(&numPixels_char);
    iftDestroyStrArray(&numNeighbors_char);
    for (int i = 0; i < graph->num_nodes; i++)
    {
        iftDestroyDblArray(&(feats[i]));
        iftDestroyStrArray(&(neighbors[i]));
        iftDestroyDict(&(nodes[i]));
    }
    iftFree(feats);
    iftFree(neighbors);
}

/* Exclui o vértice do grafo e suas arestas */
void iftDeleteVertexFLIMGraph(iftFLIMGraph *graph, int node_index)
{
    iftFLIMGraphNode *node = &graph->nodes[node_index];

    iftFLIMGraphNodeList *node_list = graph->nodes[node_index].neighbors_list_head;
    while (node_list != NULL)
    {
        // delete the node from the neighbor's list
        iftFLIMGraphNodeList *tmp = node_list->node->neighbors_list_head;
        iftFLIMGraphNodeList *prev = NULL;
        graph->nodes[node_list->node->index].numNeighbors--;
        while (tmp != NULL)
        {
            if (tmp->node->index == node_index)
            {
                if (prev == NULL)
                    tmp->node->neighbors_list_head = tmp->next;
                else
                    prev->next = tmp->next;
                tmp->next = NULL;
                tmp->node = NULL;
                iftFree(tmp);
                break;
            }
            prev = tmp;
            tmp = tmp->next;
        }

        // delete its neigbors list
        tmp = node_list;
        node_list = node_list->next;
        tmp->next = NULL;
        tmp->node = NULL;
        iftFree(tmp);
    }

    node->neighbors_list_head = NULL;
    node->neighbors_list_tail = NULL;
    iftFree(node->feats);
    node->feats = NULL;

    for (unsigned long i = node_index; i < graph->num_nodes - 1; i++)
    {
        graph->nodes[i] = graph->nodes[i + 1];
        // graph->nodes[i].index--;
    }

    graph->num_nodes--;
}

void iftDestroyFLIMGraph(iftFLIMGraph **graph)
{
    iftFLIMGraph *aux = *graph;

    if (aux != NULL)
    {
        for (int i = 0; i < aux->num_nodes; i++)
        {
            iftFLIMGraphNodeList *node_list = aux->nodes[i].neighbors_list_head;
            aux->nodes[i].neighbors_list_head = NULL;
            aux->nodes[i].neighbors_list_tail = NULL;

            while (node_list != NULL)
            {
                iftFLIMGraphNodeList *tmp = node_list;
                node_list = node_list->next;
                tmp->next = NULL;
                tmp->node = NULL;
                iftFree(tmp);
            }
            iftFree(aux->nodes[i].feats);
        }
        free(aux->nodes);
        iftFree(aux);
        *graph = NULL;
    }
}

void StatisticsFromAllSeedsGraph(iftFileSet *fs_seeds, char *dir, float *mean, float *stdev, float stdev_factor)
{
    int nseeds = 0, num_feats = 0;
    char *basename = NULL;
    char filename[200];
    iftFLIMGraph *graph = NULL;

    /*
    Para cada arquivo (imagem):
        - lê o grafo (para operar a média e desvio sobre as features)
        - lê o arquivo de sementes como iftLabeledSet
        - itera sobre as sementes, contando seu número e o somatório de suas features (cor)
    */
    for (int i = 0; i < fs_seeds->n; i++)
    {
        // lê as sementes
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds_graph.txt");

        // le o grafo
        sprintf(filename, "%s/%s.json", dir, basename);
        graph = iftReadFLIMGraph(filename);
        isValidFLIMGraph(graph);
        num_feats = graph->num_feats;

        sprintf(filename, "%s/%s.json", dir, basename);
        iftLabeledSet *S = iftReadSeedsGraph(fs_seeds->files[i]->path);

        iftFree(basename);
        while (S != NULL)
        {
            int l;
            int p = iftRemoveLabeledSet(&S, &l);
            nseeds += 1;
            for (int b = 0; b < num_feats; b++)
            {
                mean[b] += graph->nodes[p].feats[b];
            }
        }
        iftDestroyFLIMGraph(&graph);
    }

    // termina de calcular a média
    for (int b = 0; b < num_feats; b++)
    {
        mean[b] = mean[b] / nseeds;
    }

    /*
    Para cada arquivo (imagem):
        - lê o grafo (para operar a média e desvio sobre as features)
        - lê o arquivo de sementes como iftLabeledSet
        - itera sobre as sementes, contando seu número e o desvio padrão de suas features (cor)
    */
    for (int i = 0; i < fs_seeds->n; i++)
    {
        // lê as sementes
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds_graph.txt");

        // le o grafo
        sprintf(filename, "%s/%s.json", dir, basename);
        graph = iftReadFLIMGraph(filename);

        iftFree(basename);

        iftLabeledSet *S = iftReadSeedsGraph(fs_seeds->files[i]->path);
        while (S != NULL)
        {
            int l;
            int p = iftRemoveLabeledSet(&S, &l);
            for (int b = 0; b < num_feats; b++)
            {
                stdev[b] += (graph->nodes[p].feats[b] - mean[b]) * (graph->nodes[p].feats[b] - mean[b]);
            }
        }
        iftDestroyFLIMGraph(&graph);
    }

    // termina de calcualr o desvio padrão
    for (int b = 0; b < num_feats; b++)
    {
        stdev[b] = sqrtf(stdev[b] / nseeds) + stdev_factor;
    }
}


iftImage **iftGetActivations(iftFLIMGraph *graph, iftImage *labels)
{
    iftImage **activations = (iftImage **)malloc(graph->num_feats * sizeof(iftImage *));

    for (int i = 0; i < (int)(graph->num_feats); i++)
    {
        iftImage *band = iftGraphToImage(graph, labels, 255, i);

        if (iftIs3DImage(band))
        {
            activations[i] = iftGetXYSlice(band, band->zsize / 2);
        }
        else
        {
            activations[i] = iftCopyImage(band);
        }
        iftDestroyImage(&band);
    }

    return activations;
}


iftFLIMArch *iftReadFLIMArch(char *filename)
{
    iftDict *arch_dict = iftReadJson(filename);
    int nlayers = iftGetLongValFromDict("nlayers", arch_dict);
    float stdev_factor = iftGetDblValFromDict("stdev_factor", arch_dict);
    bool intrinsic_atrous = FALSE;
    if (iftDictContainKey("apply_intrinsic_atrous", arch_dict, NULL))
    {
        intrinsic_atrous = iftGetBoolValFromDict("apply_intrinsic_atrous", arch_dict);
    }
    iftFLIMArch *arch = (iftFLIMArch *)calloc(1, sizeof(iftFLIMArch));

    arch->nlayers = nlayers;
    arch->stdev_factor = stdev_factor;
    arch->layer = (iftFLIMLayer *)calloc(arch->nlayers, sizeof(iftFLIMLayer));
    arch->apply_intrinsic_atrous = intrinsic_atrous;

    for (int l = 0; l < arch->nlayers; l++)
    {

        char name[50];
        sprintf(name, "layer%d", l + 1);
        iftDict *layer_dict = iftGetDictFromDict(name, arch_dict);
        iftDict *conv_dict = iftGetDictFromDict("conv", layer_dict);
        iftIntArray *input = iftGetIntArrayFromDict("kernel_size", conv_dict);
        for (int i = 0; i < 3; i++)
            arch->layer[l].kernel_size[i] = input->val[i];
        iftFree(input);
        input = iftGetIntArrayFromDict("dilation_rate", conv_dict);
        for (int i = 0; i < 3; i++)
            arch->layer[l].dilation_rate[i] = input->val[i];
        iftFree(input);
        arch->layer[l].nkernels_per_image = iftGetLongValFromDict("nkernels_per_image", conv_dict);
        arch->layer[l].nkernels_per_marker = iftGetLongValFromDict("nkernels_per_marker", conv_dict);
        arch->layer[l].noutput_channels = iftGetLongValFromDict("noutput_channels", conv_dict);
        arch->layer[l].relu = iftGetBoolValFromDict("relu", layer_dict);
        iftDict *pool_dict = iftGetDictFromDict("pooling", layer_dict);
        arch->layer[l].pool_type = iftGetStrValFromDict("type", pool_dict);
        input = iftGetIntArrayFromDict("size", pool_dict);
        for (int i = 0; i < 3; i++)
            arch->layer[l].pool_size[i] = input->val[i];
        iftFree(input);
        arch->layer[l].pool_stride = iftGetLongValFromDict("stride", pool_dict);

        iftDestroyDict(&conv_dict);
        iftDestroyDict(&pool_dict);

        iftDestroyDict(&layer_dict);
    }

    iftDestroyDict(&arch_dict);

    return (arch);
}

void iftWriteFLIMArch(iftFLIMArch *arch, char *filename)
{
    iftDict *new_arch = iftCreateDict();

    /* Dictionary does not copy nested dicts, and therefore all nested dicts must be allocated first */
    iftDict *layer[arch->nlayers];
    iftDict *conv[arch->nlayers];
    iftDict *pooling[arch->nlayers];

    for (int l = 0; l < arch->nlayers; l++)
    {
        layer[l] = iftCreateDict();
        conv[l] = iftCreateDict();
        pooling[l] = iftCreateDict();
    }

    iftInsertIntoDict("stdev_factor", arch->stdev_factor, new_arch);
    iftInsertIntoDict("nlayers", arch->nlayers, new_arch);
    iftInsertIntoDict("apply_intrinsic_atrous", arch->apply_intrinsic_atrous, new_arch);

    for (int l = 0; l < arch->nlayers; l++)
    {

        char layername[100];
        sprintf(layername, "layer%d", l + 1);

        iftIntArray *kernel_size = iftCreateIntArray(3);
        for (int i = 0; i < 3; i++)
            kernel_size->val[i] = arch->layer[l].kernel_size[i];
        long nkernels_per_marker = arch->layer[l].nkernels_per_marker;
        iftIntArray *dilation_rate = iftCreateIntArray(3);
        for (int i = 0; i < 3; i++)
            dilation_rate->val[i] = arch->layer[l].dilation_rate[i];
        long nkernels_per_image = arch->layer[l].nkernels_per_image;
        long noutput_channels = arch->layer[l].noutput_channels;
        iftInsertIntoDict("kernel_size", kernel_size, conv[l]);
        iftInsertIntoDict("nkernels_per_marker", nkernels_per_marker, conv[l]);
        iftInsertIntoDict("dilation_rate", dilation_rate, conv[l]);
        iftInsertIntoDict("nkernels_per_image", nkernels_per_image, conv[l]);
        iftInsertIntoDict("noutput_channels", noutput_channels, conv[l]);
        iftInsertIntoDict("conv", conv[l], layer[l]);

        bool relu = TRUE;
        iftInsertIntoDict("relu", relu, layer[l]);

        char *pool_type = arch->layer[l].pool_type;
        iftIntArray *pool_size = iftCreateIntArray(3);
        for (int i = 0; i < 3; i++)
            pool_size->val[i] = arch->layer[l].pool_size[i];
        long pool_stride = arch->layer[l].pool_stride;
        iftInsertIntoDict("type", pool_type, pooling[l]);
        iftInsertIntoDict("size", pool_size, pooling[l]);
        iftInsertIntoDict("stride", pool_stride, pooling[l]);
        iftInsertIntoDict("pooling", pooling[l], layer[l]);

        iftInsertIntoDict(layername, layer[l], new_arch);
    }

    iftWriteJson(new_arch, filename);
    for (int l = 0; l < arch->nlayers; l++)
    {
        iftDestroyDict(&layer[l]);
        iftDestroyDict(&conv[l]);
        iftDestroyDict(&pooling[l]);
    }
    iftDestroyDict(&new_arch);
}

void iftDestroyFLIMArch(iftFLIMArch **arch)
{
    iftFLIMArch *aux = *arch;

    if (aux != NULL)
    {
        for (int l = 0; l < aux->nlayers; l++)
        {
            iftFree(aux->layer[l].pool_type);
        }
        iftFree(aux->layer);
        iftFree(aux);
        *arch = NULL;
    }
}

void iftFLIMLearnModel(char *orig_dir, char *markers_dir, char *param_dir, iftFLIMArch *arch)
{

    /* Set input parameters */

    iftMakeDir("tmp");

    iftFileSet *fs_orig = iftLoadFileSetFromDirOrCSV(orig_dir, 1, 1);
    iftFileSet *fs_seeds = iftLoadFileSetFromDirBySuffix(markers_dir, "-seeds.txt", 1);
    iftAdjRel *A = NULL;
    iftMImage **output = NULL;
    iftMImage *input = NULL;
    int ninput_channels = 0;
    int nimages = fs_seeds->n;
    int atrous_factor = 1;
    char *basename = NULL;
    char filename[200], ext[10];

    sprintf(ext, "%s", iftFileExt(fs_orig->files[0]->path));

    /* Generate input layer */
    for (int i = 0; i < nimages; i++)
    {
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
        sprintf(filename, "%s/%s%s", orig_dir, basename, ext);
        if (strcmp(ext, ".mimg") == 0)
        {
            input = iftReadMImage(filename);
        }
        else
        {
            input = ReadInputMImage(filename);
        }
        sprintf(filename, "tmp/%s.mimg", basename);
        iftWriteMImage(input, filename);
        iftFree(basename);
        iftDestroyMImage(&input);
    }

    /* For each layer do */

    for (int l = 0; l < arch->nlayers; l++)
    {
        basename = iftFilename(fs_seeds->files[0]->path, "-seeds.txt");
        sprintf(filename, "tmp/%s.mimg", basename);
        input = iftReadMImage(filename);
        if (l == 0)
        {
            A = iftFLIMAdjRelFromKernel(arch->layer[l], iftIs3DMImage(input));
        }
        else
        {
            if (arch->apply_intrinsic_atrous)
            {
                atrous_factor *= arch->layer[l - 1].pool_stride;
                // printf("Updating atrous factor\n");
                // fflush(stdout);
            }
            A = iftFLIMAdaptiveAdjRelFromKernel(arch->layer[l], atrous_factor, iftIs3DMImage(input));
        }
        ninput_channels = input->m;
        iftDestroyMImage(&input);

        /* Learn and save marker-based normalization parameters and kernels from each training image */

        float *mean = iftAllocFloatArray(ninput_channels);
        float *stdev = iftAllocFloatArray(ninput_channels);
        StatisticsFromAllSeeds(fs_seeds, "tmp", mean, stdev, arch->stdev_factor);
        // fflush(stdout);

        for (int i = 0; i < nimages; i++)
        {

            basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
            sprintf(filename, "tmp/%s.mimg", basename);
            input = iftReadMImage(filename);

            // printf("Processing file %s: %d of %d files\r", basename, i + 1, nimages);

            /* Apply marker-based image normalization */

            NormalizeImageByZScore(input, mean, stdev);
            sprintf(filename, "tmp/%s-norm.mimg", basename);
            iftWriteMImage(input, filename);

            /* Assign a distinct label to each marker */

            iftLabeledSet *S = iftMReadSeeds(input, fs_seeds->files[i]->path);
            iftLabeledSet *M = LabelMarkers(input, S);
            iftDestroyLabeledSet(&S);

            /* Learn and save kernel bank */

            int nsamples = iftLabeledSetSize(M);
            iftMatrix *kernelbank = LearnKernelBank(input, M, A, nsamples, arch->layer[l].nkernels_per_image, arch->layer[l].nkernels_per_marker);
            iftDestroyLabeledSet(&M);
            sprintf(filename, "tmp/%s-kernels.npy", basename);
            iftWriteMatrix(kernelbank, filename);
            iftDestroyMatrix(&kernelbank);
            iftDestroyMImage(&input);
            iftFree(basename);

            // fflush(stdout);
        }

        /* Create a consensus layer (i.e., merge kernels) and save
       final kernel bank for layer l */

        iftMatrix *kernelbank = ConsensusKernelbank(fs_seeds, "tmp", arch->layer[l].noutput_channels);

        /* BIAS: read bias array */
        char *use_bias = getenv("USE_BIAS");
        float *bias = NULL;
        if (use_bias != NULL)
        {
            bias = iftAllocFloatArray(kernelbank->ncols);
            for (int col = 0; col < kernelbank->ncols; col++)
            {
                int row = 0;
                for (int adj = 0; adj < A->n; adj++)
                {
                    for (int ch = 0; ch < ninput_channels; ch++)
                    {
                        iftMatrixElem(kernelbank, col, row) =
                            iftMatrixElem(kernelbank, col, row) / stdev[ch];
                        bias[col] -= (mean[ch] * iftMatrixElem(kernelbank, col, row));
                        row++;
                    }
                }
            }

            sprintf(filename, "%s/conv%d", param_dir, l + 1);
            WriteBias(filename, bias, kernelbank->ncols);
            iftFree(bias);
        }

        // updating number of kernels on the architecture
        arch->layer[l].noutput_channels = kernelbank->ncols;
        sprintf(filename, "%s/conv%d-kernels.npy", param_dir, l + 1);
        iftWriteMatrix(kernelbank, filename);
        iftDestroyMatrix(&kernelbank);
        sprintf(filename, "%s/conv%d", param_dir, l + 1);
        WriteMeanStdev(filename, mean, stdev, ninput_channels);
        iftFree(mean);
        iftFree(stdev);
        iftDestroyAdjRel(&A);

        /* Apply convolutional layer using the consensus kernel bank and the statistics from all markers */
        int pool_stride = arch->layer[l].pool_stride;
        arch->layer[l].pool_stride = 1;
        for (int i = 0; i < nimages; i++)
        {
            basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
            sprintf(filename, "tmp/%s.mimg", basename);
            input = iftReadMImage(filename);
            iftFree(basename);
            output = FLIMConvolutionalLayer(&input, 1, NULL, arch, l, l + 1, atrous_factor, param_dir);
            iftDestroyMImage(&input);

            iftWriteMImage(output[0], filename);
            iftDestroyMImage(&output[0]);
            iftFree(output);
        }
        arch->layer[l].pool_stride = pool_stride;
    }

    iftRemoveDir("tmp");
    iftDestroyFileSet(&fs_seeds);
    iftDestroyFileSet(&fs_orig);
}

void iftFLIMLearnModelPCA(char *orig_dir, char *markers_dir, char *param_dir, iftFLIMArch *arch)
{

    /* Set input parameters */

    iftMakeDir("tmp");

    iftFileSet *fs_orig = iftLoadFileSetFromDirOrCSV(orig_dir, 1, 1);
    iftFileSet *fs_seeds = iftLoadFileSetFromDirBySuffix(markers_dir, "-seeds.txt", 1);
    iftAdjRel *A = NULL;
    iftMImage **output = NULL;
    iftMImage *input = NULL;
    int ninput_channels = 0;
    int nimages = fs_seeds->n;
    int atrous_factor = 1;
    char *basename = NULL;
    char filename[200], ext[10];

    sprintf(ext, "%s", iftFileExt(fs_orig->files[0]->path));

    /* Generate input layer */
    for (int i = 0; i < nimages; i++)
    {
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
        sprintf(filename, "%s/%s%s", orig_dir, basename, ext);
        if (strcmp(ext, ".mimg") == 0)
        {
            input = iftReadMImage(filename);
        }
        else
        {
            input = ReadInputMImage(filename);
        }
        sprintf(filename, "tmp/%s.mimg", basename);
        iftWriteMImage(input, filename);
        iftFree(basename);
        iftDestroyMImage(&input);
    }

    /* For each layer do */

    for (int l = 0; l < arch->nlayers; l++)
    {
        basename = iftFilename(fs_seeds->files[0]->path, "-seeds.txt");
        sprintf(filename, "tmp/%s.mimg", basename);
        input = iftReadMImage(filename);
        if (l == 0)
        {
            A = iftFLIMAdjRelFromKernel(arch->layer[l], iftIs3DMImage(input));
        }
        else
        {
            if (arch->apply_intrinsic_atrous)
            {
                atrous_factor *= arch->layer[l - 1].pool_stride;
                // printf("Updating atrous factor\n");
                // fflush(stdout);
            }
            A = iftFLIMAdaptiveAdjRelFromKernel(arch->layer[l], atrous_factor, iftIs3DMImage(input));
        }
        ninput_channels = input->m;
        iftDestroyMImage(&input);

        /* Compute mean and standard deviation for marker-based
       normalization */

        float *mean = iftAllocFloatArray(ninput_channels);
        float *stdev = iftAllocFloatArray(ninput_channels);
        StatisticsFromAllSeeds(fs_seeds, "tmp", mean, stdev, arch->stdev_factor);
        /* Create and normalize the patch dataset from all seeds
           (using all training images) */

        iftDataSet *Z = iftDataSetFromAllSeeds(markers_dir, "tmp", A);
        iftNormalizeDataSetByZScoreInPlace(Z, NULL, arch->stdev_factor);

        /* Compute kernelbank for the current layer */

        int number_of_kernels = arch->layer[l].noutput_channels;
        iftMatrix *kernelbank = KernelBankByPCA(Z, &number_of_kernels);
        arch->layer[l].noutput_channels = number_of_kernels;

        // printf("\n Number of kernels=%d\n", arch->layer[l].noutput_channels);
        // fflush(stdout);

        /* BIAS: read bias array */
        char *use_bias = getenv("USE_BIAS");
        float *bias = NULL;
        if (use_bias != NULL)
        {
            bias = iftAllocFloatArray(kernelbank->ncols);
            for (int col = 0; col < kernelbank->ncols; col++)
            {
                int row = 0;
                for (int adj = 0; adj < A->n; adj++)
                {
                    for (int ch = 0; ch < ninput_channels; ch++)
                    {
                        iftMatrixElem(kernelbank, col, row) =
                            iftMatrixElem(kernelbank, col, row) / stdev[ch];
                        bias[col] -= (mean[ch] * iftMatrixElem(kernelbank, col, row));
                        row++;
                    }
                }
            }
        }

        sprintf(filename, "%s/conv%d", param_dir, l + 1);
        WriteMeanStdev(filename, mean, stdev, ninput_channels);
        iftFree(mean);
        iftFree(stdev);

        /* BIAS: write and free bias array */
        if (use_bias != NULL)
        {
            WriteBias(filename, bias, kernelbank->ncols);
            iftFree(bias);
        }

        iftDestroyAdjRel(&A);

        sprintf(filename, "%s/conv%d-kernels.npy", param_dir, l + 1);
        iftWriteMatrix(kernelbank, filename);
        iftDestroyMatrix(&kernelbank);
        iftDestroyDataSet(&Z);

        int pool_stride = arch->layer[l].pool_stride;
        arch->layer[l].pool_stride = 1;
        for (int i = 0; i < nimages; i++)
        {
            basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
            sprintf(filename, "tmp/%s.mimg", basename);
            input = iftReadMImage(filename);
            iftFree(basename);
            output = FLIMConvolutionalLayer(&input, 1, NULL, arch, l, l + 1, atrous_factor, param_dir);
            iftDestroyMImage(&input);

            iftWriteMImage(output[0], filename);
            iftDestroyMImage(&output[0]);
            iftFree(output);
        }
        arch->layer[l].pool_stride = pool_stride;
    }

    iftRemoveDir("tmp");
    iftDestroyFileSet(&fs_seeds);
    iftDestroyFileSet(&fs_orig);
}

void iftFLIMGraphLearnLayer(char *activ_dir, char *markers_dir, char *param_dir, int layer_index, iftFLIMArch *arch, char *output_dir)
{

    /* Set input parameters */
    iftMakeDir("tmp");

    iftFileSet *fs_activ = iftLoadFileSetFromDirOrCSV(activ_dir, 1, 1); // original images (for layer 1)
    iftFileSet *fs_seeds = iftLoadFileSetFromDirBySuffix(markers_dir, "-seeds_graph.txt", 1);

    int ninput_channels = 0;
    int nimages = fs_seeds->n;
    int atrous_factor = 1;
    char *basename = NULL;
    char filename[200], ext[10];

    sprintf(ext, "%s", iftFileExt(fs_activ->files[0]->path));

    /* Generate input layer */

    // read activation files and write in the temp directory
    for (int i = 0; i < nimages; i++)
    {
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds_graph.txt");

        /*************** write activation graph **************/
        sprintf(filename, "%s/%s%s", activ_dir, basename, ext);
        iftFLIMGraph *graph = iftReadFLIMGraph(filename);
        isValidFLIMGraph(graph);

        sprintf(filename, "tmp/%s.json", basename);
        iftWriteFLIMGraph(graph, filename);
        if (i == 0)
            ninput_channels = graph->num_feats;

        iftFree(basename);
        iftDestroyFLIMGraph(&graph);
    }

    /* For a specific layer do */
    /* Learn and save marker-based normalization parameters and kernels from each training image */

    // compute mean and std of the graph node seeds
    float *mean_graph = iftAllocFloatArray(ninput_channels);
    float *stdev_graph = iftAllocFloatArray(ninput_channels);

    // printf("Computing mean and std of the graph node seeds\r");
    StatisticsFromAllSeedsGraph(fs_seeds, "tmp", mean_graph, stdev_graph, arch->stdev_factor);
    // fflush(stdout);

    /*
    Para cada arquivo (grafo):
        - lê o grafo
        - lê o arquivo de sementes como iftLabeledSet: S
        - Rotula conjuntos conexos de sementes (conforme o grafo)
        - Normaliza todas as sementes conforme sua média e desvio padrão (de todas as sementes do grafo)
        - Escreve o grafo normalizado em "tmp/%s-norm.json"
        - calcula os kernels
        - escreve o kernel bank em "tmp/%s-kernels.npy"
    */
    for (int i = 0; i < nimages; i++)
    {
        // lê o arquivo de ativação de "tmp/%s.json" como iftFLIMGraph: graph
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds_graph.txt");

        // printf("Processing file %s: %d of %d files\r", basename, i + 1, nimages);
        // fflush(stdout);

        iftLabeledSet *S;
        iftLabeledSet *M;

        int nsamples;
        iftMatrix *kernelbank;

        /************ graph nodes as seeds ***************/
        sprintf(filename, "tmp/%s.json", basename);
        iftFLIMGraph *graph = iftReadFLIMGraph(filename);
        isValidFLIMGraph(graph);

        S = iftReadSeedsGraph(fs_seeds->files[i]->path);

        // LabelMarkersGraph : Dado um conjunto S de sementes em uma imagem img,
        // retorna o mesmo conjunto com rótulos iguais para sementes vizinhas.
        M = LabelMarkersGraph(graph, S);

        iftDestroyLabeledSet(&S);

        /* Apply marker-based image normalization */
        NormalizeGraphByZScore(graph, mean_graph, stdev_graph);
        sprintf(filename, "tmp/%s-norm.json", basename);

        iftWriteFLIMGraph(graph, filename);
        nsamples = iftLabeledSetSize(M); // número de conjuntos conexos de sementes

        // gera um conjunto de kernels,
        // em que cada kernel corresponde a uma semente
        // com suas N sementes vizinhas mais próximas.
        kernelbank = LearnKernelBankGraph(graph, M, arch->layer[0].kernel_size, nsamples, arch->layer[0].nkernels_per_image, arch->layer[0].nkernels_per_marker);

        iftDestroyFLIMGraph(&graph);
        iftDestroyLabeledSet(&M);

        sprintf(filename, "tmp/%s-kernels.npy", basename);
        iftWriteMatrix(kernelbank, filename);

        iftDestroyMatrix(&kernelbank);

        iftFree(basename);
    }

    /* Create a consensus layer (i.e., merge kernels) and save final kernel bank for layer l */
    // concatena as colunas dos kernels em um só kernel
    // e a quantidade de linhas é a mesma dos kernels anteriores
    // Caso a quantidade de colunas seja menor ou igual a quantidade de canais de saída da layer, apenas retorna o "merged kernel"
    // Caso contrário, se a quantidade de colunas for maior do que a de linhas, seleciona kernels com PCA. Caso contrário, com Kmeans.
    iftMatrix *kernelbank = ConsensusKernelbankGraph(fs_seeds, "tmp", arch->layer[0].noutput_channels);

    /* BIAS: read bias array */
    char *use_bias = getenv("USE_BIAS");
    float *bias = NULL;

    /************ graph nodes as seeds ***************/
    if (use_bias != NULL)
    {
        bias = iftAllocFloatArray(kernelbank->ncols); // inicia as posições com zero (calloc)

        for (int col = 0; col < kernelbank->ncols; col++)
        {
            for (int row = 0; row < kernelbank->nrows;)
            {
                for (int ch = 0; ch < ninput_channels; ch++)
                {
                    iftMatrixElem(kernelbank, col, row) =
                        iftMatrixElem(kernelbank, col, row) / stdev_graph[ch];
                    bias[col] -= (mean_graph[ch] * iftMatrixElem(kernelbank, col, row));
                    row++;
                }
            }
        }

        sprintf(filename, "%s/conv%d", param_dir, layer_index);
        WriteBias(filename, bias, kernelbank->ncols);
        iftFree(bias);
    }

    sprintf(filename, "%s/conv%d-kernels.npy", param_dir, layer_index);
    iftWriteMatrix(kernelbank, filename);
    iftDestroyMatrix(&kernelbank);
    sprintf(filename, "%s/conv%d", param_dir, layer_index);
    WriteMeanStdev(filename, mean_graph, stdev_graph, ninput_channels);
    iftFree(mean_graph);
    iftFree(stdev_graph);

    int pool_stride;
    /************ graph nodes as seeds ***************/
    /* Apply convolutional layer using the consensus kernel bank and the statistics from all markers */
    pool_stride = arch->layer[0].pool_stride;
    arch->layer[0].pool_stride = 1;
    for (int i = 0; i < nimages; i++)
    {
        iftFLIMGraph **output;
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds_graph.txt");
        sprintf(filename, "tmp/%s.json", basename);

        iftFLIMGraph *graph = iftReadFLIMGraph(filename);
        isValidFLIMGraph(graph);

        output = FLIMGraphConvolutionalLayer(&graph, 1, NULL, arch, 0, layer_index, atrous_factor, param_dir);
        isValidFLIMGraph(output[0]);
        iftDestroyFLIMGraph(&graph);

        sprintf(filename, "%s/%s.json", output_dir, basename);
        iftWriteFLIMGraph(output[0], filename);
        iftDestroyFLIMGraph(&(output[0]));
        iftFree(output);
        iftFree(basename);
    }
    arch->layer[0].pool_stride = pool_stride;

    // updating atrous factor with pooling of current layer
    if (arch->apply_intrinsic_atrous)
    {
        atrous_factor *= arch->layer[0].pool_stride;
    }
    // writing atrous factor on file
    sprintf(filename, "%s/intrinsic_atrous.txt", param_dir);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%d", atrous_factor);
    fclose(fp);

    iftRemoveDir("tmp");
    iftDestroyFileSet(&fs_seeds);
    iftDestroyFileSet(&fs_activ);
}


void iftFLIMLearnLayer(char *activ_dir, char *markers_dir, char *param_dir, int layer_index, iftFLIMArch *arch, char *output_dir)
{
    /*
        activ_dir = /tmp/tmpdir_FLIM/train/ (layer 1) or /tmp/tmpdir_FLIM/model/layerX/ (layer X+1)
        markers_dir = this->seeds_dir.toUtf8().data():/tmp/tmpdir_FLIM/seeds/
        param_dir = this->param_dir.toUtf8().data():/tmp/tmpdir_FLIM/model//
        layer_index = layer:1
        output_dir = layer_output.toUtf8().data(): /tmp/tmpdir_FLIM/model///layer1/
    */

    /* Set input parameters */
    iftMakeDir("tmp");

    iftFileSet *fs_activ = iftLoadFileSetFromDirOrCSV(activ_dir, 1, 1); // original images (for layer 1)
    iftFileSet *fs_seeds = iftLoadFileSetFromDirBySuffix(markers_dir, "-seeds.txt", 1);
    iftAdjRel *A = NULL;
    iftMImage **output = NULL;
    iftMImage *input = NULL;
    int ninput_channels = 0;
    int nimages = fs_seeds->n;
    int atrous_factor = 1;
    char *basename = NULL;
    char filename[200], ext[10];

    sprintf(ext, "%s", iftFileExt(fs_activ->files[0]->path));

    /* Generate input layer */

    // read seed file and write in a temp directory
    for (int i = 0; i < nimages; i++)
    {
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
        sprintf(filename, "%s/%s%s", activ_dir, basename, ext);
        if (strcmp(ext, ".mimg") == 0)
        {
            input = iftReadMImage(filename);
        }
        else
        {
            input = ReadInputMImage(filename);
        }
        sprintf(filename, "tmp/%s.mimg", basename);
        iftWriteMImage(input, filename);
        iftFree(basename);
        iftDestroyMImage(&input);
    }

    /* For a specific layer do */

    basename = iftFilename(fs_seeds->files[0]->path, "-seeds.txt");
    sprintf(filename, "tmp/%s.mimg", basename);
    input = iftReadMImage(filename);
    if (layer_index == 1)
    {
        A = iftFLIMAdjRelFromKernel(arch->layer[0], iftIs3DMImage(input));
    }
    else
    {
        // reading atrous factor already updated from previous layers
        sprintf(filename, "%s/intrinsic_atrous.txt", param_dir);
        char *data = iftReadFileAsString(filename);
        atrous_factor = atoi(data);
        iftFree(data);
        A = iftFLIMAdaptiveAdjRelFromKernel(arch->layer[0], atrous_factor, iftIs3DMImage(input));
    }
    ninput_channels = input->m;
    iftDestroyMImage(&input);

    /* Learn and save marker-based normalization parameters and kernels from each training image */

    float *mean = iftAllocFloatArray(ninput_channels);
    float *stdev = iftAllocFloatArray(ninput_channels);
    StatisticsFromAllSeeds(fs_seeds, "tmp", mean, stdev, arch->stdev_factor);

    /*
    Para cada arquivo (imagem):
        - lê o arquivo de sementes de "tmp/%s.mimg" como iftMImage
        - transforma de iftMImage para iftLabeledSet
        - Rotula conjuntos conexos de sementes
        - Normaliza todas as sementes conforme sua média e desvio padrão (de todas as sementes da imagem)
        - Escreve a imagem normalizada em "tmp/%s-norm.mimg"
        - calcula os kernels
        - escreve o kernel bank em "tmp/%s-kernels.npy"
    */
    for (int i = 0; i < nimages; i++)
    {

        basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
        sprintf(filename, "tmp/%s.mimg", basename);
        input = iftReadMImage(filename);

        // printf("Processing file %s: %d of %d files\r", basename, i + 1, nimages);
        // fflush(stdout);

        /* Learn mean and standard deviation */

        iftLabeledSet *S = iftMReadSeeds(input, fs_seeds->files[i]->path);
        /*
            LabelMarkers : Dado um conjunto S de sementes em uma imagem img,
            retorna o mesmo conjunto com rótulos iguais para sementes vizinhas.
        */
        iftLabeledSet *M = LabelMarkers(input, S);
        iftDestroyLabeledSet(&S);

        /* Apply marker-based image normalization */
        NormalizeImageByZScore(input, mean, stdev);
        sprintf(filename, "tmp/%s-norm.mimg", basename);
        iftWriteMImage(input, filename);

        /* Learn and save kernel bank */

        int nsamples = iftLabeledSetSize(M); // número de conjuntos conexos de sementes
        iftMatrix *kernelbank = LearnKernelBank(input, M, A, nsamples, arch->layer[0].nkernels_per_image, arch->layer[0].nkernels_per_marker);

        iftDestroyMImage(&input);
        iftDestroyLabeledSet(&M);
        sprintf(filename, "tmp/%s-kernels.npy", basename);
        iftWriteMatrix(kernelbank, filename);
        iftDestroyMatrix(&kernelbank);
        iftFree(basename);
    }

    /* Create a consensus layer (i.e., merge kernels) and save final kernel bank for layer l */

    iftMatrix *kernelbank = ConsensusKernelbank(fs_seeds, "tmp", arch->layer[0].noutput_channels);

    /* BIAS: read bias array */
    char *use_bias = getenv("USE_BIAS");
    float *bias = NULL;
    if (use_bias != NULL)
    {
        bias = iftAllocFloatArray(kernelbank->ncols);
        for (int col = 0; col < kernelbank->ncols; col++)
        {
            int row = 0;
            for (int adj = 0; adj < A->n; adj++)
            {
                for (int ch = 0; ch < ninput_channels; ch++)
                {
                    iftMatrixElem(kernelbank, col, row) =
                        iftMatrixElem(kernelbank, col, row) / stdev[ch];
                    bias[col] -= (mean[ch] * iftMatrixElem(kernelbank, col, row));
                    row++;
                }
            }
        }
        sprintf(filename, "%s/conv%d", param_dir, layer_index);
        WriteBias(filename, bias, kernelbank->ncols);
        iftFree(bias);
    }

    sprintf(filename, "%s/conv%d-kernels.npy", param_dir, layer_index);
    iftWriteMatrix(kernelbank, filename);
    iftDestroyMatrix(&kernelbank);
    sprintf(filename, "%s/conv%d", param_dir, layer_index);
    WriteMeanStdev(filename, mean, stdev, ninput_channels);
    iftFree(mean);
    iftFree(stdev);
    iftDestroyAdjRel(&A);

    /* Apply convolutional layer using the consensus kernel bank and the statistics from all markers */
    int pool_stride = arch->layer[0].pool_stride;
    arch->layer[0].pool_stride = 1;
    for (int i = 0; i < nimages; i++)
    {
        basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
        sprintf(filename, "tmp/%s.mimg", basename);
        input = iftReadMImage(filename);

        output = FLIMConvolutionalLayer(&input, 1, NULL, arch, 0, layer_index, atrous_factor, param_dir);
        iftDestroyMImage(&input);

        sprintf(filename, "%s/%s.mimg", output_dir, basename);
        iftWriteMImage(output[0], filename);
        iftDestroyMImage(&output[0]);
        iftFree(output);
        iftFree(basename);
    }
    arch->layer[0].pool_stride = pool_stride;

    // updating atrous factor with pooling of current layer
    if (arch->apply_intrinsic_atrous)
    {
        atrous_factor *= arch->layer[0].pool_stride;
    }
    // writing atrous factor on file
    sprintf(filename, "%s/intrinsic_atrous.txt", param_dir);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%d", atrous_factor);
    fclose(fp);

    iftRemoveDir("tmp");
    iftDestroyFileSet(&fs_seeds);
    iftDestroyFileSet(&fs_activ);
}

void iftFLIMExtractFeatures(char *orig_dir, char *image_list, iftFLIMArch *arch, char *param_dir,
                            char *feat_dir, char *object_dir, int device)
{
    iftFileSet *fs_images = iftLoadFileSetFromDirOrCSV(image_list, 1, 1);
    int nimages = fs_images->n;
    int batchsize = 0;
    int nbatches = 0;
    int nremaining_images = nimages;
    iftMImage *input = NULL;
    char filename[200], ext[10];
    bool batch_process = true;
    int ninput_channels = 0;
    int input_image_size = 0;
    sprintf(ext, "%s", iftFileExt(fs_images->files[0]->path));

    /* Verify if all images have the same dimension for batch processing */

    sprintf(filename, "%s/%s", orig_dir, fs_images->files[0]->path);
    if (strcmp(ext, ".mimg") == 0)
    {
        input = iftReadMImage(filename);
    }
    else
    {
        input = ReadInputMImage(filename);
    }
    int xsize = input->xsize, ysize = input->ysize, zsize = input->zsize;
    ninput_channels = input->m;
    input_image_size = input->n;
    iftDestroyMImage(&input);
    for (int i = 1; i < nimages; i++)
    {
        sprintf(filename, "%s/%s", orig_dir, fs_images->files[i]->path);
        if (strcmp(ext, ".mimg") == 0)
        {
            input = iftReadMImage(filename);
        }
        else
        {
            input = ReadInputMImage(filename);
        }
        if ((input->xsize != xsize) ||
            (input->ysize != ysize) ||
            (input->zsize != zsize))
        {
            batch_process = false;
            iftDestroyMImage(&input);
            break;
        }
    }

    /* Select batch size for GPU/CPU processing */

#ifdef IFT_GPU
    int ndevices = iftNumberOfDevices();
    if ((ndevices == 0) || (device < 0))
    {
        int first = 0, last = nimages - 1;
        FLIMExtractFeaturesPerImage(fs_images, first, last, orig_dir, arch, param_dir, feat_dir, object_dir);
    }
#endif

    if (batch_process)
    { /* process in batch */
#ifdef IFT_GPU
        if (iftStartDevice(device))
        {
            batchsize = iftFLIMBatchSizeGPU(arch, input_image_size, ninput_channels, device);
            batchsize = iftMin(batchsize, nimages);
            nbatches = nimages / batchsize;
            nremaining_images = nremaining_images - nbatches * batchsize;
            for (int batch = 0; batch < nbatches; batch++)
            {
                int first = batch * batchsize, last = first + batchsize - 1;
                // printf("Processing batch %d of %d batches\n", batch + 1, nbatches);
                // fflush(stdout);
                FLIMExtractFeaturesPerBatch(fs_images, first, last, orig_dir, arch, param_dir, feat_dir, object_dir);
            }
            int first = nimages - nremaining_images, last = nimages - 1;
            if (first <= last)
            {
                // printf("Processing remaining %d images\n", last - first + 1);
                // fflush(stdout);
                FLIMExtractFeaturesPerBatch(fs_images, first, last, orig_dir, arch, param_dir, feat_dir, object_dir);
            }
            iftStopDevice(device);
        }
        else
        {
            int first = 0, last = nimages - 1;
            FLIMExtractFeaturesPerImage(fs_images, first, last, orig_dir, arch, param_dir, feat_dir, object_dir);
        }
#else
        batchsize = iftFLIMBatchSizeCPU(arch, input_image_size, ninput_channels);
        batchsize = iftMin(batchsize, nimages);
        nbatches = nimages / batchsize;
        nremaining_images = nremaining_images - nbatches * batchsize;
        for (int batch = 0; batch < nbatches; batch++)
        {
            int first = batch * batchsize, last = first + batchsize - 1;
            // printf("Processing batch %d of %d batches\n", batch + 1, nbatches);
            // fflush(stdout);
            FLIMExtractFeaturesPerBatch(fs_images, first, last, orig_dir, arch, param_dir, feat_dir, object_dir);
        }
        int first = nimages - nremaining_images, last = nimages - 1;
        if (first <= last)
        {
            // printf("Processing remaining %d images\n", last - first + 1);
            // fflush(stdout);
            FLIMExtractFeaturesPerBatch(fs_images, first, last, orig_dir, arch, param_dir, feat_dir, object_dir);
        }
#endif
    }
    else
    { /* process per image */
        int first = 0, last = nimages - 1;
        FLIMExtractFeaturesPerImage(fs_images, first, last, orig_dir, arch, param_dir, feat_dir, object_dir);
    }

    iftDestroyFileSet(&fs_images);
}

void iftFLIMExtractFeaturesFromLayer(char *orig_dir, char *image_list, iftFLIMArch *arch, char *param_dir, int layer_index,
                                     char *feat_dir, char *object_dir, int device)
{
    iftFileSet *fs_images = iftLoadFileSetFromDirOrCSV(image_list, 1, 1);
    int nimages = fs_images->n;
    int batchsize = 0;
    int nbatches = 0;
    int nremaining_images = nimages;
    iftMImage *input = NULL;
    char filename[200], ext[10];
    bool batch_process = true;
    int ninput_channels = 0;
    int input_image_size = 0;
    sprintf(ext, "%s", iftFileExt(fs_images->files[0]->path));

    /* Verify if all images have the same dimension for batch processing */

    sprintf(filename, "%s/%s", orig_dir, fs_images->files[0]->path);
    if (strcmp(ext, ".mimg") == 0)
    {
        input = iftReadMImage(filename);
    }
    else
    {
        input = ReadInputMImage(filename);
    }

    int xsize = input->xsize, ysize = input->ysize, zsize = input->zsize;
    ninput_channels = input->m;
    input_image_size = input->n;
    iftDestroyMImage(&input);
    for (int i = 1; i < nimages; i++)
    {
        sprintf(filename, "%s/%s", orig_dir, fs_images->files[i]->path);
        if (strcmp(ext, ".mimg") == 0)
        {
            input = iftReadMImage(filename);
        }
        else
        {
            input = ReadInputMImage(filename);
        }
        if ((input->xsize != xsize) ||
            (input->ysize != ysize) ||
            (input->zsize != zsize))
        {
            batch_process = false;
            iftDestroyMImage(&input);
            break;
        }
    }

    /* Select batch size for GPU/CPU processing */

#ifdef IFT_GPU
    int ndevices = iftNumberOfDevices();
    if ((ndevices == 0) || (device < 0))
    {
        int first = 0, last = nimages - 1;
        // printf("FLIMExtractFeaturesFromLayerPerImage A \n");
        FLIMExtractFeaturesFromLayerPerImage(fs_images, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
    }
#endif

    if (batch_process)
    { /* process in batch */
#ifdef IFT_GPU
        if (iftStartDevice(device))
        {
            batchsize = iftFLIMBatchSizeGPU(arch, input_image_size, ninput_channels, device);
            batchsize = iftMin(batchsize, nimages);
            nbatches = nimages / batchsize;
            nremaining_images = nremaining_images - nbatches * batchsize;
            for (int batch = 0; batch < nbatches; batch++)
            {
                int first = batch * batchsize, last = first + batchsize - 1;
                // printf("Processing batch %d of %d batches\n", batch + 1, nbatches);
                // fflush(stdout);
                FLIMExtractFeaturesFromLayerPerBatch(fs_images, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
            }
            int first = nimages - nremaining_images, last = nimages - 1;
            if (first < last)
            {
                // printf("Processing remaining %d images\n", last - first + 1);
                // fflush(stdout);
                FLIMExtractFeaturesFromLayerPerBatch(fs_images, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
            }
            iftStopDevice(device);
        }
        else
        {
            int first = 0, last = nimages - 1;
            // printf("FLIMExtractFeaturesFromLayerPerImage B\n"):
            FLIMExtractFeaturesFromLayerPerImage(fs_images, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
        }
#else
        batchsize = iftFLIMBatchSizeCPU(arch, input_image_size, ninput_channels);
        batchsize = iftMin(batchsize, nimages);
        nbatches = nimages / batchsize;
        nremaining_images = nremaining_images - nbatches * batchsize;
        for (int batch = 0; batch < nbatches; batch++)
        {
            int first = batch * batchsize, last = first + batchsize - 1;
            // printf("Processing batch %d of %d batches\n", batch + 1, nbatches);
            // fflush(stdout);
            FLIMExtractFeaturesFromLayerPerBatch(fs_images, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
        }
        int first = nimages - nremaining_images, last = nimages - 1;
        if (first < last)
        {
            // printf("Processing remaining %d images\n", last - first + 1);
            // fflush(stdout);
            FLIMExtractFeaturesFromLayerPerBatch(fs_images, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
        }
#endif
    }
    else
    { /* process per image */
        int first = 0, last = nimages - 1;
        // printf("FLIMExtractFeaturesFromLayerPerImage C\n");
        FLIMExtractFeaturesFromLayerPerImage(fs_images, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
        // printf("end \n");
    }

    iftDestroyFileSet(&fs_images);
}

void iftFLIMGraphExtractFeaturesFromLayer(char *orig_dir, char *graph_list, iftFLIMArch *arch, char *param_dir, int layer_index,
                                          char *feat_dir, char *object_dir, int device)
{
    iftFileSet *fs_graphs = iftLoadFileSetFromDirOrCSV(graph_list, 1, 1);
    int ngraphs = fs_graphs->n;
    int batchsize = 0;
    int nbatches = 0;
    int nremaining_graphs = ngraphs;
    iftFLIMGraph *input = NULL;
    char filename[200], ext[10];
    bool batch_process = true;
    int ninput_channels = 0;
    int input_image_size = 0;

    sprintf(ext, "%s", iftFileExt(fs_graphs->files[0]->path));

    /* Verify if all graphs have the same dimension for batch processing */

    sprintf(filename, "%s/%s", orig_dir, fs_graphs->files[0]->path);
    input = iftReadFLIMGraph(filename);

    ninput_channels = input->num_feats;
    input_image_size = input->num_nodes;
    iftDestroyFLIMGraph(&input);

    for (int i = 1; i < ngraphs; i++)
    {
        sprintf(filename, "%s/%s", orig_dir, fs_graphs->files[i]->path);
        input = iftReadFLIMGraph(filename);

        if ((input_image_size != input->num_nodes))
        {
            batch_process = false;
            iftDestroyFLIMGraph(&input);
            break;
        }
    }

    /* Select batch size for GPU/CPU processing */

#ifdef IFT_GPU
    int ndevices = iftNumberOfDevices();
    if ((ndevices == 0) || (device < 0))
    {
        int first = 0, last = ngraphs - 1;
        FLIMGraphExtractFeaturesFromLayerPerImage(fs_graphs, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
    }
#endif

    if (batch_process)
    { /* process in batch */
#ifdef IFT_GPU
        if (iftStartDevice(device))
        {
            batchsize = iftFLIMBatchSizeGPU(arch, input_image_size, ninput_channels, device);
            batchsize = iftMin(batchsize, ngraphs);
            nbatches = ngraphs / batchsize;
            nremaining_graphs = nremaining_graphs - nbatches * batchsize;
            for (int batch = 0; batch < nbatches; batch++)
            {
                int first = batch * batchsize, last = first + batchsize - 1;
                // printf("Processing batch %d of %d batches\n", batch + 1, nbatches);
                // fflush(stdout);
                FLIMGraphExtractFeaturesFromLayerPerBatch(fs_graphs, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
            }
            int first = ngraphs - nremaining_graphs, last = ngraphs - 1;
            if (first < last)
            {
                // printf("Processing remaining %d graphs\n", last - first + 1);
                // fflush(stdout);
                FLIMGraphExtractFeaturesFromLayerPerBatch(fs_graphs, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
            }
            iftStopDevice(device);
        }
        else
        {
            int first = 0, last = ngraphs - 1;
            FLIMGraphExtractFeaturesFromLayerPerImage(fs_graphs, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
        }
#else
        batchsize = iftFLIMBatchSizeCPU(arch, input_image_size, ninput_channels);
        batchsize = iftMin(batchsize, ngraphs);
        nbatches = ngraphs / batchsize;
        nremaining_graphs = nremaining_graphs - nbatches * batchsize;

        for (int batch = 0; batch < nbatches; batch++)
        {
            int first = batch * batchsize, last = first + batchsize - 1;
            // printf("Processing batch %d of %d batches\n", batch + 1, nbatches);
            // fflush(stdout);
            FLIMGraphExtractFeaturesFromLayerPerBatch(fs_graphs, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
        }
        int first = ngraphs - nremaining_graphs, last = ngraphs - 1;
        if (first < last)
        {
            // printf("Processing remaining %d graphs\n", last - first + 1);
            // fflush(stdout);
            FLIMGraphExtractFeaturesFromLayerPerBatch(fs_graphs, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
        }
#endif
    }
    else
    { /* process per image */
        int first = 0, last = ngraphs - 1;
        FLIMGraphExtractFeaturesFromLayerPerImage(fs_graphs, first, last, orig_dir, arch, param_dir, layer_index, feat_dir, object_dir);
    }

    iftDestroyFileSet(&fs_graphs);
}


int iftFLIMBatchSizeCPU(iftFLIMArch *arch, int input_image_nvoxels, int input_image_nchannels)
{
    float freeMemoryMb, percMemory = 0.85;

    freeMemoryMb = (float)iftGetFreePhysicalSystemMemory() / 1024.0 / 1024.0;
    // printf("CPU: %.0fGB of free memory.\n", freeMemoryMb / 1024.0);

    float nMbytesPerDouble = sizeof(double) / 1024.0 / 1024.0;
    float ninput_channels = input_image_nchannels;
    float nMbytesInputImage = input_image_nvoxels * ninput_channels * nMbytesPerDouble;
    float max_requiredMemory = 0.0;
    for (int l = 0; l < arch->nlayers; l++)
    {
        int KS2 = 1;
        if (arch->layer[l].kernel_size[2] != 0)
            KS2 = arch->layer[l].kernel_size[2];
        float kernelsize = arch->layer[l].kernel_size[0] * arch->layer[l].kernel_size[1] * KS2;
        float nMbytesKernels = arch->layer[l].noutput_channels * kernelsize * ninput_channels * nMbytesPerDouble;
        float nMbytesOutputImage = input_image_nvoxels * arch->layer[l].noutput_channels * nMbytesPerDouble;
        float requiredMemory = nMbytesKernels + nMbytesOutputImage + nMbytesInputImage;

        if (requiredMemory > max_requiredMemory)
        {
            max_requiredMemory = requiredMemory;
        }
        ninput_channels = arch->layer[l].noutput_channels;
        nMbytesInputImage = nMbytesOutputImage / arch->layer[l].pool_stride;
        input_image_nvoxels = input_image_nvoxels / arch->layer[l].pool_stride;
    }

    int batchsize = (int)floor(percMemory * freeMemoryMb / max_requiredMemory);

    return (batchsize);
}

iftAdjRel *iftFLIMAdjRelFromKernel(iftFLIMLayer layer, bool dim3D)
{
    iftAdjRel *A;

    if (dim3D)
    {
        A = iftCuboidWithDilation(layer.kernel_size[0], layer.kernel_size[1], layer.kernel_size[2],
                                  layer.dilation_rate[0], layer.dilation_rate[1], layer.dilation_rate[2]);
    }
    else
    {
        A = iftRectangularWithDilation(layer.kernel_size[0], layer.kernel_size[1], layer.dilation_rate[0],
                                       layer.dilation_rate[1]);
    }

    return (A);
}

iftAdjRel *iftFLIMAdaptiveAdjRelFromKernel(iftFLIMLayer layer, int atrous_factor, bool dim3D)
{
    iftAdjRel *A;

    if (dim3D)
    {
        A = iftCuboidWithDilation(layer.kernel_size[0], layer.kernel_size[1], layer.kernel_size[2],
                                  layer.dilation_rate[0] * atrous_factor,
                                  layer.dilation_rate[1] * atrous_factor,
                                  layer.dilation_rate[2] * atrous_factor);
    }
    else
    {
        A = iftRectangularWithDilation(layer.kernel_size[0], layer.kernel_size[1], layer.dilation_rate[0] * atrous_factor,
                                       layer.dilation_rate[1] * atrous_factor);
    }

    return (A);
}

iftMImage *iftFLIMAtrousAveragePooling(iftMImage *mimg, int width, int height, int depth, int atrous_factor, int stride)
{

    iftAdjRel *A;

    if (iftIs3DMImage(mimg))
    {
        A = iftCuboidWithDilation(width, height, depth, atrous_factor, atrous_factor, atrous_factor);
    }
    else
    {
        A = iftRectangularWithDilation(width, height, atrous_factor, atrous_factor);
    }

    iftMImage *pool = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m);

#pragma omp parallel for
    for (int p = 0; p < mimg->n; p++)
    {
        iftVoxel u = iftMGetVoxelCoord(mimg, p);
        for (int i = 0; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);
                for (int b = 0; b < mimg->m; b++)
                    pool->val[p][b] += mimg->val[q][b];
            }
        }
        for (int b = 0; b < mimg->m; b++)
            pool->val[p][b] /= A->n;
    }

    iftDestroyAdjRel(&A);
    if (stride > 1)
    {
        iftMImage *mimgRet = iftCreateMImage(ceilf(mimg->xsize / (float)stride), ceilf(mimg->ysize / (float)stride),
                                             iftMax(ceilf(mimg->zsize / (float)stride), 1), mimg->m);
        int q = 0;
        iftVoxel u;
        for (u.z = 0; u.z < pool->zsize; u.z = u.z + stride)
        {
            for (u.y = 0; u.y < pool->ysize; u.y = u.y + stride)
            {
                for (u.x = 0; u.x < pool->xsize; u.x = u.x + stride)
                {
                    int p = iftMGetVoxelIndex(pool, u);
                    for (int b = 0; b < pool->m; b++)
                        mimgRet->val[q][b] = pool->val[p][b];
                    q++;
                }
            }
        }

        iftDestroyMImage(&pool);
        return mimgRet;
    }
    return (pool);
}

void iftFLIMGraphAtrousAveragePooling(iftMatrix *conv, iftFLIMGraph *graph, int width, int height, int depth, int atrous_factor, int stride)
{
    // a indexação dos vértices pode alterar o resultado devido ao stride

    // stride: inicia no vértice v de indice 0.
    //         O próximo vértice escolhido depende de stride.
    //         Para stride > 1, pula os vizinhos de v.
    //         A quantidade de pulos é dada por stride-1.

    // pooling: "width" e "height" contêm o tamanho do pooling.
    //          No grafo, para um dado vértice v, seleciona os width*height vizinhos de v
    //          (da mesma forma que se encontra os vizinhos para criar os kernels)
    //          e atualiza as features de v para a média das features dos vizinhos (incluindo v).

    double *dissimilarity = (double *)malloc(graph->num_nodes * sizeof(double));

    // MaxHeap para armazenar até k vizinhos mais similares.
    // Se a quantidade de vizinhos for maior que k, remove o menos similar
    // do heap (que está no topo) usando "pop" e verifica quem é o mais similar
    // para nova inserção.
    iftDHeap *maxHeap = NULL;
    maxHeap = iftCreateDHeap(graph->num_nodes, dissimilarity);
    maxHeap->removal_policy = MAXVALUE;

    // MinHeap armazena todos os vizinhos que serão incluídos na matriz de saída.
    // O minHeap fornece a ordem de inclusão na matriz de saída
    // (do mais similar para o menos similar).
    // A raíz do mineap deve sempre ser o vértice inicial,
    // pois possui maior similaridade consigo mesmo.
    iftDHeap *minHeap;
    minHeap = iftCreateDHeap((int)graph->num_nodes, dissimilarity);
    iftSet *roots = NULL;

    for (int p = 0; p < graph->num_nodes; p++)
    {
        getKNeigborsFLIMGraph(p, graph, width * height, dissimilarity, roots, maxHeap, minHeap, atrous_factor);

        // move as features dos vértices incluídos no minHeap para a matriz de saída
        int i = 0;
        iftRemoveDHeap(minHeap); // remove p from minHeap
        while (!iftEmptyDHeap(minHeap))
        {
            int node_index = iftRemoveDHeap(minHeap);
            for (int b = 0; b < graph->num_feats; b++)
            {
                iftMatrixElem(conv, b, p) += (float)(graph->nodes[node_index].feats[b]);
            }
            i++;
        }
        iftResetDHeap(minHeap);
        for (int b = 0; b < graph->num_feats; b++)
            iftMatrixElem(conv, b, p) /= (float)(width * height);
    }

    iftDestroyDHeap(&maxHeap);
    iftDestroyDHeap(&minHeap);
    iftFree(dissimilarity);
    iftDestroySet(&roots);

    if (stride > 1)
    {
        iftWarning("Stride > 1 not implemented yet", "iftFLIMGraphAtrousAveragePooling");
    }
}


void iftFLIMGraphAtrousMaxPooling(iftMatrix *conv, iftFLIMGraph *graph, int width, int height, int depth, int atrous_factor, int stride)
{

    iftSetMatrix(conv, IFT_INFINITY_FLT_NEG);

    double *dissimilarity = (double *)malloc(graph->num_nodes * sizeof(double));

    // MaxHeap para armazenar até k vizinhos mais similares.
    // Se a quantidade de vizinhos for maior que k, remove o menos similar
    // do heap (que está no topo) usando "pop" e verifica quem é o mais similar
    // para nova inserção.
    iftDHeap *maxHeap = iftCreateDHeap(graph->num_nodes, dissimilarity);
    maxHeap->removal_policy = MAXVALUE;

    // MinHeap armazena todos os vizinhos que serão incluídos na matriz de saída.
    // O minHeap fornece a ordem de inclusão na matriz de saída
    // (do mais similar para o menos similar).
    // A raíz do mineap deve sempre ser o vértice inicial,
    // pois possui maior similaridade consigo mesmo.
    iftDHeap *minHeap;
    minHeap = iftCreateDHeap((int)graph->num_nodes, dissimilarity);
    iftSet *roots = NULL;

    for (int p = 0; p < graph->num_nodes; p++)
    {

        getKNeigborsFLIMGraph(p, graph, width * height,
                              dissimilarity, roots,
                              maxHeap, minHeap, atrous_factor);

        // move as features dos vértices incluídos no minHeap para a matriz de saída
        iftRemoveDHeap(minHeap); // remove p from minHeap
        while (!iftEmptyDHeap(minHeap))
        {
            int node_index = iftRemoveDHeap(minHeap);
            for (int b = 0; b < graph->num_feats; b++)
            {
                if (iftMatrixElem(conv, b, p) < (float)(graph->nodes[node_index].feats[b]))
                    iftMatrixElem(conv, b, p) = (float)(graph->nodes[node_index].feats[b]);
            }
        }
        iftResetDHeap(minHeap);
    }
    iftDestroyDHeap(&maxHeap);
    iftDestroyDHeap(&minHeap);
    iftFree(dissimilarity);
    iftDestroySet(&roots);

    if (stride > 1)
    {
        iftWarning("Stride > 1 not implemented yet", "iftFLIMGraphMaxPooling");
    }
}


iftMImage *iftFLIMAtrousMaxPooling(iftMImage *mimg, int width, int height, int depth, int atrous_factor, int stride)
{
    iftAdjRel *A;

    if (iftIs3DMImage(mimg))
    {
        A = iftCuboidWithDilation(width, height, depth, atrous_factor, atrous_factor, atrous_factor);
    }
    else
    {
        A = iftRectangularWithDilation(width, height, atrous_factor, atrous_factor);
    }

    iftMImage *pool = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m);
    iftSetMImage(pool, IFT_INFINITY_FLT_NEG);

#pragma omp parallel for
    for (int p = 0; p < mimg->n; p++)
    {
        iftVoxel u = iftMGetVoxelCoord(mimg, p);
        for (int i = 0; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);
                for (int b = 0; b < mimg->m; b++)
                    if (mimg->val[q][b] > pool->val[p][b])
                        pool->val[p][b] = mimg->val[q][b];
            }
        }
    }

    iftDestroyAdjRel(&A);

    if (stride > 1)
    {
        iftMImage *mimgRet = iftCreateMImage(ceilf(mimg->xsize / (float)stride), ceilf(mimg->ysize / (float)stride),
                                             iftMax(ceilf(mimg->zsize / (float)stride), 1), mimg->m);
        int q = 0;
        iftVoxel u;
        for (u.z = 0; u.z < pool->zsize; u.z = u.z + stride)
        {
            for (u.y = 0; u.y < pool->ysize; u.y = u.y + stride)
            {
                for (u.x = 0; u.x < pool->xsize; u.x = u.x + stride)
                {
                    int p = iftMGetVoxelIndex(pool, u);
                    for (int b = 0; b < pool->m; b++)
                        mimgRet->val[q][b] = pool->val[p][b];
                    q++;
                }
            }
        }

        iftDestroyMImage(&pool);

        return mimgRet;
    }

    return pool;
}

iftMatrix *iftFLIMSelectKernelsManual(char *kernel_bank_path, char *selected_kernels_path)
{
    iftMatrix *input_kernels = iftReadMatrix(kernel_bank_path);
    iftDict *json = iftReadJson(selected_kernels_path);
    iftIntArray *selKernels = iftGetIntArrayFromDict("selected_kernels", json);
    int nKernels = selKernels->n;
    iftMatrix *output_kernels = iftCreateMatrix(nKernels, input_kernels->nrows);

    /* perform kernel selection */
    for (int c = 0; c < nKernels; c++)
    {
        int k = selKernels->val[c];
        for (int r = 0; r < input_kernels->nrows; r++)
        {
            iftMatrixElem(output_kernels, c, r) = iftSign(k) * iftMatrixElem(input_kernels, abs(k), r);
        }
    }

    iftDestroyMatrix(&input_kernels);
    return output_kernels;
}
