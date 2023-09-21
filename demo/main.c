#include "ift.h"
#include <time.h>
#include <omp.h>

/**
 * @brief Calculate a FLIMGraph from the original image and its superpixel labels map.
 * @author  Isabela Borlido.
 * @date    May, 2021.
 * @param num_superpixels   The number of superpixels in 'labeled_img'.
 * @param A                 An adjacency relation between pixels.
 * @param oiginal_img       The original image.
 * @param labeled_img       A map of superpixel labels.
 * @return The FLIMGraph.
 */
iftFLIMGraph* getAdjacencyGraphFromSuperpixels(unsigned long num_superpixels, iftAdjRel *A, iftMImage *oiginal_img, iftImage *labeled_img);

/**
 * @brief Insert a node in the adjacency list of another node.
 * @author  Isabela Borlido.
 * @date    May, 2021.
 * @param graph       The graph.
 * @param node_index  The index of the node to be inserted.
 * @param adj_index   The index of the node that will receive the new node.
 */
void insertFLIMGraphNode(iftFLIMGraph *graph, int node_index, int adj_index);

/**
 * @brief Compute graphs of the images in "image_dir" using DISF, 
 * save their pixel-superpixel maps in "labels_dir", 
 * and save them in "graph_dir"
 * @author  Isabela Borlido.
 * @date    May, 2021.
 * @param image_dir         Directory of the original image.
 * @param graph_dir         Directory of the graph.
 * @param num_init_seeds    Initial number of seeds in DISF.
 * @param num_superpixels   Final number of superpixels.
 */
void ImageToFLIMGraph(char *image_dir, char *labels_dir, char *graph_dir, int num_init_seeds, int num_superpixels);

/**
 * @brief Convert image seeds file to graph seeds file.
 * @author  Isabela Borlido.
 * @date    May, 2021.
 * @param num_nodes     Number of graph nodes.
 * @param labels_file   File with the superpixel labels map.
 * @param seeds_file    File with the seeds file.
 * @param output_file   File with the graph seeds file.
 */
void seedsPixelToGraph(int num_nodes, char *labels_file, char *seeds_file, char *output_file);

/**
 * @brief Convert image seeds files in 'seeds_dir' to graph seeds files in 'output_dir'.
 * Iterate over the files in 'seeds_dir' and call 'seedsPixelToGraph'.
 * @author  Isabela Borlido.
 * @date    May, 2021.
 * @param num_nodes   Number of graph nodes.
 * @param labels_dir  Directory of the superpixel labels map.
 * @param seeds_dir   Directory of the seeds file.
 * @param output_dir  Directory of the graph seeds file.
 */
void seedsPixelToGraphDir(char *graph_dir, char *labels_dir, char *seeds_dir, char *output_dir);

/**
 * @brief Convert a graph (in json format) to a .mimg file 
 * using the superpixel labels map as auxiliar.
 * @author  Isabela Borlido.
 * @date    May, 2021.
 * @param graph_dir    Directory with graph files (in json format).
 * @param labels_dir   Directory with the superpixel labels map.
 * @param output_dir   Directory to write the .mimg file.
 */
void GraphToMIMG(char *graph_dir, char *labels_dir, char *output_dir);

/**
 * @brief Given an activation directory (source) with files in .mimg or .json format, 
 * write the activations' images in destination directory with extension 'ext' (pgm or png)
 * @author  Isabela Borlido.
 * @date    May, 2021.
 * @param source        File with features. Can be in .mimg or .json format.
 * @param destination   Directory to write the images.
 * @param ext           Extension of the output images (pgm or png).
 * @param labels_dir    Directory with the superpixel labels map (pgm).
 */
void writeActivationsImg(char *source, char *destination, char *ext, char *labels_dir);

/**
 * @brief Write the layer's inputs in a csv file. Used in "ExtractFeatures" functions
 * @author  Isabela Borlido.
 * @date    May, 2021.
 * @param csv_file    File to write the csv.
 * @param files_dir   Directory of the activation files.
 */
void writeCSVFiles(char *csv_file, char *files_dir);

/**
 * @brief Selects a set of kernels from a trained model. 
 * This function updates the .npy model file, write a .json file with the selected kernels, 
 * call 'FLIMExtractFeaturesFromLayer' to update the feature kernels in 'layer' folder 
 * @param param_dir   Directory of the trained model.
 * @param orig_dir    Directory with the original images.
 * @param layer       Layer to select kernels.
 * @param seeds_dir   Directory with the seeds files.
 * @param string_kernels   String with the kernels to select.
 * @param labels_dir  Directory with the superpixel labels map.
 * @param write_img   If true, write the kernels in a .pgm file.
 */
void selectKernelsManual(char *param_dir, char *orig_dir, int layer, char *seeds_dir, char *string_kernels, char *labels_dir, bool write_img);

/**
 * @brief Train a FLIMGraph model.
 * @param layer   Layer to train.
 * @param param_dir       Directory to write the FLIMGraph model.
 * @param orig_dir        Directory with the input of layer 1 (an image or a graph in a .json file).
 * @param seeds_dir       Directory with the seeds files.
 * @param string_kernels  String with the kernels to select.
 * @param labels_dir      Directory with the superpixel labels map.
 * @param write_img       If true, write the kernels in a .pgm file.
 */
void FLIMGraphTrain(int layer, char *param_dir, char *orig_dir, char *seeds_dir);
void FLIMTrain(int layer, char *param_dir, char *orig_dir, char *seeds_dir);

void ExtractFeaturesFromLayerForTraining(char *param_dir, char *orig_dir, int layer, char *seeds_dir, char *labels_dir, bool write_img);
void FLIMGraphExtractFeaturesFromLayer(int layer, char *param_dir, char *orig_dir, char *seeds_dir, char *labels_dir, bool write_img);
void FLIMExtractFeaturesFromLayer(int layer, char *param_dir, char *orig_dir, char *seeds_dir, bool write_img);


iftFLIMGraph* getAdjacencyGraphFromSuperpixels(
    unsigned long num_superpixels,
    iftAdjRel *A,
    iftMImage *oiginal_img, iftImage *labeled_img)
{
    iftFLIMGraph *graph = (iftFLIMGraph *)calloc(1, sizeof(iftFLIMGraph));
    graph->num_nodes = num_superpixels;
    graph->num_feats = oiginal_img->m;

    graph->nodes = (iftFLIMGraphNode *)calloc(num_superpixels, sizeof(iftFLIMGraphNode));
    bool **is_adj = (bool **)calloc(num_superpixels, sizeof(bool *));

    for (unsigned long spx_index = 0; spx_index < num_superpixels; spx_index++)
    {
        graph->nodes[spx_index].feats = (double *)calloc(oiginal_img->m, sizeof(double));
        graph->nodes[spx_index].numPixels = 0;
        graph->nodes[spx_index].neighbors_list_head = NULL;
        graph->nodes[spx_index].neighbors_list_tail = NULL;
        graph->nodes[spx_index].numNeighbors = 0;
        graph->nodes[spx_index].index = spx_index;
        is_adj[spx_index] = (bool *)calloc(spx_index+1, sizeof(bool)); // maior índice armazena os vizinhos (evita preencher toda a tabela)
    }

    // for each pixel
    for (int p_index = 0; p_index < labeled_img->n; p_index++)
    {
        unsigned long label = (unsigned long)(labeled_img->val[p_index] - 1); // get its label
        graph->nodes[label].numPixels++;

        // increase feature
        for (unsigned long band_index = 0; band_index < oiginal_img->m; band_index++)
          graph->nodes[label].feats[band_index] += oiginal_img->val[p_index][band_index];

        // look at neighbors
        iftVoxel p_voxel = iftMGetVoxelCoord(oiginal_img, p_index);
        for (int i = 1; i < A->n; i++)
        {
            iftVoxel q_voxel = iftGetAdjacentVoxel(A, p_voxel, i);
            if (iftMValidVoxel(oiginal_img, q_voxel))
            {
                unsigned long q_label = (unsigned long)(labeled_img->val[iftMGetVoxelIndex(oiginal_img, q_voxel)] - 1);
                if (label > q_label && !is_adj[label][q_label])
                {
                    is_adj[label][q_label] = true;

                    insertFLIMGraphNode(graph, label, q_label);
                    insertFLIMGraphNode(graph, q_label, label);
                }
            }
        }
    }

    // for each node, compute its mean color
    for (unsigned long spx_index = 0; spx_index < num_superpixels; spx_index++)
    {
        for (unsigned long band_index = 0; band_index < oiginal_img->m; band_index++)
            graph->nodes[spx_index].feats[band_index] /= (double)(graph->nodes[spx_index].numPixels);
        free(is_adj[spx_index]);
    }

    free(is_adj);
    return graph;
}

void insertFLIMGraphNode(iftFLIMGraph *graph, int node_index, int adj_index)
{
  iftFLIMGraphNodeList *tmp = (iftFLIMGraphNodeList *)malloc(sizeof(iftFLIMGraphNodeList));
    tmp->next = NULL;
    tmp->node = &(graph->nodes[adj_index]);
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
 * Compute graphs of the images in "image_dir" using DISF, 
 * save their pixel-superpixel maps in "labels_dir", 
 * and save them in "graph_dir"
 */
void ImageToFLIMGraph(char *image_dir, char *labels_dir, char *graph_dir, int num_init_seeds, int num_superpixels)
{

  iftFileSet *fs_imgs = iftLoadFileSetFromDir(image_dir, 1);
  //printf("Processing %ld images...\n", fs_imgs->n);

  for(int i=0; i < fs_imgs->n; i++)
  {
    char *basename = NULL;
    char filename[200], ext[10], label_file[200];
    iftImage *label_img, *img, *mask;
    iftMImage *mimg;
    iftAdjRel *A;

    mask = NULL;

    sprintf(ext, "%s", iftFileExt(fs_imgs->files[i]->path));
    basename = iftFilename(fs_imgs->files[i]->path, ext);
    
    //-------------------------
    // compute iDISF
    img = iftReadImageByExt(fs_imgs->files[i]->path);

    // Init other inputs
    if (iftIsColorImage(img))
        mimg = iftImageToMImage(img, LABNorm_CSPACE);
    else
        mimg = iftImageToMImage(img, GRAYNorm_CSPACE);
    iftDestroyImage(&img);

    // 4- and 6-adjacency
    // if(iftIs3DMImage(mimg)) A = iftSpheric(sqrtf(1.0));
    // else A = iftCircular(sqrtf(1.0));
    // 8- and 24-adjacency
    if (iftIs3DMImage(mimg))
        A = iftSpheric(sqrtf(3.0));
    else
        A = iftCircular(sqrtf(2.0));

    label_img = iftDISF(mimg, A, num_init_seeds, num_superpixels, mask);

    if(labels_dir != NULL){
      sprintf(label_file, "%s/%s.pgm", labels_dir, basename);
      iftWriteImageByExt(label_img, label_file);
    }

    //-------------------------
    // Compute graph
    iftFLIMGraph *graph = getAdjacencyGraphFromSuperpixels(num_superpixels, A, mimg, label_img);

    if(graph_dir != NULL){
      sprintf(filename, "%s/%s.json", graph_dir, basename);
      iftWriteFLIMGraph(graph, filename);
    }
    
    // Clear
    iftDestroyImage(&label_img);
    iftDestroyAdjRel(&A);
    iftDestroyMImage(&mimg);
    iftDestroyFLIMGraph(&graph);
    //-------------------------
  }

}


/*
 * Convert image seeds file to graph seeds file
 */
void seedsPixelToGraph(int num_nodes, char *labels_file, char *seeds_file, char *output_file)
{
  iftImage *label_img = iftReadImageByExt(labels_file);
  iftLabeledSet *S = iftReadSeeds(label_img, seeds_file);
  iftLabeledSet *graph_seeds = NULL;

  int *foreground_count = (int*)calloc(num_nodes, sizeof(int));
  int *background_count = (int*)calloc(num_nodes, sizeof(int));

  while (S != NULL)
  {
    int l;
    int p = iftRemoveLabeledSet(&S, &l);
    
    if(l == 1) foreground_count[label_img->val[p]-1]++;
    else background_count[label_img->val[p]-1]++;
  }

  iftDestroyImage(&label_img);
  iftDestroyLabeledSet(&S);

  for (int p = 0; p < num_nodes; p++){
    if (foreground_count[p] > 0 || background_count[p] > 0){
      if (foreground_count[p] > background_count[p])
        iftInsertLabeledSet(&graph_seeds, p, 1);
      else
        iftInsertLabeledSet(&graph_seeds, p, 0);
    }
  }

  iftWriteSeedsGraph(graph_seeds, output_file);
  iftDestroyLabeledSet(&graph_seeds);
}

void seedsPixelToGraphDir(char *graph_dir, char *labels_dir, char *seeds_dir, char *output_dir)
{
  iftFileSet *fs_labels = iftLoadFileSetFromDir(labels_dir, 1);
  iftFileSet *fs_seeds = iftLoadFileSetFromDir(seeds_dir, 1);

  char graph_file[255], ext[10];
  char *basename = NULL;

  basename = iftFilename(fs_seeds->files[0]->path, "-seeds.txt");

  sprintf(graph_file, "%s/%s.json", graph_dir, basename);
  iftFLIMGraph *graph = iftReadFLIMGraph(graph_file);
  int num_nodes = graph->num_nodes;

  iftDestroyFLIMGraph(&graph);
  iftFree(basename);

  for(int i=0; i < fs_seeds->n; i++)
  {
    char seeds_file[255], output_file[255], labels_file[255];

    basename = iftFilename(fs_seeds->files[i]->path, "-seeds.txt");
    sprintf(ext, "%s", iftFileExt(fs_labels->files[i]->path));

    sprintf(seeds_file, "%s/%s-seeds.txt", seeds_dir, basename);
    sprintf(output_file, "%s/%s-seeds_graph.txt", output_dir, basename);
    sprintf(labels_file, "%s/%s%s", labels_dir, basename, ext);

    seedsPixelToGraph(num_nodes, labels_file, seeds_file, output_file);
    iftFree(basename);
    
  }
  iftDestroyFileSet(&fs_labels);
  iftDestroyFileSet(&fs_seeds);
}


/* Convert graph (json file) to image (mimg file) */
void GraphToMIMG(char *graph_dir, char *labels_dir, char *output_dir)
{
  char ext[10];
  iftFileSet *fs_graphs = iftLoadFileSetFromDir(graph_dir, 1);
  iftFileSet *fs_labels = iftLoadFileSetFromDir(labels_dir, 1);

  //printf("Processing %ld graphs...\n", fs_graphs->n);

  if(fs_graphs->n > fs_labels->n){
    printf("Error: number of graphs is greater than number of labels\n");
    exit(1);
  }

  sprintf(ext, "%s", iftFileExt(fs_labels->files[0]->path));

  for(int i=0; i < fs_graphs->n; i++)
  {
    char *basename = NULL;
    char filename[200], label_file[200];

    basename = iftFilename(fs_graphs->files[i]->path, ".json");

    sprintf(label_file, "%s/%s%s", labels_dir, basename, ext);
    
    iftFLIMGraph *graph = iftReadFLIMGraph(fs_graphs->files[i]->path);
    iftImage *label_img = iftReadImageByExt(label_file);
    iftMImage *features = iftGraphToMImage(graph, label_img);

    sprintf(filename, "%s/%s.mimg", output_dir, basename);

    iftWriteMImage(features, filename);
    iftDestroyFLIMGraph(&graph);
    iftDestroyImage(&label_img);
  }
}


/*  
 * Given an activation directory (source) with files in .mimg or .json format, 
 * write the activations' images in destination directory with extension ext (pgm or png)
 */ 
void writeActivationsImg(char *source, char *destination, char *ext, char *labels_dir)
{
  char filename[200];
  char ext_input[10];

  iftFileSet *fs_imgs = iftLoadFileSetFromDir(source, 1);
  
  for(int i=0; i < fs_imgs->n; i++)
  {
    sprintf(ext_input, "%s", iftFileExt(fs_imgs->files[i]->path));
    char *basename = iftFilename(fs_imgs->files[i]->path, ext_input);
    
    if(labels_dir != NULL)
    {
      iftFLIMGraph *graph = iftReadFLIMGraph(fs_imgs->files[i]->path);

      sprintf(filename, "%s/%s.pgm", labels_dir, basename);
      iftImage *img = iftReadImageByExt(filename);

      iftImage **activations = iftGetActivations(graph, img);
      iftDestroyImage(&img);

      int num_activations = (int)graph->num_feats;
    
      sprintf(filename, "%s/%s", destination, basename);
      iftMakeDir(filename);

      #pragma omp parallel for shared(activations, num_activations, destination, basename, ext) private(filename)
      for(int j=0; j < num_activations; j++){
        sprintf(filename, "%s/%s/kernel-%03d.%s", destination, basename, j+1, ext);
        iftWriteImageByExt(activations[j], filename);
      }

      for(int j=0; j < num_activations; j++)
        iftDestroyImage(&(activations[j]));

      iftFree(activations);
      iftDestroyFLIMGraph(&graph);

    }else{

      iftMImage *mimg = iftReadMImage(fs_imgs->files[i]->path);

      for (int j = 0; j < mimg->m; j++)
      {
        iftImage *band = iftMImageToImage(mimg,255,j);
        iftImage *band_resize = NULL;
        if (iftIs3DImage(band)){
          band_resize = iftGetXYSlice(band,band->zsize/2);
        } else {
          band_resize = iftCopyImage(band);
        }
        iftDestroyImage(&band);

        sprintf(filename, "%s/%s", destination, basename);
        iftMakeDir(filename);
        sprintf(filename, "%s/%s/kernel-%03d.%s", destination, basename, j+1, ext);
        iftWriteImageByExt(band_resize, filename);
        iftDestroyImage(&band_resize);
      }

      iftDestroyMImage(&mimg);
    }

    iftFree(basename);
  }

  iftDestroyFileSet(&fs_imgs);
}


/*
 * Write the layer's inputs in a csv file. Used in "ExtractFeatures" functions
 */
void writeCSVFiles(char *csv_file, char *files_dir)
{
  char ext[10], filename[200];
  iftFileSet *fs_imgs = iftLoadFileSetFromDir(files_dir, 1);

  FILE *fp = fopen(csv_file, "w");
  for(int i=0; i < fs_imgs->n; i++)
  {
    sprintf(ext, "%s", iftFileExt(fs_imgs->files[i]->path));

    char *basename = iftFilename(fs_imgs->files[i]->path, ext);
    sprintf(filename, "%s%s", basename,ext);

    fprintf(fp, "%s\n", filename);
    iftFree(basename);
  }
  fclose(fp);
  iftDestroyFileSet(&fs_imgs);
}


/*
 * Selects a set of kernels from a trained model
 */
void selectKernelsManual(char *param_dir, char *orig_dir, int layer, char *seeds_dir, char *string_kernels, char *labels_dir, bool write_img)
{
  char kernel_bank_path[255], manual_kernels_file[255], arch_file[255];
  
  sprintf(kernel_bank_path, "%s/conv%d-kernels.npy", param_dir, layer);
  if(!iftFileExists(kernel_bank_path)) iftError("Arquivo de kernels não existe", "selectKernelsManual");

  sprintf(manual_kernels_file, "%s/manual_kernels%d.json", param_dir, layer);
  
  iftMatrix *kernels = iftReadMatrix(kernel_bank_path);
  int nkernels = kernels->ncols;
  printf("cols: %d rows:%d nkernels:%d \n", kernels->ncols, kernels->nrows, nkernels);

  int nselected_kernels = 0;
  iftDestroyMatrix(&kernels);

  bool *selected_kernels = (bool*)calloc(nkernels, sizeof(bool));
  bool *neg_kernels = (bool*)calloc(nkernels, sizeof(bool));

  char *ptr = strtok(string_kernels, ",");

  while(ptr != NULL)
	{
    bool negative = false;
		int i = 0, min_kernel = 0, max_kernel = 0;
    if((ptr[i] < 48 || ptr[i] > 57) && ptr[i] != '-') break;

    if(ptr[i] == '-') { negative = true; i++;}

    while(ptr[i] != '\0'){
      if(ptr[i] < 48 || ptr[i] > 57) break;
      min_kernel = min_kernel*10 + (ptr[i] - '0');
      i++;
    }
    if(ptr[i] == '-') i++;
    while(ptr[i] != '\0'){
      if(ptr[i] < 48 || ptr[i] > 57) break;
      max_kernel = max_kernel*10 + (ptr[i] - '0');
      i++;
    }
    
    if(min_kernel <= nkernels){
      selected_kernels[min_kernel-1] = true;
      if(negative) neg_kernels[min_kernel-1] = true;
      nselected_kernels++;
    }
    
    min_kernel++;
    while(min_kernel <= max_kernel && min_kernel <= nkernels){
      selected_kernels[min_kernel-1] = true;
      nselected_kernels++;
      min_kernel++;
    }
    ptr = strtok(NULL, ",");
	}
  /*
	while(ptr != NULL)
	{
		//printf("%s \n", ptr);
    int i = 0, min_kernel = 0, max_kernel = 0;
    if(ptr[i] < 48 || ptr[i] > 57) break;

    while(ptr[i] != '\0' && ptr[i] != '-'){
      if(ptr[i] < 48 || ptr[i] > 57) break;
      min_kernel = min_kernel*10 + (ptr[i] - '0');
      i++;
    }
    if(ptr[i] == '-') i++;
    while(ptr[i] != '\0'){
      if(ptr[i] < 48 || ptr[i] > 57) break;
      max_kernel = max_kernel*10 + (ptr[i] - '0');
      i++;
    }
    
    if(min_kernel <= nkernels){
      selected_kernels[min_kernel] = true;
      nselected_kernels++;
    }
    
    min_kernel++;
    while(min_kernel <= max_kernel && min_kernel <= nkernels){
      selected_kernels[min_kernel] = true;
      nselected_kernels++;
      min_kernel++;
    }
    ptr = strtok(NULL, ",");
	}*/

  if(nselected_kernels == 0) iftError("selectKernelsManual", "No kernels selected\n");

  iftIntArray *arr = iftCreateIntArray(nselected_kernels);
  int j=0;
  //printf("selected_kernels: ");
  for(int i=0; i < nkernels; i++){
    if(selected_kernels[i]) { 
      arr->val[j] = i; 
      if(neg_kernels[i]) {arr->val[j] *= -1;} 
      //printf("%d ", arr->val[j]);
      j++;
    }
  }
  //printf("\n");

  iftDict *json_selected_kernels = iftCreateDict();
  iftInsertIntoDict("selected_kernels", arr, json_selected_kernels);
  iftWriteJson(json_selected_kernels, manual_kernels_file);
  iftDestroyDict(&json_selected_kernels);
  iftDestroyIntArray(&arr);
  
  free(selected_kernels);
  free(neg_kernels);

  iftMatrix *output_kernels = iftFLIMSelectKernelsManual(kernel_bank_path, manual_kernels_file);
  iftWriteMatrix(output_kernels, kernel_bank_path);
  iftDestroyMatrix(&output_kernels);

  // Updating input of next layer with selected kernels
  ExtractFeaturesFromLayerForTraining(param_dir, orig_dir, layer, seeds_dir, labels_dir, write_img);

  sprintf(arch_file, "%s/arch.json", param_dir);
  iftFLIMArch *arch = iftReadFLIMArch(arch_file);
  arch->layer[layer-1].noutput_channels = nselected_kernels;
  iftWriteFLIMArch(arch, arch_file);
  iftDestroyFLIMArch(&arch);
  return;
}


/* Train a layer */
void FLIMGraphTrain(int layer, char *param_dir, char *orig_dir, char *seeds_dir)
{
  char arch_file[255], layer_output[255], layer_input[255], arch_layer_file[255];

  if (layer == 1){
    sprintf(layer_input, "%s", orig_dir);
  }else{
    char file_kernels[255], file_mean[255], file_stdev[255];

    sprintf(file_kernels, "%s/conv%d-kernels.npy", param_dir, layer);
    sprintf(file_mean, "%s/conv%d-mean.txt", param_dir, layer);
    sprintf(file_stdev, "%s/conv%d-stdev.txt", param_dir, layer);

    if(!iftFileExists(file_kernels) && !iftFileExists(file_mean) && !iftFileExists(file_stdev))
      FLIMGraphTrain(layer-1, param_dir, orig_dir, seeds_dir); 

    sprintf(layer_input, "%s/layer%d/", param_dir, layer-1);
  }

  sprintf(layer_output, "%s/layer%d/", param_dir, layer);

  if (iftDirExists(layer_output)){
    iftRemoveDir(layer_output);
    iftMakeDir(layer_output);
  } else {
    iftMakeDir(layer_output);
  }

  // Reading architecture
  sprintf(arch_file, "%s/arch.json", param_dir);
  sprintf(arch_layer_file, "%s/arch_layer%d.json", param_dir,layer);

  // save arch on arch_layer%d.json
  iftFLIMArch *arch = iftReadFLIMArch(arch_file);
  
  iftFLIMArch *save_arch = (iftFLIMArch *)calloc(1, sizeof(iftFLIMArch)); 
  save_arch->nlayers = 1;
  save_arch->stdev_factor = arch->stdev_factor;
  save_arch->apply_intrinsic_atrous = arch->apply_intrinsic_atrous;
  save_arch->layer = (iftFLIMLayer *)calloc(1, sizeof(iftFLIMLayer));
  
  for (int i = 0; i < 3; i++){
    save_arch->layer[0].kernel_size[i] = arch->layer[layer-1].kernel_size[i];
    save_arch->layer[0].dilation_rate[i] = arch->layer[layer-1].dilation_rate[i];
    save_arch->layer[0].pool_size[i] = arch->layer[layer-1].pool_size[i];
  }
  save_arch->layer[0].nkernels_per_image = arch->layer[layer-1].nkernels_per_image;
  save_arch->layer[0].nkernels_per_marker = arch->layer[layer-1].nkernels_per_marker;
  save_arch->layer[0].noutput_channels = arch->layer[layer-1].noutput_channels;
  save_arch->layer[0].relu = arch->layer[layer-1].relu;

  save_arch->layer[0].pool_type = (char *)malloc((strlen(arch->layer[layer-1].pool_type)+1)*sizeof(char));
  strcpy(save_arch->layer[0].pool_type, arch->layer[layer-1].pool_type);
  save_arch->layer[0].pool_stride = arch->layer[layer-1].pool_stride;

  iftDestroyFLIMArch(&arch);
  iftWriteFLIMArch(save_arch, arch_layer_file);

  iftFLIMGraphLearnLayer(layer_input, seeds_dir, param_dir,layer,save_arch,layer_output);

  /* Writting architecture in case it changed */
  arch = iftReadFLIMArch(arch_file);

  arch->layer[layer-1].noutput_channels = save_arch->layer[0].noutput_channels;
  iftWriteFLIMArch(arch, arch_file);

  iftDestroyFLIMArch(&arch);
  iftDestroyFLIMArch(&save_arch);
}

void FLIMTrain(int layer, char *param_dir, char *orig_dir, char *seeds_dir)
{
  
  char arch_file[255], layer_output[255], layer_input[255], arch_layer_file[255];

  if (layer == 1){
    sprintf(layer_input, "%s", orig_dir);
  }else{
    char file_kernels[255], file_mean[255], file_stdev[255];

    sprintf(file_kernels, "%s/conv%d-kernels.npy", param_dir, layer);
    sprintf(file_mean, "%s/conv%d-mean.txt", param_dir, layer);
    sprintf(file_stdev, "%s/conv%d-stdev.txt", param_dir, layer);

    if(!iftFileExists(file_kernels) && !iftFileExists(file_mean) && !iftFileExists(file_stdev))
      FLIMTrain(layer-1, param_dir, orig_dir, seeds_dir); 

    sprintf(layer_input, "%s/layer%d/", param_dir, layer-1);
  }
  sprintf(layer_output, "%s/layer%d/", param_dir, layer);

  if (iftDirExists(layer_output)){
    iftRemoveDir(layer_output);
    iftMakeDir(layer_output);
  } else {
    iftMakeDir(layer_output);
  }

  // Reading architecture
  sprintf(arch_file, "%s/arch.json", param_dir);
  sprintf(arch_layer_file, "%s/arch_layer%d.json", param_dir,layer);

  // save arch on arch_layer%d.json
  iftFLIMArch *arch = iftReadFLIMArch(arch_file);
  
  iftFLIMArch *save_arch = (iftFLIMArch *)calloc(1, sizeof(iftFLIMArch)); 
  save_arch->nlayers = 1;
  save_arch->layer = (iftFLIMLayer *)calloc(1, sizeof(iftFLIMLayer));
  save_arch->layer[0].pool_type = (char *)malloc((strlen(arch->layer[layer-1].pool_type)+1)*sizeof(char));

  for (int i = 0; i < 3; i++){
    save_arch->layer[0].kernel_size[i] = arch->layer[layer-1].kernel_size[i];
    save_arch->layer[0].dilation_rate[i] = arch->layer[layer-1].dilation_rate[i];
    save_arch->layer[0].pool_size[i] = arch->layer[layer-1].pool_size[i];
  }

  save_arch->stdev_factor = arch->stdev_factor;
  save_arch->apply_intrinsic_atrous = arch->apply_intrinsic_atrous;
  save_arch->layer[0].nkernels_per_image = arch->layer[layer-1].nkernels_per_image;
  save_arch->layer[0].nkernels_per_marker = arch->layer[layer-1].nkernels_per_marker;
  save_arch->layer[0].noutput_channels = arch->layer[layer-1].noutput_channels;
  save_arch->layer[0].relu = arch->layer[layer-1].relu;
  strcpy(save_arch->layer[0].pool_type, arch->layer[layer-1].pool_type);
  save_arch->layer[0].pool_stride = arch->layer[layer-1].pool_stride;

  iftDestroyFLIMArch(&arch);
  iftWriteFLIMArch(save_arch, arch_layer_file);

  iftFLIMLearnLayer(layer_input, seeds_dir, param_dir,layer,save_arch,layer_output);

  /* Writting architecture in case it changed */
  arch = iftReadFLIMArch(arch_file);
  arch->layer[layer-1].noutput_channels = save_arch->layer[0].noutput_channels;
  iftWriteFLIMArch(arch, arch_file);

  iftDestroyFLIMArch(&arch);
  iftDestroyFLIMArch(&save_arch);
}


/*
 * Extract features
 */

void ExtractFeaturesFromLayerForTraining(char *param_dir, char *orig_dir, int layer, char *seeds_dir, char *labels_dir, bool write_img)
{

  char file_kernels[255], file_mean[255], file_stdev[255];
  char arch_file[255], layer_input[255], graph_list[255], activation_dir[255];

  sprintf(file_kernels, "%s/conv%d-kernels.npy", param_dir, layer);
  sprintf(file_mean, "%s/conv%d-mean.txt", param_dir, layer);
  sprintf(file_stdev, "%s/conv%d-stdev.txt", param_dir, layer);

  if(!iftFileExists(file_kernels) && !iftFileExists(file_mean) && !iftFileExists(file_stdev)){
    iftError("Arquivos de kernels, mean e stdev não existem", "ExtractFeaturesFromLayerForTraining");
  }

  sprintf(activation_dir, "%s/layer%d/", param_dir, layer);

  if (layer == 1){
    sprintf(layer_input, "%s", orig_dir);
    sprintf(graph_list, "%s/input.csv", param_dir);
  }else{
    sprintf(layer_input, "%s/layer%d/", param_dir, layer-1);
    sprintf(graph_list, "%s/layer%d.csv", param_dir, layer-1);
  }

  iftFileSet *list = iftLoadFileSetFromDir(layer_input, 1);
  iftFileSet *list_orig = iftLoadFileSetFromDir(orig_dir, 1);

  if (list->n < list_orig->n){
    if(labels_dir != NULL) FLIMGraphExtractFeaturesFromLayer(layer-1, param_dir, orig_dir, seeds_dir, labels_dir, write_img);
    else FLIMExtractFeaturesFromLayer(layer-1, param_dir, orig_dir, seeds_dir, write_img);
  }

  // Creating list of graphs
  writeCSVFiles(graph_list, layer_input);

  // Creating output dir for layer
  if (iftDirExists(activation_dir)){
    iftRemoveDir(activation_dir);
    iftMakeDir(activation_dir);
  } else {
    iftMakeDir(activation_dir);
  }

  // Reading architecture
  sprintf(arch_file, "%s/arch_layer%d.json", param_dir, layer);
  iftFLIMArch *tmp_arch = iftReadFLIMArch(arch_file);

  // Changing pooling stride to 1
  tmp_arch->layer[0].pool_stride = 1;

  if(labels_dir != NULL)
    iftFLIMGraphExtractFeaturesFromLayer(layer_input, graph_list, tmp_arch, param_dir, layer,
                                     activation_dir, NULL, -1);
  else
    iftFLIMExtractFeaturesFromLayer(layer_input, graph_list, tmp_arch, param_dir, layer,
                                     activation_dir, NULL, -1);

  iftDestroyFLIMArch(&tmp_arch);
  
  if(labels_dir != NULL)
    FLIMGraphExtractFeaturesFromLayer(layer, param_dir, orig_dir, seeds_dir, labels_dir, write_img);
  else
    FLIMExtractFeaturesFromLayer(layer, param_dir, orig_dir, seeds_dir, write_img);
  
}

void FLIMGraphExtractFeaturesFromLayer(int layer, char *param_dir, char *orig_dir, char *seeds_dir, char *labels_dir, bool write_img){
  
  char arch_file[255], layer_input[255], graph_list[255], activation_dir[255], kernels_img_dir[255];
  char file_kernels[255], file_mean[255], file_stdev[255];

  sprintf(file_kernels, "%s/conv%d-kernels.npy", param_dir, layer);
  sprintf(file_mean, "%s/conv%d-mean.txt", param_dir, layer);
  sprintf(file_stdev, "%s/conv%d-stdev.txt", param_dir, layer);

  if(!iftFileExists(file_kernels) && !iftFileExists(file_mean) && !iftFileExists(file_stdev))
    FLIMTrain(layer, param_dir, orig_dir, seeds_dir); 
  
  sprintf(activation_dir, "%s/layer%d/", param_dir, layer);

  if (layer == 1){
    sprintf(layer_input, "%s", orig_dir);
    sprintf(graph_list, "%s/input.csv", param_dir);
  }else{
    sprintf(layer_input, "%s/layer%d/", param_dir, layer-1);
    sprintf(graph_list, "%s/layer%d.csv", param_dir, layer-1);
  }

  iftFileSet *list = iftLoadFileSetFromDir(layer_input, 1);
  iftFileSet *list_orig = iftLoadFileSetFromDir(orig_dir, 1);

  if (list->n < list_orig->n){
    FLIMGraphExtractFeaturesFromLayer(layer-1, param_dir, orig_dir, seeds_dir, labels_dir, write_img);
  }

  // Creating list of graphs
  writeCSVFiles(graph_list, layer_input);
  
  // Creating output dir for layer
  if (iftDirExists(activation_dir)){
    iftRemoveDir(activation_dir);
    iftMakeDir(activation_dir);
  } else {
    iftMakeDir(activation_dir);
  }
  
  // Reading architecture
  sprintf(arch_file, "%s/arch_layer%d.json", param_dir,layer);
  iftFLIMArch *tmp_arch = iftReadFLIMArch(arch_file);

  iftFLIMGraphExtractFeaturesFromLayer(layer_input, graph_list, tmp_arch, param_dir, layer,
                                     activation_dir, NULL, -1);

  iftDestroyFLIMArch(&tmp_arch);

  // visualization

  /* Verify if layer<n> exists, which is necessary to visualize activations */
  /*char layer_output[255];
  sprintf(layer_output, "%s/layer%d/", param_dir, layer);
  if (!iftDirExists(layer_output)){
    iftMakeDir(layer_output);
  }*/

  if(write_img){
    char ext[10]; 
    sprintf(ext, "pgm");
    sprintf(kernels_img_dir, "%s/img_activations/layer%d/", param_dir, layer);
    iftRemoveDir(kernels_img_dir);
    iftMakeDir(kernels_img_dir);
    writeActivationsImg(activation_dir, kernels_img_dir, ext, labels_dir);
  }
}

void FLIMExtractFeaturesFromLayer(int layer, char *param_dir, char *orig_dir, char *seeds_dir, bool write_img){
  
  char arch_file[255], layer_input[255], graph_list[255], activation_dir[255], kernels_img_dir[255];
  char file_kernels[255], file_mean[255], file_stdev[255];

  sprintf(file_kernels, "%s/conv%d-kernels.npy", param_dir, layer);
  sprintf(file_mean, "%s/conv%d-mean.txt", param_dir, layer);
  sprintf(file_stdev, "%s/conv%d-stdev.txt", param_dir, layer);

  if(!iftFileExists(file_kernels) && !iftFileExists(file_mean) && !iftFileExists(file_stdev))
    FLIMTrain(layer, param_dir, orig_dir, seeds_dir); 

  sprintf(activation_dir, "%s/layer%d/", param_dir, layer);

  if (layer == 1){
    sprintf(layer_input, "%s", orig_dir);
    sprintf(graph_list, "%s/input.csv", param_dir);
  }else{
    sprintf(layer_input, "%s/layer%d/", param_dir, layer-1);
    sprintf(graph_list, "%s/layer%d.csv", param_dir, layer-1);
  }

  iftFileSet *list = iftLoadFileSetFromDir(layer_input, 1);
  iftFileSet *list_orig = iftLoadFileSetFromDir(orig_dir, 1);

  if (list->n < list_orig->n){
    FLIMExtractFeaturesFromLayer(layer-1, param_dir, orig_dir, seeds_dir, write_img);
  }

  // Creating list of graphs
  writeCSVFiles(graph_list, layer_input);

  // Creating output dir for layer
  if (iftDirExists(activation_dir)){
    iftRemoveDir(activation_dir);
    iftMakeDir(activation_dir);
  } else {
    iftMakeDir(activation_dir);
  }

  // Reading architecture
  sprintf(arch_file, "%s/arch_layer%d.json", param_dir,layer);
  iftFLIMArch *tmp_arch = iftReadFLIMArch(arch_file);

  iftFLIMExtractFeaturesFromLayer(layer_input, graph_list, tmp_arch, param_dir, layer,
                                     activation_dir, NULL, -1);

  iftDestroyFLIMArch(&tmp_arch);

  // visualization

  /* Verify if layer<n> exists, which is necessary to visualize activations */
  char layer_output[255];
  sprintf(layer_output, "%s/layer%d/", param_dir, layer);
  if (!iftDirExists(layer_output)){
    iftMakeDir(layer_output);
  }

  if(write_img){
    char ext[10]; 
    sprintf(ext, "pgm");
    sprintf(kernels_img_dir, "%s/img_activations/layer%d/", param_dir, layer);
    iftRemoveDir(kernels_img_dir);
    iftMakeDir(kernels_img_dir);
    writeActivationsImg(activation_dir, kernels_img_dir, ext, NULL);
  }
}



void usage()
{
  printf(" Usage: \n");
  printf("       ./flim_graph --graph <image_dir> <labels_dir> <graph_dir> [<initial_seeds> <final_spx>]\n"); // 5
  printf("       ./flim_graph --activations <source_dir> <destination_dir> [<labels_dir>] \n"); // 5
  printf("       ./flim_graph --seeds <graph_dir> <seeds_dir> <labels_dir> <out_seeds_dir> \n"); // 6
  printf("       ./flim_graph --train <orig_dir> <arch_path> <seeds_dir> <model_dir> <layers> [<labels_dir> <bool write_img>]\n"); // 8
  printf("       ./flim_graph --extract <orig_dir> <arch_path> <seeds_dir> <model_dir> <layers> [<labels_dir> <bool write_img>]\n"); // 8
  printf("       ./flim_graph --select <orig_dir> <seeds_dir> <model_dir> <layer> <kernels_list> [<labels_dir> <bool write_img>] \n"); // 8
  printf("       ./flim_graph --graph_to_mimg <graph_dir> <mimg_dir> <labels_dir>\n"); 
  printf(" Options: \n"); 
  printf("      --graph       : Compute graphs of the images in 'image_dir' using DISF, \n");
  printf("                      save their pixel-superpixel maps in 'labels_dir',\n");
  printf("                      and save the graphs in 'graph_dir'\n"); 
  printf("      --activations : Convert the .mimg or .json files in 'source_dir' to activations images, \n");
  printf("                      writting them in destination directory (destination_dir). \n");
  printf("                      For .json source files, a directory with pixel-superpixel labels map (labels_dir) is required.\n");
  printf("      --seeds       : Convert image seeds files in 'seeds_dir' directory to graph seeds files in 'output_dir' directory.\n");
  printf("      --train       : Train a FLIM or FLIMGraph model and write the output in 'model_dir'. \n");
  printf("      --extract     : Extract feature from a trained FLIM or FLIMGraph model. \n");
  printf("      --select      : Select kernels from a trained FLIM or FLIMGraph model. \n");
  printf(" Args: \n"); 
  printf("      image_dir       : Input directory with images to compute graphs.\n");
  printf("      labels_dir      : Directory with pixel-superpixel maps.\n");
  printf("      graph_dir       : Output directory for graphs.\n");
  printf("      initial_seeds   : Number of initial seeds for DISF computation.\n");
  printf("      final_spx       : Number of final superpixels for DISF computation.\n");
  printf("      source_dir      : Input folder with .mimg or .json files.\n");
  printf("      destination_dir : Output folder for activations images.\n");
  printf("      seeds_dir       : Folder with seeds files.\n");
  printf("      out_seeds_dir   : Output folder for graph-based seeds files.\n");
  printf("      model_dir       : Folder for trained model.\n");
  printf("      layers          : List of layers to train. \n");
  printf("                        Example for 3 layers: 1,2,3. Example for the 2nd layer: 2\n");
  printf("      layer           : A layer (An integer). \n");
  printf("      kernels_list    : A list of kernels to select. Example: 1,2,6,7,10\n");
  printf("      write_img       : {0,1} A boolean value to decide whether the activation files (just for visuzalization) \n");
  printf("                        will be written on the 'model_dir'/img_activations folder.\n");
  exit(1);
}

int main(int argc, const char* argv[])
{
  char param_dir[255], orig_dir[255], labels_dir[255], seeds_dir[255], arch_path[255], layers[255];
  int nLayers = 2, layer = 1;
  bool write_img = false;

  if(argc < 2) usage();
  
  if(strcmp(argv[1], "--graph") == 0 && argc > 4){
    //printf("GRAPH MODE \n");

    if(!iftDirExists(argv[2])) iftError("Input directory does not exist", "main");    
    sprintf(orig_dir, "%s", argv[2]);

    if(!iftDirExists(argv[3])) iftMakeDir(argv[3]);
    sprintf(labels_dir, "%s", argv[3]);

    if(!iftDirExists(argv[4])) iftMakeDir(argv[4]);
    sprintf(param_dir, "%s", argv[4]);

    int num_init_seeds = 8000;
    int num_superpixels = 1000;

    if(argc > 5) num_init_seeds = atoi(argv[5]);
    if(argc > 6) num_superpixels = atoi(argv[6]);

    if(num_init_seeds < 2 || num_superpixels < 2 || num_init_seeds < num_superpixels){
      usage();
      exit(1);
    }

    ImageToFLIMGraph(orig_dir, labels_dir, param_dir, num_init_seeds, num_superpixels);
    return 0;
  }

  if(strcmp(argv[1], "--activations") == 0 && (argc == 4 || argc == 5))
  {
    //printf("ACTIVATIONS MODE \n");

    //printf("checking directories.. \n");
    
    if(!iftDirExists(argv[2])) iftError("Input directory does not exist", "main");
    sprintf(orig_dir, "%s", argv[2]);
    //printf("... input OK \n");

    if(!iftDirExists(argv[3])) iftMakeDir(argv[3]);
    sprintf(param_dir, "%s", argv[3]);
    //printf("... output OK \n");

    if(argc == 4) 
    {
      writeActivationsImg(orig_dir, param_dir, ".pgm", NULL);
    }else{
      sprintf(labels_dir, "%s", argv[4]);
      //printf("... labels_dir OK \n");
      writeActivationsImg(orig_dir, param_dir, ".pgm", NULL);
    }
    
    return 0;
  }

  if(strcmp(argv[1], "--seeds") == 0 && argc == 6){
    //printf("SEEDS MODE \n");

    //printf("checking directories.. \n");
    if(!iftDirExists(argv[2])) iftError("Graph directory does not exist", "main");
    sprintf(orig_dir, "%s", argv[2]);
    //printf("... graph OK \n");

    if(!iftDirExists(argv[3])) iftError("Seeds directory does not exist", "main");
    sprintf(seeds_dir, "%s", argv[3]);
    //printf("... seeds OK \n");

    if(!iftDirExists(argv[4])) iftError("Labels directory does not exist", "main");
    sprintf(labels_dir, "%s", argv[4]);
    //printf("... labels OK \n");

    if(iftDirExists(argv[5])) iftRemoveDir(argv[5]);
    iftMakeDir(argv[5]);
    sprintf(param_dir, "%s", argv[5]);
    //printf("... seeds graph dir OK \n");

    seedsPixelToGraphDir(orig_dir, labels_dir, seeds_dir, param_dir);
    return 0;
  }

  if(strcmp(argv[1], "--select") == 0 && (argc > 6 && argc < 10))
  {
    //printf("SELECT MODE \n");

    //printf("checking directories.. \n"); // <orig_dir> <seeds_dir> <model_dir> <layer> [<labels_dir>]
    
    if(!iftDirExists(argv[2])) iftError("Input directory does not exist", "main");
    sprintf(orig_dir, "%s", argv[2]);
    //printf("... input OK \n");

    if(!iftDirExists(argv[3])) iftError("Seeds directory does not exist", "main");
    sprintf(seeds_dir, "%s", argv[3]);
    //printf("... seeds OK \n");

    if(!iftDirExists(argv[4])) iftMakeDir(argv[4]);
    sprintf(param_dir, "%s", argv[4]);
    //printf("... model OK \n");
    
    layer = atoi(argv[5]);
    //printf("... layer %d OK \n", layer);

    char *kernels_list = (char*)malloc(255*sizeof(char));
    sprintf(kernels_list, "%s", argv[6]);

    if(argc == 7 || (argc == 8 && strcmp(argv[7], "1") == 0 || strcmp(argv[7], "0") == 0)) 
    {
      if(argc > 7) write_img = atoi(argv[7]);
      selectKernelsManual(param_dir, orig_dir, layer, seeds_dir, kernels_list, NULL, write_img);
    }else{
      if(argc > 8) write_img = atoi(argv[8]);
      sprintf(labels_dir, "%s", argv[7]);
      selectKernelsManual(param_dir, orig_dir, layer, seeds_dir, kernels_list, labels_dir, write_img);
    }

    return 0;
  }

  if(strcmp(argv[1], "--train") == 0 && (argc > 6 && argc < 10)){
    
    //printf("TRAIN MODE \n");
    
    if(!iftDirExists(argv[2])) iftError("Input directory does not exist", "main");
    sprintf(orig_dir, "%s", argv[2]);
    //printf("... input OK \n");
    
    if(!iftDirExists(argv[5])) iftMakeDir(argv[5]);
    sprintf(param_dir, "%s", argv[5]);
    //printf("... model OK \n");

    if(!iftFileExists(argv[3])) iftError("Architecture file does not exist", "main");
    //printf("... architecture file OK \n");
    
    if(!iftDirExists(argv[4])) iftError("Seeds directory does not exist", "main");
    sprintf(seeds_dir, "%s", argv[4]);
    //printf("... seeds OK \n");

    //printf("Copying architecture file to %s/arch.json ... ", param_dir);
    // copy architecture file to param_dir
    iftFLIMArch *tmp_arch = iftReadFLIMArch(argv[3]);
    sprintf(arch_path, "%s/arch.json", param_dir);

    nLayers = tmp_arch->nlayers;
    int *train_layers = (int*)calloc(nLayers,sizeof(int));
    
    iftWriteFLIMArch(tmp_arch, arch_path);
    if(tmp_arch != NULL) iftDestroyFLIMArch(&tmp_arch);
    //printf("OK \n");

    char *ptr = strtok(argv[6], ",");
    while(ptr != NULL)
	  {
      int i = 0, layer = 0;
      if(ptr[i] < 48 || ptr[i] > 57) break;
      while(ptr[i] != '\0')
      {
        if(ptr[i] < 48 || ptr[i] > 57) break;
        layer = layer*10 + (ptr[i] - '0');
        i++;
      }
      if(layer > 0 && layer <= nLayers) train_layers[layer-1] = 1;
      ptr = strtok(NULL, ",");
    }

    printf("... Layers = ");
    for(int i=0; i < nLayers; i++) if(train_layers[i]) printf("%d ", i+1);
    printf("OK \n");
    
    if(argc == 7 || (argc == 8 && (strcmp(argv[7],"1") == 0 || strcmp(argv[7],"0") == 0))){ 
      if(argc == 8) write_img = atoi(argv[7]);
      
      for(int i=0; i < nLayers; i++){
        if(train_layers[i]){
          FLIMTrain(i+1, param_dir, orig_dir, seeds_dir);
          FLIMExtractFeaturesFromLayer(i+1, param_dir, orig_dir, seeds_dir, write_img);
        }
      }

    }else{
      if(argc > 8) write_img = atoi(argv[8]);
      if(!iftDirExists(argv[7])) iftError("Labels directory does not exist", "main");
      sprintf(labels_dir, "%s", argv[7]);

      //int opt = -1;
      //if(argc == 10) opt = atoi(argv[9]); 

      for(int i=0; i < nLayers; i++){
        if(train_layers[i]){
          FLIMGraphTrain(i+1, param_dir, orig_dir, seeds_dir);
          FLIMGraphExtractFeaturesFromLayer(i+1, param_dir, orig_dir, seeds_dir, labels_dir, write_img);
        }
      }
    }

    return 1;
  }
  
  if(strcmp(argv[1], "--extract") == 0 && (argc > 6 && argc < 10)){
    
    //printf("EXTRACT MODE \n");

    //printf("checking directories.. \n");
    
    if(!iftDirExists(argv[2])) iftError("Input directory does not exist", "main");
    sprintf(orig_dir, "%s", argv[2]);
    //printf("... input OK \n");
    
    if(!iftDirExists(argv[5])) iftMakeDir(argv[5]);
    sprintf(param_dir, "%s", argv[5]);
    //printf("... model OK \n");

    if(!iftFileExists(argv[3])) iftError("Architecture file does not exist", "main");
    //printf("... architecture file OK \n");
    
    if(!iftDirExists(argv[4])) iftError("Seeds directory does not exist", "main");
    sprintf(seeds_dir, "%s", argv[4]);
    //printf("... seeds OK \n");

    //printf("Copying architecture file to %s/arch.json ... ", param_dir);
    // copy architecture file to param_dir
    iftFLIMArch *tmp_arch = iftReadFLIMArch(argv[3]);
    sprintf(arch_path, "%s/arch.json", param_dir);

    nLayers = tmp_arch->nlayers;
    int *train_layers = (int*)calloc(nLayers,sizeof(int));
    
    iftWriteFLIMArch(tmp_arch, arch_path);
    if(tmp_arch != NULL) iftDestroyFLIMArch(&tmp_arch);
    //printf("OK \n");

    char *ptr = strtok(argv[6], ",");
    while(ptr != NULL)
	  {
      int i = 0, layer = 0;
      if(ptr[i] < 48 || ptr[i] > 57) break;
      while(ptr[i] != '\0')
      {
        if(ptr[i] < 48 || ptr[i] > 57) break;
        layer = layer*10 + (ptr[i] - '0');
        i++;
      }
      if(layer > 0 && layer <= nLayers) train_layers[layer-1] = 1;
      ptr = strtok(NULL, ",");
    }

    //printf("... Layers = ");
    //for(int i=0; i < nLayers; i++) if(train_layers[i]) printf("%d ", i+1);
    //printf("OK \n");

    if(argc == 7 || (argc == 8 && (strcmp(argv[7],"1") == 0 || strcmp(argv[7],"0") == 0))){ 
      if(argc == 8) write_img = atoi(argv[7]);
      //printf("write_img: %d\n", write_img);
      for(int i=0; i < nLayers; i++){
        if(train_layers[i]){
          FLIMExtractFeaturesFromLayer(i+1, param_dir, orig_dir, seeds_dir, write_img);
        }
      }
    }else{
      if(argc == 9) write_img = atoi(argv[8]);
      if(!iftDirExists(argv[7])) iftError("Labels directory does not exist", "main");
      sprintf(labels_dir, "%s", argv[7]);
      //printf("... labels OK \n");

      for(int i=0; i < nLayers; i++){
        if(train_layers[i]){
          FLIMGraphExtractFeaturesFromLayer(i+1, param_dir, orig_dir, seeds_dir, labels_dir, write_img);
        }
      }
    }
    return 0;
  }

  if(strcmp(argv[1], "--graph_to_mimg") == 0 && (argc == 5)){
    
    //printf("GRAPH TO MIMG MODE \n");
    
    if(!iftDirExists(argv[2])) iftError("Input directory does not exist", "main");
    sprintf(orig_dir, "%s", argv[2]);
    //printf("... input OK \n");

    if(!iftDirExists(argv[4])) iftError("Labels directory does not exist", "main");
    sprintf(labels_dir, "%s", argv[4]);
    //printf("... labels OK \n");

    if(!iftDirExists(argv[3])) iftMakeDir(argv[3]);
    sprintf(param_dir, "%s", argv[3]);

    GraphToMIMG(orig_dir, labels_dir, param_dir); 
    return 0;
  }

  usage();

  return 0;
}