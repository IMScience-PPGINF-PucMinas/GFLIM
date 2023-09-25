# Graph-based feature learning from image markers

In this work, we extend FLIM for general image graph modeling, allowing it for a non-strict kernel shape and taking advantage of the adjacency relation between nodes to extract feature vectors based on neighbors’ features. GFLIM is a lightweight Graph-based CNN with little supervision that generalizes the Feature Learning from Image Markers (FLIM) for superpixel graphs to explore superpixel features, providing higher level features than image pixels and reducing information redundancy. The feature extraction in GFLIM is also independent of the graph modeling and is not restricted to superpixel graphs. 

<img src="https://github.com/IsabelaBB/GFLIM/blob/main/imagesReadme/GFLIM_diagram.jpg" width=800>

GFLIM follows a similar pipeline to FLIM. For training, GFLIM estimates the kernels layer-by-layer. Also, if a layer requires retraining, only subsequent layers are affected. Likewise, one can add or remove layers without affecting the previous ones. First, the user provides foreground/background markers for a small set of training images. Then, our proposal computes a graph for each marked image using a superpixel method (e.g., DISF[1]). Afterward, a set of initial kernels is computed for each graph based on the $k$ most similar neighbors of the marked vertices. In GFLIM, representative kernels are computed using a clustering strategy (e.g., using K-means or DBSCAN). Then, the user can select the most suitable kernels by visualizing their object probability maps (saliency maps). Finally, whether another layer is desired, GFLIM extracts features from the previous layer by applying the learned kernels on the previous layer's inputs. 

Our results indicate that the proposed Graph-based FLIM, named GFLIM, not only outperforms FLIM but also produces competitive detections with deep models, even having an architecture thousands of times smaller in the number of parameters.

### Languages supported

C/C++ (Implementation)

### Requirements

- python3
- pip
- pipenv
- cmake
- CBLAS
- LAPACK
- OpenMP

### Compiling

Run make.sh to install PyIFT and compile GFLIM: `bash make.sh`

## Cite
If this work was useful for your research, please cite our paper:

```
@InProceedings{barcelos2023graph,
  title={Graph-based feature learning from image markers},
  author={Barcelos, Isabela Borlido and Jo{\~a}o, Leonardo de Melo and do Patroc{\'\i}nio Jr, Zenilton K. G. and Ewa Kijak and Falc{\~a}o, Alexandre Xavier and Guimar{\~a}es, Silvio Jamil F.},
  booktitle="26th Iberoamerican Congress on Pattern Recognition (CIARP)",
  year="2023",
  publisher="Springer International Publishing",
  address="Cham",
}
```
