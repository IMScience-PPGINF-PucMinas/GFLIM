
INCLUDESDIR = ./include

SRCDIR = ./src

OBJ = ./obj
NVCC = /usr/local/cuda/bin/nvcc

COMPILERGPU = $(NVCC)
COMPILERCPU = $(CC)

#FOR GPU support uncomment the next line (in all the Makefiles that will use GPU)
#IFT_GPU=1

ifeq ($(IFT_GPU), 1)
	SRC = $(SRC_GPU)
	INCLUDE = $(INCLUDESDIR)  -I/usr/local/cuda -I /opt/cuda -D_MWAITXINTRIN_H_INCLUDED
	FILES =  $(OBJ)/check_cuda_device.cu.o $(OBJ)/kernels.cu.o $(OBJ)/sptree.cu.o $(OBJ)/t_sne_gpu.cu.o
	MYFLAGS = -Wno-deprecated-gpu-targets -std=c++11 -lstdc++ --compiler-options '-fPIC'
else
	SRC = $(SRC_CPU)
	INCLUDE = $(INCLUDESDIR)
	FILES = $(OBJ)/sptree.o $(OBJ)/tsne.o
	MYFLAGS = -fPIC -O2 -Wall
endif

all: $(FILES)

$(OBJ)/%.o: $(SRCDIR)/%.cpp
	$(COMPILERCPU) $(MYFLAGS) -c $< -I$(INCLUDE) -o $@

$(OBJ)/%.cu.o: $(SRCDIR)/%.cu
	$(COMPILERGPU) $(MYFLAGS) -c $< -I$(INCLUDE) -o $@


clean:
	rm -f $(OBJ)/*.o;
