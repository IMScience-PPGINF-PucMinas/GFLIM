OS=$(shell uname -s)

FLAGS = -O3 -Wall -fPIC -pedantic

# Adding the OpenMP flag
# For Mac OS, add a special flag to enable the openmp.  Install: brew install libomp
ifeq ($(OS), Darwin)
	FLAGS += -Xpreprocessor -fopenmp
else
	FLAGS += -fopenmp
endif


all: obj/svm.o
	ar csr lib/libsvm.a obj/svm.o
	@echo
	@echo "libsvm.a built..."
	@echo "DONE."

obj/svm.o: src/svm.cpp
	$(CC) -I include/ $(FLAGS) -c src/svm.cpp -o obj/svm.o

clean:
	rm -f obj/*.o; rm -f lib/*.a;
