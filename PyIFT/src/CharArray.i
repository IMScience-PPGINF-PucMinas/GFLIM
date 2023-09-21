


typedef struct {
    /** Number of elements. */
    long n;
    /** Array of char values. */
    char *val;
} iftCharArray;

%extend iftCharArray {

	~iftCharArray() {
		iftCharArray* ptr = ($self);
		iftDestroyCharArray(&ptr);
	}
};

