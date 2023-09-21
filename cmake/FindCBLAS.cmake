if (CBLAS_INCLUDES AND CBLAS_LIBRARIES)
    set(CBLAS_FIND_QUIETLY TRUE)
endif (CBLAS_INCLUDES AND CBLAS_LIBRARIES)
find_path(CBLAS_INCLUDES
        NAMES
        cblas.h
        PATHS
        $ENV{CBLASDIR}/include
		/usr/local/opt/openblas/include
        ${INCLUDE_INSTALL_DIR}
        )
find_file(CBLAS_LIBRARIES
        libcblas.dylib
		libcblas.so
        libcblas.so.3
        libcblas.a
        PATHS
        /usr/lib
        /usr/lib/atlas-base
		/opt/local/lib
        /usr/local/lib
		/usr/local/opt/openblas/lib
	/usr/lib/x86_64-linux-gnu
        $ENV{CBLASDIR}/lib
        ${LIB_INSTALL_DIR}
        )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS DEFAULT_MSG
        CBLAS_INCLUDES CBLAS_LIBRARIES)
mark_as_advanced(CBLAS_INCLUDES CBLAS_LIBRARIES)
