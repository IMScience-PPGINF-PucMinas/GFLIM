if (LAPACK_LIBRARIES)
    set(LAPACK_FIND_QUIETLY TRUE)
endif (LAPACK_LIBRARIES)
find_file(LAPACK_LIBRARIES
	liblapack.dylib
        liblapack.so
        liblapack.so.3
        liblapack.a
        PATHS
        /usr/lib
	/opt/local/lib
        /usr/local/lib
	/usr/lib/atlas-base
	/usr/lib/x86_64-linux-gnu
        $ENV{LAPACKDIR}/lib
        ${LIB_INSTALL_DIR}
        )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACK DEFAULT_MSG LAPACK_LIBRARIES)
mark_as_advanced(LAPACK_LIBRARIES)
