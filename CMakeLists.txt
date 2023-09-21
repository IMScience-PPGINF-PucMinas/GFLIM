#IFT library CMAKE

cmake_minimum_required(VERSION 3.1)

project(ift)

set(CMAKE_C_STANDARD 11)

set(CMAKE_VERBOSE_MAKEFILE "ON")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/")

message(status ${CMAKE_CURRENT_LIST_DIR})

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

file(GLOB IFT_INCLUDE
        ${CMAKE_CURRENT_LIST_DIR}/include/*.h
        ${CMAKE_CURRENT_LIST_DIR}/include/ift/**/*.h
        ${CMAKE_CURRENT_LIST_DIR}/include/ift/**/**/*.h)

file(GLOB EXTERNAL
        ${CMAKE_CURRENT_LIST_DIR}/externals/*/)

foreach(lib ${EXTERNAL})
    if(EXISTS ${lib}/src)
        message(STATUS ${lib};)
        list(APPEND EXTERNAL_INCLUDE ${lib}/include)
        list(APPEND EXTERNAL_SRC ${lib}/src)
        list(APPEND EXTERNAL_LIB ${lib}/lib)
        message(STATUS ${EXTERNAL_SRC};)
    endif()
endforeach()

#if(IFT_GPU)
#    add_definitions(-DIFT_GPU=1)
#    find_package(CUDA REQUIRED)
#    include(FindCUDA)
#    include_directories(${CUDA_INCLUDE_DIRS})
#    file(GLOB IFT_CUDA_SRC ${CMAKE_CURRENT_LIST_DIR}/src/*.cu)
#endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu11 -g -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall")

if(IFT_DEBUG OR CMAKE_BUILD_TYPE STREQUAL Debug)
    add_definitions(-DIFT_DEBUG=1)
    set(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS} -fsanitize=address")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -fsanitize=address")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

find_package(Threads REQUIRED)

find_package( CBLAS REQUIRED )
find_package( LAPACK REQUIRED )
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

option (USE_LAPACKE "Use LAPACKE C interface to LAPACK" ON)
option (USE_LAPACKE "Use BLAS C interface to BLAS" ON)
option (USE_OPENMP "Use OpenMP for multithreading" ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

link_directories(${EXTERNAL_LIB})
message(Compiler: "${CMAKE_CXX_COMPILER_ID}")
message(Compiler: "${CMAKE_C_COMPILER_ID}")

file(GLOB IFT_SRC
        ${CMAKE_CURRENT_LIST_DIR}/src/*.c
        ${CMAKE_CURRENT_LIST_DIR}/src/ift/**/*.c
        ${CMAKE_CURRENT_LIST_DIR}/src/ift/**/**/*.c)
include_directories(include ${CBLAS_INCLUDES} .)

foreach(ext_src ${EXTERNAL_SRC})
    aux_source_directory(${ext_src} IFT_SRC)
endforeach()

foreach(ext_inc ${EXTERNAL_INCLUDE})
    include_directories(${ext_inc})
endforeach()


message(STATUS "External libs")
message(STATUS "\t${EXTERNAL_SRC}")

list(REMOVE_ITEM IFT_SRC "${ift_SOURCE_DIR}/src/iftNonUsedCode.c")

message(STATUS "IFT library")
foreach( iftsourcefile ${IFT_SRC} )
    message(STATUS "\t${iftsourcefile}")
endforeach(iftsourcefile ${IFT_SRC})

#add_library(ift SHARED ${IFT_SRC} ${IFT_INCLUDE})
add_library(ift STATIC ${IFT_SRC} ${IFT_INCLUDE} ${EXTERNAL_SRC} )

#if(IFT_GPU)
#    cuda_add_library(iftcuda STATIC ${IFT_CUDA_SRC} ${EXTERNAL_SRC} )
#endif()

target_link_libraries(ift ${CBLAS_LIBRARIES} ${LAPACK_LIBRARIES})

list(APPEND IFT_LIBRARIES ift)

#if(IFT_GPU)
#    list(APPEND IFT_LIBRARIES iftcuda)
#endif()

install (TARGETS ift
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

add_subdirectory(demo)
