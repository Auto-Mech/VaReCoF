project(VaReCoF CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(SLATEC REQUIRED NAMES slatec libslatec)
find_library(DL libdl.so PATHS /usr/lib64)
        
message(STATUS "Found MPI: ${MPI_LIBRARIES}")
message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
message(STATUS "Found SLATEC: ${SLATEC}")
message(STATUS "Found DL: ${DL}")

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_library(varecoflibs
    ${PROJECT_SOURCE_DIR}/src/interflux.cc
    ${PROJECT_SOURCE_DIR}/src/sjk.cc
    ${PROJECT_SOURCE_DIR}/src/integral.cc
    ${PROJECT_SOURCE_DIR}/src/raninit.cc
    ${PROJECT_SOURCE_DIR}/src/log.cc
    ${PROJECT_SOURCE_DIR}/src/run.cc
    ${PROJECT_SOURCE_DIR}/src/rfactor.cc
    ${PROJECT_SOURCE_DIR}/src/surface.cc
    ${PROJECT_SOURCE_DIR}/src/pes.cc
    ${PROJECT_SOURCE_DIR}/src/units.cc
    ${PROJECT_SOURCE_DIR}/src/molpro.cc
    ${PROJECT_SOURCE_DIR}/src/multipole.cc
    ${PROJECT_SOURCE_DIR}/src/flux.cc
    ${PROJECT_SOURCE_DIR}/src/comm.cc
    ${PROJECT_SOURCE_DIR}/src/gauss.cc
    ${PROJECT_SOURCE_DIR}/src/expression.cc
    ${PROJECT_SOURCE_DIR}/src/error.cc
    ${PROJECT_SOURCE_DIR}/src/system.cc
    ${PROJECT_SOURCE_DIR}/src/input.cc
    ${PROJECT_SOURCE_DIR}/src/rotd.cc
    ${PROJECT_SOURCE_DIR}/src/force.cc
    ${PROJECT_SOURCE_DIR}/src/random.cc
    ${PROJECT_SOURCE_DIR}/src/divsur.cc
    ${PROJECT_SOURCE_DIR}/src/math.cc)

add_executable(convert_corr ${PROJECT_SOURCE_DIR}/src/convert_corr.cc)  
add_executable(convert_multi ${PROJECT_SOURCE_DIR}/src/convert_multi.cc)  
add_executable(convolute ${PROJECT_SOURCE_DIR}/src/convolute.cc)    
add_executable(cut_multi ${PROJECT_SOURCE_DIR}/src/cut_multi.cc)     
add_executable(ej_flux ${PROJECT_SOURCE_DIR}/src/ej_flux.cc)     
add_executable(mc_flux ${PROJECT_SOURCE_DIR}/src/mc_flux.cc)     
add_executable(tst_test ${PROJECT_SOURCE_DIR}/src/tst_test.cc)       
add_executable(convert_struct ${PROJECT_SOURCE_DIR}/src/convert_struct.cc)       
add_executable(rfactor ${PROJECT_SOURCE_DIR}/src/rfactor.cc)        
add_executable(sampling ${PROJECT_SOURCE_DIR}/src/sampling.cc)      
add_executable(multi ${PROJECT_SOURCE_DIR}/src/multi.cc)       

target_link_libraries(convert_corr 
        varecoflibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(convert_multi 
        varecoflibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(convolute 
        varecoflibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(cut_multi 
        varecoflibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(ej_flux 
        varecoflibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(mc_flux
        varecoflibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(tst_test 
        varecoflibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(convert_struct
        varecoflibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(rfactor 
        varecoflibs ${MPI_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(sampling 
        varecoflibs ${MPI_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
target_link_libraries(multi 
        varecoflibs ${MPI_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} ${DL})       
#install(TARGETS convert_corr DESTINATION bin) 
#install(TARGETS convert_multi DESTINATION bin)  
#install(TARGETS convolute DESTINATION bin)
#install(TARGETS cut_multi DESTINATION bin) 
#install(TARGETS ej_flux DESTINATION bin)
#install(TARGETS mc_flux DESTINATION bin)
#install(TARGETS multi DESTINATION bin)
#install(TARGETS rfactor DESTINATION bin)  
#install(TARGETS sampling DESTINATION bin) 
#install(TARGETS convert_struct DESTINATION bin)  
#install(TARGETS tst_test DESTINATION bin)  
