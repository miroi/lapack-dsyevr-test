if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    # remove -rdynamic flag offensive for PGI Fortran
    list(REMOVE_ITEM CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-rdynamic")
endif()

if(LAPACK_FOUND)
    add_executable(dsyerv_check  dsyerv_check.F90 eispack.F)
    target_link_libraries(dsyerv_check ${MATH_LIBS})
    # copy input data into build directory
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Jz_SS_matrix.fermirp2-2  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
else()
    message(FATAL_ERROR "LAPACK library not found for the dsyerv_check !")
endif()
