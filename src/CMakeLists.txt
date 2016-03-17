if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    # remove -rdynamic flag offensive for PGI Fortran
    list(REMOVE_ITEM CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-rdynamic")
endif()

add_executable(dsyevr_check  dsyevr_check.F90 eispack.F)
target_link_libraries(dsyevr_check ${MATH_LIBS} ${EXPLICIT_LIBS})
# copy input data into build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Jz_SS_matrix.fermirp2-2  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
