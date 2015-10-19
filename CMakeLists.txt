# The name of this project is "potts".  CMakeLists files in this project can
# refer to the root source directory of the project as ${POTTS_SOURCE_DIR} and
# to the root binary directory of the project as ${POTTS_BINARY_DIR}.

cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

PROJECT(POTTS)

# Compiler and path settings:
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZATION_FLAGS}" CACHE INTERNAL "Compiler optimization flags.")
SET(EXECUTABLE_OUTPUT_PATH ${POTTS_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for executables.")
SET(LIBRARY_OUTPUT_PATH ${POTTS_BINARY_DIR}/lib CACHE INTERNAL "Single output directory for libraries.")
SET(OPTIMIZATION_FLAGS "-O3 -funroll-loops -Wall" CACHE STRING "Compiler optimization flags.")
SET(SWIG_DIR /usr/local CACHE STRING "directory containing SWIG.")


# Advanced options in the ccmake gui:

#Set up definitions and libraries:

#################################

#write a configure file 
FILE(WRITE ${POTTS_SOURCE_DIR}/src/config.h "//configure file for potts\n")
FILE(APPEND ${POTTS_SOURCE_DIR}/src/config.h "#define _POTTS_SRC_PATH \"${POTTS_SOURCE_DIR}\"\n")

#################################

#Set up source files:
AUX_SOURCE_DIRECTORY(src POTTS_SRC)

#################################

#make potts.h header file:

FILE(GLOB HEADER_FILES "${POTTS_SOURCE_DIR}/src/*.h")
STRING(REPLACE "${POTTS_SOURCE_DIR}/src/potts.h;" "" HEADER_FILES "${HEADER_FILES}")
STRING(REPLACE "${POTTS_SOURCE_DIR}/src/" "" HEADER_FILES "${HEADER_FILES}")
FILE(WRITE ${POTTS_SOURCE_DIR}/src/potts.h "#ifndef _POTTS_H\n#define _POTTS_H\n\n")
FOREACH(_FILENAME ${HEADER_FILES})
	FILE(APPEND ${POTTS_SOURCE_DIR}/src/potts.h "#include \"${_FILENAME}\"\n")
ENDFOREACH(_FILENAME ${HEADER_FILES})
FILE(APPEND ${POTTS_SOURCE_DIR}/src/potts.h "\n#endif\n\n") 

#install headers:

FILE(GLOB files "${POTTS_SOURCE_DIR}/src/*.h*")
INSTALL(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/potts)

################################

# Create a library called "potts" 

ADD_LIBRARY(potts 
	${POTTS_SRC} ${THIRD_PARTY_SRC}
)

TARGET_LINK_LIBRARIES(potts ${LIBRARIES_USED})

INSTALL(TARGETS potts
			RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
			LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)


#################################
#ADD_SUBDIRECTORY(unit_tests)
#################################
