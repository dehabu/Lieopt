cmake_minimum_required(VERSION 2.8)

PROJECT(LIEOPT)
SET(LIEOPT_VERSION 0.1.0)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LIEOPT")
SET(CPACK_PACKAGE_VENDOR "German Ros")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")

INCLUDE(CPack)

SET(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})

SET(RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

SET(CMAKE_BUILD_TYPE distribution)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-msse3 -DEIGEN_NO_DEBUG -O3")

SET(lieopt_sources_cpp
	main.cpp      
	lieopt.cpp
	optimizer.cpp
	SE3.cpp
   )

SET(lieopt_sources_hpp
	lieopt.hpp
	optimizer.hpp
	SE3.hpp
	mytime.hpp
   )

SET(alglib_sources_cpp
	alglib/alglibinternal.cpp
	alglib/alglibmisc.cpp
	alglib/ap.cpp
	alglib/linalg.cpp
	alglib/optimization.cpp
	alglib/solvers.cpp
   )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/alglib)

SET(lieopt_sources ${lieopt_sources_cpp} ${alglib_sources_cpp})

ADD_EXECUTABLE(lieopt ${lieopt_sources})


