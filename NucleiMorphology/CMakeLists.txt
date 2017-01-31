CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT( NucleiMorphology )

SET( CurrentExe "NucleiMorphology" )

SET( ND_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/Source/Main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Source/NucleiConfig.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Source/IOFunctions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Source/ITKFunctions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Source/UtilityFunctions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Source/ImageProcessing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Source/NucleiDetector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Source/NucleiFeatures.cpp
)

SET( ND_INC
  ${EXTERNAL_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}/Include/NucleiConfig.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Include/IOFunctions.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Include/ITKDeclarations.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Include/ITKFunctions.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Include/UtilityFunctions.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Include/ImageProcessing.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Include/NucleiDetector.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Include/NucleiFeatures.h
)

INCLUDE_DIRECTORIES(
  ${EXTERNAL_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include/
)

ADD_EXECUTABLE( ${CurrentExe}
  ${ND_SRC}
  ${ND_INC}
) 

TARGET_LINK_LIBRARIES( ${CurrentExe}
  ${Boost_LIBRARIES}
  ${ITK_LIBRARIES}
)
