set( EXTERNAL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

set( EXTERNAL_INC
  ${CMAKE_CURRENT_SOURCE_DIR}/ColorChannelExtraction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ColorChannelExtraction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RGBChannelExtraction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/RGBChannelExtraction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/itkCellFeatureGenerator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkCellFeatureGenerator.txx
  ${CMAKE_CURRENT_SOURCE_DIR}/itkCellForegroundExtraction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkCellForegroundExtraction.txx
  ${CMAKE_CURRENT_SOURCE_DIR}/itkCellStatistics.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkCellStatistics.txx
  ${CMAKE_CURRENT_SOURCE_DIR}/itkColorDeconvolutionImageFilter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkColorDeconvolutionImageFilter.txx
  ${CMAKE_CURRENT_SOURCE_DIR}/itkImageToVTKImageFilter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkImageToVTKImageFilter.txx
  ${CMAKE_CURRENT_SOURCE_DIR}/itkLevelSetBasedCellSegmentation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkLevelSetBasedCellSegmentation.txx
  ${CMAKE_CURRENT_SOURCE_DIR}/itkNoVerifyInputInformationGeodesicActiveContourLevelSetImageFilter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkRGBToHSIColorSpacePixelAccessor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkRGBToHSLColorSpacePixelAccessor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkRGBToHSVColorSpacePixelAccessor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkRGBToLabColorSpacePixelAccessor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkRGBToLuvColorSpacePixelAccessor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkRGBToXYZColorSpacePixelAccessor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkRGBToYUV2ColorSpacePixelAccessor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkRGBToYuvColorSpacePixelAccessor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkScalarToFractalImageFilter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkScalarToFractalImageFilter.txx
  ${CMAKE_CURRENT_SOURCE_DIR}/itkSeedExtraction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkSeedExtraction.txx
  ${CMAKE_CURRENT_SOURCE_DIR}/itkShapeLevelSetBasedCellSegmentation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkShapeLevelSetBasedCellSegmentation.txx
  ${CMAKE_CURRENT_SOURCE_DIR}/itkVoronoiBasedCellSegmentation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/itkVoronoiBasedCellSegmentation.txx
  PARENT_SCOPE
)
