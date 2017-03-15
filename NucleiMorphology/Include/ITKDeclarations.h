/*
 * IOFunctions.h
 *
 *  Created on: 20 may 2014
 *      Author: Humayun
 */
#ifndef __ITKDECLARATIONS_H__
#define __ITKDECLARATIONS_H__

#include "itkImage.h"
#include "itkImageAdaptor.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

const unsigned int Dimension = 2;

typedef bool									BoolPixelType;
typedef itk::Image< BoolPixelType, Dimension >	BoolImageType;
typedef itk::ImageFileReader< BoolImageType >	BoolReaderType;
typedef itk::ImageFileWriter< BoolImageType >	BoolWriterType;
typedef BoolImageType::Pointer					BoolImagePointer;
typedef BoolImageType::SizeType					BoolImageSizeType;
typedef BoolImageType::IndexType				BoolImageIndexType;
typedef BoolImageType::PixelType				BoolImagePixelType;
typedef BoolImageType::RegionType				BoolImageRegionType;

typedef unsigned char							CharPixelType;
typedef itk::Image< CharPixelType, Dimension >	CharImageType;
typedef itk::ImageFileReader< CharImageType >	CharReaderType;
typedef itk::ImageFileWriter< CharImageType >	CharWriterType;
typedef CharImageType::OffsetType				CharOffsetType;
typedef CharImageType::Pointer					CharImagePointer;
typedef CharImageType::SizeType					CharImageSizeType;
typedef CharImageType::IndexType				CharImageIndexType;
typedef CharImageType::PixelType				CharImagePixelType;
typedef CharImageType::RegionType				CharImageRegionType;
typedef CharImageType::SpacingType				CharImageSpacingType;

typedef unsigned short							ShortPixelType;
typedef itk::Image< ShortPixelType, Dimension >	ShortImageType;
typedef itk::ImageFileReader< ShortImageType >	ShortReaderType;
typedef itk::ImageFileWriter< ShortImageType >	ShortWriterType;
typedef ShortImageType::OffsetType				ShortOffsetType;
typedef ShortImageType::Pointer					ShortImagePointer;
typedef ShortImageType::SizeType				ShortImageSizeType;
typedef ShortImageType::IndexType				ShortImageIndexType;
typedef ShortImageType::PixelType				ShortImagePixelType;
typedef ShortImageType::RegionType				ShortImageRegionType;

typedef unsigned long							LongPixelType;
typedef itk::Image< LongPixelType, Dimension >	LongImageType;
typedef itk::ImageFileReader< LongImageType >	LongReaderType;
typedef itk::ImageFileWriter< LongImageType >	LongWriterType;
typedef LongImageType::OffsetType				LongOffsetType;
typedef LongImageType::Pointer					LongImagePointer;
typedef LongImageType::SizeType					LongImageSizeType;
typedef LongImageType::IndexType				LongImageIndexType;
typedef LongImageType::PixelType				LongImagePixelType;
typedef LongImageType::RegionType				LongImageRegionType;

typedef float									FloatPixelType;
typedef itk::Image< FloatPixelType, Dimension >	FloatImageType;
typedef itk::ImageFileReader< FloatImageType >	FloatReaderType;
typedef itk::ImageFileWriter< FloatImageType >	FloatWriterType;
typedef FloatImageType::OffsetType				FloatOffsetType;
typedef FloatImageType::Pointer					FloatImagePointer;
typedef FloatImageType::SizeType				FloatImageSizeType;
typedef FloatImageType::IndexType				FloatImageIndexType;
typedef FloatImageType::PixelType				FloatImagePixelType;
typedef FloatImageType::RegionType				FloatImageRegionType;

typedef itk::RGBPixel< unsigned char >			RGBPixelType;
typedef itk::Image< RGBPixelType, Dimension >	RGBImageType;
typedef itk::ImageFileReader< RGBImageType >	RGBReaderType;
typedef itk::ImageFileWriter< RGBImageType >	RGBWriterType;
typedef RGBImageType::OffsetType				RGBOffsetType;
typedef RGBImageType::Pointer					RGBImagePointer;
typedef RGBImageType::SizeType					RGBImageSizeType;
typedef RGBImageType::IndexType					RGBImageIndexType;
typedef RGBImageType::PixelType					RGBImagePixelType;
typedef RGBImageType::RegionType				RGBImageRegionType;

namespace itk
{
/** \class myRGBPixel
 * \brief Extends RGBPixel with operator <=
 *
 * This class overrides the <= and < operators to use Luminance as a sorting value.
 */
template< typename TComponent = unsigned short >
class myRGBPixel : public RGBPixel<TComponent>
{
public:
  typedef myRGBPixel           Self;
  typedef RGBPixel<TComponent> Superclass;
  using RGBPixel<TComponent>::operator=;
 
  bool operator<=(const Self & r) const
  {
    if (this->GetLuminance() <= r.GetLuminance())
      {
      return true;
      }
    else
      {
      return false;
      }
  }
  bool operator<(const Self & r) const
  {
    if (this->GetLuminance() < r.GetLuminance())
      {
      return true;
      }
    else
      {
      return false;
      }
  }
};
}

//Basic
#include "itkSimpleFilterWatcher.h"
#include "itkNumericSeriesFileNames.h"

#include "itkChangeInformationImageFilter.h"
typedef itk::ChangeInformationImageFilter< CharImageType >		CharChangeInformationImageFilterType;
typedef itk::ChangeInformationImageFilter< FloatImageType >		FloatChangeInformationImageFilterType;
typedef itk::ChangeInformationImageFilter< RGBImageType >		RGBChangeInformationImageFilterType;

#include "itkRegionOfInterestImageFilter.h"
typedef itk::RegionOfInterestImageFilter< CharImageType, CharImageType >		CharRegionOfInterestImageFilterType;
typedef itk::RegionOfInterestImageFilter< FloatImageType, FloatImageType >		FloatRegionOfInterestImageFilterType;
typedef itk::RegionOfInterestImageFilter< RGBImageType, RGBImageType >			RGBRegionOfInterestImageFilterType;

#include "itkRescaleIntensityImageFilter.h"
typedef itk::RescaleIntensityImageFilter< FloatImageType,FloatImageType >		FloatToFloatRescaleIntensityImageFilter;
typedef itk::RescaleIntensityImageFilter< FloatImageType, CharImageType >		FloatToCharRescaleIntensityImageFilterType;
typedef itk::RescaleIntensityImageFilter< CharImageType, FloatImageType >		CharToFloatRescaleIntensityImageFilterType;
typedef itk::RescaleIntensityImageFilter< CharImageType,CharImageType >			CharToCharRescaleIntensityImageFilter;

#include "itkImageLinearIteratorWithIndex.h"
typedef itk::ImageLinearIteratorWithIndex< CharImageType >						CharImageLinearIteratorWithIndexType;
typedef itk::ImageLinearIteratorWithIndex< RGBImageType >						RGBImageLinearIteratorWithIndexType;
typedef itk::ImageLinearConstIteratorWithIndex< CharImageType >					CharImageLinearConstIteratorWithIndexType;
typedef itk::ImageLinearConstIteratorWithIndex< RGBImageType >					RGBImageLinearConstIteratorWithIndexType;

#include "itkBinaryFillholeImageFilter.h"
typedef itk::BinaryFillholeImageFilter< CharImageType >							CharBinaryFillholeImageFilterType;
typedef itk::BinaryFillholeImageFilter< FloatImageType >						FloatBinaryFillholeImageFilterType;

#include "itkVotingBinaryHoleFillingImageFilter.h"
typedef itk::VotingBinaryHoleFillingImageFilter< CharImageType, CharImageType >	CharVotingBinaryHoleFillingImageFilterType;
typedef itk::VotingBinaryHoleFillingImageFilter< FloatImageType, FloatImageType > FloatVotingBinaryHoleFillingImageFilterType;

#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
typedef itk::VotingBinaryIterativeHoleFillingImageFilter< CharImageType >		CharVotingBinaryIterativeHoleFillingImageFilterType;
typedef itk::VotingBinaryIterativeHoleFillingImageFilter< FloatImageType >		FloatVotingBinaryIterativeHoleFillingImageFilterType;

//Morphology
#include "itkBinaryBallStructuringElement.h"
typedef itk::BinaryBallStructuringElement< CharPixelType, Dimension >										CharStructuringElementType;
typedef itk::BinaryBallStructuringElement< FloatPixelType, Dimension >										FloatStructuringElementType;

#include "itkBinaryErodeImageFilter.h"
typedef itk::BinaryErodeImageFilter< CharImageType, CharImageType, CharStructuringElementType >				CharBinaryErodeImageFilterType;
typedef itk::BinaryErodeImageFilter< FloatImageType, FloatImageType, FloatStructuringElementType >			FloatBinaryErodeImageFilterType;

#include "itkBinaryDilateImageFilter.h"
typedef itk::BinaryDilateImageFilter< CharImageType, CharImageType, CharStructuringElementType >			CharBinaryDilateImageFilterType;
typedef itk::BinaryDilateImageFilter< FloatImageType, FloatImageType, FloatStructuringElementType >			FloatBinaryDilateImageFilterType;

#include "itkBinaryMorphologicalOpeningImageFilter.h"
typedef itk::BinaryMorphologicalOpeningImageFilter< CharImageType, CharImageType, CharStructuringElementType>		CharBinaryOpeningImageFilterType;
typedef itk::BinaryMorphologicalOpeningImageFilter< FloatImageType, FloatImageType, FloatStructuringElementType>	FloatBinaryOpeningImageFilterType;

#include "itkBinaryMorphologicalClosingImageFilter.h"
typedef itk::BinaryMorphologicalClosingImageFilter< CharImageType, CharImageType, CharStructuringElementType>		CharBinaryClosingImageFilterType;
typedef itk::BinaryMorphologicalClosingImageFilter< FloatImageType, FloatImageType, FloatStructuringElementType>	FloatBinaryClosingImageFilterType;

#include "itkGrayscaleErodeImageFilter.h"
typedef itk::GrayscaleErodeImageFilter< CharImageType, CharImageType, CharStructuringElementType >			CharGrayscaleErodeImageFilterType;
typedef itk::GrayscaleErodeImageFilter< FloatImageType, FloatImageType, FloatStructuringElementType >		FloatGrayscaleErodeImageFilterType;

#include "itkGrayscaleDilateImageFilter.h"
typedef itk::GrayscaleDilateImageFilter< CharImageType, CharImageType, CharStructuringElementType >			CharGrayscaleDilateImageFilterType;
typedef itk::GrayscaleDilateImageFilter< FloatImageType, FloatImageType, FloatStructuringElementType >		FloatGrayscaleDilateImageFilterType;

#include "itkGrayscaleGeodesicErodeImageFilter.h"
typedef itk::GrayscaleGeodesicErodeImageFilter< CharImageType, CharImageType >								CharGrayscaleGeodesicErodeImageFilterType;
typedef itk::GrayscaleGeodesicErodeImageFilter< FloatImageType, FloatImageType >							FloatGrayscaleGeodesicErodeImageFilterType;

#include "itkGrayscaleGeodesicDilateImageFilter.h"
typedef itk::GrayscaleGeodesicDilateImageFilter< CharImageType, CharImageType >								CharGrayscaleGeodesicDilateImageFilterType;
typedef itk::GrayscaleGeodesicDilateImageFilter< FloatImageType, FloatImageType >							FloatGrayscaleGeodesicDilateImageFilterType;

//Segmentation
#include "itkBinaryThresholdImageFilter.h"
typedef itk::BinaryThresholdImageFilter< CharImageType, CharImageType >			CharBinaryThresholdImageFilterType;
typedef itk::BinaryThresholdImageFilter< FloatImageType, FloatImageType >		FloatBinaryThresholdImageFilterType;
typedef itk::BinaryThresholdImageFilter< FloatImageType, CharImageType >		FloatToCharThresholdingFilterType;


#include "itkConnectedComponentImageFilter.h"
typedef itk::ConnectedComponentImageFilter< CharImageType, CharImageType >		CharConnectedComponentImageFilterType;
typedef itk::ConnectedComponentImageFilter< FloatImageType, FloatImageType >	FloatConnectedComponentImageFilterType;

#include <itkVectorConfidenceConnectedImageFilter.h>
typedef itk::VectorConfidenceConnectedImageFilter< RGBImageType, CharImageType > VectorConfidenceConnectedImageFilterType;

//Labels 
typedef unsigned long																LongLabelType;

#include "itkShapeLabelObject.h"
typedef itk::ShapeLabelObject< LongLabelType, Dimension >							ShapeLabelObjectType;

#include "itkLabelMap.h"
typedef itk::LabelMap< ShapeLabelObjectType >										ShapeLabelMapType;
typedef ShapeLabelMapType::Pointer													ShapeLabelMapPointer;

#include "itkShapeLabelMapFilter.h"
typedef itk::ShapeLabelMapFilter< ShapeLabelMapType >								ShapeLabelMapFilterType;
typedef itk::ShapeOpeningLabelMapFilter< ShapeLabelMapType >						ShapeOpeningLabelMapFilterType;

#include "itkLabelMapOverlayImageFilter.h"

#include "itkLabelMapToBinaryImageFilter.h"
typedef itk::LabelMapToBinaryImageFilter< ShapeLabelMapType, CharImageType >		LabelMapToBinaryImageFilterType;

#include "itkBinaryImageToShapeLabelMapFilter.h"
typedef itk::BinaryImageToShapeLabelMapFilter< CharImageType, ShapeLabelMapType >	BinaryImageToShapeLabelMapFilterType;

#include "itkBinaryImageToLabelMapFilter.h"
typedef itk::BinaryImageToLabelMapFilter< CharImageType >							BinaryImageToLabelMapFilterType;

#include "itkLabelMapToBinaryImageFilter.h"
typedef itk::LabelMapToBinaryImageFilter< 
	BinaryImageToShapeLabelMapFilterType::OutputImageType, CharImageType >			LabelMapToBinaryImageFilterType;

#include "itkBinaryShapeOpeningImageFilter.h"
typedef itk::BinaryShapeOpeningImageFilter < CharImageType >						BinaryShapeOpeningImageFilterType;

#include "itkStatisticsLabelObject.h"
typedef itk::StatisticsLabelObject< LongLabelType, Dimension >						StatisticsLabelObjectType;
typedef itk::LabelMap< StatisticsLabelObjectType >									StatisticsLabelMapType;
typedef StatisticsLabelMapType::Pointer												StatisticsLabelMapPointer;

#include "itkStatisticsLabelMapFilter.h"
typedef itk::StatisticsLabelMapFilter< CharImageType, CharImageType >				StatisticsLabelMapFilterType;

#include "itkBinaryImageToStatisticsLabelMapFilter.h"
typedef itk::BinaryImageToStatisticsLabelMapFilter< CharImageType, CharImageType, StatisticsLabelMapType > BinaryImageToStatisticsLabelMapFilterType;

//Features
typedef itk::Neighborhood< CharPixelType, Dimension >								NeighborhoodType;

#include "itkStatisticsImageFilter.h"
typedef itk::StatisticsImageFilter< CharImageType >									StatisticsImageFilterType; 

#include "itkScalarImageToRunLengthFeaturesFilter.h"
typedef itk::Statistics::DenseFrequencyContainer2		HistogramFrequencyContainerType;
typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter <CharImageType, HistogramFrequencyContainerType> RunLengthFilterType;

#include "itkScalarImageToCooccurrenceMatrixFilter.h"
typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter < CharImageType >	ScalarImageToCooccurrenceMatrixFilterType;
typedef ScalarImageToCooccurrenceMatrixFilterType::HistogramType					HistogramType;

#include "itkHistogramToTextureFeaturesFilter.h"
typedef itk::Statistics::HistogramToTextureFeaturesFilter< HistogramType >			HistogramToTextureFeaturesFilterType;

#include "itkGaussianDistribution.h"

#endif /* __ITK_DECLARATIONS_H__ */
