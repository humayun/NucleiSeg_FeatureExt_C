/*
 * ITKFunctions.h
 *
 *  Created on: 20 May 2014
 *      Author: Humayun
 */
#ifndef __ITKFUNCTIONS_H__
#define __ITKFUNCTIONS_H__

#include "ITKDeclarations.h"

/************************************ Basic Filters and Functions **********************************************/
void CreateImage(CharImagePointer inImage, CharImageSizeType size);

#include "itkCastImageFilter.h"
typedef itk::CastImageFilter< CharImageType, FloatImageType >			CharToFloatCastImageFilterType;
typedef itk::CastImageFilter< FloatImageType, CharImageType >			FloatToCharCastImageFilterType;
typedef itk::CastImageFilter< RGBImageType, RGBImageType >				RGBToRGBCastImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP CastFilter( TInImage *inImage );

#include "itkSubtractImageFilter.h"
typedef itk::SubtractImageFilter< CharImageType, CharImageType >		CharSubtractImageFilterType;
typedef itk::SubtractImageFilter< FloatImageType, FloatImageType >		FloatSubtractImageFilterType;
template <class TImage, class TOutImageP>
TOutImageP SubtractFilter( TImage * inImage1, TImage * inImage2 );

#include "itkImageToHistogramFilter.h"
template <class TImage>
void ComputeHistogram(TImage * inImage, std::string histoFileNameCString );

#include "itkMaskedImageToHistogramFilter.h"
template <class TImage>
void ComputeHistogram(TImage * inImage, CharImagePointer mask, std::string histoFileNameCString );

#include "itkShrinkImageFilter.h"
typedef itk::ShrinkImageFilter< CharImageType, CharImageType >		CharShrinkImageFilterType;
typedef itk::ShrinkImageFilter< FloatImageType, FloatImageType >	FloatShrinkImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP ShrinkingImage( TInImage * inImage );

#include "itkInvertIntensityImageFilter.h"
typedef itk::InvertIntensityImageFilter< CharImageType >			CharInvertIntensityImageFilterType;
typedef itk::InvertIntensityImageFilter< FloatImageType >			FloatInvertIntensityImageFilterType;
template <class TImage, class TOutImageP>
TOutImageP InvertIntensityFilter( TImage* inImage, int forePixel = 255 );

#include "itkScalarToFractalImageFilter.h"
typedef itk::ScalarToFractalImageFilter< CharImageType, CharImageType >					ScalarToFractalImageFilterType;
CharImagePointer GetFractalImage( CharImagePointer inImage, CharImagePointer maskImage );

/**************************** Image Enhancmenet, Smoothing and Denoisification ************************************/

#include "itkMeanImageFilter.h"
typedef itk::MeanImageFilter< CharImageType, CharImageType >		CharMeanImageFilterType;
typedef itk::MeanImageFilter< FloatImageType, FloatImageType >		FloatMeanImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP MeanFilter( TInImage* inImage, double radius = 2 );

#include "itkMedianImageFilter.h"
typedef itk::MedianImageFilter< CharImageType, CharImageType >		CharMedianImageFilterType;
typedef itk::MedianImageFilter< FloatImageType, FloatImageType >	FloatMedianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP MedianFilter( TInImage* inImage, double radius = 2 );

#include "itkBinaryMedianImageFilter.h"
typedef itk::BinaryMedianImageFilter< CharImageType, CharImageType > BinaryMedianImageFilterType;
CharImagePointer BinaryMedianFilter( CharImagePointer inImage, double radius = 2 );

#include "itkBilateralImageFilter.h"
typedef itk::BilateralImageFilter< CharImageType, CharImageType >	CharBilateralImageFilterType;
typedef itk::BilateralImageFilter< FloatImageType, FloatImageType >	FloatBilateralImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP BilateralFilter(TInImage* inImage, double domainSigma = 1, double rangeSigma = 1 );

#include "itkCurvatureFlowImageFilter.h"
typedef itk::CurvatureFlowImageFilter< CharImageType, CharImageType >	CharCurvatureFlowImageFilterType;
typedef itk::CurvatureFlowImageFilter< FloatImageType, FloatImageType >	FloatCurvatureFlowImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP CurvatureFlowFilter(TInImage* inImage );

#include "itkDiscreteGaussianImageFilter.h"
typedef itk::DiscreteGaussianImageFilter< CharImageType, CharImageType >	CharDiscreteGaussianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP DiscreateGaussianFilter(TInImage* inImage, double variance = 1 );

#include "itkSmoothingRecursiveGaussianImageFilter.h"
typedef itk::SmoothingRecursiveGaussianImageFilter< CharImageType, CharImageType >	CharSmoothingRecursiveGaussianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP SmoothingRecursiveGaussianFilter(TInImage* inImage, double sigma = 1 );

#include "itkAnisotropicDiffusionImageFilter.h"
typedef itk::AnisotropicDiffusionImageFilter< CharImageType, CharImageType >			CharAnisotropicDiffusionImageFilterType;
typedef itk::AnisotropicDiffusionImageFilter< FloatImageType, FloatImageType >			FloatAnisotropicDiffusionImageFilterType;

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
typedef itk::CurvatureAnisotropicDiffusionImageFilter< CharImageType, CharImageType >	CharCurvatureAnisotropicDiffusionImageFilterType;
typedef itk::CurvatureAnisotropicDiffusionImageFilter< FloatImageType, FloatImageType >	FloatCurvatureAnisotropicDiffusionImageFilterType;

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP CurvatureAnisotropicDiffusionFilter( TInImage* inImage, double conductance = 9.0 );

/************************ Image Gradient, Edge Detection, Derivatives, Speed, Height ***********************************/
#include "itkRecursiveGaussianImageFilter.h"
typedef itk::RecursiveGaussianImageFilter< CharImageType, CharImageType >	CharRecursiveGaussianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP RecursiveGaussianFilter(TInImage* inImage, double direction = 0 );

#include "itkGradientAnisotropicDiffusionImageFilter.h"
typedef itk::GradientAnisotropicDiffusionImageFilter< CharImageType, CharImageType >	CharGradientAnisotropicDiffusionImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP GradientAnisotropicDiffusionFilter(TInImage* inImage, double conductance = 1.5 );

#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
typedef itk::VectorGradientAnisotropicDiffusionImageFilter< RGBImageType, RGBImageType > VectorGradientAnisotropicDiffusionImageFilterType;

#include "itkLaplacianSharpeningImageFilter.h"
typedef itk::LaplacianSharpeningImageFilter< FloatImageType, FloatImageType >		FloatLaplacianSharpeningImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP LaplacianSharpeningFilter( TInImage* inImage );

#include "itkLaplacianRecursiveGaussianImageFilter.h"
typedef itk::LaplacianRecursiveGaussianImageFilter< CharImageType, CharImageType >	CharLaplacianRecursiveGaussianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP LaplacianRecursiveGaussianFilter( TInImage* inImage );

// Speed
#include "itkSigmoidImageFilter.h"
typedef itk::SigmoidImageFilter< FloatImageType, FloatImageType >					FloatSigmoidImageFilterType;

#include "itkGradientMagnitudeImageFilter.h"
typedef itk::GradientMagnitudeImageFilter< FloatImageType, FloatImageType >			FloatGradientMagnitudeImageFilterType;

#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< CharImageType, CharImageType >	CharGradientMagnitudeRecursiveGaussianImageFilterType;
typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< FloatImageType, FloatImageType > FloatGradientMagnitudeRecursiveGaussianImageFilterType;

/***************************** Segmentation **********************************/
#include "itkScalarConnectedComponentImageFilter.h"
typedef itk::ScalarConnectedComponentImageFilter< CharImageType, CharImageType >	ScalarConnectedComponentImageFilterType;
CharImagePointer ScalarConnectedComponentFilter( CharImagePointer inImage );

#include "itkConfidenceConnectedImageFilter.h"
typedef itk::ConfidenceConnectedImageFilter< FloatImageType, FloatImageType >	ConfidenceConnectedImageFilterType;
CharImagePointer ConfidenceConnectedSegmentation( FloatImagePointer inImage, CharImageIndexType centroid, int rad = 3 );

#include "itkFastMarchingImageFilter.h"
typedef  itk::FastMarchingImageFilter< CharImageType, CharImageType >			CharFastMarchingFilterType;
typedef  itk::FastMarchingImageFilter< FloatImageType, FloatImageType >			FloatFastMarchingFilterType;
CharImagePointer FastMarchingSegmentation( CharImagePointer inImage, CharImageIndexType centroid );

#include "itkLevelSetBasedCellSegmentation.h"
typedef itk::LevelSetBasedCellSegmentation< CharImageType, CharImageType >		CharLevelSetBasedCellSegmentationType;
typedef itk::LevelSetBasedCellSegmentation< FloatImageType, FloatImageType >    FloatLevelSetBasedCellSegmentationType;
CharImagePointer LevelSetBasedCellSegmentation( CharImagePointer inImage, std::vector< CharImageIndexType > &centroids );

/******************************** *******************************************************/
#include "itkApproximateSignedDistanceMapImageFilter.h"
typedef itk::ApproximateSignedDistanceMapImageFilter< CharImageType, FloatImageType >	ApproximateSignedDistanceMapImageFilterType;
CharImagePointer ApproximateSignedDistanceMapFilter( CharImagePointer inImage );

#include "itkLabelContourImageFilter.h"
typedef itk::LabelContourImageFilter< CharImageType, CharImageType >					CharLabelContourImageFilterType;
CharImagePointer LabelContourFilter( CharImagePointer inImage );

#include "itkConnectedComponentImageFilter.h"
typedef itk::ConnectedComponentImageFilter< CharImageType, CharImageType >			CharConnectedComponentImageFilterType;
CharImagePointer ConnectedComponentWithLabelContourFilter( CharImagePointer inImage );

//Relabel region
#include "itkCustomColormapFunction.h"
typedef itk::Function::CustomColormapFunction< CharImagePixelType, RGBImagePixelType >	ColormapType;
#include "itkRelabelComponentImageFilter.h"
typedef itk::RelabelComponentImageFilter< CharImageType, CharImageType >			CharRelabelComponentImageFilterType;
#include "itkScalarToRGBColormapImageFilter.h"
typedef itk::ScalarToRGBColormapImageFilter< CharImageType, RGBImageType>			ColormapFilterType;
RGBImagePointer RelabelColormap( CharImagePointer inImage );

#include "itkScalarToRGBPixelFunctor.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkWatershedImageFilter.h"
/*
CharImagePointer WatershedSegmentation( CharImagePointer image, float lowerThreshold, float outputScaleLevel );
*/
#include "itkScalarImageKmeansImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkImageRegionIterator.h"
typedef itk::ScalarImageKmeansImageFilter< CharImageType > KMeansFilterType;
CharImagePointer KMeansImageClassifier( CharImagePointer inImage, double &meanA, double &meanB, double &meanC);

#include "itkKdTree.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkImageToListSampleAdaptor.h"
void KdTreeBasedKmeansEstimator( CharImagePointer inImage, double &meanA, double &meanB, double &meanC );

#endif // __ITKFUNCTIONS_H__