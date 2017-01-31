/*
 * ITKFunctions.cpp
 *
 *  Created on: 20 May 2014
 *      Author: Humayun
 */

#include "ITKFunctions.h"

/************************************ Basic Filters and Functions **********************************************/
void CreateImage(CharImagePointer inImage, CharImageSizeType size)
{
	CharImageRegionType region;
	CharImageIndexType start;
	start[0] = 0;
	start[1] = 0;

	region.SetSize(size);
	region.SetIndex(start);

	inImage->SetRegions(region);
	inImage->Allocate();

	//// Make a square
	//for(unsigned int r = 0; r < size[0]; r++)
	//{
	//	for(unsigned int c = 0; c < size[1]; c++)
	//	{
	//		CharImageIndexType pixelIndex;
	//		pixelIndex[0] = r;
	//		pixelIndex[1] = c;
	//		inImage->SetPixel(pixelIndex, 0);		// 255 = White image, 0 = Black image
	//	}
	//}
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP CastFilter( TInImage *inImage )
{
	typedef itk::CastImageFilter< TInImage, TOutImage> CastImageFilterType;
	typedef typename CastImageFilterType::Pointer	CastImageFilterPointer;
	CastImageFilterPointer cast = CastImageFilterType::New();
	cast->SetInput( inImage );
	cast->Update();
	return cast->GetOutput();
}

template <class TImage, class TOutImageP>
TOutImageP SubtractFilter( TImage * inImage1, TImage * inImage2 )
{
	typedef itk::SubtractImageFilter< TImage, TImage >	SubtractImageFilterType;
	typedef typename SubtractImageFilterType::Pointer SubtractImageFilterPointer;
	SubtractImageFilterPointer subtract = SubtractImageFilterType::New();
	subtract->SetInput1( inImage1 );
	subtract->SetInput2( inImage2 );
	subtract->Update();
	return subtract->GetOutput();
}

template <class TImage>
void ComputeHistogram(TImage * inImage, std::string histoFileNameCString )
{
	typedef itk::Statistics::ImageToHistogramFilter< TImage >				HistogramFilterType;
	typedef typename HistogramFilterType::Pointer							HistogramFilterPointer;
	typedef typename HistogramFilterType::HistogramMeasurementVectorType	HistogramMeasurementVectorType;
	typedef typename HistogramFilterType::HistogramSizeType					HistogramSizeType;
	typedef typename HistogramFilterType::HistogramType						HistogramType;
	typedef typename HistogramType::ConstIterator							HistogramConstIterator;

	HistogramSizeType histogramSize( 3 );
	histogramSize[0] = 2;  // number of bins for the Red   channel
	histogramSize[1] = 2;  // number of bins for the Green channel
	histogramSize[2] = 2;  // number of bins for the Blue  channel

	HistogramFilterPointer filter = HistogramFilterType::New();
	filter->SetInput( inImage );
	filter->SetAutoMinimumMaximum( true );
	filter->SetHistogramSize( histogramSize );
	filter->SetMarginalScale( 10 ); // Required (could this be set in the filter?)
	filter->Update();

	const HistogramType * histogram = filter->GetOutput();

	HistogramConstIterator histogramIterator = histogram->Begin();

	while( histogramIterator  != histogram->End() )
	{
	std::cout << "Index = " << histogram->GetIndex(histogramIterator.GetMeasurementVector())
			<< " Histogram cell center = " << histogramIterator.GetMeasurementVector()
			<< " Frequency = " << histogramIterator.GetFrequency() << std::endl;

	++histogramIterator ;
	}

	HistogramMeasurementVectorType mv(3);
	mv[0] = 255;
	mv[1] = 0;
	mv[2] = 0;
	std::cout << "Frequency = " << histogram->GetFrequency(histogram->GetIndex(mv)) << std::endl;
}

template <class TImage>
void ComputeHistogram(TImage * inImage, CharImagePointer mask, std::string histoFileNameCString )
{
	typedef itk::Statistics::MaskedImageToHistogramFilter< TImage, CharImageType >	MaskedHistogramFilterType;
	typedef typename MaskedHistogramFilterType::Pointer								MaskedHistogramFilterPointer;
	typedef typename MaskedHistogramFilterType::HistogramMeasurementVectorType		MaskedHistogramMeasurementVectorType;
	typedef typename MaskedHistogramFilterType::HistogramSizeType					MaskedHistogramSizeType;
	typedef typename MaskedHistogramFilterType::HistogramType						MaskedHistogramType;
	typedef typename MaskedHistogramType::ConstIterator								MaskedHistogramConstIterator;

	MaskedHistogramSizeType histogramSize( 3 );
	histogramSize[0] = 2;  // number of bins for the Red   channel
	histogramSize[1] = 2;  // number of bins for the Green channel
	histogramSize[2] = 2;  // number of bins for the Blue  channel

	MaskedHistogramFilterPointer filter = MaskedHistogramFilterType::New();
	filter->SetInput( inImage );
	filter->SetMaskImage( mask );
	filter->SetAutoMinimumMaximum( true );
	filter->SetHistogramSize( histogramSize );
	filter->SetMarginalScale( 10 ); // Required (could this be set in the filter?)
	filter->Update();

	const MaskedHistogramType * histogram = filter->GetOutput();

	MaskedHistogramConstIterator histogramIterator = histogram->Begin();

	while( histogramIterator  != histogram->End() )
	{
	std::cout << "Index = " << histogram->GetIndex(histogramIterator.GetMeasurementVector())
			<< " Histogram cell center = " << histogramIterator.GetMeasurementVector()
			<< " Frequency = " << histogramIterator.GetFrequency() << std::endl;

	++histogramIterator ;
	}

	MaskedHistogramMeasurementVectorType mv(3);
	mv[0] = 255;
	mv[1] = 0;
	mv[2] = 0;
	std::cout << "Frequency = " << histogram->GetFrequency(histogram->GetIndex(mv)) << std::endl;
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP ShrinkingImage( TInImage * inImage )
{
	typedef itk::ShrinkImageFilter< TInImage, TOutImage > ShrinkImageFilterType;
	typedef typename ShrinkImageFilterType::Pointer	ShrinkImageFilterPointer;
	ShrinkImageFilterPointer shrinkFilter = ShrinkImageFilterType::New();
	shrinkFilter->SetInput( inImage );
	shrinkFilter->SetShrinkFactor( 0, 2);
	shrinkFilter->SetShrinkFactor( 1, 2);
	shrinkFilter->Update();

	std::cout << "\n\nShrink ..."
		<< "\nOld Size: " << inImage->GetLargestPossibleRegion().GetSize()
		<< "\tSpacing: " << inImage->GetSpacing()
		<< "\nNew Size: " << shrinkFilter->GetOutput()->GetLargestPossibleRegion().GetSize()
		<< "\tSpacing: " << shrinkFilter->GetOutput()->GetSpacing();

	return shrinkFilter->GetOutput();
}

template <class TImage, class TOutImageP>
TOutImageP InvertIntensityFilter( TImage* inImage, int forePixel /* = 255 */ )
{
	typedef itk::InvertIntensityImageFilter< TImage > InvertIntensityImageFilterType;
	typedef typename InvertIntensityImageFilterType::Pointer InvertIntensityImageFilterPointer;
	InvertIntensityImageFilterPointer filter = InvertIntensityImageFilterType::New();
	filter->SetInput( inImage );
	filter->SetMaximum( forePixel );
	filter->Update();
	return filter->GetOutput();
}

/**************************** Image Enhancmenet, Smoothing and Denoisification ************************************/

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP MeanFilter( TInImage* inImage, double radius /* = 2 */ )
{
	typedef itk::MeanImageFilter< TInImage, TOutImage >	MeanImageFilterType;
	typedef typename MeanImageFilterType::Pointer	MeanImageFilterPointer;
	CharImageSizeType meanRadius;
	meanRadius.Fill( radius );
	MeanImageFilterPointer smoothing = MeanImageFilterType::New();
	smoothing->SetRadius( meanRadius );
	smoothing->SetInput( inImage );
	smoothing->Update();
	return smoothing->GetOutput();
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP MedianFilter( TInImage* inImage, double radius /* = 2 */ )
{
	typedef itk::MedianImageFilter< TInImage, TOutImage >	MedianImageFilterType;
	typedef typename MedianImageFilterType::Pointer MedianImageFilterPointer;
	CharImageSizeType medianRadius;
	medianRadius.Fill( radius );
	MedianImageFilterPointer smoothing = MedianImageFilterType::New();
	smoothing->SetRadius( medianRadius );
	smoothing->SetInput( inImage );
	smoothing->Update();
	return smoothing->GetOutput();
}

CharImagePointer BinaryMedianFilter( CharImagePointer inImage, double radius /* = 2 */ )
{
	CharImageSizeType medianRadius;
	medianRadius.Fill( radius );
	BinaryMedianImageFilterType::Pointer smoothing = BinaryMedianImageFilterType::New();
	smoothing->SetRadius( medianRadius );
	smoothing->SetInput( inImage );
	smoothing->Update();
	return smoothing->GetOutput();
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP BilateralFilter(TInImage* inImage, double domainSigma /* = 1 */, double rangeSigma /* = 1 */ )
{
	typedef itk::BilateralImageFilter< TInImage, TOutImage > BilateralImageFilterType;
	typedef typename BilateralImageFilterType::Pointer BilateralImageFilterPointer;
	BilateralImageFilterPointer smoothing = BilateralImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetDomainSigma( domainSigma );
	smoothing->SetRangeSigma( rangeSigma );
	smoothing->Update();
	return smoothing->GetOutput();
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP CurvatureFlowFilter(TInImage* inImage )
{
	typedef itk::CurvatureFlowImageFilter< TInImage, TOutImage > CurvatureFlowImageFilterType;
	typedef typename CurvatureFlowImageFilterType::Pointer CurvatureFlowImageFilterPointer;
	CurvatureFlowImageFilterPointer smoothing = CurvatureFlowImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetNumberOfIterations( 5 );
	smoothing->SetTimeStep( 0.125 );
	smoothing->Update();
	return smoothing->GetOutput();
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP DiscreateGaussianFilter(TInImage* inImage, double variance /* = 1 */ )
{
	typedef itk::DiscreteGaussianImageFilter< TInImage, TOutImage >	DiscreteGaussianImageFilterType;
	typedef typename DiscreteGaussianImageFilterType::Pointer DiscreteGaussianImageFilterPointer;
	DiscreteGaussianImageFilterPointer smoothing = DiscreteGaussianImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetVariance( variance );
	smoothing->Update();
	return smoothing->GetOutput();
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP SmoothingRecursiveGaussianFilter(TInImage* inImage, double sigma /* = 1 */ )
{
	typedef itk::SmoothingRecursiveGaussianImageFilter< TInImage, TOutImage >	SmoothingRecursiveGaussianImageFilterType;
	typedef typename SmoothingRecursiveGaussianImageFilterType::Pointer SmoothingRecursiveGaussianImageFilterPointer;
	SmoothingRecursiveGaussianImageFilterPointer smoothing = SmoothingRecursiveGaussianImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetSigma( sigma );
	smoothing->Update();
	return smoothing->GetOutput();
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP CurvatureAnisotropicDiffusionFilter( TInImage* inImage, double conductance /* = 9.0 */ )
{
	typedef itk::CurvatureAnisotropicDiffusionImageFilter< TInImage, TOutImage >	CurvatureAnisotropicDiffusionImageFilterType;
	typedef typename CurvatureAnisotropicDiffusionImageFilterType::Pointer CurvatureAnisotropicDiffusionImageFilterPointer;
	CurvatureAnisotropicDiffusionImageFilterPointer smoothing = CurvatureAnisotropicDiffusionImageFilterType::New();
	smoothing->SetTimeStep( 0.125 );
	smoothing->SetNumberOfIterations(  5 );
	smoothing->SetConductanceParameter( conductance );
	smoothing->SetInput( inImage );
	smoothing->Update();
	return smoothing->GetOutput();
}

/************************ Image Gradient, Edge Detection, Derivatives, Speed, Height ***********************************/

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP RecursiveGaussianFilter(TInImage* inImage, double direction /* = 0 */ )
{
	typedef itk::RecursiveGaussianImageFilter< TInImage, TOutImage > RecursiveGaussianImageFilterType;
	typedef typename RecursiveGaussianImageFilterType::Pointer RecursiveGaussianImageFilterPointer;
	RecursiveGaussianImageFilterPointer smoothing = RecursiveGaussianImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetDirection( direction );		// 0 = x-axis
	smoothing->SetSecondOrder();
	smoothing->Update();
	return smoothing->GetOutput();
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP GradientAnisotropicDiffusionFilter(TInImage* inImage, double conductance /* = 1.5 */ )
{
	typedef itk::GradientAnisotropicDiffusionImageFilter< TInImage, TOutImage >	GradientAnisotropicDiffusionImageFilterType;
	typedef typename GradientAnisotropicDiffusionImageFilterType::Pointer GradientAnisotropicDiffusionImageFilterPointer;
	GradientAnisotropicDiffusionImageFilterPointer smoothing = GradientAnisotropicDiffusionImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetConductanceParameter( conductance );
	smoothing->SetNumberOfIterations( 5 );
	smoothing->SetTimeStep( 0.125 );
	smoothing->Update();
	return smoothing->GetOutput();
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP LaplacianSharpeningFilter( TInImage* inImage )
{
	typedef itk::LaplacianSharpeningImageFilter< TInImage, TOutImage >		LaplacianSharpeningImageFilterType;
	typedef typename LaplacianSharpeningImageFilterType::Pointer LaplacianSharpeningImageFilterPointer;
	LaplacianSharpeningImageFilterPointer sharpening = LaplacianSharpeningImageFilterType::New();
	sharpening->SetInput( inImage );
	sharpening->Update();
	return sharpening->GetOutput();
}

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP LaplacianRecursiveGaussianFilter( TInImage* inImage )
{
	typedef itk::LaplacianRecursiveGaussianImageFilter< TInImage, TOutImage >	LaplacianRecursiveGaussianImageFilterType;
	typedef typename LaplacianRecursiveGaussianImageFilterType::Pointer LaplacianRecursiveGaussianImageFilterPointer;
	LaplacianRecursiveGaussianImageFilterPointer LoGFilter = LaplacianRecursiveGaussianImageFilterType::New();
	LoGFilter->SetInput( inImage );
	LoGFilter->Update();
	return LoGFilter->GetOutput();
}

/***************************** Segmentation **********************************/

CharImagePointer ScalarConnectedComponentFilter( CharImagePointer inImage )
{
	ScalarConnectedComponentImageFilterType::Pointer segment = ScalarConnectedComponentImageFilterType::New();
	segment->SetInput( inImage );
	try	{	segment->Update();	}
	catch (itk::ExceptionObject& excp)	{		std::cerr << "\n ScalarConnectedComponent Err 1:  Exception caught " << "\n" << excp << std::endl;	}
	return segment->GetOutput();
}

CharImagePointer ConfidenceConnectedSegmentation( FloatImagePointer inImage, CharImageIndexType centroid, int rad /* = 3 */ )
{
	ConfidenceConnectedImageFilterType::Pointer confidenceConnected = ConfidenceConnectedImageFilterType::New();
	confidenceConnected->SetInput( CurvatureFlowFilter<FloatImageType, FloatImageType, FloatImagePointer>( inImage ) );
	confidenceConnected->SetMultiplier( 2.5 );
	//confidenceConnected->SetTimeStep( 0.125 );
	confidenceConnected->SetNumberOfIterations( 5 );
	confidenceConnected->SetReplaceValue( 255 );
	confidenceConnected->SetInitialNeighborhoodRadius( rad );
	confidenceConnected->AddSeed( centroid );

	try	{	confidenceConnected->Update();	}
	catch (itk::ExceptionObject& excp)	{		std::cerr << "\n ConfidenceConnectedSegmentation Err 1:  Exception caught " << "\n" << excp << std::endl;	}

	return CastFilter<FloatImageType, CharImageType, CharImagePointer>( confidenceConnected->GetOutput() );
}

CharImagePointer LabelContourFilter( CharImagePointer inImage )
{
	CharLabelContourImageFilterType::Pointer labelContour = CharLabelContourImageFilterType::New();
	labelContour->SetInput( inImage );

	CharToCharRescaleIntensityImageFilter::Pointer rescaling = CharToCharRescaleIntensityImageFilter::New();
	rescaling->SetInput( labelContour->GetOutput() );
	rescaling->SetOutputMinimum( 0 );
	rescaling->SetOutputMaximum( 255 );

	rescaling->Update();

	return rescaling->GetOutput();
}

CharImagePointer ConnectedComponentWithLabelContourFilter( CharImagePointer inImage )
{
	// Label the boundary of segment region
	CharConnectedComponentImageFilterType::Pointer segmenting = CharConnectedComponentImageFilterType::New();
	segmenting->SetInput( inImage );
	segmenting->Update();

	CharLabelContourImageFilterType::Pointer labelContour = CharLabelContourImageFilterType::New();
	labelContour->SetInput( segmenting->GetOutput() );

	CharToCharRescaleIntensityImageFilter::Pointer rescaling = CharToCharRescaleIntensityImageFilter::New();
	rescaling->SetInput( labelContour->GetOutput() );
	rescaling->SetOutputMinimum( 0 );
	rescaling->SetOutputMaximum( 255 );

	rescaling->Update();

	return rescaling->GetOutput();
}

//Relabel region
RGBImagePointer RelabelColormap( CharImagePointer inImage )
{
	CharRelabelComponentImageFilterType::Pointer relabel = CharRelabelComponentImageFilterType::New();
	relabel->SetInput( inImage );

	ColormapType::Pointer largeColormap = ColormapType::New();
	// need to generate the color map for RGBColorMap
	//CreateRandomColormap( 255, largeColormap );

	ColormapFilterType::Pointer colorMapFilter = ColormapFilterType::New();

	colorMapFilter->SetInput( relabel->GetOutput() );
	colorMapFilter->SetColormap( largeColormap );
	colorMapFilter->Update();
	return colorMapFilter->GetOutput();
}
