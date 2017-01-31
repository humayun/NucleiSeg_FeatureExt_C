/*
 * ImageProcessing.cpp
 *
 *  Created on: 20 May 2014
 *      Author: Humayun
 */
#include "ImageProcessing.h"

CharImagePointer PatchExtraction( CharImagePointer inImage, CharImageIndexType centroid)
{
	CharImageSizeType imageSize = inImage->GetLargestPossibleRegion().GetSize();
	CharImageSizeType patchSize;
	patchSize.Fill( Config::getInstance()->getImagePatchSize() );
	CharImageIndexType start;

	if( centroid[0] < Config::getInstance()->getImagePatchSize()/2)
		start[0] = 0;
	else
		start[0] = centroid[0] - Config::getInstance()->getImagePatchSize()/2;
	if( start[0] + patchSize[0] > imageSize[0])
		start[0] = start[0] - (start[0] + patchSize[0] - imageSize[0]);

	if( centroid[1] < Config::getInstance()->getImagePatchSize()/2)
		start[1] = 0;
	else
		start[1] = centroid[1] - Config::getInstance()->getImagePatchSize()/2;
	if( start[1] + patchSize[1] > imageSize[1])
		start[1] = start[1] - (start[1] + patchSize[1] - imageSize[1]);

	CharImageRegionType patchRegion;
	patchRegion.SetIndex( start );
	patchRegion.SetSize( patchSize );
	CharRegionOfInterestImageFilterType::Pointer roiExtract = CharRegionOfInterestImageFilterType::New();
	roiExtract->SetRegionOfInterest( patchRegion );
	roiExtract->SetInput( inImage );
	roiExtract->Update();

	return roiExtract->GetOutput();
}

RGBImagePointer PatchExtraction( RGBImagePointer inImage, RGBImageIndexType centroid)
{
	RGBImageSizeType imageSize = inImage->GetLargestPossibleRegion().GetSize();
	RGBImageSizeType patchSize;
	patchSize.Fill( Config::getInstance()->getImagePatchSize() );
	RGBImageIndexType start;

	if( centroid[0] < Config::getInstance()->getImagePatchSize()/2)
		start[0] = 0;
	else
		start[0] = centroid[0] - Config::getInstance()->getImagePatchSize()/2;
	if( start[0] + patchSize[0] > imageSize[0])
		start[0] = start[0] - (start[0] + patchSize[0] - imageSize[0]);

	if( centroid[1] < Config::getInstance()->getImagePatchSize()/2)
		start[1] = 0;
	else
		start[1] = centroid[1] - Config::getInstance()->getImagePatchSize()/2;
	if( start[1] + patchSize[1] > imageSize[1])
		start[1] = start[1] - (start[1] + patchSize[1] - imageSize[1]);

	RGBImageRegionType patchRegion;
	patchRegion.SetIndex( start );
	patchRegion.SetSize( patchSize );
	RGBRegionOfInterestImageFilterType::Pointer roiExtract = RGBRegionOfInterestImageFilterType::New();
	roiExtract->SetRegionOfInterest( patchRegion );
	roiExtract->SetInput( inImage );
	roiExtract->Update();

	return roiExtract->GetOutput();
}

RGBImagePointer CheckandChangeInformationImage(RGBImagePointer rgbImage, CharImagePointer binImage){

	if (rgbImage->GetOrigin() == binImage->GetOrigin() || rgbImage->GetSpacing() == binImage->GetSpacing()){
		RGBChangeInformationImageFilterType::Pointer changeImageInfo = RGBChangeInformationImageFilterType::New();
		changeImageInfo->SetOutputOrigin(binImage->GetOrigin());
		changeImageInfo->ChangeOriginOn();
		changeImageInfo->SetOutputSpacing(binImage->GetSpacing());
		changeImageInfo->ChangeSpacingOn();
		changeImageInfo->SetInput(rgbImage);
		changeImageInfo->Update();
		return changeImageInfo->GetOutput();
	}
	return rgbImage;
}

RGBImagePointer DrawCircle(
		RGBImagePointer inImage,
		const std::vector< RGBImageIndexType > &centroids,
		int colorCode /* = Config::getInstance()->getOverlayColourCode() */,
		int radius /* = Config::getInstance()->getOverlayRadius() */ )
{
	RGBToRGBCastImageFilterType::Pointer copy = RGBToRGBCastImageFilterType::New();
	copy->SetInput( inImage );
	copy->Update();
	RGBImagePointer rgbImage = copy->GetOutput();

	RGBImageType::SizeType size = rgbImage->GetLargestPossibleRegion().GetSize();
	RGBImageType::PixelType pixelValue;
	RGBImageType::IndexType pixelIndex;
	switch( colorCode )
	{
	case 0:			// Red
		pixelValue.SetRed(255);
		pixelValue.SetGreen(0);
		pixelValue.SetBlue(0);
		break;
	case 1:			// Green ==> TP
		pixelValue.SetRed(0);
		pixelValue.SetGreen(255);
		pixelValue.SetBlue(0);
		break;
	case 2:			// Blue ==> FN
		pixelValue.SetRed(0);
		pixelValue.SetGreen(0);
		pixelValue.SetBlue(255);
		break;
	case 3:			// Yellow ==> FP
		pixelValue.SetRed(255);
		pixelValue.SetGreen(255);
		pixelValue.SetBlue(0);
		break;
	case 4:			// Cyan
		pixelValue.SetRed(0);
		pixelValue.SetGreen(255);
		pixelValue.SetBlue(255);
		break;
	case 5:			// Orange
		pixelValue.SetRed(255);
		pixelValue.SetGreen(128);
		pixelValue.SetBlue(0);
		break;
	default:		// Green
		pixelValue.SetRed(0);
		pixelValue.SetGreen(255);
		pixelValue.SetBlue(0);
	}
	for( int i=0; i<centroids.size(); i++ )
	{
		for( int j=0; j<360; j++ )
		{
			double angle = j*atan(1.0)/45.0;

            for( int thickness = 0; thickness < Config::getInstance()->getOverlayThickness(); ++thickness ) {
                pixelIndex[0] = (radius+thickness) * sin( angle );
                pixelIndex[1] = (radius+thickness) * cos( angle );

                pixelIndex[0] = pixelIndex[0] + centroids[i][0];
                pixelIndex[1] = pixelIndex[1] + centroids[i][1];

                if( (pixelIndex[0] >= 0 && pixelIndex[0] < size[0]) && (pixelIndex[1] >= 0 && pixelIndex[1] < size[1]) )
                    rgbImage->SetPixel( pixelIndex, pixelValue );
            }
		}
	}
	return rgbImage;
}

RGBImagePointer DrawSquare(
		RGBImagePointer inImage,
		const std::vector< RGBImageIndexType > &centroids,
		int colorCode /* = 1 */,
		int radius /* = Config::getInstance()->getOverlayRadius() */ )
{
	RGBToRGBCastImageFilterType::Pointer copy = RGBToRGBCastImageFilterType::New();
	copy->SetInput( inImage );
	copy->Update();
	RGBImagePointer rgbImage = copy->GetOutput();

	RGBImageSizeType size = rgbImage->GetLargestPossibleRegion().GetSize();
	RGBPixelType pixelValue;
	switch( colorCode )
	{
	case 0:						// Red
		pixelValue[0] = 255;
		pixelValue[1] = 0;
		pixelValue[2] = 0;
		break;
	case 1:						// Green
		pixelValue[0] = 0;
		pixelValue[1] = 255;
		pixelValue[2] = 0;
		break;
	case 2:						// Blue
		pixelValue[0] = 0;
		pixelValue[1] = 0;
		pixelValue[2] = 255;
		break;
	case 3:						// Yellow
		pixelValue[0] = 255;
		pixelValue[1] = 255;
		pixelValue[2] = 0;
		break;
	case 4:						// Cyan
		pixelValue[0] = 0;
		pixelValue[1] = 255;
		pixelValue[2] = 255;
		break;
	case 5:						// Orange
		pixelValue[0] = 255;
		pixelValue[1] = 128;
		pixelValue[2] = 0;
		break;
	default:					// Green
		pixelValue[0] = 255;
		pixelValue[1] = 255;
		pixelValue[2] = 0;
		break;
	}

	for( int i = 0; i < centroids.size(); i++ )
	{
		RGBImageIndexType startIndex = centroids[i];
		startIndex[0] -= radius;
		startIndex[1] -= radius;
		for( int r = 0; r < (2*radius+1); r++ )
		{
			RGBImageIndexType index;
			index[0] = startIndex[0] + r;
			index[1] = startIndex[1];
			for( int c = 0; c < (2*radius+1); c++ )
			{
				index[1] = startIndex[1] + c;
				if( index[0] > -1 && index[1] > -1 && index[0] < size[0] && index[1] < size[1] )
					rgbImage->SetPixel( index, pixelValue );
			}
		}
	}
	return rgbImage;
}

void ConvertIndexToVectorData(vector< CharImageIndexType > indexes1D, vector< vector< double > > &data2D)
{
	for (int i = 0; i<indexes1D.size(); i++)
	{
		vector<double> tmp;
		tmp.push_back(indexes1D[i][0]);
		tmp.push_back(indexes1D[i][1]);
		data2D.push_back(tmp);
	}
}

double ComputePDF(double mean, double variance, double value)
{
	itk::Statistics::GaussianDistribution::Pointer gaussian = itk::Statistics::GaussianDistribution::New();
	gaussian->SetMean(mean);
	gaussian->SetVariance(variance);
	return gaussian->EvaluatePDF(value);
}

void ComputeImageMeanVariance(CharImageType *inImage, double &mean, double &variance)
{
	StatisticsImageFilterType::Pointer stats = StatisticsImageFilterType::New();
	stats->SetInput( inImage );
	stats->Update();
	mean = stats->GetMean();
	variance = stats->GetVariance();
}

void ComputeImageMeanVarianceUsingMask( CharImageType *inImage, CharImageType *maskImage, double &mean, double &variance )
{
	itk::ImageRegionIterator<CharImageType> imageIterator(inImage, inImage->GetLargestPossibleRegion());
	itk::ImageRegionIterator<CharImageType> maskIterator(maskImage, maskImage->GetLargestPossibleRegion());
	double sum = 0.0, c = 0.0;
	while(!imageIterator.IsAtEnd())
	{
		if(maskIterator.Get()>0)
		{
			sum += imageIterator.Get();
			c++;
		}
		++imageIterator;
		++maskIterator;
	}
	mean = sum / c;
	imageIterator.GoToBegin();
	maskIterator.GoToBegin();
	double temp = 0;
	while(!imageIterator.IsAtEnd())
	{
		if(maskIterator.Get()>0)
			temp += pow(mean-imageIterator.Get(), 2);
		++imageIterator;
		++maskIterator;
	}
	variance = temp / c;
	double standardDeviation = sqrt( variance );
}

void ComputeCentroids( std::vector<std::vector<CharImageIndexType> > &indexes2D, std::vector<CharImageIndexType> &centroids)
{
	centroids.resize(0);
	for( int i = 0; i < indexes2D.size(); i++ )
	{
		long xSum = 0, ySum = 0;
		for( int j = 0; j < indexes2D[i].size(); j++ )
		{
			CharImageIndexType index = indexes2D[i][j];
			xSum += index[0];
			ySum += index[1];
		}
		CharImageIndexType centroid;
		centroid[0] = xSum / indexes2D[i].size();
		centroid[1] = ySum / indexes2D[i].size();
		centroids.push_back(centroid);
	}
}

CharImagePointer UpSamplingImage( CharImagePointer inImage, int UpSampleType /* = 1 */ )
	// UpSampleType,  1 for Mean, 2 for Median, 3 for Minimum, 4 for Maximum
{
	CharImagePointer upSampleImage = CharImageType::New();
	CharImageIndexType start;
	start.Fill(0);

	CharImageSizeType inImageSize = inImage->GetLargestPossibleRegion().GetSize();
	CharImageType::SizeType size;
	size[0] = inImageSize[0]/2;
	size[1] = inImageSize[1]/2;

	CharImageRegionType region( start, size );
	upSampleImage->SetRegions( region );
	upSampleImage->Allocate();
	upSampleImage->FillBuffer(0);

	for( int i = 0; i < size[0]; i++)
	{
		for ( int j = 0; j < size[1]; j++ )
		{
			CharImageIndexType ind, ind1, ind2, ind3, ind4;
			ind[0] = i;
			ind[1] = j;

			ind1[0] = i*2;
			ind1[1] = j*2;

			ind2[0] = i*2+1;
			ind2[1] = j*2;

			ind3[0] = i*2;
			ind3[1] = j*2+1;

			ind4[0] = i*2+1;
			ind4[1] = j*2+1;

			std::vector< CharImagePixelType > v;
			CharImagePixelType newValue;
			int avg=0;
			v.push_back( inImage->GetPixel( ind1 ) );
			if( ind2[0] <= inImageSize[0] )
				v.push_back( inImage->GetPixel( ind2 ) );
			if( ind3[1] <= inImageSize[1] )
				v.push_back( inImage->GetPixel( ind3 ) );
			if( ind4[0] <= inImageSize[0] && ind4[1] <= inImageSize[1] )
				v.push_back( inImage->GetPixel( ind4 ) );

			std::size_t vSize = v.size();
			std::sort( v.begin(), v.end() );

			switch( UpSampleType )
			{
			case 1:		// Upsampling based on Mean
				for(int k=0; k < vSize; k++)
					avg += v[k];
				avg /= vSize;
				newValue= avg;
				break;
			case 2:		// Upsampling based on Median
				if( vSize % 2 == 0 )
					newValue = ( v[ vSize / 2 - 1] + v[ vSize / 2 ] ) / 2;
				else
					newValue = v[ vSize / 2 ];
				break;
			case 3:		// Upsampling based on Minimum
				newValue = v[0];
				break;
			case 4:		// Upsampling based on Maximum
				newValue = v[vSize-1];
				break;
			default:
				for(int k=0; k < vSize; k++)
					avg += v[k];
				avg /= vSize;
				newValue= avg;
			}
			upSampleImage->SetPixel( ind, newValue );
		}
	}
	return upSampleImage;
}

void ComputeFeatures(CharImagePointer inImage, CharImagePointer maskImage, vector<double> &patchFeatures) {
	
	//inImage->SetOrigin(maskImage->GetOrigin());
	//inImage->SetSpacing(maskImage->GetSpacing());

	BinaryImageToStatisticsLabelMapFilterType::Pointer statisticsLabelMap = BinaryImageToStatisticsLabelMapFilterType::New();
	statisticsLabelMap->SetInput(maskImage);
	statisticsLabelMap->SetFeatureImage(inImage);
	statisticsLabelMap->Update();

	StatisticsLabelMapPointer patchLabelMap = statisticsLabelMap->GetOutput();

	long maxObjectSize = 0;
	int selectedLabelNo = 0;
	// First select the biggest object if number of objects (after segmentation) is more then one 1;
	for (unsigned int patchLabel = 1; patchLabel <= patchLabelMap->GetNumberOfLabelObjects(); patchLabel++) {
		const StatisticsLabelObjectType::Pointer patchLabelObject = patchLabelMap->GetLabelObject(patchLabel);
		if (patchLabelObject->GetPhysicalSize() > maxObjectSize) {
			selectedLabelNo = patchLabel;
			maxObjectSize = patchLabelObject->GetPhysicalSize();
		}
	}
	if (selectedLabelNo > 0 && maxObjectSize > 0) {

		const StatisticsLabelObjectType * patchLabelObject = patchLabelMap->GetLabelObject(selectedLabelNo);
		patchFeatures.push_back(boost::math::isnan(patchLabelObject->GetMean()) ? 0 : patchLabelObject->GetMean());
		patchFeatures.push_back(boost::math::isnan(patchLabelObject->GetMedian()) ? 0 : patchLabelObject->GetMedian());
		patchFeatures.push_back(boost::math::isnan(patchLabelObject->GetStandardDeviation()) ? 0 : patchLabelObject->GetStandardDeviation());
		patchFeatures.push_back(boost::math::isnan(patchLabelObject->GetKurtosis()) ? 0 : patchLabelObject->GetKurtosis());
		patchFeatures.push_back(boost::math::isnan(patchLabelObject->GetSkewness()) ? 0 : patchLabelObject->GetSkewness());
	}

	//GLCM
	NeighborhoodType CMNeighborhood;
	CMNeighborhood.SetRadius(Config::getInstance()->getCMRadius());
	const unsigned int CMCenterIndex = CMNeighborhood.GetCenterNeighborhoodIndex();
	vector < vector<double> > CMFeatures(CMCenterIndex, vector<double>(8, 0));

	CharOffsetType offset;
	for (unsigned int d = 0; d < CMCenterIndex; d++) {
		offset = CMNeighborhood.GetOffset(d);

		ScalarImageToCooccurrenceMatrixFilterType::Pointer glcmGenerator = ScalarImageToCooccurrenceMatrixFilterType::New();
		glcmGenerator->SetOffset(offset);
		glcmGenerator->SetNumberOfBinsPerAxis(Config::getInstance()->getGrayLevels());
		glcmGenerator->SetPixelValueMinMax(Config::getInstance()->getBackgroundPixel(), Config::getInstance()->getForegroundPixel());
		glcmGenerator->SetInput(inImage);

		// Computation of Regions based Texture Features 
		if (Config::getInstance()->regionFeatures()){
			glcmGenerator->SetMaskImage(maskImage);
			glcmGenerator->SetInsidePixelValue(Config::getInstance()->getForegroundPixel());
		}

		HistogramToTextureFeaturesFilterType::Pointer featureCalc = HistogramToTextureFeaturesFilterType::New();
		featureCalc->SetInput(glcmGenerator->GetOutput());
		featureCalc->Update();

		CMFeatures[d][0] =
			(boost::math::isnan(featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::Correlation))
			? 0 : featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::Correlation));
		CMFeatures[d][1] =
			(boost::math::isnan(featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::ClusterShade))
			? 0 : featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::ClusterShade));
		CMFeatures[d][2] =
			(boost::math::isnan(featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::ClusterProminence))
			? 0 : featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::ClusterProminence));
		CMFeatures[d][3] =
			(boost::math::isnan(featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::Energy))
			? 0 : featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::Energy));
		CMFeatures[d][4] =
			(boost::math::isnan(featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::Entropy))
			? 0 : featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::Entropy));
		CMFeatures[d][5] =
			(boost::math::isnan(featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::HaralickCorrelation))
			? 0 : featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::HaralickCorrelation));
		CMFeatures[d][6] =
			(boost::math::isnan(featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::Inertia))
			? 0 : featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::Inertia));
		CMFeatures[d][7] =
			(boost::math::isnan(featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::InverseDifferenceMoment))
			? 0 : featureCalc->GetFeature(HistogramToTextureFeaturesFilterType::InverseDifferenceMoment));
	}
	for (int c = 0; c < 8; c++) {
		double sum = 0;
		for (int r = 0; r < CMCenterIndex; r++)
			sum += (double)CMFeatures[r][c];
		patchFeatures.push_back(sum / CMCenterIndex);
	}

	// GLRL
	RunLengthFilterType::Pointer runLengthFilter = RunLengthFilterType::New();
	runLengthFilter->SetInput(inImage);
	runLengthFilter->SetPixelValueMinMax(Config::getInstance()->getBackgroundPixel(), Config::getInstance()->getForegroundPixel());
	runLengthFilter->SetNumberOfBinsPerAxis(Config::getInstance()->getGrayLevels());
	runLengthFilter->FastCalculationsOn();

	// Computation of Regions based Texture Features 
	if (Config::getInstance()->regionFeatures()){
		runLengthFilter->SetMaskImage(maskImage);
		runLengthFilter->SetInsidePixelValue(Config::getInstance()->getForegroundPixel());
	}

	// Set the Radius for Run-Length Matrix Calculation
	NeighborhoodType RLNeighborhood;
	RLNeighborhood.SetRadius(Config::getInstance()->getRLRadius());
	const unsigned int RLCenterIndex = RLNeighborhood.GetCenterNeighborhoodIndex();
	RunLengthFilterType::OffsetVectorPointer RLOffsets = RunLengthFilterType::OffsetVector::New();
	for (unsigned int d = 0; d < RLCenterIndex; d++)
		RLOffsets->push_back(RLNeighborhood.GetOffset(d));
	runLengthFilter->SetOffsets(RLOffsets);

	runLengthFilter->Update();

	RunLengthFilterType::FeatureValueVectorPointer means = runLengthFilter->GetFeatureMeans();
	const RunLengthFilterType::FeatureNameVector* names = runLengthFilter->GetRequestedFeatures();

	RunLengthFilterType::FeatureValueVector::ConstIterator mIt = means->Begin();
	while (mIt != means->End()) {
		patchFeatures.push_back((double)mIt.Value());
		++mIt;
	}
}

CharImagePointer SegmentationPipeLine(CharImagePointer inImage, CharImageIndexType seedPoint) {
	CharToFloatCastImageFilterType::Pointer charToFloat = CharToFloatCastImageFilterType::New();
	charToFloat->SetInput(inImage);

	// Denoise
	FloatCurvatureFlowImageFilterType::Pointer smoothing = FloatCurvatureFlowImageFilterType::New();
	smoothing->SetInput(charToFloat->GetOutput());
	smoothing->SetNumberOfIterations(5);
	smoothing->SetTimeStep(0.125);

	//CharGradientAnisotropicDiffusionImageFilterType::Pointer smoothing = CharGradientAnisotropicDiffusionImageFilterType::New();
	//smoothing->SetInput( inImage );
	//smoothing->SetConductanceParameter( 1.5 );
	//smoothing->SetNumberOfIterations( 5 );
	//smoothing->SetTimeStep( 0.125 );

	// Segment
	ConfidenceConnectedImageFilterType::Pointer confidenceConnectedFilter = ConfidenceConnectedImageFilterType::New();
	confidenceConnectedFilter->SetInitialNeighborhoodRadius(1);
	confidenceConnectedFilter->SetMultiplier(2);
	confidenceConnectedFilter->SetNumberOfIterations(10);
	confidenceConnectedFilter->SetReplaceValue(255);
	confidenceConnectedFilter->SetSeed(seedPoint);
	confidenceConnectedFilter->SetInput(smoothing->GetOutput());

	// Hole Filling
	CharVotingBinaryIterativeHoleFillingImageFilterType::Pointer holeFilling = CharVotingBinaryIterativeHoleFillingImageFilterType::New();
	CharStructuringElementType::RadiusType radius;
	radius.Fill(2);
	holeFilling->SetRadius(radius);
	holeFilling->SetMajorityThreshold(1);
	holeFilling->SetForegroundValue(Config::getInstance()->getForegroundPixel());
	holeFilling->SetBackgroundValue(Config::getInstance()->getBackgroundPixel());
	holeFilling->SetInput(CastFilter<FloatImageType, CharImageType, CharImagePointer>(confidenceConnectedFilter->GetOutput()));

	// Opening
	CharStructuringElementType openStructure;
	openStructure.SetRadius(1);
	openStructure.CreateStructuringElement();

	CharBinaryOpeningImageFilterType::Pointer opening = CharBinaryOpeningImageFilterType::New();
	opening->SetKernel(openStructure);
	opening->SetForegroundValue(Config::getInstance()->getForegroundPixel());
	opening->SetInput(holeFilling->GetOutput());

	// Closing
	CharBinaryClosingImageFilterType::Pointer closing = CharBinaryClosingImageFilterType::New();
	closing->SetKernel(openStructure);
	closing->SetForegroundValue(Config::getInstance()->getForegroundPixel());
	closing->SetInput(opening->GetOutput());

	closing->Update();

	return closing->GetOutput();
}

CharImagePointer SegmentationPipeLine(RGBImagePointer inImage, CharImageIndexType seedPoint) {
	// Denoise
	VectorGradientAnisotropicDiffusionImageFilterType::Pointer smoothing = VectorGradientAnisotropicDiffusionImageFilterType::New();
	smoothing->SetInput(inImage);

	// Segment
	VectorConfidenceConnectedImageFilterType::Pointer rgbConfidenceConnectedFilter = VectorConfidenceConnectedImageFilterType::New();
	rgbConfidenceConnectedFilter->SetInitialNeighborhoodRadius(2);
	rgbConfidenceConnectedFilter->SetMultiplier(3);
	rgbConfidenceConnectedFilter->SetNumberOfIterations(10);
	rgbConfidenceConnectedFilter->SetReplaceValue(Config::getInstance()->getForegroundPixel());
	rgbConfidenceConnectedFilter->SetSeed(seedPoint);
	rgbConfidenceConnectedFilter->SetInput(smoothing->GetOutput());

	// Hole Filling
	RGBImageSizeType holeFillingRadius;
	holeFillingRadius.Fill(2);
	CharVotingBinaryIterativeHoleFillingImageFilterType::Pointer holeFilling = CharVotingBinaryIterativeHoleFillingImageFilterType::New();
	holeFilling->SetRadius(holeFillingRadius);
	holeFilling->SetMajorityThreshold(1);
	holeFilling->SetForegroundValue(Config::getInstance()->getForegroundPixel());
	holeFilling->SetBackgroundValue(Config::getInstance()->getBackgroundPixel());
	holeFilling->SetInput(rgbConfidenceConnectedFilter->GetOutput());

	// Opening
	CharStructuringElementType openStructure;
	openStructure.SetRadius(1);
	openStructure.CreateStructuringElement();

	CharBinaryOpeningImageFilterType::Pointer opening = CharBinaryOpeningImageFilterType::New();
	opening->SetKernel(openStructure);
	opening->SetForegroundValue(Config::getInstance()->getForegroundPixel());
	opening->SetInput(holeFilling->GetOutput());

	// Closing
	CharBinaryClosingImageFilterType::Pointer closing = CharBinaryClosingImageFilterType::New();
	closing->SetKernel(openStructure);
	closing->SetForegroundValue(Config::getInstance()->getForegroundPixel());
	closing->SetInput(opening->GetOutput());

	closing->Update();

	return closing->GetOutput();
}