#include "ColorChannelExtraction.h"

CharImagePointer ColorChannelExtraction(RGBImagePointer rgbImage, std::string outputColorSpace, int colorChannel) {

	typedef itk::RescaleIntensityImageFilter<FloatImageType, CharImageType> FloatRescalerType;
	FloatRescalerType::Pointer floatRescaler = FloatRescalerType::New();
	floatRescaler->SetOutputMinimum(0);
	floatRescaler->SetOutputMaximum(255);

	if (outputColorSpace == std::string("XYZ")) {
		//std::cout << "Convert RGB to X, Y, and Z channels of XYZ color space." << std::endl;
		// Convert RGB image to XYZ color space
		typedef itk::Accessor::RGBToXYZColorSpacePixelAccessor<CharPixelType,
				FloatPixelType> RGBToXYZColorSpaceAccessorType;
		typedef itk::ImageAdaptor<RGBImageType, RGBToXYZColorSpaceAccessorType> RGBToXYZImageAdaptorType;
		RGBToXYZImageAdaptorType::Pointer rgbToXYZAdaptor =
				RGBToXYZImageAdaptorType::New();
		rgbToXYZAdaptor->SetImage(rgbImage);

		typedef itk::VectorIndexSelectionCastImageFilter<
				RGBToXYZImageAdaptorType, FloatImageType> VectorCastFilterType;
		VectorCastFilterType::Pointer vectorCastFilter =
				VectorCastFilterType::New();
		vectorCastFilter->SetInput(rgbToXYZAdaptor);

		// Selection of channel component
		vectorCastFilter->SetIndex(colorChannel);

		// Rescale float image to char image.
		floatRescaler->SetInput(vectorCastFilter->GetOutput());
		floatRescaler->Update();
	}
	if (outputColorSpace == std::string("HSI")) {
		//std::cout<<"Convert RGB to H, S, and I channels of HSI color space."<<std::endl;
		// Convert RGB image to HSI color space
		typedef itk::Accessor::RGBToHSIColorSpacePixelAccessor<CharPixelType,
				FloatPixelType> RGBToHSIColorSpacePixelAccessor;
		typedef itk::ImageAdaptor<RGBImageType, RGBToHSIColorSpacePixelAccessor> RGBToHSIImageAdaptorType;
		RGBToHSIImageAdaptorType::Pointer rgbToHSIAdaptor =
				RGBToHSIImageAdaptorType::New();
		rgbToHSIAdaptor->SetImage(rgbImage);

		typedef itk::VectorIndexSelectionCastImageFilter<
				RGBToHSIImageAdaptorType, FloatImageType> VectorCastFilterType;
		VectorCastFilterType::Pointer vectorCastFilter =
				VectorCastFilterType::New();
		vectorCastFilter->SetInput(rgbToHSIAdaptor);

		// Selection of channel component
		vectorCastFilter->SetIndex(colorChannel);

		// Rescale float image to char image.
		floatRescaler->SetInput(vectorCastFilter->GetOutput());
		floatRescaler->Update();
	}
	if (outputColorSpace == std::string("HSL")) {
		//std::cout<<"Convert RGB to H, S, and L channels of HSL color space."<<std::endl;
		// Convert RGB image to HSL color space
		typedef itk::Accessor::RGBToHSLColorSpacePixelAccessor<CharPixelType,
				FloatPixelType> RGBToHSLColorSpacePixelAccessor;
		typedef itk::ImageAdaptor<RGBImageType, RGBToHSLColorSpacePixelAccessor> RGBToHSLImageAdaptorType;
		RGBToHSLImageAdaptorType::Pointer rgbToHSLAdaptor =
				RGBToHSLImageAdaptorType::New();
		rgbToHSLAdaptor->SetImage(rgbImage);

		typedef itk::VectorIndexSelectionCastImageFilter<
				RGBToHSLImageAdaptorType, FloatImageType> VectorCastFilterType;
		VectorCastFilterType::Pointer vectorCastFilter =
				VectorCastFilterType::New();
		vectorCastFilter->SetInput(rgbToHSLAdaptor);

		// Selection of channel component
		vectorCastFilter->SetIndex(colorChannel);

		// Rescale float image to char image.
		floatRescaler->SetInput(vectorCastFilter->GetOutput());
		floatRescaler->Update();
	}
	if (outputColorSpace == std::string("HSV")) {
		//std::cout<<"Convert RGB to H, S, and V channels of HSV color space."<<std::endl;
		// Convert RGB image to HSV color space
		typedef itk::Accessor::RGBToHSVColorSpacePixelAccessor<CharPixelType,
				FloatPixelType> RGBToHSVColorSpacePixelAccessor;
		typedef itk::ImageAdaptor<RGBImageType, RGBToHSVColorSpacePixelAccessor> RGBToHSVImageAdaptorType;
		RGBToHSVImageAdaptorType::Pointer rgbToHSVAdaptor =
				RGBToHSVImageAdaptorType::New();
		rgbToHSVAdaptor->SetImage(rgbImage);

		typedef itk::VectorIndexSelectionCastImageFilter<
				RGBToHSVImageAdaptorType, FloatImageType> VectorCastFilterType;
		VectorCastFilterType::Pointer vectorCastFilter =
				VectorCastFilterType::New();
		vectorCastFilter->SetInput(rgbToHSVAdaptor);

		// Selection of channel component
		vectorCastFilter->SetIndex(colorChannel);

		// Rescale float image to char image.
		floatRescaler->SetInput(vectorCastFilter->GetOutput());
		floatRescaler->Update();
	}
	if (outputColorSpace == std::string("Lab")) {
		//std::cout<<"Convert RGB to L, a, and b channels of Lab color space."<<std::endl;
		// Convert RGB image to Lab color space
		typedef itk::Accessor::RGBToLabColorSpacePixelAccessor<CharPixelType,
				FloatPixelType> RGBToLabColorSpacePixelAccessor;
		typedef itk::ImageAdaptor<RGBImageType, RGBToLabColorSpacePixelAccessor> RGBToLabImageAdaptorType;
		RGBToLabImageAdaptorType::Pointer rgbToLabAdaptor =
				RGBToLabImageAdaptorType::New();
		rgbToLabAdaptor->SetImage(rgbImage);

		typedef itk::VectorIndexSelectionCastImageFilter<
				RGBToLabImageAdaptorType, FloatImageType> VectorCastFilterType;
		VectorCastFilterType::Pointer vectorCastFilter =
				VectorCastFilterType::New();
		vectorCastFilter->SetInput(rgbToLabAdaptor);

		// Selection of channel component
		vectorCastFilter->SetIndex(colorChannel);

		// Rescale float image to char image.
		floatRescaler->SetInput(vectorCastFilter->GetOutput());
		floatRescaler->Update();
	}
	if (outputColorSpace == std::string("Luv")) {
		//std::cout<<"Convert RGB to L, u, and v channels of Luv color space."<<std::endl;
		// Convert RGB image to Luv color space
		typedef itk::Accessor::RGBToLuvColorSpacePixelAccessor<CharPixelType,
				FloatPixelType> RGBToLuvColorSpacePixelAccessorType;
		typedef itk::ImageAdaptor<RGBImageType,
				RGBToLuvColorSpacePixelAccessorType> RGBToLuvImageAdaptorType;
		RGBToLuvImageAdaptorType::Pointer rgbToLuvAdaptor =
				RGBToLuvImageAdaptorType::New();
		rgbToLuvAdaptor->SetImage(rgbImage);

		typedef itk::VectorIndexSelectionCastImageFilter<
				RGBToLuvImageAdaptorType, FloatImageType> VectorCastFilterType;
		VectorCastFilterType::Pointer vectorCastFilter =
				VectorCastFilterType::New();
		vectorCastFilter->SetInput(rgbToLuvAdaptor);

		// Selection of channel component
		vectorCastFilter->SetIndex(colorChannel);

		// Rescale float image to char image.
		floatRescaler->SetInput(vectorCastFilter->GetOutput());
		floatRescaler->Update();
	}
	if (outputColorSpace == std::string("Yuv")) {
		//std::cout<<"Convert RGB to Y, u, and v channels of Yuv color space."<<std::endl;
		// Convert RGB image to Yuv color space
		typedef itk::Accessor::RGBToYuvColorSpacePixelAccessor<CharPixelType,
				FloatPixelType> RGBToYuvColorSpacePixelAccessorType;
		typedef itk::ImageAdaptor<RGBImageType,
				RGBToYuvColorSpacePixelAccessorType> RGBToYuvImageAdaptorType;
		RGBToYuvImageAdaptorType::Pointer rgbToYuvAdaptor =
				RGBToYuvImageAdaptorType::New();
		rgbToYuvAdaptor->SetImage(rgbImage);

		typedef itk::VectorIndexSelectionCastImageFilter<
				RGBToYuvImageAdaptorType, FloatImageType> VectorCastFilterType;
		VectorCastFilterType::Pointer vectorCastFilter =
				VectorCastFilterType::New();
		vectorCastFilter->SetInput(rgbToYuvAdaptor);

		// Selection of channel component
		vectorCastFilter->SetIndex(colorChannel);

		// Rescale float image to char image.
		floatRescaler->SetInput(vectorCastFilter->GetOutput());
		floatRescaler->Update();
	}
	if (outputColorSpace == std::string("YUV2")) {
		//std::cout<<"Convert RGB to Y, U, and V channels of YUV2 color space."<<std::endl;
		// Convert RGB image to YUV2 color space
		typedef itk::Accessor::RGBToYUV2ColorSpacePixelAccessor<CharPixelType,
				FloatPixelType> RGBToYUV2ColorSpacePixelAccessorType;
		typedef itk::ImageAdaptor<RGBImageType,
				RGBToYUV2ColorSpacePixelAccessorType> RGBToYUV2ImageAdaptorType;
		RGBToYUV2ImageAdaptorType::Pointer rgbToYUV2Adaptor =
				RGBToYUV2ImageAdaptorType::New();
		rgbToYUV2Adaptor->SetImage(rgbImage);

		typedef itk::VectorIndexSelectionCastImageFilter<
				RGBToYUV2ImageAdaptorType, FloatImageType> VectorCastFilterType;
		VectorCastFilterType::Pointer vectorCastFilter =
				VectorCastFilterType::New();
		vectorCastFilter->SetInput(rgbToYUV2Adaptor);

		// Selection of channel component
		vectorCastFilter->SetIndex(colorChannel);

		// Rescale float image to char image.
		floatRescaler->SetInput(vectorCastFilter->GetOutput());
		floatRescaler->Update();
	}

	return floatRescaler->GetOutput();
}

#define ITK_TEST_DIMENSION_MAX 6
int RegressionTestImage(std::string testImageFilename, std::string baselineImageFilename, int reportErrors, bool differences) {

	// Use the factory mechanism to read the test and baseline files and convert them to double
	typedef itk::Image<double, ITK_TEST_DIMENSION_MAX> ImageType;
	typedef itk::Image<unsigned char, ITK_TEST_DIMENSION_MAX> OutputType;
	typedef itk::Image<unsigned char, 2> DiffOutputType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	// Read the baseline file
	ReaderType::Pointer baselineReader = ReaderType::New();
	baselineReader->SetFileName(baselineImageFilename);
	try {
		baselineReader->UpdateLargestPossibleRegion();
	} catch (itk::ExceptionObject& e) {
		std::cerr << "Exception detected while reading "
				<< baselineImageFilename << " : " << e.GetDescription();
		return 1000;
	}

	// Read the file generated by the test
	ReaderType::Pointer testReader = ReaderType::New();
	testReader->SetFileName(testImageFilename);
	try {
		testReader->UpdateLargestPossibleRegion();
	} catch (itk::ExceptionObject& e) {
		std::cerr << "Exception detected while reading " << testImageFilename
				<< " : " << e.GetDescription() << std::endl;
		return 1000;
	}

	// The sizes of the baseline and test image must match
	ImageType::SizeType baselineSize;
	baselineSize =
			baselineReader->GetOutput()->GetLargestPossibleRegion().GetSize();
	ImageType::SizeType testSize;
	testSize = testReader->GetOutput()->GetLargestPossibleRegion().GetSize();

	if (baselineSize != testSize) {
		std::cerr
				<< "The size of the Baseline image and Test image do not match!"
				<< std::endl;
		std::cerr << "Baseline image: " << baselineImageFilename << " has size "
				<< baselineSize << std::endl;
		std::cerr << "Test image:     " << testImageFilename << " has size "
				<< testSize << std::endl;
		return 1;
	}

	// Now compare the two images
//	typedef itk::DifferenceImageFilter<ImageType,ImageType> DiffType;
	typedef itk::Testing::ComparisonImageFilter<ImageType, ImageType> DiffType;
	DiffType::Pointer diff = DiffType::New();
	diff->SetValidInput(baselineReader->GetOutput());
	diff->SetTestInput(testReader->GetOutput());
	diff->SetDifferenceThreshold(2.0);
	diff->UpdateLargestPossibleRegion();

	double status = diff->GetTotalDifference();

	if (reportErrors) {
		typedef itk::RescaleIntensityImageFilter<ImageType, OutputType> RescaleType;
		typedef itk::ExtractImageFilter<OutputType, DiffOutputType> ExtractType;
		typedef itk::ImageFileWriter<DiffOutputType> WriterType;
		typedef itk::ImageRegion<ITK_TEST_DIMENSION_MAX> RegionType;
		OutputType::IndexType index;
		index.Fill(0);
		OutputType::SizeType size;
		size.Fill(0);

		RescaleType::Pointer rescale = RescaleType::New();
		rescale->SetOutputMinimum(
				itk::NumericTraits<unsigned char>::NonpositiveMin());
		rescale->SetOutputMaximum(itk::NumericTraits<unsigned char>::max());
		rescale->SetInput(diff->GetOutput());
		rescale->UpdateLargestPossibleRegion();

		RegionType region;
		region.SetIndex(index);

		size = rescale->GetOutput()->GetLargestPossibleRegion().GetSize();
		for (unsigned int i = 2; i < ITK_TEST_DIMENSION_MAX; i++)
			size[i] = 0;
		region.SetSize(size);

		ExtractType::Pointer extract = ExtractType::New();
		extract->SetInput(rescale->GetOutput());
		extract->SetExtractionRegion(region);

		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(extract->GetOutput());
		if (differences) {
			// if there are discrepencies, create an diff image
			std::cout
					<< "<DartMeasurement name=\"ImageError\" type=\"numeric/double\">";
			std::cout << status;
			std::cout << "</DartMeasurement>" << std::endl;

			std::ostringstream diffName;
			diffName << testImageFilename << ".diff.png";
			try {
				rescale->SetInput(diff->GetOutput());
				rescale->Update();
			} catch (...) {
				std::cerr << "Error during rescale of " << diffName.str()
						<< std::endl;
			}
			writer->SetFileName(diffName.str().c_str());
			try {
				writer->Update();
			} catch (...) {
				std::cerr << "Error during write of " << diffName.str()
						<< std::endl;
			}

			std::cout
					<< "<DartMeasurementFile name=\"DifferenceImage\" type=\"image/png\">";
			std::cout << diffName.str();
			std::cout << "</DartMeasurementFile>" << std::endl;
		}

		std::ostringstream baseName;
		baseName << testImageFilename << ".base.png";
		try {
			rescale->SetInput(baselineReader->GetOutput());
			rescale->Update();
		} catch (...) {
			std::cerr << "Error during rescale of " << baseName.str()
					<< std::endl;
		}
		try {
			writer->SetFileName(baseName.str().c_str());
			writer->Update();
		} catch (...) {
			std::cerr << "Error during write of " << baseName.str()
					<< std::endl;
		}

		std::cout
				<< "<DartMeasurementFile name=\"BaselineImage\" type=\"image/png\">";
		std::cout << baseName.str();
		std::cout << "</DartMeasurementFile>" << std::endl;

		std::ostringstream testName;
		testName << testImageFilename << ".test.png";
		try {
			rescale->SetInput(testReader->GetOutput());
			rescale->Update();
		} catch (...) {
			std::cerr << "Error during rescale of " << testName.str()
					<< std::endl;
		}
		try {
			writer->SetFileName(testName.str().c_str());
			writer->Update();
		} catch (...) {
			std::cerr << "Error during write of " << testName.str()
					<< std::endl;
		}

		std::cout
				<< "<DartMeasurementFile name=\"TestImage\" type=\"image/png\">";
		std::cout << testName.str();
		std::cout << "</DartMeasurementFile>" << std::endl;
	}
	return (status != 0) ? 1 : 0;
}
