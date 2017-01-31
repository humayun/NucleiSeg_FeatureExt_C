/*
 * RGBChannelExtraction.h
 *
 *  Created on: 20 May 2014
 *      Author: Humayun
 */
#include "RGBChannelExtraction.h"

CharImagePointer RedChannelExtraction(RGBImagePointer rgbImage) {

	typedef itk::ImageAdaptor<RGBImageType, RedChannelPixelAccessor> ImageAdaptorType;
	ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	adaptor->SetImage(rgbImage);

	typedef itk::RescaleIntensityImageFilter<ImageAdaptorType, CharImageType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(adaptor);
	rescaler->Update();

	return rescaler->GetOutput();
}

CharImagePointer BlueChannelExtraction(RGBImagePointer rgbImage) {

	typedef itk::ImageAdaptor<RGBImageType, BlueChannelPixelAccessor> ImageAdaptorType;
	ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	adaptor->SetImage(rgbImage);

	typedef itk::RescaleIntensityImageFilter<ImageAdaptorType, CharImageType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(adaptor);
	rescaler->Update();

	return rescaler->GetOutput();
}

CharImagePointer GreenChannelExtraction(RGBImagePointer rgbImage) {

	typedef itk::ImageAdaptor<RGBImageType, GreenChannelPixelAccessor> ImageAdaptorType;
	ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	adaptor->SetImage(rgbImage);

	typedef itk::RescaleIntensityImageFilter<ImageAdaptorType, CharImageType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(adaptor);
	rescaler->Update();

	return rescaler->GetOutput();
}

CharImagePointer BlueRatioExtraction(RGBImagePointer rgbImage) {

	typedef itk::ImageAdaptor<RGBImageType, BlueRatioPixelAccessor> ImageAdaptorType;
	ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	adaptor->SetImage(rgbImage);

	typedef itk::RescaleIntensityImageFilter<ImageAdaptorType, CharImageType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(adaptor);
	rescaler->Update();

	return rescaler->GetOutput();
}

CharImagePointer BlueRatio(RGBImagePointer rgbImage) {
	RGBImageSizeType size = rgbImage->GetLargestPossibleRegion().GetSize();
	CharImagePointer charImage = CharImageType::New();
	CreateImage(charImage, size);

	RGBImageRegionType region = rgbImage->GetLargestPossibleRegion();
	RGBImageLinearConstIteratorWithIndexType it(rgbImage, region);
	CharImageLinearIteratorWithIndexType it2(charImage, region);

	for (it.GoToBegin(); !it.IsAtEnd(); it.NextLine()) {
		while (!it.IsAtEndOfLine()) {
			CharImageIndexType index = it.GetIndex();
			RGBPixelType rgbPixel = it.Get();
			CharPixelType charPixel = (100 * rgbPixel.GetBlue())
					/ (1 + rgbPixel.GetRed() + rgbPixel.GetGreen())
					* (256
							/ (1 + rgbPixel.GetRed() + rgbPixel.GetGreen()
									+ rgbPixel.GetBlue()));
			charImage->SetPixel(index, charPixel);
			++it;
		}
	}

	CharToCharRescaleIntensityImageFilter::Pointer rescaler =
			CharToCharRescaleIntensityImageFilter::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(charImage);
	rescaler->Update();

	return rescaler->GetOutput();
}

CharImagePointer HematoxylinChannelExtraction( RGBImagePointer rgbImage, int matrixType ) {

	typedef itk::ColorDeconvolutionImageFilter<RGBImageType, CharImageType> DeconvolutionType;
	DeconvolutionType::Pointer deconv = DeconvolutionType::New();
	deconv->SetInput(rgbImage);
	switch (matrixType) {
	case 1: // Antoine PSL images
		deconv->SetPSL_HEStaining();
		break;
	case 2: // Antoine NUH images
		deconv->SetNUH_HEStaining();
		break;
	default:
		deconv->SetPSL_HEStaining();
	}

	deconv->Update();
	return deconv->GetOutput1();
}

CharImagePointer EosinChannelExtraction( RGBImagePointer rgbImage, int matrixType ) {

	typedef itk::ColorDeconvolutionImageFilter<RGBImageType, CharImageType> DeconvolutionType;
	DeconvolutionType::Pointer deconv = DeconvolutionType::New();
	deconv->SetInput(rgbImage);
	switch (matrixType) {
	case 1: // Antoine PSL images
		deconv->SetPSL_HEStaining();
		break;
	case 2: // Antoine NUH images
		deconv->SetNUH_HEStaining();
		break;
	default:
		deconv->SetPSL_HEStaining();
	}

	deconv->Update();
	return deconv->GetOutput2();
}
