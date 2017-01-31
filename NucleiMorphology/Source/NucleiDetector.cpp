/*
 * NucleiDetector.cpp
 *
 *  Created on: 20 May 2014
 *      Author: Humayun
 */

#include "NucleiDetector.h"
//#include "QuickView.h"


void NucleiDetector::NucleiSegmentation(RGBImagePointer rgbImage, string imageName, path imagePath){
	m_ptrInputImage = rgbImage;
	m_strImageName = imageName;
	m_pathInputImageName = imageName;

	// Detect candidate region using selected channel
	m_ptrOutputImage = NucleiRegionDetection();
	m_ptrOutputImage->DisconnectPipeline();

	// Get all objects
	BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToLabelObjects->SetInput(m_ptrOutputImage);
	binaryImageToLabelObjects->Update();

	ShapeLabelMapPointer shapeLabelMap = binaryImageToLabelObjects->GetOutput();
	shapeLabelMap->DisconnectPipeline();

	unsigned int totalNuclei = shapeLabelMap->GetNumberOfLabelObjects();
	for (unsigned int nucleiNo = 1; nucleiNo <= totalNuclei; nucleiNo++)
	{
		BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects2 = BinaryImageToShapeLabelMapFilterType::New();
		binaryImageToLabelObjects2->SetInput(m_ptrOutputImage);
		binaryImageToLabelObjects2->Update();
		ShapeLabelMapPointer shapeLabelMap2 = binaryImageToLabelObjects2->GetOutput();
		shapeLabelMap2->DisconnectPipeline();

		const ShapeLabelObjectType * nuclei = shapeLabelMap2->GetLabelObject(nucleiNo);
		if (nuclei->GetPhysicalSize() > Config::getInstance()->getNucleiMinSize() && nuclei->GetPhysicalSize() < Config::getInstance()->getNucleiMaxSize())
		{
			for (unsigned int j = nucleiNo + 1; j <= totalNuclei; j++)
				shapeLabelMap2->RemoveLabel(shapeLabelMap2->GetLabelObject(j)->GetLabel());
			for (unsigned int j = 1; j < nucleiNo; j++)
				shapeLabelMap2->RemoveLabel(shapeLabelMap2->GetLabelObject(j)->GetLabel());

			LabelMapToBinaryImageFilterType::Pointer labelMapToBinaryImage = LabelMapToBinaryImageFilterType::New();
			labelMapToBinaryImage->SetInput(shapeLabelMap2);
			labelMapToBinaryImage->SetBackgroundValue(Config::getInstance()->getBackgroundPixel());
			labelMapToBinaryImage->SetForegroundValue(Config::getInstance()->getForegroundPixel());
			labelMapToBinaryImage->Update();

			CharImageIndexType centerIndex;
			centerIndex[0] = nuclei->GetCentroid()[0];
			centerIndex[1] = nuclei->GetCentroid()[1];
			m_vectorNucleiCentroids.push_back(centerIndex);

			CharImagePointer binaryImage = labelMapToBinaryImage->GetOutput();
			CharImagePointer binaryPatch = PatchExtraction(binaryImage, centerIndex);
			RGBImagePointer rgbPatch = PatchExtraction(rgbImage, centerIndex);

			if (Config::getInstance()->writeNucleiPatch() && imageName.length() > 0)
			{
				path newdir(imagePath);
				newdir /= imageName;
				create_directories(newdir);

				stringstream file1;
				file1 << nucleiNo << "_B" << Config::getInstance()->getImageExtension();
				path path1(newdir.make_preferred().string() / file1.str());
				CharFileWriting(binaryPatch, path1.make_preferred().string(), string("\n Char Image Writing - Exception caught\n"));

				stringstream file2;
				file2 << nucleiNo << "_RGB" << Config::getInstance()->getImageExtensionn();
				path path2(newdir.make_preferred().string() / file2.str());
				RGBFileWriting(rgbPatch, path2.make_preferred().string(), string("\n Char Image Writing - Exception caught\n"));
			}
		}
	}
}

void NucleiDetector::NucleiSegmentationFeatureExtraction(RGBImagePointer rgbImage, string name, path imagePath) {
	
	m_ptrInputImage = rgbImage;
	m_strImageName = name;
	m_pathInputImageName = imagePath;

	// Detect candidate region using selected channel
	m_ptrOutputImage = NucleiRegionDetection();
	m_ptrOutputImage->DisconnectPipeline();

	//Check the origin and spacing information of both RGB and binary image, if it is different then RGB image information is changed wrt to binary image
	m_ptrInputImage = CheckandChangeInformationImage(m_ptrInputImage, m_ptrOutputImage);

	// Get all objects
	BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToLabelObjects->SetInput( m_ptrOutputImage );
	binaryImageToLabelObjects->Update();

	ShapeLabelMapPointer shapeLabelMap = binaryImageToLabelObjects->GetOutput();
	shapeLabelMap->DisconnectPipeline();

	unsigned int totalNuclei = shapeLabelMap->GetNumberOfLabelObjects();
	for( unsigned int nucleiNo = 1; nucleiNo <= totalNuclei; nucleiNo++ )
	{
		BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects2 = BinaryImageToShapeLabelMapFilterType::New();
		binaryImageToLabelObjects2->SetInput( m_ptrOutputImage );
		binaryImageToLabelObjects2->Update();
		ShapeLabelMapPointer shapeLabelMap2 = binaryImageToLabelObjects2->GetOutput();
		shapeLabelMap2->DisconnectPipeline();

		const ShapeLabelObjectType * nuclei = shapeLabelMap2->GetLabelObject( nucleiNo );
		if (nuclei->GetPhysicalSize() > Config::getInstance()->getNucleiMinSize() && nuclei->GetPhysicalSize() <  Config::getInstance()->getNucleiMaxSize() )
		{
			for (unsigned int j = nucleiNo+1; j <= totalNuclei; j++ )
				shapeLabelMap2->RemoveLabel( shapeLabelMap2->GetLabelObject( j )->GetLabel() );
			for (unsigned int j = 1; j < nucleiNo; j++ )
				shapeLabelMap2->RemoveLabel( shapeLabelMap2->GetLabelObject( j )->GetLabel() );

			LabelMapToBinaryImageFilterType::Pointer labelMapToBinaryImage = LabelMapToBinaryImageFilterType::New();
			labelMapToBinaryImage->SetInput( shapeLabelMap2 );
			labelMapToBinaryImage->SetBackgroundValue( Config::getInstance()->getBackgroundPixel() );
			labelMapToBinaryImage->SetForegroundValue( Config::getInstance()->getForegroundPixel() );
			labelMapToBinaryImage->Update();

			CharImageIndexType centerIndex;
			centerIndex[0] = nuclei->GetCentroid()[0];
			centerIndex[1] = nuclei->GetCentroid()[1];
			m_vectorNucleiCentroids.push_back( centerIndex );

			CharImagePointer binaryImage = labelMapToBinaryImage->GetOutput();
			CharImagePointer binaryPatch = PatchExtraction( binaryImage, centerIndex );
			RGBImagePointer  rgbPatch	 = PatchExtraction( m_ptrInputImage, centerIndex );

			if( Config::getInstance()->writeNucleiPatch() && m_strImageName.length() > 0)
			{
				path newdir(m_pathInputImageName);
				newdir /= m_strImageName;
				create_directories( newdir );

				stringstream file1;
				file1 << nucleiNo << "_B" << Config::getInstance()->getImageExtension();
				path path1( newdir.make_preferred().string() / file1.str() );
				CharFileWriting( binaryPatch, path1.make_preferred().string(), string("\n Char Image Writing - Exception caught\n") );

				stringstream file2;
				file2 << nucleiNo << "_RGB" << Config::getInstance()->getImageExtensionn();
				path path2( newdir.make_preferred().string() / file2.str() );
				RGBFileWriting( rgbPatch, path2.make_preferred().string(), string("\n Char Image Writing - Exception caught\n") );
			}

			vector< double > objectFeatures;

			objectFeatures.push_back( nuclei->GetCentroid()[0] );
			objectFeatures.push_back( nuclei->GetCentroid()[1] );
			
			objectFeatures.push_back( nuclei->GetPhysicalSize() );
			objectFeatures.push_back( nuclei->GetRoundness() );
			objectFeatures.push_back( nuclei->GetElongation() );
			objectFeatures.push_back(nuclei->GetFlatness());
			objectFeatures.push_back(nuclei->GetPerimeter());
			objectFeatures.push_back( nuclei->GetEquivalentSphericalPerimeter() );
			objectFeatures.push_back(nuclei->GetEquivalentSphericalRadius());
			objectFeatures.push_back(nuclei->GetEquivalentEllipsoidDiameter()[0]);
			objectFeatures.push_back(nuclei->GetEquivalentEllipsoidDiameter()[1]);

			CharImagePointer redPatch	= RedChannelExtraction(   rgbPatch );
			CharImagePointer greenPatch = GreenChannelExtraction( rgbPatch );
			CharImagePointer bluePatch	= BlueChannelExtraction(  rgbPatch );
			CharImagePointer hsvPatch	= ColorChannelExtraction( rgbPatch, "HSV", 2 );
			CharImagePointer labPatch	= ColorChannelExtraction( rgbPatch, "Lab", 0 );
			//CharImagePointer luvPatch	= ColorChannelExtraction( rgbPatch, "Luv", 0 );
			CharImagePointer hePatch	= HematoxylinChannelExtraction( rgbPatch, Config::getInstance()->getHEMatrixType() );
			CharImagePointer brPatch	= BlueRatioExtraction(rgbPatch);

			ComputeFeatures( redPatch,   binaryPatch, objectFeatures);	// Red
			ComputeFeatures( greenPatch, binaryPatch, objectFeatures ); // Green
			ComputeFeatures( bluePatch,  binaryPatch, objectFeatures );	// Blue
			ComputeFeatures( hsvPatch,   binaryPatch, objectFeatures );	// HSV
			ComputeFeatures( labPatch,   binaryPatch, objectFeatures );	// Lab
			//ComputeFeatures( luvPatch,   binaryPatch, objectFeatures );	// Luv
			ComputeFeatures( hePatch,    binaryPatch, objectFeatures );	// H&E
			ComputeFeatures( brPatch,    binaryPatch, objectFeatures );	// BR

			for( int up = 2; up <= Config::getInstance()->getMacroPixelUpsamplingValue(); up++)
			{
				redPatch    = UpSamplingImage( redPatch, 1 );
				greenPatch  = UpSamplingImage( greenPatch, 1);
				bluePatch   = UpSamplingImage( bluePatch, 1 );
				hsvPatch    = UpSamplingImage( hsvPatch, 1 );
				labPatch    = UpSamplingImage( labPatch, 1 );
				//luvPatch    = UpSamplingImage( luvPatch, 1 );
				hePatch     = UpSamplingImage( hePatch, 1 );
				brPatch     = UpSamplingImage( brPatch, 1 );

				binaryPatch = UpSamplingImage(binaryPatch, 1);

				ComputeFeatures( redPatch,   binaryPatch, objectFeatures );	// Red
				ComputeFeatures( greenPatch, binaryPatch, objectFeatures ); // Green
				ComputeFeatures( bluePatch,  binaryPatch, objectFeatures );	// Blue
				ComputeFeatures( hsvPatch,   binaryPatch, objectFeatures );	// HSV
				ComputeFeatures( labPatch,   binaryPatch, objectFeatures );	// Lab
				//ComputeFeatures( luvPatch,   binaryPatch, objectFeatures );	// Luv
				ComputeFeatures( hePatch,    binaryPatch, objectFeatures );	// H&E
				ComputeFeatures( brPatch,    binaryPatch, objectFeatures );	// BR
			}

			m_vectorFeatures.push_back( objectFeatures );
		}
	}
}

void NucleiDetector::FeatureExtraction(RGBImagePointer rgbImage, CharImagePointer binImage, string name, path imagePath) {

	//Check the origin and spacing information of both RGB and binary image, if it is different then RGB image information is changed wrt to binary image
	rgbImage = CheckandChangeInformationImage(rgbImage, binImage);

	m_ptrInputImage = rgbImage;
	m_ptrOutputImage = binImage;
	m_strImageName = name;
	m_pathInputImageName = imagePath;

	// Get all objects
	BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToLabelObjects->SetInput(m_ptrOutputImage);
	binaryImageToLabelObjects->Update();

	ShapeLabelMapPointer shapeLabelMap = binaryImageToLabelObjects->GetOutput();
	shapeLabelMap->DisconnectPipeline();

	unsigned int totalNuclei = shapeLabelMap->GetNumberOfLabelObjects();
	for (unsigned int nucleiNo = 1; nucleiNo <= totalNuclei; nucleiNo++)
	{
		BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects2 = BinaryImageToShapeLabelMapFilterType::New();
		binaryImageToLabelObjects2->SetInput(m_ptrOutputImage);
		binaryImageToLabelObjects2->Update();
		ShapeLabelMapPointer shapeLabelMap2 = binaryImageToLabelObjects2->GetOutput();
		shapeLabelMap2->DisconnectPipeline();

		const ShapeLabelObjectType * nuclei = shapeLabelMap2->GetLabelObject(nucleiNo);
		if (nuclei->GetPhysicalSize() > Config::getInstance()->getNucleiMinSize() && nuclei->GetPhysicalSize() <  Config::getInstance()->getNucleiMaxSize())
		{
			CharImageIndexType centerIndex;
			centerIndex[0] = nuclei->GetCentroid()[0];
			centerIndex[1] = nuclei->GetCentroid()[1];
			m_vectorNucleiCentroids.push_back(centerIndex);

			for (unsigned int j = 1; j < nucleiNo; j++)
				shapeLabelMap2->RemoveLabel(shapeLabelMap2->GetLabelObject(j)->GetLabel());
			for (unsigned int j = nucleiNo + 1; j <= totalNuclei; j++)
				shapeLabelMap2->RemoveLabel(shapeLabelMap2->GetLabelObject(j)->GetLabel());

			LabelMapToBinaryImageFilterType::Pointer labelMapToBinaryImage = LabelMapToBinaryImageFilterType::New();
			labelMapToBinaryImage->SetInput(shapeLabelMap2);
			labelMapToBinaryImage->SetBackgroundValue(Config::getInstance()->getBackgroundPixel());
			labelMapToBinaryImage->SetForegroundValue(Config::getInstance()->getForegroundPixel());
			labelMapToBinaryImage->Update();

			CharImagePointer binaryImage = labelMapToBinaryImage->GetOutput();
			CharImagePointer binaryPatch = PatchExtraction(binaryImage, centerIndex);
			RGBImagePointer  rgbPatch = PatchExtraction(m_ptrInputImage, centerIndex);

			vector< double > objectFeatures;
			
			objectFeatures.push_back(nuclei->GetCentroid()[0]);
			objectFeatures.push_back(nuclei->GetCentroid()[1]);

			objectFeatures.push_back(nuclei->GetPhysicalSize());
			objectFeatures.push_back(nuclei->GetRoundness());
			objectFeatures.push_back(nuclei->GetElongation());
			objectFeatures.push_back(nuclei->GetFlatness());
			objectFeatures.push_back(nuclei->GetPerimeter());
			objectFeatures.push_back(nuclei->GetEquivalentSphericalPerimeter());
			objectFeatures.push_back(nuclei->GetEquivalentSphericalRadius());
			objectFeatures.push_back(nuclei->GetEquivalentEllipsoidDiameter()[0]);
			objectFeatures.push_back(nuclei->GetEquivalentEllipsoidDiameter()[1]);

			CharImagePointer redPatch = RedChannelExtraction(rgbPatch);
			CharImagePointer greenPatch = GreenChannelExtraction(rgbPatch);
			CharImagePointer bluePatch = BlueChannelExtraction(rgbPatch);
			CharImagePointer hsvPatch = ColorChannelExtraction(rgbPatch, "HSV", 2);
			CharImagePointer labPatch = ColorChannelExtraction(rgbPatch, "Lab", 0);
			//CharImagePointer luvPatch = ColorChannelExtraction(rgbPatch, "Luv", 0);
			CharImagePointer hePatch = HematoxylinChannelExtraction(rgbPatch, Config::getInstance()->getHEMatrixType());
			CharImagePointer brPatch = BlueRatioExtraction(rgbPatch);

			ComputeFeatures(redPatch, binaryPatch, objectFeatures);		// Red
			ComputeFeatures(greenPatch, binaryPatch, objectFeatures);	// Green
			ComputeFeatures(bluePatch, binaryPatch, objectFeatures);	// Blue
			ComputeFeatures(hsvPatch, binaryPatch, objectFeatures);		// HSV
			ComputeFeatures(labPatch, binaryPatch, objectFeatures);		// Lab
			//ComputeFeatures(luvPatch, binaryPatch, objectFeatures);		// Luv
			ComputeFeatures(hePatch, binaryPatch, objectFeatures);		// H&E
			ComputeFeatures(brPatch, binaryPatch, objectFeatures);		// BR

			for (int up = 2; up <= Config::getInstance()->getMacroPixelUpsamplingValue(); up++)
			{
				redPatch = UpSamplingImage(redPatch, 1);
				greenPatch = UpSamplingImage(greenPatch, 1);
				bluePatch = UpSamplingImage(bluePatch, 1);
				hsvPatch = UpSamplingImage(hsvPatch, 1);
				labPatch = UpSamplingImage(labPatch, 1);
				//luvPatch = UpSamplingImage(luvPatch, 1);
				hePatch = UpSamplingImage(hePatch, 1);
				brPatch = UpSamplingImage(brPatch, 1);

				binaryPatch = UpSamplingImage(binaryPatch, 1);

				ComputeFeatures(redPatch, binaryPatch, objectFeatures);	// Red
				ComputeFeatures(greenPatch, binaryPatch, objectFeatures); // Green
				ComputeFeatures(bluePatch, binaryPatch, objectFeatures); // Blue
				ComputeFeatures(hsvPatch, binaryPatch, objectFeatures);	// HSV
				ComputeFeatures(labPatch, binaryPatch, objectFeatures);	// Lab
				//ComputeFeatures(luvPatch, binaryPatch, objectFeatures);	// Luv
				ComputeFeatures(hePatch, binaryPatch, objectFeatures);	// H&E
				ComputeFeatures(brPatch, binaryPatch, objectFeatures);	// BR
			}

			m_vectorFeatures.push_back(objectFeatures);
		}
	}
}

void NucleiDetector::FeatureExtraction_LSM(RGBImagePointer rgbImage, CharImagePointer binImage, string name, path imagePath) {

	//Check the origin and spacing information of both RGB and binary image, if it is different then RGB image information is changed wrt to binary image
	rgbImage = CheckandChangeInformationImage(rgbImage, binImage);

	m_ptrInputImage = rgbImage;
	m_ptrOutputImage = binImage;
	m_strImageName = name;
	m_pathInputImageName = imagePath;

	// Get all objects
	BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToLabelObjects->SetInput(m_ptrOutputImage);
	binaryImageToLabelObjects->Update();

	ShapeLabelMapPointer shapeLabelMap = binaryImageToLabelObjects->GetOutput();
	shapeLabelMap->DisconnectPipeline();

	unsigned int totalNuclei = shapeLabelMap->GetNumberOfLabelObjects();
	for (unsigned int nucleiNo = 1; nucleiNo <= totalNuclei; nucleiNo++)
	{
		BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects2 = BinaryImageToShapeLabelMapFilterType::New();
		binaryImageToLabelObjects2->SetInput(m_ptrOutputImage);
		binaryImageToLabelObjects2->Update();
		ShapeLabelMapPointer shapeLabelMap2 = binaryImageToLabelObjects2->GetOutput();
		shapeLabelMap2->DisconnectPipeline();

		const ShapeLabelObjectType * nuclei = shapeLabelMap2->GetLabelObject(nucleiNo);
		if (nuclei->GetPhysicalSize() > Config::getInstance()->getNucleiMinSize() && nuclei->GetPhysicalSize() <  Config::getInstance()->getNucleiMaxSize())
		{
			CharImageIndexType centerIndex;
			centerIndex[0] = nuclei->GetCentroid()[0];
			centerIndex[1] = nuclei->GetCentroid()[1];
			m_vectorNucleiCentroids.push_back(centerIndex);

			for (unsigned int j = 1; j < nucleiNo; j++)
				shapeLabelMap2->RemoveLabel(shapeLabelMap2->GetLabelObject(j)->GetLabel());
			for (unsigned int j = nucleiNo + 1; j <= totalNuclei; j++)
				shapeLabelMap2->RemoveLabel(shapeLabelMap2->GetLabelObject(j)->GetLabel());

			LabelMapToBinaryImageFilterType::Pointer labelMapToBinaryImage = LabelMapToBinaryImageFilterType::New();
			labelMapToBinaryImage->SetInput(shapeLabelMap2);
			labelMapToBinaryImage->SetBackgroundValue(Config::getInstance()->getBackgroundPixel());
			labelMapToBinaryImage->SetForegroundValue(Config::getInstance()->getForegroundPixel());
			labelMapToBinaryImage->Update();

			CharImagePointer binaryImage = labelMapToBinaryImage->GetOutput();
			CharImagePointer binaryPatch = PatchExtraction(binaryImage, centerIndex);
			RGBImagePointer  rgbPatch = PatchExtraction(m_ptrInputImage, centerIndex);

			vector< double > objectFeatures;

			objectFeatures.push_back(nuclei->GetCentroid()[0]);
			objectFeatures.push_back(nuclei->GetCentroid()[1]);

			objectFeatures.push_back(nuclei->GetPhysicalSize());
			objectFeatures.push_back(nuclei->GetRoundness());
			objectFeatures.push_back(nuclei->GetElongation());
			objectFeatures.push_back(nuclei->GetFlatness());
			objectFeatures.push_back(nuclei->GetPerimeter());
			objectFeatures.push_back(nuclei->GetEquivalentSphericalPerimeter());
			objectFeatures.push_back(nuclei->GetEquivalentSphericalRadius());
			objectFeatures.push_back(nuclei->GetEquivalentEllipsoidDiameter()[0]);
			objectFeatures.push_back(nuclei->GetEquivalentEllipsoidDiameter()[1]);

			CharImagePointer greenPatch = GreenChannelExtraction(rgbPatch);

			ComputeFeatures(greenPatch, binaryPatch, objectFeatures);	// Green

			for (int up = 2; up <= Config::getInstance()->getMacroPixelUpsamplingValue(); up++)
			{
				greenPatch = UpSamplingImage(greenPatch, 1);
				binaryPatch = UpSamplingImage(binaryPatch, 1);
				ComputeFeatures(greenPatch, binaryPatch, objectFeatures); // Green
			}

			m_vectorFeatures.push_back(objectFeatures);
		}
	}
}

void NucleiDetector::FeatureExtraction_TCGA(RGBImagePointer rgbImage, CharImagePointer binImage, string name, path imagePath) {

	//Check the origin and spacing information of both RGB and binary image, if it is different then RGB image information is changed wrt to binary image
	rgbImage = CheckandChangeInformationImage(rgbImage, binImage);

	vector<string> nameToken;
	boost::split(nameToken, name, boost::is_any_of("_"));
	unsigned long frameX = 0, frameY = 0;
	if (nameToken.size() >= 2){
		frameX = string_to_double(nameToken[0]);
		frameY = string_to_double(nameToken[1]);
	}

	m_ptrInputImage = rgbImage;
	m_ptrOutputImage = binImage;
	m_strImageName = name;
	m_pathInputImageName = imagePath;

	// Get all objects
	BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToLabelObjects->SetInput(m_ptrOutputImage);
	binaryImageToLabelObjects->Update();

	ShapeLabelMapPointer shapeLabelMap = binaryImageToLabelObjects->GetOutput();
	shapeLabelMap->DisconnectPipeline();

	unsigned int totalNuclei = shapeLabelMap->GetNumberOfLabelObjects();
	for (unsigned int nucleiNo = 1; nucleiNo <= totalNuclei; nucleiNo++)
	{
		BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToLabelObjects2 = BinaryImageToShapeLabelMapFilterType::New();
		binaryImageToLabelObjects2->SetInput(m_ptrOutputImage);
		binaryImageToLabelObjects2->Update();
		ShapeLabelMapPointer shapeLabelMap2 = binaryImageToLabelObjects2->GetOutput();
		shapeLabelMap2->DisconnectPipeline();

		const ShapeLabelObjectType * nuclei = shapeLabelMap2->GetLabelObject(nucleiNo);
		if (nuclei->GetPhysicalSize() > Config::getInstance()->getNucleiMinSize() && nuclei->GetPhysicalSize() <  Config::getInstance()->getNucleiMaxSize())
		{
			CharImageIndexType centerIndex;
			centerIndex[0] = nuclei->GetCentroid()[0];
			centerIndex[1] = nuclei->GetCentroid()[1];
			m_vectorNucleiCentroids.push_back(centerIndex);

			for (unsigned int j = 1; j < nucleiNo; j++)
				shapeLabelMap2->RemoveLabel(shapeLabelMap2->GetLabelObject(j)->GetLabel());
			for (unsigned int j = nucleiNo + 1; j <= totalNuclei; j++)
				shapeLabelMap2->RemoveLabel(shapeLabelMap2->GetLabelObject(j)->GetLabel());

			LabelMapToBinaryImageFilterType::Pointer labelMapToBinaryImage = LabelMapToBinaryImageFilterType::New();
			labelMapToBinaryImage->SetInput(shapeLabelMap2);
			labelMapToBinaryImage->SetBackgroundValue(Config::getInstance()->getBackgroundPixel());
			labelMapToBinaryImage->SetForegroundValue(Config::getInstance()->getForegroundPixel());
			labelMapToBinaryImage->Update();

			CharImagePointer binaryImage = labelMapToBinaryImage->GetOutput();
			CharImagePointer binaryPatch = PatchExtraction(binaryImage, centerIndex);
			RGBImagePointer  rgbPatch = PatchExtraction(m_ptrInputImage, centerIndex);

			vector< double > objectFeatures;

			objectFeatures.push_back(frameX);
			objectFeatures.push_back(frameY);
			objectFeatures.push_back(nuclei->GetCentroid()[0]);
			objectFeatures.push_back(nuclei->GetCentroid()[1]);

			objectFeatures.push_back(nuclei->GetPhysicalSize());
			objectFeatures.push_back(nuclei->GetRoundness());
			objectFeatures.push_back(nuclei->GetElongation());
			objectFeatures.push_back(nuclei->GetFlatness());
			objectFeatures.push_back(nuclei->GetPerimeter());
			objectFeatures.push_back(nuclei->GetEquivalentSphericalPerimeter());
			objectFeatures.push_back(nuclei->GetEquivalentSphericalRadius());
			objectFeatures.push_back(nuclei->GetEquivalentEllipsoidDiameter()[0]);
			objectFeatures.push_back(nuclei->GetEquivalentEllipsoidDiameter()[1]);

			CharImagePointer redPatch = RedChannelExtraction(rgbPatch);
			CharImagePointer greenPatch = GreenChannelExtraction(rgbPatch);
			CharImagePointer bluePatch = BlueChannelExtraction(rgbPatch);
			CharImagePointer hsvPatch = ColorChannelExtraction(rgbPatch, "HSV", 2);
			CharImagePointer labPatch = ColorChannelExtraction(rgbPatch, "Lab", 0);
			//CharImagePointer luvPatch = ColorChannelExtraction(rgbPatch, "Luv", 0);
			CharImagePointer hePatch = HematoxylinChannelExtraction(rgbPatch, Config::getInstance()->getHEMatrixType());
			CharImagePointer brPatch = BlueRatioExtraction(rgbPatch);

			ComputeFeatures(redPatch, binaryPatch, objectFeatures);		// Red
			ComputeFeatures(greenPatch, binaryPatch, objectFeatures);	// Green
			ComputeFeatures(bluePatch, binaryPatch, objectFeatures);	// Blue
			ComputeFeatures(hsvPatch, binaryPatch, objectFeatures);		// HSV
			ComputeFeatures(labPatch, binaryPatch, objectFeatures);		// Lab
			//ComputeFeatures(luvPatch, binaryPatch, objectFeatures);		// Luv
			ComputeFeatures(hePatch, binaryPatch, objectFeatures);		// H&E
			ComputeFeatures(brPatch, binaryPatch, objectFeatures);		// BR

			for (int up = 2; up <= Config::getInstance()->getMacroPixelUpsamplingValue(); up++)
			{
				redPatch = UpSamplingImage(redPatch, 1);
				greenPatch = UpSamplingImage(greenPatch, 1);
				bluePatch = UpSamplingImage(bluePatch, 1);
				hsvPatch = UpSamplingImage(hsvPatch, 1);
				labPatch = UpSamplingImage(labPatch, 1);
				//luvPatch = UpSamplingImage(luvPatch, 1);
				hePatch = UpSamplingImage(hePatch, 1);
				brPatch = UpSamplingImage(brPatch, 1);

				binaryPatch = UpSamplingImage(binaryPatch, 1);

				ComputeFeatures(redPatch, binaryPatch, objectFeatures);	// Red
				ComputeFeatures(greenPatch, binaryPatch, objectFeatures); // Green
				ComputeFeatures(bluePatch, binaryPatch, objectFeatures); // Blue
				ComputeFeatures(hsvPatch, binaryPatch, objectFeatures);	// HSV
				ComputeFeatures(labPatch, binaryPatch, objectFeatures);	// Lab
				//ComputeFeatures(luvPatch, binaryPatch, objectFeatures);	// Luv
				ComputeFeatures(hePatch, binaryPatch, objectFeatures);	// H&E
				ComputeFeatures(brPatch, binaryPatch, objectFeatures);	// BR
			}

			m_vectorFeatures.push_back(objectFeatures);
		}
	}
}

CharImagePointer NucleiDetector::NucleiRegionDetection_Mosaliganti() {
	CharImagePointer inImage = ExtractChannel();
	inImage->DisconnectPipeline();

	// Cell foreground extraction
	typedef itk::CellForegroundExtraction< CharImageType, ShortImageType, FloatImageType >	CellForegroundExtractionFilterType;
	CellForegroundExtractionFilterType::Pointer cellForegroundFilter = CellForegroundExtractionFilterType::New();
	cellForegroundFilter->SetInput( 0, inImage );
	cellForegroundFilter->SetSigmaForm( 8.0 ); //2.0 ); // in real coordinates
	cellForegroundFilter->SetThresholdCellmin( Config::getInstance()->getThresholdCellMin() );
	cellForegroundFilter->SetThresholdCellmax( Config::getInstance()->getThresholdCellMax() );
	cellForegroundFilter->SetThresholdForm( 10.5 ); //0.5 );
	cellForegroundFilter->SetLargestCellRadius( 10.0 ); // in real coordinates
	float sampling[Dimension] = {2,2};//{5,5,1};
	cellForegroundFilter->SetSampling( sampling );
	cellForegroundFilter->Update();

	typedef itk::RescaleIntensityImageFilter< ShortImageType, CharImageType >	RescaleShortToCharType2;
	RescaleShortToCharType2::Pointer rescaleShortToChar2 = RescaleShortToCharType2::New();
	rescaleShortToChar2->SetInput( cellForegroundFilter->GetOutput() );
	rescaleShortToChar2->SetOutputMinimum( 0 );
	rescaleShortToChar2->SetOutputMaximum( 255 );
	rescaleShortToChar2->Update();
	typedef itk::ImageFileWriter< CharImageType > WriterType;
	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetInput( rescaleShortToChar2->GetOutput() );
	writer2->SetFileName( "foreground.png" );
	try {
		writer2->Update();
	}
	catch ( itk::ExceptionObject e ) {
		cerr << "Error: " << e << endl;
	}

	ShortImagePointer foregroundImage = cellForegroundFilter->GetOutput();
	foregroundImage->DisconnectPipeline();
//	inImage->SetSpacing( foregroundImage->GetSpacing() );
	foregroundImage->SetSpacing( inImage->GetSpacing() );

	FloatImagePointer gaussianCorrelationImage = cellForegroundFilter->GetGaussCorrImage();
	gaussianCorrelationImage->DisconnectPipeline();

	// Cell feature generator
	typedef itk::CellFeatureGenerator< FloatImageType, FloatImageType >	CellFeatureGeneratorFilterType;
	CellFeatureGeneratorFilterType::Pointer cellFeatureGeneratorFilter = CellFeatureGeneratorFilterType::New();

	typedef itk::CastImageFilter< CharImageType, FloatImageType >	CastCharToFloatFilterType;
	CastCharToFloatFilterType::Pointer castCharToFloatFilter = CastCharToFloatFilterType::New();
	castCharToFloatFilter->SetInput( inImage );
	castCharToFloatFilter->Update();
	cellFeatureGeneratorFilter->SetInput( 0, castCharToFloatFilter->GetOutput() );
	cellFeatureGeneratorFilter->SetInput( 1, gaussianCorrelationImage );

	typedef itk::CastImageFilter< ShortImageType, BoolImageType >	CastShortToBoolFilterType;
	CastShortToBoolFilterType::Pointer castShortToBoolFilter = CastShortToBoolFilterType::New();
	castShortToBoolFilter->SetInput( foregroundImage );
	castShortToBoolFilter->Update();
	BoolImageType::Pointer foregroundBoolImage = castShortToBoolFilter->GetOutput();
	cellFeatureGeneratorFilter->SetForeground( foregroundBoolImage );
	cellFeatureGeneratorFilter->SetLargestCellRadius( 4.0 );
	cellFeatureGeneratorFilter->SetSigmaCell( 0.4 );
	cellFeatureGeneratorFilter->SetSigmaCorrelation( 0.4 );
//	cellFeatureGeneratorFilter->SetSampling( sampling );
	cellFeatureGeneratorFilter->SetDistanceMapWeight( 0.7 );
	cellFeatureGeneratorFilter->SetGradientMagnitudeWeight( 0.3 );
	cellFeatureGeneratorFilter->SetGaussCorrWeight( 0.1 );

	typedef itk::RescaleIntensityImageFilter< FloatImageType, FloatImageType >	RescaleType;
	RescaleType::Pointer rescale_f = RescaleType::New();
	rescale_f->SetInput( cellFeatureGeneratorFilter->GetOutput() );
	rescale_f->SetOutputMinimum( 0 );
	rescale_f->SetOutputMaximum( 255 );
	rescale_f->Update();
	FloatImageType::Pointer cellFeatureImage = rescale_f->GetOutput();
	cellFeatureImage->DisconnectPipeline();

	// Seed extraction
	cout << "Seed extraction" << endl << flush;
	typedef itk::SeedExtraction< FloatImageType, BoolImageType > SeedExtractionFilterType;
	SeedExtractionFilterType::Pointer seedExtractionFilter = SeedExtractionFilterType::New();
	seedExtractionFilter->SetForeground( foregroundBoolImage );
	seedExtractionFilter->SetInput( 0, cellFeatureGeneratorFilter->GetDistanceMap() );
	seedExtractionFilter->SetInput( 1, gaussianCorrelationImage );
	seedExtractionFilter->SetLargestCellRadius( 2.0 ); // real coordinates
	seedExtractionFilter->Update();

	// Level set based cell segmentation
	typedef itk::LevelSetBasedCellSegmentation< FloatImageType, ShortImageType >	LevelSetBasedCellSegmentationFilterType;
	LevelSetBasedCellSegmentationFilterType::Pointer levelSetBasedCellSegmentationFilter = LevelSetBasedCellSegmentationFilterType::New();
	levelSetBasedCellSegmentationFilter->SetInput( cellFeatureImage );
	levelSetBasedCellSegmentationFilter->SetLargestCellRadius( 5.0 ); // 5.0 // in real coordinates
	levelSetBasedCellSegmentationFilter->SetSeedValue( 1.5 );
	levelSetBasedCellSegmentationFilter->SetIterations( 500 ); // max number of iterations
	levelSetBasedCellSegmentationFilter->SetPropagationScaling( 2 );
	levelSetBasedCellSegmentationFilter->SetCurvatureScaling( 1.0 ); // default value = 1.0
	levelSetBasedCellSegmentationFilter->SetAdvectionScaling( 4 );
	levelSetBasedCellSegmentationFilter->SetMaxRMSChange( 0.01 ); // used to determine when the level set solution has converged

	typedef FloatImageType::IndexType	FloatImageIndexType;
	map< float, FloatImageIndexType > seeds = seedExtractionFilter->seeds;
	map< float, FloatImageIndexType >::iterator seedsIter;
	FloatImageIndexType idx;
	seedsIter = seeds.end();
	--seedsIter;
	while( seedsIter != seeds.begin() ) {
		--seedsIter;
		float val = seedsIter->first;
	    idx = seedsIter->second;
	    levelSetBasedCellSegmentationFilter->seeds[val] = idx;
	}

	levelSetBasedCellSegmentationFilter->Update();
	ShortImagePointer levelSetImage = levelSetBasedCellSegmentationFilter->GetOutput();

	// Voronoi segmentation
	typedef itk::VoronoiBasedCellSegmentation< ShortImageType, ShortImageType >	VoronoiBasedCellSegmentationFilterType;
	VoronoiBasedCellSegmentationFilterType::Pointer voronoiBasedCellSegmentationFilter = VoronoiBasedCellSegmentationFilterType::New();
	voronoiBasedCellSegmentationFilter->SetInput( levelSetImage );
	voronoiBasedCellSegmentationFilter->SetForegroundImage( foregroundBoolImage );
	voronoiBasedCellSegmentationFilter->SetMinComponentSize( 400 );
	voronoiBasedCellSegmentationFilter->Update();

	// Cell statistics
	typedef itk::CellStatistics< ShortImageType, ShortImageType > CellStatisticsFilterType;
	CellStatisticsFilterType::Pointer cellStatisticsFilter = CellStatisticsFilterType::New();
	cellStatisticsFilter->SetInput( voronoiBasedCellSegmentationFilter->GetOutput() );
	cellStatisticsFilter->SetRawImage( castCharToFloatFilter->GetOutput() );
	cellStatisticsFilter->Update();


	typedef itk::RescaleIntensityImageFilter< ShortImageType, CharImageType >	RescaleShortToCharType;
	RescaleShortToCharType::Pointer rescaleShortToChar = RescaleShortToCharType::New();
	rescaleShortToChar->SetInput( cellStatisticsFilter->GetOutput() );
	rescaleShortToChar->SetOutputMinimum( 0 );
	rescaleShortToChar->SetOutputMaximum( 255 );
	rescaleShortToChar->Update();
	CharImageType::Pointer segmentedImage = rescaleShortToChar->GetOutput();

	typedef itk::ImageFileWriter< CharImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput( segmentedImage );
	writer->SetFileName( "toto.png" );

	try {
		writer->Update();
	}
	catch ( itk::ExceptionObject e ) {
		cerr << "Error: " << e << endl;
	}

	exit(0);
}

CharImagePointer NucleiDetector::NucleiRegionDetection()
{
	CharImagePointer inImage, smoothedImage;
	inImage = ExtractChannel();

	// Perform Grayscale Smoothing
	if( Config::getInstance()->getDetectionChannelName() == "BlueRatio" )
	{
		CharLaplacianRecursiveGaussianImageFilterType::Pointer LoGFilter = CharLaplacianRecursiveGaussianImageFilterType::New();
		LoGFilter->SetInput( inImage );
		LoGFilter->Update();
		smoothedImage = LoGFilter->GetOutput();
	}
	else
	{
		CharImageSizeType medianRadius;
		medianRadius.Fill( 2 );
		CharMedianImageFilterType::Pointer smoothing = CharMedianImageFilterType::New();
		smoothing->SetRadius( medianRadius );
		smoothing->SetInput( inImage );
		smoothing->Update();
		smoothedImage = smoothing->GetOutput();
	}

	// Perform Binary Thresholding
	CharBinaryThresholdImageFilterType::Pointer threshold = CharBinaryThresholdImageFilterType::New();
	threshold->SetLowerThreshold( Config::getInstance()->getLowerThreshold() );
	threshold->SetUpperThreshold( Config::getInstance()->getUpperThreshold() );
	threshold->SetOutsideValue( Config::getInstance()->getBackgroundPixel() );
	threshold->SetInsideValue( Config::getInstance()->getForegroundPixel() );
	threshold->SetInput( smoothedImage );
	//threshold->Update();

	// Hole Filling
	CharVotingBinaryIterativeHoleFillingImageFilterType::Pointer holeFilling = CharVotingBinaryIterativeHoleFillingImageFilterType::New();
	CharStructuringElementType::RadiusType radius;
	radius.Fill( 3 );
	holeFilling->SetRadius( radius );
	holeFilling->SetMajorityThreshold( 1 );
	holeFilling->SetForegroundValue( Config::getInstance()->getForegroundPixel() );
	holeFilling->SetBackgroundValue( Config::getInstance()->getBackgroundPixel() );
	holeFilling->SetInput( threshold->GetOutput() );
	//holeFilling->Update();

	CharStructuringElementType openStructure;
	openStructure.SetRadius( 2 );
	openStructure.CreateStructuringElement();

	// Opening
	CharBinaryOpeningImageFilterType::Pointer opening = CharBinaryOpeningImageFilterType::New();
	opening->SetKernel( openStructure );
	opening->SetForegroundValue( Config::getInstance()->getForegroundPixel() );
	opening->SetBackgroundValue( Config::getInstance()->getBackgroundPixel() );
	opening->SetInput( holeFilling->GetOutput() );
	//opening->Update();

	// Closing
	/*CharBinaryClosingImageFilterType::Pointer closing = CharBinaryClosingImageFilterType::New();
	closing->SetKernel( openStructure );
	closing->SetForegroundValue( Config::getInstance()->getForegroundPixel() );
	closing->SetInput( opening->GetOutput() );
	//closing->Update();
	*/
	// Removing small object
	BinaryShapeOpeningImageFilterType::Pointer smallObjectRemoving = BinaryShapeOpeningImageFilterType::New();
	smallObjectRemoving->SetInput(opening->GetOutput());
	smallObjectRemoving->SetForegroundValue( Config::getInstance()->getForegroundPixel() );
	smallObjectRemoving->SetBackgroundValue( Config::getInstance()->getBackgroundPixel() );
	smallObjectRemoving->SetLambda( Config::getInstance()->getNucleiMinSize() );
	smallObjectRemoving->SetReverseOrdering( false );
	smallObjectRemoving->SetFullyConnected( true );
	smallObjectRemoving->SetAttribute( "PhysicalSize" );
	//smallObjectRemoving->Update();

	// Removing large object
	BinaryShapeOpeningImageFilterType::Pointer largeObjectRemoving = BinaryShapeOpeningImageFilterType::New();
	largeObjectRemoving->SetInput( smallObjectRemoving->GetOutput() );
	largeObjectRemoving->SetForegroundValue( Config::getInstance()->getForegroundPixel() );
	largeObjectRemoving->SetBackgroundValue( Config::getInstance()->getBackgroundPixel() );
	largeObjectRemoving->SetLambda( Config::getInstance()->getNucleiMaxSize() );
	largeObjectRemoving->SetReverseOrdering( true );
	largeObjectRemoving->SetFullyConnected( true );
	largeObjectRemoving->SetAttribute( "PhysicalSize" );
	largeObjectRemoving->Update();

	return largeObjectRemoving->GetOutput();
}

CharImagePointer NucleiDetector::ExtractChannel() {
	switch( Config::getInstance()->getDetectionChannelCode() )	{
	case 0:
		return RedChannelExtraction( m_ptrInputImage );
		break;
	case 1:
		return GreenChannelExtraction( m_ptrInputImage );
		break;
	case 2:
		return BlueChannelExtraction( m_ptrInputImage );
		break;
	case 3:
		return ColorChannelExtraction( m_ptrInputImage, "HSV", 2 );
		break;
	case 4:
		return ColorChannelExtraction( m_ptrInputImage, "Lab", 0 );
		break;
	case 5:
		return HematoxylinChannelExtraction( m_ptrInputImage, Config::getInstance()->getHEMatrixType() );
		break;
	case 6:
		return BlueRatioExtraction(m_ptrInputImage);
		break;
	default:
		return RedChannelExtraction( m_ptrInputImage );
	}
}

void NucleiDetector::SetGTCentroidsFromFile(path filepath) {
	if (!exists(filepath)) {
		std::cerr << "ERROR: Ground truth centroids file " << filepath.make_preferred().string() << " does not exist." << std::endl;
		throw(std::runtime_error("Access error"));
	}

	m_pathGTCentroidsFileName = filepath;
	Read1DIndexesFromCSVFile(m_pathGTCentroidsFileName.make_preferred().string(), m_vectorGTCentroids);
}

void NucleiDetector::ComputeCandidatesLabels(bool insertLabelinFeatureVector) {
	if (m_vectorGTCentroids.size() == 0 || m_vectorNucleiCentroids.size() == 0)
		return;

	// Initialize the TP, FP and FN vector with empty vector
	m_vectorTPCentroids.resize(0);
	m_vectorFPCentroids.resize(0);
	m_vectorFNCentroids.resize(0);

	std::vector< bool > vectorTP;
	vectorTP.resize(m_vectorGTCentroids.size());
	// Initialise vectorTP to false
	for (int i = 0; i < vectorTP.size(); ++i)
		vectorTP[i] = false;

	bool candidateFound = false;
	for (int index_candidate = 0; index_candidate < m_vectorNucleiCentroids.size(); ++index_candidate) {
		candidateFound = false;
		for (int index_GT = 0; index_GT < m_vectorGTCentroids.size(); ++index_GT) {
			double GTcentroidX = m_vectorGTCentroids[index_GT][0];
			double GTcentroidY = m_vectorGTCentroids[index_GT][1];
			double diff_x = (m_vectorNucleiCentroids[index_candidate][0] - GTcentroidX) * Config::getInstance()->getScannerResolutionX();
			diff_x = diff_x * diff_x;
			double diff_y = (m_vectorNucleiCentroids[index_candidate][1] - GTcentroidY) * Config::getInstance()->getScannerResolutionY();
			diff_y = diff_y * diff_y;

			double distance = sqrt(diff_x + diff_y);
			if (distance < Config::getInstance()->getGTRadius()) {
				// Candidate center within range of circle of 8 micrometer from GT mitosis
				candidateFound = true;
				vectorTP[index_GT] = true; // Remember this GT mitosis has been matched with a candidate
				break;
			}
		}

		if (candidateFound) {
			m_vectorTPCentroids.push_back(m_vectorNucleiCentroids[index_candidate]);
			if (insertLabelinFeatureVector)
				m_vectorFeatures[index_candidate].push_back(1); // +1 for Mitosis (0)
		}
		else {
			m_vectorFPCentroids.push_back(m_vectorNucleiCentroids[index_candidate]);
			if (insertLabelinFeatureVector)
				m_vectorFeatures[index_candidate].push_back(-1); // -1 for NonMitosis (1)
		}
	}

	// All the GT mitosis that have not been matched to any candidate are FN
	for (int i = 0; i < vectorTP.size(); ++i)
	if (!vectorTP[i])
		m_vectorFNCentroids.push_back(m_vectorGTCentroids[i]);

}
