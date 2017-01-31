/*
* NucleiDetector.h
*
*  Created on: 20 May 2014
*      Author: Humayun
*/

#ifndef __NUCLEIDETECTOR_H__
#define __NUCLEIDETECTOR_H__

#include "ImageProcessing.h"
#include "IOFunctions.h"
#include "RGBChannelExtraction.h"
#include "ColorChannelExtraction.h"
#include "UtilityFunctions.h"

#include "itkCellFeatureGenerator.h"
#include "itkCellForegroundExtraction.h"
#include "itkSeedExtraction.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkLevelSetBasedCellSegmentation.h"
#include "itkVoronoiBasedCellSegmentation.h"
#include "itkCellStatistics.h"

class NucleiDetector {

public:
	NucleiDetector() { };
	~NucleiDetector() { };

	void NucleiSegmentation(RGBImagePointer rgbImage, string imageName = "", path imagePath = "");
	void NucleiSegmentationFeatureExtraction(RGBImagePointer rgbImage, string imageName = "", path imagePath = "");
	void FeatureExtraction(RGBImagePointer rgbImage, CharImagePointer binImage, string name ="", path imagePath = "");
	void FeatureExtraction_LSM(RGBImagePointer rgbImage, CharImagePointer binImage, string name = "", path imagePath = "");
	void FeatureExtraction_TCGA(RGBImagePointer rgbImage, CharImagePointer binImage, string name = "", path imagePath = "");

	void ComputeCandidatesLabels(bool insertLabelinFeatureVector = true);

	void SetGTCentroidsFromFile(path filepath);
	void SetCandidatesCentroids(vector< RGBImageIndexType > inCentroids)	{ m_vectorNucleiCentroids = inCentroids; };

	vector< CharImageIndexType >	GetGTCentroids()		{ return m_vectorGTCentroids; };
	vector< CharImageIndexType >	GetTPCentroids()		{ return m_vectorTPCentroids; };
	vector< CharImageIndexType >	GetFNCentroids()		{ return m_vectorFNCentroids; };
	vector< CharImageIndexType >	GetFPCentroids()		{ return m_vectorFPCentroids; };
	vector< RGBImageIndexType >		GetNucleiCentroids()	{ return m_vectorNucleiCentroids; };
	vector< vector< double > >		GetFeatures()			{ return m_vectorFeatures; };
	CharImagePointer				GetResultImage()		{ return m_ptrOutputImage; };

	int		GetNumberOfGT()	{ return m_vectorGTCentroids.size(); };
	int		GetNumberOfTP()	{ return m_vectorTPCentroids.size(); };
	int		GetNumberOfFN()	{ return m_vectorFNCentroids.size(); };
	int		GetNumberOfFP()	{ return m_vectorFPCentroids.size(); };

private:

	path								m_pathGTCentroidsFileName;
	path								m_pathInputImageName;
	string								m_strImageName;
	vector< RGBImageIndexType >			m_vectorNucleiCentroids;
	vector< CharImageIndexType >		m_vectorGTCentroids;
	vector< CharImageIndexType >		m_vectorTPCentroids;
	vector< CharImageIndexType >		m_vectorFNCentroids;
	vector< CharImageIndexType >		m_vectorFPCentroids;
	vector< vector< double > >			m_vectorFeatures;

	RGBImagePointer		m_ptrInputImage;
	CharImagePointer	m_ptrOutputImage;

	CharImagePointer NucleiRegionDetection();
	CharImagePointer NucleiRegionDetection_Mosaliganti();
	CharImagePointer ExtractChannel();
};

#endif /* __NUCLEIDETECTOR_H__ */
