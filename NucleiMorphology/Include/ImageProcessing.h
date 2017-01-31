/*
* ImageProcessing.h
*
*  Created on: 20 May 2014
*      Author: Humayun
*/
#ifndef __IMAGEPROCESSING_H__
#define __IMAGEPROCESSING_H__

#include "ITKFunctions.h"
#include "NucleiConfig.h"

CharImagePointer PatchExtraction( CharImagePointer inImage, CharImageIndexType centroid);

RGBImagePointer PatchExtraction( RGBImagePointer inImage, RGBImageIndexType centroid);

RGBImagePointer CheckandChangeInformationImage(RGBImagePointer rgbImage, CharImagePointer binImage);

RGBImagePointer DrawCircle( RGBImagePointer inImage,
                            const std::vector< RGBImageIndexType > &centroids,
                            int colorCode = Config::getInstance()->getOverlayColourCode(),
                            int radius = Config::getInstance()->getOverlayRadius() );

RGBImagePointer DrawSquare( RGBImagePointer inImage,
                            const std::vector< RGBImageIndexType > &centroids,
                            int colorCode = Config::getInstance()->getOverlayColourCode(),
                            int radius = Config::getInstance()->getOverlayRadius() );

void ConvertIndexToVectorData(vector< CharImageIndexType > indexes1D, vector< vector< double > > &data2D);

double ComputePDF(double mean, double variance, double value);

void ComputeImageMeanVariance( CharImageType *inImage, double &mean, double &variance );

void ComputeImageMeanVarianceUsingMask( CharImageType *inImage, CharImageType *maskImage, double &mean, double &variance );

void ComputeCentroids( std::vector<std::vector<CharImageIndexType> > &indexes2D, std::vector<CharImageIndexType> &centroids);

CharImagePointer UpSamplingImage( CharImagePointer inImage, int UpSampleType = 1 );	// UpSampleType,  1 for Mean, 2 for Median, 3 for Minimum, 4 for Maximum

void ComputeFeatures(CharImagePointer inImage, CharImagePointer maskImage, vector<double> &patchFeatures);

CharImagePointer SegmentationPipeLine(CharImagePointer inImage, CharImageIndexType seedPoint);

CharImagePointer SegmentationPipeLine(RGBImagePointer inImage, CharImageIndexType seedPoint);

#endif // __IMAGEPROCESSING_H__
