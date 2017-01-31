/*
 * IOFunctions.h
 *
 *  Created on: 20 may 2014
 *      Author: Humayun
 */
#ifndef __IOFUNCTIONS_H__
#define __IOFUNCTIONS_H__

#include "ITKDeclarations.h"
#include "NucleiConfig.h"

RGBImagePointer RGBFileReading( string fileName, string errorMsg );

CharImagePointer CharFileReading( string fileName, string errorMsg );

FloatImagePointer FloatFileReading( string fileName, string errorMsg );

void RGBFileWriting( RGBImagePointer inImage, string fileName, string errorMsg );

void CharFileWriting( CharImagePointer inImage, string fileName, string errorMsg );

void FloatFileWriting( FloatImagePointer inImage, string fileName, string errorMsg );

void Read1DIndexesFromCSVFile( string fileName, vector<CharImageIndexType>& data1D );

void Read2DIndexesFromCSVFile( string fileName, vector<vector<CharImageIndexType> >& data2D );

void Read2DDataFromCSVFile(string fileName, vector< vector< double > >& data2D );

void WriteImageToTextFile( CharImageType *inImage, string fileName);

void WriteSegmentImageIndexToCSVFile( CharImageType *inImage, string fileName);

void Write1DIndexesToCSVFile( string fileName, vector<CharImageIndexType>& indexes1D );

void Write2DIndexesToCSVFile( string fileName, vector<vector<CharImageIndexType> >& indexes2D );

void Write2DDataToCSVFile( string fileName, vector<vector<double> >& data2D );

void Write3DDataToCSVFile( string fileName, vector<vector<vector<double> > >& data3D );

void Write2DDataToArffFile( string fileName, vector<vector<double> >& data2D);

void WriteFeatureNamesToCSVFile ( string fileName );

void WriteFeaturesToCSVFile(string fileName, vector<vector< double > >& data2D);

void WriteFeaturesWithLabelToCSVFile(string fileName, vector<vector< double > >& data2D);
 
#endif // __IOFUNCTIONS_H__
