/*
 * IOFunctions.cpp
 *
 *  Created on: 20 may 2014
 *      Author: Humayun
 */
#include "IOFunctions.h"

RGBImagePointer RGBFileReading( string fileName, string errorMsg )
{
	RGBReaderType::Pointer reader = RGBReaderType::New();
	reader->SetFileName( fileName );
	try {
		reader->Update();
	}
	catch ( itk::ExceptionObject& excp )	{
		cerr << errorMsg << excp << endl;
		throw ( runtime_error("Access error") );
		return NULL;
	}
	return reader->GetOutput();
}

CharImagePointer CharFileReading( string fileName, string errorMsg )
{
	CharReaderType::Pointer reader = CharReaderType::New();
	reader->SetFileName( fileName );
	try {
		reader->Update();
	}
	catch ( itk::ExceptionObject& excp )	{
		cerr << errorMsg << excp << endl;
		throw ( runtime_error("Access error") );
		return NULL;
	}
	return reader->GetOutput();
}

FloatImagePointer FloatFileReading( string fileName, string errorMsg )
{
	FloatReaderType::Pointer reader = FloatReaderType::New();
	reader->SetFileName( fileName );
	try {
		reader->Update();
	}
	catch ( itk::ExceptionObject& excp )	{
		cerr << errorMsg << excp << endl;
		throw ( runtime_error("Access error") );
		return NULL;
	}
	return reader->GetOutput();
}

void RGBFileWriting( RGBImagePointer inImage, string fileName, string errorMsg )
{
	RGBWriterType::Pointer writer = RGBWriterType::New();
	writer->SetFileName( fileName );
	writer->SetInput( inImage );
	try {
		writer->Update();
	}
	catch ( itk::ExceptionObject& excp )	{
		cerr << errorMsg << excp << endl;
		throw ( runtime_error("Access error") );
	}
}

void CharFileWriting( CharImagePointer inImage, string fileName, string errorMsg )
{
	CharWriterType::Pointer writer = CharWriterType::New();
	writer->SetFileName( fileName );
	writer->SetInput( inImage );
	try {
		writer->Update();
	}
	catch ( itk::ExceptionObject& excp )	{
		cerr << errorMsg << excp << endl;
		throw ( runtime_error("Access error") );
	}
}

void FloatFileWriting( FloatImagePointer inImage, string fileName, string errorMsg )
{
	FloatWriterType::Pointer writer = FloatWriterType::New();
	writer->SetFileName( fileName );
	writer->SetInput( inImage );
	try {
		writer->Update();
	}
	catch ( itk::ExceptionObject& excp )	{
		cerr << errorMsg << excp << endl;
		throw ( runtime_error("Access error") );
	}
}

void Read1DIndexesFromCSVFile( string fileName, vector<CharImageIndexType>& data1D ) {
	ifstream inCSVFile( fileName.c_str() );
	string fileLine;
	int mitosisNumber = 0;
	data1D.resize( 0 );
	while( getline( inCSVFile, fileLine ) ) {
		string lineItems;
		CharImageIndexType index;
		istringstream lineStream(fileLine);
		float sumX = 0.0;
		float sumY = 0.0;
		int nbItems = 0;
		float val = 0.0;
		bool xCoord = true; // Read coordinates one by one, first coordinate to be read is X so xCoord = true
                            // second coordinate to be read is Y so xCoord = false
		while( getline( lineStream, lineItems, ',')) {
			val = atof(lineItems.c_str());
			if(xCoord) {
                sumX += val;
                xCoord = false; // Just read X coordinate ==> switch to Y coordinate
            }
			else {
                sumY += val;
                ++nbItems;
                xCoord = true; // Just read Y coordinate ==> switch to X coordinate
			}
		}

        // Check if last coordinate read from file was Y coordinate (xCoord should be true)
        if( !xCoord )
            // Last coordiante read from file was not Y coordinate ==> It was probably confidence degree
            // ==> Remove last value from sumX
            sumX -= val;

		index[(Config::getInstance()->useITKCoordinatesOrder() ? 0 : 1)] = (int)( sumX / (float)nbItems );
		index[(Config::getInstance()->useITKCoordinatesOrder() ? 1 : 0)] = (int)( sumY / (float)nbItems );
		data1D.push_back(index);
		++mitosisNumber;
	}
	inCSVFile.close();
}

void Read2DIndexesFromCSVFile( string fileName, vector<vector<CharImageIndexType> >& data2D )
{
	ifstream inCSVFile( fileName.c_str() );
	string fileLine;
	int mitosisNumber = 0;
	data2D.resize( 0 );
	while( getline( inCSVFile, fileLine ) )
	{
		data2D.resize( mitosisNumber+1 );
		string lineItems;
		int itemNumber=0;
		CharImageIndexType index;
		istringstream lineStream(fileLine);
		while( getline( lineStream, lineItems, ','))
		{
			int val = atoi(lineItems.c_str());
			if(itemNumber%2==0)
				index[0]=val;
			else
			{
				index[1]=val;
				data2D[mitosisNumber].push_back(index);
			}
			itemNumber++;
		}
		mitosisNumber++;
	}
	inCSVFile.close();
}

void Read2DDataFromCSVFile(string fileName, vector< vector< double > >& data2D )
{
	string line;
	ifstream inCSVFile( fileName.c_str() );
	while( getline( inCSVFile, line))
	{
		vector<double> tmp;
		string lineItems;
		istringstream lineStream(line);
		while( getline(lineStream, lineItems, ','))
			tmp.push_back( atof( lineItems.c_str() ) );
		data2D.push_back( tmp );
	}
	inCSVFile.close();
}

void WriteImageToTextFile( CharImageType *inImage, string fileName)
{
	ofstream file( fileName.c_str());
	CharImageRegionType region = inImage->GetLargestPossibleRegion();
	CharImageLinearConstIteratorWithIndexType it (inImage, region);
	for(it.GoToBegin(); !it.IsAtEnd(); it.NextLine())
	{
		while(!it.IsAtEndOfLine())
		{
			file << it.Get() << ",";
			++it;
		}
		file << endl;
	}
	file.close();
}

void WriteSegmentImageIndexToCSVFile( CharImageType *inImage, string fileName)
{
	ofstream file( fileName.c_str());
	CharImageRegionType region = inImage->GetLargestPossibleRegion();
	CharImageLinearConstIteratorWithIndexType it (inImage, region);
	for(it.GoToBegin(); !it.IsAtEnd(); it.NextLine())
	{
		while(!it.IsAtEndOfLine())
		{
			if(it.Get() > 0)
			{
				CharImageIndexType index = it.GetIndex();
				file << index[0] << "," << index[1];
			}
			++it;
		}
	}
	file.close();
}

void Write1DIndexesToCSVFile( string fileName, vector<CharImageIndexType>& indexes1D )
{
	ofstream file(fileName.c_str());
	for( int i=0; i<indexes1D.size(); i++)
	{
		CharImageIndexType index = indexes1D[i];
		file << index[0] << "," << index[1] << endl;
	}
	file.close();
}

void Write2DIndexesToCSVFile( string fileName, vector<vector<CharImageIndexType> >& indexes2D )
{
	ofstream file(fileName.c_str());
	for( int i = 0; i < indexes2D.size(); i++ )
	{
		for( int j = 0; j < indexes2D[i].size(); j++ )
		{
			CharImageIndexType index = indexes2D[i][j];
			file << index[0] << "," << index[1] << ",";
		}
		file << endl;
	}
	file.close();
}

void Write2DDataToCSVFile( string fileName, vector<vector<double> >& data2D )
{
	ofstream file(fileName.c_str());
	for (int i = 0; i < data2D.size(); i++){
		for( int j = 0; j < data2D[i].size(); j++ )
			file << data2D[i][j] << "," ;
		file << endl;
	}
	file.close();
}

void Write3DDataToCSVFile( string fileName, vector<vector<vector<double> > >& data3D )
{
	ofstream file(fileName.c_str());
	for( int i = 0; i < data3D.size(); i++ )
	{
		for( int j = 0; j < data3D[i].size(); j++ )
			for( int k = 0; k < data3D[i][j].size(); k++ )
				file << data3D[i][j][k] << "," ;
		file << endl;
	}
	file.close();
}

void Write2DDataToArffFile( string fileName, vector<vector<double> >& data2D)
{
	ofstream file(fileName.c_str());
	file << "@relation 'MitosisFeatures' \n\n";

	for( int i = 2; i < Config::getInstance()->getAllFeaturesNameSize(); i++)
		file << "@attribute " << Config::getInstance()->getAllFeaturesName(i) << " numeric \n";
	file << "@attribute classLabel {Mitosis,NonMitosis} \n\n@data \n";
	for (int i = 0; i < data2D.size(); i++){
		file << endl;
		int j = 2;
		for (; j < data2D[i].size()-1; j++)
			file << data2D[i][j] << ",";
		// Last column contain label of candidate i.e., +1 = Mitosis, -1 = NonMitosis 
		if (data2D[i][j] == 1)
			file << "Mitosis" << std::endl;
		else
			file << "NonMitosis" << std::endl;
	}
	file.close();
}

void WriteFeatureNamesToCSVFile ( string fileName )
{
	ofstream file( fileName.c_str() );
	for( int i = 0; i < Config::getInstance()->getAllFeaturesNameSize(); i++ )
		file << Config::getInstance()->getAllFeaturesName(i) << "," ;
	file << endl;
	file.close();
}

void WriteFeaturesToCSVFile(string fileName, vector<vector< double > >& data2D)
{
	std::ofstream file(fileName.c_str());

	// First row is name of features
	for (int i = 0; i < Config::getInstance()->getAllFeaturesNameSize(); i++)
		file << Config::getInstance()->getAllFeaturesName(i) << ", ";

	// Write features
	for (int r = 0; r < data2D.size(); r++){
		file << endl;
		for (int c = 0; c < data2D[r].size(); c++)
			file << data2D[r][c] << ",";
	}
	file.close();
}

void WriteFeaturesWithLabelToCSVFile(string fileName, vector<vector< double > >& data2D)
{
	std::ofstream file(fileName.c_str());

	// First row is name of features
	for (int i = 0; i < Config::getInstance()->getAllFeaturesNameSize(); i++)
		file << Config::getInstance()->getAllFeaturesName(i) << ", ";
	file << "ClassLabel" << std::endl;

	// Write features
	for (int r = 0; r < data2D.size(); r++){
		int c = 0;
		for (; c < data2D[r].size()-1; c++)
			file << data2D[r][c] << ",";

		// Last column contain label of candidate i.e., +1 = Mitosis, -1 = NonMitosis 
		if (data2D[r][c] == 1)
			file << "Mitosis" << std::endl;
		else
			file << "NonMitosis" << std::endl;
	}
	file.close();
}
