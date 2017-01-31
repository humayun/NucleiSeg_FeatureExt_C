/*
 * Main.cpp
 *
 *  Created on: 20 May 2014
 *      Author: Humayun
 */
//#include "QuickView.h"

#include "NucleiFeatures.h"
#include "UtilityFunctions.h"
#include "NucleiDetector.h"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

void NucleiSegmentation(void);
void NucleiSegmentationAndFeatureComputation(void);
void FeatureComputation(void);
void FeaturesComputation_TCGA(void);
void FeatureComputation_LSM(void);
void FeaturesComputation_TMA(void);
void MitosisSegFeatureComputation(void);
void MitosisFeatureComputation(void);
void Test(void);

int main(int argc, char * argv[])
{
	cout << "\n\n****************** Program start ********************* \n\n";

	/*char *gv[12];
	gv[0] = argv[1];
	gv[1] = "-i";
	gv[2] = "H:/WS/MITOS_A_TestingSet";
	gv[3] = "-f";
	gv[4] = "6";
	gv[5] = "-l";
	gv[6] = "128";
	gv[7] = "-c";
	gv[8] = "2";
	gv[9] = "-r";
	gv[10] = "2";
	//gv[11] = "-p";
	Config::getInstance()->parseCommandLine(11, gv);
	Test();
	return 1;*/

	if (argc > 1)
		Config::getInstance()->parseCommandLine(argc, argv);
	else
		Config::getInstance()->setDefaultValue();
	
	switch (Config::getInstance()->computeFeatures()){
	case 0:
		NucleiSegmentation();
		break;
	case 1:
		NucleiSegmentationAndFeatureComputation();
		break;
	case 2:
		FeatureComputation();
		break;
	case 3:
		FeatureComputation_LSM();
		break;
	case 4:
		FeaturesComputation_TCGA();
		break;
	case 5:
		MitosisSegFeatureComputation();
		break;
	case 6: 
		MitosisFeatureComputation();
		break;
	default:
		FeatureComputation();
	}
	
	cout << "\n\n****************** Program finish ******************** \n\n ";
	return 0;
}

// 0 - It takes an input folder, segment nuclei for all images in that folder
void NucleiSegmentation()
{
	cout << "\n\n\t******** Nuclei Segmentation ********" << endl;

	path directory(Config::getInstance()->getInputDirectoryOrFile());
	string directoryName = directory.make_preferred().string().substr(directory.make_preferred().string().find_last_of("/\\") + 1);

	if (!exists(directory)) {
		cerr << "ERROR: Directory " << directory.make_preferred().string() << " does not exist!" << endl;
		throw (runtime_error("Access error"));
	}
	if (!is_directory(directory)) {
		cout << "ERROR: " << directory.make_preferred().string() << " is not a directory!" << endl;
		throw (runtime_error("Access error"));
	}
	cout << "\nDirectory: " << directory.make_preferred().string() << " exists. Start reading data ... " << endl;

	// Create new Folder for storing segmentation files
	path dirSeg(directory / Config::getInstance()->getSegmentationFolder() );
	if (!exists(dirSeg.make_preferred().string())){
		if (create_directory(dirSeg))
			cout << "\nCreate new directory " << dirSeg.make_preferred().string() << " for storing segmentation files ... \n";
	}
	else
		cout << "\nDirectory " << dirSeg .make_preferred().string() << " already exist for storing segmentation files ... \n";

	// Reading images in given directory
	vector<path> pathVector;
	copy(directory_iterator(directory), directory_iterator(), back_inserter(pathVector));
	sort(pathVector.begin(), pathVector.end());
	int ct = 0;
	BOOST_FOREACH(path p, pathVector) {
		path filename = p.filename();
		std::string filenameWithoutExtension = filename.string().substr(0, filename.string().find_last_of("."));
		std::string fileExtension = filename.string().substr(filename.string().find_last_of(".") + 1);
		boost::algorithm::to_lower(fileExtension);

		if (string("bmp") == fileExtension || string("tiff") == fileExtension || string("tif") == fileExtension || string("jpeg") == fileExtension || string("png") == fileExtension){

			RGBImagePointer rgbImage = RGBFileReading(p.make_preferred().string(), string("\n Err 1: Frame File Reading Err: Exception caught\n"));
			cout << "\n" << ++ct << " - Nuclei Segmentation in image: " << p.make_preferred().string() << endl;

			// Calling Nuclei Segmentation
			NucleiDetector cd;
			cd.NucleiSegmentation(rgbImage, filenameWithoutExtension, directory);

			// Write Nuclei Centroids
			vector< CharImageIndexType > nucleiCentroid = cd.GetNucleiCentroids();
			path centroidsFile(dirSeg);
			centroidsFile /= filenameWithoutExtension + Config::getInstance()->getCentroidsFile() + Config::getInstance()->getCsvExtension();
			Write1DIndexesToCSVFile(centroidsFile.make_preferred().string(), nucleiCentroid);
			
			// Write binary segmentation image
			CharImagePointer binImage = cd.GetResultImage();
			path detectionFile(dirSeg);
			detectionFile /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtension();
			CharFileWriting(binImage, detectionFile.make_preferred().string(), string("\nError: Binary Detection image writing ..."));

			// write nuclei overlay image
			RGBImagePointer overlayImage = DrawCircle(rgbImage, nucleiCentroid, Config::getInstance()->getOverlayColourCode());
			path overlayFile(dirSeg);
			overlayFile /= filenameWithoutExtension + Config::getInstance()->getNucleiOverlayFile() + Config::getInstance()->getImageExtensionn();
			RGBFileWriting(overlayImage, overlayFile.make_preferred().string(), string("\nError: RGB Detection image writing ..."));
		}
	}
}

// 1 - It takes an input folder, segment nuclei and compute features for all images in that folder
void NucleiSegmentationAndFeatureComputation(void)
{
	cout << "\n\n\t******** Nuclei Segmentation and Feature Computation in Folder ********" << endl;

	path directory(Config::getInstance()->getInputDirectoryOrFile());
	string directoryName = directory.make_preferred().string().substr(directory.make_preferred().string().find_last_of("/\\") + 1);

	if (!exists(directory)) {
		cerr << "ERROR: Directory " << directory.make_preferred().string() << " does not exist!" << endl;
		throw (runtime_error("Access error"));
	}
	if (!is_directory(directory)) {
		cout << "ERROR: " << directory.make_preferred().string() << " is not a directory!" << endl;
		throw (runtime_error("Access error"));
	}
	cout << directory.make_preferred().string() << " directory exist. Start reading images ... " << endl;

	// Create new Folder for storing segmentation files
	path dirSeg( directory / Config::getInstance()->getSegmentationFolder() );
	if (!exists(dirSeg.make_preferred().string())){
		if (create_directory(dirSeg))
			cout << "\nCreate new directory: " << dirSeg.make_preferred().string() << " for storing segmentation files ... \n";
	}
	else
		cout << "\nDirectory: " << dirSeg.make_preferred().string() << " already exist for storing segmentation files ... \n";

	// Create new Folder for storing Features files
	string featureFolderName = directoryName;
	if (Config::getInstance()->regionFeatures())
		featureFolderName = featureFolderName + string("_C_Region_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	else
		featureFolderName = featureFolderName + string("_C_Patch_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	path dirFeatures(directory / featureFolderName);
	if (!exists(dirFeatures.make_preferred().string())){
		if (create_directory(dirFeatures))
			cout << "\nCreate new directory: " << dirFeatures.make_preferred().string() << " for storing features files ... \n";
	}
	else
		cout << "\nDirectory: " << dirFeatures.make_preferred().string() << " already exists for storing features files ... \n";

	vector< vector< double > > allFeatures;

	// Reading images in given directory
	vector<path> pathVector;
	copy(directory_iterator(directory), directory_iterator(), back_inserter(pathVector));
	sort(pathVector.begin(), pathVector.end());
	int ct = 0;
	BOOST_FOREACH(path p, pathVector) {
		path filename = p.filename();
		std::string filenameWithoutExtension = filename.string().substr(0, filename.string().find_last_of("."));
		std::string fileExtension = filename.string().substr(filename.string().find_last_of(".") + 1);
		boost::algorithm::to_lower(fileExtension);

		if (string("bmp") == fileExtension || string("tiff") == fileExtension || string("tif") == fileExtension || string("jpeg") == fileExtension || string("png") == fileExtension){

			RGBImagePointer rgbImage = RGBFileReading(p.make_preferred().string(), string("\n Err 1: RGB TMA File Reading Err: Exception caught\n"));

			cout << "\n" << ++ct << " - Nuclei Segmentation and Feature Computation in image: " << p.make_preferred().string() << endl;

			// Nuclei Segmentation and Feature Computation
			NucleiDetector cd;
			cd.NucleiSegmentationFeatureExtraction(rgbImage, filenameWithoutExtension, directory);

			// Write Binary Image
			CharImagePointer binImage = cd.GetResultImage();
			path binFileName (dirSeg);
			binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtension();
			CharFileWriting(binImage, binFileName.make_preferred().string(), "\nError: Binary File Writing ...\n");

			// Write Nuclei Features 
			vector< vector< double > > nucleiFeatures = cd.GetFeatures();
			path featuresFile(dirFeatures);
			featuresFile /= filenameWithoutExtension + Config::getInstance()->getFeaturesFile() + Config::getInstance()->getCsvExtension();
			Write2DDataToCSVFile(featuresFile.make_preferred().string(), nucleiFeatures);

			// Write Nuclei Centroids
			vector< CharImageIndexType > nucleiCentroid = cd.GetNucleiCentroids();
			path centroidFile(dirFeatures);
			centroidFile /= filenameWithoutExtension + Config::getInstance()->getCentroidsFile() + Config::getInstance()->getCsvExtension();
			Write1DIndexesToCSVFile(centroidFile.make_preferred().string(), nucleiCentroid);
			
			if (Config::getInstance()->writeDetectionImages()){
				// Write binary segmentation image
				CharImagePointer binImage = cd.GetResultImage();
				path detectionFile(dirSeg);
				detectionFile /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtension();
				CharFileWriting(binImage, detectionFile.make_preferred().string(), string("\nError: Binary Detection image writing ..."));

				// write nuclei overlay image
				RGBImagePointer overlayImage = DrawCircle(rgbImage, nucleiCentroid, Config::getInstance()->getOverlayColourCode());
				path overlayFile(dirSeg);
				overlayFile /= filenameWithoutExtension + Config::getInstance()->getNucleiOverlayFile() + Config::getInstance()->getImageExtensionn();
				RGBFileWriting(overlayImage, overlayFile.make_preferred().string(), string("\nError: RGB Detection image writing ..."));
			}

			for (int instanceNo = 0; instanceNo < nucleiFeatures.size(); instanceNo++)
				allFeatures.push_back(nucleiFeatures[instanceNo]);
		}
	}

	path allFeaturesFile(dirFeatures);
	allFeaturesFile /= featureFolderName + Config::getInstance()->getCsvExtension();
	WriteFeaturesWithLabelToCSVFile(allFeaturesFile.make_preferred().string(), allFeatures);
}

// 2 - It takes an input folder (containing Segmentation folder inside) and compute features for all images in that folder
void FeatureComputation(void)
{
	cout << "\n\n\t******** Feature Computation in Folder ********" << endl;

	path directory(Config::getInstance()->getInputDirectoryOrFile());
	string directoryName = directory.make_preferred().string().substr(directory.make_preferred().string().find_last_of("/\\") + 1);

	if (!exists(directory)) {
		cerr << "ERROR: Directory: " << directory.make_preferred().string() << " does not exist!" << endl;
		throw (runtime_error("Access error"));
	}
	if (!is_directory(directory)) {
		cout << "ERROR: " << directory.make_preferred().string() << " is not a directory!" << endl;
		throw (runtime_error("Access error"));
	}
	cout << "\nDirectory: " << directory.make_preferred().string() << " exists for reading images ... \n";

	// Check Segmentation Folder for binary/segmentation files
	path dirSeg( directory / Config::getInstance()->getSegmentationFolder() );
	if (!exists(dirSeg.make_preferred().string())){
		cerr << "\nERROR: Directory: " << dirSeg.make_preferred().string() << " does not exist for reading segmentation files! \n";
		throw (runtime_error("Access error"));
	}
	else
		cout << "\nDirectory: " << dirSeg.make_preferred().string() << " already exists for reading segmentation files ... \n";

	// Create new Folder for storing Features files
	string featureFolderName = directoryName;
	if (Config::getInstance()->regionFeatures())
		featureFolderName = featureFolderName + string("_C_Region_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	else
		featureFolderName = featureFolderName + string("_C_Patch_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	path dirFeatures(directory / featureFolderName);
	if (!exists(dirFeatures.make_preferred().string())){
		if (create_directory(dirFeatures))
			cout << "\nCreate new directory: " << dirFeatures.make_preferred().string() << " for storing features files ... \n";
	}
	else
		cout << "\nDirectory: " << dirFeatures.make_preferred().string() << " already exists for storing features files ... \n";

	vector< vector< double > > allFeatures;

	// Reading images in given directory
	vector<path> pathVector;
	copy(directory_iterator(directory), directory_iterator(), back_inserter(pathVector));
	sort(pathVector.begin(), pathVector.end());
	int ct = 0;
	BOOST_FOREACH(path p, pathVector) {
		path filename = p.filename();
		std::string filenameWithoutExtension = filename.string().substr(0, filename.string().find_last_of("."));
		std::string fileExtension = filename.string().substr(filename.string().find_last_of(".") + 1);
		boost::algorithm::to_lower(fileExtension);

		if (string("bmp") == fileExtension || string("tiff") == fileExtension || string("tif") == fileExtension || string("jpeg") == fileExtension || string("png") == fileExtension){

			// Check if feature is already computed for this image, then skin computing again
			path featuresFile(dirFeatures);
			featuresFile /= filenameWithoutExtension + Config::getInstance()->getFeaturesFile() + Config::getInstance()->getCsvExtension();
			if (exists(featuresFile.make_preferred().string()))
				continue;

			path binFileName = dirSeg;
			binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtension();
			if (!exists(binFileName.make_preferred().string())){
				binFileName = dirSeg;
				binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtensionn();
				if (!exists(binFileName.make_preferred().string())){
					binFileName = dirSeg;
					binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + std::string(".png");
					if (!exists(binFileName.make_preferred().string())){
						binFileName = dirSeg;
						binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + std::string(".jpeg");
						if (!exists(binFileName.make_preferred().string())){
							binFileName = dirSeg;
							binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + std::string(".bmp");
						}
					}
				}
			}

			RGBImagePointer rgbImage = RGBFileReading(p.make_preferred().string(), string("\n Err 1: RGB frame reading Err: Exception caught\n"));
			CharImagePointer binImage = CharFileReading(binFileName.make_preferred().string(), string("\n Err 2: Binary frame reading Err: Exception caught\n"));

			cout << "\n" << ++ct <<" - Nuclei Feature Computation in image: " << p.make_preferred().string() << endl;

			// Nuclei Feature Computation
			NucleiDetector cd;
			cd.FeatureExtraction(rgbImage, binImage, filenameWithoutExtension, directory);

			// Write Nuclei Features 
			vector< vector< double > >		nucleiFeatures = cd.GetFeatures();
			Write2DDataToCSVFile(featuresFile.make_preferred().string(), nucleiFeatures);

			// Write Nuclei Centroids
			if (Config::getInstance()->writeDetectionImages()){
				vector< CharImageIndexType > nucleiCentroid = cd.GetNucleiCentroids();
				path centroidFile(dirFeatures);
				centroidFile /= filenameWithoutExtension + Config::getInstance()->getCentroidsFile() + Config::getInstance()->getCsvExtension();
				Write1DIndexesToCSVFile(centroidFile.make_preferred().string(), nucleiCentroid);
			}

			for (int instanceNo = 0; instanceNo < nucleiFeatures.size(); instanceNo++)
				allFeatures.push_back(nucleiFeatures[instanceNo]);
		}
	}

	path allFeaturesFile(dirFeatures);
	allFeaturesFile /= featureFolderName + Config::getInstance()->getCsvExtension();
	WriteFeaturesWithLabelToCSVFile(allFeaturesFile.make_preferred().string(), allFeatures);
}

// 3 - It takes an input folder (containing Segmentation folder inside) and compute features for only Green channel
void FeatureComputation_LSM(void)
{
	cout << "\n\n\t******** Feature Computation (LSM) ********" << endl;

	path directory(Config::getInstance()->getInputDirectoryOrFile());
	string directoryName = directory.make_preferred().string().substr(directory.make_preferred().string().find_last_of("/\\") + 1);

	if (!exists(directory)) {
		cerr << "ERROR: Directory " << directory.make_preferred().string() << " does not exist!" << endl;
		throw (runtime_error("Access error"));
	}
	if (!is_directory(directory)) {
		cout << "ERROR: " << directory.make_preferred().string() << " is not a directory!" << endl;
		throw (runtime_error("Access error"));
	}
	cout << "\nDirectory: " << directory.make_preferred().string() << " exists for reading images ... \n";

	// Check Segmentation Folder for binary/segmentation files
	path dirSeg( directory / Config::getInstance()->getSegmentationFolder() );
	if (!exists(dirSeg.make_preferred().string())){
		cerr << "\nERROR: Directory: " << dirSeg.make_preferred().string() << " does not exist for reading segmentation files! \n";
		throw (runtime_error("Access error"));
	}
	else
		cout << "\nDirectory: " << dirSeg.make_preferred().string() << " exists for reading segmentation files ... \n";

	// Create new Folder for storing Features files
	string featureFolderName = directoryName;
	if (Config::getInstance()->regionFeatures())
		featureFolderName = featureFolderName + string("_C_Region_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	else
		featureFolderName = featureFolderName + string("_C_Patch_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	path dirFeatures(directory / featureFolderName);
	if (!exists(dirFeatures.make_preferred().string())){
		if (create_directory(dirFeatures))
			cout << "\nCreate new directory: " << dirFeatures.make_preferred().string() << " for storing features files ... \n";
	}
	else
		cout << "\nDirectory: " << dirFeatures.make_preferred().string() << " already exists for storing features files ... \n";

	vector< vector< double > > allFeatures;

	// Reading images in given directory
	vector<path> pathVector;
	copy(directory_iterator(directory), directory_iterator(), back_inserter(pathVector));
	sort(pathVector.begin(), pathVector.end());
	int ct = 0;
	BOOST_FOREACH(path p, pathVector) {
		path filename = p.filename();
		std::string filenameWithoutExtension = filename.string().substr(0, filename.string().find_last_of("."));
		std::string fileExtension = filename.string().substr(filename.string().find_last_of(".") + 1);
		boost::algorithm::to_lower(fileExtension);

		if (string("bmp") == fileExtension || string("tiff") == fileExtension || string("tif") == fileExtension || string("png") == fileExtension || string("jpeg") == fileExtension){

			path binFileName = dirSeg;
			binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtension();
			if (!exists(binFileName.make_preferred().string())){
				binFileName = dirSeg;
				binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtensionn();
			}

			RGBImagePointer rgbImage = RGBFileReading(p.make_preferred().string(), string("\n Err 1: RGB frame reading Err: Exception caught\n"));
			CharImagePointer binImage = CharFileReading(binFileName.make_preferred().string(), string("\n Err 2: Binary frame reading Err: Exception caught\n"));

			cout << "\n" << ++ct << " - Nuclei Feature Computation in image: " << p.make_preferred().string() << endl;

			// Nuclei Feature Computation
			NucleiDetector cd;
			cd.FeatureExtraction_LSM(rgbImage, binImage, filenameWithoutExtension, directory);

			// Write Nuclei Features 
			vector< vector< double > >		nucleiFeatures = cd.GetFeatures();
			path featuresFile(dirFeatures);
			featuresFile /= filenameWithoutExtension + Config::getInstance()->getFeaturesFile() + Config::getInstance()->getCsvExtension();
			Write2DDataToCSVFile(featuresFile.make_preferred().string(), nucleiFeatures);

			// Write Nuclei Centroids
			if (Config::getInstance()->writeDetectionImages()){
				vector< CharImageIndexType > nucleiCentroid = cd.GetNucleiCentroids();
				path centroidFile(dirFeatures);
				centroidFile /= filenameWithoutExtension + Config::getInstance()->getCentroidsFile() + Config::getInstance()->getCsvExtension();
				Write1DIndexesToCSVFile(centroidFile.make_preferred().string(), nucleiCentroid);
			}
			for (int instanceNo = 0; instanceNo < nucleiFeatures.size(); instanceNo++)
				allFeatures.push_back(nucleiFeatures[instanceNo]);

		}
	}
	path featureNamesFile(dirFeatures);
	featureNamesFile /= string("FeaturesName") + Config::getInstance()->getCsvExtension();
	WriteFeatureNamesToCSVFile(featureNamesFile.make_preferred().string());

	path allFeaturesFile(dirFeatures);
	allFeaturesFile /= featureFolderName + Config::getInstance()->getCsvExtension();
	WriteFeaturesWithLabelToCSVFile(allFeaturesFile.make_preferred().string(), allFeatures);
}

// 4 - It takes an input folder and read output.txt file in given folder (having a list of all WSI folder in the directory) and compute feature for each WSI folder
void FeaturesComputation_TCGA(void)
{
	cout << "\n\n\t******** Feature Computation (TCGA) ********" << endl;

	if (!exists(Config::getInstance()->getInputDirectoryOrFile())) {
		cerr << "ERROR: Directory/File " << Config::getInstance()->getInputDirectoryOrFile().make_preferred().string() << " does not exist." << endl;
		throw (runtime_error("Access error"));
	}

	cout << endl << Config::getInstance()->getInputDirectoryOrFile().make_preferred().string() << " directory/file exist. Start reading file ..." << endl;

	path TCGAWSIFiles(Config::getInstance()->getInputDirectoryOrFile());
	path mainTCGADir = TCGAWSIFiles.parent_path();

	ifstream TCGAOutputFile(TCGAWSIFiles.make_preferred().string().c_str());
	string WSIName;
	while (getline(TCGAOutputFile, WSIName)){
		path WSIImages(mainTCGADir);
		WSIImages /= WSIName;

		if (!exists(WSIImages.make_preferred().string())) {
			cerr << "ERROR: Directory " << WSIImages.make_preferred().string() << " does not exist." << endl;
			throw (runtime_error("Access error"));
		}
		if (!is_directory(WSIImages.make_preferred().string())) {
			cout << "ERROR: " << WSIImages.make_preferred().string() << " is not a directory!" << endl;
			throw (runtime_error("Access error"));
		}
		cout << "\nDirectory: " << WSIImages.make_preferred().string() << " exists for reading images ... \n";

		// Check Segmentation Folder for binary/segmentation files
		path dirSeg(WSIImages);
		dirSeg /= Config::getInstance()->getSegmentationFolder();
		if (!exists(dirSeg.make_preferred().string())){
			cerr << "\nERROR: Directory: " << dirSeg.make_preferred().string() << " does not exist for reading segmentation files! \n";
			throw (runtime_error("Access error"));
		}
		else
			cout << "\nDirectory: " << dirSeg.make_preferred().string() << " exists for reading segmentation files ... \n";

		// Check Features Folder for storing Features files
		string featureFolderName = WSIName;
		if (Config::getInstance()->regionFeatures())
			featureFolderName = featureFolderName + string("_C_Region_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
			+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
			+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
		else
			featureFolderName = featureFolderName + string("_C_Patch_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
			+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
			+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
		path dirFeatures( WSIImages / featureFolderName);
		if (!exists(dirFeatures.make_preferred().string())){
			if (create_directory(dirFeatures))
				cout << "\nCreate new directory: " << dirFeatures.make_preferred().string() << " for storing features files ... \n";
		}
		else
			cout << "\nDirectory: " << dirFeatures.make_preferred().string() << " already exists for storing features files ... \n";

		vector< vector< double > > allFeatures;
		vector<path> pathVector;
		copy(directory_iterator(WSIImages.make_preferred().string()), directory_iterator(), back_inserter(pathVector));
		sort(pathVector.begin(), pathVector.end());
		int ct = 0;
		BOOST_FOREACH(path p, pathVector) {
			path filename = p.filename();
			std::string filenameWithoutExtension = filename.string().substr(0, filename.string().find_last_of("."));
			std::string fileExtension = filename.string().substr(filename.string().find_last_of(".") + 1);
			boost::algorithm::to_lower(fileExtension);

			if (string("bmp") == fileExtension || string("tiff") == fileExtension || string("tif") == fileExtension || string("jpeg") == fileExtension || string("png") == fileExtension){

				path binFileName (dirSeg);
				binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtension();
				if (!exists(binFileName.make_preferred().string())){
					binFileName = dirSeg;
					binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtensionn();
				}

				RGBImagePointer rgbImage = RGBFileReading(p.make_preferred().string(), string("\n Err 1: RGB frame reading Err: Exception caught\n"));
				CharImagePointer binImage = CharFileReading(binFileName.make_preferred().string(), string("\n Err 2: Binary frame feading Err: Exception caught\n"));

				cout << "\n" << ++ct << " - Nuclei Feature Computation in image: " << p.make_preferred().string() << endl;

				// Nuclei Detection and Feature Computation
				NucleiDetector cd;
				cd.FeatureExtraction_TCGA(rgbImage, binImage, filenameWithoutExtension, WSIImages);

				// Write Nuclei Features 
				vector< vector< double > >		nucleiFeatures = cd.GetFeatures();
				path featuresFile(dirFeatures);
				featuresFile /= filenameWithoutExtension + Config::getInstance()->getFeaturesFile() + Config::getInstance()->getCsvExtension();
				Write2DDataToCSVFile(featuresFile.make_preferred().string(), nucleiFeatures);

				// Write Nuclei Centroids
				if (Config::getInstance()->writeDetectionImages()){
					vector< CharImageIndexType > nucleiCentroid = cd.GetNucleiCentroids();
					path centroidFile(dirFeatures);
					centroidFile /= filenameWithoutExtension + Config::getInstance()->getCentroidsFile() + Config::getInstance()->getCsvExtension();
					Write1DIndexesToCSVFile(centroidFile.make_preferred().string(), nucleiCentroid);
				}

				for (int instanceNo = 0; instanceNo < nucleiFeatures.size(); instanceNo++)
					allFeatures.push_back(nucleiFeatures[instanceNo]);
			}
		}

		path allFeaturesFile(mainTCGADir);
		allFeaturesFile /= featureFolderName + Config::getInstance()->getCsvExtension();
		WriteFeaturesToCSVFile(allFeaturesFile.make_preferred().string(), allFeatures);
	}
}

// 5 - It takes an input folder (containing images, ground truth for mitosis and segmentation folder inside) and compute features for each image in folder
void MitosisSegFeatureComputation(void)
{
	cout << "\n\n\t******** Mitosis Candidate Segmentation and Feature Computation ********" << endl;

	path directory(Config::getInstance()->getInputDirectoryOrFile());
	string directoryName = directory.make_preferred().string().substr(directory.make_preferred().string().find_last_of("/\\") + 1);

	if (!exists(directory)) {
		cerr << "ERROR: Directory " << directory.make_preferred().string() << " does not exist!" << endl;
		throw (runtime_error("Access error"));
	}
	if (!is_directory(directory)) {
		cout << "ERROR: " << directory.make_preferred().string() << " is not a directory!" << endl;
		throw (runtime_error("Access error"));
	}
	cout << directory.make_preferred().string() << " directory exist. Start reading images ... " << endl;

	// Create new Folder for storing segmentation files
	path dirSeg(directory);
	dirSeg /= Config::getInstance()->getSegmentationFolder();
	if (!exists(dirSeg.make_preferred().string())){
		if (create_directory(dirSeg))
			cout << "\nCreate new directory: " << dirSeg.make_preferred().string() << " for storing segmentation files ... \n";
	}
	else
		cout << "\nDirectory: " << dirSeg.make_preferred().string() << " already exist for storing segmentation files ... \n";

	// Create new Folder for storing Features files
	string featureFolderName = directoryName;
	if (Config::getInstance()->regionFeatures())
		featureFolderName = featureFolderName + string("_C_Region_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	else
		featureFolderName = featureFolderName + string("_C_Patch_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	path dirFeatures(directory / featureFolderName);
	if (!exists(dirFeatures.make_preferred().string())){
		if (create_directory(dirFeatures))
			cout << "\nCreate new directory: " << dirFeatures.make_preferred().string() << " for storing features files ... \n";
	}
	else
		cout << "\nDirectory: " << dirFeatures.make_preferred().string() << " already exists for storing features files ... \n";

	unsigned int TP = 0, FP = 0, GT = 0, FN = 0;
	vector< vector< double > > allFeatures;

	std::stringstream fileResult;
	fileResult << "CandidateDetection_Results" << Config::getInstance()->getCsvExtension();
	path CandDetResultFile(directory);
	CandDetResultFile /= fileResult.str();
	std::ofstream file(CandDetResultFile.make_preferred().string().c_str());

	// Reading images in given directory
	vector<path> pathVector;
	copy(directory_iterator(directory), directory_iterator(), back_inserter(pathVector));
	sort(pathVector.begin(), pathVector.end());
	int ct = 0;
	BOOST_FOREACH(path p, pathVector) {
		path filename = p.filename();
		std::string filenameWithoutExtension = filename.string().substr(0, filename.string().find_last_of("."));
		std::string fileExtension = filename.string().substr(filename.string().find_last_of(".") + 1);
		boost::algorithm::to_lower(fileExtension);

		if (string("bmp") != fileExtension && string("tiff") != fileExtension && string("tif") != fileExtension && string("jpeg") != fileExtension && string("png") != fileExtension)
			continue;

		RGBImagePointer rgbImage = RGBFileReading(p.make_preferred().string(), string("\n Err 1: RGB Image Reading Err: Exception caught\n"));
		cout << "\n" << ++ct << " - Mitosis candidate detection in image: " << p.make_preferred().string() << endl;

		// Compare candidates with GT mitosis and find out TP (Mitosis) and FP (NonMitosis)
		std::stringstream GTCentroidsFileName;
		GTCentroidsFileName << filenameWithoutExtension << Config::getInstance()->getCentroidsFile() << Config::getInstance()->getCsvExtension();
		path GTCentroidsFile(directory / GTCentroidsFileName.str());

		// Nuclei Segmentation and Feature Computation
		NucleiDetector cd;
		cd.NucleiSegmentationFeatureExtraction(rgbImage, filenameWithoutExtension, directory);
		cd.SetGTCentroidsFromFile(GTCentroidsFile);
		cd.ComputeCandidatesLabels(true);
		GT = GT + cd.GetNumberOfGT();
		TP = TP + cd.GetNumberOfTP();
		FN = FN + cd.GetNumberOfFN();
		FP = FP + cd.GetNumberOfFP();

		// Write Binary Image
		CharImagePointer binImage = cd.GetResultImage();
		path binFileName(dirSeg);
		binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtension();
		CharFileWriting(binImage, binFileName.make_preferred().string(), "\nError: Binary File Writing ...\n");

		// Write Nuclei Features 
		vector< vector< double > > nucleiFeatures = cd.GetFeatures();
		path featuresFile(dirFeatures);
		featuresFile /= filenameWithoutExtension + Config::getInstance()->getFeaturesFile() + Config::getInstance()->getCsvExtension();
		Write2DDataToCSVFile(featuresFile.make_preferred().string(), nucleiFeatures);

		// Write Nuclei Centroids
		vector< CharImageIndexType > nucleiCentroid = cd.GetNucleiCentroids();
		path centroidFile(dirFeatures);
		centroidFile /= filenameWithoutExtension + Config::getInstance()->getCentroidsFile() + Config::getInstance()->getCsvExtension();
		Write1DIndexesToCSVFile(centroidFile.make_preferred().string(), nucleiCentroid);

		if (Config::getInstance()->writeDetectionImages()){

			// Write binary segmentation image
			CharImagePointer binImage = cd.GetResultImage();
			path detectionFile(dirSeg);
			detectionFile /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtension();
			CharFileWriting(binImage, detectionFile.make_preferred().string(), string("\nError: Binary Detection image writing ..."));

			// write nuclei overlay image
			RGBImagePointer overlayImage = DrawCircle(rgbImage, nucleiCentroid, Config::getInstance()->getOverlayColourCode());
			path overlayFile(dirSeg);
			overlayFile /= filenameWithoutExtension + Config::getInstance()->getNucleiOverlayFile() + Config::getInstance()->getImageExtensionn();
			RGBFileWriting(overlayImage, overlayFile.make_preferred().string(), string("\nError: RGB Detection image writing ..."));
		}

		for (int instanceNo = 0; instanceNo < nucleiFeatures.size(); instanceNo++)
			allFeatures.push_back(nucleiFeatures[instanceNo]);

		std::cout << "GT = " << cd.GetNumberOfGT()
			<< "\tTP = " << cd.GetNumberOfTP()
			<< "\tFN = " << cd.GetNumberOfFN()
			<< "\tFP = " << cd.GetNumberOfFP()
			<< "\tTotal = " << nucleiFeatures.size() << std::endl;
		file << "Image_Name," << filenameWithoutExtension
			<< ",GT," << cd.GetNumberOfGT()
			<< ",TP," << cd.GetNumberOfTP()
			<< ",FN," << cd.GetNumberOfFN()
			<< ",FP," << cd.GetNumberOfFP()
			<< ",Total," << nucleiFeatures.size() << std::endl;

	}
	std::cout << "GT = " << GT << "\tTP = " << TP << "\tFN = " << FN << "\tFP = " << FP << "\tTotal = " << TP + FP << std::endl;
	file << "\n\nTotal, ,GT," << GT << ",TP," << TP << ",FN," << FN << ",FP," << FP << ",Total," << TP + FP;
	file.close();

	path allFeaturesFile(dirFeatures);
	allFeaturesFile /= featureFolderName + Config::getInstance()->getCsvExtension();
	WriteFeaturesWithLabelToCSVFile(allFeaturesFile.make_preferred().string(), allFeatures);

	path allFeaturesArff(dirFeatures);
	allFeaturesArff /= featureFolderName + Config::getInstance()->getArffExtension();
	Write2DDataToArffFile(allFeaturesArff.make_preferred().string(), allFeatures);
}

// 6 - It takes an input folder (containing images and ground truth for mitosis (%_Centroid.csv)) as input and segment and compute features for each image in folder
void MitosisFeatureComputation(void)
{
	cout << "\n\n\t******** Mitosis Candidate Feature Computation ********" << endl;

	path directory(Config::getInstance()->getInputDirectoryOrFile());
	string directoryName = directory.make_preferred().string().substr(directory.make_preferred().string().find_last_of("/\\") + 1);

	if (!exists(directory)) {
		cerr << "ERROR: Directory " << directory.make_preferred().string() << " does not exist!" << endl;
		throw (runtime_error("Access error"));
	}
	if (!is_directory(directory)) {
		cout << "ERROR: " << directory.make_preferred().string() << " is not a directory!" << endl;
		throw (runtime_error("Access error"));
	}
	cout << directory.make_preferred().string() << " directory exist. Start reading images ... " << endl;

	// Create new Folder for storing segmentation files
	path dirSeg(directory);
	dirSeg /= Config::getInstance()->getSegmentationFolder();
	if (!exists(dirSeg.make_preferred().string())){
		cerr << "\nSegmentation directory doesn't exist. : " << dirSeg.make_preferred().string() << " for storing segmentation files ... \n";
		throw (runtime_error("Exist error"));
	}
	else
		cout << "\nDirectory: " << dirSeg.make_preferred().string() << " already exist for storing segmentation files ... \n";

	// Create new Folder for storing Features files
	string featureFolderName = directoryName;
	if (Config::getInstance()->regionFeatures())
		featureFolderName = featureFolderName + string("_C_Region_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	else
		featureFolderName = featureFolderName + string("_C_Patch_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius())
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	path dirFeatures(directory / featureFolderName);
	if (!exists(dirFeatures.make_preferred().string())){
		if (create_directory(dirFeatures))
			cout << "\nCreate new directory: " << dirFeatures.make_preferred().string() << " for storing features files ... \n";
	}
	else
		cout << "\nDirectory: " << dirFeatures.make_preferred().string() << " already exists for storing features files ... \n";

	unsigned int TP = 0, FP = 0, GT = 0, FN = 0;
	vector< vector< double > > allFeatures;

	std::stringstream fileResult;
	fileResult << "CandidateDetection_Results" << Config::getInstance()->getCsvExtension();
	path CandDetResultFile(directory);
	CandDetResultFile /= fileResult.str();
	std::ofstream file(CandDetResultFile.make_preferred().string().c_str());

	// Reading images in given directory
	vector<path> pathVector;
	copy(directory_iterator(directory), directory_iterator(), back_inserter(pathVector));
	sort(pathVector.begin(), pathVector.end());
	int ct = 0;
	BOOST_FOREACH(path p, pathVector) {
		path filename = p.filename();
		std::string filenameWithoutExtension = filename.string().substr(0, filename.string().find_last_of("."));
		std::string fileExtension = filename.string().substr(filename.string().find_last_of(".") + 1);
		boost::algorithm::to_lower(fileExtension);

		if (string("bmp") != fileExtension && string("tiff") != fileExtension && string("tif") != fileExtension && string("jpeg") != fileExtension && string("png") != fileExtension)
			continue;

		path binFileName = dirSeg;
		binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtension();
		if (!exists(binFileName.make_preferred().string())){
			binFileName = dirSeg;
			binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + Config::getInstance()->getImageExtensionn();
			if (!exists(binFileName.make_preferred().string())){
				binFileName = dirSeg;
				binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + std::string(".png");
				if (!exists(binFileName.make_preferred().string())){
					binFileName = dirSeg;
					binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + std::string(".jpeg");
					if (!exists(binFileName.make_preferred().string())){
						binFileName = dirSeg;
						binFileName /= filenameWithoutExtension + Config::getInstance()->getBinaryFile() + std::string(".bmp");
					}
				}
			}
		}

		RGBImagePointer rgbImage = RGBFileReading(p.make_preferred().string(), string("\n Err 1: RGB frame reading Err: Exception caught\n"));
		CharImagePointer binImage = CharFileReading(binFileName.make_preferred().string(), string("\n Err 2: Binary frame reading Err: Exception caught\n"));

		cout << "\n" << ++ct << " - Mitosis candidate feature computation in image: " << p.make_preferred().string() << endl;

		// Compare candidates with GT mitosis and find out TP (Mitosis) and FP (NonMitosis)
		std::stringstream GTCentroidsFileName;
		GTCentroidsFileName << filenameWithoutExtension << Config::getInstance()->getCentroidsFile() << Config::getInstance()->getCsvExtension();
		path GTCentroidsFile(directory / GTCentroidsFileName.str());

		// Nuclei Feature Computation
		NucleiDetector cd;
		cd.FeatureExtraction(rgbImage, binImage, filenameWithoutExtension, directory);
		cd.SetGTCentroidsFromFile(GTCentroidsFile);
		cd.ComputeCandidatesLabels(true);
		GT = GT + cd.GetNumberOfGT();
		TP = TP + cd.GetNumberOfTP();
		FN = FN + cd.GetNumberOfFN();
		FP = FP + cd.GetNumberOfFP();

		// Write Nuclei Features 
		vector< vector< double > > nucleiFeatures = cd.GetFeatures();
		path featuresFile(dirFeatures);
		featuresFile /= filenameWithoutExtension + Config::getInstance()->getFeaturesFile() + Config::getInstance()->getCsvExtension();
		Write2DDataToCSVFile(featuresFile.make_preferred().string(), nucleiFeatures);

		// Write Nuclei Centroids
		if (Config::getInstance()->writeDetectionImages()){
			vector< CharImageIndexType > nucleiCentroid = cd.GetNucleiCentroids();
			path centroidFile(dirFeatures);
			centroidFile /= filenameWithoutExtension + Config::getInstance()->getCentroidsFile() + Config::getInstance()->getCsvExtension();
			Write1DIndexesToCSVFile(centroidFile.make_preferred().string(), nucleiCentroid);
		}

		for (int instanceNo = 0; instanceNo < nucleiFeatures.size(); instanceNo++)
			allFeatures.push_back(nucleiFeatures[instanceNo]);

		std::cout << "GT = " << cd.GetNumberOfGT()
			<< "\tTP = " << cd.GetNumberOfTP()
			<< "\tFN = " << cd.GetNumberOfFN()
			<< "\tFP = " << cd.GetNumberOfFP()
			<< "\tTotal = " << nucleiFeatures.size() << std::endl;
		file << "Image_Name," << filenameWithoutExtension
			<< ",GT," << cd.GetNumberOfGT()
			<< ",TP," << cd.GetNumberOfTP()
			<< ",FN," << cd.GetNumberOfFN()
			<< ",FP," << cd.GetNumberOfFP()
			<< ",Total," << nucleiFeatures.size() << std::endl;

	}
	std::cout << "GT = " << GT << "\tTP = " << TP << "\tFN = " << FN << "\tFP = " << FP << "\tTotal = " << TP + FP << std::endl;
	file << "\n\nTotal, ,GT," << GT << ",TP," << TP << ",FN," << FN << ",FP," << FP << ",Total," << TP + FP;
	file.close();

	path allFeaturesFile(dirFeatures);
	allFeaturesFile /= featureFolderName + Config::getInstance()->getCsvExtension();
	WriteFeaturesWithLabelToCSVFile(allFeaturesFile.make_preferred().string(), allFeatures);

	path allFeaturesArff(dirFeatures);
	allFeaturesArff /= featureFolderName + Config::getInstance()->getArffExtension();
	Write2DDataToArffFile(allFeaturesArff.make_preferred().string(), allFeatures);
}

void Test(void)
{
	cout << "\n\n\t******** Mitosis Candidate Feature Computation ********" << endl;

	path directory("H:/WS/MITOS_A_TestingSet");
	string directoryName = directory.make_preferred().string().substr(directory.make_preferred().string().find_last_of("/\\") + 1);

	if (!exists(directory)) {
		cerr << "ERROR: Directory " << directory.make_preferred().string() << " does not exist!" << endl;
		throw (runtime_error("Access error"));
	}
	if (!is_directory(directory)) {
		cout << "ERROR: " << directory.make_preferred().string() << " is not a directory!" << endl;
		throw (runtime_error("Access error"));
	}
	cout << directory.make_preferred().string() << " directory exist. Start reading images ... " << endl;

	// Create new Folder for storing segmentation files
	path dirSeg(directory);
	dirSeg /= Config::getInstance()->getSegmentationFolder();
	if (!exists(dirSeg.make_preferred().string())){
		cerr << "\nSegmentation directory doesn't exist. : " << dirSeg.make_preferred().string() << " for storing segmentation files ... \n";
		throw (runtime_error("Exist error"));
	}
	else
		cout << "\nDirectory: " << dirSeg.make_preferred().string() << " already exist for storing segmentation files ... \n";

	// Create new Folder for storing Features files
	string featureFolderName = directoryName;
	if (Config::getInstance()->regionFeatures())
		featureFolderName = featureFolderName + string("_C_Region_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		+ string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius()) 
		+ string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	else
		featureFolderName = featureFolderName + string("_C_Patch_") + boost::lexical_cast<std::string>(Config::getInstance()->getGrayLevels())
		 + string("_C") + boost::lexical_cast<std::string>(Config::getInstance()->getCMRadius()) 
		 + string("_R") + boost::lexical_cast<std::string>(Config::getInstance()->getRLRadius());
	path dirFeatures(directory / featureFolderName);
	cout << "\nFeature Directory: " << dirFeatures.make_preferred().string();
	if (!exists(dirFeatures.make_preferred().string())){
		if (create_directory(dirFeatures))
			cout << "\nCreate new directory: " << dirFeatures.make_preferred().string() << " for storing features files ... \n";
	}
	else
		cout << "\nDirectory: " << dirFeatures.make_preferred().string() << " already exists for storing features files ... \n";

	vector< vector< double > > allFeatures;

	// Reading images in given directory
	vector<path> pathVector;
	copy(directory_iterator(directory), directory_iterator(), back_inserter(pathVector));
	sort(pathVector.begin(), pathVector.end());
	int ct = 0;
	BOOST_FOREACH(path p, pathVector) {
		path filename = p.filename();
		std::string filenameWithoutExtension = filename.string().substr(0, filename.string().find_last_of("."));
		std::string fileExtension = filename.string().substr(filename.string().find_last_of(".") + 1);
		boost::algorithm::to_lower(fileExtension);

		if (string("bmp") != fileExtension && string("tiff") != fileExtension && string("tif") != fileExtension && string("jpeg") != fileExtension && string("png") != fileExtension)
			continue;

		// Write Nuclei Features */
		vector< vector< double > > nucleiFeatures;// = cd.GetFeatures();
		path featuresFile(dirFeatures);
		featuresFile /= filenameWithoutExtension + Config::getInstance()->getFeaturesFile() + Config::getInstance()->getCsvExtension();
		Read2DDataFromCSVFile(featuresFile.make_preferred().string(), nucleiFeatures);

		for (int instanceNo = 0; instanceNo < nucleiFeatures.size(); instanceNo++)
			allFeatures.push_back(nucleiFeatures[instanceNo]);


	}

	path allFeaturesFile(dirFeatures);
	allFeaturesFile /= featureFolderName + Config::getInstance()->getCsvExtension();
	WriteFeaturesWithLabelToCSVFile(allFeaturesFile.make_preferred().string(), allFeatures);

	path allFeaturesArff(dirFeatures);
	allFeaturesArff /= featureFolderName + Config::getInstance()->getArffExtension();
	Write2DDataToArffFile(allFeaturesArff.make_preferred().string(), allFeatures);
}
