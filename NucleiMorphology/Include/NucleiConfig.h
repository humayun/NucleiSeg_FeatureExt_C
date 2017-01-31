/*
 * NucleiConfig.h
 *
 *  Created on: 20 May 2014
 *      Author: Humayun
 */

#ifndef __NUCLEICONFIG_H__
#define __NUCLEICONFIG_H__

#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/exception/all.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost::filesystem;
using namespace std;

class Config {

public:
	static Config* getInstance();
	void parseCommandLine(int argc, char* argv[]);
	void setDefaultValue();
	string getVersion()		{ return m_version; };

	// Configuration file RGB section
	string getDetectionChannelName()	{ return m_channelName; };
	int    getDetectionChannelCode()	{ return m_channelCode; };
	void setDetectionChannelName(string channelName) {
		m_channelName = channelName;
		m_channelCode = getChannelNumberFromName(m_channelName);
	};

	int getHEMatrixType()				{ return m_HEMatrixType; };
	int getMacroPixelUpsamplingValue()	{ return m_macroPixelUpsamplingValue; };

	int getLowerThreshold()				{ return m_lowerThreshold; };
	int getUpperThreshold()				{ return m_upperThreshold; };
	int getThresholdCellMin()			{ return m_thresholdCellMin; };
	int getThresholdCellMax()			{ return m_thresholdCellMax; };
	int getNucleiMinSize()				{ return m_nucleiMinSize; };
	int getNucleiMaxSize()				{ return m_nucleiMaxSize; };

	// Configuration file overlay section
	int getOverlayRadius()				{ return m_overlayRadius; };
	int getOverlayThickness()			{ return m_overlayThickness; };
	int getOverlayColourCode()			{ return m_overlayColourCode; };
	bool useITKCoordinatesOrder()		{ return m_useITKCoordinatesOrder; };

	// Methods for mode
	int    computeFeatures()			{ return m_computeFeatures; };
	bool   regionFeatures()				{ return m_regionFeatures; };
	bool   writeNucleiPatch()			{ return m_writeNucleiPatch; };
	bool   writeDetectionImages()		{ return m_writeDetectionImages; };
	string getChannelName()				{ return m_channelName; };

	// Members for features
	unsigned int getAllFeaturesNameSize()	{ return m_numberOfFeatures; };
	string  getAllFeaturesName(int i)	{ return m_allFeaturesName[i]; };

	int     getImagePatchSize()			{ return m_imagePatchSize; };
	int     getForegroundPixel()		{ return m_foregroundPixel; };
	int     getBackgroundPixel()		{ return m_backgroundPixel; };

	int		getGrayLevels()				{ return m_graylevels; }
	int     getCMRadius()				{ return m_CMRadius; };
	int     getRLRadius()				{ return m_RLRadius; };

	int		getGTRadius()				{ return m_GTRadius; };
	float	getScannerResolutionX()		{ return m_scannerResolutionX; };
	float	getScannerResolutionY()		{ return m_scannerResolutionY; };

	string  getImageExtension()			{ return m_imageExtension; };
	string  getImageExtensionn()		{ return m_imageExtensionn; };
	string  getCsvExtension()			{ return m_csvExtension; };
	string  getArffExtension()			{ return m_arffExtension; };
	string  getXmlExtension()			{ return m_xmlExtension; };
	string  getFeaturesFile()			{ return m_featuresFile; };
	string  getCentroidsFile()			{ return m_centroidsFile; };
	string  getSegmentationFolder()		{ return m_segmentationFolder; };
	string  getFeaturesFolder()			{ return m_featuresFolder; };
	string  getBinaryFile()				{ return m_binaryFile; };
	string  getNucleiOverlayFile()		{ return m_nucleiOverlayFile; };

	path    getInputDirectoryOrFile()	{ return m_inputDirectoryOrFile; };


	void	setGrayLevels(int graylevel){ m_graylevels=graylevel; }
	void    setCMRadius(int CMRadius)	{ m_CMRadius=CMRadius; };
	void    setRLRadius(int RLRadius)	{ m_RLRadius=RLRadius; };

	void	setGTRadius(float GTRadius)	{ m_GTRadius=GTRadius; };
	void	setScannerResolutionX(float x)	{ m_scannerResolutionX = x; };
	void	setScannerResolutionY(float y)	{ m_scannerResolutionY = y; };

	void ImageVariablesInitialization() {
	  cout << "\nInput Directory: " << getInputDirectoryOrFile()
		<< "\nSelected Channel: " << getDetectionChannelName()
		<< "\nLow-Threshold: " << getLowerThreshold()
		<< "\tUp-Threshold: " << getUpperThreshold()
		<< "\nCandidate Min Size: " << getNucleiMinSize()
		<< "\tCandidate Max Size: " << getNucleiMaxSize()
		<< "\nMode: " << (computeFeatures() ? (regionFeatures() ? "Feature Computation (Region based)" : "Feature Computation (Patch based)") : "Nuclei Detection")
		<< (writeNucleiPatch() ? "\nWrite nuclei patch" : "")
		<< (writeDetectionImages() ? "\nWrite nuclei detection images.\n" : ".\n");
	};

private:
  Config() { // Constructor is private so that it cannot be called
    // Initialise internal members
	m_version = "1.0";
    m_channelName = "";
	m_configFilePath = boost::filesystem::path("H:/WS/NucleiDetection/configNucleiDetector.ini");

    // Configuration file RGB section
    m_channelCode = -1;
    m_HEMatrixType = -1;
    m_macroPixelUpsamplingValue = -1;

	// Different mode
	m_computeFeatures = 2;		// 0 = Nuclei Segmentation only, 1 = Nuclei Segmentation and Feature Computation, 2 = Nuclei Feature Computation, 3 = Features on LSM, 4 = Features on TCGA
	m_regionFeatures = false;
	m_writeNucleiPatch = false;
	m_writeDetectionImages = false;

	m_graylevels = 256;
	m_RLRadius = 1;
	m_CMRadius = 1;

	m_GTRadius = 8;
	m_scannerResolutionX = 0.2456;	// Aperio Resolution X
	m_scannerResolutionY = 0.2456;	// Aperio Resolution Y

    // Configuration file Nuclei Detection
    m_lowerThreshold		= -1;
    m_upperThreshold		= -1;
    m_thresholdCellMin		= -1;
    m_thresholdCellMax		= -1;
    m_nucleiMinSize			= -1;
    m_nucleiMaxSize			= -1;

    // Configuration file ovelay section
    m_overlayRadius = -1;
    m_overlayThickness = -1;
    m_overlayColourCode = -1;
    m_useITKCoordinatesOrder = true;

	m_imageExtension		= string(".tif");
	m_imageExtensionn		= string(".tiff");
	m_csvExtension			= string(".csv");
	m_arffExtension			= string(".arff");
	m_xmlExtension			= string(".xml");
	m_featuresFile			= string("_Features");
	m_centroidsFile			= string("_Centroid");
	m_nucleiOverlayFile		= string("_Overlay");
	m_segmentationFolder	= string("Segmentation");
	m_featuresFolder		= string("Features");
	m_binaryFile			= string("_Binary");
  };

  void parseConfigFile( int argc, char* argv[],
      			        boost::program_options::options_description options_command_line,
              			boost::program_options::options_description options_all );

  void initializeFeatures();
  void InitializeFeaures_All();
  void InitializeFeaures_LSM();
  void InitializeFeaures_TCGA();
  static Config *m_instance;
  boost::program_options::variables_map m_vm;
  path m_configFilePath;
  string m_version;

  // Members initialised from configuration file

  // Parameters for RGB section
  int m_channelCode;
  string m_channelName;
  int getChannelNumberFromName(string channelName);

  int m_HEMatrixType;
  int getHEMatrixTypeFromName(string name);
  int m_macroPixelUpsamplingValue;
  bool m_regionFeatures;
  int m_computeFeatures;
  bool m_writeNucleiPatch;
  bool m_writeDetectionImages;

  // Parameters for color channel
  int m_lowerThreshold;
  int m_upperThreshold;
  int m_thresholdCellMin;
  int m_thresholdCellMax;
  int m_nucleiMinSize;
  int m_nucleiMaxSize;

  int m_graylevels;
  int m_CMRadius;	// Radius of Co-occurrence Matrix 
  int m_RLRadius;	// Radius of Run-length Matrix

  float m_GTRadius;     // Radius of GT center
  float m_scannerResolutionX;
  float m_scannerResolutionY;

  static const int m_imagePatchSize = 70;    // 70*70
  static const int m_foregroundPixel = 255;
  static const int m_backgroundPixel = 0;

  // Members for features
  int m_numberOfFeatures;  //170 = 9+(7*23), 492 = 5+(3*7*23)  // Modified by the code

  // Parameters for overlay section
  int m_overlayRadius;
  int m_overlayThickness;
  int m_overlayColourCode;
  int getOverlayColourCodeFromName( string colourName );
  bool m_useITKCoordinatesOrder;


  path m_inputDirectoryOrFile;

  vector< string > m_allFeaturesName;
  vector< string > m_morphologicalFeaturesName;
  vector< string > m_intensityFeaturesName;
  vector< string > m_CMFeaturesName;
  vector< string > m_RLFeaturesName;
  vector< string > m_channelsName;

  string m_imageExtension;
  string m_imageExtensionn;
  string m_csvExtension;
  string m_arffExtension;
  string m_xmlExtension;
  string m_featuresFile;
  string m_centroidsFile;
  string m_nucleiOverlayFile;
  string m_segmentationFolder;
  string m_featuresFolder;
  string m_binaryFile;

};

#endif // __NUCLEICONFIG_H__
