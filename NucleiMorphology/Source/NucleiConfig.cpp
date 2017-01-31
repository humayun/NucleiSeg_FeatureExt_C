/*
 * NucleiConfig.cpp
 *
 *  Created on: 20 May 2014
 *      Author: Humayun
 */

#include "NucleiConfig.h"

// Global static pointer to ensure there is only one single instance of the class
Config* Config::m_instance = NULL;

Config* Config::getInstance() {
	if (!m_instance) // Make sure that only a single one instance of this class is generated
		m_instance = new Config;
	return m_instance;
}

void Config::parseCommandLine(int argc, char* argv[]) {

	// Program options
	namespace po = boost::program_options;

	// List of allowed options
	po::options_description options_general("");
	options_general.add_options()
		("help,h",										" display help message")
		("inputDirectory,i",	po::value<string>(),	" input directory (containing frames)")
		("computeFeatures,f",	po::value<int>(),		" activate feature computation (if not set will perform Nuclei Segmentation and Feature Computations)")
		("grayLevels,l",		po::value<int>(),		" number of graylevels (256, 128, 64, 32, 16, 8, 4, 2)")
		("CMRadius,c",			po::value<int>(),		" Co-occurrence Matrix Feature radius (1, 2, 4)")
		("RLRadius,r",			po::value<int>(),		" Run-Length Matrix Feature radius (1, 2, 4)")
		("regionFeatures,p",							" Patch (regionFeatures=0) or (regionFeatures=1) features")
		("write,w",										" nuclei detection frames/patches will be written");

	po::options_description options_hidden("Hidden options of configuration file");
	options_hidden.add_options()
		("RGB.nucleiDetectionChannel",	po::value<string>()->required(),	"channel name\nPossible values: BlueRatio, Red, Green, Blue, HSV, Lab, Luv, Haematoxylin")
		("RGB.HEMatrixType",			po::value<string>()->required(),	"type of H&E matrix\nPossible values: PSL, NUH, Ruitfrok")					
		("RGB.macroPixelUpsampling",	po::value<int>()->required(), 		"value for macro pixel upsampling")
		("RGB.regionFeatures",			po::value<int>()->required(),		"texture features computation on patch or region\nPossible values: Region=1, Patch=0")
		("RGB.computeFeatures",			po::value<int>()->required(),		"nucleiSegmentation (0) or featureComputation (1)")
		("RGB.writeNucleiPatch",		po::value<int>()->required(),		"write nuclei detection images (writeImages=1 or writeImages=0)")
		("RGB.writeDetectionImages",	po::value<int>()->required(),		"write nuclei detection images (writing detection images=1 or not writing detection images=0)")
		("RGB.resolutionX",				po::value<int>()->required(),		"Scanner Resolution X")
		("RGB.resolutionY",				po::value<int>()->required(),		"Scanner Resolution Y")
		("RGB.graylevels",				po::value<int>()->required(),		"No of gray levels")
		("RGB.GTRadius",				po::value<int>()->required(),		"radius size for comparing with GT centroid")
		("RGB.CMRadius",				po::value<int>()->required(),		"radius size for Co-occurrence Matrix")
		("RGB.RLRadius",				po::value<int>()->required(),		"radius size for Run-length Matrix")
		("Channel.lowerThreshold",		po::value<int>()->required(),		"lower threshold for Aperio scanner")
		("Channel.upperThreshold",		po::value<int>()->required(),		"upper threshold for Aperio scanner")
		("Channel.thresholdCellMin",	po::value<int>()->required(),		"threshold cell minimum for Aperio scanner")
		("Channel.thresholdCellMax",	po::value<int>()->required(),		"threshold cell maximum for Aperio scanner")
		("Channel.nucleiMinSize",		po::value<int>()->required(),		"nuclei min size for Aperio scanner")
		("Channel.nucleiMaxSize",		po::value<int>()->required(),		"nuclei max size for Aperio scanner")
		("overlay.radius",				po::value<int>()->required(),		"radius of overlay circles on detected nuclei")
		("overlay.thickness",			po::value<int>()->required(),		"thickness of overlay circles on detected nuclei")
		("overlay.colour",				po::value<string>()->required(),	"colour for overlay circles on detected nuclei\nPossible values: red, green, blue, yellow, cyan, orange")
		("overlay.coordinatesOrder",	po::value<string>()->required(),	"coordinates order\nPossible values: ITK, matlab");

	// All the options
	po::options_description options_all("Allowed options");
	options_all.add(options_general).add(options_hidden);

	// Command line options to be displayed to the user
	po::options_description options_command_line("Command line options");
	options_command_line.add(options_general);

	// Positional options
	po::positional_options_description pos_options;
	pos_options.add("inputDirectory", -1);

	try {
		po::store(po::parse_command_line(argc, argv, options_all), m_vm);
		po::store(po::command_line_parser(argc, argv).options(options_all).positional(pos_options).run(), m_vm);
	}
	catch (exception& e) {
		cerr << "ERROR: " << e.what() << endl << endl << "Usage: " << argv[0] << " [options] [Directory (containing frames)] " << endl;
		cerr << options_command_line << endl;
		exit(1);
	}

	// Parse config file first
	parseConfigFile(argc, argv, options_command_line, options_all);

	// Then parse command line options
	cout << "\n\nParsing command line arguments ... " << endl;

	// Help
	string entryName = string("help");
	if (m_vm.count(entryName)) {
		cout << endl << "Usage: " << argv[0] << " [options] [Directory (containing frames)] " << endl;
		cerr << options_command_line << endl;
		exit(0);
	}

	// Input input directory
	entryName = string("inputDirectory");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Missing input directory/folder" << endl;
		cerr << endl << "Usage: " << argv[0] << " [options] [Folder (containing frames)] " << endl;
		cerr << options_command_line << endl;
		exit(1);
	}
	m_inputDirectoryOrFile = m_vm[entryName].as<string>();
	cout << "\n\tInput Directory: " << m_inputDirectoryOrFile;

	// compute features
	entryName = string("computeFeatures");
	if (!m_vm.count(entryName)){
		cout << "\n\tWarning: Missing computeFeature parameter.";
		cout << "\tTake default value from config file that is: " << m_computeFeatures;
	}
	m_computeFeatures = m_vm[entryName].as<int>();

	// number of graylevels
	entryName = string("grayLevels");
	if (!m_vm.count(entryName)){
		cout << "\n\tWarning: Missing graylevels parameter.";
		cout << "\tTake default value from config file that is: " << m_graylevels;
	}
	else{
		m_graylevels = m_vm[entryName].as<int>();
		cout << "\n\tNo. of Gray Levels: " << m_graylevels;
	}

	// number of CM Radius
	entryName = string("CMRadius");
	if (!m_vm.count(entryName)){
		cout << "\n\tWarning: Missing CMRadius parameter.";
		cout << "\tTake default value from config file that is: " << m_CMRadius;
	}
	else{
		m_CMRadius = m_vm[entryName].as<int>();
		cout << "\n\tCo-occurrence Matrix Radius: " << m_CMRadius;
	}

	// number of RL Radius
	entryName = string("RLRadius");
	if (!m_vm.count(entryName)){
		cout << "\n\tWarning: Missing RLRadius parameter.";
		cout << "\tTake default value from config file that is: " << m_RLRadius;
	}
	else{
		m_RLRadius = m_vm[entryName].as<int>();
		cout << "\n\tRun-Length Matrix Radius: " << m_RLRadius;
	}

	// region or patch features
	m_regionFeatures = (0 == m_vm.count("regionFeatures") ? true : false);
	cout << (m_regionFeatures ? "\n\tRegion Features" : "\n\tPatch Features");

	// write candidate image
	m_writeDetectionImages = (0 == m_vm.count("write") ? false : true);
	cout << (m_writeNucleiPatch ? "\n\tWriting resultant images" : "\n\tNot writing resultant images");

	switch (m_computeFeatures){
	case 0:
		cout << "\n\tMode: Nuclei Segmentaion";
		break;
	case 1:
		cout << "\n\tMode: Nuclei Segmentaion and Feature Computation";
		break;
	case 2:
		cout << "\n\tMode: Nuclei Feature Computation";
		break;
	case 3:
		cout << "\n\tMode: Nuclei Feature Computation (LSM)";
		break;
	case 4:
		cout << "\n\tMode: Nuclei Feature Computation (TCGA)";
		break;
	case 5:
		cout << "\n\tMode: Mitosis Segmentation and Feature Computation";
		break;
	case 6:
		cout << "\n\tMode: Mitosis Feature Computation";
		break;
	default:
		m_computeFeatures = 2;
		cout << "\n\tMode: Nuclei Feature Computation";
	}

	initializeFeatures();
}

void Config::parseConfigFile(int argc, char* argv[], boost::program_options::options_description options_command_line, boost::program_options::options_description options_all) {

	cout << "\nReading Config file: " << m_configFilePath << " ... " << endl;

	// Program options
	namespace po = boost::program_options;

	try {
		if (!boost::filesystem::exists(	m_configFilePath.make_preferred().string())) {
			cerr << "ERROR: Configuration file '" << m_configFilePath.make_preferred().string() << "' missing"	<< endl;
			throw(runtime_error("File access error (Config::parseCommandLine)"));
		}
		ifstream config_file(m_configFilePath.make_preferred().string().c_str());
		po::store(po::parse_config_file(config_file, options_all), m_vm);
	}
	catch (exception& e) {
		cerr << "ERROR: " << e.what() << endl << "Usage: " << argv[0] << " [options] [Folder (containing frames)] " << endl;
		cerr << options_command_line << endl;
		exit(1);
	}

	po::notify(m_vm);

	// nucleiDetectionChannel
	string entryName = string("RGB.nucleiDetectionChannel");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_channelName = m_vm[entryName].as<string>();
	m_channelCode = getChannelNumberFromName(m_channelName);
	cout << "\n\tRGB.nucleiDetectionChannel: " << m_channelName;

	// HE matrix type
	entryName = string("RGB.HEMatrixType");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_HEMatrixType = getHEMatrixTypeFromName(m_vm[entryName].as<string>());
	cout << "\n\tRGB.HEMatrixType: " << m_vm[entryName].as<string>();

	// macro pixel upsampling value
	entryName = string("RGB.macroPixelUpsampling");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file "	<< m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_macroPixelUpsamplingValue = m_vm[entryName].as<int>();
	cout << "\n\tRGB.macroPixelUpsampling: " << m_macroPixelUpsamplingValue;

	// patch or region features
	entryName = string("RGB.regionFeatures");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_regionFeatures = ((m_vm[entryName].as<int>() == 0) ? false : true);
	cout << "\n\tRGB.regionFeatures: " << ((m_vm[entryName].as<int>() == 0) ? false : true);

	// nucleiDetection (0) or featureComputation (1)
	entryName = string("RGB.computeFeatures");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_computeFeatures = m_vm[entryName].as<int>();
	cout << "\n\tRGB.computeFeatures: " << m_computeFeatures;

	// write nuclei patch (writeNucleiPatch=1 or not writeNucleiPatch=0)
	entryName = string("RGB.writeNucleiPatch");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_writeNucleiPatch = ((m_vm[entryName].as<int>() == 0) ? false : true);
	cout << "\n\tRGB.writeNucleiPatch: " << ((m_vm[entryName].as<int>() == 0) ? false : true);

	// write nuclei detection images (writeDetectionImages=1 or not writeDetectionImages=0)
	entryName = string("RGB.writeDetectionImages");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_writeDetectionImages = ((m_vm[entryName].as<int>() == 0) ? false : true);
	cout << "\n\tRGB.writeDetectionImages: " << ((m_vm[entryName].as<int>() == 0) ? false : true);

	// Scanner Resolution X
	entryName = string("RGB.resolutionX");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_scannerResolutionX = m_vm[entryName].as<int>();
	m_scannerResolutionX /= 10000;
	cout << "\n\tRGB.resolutionX: " << m_scannerResolutionX;

	// Scanner Resolution Y
	entryName = string("RGB.resolutionY");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_scannerResolutionY = m_vm[entryName].as<int>();
	m_scannerResolutionY /= 10000;
	cout << "\n\tRGB.resolutionY: " << m_scannerResolutionY;

	// Number of graylevels of channel
	entryName = string("RGB.graylevels");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_graylevels = m_vm[entryName].as<int>();
	cout << "\n\tRGB.graylevels: " << m_graylevels;

	// Number of GT Radius
	entryName = string("RGB.GTRadius");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_GTRadius = m_vm[entryName].as<int>();
	cout << "\n\tRGB.GTRadius: " << m_GTRadius;

	// Radius Size for Co-occurrence Matrix
	entryName = string("RGB.CMRadius");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_CMRadius = m_vm[entryName].as<int>();
	cout << "\n\tRGB.CMRadius: " << m_CMRadius;

	// Radius Size for Run-Length Matrix
	entryName = string("RGB.RLRadius");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_RLRadius = m_vm[entryName].as<int>();
	cout << "\n\tRGB.RLRadius: " << m_RLRadius;

	/* *************** Channel *************** */
	// lower thresholds
	entryName = string("Channel.lowerThreshold");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file "	<< m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_lowerThreshold = m_vm[entryName].as<int>();
	cout << "\n\tChannel.lowerThreshold: " << m_lowerThreshold;

	// upper thresholds
	entryName = string("Channel.upperThreshold");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file "	<< m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_upperThreshold = m_vm[entryName].as<int>();
	cout << "\n\tChannel.upperThreshold: " << m_upperThreshold;

	// threshold cell min
	entryName = string("Channel.thresholdCellMin");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file "	<< m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_thresholdCellMin = m_vm[entryName].as<int>();
	cout << "\n\tChannel.thresholdCellMin: " << m_thresholdCellMin;

	// threshold cell max
	entryName = string("Channel.thresholdCellMax");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file "	<< m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_thresholdCellMax = m_vm[entryName].as<int>();
	cout << "\n\tChannel.thresholdCellMax: " << m_thresholdCellMax;

	// candidate min sizes
	entryName = string("Channel.nucleiMinSize");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_nucleiMinSize = m_vm[entryName].as<int>();
	cout << "\n\tChannel.nucleiMinSize: " << m_nucleiMinSize;

	// candidate max sizes
	entryName = string("Channel.nucleiMaxSize");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_nucleiMaxSize = m_vm[entryName].as<int>();
	cout << "\n\tChannel.nucleiMaxSize: " << m_nucleiMaxSize;

	// Overlay circle radius
	entryName = string("overlay.radius");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file "	<< m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_overlayRadius = m_vm[entryName].as<int>();
	cout << "\n\toverlay.radius: " << m_overlayRadius;

	// Overlay circle thickness
	entryName = string("overlay.thickness");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file "	<< m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_overlayThickness = m_vm[entryName].as<int>();
	cout << "\n\toverlay.thickness: " << m_overlayThickness;

	// Overlay circle colour
	entryName = string("overlay.colour");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file "	<< m_configFilePath.make_preferred().string() << endl;
		throw(logic_error("Configuration file entry missing (Config::parseCommandLine)"));
	}
	m_overlayColourCode = getOverlayColourCodeFromName(m_vm[entryName].as<string>());
	cout << "\n\toverlay.colour: " << m_overlayColourCode;

	// Coordinates order
	// ITK: (x,y)
	// matlab: (y,x)
	entryName = string("overlay.coordinatesOrder");
	if (!m_vm.count(entryName)) {
		cerr << "ERROR: Entry '" << entryName << "' missing in configuration file " << m_configFilePath.make_preferred().string() << endl;
		throw(logic_error(
			"Configuration file entry missing (Config::parseCommandLine)"));
	}
	if (m_vm[entryName].as<string>() == string("ITK"))
		m_useITKCoordinatesOrder = true;
	else if (m_vm[entryName].as<string>() == string("matlab"))
		m_useITKCoordinatesOrder = false;
	else
		throw(logic_error("Unknown value for configuration file entry overlay.coordinatesOrder (Config::parseCommandLine)"));

	cout << "\n\toverlay.coordinatesOrder: " << m_vm[entryName].as<string>();
}

int Config::getChannelNumberFromName(string channelName) {
	/*0 = Red, 1 = Green, 2 = Blue, 3 = HSV, 4 = Lab, 5 = H (H&E), 6 = BR */
	if (string("Red") == channelName)
		return 0;
	else if (string("Green") == channelName)
		return 1;
	else if (string("Blue") == channelName)
		return 2;
	else if (string("HSV") == channelName)
		return 3;
	else if (string("Lab") == channelName)
		return 4;
	else if (string("Haematoxylin") == channelName)
		return 5;
	else if (string("BlueRatio") == channelName)
		return 6;
	else {
		cerr << "ERROR: incorrect channel name '" << channelName << "'" << endl;
		throw(logic_error("Incorrect value (Config::getChannelNumberFromName)"));
	}
}

int Config::getHEMatrixTypeFromName(string name) {
	/* Antoine PSL = 1, Antoine NUH images = 2, Ruitfrok matrix = 3 */
	if (string("PSL") == name)
		return 1;
	else if (string("NUH") == name)
		return 2;
	else if (string("Ruitfrok") == name)
		return 3;
	else {
		cerr << "ERROR: incorrect classifier name '" << name << "'" << endl;
		throw(logic_error("Incorrect value (Config::getHEMatrixTypeFromName)"));
	}
}

int Config::getOverlayColourCodeFromName(string colourName) {
	/* 0 = red, 1 = green, 2 = blue, 3 = yellow, 4 = cyan, 5 = orange */
	if (string("red") == colourName)
		return 0;
	else if (string("green") == colourName)
		return 1;
	else if (string("blue") == colourName)
		return 2;
	else if (string("yellow") == colourName)
		return 3;
	else if (string("cyan") == colourName)
		return 4;
	else if (string("orange") == colourName)
		return 5;
	else {
		cerr << "ERROR: incorrect overlay colour name '" << colourName << "'" << endl;
		throw(logic_error("Incorrect value (Config::getOverlayColourCodeFromName)"));
	}
}

void Config::initializeFeatures() {

	switch (Config::getInstance()->computeFeatures()){
	case 0:
		InitializeFeaures_All();
		break;
	case 1:
		InitializeFeaures_All();
		break;
	case 2:
		InitializeFeaures_All();
		break;
	case 3:
		InitializeFeaures_LSM();
		break;
	case 4:
		InitializeFeaures_TCGA();
		break;
	default:
		InitializeFeaures_All();
	}
}

void Config::InitializeFeaures_All(void){
	// Initialize features
	m_allFeaturesName.resize(0);
	m_channelsName.resize(0);
	m_CMFeaturesName.resize(0);
	m_RLFeaturesName.resize(0);
	m_intensityFeaturesName.resize(0);
	m_morphologicalFeaturesName.resize(0);

	// Morphological Features = 9
	if( m_morphologicalFeaturesName.size() < 1 ) {
		m_morphologicalFeaturesName.push_back( "Area" );
		m_morphologicalFeaturesName.push_back( "Roundness" );
		m_morphologicalFeaturesName.push_back( "Elongation" );
		m_morphologicalFeaturesName.push_back("Flatness");
		m_morphologicalFeaturesName.push_back("Perimeter");
		m_morphologicalFeaturesName.push_back( "EqSphPerimeter" );
		m_morphologicalFeaturesName.push_back("eqSph_Radius");
		m_morphologicalFeaturesName.push_back("EqEllip_Minor");
		m_morphologicalFeaturesName.push_back("EqEllip_Major");
	}
	// Intensity Features = 5
	if( m_intensityFeaturesName.size() < 1 ) {
		m_intensityFeaturesName.push_back( "Mean" );
		m_intensityFeaturesName.push_back( "Median" );
		m_intensityFeaturesName.push_back( "StandardDev" );
		m_intensityFeaturesName.push_back( "Kurtosis" );
		m_intensityFeaturesName.push_back( "Skewness" );
	}
	// CM Features = 8
	if( m_CMFeaturesName.size() < 1 ) {
		m_CMFeaturesName.push_back( "Correlation" );
		m_CMFeaturesName.push_back( "ClusterShade" );
		m_CMFeaturesName.push_back( "ClusterProm" );
		m_CMFeaturesName.push_back( "Energy" );
		m_CMFeaturesName.push_back( "Entropy" );
		m_CMFeaturesName.push_back( "HaraCorrelation" );
		m_CMFeaturesName.push_back( "Inertia" );
		m_CMFeaturesName.push_back( "IDM" );
	}
	// RL Features = 10
	if( m_RLFeaturesName.size() < 1 ) {
		m_RLFeaturesName.push_back( "SRE" );
		m_RLFeaturesName.push_back( "LRE" );
		m_RLFeaturesName.push_back( "GLN" );
		m_RLFeaturesName.push_back( "RLN" );
		m_RLFeaturesName.push_back( "LGLRE" );
		m_RLFeaturesName.push_back( "HGLRE" );
		m_RLFeaturesName.push_back( "SRLGLE" );
		m_RLFeaturesName.push_back( "SRHGLE" );
		m_RLFeaturesName.push_back( "LRLGLE" );
		m_RLFeaturesName.push_back( "LRHGLE" );
	}
	// Channels Name = 7
	if( m_channelsName.size() < 1 ) {
		m_channelsName.push_back( "Red" );
		m_channelsName.push_back( "Ggr" );
		m_channelsName.push_back( "Blu" );
		m_channelsName.push_back( "HSV" );
		m_channelsName.push_back( "Lab" );
		//m_channelsName.push_back( "Luv" );
		m_channelsName.push_back( "HE" );
		m_channelsName.push_back( "BR" );
	}

	m_allFeaturesName.push_back("Centroid_X");
	m_allFeaturesName.push_back("Centroid_Y");
	// Copy morphological features names
	for (unsigned int i = 0; i < m_morphologicalFeaturesName.size(); i++ )
		m_allFeaturesName.push_back( m_morphologicalFeaturesName[i] );
	// Copy intensity, cooccurrence and runlength features names for each channel
	for (unsigned int i = 0; i < m_channelsName.size(); i++) {
		for (unsigned int j = 0; j < m_intensityFeaturesName.size(); j++) {
			stringstream tmp;
			tmp << m_intensityFeaturesName[j] << "_" << m_channelsName[i];
			m_allFeaturesName.push_back( tmp.str() );
		}

		for (unsigned int j = 0; j < m_CMFeaturesName.size(); j++) {
			stringstream tmp;
			tmp << m_CMFeaturesName[j] << "_" << m_channelsName[i];
			m_allFeaturesName.push_back( tmp.str() );
		}

		for (unsigned int j = 0; j < m_RLFeaturesName.size(); j++) {
			stringstream tmp;
			tmp << m_RLFeaturesName[j] << "_" << m_channelsName[i];
			m_allFeaturesName.push_back( tmp.str() );
		}
	}

	for (unsigned int k = 2; k <= getMacroPixelUpsamplingValue(); k++) {
		for (unsigned int i = 0; i < m_channelsName.size(); i++) {
			for (unsigned int j = 0; j < m_intensityFeaturesName.size(); j++) {
				stringstream tmp;
				tmp << m_intensityFeaturesName[j] << "_" << m_channelsName[i] << "_" << k;
				m_allFeaturesName.push_back(tmp.str());
			}

			for (unsigned int j = 0; j < m_CMFeaturesName.size(); j++) {
				stringstream tmp;
				tmp << m_CMFeaturesName[j] << "_" << m_channelsName[i] << "_" << k;
				m_allFeaturesName.push_back(tmp.str());
			}

			for (unsigned int j = 0; j < m_RLFeaturesName.size(); j++) {
				stringstream tmp;
				tmp << m_RLFeaturesName[j] << "_" << m_channelsName[i] << "_" << k;
				m_allFeaturesName.push_back(tmp.str());
			}
		}
	}
	m_numberOfFeatures = m_allFeaturesName.size();
	cout << "\n\t*** Number of features = " << m_numberOfFeatures << " ***" << std::endl;
}

void Config::InitializeFeaures_TCGA(void){
	// Initialize features
	m_allFeaturesName.resize(0);
	m_channelsName.resize(0);
	m_CMFeaturesName.resize(0);
	m_RLFeaturesName.resize(0);
	m_intensityFeaturesName.resize(0);
	m_morphologicalFeaturesName.resize(0);

	// Morphological Features = 9
	if (m_morphologicalFeaturesName.size() < 1) {
		m_morphologicalFeaturesName.push_back("Area");
		m_morphologicalFeaturesName.push_back("Roundness");
		m_morphologicalFeaturesName.push_back("Elongation");
		m_morphologicalFeaturesName.push_back("Flatness");
		m_morphologicalFeaturesName.push_back("Perimeter");
		m_morphologicalFeaturesName.push_back("EqSphPerimeter");
		m_morphologicalFeaturesName.push_back("eqSph_Radius");
		m_morphologicalFeaturesName.push_back("EqEllip_Minor");
		m_morphologicalFeaturesName.push_back("EqEllip_Major");
	}
	// Intensity Features = 5
	if (m_intensityFeaturesName.size() < 1) {
		m_intensityFeaturesName.push_back("Mean");
		m_intensityFeaturesName.push_back("Median");
		m_intensityFeaturesName.push_back("StandardDev");
		m_intensityFeaturesName.push_back("Kurtosis");
		m_intensityFeaturesName.push_back("Skewness");
	}
	// CM Features = 8
	if (m_CMFeaturesName.size() < 1) {
		m_CMFeaturesName.push_back("Correlation");
		m_CMFeaturesName.push_back("ClusterShade");
		m_CMFeaturesName.push_back("ClusterProm");
		m_CMFeaturesName.push_back("Energy");
		m_CMFeaturesName.push_back("Entropy");
		m_CMFeaturesName.push_back("HaraCorrelation");
		m_CMFeaturesName.push_back("Inertia");
		m_CMFeaturesName.push_back("IDM");
	}
	// RL Features = 10
	if (m_RLFeaturesName.size() < 1) {
		m_RLFeaturesName.push_back("SRE");
		m_RLFeaturesName.push_back("LRE");
		m_RLFeaturesName.push_back("GLN");
		m_RLFeaturesName.push_back("RLN");
		m_RLFeaturesName.push_back("LGLRE");
		m_RLFeaturesName.push_back("HGLRE");
		m_RLFeaturesName.push_back("SRLGLE");
		m_RLFeaturesName.push_back("SRHGLE");
		m_RLFeaturesName.push_back("LRLGLE");
		m_RLFeaturesName.push_back("LRHGLE");
	}
	// Channels Name = 7
	if (m_channelsName.size() < 1) {
		m_channelsName.push_back("Red");
		m_channelsName.push_back("Ggr");
		m_channelsName.push_back("Blu");
		m_channelsName.push_back("HSV");
		m_channelsName.push_back("Lab");
		//m_channelsName.push_back("Luv");
		m_channelsName.push_back("HE");
		m_channelsName.push_back("BR");
	}

	m_allFeaturesName.push_back("Frame_X");
	m_allFeaturesName.push_back("Frame_Y");
	m_allFeaturesName.push_back("Centroid_X");
	m_allFeaturesName.push_back("Centroid_Y");
	// Copy morphological features names
	for (unsigned int i = 0; i < m_morphologicalFeaturesName.size(); i++)
		m_allFeaturesName.push_back(m_morphologicalFeaturesName[i]);
	// Copy intensity, cooccurrence and runlength features names for each channel
	for (unsigned int i = 0; i < m_channelsName.size(); i++) {
		for (unsigned int j = 0; j < m_intensityFeaturesName.size(); j++) {
			stringstream tmp;
			tmp << m_intensityFeaturesName[j] << "_" << m_channelsName[i];
			m_allFeaturesName.push_back(tmp.str());
		}

		for (unsigned int j = 0; j < m_CMFeaturesName.size(); j++) {
			stringstream tmp;
			tmp << m_CMFeaturesName[j] << "_" << m_channelsName[i];
			m_allFeaturesName.push_back(tmp.str());
		}

		for (unsigned int j = 0; j < m_RLFeaturesName.size(); j++) {
			stringstream tmp;
			tmp << m_RLFeaturesName[j] << "_" << m_channelsName[i];
			m_allFeaturesName.push_back(tmp.str());
		}
	}

	for (unsigned int k = 2; k <= getMacroPixelUpsamplingValue(); k++) {
		for (unsigned int i = 0; i < m_channelsName.size(); i++) {
			for (unsigned int j = 0; j < m_intensityFeaturesName.size(); j++) {
				stringstream tmp;
				tmp << m_intensityFeaturesName[j] << "_" << m_channelsName[i] << "_" << k;
				m_allFeaturesName.push_back(tmp.str());
			}

			for (unsigned int j = 0; j < m_CMFeaturesName.size(); j++) {
				stringstream tmp;
				tmp << m_CMFeaturesName[j] << "_" << m_channelsName[i] << "_" << k;
				m_allFeaturesName.push_back(tmp.str());
			}

			for (unsigned int j = 0; j < m_RLFeaturesName.size(); j++) {
				stringstream tmp;
				tmp << m_RLFeaturesName[j] << "_" << m_channelsName[i] << "_" << k;
				m_allFeaturesName.push_back(tmp.str());
			}
		}
	}
	m_numberOfFeatures = m_allFeaturesName.size();
}

void Config::InitializeFeaures_LSM(void){
	// Initialize features
	m_allFeaturesName.resize(0);
	m_channelsName.resize(0);
	m_CMFeaturesName.resize(0);
	m_RLFeaturesName.resize(0);
	m_intensityFeaturesName.resize(0);
	m_morphologicalFeaturesName.resize(0);

	// Morphological Features = 5
	if (m_morphologicalFeaturesName.size() < 1) {
		m_morphologicalFeaturesName.push_back("Area");
		m_morphologicalFeaturesName.push_back("Roundness");
		m_morphologicalFeaturesName.push_back("Elongation");
		m_morphologicalFeaturesName.push_back("Flatness");
		m_morphologicalFeaturesName.push_back("Perimeter");
		m_morphologicalFeaturesName.push_back("EqSphPerimeter");
		m_morphologicalFeaturesName.push_back("eqSph_Radius");
		m_morphologicalFeaturesName.push_back("EqEllip_Minor");
		m_morphologicalFeaturesName.push_back("EqEllip_Major");
	}
	// Intensity Features = 5
	if (m_intensityFeaturesName.size() < 1) {
		m_intensityFeaturesName.push_back("Mean");
		m_intensityFeaturesName.push_back("Median");
		m_intensityFeaturesName.push_back("StandardDev");
		m_intensityFeaturesName.push_back("Kurtosis");
		m_intensityFeaturesName.push_back("Skewness");
	}
	// CM Features = 8
	if (m_CMFeaturesName.size() < 1) {
		m_CMFeaturesName.push_back("Correlation");
		m_CMFeaturesName.push_back("ClusterShade");
		m_CMFeaturesName.push_back("ClusterProm");
		m_CMFeaturesName.push_back("Energy");
		m_CMFeaturesName.push_back("Entropy");
		m_CMFeaturesName.push_back("HaraCorrelation");
		m_CMFeaturesName.push_back("Inertia");
		m_CMFeaturesName.push_back("IDM");
	}
	// RL Features = 10
	if (m_RLFeaturesName.size() < 1) {
		m_RLFeaturesName.push_back("SRE");
		m_RLFeaturesName.push_back("LRE");
		m_RLFeaturesName.push_back("GLN");
		m_RLFeaturesName.push_back("RLN");
		m_RLFeaturesName.push_back("LGLRE");
		m_RLFeaturesName.push_back("HGLRE");
		m_RLFeaturesName.push_back("SRLGLE");
		m_RLFeaturesName.push_back("SRHGLE");
		m_RLFeaturesName.push_back("LRLGLE");
		m_RLFeaturesName.push_back("LRHGLE");
	}

	m_allFeaturesName.push_back("Centroid_X");
	m_allFeaturesName.push_back("Centroid_Y");
	// Copy morphological features names
	for (unsigned int i = 0; i < m_morphologicalFeaturesName.size(); i++)
		m_allFeaturesName.push_back(m_morphologicalFeaturesName[i]);
	// Copy intensity, cooccurrence and runlength features names for each channel
	for (unsigned int j = 0; j < m_intensityFeaturesName.size(); j++) {
		stringstream tmp;
		tmp << m_intensityFeaturesName[j] << "_Gr";
		m_allFeaturesName.push_back(tmp.str());
	}

	for (unsigned int j = 0; j < m_CMFeaturesName.size(); j++) {
		stringstream tmp;
		tmp << m_CMFeaturesName[j] << "_Gr";
		m_allFeaturesName.push_back(tmp.str());
	}

	for (unsigned int j = 0; j < m_RLFeaturesName.size(); j++) {
		stringstream tmp;
		tmp << m_RLFeaturesName[j] << "_Gr";
		m_allFeaturesName.push_back(tmp.str());
	}
	m_numberOfFeatures = m_allFeaturesName.size();
}

void Config::setDefaultValue(){

	cout << "\nDefault Value reading ... \n";
	m_inputDirectoryOrFile = boost::filesystem::path("/groups/becklab/");
	m_configFilePath = boost::filesystem::path("/home/hi41/WS/NucleiDetection/configNucleiDetector.ini");

	// Configuration file RGB section
	m_channelCode = 0;
	m_channelName = "Red";
	m_HEMatrixType = 1;
	m_macroPixelUpsamplingValue = 1;
	m_scannerResolutionX = 0.2456;
	m_scannerResolutionY = 0.2456;
	m_graylevels = 256;

	// Configuration file Nuclei Detection
	m_lowerThreshold = 0;
	m_upperThreshold = 110;
	m_thresholdCellMin = 200;
	m_thresholdCellMax = 2000;
	m_nucleiMinSize = 200;
	m_nucleiMaxSize = 2000;

	// Configuration file ovelay section
	m_overlayRadius = 30;
	m_overlayThickness = 9;
	m_overlayColourCode = 3;
	m_useITKCoordinatesOrder = true;

	// Different mode
	m_computeFeatures = 2;
	m_regionFeatures = true;
	m_writeNucleiPatch = false;
	m_writeDetectionImages = false;

	initializeFeatures();
}