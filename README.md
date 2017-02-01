#Nuclei Segmentation and Feature Extraction Framework (Nuclei Morphology Analysis) Framework

This framework have different function to segment nuclei and compute nuclei features including (morphology, Intensity, Co-occurrence and Run-Length texture features). This framework is used to compute nuclei features of histopathology images of Breast for PlosOne 2014 paper (Discriminate UDH and DCIS). 

PlosOne Paper (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0114885)

This framework have different functions for computing nuclei segmentation and Feature extraction like:
1- Nuclei Segmentation only
2- Nuclei Segmentation and Feature Extraction 
3- Feature Extraction only ( The segmented nuclei images should exist in the processing folder)
4- Nuclei Segmentation on Grayscale (Flouresence Images like Maximum Intensity Projection images of 3D Lightsheet Microsocpy)
5- Nuclei Feature extraction on ROI images of TCGA WSI
6- Nuclei Segmentation and Feature Extraction of ICPR MITOS Dataset 2012
7- Nuclei Feature Extraction of ICPR MITOS Dataset 2012

To build NucleiMorphology Program, two libraries are required: (ITK and Boost).

To run the NucleiMorphology Program, you need to pass following arguments:

-i “image folder location” (Input Directory: required)
-f numeric value ranging from 0 to 7 (Execution Mode: required)
•	0 – nuclei segmentation only, 
•	1 – nuclei segmentation and feature extraction, 
•	2 – nuclei feature extraction only, 
•	3 – nuclei feature extraction for grayscale image (LSM), 
•	5 – Mitosis segmentation and feature extraction, 
•	6 – Mitosis feature Extraction only, 
•	7 – feature extraction using selected color channels (4)
-l 128 (Grayscale Levels: optional,  Default is 256; other values can be 128, 64, 32 or 16)
-c (NoOfColumns: optional, Default is 1, other values can be 2, 4, or 8, number of columns for computing texture features (Co-occurrence, Run-Length))
-r  (NoOfRows: optioanl, Default is 1, other values can be 2,4, or 8, number of rows for computing texture features (Co-occurrence, Run-Length)
-p  (Feature Computation Mode: optional,  Default is region features; other value is patch based features.)
