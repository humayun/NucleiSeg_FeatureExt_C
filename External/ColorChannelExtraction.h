#ifndef __COLORCHANNELEXTRACTION_H__
#define __COLORCHANNELEXTRACTION_H__

#include "itkTestingComparisonImageFilter.h"

#include "itkRGBToXYZColorSpacePixelAccessor.h"
#include "itkRGBToHSIColorSpacePixelAccessor.h"
#include "itkRGBToHSLColorSpacePixelAccessor.h"
#include "itkRGBToHSVColorSpacePixelAccessor.h"
#include "itkRGBToLabColorSpacePixelAccessor.h"
#include "itkRGBToLuvColorSpacePixelAccessor.h"
#include "itkRGBToYuvColorSpacePixelAccessor.h"
#include "itkRGBToYUV2ColorSpacePixelAccessor.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "ITKDeclarations.h"

#include <iostream>

CharImagePointer ColorChannelExtraction(RGBImagePointer rgbImage, std::string outputColorSpace, int colorChannel);
int RegressionTestImage (std::string testImageFilename, std::string baselineImageFilename,int reportErrors, bool differences);

#endif // __COLORCHANNELEXTRACTION_H__