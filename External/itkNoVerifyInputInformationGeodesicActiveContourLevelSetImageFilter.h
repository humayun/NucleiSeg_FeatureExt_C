#ifndef __itkNoVerifyInputInformationGeodesicActiveContourLevelSetImageFilter_h
#define __itkNoVerifyInputInformationGeodesicActiveContourLevelSetImageFilter_h

#include "itkGeodesicActiveContourLevelSetImageFilter.h"

namespace itk {

template< class TInputImage,
          class TFeatureImage,
          class TOutputPixelType = float >
class ITK_EXPORT NoVerifyInputInformationGeodesicActiveContourLevelSetImageFilter:
  public GeodesicActiveContourLevelSetImageFilter< TInputImage, TFeatureImage, TOutputPixelType >
{
public:
	/** Standard class typedefs. */
	typedef NoVerifyInputInformationGeodesicActiveContourLevelSetImageFilter							Self;
	typedef GeodesicActiveContourLevelSetImageFilter< TInputImage, TFeatureImage, TOutputPixelType >	Superclass;
	typedef SmartPointer< Self >																		Pointer;
 
	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(NoVerifyInputInformationGeodesicActiveContourLevelSetImageFilter, GeodesicActiveContourLevelSetImageFilter);

protected:
	// Override VeriyInputInformation(): no need to check if inputs have the same physical space.
	virtual void VerifyInputInformation() {}
};
}


#endif /* __itkNoVerifyInputInformationGeodesicActiveContourLevelSetImageFilter_h */