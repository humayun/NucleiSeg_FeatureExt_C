/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRGBToYuvColorSpacePixelAccessor.h,v $
  Language:  C++
  Date:      $Date: 2009-03-03 15:08:46 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRGBToYuvColorSpacePixelAccessor_h
#define __itkRGBToYuvColorSpacePixelAccessor_h


#include "itkRGBPixel.h"
#include "itkVector.h"
#include "vnl/vnl_math.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Accessor
{
/**
 * \class RGBToYuvColorSpacePixelAccessor
 * \brief Give access to a RGBPixel as if it were in Yuv Color Space as a Vector type.
 *
 * This class is intended to be used as parameter of 
 * an ImageAdaptor to make an RGBPixel image appear as being
 * an image of Vector pixel type in Yuv Color Space.
 *
 * \sa ImageAdaptor
 * \ingroup ImageAdaptors
 *
 */

template <class TInput, class TOutput>
class ITK_EXPORT RGBToYuvColorSpacePixelAccessor
{
public:
  /** Standard class typedefs. */
  typedef   RGBToYuvColorSpacePixelAccessor        Self;

 /** External typedef. It defines the external aspect
   * that this class will exhibit */
  typedef  Vector<TOutput,3>     ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data */
  typedef   RGBPixel<TInput>    InternalType;

  /** Write access to the RGBToYuvColorSpace component */
  inline void Set( InternalType & output, const ExternalType & input ) const
    { 
    // Normalize RGB values.
    double r = (double)input[0] / (double)NumericTraits<TInput>::max();
    double g = (double)input[1] / (double)NumericTraits<TInput>::max();
    double b = (double)input[2] / (double)NumericTraits<TInput>::max();
       
    double Y = (r + g + b) / 3.0;
    double u, v;
    
    if(vcl_abs(Y) > 0.05) // tolerance level above 0
      {
	  u = 3.0 * (r - b) / (2.0 * Y);
	  v = (double) vcl_sqrt(3.0) * (2.0 * g - r - b) / 2.0 * Y;	  
      }
    else
      {
	  u = 0.0;
	  v = 0.0;    
	  }

    output[0] = static_cast<TInput>(Y); // Y
    output[1] = static_cast<TInput>(u); // u
    output[2] = static_cast<TInput>(v); // v
    
    return output;
    }

  /** Read access to the RGBToYuvColorSpace component */
  inline ExternalType Get( const InternalType & input ) const
    {
    // Normalize RGB values.
    double r = (double)input[0] / (double)NumericTraits<TInput>::max();
    double g = (double)input[1] / (double)NumericTraits<TInput>::max();
    double b = (double)input[2] / (double)NumericTraits<TInput>::max();
       
    double Y = (r + g + b) / 3.0;
    double u, v;
    
    if(vcl_abs(Y) > 0.05) // tolerance level above 0
      {
	  u = 3.0 * (r - b) / 2.0 * Y;
	  v = (double) vcl_sqrt(3.0) * (2.0 * g - r - b) / 2.0 * Y;	  
      }
    else
      {
	  u = 0.0;
	  v = 0.0;    
	  }
      
    ExternalType output;
    output[0] = static_cast<TOutput>(Y); // Y
    output[1] = static_cast<TOutput>(u); // u
    output[2] = static_cast<TOutput>(v); // v
    
    return output;
    }

private:
};
  
}  // end namespace Accessor
}  // end namespace itk

#endif
