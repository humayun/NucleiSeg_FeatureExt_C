/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRGBToXYZColorSpacePixelAccessor.h,v $
  Language:  C++
  Date:      $Date: 2009-03-03 15:08:46 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRGBToXYZColorSpacePixelAccessor_h
#define __itkRGBToXYZColorSpacePixelAccessor_h


#include "itkRGBPixel.h"
#include "itkVector.h"
#include "vnl/vnl_math.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Accessor
{
/**
 * \class RGBToXYZColorSpacePixelAccessor
 * \brief Give access to a RGBPixel as if it were in XYZ Color Space as a Vector type.
 *
 * This class is intended to be used as parameter of 
 * an ImageAdaptor to make an RGBPixel image appear as being
 * an image of Vector pixel type in XYZ Color Space.
 *
 * \sa ImageAdaptor
 * \ingroup ImageAdaptors
 *
 */

template <class TInput, class TOutput>
class ITK_EXPORT RGBToXYZColorSpacePixelAccessor
{
public:
  /** Standard class typedefs. */
  typedef   RGBToXYZColorSpacePixelAccessor        Self;

 /** External typedef. It defines the external aspect
   * that this class will exhibit */
  typedef  Vector<TOutput,3>     ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data */
  typedef   RGBPixel<TInput>    InternalType;

  /** Write access to the RGBToXYZColorSpace component */
  inline void Set( InternalType & output, const ExternalType & input ) const
    {     
    // Normalize RGB values.
    double r = (double)input[0] / (double)NumericTraits<TInput>::max();
    double g = (double)input[1] / (double)NumericTraits<TInput>::max();
    double b = (double)input[2] / (double)NumericTraits<TInput>::max();
    
    if(r > 0.04045)
      r = vcl_pow(((r + 0.055) / 1.055), 2.4);
    else
      r = r / 12.92;
    
    if(g > 0.04045)
      g = vcl_pow(((g + 0.055) / 1.055), 2.4);
    else
      g = g / 12.92;
      
    if(b > 0.04045)
      b = vcl_pow(((b + 0.055) / 1.055), 2.4);
    else
      b = b / 12.92;
    
    r = r * 100.0;
    g = g * 100.0;
    b = b * 100.0;
      
    output[0] = static_cast<TInput>(r * 0.4124 + g * 0.3576 + b * 0.1805); // X
    output[1] = static_cast<TInput>(r * 0.2126 + g * 0.7152 + b * 0.0722); // Y
    output[2] = static_cast<TInput>(r * 0.0193 + g * 0.1192 + b * 0.9505); // Z
    }

  /** Read access to the RGBToXYZColorSpace component */
  inline ExternalType Get( const InternalType & input ) const
    {
    // Normalize RGB values.
    double r = (double)input[0] / (double)NumericTraits<TInput>::max();
    double g = (double)input[1] / (double)NumericTraits<TInput>::max();
    double b = (double)input[2] / (double)NumericTraits<TInput>::max();
    
    if(r > 0.04045)
      r = vcl_pow(((r + 0.055) / 1.055), 2.0);
    else
      r = r / 12.92;
    
    if(g > 0.04045)
      g = vcl_pow(((g + 0.055) / 1.055), 2.0);
    else
      g = g / 12.92;
      
    if(b > 0.04045)
      b = vcl_pow(((b + 0.055) / 1.055), 2.0);
    else
      b = b / 12.92;
    
    r = r * 100.0;
    g = g * 100.0;
    b = b * 100.0;
      
    ExternalType output;
    output[0] = static_cast<TOutput>(r * 0.4124 + g * 0.3576 + b * 0.1805); // X
    output[1] = static_cast<TOutput>(r * 0.2126 + g * 0.7152 + b * 0.0722); // Y
    output[2] = static_cast<TOutput>(r * 0.0193 + g * 0.1192 + b * 0.9505); // Z

    return output;
    }

private:
};
  
}  // end namespace Accessor
}  // end namespace itk

#endif
