/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRGBToLuvColorSpacePixelAccessor.h,v $
  Language:  C++
  Date:      $Date: 2009-03-03 15:08:46 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRGBToLuvColorSpacePixelAccessor_h
#define __itkRGBToLuvColorSpacePixelAccessor_h


#include "itkRGBPixel.h"
#include "itkVector.h"
#include "vnl/vnl_math.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Accessor
{
/**
 * \class RGBToLuvColorSpacePixelAccessor
 * \brief Give access to a RGBPixel as if it were in Luv (CIE L*u*v*) Color Space as a Vector type.
 *
 * This class is intended to be used as parameter of 
 * an ImageAdaptor to make an RGBPixel image appear as being
 * an image of Vector pixel type in Luv Color Space.
 *
 * \sa ImageAdaptor
 * \ingroup ImageAdaptors
 *
 */

template <class TInput, class TOutput>
class ITK_EXPORT RGBToLuvColorSpacePixelAccessor
{
public:
  /** Standard class typedefs. */
  typedef   RGBToLuvColorSpacePixelAccessor        Self;

 /** External typedef. It defines the external aspect
   * that this class will exhibit */
  typedef  Vector<TOutput,3>     ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data */
  typedef   RGBPixel<TInput>    InternalType;

  /** Write access to the RGBToLuvColorSpace component */
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
    
    double X = r * 0.4124 + g * 0.3576 + b * 0.1805;
    double Y = r * 0.2126 + g * 0.7152 + b * 0.0722;
    double Z = r * 0.0193 + g * 0.1192 + b * 0.9505;
   
    double u, v, l;
    if(X == 0 && Y == 0 && Z == 0)
      {
      u = 0;
      v = 0;
      }
    else
      {
      u = (4 * X) / (X + (15 * Y) + (3 * Z));
      v = (9 * Y) / (X + (15 * Y) + (3 * Z));
      }
    
    l = Y / 100.0;
    if(l > 0.008856)
      l = vcl_exp(vcl_log(l)/3.0);
    else
      l = ((7.787 * l) + (16.0 / 116.0));
    
    X = X / 95.047;
    Y = Y / 100.0;
    Z = Z / 108.883;
    
    double u2 = (4 * X) / (X + (15 * Y) + (3 * Z));
    double v2 = (9 * Y) / (X + (15 * Y) + (3 * Z));
      
    output[0] = static_cast<TInput>((116.0 * l) - 16.0); // L
    output[1] = static_cast<TInput>(13.0 * output[0] * (u - u2)); // u
    output[2] = static_cast<TInput>(13.0 * output[0] * (v - v2)); // v
    }

  /** Read access to the RGBToLuvColorSpace component */
  inline ExternalType Get( const InternalType & input ) const
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
    
    double X = r * 0.4124 + g * 0.3576 + b * 0.1805;
    double Y = r * 0.2126 + g * 0.7152 + b * 0.0722;
    double Z = r * 0.0193 + g * 0.1192 + b * 0.9505;
   
    double u, v, l;
    if(X == 0 && Y == 0 && Z == 0)
      {
      u = 0;
      v = 0;
      }
    else
      {
      u = (4 * X) / (X + (15 * Y) + (3 * Z));
      v = (9 * Y) / (X + (15 * Y) + (3 * Z));
      }
    
    l = Y / 100.0;
    if(l > 0.008856)
      l = vcl_exp(vcl_log(l)/3.0);
    else
      l = ((7.787 * l) + (16.0 / 116.0));
    
    X = X / 95.047;
    Y = Y / 100.0;
    Z = Z / 108.883;
    
    double u2 = (4 * X) / (X + (15 * Y) + (3 * Z));
    double v2 = (9 * Y) / (X + (15 * Y) + (3 * Z));
      
    ExternalType output;
    output[0] = static_cast<TOutput>((116.0 * l) - 16.0); // L         
    output[1] = static_cast<TOutput>(13.0 * output[0] * (u - u2)); // u
    output[2] = static_cast<TOutput>(13.0 * output[0] * (v - v2)); // v

    return output;
    }

private:
};
  
}  // end namespace Accessor
}  // end namespace itk

#endif
