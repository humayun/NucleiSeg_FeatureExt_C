#ifndef __itkcolorDeconvolutionImageFilter_txx
#define __itkcolorDeconvolutionImageFilter_txx

#include "itkColorDeconvolutionImageFilter.h"

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{

template< class TInput, class TOutput >
ColorDeconvolutionImageFilter<TInput, TOutput>::ColorDeconvolutionImageFilter()
{
  this->SetNumberOfRequiredOutputs(3);
  this->SetNumberOfRequiredInputs(1);

  this->SetNthOutput( 0, this->MakeOutput(0) );
  this->SetNthOutput( 1, this->MakeOutput(1) );
  this->SetNthOutput( 2, this->MakeOutput(2) );
}

template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::GenerateData()
{
	double cosx[3];
	double cosy[3];
	double cosz[3];
	double len[3];

	for (int i = 0; i < 3; i++) {
		// normalise vector length
		cosx[i] = 0.0;
		cosy[i] = 0.0;
		cosz[i] = 0.0;
		len[i] = sqrt( MODx[i]*MODx[i] + MODy[i]*MODy[i] + MODz[i]*MODz[i] );
		if (len[i] != 0.0) {
			cosx[i] = MODx[i]/len[i];
			cosy[i] = MODy[i]/len[i];
			cosz[i] = MODz[i]/len[i];
		}
	}

	// translation matrix
	if (cosx[1] == 0.0) { // 2nd colour is unspecified
		if (cosy[1] == 0.0) {
			if (cosz[1] == 0.0) {
				cosx[1] = cosz[0];
				cosy[1] = cosx[0];
				cosz[1] = cosy[0];
			}
		}
	}

	if (cosx[2] == 0.0) { // 3rd colour is unspecified
		if (cosy[2] == 0.0) {
			if (cosz[2] == 0.0) {
				if ( (cosx[0]*cosx[0] + cosx[1]*cosx[1]) > 1.0 ) {
//					std::cerr << "Colour_3 has a negative R component." << std::endl;
					cosx[2] = 0.0;
				}
				else {
					cosx[2] = sqrt( 1.0 - (cosx[0]*cosx[0]) - (cosx[1]*cosx[1]) );
				}

				if ( (cosy[0]*cosy[0] + cosy[1]*cosy[1]) > 1.0 ) {
//					std::cerr << "Colour_3 has a negative G component." << std::endl;
					cosy[2] = 0.0;
				}
				else {
					cosy[2] = sqrt( 1.0 - (cosy[0]*cosy[0]) - (cosy[1]*cosy[1]) );
				}

				if ( (cosz[0]*cosz[0] + cosz[1]*cosz[1]) > 1.0 ) {
//					std::cerr << "Colour_3 has a negative B component." << std::endl;
					cosz[2] = 0.0;
				}
				else {
					cosz[2] = sqrt( 1.0 - (cosz[0]*cosz[0]) - (cosz[1]*cosz[1]) );
				}
			}
		}
	}

	double leng = sqrt( cosx[2]*cosx[2] + cosy[2]*cosy[2] + cosz[2]*cosz[2] );

	cosx[2] = cosx[2]/leng;
	cosy[2] = cosy[2]/leng;
	cosz[2] = cosz[2]/leng;

	for (int i=0; i<3; i++) {
		if (cosx[i] == 0.0) cosx[i] = 0.001;
		if (cosy[i] == 0.0) cosy[i] = 0.001;
		if (cosz[i] == 0.0) cosz[i] = 0.001;
	}

	// matrix inversion
	double A = cosy[1] - cosx[1] * cosy[0] / cosx[0];
	double V = cosz[1] - cosx[1] * cosz[0] / cosx[0];
	double C = cosz[2] - cosy[2] * V/A + cosx[2] * (V/A * cosy[0] / cosx[0] - cosz[0] / cosx[0]);
	double q[9];
	q[2] = (-cosx[2] / cosx[0] - cosx[2] / A * cosx[1] / cosx[0] * cosy[0] / cosx[0] + cosy[2] / A * cosx[1] / cosx[0]) / C;
	q[1] = -q[2] * V / A - cosx[1] / (cosx[0] * A);
	q[0] = 1.0 / cosx[0] - q[1] * cosy[0] / cosx[0] - q[2] * cosz[0] / cosx[0];
	q[5] = (-cosy[2] / A + cosx[2] / A * cosy[0] / cosx[0]) / C;
	q[4] = -q[5] * V / A + 1.0 / A;
	q[3] = -q[4] * cosy[0] / cosx[0] - q[5] * cosz[0] / cosx[0];
	q[8] = 1.0 / C;
	q[7] = -q[8] * V / A;
	q[6] = -q[7] * cosy[0] / cosx[0] - q[8] * cosz[0] / cosx[0];


	// translate ------------------
	// Setup output 1
	typename TOutput::Pointer output1 = this->GetOutput1();
	output1->SetRegions( this->GetInput()->GetRequestedRegion() );
	output1->Allocate();
	// Setup output 2
	typename TOutput::Pointer output2 = this->GetOutput2();
	output2->SetRegions( this->GetInput()->GetRequestedRegion() );
	output2->Allocate();
	// Setup output 3
	typename TOutput::Pointer output3 = this->GetOutput3();
	output3->SetRegions( this->GetInput()->GetRequestedRegion() );
	output3->Allocate();

	typedef itk::ImageRegionConstIterator< TInput >	ConstInputIteratorType;
	typedef itk::ImageRegionIterator< TOutput >		OutputIteratorType;

	ConstInputIteratorType	inputIter	( this->GetInput(), this->GetInput()->GetRequestedRegion() );
	OutputIteratorType		output1Iter	( output1, this->GetInput()->GetRequestedRegion() );
	OutputIteratorType		output2Iter	( output2, this->GetInput()->GetRequestedRegion() );
	OutputIteratorType		output3Iter	( output3, this->GetInput()->GetRequestedRegion() );

	for( inputIter.GoToBegin(), output1Iter.GoToBegin(), output2Iter.GoToBegin(), output3Iter.GoToBegin();
			!inputIter.IsAtEnd(); ++inputIter, ++output1Iter, ++output2Iter, ++output3Iter ) {

		// log transform the RGB data
		typename TInput::PixelType pixel = inputIter.Get();
		double log255 = log( 255.0 );
		double Rlog = -((255.0*log(((double)pixel[0]+1)/255.0))/log255); // red
		double Glog = -((255.0*log(((double)pixel[1]+1)/255.0))/log255); // green
		double Blog = -((255.0*log(((double)pixel[2]+1)/255.0))/log255); // blue
		for( int i = 0; i < 3; ++i ) {
			// i=0 ==> hematoxylin / i=1 ==> eosin / i=2 ==> background
			// rescale to match original paper values
			double Rscaled = Rlog * q[i*3];
			double Gscaled = Glog * q[i*3+1];
			double Bscaled = Blog * q[i*3+2];
			double output = exp(-((Rscaled + Gscaled + Bscaled) - 255.0) * log255 / 255.0);
			if( output > 255.0 )
				output = 255.0;

			unsigned char value = (unsigned char)(floor(output+0.5));
			switch( i ) {
			case 0: // hematoxylin
				output1Iter.Set( value );
				break;
			case 1: // eosin
				output2Iter.Set( value );
				break;
			case 2: // background
				output3Iter.Set( value );
				break;
			default:
				std::cerr << "ERROR: Value " << i << " is out of bounds [0-2]" << std::endl;
				throw ( std::runtime_error("Value out of bounds") );
			}
		}
	}
}


template< class TInput, class TOutput >
DataObject::Pointer ColorDeconvolutionImageFilter<TInput, TOutput>::MakeOutput(unsigned int idx)
{
  DataObject::Pointer output;

  switch ( idx )
    {
    case 0:
      output = ( TOutput::New() ).GetPointer();
      break;
    case 1:
      output = ( TOutput::New() ).GetPointer();
      break;
    case 2:
      output = ( TOutput::New() ).GetPointer();
      break;
    default:
      std::cerr << "No output " << idx << std::endl;
      output = NULL;
      break;
    }
  return output.GetPointer();
}

template< class TInput, class TOutput >
TOutput* ColorDeconvolutionImageFilter<TInput, TOutput>::GetOutput1()
{
  return dynamic_cast< TOutput * >(
           this->ProcessObject::GetOutput(0) );
}

template< class TInput, class TOutput >
TOutput* ColorDeconvolutionImageFilter<TInput, TOutput>::GetOutput2()
{
  return dynamic_cast< TOutput * >(
           this->ProcessObject::GetOutput(1) );
}

template< class TInput, class TOutput >
TOutput* ColorDeconvolutionImageFilter<TInput, TOutput>::GetOutput3()
{
  return dynamic_cast< TOutput * >(
           this->ProcessObject::GetOutput(2) );
}



/** Set Pitié-Salpêtrière H&E staining: haematoxylin and eosin
 *  output1 = haematoxylin
 *  output2 = eosin */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetPSL_HEStaining()
{
	// GL Haematoxylin matrix
	MODx[0]= 0.8249;
	MODy[0]= 0.8040;
	MODz[0]= 0.5969;
	// GL Eosin matrix
    MODx[1]= 0.0825;	// GL E matrix 
	MODy[1]= 0.5454;
	MODz[1]= 0.3680;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set NUH H&E staining: haematoxylin and eosin
 *  output1 = haematoxylin
 *  output2 = eosin */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetNUH_HEStaining()
{
	// GL Haematoxylin matrix
	MODx[0]= 0.650;
	MODy[0]= 0.704;
	MODz[0]= 0.286;
	// GL Eosin matrix
	MODx[1]= 0.072;
	MODy[1]= 0.990;
	MODz[1]= 0.105;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set H&E staining: haematoxylin and eosin
 *  output1 = haematoxylin
 *  output2 = eosin */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetHEStaining()
{
	// GL Haematoxylin matrix
	MODx[0]= 0.644211; //0.650;
	MODy[0]= 0.716556; //0.704;
	MODz[0]= 0.266844; //0.286;
	// GL Eosin matrix
	MODx[1]= 0.092789; //0.072;
	MODy[1]= 0.954111; //0.990;
	MODz[1]= 0.283111; //0.105;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set H&E staining2: haematoxylin and eosin
 *  output1 = haematoxylin
 *  output2 = eosin */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetHE2Staining()
{
	// GL Haematoxylin matrix
	MODx[0]= 0.49015734;
	MODy[0]= 0.76897085;
	MODz[0]= 0.41040173;
	// GL Eosin matrix
	MODx[1]= 0.04615336;
	MODy[1]= 0.8420684;
	MODz[1]= 0.5373925;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set H&E and DAB staining: haematoxylin, eosin and DAB (3,3'-Diaminobenzidine)
 *  output1 = haematoxylin
 *  output2 = eosin
 *  output3 = DAB */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetHEDABStaining()
{
	// Haematoxylin matrix
	MODx[0]= 0.650;
	MODy[0]= 0.704;
	MODz[0]= 0.286;
	// Eosin matrix
	MODx[1]= 0.072;
	MODy[1]= 0.990;
	MODz[1]= 0.105;
	// DAB matrix
	MODx[2]= 0.268;
	MODy[2]= 0.570;
	MODz[2]= 0.776;
}

/** Set H and DAB staining: haematoxylin and DAB (3,3'-Diaminobenzidine)
 *  output1 = haematoxylin
 *  output2 = DAB */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetHDABStaining()
{
	// Haematoxylin matrix
	MODx[0]= 0.650;
	MODy[0]= 0.704;
	MODz[0]= 0.286;
	// DAB matrix
	MODx[1]= 0.268;
	MODy[1]= 0.570;
	MODz[1]= 0.776;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set H and AEC staining: haematoxylin and 3-amino-9-ethylcarbazole
 *  output1 = haematoxylin
 *  output2 = AEC */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetHAECStaining()
{
	// Haematoxylin matrix
	MODx[0]= 0.650;
	MODy[0]= 0.704;
	MODz[0]= 0.286;
	// AEC matrix
	MODx[1]= 0.2743;
	MODy[1]= 0.6796;
	MODz[1]= 0.6803;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set Feulgen and light green staining
 *  output1 = Feulgen
 *  output2 = light green */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetFeulgenLightGreenStaining()
{
	//Feulgen
	MODx[0]= 0.46420921;
	MODy[0]= 0.83008335;
	MODz[0]= 0.30827187;
	// light green
	MODx[1]= 0.94705542;
	MODy[1]= 0.25373821;
	MODz[1]= 0.19650764;
	// Zero matrix
	MODx[2]= 0.0; // 0.0010000
	MODy[2]= 0.0; // 0.47027777
	MODz[2]= 0.0; //0.88235928
}

/** Set Giemsa staining: methylene blue and eosin
 *  output1 = methylene blue
 *  output2 = eosin */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetGiemsaStaining()
{
	// GL methylene blue
	MODx[0]= 0.834750233;
	MODy[0]= 0.513556283;
	MODz[0]= 0.196330403;
	// GL eosin matrix
	MODx[1]= 0.092789;
	MODy[1]= 0.954111;
	MODz[1]= 0.283111;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set fast red blue DAB staining: fast red, Luxol fast blue and 3,3'-Diaminobenzidine
 *  output1 = fast red
 *  output2 = Luxol fast blue
 *  output3 = DAB */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetFastRedBlueDABStaining()
{
	// fast red
	MODx[0]= 0.21393921;
	MODy[0]= 0.85112669;
	MODz[0]= 0.47794022;
	// Luxol fast blue
	MODx[1]= 0.74890292;
	MODy[1]= 0.60624161;
	MODz[1]= 0.26731082;
	// DAB
	MODx[2]= 0.268;
	MODy[2]= 0.570;
	MODz[2]= 0.776;
}

/** Set methyl green DAB staining: methyl green and 3,3'-Diaminobenzidine
 *  output1 = methyl green
 *  output2 = DAB */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetMethylGreenDABStaining()
{
	// methyl green matrix (GL)
	MODx[0]= 0.98003;
	MODy[0]= 0.144316;
	MODz[0]= 0.133146;
	// DAB matrix
	MODx[1]= 0.268;
	MODy[1]= 0.570;
	MODz[1]= 0.776;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set AZAN Mallory staining: aniline blue, azocarmine and orange-G
 *  output1 = aniline blue
 *  output2 = azocarmine
 *  output3 = orange-G */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetAZANMalloryStaining()
{
	// GL Aniline Blue
	MODx[0]= .853033;
	MODy[0]= .508733;
	MODz[0]= .112656;
	// GL Azocarmine
	MODx[1]=0.09289875;
	MODy[1]=0.8662008;
	MODz[1]=0.49098468;
	// GL Orange-G
	MODx[2]=0.10732849;
	MODy[2]=0.36765403;
	MODz[2]=0.9237484;
}

/** Set Masson trichrome staining: methyl blue and Ponceau fuchsin
 *  output1 = methyl blue
 *  output2 = Ponceau fuchsin */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetMassonTrichromeStaining()
{
	// GL Methyl blue
	MODx[0]=0.7995107;
	MODy[0]=0.5913521;
	MODz[0]=0.10528667;
	// GL Ponceau Fuchsin has 2 hues, really this is only approximate
	MODx[1]=0.09997159;
	MODy[1]=0.73738605;
	MODz[1]=0.6680326;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set Alcian staining: Alcian blue and haematoxylin after periodic acid-Schiff (PAS)
 *  output1 = Alcian blue
 *  output2 = haematoxylin after periodic acid-Schiff (PAS) */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetAlcianStaining()
{
	// GL Alcian Blue matrix
	MODx[0]= 0.874622;
	MODy[0]= 0.457711;
	MODz[0]= 0.158256;
	// GL Haematoxylin after PAS matrix
	MODx[1]= 0.552556;
	MODy[1]= 0.7544;
	MODz[1]= 0.353744;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}

/** Set HPAS staining: haematoxylin and periodic acid-Schiff (PAS)
 *  output1 = haematoxylin
 *  output2 = periodic acid-Schiff (PAS) */
template< class TInput, class TOutput >
void ColorDeconvolutionImageFilter<TInput, TOutput>::SetHPASStaining()
{
	// GL Haematoxylin matrix
	MODx[0]= 0.644211; //0.650;
	MODy[0]= 0.716556; //0.704;
	MODz[0]= 0.266844; //0.286;
	// GL PAS matrix
	MODx[1]= 0.175411;
	MODy[1]= 0.972178;
	MODz[1]= 0.154589;
	// Zero matrix
	MODx[2]= 0.0;
	MODy[2]= 0.0;
	MODz[2]= 0.0;
}




}// end namespace

#endif /* __itkcolorDeconvolutionImageFilter_txx */
