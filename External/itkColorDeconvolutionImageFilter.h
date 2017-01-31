#ifndef __itkColorDeconvolutionImageFilter_h
#define __itkColorDeconvolutionImageFilter_h

#include "itkImageToImageFilter.h"
 
namespace itk
{
template< class TInput, class TOutput >
class ColorDeconvolutionImageFilter : public ImageToImageFilter< TInput, TOutput >
{
public:
  /** Standard class typedefs. */
  typedef ColorDeconvolutionImageFilter			Self;
  typedef ImageToImageFilter< TInput, TOutput >	Superclass;
  typedef SmartPointer< Self >					Pointer;
 
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
 
  /** Run-time type information (and related methods). */
  itkTypeMacro(ColorDeconvolutionImageFilter, ImageToImageFilter);
 
  TOutput* GetOutput1();
  TOutput* GetOutput2();
  TOutput* GetOutput3();

  /** Set Pitié-Salpêtrière H&E staining: haematoxylin and eosin
   *  output1 = haematoxylin
   *  output2 = eosin */
  void SetPSL_HEStaining();

  /** Set NUH H&E staining: haematoxylin and eosin
   *  output1 = haematoxylin
   *  output2 = eosin */
  void SetNUH_HEStaining();

  /** Set H&E staining: haematoxylin and eosin
   *  output1 = haematoxylin
   *  output2 = eosin */
  void SetHEStaining();

  /** Set H&E staining 2: haematoxylin and eosin
   *  output1 = haematoxylin
   *  output2 = eosin */
  void SetHE2Staining();
 
  /** Set H&E and DAB staining: haematoxylin, eosin and 3,3'-Diaminobenzidine
   *  output1 = haematoxylin
   *  output2 = eosin
   *  output3 = DAB */
  void SetHEDABStaining();
 
  /** Set H and DAB staining: haematoxylin and 3,3'-Diaminobenzidine
   *  output1 = haematoxylin
   *  output2 = DAB */
  void SetHDABStaining();
 
  /** Set H and AEC staining: haematoxylin and 3-amino-9-ethylcarbazole
   *  output1 = haematoxylin
   *  output2 = AEC */
  void SetHAECStaining();

  /** Set Feulgen and light green staining
   *  output1 = Feulgen
   *  output2 = light green */
  void SetFeulgenLightGreenStaining();

  /** Set Giemsa staining: methylene blue and eosin
   *  output1 = methylene blue
   *  output2 = eosin */
  void SetGiemsaStaining();

  /** Set fast red blue DAB staining: fast red, Luxol fast blue and 3,3'-Diaminobenzidine
   *  output1 = fast red
   *  output2 = Luxol fast blue
   *  output3 = DAB */
  void SetFastRedBlueDABStaining();

  /** Set methyl green DAB staining: methyl green and 3,3'-Diaminobenzidine
   *  output1 = methyl green
   *  output2 = DAB */
  void SetMethylGreenDABStaining();

  /** Set AZAN Mallory staining: aniline blue, azocarmine and orange-G
   *  output1 = aniline blue
   *  output2 = azocarmine
   *  output3 = orange-G */
  void SetAZANMalloryStaining();

  /** Set Masson trichrome staining: methyl blue and Ponceau fuchsin
   *  output1 = methyl blue
   *  output2 = Ponceau fuchsin */
  void SetMassonTrichromeStaining();

  /** Set Alcian staining: Alcian blue and haematoxylin after periodic acid-Schiff (PAS)
   *  output1 = Alcian blue
   *  output2 = haematoxylin after periodic acid-Schiff (PAS) */
  void SetAlcianStaining();

  /** Set HPAS staining: haematoxylin and periodic acid-Schiff (PAS)
   *  output1 = haematoxylin
   *  output2 = periodic acid-Schiff (PAS) */
  void SetHPASStaining();



protected:
  ColorDeconvolutionImageFilter();
  ~ColorDeconvolutionImageFilter(){}
 
  /** Does the real work. */
  virtual void GenerateData();
 
  /**  Create the Output */
  DataObject::Pointer MakeOutput(unsigned int idx);
 
private:
  ColorDeconvolutionImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  double MODx[3];
  double MODy[3];
  double MODz[3];
 
};
} //namespace ITK
 
 
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkColorDeconvolutionImageFilter.txx"
#endif

#endif /* __itkColorDeconvolutionImageFilter_h */
