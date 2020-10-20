#ifndef Mu2eUtilities_MueXSpectrum_hh
#define Mu2eUtilities_MueXSpectrum_hh

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"

// C++ includes
#include <utility>

namespace mu2e {

  double f_mu2eX_gen(double E, void *p);

  class MueXSpectrum {

    public:
      
      struct Params_t {
        double eMax;
        double mmu;
        double Emu;
        double BR;
        double me;
        double mN;
        double a0;
        double a1;
        double a2;
        double a3;
        double a4;
        double a5;
      } _par;
      
      MueXSpectrum(double maxEnergy, double bin, int RadCorrected = 0);
      
      ~MueXSpectrum(){}
   
      double getWeight(double E) const;
      double getCorrectedMueXSpectrum(double e) const ;
      double evalIntegral(double de);
      static double  f_mu2eX_gen(double E, void *p);
      void   setSpectrum   (int SpectrumType) { _spectrumType = SpectrumType; }
    
    private:

      double             _nbins;
      double             _bin;
      int                _spectrumType;// 0:delta function ; 1: rad corrected
      double             _eMax;           
      double             _me;// electron mass
      double             _integral;// over n-1 bins...
    };

} // end of namespace mu2e

#endif /* Mu2eUtilities_MueXSpectrum_hh */

  
