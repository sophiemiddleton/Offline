#ifndef __Mu2eII_ana_CrvStubPar_t__
#define __Mu2eII_ana_CrvStubPar_t__

#include "Stntuple/ana/ParBase_t.hh"

namespace Mu2eII {
class CrvStubPar_t : public ParBase_t {
public:
  float     fTime;
  float     fZ;

  float     fCorrTime;
  float     fCorrTimeTof;
  float     fCorrTimeProp;

  int       fSector;			// sector number (not type)
  int       fFirstBar;			// bar # of the first pulse
  int       fTwoEndBars;
  int       fTotalBars;
  int       fNSectors;
  int       fNDiffLSectors;
  int       fStubQN;

  float     fXCorrected;
  float     fNPePP;

  float     fTCorrAana;                 // time with all analysis-based corrections (evolving)

  float     fStubDYDZ;
  float     fStubSlopeChi2;
  float     fStubDYDZMC;
  float     fStubSlopeMCProduct;
};
}
#endif
