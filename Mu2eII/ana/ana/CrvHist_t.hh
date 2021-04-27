#ifndef __Mu2eII_ana_CrvHist_t_hh
#define __Mu2eII_ana_CrvHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace Mu2eII {

  struct CrvClusterHist_t {
    TH1F*    fSector;
    TH1F*    fFirstBar;                           // bar # of the first pulse
    TH1F*    fNPulses;
    TH1F*    fNPe;                                // N(PE) - apparently, the sum
    TH1F*    fNPePP;                              // N(PE) per pulse
    TH1F*    fStartTime;
    TH1F*    fEndTime;
    TH1F*    fWidth;
    TH2F*    fXVsZ;
    TH2F*    fYVsZ;
    TH1F*    fCorrTime;
    TH1F*    fBarsOneEnd;
    TH1F*    fCrvPropdT;
    TH1F*    fNSectors;
    TH1F*    fBarsTwoEnd;
    TH1F*    fNDiffLSectors;
    TH1F*    fStubSlope;
    TH1F*    fStubSlopeChi2;
    TH1F*    fStubSlopeDelta;
    TH1F*    fStubQN;
    TH1F*    fStubSlopeMCProduct;  

  };

  struct CrvPulseHist_t {
    TH1F*    fNPe;
    TH1F*    fNPeHeight;
    TH1F*    fNDigis;
    TH1F*    fBar;
    TH1F*    fSipm;
    TH1F*    fTime;
    TH1F*    fHeight;
    TH1F*    fWidth;
    TH1F*    fChi2;
    TH1F*    fLeTime;
    TH1F*    fDt;
  };

  struct CrvCoincidenceHist_t {
    TH1F*    fSectorType;
    TH1F*    fNPulses;
  };

}
#endif
