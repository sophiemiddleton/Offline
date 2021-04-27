#ifndef __Mu2eII_ana_ClusterHist_t_hh
#define __Mu2eII_ana_ClusterHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace Mu2eII {

  struct ClusterHist_t {
    TH1F*    fDiskID;
    TH1F*    fEnergy;
    TH1F*    fEnergyDiff;
    TH1F*    fT0;
    TH1F*    fRow;
    TH1F*    fCol;
    TH1F*    fX;
    TH1F*    fY;
    TH1F*    fZ;
    TH1F*    fR;
    TH1F*    fNCr0;			// all clustered
    TH1F*    fNCr1;			// above 1MeV
    TH1F*    fYMean;
    TH1F*    fZMean;
    TH1F*    fSigY;
    TH1F*    fSigZ;
    TH1F*    fSigR;
    TH1F*    fFrE1;
    TH1F*    fFrE2;
    TH1F*    fSigE1;
    TH1F*    fSigE2;
  };
}
#endif
