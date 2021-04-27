#ifndef __Mu2eII_ana_SimpHist_t_hh
#define __Mu2eII_ana_SimpHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace Mu2eII {

  struct SimpHist_t {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fNStrawHits;               // 
    TH1F*    fMomTargetEnd;             // 
    TH1F*    fMomTrackerFront;          // 
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fCosTh;			// 
  };

}
#endif
