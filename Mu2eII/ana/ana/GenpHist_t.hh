#ifndef __Mu2eII_ana_GenpHist_t_hh
#define __Mu2eII_ana_GenpHist_t_hh

namespace Mu2eII {

  struct GenpHist_t {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fCosTh;			// 
  };
					// histograms for the simulated CE
}
#endif
