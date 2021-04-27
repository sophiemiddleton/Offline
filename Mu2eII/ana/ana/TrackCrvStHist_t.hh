#ifndef __Mu2eII_ana_TrackCrvStHist_t_hh
#define __Mu2eII_ana_TrackCrvStHist_t_hh

namespace Mu2eII {

  struct TrackCrvStHist_t {
    TH1F*    fDt;			// track-TC time difference 
    TH1F*    fDt2;			// track-TC time difference 
    TH2F*    fDtVsZCrv;
    TH2F*    fDt2VsZCrv;
  };
}
#endif
