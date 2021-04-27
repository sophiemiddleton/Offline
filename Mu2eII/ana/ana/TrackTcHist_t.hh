#ifndef __Mu2eII_ana_TrackTcHist_t_hh
#define __Mu2eII_ana_TrackTcHist_t_hh

namespace Mu2eII {

  struct TrackTcHist_t {
    TH1F*    fDt;			// track-TC time difference 
    TH1F*    fClusterZ;
    TH1F*    fUeDt;
    TH1F*    fUeClusterZ;
    TH1F*    fDeUeDt;
  };
}
#endif
