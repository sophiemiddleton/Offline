#ifndef __Mu2eII_ana_TrackTcPar_t__
#define __Mu2eII_ana_TrackTcPar_t__

#include "Stntuple/ana/ParBase_t.hh"

class TStnTrackID;

namespace Mu2eII {
  class TrackTcPar_t : public ParBase_t {
  public:
    float fDt; 
    float fClusterZ; 
    float fUeDt; 
    float fUeClusterZ; 
    float fDeUeDt; 
  };
}
#endif
