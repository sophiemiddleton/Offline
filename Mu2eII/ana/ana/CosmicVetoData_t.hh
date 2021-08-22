///////////////////////////////////////////////////////////////////////////////
// structure to facilitate passing multiple data blocks with the event data around
///////////////////////////////////////////////////////////////////////////////

namespace Mu2eII {

  struct CosmicVetoData_t {
    TStnTimeClusterBlock* fTCFinderBlockUe;  // dont immediately need the De one 
    TStnHelixBlock*       fHelixBlockDe;
    TStnHelixBlock*       fHelixBlockUe;
    TStnTrackBlock*       fTrackBlockDe;
    TrackPar_t*           fTrackParDe;
    TStnTrackBlock*       fTrackBlockUe;
    TStnClusterBlock*     fClusterBlock;
  };

};
