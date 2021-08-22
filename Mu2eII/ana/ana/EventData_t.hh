///////////////////////////////////////////////////////////////////////////////
// structure to facilitate passing multiple data blocks with the event data around
///////////////////////////////////////////////////////////////////////////////

namespace Mu2eII {

  struct CosmicsVetoData_t {
    TStnTimeClusterBlock* fTCFinderBlockUe;
    TStnHelixBlock*       fHelixBlockDe;
    TStnHelixBlock*       fHelixBlockUe;
    TStnTrackBlock*       fTrackBlockDe;
    TStnTrackBlock*       fTrackBlockUe;
    TStnClusterBlock*     fClusterBlock;
  };

};
