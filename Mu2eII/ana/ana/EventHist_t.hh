#ifndef __Mu2eII_ana_EventHist_t_hh
#define __Mu2eII_ana_EventHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace Mu2eII {

  struct EventHist_t {
    TH1F*    fEventWeight[2];           // MC truth
    TH1F*    fEventE;                   // MC truth relevant event energy
    TH1F*    fInstLumi[3];              // MC truth lumi info: 0: nominal 1: undo batch weight 2: opposite batch weight
    TH1F*    fBatchWeight[2];           // MC truth
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;

    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNshCE;

    TH1F*    fMcMom;
    TH1D*    fDioMom;
    TH1F*    fMcCosTh;
    TH1F*    fNHelicesDe;
    TH1F*    fNHelicesUe;
    TH1F*    fNTracksDe;
    TH1F*    fNTracksUe;
    TH1F*    fNShTot [2];
    TH1F*    fNGoodSH;
    TH1F*    fDtClT;
    TH1F*    fDtClS;
    TH1F*    fSHTime;
    TH1F*    fNHyp;
    TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    TH1F*    fNGenp;                    // N(particles in GENP block)

    TH1F*    fNClusters;
    TH1F*    fEClMax;			// energy of the first (highest) reconstructed cluster
    TH1F*    fTClMax;			// time   of the first (highest) reconstructed cluster
    TH1F*    fDp;                       // P(TrkPatRec)-P(CalPatRec)
    TH1F*    fWeight;			// weight, need with statistics
    TH1F*    fGMom;			// photon momentum
    TH1F*    fGMomRMC;                  // photon momentum, RMC weighted

    TH1F*    fNCrvClusters;
    TH1F*    fNCrvCoincidences[2];
    TH1F*    fNCrvPulses[2];
    TH1F*    fTimeClusterDt[11];            // time difference between time clusters in events with 2 time clusters
    TH1F*    fAbsTimeClusterDt;         // absolute value of fTimeClusterDt
    TH1F*    fAbsTimeClusterDt2;        // zoomed absolute value of fTimeClusterDt
    TH1F*    fNTimeClusters;            // number of time clusters
    TH1F*    fNEffTimeClusters;         // effective number of time clusters
    TH1F*    fCrvCutFlow;               // series of event cuts to determine outliers
    TH1F*    fTimeClusterVeto;          // 
    TH1F*    fCosmicVeto;               // overall veto =0 is events passes
    
    TH1F*    fDtUe;
    TH1F*    fDpUe;
    TH1F*    fDchiUe;
    TH1F*    fUeGate;
    TH1F*    fDtUe_goodUe;
    TH1F*    fDtUe_goodDe;
    TH1F*    fSameLegDp;
    TH1F*    fSameLegDchi2;
    TH1F*    fDiffLegDp;
    TH1F*    fDiffLegDchi2;
    TH1F*    fNGoodTracksTotal;
    
  };
}
#endif
