#ifndef __Mu2eII_ana_EventPar_t__
#define __Mu2eII_ana_EventPar_t__

namespace Mu2eII {
					// cosmics rejection
  enum {
    kNTrkDeVetoBit     = 0x0001, // N(De tracks) > 1
    kNTrkUeVetoBit     = 0x0002, // N(Ue tracks) > 1
    kTrkUeDtVetoBit    = 0x0004, // Dt(De - Ue) > 50 ns
    kTrkTcVetoBit      = 0x0008, // Dt(De - Ue) > 50 ns
    kNHelDeVetoBit     = 0x0010, // N(De helices) > 1
    kNHelUeVetoBit     = 0x0020, // N(Ue helices) > 1
    kCaloInTimeVetoBit = 0x0040, // too energetic cluster in the calorimeter
    kCaloEarlyVetoBit  = 0x0080, // too energetic cluster in the calorimeter
    kCrvStubVetoBit    = 0x0100  // CRV stub close in time
  };

  struct EventPar_t {
    int           fNGenp;
    TGenParticle* fParticle;                // generator "signal" particle
    float         fGenE;		    // generator (signal) particle energy

    int           fNTracksDe;               // might need two for different reasons: DAR<->PAR, ELE<->MUO
    int           fNTracksUe;               // number of reconstructed UE tracks

    int           fNGoodTracks;
    int           fNStrawHits;

    int           fNHelicesDe;              // number of reconstructed downstream helices
    int           fNHelicesUe;              // number of reconstructed upstream   helices

    int           fNClusters;               // calorimeter clusters

    int           fNCrvClusters;
    int           fNCrvCoincidences;
    int           fNCrvPulses;
    
    int           fNEleCandidates_BOX;      // N(electron candidates), BOX cuts
    int           fNEleCandidates_MVA;      // N(electron candidates), MVA cuts

    float         fTimeClusterDt;           // difference in time between time clusters in events with 2 TCs
    float         fAbsTimeClusterDt;        // abs difference in time between time clusters in events with 2 TCs
    int           fNTimeClusters;           // number of time clusters present in each event
    int           fNEffTimeClusters;        // effective number of time clusters in each event
    int           fNTCIndex[2];             // holds indices of 'good' time clusters, with no duplicates
    int           fCutCounter[4];           // counter to find outlier events

    int           fCosmicVeto;              // overall veto flag: 10*fNHelices + fTimeClusterVeto
    int           fTimeClusterVeto;         // 

    int           fCandidate_BOX;           // if 1, passes analysis cuts = event candidate
    int           fCandidate_MVA;           // if 1, passes analysis cuts = event candidate

    int           fVetoedCandidate_BOX;     // if 1, passes analysis cuts but fails CRV veto
    int           fVetoedCandidate_MVA;     // if 1, passes analysis cuts but fails CRV veto
    int           fTCType;

    float         fInstLum;
					    // different weights, one per event
    double        fOneBatchWeight;
    double        fTwoBatchWeight;
    double        fLumiWt;                  // legacy reweighting factor: either mode=2: Two/One, mode=1: One/Two

    double        fDioLOWt;
    double        fDioLLWt;

    double        fDtUe;
    double        fDpUe;
    double        fDchiUe;
    double        fUeGate;
    double        fDtUe_goodUe;
    double        fDtUe_goodDe;
    double        fSameLegDp;
    double        fDiffLegDp;
    double        fSameLegDchi2;
    double        fDiffLegDchi2;
    double        fNGoodTracksTotal;
  };
}
#endif
