///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TCosmicAnaModule_hh
#define Stntuple_ana_TCosmicAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

#include "Mu2eII/ana/SimPar_t.hh"
#include "Mu2eII/ana/TrackPar_t.hh"
#include "Mu2eII/ana/TrackHist_t.hh"
#include "Mu2eII/ana/CrvStubPar_t.hh"

#include "Mu2eII/ana/TAnaModule.hh"

namespace Mu2eII {
class TCosmicAnaModule: public TAnaModule {
public:
  enum { kMaxCrvStubs             =   100 };
  enum { kNTrackPar               =    20 };

  enum { kNEventHistSets          =   100 };
  enum { kNTrackHistSets          = 10000 };
  enum { kNTrackTcHistSets        =   100 };
  enum { kNTrackCrvStHistSets     =   100 };
  enum { kNGenpHistSets           =   100 };
  enum { kNSimpHistSets           =   100 };
  enum { kNCrvClusterHistSets     =   100 };
  enum { kNCrvPulseHistSets       =   100 };
  enum { kNCrvCoincidenceHistSets =   100 };

  struct Hist_t {
    EventHist_t*                   fEvent         [kNEventHistSets];
    Mu2eII::TrackHist_t*           fTrack         [kNTrackHistSets];
    Mu2eII::TrackCrvStHist_t*      fTrackCrvSt    [kNTrackCrvStHistSets];
    Mu2eII::TrackTcHist_t*         fTrackTc       [kNTrackTcHistSets];
    Mu2eII::GenpHist_t*            fGenp          [kNGenpHistSets];
    Mu2eII::SimpHist_t*            fSimp          [kNSimpHistSets];
    Mu2eII::CrvClusterHist_t*      fCrvCluster    [kNCrvClusterHistSets];
    Mu2eII::CrvCoincidenceHist_t*  fCrvCoincidence[kNCrvCoincidenceHistSets];
    Mu2eII::CrvPulseHist_t*        fCrvPulse      [kNCrvPulseHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*       fTrackBlockDe;
  TStnTrackBlock*       fTrackBlockUe;

  TStnHelixBlock*       fHelixBlockDe;
  TStnHelixBlock*       fHelixBlockUe;  // need both

  TStnTrackSeedBlock*   fTrackSeedBlockDe;
  TStnTrackSeedBlock*   fTrackSeedBlockUe;

  TStnClusterBlock*     fClusterBlock;
  TGenpBlock*           fGenpBlock;
  TSimpBlock*           fSimpBlock;
  TStepPointMCBlock*    fSpmcBlockVDet;

  TStnTimeClusterBlock* fTimeClusterBlockDe;   // merged DE, helices only
  TStnTimeClusterBlock* fTimeClusterBlockUe;   // merged UE, helices only

  TStnTimeClusterBlock* fTCFinderBlockDe;      // non-merged ones, De
  TStnTimeClusterBlock* fTCFinderBlockUe;      // non-merged ones, Ue
  TStnTimeClusterBlock* fCTPFinderBlock;       // non-merged ones, CalTimePeakFinder, De only

  TCrvClusterBlock*     fCrvClusterBlock;
  TCrvPulseBlock*       fCrvPulseBlock;

  CosmicVetoData_t      fCosmicVetoData;      // holder for non-CRV data blocks used to veto cosmics

					// additional track parameters (assume ntracks < 20)
  Mu2eII::EventPar_t    fEvtPar;        // defined in TAnaModule.hh
  Mu2eII::CrvStubPar_t  fCRVStubPar[kMaxCrvStubs];
  Mu2eII::TrackPar_t    fTrackParDe[kNTrackPar];
  Mu2eII::TrackPar_t    fTrackParUe[kNTrackPar];
  Mu2eII::SimPar_t      fSimPar;		// additional parameters of the simulated MC particle
					// histograms filled
  Hist_t                fHist;
					// cut values
  double                fPtMin;

  TGenParticle*         fParticle;		// electron or muon

  TSimParticle*         fSimp;
  double                fEleE;		// electron energy

  int                   fNGoodTracks;
  int                   fNGenp;		// N(generated particles)

  TString               fTrackBlockNameDe;
  TString               fTrackBlockNameUe;

  int                   fUseAllPulses;

  int                   fNGoodTracks_BOX;
  int                   fNGoodTracks_MVA;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TCosmicAnaModule(const char* name="Mu2eII_CosmicAna", const char* title="CosmicAna");
  ~TCosmicAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnTrackBlock*    GetTrackBlockDe() { return fTrackBlockDe; }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void               SetTrackBlockName(int Flag, const char* Name) { 
    if      (Flag == 0) fTrackBlockNameDe = Name ; 
    else if (Flag == 1) fTrackBlockNameUe = Name ; 
  }

  void               SetUseAllPulses(int Flag) { fUseAllPulses = Flag; }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookTrkCrvHistograms();
  void    BookHistograms();

  void    FillrkCrvHistograms();
  void    FillHistograms();

  void    Debug();
  void    PrintTrack(TStnTrack* Track, Mu2eII::TrackPar_t* Tp, Option_t* Option) const ;
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

};
}
#endif
