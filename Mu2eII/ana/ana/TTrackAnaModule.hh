///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TTrackAnaModule_hh
#define Stntuple_ana_TTrackAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TString.h"

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

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/alg/TStnTrackID.hh"

#include "Mu2eII/ana/SimPar_t.hh"
#include "Mu2eII/ana/TrackPar_t.hh"
#include "Mu2eII/ana/EventPar_t.hh"

#include "Mu2eII/ana/TrackHist_t.hh"

#include "Mu2eII/ana/TAnaModule.hh"

namespace Mu2eII {

class TTrackAnaModule: public Mu2eII::TAnaModule {
public:

  enum { kNTrackPar        =    20 };

  enum { kNEventHistSets   =   100 };
  enum { kNTrackHistSets   = 10000 };
  enum { kNGenpHistSets    =   100 };
  enum { kNSimpHistSets    =   100 };
  enum { kNTrackIDHistSets =    10 };

  struct Hist_t {
    Mu2eII::EventHist_t*  fEvent     [kNEventHistSets  ];
    Mu2eII::TrackHist_t*  fTrack     [kNTrackHistSets  ];
    Mu2eII::GenpHist_t*   fGenp      [kNGenpHistSets   ];
    Mu2eII::SimpHist_t*   fSimp      [kNSimpHistSets   ];
    TStnTrackID::Hist_t*  fTrackID   [kNTrackIDHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*       fTrackBlockDe;
  TStnTrackBlock*       fTrackBlockUe;

  TStnHelixBlock*       fHelixBlockDe;
  TStnHelixBlock*       fHelixBlockUe;

  TStnTimeClusterBlock* fTCFinderBlockUe; // logically, need only UE to deal with cosmics

  TStnClusterBlock*     fClusterBlock;

  TGenpBlock*           fGenpBlock;
  TSimpBlock*           fSimpBlock;
  TStepPointMCBlock*    fSpmcBlockVDet;

  Mu2eII::CosmicVetoData_t fCosmicVetoData;

					// additional track parameters (assume ntracks < 20)
  Mu2eII::TrackPar_t  fTrackPar[kNTrackPar];
  Mu2eII::SimPar_t    fSimPar;		// additional parameters of the simulated MC particle
					// histograms filled
  Hist_t              fHist;
					// cut values
    //  double              fPtMin;

  TGenParticle*       fParticle;	// electron or muon

  TSimParticle*       fSimp;

  int                 fNGoodTracks;
  int                 fNMatchedTracks;
  int                 fNGenp;		// N(generated particles)

  int                 fNHyp;
  int                 fBestHyp[10];
  int                 fFillDioHist;
					// fTrackNumber[i]: track number, 
					// corresponding to OBSP particle #i
					// or -1
    //  TStnArrayI          fTrackNumber;

  TStnTrack*          fTrack;
    //  int                 fBestID;

  TString             fTrackBlockNameDe;
  TString             fTrackBlockNameUe;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackAnaModule(const char* name="Mu2eII_TrackAna", const char* title="TrackAna");
  ~TTrackAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnTrackBlock*    GetTrackBlockDe() { return fTrackBlockDe; }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void               SetTrackBlockNameDe(const char* Name) { 
    fTrackBlockNameDe = Name ; 
    printf("%s::SetTrackBlockNameDe : track block name set to: %s\n",GetName(),fTrackBlockNameDe.Data());
  }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int     BeginJob();
    //  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookHistograms();
  void    FillHistograms();

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

};
}
#endif
