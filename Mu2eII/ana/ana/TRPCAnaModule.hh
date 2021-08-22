///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TRPCAnaModule_hh
#define Stntuple_ana_TRPCAnaModule_hh

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

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/alg/TStnTrackID.hh"

#include "Mu2eII/ana/SimPar_t.hh"
#include "Mu2eII/ana/TrackPar_t.hh"
#include "Mu2eII/ana/TrackHist_t.hh"

#include "Mu2eII/ana/TAnaModule.hh"

namespace Mu2eII {
class TRPCAnaModule: public TAnaModule {
public:
  enum { kNTrackPar        =    20 };

  enum { kNEventHistSets   =   100 };
  enum { kNTrackHistSets   = 10000 };
  enum { kNGenpHistSets    =   100 };
  enum { kNSimpHistSets    =   100 };

  struct Hist_t {
    EventHist_t*          fEvent     [kNEventHistSets];
    Mu2eII::TrackHist_t*  fTrack     [kNTrackHistSets];
    GenpHist_t*           fGenp      [kNGenpHistSets];
    SimpHist_t*           fSimp      [kNSimpHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*     fTrackBlock;
  TStnClusterBlock*   fClusterBlock;
  TGenpBlock*         fGenpBlock;
  TSimpBlock*         fSimpBlock;
  TStepPointMCBlock*  fSpmcBlockVDet;

					// additional track parameters (assume ntracks < 20)
  Mu2eII::TrackPar_t  fTrackPar[kNTrackPar];
  Mu2eII::SimPar_t    fSimPar;		// additional parameters of the simulated MC particle
					// histograms filled
  Hist_t              fHist;
					// cut values
  double              fPtMin;

  TGenParticle*       fParticle;		// electron or muon

  TSimParticle*       fSimp;
  double              fEleE;		// electron energy

  int                 fNGoodTracks;
  int                 fNMatchedTracks;
  int                 fNGenp;		// N(generated particles)

  int                 fNHyp;
  int                 fBestHyp[10];
  int                 fFillDioHist;
					// fTrackNumber[i]: track number, 
					// corresponding to OBSP particle #i
					// or -1
  TStnArrayI          fTrackNumber;

  TStnTrack*          fTrack;

  double              fTau;            // capture time / pion lifetime
  double              fSurvProb;       // simulated captured pion survival probability 

  int                 fBestID;

  TString             fTrackBlockName;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TRPCAnaModule(const char* name="Mu2eII_RPCAna", const char* title="RPCAna");
  ~TRPCAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnTrackBlock*    GetTrackBlock  () { return fTrackBlock;   }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void               SetTrackBlockName(const char* Name) { fTrackBlockName = Name ; }
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
  void    PrintTrack(TStnTrack* Track, Mu2eII::TrackPar_t* Tp, Option_t* Option) const ;
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

};
}
#endif
