///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TCrvAnaModule_hh
#define Stntuple_ana_TCrvAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Mu2eII/ana/SimPar_t.hh"
#include "Mu2eII/ana/CrvStubPar_t.hh"

#include "Mu2eII/ana/TAnaModule.hh"

namespace Mu2eII {
class TCrvAnaModule: public TAnaModule {
public:
  enum { kMaxCrvStubs             =   100 };

  enum { kNEventHistSets          =   100 };
  enum { kNGenpHistSets           =   100 };
  enum { kNSimpHistSets           =   100 };
  enum { kNCrvClusterHistSets     =   200 };
  enum { kNCrvPulseHistSets       =   200 };
  enum { kNCrvCoincidenceHistSets =   200 };

  struct Hist_t {
    EventHist_t*                   fEvent         [kNEventHistSets];
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
  TGenpBlock*           fGenpBlock;
  TSimpBlock*           fSimpBlock;
  TStepPointMCBlock*    fSpmcBlockVDet;
  TCrvClusterBlock*     fCrvClusterBlock;
  TCrvPulseBlock*       fCrvPulseBlock;

					// additional track parameters (assume ntracks < 20)
  Mu2eII::EventPar_t    fEvtPar;        // defined in TAnaModule.hh
  Mu2eII::CrvStubPar_t  fCrvStubPar[kMaxCrvStubs];
  Mu2eII::SimPar_t      fSimPar;		// additional parameters of the simulated MC particle
					// histograms filled
  Hist_t                fHist;
					// cut values
  TGenParticle*         fParticle;		// electron or muon

  TSimParticle*         fSimp;
  double                fEleE;		// electron energy

  int                   fNGenp;		// N(generated particles)

  int                   fUseAllPulses;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TCrvAnaModule(const char* name="Mu2eII_CrvAna", const char* title="CrvAna");
  ~TCrvAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void     SetUseAllPulses(int Flag) { fUseAllPulses = Flag; }
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
