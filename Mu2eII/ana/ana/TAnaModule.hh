///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __Mu2eII_ana_TAnaModule_hh__
#define __Mu2eII_ana_TAnaModule_hh__

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"

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
#include "Stntuple/alg/TStntuple.hh"

#include "Mu2eII/ana/SimPar_t.hh"
#include "Mu2eII/ana/TrackPar_t.hh"
#include "Mu2eII/ana/TrackTcPar_t.hh"
#include "Mu2eII/ana/EventPar_t.hh"

#include "Mu2eII/ana/SimpHist_t.hh"
#include "Mu2eII/ana/CrvHist_t.hh"
#include "Mu2eII/ana/ClusterHist_t.hh"
#include "Mu2eII/ana/TrackCrvStHist_t.hh"
#include "Mu2eII/ana/TrackHist_t.hh"
#include "Mu2eII/ana/TrackTcHist_t.hh"
#include "Mu2eII/ana/EventHist_t.hh"
#include "Mu2eII/ana/GenpHist_t.hh"
#include "Mu2eII/ana/CrvStubPar_t.hh"

#include "Mu2eII/ana/mva_data.hh"

#include "Mu2eII/ana/CosmicVetoData_t.hh"

namespace Mu2eII {
class TAnaModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
  int                 fNTrkID;              // counts fTrackID_BOX + fTrackID_MVA 
  TStnTrackID*        fTrackID_BOX;
  TStnTrackID*        fTrackID_MVA;
  TStnTrackID*        fTrackID_EP_BOX;      // specifically mu- --> e+

  double              fMbTime;              // microbunch time
  double              fMinT0;
  int                 fApplyCorr;           // default=1

  int                 fPDGCode;             // 
  int                 fMCProcessCode;       // 
  double              fEventWeight;         // weight applied when filling histograms
  int                 fBatchMode;	    // use for reweighting
  double              fEleE;                // relevant event generated energy - to migrate to fEvtPar **FIXME**

  Mu2eII::EventPar_t  fEvtPar;
  TStntuple*          fStnt;                   // STNTUPLE singleton
//-----------------------------------------------------------------------------
// MVA-based track quality (TRQ) and particle ID (PID) classifiers
//-----------------------------------------------------------------------------
  int                 fUseTrqMVA;
  mva_data*           fTrqMVA[2];              // TRQ 

  int                 fUsePidMVA;
  mva_data*           fPidMVA;                 // PID

  int                 fDebugLevel;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TAnaModule(const char* name="Mu2eII_Ana", const char* title="Ana");
  ~TAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void   SetMinT0        (double T0   ) { fMinT0         = T0  ;  }
  void   SetPDGCode      (int    Code ) { fPDGCode       = Code;  }
  void   SetMCProcessCode(int    Code ) { fMCProcessCode = Code;  }
  void   SetApplyCorr    (int    Flag ) { fApplyCorr     = Flag;  }
  void   SetDebugLevel   (int    Level) { fDebugLevel    = Level; }
//-----------------------------------------------------------------------------
// TRQ MVA Training Codes: 
//
// 0060 : PAR dPf > 0.60
// 0070 : PAR dPf > 0.70
// 1060 : DAR dPf > 0.60
// 1070 : DAR dPf > 0.70
//-----------------------------------------------------------------------------
  void                SetTrqMVA      (const char* Dataset, int MvaTrainingCode);
//-----------------------------------------------------------------------------
// PID MVA
//-----------------------------------------------------------------------------
  void                SetPidMVA      (const char* Dataset, int MvaTrainingCode);
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int     BeginJob();
  virtual int     BeginRun();
  // virtual int     Event   (int ientry);
  //  virtual int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  double  BatchModeWeight(float lumi, int mode);

  void    BookClusterHistograms   (Mu2eII::ClusterHist_t*    Hist, const char* Folder);
  void    BookCrvClusterHistograms(Mu2eII::CrvClusterHist_t* Hist, const char* Folder);
  void    BookCrvPulseHistograms  (Mu2eII::CrvPulseHist_t*   Hist, const char* Folder);
  void    BookGenpHistograms      (Mu2eII::GenpHist_t*       Hist, const char* Folder);
  void    BookEventHistograms     (Mu2eII::EventHist_t*      Hist, const char* Folder);
  void    BookSimpHistograms      (Mu2eII::SimpHist_t*       Hist, const char* Folder);
  void    BookTrackHistograms     (Mu2eII::TrackHist_t*      Hist, const char* Folder);
  void    BookTrackIDHistograms   (TStnTrackID::Hist_t*      Hist, const char* Folder);
  void    BookTrackTcHistograms   (Mu2eII::TrackTcHist_t*    Hist, const char* Folder);
  void    BookTrackCrvStHistograms(Mu2eII::TrackCrvStHist_t* Hist, const char* Folder);

  void    FillCrvClusterHistograms(Mu2eII::CrvClusterHist_t*  Hist, TCrvCoincidenceCluster* CrvCl,
				   Mu2eII::CrvStubPar_t*      CrvPar);
  void    FillCrvPulseHistograms  (Mu2eII::CrvPulseHist_t*    Hist, TCrvRecoPulse* CrvCl);
  void    FillClusterHistograms   (Mu2eII::ClusterHist_t* Hist, TStnCluster*  Cluster, double Weight = 1.);

  void    FillEventHistograms     (Mu2eII::EventHist_t*  Hist, Mu2eII::EventPar_t*  Evtpar  );

  void    FillGenpHistograms      (Mu2eII::GenpHist_t*   Hist, TGenParticle* Genp);
  void    FillSimpHistograms      (Mu2eII::SimpHist_t*   Hist, TSimParticle* Simp);

  void    FillTrackHistograms     (Mu2eII::TrackHist_t* Hist, 
				   TStnTrack*           Trk, 
				   Mu2eII::TrackPar_t*  Tp, 
				   Mu2eII::SimPar_t*    SimPar,
				   double               Weight = 1.);

  void    FillTrackTcHistograms   (Mu2eII::TrackTcHist_t* Hist,
				   TrackTcPar_t*          TTc);

  void    FillTrackCrvStHistograms(Mu2eII::TrackCrvStHist_t* Hist    ,
				   TrackPar_t*               TrackPar,
				   CrvStubPar_t*             CrvStPar);

  int     InitCrvStubPar(TCrvClusterBlock*      ClusterBlock, 
			 Mu2eII::CrvStubPar_t*  CrvStubPar  , 
			 TSimParticle*          SimPar=0);

  int     InitTrackPar(TStnTrackBlock*      TrackBlock  , 
		       TStnClusterBlock*    ClusterBlock, 
		       Mu2eII::TrackPar_t*  TrackPar    ,
		       Mu2eII::SimPar_t*    SimPar      );

					// veto cosmics based on the tracker+caloriemter

  int     NonCrvCosmicVeto(Mu2eII::CosmicVetoData_t* Data, Mu2eII::EventPar_t* EvtPar);

  void    PrintTrack(Mu2eII::TrackPar_t* Tp, Option_t* Option) const ;

  //  ClassDef(TAnaModule,0)
};
}
#endif
