///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Mu2eII_ana_TEmuAnaModule_hh
#define Mu2eII_ana_TEmuAnaModule_hh

#include "Mu2eII/ana/TAnaModule.hh"

#include "Mu2eII/ana/TPidMva.hh"
#include "Mu2eII/ana/mva_data.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

namespace Mu2eII {
class TEmuAnaModule: public Mu2eII::TAnaModule {
public:
					//  parameters
  // struct DTrackHist_t : public HistBase_t {
  //   TH1F*    fDp;			// difference between the two track momenta
  //   TH1F*    fRMomErr10;
  //   TH1F*    fDT0;                      // shift in T0 
  // };
//-----------------------------------------------------------------------------
//  fTrackHist[  0]: all tracks
//  fTrackHist[100]: Set C tracks
//-----------------------------------------------------------------------------
  enum { kELE = 0, kMUO = 1 };

  enum { kNTrackPar          =    10 };

  enum { kNEventHistSets     =   100 };
  enum { kNClusterHistSets   =   100 };
  enum { kNTrackSeedHistSets =   100 };
  enum { kNTrackHistSets     = 10000 };
  enum { kNDTrackHistSets    =   100 };
  enum { kNSimpHistSets      =   100 };

  struct Hist_t {
    ClusterHist_t*   fCluster  [kNClusterHistSets];
    EventHist_t*     fEvent    [kNEventHistSets];
    TrackHist_t*     fTrack    [kNTrackHistSets];
    //    TrackSeedHist_t* fTrackSeed[kNTrackSeedHistSets];
    //    DTrackHist_t*    fDTrack   [kNDTrackHistSets];
  };


  struct Cut_t {
    double fXMin;
    double fXMax;
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*        fTrackBlock[2];	// [1]: DeDar [1]:DmuDar
  TStnClusterBlock*      fClusterBlock;
  TGenpBlock*            fGenpBlock;
  TSimpBlock*            fSimpBlock;
  TStepPointMCBlock*     fSpmcBlockVDet;

  TString                fTrackBlockName[2];

  TStnTrackID*           fTrackID[10];            // just in case, need only one
					
  TrackPar_t             fTrackPar[2][kNTrackPar];	// additional track parameters (assume ntracks < 10)
  SimPar_t               fSimPar;		        // additional parameters of the simulated MC particle
  Hist_t                 fHist;		        // histograms filled

  int                    fNID;
  int                    fBestID[2];

  int                    fNHelices;         // 
  int                    fNTrackSeeds;	// N reconstructed track seeds in the event
  int                    fNGoodTrackSeeds;  // N track seeds passing quality cuts
  int                    fNTracks    [2];	// 0:TrkPatRec 1:CalPatRec
  int                    fNGoodTracks[2];
  int                    fNSimp;		// N(simulated particles)

  TStnCluster*           fCluster;
  int                    fNClusters;
//-----------------------------------------------------------------------------
//  MVA
//-----------------------------------------------------------------------------
  int                    fWritePidMvaTree   ;       // =0:PAR   1:DAR  (default=-1)

  TFile*                 fPidMvaFile;
  TTree*                 fPidMvaTree;

  TPidMvaTrainData_t     fPidMvaData;
  TPidMvaTrainBranches_t fPidMvaBranch;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TEmuAnaModule(const char* name="Mu2eII_EmuAna", const char* title="SU2020 EmuAna");
  ~TEmuAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist      () { return &fHist;      }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    SetTrackBlockName(int I, const char* Name) { fTrackBlockName[I] = Name; }

  void    SetWriteMvaTree (int Flag) { fWritePidMvaTree = Flag; }
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

  ClassDef(Mu2eII::TEmuAnaModule,0)
};
}
#endif
