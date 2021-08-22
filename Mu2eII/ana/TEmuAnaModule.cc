//////////////////////////////////////////////////////////////////////////////
// 2020-12-26 P.Murat
// DAR fits only, compare downstream-only reco for electron and muon events 
//
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : events with E/P > 0.7 - learn how to reject mu-->e decays
//  4  : events with MVA PID > 0.5 - learn about the PID filures
//
// call: "temu_ana(3)
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/base/TStnDataset.hh"
#include "Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

#include <string>
#include "math.h"

#include "Mu2eII/ana/TEmuAnaModule.hh"
ClassImp(Mu2eII::TEmuAnaModule)

namespace Mu2eII {
//-----------------------------------------------------------------------------
TEmuAnaModule::TEmuAnaModule(const char* name, const char* title): TAnaModule(name,title) {

  fTrackBlockName[0] = "TrackBlockDarDe";
  fTrackBlockName[1] = "TrackBlockDarDmu";
//-----------------------------------------------------------------------------
// TrackID[0]     : Michael's cuts
// fTrackID[1-19] : cut on MVA-based TrkQual with the step of 0.05
//-----------------------------------------------------------------------------
  fBestID[0] = 0;                       // default best ID word for electrons
  fBestID[1] = 0;                       // default best ID word for muons

  //  fNPidMVA   = 1; 			// PID MVA, just one

  fNID       = 1;                      // fNID has to be < TAnaModule::kMaxNTrackID

  for (int i=0; i<fNID; i++) {
    fTrackID[i] = new TStnTrackID();

    fTrackID[i]->SetMinTanDip (  0.5);
    fTrackID[i]->SetMaxTanDip (  1.0);
    fTrackID[i]->SetMinD0     (-100.);
    fTrackID[i]->SetMaxD0     ( 100.);
    fTrackID[i]->SetMinTrkQual(  0.8);   // assuming I know what I'm doind

    int mask = TStnTrackID::kT0Bit | TStnTrackID::kTanDipBit | TStnTrackID::kD0Bit | TStnTrackID::kTrkQualBit  ;

    fTrackID[0]->SetUseMask(mask);
  }

  SetPidMVA(  "ele00s61b0",1000);

  fWritePidMvaTree = 0;
}

//-----------------------------------------------------------------------------
TEmuAnaModule::~TEmuAnaModule() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TEmuAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName[0].Data(), "TStnTrackBlock"   , &fTrackBlock[0] );
  RegisterDataBlock(fTrackBlockName[1].Data(), "TStnTrackBlock"   , &fTrackBlock[1] );
  RegisterDataBlock("ClusterBlock"           , "TStnClusterBlock" , &fClusterBlock  );
  RegisterDataBlock("SimpBlock"              , "TSimpBlock"       , &fSimpBlock     );
  RegisterDataBlock("GenpBlock"              , "TGenpBlock"       , &fGenpBlock     );
  RegisterDataBlock("SpmcBlockVDet"          , "TStepPointMCBlock", &fSpmcBlockVDet );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// init track ID pointers in TrackPar - do it just once
//  .. assume TrkQual to be precalculated, do only beed to cut on
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    for (int id=0; id<fNID; id++) {
      for (int ip=0; ip<kNTrackPar; ip++) {
	TrackPar_t* tp = &fTrackPar[i][ip];
	tp->fFitType     = 1;                           // both blocks use DAR fits, kDAR=1
	tp->fTrqMvaIndex = 0;                           // both blocks use the same fits, same track quality MVA
	tp->fTrackID[id] = fTrackID[id];
	if (fUseTrqMVA) {
//-----------------------------------------------------------------------------
// in case of on-the-fly calculation store result in TStnTrack::fTmp[0] 
// but do not change the cut value on the MVA output
// in principle, different track blocks could use different track ID definitions
//-----------------------------------------------------------------------------
	  if (fTrqMVA[i]) {
	    tp->fTrackID[id]->SetLocTrkQual(0);
	  }
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// fWriteTmvaTree =  -1 : do not write
//                >=  0 : write PID MVA training tree
//-----------------------------------------------------------------------------
  if (fWritePidMvaTree != 0) {
    TDirectory* dir = gDirectory;

    const char* dsname = GetAna()->GetInputModule()->GetDataset(0)->GetName();

    fPidMvaFile  = new TFile(Form("%s.tmva_training_%04i.root",dsname,1000*fWritePidMvaTree),"recreate");
    fPidMvaTree  = new TTree("tmva_training_tree","TMVA Training Tree");

    fPidMvaBranch.fP          = fPidMvaTree->Branch("p"       ,&fPidMvaData.fP         ,"F");
    fPidMvaBranch.fEcl        = fPidMvaTree->Branch("ecl"     ,&fPidMvaData.fEcl       ,"F");
    fPidMvaBranch.fNCrystals  = fPidMvaTree->Branch("ncr"     ,&fPidMvaData.fNCrystals ,"F");
    fPidMvaBranch.fSeedFr     = fPidMvaTree->Branch("seedfr"  ,&fPidMvaData.fSeedFr    ,"F");

    fPidMvaBranch.fEleTchDt   = fPidMvaTree->Branch("ele_dt"  ,&fPidMvaData.fEleTchDt  ,"F");
    fPidMvaBranch.fEleTchDz   = fPidMvaTree->Branch("ele_dz"  ,&fPidMvaData.fEleTchDz  ,"F");
    fPidMvaBranch.fEleTchDr   = fPidMvaTree->Branch("ele_dr"  ,&fPidMvaData.fEleTchDr  ,"F");
    fPidMvaBranch.fElePath    = fPidMvaTree->Branch("ele_path",&fPidMvaData.fElePath   ,"F");

    fPidMvaBranch.fMuoTchDt   = fPidMvaTree->Branch("muo_dt"  ,&fPidMvaData.fMuoTchDt  ,"F");
    fPidMvaBranch.fMuoTchDz   = fPidMvaTree->Branch("muo_dz"  ,&fPidMvaData.fMuoTchDz  ,"F");
    fPidMvaBranch.fMuoTchDr   = fPidMvaTree->Branch("muo_dr"  ,&fPidMvaData.fMuoTchDr  ,"F");
    fPidMvaBranch.fMuoPath    = fPidMvaTree->Branch("muo_path",&fPidMvaData.fMuoPath   ,"F");

    dir->cd();
  }

  TAnaModule::BeginJob();

  return 0;
}


//_____________________________________________________________________________
void TEmuAnaModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();

  TH1::SetDefaultSumw2(kTRUE);	

  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book cluster histograms
//-----------------------------------------------------------------------------
  int book_cluster_histset[kNClusterHistSets];
  for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

  book_cluster_histset[ 0] = 1;		// all clusters
  book_cluster_histset[ 1] = 1;		// first disk
  book_cluster_histset[ 2] = 1;         // second disk

  for (int i=0; i<kNClusterHistSets; i++) {
    if (book_cluster_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCluster[i] = new ClusterHist_t;
      BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  TString*  track_selection   [kNTrackHistSets];

  for (int i=0; i<kNTrackHistSets; i++) { track_selection[i] = NULL; }

  track_selection[  0] = new TString("good [0] tracks");
  track_selection[  1] = new TString("good [1] tracks in events with no good [0] tracks TrkQual>0.3");
  track_selection[  2] = new TString("good [1] tracks in events with no good [0] tracks TrkQual>0.2");

  track_selection[100] = new TString("[0] all tracks");
  track_selection[101] = new TString("[0] BestTrackID");

  track_selection[102] = new TString("[0] BestTrackID+cluster");
  track_selection[103] = new TString("[0] BestTrackID+cluster, first disk");
  track_selection[104] = new TString("[0] BestTrackID+cluster, second disk");

  track_selection[105] = new TString("[0] BestTrackID+cluster+(p>100)");
  track_selection[106] = new TString("[0] BestTrackID+cluster+(p>100)+(|dr|<100)");
  track_selection[107] = new TString("[0] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)" );
  track_selection[108] = new TString("[0] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)");
  track_selection[109] = new TString("[0] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)");
  track_selection[110] = new TString("[0] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(pid<0.5)");
  track_selection[111] = new TString("[0] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(pid<0.5)");

  track_selection[115] = new TString("[0] BestTrackID+cluster+(87.5<p<92.5)");
  track_selection[116] = new TString("[0] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)");
  track_selection[117] = new TString("[0] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)" );
  track_selection[118] = new TString("[0] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)+(-50<dz<250)");
  track_selection[119] = new TString("[0] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)");
  track_selection[120] = new TString("[0] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(pid<0.5)");
  track_selection[121] = new TString("[0] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(pid<0.5)");
							    
  track_selection[200] = new TString("[1] all tracks");
  track_selection[201] = new TString("[1] BestTrackID");

  track_selection[202] = new TString("[1] BestTrackID, cluster");
  track_selection[203] = new TString("[1] BestTrackID, cluster, first disk");
  track_selection[204] = new TString("[1] BestTrackID, cluster, second disk");

  track_selection[205] = new TString("[1] BestTrackID+cluster+(p>100)");
  track_selection[206] = new TString("[1] BestTrackID+cluster+(p>100)+(|dr|<100)");
  track_selection[207] = new TString("[1] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)" );
  track_selection[208] = new TString("[1] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)");
  track_selection[209] = new TString("[1] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)");
  track_selection[210] = new TString("[1] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(pid<0.5)");
  track_selection[211] = new TString("[1] BestTrackID+cluster+(p>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(pid<0.5)");

  track_selection[215] = new TString("[1] BestTrackID+cluster+(87.5<p<92.5)");
  track_selection[216] = new TString("[1] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)");
  track_selection[217] = new TString("[1] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)" );
  track_selection[218] = new TString("[1] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)+(-50<dz<250)");
  track_selection[219] = new TString("[1] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)");
  track_selection[220] = new TString("[1] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(pid<0.5)");
  track_selection[221] = new TString("[1] BestTrackID+cluster+(87.5<p<92.5)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(pid<0.5)");
							    
  const char* folder_title;
  for (int i=0; i<kNTrackHistSets; i++) {
    if (track_selection[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title = folder_name;
      folder_title = track_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  }

}

//_____________________________________________________________________________

void TEmuAnaModule::FillHistograms() {

  TStnTrack*   trk;
  TrackPar_t*  tp;
  int          ihist;

  int bid_par = fBestID[kELE];
  int bid_dar = fBestID[kMUO];
//-----------------------------------------------------------------------------
// event histograms
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);
//-----------------------------------------------------------------------------
// what does kMUO add ?
// TRK_0 : kELE tracks, BEST_ID
//-----------------------------------------------------------------------------
  if ((fTrackBlock[kELE]->NTracks() > 0) && (fTrackPar[kELE][0].fIDWord[bid_par] == 0)) {
    TStnTrack* trk = fTrackBlock[kELE]->Track(0);
    TrackPar_t* tp = &fTrackPar[kELE][0];
    FillTrackHistograms(fHist.fTrack[0],trk,tp,&fSimPar);
  }
  else {
    if (fTrackBlock[kMUO]->NTracks() > 0) {
      if (fTrackPar[kMUO][0].fIDWord[bid_dar] == 0) {
//-----------------------------------------------------------------------------
// TRK_1: have KDAR track TrkQual > 0.3 and no KPAR TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[kMUO]->Track(0);
	TrackPar_t* tp = &fTrackPar[kMUO][0];
	FillTrackHistograms(fHist.fTrack[1],trk,tp,&fSimPar);
      }
      if (fTrackPar[1][0].fIDWord[2] == 0) {
//-----------------------------------------------------------------------------
// TRK_2: have KDAR track TrkQual > 0.2 and no KPAR TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[1]->Track(0);
	TrackPar_t* tp = &fTrackPar[1][0];
	FillTrackHistograms(fHist.fTrack[2],trk,tp,&fSimPar);
      }
    }
  }
//-----------------------------------------------------------------------------
// KPAR and KDAR histograms, inclusive, ihist defines the offset
// i=0:KPAR, i=1:KDAR
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    //n_setc_tracks[i] = 0;
    ihist            = 100*(i+1);
    int best_id      = fBestID[i];

    for (int itrk=0; itrk<fNTracks[i]; itrk++) {
      trk = fTrackBlock[i]->Track(itrk);
      tp  = fTrackPar  [i]+itrk;
//-----------------------------------------------------------------------------
// TRK_100, TRK_200: all reconstructed tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[ihist+0],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// IHIST+ 1: BestID 
// IHIST+ 2: BestID+cluster
// IHIST+ 3: BestID+cluster+(1st disk)
// IHIST+ 4: BestID+cluster+(2nd disk)
// IHIST+ 5: BestID+(P>100)
// IHIST+ 6: BestID+(P>100)+(|dr|<100)
// IHIST+ 7: BestID+(P>100)+(|dr|<100)+(|dt|<10)
// IHIST+ 8: BestID+(P>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)
// IHIST+ 9: BestID+(P>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)
// IHIST+10: BestID+(P>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(PID>0.5)
// IHIST+11: BestID+(P>100)+(|dr|<100)+(|dt|<10)+(-50<dz<250)+(E/P<1.2)+(PID<0.5)
//-----------------------------------------------------------------------------
      if (tp->fIDWord[best_id] == 0) {
	FillTrackHistograms(fHist.fTrack[ihist+1],trk,tp,&fSimPar);

	if (tp->fDiskID >= 0) {
	  FillTrackHistograms(fHist.fTrack[ihist+2],trk,tp,&fSimPar);
	  if      (tp->fDiskID == 0) FillTrackHistograms(fHist.fTrack[ihist+3],trk,tp,&fSimPar);
	  else if (tp->fDiskID == 1) FillTrackHistograms(fHist.fTrack[ihist+4],trk,tp,&fSimPar);
	}
	
	if (tp->fP > 100) {
	  FillTrackHistograms(fHist.fTrack[ihist+5],trk,tp,&fSimPar);
	  if (fabs(tp->fTchDr) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+6],trk,tp,&fSimPar);
	    if (fabs(tp->fTchDt) < 10) {
	      FillTrackHistograms(fHist.fTrack[ihist+7],trk,tp,&fSimPar);
	      if ((tp->fTchDz > -50) && (tp->fTchDz < 250)) {
		FillTrackHistograms(fHist.fTrack[ihist+8],trk,tp,&fSimPar);
		if (tp->fEp < 1.2) {
		  FillTrackHistograms(fHist.fTrack[ihist+9],trk,tp,&fSimPar);
		  if (tp->fPidMvaOut[0] > 0.5) FillTrackHistograms(fHist.fTrack[ihist+10],trk,tp,&fSimPar);
		  else                         FillTrackHistograms(fHist.fTrack[ihist+11],trk,tp,&fSimPar);
		}
	      }
	    }
	  }
	}

	if ((tp->fP > 87.5) && ( tp->fP < 92.5)) {
	  FillTrackHistograms(fHist.fTrack[ihist+15],trk,tp,&fSimPar);
	  if (fabs(tp->fTchDr) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+16],trk,tp,&fSimPar);
	    if (fabs(tp->fTchDt) < 10) {
	      FillTrackHistograms(fHist.fTrack[ihist+17],trk,tp,&fSimPar);
	      if ((tp->fTchDz > -50) && (tp->fTchDz < 250)) {
		FillTrackHistograms(fHist.fTrack[ihist+18],trk,tp,&fSimPar);
		if (tp->fEp < 1.2) {
		  FillTrackHistograms(fHist.fTrack[ihist+19],trk,tp,&fSimPar);
		  if (tp->fPidMvaOut[0] > 0.5) FillTrackHistograms(fHist.fTrack[ihist+20],trk,tp,&fSimPar);
		  else                         FillTrackHistograms(fHist.fTrack[ihist+21],trk,tp,&fSimPar);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TEmuAnaModule::Event(int ientry) {

  fTrackBlock[0]->GetEntry(ientry);
  fTrackBlock[1]->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSpmcBlockVDet->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to
// be changed
//-----------------------------------------------------------------------------
  fEvtPar.fDioLOWt          = 1.;
  fEvtPar.fDioLLWt          = 1.;
  fEvtPar.fNCrvClusters     = -1;
  fEvtPar.fNCrvPulses       = -1;
  fEvtPar.fNCrvCoincidences = -1;
  fEvtPar.fNGenp            = fGenpBlock->NParticles();
  fEvtPar.fParticle         = NULL;

  fNSimp       = fSimpBlock->NParticles();
  fNClusters   = fClusterBlock->NClusters();

  fCluster = NULL;
  // fEClMax  = -1;
  // fTClMax  = -1;
  if (fNClusters > 0) {
    fCluster = fClusterBlock->Cluster(0);
    // fEClMax  = fCluster->Energy();
    // fTClMax  = fCluster->Time  ();
  }
//-----------------------------------------------------------------------------
// MC generator info
//-----------------------------------------------------------------------------
  for (int i=fEvtPar.fNGenp-1; i>=0; i--) {
    TGenParticle* genp = fGenpBlock->Particle(i);
    int pdg_code       = genp->GetPdgCode();
    int process_code   = genp->GetStatusCode();
    if ((abs(pdg_code) == fPDGCode) && (process_code == fMCProcessCode)) {
      fEvtPar.fParticle = genp;
      break;
    }
  }

  //  fWeight  =  1.;
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimPar in future
//-----------------------------------------------------------------------------
  fSimPar.fParticle = fSimpBlock->Particle(0);
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = NULL;
//-----------------------------------------------------------------------------
// virtual detectors - for fSimPar need parameters at the tracker front
//-----------------------------------------------------------------------------
  int nsteps = fSpmcBlockVDet->NStepPoints();
  
  for (int i=0; i<nsteps; i++) {
    TStepPointMC* step = fSpmcBlockVDet->StepPointMC(i);
    if (step->PDGCode() == fSimPar.fParticle->PDGCode()) {
      if ((step->VolumeID() == 13) || (step->VolumeID() == 14)) {
	fSimPar.fTFront = step;
      }
      else if ((step->VolumeID() == 11) || (step->VolumeID() == 12)) {
	fSimPar.fTMid = step;
      }
    }
  }
//-----------------------------------------------------------------------------
// initialize additional track parameters
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    fNTracks    [i] = fTrackBlock[i]->NTracks();
    fNGoodTracks[i] = 0;

    for (int it=0; it<fNTracks[i]; it++) {
      TrackPar_t* tp = &fTrackPar[i][it];
      tp->fDioLOWt   = fEvtPar.fDioLOWt;
      tp->fDioLLWt   = fEvtPar.fDioLLWt;
    }

    InitTrackPar(fTrackBlock[i],fClusterBlock,fTrackPar[i],&fSimPar);
  }
//-----------------------------------------------------------------------------
// fill histograms
//-----------------------------------------------------------------------------
  FillHistograms();
//-----------------------------------------------------------------------------
// when writing an ntuple for PID MVA training, use the first track
// fWriteTmvaTree = algorithm (0 or 1)
// require fits with both mass hypotheses to succeed and both tracks to have clusters
// also reject muon decays in flight : zmax > 10000 should do that
//-----------------------------------------------------------------------------
  if ((fWritePidMvaTree != 0) && (fNTracks[0] == 1)) {

    TrackPar_t* tpe = &fTrackPar[0][0];   // electron reco, first track
    TrackPar_t* tpm = &fTrackPar[1][0];   // muon reco    , first track

    float zmax = fSimPar.fParticle->EndPos()->Z();

    if ((zmax > 10000) && (tpe->fIDWord[fBestID[0]] == 0) && (tpe->fEcl > 0)) {

      fPidMvaData.fP          = tpe->fP;
      fPidMvaData.fEcl        = tpe->fEcl;
      fPidMvaData.fNCrystals  = tpe->fCluster->NCrystals();
      fPidMvaData.fSeedFr     = tpe->fCluster->SeedFr();
    
      fPidMvaData.fEleTchDt   = tpe->fTchDt;
      fPidMvaData.fEleTchDz   = tpe->fTchDz;
      fPidMvaData.fEleTchDr   = tpe->fTchDr;
      fPidMvaData.fElePath    = tpe->fPath;
    
      if (fNTracks[1] > 0) {
	fPidMvaData.fMuoTchDt   = tpm->fTchDt;
	fPidMvaData.fMuoTchDz   = tpm->fTchDz;
	fPidMvaData.fMuoTchDr   = tpm->fTchDr;
	fPidMvaData.fMuoPath    = tpm->fPath;
      }
      else {
	fPidMvaData.fMuoTchDt   = -999.; // tpm->fTchDt;
	fPidMvaData.fMuoTchDz   = -999.; // tpm->fTchDz;
	fPidMvaData.fMuoTchDr   = -999.; // tpm->fTchDr;
	fPidMvaData.fMuoPath    = -999.; // tpm->fPath;
      }
    
      fPidMvaTree->Fill();
    }
  }

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
// looking mostly at the KDAR tracks
//-----------------------------------------------------------------------------
// diagnostics is related to KDAR
//-----------------------------------------------------------------------------
void TEmuAnaModule::Debug() {

  // TStnTrack* trk;
  // TrackPar_t* tp(NULL);
  //  char        text[500];
//-----------------------------------------------------------------------------
// bit 0: All Events
//-----------------------------------------------------------------------------
  if (GetDebugBit(0) == 1) {
    GetHeaderBlock()->Print(Form("TEmuAnaModule :bit000:"));
    printf("KPAR:\n");
    fTrackBlock[kELE]->Print();

    for (int i=0; i<fNTracks[kELE]; i++) { 
      PrintTrack(&fTrackPar[kELE][i],"data");
    }

    printf("KDAR:\n");
    fTrackBlock[kMUO]->Print();

    for (int i=0; i<fNTracks[kMUO]; i++) { 
      PrintTrack(&fTrackPar[kMUO][i],"data");
    }
  }

  if (GetDebugBit(3) == 1) {
//-----------------------------------------------------------------------------
// bit 3: E/P > 0.7
//-----------------------------------------------------------------------------
    TStnTrackBlock* tb = fTrackBlock[kELE];
    int ntrk = tb->NTracks();
    for (int itrk=0; itrk<ntrk; itrk++) {
      TrackPar_t* tp = &fTrackPar[kELE][itrk];
      if ((tp->fIDWord[fBestID[kELE]] == 0) && (tp->fEp > 0.7)) {
	GetHeaderBlock()->Print(Form("TEmuAnaModule bit:003: tp->fEP = %10.3f",tp->fEp));
      }
    }
  }
  if (GetDebugBit(4) == 1) {
//-----------------------------------------------------------------------------
// bit 3: E/P > 0.7
//-----------------------------------------------------------------------------
    TStnTrackBlock* tb = fTrackBlock[kELE];
    int ntrk = tb->NTracks();
    for (int itrk=0; itrk<ntrk; itrk++) {
      TrackPar_t* tp = &fTrackPar[kELE][itrk];
      if ((tp->fIDWord[fBestID[kELE]] == 0) && (tp->fPidMvaOut[0] > 0.5)) {
	GetHeaderBlock()->Print(Form("TEmuAnaModule bit:004: tp->fPidMvaOut[0] = %10.3f",tp->fPidMvaOut[0]));
      }
    }
  }
}

//-----------------------------------------------------------------------------
int TEmuAnaModule::EndJob() {

  if (fWritePidMvaTree != 0) {
    printf("[%s::EndJob] Writing output MVA training file %s\n",GetName(),fPidMvaFile->GetName());

    fPidMvaFile->Write();

    delete fPidMvaFile;
    fPidMvaFile  = 0;
    fPidMvaTree  = 0;
  }
  
  TAnaModule::EndJob();
  return 0;
}

//_____________________________________________________________________________
void TEmuAnaModule::Test001() {
}


}
