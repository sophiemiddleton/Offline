//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
// 
//  3  : events with set C tracks and 70mm < |dx|  < 90 mm
//  4  : events with DpF > 1 MeV : obviously, misreconstructed ones
//  5  : events with N(tracks) > 1
//  6  : events trk_41 with 0.8< E/P < 1.1 - tracks missed by CalPatRec
//  7  : events (muo) with LogLHRCal >   20
//  8  : events (ele) with LogLHRCal < - 20
//  9  : events (muo) with 0.42 < E/P < 0.46
// 10  : events (muo) with Set C track with ECL > 80 MeV
// 28  : Set C DEM tracks with E/P > 1.1
// 29  : TRK_19 (Set C DEM tracks with a cluster) and LLHR(cal) < 0
// 31  : EVT_6 events with ce_costh > 0.8 
// 32  : TRK_1 events with chi2tcm > 100. 
// 33  : DU < -80mm - study edge effects
// 34  : EVT_7: events with E_CL > 60 and no tracks (makes sense only for single CE events)
// 35  : TRK_1: events with P > 106 MeV/c - misreconstruction
// 36  : TRK_23 events with P < 80: odd misidentified muons - turned out to be DIO electrons
// 37  : TRK_26 LLHR_CAL > 5
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
#include "Mu2eII/ana/TPbarAnaModule.hh"

// ClassImp(TPbarAnaModule)

namespace Mu2eII {

//-----------------------------------------------------------------------------
TPbarAnaModule::TPbarAnaModule(const char* name, const char* title):
  TAnaModule(name,title)
{
  fPtMin  = 1.;
  fTrackNumber.Set(100);
//-----------------------------------------------------------------------------
// MC truth is defined in TAnaModule
//-----------------------------------------------------------------------------
  fTrackBlockName = "TrackBlockDar";
  fBestID         = 1;                   // best ID : MVA, not used here
  fNTrkID         = 2;
}

//-----------------------------------------------------------------------------
TPbarAnaModule::~TPbarAnaModule() {
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TPbarAnaModule::BeginJob() {

  TAnaModule::BeginJob();

//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data(),"TStnTrackBlock"      ,&fTrackBlock     );
  RegisterDataBlock("ClusterBlock"        ,"TStnClusterBlock"    ,&fClusterBlock   );
  RegisterDataBlock("GenpBlock"           ,"TGenpBlock"          ,&fGenpBlock      );
  RegisterDataBlock("SimpBlock"           ,"TSimpBlock"          ,&fSimpBlock      );
  RegisterDataBlock("SpmcBlockVDet"       , "TStepPointMCBlock"  ,&fSpmcBlockVDet );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// init track ID pointers in TrackPar - do it just once
//-----------------------------------------------------------------------------
  for (int i=0; i<kNTrackPar; i++) {
    TrackPar_t* tp = &fTrackPar[i];

    tp->fFitType     = 1;                                // assume DAR, PAR needs to be set specially
    tp->fTrqMvaIndex = 0;                                // index of the TRQ MVA used by this block

    tp->fDioLOWt     = 1.;                               // no DIO here, avoid confusions 
    tp->fDioLLWt     = 1.;
    
    tp->fTrackID[0]  = fTrackID_BOX;
    tp->fTrackID[1]  = fTrackID_MVA;
  }
  return 0;
}

//_____________________________________________________________________________
void TPbarAnaModule::BookHistograms() {

  //  char name [200];
  //  char title[200];
  
  TFolder*    fol;
  TFolder*    hist_folder;
  char        folder_name[200];
  const char* folder_title;
  
  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
  
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  TString*    event_selection [kNEventHistSets];

  for (int i=0; i<kNEventHistSets; i++) event_selection[i] = 0;

  event_selection[ 0] = new TString("all events");
  event_selection[ 1] = new TString("events with a reconstructed track");

  for (int i=0; i<kNEventHistSets; i++) {
    if (event_selection[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title    = event_selection[i]->Data();
      if (! fol) fol  = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book simp histograms
//-----------------------------------------------------------------------------
  TString*  simp_selection [kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) simp_selection[i] = 0;

  simp_selection[ 0] = new TString("all events");

  for (int i=0; i<kNSimpHistSets; i++) {
    if (simp_selection[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title   = simp_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  TString* track_selection[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) track_selection[i] = 0;

  track_selection[   0] = new TString("all tracks");
  track_selection[   1] = new TString("negative tracks");
  track_selection[   2] = new TString("positive tracks");

  track_selection[1000] = new TString("e- tracks passing BOX TRQ + PID");
  track_selection[1004] = new TString("e- tracks passing BOX TRQ + PID batch 1 wt");
  track_selection[1005] = new TString("e- tracks passing BOX TRQ + PID batch 2 wt");
  track_selection[1010] = new TString("e- tracks passing BOX TRQ + PID, T>700");
  track_selection[1014] = new TString("e- tracks passing BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[1015] = new TString("e- tracks passing BOX TRQ + PID, T>700 batch 2 wt");

  track_selection[1500] = new TString("all mu- tracks passing analysis_box cuts");
  track_selection[1501] = new TString("all mu- tracks passing analysis_box cuts, T0>700");

  track_selection[2000] = new TString("e- tracks passing MVA TRQ + PID");
  track_selection[2001] = new TString("e- tracks passing MVA TRQ");
  track_selection[2004] = new TString("e- tracks passing MVA TRQ + PID batch 1 wt");
  track_selection[2005] = new TString("e- tracks passing MVA TRQ + PID batch 2 wt");
  track_selection[2010] = new TString("e- tracks passing MVA TRQ + PID, T>700");
  track_selection[2014] = new TString("e- tracks passing MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[2015] = new TString("e- tracks passing MVA TRQ + PID, T>700 batch 2 wt");

  track_selection[2020] = new TString("e- CD3 cut MVA TRQ + PID");
  track_selection[2021] = new TString("e- CD3 cut MVA TRQ + PID, batch 1 wt");
  track_selection[2022] = new TString("e- CD3 cut MVA TRQ + PID, batch 2 wt");

  track_selection[2024] = new TString("e- old CD3 cut MVA TRQ + PID");
  track_selection[2025] = new TString("e- old CD3 cut MVA TRQ + PID, batch 1 wt");
  track_selection[2026] = new TString("e- old CD3 cut MVA TRQ + PID, batch 2 wt");

  track_selection[2500] = new TString("all mu- tracks passing analysis_MVA cuts");
  track_selection[2501] = new TString("all mu- tracks passing analysis_MVA cuts, T0>700");

  track_selection[3000] = new TString("e+ tracks passing BOX TRQ + PID");
  track_selection[3004] = new TString("e+ tracks passing BOX TRQ + PID batch 1 wt");
  track_selection[3005] = new TString("e+ tracks passing BOX TRQ + PID batch 2 wt");
  track_selection[3010] = new TString("e+ tracks passing BOX TRQ + PID, T>700");
  track_selection[3014] = new TString("e+ tracks passing BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[3015] = new TString("e+ tracks passing BOX TRQ + PID, T>700 batch 2 wt");

  track_selection[3500] = new TString("all mu+ tracks passing analysis_box cuts");
  track_selection[3501] = new TString("all mu+ tracks passing analysis_box cuts, T0>700");

  track_selection[4000] = new TString("e+ tracks passing MVA TRQ + PID");
  track_selection[4001] = new TString("e+ tracks passing MVA TRQ");
  track_selection[4004] = new TString("e+ tracks passing MVA TRQ + PID batch 1 wt");
  track_selection[4005] = new TString("e+ tracks passing MVA TRQ + PID batch 2 wt");
  track_selection[4010] = new TString("e+ tracks passing MVA TRQ + PID, T>700");
  track_selection[4014] = new TString("e+ tracks passing MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[4015] = new TString("e+ tracks passing MVA TRQ + PID, T>700 batch 2 wt");

  track_selection[4020] = new TString("e+ CD3 cut MVA TRQ + PID");
  track_selection[4021] = new TString("e+ CD3 cut MVA TRQ + PID, batch 1 wt");
  track_selection[4022] = new TString("e+ CD3 cut MVA TRQ + PID, batch 2 wt");

  track_selection[4500] = new TString("all mu+ tracks passing analysis_MVA cuts");
  track_selection[4501] = new TString("all mu+ tracks passing analysis_MVA cuts, T0>700");

  for (int i=0; i<kNTrackHistSets; i++) {
    if (track_selection[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title    = track_selection[i]->Data();
      if (! fol) fol  = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fTrack[i] = new Mu2eII::TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book Genp histograms
//-----------------------------------------------------------------------------
  TString* genp_selection[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) genp_selection[i] = 0;

  genp_selection[0] = new TString("all events");

  for (int i=0; i<kNGenpHistSets; i++) {
    if (genp_selection[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title   = genp_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }
}

//_____________________________________________________________________________
void TPbarAnaModule::FillHistograms() {

  double wt_b1(fEventWeight), wt_b2(fEventWeight);
  if(fBatchMode == 2)
    wt_b1 *= fEvtPar.fOneBatchWeight / fEvtPar.fTwoBatchWeight;
  if(fBatchMode == 1)
    wt_b2 *= fEvtPar.fTwoBatchWeight / fEvtPar.fOneBatchWeight;

//-----------------------------------------------------------------------------
// 1. fill event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);

  if (fEvtPar.fNTracksDe> 0) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
//-----------------------------------------------------------------------------
// 2. fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fEvtPar.fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
  }
//-----------------------------------------------------------------------------
// 3. Simp histograms
//-----------------------------------------------------------------------------
  if (fSimPar.fParticle) {
    FillSimpHistograms(fHist.fSimp[0],fSimPar.fParticle);
  }
//-----------------------------------------------------------------------------
// 4. track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  Mu2eII::TrackPar_t*  tp;

  for (int i=0; i<fEvtPar.fNTracksDe; ++i ) {
    trk = fTrackBlock->Track(i);
    tp  = fTrackPar+i;

    FillTrackHistograms(fHist.fTrack[0],trk,tp,&fSimPar);
    if (trk->Charge() < 0) FillTrackHistograms(fHist.fTrack[1],trk,tp,&fSimPar);
    else FillTrackHistograms(fHist.fTrack[2],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// good track, BOX cuts + PID
//-----------------------------------------------------------------------------
    if ((tp->fIDWord[0] == 0) && (tp->fPidMvaOut[0] > 0.5)) { // "good" track, MVA trained on DAR tracks
      if (trk->Charge() < 0) { 
	FillTrackHistograms(fHist.fTrack[1000],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[1004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[1005],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[1010],trk,tp,&fSimPar);
	  FillTrackHistograms(fHist.fTrack[1014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[1015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	}
      }
      else {
	FillTrackHistograms(fHist.fTrack[3000],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[3004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[3005],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[3010],trk,tp,&fSimPar);
	  FillTrackHistograms(fHist.fTrack[3014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[3015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	}
      }

    }
    else if(tp->fIDWord[0] == 0) {
      // muons
      if (trk->Charge() < 0) { 
	// mu-
	FillTrackHistograms(fHist.fTrack[1500],trk,tp,&fSimPar);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[1501],trk,tp,&fSimPar);
      }
      else {
	// mu+
	FillTrackHistograms(fHist.fTrack[3500],trk,tp,&fSimPar);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[3501],trk,tp,&fSimPar);
      }
    }
//-----------------------------------------------------------------------------
// good track
//-----------------------------------------------------------------------------
    if (tp->fIDWord[1] == 0) { // "good" track, MVA trained on DAR tracks
      if (trk->Charge() < 0) FillTrackHistograms(fHist.fTrack[2001],trk,tp,&fSimPar);
      else FillTrackHistograms(fHist.fTrack[4001],trk,tp,&fSimPar);
    }
//-----------------------------------------------------------------------------
// good track, MVA cuts + PID
//-----------------------------------------------------------------------------
    if ((tp->fIDWord[1] == 0) && (tp->fPidMvaOut[0] > 0.5)) { // "good" track, MVA trained on DAR tracks
      if (trk->Charge() < 0) { 
	FillTrackHistograms(fHist.fTrack[2000],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[2004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[2005],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[2010],trk,tp,&fSimPar);
	  FillTrackHistograms(fHist.fTrack[2014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[2015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	  //
	  if (trk->P()>103.85 && trk->P()<104.9){
	    FillTrackHistograms(fHist.fTrack[2020],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[2021],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	    FillTrackHistograms(fHist.fTrack[2022],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	  }
	  if (trk->P()>103.85 && trk->P()<105.1){
	    FillTrackHistograms(fHist.fTrack[2024],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[2025],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	    FillTrackHistograms(fHist.fTrack[2026],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	  }
	}
      }
      else {
	FillTrackHistograms(fHist.fTrack[4000],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[4004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[4005],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[4010],trk,tp,&fSimPar);
	  FillTrackHistograms(fHist.fTrack[4014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[4015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	  //
	  if (trk->P()>90.5 && trk->P()<92.5){
	    FillTrackHistograms(fHist.fTrack[4020],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[4021],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	    FillTrackHistograms(fHist.fTrack[4022],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	  }
	}
      }
    }
    else if(tp->fIDWord[1] == 0) {
      // muons
      if (trk->Charge() < 0) { 
	// mu-
	FillTrackHistograms(fHist.fTrack[2500],trk,tp,&fSimPar);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[2501],trk,tp,&fSimPar);
      }
      else {
	// mu+
	FillTrackHistograms(fHist.fTrack[4500],trk,tp,&fSimPar);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[4501],trk,tp,&fSimPar);
      }
    }
  } //end track loop
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TPbarAnaModule::Event(int ientry) {

  //  double                xs, p;
  //  TEmuLogLH::PidData_t  dat;
  //  TStnTrack*            track;
  //  int                   id_word;
  TLorentzVector        mom;

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock  ->GetEntry(ientry);
  //  fTriggerBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
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
  fEvtPar.fNStrawHits = GetHeaderBlock()->fNStrawHits;
  fEvtPar.fInstLum = GetHeaderBlock()->fInstLum;
  fEvtPar.fOneBatchWeight = BatchModeWeight(fEvtPar.fInstLum, 1); //1 batch mode
  fEvtPar.fTwoBatchWeight = BatchModeWeight(fEvtPar.fInstLum, 2); //2 batch mode
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
					// may want to revisit the definition of fSimp

  fSimp             = fSimpBlock->Particle(0);
  fSimPar.fParticle = fSimp;
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = fEvtPar.fParticle;

  if (fSimPar.fGenp) fSimPar.fEleE = fSimPar.fGenp->Energy();
  else               fSimPar.fEleE = -1;
//-----------------------------------------------------------------------------
// virtual detectors - for fSimp need parameters at the tracker front
//-----------------------------------------------------------------------------
  int nsteps = fSpmcBlockVDet->NStepPoints();
  
  for (int i=0; i<nsteps; i++) {
    TStepPointMC* step = fSpmcBlockVDet->StepPointMC(i);
    if (step->PDGCode() == fSimp->fPdgCode) {
      if ((step->VolumeID() == 13) || (step->VolumeID() == 14)) {
	fSimPar.fTFront = step;
      }
      else if ((step->VolumeID() == 11) || (step->VolumeID() == 12)) {
	fSimPar.fTMid = step;
      }
    }
  }

  fEvtPar.fNTracksDe   = fTrackBlock->NTracks();
  fEvtPar.fNGoodTracks = 0;

  if (fEvtPar.fNTracksDe == 0) fTrack = 0;
  else                       fTrack = fTrackBlock->Track(0);


  for (int i=0; i<fEvtPar.fNTracksDe; i++) {
    TrackPar_t*   tp = fTrackPar+i;
    if (fTrackBlockName.Index("TrackBlockPar") == 0) tp->fFitType = 0;  // this better be made more error-proof
    if (fTrackBlockName.Index("TrackBlockDar") == 0) tp->fFitType = 1;
  }

  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,&fSimPar);
  fEventWeight = fGenpBlock->Weight();

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TPbarAnaModule::Debug() {

  TStnTrack* trk;
  int ntrk = fTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock->Track(itrk);
    //    tp  = &fTrackPar[itrk];
//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (trk->fIDWord == 0) {
	TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	if ((vr && (fabs(vr->fDx) > 70) && (fabs(vr->fDx) < 90))) {
	  GetHeaderBlock()->Print(Form("large DX: %f",vr->fDx));
	}
      }
    }
  }
}

//_____________________________________________________________________________
int TPbarAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TPbarAnaModule::Test001() {
}
}
