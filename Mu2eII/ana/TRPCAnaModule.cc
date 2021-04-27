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
#include "Mu2eII/ana/TRPCAnaModule.hh"

// ClassImp(TRPCAnaModule)

namespace Mu2eII {

//-----------------------------------------------------------------------------
TRPCAnaModule::TRPCAnaModule(const char* name, const char* title):
  TAnaModule(name,title)
{
  fPtMin  = 1.;
  fTrackNumber.Set(100);
//-----------------------------------------------------------------------------
// MC truth is defined in TAnaModule
//-----------------------------------------------------------------------------
  fBestID          = 1;                   // best ID word - not used here
  fTrackBlockName  = "TrackBlockDar";
  fNTrkID          = 2;
}

//-----------------------------------------------------------------------------
TRPCAnaModule::~TRPCAnaModule() {
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TRPCAnaModule::BeginJob() {

  TAnaModule::BeginJob();

//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data(),"TStnTrackBlock"      ,&fTrackBlock     );
  RegisterDataBlock("ClusterBlock"        , "TStnClusterBlock"   ,&fClusterBlock   );
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
void TRPCAnaModule::BookHistograms() {

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

  track_selection[   0] = new TString("all e- tracks");
  track_selection[   1] = new TString("all e- tracks, w=SurvProb");

  track_selection[ 100] = new TString("all e- tracks BOX TRQ");
  track_selection[ 101] = new TString("all e- tracks BOX TRQ w=SurvProb");
  track_selection[ 102] = new TString("all e- tracks BOX TRQ w=SurvProb P>80");

  track_selection[ 200] = new TString("all e- tracks MVA TRQ");
  track_selection[ 201] = new TString("all e- tracks MVA TRQ w=SurvProb");
  track_selection[ 202] = new TString("all e- tracks MVA TRQ w=SurvProb P>80");

  track_selection[ 300] = new TString("all e+ tracks");
  track_selection[ 301] = new TString("all e+ tracks w=SurvProb");
  track_selection[ 302] = new TString("all e+ tracks w=SurvProb P>80");

  track_selection[ 400] = new TString("all e+ tracks BOX TRQ");
  track_selection[ 401] = new TString("all e+ tracks BOX TRQ w=SurvProb");
  track_selection[ 402] = new TString("all e+ tracks NOX TRQ w=SurvProb P>80");

  track_selection[ 500] = new TString("all e+ tracks MVA TRQ");
  track_selection[ 501] = new TString("all e+ tracks MVA TRQ w=SurvProb");
  track_selection[ 502] = new TString("all e+ tracks NOX TRQ w=SurvProb P>80");

  track_selection[1000] = new TString("e- tracks BOX TRQ + PID");
  track_selection[1002] = new TString("e- tracks BOX TRQ w=SurvProb P>80");

  track_selection[1004] = new TString("e- tracks BOX TRQ + PID batch 1 wt");
  track_selection[1005] = new TString("e- tracks BOX TRQ + PID batch 2 wt");
  track_selection[1006] = new TString("all e- 2 batch weighted tracks BOX, T>700, 103.85<mom<105.1");

  track_selection[1010] = new TString("e- tracks BOX TRQ + PID, T>700");
  track_selection[1012] = new TString("e- tracks BOX TRQ w=SurvProb P>80");
  track_selection[1014] = new TString("e- tracks BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[1015] = new TString("e- tracks BOX TRQ + PID, T>700 batch 2 wt");

  track_selection[2000] = new TString("e- tracks MVA TRQ + PID");
  track_selection[2002] = new TString("e- tracks MVA TRQ w=SurvProb P>80");
  track_selection[2004] = new TString("e- tracks MVA TRQ + PID batch 1 wt");
  track_selection[2005] = new TString("e- tracks MVA TRQ + PID batch 2 wt");
  track_selection[2006] = new TString("all e- 2 batch weighted tracks MVA TRQ, T>700, 103.85<mom<105.1");

  track_selection[2010] = new TString("e- tracks MVA TRQ + PID, T>700");
  track_selection[2012] = new TString("e- tracks MVA TRQ + PID w=SurvProb P>80");

  track_selection[2014] = new TString("e- tracks MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[2015] = new TString("e- tracks MVA TRQ + PID, T>700 batch 2 wt");

  track_selection[3000] = new TString("e+ tracks BOX TRQ + PID");
  track_selection[3004] = new TString("e+ tracks BOX TRQ + PID batch 1 wt");
  track_selection[3005] = new TString("e+ tracks BOX TRQ + PID batch 2 wt");
  track_selection[3006] = new TString("all e+ 2 batch weighted tracks passing BOX cuts, T>700, 90.85 MeV/c < p < 92.1 MeV/c");
  track_selection[3010] = new TString("e+ tracks BOX TRQ + PID, T>700");
  track_selection[3014] = new TString("e+ tracks BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[3015] = new TString("e+ tracks BOX TRQ + PID, T>700 batch 2 wt");

  track_selection[4000] = new TString("e+ tracks MVA TRQ + PID");
  track_selection[4002] = new TString("e+ tracks MVA TRQ w=SurvProb P>80");
  track_selection[4004] = new TString("e+ tracks MVA TRQ + PID batch 1 wt");
  track_selection[4005] = new TString("e+ tracks MVA TRQ + PID batch 2 wt");
  track_selection[4006] = new TString("all e+ 2 batch weighted tracks passing MVA cuts, T>700, 90.85 MeV/c < p < 92.1 MeV/c");

  track_selection[4010] = new TString("e+ tracks MVA TRQ + PID, T>700");
  track_selection[4012] = new TString("e+ tracks MVA TRQ + PID w=SurvProb P>80");
  track_selection[4014] = new TString("e+ tracks MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[4015] = new TString("e+ tracks MVA TRQ + PID, T>700 batch 2 wt");

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
void TRPCAnaModule::FillHistograms() {

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

    if (trk->Charge() < 0) {
//-----------------------------------------------------------------------------
// electrons
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[0],trk,tp,&fSimPar);
      FillTrackHistograms(fHist.fTrack[1],trk,tp,&fSimPar,fSurvProb);

      if ((tp->fIDWord[0] == 0) && (tp->fPidMvaOut[0] > 0.5)) {  // "good" track, BOX cuts, PID MVA trained on DAR tracks

	FillTrackHistograms(fHist.fTrack[ 100],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[ 101],trk,tp,&fSimPar,fSurvProb);
//-----------------------------------------------------------------------------
// a lazy attempt to get rid of pileup tracks and have a meaningful timing distribution
//-----------------------------------------------------------------------------
	if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[102],trk,tp,&fSimPar,fSurvProb);

	FillTrackHistograms(fHist.fTrack[1000],trk,tp,&fSimPar,fSurvProb);
	if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[1002],trk,tp,&fSimPar,fSurvProb);

	FillTrackHistograms(fHist.fTrack[1004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[1005],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	if (trk->T0() > 700.  and trk->P()<105.1 and trk->P()>103.85) FillTrackHistograms(fHist.fTrack[1006],trk,tp,&fSimPar,fSurvProb);
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[1010],trk,tp,&fSimPar,fSurvProb);
	  if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[1012],trk,tp,&fSimPar,fSurvProb);

	  FillTrackHistograms(fHist.fTrack[1014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[1015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	}
      } //end BOX ID

      if ((tp->fIDWord[1] == 0) &&  (tp->fPidMvaOut[0] > 0.5)) { // "good" track, MVA cuts trained on DAR tracks

	FillTrackHistograms(fHist.fTrack[200],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[201],trk,tp,&fSimPar,fSurvProb);
//-----------------------------------------------------------------------------
// a lazy attempt to get rid of pileup tracks and have a meaningful timing distribution
//-----------------------------------------------------------------------------
	if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[202],trk,tp,&fSimPar,fSurvProb);

	FillTrackHistograms(fHist.fTrack[2000],trk,tp,&fSimPar,fSurvProb);
	if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[2002],trk,tp,&fSimPar,fSurvProb);

	FillTrackHistograms(fHist.fTrack[2004],trk,tp,&fSimPar,wt_b1);   // batch re-weighting
	FillTrackHistograms(fHist.fTrack[2005],trk,tp,&fSimPar,wt_b2);   // batch re-weighting

	if (trk->T0() > 700. and trk->P()<105.1 and trk->P()>103.85 ) FillTrackHistograms(fHist.fTrack[2006],trk,tp,&fSimPar,fSurvProb);

	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[2010],trk,tp,&fSimPar,fSurvProb);
	  if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[2012],trk,tp,&fSimPar,fSurvProb);

	  FillTrackHistograms(fHist.fTrack[2014],trk,tp,&fSimPar,wt_b1); // batch re-weighting
	  FillTrackHistograms(fHist.fTrack[2015],trk,tp,&fSimPar,wt_b2); // batch re-weighting
	}
      } //end MVA ID
    } //end negative track charge
    else {
//-----------------------------------------------------------------------------
// positrons
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[300],trk,tp,&fSimPar);
      FillTrackHistograms(fHist.fTrack[301],trk,tp,&fSimPar,fSurvProb);

      if ((tp->fIDWord[0] == 0) && (tp->fPidMvaOut[0] > 0.5)) {  // "good" track, BOX cuts, PID MVA trained on DAR tracks

	FillTrackHistograms(fHist.fTrack[ 400],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[ 401],trk,tp,&fSimPar,fSurvProb);
//-----------------------------------------------------------------------------
// a lazy attempt to get rid of pileup tracks
//-----------------------------------------------------------------------------
	if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[402],trk,tp,&fSimPar,fSurvProb);

	FillTrackHistograms(fHist.fTrack[3000],trk,tp,&fSimPar,fSurvProb);

	FillTrackHistograms(fHist.fTrack[3004],trk,tp,&fSimPar,wt_b1); // batch re-weighting
	FillTrackHistograms(fHist.fTrack[3005],trk,tp,&fSimPar,wt_b2); // batch re-weighting
	if (trk->T0() > 700. and trk->P()< 92.1  and trk->P()>90.85 ) FillTrackHistograms(fHist.fTrack[3006],trk,tp,&fSimPar,fSurvProb);
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[3010],trk,tp,&fSimPar,fSurvProb);

	  FillTrackHistograms(fHist.fTrack[3014],trk,tp,&fSimPar,wt_b1); // batch re-weighting
	  FillTrackHistograms(fHist.fTrack[3015],trk,tp,&fSimPar,wt_b2); // batch re-weighting
	}
      } //end BOX ID

      if ((tp->fIDWord[1] == 0) &&  (tp->fPidMvaOut[0] > 0.5)) { // "good" track, MVA cuts trained on DAR tracks

	FillTrackHistograms(fHist.fTrack[ 500],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[ 501],trk,tp,&fSimPar,fSurvProb);
	if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[502],trk,tp,&fSimPar,fSurvProb);

	FillTrackHistograms(fHist.fTrack[4000],trk,tp,&fSimPar,fSurvProb);
	if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[2002],trk,tp,&fSimPar,fSurvProb);

	FillTrackHistograms(fHist.fTrack[4004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[4005],trk,tp,&fSimPar,wt_b2); //batch re-weighting

	if (trk->T0() > 700. and trk->P()< 92.1  and trk->P()>90.85  ) FillTrackHistograms(fHist.fTrack[4006],trk,tp,&fSimPar,fSurvProb);

	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[4010],trk,tp,&fSimPar,fSurvProb);
	  if (tp->fP > 80.) FillTrackHistograms(fHist.fTrack[4012],trk,tp,&fSimPar,fSurvProb);

	  FillTrackHistograms(fHist.fTrack[4014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[4015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	}
      } //end MVA ID
    }//end positive tracks
  }//end track loop
 
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TRPCAnaModule::Event(int ientry) {

  //  double                xs, p;
  //  TEmuLogLH::PidData_t  dat;
  //  TStnTrack*            track;
  //  int                   id_word;
  TLorentzVector        mom;

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock  ->GetEntry(ientry);
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
//-----------------------------------------------------------------------------
// for RPC, the weight saved by the generator is supposed to be the pion survival probability
//-----------------------------------------------------------------------------
  fSurvProb         = fGenpBlock->Weight();
  fEventWeight      = fSurvProb;
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
  else                         fTrack = fTrackBlock->Track(0);

  for (int i=0; i<fEvtPar.fNTracksDe; i++) {
    TrackPar_t*   tp = fTrackPar+i;

    if (fTrackBlockName.Index("TrackBlockPar") == 0) tp->fFitType = 0;  // this better be made more error-proof
    if (fTrackBlockName.Index("TrackBlockDar") == 0) tp->fFitType = 1;
  }

  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,&fSimPar);
//-----------------------------------------------------------------------------
// RPC-specific initializations, assume running on RPC dataset
//-----------------------------------------------------------------------------
  for (int i=0; i<fEvtPar.fNTracksDe; i++) {
    fTrackPar[i].fRPCTimeWt = fGenpBlock->Weight();
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TRPCAnaModule::Debug() {

  int ntrk = fTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    TStnTrack*          trk = fTrackBlock->Track(itrk);
    Mu2eII::TrackPar_t* tp  = &fTrackPar[itrk];
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
    
    if (GetDebugBit(4) == 1) {
      if ((tp->fIDWord[1] == 0) and (tp->fPidMvaOut[0] > 0.5)) {
	if ((trk->T0() > 700.) and (tp->fP > 70.) and (tp->fP < 73.5)) {
	  GetHeaderBlock()->Print(Form("TRPCAnaModule::Debug::bit_004: P=%8.3e fSurvProb:%12.5e",tp->fP,fSurvProb));
	}
      }
    }
  }
}

//_____________________________________________________________________________
int TRPCAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TRPCAnaModule::Test001() {
}
}
