//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
// bit 000  : all events
// bit 001  : passed events
// bit 002  : rejected events
// ----------------------------- 
// bit 003  : events with set C tracks and 70mm < |dx|  < 90 mm
// bit 004  : events with tp->fIDWord[1] == 0 , tp->fDpf > 3 MeV (look at the resolution tail)
// bit 051  : events with tp->fIDWord[1] == 0 and fEvtPar.fCosmicVeto != 0
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
#include "Mu2eII/ana/TTrackAnaModule.hh"

// ClassImp(TTrackAnaModule)

namespace Mu2eII {

//-----------------------------------------------------------------------------
TTrackAnaModule::TTrackAnaModule(const char* name, const char* title):
  TAnaModule(name,title)
{
  // fPtMin          = 1.;
  // fTrackNumber.Set(100);
//-----------------------------------------------------------------------------
// fTrackID[0] reserved for box cuts, [1] for MVA; fBestID is not used in this module
//-----------------------------------------------------------------------------
//  fBestID         = 1;                   // best ID word for TRQ MVA

  fTrackBlockNameDe   = "TrackBlockDarDe";
  fTrackBlockNameUe   = "TrackBlockDarUe";
                                         // best ID word for TRQ MVA
  TAnaModule::fNTrkID = 2;
}

//-----------------------------------------------------------------------------
TTrackAnaModule::~TTrackAnaModule() {
  // delete fTrackBlockName;
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrackAnaModule::BeginJob() {

  TAnaModule::BeginJob();
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockNameDe.Data(), "TStnTrackBlock"      , &fTrackBlockDe   );
  RegisterDataBlock(fTrackBlockNameUe.Data(), "TStnTrackBlock"      , &fTrackBlockUe   );
  RegisterDataBlock("HelixBlockDe"          , "TStnHelixBlock"      , &fHelixBlockDe   );
  RegisterDataBlock("HelixBlockUe"          , "TStnHelixBlock"      , &fHelixBlockUe   );
  RegisterDataBlock("TCFinderBlockUe"       , "TStnTimeClusterBlock", &fTCFinderBlockUe);

  RegisterDataBlock("ClusterBlock"          , "TStnClusterBlock" , &fClusterBlock );

  RegisterDataBlock("GenpBlock"             , "TGenpBlock"       , &fGenpBlock    );
  RegisterDataBlock("SimpBlock"             , "TSimpBlock"       , &fSimpBlock    );
  RegisterDataBlock("SpmcBlockVDet"         , "TStepPointMCBlock", &fSpmcBlockVDet);

					// initialize structure for non-CRV cosmics tagging

  fCosmicVetoData.fTrackBlockDe    = fTrackBlockDe;
  fCosmicVetoData.fTrackParDe      = fTrackPar;         // defined for De only
  fCosmicVetoData.fTrackBlockUe    = fTrackBlockUe;
  fCosmicVetoData.fClusterBlock    = fClusterBlock;
  fCosmicVetoData.fHelixBlockDe    = fHelixBlockDe;
  fCosmicVetoData.fHelixBlockUe    = fHelixBlockUe;
  fCosmicVetoData.fTCFinderBlockUe = fTCFinderBlockUe;
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

    tp->fTrackID[0]  = fTrackID_BOX;                     // these poiters need to be set just once
    tp->fTrackID[1]  = fTrackID_MVA;
  }

  return 0;
}
//-----------------------------------------------------------------------------
void TTrackAnaModule::BookHistograms() {
  
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
  track_selection[   1] = new TString("all tracks DIO LO");
  track_selection[   2] = new TString("all tracks DIO LL");
  track_selection[  10] = new TString("all tracks (T0>700)");
  track_selection[  11] = new TString("all tracks (T0>700) DIO LO");
  track_selection[  12] = new TString("all tracks (T0>700) DIO LL");

  track_selection[ 100] = new TString("all tracks passing box cuts ");
  track_selection[ 101] = new TString("all tracks passing box cuts  DIOLO");
  track_selection[ 102] = new TString("all tracks passing box cuts  DIOLL");
  track_selection[ 103] = new TString("all tracks passing BOX cuts + (T>700)");
  track_selection[ 104] = new TString("all tracks passing BOX cuts + (T>700) DIOLL");

  track_selection[ 110] = new TString("all tracks passing box cuts (Michael) 1st disk");
  track_selection[ 111] = new TString("all tracks passing box cuts (Michael) 2nd disk");

  track_selection[ 119] = new TString("all tracks passing BOX TID + pre-PID");

  track_selection[ 200] = new TString("all tracks passing MVA TID");
  track_selection[ 201] = new TString("all tracks passing MVA TID DIOLO");
  track_selection[ 202] = new TString("all tracks passing MVA TID DIOLL");
  track_selection[ 203] = new TString("all tracks passing MVA TID lumiwt");
  track_selection[ 204] = new TString("all tracks passing MVA TID DIOLL*lumiwt");
  track_selection[ 205] = new TString("all tracks passing MVA TID + (T>700)");
  track_selection[ 206] = new TString("all tracks passing MVA TID + (T>700) lumiwt");
  track_selection[ 207] = new TString("all tracks passing MVA TID + (T>700) DIOLL");
  track_selection[ 208] = new TString("all tracks passing MVA TID + (T>700) DIOLL*lumiwt");

  track_selection[ 210] = new TString("all tracks passing TRQ MVA  1st disk");
  track_selection[ 211] = new TString("all tracks passing TRQ MVA  2nd disk");

  track_selection[ 219] = new TString("all tracks passing MVA TID + pre-PID");

  track_selection[ 220] = new TString("all tracks passing TRQ+PID MVA cuts");
  track_selection[ 222] = new TString("all tracks passing TRQ+PID MVA cuts w/DIO LL");

  track_selection[ 230] = new TString("all tracks passing TRQ+PID MVA cuts+(T0>700)");
  track_selection[ 232] = new TString("all tracks passing TRQ+PID MVA cuts+(T0>700) w/DIO LL");

  track_selection[ 240] = new TString("all tracks passing TRQ+PID MVA cuts + (T0>700) + CR cuts");
  track_selection[ 242] = new TString("all tracks passing TRQ+PID MVA cuts + (T0>700) + CR cuts w/DIO LL");

  track_selection[1000] = new TString("e- tracks passing BOX TRQ + PID");
  track_selection[1001] = new TString("e- tracks passing BOX TRQ + PID *(Legacy lumiwt)");
  track_selection[1002] = new TString("e- tracks passing BOX TRQ + PID DIOLL");
  track_selection[1003] = new TString("e- tracks passing BOX TRQ + PID DIOLL*(Legacy lumiwt)");
  track_selection[1004] = new TString("e- tracks passing BOX TRQ + PID batch 1 wt");
  track_selection[1005] = new TString("e- tracks passing BOX TRQ + PID batch 2 wt");
  track_selection[1006] = new TString("e- tracks passing BOX TRQ + PID DIOLL*batch 1 wt");
  track_selection[1007] = new TString("e- tracks passing BOX TRQ + PID DIOLL*batch 2 wt");
  track_selection[1010] = new TString("e- tracks passing BOX TRQ + PID (T>700)");
  track_selection[1011] = new TString("e- tracks passing BOX TRQ + PID (T>700)*(Legacy lumiwt)");
  track_selection[1012] = new TString("e- tracks passing BOX TRQ + PID (T>700) DIOLL");
  track_selection[1013] = new TString("e- tracks passing BOX TRQ + PID (T>700) DIOLL*(Legacy lumiwt)");
  track_selection[1014] = new TString("e- tracks passing BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[1015] = new TString("e- tracks passing BOX TRQ + PID, T>700 batch 2 wt");
  track_selection[1016] = new TString("e- tracks passing BOX TRQ + PID, T>700 DIOLL*batch 1 wt");
  track_selection[1017] = new TString("e- tracks passing BOX TRQ + PID, T>700 DIOLL*batch 2 wt");

  track_selection[2000] = new TString("e- tracks passing MVA TRQ + PID");
  track_selection[2001] = new TString("e- tracks passing MVA TRQ + PID*(Legacy lumiwt)");
  track_selection[2002] = new TString("e- tracks passing MVA TRQ + PID DIOLL");
  track_selection[2003] = new TString("e- tracks passing MVA TRQ + PID DIOLL*(Legacy lumiwt)");
  track_selection[2004] = new TString("e- tracks passing MVA TRQ + PID batch 1 wt");
  track_selection[2005] = new TString("e- tracks passing MVA TRQ + PID batch 2 wt");
  track_selection[2006] = new TString("e- tracks passing MVA TRQ + PID DIOLL*batch 1 wt");
  track_selection[2007] = new TString("e- tracks passing MVA TRQ + PID DIOLL*batch 2 wt");
  track_selection[2010] = new TString("e- tracks passing MVA TRQ + PID (T>700)");
  track_selection[2011] = new TString("e- tracks passing MVA TRQ + PID (T>700)*(Legacy lumiwt)");
  track_selection[2012] = new TString("e- tracks passing MVA TRQ + PID (T>700) DIOLL");
  track_selection[2013] = new TString("e- tracks passing MVA TRQ + PID (T>700) DIOLL*(Legacy lumiwt)");
  track_selection[2014] = new TString("e- tracks passing MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[2015] = new TString("e- tracks passing MVA TRQ + PID, T>700 batch 2 wt");
  track_selection[2016] = new TString("e- tracks passing MVA TRQ + PID, T>700 DIOLL*batch 1 wt");
  track_selection[2017] = new TString("e- tracks passing MVA TRQ + PID, T>700 DIOLL*batch 2 wt");

  track_selection[2020] = new TString("e- tracks passing MVA TRQ + PID (T>700) + non-CRV cosmic veto");

  track_selection[3000] = new TString("e+ tracks passing BOX TRQ + PID");
  track_selection[3001] = new TString("e+ tracks passing BOX TRQ + PID wt_b2");
  track_selection[3002] = new TString("e+ tracks passing BOX TRQ + PID DIOLL");
  track_selection[3003] = new TString("e+ tracks passing BOX TRQ + PID DIOLL*(Legacy lumiwt)");
  track_selection[3004] = new TString("e+ tracks passing BOX TRQ + PID batch 1 wt");
  track_selection[3005] = new TString("e+ tracks passing BOX TRQ + PID batch 2 wt");
  track_selection[3006] = new TString("e+ tracks passing BOX TRQ + PID DIOLL*batch 1 wt");
  track_selection[3007] = new TString("e+ tracks passing BOX TRQ + PID DIOLL*batch 2 wt");
  track_selection[3010] = new TString("e+ tracks passing BOX TRQ + PID (T>700)");
  track_selection[3011] = new TString("e+ tracks passing BOX TRQ + PID (T>700)*(Legacy lumiwt)");
  track_selection[3012] = new TString("e+ tracks passing BOX TRQ + PID (T>700) DIOLL");
  track_selection[3013] = new TString("e+ tracks passing BOX TRQ + PID (T>700) DIOLL*(Legacy lumiwt)");
  track_selection[3014] = new TString("e+ tracks passing BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[3015] = new TString("e+ tracks passing BOX TRQ + PID, T>700 batch 2 wt");
  track_selection[3016] = new TString("e+ tracks passing BOX TRQ + PID, T>700 DIOLL*batch 1 wt");
  track_selection[3017] = new TString("e+ tracks passing BOX TRQ + PID, T>700 DIOLL*batch 2 wt");

  track_selection[4000] = new TString("e+ tracks passing MVA TRQ + PID");
  track_selection[4001] = new TString("e+ tracks passing MVA TRQ + PID*(Legacy lumiwt)");
  track_selection[4002] = new TString("e+ tracks passing MVA TRQ + PID DIOLL");
  track_selection[4003] = new TString("e+ tracks passing MVA TRQ + PID DIOLL*(Legacy lumiwt)");
  track_selection[4004] = new TString("e+ tracks passing MVA TRQ + PID batch 1 wt");
  track_selection[4005] = new TString("e+ tracks passing MVA TRQ + PID batch 2 wt");
  track_selection[4006] = new TString("e- tracks passing MVA TRQ + PID DIOLL*batch 1 wt");
  track_selection[4007] = new TString("e- tracks passing MVA TRQ + PID DIOLL*batch 2 wt");
  track_selection[4010] = new TString("e+ tracks passing MVA TRQ + PID (T>700)");
  track_selection[4011] = new TString("e+ tracks passing MVA TRQ + PID (T>700)*(Legacy lumiwt)");
  track_selection[4012] = new TString("e+ tracks passing MVA TRQ + PID (T>700) DIOLL");
  track_selection[4013] = new TString("e+ tracks passing MVA TRQ + PID (T>700) DIOLL*(Legacy lumiwt)");
  track_selection[4014] = new TString("e+ tracks passing MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[4015] = new TString("e+ tracks passing MVA TRQ + PID, T>700 batch 2 wt");
  track_selection[4016] = new TString("e+ tracks passing MVA TRQ + PID, T>700 DIOLL*batch 1 wt");
  track_selection[4017] = new TString("e+ tracks passing MVA TRQ + PID, T>700 DIOLL*batch 2 wt");

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
// book track ID histograms
//-----------------------------------------------------------------------------
  TString* track_id_selection[kNTrackIDHistSets];
  for (int i=0; i<kNTrackIDHistSets; i++) track_id_selection[i] = 0;

  track_id_selection[  0] = new TString("e- tracks, 103.85 <P< 105.1, BOX");
  track_id_selection[  1] = new TString("e- tracks, 103.85 <P< 105.1, MVA");
  track_id_selection[  2] = new TString("e+ tracks, 103.85 <P< 105.1, BOX");
  track_id_selection[  3] = new TString("e+ tracks, 103.85 <P< 105.1, MVA");

  track_id_selection[  4] = new TString("e- tracks, 90.85<P<92.1, BOX");
  track_id_selection[  5] = new TString("e- tracks, 90.85<P<92.1, MVA");
  track_id_selection[  6] = new TString("e+ tracks, 90.85<P<92.1, BOX");
  track_id_selection[  7] = new TString("e+ tracks, 90.85<P<92.1, MVA");

  for (int i=0; i<kNTrackIDHistSets; i++) {
    if (track_id_selection[i] != 0) {
      sprintf(folder_name,"tid_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title   = track_id_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fTrackID[i] = new TStnTrackID::Hist_t;
      BookTrackIDHistograms(fHist.fTrackID[i],Form("Hist/%s",folder_name));
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

void TTrackAnaModule::FillHistograms() {

  double legacy_lumiwt = fEvtPar.fLumiWt; //legacy lumi re-weighting

  double wt_b1(fEventWeight), wt_b2(fEventWeight); //weights for one and two batch mode luminosity respectively
  if(fBatchMode == 2) //sample generated using two batch mode
    wt_b1 *= fEvtPar.fOneBatchWeight / fEvtPar.fTwoBatchWeight;
  if(fBatchMode == 1) //sample generated using one batch mode
    wt_b2 *= fEvtPar.fTwoBatchWeight / fEvtPar.fOneBatchWeight;

//-----------------------------------------------------------------------------
// 1. fill event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);

  if (fEvtPar.fNTracksDe > 0) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
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
  for (int i=0; i<fEvtPar.fNTracksDe; ++i ) {
    TStnTrack* trk          = fTrackBlockDe->Track(i);
    Mu2eII::TrackPar_t* tp  = fTrackPar+i;

    FillTrackHistograms(fHist.fTrack[0],trk,tp,&fSimPar);
    FillTrackHistograms(fHist.fTrack[1],trk,tp,&fSimPar,tp->fDioLOWt);
    FillTrackHistograms(fHist.fTrack[2],trk,tp,&fSimPar,tp->fDioLLWt);

    if (trk->T0() > 700.) {
      FillTrackHistograms(fHist.fTrack[10],trk,tp,&fSimPar);
      FillTrackHistograms(fHist.fTrack[11],trk,tp,&fSimPar,tp->fDioLOWt);
      FillTrackHistograms(fHist.fTrack[12],trk,tp,&fSimPar,tp->fDioLLWt);
    }
//-----------------------------------------------------------------------------
// good track, BOX cuts
//-----------------------------------------------------------------------------
    if (tp->fIDWord[0] == 0) {
      FillTrackHistograms(fHist.fTrack[ 100],trk,tp,&fSimPar);
      FillTrackHistograms(fHist.fTrack[ 101],trk,tp,&fSimPar,tp->fDioLOWt);
      FillTrackHistograms(fHist.fTrack[ 102],trk,tp,&fSimPar,tp->fDioLLWt);

      if (trk->T0() > 700.) { 
	FillTrackHistograms(fHist.fTrack[ 103],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[ 104],trk,tp,&fSimPar,tp->fDioLLWt);
      }

      if (tp->fDiskID == 0) FillTrackHistograms(fHist.fTrack[110],trk,tp,&fSimPar);
      if (tp->fDiskID == 1) FillTrackHistograms(fHist.fTrack[111],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// particle ID: the PID ANN has been trained with DAR tracks
//-----------------------------------------------------------------------------
      if ((fabs(tp->fTchDr) <  100) && (fabs(tp->fTchDt) < 10.) && 
          (tp->fTchDz       >  -50) && (tp->fTchDz       < 250) && (tp->fEp < 1.2)) {
	
	FillTrackHistograms(fHist.fTrack[ 119],trk,tp,&fSimPar);

	if (tp->fPidMvaOut[0] > 0.5) {
	  if (trk->Charge() < 0) {
	    FillTrackHistograms(fHist.fTrack[1000],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[1001],trk,tp,&fSimPar,legacy_lumiwt);
	    FillTrackHistograms(fHist.fTrack[1002],trk,tp,&fSimPar,tp->fDioLLWt);
	    FillTrackHistograms(fHist.fTrack[1003],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
	    FillTrackHistograms(fHist.fTrack[1004],trk,tp,&fSimPar,wt_b1);
	    FillTrackHistograms(fHist.fTrack[1005],trk,tp,&fSimPar,wt_b2);
	    FillTrackHistograms(fHist.fTrack[1006],trk,tp,&fSimPar,wt_b1*tp->fDioLLWt);
	    FillTrackHistograms(fHist.fTrack[1007],trk,tp,&fSimPar,wt_b2*tp->fDioLLWt);

	    if (trk->T0() > 700.  ) {
	      FillTrackHistograms(fHist.fTrack[1010],trk,tp,&fSimPar);
	      FillTrackHistograms(fHist.fTrack[1011],trk,tp,&fSimPar,legacy_lumiwt);
	      FillTrackHistograms(fHist.fTrack[1012],trk,tp,&fSimPar,tp->fDioLLWt);
	      FillTrackHistograms(fHist.fTrack[1013],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
	      FillTrackHistograms(fHist.fTrack[1014],trk,tp,&fSimPar,wt_b1);
	      FillTrackHistograms(fHist.fTrack[1015],trk,tp,&fSimPar,wt_b2);
	      FillTrackHistograms(fHist.fTrack[1016],trk,tp,&fSimPar,wt_b1*tp->fDioLLWt);
	      FillTrackHistograms(fHist.fTrack[1017],trk,tp,&fSimPar,wt_b2*tp->fDioLLWt);
	    }
	  } 
	  else {
	    FillTrackHistograms(fHist.fTrack[3000],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[3001],trk,tp,&fSimPar,legacy_lumiwt);
	    FillTrackHistograms(fHist.fTrack[3002],trk,tp,&fSimPar,tp->fDioLLWt);
	    FillTrackHistograms(fHist.fTrack[3003],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
	    FillTrackHistograms(fHist.fTrack[3004],trk,tp,&fSimPar,wt_b1);
	    FillTrackHistograms(fHist.fTrack[3005],trk,tp,&fSimPar,wt_b2);
	    FillTrackHistograms(fHist.fTrack[3006],trk,tp,&fSimPar,wt_b1*tp->fDioLLWt);
	    FillTrackHistograms(fHist.fTrack[3007],trk,tp,&fSimPar,wt_b2*tp->fDioLLWt);

	    if (trk->T0() > 700.) {
	      FillTrackHistograms(fHist.fTrack[3010],trk,tp,&fSimPar);
	      FillTrackHistograms(fHist.fTrack[3011],trk,tp,&fSimPar,legacy_lumiwt);
	      FillTrackHistograms(fHist.fTrack[3012],trk,tp,&fSimPar,tp->fDioLLWt);
	      FillTrackHistograms(fHist.fTrack[3013],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
	      FillTrackHistograms(fHist.fTrack[3014],trk,tp,&fSimPar,wt_b1);
	      FillTrackHistograms(fHist.fTrack[3015],trk,tp,&fSimPar,wt_b2);
	      FillTrackHistograms(fHist.fTrack[3016],trk,tp,&fSimPar,wt_b1*tp->fDioLLWt);
	      FillTrackHistograms(fHist.fTrack[3017],trk,tp,&fSimPar,wt_b2*tp->fDioLLWt);
	    }
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// good track, MVA cuts (PAR:Dave/Andy or DAR:on_the_fly)
//-----------------------------------------------------------------------------
    if (tp->fIDWord[1] == 0) {
      FillTrackHistograms(fHist.fTrack[ 200],trk,tp,&fSimPar);
      FillTrackHistograms(fHist.fTrack[ 201],trk,tp,&fSimPar,tp->fDioLOWt);
      FillTrackHistograms(fHist.fTrack[ 202],trk,tp,&fSimPar,tp->fDioLLWt);
      FillTrackHistograms(fHist.fTrack[ 203],trk,tp,&fSimPar,legacy_lumiwt);
      FillTrackHistograms(fHist.fTrack[ 204],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
      
      if (trk->T0() > 700.) { 
	FillTrackHistograms(fHist.fTrack[ 205],trk,tp,&fSimPar);
	FillTrackHistograms(fHist.fTrack[ 206],trk,tp,&fSimPar,legacy_lumiwt);
	FillTrackHistograms(fHist.fTrack[ 207],trk,tp,&fSimPar,tp->fDioLLWt);
	FillTrackHistograms(fHist.fTrack[ 208],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
      }

      if (tp->fDiskID == 0) FillTrackHistograms(fHist.fTrack[210],trk,tp,&fSimPar);
      if (tp->fDiskID == 1) FillTrackHistograms(fHist.fTrack[211],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// particle ID: the PID ANN has been trained with DAR tracks
// prototype the preselection cuts - ere they need to be listed explicitly
//-----------------------------------------------------------------------------
      if ((fabs(tp->fTchDr) <  100) and (fabs(tp->fTchDt) <  10.) and
          (tp->fTchDz       >  -50) and (tp->fTchDz       <  250) and 
	  (tp->fEp          >    0) and (tp->fEp          < 1.05)     ) {

	FillTrackHistograms(fHist.fTrack[ 219],trk,tp,&fSimPar);
				       
	if (tp->fPidMvaOut[0] > 0.5) {
	  FillTrackHistograms(fHist.fTrack[ 220],trk,tp,&fSimPar);
	  FillTrackHistograms(fHist.fTrack[ 222],trk,tp,&fSimPar,tp->fDioLLWt);
	  if (trk->T0() > 700.) {
	    FillTrackHistograms(fHist.fTrack[ 230],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[ 232],trk,tp,&fSimPar,tp->fDioLLWt);
	  }
//-----------------------------------------------------------------------------
// add explicit timing cut T > 700 
// 2021-01-01: some work still remains on cosmics rejection
//-----------------------------------------------------------------------------
	  if (trk->Charge() < 0) {
	    FillTrackHistograms(fHist.fTrack[2000],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[2001],trk,tp,&fSimPar,legacy_lumiwt);
	    FillTrackHistograms(fHist.fTrack[2002],trk,tp,&fSimPar,tp->fDioLLWt);
	    FillTrackHistograms(fHist.fTrack[2003],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
	    FillTrackHistograms(fHist.fTrack[2004],trk,tp,&fSimPar,wt_b1);
	    FillTrackHistograms(fHist.fTrack[2005],trk,tp,&fSimPar,wt_b2);
	    FillTrackHistograms(fHist.fTrack[2006],trk,tp,&fSimPar,wt_b1*tp->fDioLLWt);
	    FillTrackHistograms(fHist.fTrack[2007],trk,tp,&fSimPar,wt_b2*tp->fDioLLWt);

	    if (trk->T0() > 700.  ) {
	      FillTrackHistograms(fHist.fTrack[2010],trk,tp,&fSimPar);
	      FillTrackHistograms(fHist.fTrack[2011],trk,tp,&fSimPar,legacy_lumiwt);
	      FillTrackHistograms(fHist.fTrack[2012],trk,tp,&fSimPar,tp->fDioLLWt);
	      FillTrackHistograms(fHist.fTrack[2013],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
	      FillTrackHistograms(fHist.fTrack[2014],trk,tp,&fSimPar,wt_b1);
	      FillTrackHistograms(fHist.fTrack[2015],trk,tp,&fSimPar,wt_b2);
	      FillTrackHistograms(fHist.fTrack[2016],trk,tp,&fSimPar,wt_b1*tp->fDioLLWt);
	      FillTrackHistograms(fHist.fTrack[2017],trk,tp,&fSimPar,wt_b2*tp->fDioLLWt);

	      if (fEvtPar.fCosmicVeto == 0) {
		FillTrackHistograms(fHist.fTrack[2020],trk,tp,&fSimPar);
	      }
	    }
	  } 
	  else {
	    FillTrackHistograms(fHist.fTrack[4000],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[4001],trk,tp,&fSimPar,legacy_lumiwt);
	    FillTrackHistograms(fHist.fTrack[4002],trk,tp,&fSimPar,tp->fDioLLWt);
	    FillTrackHistograms(fHist.fTrack[4003],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
	    FillTrackHistograms(fHist.fTrack[4004],trk,tp,&fSimPar,wt_b1);
	    FillTrackHistograms(fHist.fTrack[4005],trk,tp,&fSimPar,wt_b2);
	    FillTrackHistograms(fHist.fTrack[4006],trk,tp,&fSimPar,wt_b1*tp->fDioLLWt);
	    FillTrackHistograms(fHist.fTrack[4007],trk,tp,&fSimPar,wt_b2*tp->fDioLLWt);

	    if (trk->T0() > 700.) {
	      FillTrackHistograms(fHist.fTrack[4010],trk,tp,&fSimPar);
	      FillTrackHistograms(fHist.fTrack[4011],trk,tp,&fSimPar,legacy_lumiwt);
	      FillTrackHistograms(fHist.fTrack[4012],trk,tp,&fSimPar,tp->fDioLLWt);
	      FillTrackHistograms(fHist.fTrack[4013],trk,tp,&fSimPar,tp->fDioLLWt*legacy_lumiwt);
	      FillTrackHistograms(fHist.fTrack[4014],trk,tp,&fSimPar,wt_b1);
	      FillTrackHistograms(fHist.fTrack[4015],trk,tp,&fSimPar,wt_b2);
	      FillTrackHistograms(fHist.fTrack[4016],trk,tp,&fSimPar,wt_b1*tp->fDioLLWt);
	      FillTrackHistograms(fHist.fTrack[4017],trk,tp,&fSimPar,wt_b2*tp->fDioLLWt);
	    }
	  }
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// track ID histograms: define ID to a 2 MeV momentum window and T>700 ns
// (this is just a benchmark)
//-----------------------------------------------------------------------------
  for (int i=0; i<fEvtPar.fNTracksDe; ++i) {
    TStnTrack* trk          = fTrackBlockDe->Track(i);
    Mu2eII::TrackPar_t* tp  = fTrackPar+i;

    if ((tp->fP > 103.85) && (tp->fP < 105.1) && (trk->T0() > 700)) { 
      if (trk->Charge() < 0) { 
	fTrackID_BOX->FillHistograms(fHist.fTrackID[0],trk,1);
	fTrackID_MVA->FillHistograms(fHist.fTrackID[1],trk,1);
      }
      else {
	fTrackID_BOX->FillHistograms(fHist.fTrackID[2],trk,1);
	fTrackID_MVA->FillHistograms(fHist.fTrackID[3],trk,1);
      }
    }

    if ((tp->fP > 90.85) && (tp->fP < 92.1) && (trk->T0() > 700)) { 
      if (trk->Charge() < 0) { 
	fTrackID_BOX->FillHistograms(fHist.fTrackID[4],trk,1);
	fTrackID_MVA->FillHistograms(fHist.fTrackID[5],trk,1);
      }
      else {
	fTrackID_BOX->FillHistograms(fHist.fTrackID[6],trk,1);
	fTrackID_MVA->FillHistograms(fHist.fTrackID[7],trk,1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TTrackAnaModule::Event(int ientry) {

  fTrackBlockDe   ->GetEntry(ientry);
  fTrackBlockUe   ->GetEntry(ientry);
  fHelixBlockDe   ->GetEntry(ientry);
  fHelixBlockUe   ->GetEntry(ientry);
  fTCFinderBlockUe->GetEntry(ientry);

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
  fEvtPar.fNGenp            = fGenpBlock->NParticles();
  fEvtPar.fParticle         = NULL;
  fEvtPar.fGenE             = -1;

  fEvtPar.fNStrawHits       = GetHeaderBlock()->fNStrawHits;
  fEvtPar.fInstLum          = GetHeaderBlock()->fInstLum;

  fEvtPar.fNTracksDe        = fTrackBlockDe->NTracks();
  fEvtPar.fNTracksUe        = fTrackBlockUe->NTracks();
  fEvtPar.fNHelicesDe       = fHelixBlockDe->NHelices();
  fEvtPar.fNHelicesUe       = fHelixBlockUe->NHelices();
  fEvtPar.fNGoodTracks      = 0;

  fEvtPar.fNCrvClusters     = -1;
  fEvtPar.fNCrvPulses       = -1;
  fEvtPar.fNCrvCoincidences = -1;

  fEvtPar.fOneBatchWeight   = BatchModeWeight(fEvtPar.fInstLum, 1); //1 batch mode
  fEvtPar.fTwoBatchWeight   = BatchModeWeight(fEvtPar.fInstLum, 2); //2 batch mode

  fEvtPar.fLumiWt           = 1.; //legacy lumi re-weighting
  if      (fBatchMode == 1) fEvtPar.fLumiWt = fEvtPar.fOneBatchWeight/fEvtPar.fTwoBatchWeight;
  else if (fBatchMode == 2) fEvtPar.fLumiWt = fEvtPar.fTwoBatchWeight/fEvtPar.fOneBatchWeight;
//-----------------------------------------------------------------------------
// MC generator info
//-----------------------------------------------------------------------------
  TLorentzVector        mom;
  for (int i=fEvtPar.fNGenp-1; i>=0; i--) {
    TGenParticle* genp = fGenpBlock->Particle(i);
    int pdg_code       = genp->GetPdgCode();
    int process_code   = genp->GetStatusCode();
    // printf("%d\t%d\t%d\t%d\n",pdg_code,fPDGCode,process_code,fMCProcessCode);
    if ((pdg_code == fPDGCode) && (process_code == fMCProcessCode)) {
      fEvtPar.fParticle = genp;
      genp->Momentum(mom);
      fEvtPar.fGenE     = mom.Energy();
      break;
    }
  }
//-----------------------------------------------------------------------------
// cache SimP parameters - overlaps with EvtPar, so may want to revisit in future
//-----------------------------------------------------------------------------
  fSimp             = fSimpBlock->Particle(0);
  fSimPar.fParticle = fSimp;
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = fEvtPar.fParticle;
  if (fSimPar.fGenp) fSimPar.fEleE = fSimPar.fGenp->Energy(); 
  else               fSimPar.fEleE = -1;
//-----------------------------------------------------------------------------
// calculate DIO weights once per event
//-----------------------------------------------------------------------------
  if (fEvtPar.fGenE > 0) {
    fEvtPar.fDioLOWt = TStntuple::DioWeightAl   (fEvtPar.fGenE);
    fEvtPar.fDioLLWt = TStntuple::DioWeightAl_LL(fEvtPar.fGenE);
  }
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

  if (fEvtPar.fNTracksDe == 0) fTrack = 0;
  else                         fTrack = fTrackBlockDe->Track(0);
//-----------------------------------------------------------------------------
// add on-the fly MVA calculation to the list
// pre-initialization of TrackPar
//-----------------------------------------------------------------------------
  for (int i=0; i<fEvtPar.fNTracksDe; i++) {
    TrackPar_t*  tp = fTrackPar+i;

    tp->fDioLOWt    = fEvtPar.fDioLOWt;
    tp->fDioLLWt    = fEvtPar.fDioLLWt;
    
    if (fTrackBlockNameDe.Index("TrackBlockPar") == 0) tp->fFitType = 0;  // this better be made more error-proof
    if (fTrackBlockNameDe.Index("TrackBlockDar") == 0) tp->fFitType = 1;
  }

  InitTrackPar(fTrackBlockDe,fClusterBlock,fTrackPar,&fSimPar);

  NonCrvCosmicVeto(&fCosmicVetoData,&fEvtPar);

  // P.M.: not sure what this was for
  fEventWeight = 1.;                                          // fEvtPar.fDioLLWt;

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTrackAnaModule::Debug() {

  int ntrk = fTrackBlockDe->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    TStnTrack*  trk = fTrackBlockDe->Track(itrk);
    TrackPar_t* tp  = fTrackPar+itrk;
//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (trk->fIDWord == 0) {
	TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	if ((vr && (fabs(vr->fDx) > 70) and (fabs(vr->fDx) < 90))) {
	  GetHeaderBlock()->Print(Form("%s::bit_003: large DX: %f",GetName(),vr->fDx));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 4  : events with tp->fIDWord[1] == 0 , tp->fDpf > 3 MeV (look at the resolution tail)
//-----------------------------------------------------------------------------
    if (GetDebugBit(4) == 1) {
      if ((tp->fIDWord[1] == 0) and (tp->fP > 104.5) and (tp->fDpF > 3.)) {
	GetHeaderBlock()->Print(Form("%s::bit_004: tp->fIDWord[1] = 0, tp->fP = %10.3f tp->fDpf = %8.3f ttk->T0() = %10.3f",
				     GetName(),tp->fP, tp->fDpF, trk->T0()));
      }
    }
//-----------------------------------------------------------------------------
// bit 051  : events with fEvtPar.fCosmicVeto != 0
//-----------------------------------------------------------------------------
    if (GetDebugBit(51) == 1) {
      if ((tp->fIDWord[1] == 0) and (fEvtPar.fCosmicVeto != 0)) {
	GetHeaderBlock()->Print(Form("%s::bit_051: tid=0,fEvtPar.fCosmicVeto = 0x%08x",
				     GetName(),fEvtPar.fCosmicVeto));
      }
    }
  }

}

//_____________________________________________________________________________
int TTrackAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTrackAnaModule::Test001() {
}
}
