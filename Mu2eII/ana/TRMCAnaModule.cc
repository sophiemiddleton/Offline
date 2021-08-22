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
#include "Mu2eII/ana/TRMCAnaModule.hh"

// ClassImp(TRMCAnaModule)

namespace Mu2eII {

//-----------------------------------------------------------------------------
TRMCAnaModule::TRMCAnaModule(const char* name, const char* title):
  TAnaModule(name,title)
{
  fPtMin  = 1.;
  fTrackNumber.Set(100);
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
//-----------------------------------------------------------------------------
  fKMax           = 90.1;
  fKinematicLimit = 101.853;
  fFlatInput      = 1;
  fInternalRMC    = 0;
  fSpectrum       = 0;
  fSpectrumParam[0] = 0.;
  fSpectrumParam[1] = 0.;
  fMCProcessCode = 41;
  fPDGCode = 22;
  
  fTrackBlockName = "TrackBlockDar";
  fNTrkID         = 2;
}

//-----------------------------------------------------------------------------
TRMCAnaModule::~TRMCAnaModule() {
}


//_____________________________________________________________________________
void TRMCAnaModule::BookHistograms() {

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
  event_selection[ 2] = new TString("1 batch mode weighted");
  event_selection[ 3] = new TString("2 batch mode weighted");

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
  track_selection[   1] = new TString("all tracks unweighted");
  track_selection[   2] = new TString("all tracks, 1 batch mode weighted");
  track_selection[   3] = new TString("all tracks, 2 batch mode weighted");

  track_selection[1000] = new TString("e- tracks passing BOX TRQ + PID");
  track_selection[1004] = new TString("e- tracks passing BOX TRQ + PID batch 1 wt");
  track_selection[1005] = new TString("e- tracks passing BOX TRQ + PID batch 2 wt");
  track_selection[1010] = new TString("e- tracks passing BOX TRQ + PID, T>700");
  track_selection[1014] = new TString("e- tracks passing BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[1015] = new TString("e- tracks passing BOX TRQ + PID, T>700 batch 2 wt");

  track_selection[2000] = new TString("e- tracks passing MVA TRQ + PID");
  track_selection[2004] = new TString("e- tracks passing MVA TRQ + PID batch 1 wt");
  track_selection[2005] = new TString("e- tracks passing MVA TRQ + PID batch 2 wt");
  track_selection[2010] = new TString("e- tracks passing MVA TRQ + PID, T>700");
  track_selection[2014] = new TString("e- tracks passing MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[2015] = new TString("e- tracks passing MVA TRQ + PID, T>700 batch 2 wt");

  track_selection[3000] = new TString("e+ tracks passing BOX TRQ + PID");
  track_selection[3004] = new TString("e+ tracks passing BOX TRQ + PID batch 1 wt");
  track_selection[3005] = new TString("e+ tracks passing BOX TRQ + PID batch 2 wt");
  track_selection[3010] = new TString("e+ tracks passing BOX TRQ + PID, T>700");
  track_selection[3014] = new TString("e+ tracks passing BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[3015] = new TString("e+ tracks passing BOX TRQ + PID, T>700 batch 2 wt");

  track_selection[4000] = new TString("e+ tracks passing MVA TRQ + PID");
  track_selection[4004] = new TString("e+ tracks passing MVA TRQ + PID batch 1 wt");
  track_selection[4005] = new TString("e+ tracks passing MVA TRQ + PID batch 2 wt");
  track_selection[4010] = new TString("e+ tracks passing MVA TRQ + PID, T>700");
  track_selection[4014] = new TString("e+ tracks passing MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[4015] = new TString("e+ tracks passing MVA TRQ + PID, T>700 batch 2 wt");

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
//-----------------------------------------------------------------------------
// book cluster histograms
//-----------------------------------------------------------------------------
  TString*  cluster_selection [kNClusterHistSets];
  for (int i=0; i<kNClusterHistSets; i++) cluster_selection[i] = 0;

  cluster_selection[ 0] = new TString("all clusters");
  cluster_selection[ 1] = new TString("at least 1 track");
  cluster_selection[ 2] = new TString("0 tracks");
  cluster_selection[ 3] = new TString("Cluster E > 10 MeV");
  cluster_selection[ 4] = new TString("Cluster E > 55 MeV");
  cluster_selection[ 5] = new TString("Cluster E > 80 MeV");
  cluster_selection[ 6] = new TString("Cluster in Disk I");
  cluster_selection[ 7] = new TString("Cluster in Disk II");
  cluster_selection[ 8] = new TString("Cluster R > 500 mm");
  cluster_selection[ 9] = new TString("Cluster R > 500 mm, E > 55 MeV");
  cluster_selection[10] = new TString("Cluster R > 500 mm, Disk I");
  cluster_selection[11] = new TString("Cluster R > 500 mm, E > 55 MeV, Disk I");

  for (int i=0; i<kNClusterHistSets; i++) {
    if (cluster_selection[i] != 0) {
      sprintf(folder_name,"cls_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title   = cluster_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fCluster[i] = new ClusterHist_t;
      BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TRMCAnaModule::BeginJob() {

  TAnaModule::BeginJob();
  
  if(fInternalRMC && fIntSpectrum == 0 && fMCProcessCode != 28) {fMCProcessCode = 42; fPDGCode = 11;}
  if(fInternalRMC && fIntSpectrum == 1) {fMCProcessCode = 28; fPDGCode = 11;}
  if(fMCProcessCode == 45) {fPDGCode = 11;}
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  std::cout << "TRMCAnaModule::" << __func__ << ": Using track block name "
	    << fTrackBlockName.Data() << std::endl;
  RegisterDataBlock(fTrackBlockName.Data(),"TStnTrackBlock"      ,&fTrackBlock     );
  RegisterDataBlock("ClusterBlock"        ,"TStnClusterBlock"    ,&fClusterBlock   );
  RegisterDataBlock("GenpBlock"           ,"TGenpBlock"          ,&fGenpBlock      );
  RegisterDataBlock("SimpBlock"           ,"TSimpBlock"          ,&fSimpBlock      );
  RegisterDataBlock("SpmcBlockVDet"       ,"TStepPointMCBlock"   ,&fSpmcBlockVDet );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

//-----------------------------------------------------------------------------
// initialize spectrum information
//-----------------------------------------------------------------------------
  if(fInternalRMC && fMCProcessCode == 42 && fIntSpectrum == 0) //use external RMC spectrum
    fRMCSpectra = new RMCSpectra(fKMax, fKinematicLimit, fSpectrum, 0, fIntSpectrum);
  else
    fRMCSpectra = new RMCSpectra(fKMax, fKinematicLimit, fSpectrum, fInternalRMC, fIntSpectrum);
  fRMCSpectra->verbose_ = fVerbose;
  fRMCSpectra->InitializeSpectrum();
  
  if(fVerbose > 0) std::cout << "TRMCAnaModule::" << __func__ << ": Using fInternalRMC = " << fInternalRMC
			     << " fIntSpectrum = " << fIntSpectrum << " fExtSpectrum = " << fSpectrum 
			     << " fBatchMode = " << fBatchMode << " fMCProcessCode = " << fMCProcessCode
			     << std::endl;
  return 0;
}


//_____________________________________________________________________________
void TRMCAnaModule::FillHistograms() {

  double wt_b1(fEventWeight), wt_b2(fEventWeight);
  if(fBatchMode == 2)
    wt_b1 *= fEvtPar.fOneBatchWeight / fEvtPar.fTwoBatchWeight;
  if(fBatchMode == 1)
    wt_b2 *= fEvtPar.fTwoBatchWeight / fEvtPar.fOneBatchWeight;

//-----------------------------------------------------------------------------
// 1. fill event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);

  if (fEvtPar.fNTracksDe > 0) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
  double wt_prev = fEventWeight;
  fEventWeight = wt_b1;
  FillEventHistograms(fHist.fEvent[2],&fEvtPar);
  fEventWeight = wt_b2;
  FillEventHistograms(fHist.fEvent[3],&fEvtPar);
  fEventWeight = wt_prev;

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

    FillTrackHistograms(fHist.fTrack[   0],trk,tp,&fSimPar, tp->fTotWt);
    FillTrackHistograms(fHist.fTrack[   1],trk,tp,&fSimPar); //unweighted tracks
    FillTrackHistograms(fHist.fTrack[   2],trk,tp,&fSimPar, wt_b1);
    FillTrackHistograms(fHist.fTrack[   3],trk,tp,&fSimPar, wt_b2);
//-----------------------------------------------------------------------------
// PID preselection cuts are applied in TAnaModule.cc, if they fail, the score is < 0
//-----------------------------------------------------------------------------
    if ((tp->fIDWord[0] == 0) && (tp->fPidMvaOut[0] > 0.5)) { // "good" track, MVA trained on DAR tracks
      if (trk->Charge() < 0) { 
	FillTrackHistograms(fHist.fTrack[1000],trk,tp,&fSimPar,tp->fTotWt);
	FillTrackHistograms(fHist.fTrack[1004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[1005],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[1010],trk,tp,&fSimPar,tp->fTotWt);
	  FillTrackHistograms(fHist.fTrack[1014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[1015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	}
      }
      else {
	FillTrackHistograms(fHist.fTrack[3000],trk,tp,&fSimPar,tp->fTotWt);
	FillTrackHistograms(fHist.fTrack[3004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[3005],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[3010],trk,tp,&fSimPar,tp->fTotWt);
	  FillTrackHistograms(fHist.fTrack[3014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[3015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	}
      }
    }

    if ((tp->fIDWord[1] == 0) && (tp->fPidMvaOut[0] > 0.5)) { // "good" track, MVA trained on DAR tracks
      if (trk->Charge() < 0) {
	FillTrackHistograms(fHist.fTrack[2000],trk,tp,&fSimPar,tp->fTotWt);
	FillTrackHistograms(fHist.fTrack[2004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[2005],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[2010],trk,tp,&fSimPar,tp->fTotWt);
	  FillTrackHistograms(fHist.fTrack[2014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[2015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	}
      }
      else {
	FillTrackHistograms(fHist.fTrack[4000],trk,tp,&fSimPar,tp->fTotWt);
	FillTrackHistograms(fHist.fTrack[4004],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	FillTrackHistograms(fHist.fTrack[4005],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	if (trk->T0() > 700.  ) {
	  FillTrackHistograms(fHist.fTrack[4010],trk,tp,&fSimPar,tp->fTotWt);
	  FillTrackHistograms(fHist.fTrack[4014],trk,tp,&fSimPar,wt_b1); //batch re-weighting
	  FillTrackHistograms(fHist.fTrack[4015],trk,tp,&fSimPar,wt_b2); //batch re-weighting
	}
      }
    }
  } //end track loop

//-----------------------------------------------------------------------------
// 5. cluster histograms 
//-----------------------------------------------------------------------------
  TStnCluster*  cl;
  int           id;
  for (int i=0; i<fEvtPar.fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);
    id = cl->DiskID();
    FillClusterHistograms(fHist.fCluster[0],cl,fEventWeight);
    
    if (fEvtPar.fNTracksDe> 0) FillClusterHistograms(fHist.fCluster[1],cl,fEventWeight);
    else                       FillClusterHistograms(fHist.fCluster[2],cl,fEventWeight);
    if (cl->Energy()    > 10.) FillClusterHistograms(fHist.fCluster[3],cl,fEventWeight);
    if (cl->Energy()    > 55.) FillClusterHistograms(fHist.fCluster[4],cl,fEventWeight);
    if (cl->Energy()    > 80.) FillClusterHistograms(fHist.fCluster[5],cl,fEventWeight);
    
    if      (id == 0         ) FillClusterHistograms(fHist.fCluster[6],cl,fEventWeight);
    else if (id == 1         ) FillClusterHistograms(fHist.fCluster[7],cl,fEventWeight);
    
    double cl_x   = cl->fX;
    double cl_y   = cl->fY;
    double cl_r   = sqrt(cl_x*cl_x+cl_y*cl_y);
    if(cl_r > 500.) FillClusterHistograms(fHist.fCluster[8],cl,fEventWeight);
    if(cl_r > 500. && cl->Energy() > 55.) FillClusterHistograms(fHist.fCluster[9],cl,fEventWeight);
    if(id == 0) {
      if(cl_r > 500.) FillClusterHistograms(fHist.fCluster[10],cl,fEventWeight);
      if(cl_r > 500. && cl->Energy() > 55.) FillClusterHistograms(fHist.fCluster[11],cl,fEventWeight);
    } 
  }
}

double TRMCAnaModule::RMCWeight(double energy) {
  if(energy < 0. || energy > fKinematicLimit) return 0.;
  double weight = 1.;
  if(fSpectrum == kClosure && !(fInternalRMC && fIntSpectrum != 0)) {
    weight = (energy > fKMax) ? 0. : TStntuple::RMC_ClosureAppxWeight(energy, fKMax);
  } else {    
    weight = fRMCSpectra->Weight(energy);
  }
  //   if(fSpectrum == kClosureFlat) {
  //   weight = (energy < fSpectrumParam[1]) ? std::max(fSpectrumParam[0]/(fSpectrumParam[1]-fKMax),  (float) TStntuple::RMC_ClosureAppxWeight(energy, fKMax)) : 0.;
  // } else if(fSpectrum == kClosureExp) {
  //   weight = (energy < fSpectrumParam[0]) ? TStntuple::RMC_ClosureAppxWeight(energy, fKMax) : 
  //     TStntuple::RMC_ClosureAppxWeight(fSpectrumParam[0], fKMax)*std::exp((energy - fSpectrumParam[0])*fSpectrumParam[1]); //match at intersection point
  // } else if(fSpectrum == kClosureTransition) {
  //   weight = (abs(energy - fSpectrumParam[1]) < 0.1) ? fSpectrumParam[0] : TStntuple::RMC_ClosureAppxWeight(energy, fKMax);
  // } else {
  //   std::cout << "ERROR! Unknown RMC spectrum being used: " << fSpectrum << std::endl;
  // }
  return weight;
}

void TRMCAnaModule::SetSpectrum(int spectrum) {
  fSpectrum = spectrum;
}
//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TRMCAnaModule::Event(int ientry) {

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

  fEvtPar.fNGenp    = fGenpBlock->NParticles();
  fEvtPar.fNStrawHits = GetHeaderBlock()->fNStrawHits;
  fEvtPar.fInstLum = GetHeaderBlock()->fInstLum;
  fEvtPar.fOneBatchWeight = BatchModeWeight(fEvtPar.fInstLum, 1); //1 batch mode
  fEvtPar.fTwoBatchWeight = BatchModeWeight(fEvtPar.fInstLum, 2); //2 batch mode
//-----------------------------------------------------------------------------
// MC generator info
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  int           pdg_code, generator_code;

  fEvtPar.fParticle = NULL;
  double k = fGenpBlock->GenEnergy();
  bool replaceK = false;
  if(k < 0.) {k = 0.; replaceK = true;} //< 0 means not saved
  for (int i=fEvtPar.fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if(replaceK && generator_code == fMCProcessCode) {
      TLorentzVector mom;
      genp->Momentum(mom);
      k += mom.P();
    }
    if ((abs(pdg_code) == fPDGCode) && (generator_code == fMCProcessCode)) {
      fEvtPar.fParticle = genp;
      if(!(fMCProcessCode == 42 || fMCProcessCode == 45)) break; //internal RMC/gammaPairProduction need to look at both particles
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
  fEleE = fSimPar.fEleE;
  // if(fSimPar.fGenp) fSimPar.fGenp->Print();
  // std::cout << "Event energy: " << fEleE << std::endl;
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

  fEvtPar.fNClusters   = fClusterBlock->NClusters();

  fEvtPar.fNTracksDe   = fTrackBlock->NTracks();
  fEvtPar.fNGoodTracks = 0;

  if (fEvtPar.fNTracksDe == 0) fTrack = 0;
  else                         fTrack = fTrackBlock->Track(0);


  for (int i=0; i<fEvtPar.fNTracksDe; i++) {
    TrackPar_t*   tp = fTrackPar+i;

    tp->fTrackID[0] = TAnaModule::fTrackID_BOX;
    tp->fTrackID[1] = TAnaModule::fTrackID_MVA;
    tp->fTrackID[2] = TRMCAnaModule::fTrackID_RMC_BOX;

    tp->fDioLOWt    = fEvtPar.fDioLOWt;
    tp->fDioLLWt    = fEvtPar.fDioLLWt;

    if (fTrackBlockName.Index("TrackBlockPar") == 0) tp->fFitType = 0;  // this better be made more error-proof
    if (fTrackBlockName.Index("TrackBlockDar") == 0) tp->fFitType = 1;
  }
  
  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,&fSimPar);
//-----------------------------------------------------------------------------
// RMC-specific initializations, assume running on a RMC dataset, generated 
// flat in the photon energy
//-----------------------------------------------------------------------------
  fEventWeight = (fFlatInput) ? RMCWeight(k) : 1.;
  fEventWeight *= fGenpBlock->Weight();
  for (int i=0; i<fEvtPar.fNTracksDe; i++) {
    fTrackPar[i].fRMCEnergyWt = (fFlatInput) ? RMCWeight(k) : 1.;
    fTrackPar[i].fTotWt       = fTrackPar[i].fRMCEnergyWt*fGenpBlock->Weight();
  }
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TRMCAnaModule::Debug() {

  TStnTrack* trk;
  Mu2eII::TrackPar_t*  tp;
  int ntrk = fTrackBlock->NTracks();
  if(GetDebugBit(0) == 1) {
    GetHeaderBlock()->Print(Form("Event weight = %.3e", fEventWeight));
  }

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock->Track(itrk);
    tp  = &fTrackPar[itrk];
//-----------------------------------------------------------------------------
// bit 1: All tracks
//-----------------------------------------------------------------------------
    if (GetDebugBit(1) == 1) {
      GetHeaderBlock()->Print(Form("Track: id_word[0] = %i (%i), fNTrkID = %i",tp->fIDWord[0],tp->fTrackID[0]->IDWord(trk), fNTrkID));
    }
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
int TRMCAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TRMCAnaModule::Test001() {
}
}
