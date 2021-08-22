//////////////////////////////////////////////////////////////////////////////
// CRY3 : all tracks have T0 > 700...
// use of debug bits:
//
// bit 003:  
// bit 004:  
// bit 005: CRV stubs sector==0 and x>-2000.
// bit 006: CRV pulses with charge < 9 PE
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
#include "Mu2eII/ana/TCrvAnaModule.hh"

// ClassImp(TCrvAnaModule)

namespace Mu2eII {

//-----------------------------------------------------------------------------
TCrvAnaModule::TCrvAnaModule(const char* name, const char* title):
  TAnaModule(name,title)
{
  //  fPtMin           = 1.;
  // fTrackNumber.Set(100);
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
//-----------------------------------------------------------------------------
  fPDGCode         = 11;
  fMCProcessCode   =  2;                  // conversionGun, 28:StoppedParticleReactionGun
  //  fBestID          = 0;                   // best ID word

  // fTrackBlockName  = "TrackBlockDar";
  // fNTrkID          = 2;                   // keep comparing the box and MVA cuts
}

//-----------------------------------------------------------------------------
TCrvAnaModule::~TCrvAnaModule() {
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TCrvAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  // RegisterDataBlock(fTrackBlockName.Data(), "TStnTrackBlock"   , &fTrackBlock     );
  // RegisterDataBlock("ClusterBlock"        , "TStnClusterBlock" , &fClusterBlock   );
  // RegisterDataBlock("TimeClusterBlockDe"  , "TStnTimeClusterBlock", &fTimeClusterBlock);
  // RegisterDataBlock("HelixBlock"          , "TStnHelixBlock"   , &fHelixBlock     );

  RegisterDataBlock("GenpBlock"           , "TGenpBlock"       , &fGenpBlock      );
  RegisterDataBlock("SimpBlock"           , "TSimpBlock"       , &fSimpBlock      );
  RegisterDataBlock("SpmcBlockVDet"       , "TStepPointMCBlock", &fSpmcBlockVDet  );
  RegisterDataBlock("CrvClusterBlock"     , "TCrvClusterBlock" , &fCrvClusterBlock);
  RegisterDataBlock("CrvPulseBlock"       , "TCrvPulseBlock"   , &fCrvPulseBlock  );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//_____________________________________________________________________________
void TCrvAnaModule::BookHistograms() {

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
  TString* event_selection[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) event_selection[i] = 0;

  event_selection[ 0] = new TString("all events");

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
  TString* simp_selection[kNSimpHistSets];
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
// book Genp histograms
//-----------------------------------------------------------------------------
  TString* genp_selection[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) genp_selection[i] = 0;

  genp_selection[0] = new TString("all particles");

  for (int i=0; i<kNGenpHistSets; i++) {
    if (genp_selection[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title    = genp_selection[i]->Data();
      if (! fol ) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fGenp[i]  = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book CRV pulse histograms
//-----------------------------------------------------------------------------
  TString* crvp_selection[kNCrvPulseHistSets];
  for (int i=0; i<kNCrvPulseHistSets; i++) crvp_selection[i] = 0;

  crvp_selection[0] = new TString("all");
  crvp_selection[1] = new TString("from clusters");

  for (int i=0; i<kNCrvPulseHistSets; i++) {
    if (crvp_selection[i] != 0) {
      sprintf(folder_name,"crvp_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title       = crvp_selection[i]->Data();
      fol                = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fCrvPulse[i] = new CrvPulseHist_t;
      BookCrvPulseHistograms(fHist.fCrvPulse[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book CRV cluster histograms
//-----------------------------------------------------------------------------
  TString* crvc_selection[kNCrvClusterHistSets];
  for (int i=0; i<kNCrvClusterHistSets; i++) crvc_selection[i] = nullptr;

  crvc_selection[  0]   = new TString("all particles");
  crvc_selection[100] = new TString("CRV sector= 0");
  crvc_selection[101] = new TString("CRV sector= 1");
  crvc_selection[102] = new TString("CRV sector= 2");
  crvc_selection[103] = new TString("CRV sector= 3");
  crvc_selection[104] = new TString("CRV sector= 4");
  crvc_selection[105] = new TString("CRV sector= 5");
  crvc_selection[106] = new TString("CRV sector= 6");
  crvc_selection[107] = new TString("CRV sector= 7");
  crvc_selection[108] = new TString("CRV sector= 8");
  crvc_selection[109] = new TString("CRV sector= 9");
  crvc_selection[110] = new TString("CRV sector=10");
  crvc_selection[111] = new TString("CRV sector=11");
  crvc_selection[112] = new TString("CRV sector=12");
  crvc_selection[113] = new TString("CRV sector=13");
  crvc_selection[114] = new TString("CRV sector=14");
  crvc_selection[115] = new TString("CRV sector=15");
  crvc_selection[116] = new TString("CRV sector=16");
  crvc_selection[117] = new TString("CRV sector=17");
  crvc_selection[118] = new TString("CRV sector=18");
  crvc_selection[119] = new TString("CRV sector=19");
  crvc_selection[120] = new TString("CRV sector=20");
  crvc_selection[121] = new TString("CRV sector=21");

  for (int i=0; i<kNCrvClusterHistSets; i++) {
    if (crvc_selection[i] != 0) {
      sprintf(folder_name,"crvc_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title         = crvc_selection[i]->Data();
      if (! fol) fol       = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fCrvCluster[i] = new CrvClusterHist_t;
      BookCrvClusterHistograms(fHist.fCrvCluster[i],Form("Hist/%s",folder_name));
    }
  }
}

//_____________________________________________________________________________
int TCrvAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}

//_____________________________________________________________________________
void TCrvAnaModule::FillHistograms() {

  double wt_b1(fEventWeight), wt_b2(fEventWeight);

  if(fBatchMode == 2) wt_b1 *= fEvtPar.fOneBatchWeight / fEvtPar.fTwoBatchWeight;
  if(fBatchMode == 1) wt_b2 *= fEvtPar.fTwoBatchWeight / fEvtPar.fOneBatchWeight;
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
// 4. CRV pulse histograms
//  CRVP_0 : all pulses
//-----------------------------------------------------------------------------
  int ncrvp = fCrvPulseBlock->NPulses();

  for (int i=0; i<ncrvp; i++) {
    TCrvRecoPulse* crp = fCrvPulseBlock->Pulse(i);
    FillCrvPulseHistograms(fHist.fCrvPulse[0],crp);
  }
//-----------------------------------------------------------------------------
// 5. CRV cluster histograms
// CRVC_0: all clusters
//-----------------------------------------------------------------------------
  int nstubs = fCrvClusterBlock->NClusters();

  for (int is=0; is<nstubs; is++) {
    TCrvCoincidenceCluster* stub = fCrvClusterBlock->Cluster(is);
    CrvStubPar_t* sp = fCrvStubPar+is;

    FillCrvClusterHistograms(fHist.fCrvCluster[0],stub,sp);

    int sector = sp->fSector;
    FillCrvClusterHistograms(fHist.fCrvCluster[100+sector],stub,sp);

    int np = fCrvClusterBlock->NClusterPulses(is);
    for (int ip=0; ip<np; ip++) {
      int loc = fCrvClusterBlock->ClusterPulseIndex(is,ip);
      TCrvRecoPulse* crp = fCrvClusterBlock->Pulse(loc);
//-----------------------------------------------------------------------------
//  CRVP_1 : pulses from clusters
//-----------------------------------------------------------------------------
      FillCrvPulseHistograms(fHist.fCrvPulse[1],crp);
    }
  }
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TCrvAnaModule::Event(int ientry) {

  fEvtPar.fCutCounter[0] = 0;
  fEvtPar.fCutCounter[1] = 0;
  fEvtPar.fCutCounter[2] = 0;
  fEvtPar.fCutCounter[3] = 0;

  TLorentzVector        mom;

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fSpmcBlockVDet->GetEntry(ientry);

  fCrvClusterBlock->GetEntry(ientry);
  fCrvPulseBlock->GetEntry(ientry);

  fEventWeight              = 1.;
  fEvtPar.fDioLOWt          = 1.;
  fEvtPar.fDioLLWt          = 1.;

  fEvtPar.fNCrvClusters     = fCrvClusterBlock->NClusters();
  fEvtPar.fNCrvPulses       = fCrvPulseBlock->NPulses();
  fEvtPar.fNCrvCoincidences = fCrvPulseBlock->NCoincidences();
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fEvtPar.fNGenp          = fGenpBlock->NParticles();
  fEvtPar.fNStrawHits     = GetHeaderBlock()->fNStrawHits;
  fEvtPar.fInstLum        = GetHeaderBlock()->fInstLum;
  fEvtPar.fOneBatchWeight = BatchModeWeight(fEvtPar.fInstLum, 1); //1 batch mode
  fEvtPar.fTwoBatchWeight = BatchModeWeight(fEvtPar.fInstLum, 2); //2 batch mode
//-----------------------------------------------------------------------------
// MC generator info
//-----------------------------------------------------------------------------
  fEvtPar.fParticle = NULL;
  for (int i=fEvtPar.fNGenp-1; i>=0; i--) {
    TGenParticle* genp = fGenpBlock->Particle(i);
    int pdg_code       = genp->GetPdgCode();
    int generator_code = genp->GetStatusCode();
    if ((abs(pdg_code) == fPDGCode) && (generator_code == fMCProcessCode)) {
      fEvtPar.fParticle = genp;
      break;
    }
  }
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimp and remove overlaps with fEvtPar
//-----------------------------------------------------------------------------
  fSimPar.fParticle = fSimpBlock->Particle(0);
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
  
  if (fSimPar.fParticle) {
    for (int i=0; i<nsteps; i++) {
      TStepPointMC* step = fSpmcBlockVDet->StepPointMC(i);
      if (step->PDGCode() == fSimPar.fParticle->fPdgCode) {
	if ((step->VolumeID() == 13) || (step->VolumeID() == 14)) {
	  fSimPar.fTFront = step;
	}
	else if ((step->VolumeID() == 11) || (step->VolumeID() == 12)) {
	  fSimPar.fTMid = step;
	}
      }
    }
  }

  fEvtPar.fNTracksDe   = 0;
  fEvtPar.fNGoodTracks = 0;

  InitCrvStubPar(fCrvClusterBlock, fCrvStubPar);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TCrvAnaModule::Debug() {

  if (GetDebugBit(3) == 1) {
    if (fEvtPar.fNCrvClusters > 0) { 
      GetHeaderBlock()->Print(Form("bit_003: N(CRV clusters): %3i",fEvtPar.fNCrvClusters));
    }
  }

  for (int i=0; i<fEvtPar.fNCrvClusters; i++) {
    TCrvCoincidenceCluster* crvcc = fCrvClusterBlock->Cluster(i);
    CrvStubPar_t* sp              = fCrvStubPar+i;
    const TVector3* pos           = crvcc->Position();

    if (GetDebugBit(4) == 1) {
      if ((fabs(pos->X()) < 1000) and (pos->Z() > 0) and (pos->Z() < 400)) {
	GetHeaderBlock()->Print(Form("bit_004: CRV cluster X,Y,Z = %11.4e %11.4e %11.4e",pos->X(),pos->Y(),pos->Z()));
      }
    }

    if (GetDebugBit(5) == 1) {
      if ((sp->fSector == 0) and (pos->X() > -2000.)) {
	GetHeaderBlock()->Print(Form("bit_005: CRV cluster sector == 0 and x = %10.3e",pos->X()));
      }
    }
  }

  for (int i=0; i<fEvtPar.fNCrvPulses; i++) {
    TCrvRecoPulse* p = fCrvPulseBlock->Pulse(i);
    if (GetDebugBit(6) == 1) {
      if (p->NPe() < 9) {
	GetHeaderBlock()->Print(Form("bit_006: CRV pulse with N(PE) = %3i",p->NPe()));
      }
    }
  }
}

//_____________________________________________________________________________
int TCrvAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TCrvAnaModule::Test001() {
}
}
