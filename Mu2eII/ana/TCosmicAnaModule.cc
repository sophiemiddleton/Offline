//////////////////////////////////////////////////////////////////////////////
// CRY3 : all tracks have T0 > 700...
// use of debug bits:
//
// bit   3: print events with good tracks and no CRV coincidences (not even clusters)
// bit   4: print events with T(track) - T(CRV start) < -100
// bit   5: print events with high dT, print NCrvStubs, t0(trk), t(startCRV), dPf
// bit   6: print events with dpf < -2, print t0 of track
// bit   7 and 8 are also in use
// bit   9: print events with dT > +100, print NCrvStubs, t0(trk), t(startCRV), dPf
// bit  10: print events with more than 2 time clusters
// bit  11: used 
// bit  12: used 
// bit  13: used 
// bit  14: used 
// bit  15: figure out structures in fDtCRV...
// bit  18: potentially, events electrons/positrons and lost downstream leg
// bit  19: lower leg of Ralf's V3059
// bit  20: fNCrvClusters == 0, tp->fIDWord[1]=0, pid_ele
// bit  21: upper leg of Ralf's V3059 constellation
// bit  22: print events with good tracks and no CRV stub candidates (evt_4)
// bit_023: upper tail of the DT corrected distribution - track of any sign
// bit_024: CRV background , region 1
// bit_025: CRV background , region 2
// bit_026: CRV background , region 3
// bit_027: CRV background , region 4
// bit_028: CRV background , region 5
// bit_029: CRV background , all regions, dt > 100, or dt < -40
// bit_035: fEvtPar.fCandidate_MVA
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
#include "Mu2eII/ana/TCosmicAnaModule.hh"

// ClassImp(TCosmicAnaModule)

namespace Mu2eII {

//-----------------------------------------------------------------------------
TCosmicAnaModule::TCosmicAnaModule(const char* name, const char* title):
  TAnaModule(name,title)
{
  fPtMin           = 1.;
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
//-----------------------------------------------------------------------------
  fPDGCode         = 11;
  fMCProcessCode   =  2;                  // conversionGun, 28:StoppedParticleReactionGun
  //  fBestID          = 0;                   // best ID word

  fTrackBlockNameDe = "TrackBlockDarDe";
  fTrackBlockNameUe = "TrackBlockDarUe";
  fNTrkID           = 2;                   // keep comparing the box and MVA cuts
//-----------------------------------------------------------------------------
// assuming a track has an associated cluster , the fit error on the reconstructed T0 
// should be small, require T0Err < 0.9 ns
// if not, the cluster has been dropped - don't want that
//-----------------------------------------------------------------------------
  int mask = TStnTrackID::kTrkQualBit | TStnTrackID::kD0Bit | TStnTrackID::kTanDipBit | TStnTrackID::kT0Bit;
  mask     = mask | TStnTrackID::kT0ErrBit ;
  fTrackID_MVA->SetUseMask(mask);
  fTrackID_MVA->SetMaxT0Err(0.9);
}

//-----------------------------------------------------------------------------
TCosmicAnaModule::~TCosmicAnaModule() {
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TCosmicAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockNameDe.Data(), "TStnTrackBlock"      , &fTrackBlockDe    );
  RegisterDataBlock(fTrackBlockNameUe.Data(), "TStnTrackBlockUe"    , &fTrackBlockUe    );
  RegisterDataBlock("ClusterBlock"          , "TStnClusterBlock"    , &fClusterBlock    );

  RegisterDataBlock("TimeClusterBlockDe"    , "TStnTimeClusterBlock", &fTimeClusterBlockDe);
  RegisterDataBlock("TimeClusterBlockUe"    , "TStnTimeClusterBlock", &fTimeClusterBlockUe);

  RegisterDataBlock("TCFinderBlockDe"       , "TStnTimeClusterBlock", &fTCFinderBlockDe);
  RegisterDataBlock("TCFinderBlockUe"       , "TStnTimeClusterBlock", &fTCFinderBlockUe);
  RegisterDataBlock("CTPFinderBlock"        , "TStnTimeClusterBlock", &fCTPFinderBlock );

  RegisterDataBlock("GenpBlock"             , "TGenpBlock"          , &fGenpBlock       );
  RegisterDataBlock("SimpBlock"             , "TSimpBlock"          , &fSimpBlock       );
  RegisterDataBlock("SpmcBlockVDet"         , "TStepPointMCBlock"   , &fSpmcBlockVDet   );
  RegisterDataBlock("CrvClusterBlock"       , "TCrvClusterBlock"    , &fCrvClusterBlock );
  RegisterDataBlock("CrvPulseBlock"         , "TCrvPulseBlock"      , &fCrvPulseBlock   );

  RegisterDataBlock("HelixBlockDe"          , "TStnHelixBlock"      , &fHelixBlockDe    );
  RegisterDataBlock("HelixBlockUe"          , "TStnHelixBlock"      , &fHelixBlockUe    );
  RegisterDataBlock("TrackSeedBlockDe"      , "TStnTrackSeedBlock"  , &fTrackSeedBlockDe);
  RegisterDataBlock("TrackSeedBlockUe"      , "TStnTrackSeedBlock"  , &fTrackSeedBlockUe);
//-----------------------------------------------------------------------------
// cache pointers to non-CRV data blocks needed for cosmic rejection, do that once
//-----------------------------------------------------------------------------
  fCosmicVetoData.fTrackBlockDe    = fTrackBlockDe;
  fCosmicVetoData.fTrackParDe      = fTrackParDe;
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
    TrackPar_t* tp = &fTrackParDe[i];

    tp->fFitType     = 1;                                // assume DAR, PAR needs to be set specially
    tp->fTrqMvaIndex = 0;                                // index of the TRQ MVA used by this block

    tp->fTrackID[0]  = fTrackID_BOX;                     // these poiters need to be set just once
    tp->fTrackID[1]  = fTrackID_MVA;

    tp = &fTrackParUe[i];

    tp->fFitType     = 1;                                // assume DAR, PAR needs to be set specially
    tp->fTrqMvaIndex = 0;                                // index of the TRQ MVA used by this block

    tp->fTrackID[0]  = fTrackID_BOX;                     // these poiters need to be set just once
    tp->fTrackID[1]  = fTrackID_MVA;
  }

  return 0;
}

//_____________________________________________________________________________
void TCosmicAnaModule::BookHistograms() {

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
  event_selection[ 1] = new TString("events with a reconstructed track");
  event_selection[ 2] = new TString("events with a reconstructed track passing BOX cuts");
  event_selection[ 3] = new TString("events with a reconstructed track passing MVA cut");
  event_selection[ 4] = new TString("events with a reconstructed track passing BOX cut and no clusters");
  event_selection[ 5] = new TString("event candidates, no momentum cut");
  event_selection[ 6] = new TString("event candidates outside of CRV dead time window, but good timing");

  event_selection[11] = new TString("events passing the cosmic veto cut");
  event_selection[12] = new TString("events failing the cosmic veto cut");

  event_selection[13] = new TString("candidate events, BOX cuts");
  event_selection[14] = new TString("events failing the candidate BOX cuts");

  event_selection[15] = new TString("candidate events, MVA cuts");
  event_selection[16] = new TString("events failing the candidate MVA cuts");

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
// book track histograms
//-----------------------------------------------------------------------------
  TString* track_selection[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) track_selection[i] = 0;

  track_selection[   0] = new TString("all tracks");

  track_selection[  10] = new TString("all e- tracks");
  track_selection[  11] = new TString("all e- tracks   0 <  Ep < 0.6");
  track_selection[  12] = new TString("all e- tracks 0.6 <= Ep < 1.05");
  track_selection[  13] = new TString("all e- tracks 1.05 <= Ep");

  track_selection[  20] = new TString("all e+ tracks");
  track_selection[  21] = new TString("all e+ tracks   0 <= Ep < 0.6");
  track_selection[  22] = new TString("all e+ tracks 0.6 <= Ep < 1.0.5");
  track_selection[  23] = new TString("all e+ tracks 1.05 <= Ep");

  track_selection[ 100] = new TString("all tracks passing BOX cuts (Michael)");

  track_selection[ 110] = new TString("all e- tracks passing box cuts");
  track_selection[ 111] = new TString("all e- tracks passing box cuts   0 <  Ep < 0.6");
  track_selection[ 112] = new TString("all e- tracks passing box cuts 0.6 <= Ep < 1.05");
  track_selection[ 113] = new TString("all e- tracks passing box cuts 1.05 <= Ep");
  track_selection[ 114] = new TString("all e- tracks passing box cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0");
  track_selection[ 115] = new TString("all e- tracks passing box cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 and NCRVCl == 0");
  track_selection[ 116] = new TString("all e- tracks passing box cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 and 100 <= P < 110");

  track_selection[ 120] = new TString("all e+ tracks passing box cuts");
  track_selection[ 121] = new TString("all e+ tracks passing box cuts   0  <  Ep < 0.6");
  track_selection[ 122] = new TString("all e+ tracks passing box cuts 0.6  <= Ep < 1.05");
  track_selection[ 123] = new TString("all e+ tracks passing box cuts 1.05 <= Ep");
  track_selection[ 124] = new TString("all e+ tracks passing box cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0");
  track_selection[ 125] = new TString("all e+ tracks passing box cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 and NCRVCl == 0");
  track_selection[ 126] = new TString("all e+ tracks passing box cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 and 100 <= P < 110");

  track_selection[ 200] = new TString("all tracks passing TRQ MVA cuts");
  track_selection[ 201] = new TString("all tracks passing TRQ MVA cuts, TCVeto==0");
  track_selection[ 202] = new TString("all tracks passing TRQ MVA cuts, TCVeto>=1");
  track_selection[ 203] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid");
  track_selection[ 204] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ~ele_pid");
  track_selection[ 205] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ep>1.05");
  track_selection[ 206] = new TString("all tracks passing TRQ MVA, PIDe");
  track_selection[ 207] = new TString("all tracks passing TRQ MVA, PIDe, non_crv_veto = 0");
  track_selection[ 208] = new TString("all tracks passing TRQ MVA, PIDe, NCRVStubs>0,crv_veto=0");

  track_selection[ 209] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, NTracksDe<2, NTracksUe<2, at least 1 U-TC, dt(track-timeClusterUe) not in [50,200]");
  track_selection[ 210] = new TString("all e- tracks passing TRQ MVA cuts");
  track_selection[ 211] = new TString("all e- tracks passing TRQ MVA cuts   0  <  Ep < 0.6");
  track_selection[ 212] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05");
  track_selection[ 213] = new TString("all e- tracks passing TRQ MVA cuts 1.05 <= Ep");
  track_selection[ 214] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0");
  track_selection[ 215] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 and NCRVCl == 0");
  track_selection[ 216] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 and 100 <= P < 110");
  track_selection[ 217] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 with 1 time cluster");
  track_selection[ 218] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 with more than 1 time cluster");

  // duplicate of trk_214
  // track_selection[ 219] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 with corrected dT");

  track_selection[ 220] = new TString("all e+ tracks passing TRQ MVA cuts");
  track_selection[ 221] = new TString("all e+ tracks passing TRQ MVA cuts   0  <  Ep < 0.6");
  track_selection[ 222] = new TString("all e+ tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05");
  track_selection[ 223] = new TString("all e+ tracks passing TRQ MVA cuts 1.05 <= Ep");
  track_selection[ 224] = new TString("all e+ tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0");

  track_selection[ 225] = new TString("all e+ tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 and NCRVCl == 0");
  track_selection[ 226] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 and 100 <= P < 110");

  track_selection[ 240] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0, with PDG e-");
  track_selection[ 241] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0, with PDG e+");
  track_selection[ 242] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0, with PDG #mu-");
  track_selection[ 243] = new TString("all e- tracks passing TRQ MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0, with PDG #mu+");

  track_selection[ 262] = new TString("all e- tracks passing MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 with 1 effective time cluster, corrected");
  track_selection[ 263] = new TString("all e- tracks passing MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 with more than 1 effective time cluster, corrected");
  track_selection[ 264] = new TString("all e- tracks passing MVA cuts 0.6  <= Ep < 1.05 and -3.5 < dt < 1.0 outside of CRV window, but good timing, corrected");
  // these are the folders of track histograms filled by CRV sector number
  track_selection[ 265] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 0");
  track_selection[ 266] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 1");
  track_selection[ 267] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 2");
  track_selection[ 268] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 3");
  track_selection[ 269] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 4");
  track_selection[ 270] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 5");
  track_selection[ 271] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 6");
  track_selection[ 272] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 7");
  track_selection[ 273] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 8");
  track_selection[ 274] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 9");
  track_selection[ 275] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 10");
  track_selection[ 276] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 11");
  track_selection[ 277] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 12");
  track_selection[ 278] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 13");
  track_selection[ 279] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 14");
  track_selection[ 280] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 15");
  track_selection[ 281] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 16");
  track_selection[ 282] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 17");
  track_selection[ 283] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 18");
  track_selection[ 284] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 19");
  track_selection[ 285] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 20");
  track_selection[ 286] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV sector = 21");
  track_selection[ 287] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV stub slope > 0");
  track_selection[ 288] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV stub slope < 0");
  track_selection[ 289] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, in sector 0 with CRV stub slope > 0");
  track_selection[ 290] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, in sector 0 with CRV stub slope < 0");
  track_selection[ 291] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, in sector 10 with CRV stub slope > 0");
  track_selection[ 292] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, in sector 10 with CRV stub slope < 0");
  track_selection[ 293] = new TString("all tracks passing TRQ MVA cuts, CrvVeto, ele_pid, with CRV stub slope * slopeMC product < 0");
  track_selection[ 301] = new TString("all tracks passing MVA cuts dTCRV < 50 ZCRV>12000");  // investigation, so far - unsuccessfull P.M.
  track_selection[ 302] = new TString("all tracks passing MVA cuts dTCRV > 50 ZCRV>12000");  // investigation, so far - unsuccessfull  P.M.

  //upstream tracks
  track_selection[ 500] = new TString("all upstream tracks");
  track_selection[ 501] = new TString("all upstream tracks passing MVA cuts");
  track_selection[ 502] = new TString("all upstream tracks NOT passing MVA cuts");

  // electrons
  track_selection[1000] = new TString("e- tracks passing BOX TRQ + PID");
  track_selection[1004] = new TString("e- tracks passing BOX TRQ + PID batch 1 wt");
  track_selection[1005] = new TString("e- tracks passing BOX TRQ + PID batch 2 wt");
  track_selection[1006] = new TString("e- tracks passing BOX TRQ + PID, no veto batch 1 wt");
  track_selection[1007] = new TString("e- tracks passing BOX TRQ + PID, no veto batch 2 wt");
  track_selection[1010] = new TString("e- tracks passing BOX TRQ + PID, T>700");
  track_selection[1014] = new TString("e- tracks passing BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[1015] = new TString("e- tracks passing BOX TRQ + PID, T>700 batch 2 wt");
  track_selection[1016] = new TString("e- tracks passing BOX TRQ + PID, T>700, no veto batch 1 wt");
  track_selection[1017] = new TString("e- tracks passing BOX TRQ + PID, T>700, no veto batch 2 wt");
  track_selection[1020] = new TString("e- tracks passing BOX TRQ + N(CRVStubs=0)");
  track_selection[1021] = new TString("e- tracks passing BOX TRQ + N(CRVStubs=0) + PIDe");
  track_selection[1022] = new TString("e- tracks passing BOX TRQ + N(CRVStubs=0) + PIDm");
  track_selection[1023] = new TString("e- tracks passing BOX TRQ + N(CRVStubs=0) + T>700");
  track_selection[1024] = new TString("e- tracks passing BOX TRQ + N(CRVStubs=0) + T>700 + PIDe");
  track_selection[1025] = new TString("e- tracks passing BOX TRQ + N(CRVStubs=0) + T>700 + PIDm");
  track_selection[1026] = new TString("e- tracks passing MVA TRQ + N(CRVStubs=0) + T>700 + PIDe+ 100<P<110");

  track_selection[2000] = new TString("e- tracks passing MVA TRQ + PID");
  track_selection[2004] = new TString("e- tracks passing MVA TRQ + PID batch 1 wt");
  track_selection[2005] = new TString("e- tracks passing MVA TRQ + PID batch 2 wt");
  track_selection[2006] = new TString("e- tracks passing MVA TRQ + PID, no veto batch 1 wt");
  track_selection[2007] = new TString("e- tracks passing MVA TRQ + PID, no veto batch 2 wt");
  track_selection[2010] = new TString("e- tracks passing MVA TRQ + PID, T>700");
  track_selection[2014] = new TString("e- tracks passing MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[2015] = new TString("e- tracks passing MVA TRQ + PID, T>700 batch 2 wt");
  track_selection[2016] = new TString("e- tracks passing MVA TRQ + PID, T>700, no veto batch 1 wt");
  track_selection[2017] = new TString("e- tracks passing MVA TRQ + PID, T>700, no veto batch 2 wt");
  track_selection[2020] = new TString("e- tracks passing MVA TRQ + N(CRVStubs=0)");
  track_selection[2021] = new TString("e- tracks passing MVA TRQ + N(CRVStubs=0) + PIDe");
  track_selection[2022] = new TString("e- tracks passing MVA TRQ + N(CRVStubs=0) + PIDm");
  track_selection[2023] = new TString("e- tracks passing MVA TRQ + N(CRVStubs=0) + T>700");
  track_selection[2024] = new TString("e- tracks passing MVA TRQ + N(CRVStubs=0) + T>700 + PIDe");
  track_selection[2025] = new TString("e- tracks passing MVA TRQ + N(CRVStubs=0) + T>700 + PIDm");
  track_selection[2026] = new TString("e- tracks passing MVA TRQ + N(CRVStubs=0) + T>700 + PIDe+ 100<P<110");

  track_selection[3000] = new TString("e+ tracks passing BOX TRQ + PID");
  track_selection[3004] = new TString("e+ tracks passing BOX TRQ + PID batch 1 wt");
  track_selection[3005] = new TString("e+ tracks passing BOX TRQ + PID batch 2 wt");
  track_selection[3006] = new TString("e+ tracks passing BOX TRQ + PID, no veto batch 1 wt");
  track_selection[3007] = new TString("e+ tracks passing BOX TRQ + PID, no veto batch 2 wt");
  track_selection[3010] = new TString("e+ tracks passing BOX TRQ + PID, T>700");
  track_selection[3014] = new TString("e+ tracks passing BOX TRQ + PID, T>700 batch 1 wt");
  track_selection[3015] = new TString("e+ tracks passing BOX TRQ + PID, T>700 batch 2 wt");
  track_selection[3016] = new TString("e+ tracks passing BOX TRQ + PID, T>700, no veto batch 1 wt");
  track_selection[3017] = new TString("e+ tracks passing BOX TRQ + PID, T>700, no veto batch 2 wt");
  track_selection[3020] = new TString("e+ tracks passing BOX TRQ + PID N(CRVStubs=0)");
  track_selection[3021] = new TString("e+ tracks passing BOX TRQ + N(CRVStubs=0) + PIDe");
  track_selection[3022] = new TString("e+ tracks passing BOX TRQ + N(CRVStubs=0) + PIDm");
  track_selection[3023] = new TString("e+ tracks passing BOX TRQ + N(CRVStubs=0) + T>700");
  track_selection[3024] = new TString("e+ tracks passing BOX TRQ + N(CRVStubs=0) + T>700 + PIDe");
  track_selection[3025] = new TString("e+ tracks passing BOX TRQ + N(CRVStubs=0) + T>700 + PIDm");
  track_selection[3026] = new TString("e+ tracks passing BOX TRQ + N(CRVStubs=0) + T>700 + PIDe+ 100<P<110");

  track_selection[4000] = new TString("e+ tracks passing MVA TRQ + PID");
  track_selection[4004] = new TString("e+ tracks passing MVA TRQ + PID batch 1 wt");
  track_selection[4005] = new TString("e+ tracks passing MVA TRQ + PID batch 2 wt");
  track_selection[4006] = new TString("e+ tracks passing MVA TRQ + PID, no veto batch 1 wt");
  track_selection[4007] = new TString("e+ tracks passing MVA TRQ + PID, no veto batch 2 wt");
  track_selection[4010] = new TString("e+ tracks passing MVA TRQ + PID, T>700");
  track_selection[4014] = new TString("e+ tracks passing MVA TRQ + PID, T>700 batch 1 wt");
  track_selection[4015] = new TString("e+ tracks passing MVA TRQ + PID, T>700 batch 2 wt");
  track_selection[4016] = new TString("e+ tracks passing MVA TRQ + PID, T>700, no veto batch 1 wt");
  track_selection[4017] = new TString("e+ tracks passing MVA TRQ + PID, T>700, no veto batch 2 wt");
  track_selection[4020] = new TString("e+ tracks passing MVA TRQ + PID N(CRVStubs=0)");
  track_selection[4021] = new TString("e+ tracks passing MVA TRQ + N(CRVStubs=0) + PIDe");
  track_selection[4022] = new TString("e+ tracks passing MVA TRQ + N(CRVStubs=0) + PIDm");
  track_selection[4023] = new TString("e+ tracks passing MVA TRQ + N(CRVStubs=0) + T>700");
  track_selection[4024] = new TString("e+ tracks passing MVA TRQ + N(CRVStubs=0) + T>700 + PIDe");
  track_selection[4025] = new TString("e+ tracks passing MVA TRQ + N(CRVStubs=0) + T>700 + PIDm");
  track_selection[4026] = new TString("e+ tracks passing MVA TRQ + N(CRVStubs=0) + T>700 + PIDe+ 100<P<110");

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
// book track-time cluster histograms
//-----------------------------------------------------------------------------
  TString* track_tc_selection[kNTrackTcHistSets];
  for (int i=0; i<kNTrackTcHistSets; i++) track_tc_selection[i] = 0;

  track_tc_selection[0] = new TString("all tracks");
  track_tc_selection[1] = new TString("tracks fTrackID[0]=0");

  for (int i=0; i<kNTrackTcHistSets; i++) {
    if (track_tc_selection[i] != 0) {
      sprintf(folder_name,"ttc_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title    = track_tc_selection[i]->Data();
      if (! fol ) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fTrackTc[i]  = new TrackTcHist_t;
      BookTrackTcHistograms(fHist.fTrackTc[i],Form("Hist/%s",folder_name));
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

  crvp_selection[0] = new TString("all pulses");

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

  crvc_selection[0] = new TString("all particles");
  crvc_selection[1] = new TString("sector 10");
  crvc_selection[2] = new TString("stub dYdZ.qn() = 2");
  crvc_selection[3] = new TString("stub dYdZ.qn() = 4");
  crvc_selection[4] = new TString("stub slope * slope MC product < 0");

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
//-----------------------------------------------------------------------------
// book track-CRV cluster histograms
//-----------------------------------------------------------------------------
  TString* trk_crvst_selection[kNTrackCrvStHistSets];
  for (int i=0; i<kNTrackCrvStHistSets; i++) trk_crvst_selection[i] = nullptr;

  trk_crvst_selection[0] = new TString("all TID De tracks");
  trk_crvst_selection[1] = new TString("all TID+PIDe De tracks");
  trk_crvst_selection[2] = new TString("all TID+PIDm De tracks");

  for (int i=0; i<kNTrackCrvStHistSets; i++) {
    if (trk_crvst_selection[i] != 0) {
      sprintf(folder_name,"tcrs_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title         = trk_crvst_selection[i]->Data();
      if (! fol) fol       = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fTrackCrvSt[i] = new TrackCrvStHist_t;
      BookTrackCrvStHistograms(fHist.fTrackCrvSt[i],Form("Hist/%s",folder_name));
    }
  }
}

//_____________________________________________________________________________
int TCosmicAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}

//_____________________________________________________________________________
void TCosmicAnaModule::FillHistograms() {

  double wt_b1(fEventWeight), wt_b2(fEventWeight);

  if(fBatchMode == 2) wt_b1 *= fEvtPar.fOneBatchWeight / fEvtPar.fTwoBatchWeight;
  if(fBatchMode == 1) wt_b2 *= fEvtPar.fTwoBatchWeight / fEvtPar.fOneBatchWeight;
//-----------------------------------------------------------------------------
// 1. fill event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);

  if (fEvtPar.fNTracksDe > 0) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
  if (fNGoodTracks_BOX   > 0) FillEventHistograms(fHist.fEvent[2],&fEvtPar);
  if (fNGoodTracks_MVA   > 0) FillEventHistograms(fHist.fEvent[3],&fEvtPar);

  if ((fNGoodTracks_BOX > 0) and (fEvtPar.fNCrvClusters == 0)) {
     FillEventHistograms(fHist.fEvent[4],&fEvtPar);
  }

  if (fEvtPar.fCandidate_BOX    ) FillEventHistograms(fHist.fEvent[5],&fEvtPar);
  if (fEvtPar.fCutCounter[3] > 0) FillEventHistograms(fHist.fEvent[6],&fEvtPar);
//-----------------------------------------------------------------------------
// EVT_11 : events passing the cosmic veto
// EVT_12 : events failing the cosmic veto
//-----------------------------------------------------------------------------
  if (fEvtPar.fCosmicVeto   == 0) FillEventHistograms(fHist.fEvent[11],&fEvtPar);
  else                            FillEventHistograms(fHist.fEvent[12],&fEvtPar);
//-----------------------------------------------------------------------------
// EVT_13 : candidate events, BOX cuts
// EVT_14 : the rest
// EVT_15 : candidate events, MVA cuts
// EVT_16 : the rest
//-----------------------------------------------------------------------------
  if (fEvtPar.fCandidate_BOX == 1) FillEventHistograms(fHist.fEvent[13],&fEvtPar);
  else                             FillEventHistograms(fHist.fEvent[14],&fEvtPar);

  if (fEvtPar.fCandidate_MVA == 1) FillEventHistograms(fHist.fEvent[15],&fEvtPar);
  else                             FillEventHistograms(fHist.fEvent[16],&fEvtPar);
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
    TStnTrack* trk         = fTrackBlockDe->Track(i);
    Mu2eII::TrackPar_t* tp = fTrackParDe+i;

    double p     = trk->P();

    bool low_ep  = (tp->fEp >  0  ) and (tp->fEp <  0.6);
    bool good_ep = (tp->fEp >= 0.6) and (tp->fEp < 1.05);
    bool high_ep = (tp->fEp >= 1.05);

    bool good_dt = (tp->fDt >=-3.5) and (tp->fDt < 1.0 );

    bool pid_ele = (tp->fPidMvaOut[0] >= 0.5);            // includes the pre-PID cuts
//-----------------------------------------------------------------------------
// all tracks
//-----------------------------------------------------------------------------
    FillTrackHistograms(fHist.fTrack[  0],trk,tp,&fSimPar);

    if (trk->Charge() < 0) { 
      FillTrackHistograms(fHist.fTrack[  10],trk,tp,&fSimPar);
      if (low_ep ) FillTrackHistograms(fHist.fTrack[  11],trk,tp,&fSimPar);
      if (good_ep) FillTrackHistograms(fHist.fTrack[  12],trk,tp,&fSimPar);
      if (high_ep) FillTrackHistograms(fHist.fTrack[  13],trk,tp,&fSimPar);
    }
    if (trk->Charge() > 0) {
      FillTrackHistograms(fHist.fTrack[  20],trk,tp,&fSimPar);
      if (low_ep ) FillTrackHistograms(fHist.fTrack[  21],trk,tp,&fSimPar);
      if (good_ep) FillTrackHistograms(fHist.fTrack[  22],trk,tp,&fSimPar);
      if (high_ep) FillTrackHistograms(fHist.fTrack[  23],trk,tp,&fSimPar);
    }
//-----------------------------------------------------------------------------
// "good" tracks, BOX cuts
//-----------------------------------------------------------------------------
    if (tp->fIDWord[0] == 0) {
      FillTrackHistograms(fHist.fTrack[100],trk,tp,&fSimPar);
      if (trk->Charge() < 0) {

	FillTrackHistograms(fHist.fTrack[110],trk,tp,&fSimPar);

	if (low_ep ) FillTrackHistograms(fHist.fTrack[ 111],trk,tp,&fSimPar);
	if (good_ep) FillTrackHistograms(fHist.fTrack[ 112],trk,tp,&fSimPar);

	if (high_ep) FillTrackHistograms(fHist.fTrack[ 113],trk,tp,&fSimPar);

	if (good_ep and good_dt) {
	  FillTrackHistograms(fHist.fTrack[ 114],trk,tp,&fSimPar);
	  if (fEvtPar.fNCrvClusters == 0) {
	    FillTrackHistograms(fHist.fTrack[ 115],trk,tp,&fSimPar);
	  }

	  if ((p >= 100) and (p < 110)) FillTrackHistograms(fHist.fTrack[ 116],trk,tp,&fSimPar);
	}

	if (fEvtPar.fCandidate_BOX) { 
	  FillTrackHistograms(fHist.fTrack[1000],trk,tp,&fSimPar);
	  FillTrackHistograms(fHist.fTrack[1004],trk,tp,&fSimPar, wt_b1);
	  FillTrackHistograms(fHist.fTrack[1005],trk,tp,&fSimPar, wt_b2);
	  if (trk->T0() > 700.  ) {
	    FillTrackHistograms(fHist.fTrack[1010],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[1014],trk,tp,&fSimPar, wt_b1);
	    FillTrackHistograms(fHist.fTrack[1015],trk,tp,&fSimPar, wt_b2);
	  }
	}
//-----------------------------------------------------------------------------
// passes candidate without veto requirement
//-----------------------------------------------------------------------------
	if (fEvtPar.fVetoedCandidate_BOX || fEvtPar.fCandidate_BOX) {
	  FillTrackHistograms(fHist.fTrack[1006],trk,tp,&fSimPar, wt_b1);
	  FillTrackHistograms(fHist.fTrack[1007],trk,tp,&fSimPar, wt_b2);
	  if (trk->T0() > 700.  ) {
	    FillTrackHistograms(fHist.fTrack[1016],trk,tp,&fSimPar, wt_b1);
	    FillTrackHistograms(fHist.fTrack[1017],trk,tp,&fSimPar, wt_b2);
	  }
	}
	if (fEvtPar.fNCrvClusters == 0) {
	  FillTrackHistograms(fHist.fTrack[1020],trk,tp,&fSimPar);
	  if (pid_ele) FillTrackHistograms(fHist.fTrack[1021],trk,tp,&fSimPar);
	  else         FillTrackHistograms(fHist.fTrack[1022],trk,tp,&fSimPar);
	  if (trk->T0() > 700.  ) { 
	    FillTrackHistograms(fHist.fTrack[1023],trk,tp,&fSimPar);
	    if (pid_ele) FillTrackHistograms(fHist.fTrack[1024],trk,tp,&fSimPar);
	    else         FillTrackHistograms(fHist.fTrack[1025],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// finally, checking the ~ signal band - mostly to see what the PID code is doing
// and limit the number of event candidates to look at
//-----------------------------------------------------------------------------
	    if (pid_ele and (tp->fP > 100) and (tp->fP < 110)) {
	      FillTrackHistograms(fHist.fTrack[1026],trk,tp,&fSimPar);
	    }
	  }
	}
      }
//-----------------------------------------------------------------------------
// positive tracks
//-----------------------------------------------------------------------------
      if (trk->Charge() > 0) {
	FillTrackHistograms(fHist.fTrack[120],trk,tp,&fSimPar);
	if (low_ep ) FillTrackHistograms(fHist.fTrack[ 121],trk,tp,&fSimPar);
	if (good_ep) FillTrackHistograms(fHist.fTrack[ 122],trk,tp,&fSimPar);
	if (high_ep) FillTrackHistograms(fHist.fTrack[ 123],trk,tp,&fSimPar);

	if (good_ep and good_dt) {
	  FillTrackHistograms(fHist.fTrack[124],trk,tp,&fSimPar);
	  if (fEvtPar.fNCrvClusters == 0) {
	    FillTrackHistograms(fHist.fTrack[125],trk,tp,&fSimPar);
	  }
	  if ((p >= 100) and (p < 110)) FillTrackHistograms(fHist.fTrack[ 126],trk,tp,&fSimPar);
	}

	if (fEvtPar.fCandidate_BOX) { 
	  FillTrackHistograms(fHist.fTrack[3000],trk,tp,&fSimPar);
	  FillTrackHistograms(fHist.fTrack[3004],trk,tp,&fSimPar, wt_b1);
	  FillTrackHistograms(fHist.fTrack[3005],trk,tp,&fSimPar, wt_b2);
	  if (trk->T0() > 700.  ) {
	    FillTrackHistograms(fHist.fTrack[3010],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[3014],trk,tp,&fSimPar, wt_b1);
	    FillTrackHistograms(fHist.fTrack[3015],trk,tp,&fSimPar, wt_b2);
	  }
	}
	// passes candidate without veto requirement
	if (fEvtPar.fVetoedCandidate_BOX || fEvtPar.fCandidate_BOX) {
	  FillTrackHistograms(fHist.fTrack[3006],trk,tp,&fSimPar, wt_b1);
	  FillTrackHistograms(fHist.fTrack[3007],trk,tp,&fSimPar, wt_b2);
	  if (trk->T0() > 700.  ) {
	    FillTrackHistograms(fHist.fTrack[3016],trk,tp,&fSimPar, wt_b1);
	    FillTrackHistograms(fHist.fTrack[3017],trk,tp,&fSimPar, wt_b2);
	  }
	}

	if (fEvtPar.fNCrvClusters == 0) {
	  FillTrackHistograms(fHist.fTrack[3020],trk,tp,&fSimPar);
	  if (pid_ele) FillTrackHistograms(fHist.fTrack[3021],trk,tp,&fSimPar);
	  else         FillTrackHistograms(fHist.fTrack[3022],trk,tp,&fSimPar);
	  if (trk->T0() > 700.  ) {
	    FillTrackHistograms(fHist.fTrack[3023],trk,tp,&fSimPar);
	    if (pid_ele) FillTrackHistograms(fHist.fTrack[3024],trk,tp,&fSimPar);
	    else         FillTrackHistograms(fHist.fTrack[3025],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// finally, checking the ~ signal band - mostly to see what the PID code is doing
// and limit the number of event candidates to look at
//-----------------------------------------------------------------------------
	    if (pid_ele and (tp->fP > 100) and (tp->fP < 110)) {
	      FillTrackHistograms(fHist.fTrack[3026],trk,tp,&fSimPar);
	    }
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// "good" tracks, MVA cuts
//-----------------------------------------------------------------------------
    if (tp->fIDWord[1] == 0) {

      FillTrackHistograms(fHist.fTrack[200],trk,tp,&fSimPar);

      if (fEvtPar.fCosmicVeto == 0) FillTrackHistograms(fHist.fTrack[201],trk,tp,&fSimPar);
      else                          FillTrackHistograms(fHist.fTrack[202],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// prototype ultimate cuts, both charges !
//-----------------------------------------------------------------------------
      if (fEvtPar.fCosmicVeto == 0) {
	if (pid_ele) {
	  FillTrackHistograms(fHist.fTrack[ 203],trk,tp,&fSimPar);
	  if ( fTCFinderBlockUe->NTimeClusters()>0){
	    FillTrackHistograms(fHist.fTrack[ 209],trk,tp,&fSimPar);
	    if ( (tp->fZCRV > 12000.) && (tp->fDtCRVCorrTof<-70)) {
	      //FillTrackHistograms(fHist.fTrack[ 209],trk,tp,&fSimPar);
	      if (GetDebugBit(37) == 1) {
		GetHeaderBlock()->Print(Form("bit_037: NTracksDe = %d NTracksUe = %d", fTCFinderBlockDe->NTimeClusters(), fTCFinderBlockUe->NTimeClusters()));
	      }
	    }
	    
	  }else {
	    if ( (tp->fZCRV > 12000.) && (tp->fDtCRVCorrTof<-70)) {
	      if (GetDebugBit(38) == 1) {
		GetHeaderBlock()->Print("bit_038: event with no TimeClusterUe and NTracksDe<2");
	      }
	    }
	  }	  
	}
	else         FillTrackHistograms(fHist.fTrack[ 204],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// TRK_205: events with E/P > 1.05
//-----------------------------------------------------------------------------
	if (high_ep) FillTrackHistograms(fHist.fTrack[ 205],trk,tp,&fSimPar);
      }

      int crv_veto     = (fEvtPar.fCosmicVeto &  Mu2eII::kCrvStubVetoBit);
      int non_crv_veto = (fEvtPar.fCosmicVeto & ~Mu2eII::kCrvStubVetoBit);

      if (pid_ele) {
	FillTrackHistograms(fHist.fTrack[ 206],trk,tp,&fSimPar);
	if (non_crv_veto == 0)                               FillTrackHistograms(fHist.fTrack[ 207],trk,tp,&fSimPar);
	if ((fEvtPar.fNCrvClusters > 0) and (crv_veto == 0)) FillTrackHistograms(fHist.fTrack[ 208],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// negative tracks
//-----------------------------------------------------------------------------
      if (trk->Charge() < 0) {
	FillTrackHistograms(fHist.fTrack[ 210],trk,tp,&fSimPar);

	if (pid_ele) FillTrackHistograms(fHist.fTrack[ 211],trk,tp,&fSimPar);
	else         FillTrackHistograms(fHist.fTrack[ 212],trk,tp,&fSimPar);
	if (high_ep) FillTrackHistograms(fHist.fTrack[ 213],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// poor man's PID - will go away ? - may be worth keeping around, as these cuts are easy to understand
//-----------------------------------------------------------------------------
	if (good_ep and good_dt) {
	  FillTrackHistograms(fHist.fTrack[214],trk,tp,&fSimPar);
	  
	  if (fEvtPar.fNCrvClusters == 0 ) FillTrackHistograms(fHist.fTrack[215],trk,tp,&fSimPar);

	  if ((p >= 100) and (p < 110)    ) FillTrackHistograms(fHist.fTrack[216],trk,tp,&fSimPar);
	  
	  if (fEvtPar.fNTimeClusters == 1) FillTrackHistograms(fHist.fTrack[217],trk,tp,&fSimPar);
	  else                             FillTrackHistograms(fHist.fTrack[218],trk,tp,&fSimPar);
	  
	  if      (trk->fPdgCode ==  11  ) FillTrackHistograms(fHist.fTrack[240],trk,tp,&fSimPar);
	  else if (trk->fPdgCode == -11  ) FillTrackHistograms(fHist.fTrack[241],trk,tp,&fSimPar);
	  else if (trk->fPdgCode ==  13  ) FillTrackHistograms(fHist.fTrack[242],trk,tp,&fSimPar);
	  else if (trk->fPdgCode == -13  ) FillTrackHistograms(fHist.fTrack[243],trk,tp,&fSimPar);
	}
//-----------------------------------------------------------------------------
// negative tracks, poor man's PID + 10 MeV momentum window
//-----------------------------------------------------------------------------
	if ( (good_ep and good_dt) and (trk->fP0 > 100.) and (trk->fP0 < 110.)) {
       
	  if (fEvtPar.fNEffTimeClusters == 1) FillTrackHistograms(fHist.fTrack[262],trk,tp,&fSimPar);
	  if (fEvtPar.fNEffTimeClusters  > 1) FillTrackHistograms(fHist.fTrack[263],trk,tp,&fSimPar);
	  
	  if ((tp->fDtCRVCorr > 30.) || (tp->fDtCRVCorr < -50.)) {
	    if ((fEvtPar.fNEffTimeClusters == 1) || (fEvtPar.fAbsTimeClusterDt > 200.)) {
	      FillTrackHistograms(fHist.fTrack[264],trk,tp,&fSimPar);
	      if (GetDebugBit(14) == 1) {
		if (fEvtPar.fCutCounter[3] == 1) {
		  GetHeaderBlock()->Print(Form("I am an outlier! ntc = %3i, nCRVStubs = %3i, t0(trk) = %6.3f, dT = %6.3f", 
					       fTimeClusterBlockDe->NTimeClusters(), fCrvClusterBlock->NClusters(), 
					       trk->T0(), tp->fDtCRVCorr));
		}
	      }
	    }
	  }
	}
//-----------------------------------------------------------------------------
// MVA-based track selection
//-----------------------------------------------------------------------------
	if (fEvtPar.fCandidate_MVA == 1) { 
	  FillTrackHistograms(fHist.fTrack[2000],trk,tp,&fSimPar);
	  FillTrackHistograms(fHist.fTrack[2004],trk,tp,&fSimPar, wt_b1);
	  FillTrackHistograms(fHist.fTrack[2005],trk,tp,&fSimPar, wt_b2);
	  if (trk->T0() > 700.  ) {
	    FillTrackHistograms(fHist.fTrack[2010],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[2014],trk,tp,&fSimPar, wt_b1);
	    FillTrackHistograms(fHist.fTrack[2015],trk,tp,&fSimPar, wt_b2);
	  }
	}
//-----------------------------------------------------------------------------
// passes candidate without veto requirement
//-----------------------------------------------------------------------------
	if (fEvtPar.fVetoedCandidate_MVA || fEvtPar.fCandidate_MVA) {
	  FillTrackHistograms(fHist.fTrack[2006],trk,tp,&fSimPar, wt_b1);
	  FillTrackHistograms(fHist.fTrack[2007],trk,tp,&fSimPar, wt_b2);
	  if (trk->T0() > 700.  ) {
	    FillTrackHistograms(fHist.fTrack[2016],trk,tp,&fSimPar, wt_b1);
	    FillTrackHistograms(fHist.fTrack[2017],trk,tp,&fSimPar, wt_b2);
	  }
	}
	
	if (fEvtPar.fNCrvClusters == 0) {
	  FillTrackHistograms(fHist.fTrack[2020],trk,tp,&fSimPar);
	  if (pid_ele) FillTrackHistograms(fHist.fTrack[2021],trk,tp,&fSimPar);
	  else         FillTrackHistograms(fHist.fTrack[2022],trk,tp,&fSimPar);
	  if (trk->T0() > 700.  ) { 
	    FillTrackHistograms(fHist.fTrack[2023],trk,tp,&fSimPar);
	    if (pid_ele) FillTrackHistograms(fHist.fTrack[2024],trk,tp,&fSimPar);
	    else         FillTrackHistograms(fHist.fTrack[2025],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// finally, checking the ~ signal band - mostly to see what the PID code is doing
// and limit the number of event candidates to look at
//-----------------------------------------------------------------------------
	    if (pid_ele and (tp->fP > 100) and (tp->fP < 110)) {
	      FillTrackHistograms(fHist.fTrack[2026],trk,tp,&fSimPar);
	      //	      GetHeaderBlock()->Print(Form("bit_2026: IDWord[1]=0:pid_ele:CosmicVeto=0 p = %10.3f",tp->fP*trk->Charge()));
	    }
	  }
	}
      }
//-----------------------------------------------------------------------------
// positive tracks
//-----------------------------------------------------------------------------
      if (trk->Charge() > 0) {
	FillTrackHistograms(fHist.fTrack[220],trk,tp,&fSimPar);
	if (low_ep ) FillTrackHistograms(fHist.fTrack[ 221],trk,tp,&fSimPar);
	if (good_ep) FillTrackHistograms(fHist.fTrack[ 222],trk,tp,&fSimPar);
	if (high_ep) FillTrackHistograms(fHist.fTrack[ 223],trk,tp,&fSimPar);

	if (good_ep and good_dt)   	          FillTrackHistograms(fHist.fTrack[224],trk,tp,&fSimPar);
	
	if (pid_ele) {
	  if (fEvtPar.fNCrvClusters == 0) FillTrackHistograms(fHist.fTrack[225],trk,tp,&fSimPar);
	  if ((p >= 100) and (p < 110)  ) FillTrackHistograms(fHist.fTrack[226],trk,tp,&fSimPar);
	}
	
	if (fEvtPar.fCandidate_MVA == 1) { 
	  FillTrackHistograms(fHist.fTrack[4000],trk,tp,&fSimPar);
	  FillTrackHistograms(fHist.fTrack[4004],trk,tp,&fSimPar, wt_b1);
	  FillTrackHistograms(fHist.fTrack[4005],trk,tp,&fSimPar, wt_b2);
	  if (trk->T0() > 700.  ) {
	    FillTrackHistograms(fHist.fTrack[4010],trk,tp,&fSimPar);
	    FillTrackHistograms(fHist.fTrack[4014],trk,tp,&fSimPar, wt_b1);
	    FillTrackHistograms(fHist.fTrack[4015],trk,tp,&fSimPar, wt_b2);
	  }
	}
//-----------------------------------------------------------------------------
// passes candidate without veto requirement
//-----------------------------------------------------------------------------
	if (fEvtPar.fVetoedCandidate_MVA || fEvtPar.fCandidate_MVA) {
	  FillTrackHistograms(fHist.fTrack[4006],trk,tp,&fSimPar, wt_b1);
	  FillTrackHistograms(fHist.fTrack[4007],trk,tp,&fSimPar, wt_b2);
	  if (trk->T0() > 700.  ) {
	    FillTrackHistograms(fHist.fTrack[4016],trk,tp,&fSimPar, wt_b1);
	    FillTrackHistograms(fHist.fTrack[4017],trk,tp,&fSimPar, wt_b2);
	  }
	}

	if (fEvtPar.fNCrvClusters == 0) {
	  FillTrackHistograms(fHist.fTrack[4020],trk,tp,&fSimPar);
	  if (pid_ele) FillTrackHistograms(fHist.fTrack[4021],trk,tp,&fSimPar);
	  else         FillTrackHistograms(fHist.fTrack[4022],trk,tp,&fSimPar);
	  if (trk->T0() > 700.  ) {
            FillTrackHistograms(fHist.fTrack[4023],trk,tp,&fSimPar);
	    if (pid_ele) FillTrackHistograms(fHist.fTrack[4024],trk,tp,&fSimPar);
	    else         FillTrackHistograms(fHist.fTrack[4025],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// finally, checking the ~ signal band - mostly to see what the PID code is doing
// and limit the number of event candidates to look at
//-----------------------------------------------------------------------------
	    if (pid_ele and (tp->fP > 100) and (tp->fP < 110)) {
	      FillTrackHistograms(fHist.fTrack[4026],trk,tp,&fSimPar);
	      //	      GetHeaderBlock()->Print(Form("bit_4026: IDWord[1]=0:pid_ele:CosmicVeto=0 p = %10.3f",tp->fP*trk->Charge()));
	    }
          }
	}
      }
//-----------------------------------------------------------------------------
// back to good MVA tracks, both charges
// TRK_265-286: by sector
//-----------------------------------------------------------------------------
      if ((fEvtPar.fCosmicVeto == 0) and pid_ele) {
	if (tp->fCRVSector ==  0) FillTrackHistograms(fHist.fTrack[265],trk,tp,&fSimPar);
	if (tp->fCRVSector ==  1) FillTrackHistograms(fHist.fTrack[266],trk,tp,&fSimPar);
	if (tp->fCRVSector ==  2) FillTrackHistograms(fHist.fTrack[267],trk,tp,&fSimPar);
	if (tp->fCRVSector ==  3) FillTrackHistograms(fHist.fTrack[268],trk,tp,&fSimPar);
	if (tp->fCRVSector ==  4) FillTrackHistograms(fHist.fTrack[269],trk,tp,&fSimPar);
	if (tp->fCRVSector ==  5) FillTrackHistograms(fHist.fTrack[270],trk,tp,&fSimPar);
	if (tp->fCRVSector ==  6) FillTrackHistograms(fHist.fTrack[271],trk,tp,&fSimPar);
	if (tp->fCRVSector ==  7) FillTrackHistograms(fHist.fTrack[272],trk,tp,&fSimPar);
	if (tp->fCRVSector ==  8) FillTrackHistograms(fHist.fTrack[273],trk,tp,&fSimPar);
	if (tp->fCRVSector ==  9) FillTrackHistograms(fHist.fTrack[274],trk,tp,&fSimPar);
	if (tp->fCRVSector == 10) FillTrackHistograms(fHist.fTrack[275],trk,tp,&fSimPar);
	if (tp->fCRVSector == 11) FillTrackHistograms(fHist.fTrack[276],trk,tp,&fSimPar);
	if (tp->fCRVSector == 12) FillTrackHistograms(fHist.fTrack[277],trk,tp,&fSimPar);
	if (tp->fCRVSector == 13) FillTrackHistograms(fHist.fTrack[278],trk,tp,&fSimPar);
	if (tp->fCRVSector == 14) FillTrackHistograms(fHist.fTrack[279],trk,tp,&fSimPar);
	if (tp->fCRVSector == 15) FillTrackHistograms(fHist.fTrack[280],trk,tp,&fSimPar);
	if (tp->fCRVSector == 16) FillTrackHistograms(fHist.fTrack[281],trk,tp,&fSimPar);
	if (tp->fCRVSector == 17) FillTrackHistograms(fHist.fTrack[282],trk,tp,&fSimPar);
	if (tp->fCRVSector == 18) FillTrackHistograms(fHist.fTrack[283],trk,tp,&fSimPar);
	if (tp->fCRVSector == 19) FillTrackHistograms(fHist.fTrack[284],trk,tp,&fSimPar);
	if (tp->fCRVSector == 20) FillTrackHistograms(fHist.fTrack[285],trk,tp,&fSimPar);
	if (tp->fCRVSector == 21) FillTrackHistograms(fHist.fTrack[286],trk,tp,&fSimPar);
	if (tp->fCRVStubSlope > 0.) FillTrackHistograms(fHist.fTrack[287],trk,tp,&fSimPar);
	if (tp->fCRVStubSlope < 0.) FillTrackHistograms(fHist.fTrack[288],trk,tp,&fSimPar);
	if ((tp->fCRVSector == 0) && (tp->fCRVStubSlope > 0.)) FillTrackHistograms(fHist.fTrack[289],trk,tp,&fSimPar);
	if ((tp->fCRVSector == 0) && (tp->fCRVStubSlope < 0.)) FillTrackHistograms(fHist.fTrack[290],trk,tp,&fSimPar);
	if ((tp->fCRVSector == 10) && (tp->fCRVStubSlope > 0.)) FillTrackHistograms(fHist.fTrack[291],trk,tp,&fSimPar);
	if ((tp->fCRVSector == 10) && (tp->fCRVStubSlope < 0.)) FillTrackHistograms(fHist.fTrack[292],trk,tp,&fSimPar);
	if (tp->fCRVStubSlopeMCProduct < 0.) FillTrackHistograms(fHist.fTrack[293],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// TRK_301, TRK_302: debug tp->fDtCRV splitting - likely electrons and muons, to be removed
//-----------------------------------------------------------------------------
      if ((fEvtPar.fNCrvClusters > 0) and (tp->fDtCRV < 50) and (tp->fZCRV > 12000)) {
	FillTrackHistograms(fHist.fTrack[301],trk,tp,&fSimPar);
      }

      if ((fEvtPar.fNCrvClusters > 0) and (tp->fDtCRV > 50) and (tp->fZCRV > 12000)) {
	FillTrackHistograms(fHist.fTrack[302],trk,tp,&fSimPar);
      }
    }
  }


//-----------------------------------------------------------------------------
// fill upstream-only track blocks
//-----------------------------------------------------------------------------

  for (int i=0; i<fTrackBlockUe->NTracks(); ++i ) {
    TStnTrack* trk         = fTrackBlockUe->Track(i);
    Mu2eII::TrackPar_t* tp = fTrackParUe+i;

    FillTrackHistograms(fHist.fTrack[500],trk,tp,&fSimPar);
    if (tp->fIDWord[1] == 0) FillTrackHistograms(fHist.fTrack[501],trk,tp,&fSimPar);
    else FillTrackHistograms(fHist.fTrack[502],trk,tp,&fSimPar);
  }


//-----------------------------------------------------------------------------
// track-timeCluster timing
//-----------------------------------------------------------------------------
  Mu2eII::TrackTcPar_t ttc;

  for (int it=0; it<fEvtPar.fNTracksDe; it++) {
    TStnTrack* trk         = fTrackBlockDe->Track(it);
    Mu2eII::TrackPar_t* tp = fTrackParDe+it;
    bool pid_ele = (tp->fPidMvaOut[0] >= 0.5); 
    if ((!pid_ele) || (tp->fIDWord[1] != 0)){
      continue;
    }
    if (fTCFinderBlockUe->fNTimeClusters>0){
      TStnTimeCluster* tc = fTCFinderBlockUe->TimeCluster(0);
      ttc.fDt = trk->T0()-tc->T0();  //d. track with 0th u. tc, bouncing identifier
      ttc.fClusterZ = tc->ClusterZ();
      FillTrackTcHistograms(fHist.fTrackTc[0],&ttc);
    }
    // for (int ic=0; ic<fEvtPar.fNTimeClusters; ic++) {
    //   TStnTimeCluster* tc = fTimeClusterBlockDe->TimeCluster(ic);
    //   ttc.fDt = trk->T0()-tc->T0();
    //   FillTrackTcHistograms(fHist.fTrackTc[0],&ttc);
    //   if (tp->fIDWord[1] == 0) FillTrackTcHistograms(fHist.fTrackTc[1],&ttc);
    // }
  }
//-----------------------------------------------------------------------------
// track-CRV stub pairs
//-----------------------------------------------------------------------------
  for (int it=0; it<fEvtPar.fNTracksDe; it++) {
    // TStnTrack* trk         = fTrackBlockDe->Track(it);
    Mu2eII::TrackPar_t* tp = fTrackParDe+it;
    if (tp->fIDWord[1] == 0) {
      bool        pid_ele  = (tp->fPidMvaOut[0] >= 0.5);            // includes the pre-PID cuts
      for (int ic=0; ic<fEvtPar.fNCrvClusters; ic++) {
	CrvStubPar_t* crv_sp = fCRVStubPar+ic;
	FillTrackCrvStHistograms(fHist.fTrackCrvSt[0],tp,crv_sp);
	if   (pid_ele) FillTrackCrvStHistograms(fHist.fTrackCrvSt[1],tp,crv_sp);
	else           FillTrackCrvStHistograms(fHist.fTrackCrvSt[2],tp,crv_sp);
      }
    }
  }
//-----------------------------------------------------------------------------
// 5. CRV cluster histograms
//-----------------------------------------------------------------------------
  int ncrvc = fCrvClusterBlock->NClusters();

  for (int i=0; i<ncrvc; i++) {
    TCrvCoincidenceCluster* ccc = fCrvClusterBlock->Cluster(i);
    FillCrvClusterHistograms(fHist.fCrvCluster[0],ccc,fCRVStubPar);
    if (fCRVStubPar[0].fSector == 10) FillCrvClusterHistograms(fHist.fCrvCluster[1],ccc,fCRVStubPar);
    if (fCRVStubPar[0].fStubQN == 2) FillCrvClusterHistograms(fHist.fCrvCluster[2],ccc,fCRVStubPar);
    if (fCRVStubPar[0].fStubQN == 4) FillCrvClusterHistograms(fHist.fCrvCluster[3],ccc,fCRVStubPar);
    if (fCRVStubPar[0].fStubSlopeMCProduct < 0.) FillCrvClusterHistograms(fHist.fCrvCluster[4],ccc,fCRVStubPar);

  }
//-----------------------------------------------------------------------------
// 6. CRV pulse histograms
//-----------------------------------------------------------------------------
  int ncrvp = fCrvPulseBlock->NPulses();

  for (int i=0; i<ncrvp; i++) {
    TCrvRecoPulse* p = fCrvPulseBlock->Pulse(i);
    FillCrvPulseHistograms(fHist.fCrvPulse[0],p);
  }
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TCosmicAnaModule::Event(int ientry) {

  fEvtPar.fCutCounter[0] = 0;
  fEvtPar.fCutCounter[1] = 0;
  fEvtPar.fCutCounter[2] = 0;
  fEvtPar.fCutCounter[3] = 0;

  TLorentzVector        mom;

  fHelixBlockDe->GetEntry(ientry);
  fHelixBlockUe->GetEntry(ientry);

  fTrackSeedBlockDe->GetEntry(ientry);
  fTrackSeedBlockUe->GetEntry(ientry);

  fTrackBlockDe->GetEntry(ientry);
  fTrackBlockUe->GetEntry(ientry);

  fTimeClusterBlockDe->GetEntry(ientry);
  fTimeClusterBlockUe->GetEntry(ientry);

  fTCFinderBlockDe->GetEntry(ientry);
  fTCFinderBlockUe->GetEntry(ientry);
  fCTPFinderBlock->GetEntry(ientry);

  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fSpmcBlockVDet->GetEntry(ientry);

  fCrvClusterBlock->GetEntry(ientry);
  fCrvPulseBlock->GetEntry(ientry);

  fEventWeight              = 1.;
  fEvtPar.fDioLOWt          = 1.;
  fEvtPar.fDioLLWt          = 1.;

  fEvtPar.fNCrvClusters     = fCrvClusterBlock->NClusters();
  fEvtPar.fNCrvPulses       = fCrvPulseBlock->NPulses();
  fEvtPar.fNCrvCoincidences = fCrvPulseBlock->NCoincidences();

  fEvtPar.fNHelicesDe       = fHelixBlockDe->NHelices();
  fEvtPar.fNHelicesUe       = fHelixBlockUe->NHelices();

  fEvtPar.fTCType           = -10;

  int ntc                   = fTCFinderBlockDe->NTimeClusters();
  fEvtPar.fNTimeClusters    = ntc;
//-----------------------------------------------------------------------------
// calculate effective number of time clusters: check for dt between clusters, nsh and nch, and matching cluster energies
//-----------------------------------------------------------------------------
  float tcdt = 0;
  int   effntc = 0;
  if (ntc > 1){
    std::vector<int> dups;
    for(int i=0; i < ntc; i++){
      TStnTimeCluster* tci = fTCFinderBlockDe->TimeCluster(i);
      for(int j=i+1; j < ntc; j++){
	TStnTimeCluster* tcj = fTCFinderBlockDe->TimeCluster(j);
	tcdt = abs(tci->T0() - tcj->T0());
	if(std::find(dups.begin(), dups.end(), j) != dups.end()){
	  continue;
	}
	if(tcdt < 40.){
	  dups.push_back(j);
	}
      }
    }
    effntc = ntc - int(dups.size());
    fEvtPar.fNEffTimeClusters = effntc;
    fEvtPar.fNTCIndex[0] = -9999;
    fEvtPar.fNTCIndex[1] = -9999;
    for(int i=0; i < ntc; i++){
      if(std::find(dups.begin(), dups.end(), i) != dups.end()){
	continue;
      }
      if(fEvtPar.fNTCIndex[0] < 0.){
	fEvtPar.fNTCIndex[0] = i;
	continue;
      }
      if(fEvtPar.fNTCIndex[1] < 0.){
	fEvtPar.fNTCIndex[1] = i;
      }
    }
  }
  else{
    fEvtPar.fNEffTimeClusters = ntc;
  }
  
  if (fTrackBlockDe->NTracks() > 0){
    if (fEvtPar.fNEffTimeClusters > 1){
      TStnTimeCluster* tc0 = fTCFinderBlockDe->TimeCluster(fEvtPar.fNTCIndex[0]);
      TStnTimeCluster* tc1 = fTCFinderBlockDe->TimeCluster(fEvtPar.fNTCIndex[1]);
      fEvtPar.fTimeClusterDt    = tc1->T0() - tc0->T0();
      fEvtPar.fAbsTimeClusterDt = abs(tc1->T0() - tc0->T0());
      if ((tc0->ClusterZ() > 1865.) && (tc0->ClusterZ() < 1866.)) { // first TC comes from first cal disk
	if ((tc1->ClusterZ() > 1865.) && (tc1->ClusterZ() < 1866.)) fEvtPar.fTCType = 0; // second TC comes from first cal disk
	if ((tc1->ClusterZ() > 2565.) && (tc1->ClusterZ() < 2566.)) fEvtPar.fTCType = 1;  // second TC comes from second cal disk
	if ((tc1->ClusterZ() > -0.1) && (tc1->ClusterZ() < 0.1)) fEvtPar.fTCType = 2; // second TC has no cal cluster
      }
      else if ((tc0->ClusterZ() > 2565.) && (tc0->ClusterZ() < 2566.)){ // first TC comes from second cal disk
	if ((tc1->ClusterZ() > 1865.) && (tc1->ClusterZ() < 1866.)) fEvtPar.fTCType = 3; // second TC comes from first cal disk
	if ((tc1->ClusterZ() > 2565.) && (tc1->ClusterZ() < 2566.)) fEvtPar.fTCType = 4; // second TC comes from second cal disk
	if ((tc1->ClusterZ() > -0.1) && (tc1->ClusterZ() < 0.1)) fEvtPar.fTCType = 5; // second TC has no cal cluster
      }
      else if ((tc0->ClusterZ() > -0.1) && (tc0->ClusterZ() < 0.1)){ // first TC has no cal cluster
	if ((tc1->ClusterZ() > 1865.) && (tc1->ClusterZ() < 1866.)) fEvtPar.fTCType = 6; // second TC comes from first cal disk
	if ((tc1->ClusterZ() > 2565.) && (tc1->ClusterZ() < 2566.)) fEvtPar.fTCType = 7; // second TC comes from second cal disk
	if ((tc1->ClusterZ() > -0.1) && (tc1->ClusterZ() < 0.1)) fEvtPar.fTCType = 8; // second TC has no cal cluster
      }
      if (fEvtPar.fNEffTimeClusters > 2) {
	fEvtPar.fTimeClusterDt     = -1.e6;
	fEvtPar.fAbsTimeClusterDt  = -1.e6;
      }
    }
    else {
      fEvtPar.fTimeClusterDt       = 1.e6;
      fEvtPar.fAbsTimeClusterDt    = 1.e6;
    }
  }
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
    if ((abs(pdg_code) == fPDGCode) and (generator_code == fMCProcessCode)) {
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
  else {
    GetHeaderBlock()->Print("ERROR: fSimPar.fParticle = NULL");
  }

  fEvtPar.fNTracksDe   = fTrackBlockDe->NTracks();
  fEvtPar.fNTracksUe   = fTrackBlockUe->NTracks();
  fEvtPar.fNClusters   = fClusterBlock->NClusters();

  fEvtPar.fNGoodTracks = 0;

  fEvtPar.fNGoodTracksTotal = 0;

  InitTrackPar(fTrackBlockDe,fClusterBlock,fTrackParDe,&fSimPar);
  InitTrackPar(fTrackBlockUe,fClusterBlock,fTrackParUe,&fSimPar);

  for (int i=0; i<fEvtPar.fNTracksDe; i++) {
    TrackPar_t* tp = fTrackParDe+i;
    if (tp->fIDWord[1] != 0) {
      continue;
    }
    fEvtPar.fNGoodTracksTotal++;
  }
  for (int i1=0; i1<fEvtPar.fNTracksUe; i1++) {
    TrackPar_t* tp = fTrackParUe+i1;
    if (tp->fIDWord[1] != 0) {
      continue;
    }
    fEvtPar.fNGoodTracksTotal++;
  }

  InitCrvStubPar(fCrvClusterBlock, fCRVStubPar, fSimPar.fParticle);
//-----------------------------------------------------------------------------
// now, for each track compare its T0 to the timing of the CRV stub candidate(s),
// if those are present
//-----------------------------------------------------------------------------
  fEvtPar.fCandidate_BOX       = 0;
  fEvtPar.fCandidate_MVA       = 0;
  fEvtPar.fVetoedCandidate_BOX = 0;
  fEvtPar.fVetoedCandidate_MVA = 0;

  fEvtPar.fNEleCandidates_BOX  = 0;
  fEvtPar.fNEleCandidates_MVA  = 0;

  fNGoodTracks_BOX       = 0;
  fNGoodTracks_MVA       = 0;

//--------------------------------------------------------------------------------
// compare De and Ue Tracks
//--------------------------------------------------------------------------------
  for (int i=0; i< fEvtPar.fNTracksDe; i++) {
    TStnTrack*  trk = fTrackBlockDe->Track(i);
    TrackPar_t* tp  = fTrackParDe+i;
    float       td  = trk->T0();
    float       tdchi2 = trk->Chi2Dof();
    float       tdp = trk->fP0;
    if ( (tp->fIDWord[1] != 0)) {
      continue;
    }
    for (int i1=0; i1<fEvtPar.fNTracksUe; i1++) {
      TStnTrack* trk_ue = fTrackBlockUe->Track(i1);
      // TrackPar_t* tp_ue = fTrackParUe+i1;
      // if (tp_ue->fIDWord[1] != 0) {
      // 	continue;
      // }
      float       tu    = trk_ue->T0();
      float       dt    = td - tu;
      float      tuchi2 = trk_ue->Chi2Dof();
      float       tup   = trk_ue->fP0;
      fEvtPar.fDchiUe = tdchi2 - tuchi2;
      fEvtPar.fDtUe   = dt;
      fEvtPar.fDpUe   = tdp - tup;
      tp->fUeT0       = tu;
      tp->fUeP0       = tup;
      tp->fDtUe       = dt;
      tp->fDpUe       = tdp - tup;
      if (trk->fClusterE > 0.) {
	// if (tdchi2 < tuchi2) { //look for 'good downstream' tracks, d. track has better chi2 and see if it has u. track ~100ns earlier 
	if ((dt > 50.) && (dt < 200.)){
	  tp->fBounceDt   = dt;
	  tp->fBounceDp   = tdp - tup;
	  tp->fBounceDchi = tdchi2 - tuchi2;
	}
	//}
      }
      
      if (dt > 50.) {
        fEvtPar.fUeGate = dt;// will give us a time gate to veto bouncing particles using tracks of D.S. track is reco'd as D.S and U.S. track is reco'd as U.S.
      }
      else if (dt < 50.) {
	fEvtPar.fSameLegDchi2 = tdchi2 - tuchi2;
	fEvtPar.fSameLegDp    = tdp - tup;
	// compare chi2 to find if track is D.S. or U.S. in the case of track pairs close in time
	// removes fake D.S. tracks in case of misreco'd U.S. as D.S. and we have an instance of both U.S. and D.S.
	if (tuchi2 < tdchi2) {
	  fEvtPar.fDtUe_goodUe = dt;
	  //  fEvtPar.fCosmicVeto |= 0x16;
	  //break;
	}
	if (tuchi2 > tdchi2) {
	  fEvtPar.fDtUe_goodDe = dt;
	}
      }
      else if (dt > 50.) {
	fEvtPar.fDiffLegDchi2 = tdchi2 - tuchi2;
	fEvtPar.fDiffLegDp    = tdp - tup;
      }
    }
  }
//-----------------------------------------------------------------------------
// evalaute CRV corrected time
//-----------------------------------------------------------------------------
  for (int i=0; i< fEvtPar.fNTracksDe; i++) {
    TStnTrack*  trk = fTrackBlockDe->Track(i);
    TrackPar_t* tp  = fTrackParDe+i;
//-----------------------------------------------------------------------------
// find the closest CRV stub candidate
//-----------------------------------------------------------------------------
    for (int icl=0; icl<fEvtPar.fNCrvClusters; icl++) {
      if (icl >= kMaxCrvStubs) {
	GetHeaderBlock()->Print(Form("*** ERROR TCosmicAnaModule::Event::001: N(CRV) stubs > %i\n",kMaxCrvStubs));
	break;
      }
      TCrvCoincidenceCluster* crvcc  = fCrvClusterBlock->Cluster(icl);
      Mu2eII::CrvStubPar_t*   crv_sp = fCRVStubPar+icl;

//-----------------------------------------------------------------------------
// introduce analysis-based corrections
//-----------------------------------------------------------------------------
      float dt1 = trk->T0()-crv_sp->fCorrTimeProp; // crvcc->StartTime();
      float dt2 = dt1;
      float z   = crvcc->Position()->Z();

      if (z < 2000) {
//------------------------------------------------------------------------------
// region 1: Z< 2000 : more or less unknown
//------------------------------------------------------------------------------
	dt2 = dt2-28; 
      }
      else if ((z >=2000) and (z < 5000)) {
//------------------------------------------------------------------------------
// region 2: 2000 <= Z < 6500 : correction is, effectively, the same
//------------------------------------------------------------------------------
	dt2 = dt2 - 30;
      }
      else if ((z >= 5000) and (dt1 >= -1025)) {
//------------------------------------------------------------------------------
// region 3: Z > 6500, dt1 > 25 : "upstream decays in flight" 
// region 4: downstream DIF
//------------------------------------------------------------------------------
	dt2 = dt1-4.49-5.64e-3*z ; // fit: 4.49 + 0.00564*z

	if (z < 11000) {
	  double ddt2 = dt2 - 5.972 + 0.007676*z; //"rotate" downstream DIF back
	  if (fabs(ddt2) < 20) {
	    dt2 = ddt2;
	  } 
	}
      }

      if ((z > 10000.) and (dt1 < 25) and (dt1 >= -75)) {
//------------------------------------------------------------------------------
// region 5: reconstruction bugs
// interactions in the calorimeter coming through because of the 
// upstream TC finder using ANN...
//-----------------------------------------------------------------------------
	dt2 = dt1+13 ;
      }

      if (fabs(dt2) < fabs(tp->fDtCRV2)) {
    	tp->fDtCRV         = dt1;
	tp->fDtCRV2        = dt2;
    	tp->fZCRV          = crvcc->Position()->Z();
    	tp->fDtCRVCorr     = trk->T0()-crv_sp->fCorrTime;
	tp->fDtCRVCorrProp = trk->T0()-crv_sp->fCorrTimeProp;
	tp->fDtCRVCorrTof  = dt1-crv_sp->fCorrTimeTof;
	tp->fXCorrected    = crv_sp->fXCorrected;
      	tp->fCRVSector     = crv_sp->fSector;
	tp->fTwoEndBars    = crv_sp->fTwoEndBars;
	tp->fUeDtCRV       = tp->fUeT0-crvcc->StartTime();
	tp->fCRVStubSlope  = crv_sp->fStubDYDZ;
	tp->fCRVStubSlopeMCProduct = crv_sp->fStubSlopeMCProduct;
      }

      if (GetDebugBit(16) == 1) {
	if ((fEvtPar.fNCrvClusters > 0) and (tp->fIDWord[1] == 0) and (tp->fDtCRVCorr < -80.)) {
	  GetHeaderBlock()->Print(Form(" bit_016: N(CRV): %4i T(track)-T(CRV) = %10.3f Tp->fZCRV = %10.3f Tp->fDtCRV = %10.3f",
				       fEvtPar.fNCrvClusters,tp->fDtCRV,tp->fZCRV,tp->fDtCRV));
	}
      }
    }

    if (tp->fIDWord[0] == 0) fNGoodTracks_BOX++;
    if (tp->fIDWord[1] == 0) fNGoodTracks_MVA++;

    if ((tp->fIDWord[0] == 0) and (tp->fPidMvaOut[0] > 0.5)) fEvtPar.fNEleCandidates_BOX += 1; 
    if ((tp->fIDWord[1] == 0) and (tp->fPidMvaOut[0] > 0.5)) fEvtPar.fNEleCandidates_MVA += 1; 
  }
//-----------------------------------------------------------------------------
// vetoing cosmics
// 1. use tracker+calorimeter information
//-----------------------------------------------------------------------------
  fEvtPar.fCosmicVeto  = 0;
  NonCrvCosmicVeto(&fCosmicVetoData,&fEvtPar);

//-----------------------------------------------------------------------------
// 2. handle CRV stubs .. P.M. set to 60 ns
//-----------------------------------------------------------------------------
  for (int i=0; i<fEvtPar.fNTracksDe; i++) {
    TrackPar_t* tp  = fTrackParDe+i;
//-----------------------------------------------------------------------------
// window is asymmetric, as is the nature of the distribution
//-----------------------------------------------------------------------------
    if ((tp->fDtCRV2 > -50.) and (tp->fDtCRV2 < 80)) {
      fEvtPar.fCosmicVeto |= Mu2eII::kCrvStubVetoBit;
    }
  }
  //  printf("fEvtPar.fCosmicVeto = %08x\n",fEvtPar.fCosmicVeto);

  if ((fEvtPar.fNEleCandidates_BOX == 1) and (fEvtPar.fCosmicVeto == 0)) fEvtPar.fCandidate_BOX       = 1;
  else                                                                   fEvtPar.fVetoedCandidate_BOX = 1;

  if ((fEvtPar.fNEleCandidates_MVA == 1) and (fEvtPar.fCosmicVeto == 0)) fEvtPar.fCandidate_MVA       = 1;
  else                                                                   fEvtPar.fVetoedCandidate_MVA = 1;

  if (fEvtPar.fNTracksDe > 0) fEvtPar.fCutCounter[0] += 1;

  bool  conditions[3] = {false};

  for (int i=0; i<fEvtPar.fNTracksDe; ++i ) {
    TStnTrack*          trk = fTrackBlockDe->Track(i);
    Mu2eII::TrackPar_t* tp  = fTrackParDe+i; 
    
    //HERE GOES BUNCH OF IF SELECTIONS...
    int good_ep = (tp->fEp >= 0.6) and (tp->fEp < 1.05);
    int good_dt = (tp->fDt >=-3.5) and (tp->fDt < 1.0 );

    if ((tp->fIDWord[1] == 0) and (trk->Charge() < 0) and (good_ep and good_dt) and (trk->fP0 > 100.) and (trk->fP0 < 110.)) {
      if (conditions[0] == false){
	fEvtPar.fCutCounter[1] += 1;
	conditions[0]           = true;
      }

      if ((tp->fDtCRVCorr > 30.) || (tp->fDtCRVCorr < -50.)){
	if(conditions[1] == false){
	  fEvtPar.fCutCounter[2] += 1;
	  conditions[1]           = true;
	}
	if ((fEvtPar.fNEffTimeClusters == 1) || (fEvtPar.fAbsTimeClusterDt > 200.)){
	  if (conditions[2] == false){
	    fEvtPar.fCutCounter[3] += 1;
	    conditions[2]           = true;
	  }
	}
      }
    }
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TCosmicAnaModule::Debug() {

  if (GetDebugBit(3) == 1) {
    if ((fNGoodTracks_BOX > 0) and (fEvtPar.fNCrvClusters == 0)) {
      GetHeaderBlock()->Print(Form("good track, no CRV coincidences"));
    }
  }

  if (GetDebugBit(11) == 1) {
    GetHeaderBlock()->Print(Form("NTracks from track block = %d",fTrackBlockDe->NTracks()));
  }

  if (GetDebugBit(12) == 1) {
    if (fEvtPar.fAbsTimeClusterDt == 0.) {
      GetHeaderBlock()->Print(Form("event with 2 time clusters, TCDT = 0; ntc = %3i, nCRVStubs = %3i, tcdt = %8f", 
				   fTimeClusterBlockDe->NTimeClusters(), fCrvClusterBlock->NClusters(), fEvtPar.fAbsTimeClusterDt));
    }
  }

  int ntc = fTimeClusterBlockDe->NTimeClusters();

  for (int i=0; i< fEvtPar.fNTracksDe; i++) {
    TStnTrack*  trk     = fTrackBlockDe->Track(i);
    TrackPar_t* tp      = fTrackParDe+i;
    bool        pid_ele = (tp->fPidMvaOut[0] >= 0.5);            // includes the pre-PID cuts

    if (GetDebugBit(4) == 1) {
      if ((fEvtPar.fNCrvClusters > 0) and (tp->fIDWord[1] == 0) and (tp->fDtCRV < -100)) {
	GetHeaderBlock()->Print(Form("bit_004: N(CRV): %4i T(track)-T(CRV) = %10.3f",
				     fEvtPar.fNCrvClusters,tp->fDtCRV));
      }
    }

    if (GetDebugBit(5) == 1) {
      if (tp->fDtCRV < -40.) {
	GetHeaderBlock()->Print(Form("event with high abs Dt, dT = %6.3f, NCRVStubs = %3i, t0(trk) = %6.3f, DPf = %6.3f",
				     tp->fDtCRV,fCrvClusterBlock->NClusters(),trk->T0(),tp->fDpF));
      }
    }

    if (GetDebugBit(6) == 1) {
      if (tp->fDpF < -2.) {
	GetHeaderBlock()->Print(Form("event on tail of DPf dist, DPf = %6.3f, t0 of track = %6.3f",tp->fDpF,trk->T0()));
      }
    }

    if (GetDebugBit(9) == 1) {
      if (tp->fDtCRV > 100.) {
	GetHeaderBlock()->Print(Form("event with high positive Dt, dT = %6.3f, NCRVStubs = %3i, t0(trk) = %6.3f, DPf = %6.3f",
				     tp->fDtCRV,fCrvClusterBlock->NClusters(),trk->T0(),tp->fDpF));
      }
    }

    if (GetDebugBit(10) == 1) {
      if (ntc > 2) {
	GetHeaderBlock()->Print(Form("event with > 2 time clusters, ntc = %3i, nCRVStubs = %3i, t0(trk) = %6.3f, dT = %6.3f, DPf = %6.3f", 
				     fTimeClusterBlockDe->NTimeClusters(), fCrvClusterBlock->NClusters(), trk->T0(), tp->fDtCRV, tp->fDpF));
      }
    }

    if (GetDebugBit(13) == 1) {
      if(fEvtPar.fNCrvClusters > 0) {
      GetHeaderBlock()->Print(Form("print for prop time correction, dT = %6.3f, NCRVStubs = %3i, t0(trk) = %6.3f, propCorrdT = %6.3f",
				   tp->fDtCRV,fCrvClusterBlock->NClusters(),trk->T0(),
				   fCRVStubPar[0].fCorrTimeProp - fCrvClusterBlock->Cluster(0)->StartTime()));
      }
    }

    if (GetDebugBit(15) == 1) {
      if ((fEvtPar.fNCrvClusters > 0) and (tp->fIDWord[1] == 0) and (tp->fDtCRV < 50) and (tp->fZCRV > 12000)) {
	GetHeaderBlock()->Print(Form(" bit_015: N(CRV): %4i T(track)-T(CRV) = %10.3f Tp->fZCRV = %10.3f",
				     fEvtPar.fNCrvClusters,tp->fDtCRV,tp->fZCRV));
      }
    }

    if (GetDebugBit(16) == 1) {
      if ((fEvtPar.fNCrvClusters > 0) and (tp->fIDWord[1] == 0) and (tp->fDtCRVCorr < -80.)) {
	GetHeaderBlock()->Print(Form(" bit_016: N(CRV): %4i T(track)-T(CRV) = %10.3f Tp->fZCRV = %10.3f",
				     fEvtPar.fNCrvClusters,tp->fDtCRV,tp->fZCRV));
      }
    }

    if (GetDebugBit(17) == 1) {
      if ((fEvtPar.fNTimeClusters > 1) and (tp->fIDWord[1] == 0)){
	GetHeaderBlock()->Print(Form("bit_017: good track, N(time clusters): %2i", fEvtPar.fNTimeClusters));
      }
    }

    if (GetDebugBit(18) == 1) {
//-----------------------------------------------------------------------------
// bit 18: potentially, events electrons/positrons and lost downstream leg
//-----------------------------------------------------------------------------
      bool pid_ele = (tp->fPidMvaOut[0] >= 0.5); // includes the pre-PID cuts

      if ((tp->fIDWord[1] == 0) and (fEvtPar.fCosmicVeto == 0) and pid_ele) {
	if ((tp->fDtCRVCorr < -50) and (tp->fZCRV > 10000)) { 
	  GetHeaderBlock()->Print(Form("bit_018: good track, N(time clusters): %2i tp->fDtCRVCorr = %8.3f tp->fZCRV: %10.3f", 
				       fEvtPar.fNTimeClusters,tp->fDtCRVCorr,tp->fZCRV));
	}
      }
    }

//-----------------------------------------------------------------------------
// bit 19: lower leg of Ralf's V3059
//-----------------------------------------------------------------------------
    if (GetDebugBit(19) == 1) {
      if ((tp->fIDWord[1] == 0) and (fEvtPar.fCosmicVeto == 0)) {
	if ((tp->fDtCRV < 0) and (tp->fZCRV > 7000) and (tp->fZCRV < 10000)) { 
	  GetHeaderBlock()->Print(Form("bit_019: good track, N(time clusters): %2i tp->fDtCRV = %8.3f tp->fZCRV: %10.3f", 
				       fEvtPar.fNTimeClusters,tp->fDtCRV,tp->fZCRV));
	}
      }
    }

    if (GetDebugBit(20) == 1) {
      if (fEvtPar.fNCrvClusters == 0) {
	if ((tp->fIDWord[1] == 0) and (pid_ele)) {
	  GetHeaderBlock()->Print(Form("bit_020: CANDIDATE EVENT : N(CRV stub candidates)=0 IDWord[1]=0 pid_ele=%10.3f", 
				       tp->fPidMvaOut[0]));
	}
      }
    }

    if (GetDebugBit(21) == 1) {
      if ((tp->fIDWord[1] == 0) and (tp->fDtCRV > 0) and (tp->fZCRV > 10000)) {
	GetHeaderBlock()->Print(Form("bit_021: IDWord[1]=0 tp->DtCrv = %10.3f ZCrv = %10.3f",
				     tp->fDtCRV,tp->fZCRV));
      }
    }

    if (GetDebugBit(22) == 1) {
      if ((tp->fIDWord[1] == 0) and (fEvtPar.fNCrvClusters == 0)) { 
	GetHeaderBlock()->Print(Form("bit_022: IDWord[1]=0 fEvtPar.fNCrvClusters=0"));
      }
    }
//-----------------------------------------------------------------------------
// bit_023: upper tail of the DT corrected distribution - track of any sign
//-----------------------------------------------------------------------------
    if (GetDebugBit(23) == 1) {
      if ((tp->fIDWord[1] == 0) and (tp->fDtCRVCorr > 80.)) {
	GetHeaderBlock()->Print(Form("bit_023: IDWord[1]=0 tp->fDtCRVCorr=%10.3f",tp->fDtCRVCorr));
      }
    }
//-----------------------------------------------------------------------------
// bit_024: CRV background region 1 
//-----------------------------------------------------------------------------
    if (GetDebugBit(24) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele and (fEvtPar.fNCrvClusters > 0) and (fEvtPar.fCosmicVeto == 0)) {
	if ((tp->fDtCRVCorrProp < 0) and (tp->fZCRV > 11000)) {
	  GetHeaderBlock()->Print(Form("bit_024: tid:pid:crv_veto: tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV: %10.3f %10.3f %10.3f",
				       tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit_025: CRV background region 2
//-----------------------------------------------------------------------------
    if (GetDebugBit(25) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele and (fEvtPar.fNCrvClusters > 0) and (fEvtPar.fCosmicVeto == 0)) {
	if ((tp->fDtCRVCorrProp > 50) and (tp->fZCRV > 7000)) {
	  GetHeaderBlock()->Print(Form("bit_025: tid:pid:crv_veto: tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV: %10.3f %10.3f %10.3f",
				       tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit_026: CRV background, region 3
//-----------------------------------------------------------------------------
    if (GetDebugBit(26) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele and (fEvtPar.fNCrvClusters > 0) and (fEvtPar.fCosmicVeto == 0)) {
	if ((tp->fDtCRVCorrProp < 0) and (tp->fZCRV > 7000) and (tp->fZCRV < 10000)) {
	  GetHeaderBlock()->Print(Form("bit_026: tid:pid:crv_veto: tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV: %10.3f %10.3f %10.3f",
				       tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit_027: CRV background, region 4
//-----------------------------------------------------------------------------
    if (GetDebugBit(27) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele and (fEvtPar.fNCrvClusters > 0) and (fEvtPar.fCosmicVeto == 0)) {
	if (tp->fZCRV < 2000) {
	  GetHeaderBlock()->Print(Form("bit_027: tid:pid:crv_veto: tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV: %10.3f %10.3f %10.3f",
				       tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit_028: CRV background, region 5
//-----------------------------------------------------------------------------
    if (GetDebugBit(28) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele and (fEvtPar.fNCrvClusters > 0) and (fEvtPar.fCosmicVeto == 0)) {
	if ((tp->fZCRV > 3000) and (tp->fZCRV < 4000)) {
	  GetHeaderBlock()->Print(Form("bit_028: tid:pid:crv_veto: tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV: %10.3f %10.3f %10.3f",
				       tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit_029: CRV background, just large delta T
//-----------------------------------------------------------------------------
    if (GetDebugBit(29) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele and (fEvtPar.fNCrvClusters > 0) and (fEvtPar.fCosmicVeto == 0)) {
	if ((tp->fDtCRVCorrProp > 100) or (tp->fDtCRVCorrProp < -40)) {
	  GetHeaderBlock()->Print(Form("bit_029: tid:pid:crv_veto: tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV: %10.3f %10.3f %10.3f",
				       tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit_030: CRV background, large track-CRV deltaT : |deltaT| > 200
//-----------------------------------------------------------------------------
    if (GetDebugBit(30) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele) {
	int ncrv = fEvtPar.fNCrvClusters;
	for (int ic=0; ic<ncrv; ic++) {
	  CrvStubPar_t* crv_sp = fCRVStubPar + ic;
	  float dt2 = trk->T0()-crv_sp->fCorrTimeProp;
	  if (fabs(dt2) > 200) {
	    GetHeaderBlock()->Print(Form("bit_030: tid:pid_ele: DtCRV2 = %10.3f",dt2));
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// bit_031: CRV background, just large dT2
//-----------------------------------------------------------------------------
    if (GetDebugBit(31) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele and (fEvtPar.fNCrvClusters > 0) and (fEvtPar.fCosmicVeto == 0)) {
	if ((tp->fDtCRV2 > 50) or (tp->fDtCRV2 < -50)) {
	  GetHeaderBlock()->Print(Form("bit_031: tid:pid:crv_veto: Dt,dt_corr_prop, Dt2,tp->fZCRV: %10.3f %10.3f %10.3f %10.3f",
				       tp->fDtCRV,tp->fDtCRVCorrProp,tp->fDtCRV2,tp->fZCRV));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit_32: trk_203:
//-----------------------------------------------------------------------------
    if (GetDebugBit(32) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele and (fEvtPar.fCosmicVeto == 0)) {
	GetHeaderBlock()->Print(Form("bit_032: tid:pid:crv_veto: Dt,dt_corr_prop, Dt2,tp->fZCRV: %10.3f %10.3f %10.3f %10.3f",
				       tp->fDtCRV,tp->fDtCRVCorrProp,tp->fDtCRV2,tp->fZCRV));
      }
    }
//-----------------------------------------------------------------------------
// bit_33: trk_203 + e- in 100-105 MeV/c
//-----------------------------------------------------------------------------
    if (GetDebugBit(33) == 1) {
      if ((tp->fIDWord[1] == 0) and pid_ele and (fEvtPar.fCosmicVeto == 0) and 
	  (trk->Charge() < 0) and (tp->fP > 100) and (tp->fP < 105)) {
	GetHeaderBlock()->Print(Form("bit_033: p:tid:pid:crv_veto: Dt,dt_corr_prop, Dt2,tp->fZCRV: %10.3f %10.3f %10.3f %10.3f %10.3f",
				     tp->fP*trk->fCharge,tp->fDtCRV,tp->fDtCRVCorrProp,tp->fDtCRV2,tp->fZCRV));
      }
    }
//-----------------------------------------------------------------------------
// bit_34: trk_206
//-----------------------------------------------------------------------------
    if (GetDebugBit(34) == 1) {
      int non_crv_veto = (fEvtPar.fCosmicVeto & ~Mu2eII::kCrvStubVetoBit);
      float dt         = tp->fDtCRVCorrProp; 
      if ((tp->fIDWord[1] == 0) and pid_ele and ((dt > 150) or (dt < - 30)) and (non_crv_veto == 0)) { 
	
	GetHeaderBlock()->Print(Form("bit_034: p:tid:pid_ele: dt,dt_corr_prop, zcrv: %10.3f %10.3f %10.3f",
				     tp->fDtCRV,tp->fDtCRVCorrProp,tp->fZCRV));
      }
    }
//-----------------------------------------------------------------------------
// bit_35: fEvtPar.fCandidate_MVA
//-----------------------------------------------------------------------------
    if (GetDebugBit(35) == 1) {
      if ((trk->Charge() < 0) and fEvtPar.fCandidate_MVA) {
	GetHeaderBlock()->Print(Form("bit_035: p:tid:pid_ele: dt,dt2,dt_corr_prop, zcrv: %8.1f %8.1f %8.1f %8.1f",
				     tp->fDtCRV,tp->fDtCRV2,tp->fDtCRVCorrProp,tp->fZCRV));
      }
    }
//--------------------------------------------------------------------------------
// bit_38: good De track with |dtCRV2| > 100
//--------------------------------------------------------------------------------
    if (GetDebugBit(38) == 1) {
       if ((tp->fIDWord[1] == 0) and(std::fabs(tp->fDtCRV2) > 100) ) {
	GetHeaderBlock()->Print(Form("bit_038: p:tid:pid:crv_veto: Dt,dt_corr_prop, Dt2,tp->fZCRV: %10.3f %10.3f %10.3f %10.3f %10.3f",
				     tp->fP*trk->fCharge,tp->fDtCRV,tp->fDtCRVCorrProp,tp->fDtCRV2,tp->fZCRV));
      }
    }
  }
}

//_____________________________________________________________________________
int TCosmicAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TCosmicAnaModule::Test001() {
}
}
