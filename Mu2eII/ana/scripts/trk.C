///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "Mu2eII/ana/scripts/modules.hh"

def_name Mu2eII_trk_0010("Mu2eII_track_ana");
def_name Mu2eII_trk_0011("Mu2eII_track_ana_dar");
def_name Mu2eII_trk_0012("Mu2eII_mva_test");
def_name Mu2eII_trk_0002("Mu2eII_cosmic_ana");
def_name Mu2eII_trk_0003("Mu2eII_pbar_ana");
def_name Mu2eII_trk_0004("Mu2eII_rmc_ana");
def_name Mu2eII_trk_0005("Mu2eII_rpc_ana");
def_name Mu2eII_trk_0006("Mu2eII_conv_ana");
//-----------------------------------------------------------------------------
// configure track analysis module
// 
// RunningMode = 1000*track_type+100*channel+10*UseTrqMVA + BatchMode ; 
//
// track_type = 0: PAR
//              1: DAR (default, interface to PAR not validated)
//
// channel    = 0: mu- --> e- (105 MeV/c)
//              1: mu- --> e+ (92 MeV/c)
//
// UseTrqMVA  = 0
//            = 1 use the ANN trained for DAR tracks
// 
// BatchMode  = 1 generated with 1 batch mode luminosity
//            = 2 generated with 2 batch mode luminosity
//            = 0 generated without pileup
//
// examples:  
// ---------
// RunninMode = 1010 : DAR tracks, 105 MeV e- training, no lumi reweighing
//              1111 : DAR tarcks,  92 MeV e+ training, lumi reweighing from 1-batch mode to 2-batch mode
//                     for 2-batch mode histogram sets
//                      (may break when running on single particle dataset with lumi=0)
//
// default track block: "TrackBlockDar"
// for cataloged datasets, the MCProcessCode and PDGCode are taken from the dataset catalogs
// for RunningMode < 10, have to use "TrackBlockPar"
// dioTail:7, conversions:43 
// for cosmic reco path, the track block name "TrackBlockDarDe" and "TrackBlockDarUe"
//-----------------------------------------------------------------------------
void Mu2eII_track_ana(const char* TrackBlockName = nullptr, int RunningMode = 1000, 
		      int PDGCode = 0, int MCProcessCode = -1, int DebugBit = -1) {

  Mu2eII::m_trk = (Mu2eII::TTrackAnaModule*) g.x->AddModule("Mu2eII::TTrackAnaModule",0);  

  int channel      = (RunningMode % 1000) / 100;            // 0:105 MeV, 1:92 MeV
  int use_trq_mva  = (RunningMode %  100) / 10;
  int batch_mode   =  RunningMode %   10;

  Mu2eII::m_trk->fBatchMode = batch_mode;

  if (use_trq_mva > 0) { 
    if      (channel == 0) Mu2eII::m_trk->SetTrqMVA("fele2s51b1",1070);  // signal: 105 MeV e-
    else if (channel == 1) Mu2eII::m_trk->SetTrqMVA("fpos2s51b1",1170);  // signal:  92 MeV e+
  }

  // PID MVA is always used

  if      (channel == 0) Mu2eII::m_trk->SetPidMVA("ele00s61b0",1000);
  else if (channel == 1) Mu2eII::m_trk->SetPidMVA("ele01s51b0",1100);

  Mu2eII::m_trk->SetPDGCode      (PDGCode);
  Mu2eII::m_trk->SetMCProcessCode(MCProcessCode);

  if (DebugBit >= 0) {
    Mu2eII::m_trk->SetDebugBit(DebugBit,1);
  }

  if (TrackBlockName != nullptr) {
    printf("set track block name:%s\n",TrackBlockName);
    Mu2eII::m_trk->SetTrackBlockNameDe(TrackBlockName);
    printf("--- done\n");
  }
}

//-----------------------------------------------------------------------------
// track block name fixed to TrackBlockDar
// if use_trq_mva is set to zero, use default calculation
//-----------------------------------------------------------------------------
void Mu2eII_track_ana_dar(int RunningMode = 1000, int PDGCode = 0, int MCProcessCode = -1) {

  Mu2eII::m_trk = (Mu2eII::TTrackAnaModule*) g.x->AddModule("Mu2eII::TTrackAnaModule",0);  

  int channel      = (RunningMode % 1000) / 100;            // 0:105 MeV, 1:92 MeV
  int use_trq_mva  = (RunningMode %  100) / 10;
  int batch_mode   =  RunningMode %   10;

  Mu2eII::m_trk->fBatchMode = batch_mode;

  if (use_trq_mva > 0) { 
    if      (channel == 0) Mu2eII::m_trk->SetTrqMVA("fele2s51b1",1070);  // signal: 105 MeV e-
    else if (channel == 1) Mu2eII::m_trk->SetTrqMVA("fpos2s51b1",1170);  // signal:  92 MeV e+
  }

  // PID MVA is always used

  if      (channel == 0) Mu2eII::m_trk->SetPidMVA("ele00s61b0",1000);
  else if (channel == 1) Mu2eII::m_trk->SetPidMVA("ele01s51b0",1100);

  Mu2eII::m_trk->SetPDGCode      (PDGCode);
  Mu2eII::m_trk->SetMCProcessCode(MCProcessCode);
}

//-----------------------------------------------------------------------------
// batch mode always 1
//-----------------------------------------------------------------------------
void Mu2eII_mva_test(const char* TrackBlockName = nullptr     , 
		     const char* MVATrainDS     = "fele2s51b1", 
		     int         MVATrainCode   = -1          , 
		     int         PDGCode        = 0           , 
		     int         MCProcessCode  = -1          ) {

  Mu2eII::m_trk = (Mu2eII::TTrackAnaModule*) g.x->AddModule("Mu2eII::TTrackAnaModule",0);  

  if (PDGCode       != 0) Mu2eII::m_trk->SetPDGCode      (PDGCode);
  if (MCProcessCode >= 0) Mu2eII::m_trk->SetMCProcessCode(MCProcessCode);

  if (TrackBlockName != nullptr) {
    Mu2eII::m_trk->SetTrackBlockNameDe(TrackBlockName);
  }

  if (MVATrainCode > 0) Mu2eII::m_trk->SetTrqMVA(MVATrainDS,MVATrainCode);

}

void Mu2eII_conv_ana(const char* TrackBlockName = nullptr, int RunningMode = 1000, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure CE analysis module
//-----------------------------------------------------------------------------
  Mu2eII::m_cnv = (Mu2eII::TConvAnaModule*) g.x->AddModule("Mu2eII::TConvAnaModule",0);  

  int channel      = (RunningMode % 1000) / 100;            // 0:105 MeV, 1:92 MeV
  int use_trq_mva  = (RunningMode %  100) / 10;
  int batch_mode   =  RunningMode %   10;

  Mu2eII::m_cnv->fBatchMode = batch_mode;

  if (use_trq_mva > 0) { 
    if      (channel == 0) Mu2eII::m_cnv->SetTrqMVA("fele2s51b1",1070);  // signal: 105 MeV e-
    else if (channel == 1) Mu2eII::m_cnv->SetTrqMVA("fpos2s51b1",1170);  // signal:  92 MeV e+
  }

  if (TrackBlockName != nullptr) Mu2eII::m_cnv->SetTrackBlockName(TrackBlockName);

  // PID MVA is always used

  if      (channel == 0) Mu2eII::m_cnv->SetPidMVA("ele00s61b0",1000);
  else if (channel == 1) Mu2eII::m_cnv->SetPidMVA("ele01s51b0",1100);

  if (DebugBit > 0) {
    Mu2eII::m_cnv->SetDebugBit(DebugBit,1);
  }
}

void Mu2eII_cosmic_ana(const char* TrackBlockDeName = nullptr, int RunningMode = 1000, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure cosmics module, no pileup
// cosmics module reads several track blocks - at the very least, De and Ue
//-----------------------------------------------------------------------------
  Mu2eII::m_cos = (Mu2eII::TCosmicAnaModule*) g.x->AddModule("Mu2eII::TCosmicAnaModule",0);  

  if (TrackBlockDeName != nullptr) Mu2eII::m_cos->SetTrackBlockName(0,TrackBlockDeName);

  int channel      = (RunningMode % 1000) / 100;            // 0:105 MeV, 1:92 MeV
  int use_trq_mva  = (RunningMode %  100) / 10;
  int batch_mode   =  RunningMode %   10;

  if (use_trq_mva > 0) { 
    if      (channel == 0) Mu2eII::m_cos->SetTrqMVA("fele2s51b1",1070);  // signal: 105 MeV e-
    else if (channel == 1) Mu2eII::m_cos->SetTrqMVA("fpos2s51b1",1170);  // signal:  92 MeV e+
  }

  // PID MVA is always used

  if      (channel == 0) Mu2eII::m_cos->SetPidMVA("ele00s61b0",1000);
  else if (channel == 1) Mu2eII::m_cos->SetPidMVA("ele01s51b0",1100);

  if (DebugBit > 0) {
    Mu2eII::m_cos->SetDebugBit(DebugBit,1);
  }
}

//-----------------------------------------------------------------------------
// configure pbar analysis module, no pileup
//-----------------------------------------------------------------------------
void Mu2eII_pbar_ana(const char* TrackBlockName = nullptr, int RunningMode = 1000, int DebugBit = -1) {

  Mu2eII::m_pbr = (Mu2eII::TPbarAnaModule*) g.x->AddModule("Mu2eII::TPbarAnaModule",0);  

  if (TrackBlockName != nullptr) Mu2eII::m_pbr->SetTrackBlockName(TrackBlockName);

  int channel      = (RunningMode % 1000) / 100;            // 0:105 MeV, 1:92 MeV
  int use_trq_mva  = (RunningMode %  100) / 10;
  int batch_mode   =  RunningMode %   10;

  if (use_trq_mva > 0) { 
    if      (channel == 0) Mu2eII::m_pbr->SetTrqMVA("fele2s51b1",1070);  // signal: 105 MeV e-
    else if (channel == 1) Mu2eII::m_pbr->SetTrqMVA("fpos2s51b1",1170);  // signal:  92 MeV e+
  }

  // PID MVA is always used

  if      (channel == 0) Mu2eII::m_pbr->SetPidMVA("ele00s61b0",1000);
  else if (channel == 1) Mu2eII::m_pbr->SetPidMVA("ele01s51b0",1100);

  if (DebugBit > 0) Mu2eII::m_pbr->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
// configure RMC analysis module
//
// RunningMode = 1000*flat input + 100*charge_mode + 10*UseTrqMVA + BatchMode, UseTrqMVA!=0 tells to use the ANN trained for DAR tracks
// examples:  
// ---------
// RunninMode = 0100: physics spectrum input, positrons, 105 MeV e- training, no lumi reweighing
//              1111: flat spectrum input, 92 MeV e+ training, lumi reweighing from 1-batch mode to 2-batch mode 
//                      (may break when running on single particle dataset with lumi=0)
// ---------
//
// Spectrum = 100*is internal + 10*internal version + external version
// examples:  
// ---------
// Spectrum = 0  : external RMC version 0 (Closure Approx.)
//            120: internal RMC version 2 (Plestid+Hill) external version 0 (Closure Approx.)
//
//-----------------------------------------------------------------------------
void Mu2eII_rmc_ana(const char* TrackBlockName = nullptr, double kmax = 90.1, int RunningMode = 1111,
		    int Spectrum = 000, int MCProcessCode = -1, int DebugBit = -1) {
  //Parse running mode options
  int FlatInput = RunningMode / 1000;
  RunningMode = RunningMode % 1000;
  int ChargeMode = RunningMode / 100;
  RunningMode = RunningMode % 100;
  int Training = RunningMode / 10;
  RunningMode = RunningMode % 10;
  int BatchMode = RunningMode;

  //Parse spectrum options
  bool negative = Spectrum < 0;
  Spectrum = abs(Spectrum);
  int Internal = Spectrum / 100;
  Spectrum = Spectrum % 100;
  int InternalVersion = Spectrum / 10;
  Spectrum = Spectrum % 10;
  int ExternalVersion = (1-2*negative)*Spectrum - negative; //if negative, make -(Spectrum) - 1

  Mu2eII::m_rmc = (Mu2eII::TRMCAnaModule*) g.x->AddModule("Mu2eII::TRMCAnaModule",0);  
  
  if (Training > 0) { 
    if   (ChargeMode == 0) Mu2eII::m_rmc->SetTrqMVA("fele2s51b1",1070);  // signal: 105 MeV e-
    else                   Mu2eII::m_rmc->SetTrqMVA("fpos2s51b1",1170);  // signal:  92 MeV e+
  }

  // PID MVA is always used

  if      (ChargeMode == 0) Mu2eII::m_rmc->SetPidMVA("ele00s61b0",1000);
  else if (ChargeMode == 1) Mu2eII::m_rmc->SetPidMVA("ele01s51b0",1100);

  if (TrackBlockName != nullptr) {
    Mu2eII::m_rmc->SetTrackBlockName(TrackBlockName);
  }

  Mu2eII::m_rmc->fBatchMode = BatchMode;

  if (DebugBit >= 0) {
    Mu2eII::m_rmc->SetDebugBit(DebugBit,1);
  }
  Mu2eII::m_rmc->fKMax = kmax;
  Mu2eII::m_rmc->fFlatInput = FlatInput;
  Mu2eII::m_rmc->fInternalRMC = Internal;
  Mu2eII::m_rmc->SetSpectrum(ExternalVersion);
  Mu2eII::m_rmc->fIntSpectrum = InternalVersion;
  Mu2eII::m_rmc->fVerbose = 1;
  if(MCProcessCode > 0) Mu2eII::m_rmc->fMCProcessCode = MCProcessCode;
}


//-----------------------------------------------------------------------------
// configure RPC analysis module
//-----------------------------------------------------------------------------
void Mu2eII_rpc_ana(const char* TrackBlockName = nullptr, int RunningMode = 1000, int DebugBit = -1) {

  Mu2eII::m_rpc = (Mu2eII::TRPCAnaModule*) g.x->AddModule("Mu2eII::TRPCAnaModule",0);  

  if (TrackBlockName != nullptr) Mu2eII::m_rpc->SetTrackBlockName(TrackBlockName);

  int channel      = (RunningMode % 1000) / 100;            // 0:105 MeV, 1:92 MeV
  int use_trq_mva  = (RunningMode %  100) / 10;
  int batch_mode   =  RunningMode %   10;

  Mu2eII::m_rpc->fBatchMode = batch_mode;

  if (use_trq_mva > 0) { 
    if      (channel == 0) Mu2eII::m_rpc->SetTrqMVA("fele2s51b1",1070);  // signal: 105 MeV e-
    else if (channel == 1) Mu2eII::m_rpc->SetTrqMVA("fpos2s51b1",1170);  // signal:  92 MeV e+
  }

  // PID MVA is always used

  if      (channel == 0) Mu2eII::m_rpc->SetPidMVA("ele00s61b0",1000);
  else if (channel == 1) Mu2eII::m_rpc->SetPidMVA("ele01s51b0",1100);

  if (DebugBit > 0) {
    Mu2eII::m_rpc->SetDebugBit(DebugBit,1);
  }
}
