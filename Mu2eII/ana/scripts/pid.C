///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "Mu2eII/ana/scripts/modules.hh"

def_name Mu2eII_pid_001("Mu2eII_pid_ana");
def_name Mu2eII_pid_002("Mu2eII_pid_ana_write_mva_tree");
//-----------------------------------------------------------------------------
// use MVA-based PID , assume DAR tracks
// - TrackBlockName and MVA training code need to be in sync - training depends on the track reco
// - if MinTrq < 0, the default stored in the training results is used
// defaults: TrackBlockName[0]: "TrackBlockDarDe"
//           TrackBlockName[1]: "TrackBlockDarDmu"
//-----------------------------------------------------------------------------
void  Mu2eII_pid_ana(const char* TrackBlockName=nullptr, int RunningMode = 1000, int PDGCode = 0, int MCProcessCode = -1, int DebugBit = -1) {
  Mu2eII::m_emu = (Mu2eII::TEmuAnaModule*) g.x->AddModule("Mu2eII::TEmuAnaModule",0);  

  if (TrackBlockName) Mu2eII::m_emu->SetTrackBlockName(0,TrackBlockName);

  int channel      = (RunningMode % 1000) / 100;            // 0:105 MeV, 1:92 MeV
  int use_trq_mva  = (RunningMode %  100) / 10;
  int batch_mode   =  RunningMode %   10;

  Mu2eII::m_emu->fBatchMode = batch_mode;

  
  if (use_trq_mva > 0) { 
    if      (channel == 0) Mu2eII::m_emu->SetTrqMVA("fele2s51b1",1070);  // signal: 105 MeV e-
    else if (channel == 1) Mu2eII::m_emu->SetTrqMVA("fpos2s51b1",1170);  // signal:  92 MeV e+
  }

  // PID MVA is always used

  if      (channel == 0) Mu2eII::m_emu->SetPidMVA("ele00s61b0",1000);
  else if (channel == 1) Mu2eII::m_emu->SetPidMVA("ele01s51b0",1100);


  if (DebugBit >= 0) Mu2eII::m_emu->SetDebugBit(DebugBit,1);

  Mu2eII::m_emu->SetPDGCode      (PDGCode);
  Mu2eII::m_emu->SetMCProcessCode(MCProcessCode);

  if (TrackBlockName != nullptr) {
    printf("set track block name:%s\n",TrackBlockName);
    Mu2eII::m_emu->SetTrackBlockName(0,TrackBlockName);
    printf("--- done\n");
  }
}

//-----------------------------------------------------------------------------
void  Mu2eII_pid_ana_write_mva_tree(const char* TrackBlockName=nullptr, int RunningMode = 1000, 
				    int PDGCode = 0, int MCProcessCode = -1, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  Mu2eII::m_emu = (Mu2eII::TEmuAnaModule*) g.x->AddModule("Mu2eII::TEmuAnaModule",0);  
  Mu2eII::m_emu->SetWriteMvaTree(1);

  int channel      = (RunningMode % 1000) / 100;            // 0:105 MeV, 1:92 MeV
  int use_trq_mva  = (RunningMode %  100) / 10;
  int batch_mode   =  RunningMode %   10;

  Mu2eII::m_emu->fBatchMode = batch_mode;

  if (use_trq_mva > 0) { 
    if      (channel == 0) Mu2eII::m_emu->SetTrqMVA("fele2s51b1",1070);  // signal: 105 MeV e-
    else if (channel == 1) Mu2eII::m_emu->SetTrqMVA("fpos2s51b1",1170);  // signal:  92 MeV e+
  }

  // PID MVA is always used

  if      (channel == 0) Mu2eII::m_emu->SetPidMVA("ele00s61b0",1000);
  else if (channel == 1) Mu2eII::m_emu->SetPidMVA("ele01s51b0",1100);


  if (DebugBit >= 0) Mu2eII::m_emu->SetDebugBit(DebugBit,1);

  Mu2eII::m_emu->SetPDGCode      (PDGCode);
  Mu2eII::m_emu->SetMCProcessCode(MCProcessCode);

  if (TrackBlockName != nullptr) {
    printf("set track block name:%s\n",TrackBlockName);
    Mu2eII::m_emu->SetTrackBlockName(0,TrackBlockName);
    printf("--- done\n");
  }

  if (DebugBit >= 0) Mu2eII::m_emu->SetDebugBit(DebugBit,1);
}


