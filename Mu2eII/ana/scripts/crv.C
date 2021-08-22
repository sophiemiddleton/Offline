///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "Mu2eII/ana/scripts/modules.hh"

def_name Mu2eII_crv_001("Mu2eII_crv_ana");
//-----------------------------------------------------------------------------
// beam-induced CRV noise
//-----------------------------------------------------------------------------
void  Mu2eII_crv_ana(int DebugBit = -1) {
  Mu2eII::m_crv = (Mu2eII::TCrvAnaModule*) g.x->AddModule("Mu2eII::TCrvAnaModule",0);  

  if (DebugBit >= 0) Mu2eII::m_crv->SetDebugBit(DebugBit,1);
}
