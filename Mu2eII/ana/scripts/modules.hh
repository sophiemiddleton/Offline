#ifndef Mu2eII_ana_scripts_modules_hh
#define Mu2eII_ana_scripts_modules_hh

#include "Mu2eII/ana/TTrackAnaModule.hh"
#include "Mu2eII/ana/TConvAnaModule.hh"
#include "Mu2eII/ana/TCosmicAnaModule.hh"
#include "Mu2eII/ana/TCrvAnaModule.hh"
#include "Mu2eII/ana/TEmuAnaModule.hh"
#include "Mu2eII/ana/TPbarAnaModule.hh"
#include "Mu2eII/ana/TRMCAnaModule.hh"
#include "Mu2eII/ana/TRPCAnaModule.hh"

namespace Mu2eII {
  Mu2eII::TConvAnaModule*          m_cnv   = nullptr;
  Mu2eII::TCosmicAnaModule*        m_cos   = nullptr;
  Mu2eII::TCrvAnaModule*           m_crv   = nullptr;
  Mu2eII::TEmuAnaModule*           m_emu   = nullptr;
  Mu2eII::TTrackAnaModule*         m_trk   = nullptr;
  Mu2eII::TPbarAnaModule*          m_pbr   = nullptr;
  Mu2eII::TRMCAnaModule*           m_rmc   = nullptr;
  Mu2eII::TRPCAnaModule*           m_rpc   = nullptr;
};

#endif
