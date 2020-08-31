/*
// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

// cetlib includes
#include "cetlib_except/exception.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/PionCaptureSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "MCDataProducts/inc/FixedTimeMap.hh"
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "DataProducts/inc/EventWindowMarker.hh"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"

namespace mu2e {

  //================================================================
  class MotherMomAna : public art::EDAnalyzer {
    fhicl::ParameterSet psphys_;
    RootTreeSampler<IO::StoppedParticleTauNormF> stops_;

    TTree*  _Ntup;
    Float_t _PionMom;

  public:
    explicit MotherMomAna(const fhicl::ParameterSet& pset);
    virtual void beginRun(art::Run&   r) override;
    virtual void analyze(art::Event& event);
  };

  //================================================================
  MotherMomAna::MotherMomAna(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer{pset}
    , stops_(eng_, pset.get<fhicl::ParameterSet>("pionStops"))
    {
      
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "MotherMomAna" );
      _Ntup  = tfs->make<TTree>("GenTree", "GenTree");  
      _Ntup->Branch("PionMom", &_PionMom , "PionMom/F");
        
    }

  //================================================================

  void MotherMomAna::beginRun(art::Run& run) {
   
  }

  void MotherMomAna::analyze(art::Event& event) {

      IO::StoppedParticleTauNormF stop;
      double PionMOm = stop.pt;
      _Ntup->Fill();
    }


  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MotherMomAna);*/
