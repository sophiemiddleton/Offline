//
// Empty digi maker

#include "CRVResponse/inc/MakeCrvRecoPulses.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CrvParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulseCollection.hh"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <TH1D.h>

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TMath.h>

namespace mu2e 
{
  class DigiMaker : public art::EDProducer 
  {

    public:
    explicit DigiMaker(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();
    private:
    std::string _crvDigiModuleLabel;
    
  };

  DigiMaker::DigiMaker(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _crvDigiModuleLabel(pset.get<std::string>("crvDigiModuleLabel"))
  {
    produces<CrvDigiCollection>();
  }

  void DigiMaker::beginJob()
  {
  }

  void DigiMaker::endJob()
  {
  }

  void DigiMaker::beginRun(art::Run &run)
  {
  }

  void DigiMaker::produce(art::Event& event) 
  {
    
    std::unique_ptr<CrvDigiCollection> crvDigiCollection2(new CrvDigiCollection);
    
    event.put(std::move(crvDigiCollection2));
    
  } // end produce

} // end namespace mu2e

using mu2e::DigiMaker;
DEFINE_ART_MODULE(DigiMaker)
