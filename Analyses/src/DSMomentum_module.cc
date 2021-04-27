
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "art/Framework/Principal/Provenance.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectory.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

namespace mu2e {

  class DSMomentum : public art::EDAnalyzer {
  public:

    typedef SimParticleCollection::key_type key_type;

    explicit DSMomentum(fhicl::ParameterSet const& pset) :
      art::EDAnalyzer(pset),
      _nAnalyzed(0),
      _stepModuleLabel(pset.get<std::string>("stepModuleLabel", "virtualdetector"))
    {
    }

    virtual ~DSMomentum() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);

    void analyze(const art::Event& e);

  private:
    int _nAnalyzed;

    TNtuple* _ntpstep;

    std::string _stepModuleLabel;
  };

  void DSMomentum::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
                           
  _ntpstep = tfs->make<TNtuple>( "ntpstep", "StepPoint ntuple",
                                  "EndMomx:EndMomy:EndMomz:EndPTot:PosX:PosY:PosZ:PDG:VolID"
                                  );
  }
  

  void DSMomentum::beginRun(art::Run const& run){

  }

  void DSMomentum::analyze(const art::Event& event) {
        //bool first = true;
        //bool inVolume = false;
        ++_nAnalyzed;

        // ntuple buffer.
        float nstep[_ntpstep->GetNvar()];

        art::Handle<StepPointMCCollection> steps;
        event.getByLabel(_stepModuleLabel, steps);
        
        
         const StepPointMCCollection& stepPC = *steps;
        
        for (const auto& step : stepPC) {
               
                nstep[0] = step.momentum().x();
                nstep[1] = step.momentum().y();
                nstep[2] = step.momentum().z();
                nstep[3] = sqrt(step.momentum().x()*step.momentum().x() + step.momentum().y()*step.momentum().y() +step.momentum().z()*step.momentum().z());
                nstep[4] = step.position().x();
                nstep[5] = step.position().y();
                nstep[6] = step.position().z();
                _ntpstep->Fill(nstep);
                 
                 art::Ptr<SimParticle> sim =  step.simParticle();
                 nstep[7] = sim->pdgId();
                 nstep[8] = step.volumeId();
             
                
        }
  }

}  // end namespace mu2e

using mu2e::DSMomentum;
DEFINE_ART_MODULE(DSMomentum);
