// Purpose: Event filter for DIO simulations
// author: S Middleton  2024
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include <iostream>
#include <string>

#include "TTree.h"
using namespace std;
namespace mu2e {

  class PionFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> SimToken{Name("SimParticleCollection"),Comment("")};
        fhicl::Atom<art::InputTag>SimTag{Name("SimParticleCollection"),Comment("SimTag")};
        fhicl::Atom<double> tmin{Name("tmin"),0};
        fhicl::Atom<bool> isNull{Name("isNull"),true};
      };
      explicit PionFilter(const art::EDFilter::Table<Config>& config);
      virtual bool filter(art::Event& event) override;

    private:
      art::InputTag _SimToken;
      const SimParticleCollection* _SimCol;
      double tmin_;
      bool isNull_;
      TTree* genTree;
      Float_t _endtime;
      Float_t _starttime;
  };

  PionFilter::PionFilter(const art::EDFilter::Table<Config>& config) :
     EDFilter{config}
    , _SimToken(config().SimToken())
    , tmin_{config().tmin()}
    , isNull_{config().isNull()}
  {
      art::ServiceHandle<art::TFileService> tfs;
      genTree  = tfs->make<TTree>("GenAna", "GenAna");
      genTree->Branch("endtime", &_endtime, "endtime/F");
      genTree->Branch("starttime", &_starttime, "starttime/F");
  }

  bool PionFilter::filter(art::Event& evt) {
    if(isNull_) return true;
    bool passed = false;

     std::vector<art::Handle<SimParticleCollection>> vah = evt.getMany<SimParticleCollection>();
      int n = 0;
      for (auto const& ah : vah) { //always one collection
        for(const auto& aParticle : *ah){
          
          art::Ptr<SimParticle> pp(ah, aParticle.first.asUint());
          _endtime = pp->endGlobalTime() ;
          _starttime = pp->startGlobalTime() ;
          if( (abs(pp->pdgId())  == 211 and _endtime > tmin_ )){
            passed = true; 
          }
          if(n==2 and (abs(pp->pdgId())  == 211)){ //geneolgy means we want the last pion (the stop)
            genTree->Fill();
          }
          n++;
        }
      }
    return passed;
  }
}

using mu2e::PionFilter;
DEFINE_ART_MODULE(PionFilter)
