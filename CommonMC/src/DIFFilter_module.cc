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

  class DIFFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> SimToken{Name("SimParticleCollection"),Comment("")};
        fhicl::Atom<bool> makeplots{Name("makeplots"),false};
        fhicl::Atom<bool> isNull{Name("isNull"),true};
      };
      explicit DIFFilter(const art::EDFilter::Table<Config>& config);
      virtual bool filter(art::Event& event) override;

    private:
	  	art::InputTag _SimToken;
	  	const SimParticleCollection* _SimCol;
      bool makeplots_;
      bool isNull_;
      TTree* genTree;
      Float_t _maxr;
      Float_t _momT;
      Float_t _posT;
      Float_t _cosTheta;
      Float_t _time;
  };

  DIFFilter::DIFFilter(const art::EDFilter::Table<Config>& config) :
     EDFilter{config}
    , _SimToken(config().SimToken())
    , makeplots_{config().makeplots()}
    , isNull_{config().isNull()}
  {
    if(makeplots_){
      art::ServiceHandle<art::TFileService> tfs;
      genTree  = tfs->make<TTree>("GenAna", "GenAna");
      genTree->Branch("maxr", &_maxr, "maxr/F");   
      genTree->Branch("momT", &_momT, "momT/F"); 
      genTree->Branch("posT", &_posT, "posT/F");
      genTree->Branch("cosTheta", &_cosTheta, "cosTheta/F");
      genTree->Branch("time", &_time, "time/F");
    }
  }

  bool DIFFilter::filter(art::Event& event) {
    if(isNull_) return true;
    bool passed = false;
   
   /*std::cout<<"====================================="<<std::endl;

    // get all instances of products of type T
      std::vector<art::Handle<SimParticleCollection>> vah = event.getMany<SimParticleCollection>();
      
      // loop over the list of instances of products of this type
      for (auto const& ah : vah) {
          const art::Provenance* prov = ah.provenance();
          
          std::string fcn = prov->friendlyClassName();
          std::string modn = prov->moduleLabel();
          std::string instn = prov->processName();
          std::cout<<"extracting name =  "<<fcn<<" "<<modn<<" "<<instn<<std::endl; 
          std::cout<<"with type =  "<<typeid(prov).name()<<std::endl;
          auto _SimCol = ah.product();
          std::cout<<"sim size "<<(_SimCol->size()) <<std::endl;
          if(modn == "beamResampler"){
            if(_SimCol->size() !=0){
              for ( SimParticleCollection::const_iterator i=_SimCol->begin(); i!=_SimCol->end(); ++i ){
                SimParticle const& sim = i->second;
                std::cout<<" PDG "<<sim.pdgId()<<" creation "<< sim.creationCode()<<std::endl;
                if(makeplots_){ genTree->Fill();}
              }
            }
          }   
    }*/

    //------------SimParticles-------------//
    /*auto sH = event.getValidHandle<mu2e::SimParticleCollection>(_SimToken);
    _SimCol = sH.product();

    for ( SimParticleCollection::const_iterator i=_SimCol->begin(); i!=_SimCol->end(); ++i ){
      SimParticle const& sim = i->second;
      if(sim.pdgId() == 11 and sim.creationCode()== 14 and sim.parent()->pdgId() == 13){ passed = true;}
      std::cout<<" PDG "<<sim.pdgId()<<" creation "<< sim.creationCode()<<" parent "<< sim.parent()->pdgId()<<std::endl;
      if(makeplots_){ genTree->Fill();}
    }*/
    return passed;
  }
}

using mu2e::DIFFilter;
DEFINE_ART_MODULE(DIFFilter)
