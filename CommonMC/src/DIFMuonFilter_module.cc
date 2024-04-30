// Purpose: Event filter for DIF simulations
// author: S Middleton  2024
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
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

  class DIFMuonFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> SimToken{Name("SimParticleCollection"),Comment("tag for Sim collection")};
        fhicl::Atom<double> max_mom{Name("max_mom"),1e7};
        fhicl::Atom<double> min_mom{Name("min_mom"),0};
        fhicl::Atom<double> max_pos{Name("max_pos"),1e7};
        fhicl::Atom<double> min_pos{Name("min_pos"),0};
        fhicl::Atom<bool> makeplots{Name("makeplots"),false};
        fhicl::Atom<bool> isNull{Name("isNull"),true};
      };
      explicit DIFMuonFilter(const art::EDFilter::Table<Config>& config);
      virtual bool filter(art::Event& event) override;

    private:
      art::InputTag _SimToken;
	    const SimParticleCollection* _SimCol;
	  
      double max_mom_;
      double min_mom_;
      double max_pos_;
      double min_pos_;
      bool makeplots_;
      bool isNull_;
      TTree* genTree;
      Int_t _pdg;
      Int_t _stopCode;
      Float_t _momT;
      Float_t _posT;
      Float_t _cosTheta;
      Float_t _time;
  };

  DIFMuonFilter::DIFMuonFilter(const art::EDFilter::Table<Config>& config) :
     EDFilter{config}
    , _SimToken(config().SimToken())
    , max_mom_(config().max_mom())
    , min_mom_(config().min_mom())
    , max_pos_(config().max_pos())
    , min_pos_(config().min_pos())
    , makeplots_{config().makeplots()}
    , isNull_{config().isNull()}
  {
    if(makeplots_){
      art::ServiceHandle<art::TFileService> tfs;
      genTree  = tfs->make<TTree>("GenAna", "GenAna");
      genTree->Branch("pdg", &_pdg, "pdg/I");
      genTree->Branch("stopCode", &_stopCode, "stopCode/I");
      genTree->Branch("momT", &_momT, "momT/F"); 
      genTree->Branch("posT", &_posT, "posT/F");
      genTree->Branch("cosTheta", &_cosTheta, "cosTheta/F");
      genTree->Branch("time", &_time, "time/F");
    }
  }

  bool DIFMuonFilter::filter(art::Event& event) {
    if(isNull_) return true;
    bool passed = true;
    //std::cout<<"====================================="<<std::endl;

    // get all instances of products of type T
      std::vector<art::Handle<SimParticleCollection>> vah = event.getMany<SimParticleCollection>();
      
      /*std::vector<art::Handle<StepPointMCCollection>> vab = event.getMany<StepPointMCCollection>();
      for (auto const& ah : vab) {
          const art::Provenance* prov = ah.provenance();
          
          std::string fcn = prov->friendlyClassName();
          std::string modn = prov->moduleLabel();
          std::string instn = prov->processName();
          std::cout<<"extracting name =  "<<fcn<<" "<<modn<<" "<<instn<<std::endl; 
          std::cout<<"with type =  "<<typeid(prov).name()<<std::endl;
      }*/
      // loop over the list of instances of products of this type
      for (auto const& ah : vah) {
          const art::Provenance* prov = ah.provenance();
          
          std::string fcn = prov->friendlyClassName();
          std::string modn = prov->moduleLabel();
          std::string instn = prov->processName();
          //std::cout<<"extracting name =  "<<fcn<<" "<<modn<<" "<<instn<<std::endl; 
          //std::cout<<"with type =  "<<typeid(prov).name()<<std::endl;
          auto _SimCol = ah.product();
          //std::cout<<"sim size "<<(_SimCol->size()) <<std::endl;
          if(modn == "beamResampler"){
            if(_SimCol->size() !=0){
              GeomHandle<DetectorSystem> det;
              for ( SimParticleCollection::const_iterator i=_SimCol->begin(); i!=_SimCol->end(); ++i ){
                SimParticle const& sim = i->second;
                //std::cout<<" PDG "<<sim.pdgId()<<" stop "<< sim.stoppingCode()<<std::endl;
                //std::cout<<"end mom "<<sim.endMomentum().rho()<<" pos "<< det->toDetector(sim.endPosition()).rho() <<std::endl;
                if(sim.pdgId() == 13 and  sim.endMomentum().rho() < 25 and det->toDetector(sim.endPosition()).rho() < 125 and det->toDetector(sim.endPosition()).rho() > 75){// and sim.endMomentum().rho() < 17 and det->toDetector(sim.endPosition()).rho() < 125 and det->toDetector(sim.endPosition()).rho() > 75){ 
                passed = true; }
                if(makeplots_ ){ 
                  _pdg = sim.pdgId();
                  _stopCode = sim.stoppingCode();
                  _momT = sim.endMomentum().rho();
                  _posT = det->toDetector(sim.endPosition()).rho();
                  _cosTheta = cos(atan2(_momT,sim.startMomentum().z()));
                  _time = sim.startGlobalTime();
                  genTree->Fill();
                
                }
              }
            }
          }   
    }
    return passed;
  }
}

using mu2e::DIFMuonFilter;
DEFINE_ART_MODULE(DIFMuonFilter)
