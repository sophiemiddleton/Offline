// Module to plot generator output for any generator
// S. Middleton 2021

#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>
#include <vector>

#include "TTree.h"
using namespace std;
namespace mu2e {

  class RPCAna : public art::EDAnalyzer {

     public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> SimToken{Name("SimParticleCollection"),Comment("")};
        fhicl::Atom<int> pdg{Name("pdg"),Comment("")};
      };
      typedef art::EDAnalyzer::Table<Config> Parameters;

      explicit RPCAna(const Parameters& conf);
      virtual ~RPCAna() {}
      virtual void beginJob();
      virtual void endJob();
      virtual void analyze(const art::Event& e) override;

     private:

      art::InputTag _SimToken;
      int _pdg;
      const SimParticleCollection* _SimCol;

      TTree* genTree;
      Float_t _startmom;
      Float_t _startCode;
  };

  RPCAna::RPCAna(const Parameters& conf):
  art::EDAnalyzer(conf)
    , _SimToken(conf().SimToken())
    , _pdg(conf().pdg())
  {}

  void RPCAna::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
      genTree  = tfs->make<TTree>("GenAna", "GenAna");
      genTree->Branch("startmom", &_startmom, "startmom/F");
      genTree->Branch("startCode", &_startCode, "startCode/F");
  }

  void RPCAna::analyze(const art::Event& evt) {
     std::vector<art::Handle<SimParticleCollection>> vah = evt.getMany<SimParticleCollection>();
      for (auto const& ah : vah) { //always one collection
        for(const auto& aParticle : *ah){
          
          art::Ptr<SimParticle> pp(ah, aParticle.first.asUint());
          _startCode = pp->creationCode();
          _startmom = sqrt(pp->startMomXYZT().x()*pp->startMomXYZT().x() + pp->startMomXYZT().y()*pp->startMomXYZT().y() + pp->startMomXYZT().z()*pp->startMomXYZT().z());

          if( ((pp->pdgId())  == _pdg and _startCode == 179 )){
            genTree->Fill();
          }
        }
      }
    //return passed;
  }
  
 void RPCAna::endJob(){}
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::RPCAna)

