/*//Author: S Middleton 
//Date: March 2020
//Purpose : To help analyze which features if the tracks belong to high momentum tracks

#include "DataProducts/inc/PDGCode.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TNtuple.h"
#include "Math/VectorUtil.h"
using namespace std;
using namespace mu2e;

namespace mu2e {

  class CosmicTrackSelector : public art::EDFilter {
  public:

    explicit CosmicTrackSelector(fhicl::ParameterSet const& pset);
    virtual bool filter  (art::Event& event) override;

  private:
    art::InputTag _strawHitsTag, _strawDigisTag, _strawDigiMCsTag;
    art::InputTag _panelHitsTag;
    int           _diagLevel;
    struct Cuts {
      double   pmin;
      double   pmax;
      unsigned minStrawDigis;
      unsigned minPlanes;
      int      minBackground;
      int      maxBackground;
      Cuts ( fhicl::ParameterSet const& pset):
        pmin(pset.get<double>("pmin")),
	pmax(pset.get<double>("pmax")),
        minStrawDigis(pset.get<unsigned>("minStrawDigis")),
        minPlanes(pset.get<unsigned>("minPlanes")),
        minBackground(pset.get<int>("minBackground")),
        maxBackground(pset.get<int>("maxBackground")){}
    };
    Cuts _cuts;
    TTree* _cosmic_track_selector;

   Int_t _evt;
   //CosmicTrack Features:
    Float_t _TrueA1;
    Float_t _TrueB1;
    Float_t _TrueA0;
    Float_t _TrueB0;
    Float_t _DOCA;
    Float_t _Chi2;
    

    Float_t _phiMC;
    Float_t _thetaMC;
    
    Float_t _NStrawHits;
    Float_t _NPanelHits;
  
    Float_t _NPlanes;
    Float_t _NPanels;
    Float_t _NStraws;

    Float_t _DigisPerMuon;
  
    Float_t _MomentumPre;
    Float_t _MomentumPost;

};
}
namespace {
  class DigisBySim {

  public:

    struct TrkMCInfo{

      std::vector<int> digi_indices;
      double p;

      TrkMCInfo():digi_indices(),p(0){}

      void addIndex( size_t i, double pStep ){
        p = std::max(p,pStep);
        digi_indices.push_back(i);
      }

      int at( size_t i ) const{
        return digi_indices.at(i);
      }

      size_t size() const{
        return digi_indices.size();
      }

    };
    typedef std::map<art::Ptr<mu2e::SimParticle>,TrkMCInfo > impl_type;

    DigisBySim( mu2e::StrawDigiMCCollection const& digimcs ):
      _bySim(){

      int n{-1};
      for ( auto const& digimc : digimcs ){
        auto const& sim_cal = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
        auto const& sim_hv  = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
        if ( sim_cal == sim_hv ){
          double p = sqrt(digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().mag2());
          _bySim[digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle()].addIndex(++n,p);
          ++_sum;
        }else{
          ++_nbad;
        }
      }
    }

    impl_type const& bySim() const { return _bySim; }
    int sum()  const { return _sum; }
    int nBad() const { return _nbad; }
    int nTracks() const  { return _bySim.size(); }

    int maxDigis() const {
      int maxD{0};
      for ( auto const& i : _bySim ){
        maxD = std::max( maxD, int(i.second.size()) );
      }
      return maxD;
    }

    int nGoodTracks ( int minDigis ){
      int ngood{0};
      for ( auto const& i : _bySim ){
        if ( int(i.second.digi_indices.size()) >= minDigis ){
          ++ngood;
        }
      }
      return ngood;
    }

    int nMuons () const{
      int nMu{0};
      for ( auto const& i : _bySim ){
        mu2e::SimParticle const& sim = *i.first;
        if ( std::abs(sim.pdgId()) == mu2e::PDGCode::mu_minus ){
          ++nMu;
        }
      }
      return nMu;
    }

  private:
    impl_type  _bySim;
    int        _sum  = 0; // Sum of all digis over all SimParticles
    int        _nbad = 0; // Number of digis rejected since opposite ends
    
  };

}

mu2e::CosmicTrackSelector::CosmicTrackSelector(fhicl::ParameterSet const& pset):
  art::EDFilter{pset},
  _strawHitsTag(pset.get<std::string>("strawHitsTag")),
  _strawDigisTag(pset.get<std::string>("strawDigisTag")),
  _strawDigiMCsTag(pset.get<std::string>("strawDigiMCsTag")),
  _panelHitsTag(pset.get<std::string>("panelHitsTag")),
  _diagLevel(pset.get<int>("diagLevel",0)),
  _cuts(pset.get<fhicl::ParameterSet>("filterCuts")){

  art::ServiceHandle<art::TFileService> tfs;
  _cosmic_track_selector->Branch("evt", &_evt,"evt/I");  
  _cosmic_track_selector->Branch("phiMC" ,&_phiMC,"phiMC/F"); 
  _cosmic_track_selector->Branch("thetaMC" ,&_thetaMC,"thetaMC/F"); 
  _cosmic_track_selector->Branch("NStrawHits" ,&_NStrawHits,"NStrawHits/F"); 
  _cosmic_track_selector->Branch("NPanelHits" ,&_NPanelHits,"NPanelHits/F"); 
  _cosmic_track_selector->Branch("NPlanes" ,&_NPlanes,"NPlanes/F"); 
  _cosmic_track_selector->Branch("NPanels" ,&_NPanels,"NPanels/F"); 
  _cosmic_track_selector->Branch("NStraws" ,&_NStraws,"NStraws/F"); 
  
  _cosmic_track_selector->Branch("NDigisPerMuon" ,&_DigisPerMuon,"NDigisPerMuon/F");
  _cosmic_track_selector->Branch("MomentumPre" ,&_MomentumPre,"MomentumPre/F");
  _cosmic_track_selector->Branch("MomentumPost" ,&_MomentumPost,"MomentumPost/F");



}

bool mu2e::CosmicTrackSelector::filter(art::Event& event) {

  bool retval{false};
  auto strawDigis   = event.getValidHandle<StrawDigiCollection>( _strawDigisTag );
  auto strawDigiMCs = event.getValidHandle<StrawDigiMCCollection>( _strawDigiMCsTag );
  auto strawHits    = event.getValidHandle<StrawHitCollection>( _strawHitsTag );
  auto comboHits    = event.getValidHandle<ComboHitCollection>( _strawHitsTag );
  auto panelHits    = event.getValidHandle<ComboHitCollection>( _panelHitsTag );
  
  _NStrawHits = (strawHits->size());
  _NPanelHits = (panelHits->size());
  
  DigisBySim digisBySim( *strawDigiMCs );
  
  for ( auto const& trkinfo : digisBySim.bySim() ){
    auto const& sim  = *trkinfo.first;
    const int nBackground = strawDigis->size() - trkinfo.second.digi_indices.size();
    const double pStart   = sim.startMomentum().vect().mag();
    const double p        = trkinfo.second.p;
    const double pDelta   = pStart-p;
    const bool isMuon = (std::abs(sim.pdgId()) == PDGCode::mu_minus);
    XYZVec momStart(sim.startMomentum().vect().x(),sim.startMomentum().vect().y(), sim.startMomentum().vect().z());
    
    double mag = sqrt((sim.startMomentum().vect().x()*sim.startMomentum().vect().x())+(sim.startMomentum().vect().x()*sim.startMomentum().vect().x())+(sim.startMomentum().vect().x()*sim.startMomentum().vect().x()));
    
    const double theta = acos(sim.startMomentum().vect().z()/mag);
    const double phi = atan(sim.startMomentum().vect().y()/sim.startMomentum().vect().x());
    
    if ( isMuon ) {
      _DigisPerMuon ( trkinfo.second.digi_indices.size() );
     
      _phiMC = phi;
      _thetaMC = theta_start;
    }

    set<int> planes;
    for ( int i : trkinfo.second.digi_indices ){
      auto const& digi = strawDigis->at(i);
      planes.insert( digi.strawId().getPlane() );
      
    }
    
    if ( isMuon && p  >= _cuts.pmin  &&  p  <= _cuts.pmax &&
         trkinfo.second.digi_indices.size() >= _cuts.minStrawDigis
         ){
      
      _MomentumPre = (p);//Pre Cuts

    for ( auto const& trkinfo : digisBySim.bySim() ){
    auto const& sim  = *trkinfo.first;
    const double p        = trkinfo.second.p;
    const bool isMuon = (std::abs(sim.pdgId()) == PDGCode::mu_minus);

    double p;
    int _sum, _nbad;
    for ( auto const& digimc : strawDigiMCs ){
        auto const& sim_cal = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
        auto const& sim_hv  = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
        if ( sim_cal == sim_hv ){
          double p = sqrt(digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().mag2());
          ++_sum;
        }else{
          ++_nbad;
        }

  if ( isMuon && p  >= _cuts.pmin  &&  p  <= _cuts.pmax &&
         trkinfo.second.digi_indices.size() >= _cuts.minStrawDigis){
         _MomentumPost = (p);

   }

}
}
}
}
}




using mu2e::CosmicTrackSelector;
DEFINE_ART_MODULE(CosmicTrackSelector);*/
