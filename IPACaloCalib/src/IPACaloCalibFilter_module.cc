
// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrackClusterMatch.hh"
// Mu2e includes.
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"

// BaBar Kalman filter includes
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/BaBar/BaBar.hh"

// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// Root includes
//#include "TNtuple.h"
#include "TTree.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <vector>

using namespace std;

namespace mu2e {

  class IPACaloCalibFilter : public art::EDFilter {

  public:
	struct Config {
		     using Name=fhicl::Name;
		     using Comment=fhicl::Comment;
	 	     fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
		     fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
		     fhicl::Atom<art::InputTag> kalrepTag{Name("KalRepPtrCollection"),Comment("outcome of Kalman filter (for tracker momentum info)")};
		     fhicl::Atom<art::InputTag> tcmatchTag{Name("TrackClusterMatchCollection"), Comment("track calo match"), "TrackCaloMatching"};
          fhicl::Atom<int> minTrackerHits {Name("minTrackerHits"), Comment("minimum number of straw hits in tracker "),25};//From Pasha's Study
          fhicl::Atom<int> minCrystalHits{Name("minCrystalHits"), Comment("minimum number of crystal hits "), 4};
		      fhicl::Atom<float> minClusterEDep {Name("minClusterEDep"), Comment("minimum amount of energy deposited "),10};
          fhicl::Atom<float>  maxTrackChi2 {Name("maxTrackChi2"), Comment("minimum allowed chi2 "),3};
	    	};
  
    typedef art::EDFilter::Table<Config> Parameters;

		explicit IPACaloCalibFilter(const Parameters& conf);
    virtual bool filter  (art::Event& event) override;
		virtual ~IPACaloCalibFilter() {}
   
  private:
    Config _conf;
		int _diagLevel;
		int _mcdiag;

		art::InputTag _kalrepTag;
		art::InputTag _tcmatchTag;
   
		const KalRepPtrCollection* _kalrepcol;		
		const TrackClusterMatchCollection* _tcmatchcol;

    int _minTrackerHits, _minCrystalHits;
    float _minClusterEDep, _maxTrackChi2, _momCut;
    
    bool passCH = false;
    bool passEdep = false;
    bool passChi2ndf = false;
    bool passMomCut  = false;
    bool passtrackerHit = false;
    bool findData(const art::Event& evt);
  };

  IPACaloCalibFilter::IPACaloCalibFilter(const Parameters& conf):
    art::EDFilter{conf},
    _diagLevel(conf().diagLevel()),
		_mcdiag(conf().mcdiag()),
		_kalrepTag(conf().kalrepTag()),		
		_tcmatchTag(conf().tcmatchTag()),
    _minTrackerHits(conf().minTrackerHits()),
    _minCrystalHits(conf().minCrystalHits()),
    _minClusterEDep(conf().minClusterEDep()),
    _maxTrackChi2(conf().maxTrackChi2())
{}

  bool IPACaloCalibFilter::filter(art::Event& event) {
    std::cout<<"[In IPACaloCalib::filter() ] ... "<<std::endl;
    if(!findData(event)) 
		  throw cet::exception("RECO")<<"No data in  event"<< endl; 

	  art::ServiceHandle<mu2e::GeometryService>   geom;
	
    mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
    mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;
   
    Hep3Vector vd_tt_back = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_Back));
    double     Z      = vd_tt_back.z();

    for(unsigned int i=0;i<_kalrepcol->size();i++){
      art::Ptr<KalRep> const& ptr = _kalrepcol->at(i);
      const KalRep* TrackKrep = ptr.get();
     
      double  ds(10.), s0, s1, s2, z0, z1, z2, dzds, sz, sz1, z01;
      const TrkHitVector* hots = &TrackKrep->hitVector();
      int nh = hots->size();
      if(nh > _minTrackerHits){
        passtrackerHit = true;
      }
      const TrkHit *first(nullptr), *last(nullptr);

      for (int ih=0; ih<nh; ++ih) {
        const TrkHit* hit = hots->at(ih);
        if (hit  != nullptr) {
          if (first == nullptr) first = hit;
          last = hit;
          }
    }

    s1 = first->fltLen();
    s2 = last ->fltLen();

    z1     = TrackKrep->position(s1).z();
    z2     = TrackKrep->position(s2).z();

    dzds   = (z2-z1)/(s2-s1);

    if (fabs(Z-z1) > fabs(Z-z2)) {
      z0 = z2;
      s0 = s2;
    }
    else {
      z0 = z1;
      s0 = s1;
    }

    sz    = s0+(Z-z0)/dzds;

    z0     = TrackKrep->position(sz).z();     // z0 has to be close to Z(TT_FrontPA)
    z01    = TrackKrep->position(sz+ds).z();

    dzds   = (z01-z0)/ds;
    sz1    = sz+(Z-z0)/dzds;	          // should be good enough

    double EndMom= TrackKrep->momentum(sz1).mag();
    if(EndMom > _momCut) passMomCut = true;
    //enact tracker cuts:
    if(TrackKrep->chisq()/TrackKrep->nActive() < _maxTrackChi2) {
      std::cout<<"passes chi2 "<<std::endl;
      passChi2ndf = true;
    }
    if(TrackKrep->chisq()/TrackKrep->nActive() > _maxTrackChi2) return false;
    //============================= Add Matches ===============================//
    if(_tcmatchcol->size() ==0) continue;
    for(unsigned int c=0;c<_tcmatchcol->size();c++){
        TrackClusterMatch const& tcm = (*_tcmatchcol)[c];
        const TrkCaloIntersect* extrk = tcm.textrapol();
        const KalRep* Krep  = extrk->trk().get();
        if (Krep == TrackKrep) {
          const mu2e::CaloCluster* cl = tcm.caloCluster();
            if(cl->size() > _minCrystalHits) { 
              std::cout<<"passes cl size "<<cl->size()<<std::endl;
              passCH = true;
          }
          if(cl->size() < _minCrystalHits) return false;
          if(cl->energyDep() > _minClusterEDep) {
            passEdep = true;
            std::cout<<"passes cl edep "<<cl->energyDep()<<std::endl;
          }
          if(cl->energyDep() < _minClusterEDep) return false;
      }
    }
  }
  if(passCH==true && passtrackerHit == true && passEdep==true &&  passChi2ndf==true ){ 
    std::cout<<"passes ALL"<<std::endl;
    return true;
  }
  else{ return false; }
 }

  

  bool IPACaloCalibFilter::findData(const art::Event& evt){

	
	  _tcmatchcol=0;
	  _kalrepcol = 0;

	  auto kalrep = evt.getValidHandle<KalRepPtrCollection>(_kalrepTag);
	  _kalrepcol =kalrep.product();
	  auto tcmatch = evt.getValidHandle<TrackClusterMatchCollection>(_tcmatchTag);
	  _tcmatchcol = tcmatch.product();
    

	  return _kalrepcol!=0 and _tcmatchcol!=0;
 }


}

using mu2e::IPACaloCalibFilter;
DEFINE_ART_MODULE(IPACaloCalibFilter);
