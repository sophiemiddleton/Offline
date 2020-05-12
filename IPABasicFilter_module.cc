
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

  class IPABasicFilter : public art::EDFilter {

  public:
	struct Config {
		     using Name=fhicl::Name;
		     using Comment=fhicl::Comment;
	 	     fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
		     fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
		     fhicl::Atom<art::InputTag> tcmatchTag{Name("TrackClusterMatchCollection"), Comment("track calo match"), "TrackCaloMatching"};
          fhicl::Atom<int> minCrystalHits{Name("minCrystalHits"), Comment("minimum number of crystal hits "), 2};
         
	    	};
  
    typedef art::EDFilter::Table<Config> Parameters;

		explicit IPABasicFilter(const Parameters& conf);
    virtual bool filter  (art::Event& event) override;
		virtual ~IPABasicFilter() {}
   
  private:
    Config _conf;
		int _diagLevel;
		int _mcdiag;

		art::InputTag _tcmatchTag;
   		
		const TrackClusterMatchCollection* _tcmatchcol;

    int _minTrackerHits, _minCrystalHits, _maxCrystalHits;
    float _minClusterEDep, _maxTrackChi2, _momCut, _maxDt, _maxEoP;
    
    bool passCH = false;
    bool passEdep = false;
    bool passChi2ndf = false;
    bool passMomCut  = false;
    bool passtrackerHit = false;
    bool findData(const art::Event& evt);
  };

  IPABasicFilter::IPABasicFilter(const Parameters& conf):
    art::EDFilter{conf},
    _diagLevel(conf().diagLevel()),
		_mcdiag(conf().mcdiag()),
		_tcmatchTag(conf().tcmatchTag()),
    _minCrystalHits(conf().minCrystalHits()),
 
{}

  bool IPABasicFilter::filter(art::Event& event) {
    std::cout<<"[In IPACaloCalib::filter() ] ... "<<std::endl;
    if(!findData(event)) 
		  throw cet::exception("RECO")<<"No data in  event"<< endl; 

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
          if(cl->size() < _minCrystalHits ) return false;
          
      }
    }
  }
  if(passCH==true ){ 
    std::cout<<"passes ALL"<<std::endl;
    return true;
  }
  else{ return false; }
 }

  

  bool IPABasicFilter::findData(const art::Event& evt){

	
	  _tcmatchcol=0;
	 
	  auto tcmatch = evt.getValidHandle<TrackClusterMatchCollection>(_tcmatchTag);
	  _tcmatchcol = tcmatch.product();
    

	  return _kalrepcol!=0 and _tcmatchcol!=0;
 }


}

using mu2e::IPABasicFilter;
DEFINE_ART_MODULE(IPABasicFilter);
