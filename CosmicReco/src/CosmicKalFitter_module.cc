// Author: S Middleton
// Purpose: Use Kalman BTrk for the cosmic fit
// Date: Dec 2019
#include "CosmicReco/inc/CosmicKalFit.hh"
#include "CosmicReco/inc/CosmicKalFitData.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "CosmicReco/inc/CosmicTrackFit.hh"
#include "CosmicReco/inc/CosmicTrackFinderData.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "CosmicReco/inc/CosmicTrkMomCalc.hh"
#include "CosmicReco/inc/CosmicTrkUtils.hh"

// ART:
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "GeneralUtilities/inc/Angles.hh"
#include "art/Utilities/make_tool.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/Mu2eMaterial.hh"
#include "TrackerConditions/inc/Mu2eDetector.hh"
// art
#include "canvas/Persistency/Common/Ptr.h"
// MC data
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
//MU2E:
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
//utils:
#include "Mu2eUtilities/inc/ParametricFit.hh"
//For Drift:
#include "TrkReco/inc/PanelAmbigResolver.hh"
#include "TrkReco/inc/PanelStateIterator.hh"
#include "TrkReco/inc/TrkFaceData.hh"
// Mu2e BaBar
#include "BTrkData/inc/TrkStrawHit.hh"
// BaBar
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
//#include "BTrk/TrkBase/TrkCosmicUtils.hh" //TODO
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"

//CLHEP:
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
//C++:
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

using namespace std;
using namespace ROOT::Math::VectorUtil;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;
namespace{
    //create a compare struct to allow height ordering of the hits in an event. 
    struct ycomp : public std::binary_function<mu2e::ComboHit,mu2e::ComboHit,bool> {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1._pos.y() > p2._pos.y(); }
  };
}
namespace mu2e{

  class CosmicKalFitter : public art::EDProducer {
  public:
    struct Config{
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<int> diag{Name("diag"), Comment("set on for MC info"),2};
	      fhicl::Atom<int> debug{Name("debugLevel"), Comment("set to 1 for debug prints"),1};
              fhicl::Atom<art::InputTag> seedToken{Name("CosmicTrackSeedCollection"),Comment("tag for combo hit collection")};
              fhicl::Atom<unsigned> minnhits{Name("minnhits"), Comment("min number of trkstrawhits"), 10};
	      fhicl::Table<CosmicKalFit::Config> kfit{Name("CosmicKalFit"), Comment("fit")};
	      fhicl::Atom<bool> saveall{Name("saveall"), Comment("save all"), false);
    };
    typedef art::EDProducer::Table<Config> Parameters;
    explicit CosmicKalFitter(const Parameters& conf);
    virtual ~CosmicKalFitter();
    virtual void beginJob() override;
    virtual void beginRun(art::Run& run) override ;
    virtual void produce(art::Event& event ) override;
    
  private:
    
    Config _conf;

    int 			        _diag, _debug;
    art::InputTag 		        _seedToken;
    unsigned				_minnhits;
    std::vector<double> 		_perr;
    double _amsign; // cached sign of angular momentum WRT the z axis 
    double _bz000;        // sign of the magnetic field at (0,0,0)
    HepSymMatrix _hcovar; // cache of parameter error covariance matrix

    CosmicKalFit     _kfit;

    TrkParticle      _tpart; //TODO --> mu+/-
    StrawHitFlag     _dontuseflag;
    TrkFitDirection  _fdir;
    double upz, downz; //TODO - may not need
    bool _saveall;
    CosmicKalFitData               _kalResult;
    CosmicTrackSeedCollection	     *_seedcol;
    ProditionsHandle<StrawResponse> _strawResponse_h; 
    ProditionsHandle<Mu2eMaterial> _mu2eMaterial_h;
    ProditionsHandle<Mu2eDetector> _mu2eDetector_h;

};


 CosmicKalFitter::CosmicKalFitter(const Parameters& conf) :
   art::EDProducer(conf),
   	_diag (conf().diag()),
	_debug  (conf().debug()),
	_seedToken (conf().seedToken()),
        _minnhits (conf().minnhits()),
	_kfit (conf().kfit()),
	_saveall(conf().saveall())
{
	    consumes<CosmicTrackSeedCollection>(_seedToken);
	    produces<CosmicKalSeedCollection>();
	     if(_perr.size() != HelixTraj::NHLXPRM)
      		throw cet::exception("RECO")<<"mu2e::CosmicKalFit: parameter error vector has wrong size"<<std::endl;
    	     _hcovar = HepSymMatrix(CosmicLineTraj::NHLXPRM,0);
    		for(size_t ipar = 0; ipar < CosmicLineTraj::NHLXPRM; ++ipar){
      			_hcovar(ipar+1,ipar+1) = _perr[ipar]*_perr[ipar];   
 }
}
  CosmicKalFitter::~CosmicKalFitter(){}

  void CosmicKalFitter::beginJob() {
   art::ServiceHandle<art::TFileService> tfs;
    
  }

  void CosmicKalFitter::beginRun(art::Run& run) {
   mu2e::GeomHandle<mu2e::Tracker> th;
   mu2e::GeomHandle<BFieldManager> bfmgr;
   mu2e::GeomHandle<DetectorSystem> det;
   mu2e::GeomHandle<mu2e::Tracker> th;
 
   const Tracker* tracker = th.get();
   _kalResult.run = &run;
   _kfit.setTracker  (tracker);
   _mu2eMaterial_h.get(run.id());

   // change coordinates to mu2e
    CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
    CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
    _bz000    = field.z();
   
   }

  void CosmicKalFitter::produce(art::Event& event ) {
     int _iev=event.id().event();

     auto _srep = _strawResponse_h.getPtr(event.id());
     auto detmodel = _mu2eDetector_h.getPtr(event.id());

     unique_ptr<CosmicKalSeedCollection> kal_col(new CosmicKalSeedCollection());
     
     auto const& seedH = event.getValidHandle<CosmicTrackSeedCollection>(_seedToken);
     const CosmicTrackSeedCollection& seedcol(*seedH);

     _kalResult.event   = &event;
     _kalResult.seedcol  = seedcol; 
     _kalResult.fdir = _fdir;
     
     for (size_t index=0;index< seedcol.size();++index) {
       
       CosmicKalSeed kalseed;
       _kalResult._kalseed            = kalseed;
      //_kalResult.kalseed._chcol    =&
      _kalResult._kalseed._status.merge(TrkFitFlag::kalmanOK);

      CosmicTrackSeed sts =(*seedcol)[ist];
      TrkParticle tpart(_tpart);
      TrkParticle::type t = (TrkParticle::type)(-(int) _tpart.particleType());//TODO - ubnderstand how to pass this in in new fhicil format
      tpart = TrkPatricle t;
      if (sts.track().converged == true ) {
	      std::vector<CosmciKalSeed>   kal_vec;
	      CosmicKalFitData tmpResult(_kalResult);
	      HepVector hpvec(CosmicLineTraj::NHLXPRM);
              CosmicTrack2Traj();//IMPORTANT: This is were the switch happens!!
	      CosmicLineTraj cosmictraj(hpvec,_hcovar);
	      TimeCluster tclust;
	      tclust._t0 = sts._t0;
	      for(uint16_t ihit=0;ihit < sts.hits().size(); ++ihit){
	  	ComboHit const& ch = sts.hits()[ihit];
	  	sts.hits().fillStrawHitIndices(event,ihit,tclust._strawHitIdxs);
		}
	       TrkDef seeddef(tclust,cosmictraj,tpart,_fdir); //FIXME - is this needed?
	       const CosmicTraj* traj = &seeddef.cosmic();
	       double flt0  = cosmictraj->zFlight(0.0); //FIXME - does this work?
	       double mom   = CosmicTrkMomCalc::vecMom(*cosmictraj, _kfit.bField(), flt0).mag(); //TODO check
		double vflt  = seeddef.particle().beta(mom)*CLHEP::c_light;
	        double  cosmict0 = sts.t0().t0();
		CosmicKalSeed kf(tpart,_fdir, sts.t0(), flt0, sts.status());
	        auto cosH = event.getValidHandle(_hsToken);
	        kf._cosmicseed = art::Ptr<CosmicTrackSeed>(cosH,iseed);
	        std::vector<TrkStrawHitSeed>const trkseeds = sts.trkstrawhits();
              	for(auto const& ths : trkseeds ){
      	             size_t index = ths.index();
		     const ComboHit& strawhit(sts._straw_chits.at(index));
		     const Straw& straw = _tracker->getStraw(strawhit.strawId());
                     double fltlen = cosmictraj->zFlight(straw.getMidPoint().z());
		     double propTime = (fltlen-flt0)/vflt;
		     TrkStrawHit* trkhit = new TrkStrawHit(_srep, strawhit, straw, ths.index(),ths.t0(),fltlen, 5.,1.); //maxpull set to 5, time weight set to 1.
		}  
	}//end "is converged"
				
	if(kf._hits.size() >= _minnhits) kf._status.merge(TrkFitFlag::hitsOK);

	if(!cosmictraj !=0){
		KalSegment kseg;
		BbrVectorErr = momerr;
		CosmicTrkUtils::fillSegment(*cosmictraj,momerr,0.0, kseg);
		kf._segments.push_back(kseg);
	}else {
		throw cet::exception("RECO")<<"mu2e::CosmicKalFitter: Cant extract comsic traj from seed "<<endl;
	}

	_tmpResult.kalSeed = &kf;		      
	_kfit.MakeTrack( _srep, detmodel, tmpResult);  
	if(_tmpResult.krep != 0 && (_tmpResult.krep->fitStatus().success() || _saveall)){ 
	 
          //copied from KalSeedFit: create seed from the fit
	  KalSeed kseed(_tmpResult.krep->particleType(),_fdir,_tmpResult.krep->t0(),_tmpResult.krep->flt0(),kf.status());
	  kseed._status.merge(_ksf);

          auto cosH = event.getValidHandle(_seedToken);
	  kseed._cosmicSeed = art::Ptr<CosmicTrackSeed>(cosH,index);
          //extract hits from the rep and put into CosmicKalSeed
          CosmicTrkUtils::fillStrawHitSeeds(_tmpResult.krep,*_chcol,kseed._hits);
	  if(_tmpResult.krep->fitStatus().success())kseed._status.merge(TrkFitFlag::kalmanOK);
	  if(_tmpResult.krep->fitStatus().success()==1)kseed._status.merge(TrkFitFlag::kalmanConverged);
	  if(kseed._hits.size() >= _minnhits)kseed._status.merge(TrkFitFlag::hitsOK);
	  kseed._chisq = _tmpResult.krep->chisq();
	  kseed._fitcon = _tmpResult.krep->chisqConsistency().significanceLevel();
	
	  double locflt;
	  const CosmicLineTraj* htraj = dynamic_cast<const CosmicLineTraj*>(_tmpResult.krep->localTrajectory(_tmpResult.krep->flt0(),locflt));
	 
	  if(htraj != 0){
	    KalSegment kseg;
	    BbrVectorErr momerr = _tmpResult.krep->momentumErr(_tmpResult.krep->flt0());
	    CosmicTrkUtils::fillSegment(*htraj,momerr,locflt-_tmpResult.krep->flt0(),kseg);
	    double upflt(0.0), downflt(0.0);
	    //TrkHelixUtils::findZFltlen(*htraj,_upz,upflt);//TODO
	    //TrkHelixUtils::findZFltlen(*htraj,_downz,downflt);//TODO
	    if(_fdir == TrkFitDirection::downstream){
		kseg._fmin = upflt;
		kseg.fmax = downflt;
	    }else{
		kseg._fmax = upflt;
		kseg._fmin = downflt;
            }
		kseed._segments.push_bck(kseg);
	    }
	    tmpResult._kalseed = kseed;
	    kal_vec.push_back(kseed);
	    CosmicKalSeedCollection* col = kal_col.get();
	    iff (kal_vec.size() == 0)     continue;
	      col->push_back(kseed);  
			          
              }
	    tmpResult.deleteTrack();
        }
    }
 
  event.put(std::move(kal_col));    
  } //end produce

}//end namespace

using mu2e::CosmicKalFitter;
DEFINE_ART_MODULE(CosmicKalFitter);
