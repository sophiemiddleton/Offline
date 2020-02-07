// Author: S Middleton
// Purpose: Use Kalman BTrk for the cosmic fit
// Date: Dec 2019
#include "CosmicReco/inc/CosmicKalFit.hh"
#include "CosmicReco/inc/CosmicKalFitData.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "CosmicReco/inc/CosmicTrackFit.hh"
#include "CosmicReco/inc/CosmicTrackFinderData.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
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
#include "CosmicReco/inc/CosmicTrkDef.hh"
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
			explicit CosmicKalFitter(fhicl::ParameterSet const&);
			virtual ~CosmicKalFitter();
			virtual void beginJob() override;
			virtual void beginRun(art::Run& run) override ;
			
		private:
			void produce(art::Event& event ) override;
			int 			        _debug,_diag;
			unsigned _iev;
  			int _printfreq;
			art::ProductToken<CosmicTrackSeedCollection> const  _seedToken;
			art::ProductToken<ComboHitCollection> const  _chToken;
			unsigned				_minnhits;
			std::vector<double> 		_perr;// diagonal parameter errors to use in the fit
			
			double _bz000;        // sign of the magnetic field at (0,0,0)
			HepSymMatrix _hcovar; // cache of parameter error covariance matrix
			bool _saveall;
			
			TrkParticle      _tpart; //TODO --> mu+/-
			CosmicKalFit     _kfit;

			StrawHitFlag     _dontuseflag;
			TrkFitDirection  _fdir;
			double _upz, _downz; 
			
			CosmicKalFitData               _kalResult;
			const CosmicTrackSeedCollection*  _seedcol;
			const ComboHitCollection *_chcol;
			ProditionsHandle<StrawResponse> _strawResponse_h; 
			ProditionsHandle<Mu2eMaterial> _mu2eMaterial_h;
			ProditionsHandle<Mu2eDetector> _mu2eDetector_h;
			const Tracker* tracker;     // straw tracker geometry
			bool findData(const art::Event& e);
			
	};


	CosmicKalFitter::CosmicKalFitter(fhicl::ParameterSet const& pset) :
    	art::EDProducer{pset},
    	_debug(pset.get<int>("debugLevel", 1)),
    	_diag(pset.get<int>("diagLevel",1)),
	_printfreq(pset.get<int>("printFrequency",101)),	
	_seedToken{consumes<CosmicTrackSeedCollection> (pset.get<art::InputTag>("CosmicTrackSeedCollection"))},
	_chToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
	_minnhits(pset.get<unsigned>("minnhits",1)),
	_perr(pset.get<vector<double> >("ParameterErrors")),
	_saveall(pset.get<bool>("saveall", false)),
	_tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::mu_minus))),//TODO+/-
	_kfit(pset.get<fhicl::ParameterSet>("CosmicKalFit",fhicl::ParameterSet())),
	_upz(pset.get<double>("UpstreamZ",-1500)),
    	_downz(pset.get<double>("DownstreamZ",1500))
	{
		//consumes<CosmicTrackSeedCollection>(_seedToken);
		produces<CosmicKalSeedCollection>();
		if(_perr.size() != CosmicLineTraj::NHLXPRM){
			throw cet::exception("RECO")<<"mu2e::CosmicKalFit: parameter error vector has wrong size"<<std::endl;
		}
		_hcovar = HepSymMatrix(CosmicLineTraj::NHLXPRM,0);
		for(size_t ipar = 0; ipar < CosmicLineTraj::NHLXPRM; ++ipar){
			_hcovar(ipar+1,ipar+1) = _perr[ipar]*_perr[ipar];//sets it to all 1's?   
		}
	}

	CosmicKalFitter::~CosmicKalFitter(){}

  	void CosmicKalFitter::beginJob() {
		cout<<"Beginning the Kal Fit module "<<endl;
   		art::ServiceHandle<art::TFileService> tfs;
  	}

  	void CosmicKalFitter::beginRun(art::Run& run) {
		cout<<"Beginning Kal Fit run ... "<<endl;
		mu2e::GeomHandle<mu2e::Tracker> th;
		mu2e::GeomHandle<BFieldManager> bfmgr;
		mu2e::GeomHandle<DetectorSystem> det;
		const Tracker* tracker = th.get();
		_kfit.setTracker  (tracker);
		_mu2eMaterial_h.get(run.id());

		// change coordinates to mu2e
		CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
		CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
		CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
		_bz000    = field.z();
		

   	}

  	void CosmicKalFitter::produce(art::Event& event ) {
		std::cout<<"CosmicKal produce "<<std::endl;
		auto _srep = _strawResponse_h.getPtr(event.id());
		auto detmodel = _mu2eDetector_h.getPtr(event.id());
		unique_ptr<CosmicKalSeedCollection> kal_col(new CosmicKalSeedCollection());
std::cout<<"ks ptr"<<std::endl;
		_iev=event.id().event();
		
		if(_debug > 0 && (_iev%_printfreq)==0)cout<<"CosmicKalFit: event="<<_iev<<endl;
		
		if(!findData(event)){
      			throw cet::exception("RECO")<<"mu2e::CosmicKalFit: data missing or incomplete"<< endl;
		
    		}
		
		_kalResult.event   = &event;
		_kalResult.seedcol  = _seedcol; 
		_kalResult.fdir = _fdir; //TODO - we dont care about fitdirection for Cosmics?
		_kalResult.chcol  = _chcol; 
		
		for (size_t index=0;index< _kalResult.seedcol->size();++index) {
			
			CosmicTrackSeed const& sts(_kalResult.seedcol->at(index));
			ComboHitCollection _strawCHcol = sts._straw_chits ;

			TrkParticle tpart(_tpart); //TODO we want either +/-13
			TrkParticle::type t = (TrkParticle::type) fabs(((-(int) _tpart.particleType())));
			tpart = TrkParticle(t);

      			if (sts.track().converged == true ) {
				if(_diag> 1){
					std::cout<<"[CosmicKalFitter::produce] seed has converged seed track: "<<index<<std::endl;
				}
				std::vector<CosmicKalSeed>   kal_vec;
				CosmicKalFitData _tmpResult(_kalResult);
				HepVector hpvec(CosmicLineTraj::NHLXPRM);
				float AM = sts.track().AMSIGN();
				if(CosmicTrkUtils::CosmicTrack2Traj(sts.track(), hpvec,AM))				{					
					CosmicLineTraj cosmictraj(hpvec,_hcovar);
					if(_debug > 1) {
		  				cout << "Seed Fit HelixTraj parameters " << cosmictraj.parameters()->parameter()
		       				<< "and covariance " << cosmictraj.parameters()->covariance() <<  endl;
					}
					
					
					TimeCluster tclust;
					tclust._t0 = sts._t0;
					for(uint16_t ihit=0;ihit < sts.hits().size(); ++ihit){
						sts.hits().fillStrawHitIndices(event,ihit,tclust._strawHitIdxs);
					}	
					CosmicTrkDef seeddef(tclust,cosmictraj,tpart,_fdir); 
					const CosmicLineTraj* traj = &seeddef.cosmic();
					double flt0  = traj->zFlight(0.0,traj->z0()); 
					double mom   = TrkMomCalculator::vecMom(*traj, _kfit.bField(), flt0).mag(); //FIXME - mom never set....
					double vflt  = seeddef.particle().beta(mom)*CLHEP::c_light;
					double  cosmict0 = sts.t0().t0();
					cout<<"vflt "<<vflt<<" mom "<<mom<<endl;
					CosmicKalSeed kf(tpart,_fdir, sts.t0(), flt0, sts.status());
					auto cosH = event.getValidHandle(_seedToken);
					kf._cosmicseed = art::Ptr<CosmicTrackSeed>(cosH,index);
					
					int nsh = sts._shits().size();
					
					for (int i=0; i< nsh; ++i){
	  					size_t i_straw   = sts._shits().at(i);
          					
						//const ComboHit& strawhit(_tmpResult.chcol->at(i_straw));
						//const Straw& straw = tracker->getStraw(strawhit.strawId());
	 					double  fltlen  = traj->zFlight(1.0, traj->z0());//straw.getMidPoint().z(), traj->z0());
						
	  					double  propTime = (fltlen-flt0)/0.062;

	 					TrkStrawHitSeed tshs;
						  tshs._index  = i_straw;
						  tshs._t0     = TrkT0(cosmict0 + propTime, sts._t0.t0Err());
						  tshs._trklen = fltlen; 
						  kf._hits.push_back(tshs);
					}
					if(kf._hits.size() >= _minnhits) kf._status.merge(TrkFitFlag::hitsOK);
					
					if(traj !=0){
						KalSegment kseg;
						BbrVectorErr momerr;
						CosmicTrkUtils::fillSegment(*traj,momerr,0.0, kseg);
						kf._segments.push_back(kseg);
						if(_diag> 1){
							std::cout<<"[CosmicKalFitter::produce] filled segments: "<<index<<std::endl;
						}
					}else {
						throw cet::exception("RECO")<<"mu2e::CosmicKalFitter: Cant extract comsic traj from seed "<<endl;
					}

					_tmpResult.cosmicKalSeed = &kf;	
					if(_diag> 1){
						std::cout<<"[CosmicKalFitter::produce] producing track: "<<index<<std::endl;
				}	      
					_kfit.MakeTrack( _srep, detmodel, _tmpResult);  
					if(_tmpResult.krep != 0 && (_tmpResult.krep->fitStatus().success() || _saveall)){ 
						cout<<"has Rep "<<endl;
						//copied from KalSeedFit: create seed from the fit
						CosmicKalSeed kseed(_tmpResult.krep->particleType(),_fdir,_tmpResult.krep->t0(),_tmpResult.krep->flt0(),kf.status());
						kseed._status.merge(TrkFitFlag::kalmanOK);

						auto cosH = event.getValidHandle(_seedToken);
						kseed._cosmicseed = art::Ptr<CosmicTrackSeed>(cosH,index);
						cout<<"fill SH's "<<endl;	CosmicTrkUtils::fillStrawHitSeeds(_tmpResult.krep,_strawCHcol,kseed._hits);
						if(_tmpResult.krep->fitStatus().success())kseed._status.merge(TrkFitFlag::kalmanOK);
						if(_tmpResult.krep->fitStatus().success()==1)kseed._status.merge(TrkFitFlag::kalmanConverged);
						if(kseed._hits.size() >= _minnhits)kseed._status.merge(TrkFitFlag::hitsOK);
						kseed._chisq = _tmpResult.krep->chisq();
						kseed._fitcon = _tmpResult.krep->chisqConsistency().significanceLevel();
						cout<<"test point 1"<<endl;
						double locflt;
						const CosmicLineTraj* htraj = dynamic_cast<const CosmicLineTraj*>(_tmpResult.krep->localTrajectory(_tmpResult.krep->flt0(),locflt));

						 if(htraj != 0){
						    	KalSegment kseg;
						    	BbrVectorErr momerr = _tmpResult.krep->momentumErr(_tmpResult.krep->flt0());
						    cout<<"test point 2"<<endl;	CosmicTrkUtils::fillSegment(*htraj,momerr,locflt-_tmpResult.krep->flt0(),kseg);
						    	double upflt(0.0), downflt(0.0);

				  		    	static const unsigned maxniter(10);
							double dz;
							unsigned niterUp = 0;
							unsigned niterDown =0;
							double tol =0.1;
							bool Up = true;
							bool Down = true;
							do {
								double ztraj = traj->position(upflt).z();
								dz = _upz - ztraj;
								upflt += dz/traj->direction(upflt).z();
								niterUp++;
								if( fabs(dz) > tol) Up = false;
							}while(fabs(dz) > tol && niterUp < maxniter && Up );
						cout<<"test point 3"<<endl;
							//This is : TrkHelixUtils::findZFltlen(*htraj,_upz,upflt);
						    	do {
								double ztraj = traj->position(downflt).z();
								dz = _downz - ztraj;
								downflt += dz/traj->direction(downflt).z();
								niterDown++;
								if( fabs(dz) > tol) Down = false;
							} while(fabs(dz) > tol && niterDown < maxniter && Down);
							cout<<"test point 4"<<endl;
							//TrkHelixUtils::findZFltlen(*htraj,_downz,downflt);*/
							if(_fdir == TrkFitDirection::downstream){
								kseg._fmin = upflt;
								kseg._fmax = downflt;
							}else{
								kseg._fmax = upflt;
								kseg._fmin = downflt;
							}
							kseed._segments.push_back(kseg);
			    				}
							cout<<"test point 4"<<endl;
							kal_vec.push_back(kseed);
							CosmicKalSeedCollection* col = kal_col.get();
							if (kal_vec.size() == 0)     continue;
						      	col->push_back(kseed);  
			    		}
		      			cout<<"test point 5"<<endl;
		    			_tmpResult.deleteTrack();
					}
        		}
    		}	
   		
  		event.put(std::move(kal_col));    
  	} 

	bool CosmicKalFitter::findData(const art::Event& evt){
		_seedcol = 0;
		_chcol   = 0;
		auto cossH = evt.getValidHandle(_seedToken);
		_seedcol = cossH.product();
		auto shH = evt.getValidHandle(_chToken);
   		_chcol = shH.product();
		return _seedcol != 0 and _chcol!=0 ;
	}

}
using mu2e::CosmicKalFitter;
DEFINE_ART_MODULE(CosmicKalFitter);
