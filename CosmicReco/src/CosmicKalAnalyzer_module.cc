//Author: S Middleton
//Date: March 2020
//Purpose: Analyzer for Kalman Cosmics

#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "CosmicReco/inc/CosmicTrackFit.hh"
#include "CosmicReco/inc/CosmicTrackFinderData.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "CosmicReco/inc/CosmicTrackMCInfo.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"

//Cosmic Kalman:
#include "CosmicReco/inc/CosmicKalFit.hh"
#include "CosmicReco/inc/CosmicKalFitData.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"

//Mu2e Data Prods:
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/XYZVec.hh"

//Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
#include "CosmicReco/inc/DriftFitUtils.hh"
#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "TrackerConditions/inc/StrawResponse.hh"

// Mu2e diagnostics
#include "TrkDiag/inc/ComboHitInfo.hh"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

// ROOT incldues
#include "TStyle.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TProfile.h"

//Geom
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"

using namespace std; 

namespace mu2e 
{
	class CosmicKalAnalyzer : public art::EDAnalyzer {
	public:
	struct Config{
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<int> diag{Name("diagLevel"), Comment("set to 1 for info"),1};
	      fhicl::Atom<bool> mcdiag{Name("mcdiag"), Comment("set on for MC info"),true};
	      fhicl::Atom<art::InputTag> kaltag{Name("CosmicKalSeedCollection"),Comment("tag for cosmci track seed collection")};
	      fhicl::Atom<art::InputTag> mcdigistag{Name("StrawDigiMCCollection"),Comment("StrawDigi collection tag"),"makeSD"};
	      fhicl::Table<SimParticleTimeOffset::Config> toff{Name("TimeOffsets"), Comment("Sim particle time offset ")};
	    };
	typedef art::EDAnalyzer::Table<Config> Parameters;

	explicit CosmicKalAnalyzer(const Parameters& conf);
	virtual ~CosmicKalAnalyzer();
	virtual void beginJob() override;
	virtual void beginRun(const art::Run& r) override;
	virtual void analyze(const art::Event& e) override;
	
    private: 
      
	Config _conf;

	int  _diag;
	bool _mcdiag;
	std::ofstream outputfile;
	art::InputTag   _kaltag;
	art::InputTag   _mcdigistag; //MC Digis
	SimParticleTimeOffset _toff;

	const CosmicKalSeedCollection* _kalcol;
	const StrawDigiMCCollection* _mcdigis;
	CosmicTrackMCInfo trueinfo;

	TTree* _cosmic_tree;
	
	Float_t _Particle;
	Float_t _Direction;
	Float_t _Chi2;
	Float_t _T0;
  Float_t _ChiCon;
  Float_t _d0;
  Float_t _z0;
  Float_t _phi0;
  Float_t _theta;

	Int_t _evt; 
	int _nused;
        Bool_t _KalmanConverged;

	ProditionsHandle<Tracker> _alignedTracker_h;
	ProditionsHandle<StrawResponse> _strawResponse_h; 

	Int_t _strawid; 
	vector<ComboHitInfoMC> _chinfomc;

	CosmicTrackMCInfo FitMC(const StrawDigiMCCollection*& _mcdigis);
	CosmicTrackMCInfo FillDriftMC(ComboHit const& chi, double reco_ambig, CosmicTrackMCInfo info, double t0, const Tracker* tracker);
      	bool findData(const art::Event& evt);

    };

    CosmicKalAnalyzer::CosmicKalAnalyzer(const Parameters& conf) :
	art::EDAnalyzer(conf),
	_diag (conf().diag()),
	_mcdiag (conf().mcdiag()),
	_kaltag (conf().kaltag()),
	_mcdigistag (conf().mcdigistag()),
	_toff (conf().toff())
	{
      		if(_mcdiag){
			for (auto const& tag : conf().toff().inputs()) {
				consumes<SimParticleTimeMap>(tag);
			}
		}
       }

    CosmicKalAnalyzer::~CosmicKalAnalyzer(){}

    void CosmicKalAnalyzer::beginJob() {
    
	if(_diag > 0){

		art::ServiceHandle<art::TFileService> tfs;
		_cosmic_tree=tfs->make<TTree>("cosmic_tree"," Diagnostics for Cosmic Track Fitting");
		_cosmic_tree->Branch("nused",  &_nused ,   "nused/I");
		_cosmic_tree->Branch("evt",&_evt,"evt/I");  
		_cosmic_tree->Branch("Patricle", &_Particle,"Particle/F");
		_cosmic_tree->Branch("Direction", &_Direction,"Direction/F");
		_cosmic_tree->Branch("Chi2", &_Chi2, "Chi2/F");
		_cosmic_tree->Branch("T0", &_T0, "T0/F");
    _cosmic_tree->Branch("ChiCon", &_ChiCon, "ChiCon/F");
		_cosmic_tree->Branch("d0", &_d0, "d0/F");
    _cosmic_tree->Branch("z0", &_z0, "z0/F");
    _cosmic_tree->Branch("phi0", &_phi0, "phi0/F");
    _cosmic_tree->Branch("theta", &_theta, "theta/F");
	}
      }

  void CosmicKalAnalyzer::beginRun(const art::Run& run){}

  void CosmicKalAnalyzer::analyze(const art::Event& event) {

	  _evt = event.id().event();  

	  if(!findData(event)) { throw cet::exception("RECO")<<"No Time Clusters in event"<< endl; }

    for(size_t ist = 0;ist < _kalcol->size(); ++ist){

		  CosmicKalSeed sts =(*_kalcol)[ist];
		  double t0 = sts._t0.t0();
		  TrkFitFlag const& status = sts._status;
      const CosmicLineTraj *traj = sts._traj;
		  if (status.hasAllProperties(TrkFitFlag::kalmanConverged) ){ 
		    //Particle = sts._tpart;
		    //Direction = sts._fdir;
		    _Chi2 = sts._chisq;
        _ChiCon = sts._fitcon;
		    _T0 = t0;
        _d0 = traj->d0();
        _z0 = traj->z0();
        _phi0 = traj->phi0();
        _theta = traj->theta();

		    _nused ++;
		
	
	
 		    _cosmic_tree->Fill();
      }
  }
      
}

  bool CosmicKalAnalyzer::findData(const art::Event& evt){
	_kalcol = 0; 
	auto stH = evt.getValidHandle<CosmicKalSeedCollection>(_kaltag);
	_kalcol =stH.product();
	if(_mcdiag){
		_mcdigis=0;
		auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
		_mcdigis = mcdH.product();
		_toff.updateMap(evt);
	}
	return _kalcol !=0 && (_mcdigis != 0 || !_mcdiag);
       }

}

using mu2e::CosmicKalAnalyzer;
DEFINE_ART_MODULE(CosmicKalAnalyzer);
