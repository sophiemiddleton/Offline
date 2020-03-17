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
	      fhicl::Atom<art::InputTag> costag{Name("CosmicKalSeedCollection"),Comment("tag for cosmci track seed collection")};
	      fhicl::Atom<art::InputTag> mcdigistag{Name("StrawDigiMCCollection"),Comment("StrawDigi collection tag"),"makeSD"};
	      fhicl::Table<SimParticleTimeOffset::Config> toff{Name("TimeOffsets"), Comment("Sim particle time offset ")};
	    };
	typedef art::EDAnalyzer::Table<Config> Parameters;

	explicit CosmicKalAnalyzer(const Parameters& conf);
	virtual ~CosmicKalAnalyzer();
	virtual void beginJob() override;
	virtual void beginRun(const art::Run& r) override;
	virtual void analyze(const art::Event& e) override;
	virtual void endJob() override;
    private: 
      
	Config _conf;

	int  _diag;
	bool _mcdiag;
	std::ofstream outputfile;
	art::InputTag   _costag;
	art::InputTag   _mcdigistag; //MC Digis
	SimParticleTimeOffset _toff;

	const CosmicKalSeedCollection* _coscol;
	const StrawDigiMCCollection* _mcdigis;
	CosmicTrackMCInfo trueinfo;

	TTree* _cosmic_tree;
	
	Float_t _Particle;
	Float_t _Direction;
	Float_t _Chi2;
	Float_t _T0;

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
	_costag (conf().costag()),
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
		
	}
      }

    void CosmicKalAnalyzer::beginRun(const art::Run& run){}

    void CosmicKalAnalyzer::analyze(const art::Event& event) {

	_evt = event.id().event();  

	if(!findData(event)) 
		throw cet::exception("RECO")<<"No Time Clusters in event"<< endl; 

        for(size_t ist = 0;ist < _coscol->size(); ++ist){

		CosmicKalSeed sts =(*_coscol)[ist];
		double t0 = sts._t0.t0();
		TrkFitFlag const& status = sts._status;

		if (!status.hasAllProperties(TrkFitFlag::helixOK) ){ continue; }
		//Particle = sts._tpart;
		//Direction = sts._fdir;
		_Chi2 = sts._chisq;
		_T0 = t0;
		if (status.hasAllProperties(TrkFitFlag::kalmanConverged) ){ 
			_nused ++;
		}
	}
	
 		_cosmic_tree->Fill();
}
      



    void CosmicKalAnalyzer::endJob() {}

    CosmicTrackMCInfo CosmicKalAnalyzer::FitMC(const StrawDigiMCCollection*& _mcdigis){	
	
	::BuildLinearFitMatrixSums S; 
	CosmicTrackMCInfo TrackTrueInfo;

	StrawDigiMC hitP1; 
	StrawDigiMC first = (*_mcdigis)[0];

	auto const& spmcp0= first.earlyStrawGasStep();
	XYZVec pos0(spmcp0->position().x(), spmcp0->position().y(), spmcp0->position().z());
	XYZVec dir0(spmcp0->momentum().x(), spmcp0->momentum().y(), spmcp0->momentum().z());
	
        for(size_t ich = 0;ich < _mcdigis->size(); ++ich){
            hitP1 = (*_mcdigis)[ich];
	    
	    auto const& spmcp= hitP1.earlyStrawGasStep();
            XYZVec posN(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());
            
            XYZVec ZPrime = (spmcp->momentum().Unit());
            
            TrackAxes TrueAxes = ParametricFit::GetTrackAxes(ZPrime);
            TrackTrueInfo.TrueTrackCoordSystem = (TrueAxes);
	    
            XYZVec point(posN.x(), posN.y(), posN.z());
            XYZVec X(1,0,0);
            XYZVec Y(0,1,0);
            XYZVec Z(0,0,1);
            S.addPoint( point, X,Y,Z, 1,1);
            
        }   
    
	TrackParams RawTrueParams(S.GetAlphaX()[0][0], S.GetAlphaX()[1][0], S.GetAlphaY()[0][0], S.GetAlphaY()[1][0]);

	XYZVec TruePos(S.GetAlphaX()[0][0], S.GetAlphaY()[0][0], 0);

	XYZVec TrueDir(S.GetAlphaX()[1][0], S.GetAlphaY()[1][0], 1);
	TrueDir = TrueDir.Unit();
	TrueDir = TrueDir/TrueDir.Z();

	pos0.SetX(pos0.X()-(dir0.X()*pos0.Z()/dir0.Z()));
	pos0.SetY(pos0.Y()-(dir0.Y()*pos0.Z()/dir0.Z()));
	pos0.SetZ(pos0.Z()-(dir0.Z()*pos0.Z()/dir0.Z()));
	dir0 = dir0/dir0.Z();

	TrackEquation TrueTrack(pos0, dir0);

	TrackTrueInfo.TrueFitEquation = (TrueTrack);
	TrackTrueInfo.TruePhi =(atan(TrueDir.y()/TrueDir.x()));
	TrackTrueInfo.TrueTheta = (acos(TrueDir.x()/sqrt(TrueDir.Mag2())));


	int n{-1};
	for(size_t ich = 0;ich < _mcdigis->size(); ++ich){
		StrawDigiMC digimc= (*_mcdigis)[ich];
		auto const& sim_cal = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
		auto const& sim_hv  = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
		if ( sim_cal == sim_hv ){
		  
		  if (n == 0){
			TrackTrueInfo.TrueMomentum = sqrt(digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().mag2());;
		    	TrackTrueInfo.TrueThetaSIM = acos(digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().z()/TrackTrueInfo.TrueMomentum );
		    	TrackTrueInfo.TruePhiSIM = atan(digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().y()/digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().x());
		 }
		 n++;
      }
    }
     
    return TrackTrueInfo;
   }

    CosmicTrackMCInfo CosmicKalAnalyzer::FillDriftMC(ComboHit const& chit, double RecoAmbig, CosmicTrackMCInfo info, double t0,  const Tracker* tracker){
	double true_doca = DriftFitUtils::GetTestDOCA(chit, info.TrueFitEquation.Pos.X(), info.TrueFitEquation.Dir.X(), info.TrueFitEquation.Pos.Y(),info.TrueFitEquation.Dir.Y(),  tracker);
	double trueambig = DriftFitUtils::GetAmbig(chit, info.TrueFitEquation.Pos.X(), info.TrueFitEquation.Dir.X(), info.TrueFitEquation.Pos.Y(),info.TrueFitEquation.Dir.Y(),  tracker);
	info.Ambig.push_back(trueambig);
	info.TrueDOCA.push_back(true_doca);

	return info;
    }

    bool CosmicKalAnalyzer::findData(const art::Event& evt){
	_coscol = 0; 
	auto stH = evt.getValidHandle<CosmicKalSeedCollection>(_costag);
	_coscol =stH.product();
	if(_mcdiag){
		_mcdigis=0;
		auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
		_mcdigis = mcdH.product();
		_toff.updateMap(evt);
	}
	return _coscol !=0 && (_mcdigis != 0 || !_mcdiag);
       }

}

using mu2e::CosmicKalAnalyzer;
DEFINE_ART_MODULE(CosmicKalAnalyzer);
