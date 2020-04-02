// Author : S. middleton 
// Date: March 2020
// Purpose: IPA Analysis
#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"

#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

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

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

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

// ROOT incldues
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"

#include "Rtypes.h"
#include "TApplication.h"
#include "TArc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TNtuple.h"

#include "TStyle.h"
#include "TText.h"
#include "TRotMatrix.h"
#include "TColor.h"
#include "TLorentzVector.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;

#define me 0.511 //MeV

namespace mu2e {

  class IPAMatchAna : public art::EDAnalyzer {
     public:
	struct Config {
		     using Name=fhicl::Name;
		     using Comment=fhicl::Comment;
	 	     fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
		     fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
		     fhicl::Atom<art::InputTag> kalrepTag{Name("KalRepPtrCollection"),Comment("outcome of Kalman filter (for tracker momentum info)")};
		     fhicl::Atom<art::InputTag> tcmatchTag{Name("TrackClusterMatchCollection"), Comment("track calo match"), "TrackCaloMatching"};
          fhicl::Atom<int> minTrackerHits {Name("mintrackerhits"), Comment("minimum number of straw hits in tracker "),25};//From Pasha's Study
          fhicl::Atom<int> minCrystalHits{Name("mincrystalhits"), Comment("minimum number of crystal hits "), 2};
		      fhicl::Atom<float> minClusterEDep {Name("minclusterEdep"), Comment("minimum amount of energy deposited "),2};
          fhicl::Atom<float>  minTrackChi2 {Name("maxtrackchi2"), Comment("minimum allowed chi2 "),2.5};
	    	};
		typedef art::EDAnalyzer::Table<Config> Parameters;

		explicit IPAMatchAna(const Parameters& conf);
		virtual ~IPAMatchAna() {}

    
		virtual void beginJob();
		virtual void endJob();
		virtual void analyze(const art::Event& e) override;

     	private:
		std::ofstream outputfile;
		Config _conf;
		int _diagLevel;
		int _mcdiag;

		art::InputTag _kalrepTag;
		art::InputTag _tcmatchTag;
    
		
		const KalRepPtrCollection* _kalrepcol;		
		const TrackClusterMatchCollection* _tcmatchcol;

		TTree* _Ntup;

		Int_t   _nEvents = 0;
		Int_t _evt, _run, _nTracks, _nMatches, _nTrackMatched;
		Float_t _TrackT0, _TrackT0Err, _TrackMom, _MaxEoP, _TrackBackTime , _TrackBackOmega ,_TrackBackD0 , _TrackBackZ0, _TrackBackPhi0, _TrackBackTanDip, _TrackChi2, _TrackChi2DOF, _TrackCosTheta, _TrackEnergy, _TrackEoP;

		Float_t _matchChi2, _matchEDep, _matchPosXCl, _matchPosYCl, _matchPosZCl, _matchPathLen, _matchR, _matchDt, _matchPosXtrk, _matchPosYtrk, _matchPosZtrk,_matchTtrk, _matchClusterSize, _CaloEoP, _matchcrySum;

    int _minTrackerHits, _minCrystalHits;
    float _minClusterEDep, _minTrackChi2;

    Float_t _EoP_diff, _E_diff;
    bool passCH = false;
    bool passEdep = false;
		bool findData(const art::Event& evt);

	};

   IPAMatchAna::IPAMatchAna(const Parameters& conf):
		art::EDAnalyzer(conf),
		_diagLevel(conf().diagLevel()),
		_mcdiag(conf().mcdiag()),
		_kalrepTag(conf().kalrepTag()),		
		_tcmatchTag(conf().tcmatchTag()),
    _minTrackerHits(conf().minTrackerHits()),
    _minCrystalHits(conf().minCrystalHits()),
    _minClusterEDep(conf().minClusterEDep()),
    _minTrackChi2(conf().minTrackChi2())
    {}

  void IPAMatchAna::beginJob(){
	//std::cout<<"[In BeginJob()] Beginning ..."<<std::endl;
  art::ServiceHandle<art::TFileService> tfs;
  _Ntup  = tfs->make<TTree>("CaloMatchAna", "CaloMatchAna");
  _Ntup->Branch("evt",          	&_evt ,        "evt/I");
  _Ntup->Branch("run",          	&_run ,        "run/I");

  _Ntup->Branch("nTracks",	&_nTracks,		"nTracks/I");
  _Ntup->Branch("TrackT0", 	&_TrackT0, 		"TrackT0/F");
  _Ntup->Branch("TrackT0Err", 	&_TrackT0Err,		"TrackT0Err/F");
  _Ntup->Branch("TrackBackTime", 	&_TrackBackTime,   	"TrackBackTime/F"); 
  _Ntup->Branch("TrackBackOmega", &_TrackBackOmega,	"TrackBackOmega/F");
  _Ntup->Branch("TrackBackD0",  	&_TrackBackD0, 		"TrackBackD0/F"); 	
  _Ntup->Branch("TrackBackZ0", 	&_TrackBackZ0, 		"TrackBackZ0/F");
  _Ntup->Branch("TrackBackPhi0",	&_TrackBackPhi0,	"TrackBackPhi0/F");
  _Ntup->Branch("TrackBackTanDip",&_TrackBackTanDip,	"TrackBackTanDip/F");
  _Ntup->Branch("TrackChi2",	&_TrackChi2,		"TrackChi2/F");
  _Ntup->Branch("TrackChi2DOF", 	&_TrackChi2DOF, 	"TrackCho2DOF/F");
  _Ntup->Branch("TrackCosTheta",  &_TrackCosTheta,	"TrackCosTheta/F");
  _Ntup->Branch("TrackEnergy",  	&_TrackEnergy,		"TrackEnergy/F");
  _Ntup->Branch("TrackMom",	&_TrackMom, 		"TrackMom/F");
  _Ntup->Branch("TrackEoP",		&_TrackEoP,			"TrackEoP/F");

  _Ntup->Branch("nMatches",	&_nMatches,		"nMatches/I");	
  _Ntup->Branch("matchPosXtrk",	&_matchPosXtrk,		"matchPosXtrk/F");
  _Ntup->Branch("matchPosYtrk",	&_matchPosYtrk,		"matchPosYtrk/F");
  _Ntup->Branch("matchPosZtrk",	&_matchPosZtrk,		"matchPosZtrk/F");
  _Ntup->Branch("matchTtrk",	&_matchTtrk,		"matchTtrk/F");
  _Ntup->Branch("matchChi2", 	&_matchChi2, 		"matchChi/F");
  _Ntup->Branch("matchEDep", 	&_matchEDep, 		"matchEDep/F");
  _Ntup->Branch("matchPosXcl", 	&_matchPosXCl, 		"matchPosXCl/F");
  _Ntup->Branch("matchPosYCl", 	&_matchPosYCl, 		"matchPosYCl/F");
  _Ntup->Branch("matchPosZCl", 	&_matchPosZCl, 		"matchPosZCl/F"); 
  _Ntup->Branch("matchPathLen", 	&_matchPathLen,	        "matchPathLen/F");
  _Ntup->Branch("matchR",		&_matchR, 		"matchR/F");
  _Ntup->Branch("matchDt",	&_matchDt,		"matchDt/F");
  _Ntup->Branch("matchClusterSize",	&_matchClusterSize,		"matchClusterSize/F");
  _Ntup->Branch("CaloEoP", 	&_CaloEoP, 		"CaloEoP/F");
  _Ntup->Branch("matchcrySum", 	&_matchcrySum, 		"_matchcrySum/F");

  _Ntup->Branch("EoP_diff",	&_EoP_diff,		"_EoP_diff/F");
  _Ntup->Branch("E_diff",	&_E_diff,		"_E_diff/F");
  
  outputfile.open("IPAAnaMatchedTrackClusters.csv");
        outputfile<<"event,run,size,CaloEoP,CaloE,TrackerEoP,TrackerE,P"<<std::endl;
  }


  void IPAMatchAna::analyze(const art::Event& event) {
	//std::cout<<"[In Analyze()] Beginning ..."<<std::endl;
	_evt = event.id().event();
	_run = event.run();

	if(!findData(event)) 
		throw cet::exception("RECO")<<"No data in  event"<< endl; 

	art::ServiceHandle<mu2e::GeometryService>   geom;
	const mu2e::Calorimeter* bc(nullptr);
	if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
    		mu2e::GeomHandle<mu2e::DiskCalorimeter> h;
    		bc = (const mu2e::Calorimeter*) h.get();
  	}
  mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
  mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;
  _nTracks = 0;
  _nMatches =0;

  Hep3Vector vd_tt_back = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_Back));
  double     Z      = vd_tt_back.z();

  for(unsigned int i=0;i<_kalrepcol->size();i++){
    int iv = 0;
    art::Ptr<KalRep> const& ptr = _kalrepcol->at(i);
    const KalRep* TrackKrep = ptr.get();
    const CaloCluster* ClosestCluster(nullptr);
    double best_chi2_match(1.e6); //high number to start
    double  ds(10.), s0, s1, s2, z0, z1, z2, dzds, sz, sz1, z01;
    const TrkHitVector* hots = &TrackKrep->hitVector();
    int nh = hots->size();
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

    double EndMom= TrackKrep->momentum(sz1).mag();//TODO
    //enact tracker cuts:
    //if(TrackKrep->chisq()/TrackKrep->nActive() < minTrackChi2) continue
    //minTrackerHits
    //=========================== Add Tracks =============================//
    _TrackT0 = TrackKrep->t0().t0();
    _TrackT0Err = TrackKrep->t0().t0Err();
    _TrackMom = TrackKrep->momentum(sz1).mag();
    _TrackBackTime =   TrackKrep->arrivalTime(sz1);
    HelixParams helx  = TrackKrep->helix(sz1);
    _TrackBackOmega       = helx.omega(); 
    _TrackBackD0       = helx.d0();
    _TrackBackZ0       = helx.z0();
    _TrackBackPhi0     = helx.phi0();
    _TrackBackTanDip   = helx.tanDip(); 
    _TrackEnergy = sqrt(TrackKrep->momentum(sz1).mag()*TrackKrep->momentum(sz1).mag() + me*me);
    double entlen         = std::min(s1,s2);
    CLHEP::Hep3Vector fitmom = TrackKrep->momentum(entlen);
    TLorentzVector  Momentum(fitmom.x(),fitmom.y(),fitmom.z(),0.511);
    _TrackCosTheta = Momentum.CosTheta();
    _TrackChi2 =TrackKrep->chisq();
    _TrackChi2DOF= TrackKrep->chisq()/TrackKrep->nActive();
    _nTracks ++;
    //============================= Add Matches ===============================//
    if(_tcmatchcol->size() ==0) continue;
    for(unsigned int c=0;c<_tcmatchcol->size();c++){

      TrackClusterMatch const& tcm = (*_tcmatchcol)[c];
      const TrkCaloIntersect* extrk = tcm.textrapol();
      const KalRep* Krep  = extrk->trk().get();
      if (Krep == TrackKrep) {
      const mu2e::CaloCluster* cl = tcm.caloCluster();
      
      // enact cluster cuts:
      if(cl->size() < _minCrystalHits) continue;
      passCH = true;
      if(cl->energyDep() < _minClusterEDep) continue;
      passEdep = true;
      _matchClusterSize = cl->size();
      _matchEDep = cl->energyDep();

      _matchcrySum = 0;
      for(unsigned i =0 ; i< cl->caloCrystalHitsPtrVector().size();i++){
			   art::Ptr< CaloCrystalHit>  cry=cl->caloCrystalHitsPtrVector()[i] ;
        _matchcrySum += cry->energyDep();
      }
      iv   = cl->diskId();
      CLHEP::Hep3Vector x1   = bc->geomUtil().mu2eToDisk(iv,cl->cog3Vector());

    if ((ClosestCluster == nullptr) || (tcm.chi2() < best_chi2_match )) {
      ClosestCluster = cl;
	    best_chi2_match    = tcm.chi2();
    }
    _matchPosXtrk = tcm.xtrk();
    _matchPosYtrk = tcm.ytrk();
    _matchPosZtrk = tcm.ztrk();
    _matchTtrk	 = tcm.ttrk();
    _matchChi2 = tcm.chi2();
    _matchPosXCl =x1.x();
    _matchPosYCl =x1.y();
    _matchPosZCl = x1.z();
    _matchPathLen = tcm.ds();
    _matchDt = tcm.dt();
    _matchR = sqrt(x1.x()*x1.x() + x1.y()*x1.y());
    
    _nMatches++;
    }
    }

    //double Ep=ClosestCluster->energyDep();
    if(EndMom!=0) { 
      _CaloEoP = ClosestCluster->energyDep()/EndMom;
      _TrackEoP = _TrackEnergy/EndMom;
      _EoP_diff = _CaloEoP - _TrackEoP;
      _E_diff = ClosestCluster->energyDep() - _TrackEnergy;
    }
    outputfile<<_evt<<","<<_run<<","<<_matchClusterSize<<","<<_CaloEoP<<","<<ClosestCluster->energyDep()<<","<<_TrackEoP<<","<<_TrackEnergy<<","<<_TrackMom<<std::endl;

    _nTrackMatched++;
  }	
  
    _Ntup->Fill();
    _nEvents++;

}



bool IPAMatchAna::findData(const art::Event& evt){

	
	_tcmatchcol=0;
	_kalrepcol = 0;

	auto kalrep = evt.getValidHandle<KalRepPtrCollection>(_kalrepTag);
	_kalrepcol =kalrep.product();
	auto tcmatch = evt.getValidHandle<TrackClusterMatchCollection>(_tcmatchTag);
	_tcmatchcol = tcmatch.product();
  

	return _kalrepcol!=0 and _tcmatchcol!=0;
       }

 void IPAMatchAna::endJob(){
	
  } 

}

DEFINE_ART_MODULE(mu2e::IPAMatchAna);
