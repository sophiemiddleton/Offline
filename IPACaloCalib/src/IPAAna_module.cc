/*
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
//#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

//#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"

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

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <fstream>

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

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class IPAAna : public art::EDAnalyzer {
     public:
	struct Config {
		     using Name=fhicl::Name;
		     using Comment=fhicl::Comment;
	 	     fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
		     fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
		     fhicl::Atom<art::InputTag> calocrysTag{Name("CaloCrystalHitCollection"),Comment("cal reco crystal hit info")};
		     fhicl::Atom<art::InputTag> caloclusterTag{Name("CaloClusterCollection"),Comment("cal reco cluster info")};
		     fhicl::Atom<art::InputTag> genTag{Name("GenParticleCollection"), Comment("gen particle info")};
		     fhicl::Atom<art::InputTag> kalrepTag{Name("KalRepPtrCollection"),Comment("outcome of Kalman filter (for tracker momentum info)")};
		     fhicl::Atom<art::InputTag> tcmatchTag{Name("TrackClusterMatchCollection"), Comment("track calo match"), "TrackCaloMatching"};
		     
	    	};
		typedef art::EDAnalyzer::Table<Config> Parameters;

		explicit IPAAna(const Parameters& conf);
		virtual ~IPAAna() {}

		virtual void beginJob();
		virtual void endJob();
		virtual void analyze(const art::Event& e) override;

     	private:
		std::ofstream outputfile;
		Config _conf;
		int _diagLevel;
		int _mcdiag;

		art::InputTag _calocrysTag;
		art::InputTag _caloclusterTag;
		art::InputTag _genTag;
		art::InputTag _kalrepTag;
		art::InputTag _tcmatchTag;

		const KalRepPtrCollection* _kalrepcol;		
		const TrackClusterMatchCollection* _tcmatchcol;
		const CaloCrystalHitCollection*  _calcryhitcol;	
		const CaloClusterCollection* _calclustercol;
		const GenParticleCollection *_gencol;

		TTree* _Ntup;

		Int_t   _nEvents = 0;
		Int_t _evt, _run, _nHits, _nClusters, _nTracks, _nMatches, _nTrackMatched, _nGen;
		
		Int_t _genPdgId, _genCrCode;

		Float_t _genmomX,_genmomY, _genmomZ, _genStartX,  _genStartY, _genStartZ, _genStartT;

		Int_t   _cryId[8129], _crySectionId[8129], _crySimIdx[8129], _crySimLen[8129];
		Float_t _cryTime[8129], _cryEdep[8129],_cryDose[8129], _cryPosX[8129], _cryPosY[8129], _cryPosZ[8129], _cryLeak[8129], _cryTotE[8129], _cryTotSum[8129], _cryTotEErr[8129], _cryRadius[8129],  _cryMaxR[8129];

		Int_t   _clusterId, _clusterdiskId;
		Float_t _clustertime, _clustertimeErr, _clusterEdep, _clusterEDepErr, _clusterangle, _clustercog3VectorX, _clustercog3VectorY, _clustercog3VectorZ,_clusterR,_clusterNHits, _clustermaxECrystal, _clusterindexMaxECrystal, _clusterERatio;

                
		Int_t _evt, _run, _nTracks, _nMatches, _nTrackMatched;
		Float_t _TrackT0, _TrackT0Err, _TrackMom, _MaxEoP, _EoP, _TrackBackTime , _TrackBackOmega ,_TrackBackD0 , _TrackBackZ0, _TrackBackPhi0, _TrackBackTanDip, _TrackChi2, _TrackChi2DOF, _TrackCosTheta, _TrackEnergy, _TrackEoP;

		Float_t _matchChi2, _matchEDep, _matchPosXCl, _matchPosYCl, _matchPosZCl, _matchPathLen, _matchR, _matchDt, _matchPosXtrk, _matchPosYtrk, _matchPosZtrk,_matchTtrk;

		bool findData(const art::Event& evt);

	};

   IPAAna::IPAAna(const Parameters& conf):
		art::EDAnalyzer(conf),
		_diagLevel(conf().diagLevel()),
		_mcdiag(conf().mcdiag()),
		_calocrysTag(conf().calocrysTag()),
		_caloclusterTag(conf().caloclusterTag()),
		_genTag(conf().genTag()),
		_kalrepTag(conf().kalrepTag()),		
		_tcmatchTag(conf().tcmatchTag()){}
	

  void IPAAna::beginJob(){
	//std:cout<<"[In BeginJob()] Beginning ..."<<std::endl;
	art::ServiceHandle<art::TFileService> tfs;
	_Ntup  = tfs->make<TTree>("CaloCalibAna", "CaloCalibAna");
	_Ntup->Branch("evt",          	&_evt ,        "evt/I");
	_Ntup->Branch("run",          	&_run ,        "run/I");

	_Ntup->Branch("nGen",         	&_nGen ,        "nGen/I");
	_Ntup->Branch("genId",        	&_genPdgId,     "genId/I");
	_Ntup->Branch("genCrCode",    	&_genCrCode,    "genCrCode/I");
	_Ntup->Branch("genMomX",      	&_genmomX,      "genMomX/F");
	_Ntup->Branch("genMomY",      	&_genmomY,      "genMomY/F");
	_Ntup->Branch("genMomZ",      	&_genmomZ,      "genMomZ/F");
	_Ntup->Branch("genStartX",    	&_genStartX,    "genStartX/F");
	_Ntup->Branch("genStartY",    	&_genStartY,    "genStartY/F");
	_Ntup->Branch("genStartZ",    	&_genStartZ,    "genStartZ/F");
	_Ntup->Branch("genStartT",    	&_genStartT,    "genStartT/F");

	_Ntup->Branch("nCry",         	&_nHits ,       "nCry/I");
	_Ntup->Branch("cryId",        	&_cryId ,       "cryId[nCry]/I");
	_Ntup->Branch("crySectionId", 	&_crySectionId, "crySectionId[nCry]/I");
	_Ntup->Branch("cryPosX",      	&_cryPosX ,     "cryPosX[nCry]/F");
	_Ntup->Branch("cryPosY",      	&_cryPosY ,     "cryPosY[nCry]/F");
	_Ntup->Branch("cryPosZ",      	&_cryPosZ ,     "cryPosZ[nCry]/F");
	_Ntup->Branch("cryEdep",      	&_cryEdep ,     "cryEdep[nCry]/F");
	_Ntup->Branch("cryTime",      	&_cryTime ,     "cryTime[nCry]/F");
	_Ntup->Branch("cryDose",      	&_cryDose ,     "cryDose[nCry]/F");
	_Ntup->Branch("cryRadius",	&_cryRadius,	"cryRadius[nCry]/F");
	
	_Ntup->Branch("nClu",         	&_nClusters ,       	"nClu/I");
	_Ntup->Branch("clustertime", 	&_clustertime,		"clustertime/F");
	_Ntup->Branch("clustertimeErr", &_clustertimeErr,	"clustertimeErr/F");
	_Ntup->Branch("clusterEdep", 	&_clusterEdep,		"clusterEdep/F");
	_Ntup->Branch("clusterEDepErr", &_clusterEDepErr, 	"clusterEdepErr/F");
	_Ntup->Branch("clusterangle", 	&_clusterangle, 	"clusterangle/F");
	_Ntup->Branch("clusterPosX",	&_clustercog3VectorX, 	"clustercog3VectorX/F");
	_Ntup->Branch("clusterPosY",	&_clustercog3VectorY, 	"clustercog3VectorY/F");
	_Ntup->Branch("clusterPosZ",	&_clustercog3VectorZ, 	"clustercog3VectorZ/F");
	_Ntup->Branch("clusterR", 	&_clusterR,		"clusterR/F");
	_Ntup->Branch("clusterNHits", 	&_clusterNHits,		"clusterNHits/F");

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

	_Ntup->Branch("nMatches",	&_nMatches,		"nMatches/I");	
	_Ntup->Branch("matchPosXtrk",	&_matchPosXtrk,		"matchPosXtrk/F");
	_Ntup->Branch("matchPosYtrk",	&_matchPosYtrk,		"matchPosYtrk/F");
	_Ntup->Branch("matchPosZtrk",	&_matchPosZtrk,		"matchPosZtrk/F");
	_Ntup->Branch("matchTtrk",	&_matchTtrk,		"matchTtrk/F");
	_Ntup->Branch("matchChi2", 	&_matchChi2, 		"matchChi/F");
	_Ntup->Branch("matchEDep", 	&_matchEDep, 		"matchEDep/F");
	_Ntup->Branch("matchPosXcl", 	&_matchPosXCl, 		"matchPosXCl/F");
	_Ntup->Branch("matchPosYCl", 	&_matchPosYCl, 		"matchPosYCl/F");
	_Ntup->Branch("matchPosZCl", 	&_matchPosZCl, 		"matchPosZClF"); 
	_Ntup->Branch("matchPathLen", 	&_matchPathLen,	        "matchPathLen/F");
	_Ntup->Branch("matchR",		&_matchR, 		"matchR/F");
	_Ntup->Branch("matchDt",	&_matchDt,		"matchDt/F");

	_Ntup->Branch("trackMom",	&_TrackMom, 		"trackMom/F");
	_Ntup->Branch("EoP",		&_TrackEoP,			"EoP/F");
	outputfile.open("IPAAna.csv");
        outputfile<<event,run,cryId,cryE,trackerE<<std::endl;s
  }


  void IPAAna::analyze(const art::Event& event) {
	//std:cout<<"[In Analyze()] Beginning ..."<<std::endl;
	_evt = event.id().event();
	_run = event.run();

	if(!findData(event)) 
		throw cet::exception("RECO")<<"No data in  event"<< endl; 

	//std:cout<<"[In Analyze()] Found Data ..."<<std::endl;

      	art::ServiceHandle<GeometryService> geom;
      	if( ! geom->hasElement<Calorimeter>() ) return;
      	Calorimeter const & cal = *(GeomHandle<Calorimeter>());
	
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
			_matchEDep = cl->energyDep();
			_matchPosXCl =x1.x();
			_matchPosYCl =x1.y();
		        _matchPosZCl = x1.z();
		        _matchPathLen = tcm.ds();
			_matchDt = tcm.dt();
			_matchR = sqrt(x1.x()*x1.x() + x1.y()*x1.y());
		       _nMatches++;
			}
		}
      
	     double Ep=ClosestCluster->energyDep();
	     if(EndMom!=0) _TrackEoP = Ep/EndMom;

//================== GenPartile Info======================//
	//std:cout<<"[In Analyze()] Getting GenInfo ..."<<std::endl;
	_nGen = _gencol->size();
	for (unsigned int i=0; i <_gencol->size(); ++i)
	{
		GenParticle const& gen = (*_gencol)[i];
		_genPdgId   = gen.pdgId();
		_genCrCode  = gen.generatorId().id();
		_genmomX    = gen.momentum().vect().x();
		_genmomY    = gen.momentum().vect().y();
		_genmomZ    = gen.momentum().vect().z();
		_genStartX  = gen.position().x()+ 3904;
		_genStartY  = gen.position().y();
		_genStartZ  = gen.position().z();
		_genStartT  = gen.time();
	} 


//=======================Get Cluster Info ==============//
	//std:cout<<"[In Analyze()] Getting Cluster Info..."<<std::endl;
        for (unsigned int tclu=0; tclu<_calclustercol->size();++tclu){
		CaloCluster const& cluster = (*_calclustercol)[tclu];
		const CaloCluster::CaloCrystalHitPtrVector caloClusterHits = 		
						 cluster.caloCrystalHitsPtrVector();

		CLHEP::Hep3Vector crystalPos   	= cal.geomUtil().mu2eToDiskFF(cluster.diskId(), 	 						cluster.cog3Vector()); 
		_clusterNHits 			=  caloClusterHits.size();
		_clusterdiskId  		= cluster.diskId();
		_clustertime      		= cluster.time();
		_clusterEdep      		= cluster.energyDep(); 
		_clusterEDepErr     		= cluster.energyDepErr();
		_clusterangle      		= cluster.angle();
		_clustercog3VectorX      	= cluster.cog3Vector().x();
		_clustercog3VectorY      	= cluster.cog3Vector().y();
		_clustercog3VectorZ      	= cluster.cog3Vector().z();
		_clusterR  			= sqrt(cluster.cog3Vector().x()* cluster.cog3Vector().x()
								+ cluster.cog3Vector().y()*cluster.cog3Vector().y());
		outputfile<<_evt<<","<<_run<<","<<_calcryhitcol->size()<<std::endl;
        	_nClusters++;
	}
		
//=====================Crystal Hits Info =======================//
	//std:cout<<"[In Analyze()] Getting Crystal Info..."<<std::endl;
	_nHits = _calcryhitcol->size();
	for (unsigned int ic=0; ic<_calcryhitcol->size();++ic) 
	{	   
		   CaloCrystalHit const& hit      = (*_calcryhitcol)[ic];
		   int diskId                     = cal.crystal(hit.id()).diskId();
		   CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position());  
		   _cryId[ic] 	 	= hit.id();
		   _cryTime[ic]      	= hit.time();
		   _cryEdep[ic]      	= hit.energyDep();
		   _cryTotE[ic] 	= hit.energyDepTot();
		   _cryTotEErr[ic] 	= hit.energyDepTotErr();
		   _cryPosX[ic]     	= crystalPos.x();
		   _cryPosY[ic]      	= crystalPos.y();
		   _cryPosZ[ic]      	= crystalPos.z();
		   _cryRadius[ic] 	= sqrt(crystalPos.x()*crystalPos.x() +  
						    crystalPos.y()*crystalPos.y());
	}
	_Ntup->Fill();
	_nEvents++;
}


bool IPAAna::findData(const art::Event& evt){

	_calcryhitcol =0;
	_calclustercol=0;
	_gencol=0;
	
	auto genpart = evt.getValidHandle<GenParticleCollection>(_genTag);
	_gencol = genpart.product();
	auto cryhit = evt.getValidHandle<CaloCrystalHitCollection>(_calocrysTag);
	_calcryhitcol =cryhit.product();
	auto cluster= evt.getValidHandle<CaloClusterCollection>(_caloclusterTag);
	_calclustercol =cluster.product();
	
        
	return  _gencol!=0 and _calcryhitcol!=0 && _calclustercol !=0;
       }

 void IPAAna::endJob(){} 

}

DEFINE_ART_MODULE(mu2e::IPAAna);*/
