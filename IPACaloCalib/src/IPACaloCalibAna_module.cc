
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

  class IPACaloCalibAna : public art::EDAnalyzer {
     
     public:

     struct Config {
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
 	      fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
	     fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
	     //fhicl::Atom<art::InputTag> calohitsTag{Name("CaloHitCollection"), Comment("readout")};
	     fhicl::Atom<art::InputTag> calocrysTag{Name("CaloCrystalHitCollection"),Comment("cal reco crystal hit info")};
	     fhicl::Atom<art::InputTag> kalrepTag{Name("KalRepPtrCollection"),Comment("outcome of Kalman filter (for tracker momentum info)")};
	     //fhicl::Atom<art::InputTag> calohitMCcrysTag{Name("CaloHitMCTruthCollection"),Comment("cal hit mc truth info")};
	     fhicl::Atom<art::InputTag> caloclusterTag{Name("CaloClusterCollection"),Comment("cal reco cluster info")};
	     fhicl::Atom<art::InputTag> tcmatchTag{Name("TrackClusterMatchCollection"), Comment("track calo match"), "TrackCaloMatching"};
	    fhicl::Atom<art::InputTag> genTag{Name("GenParticleCollection"), Comment("gen particle info")};
	     //fhicl::Atom<art::InputTag> caloSimPartMCTag{Name("CaloHitSimPartMCCollection"), Comment("cal hit sim particle mc info")};
    };
       typedef art::EDAnalyzer::Table<Config> Parameters;

       explicit IPACaloCalibAna(const Parameters& conf);
       virtual ~IPACaloCalibAna() {}
       
      
       virtual void beginJob();
       virtual void endJob();
	virtual void analyze(const art::Event& e) override;

     private:
       Config _conf;
       int _diagLevel;
       int _nProcess;
       int _mcdiag;
       
       art::InputTag _calohitsTag;
       art::InputTag _calocrysTag;
       art::InputTag _kalrepTag;
       art::InputTag _calohitMCcrysTag;
       art::InputTag _caloclusterTag;
       art::InputTag _caloSimPartMCTag;
       art::InputTag _tcmatchTag;
       art::InputTag _genTag;
       
       const CaloHitCollection* _calhitcol;
       const CaloCrystalHitCollection*  _calcryhitcol;
       const KalRepPtrCollection* _kalrepcol;
       const CaloHitMCTruthCollection* _calhitMCtruecol;
       const CaloClusterCollection* _calclustercol;
       const CaloHitSimPartMCCollection* _simpartcol;
       const TrackClusterMatchCollection* _tcmatchcol;
       const GenParticleCollection *_gencol;
       
       //TCanvas* canvas_ = nullptr;
       TH2F* _hviewxy = nullptr;
       TH2F* _hviewxz = nullptr;
       TH1F* _hfit = nullptr;
       TTree* _Ntup;
       std::vector<double> calibration_constants;
       std::vector<int> crystallist;

       double MostProbEoP;
       
       int   _nEvents =0;
       int _evt, _run, _nHits, _nClusters, _nTracks, _nSim, _nMatches, _nTrackMatched;
       int   _nGen,_genPdgId[163840],_genCrCode[163840];

       float _genmomX[163840],_genmomY[163840],_genmomZ[163840], _genStartX[163840],_genStartY[163840],_genStartZ[163840],_genStartT[163840];

       int   _cryId[1638400],_crySectionId[1638400],_crySimIdx[1638400],_crySimLen[1638400];
       float _cryTime[1638400],_cryEdep[1638400],_cryDose[1638400],_cryPosX[1638400],_cryPosY[1638400],_cryPosZ[1638400],_cryLeak[1638400], _cryTotE[1638400],_cryTotSum[1638400], _cryTotEErr[1638400], _cryRadius[1638400],_cryMaxEoP[1638400],_cryMaxR[1638400];

       int   _clusterId[163840],_clusterdiskId[163840];
       float _clustertime[163840], _clustertimeErr[163840], _clusterEdep[163840], _clusterEDepErr[163840], _clusterangle[163840], _clustercog3VectorX[163840], _clustercog3VectorY[163840], _clustercog3VectorZ[163840],_clusterR[163840],_clusterNHits[163840], _clustermaxECrystal[163840],	_clusterindexMaxECrystal[163840], _clusterERatio[163840];
    
       float _TrackT0[163840], _TrackT0Err[163840], _TrackMom[163840], _MaxEoP[163840], _EoP[163840], _TrackBackTime[163840] , _TrackBackOmega[163840] ,_TrackBackD0[163840] , _TrackBackZ0[163840], _TrackBackPhi0[163840], _TrackBackTanDip[163840], _TrackChi2[163840], _TrackChi2DOF[163840], _TrackCosTheta[163840];
      
	float _matchChi2[163840], _matchEDep[163840], _matchPosXCl[163840], _matchPosYCl[163840], _matchPosZCl[163840], _matchPathLen[163840], _matchR[163840], _matchDt[163840], _matchPosXtrk[163840], _matchPosYtrk[163840], _matchPosZtrk[163840],_matchTtrk[163840];
    
      int NCrystals = 674;

       bool findData(const art::Event& evt);
      
       double GetMostProbEoP(TH1F* h);
       void GetGenPartInfo(const art::Event& evt);
       void GetKalInfo(const art::Event& evt);
       
       void FillXY();
       void fitEoP(TH1F* h);
};


  IPACaloCalibAna::IPACaloCalibAna(const Parameters& conf):
    art::EDAnalyzer(conf),
    _diagLevel(conf().diagLevel()),
    _mcdiag(conf().mcdiag()),
    //_calohitsTag(conf().calohitsTag()),
    _calocrysTag(conf().calocrysTag()),
    _kalrepTag(conf().kalrepTag()),
    //_calohitMCcrysTag(conf().calohitMCcrysTag()),
    _caloclusterTag(conf().caloclusterTag()),
    //_caloSimPartMCTag(conf().caloSimPartMCTag())
    _tcmatchTag(conf().tcmatchTag()),
    _genTag(conf().genTag())
   //_pidTag(conf().pidTag())
{}

  void IPACaloCalibAna::beginJob(){
    for(int i=0;i<NCrystals;i++){
		crystallist.push_back(0);
	}
    art::ServiceHandle<art::TFileService> tfs;
    _Ntup  = tfs->make<TTree>("CaloCalibAna", "CaloCalibAna");
    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");

    _Ntup->Branch("nGen",         &_nGen ,        "nGen/I");
    _Ntup->Branch("genId",        &_genPdgId,     "genId[nGen]/I");
    _Ntup->Branch("genCrCode",    &_genCrCode,    "genCrCode[nGen]/I");
    _Ntup->Branch("genMomX",      &_genmomX,      "genMomX[nGen]/F");
    _Ntup->Branch("genMomY",      &_genmomY,      "genMomY[nGen]/F");
    _Ntup->Branch("genMomZ",      &_genmomZ,      "genMomZ[nGen]/F");
    _Ntup->Branch("genStartX",    &_genStartX,    "genStartX[nGen]/F");
    _Ntup->Branch("genStartY",    &_genStartY,    "genStartY[nGen]/F");
    _Ntup->Branch("genStartZ",    &_genStartZ,    "genStartZ[nGen]/F");
    _Ntup->Branch("genStartT",    &_genStartT,    "genStartT[nGen]/F");

    _Ntup->Branch("nCry",         &_nHits ,       "nCry/I");
    _Ntup->Branch("nEvents",         &_nEvents ,       "nEvents/I");
    _Ntup->Branch("cryId",        &_cryId ,       "cryId[nCry]/I");
    _Ntup->Branch("crySectionId", &_crySectionId, "crySectionId[nCry]/I");
    _Ntup->Branch("cryPosX",      &_cryPosX ,     "cryPosX[nCry]/F");
    _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosY[nCry]/F");
    _Ntup->Branch("cryPosZ",      &_cryPosZ ,     "cryPosZ[nCry]/F");
    _Ntup->Branch("cryEdep",      &_cryEdep ,     "cryEdep[nCry]/F");
    _Ntup->Branch("cryTime",      &_cryTime ,     "cryTime[nCry]/F");
    _Ntup->Branch("cryDose",      &_cryDose ,     "cryDose[nCry]/F");
    _Ntup->Branch("cryTotE",      &_cryTotE ,     "cryTotE[nCry]/F");
    _Ntup->Branch("cryTotSum",    &_cryTotSum,    "cryTotSum/F");
    _Ntup->Branch("cryTotEErr",   &_cryTotEErr,	  "cryTotEErr[nCry]/F");   
    _Ntup->Branch("cryRadius",    &_cryRadius,	  "cryRadius[nCry]/F");
    _Ntup->Branch("cryMaxEoP",    &_cryMaxEoP,    "cryMaxEoP[nEvents]/F");
    _Ntup->Branch("cryMaxR",	  &_cryMaxR,      "cryMaxR[nEvents]/F");

    _Ntup->Branch("nClu",         	&_nClusters ,       	"nClu/I");
    _Ntup->Branch("clustertime", 	&_clustertime,		"clustertime[nClu]/F");
    _Ntup->Branch("clustertimeErr", 	&_clustertimeErr,	"clustertimeErr[nClu]/F");
    _Ntup->Branch("clusterEdep", 	&_clusterEdep,		"clusterEdep[nClu]/F");
    _Ntup->Branch("clusterEDepErr", 	&_clusterEDepErr, 	"clusterEdepErr[nClu]/F");
    _Ntup->Branch("clusterangle", 	&_clusterangle, 	"clusterangle[nClu]/F");
    _Ntup->Branch("clusterPosX",	&_clustercog3VectorX, 	"clustercog3VectorX[nClu]/F");
    _Ntup->Branch("clusterPosY",	&_clustercog3VectorY, 	"clustercog3VectorY[nClu]/F");
    _Ntup->Branch("clusterPosZ",	&_clustercog3VectorZ, 	"clustercog3VectorZ[nClu]/F");
    _Ntup->Branch("clusterR", 		&_clusterR,		"clusterR[nClu]/F");
    _Ntup->Branch("clusterNHits", 		&_clusterNHits,		"clusterNHits[nClu]/F");
    _Ntup->Branch("clustermaxECrystsl", 		&_clustermaxECrystal,		"clustermaxECrystal[nClu]/F");
    _Ntup->Branch("clusterindexMaxECrystal", 		&_clusterindexMaxECrystal, "clusterindexMaxECrystal[nClu]/F");
	_Ntup->Branch("clusterERatio",    &_clusterERatio, "clusterERatio[nClu]/F");
    
    _Ntup->Branch("nTracks",		&_nTracks,		"nTracks/I");
    _Ntup->Branch("TrackT0", 		&_TrackT0, 		"TrackT0[nTracks]/F");
    _Ntup->Branch("TrackT0Err", 	&_TrackT0Err,		"TrackT0Err[nTracks]/F");
    _Ntup->Branch("TrackBackTime", 	&_TrackBackTime,   	"TrackBackTime[nTracks]/F"); 
    _Ntup->Branch("TrackBackOmega", 	&_TrackBackOmega,	"TrackBackOmega[nTracks]/F");
    _Ntup->Branch("TrackBackD0",  	&_TrackBackD0, 		"TrackBackD0[nTracks]/F"); 	
    _Ntup->Branch("TrackBackZ0", 	&_TrackBackZ0, 		"TrackBackZ0[nTracks]/F");
    _Ntup->Branch("TrackBackPhi0",	&_TrackBackPhi0,	"TrackBackPhi0[nTracks]/F");
    _Ntup->Branch("TrackBackTanDip",	&_TrackBackTanDip,	"TrackBackTanDip[nTracks]/F");
    _Ntup->Branch("TrackChi2"	,	&_TrackChi2,		"TrackChi2[nTracks]/F");
    _Ntup->Branch("TrackChi2DOF", 	&_TrackChi2DOF, 	"TrackCho2DOF[nTracks]/F");
    _Ntup->Branch("TrackCosTheta",       &_TrackCosTheta,	"TrackCosTheta[nTracks]/F");

    _Ntup->Branch("nMatches",		&_nMatches,		"nMatches/I");	
    _Ntup->Branch("matchPosXtrk",	&_matchPosXtrk,		"matchPosXtrk[nMatches]/F");
    _Ntup->Branch("matchPosYtrk",	&_matchPosYtrk,		"matchPosYtrk[nMatches]/F");
    _Ntup->Branch("matchPosZtrk",	&_matchPosZtrk,		"matchPosZtrk[nMatches]/F");
    _Ntup->Branch("matchTtrk",	        &_matchTtrk,		"matchTtrk[nMatches]/F");
    _Ntup->Branch("matchChi2", 		&_matchChi2, 		"matchChi[nMatches]/F");
    _Ntup->Branch("matchEDep", 		&_matchEDep, 		"matchEDep[nMatches]/F");
    _Ntup->Branch("matchPosXcl", 	&_matchPosXCl, 		"matchPosXCl[nMatches]/F");
    _Ntup->Branch("matchPosYCl", 	&_matchPosYCl, 		"matchPosYCl[nMatches]/F");
    _Ntup->Branch("matchPosZCl", 	&_matchPosZCl, 		"matchPosZCl[nMatches]F"); 
    _Ntup->Branch("matchPathLen", 	&_matchPathLen,	        "matchPathLen[nMatches]/F");
    _Ntup->Branch("matchR",		&_matchR, 		"matchR[nMatches]/F");
    _Ntup->Branch("matchDt",		&_matchDt,		"matchDt[nMatches]/F");
   
    _Ntup->Branch("trackMom",		&_TrackMom, 		"trackMom[nTracks]/F");
    _Ntup->Branch("EoP",		 &_EoP,			"EoP[nTracks]/F");
    _Ntup->Branch("MaxEoP",		&_MaxEoP, 		"MaxEoP[nTracks]/F");
  
    _hviewxy = tfs->make<TH2F>("hxy", "hxy",  350,-700,700,350,-700,700  );
   
    _hfit = tfs->make<TH1F>("EoP","EoP", 50, 0, 1);
    //canvas_ = tfs->make<TCanvas>("canvas XY", "XY" ,1300,800);
  }


  void IPACaloCalibAna::analyze(const art::Event& event) {
	
       _evt = event.id().event();
       _run = event.run();
    
       if(!findData(event)) 
      		throw cet::exception("RECO")<<"No data in  event"<< endl; 
       
      ++_nProcess;
      if (_nProcess%1000==0) std::cout<<"Processing event "<<_nProcess<<std::endl;
      
      GetKalInfo(event);
      GetGenPartInfo(event);

      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());
  //Utility to match  cloHits with MCtruth, simParticles and StepPoints
      //CaloHitMCNavigator caloHitNavigator(caloHits, caloHitsMCTruth, caloHitSimPartMC);
   
	//TH1F *histEoP = new TH1F("histEoP", "histEoP",100,0,1);//local hist for max bins
       
        for (unsigned int tclu=0; tclu<_calclustercol->size();++tclu){
	   CaloCluster const& cluster = (*_calclustercol)[tclu];
	   _clusterdiskId[_nClusters]  = cluster.diskId();
           CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(cluster.diskId(),cluster.cog3Vector()); 
	   _clustertime[_nClusters]      = cluster.time();
	   _clusterEdep[_nClusters]      = cluster.energyDep(); 
	   _clusterEDepErr[_nClusters]      = cluster.energyDepErr();
           _clusterangle[_nClusters]      = cluster.angle();
	   _clustercog3VectorX[_nClusters]      = cluster.cog3Vector().x();
	   _clustercog3VectorY[_nClusters]      = cluster.cog3Vector().y();
	   _clustercog3VectorZ[_nClusters]      = cluster.cog3Vector().z();
	   _clusterR[_nClusters]  		=sqrt(cluster.cog3Vector().x()*cluster.cog3Vector().x()+cluster.cog3Vector().y()*cluster.cog3Vector().y());
        const CaloCluster::CaloCrystalHitPtrVector caloClusterHits = cluster.caloCrystalHitsPtrVector();
       _clusterNHits[_nClusters] =  caloClusterHits.size();
       double    maxECrystal(0);
       for (unsigned int j=0; j<caloClusterHits.size(); ++j){
	       const CaloCrystalHit* crystalHit = &(*caloClusterHits.at(j));
	       double   crystalEnergy = crystalHit->energyDep();
	       if (crystalEnergy > maxECrystal) {
	           _clustermaxECrystal[_nClusters]        = crystalEnergy;
	           _clusterindexMaxECrystal[_nClusters]   = j;
		   maxECrystal = crystalEnergy;
	       }
        }
	if (maxECrystal !=0){
        	_clusterERatio[_nClusters]   = cluster.energyDep()/maxECrystal;
	}
        _nClusters++;
	}

	
        _cryTotSum[_evt] =0;
	double maxCrysEoP = 0;
	double maxCrysR = 0;
	_nHits = _calcryhitcol->size();
	for (unsigned int ic=0; ic<_calcryhitcol->size();++ic) 
       {	   
	   CaloCrystalHit const& hit      = (*_calcryhitcol)[ic];
	   int diskId                     = cal.crystal(hit.id()).diskId();
           CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position());  //in disk FF frame
           //CaloHit const& caloHit         = *(hit.readouts().at(0));
	   //CaloHitSimPartMC const& hitSim = caloHitNavigator.sim(caloHit);
           //int nPartInside                = hitSim.simParticles().size();
	   _cryId[ic] 	 =     hit.id();
	   _cryTime[ic]      = hit.time();
	   _cryEdep[ic]      = hit.energyDep();
	   _cryTotSum[_nEvents] 	+= hit.energyDep();
	   _cryTotE[ic] 	= hit.energyDepTot();
	   _cryTotEErr[ic] =hit.energyDepTotErr();
	   _cryPosX[ic]      = crystalPos.x();
	   _cryPosY[ic]      = crystalPos.y();
	   _cryPosZ[ic]      = crystalPos.z();
	   _hviewxy->Fill(crystalPos.x(),crystalPos.y(),hit.energyDep());
	   //histEoP->Fill(hit.energyDep());
	   _cryRadius[ic] 	  = sqrt(crystalPos.x()*crystalPos.x() + crystalPos.y()*crystalPos.y());
	   //int maxBin = histEoP->GetMaximumBin();
	   //double maxEoP = histEoP->GetBinCenter(maxBin);
	   //if (maxEoP > maxCrysEoP){
		//maxCrysEoP=maxEoP;
		//maxCrysR = sqrt(crystalPos.x()*crystalPos.x() + crystalPos.y()*crystalPos.y());
	   //}
	   crystallist[hit.id()]+=1;
	   
	}
	if(maxCrysEoP!=0 and maxCrysR!=0){
		_cryMaxEoP[_nEvents] = maxCrysEoP;
		_cryMaxR[_nEvents] = maxCrysR;
	}
       //_MaxEoP[_evt] = GetMostProbEoP(histEoP);
	
       _Ntup->Fill();
       _nEvents++;
}


void IPACaloCalibAna::GetGenPartInfo(const art::Event& evt){
	
	_nGen = _gencol->size();
       for (unsigned int i=0; i <_gencol->size(); ++i)
       {
           GenParticle const& gen = (*_gencol)[i];
	   _genPdgId[i]   = gen.pdgId();
	   _genCrCode[i]  = gen.generatorId().id();
	   _genmomX[i]    = gen.momentum().vect().x();
	   _genmomY[i]    = gen.momentum().vect().y();
	   _genmomZ[i]    = gen.momentum().vect().z();
	   _genStartX[i]  = gen.position().x()+ 3904;
	   _genStartY[i]  = gen.position().y();
	   _genStartZ[i]  = gen.position().z();
	   _genStartT[i]  = gen.time();
       } 
}

void IPACaloCalibAna::GetKalInfo(const art::Event& evt){
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

		_TrackT0[_nTracks] = TrackKrep->t0().t0();
    		_TrackT0Err[_nTracks] = TrackKrep->t0().t0Err();
		_TrackMom[_nTracks] = TrackKrep->momentum(sz1).mag();
		_TrackBackTime[_nTracks] =   TrackKrep->arrivalTime(sz1);
		 HelixParams helx  = TrackKrep->helix(sz1);
    		_TrackBackOmega[_nTracks]       = helx.omega(); 
    		_TrackBackD0[_nTracks]       = helx.d0();
    		_TrackBackZ0[_nTracks]       = helx.z0();
    		_TrackBackPhi0[_nTracks]     = helx.phi0();
    		_TrackBackTanDip[_nTracks]   = helx.tanDip(); 

		double entlen         = std::min(s1,s2);
		CLHEP::Hep3Vector fitmom = TrackKrep->momentum(entlen);
		TLorentzVector  Momentum(fitmom.x(),fitmom.y(),fitmom.z(),0.511);
		_TrackCosTheta[_nTracks] = Momentum.CosTheta();
    		_TrackChi2[_nTracks] =TrackKrep->chisq();
		_TrackChi2DOF[_nTracks]= TrackKrep->chisq()/TrackKrep->nActive();
		_nTracks ++;
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
			_matchPosXtrk[_nMatches] = tcm.xtrk();
			_matchPosYtrk[_nMatches] = tcm.ytrk();
			_matchPosZtrk[_nMatches] = tcm.ztrk();
			_matchTtrk[_nMatches]	 = tcm.ttrk();
			_matchChi2[_nMatches] = tcm.chi2();
			_matchEDep[_nMatches] = cl->energyDep();
			_matchPosXCl[_nMatches] =x1.x();
			_matchPosYCl[_nMatches] =x1.y();
		        _matchPosZCl[_nMatches] = x1.z();
		        _matchPathLen[_nMatches] = tcm.ds();
			_matchDt[_nMatches] = tcm.dt();
			_matchR[_nMatches] = sqrt(x1.x()*x1.x() + x1.y()*x1.y());
		       _nMatches++;
			}
		}
      
    
     double Ep=ClosestCluster->energyDep();
    _EoP[_nTracks] = Ep/EndMom;
    _hfit->Fill(Ep/EndMom);
    _nTrackMatched++;
	}

}

bool IPACaloCalibAna::findData(const art::Event& evt){

	_calhitcol =0;
	_calcryhitcol =0;
	_calclustercol=0;
	_tcmatchcol=0;
	_gencol=0;
	
	auto genpart = evt.getValidHandle<GenParticleCollection>(_genTag);
	_gencol = genpart.product();
	auto cryhit = evt.getValidHandle<CaloCrystalHitCollection>(_calocrysTag);
	_calcryhitcol =cryhit.product();
	auto cluster= evt.getValidHandle<CaloClusterCollection>(_caloclusterTag);
	_calclustercol =cluster.product();
        auto kalrep = evt.getValidHandle<KalRepPtrCollection>(_kalrepTag);
	_kalrepcol =kalrep.product();
	auto tcmatch = evt.getValidHandle<TrackClusterMatchCollection>(_tcmatchTag);
	_tcmatchcol = tcmatch.product();
        if(_mcdiag){
	    
	   _calhitMCtruecol=0;
	   _simpartcol = 0;
           auto simpar= evt.getValidHandle<CaloHitSimPartMCCollection>(_caloSimPartMCTag);
	   _simpartcol =simpar.product();
	   auto calhitMC = evt.getValidHandle<CaloHitMCTruthCollection>(_calohitMCcrysTag);
	   _calhitMCtruecol =calhitMC.product();
          
        }
	return  _gencol!=0 && _kalrepcol!=0 && _calcryhitcol!=0 && _calclustercol !=0 && (_simpartcol != 0 || !_mcdiag) &&_tcmatchcol!=0;
       }

 void IPACaloCalibAna::endJob(){
	//MostProbEoP= GetMostProbEoP(_hfit);
	//fitEoP(_hfit);
	//FillXY();
	
  } 

/*
void IPACaloCalibAna::FillXY(){
	TBox box;
	TArc arc;
	auto xyplot = canvas_->DrawFrame(-1000,-1000, 1000,1000);
	xyplot->GetYaxis()->SetTitleOffset(1.25);
	xyplot->SetTitle( "XY; X(mm);Y(mm)");
	art::ServiceHandle<GeometryService> geom;
        if( ! geom->hasElement<Calorimeter>() ) return;
        Calorimeter const & cal = *(GeomHandle<Calorimeter>());
	canvas_->SetTitle("foo title");
	Disk const & disk =  cal.disk(1);
	double disk_outR = disk.outerRadius();
	double disk_inR= disk.innerRadius();
	arc.SetFillStyle(0);
      	arc.DrawArc(0.,0., disk_inR);
      	arc.DrawArc(0.,0., disk_outR);
	for(int i=0;i<NCrystals;i++){
	        Crystal const &crystal = cal.crystal(i);
	   	double crystalXLen = crystal.size().x();
		double crystalYLen = crystal.size().y();
		CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(1,crystal.position());
	  	if(crystallist[i] ==0) box.SetFillColor(0);
		if(crystallist[i]>0 and crystallist[i]<11) box.SetFillColor(kViolet-10+crystallist[i]);
		if(crystallist[i]>10 and crystallist[i]<21)box.SetFillColor(kBlue-20+crystallist[i]);
		if(crystallist[i]>20 and crystallist[i]<31)box.SetFillColor(kAzure-30+crystallist[i]);
		if(crystallist[i]>31 and crystallist[i]<41) box.SetFillColor(kAzure-30+crystallist[i]);
		if(crystallist[i]>40 and crystallist[i]<51)box.SetFillColor(kTeal-50+crystallist[i]);
		if(crystallist[i]>50 and crystallist[i]<61)box.SetFillColor(kTeal-50+crystallist[i]);
		if(crystallist[i]>60 and crystallist[i]<71) box.SetFillColor(kGreen-70+crystallist[i]);
		if(crystallist[i]>70 and crystallist[i]<81)box.SetFillColor(kGreen-70+crystallist[i]);
		if(crystallist[i]>80 and crystallist[i]<91)box.SetFillColor(kYellow-90+crystallist[i]);
		if(crystallist[i]>90 and crystallist[i]<101)box.SetFillColor(kYellow-90+crystallist[i]);
		if(crystallist[i]>100 and crystallist[i]<111)box.SetFillColor(kOrange-110+crystallist[i]);
		if(crystallist[i]>110 and crystallist[i]<121)box.SetFillColor(kOrange-110+crystallist[i]);
		if(crystallist[i]>100) box.SetFillColor(kRed+4);
		box.DrawBox(crystalPos.x()-crystalXLen/2, crystalPos.y()-crystalYLen/2,crystalPos.x()+crystalXLen/2, crystalPos.y()+crystalYLen/2);
		if(crystallist[i] !=0){
			TLatex latex;
			stringstream crys;
                	crys<<crystallist[i];
                	const char* str_crys = crys.str().c_str();
		   	latex.SetTextSize(0.02);
		   	latex.DrawLatex(crystalPos.x()-crystalXLen/2, crystalPos.y()-crystalYLen/2,str_crys);
		   }
	}
	canvas_->Update();
	canvas_->SaveAs("AllCrystalsInAllEvents.root");
	
}

*/
 
double IPACaloCalibAna::GetMostProbEoP(TH1F* h){
	int maxBin = h->GetMaximumBin();
	double maxEdep = h->GetBinCenter(maxBin);
	return maxEdep;
}

void IPACaloCalibAna::fitEoP(TH1F* h){
 	TF1 *_fitF = new TF1("fitF","[4]*ROOT::Math::crystalball_function(x, [3],[2],[1],[0])",0.6,1);
	//h->Scale(1/h->Integral());
	int bin = h->GetMaximumBin();
	int content = h->GetBinContent(bin);
        _fitF->SetParameter(4,content);//N
	_fitF->SetParameter(0,0.8);//mean
	_fitF->SetParameter(1,0.12);//sigma
	_fitF->SetParameter(2,2);//n
	_fitF->SetParameter(3,1);//alpha
	h->Fit(_fitF);//(alpha, n sigma, mu)
	h->SaveAs("fit.root");

}

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::IPACaloCalibAna);
