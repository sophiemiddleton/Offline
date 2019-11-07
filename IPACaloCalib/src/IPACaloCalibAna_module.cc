
#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

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

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TDirectory.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2.h"

// Mu2e includes.
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"

// BaBar Kalman filter includes
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"


#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>

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
     fhicl::Atom<std::string> g4ModuleLabel{Name("g4ModuleLabel"), Comment("g4 label"),"g4run"};
     fhicl::Atom<std::string> generatorModuleLabel{Name("generatorModuleLabel"),Comment("gen label"),"generate"};
     fhicl::Atom<std::string> caloReadoutModuleLabel{Name("caloReadoutModuleLabel"), Comment("readout"),"CaloReadoutHitsMaker"};
     fhicl::Atom<std::string> caloCrystalModuleLabel{Name("caloCrystalModuleLabel"),Comment("cal label"),"CaloCrystalHitsMaker"};
	fhicl::Atom<std::string> caloHitMCCrystalPtrLabel{Name("caloHitMCCrystalPtrLabel"),Comment("cal label"),"CaloHitMCCrystalPtr"};
    fhicl::Atom<std::string> caloClusterPtrLabel{Name("caloClusterModuleLabel"),Comment("caloClusterModuleLabel"),"CaloClusterMakerNew"};
    fhicl::Atom<std::string> virtualDetectorLabel{Name("virtualDetectorLabel"),Comment("virtualDetectorLabel"),"virtualdetector"};
    fhicl::Atom<std::string> trkPatRecModuleLabel{Name("trkPatRecModuleLabel"),Comment("trkPatRecModuleLabel"),"TPRDownstreameMinus"};
    };
       typedef art::Ptr<StepPointMC> StepPtr;
       typedef std::vector<StepPtr>  StepPtrs;
       typedef std::map<int,StepPtrs > HitMap;
       typedef art::EDAnalyzer::Table<Config> Parameters;

       explicit IPACaloCalibAna(const Parameters& conf);
       virtual ~IPACaloCalibAna() {}

       virtual void beginJob();
       virtual void endJob();

       virtual void analyze(const art::Event& e) override;

     private:
       Config _conf;
       typedef std::vector<art::Handle<StepPointMCCollection>> HandleVector;
       typedef art::Ptr<CaloCrystalHit> CaloCrystalHitPtr;
       typedef art::Ptr<SimParticle> SimParticlePtr;


       int _diagLevel;
       int _nProcess;

       std::string _g4ModuleLabel;
       std::string _generatorModuleLabel;
       std::string _caloReadoutModuleLabel;
       std::string _caloCrystalModuleLabel;
       std::string _caloHitMCCrystalPtrLabel;
       std::string _caloClusterModuleLabel;
       std::string _caloClusterAlgorithm;
       std::string _caloClusterSeeding;
       const std::string _producerlabel;
       std::string _virtualDetectorLabel;
       std::string _trkPatRecModuleLabel;
       std::string _instancelabel;

       TrkParticle _tpart;
       TrkFitDirection _fdir;
       
       TH2F* _hviewxy;
       TH2F* _hviewxz;
       
       TTree* _Ntup;

       int   _evt,_run;
       int   _nHits;
	int   _nSim;
       int   _nGen,_genPdgId[16384],_genCrCode[16384];
	float _genmomX[16384],_genmomY[16384],_genmomZ[16384],_genStartX[16384],_genStartY[16384],_genStartZ[16384],_genStartT[16384];
       float _cryTime[16384],_cryEdep[16384],_cryDose[16384],_cryPosX[16384],_cryPosY[16384],_cryPosZ[16384],_cryLeak[16384];
       float _trueEdep[16384], _truetime[16384], trueCalX[16384] , trueCalY[16384] , trueCalZ[16384];
  };


  IPACaloCalibAna::IPACaloCalibAna(const Parameters& conf):
    art::EDAnalyzer(conf),
    _diagLevel(conf().diagLevel()),
    _g4ModuleLabel(conf().g4ModuleLabel()),
    _generatorModuleLabel(conf().generatorModuleLabel()),
    _caloReadoutModuleLabel(conf().caloReadoutModuleLabel()),
    _caloCrystalModuleLabel(conf().caloCrystalModuleLabel()),
    _caloHitMCCrystalPtrLabel(conf().caloHitMCCrystalPtrLabel()),
    _caloClusterModuleLabel(conf().caloClusterPtrLabel()),
    _virtualDetectorLabel(conf().virtualDetectorLabel()),
    _trkPatRecModuleLabel(conf().trkPatRecModuleLabel())
{}

  void IPACaloCalibAna::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _Ntup  = tfs->make<TTree>("CaloCalibAna", "CaloCalibAna");



    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");

    _Ntup->Branch("nGen",         &_nGen ,        "nGen/I");
    _Ntup->Branch("genId",        &_genPdgId,     "genId[nGen]/I");
    _Ntup->Branch("genCrCode",    &_genCrCode,    "genCrCode[nGen]/I");
    _Ntup->Branch("genMomX",      &_genmomX,      "genMomX[nGen]/F");
    _Ntup->Branch("genMomY",      &_genmomY,      "genMomX[nGen]/F");
    _Ntup->Branch("genMomZ",      &_genmomZ,      "genMomX[nGen]/F");
    _Ntup->Branch("genStartX",    &_genStartX,    "genStartX[nGen]/F");
    _Ntup->Branch("genStartY",    &_genStartY,    "genStartY[nGen]/F");
    _Ntup->Branch("genStartZ",    &_genStartZ,    "genStartZ[nGen]/F");
    _Ntup->Branch("genStartT",    &_genStartT,    "genStartT[nGen]/F");

  }

	

  void IPACaloCalibAna::endJob(){
  }




  void IPACaloCalibAna::analyze(const art::Event& event) {

      ++_nProcess;
      if (_nProcess%1000==0) std::cout<<"Processing event "<<_nProcess<<std::endl;
      
      
      //Get handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());
  
	/*
      //Get generated particles
      art::Handle<GenParticleCollection> gensHandle;
      event.getByLabel(_generatorModuleLabel, gensHandle);
      GenParticleCollection const& genParticles(*gensHandle);
*/
      //Get calorimeter readout hits (2 readout / crystal as of today)
      art::Handle<CaloHitCollection> caloHitsHandle;
      event.getByLabel(_caloReadoutModuleLabel, caloHitsHandle);
      CaloHitCollection const& caloHits(*caloHitsHandle);
      
      //Get calorimeter readout hits MC level - energy/time/type
      art::Handle<CaloHitMCTruthCollection> caloHitMCTruthHandle;
      event.getByLabel(_caloReadoutModuleLabel, caloHitMCTruthHandle);
      CaloHitMCTruthCollection const& caloHitsMCTruth(*caloHitMCTruthHandle);
      
      
      //Get simParticles and stepPointMC summary for crystal readout hits
      art::Handle<CaloHitSimPartMCCollection> caloHitSimMCHandle;
      event.getByLabel(_caloReadoutModuleLabel, caloHitSimMCHandle);
      CaloHitSimPartMCCollection const& caloHitSimPartMC(*caloHitSimMCHandle);
            
      //Get calo crystal hits (average from readouts)
      art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
      event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
      CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);

   
      //Utility to match  cloHits with MCtruth, simParticles and StepPoints
      CaloHitMCNavigator caloHitNavigator(caloHits, caloHitsMCTruth, caloHitSimPartMC);
   //--------------------------  Do generated particles --------------------------------

       _evt = event.id().event();
       _run = event.run();
              /*
       _nGen = genParticles.size();
	if(_nGen !=0){
	       for (unsigned int i=0; i < genParticles.size(); ++i)
	       {
		   GenParticle const* gen = &genParticles[i];
		   _genPdgId[i]   = gen->pdgId();
		   _genCrCode[i]  = gen->generatorId().id();
		   _genmomX[i]    = gen->momentum().vect().x();
		   _genmomY[i]    = gen->momentum().vect().y();
		   _genmomZ[i]    = gen->momentum().vect().z();
		   _genStartX[i]  = gen->position().x()+ 3904;
		   _genStartY[i]  = gen->position().y();
		   _genStartZ[i]  = gen->position().z();
		   _genStartT[i]  = gen->time();
	       } 
      }
      
      */
       //--------------------------  Do calorimeter hits --------------------------------
      
       _nHits = _nSim = 0;
    
	for (unsigned int ic=0; ic<caloCrystalHits.size();++ic) 
       {	   
	   CaloCrystalHit const& hit      = caloCrystalHits.at(ic);
	   int diskId                     = cal.crystal(hit.id()).diskId();
           CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position());  //in disk FF frame
           //CaloHit const& caloHit         = *(hit.readouts().at(0));
	   //CaloHitSimPartMC const& hitSim = caloHitNavigator.sim(caloHit);
           //int nPartInside                = hitSim.simParticles().size();
	   _cryTime[_nHits]      = hit.time();
	   _cryEdep[_nHits]      = hit.energyDep();
	  
	   _cryPosX[_nHits]      = crystalPos.x();
	   _cryPosY[_nHits]      = crystalPos.y();
	   _cryPosZ[_nHits]      = crystalPos.z();
	   _nHits++;
	   
	}
  	_Ntup->Fill();
  

  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::IPACaloCalibAna);
