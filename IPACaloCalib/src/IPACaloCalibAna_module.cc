
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
	     fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("diag level"),0};
	     //fhicl::Atom<art::InputTag> calohitsTag{Name("CaloHitCollection"), Comment("readout")};
	     fhicl::Atom<art::InputTag> calocrysTag{Name("CaloCrystalHitCollection"),Comment("cal reco crystal hit info")};
	     //fhicl::Atom<art::InputTag> calohitMCcrysTag{Name("CaloHitMCTruthCollection"),Comment("cal hit mc truth info")};
	     fhicl::Atom<art::InputTag> caloclusterTag{Name("CaloClusterCollection"),Comment("cal reco cluster info")};
	     //fhicl::Atom<art::InputTag> caloSimPartMCTag{Name("CaloHitSimPartMCCollection"), Comment("cal hit sim particle mc info")};
    };
       typedef art::EDAnalyzer::Table<Config> Parameters;

       typedef art::Ptr<StepPointMC> StepPtr;
       typedef std::vector<StepPtr>  StepPtrs;
       typedef std::map<int,StepPtrs > HitMap;
       
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
	int _mcdiag;
       //std::string _g4ModuleLabel;
       //std::string _generatorModuleLabel;
       art::InputTag _calohitsTag;
       art::InputTag _calocrysTag;
       art::InputTag _calohitMCcrysTag;
       art::InputTag _caloclusterTag;
       art::InputTag _caloSimPartMCTag;
       const CaloHitCollection* _calhitcol;
       const CaloCrystalHitCollection*  _calcryhitcol;
       const CaloHitMCTruthCollection* _calhitMCtruecol;
       const CaloClusterCollection* _calclustercol;
       const CaloHitSimPartMCCollection* _simpartcol;
       TrkParticle _tpart;
       TrkFitDirection _fdir;
       
       TH2F* _hviewxy;
       TH2F* _hviewxz;
       
       TTree* _Ntup;

       int   _evt,_run;
       int   _nHits;
       int   _nSim;
       int   _nGen;
       int   _cryId[16384],_crySectionId[16384],_crySimIdx[16384],_crySimLen[16384];
       float _cryTime[16384],_cryEdep[16384],_cryDose[16384],_cryPosX[16384],_cryPosY[16384],_cryPosZ[16384],_cryLeak[16384];
       bool findData(const art::Event& evt);
 };


  IPACaloCalibAna::IPACaloCalibAna(const Parameters& conf):
    art::EDAnalyzer(conf),
    _diagLevel(conf().diagLevel()),
    _mcdiag(conf().mcdiag()),
    //_calohitsTag(conf().calohitsTag()),
    _calocrysTag(conf().calocrysTag()),
    //_calohitMCcrysTag(conf().calohitMCcrysTag()),
    _caloclusterTag(conf().caloclusterTag())//,
    //_caloSimPartMCTag(conf().caloSimPartMCTag())
{}

  void IPACaloCalibAna::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    _Ntup  = tfs->make<TTree>("CaloCalibAna", "CaloCalibAna");
    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");
    _Ntup->Branch("nCry",         &_nHits ,       "nCry/I");
    _Ntup->Branch("cryId",        &_cryId ,       "cryId[nCry]/I");
    _Ntup->Branch("crySectionId", &_crySectionId, "crySectionId[nCry]/I");
    _Ntup->Branch("cryPosX",      &_cryPosX ,     "cryPosX[nCry]/F");
    _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosY[nCry]/F");
    _Ntup->Branch("cryPosZ",      &_cryPosZ ,     "cryPosZ[nCry]/F");
    _Ntup->Branch("cryEdep",      &_cryEdep ,     "cryEdep[nCry]/F");
    _Ntup->Branch("cryTime",      &_cryTime ,     "cryTime[nCry]/F");
    _Ntup->Branch("cryDose",      &_cryDose ,     "cryDose[nCry]/F");
    _Ntup->Branch("crySimIdx",    &_crySimIdx ,   "crySimIdx[nCry]/I");
    _Ntup->Branch("crySimLen",    &_crySimLen ,   "crySimLen[nCry]/I");

    _hviewxy = tfs->make<TH2F>("viewxy","StepPoint MC hits dist xy",100,-20.,20.,100,-20.,20.);
    _hviewxz = tfs->make<TH2F>("viewxz","StepPoint MC hits dist xz",230,-10.,220.,100,-20.,20.);
   

  }

  void IPACaloCalibAna::endJob(){
  }




  void IPACaloCalibAna::analyze(const art::Event& event) {
	
       _evt = event.id().event();
       _run = event.run();
       if(!findData(event)) // find data
      		throw cet::exception("RECO")<<"No data in  event"<< endl; 
           
      ++_nProcess;
      if (_nProcess%1000==0) std::cout<<"Processing event "<<_nProcess<<std::endl;
      
      
      //Get handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());
  //Utility to match  cloHits with MCtruth, simParticles and StepPoints
      //CaloHitMCNavigator caloHitNavigator(caloHits, caloHitsMCTruth, caloHitSimPartMC);
   //--------------------------  Do generated particles --------------------------------
 
       //--------------------------  Do calorimeter hits --------------------------------
      
       _nHits = _nSim = 0;
    
	for (unsigned int ic=0; ic<_calcryhitcol->size();++ic) 
       {	   
	   CaloCrystalHit const& hit      = (*_calcryhitcol)[ic];
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

bool IPACaloCalibAna::findData(const art::Event& evt){

	_calhitcol =0;
	_calcryhitcol =0;
	_calclustercol=0;
	//auto calhit = evt.getValidHandle<CaloHitCollection>(_calohitsTag);
	//_calhitcol = calhit.product();
	auto cryhit = evt.getValidHandle<CaloCrystalHitCollection>(_calocrysTag);
	_calcryhitcol =cryhit.product();
	auto cluster= evt.getValidHandle<CaloClusterCollection>(_caloclusterTag);
	_calclustercol =cluster.product();
        if(_mcdiag){
	    
	   _calhitMCtruecol=0;
	   _simpartcol = 0;
           auto simpar= evt.getValidHandle<CaloHitSimPartMCCollection>(_caloSimPartMCTag);
	   _simpartcol =simpar.product();
	   auto calhitMC = evt.getValidHandle<CaloHitMCTruthCollection>(_calohitMCcrysTag);
	   _calhitMCtruecol =calhitMC.product();
           //_toff.updateMap(evt);
        }
	return _calcryhitcol!=0 && _calclustercol !=0 && (_simpartcol != 0 || !_mcdiag);//_calhitcol != 0 && 
       }

//}

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::IPACaloCalibAna);
