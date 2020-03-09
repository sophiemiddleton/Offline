
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

#include <cmath>
#include <iostream>
#include <fstream>
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

  class IPACaloCrystalHitsAna : public art::EDAnalyzer {
     public:
	struct Config {
		     using Name=fhicl::Name;
		     using Comment=fhicl::Comment;
	 	     fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
		     fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
		     fhicl::Atom<art::InputTag> calocrysTag{Name("CaloCrystalHitCollection"),Comment("cal reco crystal hit info")};
		  
	    	};
		typedef art::EDAnalyzer::Table<Config> Parameters;

		explicit IPACaloCrystalHitsAna(const Parameters& conf);
		virtual ~IPACaloCrystalHitsAna() {}


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

		const CaloCrystalHitCollection*  _calcryhitcol;	
		const CaloClusterCollection* _calclustercol;
		const GenParticleCollection *_gencol;

		TTree* _Ntup;

		Int_t   _nEvents = 0;
		Int_t _evt, _run, _nHits;
		
		Int_t   _cryId[16384], _crySectionId[16384], _crySimIdx[16384], _crySimLen[16384];
		Float_t _cryTime[16384], _cryEdep[16384],_cryDose[16384], _cryPosX[16384], _cryPosY[16384], _cryPosZ[16384], _cryLeak[16384], _cryTotE[16384], _cryTotSum[16384], _cryTotEErr[16384], _cryRadius[16384],  _cryMaxR[16384];

		
		bool findData(const art::Event& evt);

	};

   IPACaloCrystalHitsAna::IPACaloCrystalHitsAna(const Parameters& conf):
		art::EDAnalyzer(conf),
		_diagLevel(conf().diagLevel()),
		_mcdiag(conf().mcdiag()),
		_calocrysTag(conf().calocrysTag())	
	{}

  void IPACaloCrystalHitsAna::beginJob(){
	//std::cout<<"[In BeginJob()] Beginning ..."<<std::endl;
	art::ServiceHandle<art::TFileService> tfs;
	_Ntup  = tfs->make<TTree>("CaloCalibAna", "CaloCalibAna");
	_Ntup->Branch("evt",          	&_evt ,        "evt/I");
	_Ntup->Branch("run",          	&_run ,        "run/I");

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
	outputfile.open("IPAAnaCaloCrystalHits.csv");
	
  }


  void IPACaloCrystalHitsAna::analyze(const art::Event& event) {
	//std::cout<<"[In Analyze()] Beginning ..."<<std::endl;
	_evt = event.id().event();
	_run = event.run();

	if(!findData(event)) 
		throw cet::exception("RECO")<<"No data in  event"<< endl; 

	//std::cout<<"[In Analyze()] Found Data ..."<<std::endl;

      	art::ServiceHandle<GeometryService> geom;
      	if( ! geom->hasElement<Calorimeter>() ) return;
      	Calorimeter const & cal = *(GeomHandle<Calorimeter>());	
//=====================Crystal Hits Info =======================//
	//std::cout<<"[In Analyze()] Getting Crystal Info..."<<std::endl;
	_nHits = _calcryhitcol->size();
	for (unsigned int ic=0; ic<_calcryhitcol->size();++ic) 
	{	   
		   CaloCrystalHit const& hit      = (*_calcryhitcol)[ic];
		   int diskId                     = cal.crystal(hit.id()).diskId();
		   CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position());  
		   _cryId[ic] 	 	= hit.id();
		   _cryTime[ic]       	= hit.time();
		   _cryEdep[ic]       	= hit.energyDep();
		   _cryTotE[ic]  		= hit.energyDepTot();
		   _cryTotEErr[ic]  		= hit.energyDepTotErr();
		   _cryPosX[ic]      	= crystalPos.x();
		   _cryPosY[ic]       	= crystalPos.y();
		   _cryPosZ[ic]       	= crystalPos.z();
		   _cryRadius[ic]  	 	= sqrt(crystalPos.x()*crystalPos.x() +  
						    crystalPos.y()*crystalPos.y());
		   outputfile<<_evt<<","<<_run<<","<<hit.id()<<","<<hit.energyDep()<<std::endl;
	}
	_Ntup->Fill();
	_nEvents++;
}


bool IPACaloCrystalHitsAna::findData(const art::Event& evt){
	_calcryhitcol =0;
	auto cryhit = evt.getValidHandle<CaloCrystalHitCollection>(_calocrysTag);
	_calcryhitcol =cryhit.product();
	return   _calcryhitcol!=0;
}

 void IPACaloCrystalHitsAna::endJob(){} 

}

DEFINE_ART_MODULE(mu2e::IPACaloCrystalHitsAna);
