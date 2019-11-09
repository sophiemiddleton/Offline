
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

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"

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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


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
// ROOT incldues
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"
#include "TH3D.h"
#include "Rtypes.h"
#include "TApplication.h"
#include "TArc.h"
#include "TTUBE.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TText.h"
#include "TRotMatrix.h"

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
	     fhicl::Atom<art::InputTag> kalrepTag{Name("KalRepPtrCollection"),Comment("outcome of Kalman filter (for tracker momentum info)")};
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
       art::InputTag _kalrepTag;
       art::InputTag _calohitMCcrysTag;
       art::InputTag _caloclusterTag;
       art::InputTag _caloSimPartMCTag;

       const CaloHitCollection* _calhitcol;
       const CaloCrystalHitCollection*  _calcryhitcol;
       const KalRepPtrCollection* _kalrepcol;
       const CaloHitMCTruthCollection* _calhitMCtruecol;
       const CaloClusterCollection* _calclustercol;
       const CaloHitSimPartMCCollection* _simpartcol;
       TrkParticle _tpart;
       TrkFitDirection _fdir;
       
       TH2F* _hviewxy;
       TH2F* _hviewxz;
       TH1F* _hfit;
       TTree* _Ntup;
       
       int   _evt,_run;
       int   _nEvents, _nHits, _nClusters, _nSim, _nGen;

//Crystel:
       int   _cryId[16384],_crySectionId[16384],_crySimIdx[16384],_crySimLen[16384];
       float _cryTime[16384],_cryEdep[16384],_cryDose[16384],_cryPosX[16384],_cryPosY[16384],_cryPosZ[16384],_cryLeak[16384];
	float _cryTotE[16384];

//Cluster:
       int   _clusterId[16384],_clusterdiskId[16384];
       float _clustertime[16384], _clustertimeErr[16384], _clusterEdep[16384], _clusterEDepErr[16384], _clusterangle[16384], _clustercog3VectorX[16384], _clustercog3VectorY[16384], _clustercog3VectorZ[16384];
	float _EoP[16384];
//Simparticle:
       float _simPartEdep[16384] , _simPartTime[16384], _simPartPosX[16384], _simPartPosY[16384], _simPartPosZ[16384];

//KalFilter:
       float _TrackMom[16382];
       bool findData(const art::Event& evt);
       double CrystalBall(double x, double alpha, double n, double sigma, double mean);
       double CrystalBall(const double *x, const double *p);
       double GetMax(TH1F* h);
       double GetTrackerMom(const art::Event& evt);
};


  IPACaloCalibAna::IPACaloCalibAna(const Parameters& conf):
    art::EDAnalyzer(conf),
    _diagLevel(conf().diagLevel()),
    _mcdiag(conf().mcdiag()),
    //_calohitsTag(conf().calohitsTag()),
    _calocrysTag(conf().calocrysTag()),
    _kalrepTag(conf().kalrepTag()),
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
    _Ntup->Branch("cryTotE",    &_cryTotE ,   "cryTotE[nCry]/F");
    
    _Ntup->Branch("nClu",         &_nClusters ,       "nClu/I");
    _Ntup->Branch("clustertime", &_clustertime,"clustertime[nClu]/F");
    _Ntup->Branch("clustertimeErr", &_clustertimeErr,"clustertimeErr[nClu]/F");
    _Ntup->Branch("clusterEdep", &_clusterEdep,"clusterEdep[nClu]/F");
    _Ntup->Branch("clusterEDepErr", &_clusterEDepErr, "clusterEdepErr[nClu]/F");
    _Ntup->Branch("clusterangle", &_clusterangle, "clusterangle[nClu]/F");
    _Ntup->Branch("clusterPosX",&_clustercog3VectorX, "clustercog3VectorX[nClu]/F");
    _Ntup->Branch("clusterPosY",&_clustercog3VectorY, "clustercog3VectorY[nClu]/F");
    _Ntup->Branch("clusterPosZ",&_clustercog3VectorZ, "clustercog3VectorZ[nClu]/F");

    _Ntup->Branch("trackMom",&_TrackMom, "trackMom/F");
    _Ntup->Branch("EoP", &_EoP,"EoP/F");
    _hviewxy = tfs->make<TH2F>("viewxy","StepPoint MC hits dist xy",100,-20.,20.,100,-20.,20.);
    _hviewxz = tfs->make<TH2F>("viewxz","StepPoint MC hits dist xz",230,-10.,220.,100,-20.,20.);
    _hfit = tfs->make<TH1F>("fitEdep","FitEDep", 100, 0, 100);
   
  }

  void IPACaloCalibAna::endJob(){
	std::cout<<"Max bin"<<GetMax(_hfit)<<std::endl;
  }




double IPACaloCalibAna::CrystalBall(double x, double alpha, double n, double sigma, double mean) {
  // evaluate the crystal ball function
  if (sigma < 0.)     return 0.;
  double z = (x - mean)/sigma; 
  if (alpha < 0) z = -z; 
  double abs_alpha = std::abs(alpha);
  
  if (z  > - abs_alpha)
    return std::exp(- 0.5 * z * z);
  else {
   
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    return AA * std::pow(arg,n);
  }
}

double IPACaloCalibAna::CrystalBall(const double *x, const double *p) {
  return (p[0] * CrystalBall(x[0], p[3], p[4], p[2], p[1]));
}

  void IPACaloCalibAna::analyze(const art::Event& event) {
     /*
 	TF1 *_fitF = new TF1("fitF","CrystalBall");
        _fitF->SetParameter(0,1);
	_fitF->SetParameter(1,1);
	_fitF->SetParameter(2,1);
	_fitF->SetParameter(3,1);
	_fitF->SetParameter(4,1);
*/
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
       double p =GetTrackerMom(event);
	_TrackMom[_evt] = p;
       _nHits = _nSim = 0;
	
        std::cout<<"Number of Clusters"<< _calclustercol->size()<<std::endl;
        for (unsigned int tclu=0; tclu<_calclustercol->size();++tclu){
	   CaloCluster const& cluster = (*_calclustercol)[tclu];
	   //std::cout<<"crystal hits in cluster "<<cluster.caloCrystalHitsPtrVector().size()<<std::endl;

	   _clusterdiskId[_nClusters]                    = cluster.diskId();
           //CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(_clusterdiskId,cal.crystal(cluster.cog3Vector()); 
	   _clustertime[_nClusters]      = cluster.time();
	   _clusterEdep[_nClusters]      = cluster.energyDep();
	   _clusterEDepErr[_nClusters]      = cluster.energyDepErr();
           _clusterangle[_nClusters]      = cluster.angle();
	   _clustercog3VectorX[_nHits]      = cluster.cog3Vector().x();
	   _clustercog3VectorY[_nHits]      = cluster.cog3Vector().y();
	   _clustercog3VectorZ[_nHits]      = cluster.cog3Vector().z();
	   _EoP[_nHits] = cluster.energyDep()/p;
	   _nClusters++;
	}
	std::cout<<"Crystel hit Collection Size "<<_calcryhitcol->size()<<std::endl;
        _cryTotE[_evt] =0;
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
	   _cryTotE[_evt] 	 +=hit.energyDep();
          
	   _cryPosX[_nHits]      = crystalPos.x();
	   _cryPosY[_nHits]      = crystalPos.y();
	   _cryPosZ[_nHits]      = crystalPos.z();
	   _nHits++;
	  
	}
	
	if(_cryTotE[_evt] !=0){
	    _hfit->Fill(_cryTotE[_evt]);
         }

	//_hfit->Fit("fitF");
       
       _Ntup->Fill();
       
       
  }
double IPACaloCalibAna::GetMax(TH1F* h){
	std::cout<<"getting maximum bin "<<std::endl;
	/*int max_bin =0;
	double Xmax_bin = 0;
	for(int i=1;i<h->GetNbinsX();h++){
		
		if(h->GetBinContent(i) > max_bin){
			max_bin = i;
			Xmax_bin = h->GetBinCenter(i);
		}
 
	}
	
	return Xmax_bin;*/
	return h->GetMaximumBin();
}



double IPACaloCalibAna::GetTrackerMom(const art::Event& evt){

	art::ServiceHandle<mu2e::GeometryService>   geom;
  	mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
  	mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;

	double EndMom =0;
	_kalrepcol =0;
	auto kal = evt.getValidHandle<KalRepPtrCollection>(_kalrepTag);
	_kalrepcol =kal.product();
	Hep3Vector vd_tt_back = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_Back));
    	double     Z      = vd_tt_back.z();
    	for(unsigned int i=0;i<_kalrepcol->size();i++){
		art::Ptr<KalRep> const& ptr = _kalrepcol->at(i);
		const KalRep* Krep = ptr.get();
	    	double  ds(10.), s0, s1, s2, z0, z1, z2, dzds, sz, sz1, z01;

		const TrkHitVector* hots = &Krep->hitVector();
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

		z1     = Krep->position(s1).z();
		z2     = Krep->position(s2).z();

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

		z0     = Krep->position(sz).z();     // z0 has to be close to Z(TT_FrontPA)
		z01    = Krep->position(sz+ds).z();

		dzds   = (z01-z0)/ds;
		sz1    = sz+(Z-z0)/dzds;	          // should be good enough
	        EndMom= Krep->momentum(sz1).mag();
		
	}

	return EndMom;
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
        auto kalrep = evt.getValidHandle<KalRepPtrCollection>(_kalrepTag);
	_kalrepcol =kalrep.product();
        if(_mcdiag){
	    
	   _calhitMCtruecol=0;
	   _simpartcol = 0;
           auto simpar= evt.getValidHandle<CaloHitSimPartMCCollection>(_caloSimPartMCTag);
	   _simpartcol =simpar.product();
	   auto calhitMC = evt.getValidHandle<CaloHitMCTruthCollection>(_calohitMCcrysTag);
	   _calhitMCtruecol =calhitMC.product();
           //_toff.updateMap(evt);
        }
	return _kalrepcol!=0 && _calcryhitcol!=0 && _calclustercol !=0 && (_simpartcol != 0 || !_mcdiag);//_calhitcol != 0 && 
       }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::IPACaloCalibAna);
