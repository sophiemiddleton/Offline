//
// An EDAnalyzer module that reads the Trigger Info 
//
// $Id:  $
// $Author:  $
// $Date:  $
//
// Original author G. Pezzullo
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// #include "canvas/Utilities/InputTag.h"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"

//Conditions
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

//Dataproducts
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloTrigSeed.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TriggerAlg.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "DataProducts/inc/XYZVec.hh"

//MC dataproducts
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "TrkDiag/inc/TrkMCTools.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"


//ROOT
#include "TH1F.h"
#include "TH2F.h"

#include <cmath>
// #include <iostream>
#include <string>
// #include <map>
#include <vector>



namespace mu2e {


  class ReadTriggerInfo : public art::EDAnalyzer {

  public:

    enum {
      kNTrigInfo     = 20,
      kNTrackTrig    = 10,
      kNTrackTrigVar = 30,
      kNHelixTrig    = 10,
      kNHelixTrigVar = 30,
      kNCaloCalib    = 5,
      kNCaloCalibVar = 5,
      kNCaloOnly     = 5,
      kNCaloOnlyVar  = 5,
      kNOcc          = 20,
      kNOccVar       = 10
    };

    struct  trigInfo_ {
      int           counts;
      int           exclusive_counts;
      std::string   label;
      
      trigInfo_ ():counts(0), exclusive_counts(0){}
    };

    struct  summaryInfoHist_  {
      TH1F *_hTrigInfo  [kNTrigInfo];
      TH2F *_h2DTrigInfo[kNTrigInfo];
      TH1F *_hTrigBDW   [kNTrigInfo];
      
      summaryInfoHist_() {
	for (int i=0; i<kNTrigInfo; ++i){ 
	  _hTrigInfo  [i] = NULL;
	  _h2DTrigInfo[i] = NULL;
	  _hTrigBDW   [i] = NULL;
	}
      }
    };
    
    struct  trackInfoHist_    {
      TH1F *_hTrkInfo [kNTrackTrig][kNTrackTrigVar];

      trackInfoHist_ (){
	for (int i=0; i<kNTrackTrig; ++i){ 
	  for (int j=0; j<kNTrackTrigVar; ++j){
	    _hTrkInfo  [i][j] = NULL;
	  }
	}     
      }
    };
    
    struct  helixInfoHist_    {
      TH1F *_hHelInfo [kNHelixTrig][kNHelixTrigVar];
      
      helixInfoHist_(){
	for (int i=0; i<kNHelixTrig; ++i){ 
	  for (int j=0; j<kNHelixTrigVar; ++j){
	    _hHelInfo  [i][j] = NULL;
	  }
	}	
      }
    };  

    struct  caloTrigSeedHist_ {
      TH1F *_hCaloOnlyInfo [kNCaloOnly][kNCaloOnlyVar];
      
      caloTrigSeedHist_(){
	for (int i=0; i<kNCaloOnly; ++i){ 
	  for (int j=0; j<kNCaloOnlyVar; ++j){
	    _hCaloOnlyInfo  [i][j] = NULL;
	  }
	}	
      }
    };
    
    struct  caloCalibrationHist_ {
      TH1F *_hCaloCalibInfo[kNCaloCalib][kNCaloCalibVar];
      
      caloCalibrationHist_ (){
	for (int i=0; i<kNCaloCalib; ++i){ 
	  for (int j=0; j<kNCaloCalibVar; ++j){
	    _hCaloCalibInfo  [i][j] = NULL;
	  }
	}
      }
   };
    
    struct  occupancyHist_       {
      TH1F *_hOccInfo  [kNOcc][kNOccVar];
      TH2F *_h2DOccInfo[kNOcc][kNOccVar];
      
      occupancyHist_ (){
	for (int i=0; i<kNOcc; ++i){ 
	  for (int j=0; j<kNOccVar; ++j){
	    _hOccInfo    [i][j] = NULL;
	    _h2DOccInfo  [i][j] = NULL;
	  }
	}
      }
    };
      
    explicit ReadTriggerInfo(fhicl::ParameterSet const& pset);
    virtual ~ReadTriggerInfo() { }

    virtual void beginJob();
    virtual void endJob();
    virtual void endSubRun(const art::SubRun& sr);

    // This is called for each event.
    virtual void analyze(const art::Event& e);
    virtual void beginRun(const art::Run & run);

    void     bookHistograms           ();
    void     bookTrigInfoHist         (art::ServiceHandle<art::TFileService> & Tfs, summaryInfoHist_       &Hist);
    void     bookTrackInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, trackInfoHist_         &Hist);
    void     bookHelixInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, helixInfoHist_         &Hist);
    void     bookCaloTrigSeedInfoHist (art::ServiceHandle<art::TFileService> & Tfs, caloTrigSeedHist_      &Hist);
    void     bookCaloCalibInfoHist    (art::ServiceHandle<art::TFileService> & Tfs, caloCalibrationHist_   &Hist);
    void     bookOccupancyInfoHist    (art::ServiceHandle<art::TFileService> & Tfs, occupancyHist_         &Hist);


    void     findTrigIndex            (std::vector<trigInfo_> Vec, std::string ModuleLabel, int &Index);
    void     fillTrackTrigInfo        (int TrkTrigIndex  , const KalSeed*   KSeed, trackInfoHist_         &Hist);
    void     fillHelixTrigInfo        (int HelTrigIndex  , const HelixSeed* HSeed, helixInfoHist_         &Hist);
    void     fillCaloTrigSeedInfo     (int CTrigSeedIndex, const CaloTrigSeed*HCl, caloTrigSeedHist_      &Hist);
    void     fillCaloCalibTrigInfo    (int ClCalibIndex  , const CaloCluster* HCl, caloCalibrationHist_   &Hist);
    void     fillOccupancyInfo        (int Index         , const StrawDigiCollection*SDCol, const CaloDigiCollection*CDCol, occupancyHist_   &Hist);

    void     findCorrelatedEvents (std::vector<string>& VecLabels, double &NCorrelated);
    void     evalTriggerRate      ();

  private:

    int             _diagLevel;
    size_t          _nMaxTrig;    
    int             _nTrackTrig;
    int             _nCaloTrig;
    int             _nCaloCalibTrig;
    art::InputTag   _trigAlgTag;
    art::InputTag   _sdMCTag;
    art::InputTag   _sdTag;
    art::InputTag   _chTag;    
    art::InputTag   _cdTag;
    art::InputTag   _evtWeightTag;
    double          _duty_cycle;

    int             _nProcess;
    int             _numEvents;    
    double          _bz0;
    double          _nPOT;
    
    std::vector<trigInfo_>    _trigAll;	     
    std::vector<trigInfo_>    _trigFinal;    
    std::vector<trigInfo_>    _trigCaloOnly; 
    std::vector<trigInfo_>    _trigCaloCalib;
    std::vector<trigInfo_>    _trigTrack;    
    std::vector<trigInfo_>    _trigHelix;    
    std::vector<trigInfo_>    _trigEvtPS;    
    
    summaryInfoHist_          _sumHist;
    trackInfoHist_            _trkHist;
    helixInfoHist_            _helHist;

    caloTrigSeedHist_         _caloTSeedHist;
    caloCalibrationHist_      _caloCalibHist;
    occupancyHist_            _occupancyHist;

    //the following pointer is needed to navigate the MC truth info of the strawHits
    const mu2e::StrawDigiMCCollection* _mcdigis;
    const mu2e::ComboHitCollection*    _chcol;
    const art::Event*                  _event;
  };


  ReadTriggerInfo::ReadTriggerInfo(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset), 
    _diagLevel     (pset.get<int>   ("diagLevel", 0)),
    _nMaxTrig      (pset.get<size_t>("nFilters", 70)),
    _nTrackTrig    (pset.get<size_t>("nTrackTriggers", 4)),
    _nCaloTrig     (pset.get<size_t>("nCaloTriggers", 4)),
    _nCaloCalibTrig(pset.get<size_t>("nCaloCalibTriggers", 4)),
    _trigAlgTag    (pset.get<art::InputTag>("triggerAlgMerger"     , "triggerInfoMerger")),
    _sdMCTag       (pset.get<art::InputTag>("strawDigiMCCollection", "compressDigiMCs")),
    _sdTag         (pset.get<art::InputTag>("strawDigiCollection"  , "makeSD")),
    _chTag         (pset.get<art::InputTag>("comboHitCollection"   , "TTmakeSH")),
    _cdTag         (pset.get<art::InputTag>("caloDigiCollection"   , "CaloDigiFromShower")),
    _evtWeightTag  (pset.get<art::InputTag>("protonBunchIntensity" , "protonBunchIntensity")),
    _duty_cycle    (pset.get<float> ("dutyCycle", 1.)),
    _nProcess(0), 
    _numEvents(0)
  {
    _trigAll.      resize(_nMaxTrig);	     
    _trigFinal.    resize(_nMaxTrig);    
    _trigCaloOnly. resize(_nMaxTrig); 
    _trigCaloCalib.resize(_nMaxTrig);
    _trigTrack.    resize(_nMaxTrig);    
    _trigHelix.    resize(_nMaxTrig);    
    _trigEvtPS.    resize(_nMaxTrig);    
  }
  
  
  void     ReadTriggerInfo::bookHistograms           (){
    art::ServiceHandle<art::TFileService> tfs;
    
    bookTrigInfoHist(tfs, _sumHist);

    bookTrackInfoHist(tfs, _trkHist);
    
    bookHelixInfoHist(tfs, _helHist);
    
    bookCaloTrigSeedInfoHist(tfs, _caloTSeedHist);
    
    bookCaloCalibInfoHist(tfs, _caloCalibHist);

    bookOccupancyInfoHist(tfs, _occupancyHist);
  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookTrigInfoHist         (art::ServiceHandle<art::TFileService> & Tfs, summaryInfoHist_       &Hist){
    art::TFileDirectory trigInfoDir = Tfs->mkdir("trigInfo");

    Hist._hTrigInfo[0]   = trigInfoDir.make<TH1F>("hTrigInfo_global"    , "Global Trigger rejection"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    Hist._hTrigInfo[1]   = trigInfoDir.make<TH1F>("hTrigInfo_track"     , "Calo-only Triggers rejection"               , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    Hist._hTrigInfo[2]   = trigInfoDir.make<TH1F>("hTrigInfo_calo"      , "Track Triggers rejection"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    Hist._hTrigInfo[3]   = trigInfoDir.make<TH1F>("hTrigInfo_evtPS"     , "Event prescaler Trigger bits distribution"  , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    Hist._hTrigInfo[4]   = trigInfoDir.make<TH1F>("hTrigInfo_helix"     , "HelixSeed Triggers rejection"               , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    Hist._hTrigInfo[5]   = trigInfoDir.make<TH1F>("hTrigInfo_caloCalib" , "Calo Calibration rejection"                 , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    Hist._hTrigInfo[6]   = trigInfoDir.make<TH1F>("hTrigInfo_final"     , "Global Trigger rejection of the paths"      , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       

    Hist._hTrigInfo[10]  = trigInfoDir.make<TH1F>("hTrigInfo_unique_all", "Events found only by each Trig path"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    Hist._hTrigInfo[11]  = trigInfoDir.make<TH1F>("hTrigInfo_unique"    , "Events found only by each Trig path"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       

    Hist._hTrigInfo[15]  = trigInfoDir.make<TH1F>("hTrigInfo_paths"     , "Rejection of all the Trigger paths"         , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       


    Hist._h2DTrigInfo[0] = trigInfoDir.make<TH2F>("h2DTrigInfo_map_all" , "Trigger correlation map from all filters"   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5), (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    Hist._h2DTrigInfo[1] = trigInfoDir.make<TH2F>("h2DTrigInfo_map"     , "Trigger correlation map"                    , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5), (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));   

    art::TFileDirectory trigBDWDir = Tfs->mkdir("trigBDW");

    Hist._hTrigBDW[0]   = trigBDWDir.make<TH1F>("hTrigBDW_global"    , "Trigger bandwidth; ; rate [Hz]"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    Hist._hTrigBDW[1]   = trigBDWDir.make<TH1F>("hTrigBDW_cumulative", "Cumulative Trigger bandwidth; ; rate [Hz]"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    
  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookTrackInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, trackInfoHist_         &Hist){
    for (int i=0; i<_nTrackTrig; ++i){
      art::TFileDirectory trkInfoDir  = Tfs->mkdir(Form("trk_%i",i));
      Hist._hTrkInfo[i][0] = trkInfoDir.make<TH1F>(Form("hP_%i"    , i), "Track Momentum; p[MeV/c]", 400, 0, 200);
      Hist._hTrkInfo[i][1] = trkInfoDir.make<TH1F>(Form("hPt_%i"   , i), "Track Pt; p_{t} [MeV/c]", 400, 0, 200);
      Hist._hTrkInfo[i][2] = trkInfoDir.make<TH1F>(Form("hNSh_%i"  , i), "N-StrawHits; nStrawHits", 101, -0.5, 100.5);
      Hist._hTrkInfo[i][3] = trkInfoDir.make<TH1F>(Form("hD0_%i"   , i), "Track impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hTrkInfo[i][4] = trkInfoDir.make<TH1F>(Form("hChi2d_%i", i), "Track #chi^{2}/ndof;#chi^{2}/ndof", 100, 0, 50);
      Hist._hTrkInfo[i][5] = trkInfoDir.make<TH1F>(Form("hClE_%i"  , i), "calorimeter Cluster energy; E [MeV]", 240, 0, 120);
      
      Hist._hTrkInfo[i][10] = trkInfoDir.make<TH1F>(Form("hPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hTrkInfo[i][11] = trkInfoDir.make<TH1F>(Form("hPtMC_%i", i), "MC Track Pt @ tracker front; p_{t} [MeV/c]" , 400, 0, 200);
      Hist._hTrkInfo[i][12] = trkInfoDir.make<TH1F>(Form("hPzMC_%i", i), "MC Track Pt @ tracker front; p_{t} [MeV/c]" , 400, 0, 200);
      Hist._hTrkInfo[i][13] = trkInfoDir.make<TH1F>(Form("hDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{trk} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hTrkInfo[i][14] = trkInfoDir.make<TH1F>(Form("hDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{trk} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hTrkInfo[i][15] = trkInfoDir.make<TH1F>(Form("hDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{trk} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hTrkInfo[i][16] = trkInfoDir.make<TH1F>(Form("hPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hTrkInfo[i][17] = trkInfoDir.make<TH1F>(Form("hGenZ_%i" , i), "z origin; z-origin [mm]"              , 300,   0,   15000);
      Hist._hTrkInfo[i][18] = trkInfoDir.make<TH1F>(Form("hGenR_%i" , i), "radial position origin; r-origin [mm]", 500,   0,   5000);
      Hist._hTrkInfo[i][19] = trkInfoDir.make<TH1F>(Form("hPDGM_%i" , i), "PDG Mother Id; PdgId-mother"          , 2253,   -30.5,   2222.5);
      Hist._hTrkInfo[i][20] = trkInfoDir.make<TH1F>(Form("hEMC_%i"  , i), "MC Energy; E_{MC} [MeV]"              , 400,   0,   200);


  
    }
  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookHelixInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, helixInfoHist_         &Hist){
    for (int i=0; i<_nTrackTrig; ++i){
      art::TFileDirectory helInfoDir  = Tfs->mkdir(Form("helix_%i", i));
      Hist._hHelInfo[i][0] = helInfoDir.make<TH1F>(Form("hP_%i"        , i), "Helix Momentum; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][1] = helInfoDir.make<TH1F>(Form("hPt_%i"       , i), "Helix Pt; p_{t} [MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][2] = helInfoDir.make<TH1F>(Form("hNSh_%i"      , i), "N-StrawHits; nStrawHits", 101, -0.5, 100.5);
      Hist._hHelInfo[i][3] = helInfoDir.make<TH1F>(Form("hD0_%i"       , i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][4] = helInfoDir.make<TH1F>(Form("hChi2dXY_%i"  , i), "Helix #chi^{2}_{xy}/ndof;#chi^{2}_{xy}/ndof"      , 100, 0, 50);
      Hist._hHelInfo[i][5] = helInfoDir.make<TH1F>(Form("hChi2dZPhi_%i", i), "Helix #chi^{2}_{z#phi}/ndof;#chi^{2}_{z#phi}/ndof", 100, 0, 50);
      Hist._hHelInfo[i][6] = helInfoDir.make<TH1F>(Form("hClE_%i"      , i), "calorimeter Cluster energy; E [MeV]", 240, 0, 120);
      Hist._hHelInfo[i][7] = helInfoDir.make<TH1F>(Form("hLambda_%i"   , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);

        
      Hist._hHelInfo[i][10] = helInfoDir.make<TH1F>(Form("hPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][11] = helInfoDir.make<TH1F>(Form("hPtMC_%i", i), "MC Track Pt @ tracker front; p_{t} [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][12] = helInfoDir.make<TH1F>(Form("hPzMC_%i", i), "MC Track Pt @ tracker front; p_{z} [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][13] = helInfoDir.make<TH1F>(Form("hDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][14] = helInfoDir.make<TH1F>(Form("hDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][15] = helInfoDir.make<TH1F>(Form("hDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][16] = helInfoDir.make<TH1F>(Form("hPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][17] = helInfoDir.make<TH1F>(Form("hGenZ_%i" , i), "z origin; z-origin [mm]"              , 300,   0,   15000);
      Hist._hHelInfo[i][18] = helInfoDir.make<TH1F>(Form("hGenR_%i" , i), "radial position origin; r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][19] = helInfoDir.make<TH1F>(Form("hPDGM_%i" , i), "PDG Mother Id; PdgId-mother"          , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][20] = helInfoDir.make<TH1F>(Form("hEMC_%i"  , i), "MC Energy; E_{MC} [MeV]"              , 400,   0,   200);

    }
  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookCaloTrigSeedInfoHist (art::ServiceHandle<art::TFileService> & Tfs, caloTrigSeedHist_      &Hist){
    for (int i=0; i<_nCaloTrig; ++i){
      art::TFileDirectory caloInfoDir = Tfs->mkdir(Form("caloOnly_%i",i));
      Hist._hCaloOnlyInfo[i][0] = caloInfoDir.make<TH1F>(Form("hEPeak_%i"   , i), "peak energy; E[MeV]"        , 400, 0, 200);
      Hist._hCaloOnlyInfo[i][1] = caloInfoDir.make<TH1F>(Form("hR1Max1_%i"   , i), "ring1 max; ring1max [MeV]" , 400, 0, 200);
      Hist._hCaloOnlyInfo[i][2] = caloInfoDir.make<TH1F>(Form("hR1Max2_%i"   , i), "ring1 max; ring1max2 [MeV]", 400, 0, 200);    
    }
  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookCaloCalibInfoHist    (art::ServiceHandle<art::TFileService> & Tfs, caloCalibrationHist_   &Hist){      
    for (int i=0; i<_nCaloCalibTrig; ++i){
      art::TFileDirectory caloCalibInfoDir = Tfs->mkdir(Form("caloCalib_%i",i));
      Hist._hCaloCalibInfo[i][0] = caloCalibInfoDir.make<TH1F>(Form("hE_%i"   , i), "Cluster energy; E[MeV]", 800, 0, 800);
      Hist._hCaloCalibInfo[i][1] = caloCalibInfoDir.make<TH1F>(Form("hN_%i"   , i), "Cluster size; nCrystalHits", 101, -0.5, 100.5);
    }    
  }

  //--------------------------------------------------------------------------------//
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookOccupancyInfoHist         (art::ServiceHandle<art::TFileService> & Tfs, occupancyHist_       &Hist){
    
    for (int i=0; i<_nTrackTrig; ++i){
      art::TFileDirectory occInfoDir = Tfs->mkdir(Form("occInfoTrk_%i", i));
      Hist._hOccInfo  [i][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,i),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, 1e6, 4e8);
      
      Hist._h2DOccInfo[i][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,i),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, 1e6, 4e8, 5000, 0., 20000.);
      Hist._h2DOccInfo[i][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,i),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, 1e6, 4e8, 5000, 0., 20000.);
    }
    
    for (int i=_nTrackTrig; i<_nTrackTrig*2; ++i){
      art::TFileDirectory occInfoDir = Tfs->mkdir(Form("occInfoHel_%i", i));
      Hist._hOccInfo  [i][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,i),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, 1e6, 4e8);
      
      Hist._h2DOccInfo[i][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,i),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, 1e6, 4e8, 5000, 0., 20000.);
      Hist._h2DOccInfo[i][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,i),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, 1e6, 4e8, 5000, 0., 20000.);
    }
    
    for (int i=_nTrackTrig*2; i<_nTrackTrig*2+_nCaloTrig; ++i){
      art::TFileDirectory occInfoDir = Tfs->mkdir(Form("occInfoCaloTrig_%i", i));
      Hist._hOccInfo  [i][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,i),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, 1e6, 4e8);
      			    
      Hist._h2DOccInfo[i][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,i),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, 1e6, 4e8, 5000, 0., 20000.);
      Hist._h2DOccInfo[i][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,i),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, 1e6, 4e8, 5000, 0., 20000.);
    }
    
     int    index_last = _nTrackTrig+_nCaloTrig;
     art::TFileDirectory occInfoDir = Tfs->mkdir("occInfoGeneral");
     Hist._hOccInfo  [index_last][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,index_last),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, 1e6, 4e8);
      
     Hist._h2DOccInfo[index_last][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,index_last),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, 1e6, 4e8, 5000, 0., 20000.);
     Hist._h2DOccInfo[index_last][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,index_last),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, 1e6, 4e8, 5000, 0., 20000.);
 
	
    
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::beginJob(){

    bookHistograms();
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::endJob(){

    //set hitograms' titles

    //Helix
    for (int i=0; i<_nTrackTrig; ++i){
      for (int j=0; j<kNHelixTrigVar; ++j){
	if (_helHist._hHelInfo[i][j] == NULL)    continue;
	string title = _trigHelix[i].label +": "+ _helHist._hHelInfo[i][j]->GetTitle();
	_helHist._hHelInfo[i][j]->SetTitle(title.c_str());
      }
    }
    
    //Tracks
    for (int i=0; i<_nTrackTrig; ++i){
      for (int j=0; j<kNTrackTrigVar; ++j){
	if (_trkHist._hTrkInfo[i][j] == NULL)    continue;
	string title = _trigTrack[i].label +": "+ _trkHist._hTrkInfo[i][j]->GetTitle();
	_trkHist._hTrkInfo[i][j]->SetTitle(title.c_str());
      }
    }

 
    //occupancy
    //tracks
    for (int i=0; i<_nTrackTrig; ++i){
      for (int j=0; j<kNOccVar; ++j){
	if (_occupancyHist._hOccInfo[i][j] != NULL)    {
	  string title = _trigTrack[i].label +": "+ _occupancyHist._hOccInfo[i][j]->GetTitle();
	  _occupancyHist._hOccInfo[i][j]->SetTitle(title.c_str());
	}
	if (_occupancyHist._h2DOccInfo[i][j] != NULL)    {
	  string title = _trigTrack[i].label +": "+ _occupancyHist._h2DOccInfo[i][j]->GetTitle();
	  _occupancyHist._h2DOccInfo[i][j]->SetTitle(title.c_str());
	}
      }
    }
    //helix
    for (int i=_nTrackTrig; i<_nTrackTrig*2; ++i){
      for (int j=0; j<kNOccVar; ++j){
	if (_occupancyHist._hOccInfo[i][j] != NULL)    {
	  string title = _trigHelix[i].label +": "+ _occupancyHist._hOccInfo[i][j]->GetTitle();
	  _occupancyHist._hOccInfo[i][j]->SetTitle(title.c_str());
	}
	if (_occupancyHist._h2DOccInfo[i][j] != NULL)    {
	  string title = _trigHelix[i].label +": "+ _occupancyHist._h2DOccInfo[i][j]->GetTitle();
	  _occupancyHist._h2DOccInfo[i][j]->SetTitle(title.c_str());
	}
      }
    }
    //calo trig
    for (int i=_nTrackTrig*2; i<_nTrackTrig*2+_nCaloTrig; ++i){
      for (int j=0; j<kNOccVar; ++j){
	if (_occupancyHist._hOccInfo[i][j] != NULL)    {
	  string title = _trigCaloOnly[i].label +": "+ _occupancyHist._hOccInfo[i][j]->GetTitle();
	  _occupancyHist._hOccInfo[i][j]->SetTitle(title.c_str());
	}
	if (_occupancyHist._h2DOccInfo[i][j] != NULL)    {
	  string title = _trigCaloOnly[i].label +": "+ _occupancyHist._h2DOccInfo[i][j]->GetTitle();
	  _occupancyHist._h2DOccInfo[i][j]->SetTitle(title.c_str());
	}
      }   
    }


    int    indexTrigInfo11(0);

    //fill the histograms
    for (size_t i=0; i<_trigAll.size(); ++i ){
      _sumHist._hTrigInfo  [0]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());
      _sumHist._h2DTrigInfo[0]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());

      if (_trigAll[i].counts > 0) {
	_sumHist._hTrigInfo[0]->SetBinContent(i+1, _nProcess/_trigAll[i].counts);
	for (size_t j=0; j<_trigAll.size(); ++j ){
	  _sumHist._h2DTrigInfo[0]->GetYaxis()->SetBinLabel(j+1, _trigAll[j].label.c_str());
	}
      }

      _sumHist._hTrigInfo[1]->GetXaxis()->SetBinLabel(i+1, _trigTrack[i].label.c_str());
      if (_trigTrack[i].counts > 0) _sumHist._hTrigInfo[1]->SetBinContent(i+1, _nProcess/_trigTrack[i].counts);

      _sumHist._hTrigInfo[2]->GetXaxis()->SetBinLabel(i+1, _trigCaloOnly[i].label.c_str());
      if (_trigCaloOnly[i].counts > 0) _sumHist._hTrigInfo[2]->SetBinContent(i+1, _nProcess/_trigCaloOnly[i].counts);

      _sumHist._hTrigInfo[3]->GetXaxis()->SetBinLabel(i+1, _trigEvtPS[i].label.c_str());
      if (_trigEvtPS[i].counts > 0) _sumHist._hTrigInfo[3]->SetBinContent(i+1, _trigEvtPS[i].counts);

      _sumHist._hTrigInfo[4]->GetXaxis()->SetBinLabel(i+1, _trigHelix[i].label.c_str());
      if (_trigHelix[i].counts > 0) _sumHist._hTrigInfo[4]->SetBinContent(i+1, _trigHelix[i].counts);

      _sumHist._hTrigInfo[5]->GetXaxis()->SetBinLabel(i+1, _trigCaloCalib[i].label.c_str());
      if (_trigCaloCalib[i].counts > 0) _sumHist._hTrigInfo[5]->SetBinContent(i+1, _trigCaloCalib[i].counts);

      if (_trigFinal[i].counts > 0) {
	_sumHist._hTrigInfo  [6]->GetXaxis()->SetBinLabel(i+1, _trigFinal[i].label.c_str());
	_sumHist._hTrigInfo  [6]->SetBinContent(i+1, _nProcess/_trigFinal[i].counts);
      }

      //fill  the histograms that shows how many events were found exclusively by each trigger path
      _sumHist._hTrigInfo  [10]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());
      double    content_trigInfo11 = _sumHist._hTrigInfo [10]->GetBinContent(i+1);
      if (content_trigInfo11>0){
	_sumHist._hTrigInfo  [11]->GetXaxis()->SetBinLabel(indexTrigInfo11 +1, _trigAll[i].label.c_str());
	_sumHist._hTrigInfo  [11]->SetBinContent(indexTrigInfo11 +1, content_trigInfo11);
	++indexTrigInfo11;
      }

    }

    //now let's filter the 2D correlation histogram with only those that actually triggered at least one event
    int                nbinsx = _sumHist._h2DTrigInfo[0]->GetNbinsX();
    int                nbinsy = _sumHist._h2DTrigInfo[0]->GetNbinsY();
    std::vector<int>   binsToSkip;

    for (int i=0; i<nbinsx; ++i){
      bool used(false);

      for (int j=0; j<nbinsy; ++j){
	if (_sumHist._h2DTrigInfo[0]->GetBinContent(i+1, j+1) > 0) {
	  used = true;
	  break;
	}
      }
      if (!used) binsToSkip.push_back(i);
    }

    int   index_x(0);
    for (int i=0; i<nbinsx; ++i){
      int    counts = std::count(binsToSkip.begin(), binsToSkip.end(), i);
      if (counts >= 1)       continue;
      //set the label
      _sumHist._h2DTrigInfo[1]->GetXaxis()->SetBinLabel(index_x+1, _trigAll[i].label.c_str());

      int    index_y(0);

      for (int j=0; j<nbinsy; ++j){
	counts = std::count(binsToSkip.begin(), binsToSkip.end(), j);
	if (counts >= 1)       continue;
	double  content =  _sumHist._h2DTrigInfo[0]->GetBinContent(i+1, j+1);
	_sumHist._h2DTrigInfo[1]->SetBinContent(index_x+1, index_y+1, content);
	
	//set the label
	if (index_x == 0){
	  _sumHist._h2DTrigInfo[1]->GetYaxis()->SetBinLabel(index_y+1, _trigAll[j].label.c_str());
	}

	++index_y;
      }
      ++index_x;
    }
    
    // now evaluate the bandwidth
    // NOTE: "evalTriggerrate" re-order the vectors _trigFinal
    evalTriggerRate();
  }
  
  //--------------------------------------------------------------------------------  
  void   ReadTriggerInfo::evalTriggerRate        (){
    //order the array with the filter used at the end of each path
    std::sort(_trigFinal.begin(), _trigFinal.end(), [](const auto a, const auto b) {return a.counts < b.counts; });
    
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    double    mbtime         = accPar->deBuncherPeriod;
    double    mean_mb_rate   = 1./(mbtime/CLHEP::s)*_duty_cycle;

    bool      isFirst(true);
    int       index(0);
   
    std::vector<string>    labels_by_rate;

    for (size_t i=0; i< _trigFinal.size(); ++i){
      double  nEvents = (double)_trigFinal[i].counts;
      if ( nEvents <= 1e-3)                 continue;

      labels_by_rate.push_back(_trigFinal[i].label);

      double  eff   = nEvents/(double)_nProcess;
      double  rate  = mean_mb_rate*eff;
      _sumHist._hTrigBDW[0]->GetXaxis()->SetBinLabel(index+1, _trigFinal[i].label.c_str());
      _sumHist._hTrigBDW[1]->GetXaxis()->SetBinLabel(index+1, _trigFinal[i].label.c_str());
      _sumHist._hTrigBDW[0]->SetBinContent(index+1, rate);
      
      if (isFirst) {
      	_sumHist._hTrigBDW[1]->SetBinContent(index+1, rate);
      	//	rate_ref  = rate;
      	isFirst = false;
      }else{
	double    nCorrelated(0);
	findCorrelatedEvents(labels_by_rate, nCorrelated);

	rate = _sumHist._hTrigBDW[1]->GetBinContent(index) + (nEvents-nCorrelated)/(double)_nProcess*mean_mb_rate;
      	_sumHist._hTrigBDW[1]->SetBinContent(index+1, rate);
      }

      ++index;
    }

  }



  void   ReadTriggerInfo::findCorrelatedEvents(std::vector<string>& VecLabels, double &NCorrelated){
    
    NCorrelated = 0;
    
    const char* label_ref = VecLabels.at(VecLabels.size()-1).c_str();

    //    char* label(0);
    
    int        nLabels = VecLabels.size() -1;
    int        nbins   = _sumHist._h2DTrigInfo[1]->GetNbinsX();
    for(int i=0; i<nbins; ++i){
      const char* label = _sumHist._h2DTrigInfo[1]->GetXaxis()->GetBinLabel(i+1);
      if (std::strcmp(label_ref, label) != 0)      continue;
      
      for (int k=0; k<nLabels; ++k){
	label_ref = VecLabels.at(k).c_str();
	for (int j=0; j<nbins; ++j){
	  if (j == i)      break;
	  label =   _sumHist._h2DTrigInfo[1]->GetYaxis()->GetBinLabel(j+1);
	  if (std::strcmp(label_ref, label) != 0)        continue;
	  NCorrelated += _sumHist._h2DTrigInfo[1]->GetBinContent(i+1, j+1);
	}
      }
      break;
    }
    
    
    
  }

  //================================================================
  void   ReadTriggerInfo::beginRun(const art::Run & run){
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();
  }

  void ReadTriggerInfo::endSubRun(const art::SubRun& sr){}

  void ReadTriggerInfo::analyze(const art::Event& event) {

    ++_nProcess;

    //get the number of POT
    _nPOT  = -1.;
    art::Handle<ProtonBunchIntensity> evtWeightH;
    event.getByLabel(_evtWeightTag, evtWeightH);
    if (evtWeightH.isValid()){
      _nPOT  = (double)evtWeightH->intensity();
    }

    std::vector<art::Handle<TriggerInfo> > hTrigInfoVec;
    event.getManyByType(hTrigInfoVec);

    //get the TriggerAlg from the Event
    art::Handle<mu2e::TriggerAlg> trigAlgH;
    event.getByLabel(_trigAlgTag, trigAlgH);
    const mu2e::TriggerAlg*       trigAlg(0);
    if (trigAlgH.isValid()){
      trigAlg = trigAlgH.product();
    }

    //get the strawDigiMC truth if present
    art::Handle<mu2e::StrawDigiMCCollection> mcdH;
    event.getByLabel(_sdMCTag, mcdH);
    if (mcdH.isValid()) {
      _mcdigis = mcdH.product();
      _event   = &event;
    }else {
      _mcdigis = NULL;
    }

    //get the StrawDigi Collection
    art::Handle<mu2e::StrawDigiCollection> sdH;
    event.getByLabel(_sdTag, sdH);
    const StrawDigiCollection* sdCol(0);
    if (sdH.isValid()) {
      sdCol = sdH.product();
    }

    //get the ComboHitCollection
    art::Handle<mu2e::ComboHitCollection> chH;
    event.getByLabel(_chTag, chH);
    if (chH.isValid()) {
      _chcol = chH.product();
    }else {
      _chcol = NULL;
    }

    //get the CaloDigi Collection
    art::Handle<mu2e::CaloDigiCollection> cdH;
    event.getByLabel(_cdTag, cdH);
    const CaloDigiCollection* cdCol(0);
    if (cdH.isValid()) {
      cdCol = cdH.product();
    }

    if (trigAlg != 0) {
      if (trigAlg->hasAnyProperty(TriggerAlg::caloCalibCosmic)) _sumHist._hTrigInfo[15]->Fill(0.);
      if (trigAlg->hasAnyProperty(TriggerAlg::caloMVACE))       _sumHist._hTrigInfo[15]->Fill(1);
      if (trigAlg->hasAnyProperty(TriggerAlg::tprSeedDeM))      _sumHist._hTrigInfo[15]->Fill(2.);
      if (trigAlg->hasAnyProperty(TriggerAlg::tprSeedDeP))      _sumHist._hTrigInfo[15]->Fill(3.);
      if (trigAlg->hasAnyProperty(TriggerAlg::cprSeedDeM))      _sumHist._hTrigInfo[15]->Fill(4.);
      if (trigAlg->hasAnyProperty(TriggerAlg::cprSeedDeP))      _sumHist._hTrigInfo[15]->Fill(5.);
    }

    
    //fill the general occupancy histogram
    fillOccupancyInfo   (_nTrackTrig+_nCaloTrig, sdCol, cdCol, _occupancyHist);

    art::Handle<TriggerInfo>       hTrigInfo;
    TriggerFlag                    prescalerFlag       = TriggerFlag::prescaleRandom;
    TriggerFlag                    trackFlag           = TriggerFlag::track;
    TriggerFlag                    helixFlag           = TriggerFlag::helix;
    TriggerFlag                    caloFlag            = TriggerFlag::caloCluster;
    TriggerFlag                    caloCalibFlag       = TriggerFlag::caloCalib;
    TriggerFlag                    caloTrigSeedFlag    = TriggerFlag::caloTrigSeed;
    TriggerFlag                    caloOrTrackFlag     = trackFlag; caloOrTrackFlag.merge(caloFlag); caloOrTrackFlag.merge(caloCalibFlag); caloOrTrackFlag.merge(caloTrigSeedFlag);// caloOrTrackFlag.merge(helixFlag);
    
    std::vector<int>   trigFlagAll_index, trigFlag_index;
    
    for (size_t i=0; i<hTrigInfoVec.size(); ++i){
      hTrigInfo = hTrigInfoVec.at(i);
      if (!hTrigInfo.isValid())         continue;
      const TriggerInfo* trigInfo  = hTrigInfo.product();
      const TriggerFlag  flag      = trigInfo->triggerBits();

      std::string    moduleLabel   = hTrigInfo.provenance()->moduleLabel();
      int            index_all(0);         
      int            index(0);         
      
      //fill the Global Trigger bits info
      findTrigIndex(_trigAll, moduleLabel, index_all);
      _trigAll[index_all].label  = moduleLabel;
      //      if ( flag.hasAnyProperty(caloOrTrackFlag)){ 
	//      }

      findTrigIndex(_trigFinal, moduleLabel, index);
      if ( flag.hasAnyProperty(caloOrTrackFlag)){ 
	_trigFinal[index].label  = moduleLabel;
	_trigFinal[index].counts = _trigFinal[index].counts + 1;
	_trigAll[index_all].counts = _trigAll[index_all].counts + 1;
	trigFlagAll_index.push_back(index_all);
      }
      //fill the Calo-Only Trigger bits info
      findTrigIndex(_trigCaloOnly, moduleLabel, index);
      if ( flag.hasAnyProperty(caloFlag) || flag.hasAnyProperty(caloTrigSeedFlag)){ 
	_trigCaloOnly[index].label  = moduleLabel;
	_trigCaloOnly[index].counts = _trigCaloOnly[index].counts + 1;
	const CaloTrigSeed*clseed = trigInfo->caloTrigSeed().get();
	if(clseed) {
	  fillCaloTrigSeedInfo(index, clseed, _caloTSeedHist);
	  fillOccupancyInfo   (_nTrackTrig*2+index, sdCol, cdCol, _occupancyHist);
	}
	trigFlag_index.push_back(index_all);
      }

      //fill the CaloCalib Trigger bits info
      findTrigIndex(_trigCaloCalib, moduleLabel, index);
      if ( flag.hasAnyProperty(caloCalibFlag)){ 
	_trigCaloCalib[index].label  = moduleLabel;
	_trigCaloCalib[index].counts = _trigCaloCalib[index].counts + 1;
	const CaloCluster*cluster = trigInfo->caloCluster().get();
	if(cluster) fillCaloCalibTrigInfo(index, cluster, _caloCalibHist);
	trigFlag_index.push_back(index_all);
      }
      
      //fill the Track Trigger bits info
      findTrigIndex(_trigTrack, moduleLabel, index);
      if ( flag.hasAnyProperty(trackFlag)){ 
	_trigTrack[index].label  = moduleLabel;
	_trigTrack[index].counts = _trigTrack[index].counts + 1;
	const KalSeed*kseed = trigInfo->track().get();
	if(kseed) {
	  fillTrackTrigInfo(index, kseed, _trkHist);
	  fillOccupancyInfo(index, sdCol, cdCol, _occupancyHist);
	}
	trigFlag_index.push_back(index_all);
      }
       //fill the Helix Trigger bits info
      findTrigIndex(_trigHelix, moduleLabel, index);
      if ( flag.hasAnyProperty(helixFlag)){ 
	_trigHelix[index].label  = moduleLabel;
	_trigHelix[index].counts = _trigHelix[index].counts + 1;
	const HelixSeed*hseed = trigInfo->helix().get();
	if(hseed) {
	  fillHelixTrigInfo(index, hseed, _helHist);
	  fillOccupancyInfo(_nTrackTrig+index, sdCol, cdCol, _occupancyHist);
	}
      }
      
      //fill the Event-Prescaler Trigger bits info
      findTrigIndex(_trigEvtPS, moduleLabel, index);
      if ( flag.hasAnyProperty(prescalerFlag)){ 
	_trigEvtPS[index].label  = moduleLabel;
	_trigEvtPS[index].counts = _trigEvtPS[index].counts + 1;
      }
      
      
    }//end loop over the TriggerInfo Handles

    //now fill the correlation matrix
    for (size_t i=0; i<trigFlagAll_index.size(); ++i){
      for (size_t j=0; j<trigFlagAll_index.size(); ++j){
	_sumHist._h2DTrigInfo[0]->Fill(trigFlagAll_index.at(i), trigFlagAll_index.at(j));
      }
    }
    
    if (trigFlagAll_index.size() == 1) _sumHist._hTrigInfo[10]->Fill(trigFlagAll_index.at(0));

  }
  
  void   ReadTriggerInfo::findTrigIndex(std::vector<trigInfo_> Vec, std::string ModuleLabel, int &Index){
    //reset the index value
    Index = 0;
    for (size_t i=0; i<Vec.size(); ++i){
      if (Vec[i].label == ModuleLabel) { 
	Index = i;
	break;
      }else if (Vec[i].label != ""){
	Index = i+1;
      }
    }
  }

  void   ReadTriggerInfo::fillTrackTrigInfo(int TrkTrigIndex, const KalSeed*KSeed, trackInfoHist_   &Hist){
    GlobalConstantsHandle<ParticleDataTable> pdt;

    int                nsh = (int)KSeed->hits().size();
    KalSegment const& fseg = KSeed->segments().front();

    double     ndof  = std::max(1.0,nsh - 5.0);
    double     p     = fseg.mom();
    double     chi2d = KSeed->chisquared()/ndof;
    double     pt    = p*std::cos(std::atan(fseg.helix().tanDip()));
    double     d0    = fseg.helix().d0();
    double     clE(-1.);
    if (KSeed->caloCluster()) clE = KSeed->caloCluster()->energyDep();

    Hist._hTrkInfo[TrkTrigIndex][0]->Fill(p);
    Hist._hTrkInfo[TrkTrigIndex][1]->Fill(pt);
    Hist._hTrkInfo[TrkTrigIndex][2]->Fill(nsh);
    Hist._hTrkInfo[TrkTrigIndex][3]->Fill(d0);
    Hist._hTrkInfo[TrkTrigIndex][4]->Fill(chi2d);
    Hist._hTrkInfo[TrkTrigIndex][5]->Fill(clE);
    
    //add the MC info if available
    if (_mcdigis) {
      const mu2e::ComboHit*    hit(0), *hit_0(0);
      hit_0     = &_chcol->at(0);

      int                      loc(-1);
      std::vector<int>         hits_simp_id, hits_simp_index, hits_simp_z;

      for (int j=0; j<nsh; ++j){
	int  hitIndex  = int(KSeed->hits().at(j).index());
	hit            = &_chcol->at(hitIndex);
	loc            = hit - hit_0;
	const mu2e::StepPointMC* step(0);
	const mu2e::StrawDigiMC* sdmc = &_mcdigis->at(loc);
	if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = sdmc->stepPointMC(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = sdmc->stepPointMC(mu2e::StrawEnd::hv ).get();
	}
	
	if (step) {
	  art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	  int sim_id        = simptr->id().asInt();

	  hits_simp_id.push_back   (sim_id); 
	  hits_simp_index.push_back(loc);
	  hits_simp_z.push_back(step->position().z());
	}
      }//end loop over the hits
    
      int     max(0), mostvalueindex(-1);//, mostvalue= hits_simp_id[0];
      float   dz_most(1e4);
      for (int k=0; k<(int)hits_simp_id.size(); ++k){
	int co = (int)std::count(hits_simp_id.begin(), hits_simp_id.end(), hits_simp_id[k]);
	if ( (co>0) &&  (co>max)) {
	  float  dz      = std::fabs(hits_simp_z[k]);
	  if (dz < dz_most){
	    dz_most        = dz;
	    max            = co;
	    //	    mostvalue      = hits_simp_id[k];
	    mostvalueindex = hits_simp_index[k];
	  }
	}
      }
    
      //finally, get the info of the first StrawDigi
      const mu2e::StepPointMC* step(0);
      const mu2e::StrawDigiMC* sdmc = &_mcdigis->at(mostvalueindex);
      if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	step = sdmc->stepPointMC(mu2e::StrawEnd::cal).get();
      }
      else {
	step = sdmc->stepPointMC(mu2e::StrawEnd::hv ).get();
      }
      
      const mu2e::SimParticle * sim (0);

      if (step) {
	art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	int     pdg   = simptr->pdgId();
	art::Ptr<mu2e::SimParticle> mother = simptr;

	while(mother->hasParent()) mother = mother->parent();
	sim = mother.operator ->();
	int      pdgM   = sim->pdgId();
	double   pXMC   = simptr->startMomentum().x();
	double   pYMC   = simptr->startMomentum().y();
	double   pZMC   = simptr->startMomentum().z();
	double   mass(-1.);//  = part->Mass();
	double   energy(-1.);// = sqrt(px*px+py*py+pz*pz+mass*mass);
	mass   = pdt->particle(pdg).ref().mass();
	energy = sqrt(pXMC*pXMC+pYMC*pYMC+pZMC*pZMC+mass*mass);
      
	double   pTMC   = sqrt(pXMC*pXMC + pYMC*pYMC);
	double   pMC    = sqrt(pZMC*pZMC + pTMC*pTMC);
      
	const CLHEP::Hep3Vector* sp = &simptr->startPosition();
	XYZVec origin;
	origin.SetX(sp->x()+3904);
	origin.SetY(sp->y());
	origin.SetZ(sp->z());
	double origin_r = sqrt(origin.x()*origin.x() + origin.y()*origin.y());
	double pz     = sqrt(p*p - pt*pt);

	//now fill the MC histograms
	Hist._hTrkInfo[TrkTrigIndex][10]->Fill(pMC);
	Hist._hTrkInfo[TrkTrigIndex][11]->Fill(pTMC);
	Hist._hTrkInfo[TrkTrigIndex][12]->Fill(pZMC);
	Hist._hTrkInfo[TrkTrigIndex][13]->Fill(p - pMC);
	Hist._hTrkInfo[TrkTrigIndex][14]->Fill(pt - pTMC);
	Hist._hTrkInfo[TrkTrigIndex][15]->Fill(pz - pZMC);
	Hist._hTrkInfo[TrkTrigIndex][16]->Fill(pdg);
	Hist._hTrkInfo[TrkTrigIndex][17]->Fill(origin.z());
	Hist._hTrkInfo[TrkTrigIndex][18]->Fill(origin_r);
	Hist._hTrkInfo[TrkTrigIndex][19]->Fill(pdgM);
	Hist._hTrkInfo[TrkTrigIndex][20]->Fill(energy);
      }
    }
  }

  
  void   ReadTriggerInfo::fillHelixTrigInfo(int HelTrigIndex, const HelixSeed*HSeed, helixInfoHist_  &Hist){
    GlobalConstantsHandle<ParticleDataTable> pdt;

    int        nch       = (int)HSeed->hits().size();
    int        nsh(0);
    for (int i=0; i<nch; ++i) {
      nsh += HSeed->hits().at(i).nStrawHits();
    }
    float      mm2MeV    = (3./10.)*_bz0;

    double     p         = HSeed->helix().momentum()*mm2MeV;
    double     chi2dZPhi = HSeed->helix().chi2dZPhi();
    double     chi2dXY   = HSeed->helix().chi2dXY();
    double     pt        = HSeed->helix().radius()*mm2MeV;
    double     d0        = HSeed->helix().rcent() - HSeed->helix().radius();
    double     clE(-1.);
    double     lambda    = fabs(HSeed->helix().lambda());
    if (HSeed->caloCluster()) clE = HSeed->caloCluster()->energyDep();

    Hist._hHelInfo[HelTrigIndex][0]->Fill(p);
    Hist._hHelInfo[HelTrigIndex][1]->Fill(pt);
    Hist._hHelInfo[HelTrigIndex][2]->Fill(nsh);
    Hist._hHelInfo[HelTrigIndex][3]->Fill(d0);
    Hist._hHelInfo[HelTrigIndex][4]->Fill(chi2dXY);
    Hist._hHelInfo[HelTrigIndex][5]->Fill(chi2dZPhi);
    Hist._hHelInfo[HelTrigIndex][6]->Fill(clE);
    Hist._hHelInfo[HelTrigIndex][7]->Fill(lambda);

     //add the MC info if available
    if (_mcdigis) {
      //      const mu2e::ComboHit*    hit(0);
      std::vector<int>         hits_simp_id, hits_simp_index, hits_simp_z;

      for (int j=0; j<nch; ++j){
	std::vector<StrawDigiIndex> shids;
	HSeed->hits().fillStrawDigiIndices((*_event),j,shids);	
	//	hit            = &HSeed->hits().at(j);

	for (size_t k=0; k<shids.size(); ++k) {
	  const mu2e::StrawDigiMC* sdmc = &_mcdigis->at(shids[k]);
	  art::Ptr<mu2e::StepPointMC>  spmcp;
	  mu2e::TrkMCTools::stepPoint(spmcp,*sdmc);
	  const mu2e::StepPointMC* step = spmcp.get();
	  if (step) {
	    art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	    int sim_id        = simptr->id().asInt();
	    float   dz        = step->position().z();// - trackerZ0;
	    hits_simp_id.push_back   (sim_id); 
	    hits_simp_index.push_back(shids[k]);
	    hits_simp_z.push_back(dz);
	    break;
	  }
	}
      }//end loop over the hits
    
      int     max(0), mostvalueindex(-1);//, mostvalue= hits_simp_id[0];
      float   dz_most(1e4);
      for (int k=0; k<(int)hits_simp_id.size(); ++k){
	int co = (int)std::count(hits_simp_id.begin(), hits_simp_id.end(), hits_simp_id[k]);
	if ( (co>0) &&  (co>max)) {
	  float  dz      = std::fabs(hits_simp_z[k]);
	  if (dz < dz_most){
	    dz_most        = dz;
	    max            = co;
	    //	    mostvalue      = hits_simp_id[k];
	    mostvalueindex = hits_simp_index[k];
	  }
	}
      }
    
      //finally, get the info of the first StrawDigi
      const mu2e::StrawDigiMC* sdmc = &_mcdigis->at(mostvalueindex);
      art::Ptr<mu2e::StepPointMC>  spmcp;
      mu2e::TrkMCTools::stepPoint(spmcp,*sdmc);
      const mu2e::StepPointMC*     step = spmcp.get();  
      const mu2e::SimParticle *    sim (0);

      if (step) {
	art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	int     pdg   = simptr->pdgId();
	art::Ptr<mu2e::SimParticle> mother = simptr;

	while(mother->hasParent()) mother = mother->parent();
	sim = mother.operator ->();
	int      pdgM   = sim->pdgId();
	double   pXMC   = step->momentum().x();
	double   pYMC   = step->momentum().y();
	double   pZMC   = step->momentum().z();
	double   mass(-1.);//  = part->Mass();
	double   energy(-1.);// = sqrt(px*px+py*py+pz*pz+mass*mass);
	mass   = pdt->particle(pdg).ref().mass();
	energy = sqrt(pXMC*pXMC+pYMC*pYMC+pZMC*pZMC+mass*mass);
      
	double   pTMC   = sqrt(pXMC*pXMC + pYMC*pYMC);
	double   pMC    = sqrt(pZMC*pZMC + pTMC*pTMC);
      
	const CLHEP::Hep3Vector* sp = &simptr->startPosition();
	XYZVec origin;
	origin.SetX(sp->x()+3904);
	origin.SetY(sp->y());
	origin.SetZ(sp->z());
	double origin_r = sqrt(origin.x()*origin.x() + origin.y()*origin.y());
	// trackSeed->fOrigin1.SetXYZT(sp->x(),sp->y(),sp->z(),simptr->startGlobalTime());
	double pz     = sqrt(p*p - pt*pt);

	//now fill the MC histograms
	Hist._hHelInfo[HelTrigIndex][10]->Fill(pMC);
	Hist._hHelInfo[HelTrigIndex][11]->Fill(pTMC);
	Hist._hHelInfo[HelTrigIndex][12]->Fill(pZMC);
	Hist._hHelInfo[HelTrigIndex][13]->Fill(p - pMC);
	Hist._hHelInfo[HelTrigIndex][14]->Fill(pt - pTMC);
	Hist._hHelInfo[HelTrigIndex][15]->Fill(pz - pZMC);
	Hist._hHelInfo[HelTrigIndex][16]->Fill(pdg);
	Hist._hHelInfo[HelTrigIndex][17]->Fill(origin.z());
	Hist._hHelInfo[HelTrigIndex][18]->Fill(origin_r);
	Hist._hHelInfo[HelTrigIndex][19]->Fill(pdgM);
	Hist._hHelInfo[HelTrigIndex][20]->Fill(energy);
      }
    }
  }
  //--------------------------------------------------------------------------------
  
  void   ReadTriggerInfo::fillCaloCalibTrigInfo(int ClCalibIndex, const CaloCluster*HCl, caloCalibrationHist_   &Hist){
    int        clsize    = HCl->size();
    double     energy    = HCl->energyDep();

    Hist._hCaloCalibInfo[ClCalibIndex][0]->Fill(energy);
    Hist._hCaloCalibInfo[ClCalibIndex][1]->Fill(clsize);
  }
  //--------------------------------------------------------------------------------

  void   ReadTriggerInfo::fillCaloTrigSeedInfo(int Index, const CaloTrigSeed*HCl, caloTrigSeedHist_      &Hist){
    Hist._hCaloOnlyInfo[Index][0]->Fill(HCl->epeak());
    Hist._hCaloOnlyInfo[Index][1]->Fill(HCl->ring1max());
    Hist._hCaloOnlyInfo[Index][2]->Fill(HCl->ring1max2());
  }
  //--------------------------------------------------------------------------------
  
  void   ReadTriggerInfo::fillOccupancyInfo(int Index         , const StrawDigiCollection*SDCol, const CaloDigiCollection*CDCol, occupancyHist_   &Hist){
    if (_nPOT < 0)          return;
    int   nSD(-1), nCD(-1);
    if (SDCol) nSD = SDCol->size();
    if (CDCol) nCD = CDCol->size();
    
    Hist._hOccInfo  [Index][0]->Fill(_nPOT);
    			    
    Hist._h2DOccInfo[Index][0]->Fill(_nPOT, nSD);
    Hist._h2DOccInfo[Index][1]->Fill(_nPOT, nCD);
  }

  
  
}  

DEFINE_ART_MODULE(mu2e::ReadTriggerInfo);