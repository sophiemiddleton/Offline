//Author: SMiddleton 
//Date: Jan 2020
//Purpose: To make TEVe based event displays in Offline environment
//This is the first stages of this development 
// ... libCore
#include <TApplication.h>
#include <TString.h>
#include <TSystem.h>
#include <TList.h>
#include <TObjArray.h>
#include <Rtypes.h>

// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGString.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
#include<TGLViewer.h>
#include <TGMsgBox.h>
// ... libGeom
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TGeoNode.h>
#include <TGeoPhysicalNode.h>
// ... libEG
#include <TParticle.h>
// ... libRGL
#include <TGLViewer.h>
// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEvePointSet.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveStraightLineSet.h>

//#include <sstream>
#include "fstream"

//TEveEventDisplay Headers:
#include  "TEveEventDisplay/src/dict_classes/NavState.h"
#include  "TEveEventDisplay/src/dict_classes/EvtDisplayUtils.h"

// Mu2e Utilities
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"

//Mu2e Tracker Geom:
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrkDiag/inc/ComboHitInfo.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

//Collections:
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
// Mu2e diagnostics
using namespace std;
using namespace mu2e;

void setRecursiveColorTransp(TGeoVolume *vol, Int_t color, Int_t transp)
  {
     if(color>=0)vol->SetLineColor(color);
     if(transp>=0)vol->SetTransparency(transp);
     Int_t nd = vol->GetNdaughters();
     for (Int_t i=0; i<nd; i++) {
        setRecursiveColorTransp(vol->GetNode(i)->GetVolume(), color, transp);
     }
  }

namespace mu2e 
{
  class TEveEventDisplay : public art::EDAnalyzer {
    public:

struct Config{
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<int> mcdiag{Name("mcdiag"), Comment("set on for MC info"),2};
	      //fhicl::Atom<art::InputTag> modulelabel_{Name("modulelabel"), Comment("")};
	      
	      //fhicl::Atom<art::InputTag> strawStepsTag{Name("strawStepsTag"), Comment("")};
	      fhicl::Atom<art::InputTag> strawDigisTag{Name("strawDigisTag"),Comment("strawDigiTag")};
	      //fhicl::Atom<art::InputTag> strawHitsTag{Name("strawHitsTag"),Comment("straw hit tag")};
	      //fhicl::Atom<std::string> generatorModuleLabel{Name("generatorModuleLabel"), Comment("")};
   	      fhicl::Atom<std::string> g4ModuleLabel{Name("g4ModuleLabel"), Comment("")};
	      //fhicl::Atom<std::string> hitMakerModuleLabel{Name("hitMakerModuleLabel"), Comment("")};
             // fhicl::Atom<std::string> trackerStepPoints{Name("trackerStepPoints"), Comment("")};
     	      fhicl::Atom<double> minEnergyDep{Name("minEnergyDep"), Comment("choose minium energy"), 50};
	       fhicl::Atom<int> minHits{Name("minHits"), Comment(""), 2};
	       fhicl::Atom<bool> doDisplay{Name("doDisplay"), Comment(""), true};
	       fhicl::Atom<bool> clickToAdvance{Name("clickToAdvance"), Comment(""), true}; 
               fhicl::Atom<bool> showEvent{Name("showEvent"), Comment(""),true};     
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;
    explicit TEveEventDisplay(const Parameters& conf);
     virtual ~TEveEventDisplay();
     virtual void beginJob();
     virtual void beginRun(const art::Run& run);
     void analyze(const art::Event& e);
     virtual void endJob();
 private:
     Config _conf;
     int _mcdiag;
     Int_t _evt; 
     //std::string moduleLabel_;
     const StrawDigiCollection* _stcol;
     //art::InputTag strawStepsTag_;
     art::InputTag strawDigisTag_;
     //art::InputTag strawHitsTag_;
    
     //std::string generatorModuleLabel_;
     std::string g4ModuleLabel_;
     //std::string hitMakerModuleLabel_;

     // Name of the tracker StepPoint collection
     //std::string trackerStepPoints_;

     // Cuts used inside SimParticleWithHits:
     //  - drop hits with too little energy deposited.
     //  - drop SimParticles with too few hits.
      double minEnergyDep_;
      size_t minHits_;

      bool doDisplay_;
      bool clickToAdvance_;
      bool showEvent_;
      TApplication* application_;
      TDirectory*   directory_ = nullptr;
      bool            drawGenTracks_;
      bool            drawHits_;
      Double_t        hitMarkerSize_;
      Double_t        trkMaxR_;
      Double_t        trkMaxZ_;
      Double_t        trkMaxStepSize_;
      Double_t        camRotateCenterH_;
      Double_t        camRotateCenterV_;
      Double_t        camDollyDelta_;

      TEveGeoShape* fSimpleGeom;

      TEveViewer *fXYView;
      TEveViewer *fRZView;
      TEveProjectionManager *fXYMgr;
      TEveProjectionManager *fRZMgr;
      TEveScene *fDetXYScene;
      TEveScene *fDetRZScene;
      TEveScene *fEvtXYScene;
      TEveScene *fEvtRZScene;

      TGTextEntry      *fTeRun,*fTeEvt;
      TGLabel          *fTlRun,*fTlEvt;
   
      TEveTrackList *fTrackList;
      TEveElementList *fHitsList;
      
      EvtDisplayUtils *visutil_ = new EvtDisplayUtils();
      
      bool foundEvent = false;
      void MakeNavPanel();
      void InsideDS( TGeoNode * node, bool inDSVac );
      void hideTop(TGeoNode* node);
      void hideNodesByName(TGeoNode* node, const std::string& str,bool onOff) ;
      void hideNodesByMaterial(TGeoNode* node, const std::string& mat, bool onOff);
      void hideBuilding(TGeoNode* node);
      void AddTrack(const art::Event& event, std::vector<int> PDGToAdd);
      void AddHits(const art::Event& event);
      bool FindData(const art::Event& event);
};

TEveEventDisplay::TEveEventDisplay(const Parameters& conf) :
	art::EDAnalyzer(conf),
	_mcdiag	(conf().mcdiag()),
        //strawStepsTag_(conf().strawStepsTag()),
	strawDigisTag_(conf().strawDigisTag()),
        //strawHitsTag_(conf().strawHitsTag()),
        //generatorModuleLabel_(conf().generatorModuleLabel()),
        g4ModuleLabel_(conf().g4ModuleLabel()),
        //hitMakerModuleLabel_(conf().hitMakerModuleLabel()),
        //trackerStepPoints_(conf().trackerStepPoints()),
        minEnergyDep_(conf().minEnergyDep()),
        minHits_(conf().minHits()),
	doDisplay_(conf().doDisplay()),
        clickToAdvance_(conf().clickToAdvance()),
        showEvent_(conf().showEvent()){
		visutil_ = new EvtDisplayUtils();
	}


TEveEventDisplay::~TEveEventDisplay(){}

/*-------Create Control Panel For Event Navigation----""*/
void TEveEventDisplay::MakeNavPanel()
{
  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft); 

  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("EVT NAV");
  frmMain->SetCleanup(kDeepCleanup);
  
  TGHorizontalFrame* navFrame = new TGHorizontalFrame(frmMain);
  TGVerticalFrame* evtidFrame = new TGVerticalFrame(frmMain);
  {
    
    TGPictureButton* b = 0;
    // Create back button and connect to "PrevEvent" in utils. The Buttons are images:
    b = new TGPictureButton(navFrame, gClient->GetPicture("TEveEventDisplay/src/Icons/back.png"));
    //Add button to frame:
    navFrame->AddFrame(b);
    b->Connect("Clicked()", "mu2e::EvtDisplayUtils", visutil_, "PrevEvent()");
    
    //Create forward button and connect to "NextEvent" utils
    b = new TGPictureButton(navFrame, gClient->GetPicture("TEveEventDisplay/src/Icons/forward.png"));
    navFrame->AddFrame(b);
    b->Connect("Clicked()", "mu2e::EvtDisplayUtils", visutil_, "NextEvent()");

    //Create run num text entry widget and connect to "GotoEvent" rcvr in visutils
    TGHorizontalFrame* runFrame = new TGHorizontalFrame(evtidFrame);
    fTlRun = new TGLabel(runFrame,"Run Number");
    fTlRun->SetTextJustify(kTextLeft);
    fTlRun->SetMargins(5,5,5,0);
    runFrame->AddFrame(fTlRun);
    
    fTeRun = new TGTextEntry(runFrame, visutil_->fTbRun = new TGTextBuffer(5), 1);
    visutil_->fTbRun->AddText(0, "1");
    fTeRun->Connect("ReturnPressed()","mu2e::EvtDisplayUtils", visutil_,"GotoEvent()");
    runFrame->AddFrame(fTeRun,new TGLayoutHints(kLHintsExpandX));

    //Create evt num text entry widget and connect to "GotoEvent" rcvr in visutils
    TGHorizontalFrame* evtFrame = new TGHorizontalFrame(evtidFrame);
    fTlEvt = new TGLabel(evtFrame,"Evt Number");
    fTlEvt->SetTextJustify(kTextLeft);
    fTlEvt->SetMargins(5,5,5,0);
    evtFrame->AddFrame(fTlEvt);

    fTeEvt = new TGTextEntry(evtFrame, visutil_->fTbEvt = new TGTextBuffer(5), 1);
    visutil_->fTbEvt->AddText(0, "1");
    fTeEvt->Connect("ReturnPressed()","mu2e::EvtDisplayUtils", visutil_,"GotoEvent()");
    evtFrame->AddFrame(fTeEvt,new TGLayoutHints(kLHintsExpandX));

    // Add horizontal run & event number subframes to vertical evtidFrame
    evtidFrame->AddFrame(runFrame,new TGLayoutHints(kLHintsExpandX));
    evtidFrame->AddFrame(evtFrame,new TGLayoutHints(kLHintsExpandX));

    // Add navFrame and evtidFrame to MainFrame
    frmMain->AddFrame(navFrame);
    TGHorizontal3DLine *separator = new TGHorizontal3DLine(frmMain);
    frmMain->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));
    frmMain->AddFrame(evtidFrame);

    frmMain->MapSubwindows();
    frmMain->Resize();
    frmMain->MapWindow();
    cout<<" finsihing panel making "<<endl;
    browser->StopEmbedding();
    browser->SetTabTitle("Event Nav", 0);
  }
}

void TEveEventDisplay::beginJob(){
  cout<<"Beginning Job ... "<<endl;
/*
  directory_ = gDirectory;
  // Create application environment:
   if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      application_ = new TApplication( "noapplication", &tmp_argc, tmp_argv );
   }
  // Initialize global Eve application manager (return gEve)
  TEveManager::Create();
  
  // Create detector and event scenes for ortho views
  fDetXYScene = gEve->SpawnNewScene("Det XY Scene", "");
  fDetRZScene = gEve->SpawnNewScene("Det RZ Scene", "");
  fEvtXYScene = gEve->SpawnNewScene("Evt XY Scene", "");
  fEvtRZScene = gEve->SpawnNewScene("Evt RZ Scene", "");

  // Create XY/RZ projection mgrs, draw projected axes, & add them to scenes
  fXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
  TEveProjectionAxes* axes_xy = new TEveProjectionAxes(fXYMgr);
  fDetXYScene->AddElement(axes_xy);

  fRZMgr = new TEveProjectionManager(TEveProjection::kPT_RhoZ);
  TEveProjectionAxes* axes_rz = new TEveProjectionAxes(fRZMgr);
  fDetRZScene->AddElement(axes_rz);

  // Create side-by-side ortho XY & RZ views in new tab & add det/evt scenes
  TEveWindowSlot *slot = 0;
  TEveWindowPack *pack = 0;

  slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
  pack = slot->MakePack();
  pack->SetElementName("Ortho Views");
  pack->SetHorizontal();
  pack->SetShowTitleBar(kFALSE);

  pack->NewSlot()->MakeCurrent();
  fXYView = gEve->SpawnNewViewer("XY View", "");
  fXYView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  fXYView->AddScene(fDetXYScene);
  fXYView->AddScene(fEvtXYScene);

  pack->NewSlot()->MakeCurrent();
  fRZView = gEve->SpawnNewViewer("RZ View", "");
  fRZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  fRZView->AddScene(fDetRZScene);
  fRZView->AddScene(fEvtRZScene);

  gEve->GetBrowser()->GetTabRight()->SetTab(0);

  // Create navigation panel
  MakeNavPanel();
  cout<<"adding event "<<endl;
  // Add new Eve event into the "Event" scene and make it the current event
  // (Subsequent elements added using "AddElements" will be added to this event)
  gEve->AddEvent(new TEveEventManager("Event", "Toy Detector Event"));

  TGLViewer *glv = gEve->GetDefaultGLViewer();
  glv->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
  glv->CurrentCamera().RotateRad(camRotateCenterH_,camRotateCenterV_);
  glv->CurrentCamera().Dolly(camDollyDelta_,kFALSE,kFALSE);
*/
}


void TEveEventDisplay::beginRun(const art::Run& run){
  /*
  if(gGeoManager){
    gGeoManager->GetListOfNodes()->Delete();
    gGeoManager->GetListOfVolumes()->Delete();
    gGeoManager->GetListOfShapes()->Delete();
  }
  gEve->GetGlobalScene()->DestroyElements();
  fDetXYScene->DestroyElements();
  fDetRZScene->DestroyElements();
  */
  cout<<"importing GDML "<<endl;
  // Import the GDML of entire Mu2e Geometry
  TGeoManager* geom = new TGeoManager("geom","Geom");
  geom->TGeoManager::SetNavigatorsLock(kFALSE);
  geom->TGeoManager::ClearNavigators();

  geom = geom->TGeoManager::Import("TEveEventDisplay/src/mu2e.gdml");
  cout<<"Getting top volume "<<endl;
  //Get Top Volume
  TGeoVolume* topvol = geom->GetTopVolume();
  
  //Set Top Volume for gGeoManager:
  cout<<"setting top "<<endl;
  gGeoManager->SetTopVolume(topvol);
  gGeoManager->SetTopVisible(kTRUE);
  int nn = gGeoManager->GetNNodes();
  printf("nodes in geom = %d\n",nn);
  //Get Top Node:
  TGeoNode* topnode = gGeoManager->GetTopNode();
  TEveGeoTopNode* etopnode = new TEveGeoTopNode(gGeoManager, topnode);
  etopnode->SetVisLevel(4);
  etopnode->GetNode()->GetVolume()->SetVisibility(kFALSE);
  //Set colours to allow transparency:
  setRecursiveColorTransp(etopnode->GetNode()->GetVolume(), kWhite-10,70);
  //This is a basic function to allow us to just see tracker and calo, it needs fixing:
  cout<<"extracting"<<endl;
  //InsideDS( topnode, false );
  //hideBuilding(topnode);
  //hideTop(topnode);
 
  //Add static detector geometry to global scene
  gEve->AddGlobalElement(etopnode);
  
}

void TEveEventDisplay::InsideDS( TGeoNode * node, bool inDSVac ){
  cout<<"In side DS"<<endl;
  std::string _name = (node->GetVolume()->GetName());
  if ( node->GetMotherVolume() ) {
    std::string motherName(node->GetMotherVolume()->GetName());
    if ( motherName == "DS2Vacuum" || motherName == "DS3Vacuum" ){
      inDSVac = true;
    }
  }

  if ( inDSVac && _name.find("Virtual") != 0 ) {
    node->SetVisibility(kTRUE);
  } else{
    node->SetVisibility(kFALSE);
  }

  // Descend into each daughter TGeoNode.
  int ndau = node->GetNdaughters();
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * dau = node->GetDaughter(i);
    InsideDS( dau, inDSVac );
  }

}

void TEveEventDisplay::analyze(const art::Event& event){
 
  _evt = event.id().event();
  if(showEvent_ ){
  	FindData(event);
  }
  std::ostringstream sstr;
  sstr << event.id().run();
  visutil_->fTbRun->Clear();
  visutil_->fTbRun->AddText(0,sstr.str().c_str());
  gClient->NeedRedraw(fTeRun);

  sstr.str("");
  sstr << event.id().event();
  visutil_->fTbEvt->Clear();
  visutil_->fTbEvt->AddText(0,sstr.str().c_str());
  gClient->NeedRedraw(fTeEvt);
  
  // Delete visualization structures associated with previous event
  gEve->GetViewers()->DeleteAnnotations();
  gEve->GetCurrentEvent()->DestroyElements();

 /*
  AddHits(event);
  AddTracks(event, [11,13]);
  // Import event into ortho views and apply projections
  TEveElement* currevt = gEve->GetCurrentEvent();

  fEvtXYScene->DestroyElements();
  fXYMgr->ImportElements(currevt, fEvtXYScene);

  fEvtRZScene->DestroyElements();
  fRZMgr->ImportElements(currevt, fEvtRZScene);
*/
 
} 

void TEveEventDisplay::hideTop(TGeoNode* node) {
cout<<"hitde top "<<endl;
  TString name = node->GetName();
  if(name.Index("Shield")>0) {
    std::cout << name << " " <<  name.Index("mBox_") << std::endl;
  }
  bool test = false;

  // from gg1
  if(name.Index("mBox_45_")>=0) test = true;
  if(name.Index("mBox_46_")>=0) test = true;
  if(name.Index("mBox_47_")>=0) test = true;
  if(name.Index("mBox_48_")>=0) test = true;
  if(name.Index("mBox_49_")>=0) test = true;
  if(name.Index("mBox_74_")>=0) test = true;

  if(test) {
    std::cout << "turning off " << name << std::endl;
    node->SetVisibility(false);
  }


  /*
  std::string str("ExtShield");
  if ( name.find(str) != std::string::npos ){
    TGeoVolume* gv = node->GetVolume(); 
    printf("vol %ul\n",gv);
    TGeoShape* gs = gv->GetShape();
    printf("ExtShield shape %s\n",gs->IsA()->GetName());
    Double_t dx,dy,dz;
  }
  */
  
  // Descend recursively into each daughter TGeoNode.
  int ndau = node->GetNdaughters();
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * dau = node->GetDaughter(i);
    hideTop( dau );
  }

}

void TEveEventDisplay::hideNodesByName(TGeoNode* node, const std::string& str,
				     bool onOff) {

  std::string name(node->GetName());
  if ( name.find(str) != std::string::npos ){
    node->SetVisibility(onOff);
    std::cout <<"hiding "<< name << std::endl;
  }

  // Descend recursively into each daughter TGeoNode.
  int ndau = node->GetNdaughters();
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * dau = node->GetDaughter(i);
    hideNodesByName( dau, str, onOff);
  }

}

void TEveEventDisplay::hideNodesByMaterial(TGeoNode* node, 
					 const std::string& mat, bool onOff) {

  std::string material(node->GetVolume()->GetMaterial()->GetName());
  if ( material.find(mat) != std::string::npos ) node->SetVisibility(onOff);

  // Descend recursively into each daughter TGeoNode.
  int ndau = node->GetNdaughters();
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * dau = node->GetDaughter(i);
    hideNodesByMaterial( dau, mat, onOff);
  }

}

void TEveEventDisplay::hideBuilding(TGeoNode* node) {

  static std::vector <std::string> substrings  { "Ceiling",
      "backfill", "dirt", "concrete", "VirtualDetector",
      "pipeType","CRSAluminium","CRV","CRS", "ExtShield", "PSShield"};
  for(auto& i: substrings) hideNodesByName(node,i,kFALSE);

  static std::vector <std::string> materials { "MBOverburden", "CONCRETE"};
  for(auto& i: materials) hideNodesByMaterial(node,i,kFALSE);


}

void TEveEventDisplay::AddTrack(const art::Event& event, std::vector<int> PDGToAdd){
    /*if (drawGenTracks_) {
    auto gens = event.getValidHandle<GenParticleCollection>(gensTag_);

    if (fTrackList == 0) {
      fTrackList = new TEveTrackList("Tracks");
      fTrackList->SetLineWidth(4);
      fTrackList->IncDenyDestroy(); 
    }
    else {
      fTrackList->DestroyElements();         
    }

    TEveTrackPropagator* trkProp = fTrackList->GetPropagator();
    trkProp->SetMagField(-geom_->bz()*1000.);
    trkProp->SetMaxR(trkMaxR_);
    trkProp->SetMaxZ(trkMaxZ_);
    trkProp->SetMaxStep(trkMaxStepSize_);

    int mcindex=-1;
    for ( auto const& gen: *gens){
      mcindex++;
      if ( gen.hasChildren() ) continue;
      TParticle mcpart;
      mcpart.SetMomentum(gen.momentum().px(),gen.momentum().py(),gen.momentum().pz(),gen.momentum().e());
      mcpart.SetProductionVertex(gen.position().x()*0.1,gen.position().y()*0.1,gen.position().z()*0.1,0.);
      mcpart.SetPdgCode(gen.pdgId());
      TEveTrack* track = new TEveTrack(&mcpart,mcindex,trkProp);
      track->SetIndex(0);
      track->SetStdTitle();
      track->SetAttLineAttMarker(fTrackList);
      if ( gen.pdgId() == 11 ){
        track->SetMainColor(kGreen);
      } else {
        track->SetMainColor(kViolet+1);
      }
      fTrackList->AddElement(track);
    }
    fTrackList->MakeTracks();
    gEve->AddElement(fTrackList);
  }*/
}

void TEveEventDisplay::AddHits(const art::Event& event){
    
  /* // Draw the detector hits
   if (fHitsList == 0) {
      fHitsList = new TEveElementList("Hits");
      fHitsList->IncDenyDestroy();     
    }
    else {
      fHitsList->DestroyElements();  
    }

    TEveElementList* ElectronHitsList  = new TEveElementList("Electron Hits");
    TEveElementList* BkgHitsList = new TEveElementList("Bkg Hits");
    int ie = 0;
    int ibkg = 0;
    drawHit("e",kGreen,hitMarkerSize_,ie++,hit,ElectronHitsList);
    drawHit("Bkg",kViolet+1,hitMarkerSize_,ibkg++,hit,BkgHitsList);
        
    fHitsList->AddElement(ElectronHitsList);  
    fHitsList->AddElement(BkgHitsList);  
    gEve->AddElement(fHitsList);
  */
}


bool TEveEventDisplay::FindData(const art::Event& evt){
	_stcol = 0; 
        cout<<"finding data "<<endl;
	auto chH = evt.getValidHandle<mu2e::StrawDigiCollection>(strawDigisTag_);
	_stcol = chH.product();
	foundEvent = true;
	return _stcol != 0;
       }

void TEveEventDisplay::endJob(){
	if(!foundEvent){
		char msg[300];
		sprintf(msg, "Reached end of file but #%i has not been found", true);
	        new TGMsgBox(gClient->GetRoot(), gClient->GetRoot(), "Event Not Found", msg, kMBIconExclamation,kMBOk);
	}

}  
	
}

using mu2e::TEveEventDisplay;
DEFINE_ART_MODULE(TEveEventDisplay);
