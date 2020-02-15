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
#include <TPolyLine3D.h>
// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGString.h>
#include <TGLabel.h>
#include <TGIcon.h>
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
#include  "TEveEventDisplay/src/dict_classes/Geom_Interface.h"

// Mu2e Utilities
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/WorldG4Maker.hh"
#include "GeometryService/inc/Mu2eCoordTransform.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
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
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
// Mu2e diagnostics
using namespace std;
using namespace mu2e;

void DrawHit(const std::string &pstr, Int_t mColor, Int_t mSize, Int_t n, CLHEP::Hep3Vector HitPos, TEveElementList *list, Geom_Interface *g) //FIXME - should be combo hit here , t should take in the CLHEPvec
  {
	
	g->GetTrackerCenter();
        CLHEP::Hep3Vector pointInMu2e = g->PointToTracker(HitPos);
	std::string hstr=" hit %d";
	std::string dstr=" hit# %d\nLayer: %d";
	std::string strlst=pstr+hstr;
	std::string strlab=pstr+dstr;

	TEvePointSet* h = new TEvePointSet(Form(strlst.c_str(),n));
	h->SetTitle(Form(strlab.c_str(),n,hstr));

	cout<<"in mu2e : "<<n<<" "<<pointInMu2e.x()/10<<" "<<pointInMu2e.y()/10<<" "<<pointInMu2e.z()/10<<endl;
	h->SetMarkerColor(kRed);
	h->SetMarkerSize(mSize);

	h->SetNextPoint(pointInMu2e.x()/10, pointInMu2e.y()/10, pointInMu2e.z()/10); //as GDML is in cm (I think...and hope)
	h->SetMarkerColor(mColor);
	h->SetMarkerSize(mSize);
	list->AddElement(h);
  }


void DrawCluster(const std::string &pstr, Int_t mColor, Int_t mSize, Int_t n, CLHEP::Hep3Vector ClusterPos, int nDisk, TEveElementList *list, Geom_Interface *g) //FIXME - should be combo hit here , t should take in the CLHEPvec
  {
	
	g->GetCaloCenter(nDisk);
        CLHEP::Hep3Vector pointInMu2e = g->PointToCalo(ClusterPos, nDisk);
	std::string hstr=" cluster %d";
	std::string dstr=" cluster# %d\nLayer: %d";
	std::string strlst=pstr+hstr;
	std::string strlab=pstr+dstr;

	TEvePointSet* h = new TEvePointSet(Form(strlst.c_str(),n));
	h->SetTitle(Form(strlab.c_str(),n,hstr));

	cout<<"in mu2e : "<<n<<" "<<pointInMu2e.x()/10<<" "<<pointInMu2e.y()/10<<" "<<pointInMu2e.z()/10<<endl;
	h->SetMarkerColor(kRed);
	h->SetMarkerSize(mSize);

	h->SetNextPoint(pointInMu2e.x()/10, pointInMu2e.y()/10, pointInMu2e.z()/10); //as GDML is in cm (I think...and hope)
	h->SetMarkerColor(mColor);
	h->SetMarkerSize(mSize);
	list->AddElement(h);
  }


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
			fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("for info"),0};
			fhicl::Atom<art::InputTag>chTag{Name("ComboHitCollection"),Comment("chTag")};
			fhicl::Atom<art::InputTag>gensTag{Name("GenParticleCollection"),Comment("gensTag")};
			fhicl::Atom<art::InputTag>strawdigiTag{Name("StrawDigiCollection"),Comment("strawdigiTag")};
			fhicl::Atom<art::InputTag>crvdigiTag{Name("CrvDigiCollection"),Comment("crvTag")};
			fhicl::Atom<art::InputTag>cosmicTag{Name("CosmicTrackSeedCollection"),Comment("cosmicTag")};
			fhicl::Atom<art::InputTag>cluTag{Name("CaloClusterCollection"),Comment("cluTag")};
			fhicl::Atom<std::string> g4ModuleLabel{Name("g4ModuleLabel"), Comment("")};
			fhicl::Atom<double> minEnergyDep{Name("minEnergyDep"), Comment("choose minium energy"), 50};
			fhicl::Atom<int> minHits{Name("minHits"), Comment(""), 2};
			fhicl::Atom<bool> doDisplay{Name("doDisplay"), Comment(""), true};
			fhicl::Atom<bool> clickToAdvance{Name("clickToAdvance"), Comment(""), true}; 
			fhicl::Atom<bool> showEvent{Name("showEvent"), Comment(""),true};     
			fhicl::Atom<bool> addHits{Name("addHits"), Comment("set to add the hits"),false};
			fhicl::Atom<bool> addTracks{Name("addTracks"), Comment("set to add tracks"),false};
			fhicl::Atom<bool> addClusters{Name("addClusters"), Comment("set to add calo lusters"),false};
			fhicl::Atom<bool> addCrvHits{Name("addCrvHits"), Comment("set to add crv hits"),false};	
			fhicl::Atom<bool> addCosmicSeedFit{Name("addCosmicSeedFit"), Comment("for fitted cosmic track"), false};
			fhicl::Atom<bool> isCosmic{Name("isCosmic"), Comment("flag for cosmic track v helix track"), false};	
	    };
		typedef art::EDAnalyzer::Table<Config> Parameters;
		explicit TEveEventDisplay(const Parameters& conf);
		virtual ~TEveEventDisplay();
		virtual void beginJob() override;
		virtual void beginRun(const art::Run& run) override;
		virtual void analyze(const art::Event& e);
		virtual void endJob() override;
	private:
		     Config _conf;
		     int _diagLevel;
		     Int_t _evt; 
		     //std::string moduleLabel_;
		     const StrawDigiCollection* _stcol;
		     const ComboHitCollection* _chcol;
		     const StrawDigiCollection* _strawdigicol;
		     const CrvDigiCollection* _crvdigicol;
		     const CosmicTrackSeedCollection* _cosmiccol;
		     const GenParticleCollection* _gencol;
		     const CaloClusterCollection* _clustercol;
		     art::InputTag chTag_;
		     art::InputTag gensTag_;
		     art::InputTag strawdigiTag_;
		     art::InputTag crvdigiTag_;
		     art::InputTag cosmicTag_;
		     art::InputTag cluTag_;
		     std::string g4ModuleLabel_;
		     //std::string hitMakerModuleLabel_;

		     // Name of the tracker StepPoint collection
		     //std::string trackerStepPoints_;


		      double minEnergyDep_;
		      size_t minHits_;
		      bool isFirstEvent = true;
		      bool doDisplay_;
		      bool clickToAdvance_;
		      bool showEvent_;
		      TApplication* application_;
		      TDirectory*   directory_ = nullptr;
		      
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
		   
		     // GeomHandle<mu2e::BFieldManager> bfmgr;
		     
		      TEveTrackList *fTrackList;
		      TEveElementList *fHitsList;
		      TEveElementList *fClusters;
		      
		      bool addHits_, addTracks_, addClusters_, addCrvHits_, addCosmicSeedFit_, isCosmic_;
		     
		      EvtDisplayUtils *visutil_ = new EvtDisplayUtils();
		      Geom_Interface *gdml_geom	=new Geom_Interface(); 
		      //Particle_Interface *particle_info = new Particle_Interface();

		      TGeoManager* geom = new TGeoManager("geom","Geom");
		      std::vector<CLHEP::Hep3Vector> GDMLt;

		  
		      bool foundEvent = false;
		      void MakeNavPanel();
		      void Heirarchy(TGeoNode *node, std::vector<CLHEP::Hep3Vector>& t);
		      void InsideDS( TGeoNode * node, bool inDSVac );
		      void hideTop(TGeoNode* node);
		      void hideNodesByName(TGeoNode* node, const std::string& str,bool onOff) ;
		      void hideNodesByMaterial(TGeoNode* node, const std::string& mat, bool onOff);
		      void hideBuilding(TGeoNode* node);
		      void AddCosmicTrack(const art::Event& event);
		      void AddHelicalTrack(const art::Event& event, mu2e::BFieldManager const& fm);
		      void AddHits(const art::Event& event);
		      void AddCaloCluster(const art::Event& event);
		      void AddCrvHits(const art::Event& event);
		      bool FindData(const art::Event& event);
	};

TEveEventDisplay::TEveEventDisplay(const Parameters& conf) :
	art::EDAnalyzer(conf),
	_diagLevel(conf().diagLevel()),
	chTag_(conf().chTag()),
	gensTag_(conf().gensTag()),
	strawdigiTag_(conf().strawdigiTag()),
	crvdigiTag_(conf().crvdigiTag()),
	g4ModuleLabel_(conf().g4ModuleLabel()),
	minEnergyDep_(conf().minEnergyDep()),
	minHits_(conf().minHits()),
	doDisplay_(conf().doDisplay()),
	clickToAdvance_(conf().clickToAdvance()),
	showEvent_(conf().showEvent()),
	addHits_(conf().addHits()),
	addTracks_(conf().addTracks()),
	addClusters_(conf().addClusters()),
	addCrvHits_(conf().addCrvHits()),
	addCosmicSeedFit_(conf().addCosmicSeedFit()),
	isCosmic_(conf().isCosmic()){
		visutil_ = new EvtDisplayUtils();
	}


TEveEventDisplay::~TEveEventDisplay(){}

/*-------Create Control Panel For Event Navigation----""*/
void TEveEventDisplay::MakeNavPanel()
{
	TEveBrowser* browser = gEve->GetBrowser();
	browser->StartEmbedding(TRootBrowser::kLeft); // insert nav frame as new tab in left pane

	TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
	frmMain->SetWindowName("EVT NAV");
	frmMain->SetCleanup(kDeepCleanup);

	TGHorizontalFrame* navFrame = new TGHorizontalFrame(frmMain);
	TGVerticalFrame* evtidFrame = new TGVerticalFrame(frmMain);
	{
		TString icondir(TString::Format("%s/icons/", gSystem->Getenv("ROOTSYS")) );
		TGPictureButton* b = 0;

		// ... Create back button and connect to "PrevEvent" rcvr in visutils
		b = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoBack.gif"));
		navFrame->AddFrame(b);
		b->Connect("Clicked()", "mu2e::EvtDisplayUtils", visutil_, "PrevEvent()");

		// ... Create forward button and connect to "NextEvent" rcvr in visutils
		b = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoForward.gif"));
		navFrame->AddFrame(b);
		b->Connect("Clicked()", "mu2e::EvtDisplayUtils", visutil_, "NextEvent()");

		// ... Create run num text entry widget and connect to "GotoEvent" rcvr in visutils
		TGHorizontalFrame* runoFrame = new TGHorizontalFrame(evtidFrame);
		fTlRun = new TGLabel(runoFrame,"Run Number");
		fTlRun->SetTextJustify(kTextLeft);
		fTlRun->SetMargins(5,5,5,0);
		runoFrame->AddFrame(fTlRun);
    
		fTeRun = new TGTextEntry(runoFrame, visutil_->fTbRun = new TGTextBuffer(5), 1);
		visutil_->fTbRun->AddText(0, "1");
		fTeRun->Connect("ReturnPressed()","mu2e::EvtDisplayUtils", visutil_,"GotoEvent()");
		runoFrame->AddFrame(fTeRun,new TGLayoutHints(kLHintsExpandX));

		// ... Create evt num text entry widget and connect to "GotoEvent" rcvr in visutils
		TGHorizontalFrame* evnoFrame = new TGHorizontalFrame(evtidFrame);
		fTlEvt = new TGLabel(evnoFrame,"Evt Number");
		fTlEvt->SetTextJustify(kTextLeft);
		fTlEvt->SetMargins(5,5,5,0);
		evnoFrame->AddFrame(fTlEvt);

		//Add logo
		std::string logoFile = "TEveEventDisplay/src/Icons/mu2e_logo_oval.png";
		const TGPicture *logo = gClient->GetPicture(logoFile.c_str());
		TGIcon *icon = new TGIcon(navFrame,logo,50,50);
		navFrame->AddFrame(icon,new TGLayoutHints(kLHintsLeft,20,0,0,0));

		fTeEvt = new TGTextEntry(evnoFrame, visutil_->fTbEvt = new TGTextBuffer(5), 1);
		visutil_->fTbEvt->AddText(0, "1");
		fTeEvt->Connect("ReturnPressed()","mu2e::EvtDisplayUtils", visutil_,"GotoEvent()");
		evnoFrame->AddFrame(fTeEvt,new TGLayoutHints(kLHintsExpandX));

		// ... Add horizontal run & event number subframes to vertical evtidFrame
		evtidFrame->AddFrame(runoFrame,new TGLayoutHints(kLHintsExpandX));
		evtidFrame->AddFrame(evnoFrame,new TGLayoutHints(kLHintsExpandX));

		// ... Add navFrame and evtidFrame to MainFrame
		frmMain->AddFrame(navFrame);
		TGHorizontal3DLine *separator = new TGHorizontal3DLine(frmMain);
		frmMain->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));
		frmMain->AddFrame(evtidFrame);

		frmMain->MapSubwindows();
		frmMain->Resize();
		frmMain->MapWindow();

		browser->StopEmbedding();
		browser->SetTabTitle("Event Nav", 0);
 	 }
}

void TEveEventDisplay::beginJob(){
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

	MakeNavPanel();

	gEve->AddEvent(new TEveEventManager("Event", "Empty Event"));

	TGLViewer *glv = gEve->GetDefaultGLViewer();
	glv->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
	glv->CurrentCamera().RotateRad(camRotateCenterH_,camRotateCenterV_);
	glv->CurrentCamera().Dolly(camDollyDelta_,kFALSE,kFALSE);
}


void TEveEventDisplay::beginRun(const art::Run& run){

	if(gGeoManager){
		gGeoManager->GetListOfNodes()->Delete();
		gGeoManager->GetListOfVolumes()->Delete();
		gGeoManager->GetListOfShapes()->Delete();
	}
	gEve->GetGlobalScene()->DestroyElements();
	fDetXYScene->DestroyElements();
	fDetRZScene->DestroyElements();

	// Import the GDML of entire Mu2e Geometry
	geom = gdml_geom->Geom_Interface::getGeom("TEveEventDisplay/src/fix.gdml");
	gdml_geom->GetTrackerCenter();

	//Get Top Volume
	TGeoVolume* topvol = geom->GetTopVolume();

	//Set Top Volume for gGeoManager:
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


	hideBuilding(topnode);
	hideTop(topnode);
	InsideDS( topnode, false );

	Heirarchy(topnode, GDMLt);
	for(auto const& transformation : GDMLt){
		cout<<"Transformation "<<transformation.x()<<" "<<transformation.y()<<" "<<transformation.z()<<endl;
  	}
  	//Add static detector geometry to global scene
  	gEve->AddGlobalElement(etopnode);
  
}


void TEveEventDisplay::Heirarchy( TGeoNode * node, std::vector<CLHEP::Hep3Vector> &TransformList ){
  std::string _name = (node->GetVolume()->GetName());
  if( _name == "HallAir") {
	cout<<"HallAir Origin IS "<<node->GetMotherVolume()->GetName();
        TGeoVolume *vol = node->GetVolume();
	TGeoBBox *shape = (TGeoBBox*)vol->GetShape();
	Double_t master[3];
	const Double_t *local = shape->GetOrigin();
	if(shape!=NULL){
		gGeoManager->LocalToMaster(local,master);
		CLHEP::Hep3Vector hallToworld(master[0], master[1], master[2]);
		TransformList.push_back(hallToworld);
        }
  }
  if( _name == "DS3Vacuum") {
	cout<<"DS3 Origin IS "<<node->GetMotherVolume()->GetName();
        TGeoVolume *vol = node->GetVolume();
	TGeoBBox *shape = (TGeoBBox*)vol->GetShape();
	Double_t master[3];
	const Double_t *local = shape->GetOrigin();
	if(shape!=NULL){
		gGeoManager->LocalToMaster(local,master);
		CLHEP::Hep3Vector DSTohall(master[0], master[1], master[2]);
		TransformList.push_back(DSTohall);
        }
 }
  if( _name == "TrackerMother") {
	cout<<"Tracker Origin IS "<<node->GetMotherVolume()->GetName();
        TGeoVolume *vol = node->GetVolume();
	TGeoBBox *shape = (TGeoBBox*)vol->GetShape();
	Double_t master[3];
	const Double_t *local = shape->GetOrigin();
	if(shape!=NULL){
		gGeoManager->LocalToMaster(local,master);
		CLHEP::Hep3Vector TrackerToDS(master[0], master[1], master[2]);
		TransformList.push_back(TrackerToDS);
        }

  }

  // Descend into each daughter TGeoNode.
  int ndau = node->GetNdaughters();
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * dau = node->GetDaughter(i);
    Heirarchy(dau, TransformList);
  }
  
}

void TEveEventDisplay::InsideDS( TGeoNode * node, bool inDSVac ){
	std::string _name = (node->GetVolume()->GetName());
	if ( node->GetMotherVolume() ) {
		std::string motherName(node->GetMotherVolume()->GetName());
		if ( motherName == "DS2Vacuum" || motherName == "DS3Vacuum" ){
		inDSVac = true;
		}
	}
	if ( inDSVac && _name.find("VirtualDetector_TT_Mid") != 0 ) {
		node->SetVisibility(kTRUE);
	} else{
		node->SetVisibility(kFALSE);
	}
	int ndau = node->GetNdaughters();
	for ( int i=0; i<ndau; ++i ){
		TGeoNode * dau = node->GetDaughter(i);
		InsideDS( dau, inDSVac );
  	}

}

void TEveEventDisplay::analyze(const art::Event& event){
	if(_diagLevel > 0){
		cout<<"Analyzing Event :"<<event.id()<<endl;
	}
	GeomHandle<mu2e::BFieldManager> bfmgr;
	_evt = event.id().event();
	if(showEvent_ ){
		if(addHits_) AddHits(event);
		if(addTracks_ ) AddHelicalTrack(event, *bfmgr);
		//if(addCosmicSeedFit_ and isCosmic_) AddCosmicTrack(event);
		//if(AddCrvHits_)
		
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

	// Delete visualization structures associated with previous event
	if(!isFirstEvent){
		gEve->GetViewers()->DeleteAnnotations();
		gEve->GetCurrentEvent()->DestroyElements();
	}
	// Import event into ortho views and apply projections
	TEveElement* currevt = gEve->GetCurrentEvent();

	fEvtXYScene->DestroyElements();
	fXYMgr->ImportElements(currevt, fEvtXYScene);

	fEvtRZScene->DestroyElements();
	fRZMgr->ImportElements(currevt, fEvtRZScene);

	geom->Draw("ogl");
	gPad->WaitPrimitive();
	isFirstEvent = false;
} 

void TEveEventDisplay::hideTop(TGeoNode* node) {
  	TString name = node->GetName();
  	if(_diagLevel > 0 and name.Index("Shield")>0) {
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
		if(_diagLevel > 0) std::cout <<"hiding "<< name << std::endl;
	}
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

void TEveEventDisplay::AddCosmicTrack(const art::Event& event){
	auto cosH = event.getValidHandle<mu2e::CosmicTrackSeedCollection>(cosmicTag_);
	_cosmiccol = cosH.product();
	if(!_cosmiccol->empty()){
		TEveStraightLineSet *CosmicTrackList = new TEveStraightLineSet();
		for(size_t ist = 0;ist < _cosmiccol->size(); ++ist){
			CosmicTrackSeed sts =(*_cosmiccol)[ist];
			CosmicTrack st = sts._track;
			
			CosmicTrackList->SetLineColor(kGreen);
			Float_t tz1 = -150;
			Float_t tz2 = 150;
			Float_t tx1 = st.InitParams.A0  + st.InitParams.A1*tz1;
			Float_t tx2 = st.InitParams.A0  + st.InitParams.A1*tz2;
			Float_t ty1 = st.InitParams.B0  + st.InitParams.B1*tz1;
			Float_t ty2 = st.InitParams.B0  + st.InitParams.B1*tz2; 	
			CosmicTrackList->AddLine(tx1, ty1, tz1, tx2, ty2, tz2);
		
			cout<<st.InitParams.A0<<"track "<<endl;
			gEve->AddElement(CosmicTrackList);
		    	gEve->Redraw3D();
		}
	}
}


void TEveEventDisplay::AddHelicalTrack(const art::Event& event, mu2e::BFieldManager const& bf){
	auto genH = event.getValidHandle<GenParticleCollection>(gensTag_);
	_gencol = genH.product();
	if (fTrackList == 0) {
		fTrackList = new TEveTrackList("Tracks");
		fTrackList->SetLineWidth(4);
		fTrackList->IncDenyDestroy(); 
	}
	else {
		fTrackList->DestroyElements();         
	}

	int mcindex=-1;
	for ( auto const& gen: *_gencol){
		TEveTrackPropagator* trkProp = fTrackList->GetPropagator();
		//if(!isCosmic_) CLHEP::Hep3Vector field = bf.getBField(CLHEP::Hep3Vector(gen.position().x(), gen.position().y(), gen.position().z()));
		//if(isCosmic_) CLHEP::Hep3Vector field(CLHEP::Hep3Vector(0,0,0)); 
		trkProp->SetMagField(0);//-1*1000.);
		trkProp->SetMaxR(trkMaxR_);
		trkProp->SetMaxZ(trkMaxZ_);
		trkProp->SetMaxStep(trkMaxStepSize_);
		mcindex++;
		      //if ( gen.hasChildren() ) continue;
		TParticle mcpart;
		mcpart.SetMomentum(gen.momentum().px(),
			gen.momentum().py(),	
			gen.momentum().pz(),
			gen.momentum().e());
		 	mcpart.SetProductionVertex(gen.position().x(), 		
						gen.position().y(), gen.position().z(),0.);
		      mcpart.SetPdgCode(gen.pdgId());
		      TEveTrack* track = new TEveTrack(&mcpart,mcindex,trkProp);
		      track->SetIndex(0);
		      track->SetStdTitle();
		      track->SetAttLineAttMarker(fTrackList);
			
		      if ( abs(gen.pdgId()) == 11 ){
		        track->SetMainColor(kRed);
		      }  if (abs(gen.pdgId()) == 13 ){
			track->SetMainColor(kGreen);
		      } else {
			track->SetMainColor(kBlue);
		      }
		      fTrackList->AddElement(track);
		    }
		    fTrackList->MakeTracks();
	            fTrackList->SetLineWidth(10);
		    gEve->AddElement(fTrackList);
		    gEve->Redraw3D();
}



void TEveEventDisplay::AddHits(const art::Event& event){
   
	auto chH = event.getValidHandle<mu2e::ComboHitCollection>(chTag_);
	_chcol = chH.product(); //this should be any collection eventually

	if (fHitsList == 0) {
		fHitsList = new TEveElementList("Hits");
		fHitsList->IncDenyDestroy();     
	}
	else {
		fHitsList->DestroyElements();  
	}
	TEveElementList* HitsList  = new TEveElementList("Combo Hits");
	if(!_chcol->empty()){
	    cout<<"Number of Hits "<<_chcol->size()<<endl;
	    for(unsigned ih = 0 ; ih < _chcol->size() ; ih ++){
		    ComboHit const& hit = (*_chcol)[ih];
		    CLHEP::Hep3Vector HitPos(hit.pos().x(), hit.pos().y(), hit.pos().z());
		    if(ih%100 ==0 and !isCosmic_) DrawHit("ComboHits",kRed, 2, ih, HitPos, HitsList, gdml_geom);
		    else DrawHit("ComboHits",kRed, 1, ih, HitPos, HitsList, gdml_geom);
		    fHitsList->AddElement(HitsList);  
		    gEve->AddElement(fHitsList);
		    gEve->Redraw3D();
	    }
	}

	
}

void TEveEventDisplay::AddCaloCluster(const art::Event& event){
	 auto cluH = event.getValidHandle<mu2e::CaloClusterCollection>(cluTag_);
   	_clustercol = cluH.product(); 
	cout<<"has clusters "<<endl;
	cout<<_clustercol->size()<<endl;
	if (fClusters == 0) {
		fClusters= new TEveElementList("Hits");
		fClusters->IncDenyDestroy();     
	}
	else {
		fClusters->DestroyElements();  
	}
	
	if(!_clustercol->empty()){
		cout<<"in cluster loop "<<endl;
		TEveElementList* ClusterList  = new TEveElementList("Combo Hits");
		for(unsigned c = 0; c < _clustercol->size(); c++){
			 CaloCluster const& cluster = (*_clustercol)[c];
			 const CLHEP::Hep3Vector&      ClusterCOG = cluster.cog3Vector();
			 int nDisk = cluster.diskId();
			 DrawCluster("CaloCluster",kRed, 2, c, ClusterCOG, nDisk, ClusterList, gdml_geom);
			 fClusters->AddElement(ClusterList);  
			 gEve->AddElement(fClusters);
			 gEve->Redraw3D(); 
		}
  	}

}



bool TEveEventDisplay::FindData(const art::Event& evt){
	_chcol = 0; 
        auto chH = evt.getValidHandle<mu2e::ComboHitCollection>(chTag_);
	_chcol = chH.product();
	foundEvent = true;
	return _chcol != 0;
  }

void TEveEventDisplay::endJob(){
	if(foundEvent){
		char msg[300];
		sprintf(msg, "Reached end of file but #%i has not been found", true);
	        new TGMsgBox(gClient->GetRoot(), gClient->GetRoot(), "Event Not Found", msg, kMBIconExclamation,kMBOk);
	}

}  
	

}
using mu2e::TEveEventDisplay;
DEFINE_ART_MODULE(TEveEventDisplay);
