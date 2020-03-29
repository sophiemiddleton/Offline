
#include <TObject.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGIcon.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGString.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
#include<TGLViewer.h>
#include <TGMsgBox.h>

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

#include "TEveEventDisplay/src/dict_classes/NavPanel.h"

namespace fhicl
{
  class ParameterSet;
}

using namespace mu2e;
namespace mu2e{

	NavPanel::NavPanel() : TGMainFrame(gClient->GetRoot(), 320, 320) {}

    void NavPanel::MakeNavPanel(EvtDisplayUtils *visutil_)
    {
        gClient->GetRoot();
	    TEveBrowser* browser = gEve->GetBrowser();
        FontStruct_t buttonfont = gClient->GetFontByName("-*-helvetica-medium-r-*-*-8-*-*-*-*-*-iso8859-1");
        GCValues_t gval;
        gval.fMask = kGCForeground | kGCFont;
        gval.fFont = gVirtualX->GetFontHandle(buttonfont);
        gClient->GetColorByName("black", gval.fForeground);
       
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

		    // ... Create forward button and connect to "Exit" rcvr in visutils
		    TGTextButton *fin = new TGTextButton(navFrame,"&Exit","gApplication->Terminate(0)");
		    navFrame->AddFrame(fin);

		    
		    //....Add in check list
		    TGGroupFrame *options = new TGGroupFrame(navFrame, "Options", kVerticalFrame);
		    options->SetTitlePos(TGGroupFrame::kLeft);
		    
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


/*
    Bool_t NavPanel::ProcessMessage(Long_t msg, Long_t param1, Long_t param2){
      switch (GET_MSG(msg))
      {
        
        case kC_COMMAND:
          switch (GET_SUBMSG(msg))
          {
            case kCM_BUTTON: 
                 if(param1==1111)
                 {
                   gApplication->Terminate();
                 }
                 if(param1==1112)
                 {
                   
                   fillEvent();
                 }
                 
                 if(param1==1102)
                 {
                   //_eventToFind=atoi(_eventToFindField->GetText());
                  // _findEvent=true;
                   
                   gApplication->Terminate();
                 }
                break;
      }
      break;
  }
  return kTRUE;
}


    void NavPanel::setEvent(const art::Event& event, bool firstLoop)
    {
      _eventNumber=event.id().event();
      _subrunNumber=event.id().subRun();
      _runNumber=event.id().run();
      //TODO - call event drawing functions from collection interface here!
       gApplication->Run(true);
    }
*/
}
