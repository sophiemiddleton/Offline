
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

    NavPanel::NavPanel(const TGWindow* p, UInt_t w, UInt_t h, fhicl::ParameterSet const &pset) : 
        TGMainFrame(p, w, h)
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
        //_mainCanvas = new TRootEmbeddedCanvas("EventDisplayCanvas",frmMain,1000,600);
	    frmMain->SetWindowName("EVT NAV");
	    frmMain->SetCleanup(kDeepCleanup);

	    TGHorizontalFrame* navFrame = new TGHorizontalFrame(frmMain);
	    TGVerticalFrame* evtidFrame = new TGVerticalFrame(frmMain);
	    {
		    TString icondir(TString::Format("%s/icons/", gSystem->Getenv("ROOTSYS")) );

		    // ... Create forward button and connect to "Exit" rcvr in visutils
		    TGTextButton *fin = new TGTextButton(navFrame,"&Exit","gApplication->Terminate(0)");
		    navFrame->AddFrame(fin);

		    
		    //....Add in check list
		    TGGroupFrame *options = new TGGroupFrame(navFrame, "Options", kVerticalFrame);
		    options->SetTitlePos(TGGroupFrame::kLeft);
            navFrame->AddFrame(options);

            TGTextButton *nextButton         = new TGTextButton(navFrame, "&Next", 1111);
            TGTextButton *quitButton         = new TGTextButton(navFrame, "&Quit", 1001);
            navFrame->AddFrame(quitButton, new TGLayoutHints(kLHintsLeft,3,0,3,0));
            navFrame->AddFrame(nextButton, new TGLayoutHints(kLHintsLeft,3,0,3,0));
                       
            quitButton->Associate(this);
            nextButton->Associate(this);
		    
		    // ... Create run num text entry widget and connect to "GotoEvent" rcvr in visutils
		    TGHorizontalFrame* runoFrame = new TGHorizontalFrame(evtidFrame);
		    fTlRun = new TGLabel(runoFrame,"Run Number");
		    fTlRun->SetTextJustify(kTextLeft);
		    fTlRun->SetMargins(5,5,5,0);
		    runoFrame->AddFrame(fTlRun);
        
		    fTeRun = new TGTextEntry(runoFrame, _runNumber = new TGTextBuffer(5), 1);
		    _runNumber->AddText(0, "1");
		    //fTeRun->Connect("ReturnPressed()","mu2e::EvtDisplayUtils", visutil_,"GotoEvent()");
		    runoFrame->AddFrame(fTeRun,new TGLayoutHints(kLHintsExpandX));

		    // ... Create evt num text entry widget and connect to "GotoEvent" rcvr in visutils
		    TGHorizontalFrame* evnoFrame = new TGHorizontalFrame(evtidFrame);
		    fTlEvt = new TGLabel(evnoFrame,"Evt Number");
		    fTlEvt->SetTextJustify(kTextLeft);
		    fTlEvt->SetMargins(5,5,5,0);
		    evnoFrame->AddFrame(fTlEvt);

             fTeEvt = new TGTextEntry(evnoFrame, _eventNumber = new TGTextBuffer(5), 1);
		    _eventNumber->AddText(0, "1");
		   /// fTeEvt->Connect("ReturnPressed()","mu2e::EvtDisplayUtils", visutil_,"GotoEvent()");
		    evnoFrame->AddFrame(fTeEvt,new TGLayoutHints(kLHintsExpandX));

		    //Add logo
		    std::string logoFile = "TEveEventDisplay/src/Icons/mu2e_logo_oval.png";
		    const TGPicture *logo = gClient->GetPicture(logoFile.c_str());
		    TGIcon *icon = new TGIcon(navFrame,logo,50,50);
		    navFrame->AddFrame(icon,new TGLayoutHints(kLHintsLeft,20,0,0,0));

		   

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
        //syntax for Format from http://root.cern.ch/phpBB3/viewtopic.php?t=8700
        gPad->AddExec("keyboardInput",TString::Format("((mu2e::TEveEventDisplay)%p)->keyboardInput()",this));
    }

    //void NavPanel::FillGeomView()
    //{}

    Bool_t NavPanel::ProcessMessage(Long_t msg, Long_t param1, Long_t param2){
      switch (GET_MSG(msg))
      {
        
        case kC_COMMAND:
          switch (GET_SUBMSG(msg))
          {
            case kCM_BUTTON: 
                 if(param1==1111)
                 {
                   std::cout<<"Next"<<std::endl;
                 }
                 if(param1==1001)
                 {
                   std::cout<<"Exit"<<std::endl;
                 }
                 
                break;
      }
      break;
  }
  return kTRUE;
}


    void NavPanel::setEvent(const art::Event& event, bool firstLoop)
    {
      _event=event.id().event();
      _subrun=event.id().subRun();
      _run=event.id().run();
      //TODO - call event drawing functions from collection interface here!
       fillEvent(firstLoop);
       //gApplication->Run(true);
    }

    void NavPanel::fillEvent(bool firstLoop)
    {
       // _findEvent=false;
       
        std::string eventInfoText;
        eventInfoText=Form("Event #: %i",_event);
        if(_eventNumberText==nullptr) 
        {
            _eventNumberText = new TText(0.6,-0.8, eventInfoText.c_str());
            _eventNumberText->SetTextColor(5);
            _eventNumberText->SetTextSize(0.025);
            _eventNumberText->Draw("same");
        }
        else _eventNumberText->SetTitle(eventInfoText.c_str());
        eventInfoText=Form("Sub Run #: %i",_subrun);
        if(_subrunNumberText==nullptr)
        {
            _subrunNumberText = new TText(0.6,-0.75,eventInfoText.c_str());
            _subrunNumberText->SetTextColor(5);
            _subrunNumberText->SetTextSize(0.025);
            _subrunNumberText->Draw("same");
        }
        else _subrunNumberText->SetTitle(eventInfoText.c_str());
        eventInfoText=Form("Run #: %i",_run);
          if(_runNumberText==nullptr)
          {
            _runNumberText = new TText(0.6,-0.7,eventInfoText.c_str());
            _runNumberText->SetTextColor(5);
            _runNumberText->SetTextSize(0.025);
            _runNumberText->Draw("same");
          }
        else _runNumberText->SetTitle(eventInfoText.c_str());

       //Collections Called Here

        this->Layout();

}

    int NavPanel::getEventToFind(bool &findEvent) const
    {
      findEvent=_findEvent;
      return _eventToFind;
    }

    bool NavPanel::isClosed() const
    {
        return _isClosed;
    }

}
