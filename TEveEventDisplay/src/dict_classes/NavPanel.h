#ifndef NavPanel_h
#define NavPanel_h

#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TText.h>
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "ConfigTools/inc/SimpleConfig.hh"

class TBox;
class TGTextEntry;
class TPad;
class TGCanvas;
class TRootEmbeddedCanvas;
class TGTextButton;
class TText;
namespace mu2e{
	class NavPanel : public TGMainFrame {
                public:
		        #ifndef __CINT__
                 NavPanel();
                 NavPanel(const NavPanel &);
                 NavPanel& operator=(const NavPanel &);

                 NavPanel(const TGWindow* p, UInt_t w, UInt_t h, fhicl::ParameterSet const &pset);
            
		         virtual ~NavPanel(){};

                 Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2);
                 void  setEvent(const art::Event& event, bool firstLoop=false);
                 void  fillEvent(bool firstLoop=false);
                 bool  isClosed() const;
                 int   getEventToFind(bool &findEvent) const;
                 #endif
                 TGTextEntry      *fTeRun,*fTeEvt;
		         TGLabel          *fTlRun,*fTlEvt;
                 TGTextBuffer *_eventNumber, *_subrunNumber, *_runNumber;
                 int  _eventToFind = 0;
                 bool  _showBuilding = false;
                 bool _showCRV=false;
                 bool _showDSOnly = true;
                 bool _showTracker = true;
                 bool _showCalo = true;
                 bool _isClosed = false;
                 bool _findEvent = true;
                 //TGHorizontalFrame* navFrame;
	             //TGVerticalFrame* evtidFrame;
                 //TRootEmbeddedCanvas *_mainCanvas;
                 //TGCanvas  *_infoCanvas;
                 //TPad *_mainPad, *_infoPad;
                 TText  *_eventNumberText, *_subrunNumberText, *_runNumberText;
                 int _event, _subrun, _run;
		         ClassDef(NavPanel,0);

	}; //end class def

}//end namespace mu2e

#endif /*NavPanel.h*/
