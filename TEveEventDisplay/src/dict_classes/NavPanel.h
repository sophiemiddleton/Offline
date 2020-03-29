#ifndef NavPanel_h
#define NavPanel_h

#include <TGLabel.h>
#include <TGTextEntry.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "ConfigTools/inc/SimpleConfig.hh"

#include "TEveEventDisplay/src/dict_classes/EvtDisplayUtils.h"

namespace mu2e{
	class NavPanel : public TGMainFrame {
            
               public:
		         #ifndef __CINT__
		         explicit NavPanel();
		         virtual ~NavPanel(){};
		         void MakeNavPanel(EvtDisplayUtils *_visutil);
                 #endif
                 TGTextEntry      *fTeRun,*fTeEvt;
		         TGLabel          *fTlRun,*fTlEvt;
                 int _eventNumber, _subrunNumber, _runNumber;
		         ClassDef(NavPanel,0);

	}; //end class def

}//end namespace mu2e

#endif /*NavPanel.h*/
