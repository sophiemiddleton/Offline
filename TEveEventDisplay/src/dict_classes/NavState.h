#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#ifndef NavState_h
#define NavState_h

namespace mu2e{

	enum nav_states_ {
	    kNEXT_EVENT,
	    kPREV_EVENT,
	    kRELOAD_EVENT,
	    kGOTO_EVENT
	};

	class NavState {
        
	  public:
	    static int  Which();
	    static void Set(int which);
	    static void SetArtTarget(art::Run &run, art::Event& event);
	    static void SetTarget(int run, int event);
	    static int  TargetRun();
	    static int  TargetEvent();
	    art::Run *TargetArtRun();
	    art::Event *TargetArtEvent();

	    virtual ~NavState() {};
       
	  private:
	    NavState() {};
	  
	};
}
#endif

