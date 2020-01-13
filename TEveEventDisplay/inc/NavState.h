#ifndef TEveEventDisplay_inc_NavState_h
#define TEveEventDisplay_inc_NavState_h


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
	    static void SetTarget(int run, int event);
	    static int  TargetRun();
	    static int  TargetEvent();
	    virtual ~NavState() {};
       
	  private:
	    NavState() {};
	  
	};

#endif

