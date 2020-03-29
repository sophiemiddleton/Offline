//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//

#include <TObject.h>
#include <TApplication.h>
#include <TGTextBuffer.h>
#include <iostream>

#include "TROOT.h"
#include "TApplication.h"
#include "TEveEventDisplay/src/dict_classes/NavState.h"
using namespace mu2e;

  static int gsNavState    = 0;
  static int gsTargetRun   = 0;
  static int gsTargetEvent = 0;
 
  //......................................................................

  int NavState::Which() { return gsNavState; }

  //......................................................................

  void NavState::Set(int which)
  {
    
    if(gsNavState != kSEQUENTIAL_ONLY)
        gsNavState = which;
    else gROOT->GetApplication()->Terminate();
  }

  //......................................................................

  void NavState::SetTarget(int run, int event)
  {
    gsTargetRun = run;
    gsTargetEvent = event;
  }

 
  int NavState::TargetRun() { return gsTargetRun; }

  int NavState::TargetEvent() { return gsTargetEvent; }

  void NavState::Print(){

        std::cout<<"The print button has been pressed on screen"<<std::endl;
   }

 /*void NavState::SetArtTarget(art::Run& run, art::Event& event)
  {
    runID= run;
    eventID = event;
  }
  art::Run *NavState::TargetArtRun() { return runID; }

  art::Event *NavState::TargetArtEvent() { return eventID; }

 void NavState::PreviousEvent() {

    GotoEvent(eventNum_ - 1);
   }

  void NavState::NextEvent() {

    GotoEvent(eventNum_ + 1);
  }


 void NavState::GotoEvent(int event) {

        eventNum_ = event;
        std::cout << "[ In NavState::GoToEvent() ] : Loading event " << event << "... " << std::flush;
        
    }
*/

