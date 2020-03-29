#ifndef EvtDisplayUtils_h
#define EvtDisplayUtils_h

#include <TObject.h>
#include <TApplication.h>
#include <TGTextBuffer.h>
#include <iostream>
#include "TROOT.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"


namespace mu2e {

 class EvtDisplayUtils  : public TObject
  {
#ifndef __CINT__
    public:
        explicit EvtDisplayUtils();
        void PrevEvent();
        void NextEvent();
        void GotoEvent();
        void Test();
     

      TGTextBuffer *fTbRun;
      TGTextBuffer *fTbEvt;
 
        
      virtual ~EvtDisplayUtils() {}
#endif
     ClassDef(EvtDisplayUtils,0);
  };
}
#endif /*TEveEventDisplay_inc_dict_classes_EvtDisplayUtils */
