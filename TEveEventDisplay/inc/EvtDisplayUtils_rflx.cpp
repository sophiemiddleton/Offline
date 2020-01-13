// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TEveEventDisplaydIincdIEvtDisplayUtils_rflx
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TEveEventDisplay/inc/EvtDisplayUtils.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_EvtDisplayUtils(void *p = 0);
   static void *newArray_EvtDisplayUtils(Long_t size, void *p);
   static void delete_EvtDisplayUtils(void *p);
   static void deleteArray_EvtDisplayUtils(void *p);
   static void destruct_EvtDisplayUtils(void *p);
   static void streamer_EvtDisplayUtils(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EvtDisplayUtils*)
   {
      ::EvtDisplayUtils *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EvtDisplayUtils >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EvtDisplayUtils", ::EvtDisplayUtils::Class_Version(), "", 12,
                  typeid(::EvtDisplayUtils), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EvtDisplayUtils::Dictionary, isa_proxy, 16,
                  sizeof(::EvtDisplayUtils) );
      instance.SetNew(&new_EvtDisplayUtils);
      instance.SetNewArray(&newArray_EvtDisplayUtils);
      instance.SetDelete(&delete_EvtDisplayUtils);
      instance.SetDeleteArray(&deleteArray_EvtDisplayUtils);
      instance.SetDestructor(&destruct_EvtDisplayUtils);
      instance.SetStreamerFunc(&streamer_EvtDisplayUtils);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EvtDisplayUtils*)
   {
      return GenerateInitInstanceLocal((::EvtDisplayUtils*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EvtDisplayUtils*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr EvtDisplayUtils::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EvtDisplayUtils::Class_Name()
{
   return "EvtDisplayUtils";
}

//______________________________________________________________________________
const char *EvtDisplayUtils::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EvtDisplayUtils*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EvtDisplayUtils::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EvtDisplayUtils*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EvtDisplayUtils::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EvtDisplayUtils*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EvtDisplayUtils::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EvtDisplayUtils*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void EvtDisplayUtils::Streamer(TBuffer &R__b)
{
   // Stream an object of class EvtDisplayUtils.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EvtDisplayUtils::Class(),this);
   } else {
      R__b.WriteClassBuffer(EvtDisplayUtils::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EvtDisplayUtils(void *p) {
      return  p ? new(p) ::EvtDisplayUtils : new ::EvtDisplayUtils;
   }
   static void *newArray_EvtDisplayUtils(Long_t nElements, void *p) {
      return p ? new(p) ::EvtDisplayUtils[nElements] : new ::EvtDisplayUtils[nElements];
   }
   // Wrapper around operator delete
   static void delete_EvtDisplayUtils(void *p) {
      delete ((::EvtDisplayUtils*)p);
   }
   static void deleteArray_EvtDisplayUtils(void *p) {
      delete [] ((::EvtDisplayUtils*)p);
   }
   static void destruct_EvtDisplayUtils(void *p) {
      typedef ::EvtDisplayUtils current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_EvtDisplayUtils(TBuffer &buf, void *obj) {
      ((::EvtDisplayUtils*)obj)->::EvtDisplayUtils::Streamer(buf);
   }
} // end of namespace ROOT for class ::EvtDisplayUtils

namespace {
  void TriggerDictionaryInitialization_lib_mu2e_TEveEventDisplay_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
"/cvmfs/mu2e.opensciencegrid.org/artexternals/root/v6_18_04c/Linux64bit+3.10-2.17-e19-prof/include",
"/home/sophie/Offlinesophie/MDC/Offline/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "lib_mu2e_TEveEventDisplay dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class EvtDisplayUtils;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "lib_mu2e_TEveEventDisplay dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#ifndef TEveEventDisplay_inc_EvtDisplayUtils_h
#define TEveEventDisplay_inc_EvtDisplayUtils_h

#include <TObject.h>
#include <TApplication.h>
#include <TGTextBuffer.h>
#include <iostream>
#include "TROOT.h"

 class EvtDisplayUtils  : public TObject
  {
#ifndef __CINT__
    public:
#endif
    explicit EvtDisplayUtils();
    public:
#ifndef __CINT__
      void PrevEvent();
      void NextEvent();
      void GotoEvent();
     
      TGTextBuffer *fTbRun;
      TGTextBuffer *fTbEvt;

      virtual ~EvtDisplayUtils() {}
#endif
     ClassDef(EvtDisplayUtils,0);
  };

#endif /*TEveEventDisplay_inc_EvtDisplayUtils */

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"EvtDisplayUtils", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("lib_mu2e_TEveEventDisplay",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_lib_mu2e_TEveEventDisplay_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_lib_mu2e_TEveEventDisplay_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_lib_mu2e_TEveEventDisplay() {
  TriggerDictionaryInitialization_lib_mu2e_TEveEventDisplay_Impl();
}
