// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME EventDict

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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "../Source/Core/LMCEvent.hh"

// Header files passed via #pragma extra_include

namespace locust {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *locust_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("locust", 0 /*version*/, "../Source/Core/LMCTrack.hh", 14,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &locust_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *locust_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static void *new_locustcLcLEvent(void *p = 0);
   static void *newArray_locustcLcLEvent(Long_t size, void *p);
   static void delete_locustcLcLEvent(void *p);
   static void deleteArray_locustcLcLEvent(void *p);
   static void destruct_locustcLcLEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::locust::Event*)
   {
      ::locust::Event *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::locust::Event >(0);
      static ::ROOT::TGenericClassInfo 
         instance("locust::Event", ::locust::Event::Class_Version(), "../Source/Core/LMCEvent.hh", 24,
                  typeid(::locust::Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::locust::Event::Dictionary, isa_proxy, 4,
                  sizeof(::locust::Event) );
      instance.SetNew(&new_locustcLcLEvent);
      instance.SetNewArray(&newArray_locustcLcLEvent);
      instance.SetDelete(&delete_locustcLcLEvent);
      instance.SetDeleteArray(&deleteArray_locustcLcLEvent);
      instance.SetDestructor(&destruct_locustcLcLEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::locust::Event*)
   {
      return GenerateInitInstanceLocal((::locust::Event*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::locust::Event*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace locust {
//______________________________________________________________________________
atomic_TClass_ptr Event::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Event::Class_Name()
{
   return "locust::Event";
}

//______________________________________________________________________________
const char *Event::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::locust::Event*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Event::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::locust::Event*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Event::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::locust::Event*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Event::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::locust::Event*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace locust
namespace locust {
//______________________________________________________________________________
void Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class locust::Event.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(locust::Event::Class(),this);
   } else {
      R__b.WriteClassBuffer(locust::Event::Class(),this);
   }
}

} // namespace locust
namespace ROOT {
   // Wrappers around operator new
   static void *new_locustcLcLEvent(void *p) {
      return  p ? new(p) ::locust::Event : new ::locust::Event;
   }
   static void *newArray_locustcLcLEvent(Long_t nElements, void *p) {
      return p ? new(p) ::locust::Event[nElements] : new ::locust::Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_locustcLcLEvent(void *p) {
      delete ((::locust::Event*)p);
   }
   static void deleteArray_locustcLcLEvent(void *p) {
      delete [] ((::locust::Event*)p);
   }
   static void destruct_locustcLcLEvent(void *p) {
      typedef ::locust::Event current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::locust::Event

namespace ROOT {
   static TClass *vectorlElocustcLcLTrackgR_Dictionary();
   static void vectorlElocustcLcLTrackgR_TClassManip(TClass*);
   static void *new_vectorlElocustcLcLTrackgR(void *p = 0);
   static void *newArray_vectorlElocustcLcLTrackgR(Long_t size, void *p);
   static void delete_vectorlElocustcLcLTrackgR(void *p);
   static void deleteArray_vectorlElocustcLcLTrackgR(void *p);
   static void destruct_vectorlElocustcLcLTrackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<locust::Track>*)
   {
      vector<locust::Track> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<locust::Track>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<locust::Track>", -2, "vector", 216,
                  typeid(vector<locust::Track>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElocustcLcLTrackgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<locust::Track>) );
      instance.SetNew(&new_vectorlElocustcLcLTrackgR);
      instance.SetNewArray(&newArray_vectorlElocustcLcLTrackgR);
      instance.SetDelete(&delete_vectorlElocustcLcLTrackgR);
      instance.SetDeleteArray(&deleteArray_vectorlElocustcLcLTrackgR);
      instance.SetDestructor(&destruct_vectorlElocustcLcLTrackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<locust::Track> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<locust::Track>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElocustcLcLTrackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<locust::Track>*)0x0)->GetClass();
      vectorlElocustcLcLTrackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElocustcLcLTrackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElocustcLcLTrackgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<locust::Track> : new vector<locust::Track>;
   }
   static void *newArray_vectorlElocustcLcLTrackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<locust::Track>[nElements] : new vector<locust::Track>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElocustcLcLTrackgR(void *p) {
      delete ((vector<locust::Track>*)p);
   }
   static void deleteArray_vectorlElocustcLcLTrackgR(void *p) {
      delete [] ((vector<locust::Track>*)p);
   }
   static void destruct_vectorlElocustcLcLTrackgR(void *p) {
      typedef vector<locust::Track> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<locust::Track>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 216,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_EventDict_Impl() {
    static const char* headers[] = {
"../Source/Core/LMCEvent.hh",
0
    };
    static const char* includePaths[] = {
"/usr/local/include",
"/home/penny/locust_mc/cbuild/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "EventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace locust{class __attribute__((annotate(R"ATTRDUMP(Root syntax.)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$../Source/Core/LMCEvent.hh")))  Event;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "EventDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "../Source/Core/LMCEvent.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"locust::Event", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("EventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_EventDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_EventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_EventDict() {
  TriggerDictionaryInitialization_EventDict_Impl();
}
