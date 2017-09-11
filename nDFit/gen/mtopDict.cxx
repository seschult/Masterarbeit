// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME gendImtopDict

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
#include "config.h"
#include "Code/mtop_fit.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_mtop_fit(void *p = 0);
   static void *newArray_mtop_fit(Long_t size, void *p);
   static void delete_mtop_fit(void *p);
   static void deleteArray_mtop_fit(void *p);
   static void destruct_mtop_fit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::mtop_fit*)
   {
      ::mtop_fit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::mtop_fit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("mtop_fit", ::mtop_fit::Class_Version(), "Code/mtop_fit.h", 6,
                  typeid(::mtop_fit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::mtop_fit::Dictionary, isa_proxy, 4,
                  sizeof(::mtop_fit) );
      instance.SetNew(&new_mtop_fit);
      instance.SetNewArray(&newArray_mtop_fit);
      instance.SetDelete(&delete_mtop_fit);
      instance.SetDeleteArray(&deleteArray_mtop_fit);
      instance.SetDestructor(&destruct_mtop_fit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::mtop_fit*)
   {
      return GenerateInitInstanceLocal((::mtop_fit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::mtop_fit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr mtop_fit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *mtop_fit::Class_Name()
{
   return "mtop_fit";
}

//______________________________________________________________________________
const char *mtop_fit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::mtop_fit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int mtop_fit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::mtop_fit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *mtop_fit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::mtop_fit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *mtop_fit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::mtop_fit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void mtop_fit::Streamer(TBuffer &R__b)
{
   // Stream an object of class mtop_fit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(mtop_fit::Class(),this);
   } else {
      R__b.WriteClassBuffer(mtop_fit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_mtop_fit(void *p) {
      return  p ? new(p) ::mtop_fit : new ::mtop_fit;
   }
   static void *newArray_mtop_fit(Long_t nElements, void *p) {
      return p ? new(p) ::mtop_fit[nElements] : new ::mtop_fit[nElements];
   }
   // Wrapper around operator delete
   static void delete_mtop_fit(void *p) {
      delete ((::mtop_fit*)p);
   }
   static void deleteArray_mtop_fit(void *p) {
      delete [] ((::mtop_fit*)p);
   }
   static void destruct_mtop_fit(void *p) {
      typedef ::mtop_fit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::mtop_fit

namespace {
  void TriggerDictionaryInitialization_mtopDict_Impl() {
    static const char* headers[] = {
"config.h",
"Code/mtop_fit.h",
0
    };
    static const char* includePaths[] = {
"lib",
"./lib",
"src",
"./src",
"Tool",
"/home/sebastian/root/include",
"/home/sebastian/root/include",
"/home/sebastian/Topmass13TeV.git/trunk/CodeFactory/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "mtopDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$Code/mtop_fit.h")))  mtop_fit;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "mtopDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "config.h"
#include "Code/mtop_fit.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"mtop_fit", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("mtopDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_mtopDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_mtopDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_mtopDict() {
  TriggerDictionaryInitialization_mtopDict_Impl();
}
