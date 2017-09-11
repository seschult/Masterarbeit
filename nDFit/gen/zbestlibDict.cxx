// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME gendIzbestlibDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
#include "src/ZBestNumber.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_ZBestNumber(void *p = 0);
   static void *newArray_ZBestNumber(Long_t size, void *p);
   static void delete_ZBestNumber(void *p);
   static void deleteArray_ZBestNumber(void *p);
   static void destruct_ZBestNumber(void *p);
   static void streamer_ZBestNumber(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ZBestNumber*)
   {
      ::ZBestNumber *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ZBestNumber >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ZBestNumber", ::ZBestNumber::Class_Version(), "ZBestNumber.h", 38,
                  typeid(::ZBestNumber), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ZBestNumber::Dictionary, isa_proxy, 16,
                  sizeof(::ZBestNumber) );
      instance.SetNew(&new_ZBestNumber);
      instance.SetNewArray(&newArray_ZBestNumber);
      instance.SetDelete(&delete_ZBestNumber);
      instance.SetDeleteArray(&deleteArray_ZBestNumber);
      instance.SetDestructor(&destruct_ZBestNumber);
      instance.SetStreamerFunc(&streamer_ZBestNumber);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ZBestNumber*)
   {
      return GenerateInitInstanceLocal((::ZBestNumber*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ZBestNumber*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEstringcOZBestNumbergR_Dictionary();
   static void pairlEstringcOZBestNumbergR_TClassManip(TClass*);
   static void *new_pairlEstringcOZBestNumbergR(void *p = 0);
   static void *newArray_pairlEstringcOZBestNumbergR(Long_t size, void *p);
   static void delete_pairlEstringcOZBestNumbergR(void *p);
   static void deleteArray_pairlEstringcOZBestNumbergR(void *p);
   static void destruct_pairlEstringcOZBestNumbergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,ZBestNumber>*)
   {
      pair<string,ZBestNumber> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,ZBestNumber>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,ZBestNumber>", "string", 96,
                  typeid(pair<string,ZBestNumber>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOZBestNumbergR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,ZBestNumber>) );
      instance.SetNew(&new_pairlEstringcOZBestNumbergR);
      instance.SetNewArray(&newArray_pairlEstringcOZBestNumbergR);
      instance.SetDelete(&delete_pairlEstringcOZBestNumbergR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOZBestNumbergR);
      instance.SetDestructor(&destruct_pairlEstringcOZBestNumbergR);

      ::ROOT::AddClassAlternate("pair<string,ZBestNumber>","pair<std::string,ZBestNumber>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const pair<string,ZBestNumber>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOZBestNumbergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const pair<string,ZBestNumber>*)0x0)->GetClass();
      pairlEstringcOZBestNumbergR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOZBestNumbergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *_Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR_Dictionary();
   static void _Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR_TClassManip(TClass*);
   static void *new__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(void *p = 0);
   static void *newArray__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(Long_t size, void *p);
   static void delete__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(void *p);
   static void deleteArray__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(void *p);
   static void destruct__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::_Rb_tree_iterator<pair<const string,ZBestNumber> >*)
   {
      ::_Rb_tree_iterator<pair<const string,ZBestNumber> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::_Rb_tree_iterator<pair<const string,ZBestNumber> >));
      static ::ROOT::TGenericClassInfo 
         instance("_Rb_tree_iterator<pair<const string,ZBestNumber> >", "set", 172,
                  typeid(::_Rb_tree_iterator<pair<const string,ZBestNumber> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &_Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR_Dictionary, isa_proxy, 1,
                  sizeof(::_Rb_tree_iterator<pair<const string,ZBestNumber> >) );
      instance.SetNew(&new__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR);
      instance.SetNewArray(&newArray__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR);
      instance.SetDelete(&delete__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR);
      instance.SetDeleteArray(&deleteArray__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR);
      instance.SetDestructor(&destruct__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR);

      ::ROOT::AddClassAlternate("_Rb_tree_iterator<pair<const string,ZBestNumber> >","map<std::string,ZBestNumber>::iterator");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::_Rb_tree_iterator<pair<const string,ZBestNumber> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *_Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::_Rb_tree_iterator<pair<const string,ZBestNumber> >*)0x0)->GetClass();
      _Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void _Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *__gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR_Dictionary();
   static void __gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR_TClassManip(TClass*);
   static void *new___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(void *p = 0);
   static void *newArray___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(Long_t size, void *p);
   static void delete___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(void *p);
   static void deleteArray___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(void *p);
   static void destruct___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >*)
   {
      ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >));
      static ::ROOT::TGenericClassInfo 
         instance("__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >", "string", 709,
                  typeid(::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &__gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR_Dictionary, isa_proxy, 1,
                  sizeof(::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >) );
      instance.SetNew(&new___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR);
      instance.SetNewArray(&newArray___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR);
      instance.SetDelete(&delete___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR);
      instance.SetDeleteArray(&deleteArray___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR);
      instance.SetDestructor(&destruct___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR);

      ::ROOT::AddClassAlternate("__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >","vector<ZBestNumber>::iterator");
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >*)
   {
      return GenerateInitInstanceLocal((::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *__gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >*)0x0)->GetClass();
      __gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void __gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ZBestNumber::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ZBestNumber::Class_Name()
{
   return "ZBestNumber";
}

//______________________________________________________________________________
const char *ZBestNumber::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ZBestNumber*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ZBestNumber::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ZBestNumber*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ZBestNumber::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ZBestNumber*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ZBestNumber::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ZBestNumber*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void ZBestNumber::Streamer(TBuffer &R__b)
{
   // Stream an object of class ZBestNumber.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ZBestNumber(void *p) {
      return  p ? new(p) ::ZBestNumber : new ::ZBestNumber;
   }
   static void *newArray_ZBestNumber(Long_t nElements, void *p) {
      return p ? new(p) ::ZBestNumber[nElements] : new ::ZBestNumber[nElements];
   }
   // Wrapper around operator delete
   static void delete_ZBestNumber(void *p) {
      delete ((::ZBestNumber*)p);
   }
   static void deleteArray_ZBestNumber(void *p) {
      delete [] ((::ZBestNumber*)p);
   }
   static void destruct_ZBestNumber(void *p) {
      typedef ::ZBestNumber current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_ZBestNumber(TBuffer &buf, void *obj) {
      ((::ZBestNumber*)obj)->::ZBestNumber::Streamer(buf);
   }
} // end of namespace ROOT for class ::ZBestNumber

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOZBestNumbergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,ZBestNumber> : new pair<string,ZBestNumber>;
   }
   static void *newArray_pairlEstringcOZBestNumbergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,ZBestNumber>[nElements] : new pair<string,ZBestNumber>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOZBestNumbergR(void *p) {
      delete ((pair<string,ZBestNumber>*)p);
   }
   static void deleteArray_pairlEstringcOZBestNumbergR(void *p) {
      delete [] ((pair<string,ZBestNumber>*)p);
   }
   static void destruct_pairlEstringcOZBestNumbergR(void *p) {
      typedef pair<string,ZBestNumber> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class pair<string,ZBestNumber>

namespace ROOT {
   // Wrappers around operator new
   static void *new__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::_Rb_tree_iterator<pair<const string,ZBestNumber> > : new ::_Rb_tree_iterator<pair<const string,ZBestNumber> >;
   }
   static void *newArray__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::_Rb_tree_iterator<pair<const string,ZBestNumber> >[nElements] : new ::_Rb_tree_iterator<pair<const string,ZBestNumber> >[nElements];
   }
   // Wrapper around operator delete
   static void delete__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(void *p) {
      delete ((::_Rb_tree_iterator<pair<const string,ZBestNumber> >*)p);
   }
   static void deleteArray__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(void *p) {
      delete [] ((::_Rb_tree_iterator<pair<const string,ZBestNumber> >*)p);
   }
   static void destruct__Rb_tree_iteratorlEpairlEconstsPstringcOZBestNumbergRsPgR(void *p) {
      typedef ::_Rb_tree_iterator<pair<const string,ZBestNumber> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::_Rb_tree_iterator<pair<const string,ZBestNumber> >

namespace ROOT {
   // Wrappers around operator new
   static void *new___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> > : new ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >;
   }
   static void *newArray___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >[nElements] : new ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >[nElements];
   }
   // Wrapper around operator delete
   static void delete___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(void *p) {
      delete ((::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >*)p);
   }
   static void deleteArray___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(void *p) {
      delete [] ((::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >*)p);
   }
   static void destruct___gnu_cxxcLcL__normal_iteratorlEZBestNumbermUcOvectorlEZBestNumbergRsPgR(void *p) {
      typedef ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >

namespace ROOT {
   static TClass *vectorlEZBestNumbergR_Dictionary();
   static void vectorlEZBestNumbergR_TClassManip(TClass*);
   static void *new_vectorlEZBestNumbergR(void *p = 0);
   static void *newArray_vectorlEZBestNumbergR(Long_t size, void *p);
   static void delete_vectorlEZBestNumbergR(void *p);
   static void deleteArray_vectorlEZBestNumbergR(void *p);
   static void destruct_vectorlEZBestNumbergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ZBestNumber>*)
   {
      vector<ZBestNumber> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ZBestNumber>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ZBestNumber>", -2, "vector", 214,
                  typeid(vector<ZBestNumber>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEZBestNumbergR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ZBestNumber>) );
      instance.SetNew(&new_vectorlEZBestNumbergR);
      instance.SetNewArray(&newArray_vectorlEZBestNumbergR);
      instance.SetDelete(&delete_vectorlEZBestNumbergR);
      instance.SetDeleteArray(&deleteArray_vectorlEZBestNumbergR);
      instance.SetDestructor(&destruct_vectorlEZBestNumbergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ZBestNumber> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<ZBestNumber>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEZBestNumbergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<ZBestNumber>*)0x0)->GetClass();
      vectorlEZBestNumbergR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEZBestNumbergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEZBestNumbergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ZBestNumber> : new vector<ZBestNumber>;
   }
   static void *newArray_vectorlEZBestNumbergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ZBestNumber>[nElements] : new vector<ZBestNumber>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEZBestNumbergR(void *p) {
      delete ((vector<ZBestNumber>*)p);
   }
   static void deleteArray_vectorlEZBestNumbergR(void *p) {
      delete [] ((vector<ZBestNumber>*)p);
   }
   static void destruct_vectorlEZBestNumbergR(void *p) {
      typedef vector<ZBestNumber> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<ZBestNumber>

namespace ROOT {
   static TClass *maplEstringcOZBestNumbergR_Dictionary();
   static void maplEstringcOZBestNumbergR_TClassManip(TClass*);
   static void *new_maplEstringcOZBestNumbergR(void *p = 0);
   static void *newArray_maplEstringcOZBestNumbergR(Long_t size, void *p);
   static void delete_maplEstringcOZBestNumbergR(void *p);
   static void deleteArray_maplEstringcOZBestNumbergR(void *p);
   static void destruct_maplEstringcOZBestNumbergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,ZBestNumber>*)
   {
      map<string,ZBestNumber> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,ZBestNumber>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,ZBestNumber>", -2, "map", 96,
                  typeid(map<string,ZBestNumber>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOZBestNumbergR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,ZBestNumber>) );
      instance.SetNew(&new_maplEstringcOZBestNumbergR);
      instance.SetNewArray(&newArray_maplEstringcOZBestNumbergR);
      instance.SetDelete(&delete_maplEstringcOZBestNumbergR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOZBestNumbergR);
      instance.SetDestructor(&destruct_maplEstringcOZBestNumbergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,ZBestNumber> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<string,ZBestNumber>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOZBestNumbergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,ZBestNumber>*)0x0)->GetClass();
      maplEstringcOZBestNumbergR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOZBestNumbergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOZBestNumbergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,ZBestNumber> : new map<string,ZBestNumber>;
   }
   static void *newArray_maplEstringcOZBestNumbergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,ZBestNumber>[nElements] : new map<string,ZBestNumber>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOZBestNumbergR(void *p) {
      delete ((map<string,ZBestNumber>*)p);
   }
   static void deleteArray_maplEstringcOZBestNumbergR(void *p) {
      delete [] ((map<string,ZBestNumber>*)p);
   }
   static void destruct_maplEstringcOZBestNumbergR(void *p) {
      typedef map<string,ZBestNumber> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,ZBestNumber>

namespace {
  void TriggerDictionaryInitialization_zbestlibDict_Impl() {
    static const char* headers[] = {
"config.h",
"src/ZBestNumber.h",
0
    };
    static const char* includePaths[] = {
"lib",
"./lib",
"src",
"./src",
"Tool",
"/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.06.02-x86_64-slc6-gcc49-opt/include",
"/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.06.02-x86_64-slc6-gcc49-opt/include",
"/ptmp/mpp/sschulte/Topmass13TeV.git/trunk/CodeFactory/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "zbestlibDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$src/ZBestNumber.h")))  ZBestNumber;
namespace std{template <class _CharT> struct __attribute__((annotate("$clingAutoload$string")))  char_traits;
}
namespace std{template <typename > class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "zbestlibDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "config.h"
#include "src/ZBestNumber.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"ZBestNumber", payloadCode, "@",
"__gnu_cxx::__normal_iterator<ZBestNumber*,vector<ZBestNumber> >", payloadCode, "@",
"vector<ZBestNumber>::iterator", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("zbestlibDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_zbestlibDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_zbestlibDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_zbestlibDict() {
  TriggerDictionaryInitialization_zbestlibDict_Impl();
}
