#define private public
#define protected public
/** \class ExRootAnalysisLinkDef
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  $Date: 2006/12/19 19:18:46 $
 *  $Revision: 1.5 $
 *
 *  
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"

#include "ExRootAnalysis/ExRootFactory.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TSortableObject+;
#pragma link C++ class TRootWeight;
#pragma link C++ class TRootLHEFEvent+;
#pragma link C++ class TRootLHEFParticle+;
#pragma link C++ class TRootGenEvent+;
#pragma link C++ class TRootGenParticle+;
#pragma link C++ class TRootGenJet+;
#pragma link C++ class TRootEvent+;
#pragma link C++ class TRootMissingET+;
#pragma link C++ class TRootPhoton+;
#pragma link C++ class TRootElectron+;
#pragma link C++ class TRootMuon+;
#pragma link C++ class TRootTau+;
#pragma link C++ class TRootJet+;

#pragma link C++ class ExRootTreeReader+;
#pragma link C++ class ExRootTreeBranch+;
#pragma link C++ class ExRootTreeWriter+;
#pragma link C++ class ExRootResult+;
#pragma link C++ class ExRootClassifier+;
#pragma link C++ class ExRootFilter+;

#pragma link C++ class ExRootFactory+;

#pragma link C++ function HistStyle;
#pragma link C++ function FillChain;

#endif

// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME tmpdIsrcdIExRootAnalysisDict

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

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_TSortableObject(void *p);
   static void deleteArray_TSortableObject(void *p);
   static void destruct_TSortableObject(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TSortableObject*)
   {
      ::TSortableObject *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TSortableObject >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TSortableObject", ::TSortableObject::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 37,
                  typeid(::TSortableObject), DefineBehavior(ptr, ptr),
                  &::TSortableObject::Dictionary, isa_proxy, 4,
                  sizeof(::TSortableObject) );
      instance.SetDelete(&delete_TSortableObject);
      instance.SetDeleteArray(&deleteArray_TSortableObject);
      instance.SetDestructor(&destruct_TSortableObject);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TSortableObject*)
   {
      return GenerateInitInstanceLocal((::TSortableObject*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TSortableObject*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootWeight(void *p = 0);
   static void *newArray_TRootWeight(Long_t size, void *p);
   static void delete_TRootWeight(void *p);
   static void deleteArray_TRootWeight(void *p);
   static void destruct_TRootWeight(void *p);
   static void streamer_TRootWeight(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootWeight*)
   {
      ::TRootWeight *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootWeight >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootWeight", ::TRootWeight::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 51,
                  typeid(::TRootWeight), DefineBehavior(ptr, ptr),
                  &::TRootWeight::Dictionary, isa_proxy, 16,
                  sizeof(::TRootWeight) );
      instance.SetNew(&new_TRootWeight);
      instance.SetNewArray(&newArray_TRootWeight);
      instance.SetDelete(&delete_TRootWeight);
      instance.SetDeleteArray(&deleteArray_TRootWeight);
      instance.SetDestructor(&destruct_TRootWeight);
      instance.SetStreamerFunc(&streamer_TRootWeight);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootWeight*)
   {
      return GenerateInitInstanceLocal((::TRootWeight*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootWeight*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootLHEFEvent(void *p = 0);
   static void *newArray_TRootLHEFEvent(Long_t size, void *p);
   static void delete_TRootLHEFEvent(void *p);
   static void deleteArray_TRootLHEFEvent(void *p);
   static void destruct_TRootLHEFEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootLHEFEvent*)
   {
      ::TRootLHEFEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootLHEFEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootLHEFEvent", ::TRootLHEFEvent::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 61,
                  typeid(::TRootLHEFEvent), DefineBehavior(ptr, ptr),
                  &::TRootLHEFEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TRootLHEFEvent) );
      instance.SetNew(&new_TRootLHEFEvent);
      instance.SetNewArray(&newArray_TRootLHEFEvent);
      instance.SetDelete(&delete_TRootLHEFEvent);
      instance.SetDeleteArray(&deleteArray_TRootLHEFEvent);
      instance.SetDestructor(&destruct_TRootLHEFEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootLHEFEvent*)
   {
      return GenerateInitInstanceLocal((::TRootLHEFEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootLHEFEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootLHEFParticle(void *p = 0);
   static void *newArray_TRootLHEFParticle(Long_t size, void *p);
   static void delete_TRootLHEFParticle(void *p);
   static void deleteArray_TRootLHEFParticle(void *p);
   static void destruct_TRootLHEFParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootLHEFParticle*)
   {
      ::TRootLHEFParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootLHEFParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootLHEFParticle", ::TRootLHEFParticle::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 80,
                  typeid(::TRootLHEFParticle), DefineBehavior(ptr, ptr),
                  &::TRootLHEFParticle::Dictionary, isa_proxy, 4,
                  sizeof(::TRootLHEFParticle) );
      instance.SetNew(&new_TRootLHEFParticle);
      instance.SetNewArray(&newArray_TRootLHEFParticle);
      instance.SetDelete(&delete_TRootLHEFParticle);
      instance.SetDeleteArray(&deleteArray_TRootLHEFParticle);
      instance.SetDestructor(&destruct_TRootLHEFParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootLHEFParticle*)
   {
      return GenerateInitInstanceLocal((::TRootLHEFParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootLHEFParticle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootGenEvent(void *p = 0);
   static void *newArray_TRootGenEvent(Long_t size, void *p);
   static void delete_TRootGenEvent(void *p);
   static void deleteArray_TRootGenEvent(void *p);
   static void destruct_TRootGenEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootGenEvent*)
   {
      ::TRootGenEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootGenEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootGenEvent", ::TRootGenEvent::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 119,
                  typeid(::TRootGenEvent), DefineBehavior(ptr, ptr),
                  &::TRootGenEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TRootGenEvent) );
      instance.SetNew(&new_TRootGenEvent);
      instance.SetNewArray(&newArray_TRootGenEvent);
      instance.SetDelete(&delete_TRootGenEvent);
      instance.SetDeleteArray(&deleteArray_TRootGenEvent);
      instance.SetDestructor(&destruct_TRootGenEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootGenEvent*)
   {
      return GenerateInitInstanceLocal((::TRootGenEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootGenEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootGenParticle(void *p = 0);
   static void *newArray_TRootGenParticle(Long_t size, void *p);
   static void delete_TRootGenParticle(void *p);
   static void deleteArray_TRootGenParticle(void *p);
   static void destruct_TRootGenParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootGenParticle*)
   {
      ::TRootGenParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootGenParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootGenParticle", ::TRootGenParticle::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 130,
                  typeid(::TRootGenParticle), DefineBehavior(ptr, ptr),
                  &::TRootGenParticle::Dictionary, isa_proxy, 4,
                  sizeof(::TRootGenParticle) );
      instance.SetNew(&new_TRootGenParticle);
      instance.SetNewArray(&newArray_TRootGenParticle);
      instance.SetDelete(&delete_TRootGenParticle);
      instance.SetDeleteArray(&deleteArray_TRootGenParticle);
      instance.SetDestructor(&destruct_TRootGenParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootGenParticle*)
   {
      return GenerateInitInstanceLocal((::TRootGenParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootGenParticle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootGenJet(void *p = 0);
   static void *newArray_TRootGenJet(Long_t size, void *p);
   static void delete_TRootGenJet(void *p);
   static void deleteArray_TRootGenJet(void *p);
   static void destruct_TRootGenJet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootGenJet*)
   {
      ::TRootGenJet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootGenJet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootGenJet", ::TRootGenJet::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 164,
                  typeid(::TRootGenJet), DefineBehavior(ptr, ptr),
                  &::TRootGenJet::Dictionary, isa_proxy, 4,
                  sizeof(::TRootGenJet) );
      instance.SetNew(&new_TRootGenJet);
      instance.SetNewArray(&newArray_TRootGenJet);
      instance.SetDelete(&delete_TRootGenJet);
      instance.SetDeleteArray(&deleteArray_TRootGenJet);
      instance.SetDestructor(&destruct_TRootGenJet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootGenJet*)
   {
      return GenerateInitInstanceLocal((::TRootGenJet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootGenJet*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootEvent(void *p = 0);
   static void *newArray_TRootEvent(Long_t size, void *p);
   static void delete_TRootEvent(void *p);
   static void deleteArray_TRootEvent(void *p);
   static void destruct_TRootEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootEvent*)
   {
      ::TRootEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootEvent", ::TRootEvent::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 189,
                  typeid(::TRootEvent), DefineBehavior(ptr, ptr),
                  &::TRootEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TRootEvent) );
      instance.SetNew(&new_TRootEvent);
      instance.SetNewArray(&newArray_TRootEvent);
      instance.SetDelete(&delete_TRootEvent);
      instance.SetDeleteArray(&deleteArray_TRootEvent);
      instance.SetDestructor(&destruct_TRootEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootEvent*)
   {
      return GenerateInitInstanceLocal((::TRootEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootMissingET(void *p = 0);
   static void *newArray_TRootMissingET(Long_t size, void *p);
   static void delete_TRootMissingET(void *p);
   static void deleteArray_TRootMissingET(void *p);
   static void destruct_TRootMissingET(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootMissingET*)
   {
      ::TRootMissingET *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootMissingET >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootMissingET", ::TRootMissingET::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 201,
                  typeid(::TRootMissingET), DefineBehavior(ptr, ptr),
                  &::TRootMissingET::Dictionary, isa_proxy, 4,
                  sizeof(::TRootMissingET) );
      instance.SetNew(&new_TRootMissingET);
      instance.SetNewArray(&newArray_TRootMissingET);
      instance.SetDelete(&delete_TRootMissingET);
      instance.SetDeleteArray(&deleteArray_TRootMissingET);
      instance.SetDestructor(&destruct_TRootMissingET);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootMissingET*)
   {
      return GenerateInitInstanceLocal((::TRootMissingET*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootMissingET*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootPhoton(void *p = 0);
   static void *newArray_TRootPhoton(Long_t size, void *p);
   static void delete_TRootPhoton(void *p);
   static void deleteArray_TRootPhoton(void *p);
   static void destruct_TRootPhoton(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootPhoton*)
   {
      ::TRootPhoton *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootPhoton >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootPhoton", ::TRootPhoton::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 212,
                  typeid(::TRootPhoton), DefineBehavior(ptr, ptr),
                  &::TRootPhoton::Dictionary, isa_proxy, 4,
                  sizeof(::TRootPhoton) );
      instance.SetNew(&new_TRootPhoton);
      instance.SetNewArray(&newArray_TRootPhoton);
      instance.SetDelete(&delete_TRootPhoton);
      instance.SetDeleteArray(&deleteArray_TRootPhoton);
      instance.SetDestructor(&destruct_TRootPhoton);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootPhoton*)
   {
      return GenerateInitInstanceLocal((::TRootPhoton*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootPhoton*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootElectron(void *p = 0);
   static void *newArray_TRootElectron(Long_t size, void *p);
   static void delete_TRootElectron(void *p);
   static void deleteArray_TRootElectron(void *p);
   static void destruct_TRootElectron(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootElectron*)
   {
      ::TRootElectron *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootElectron >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootElectron", ::TRootElectron::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 231,
                  typeid(::TRootElectron), DefineBehavior(ptr, ptr),
                  &::TRootElectron::Dictionary, isa_proxy, 4,
                  sizeof(::TRootElectron) );
      instance.SetNew(&new_TRootElectron);
      instance.SetNewArray(&newArray_TRootElectron);
      instance.SetDelete(&delete_TRootElectron);
      instance.SetDeleteArray(&deleteArray_TRootElectron);
      instance.SetDestructor(&destruct_TRootElectron);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootElectron*)
   {
      return GenerateInitInstanceLocal((::TRootElectron*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootElectron*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootMuon(void *p = 0);
   static void *newArray_TRootMuon(Long_t size, void *p);
   static void delete_TRootMuon(void *p);
   static void deleteArray_TRootMuon(void *p);
   static void destruct_TRootMuon(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootMuon*)
   {
      ::TRootMuon *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootMuon >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootMuon", ::TRootMuon::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 254,
                  typeid(::TRootMuon), DefineBehavior(ptr, ptr),
                  &::TRootMuon::Dictionary, isa_proxy, 4,
                  sizeof(::TRootMuon) );
      instance.SetNew(&new_TRootMuon);
      instance.SetNewArray(&newArray_TRootMuon);
      instance.SetDelete(&delete_TRootMuon);
      instance.SetDeleteArray(&deleteArray_TRootMuon);
      instance.SetDestructor(&destruct_TRootMuon);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootMuon*)
   {
      return GenerateInitInstanceLocal((::TRootMuon*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootMuon*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootTau(void *p = 0);
   static void *newArray_TRootTau(Long_t size, void *p);
   static void delete_TRootTau(void *p);
   static void deleteArray_TRootTau(void *p);
   static void destruct_TRootTau(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootTau*)
   {
      ::TRootTau *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootTau >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootTau", ::TRootTau::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 282,
                  typeid(::TRootTau), DefineBehavior(ptr, ptr),
                  &::TRootTau::Dictionary, isa_proxy, 4,
                  sizeof(::TRootTau) );
      instance.SetNew(&new_TRootTau);
      instance.SetNewArray(&newArray_TRootTau);
      instance.SetDelete(&delete_TRootTau);
      instance.SetDeleteArray(&deleteArray_TRootTau);
      instance.SetDestructor(&destruct_TRootTau);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootTau*)
   {
      return GenerateInitInstanceLocal((::TRootTau*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootTau*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRootJet(void *p = 0);
   static void *newArray_TRootJet(Long_t size, void *p);
   static void delete_TRootJet(void *p);
   static void deleteArray_TRootJet(void *p);
   static void destruct_TRootJet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootJet*)
   {
      ::TRootJet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootJet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRootJet", ::TRootJet::Class_Version(), "ExRootAnalysis/ExRootClasses.h", 305,
                  typeid(::TRootJet), DefineBehavior(ptr, ptr),
                  &::TRootJet::Dictionary, isa_proxy, 4,
                  sizeof(::TRootJet) );
      instance.SetNew(&new_TRootJet);
      instance.SetNewArray(&newArray_TRootJet);
      instance.SetDelete(&delete_TRootJet);
      instance.SetDeleteArray(&deleteArray_TRootJet);
      instance.SetDestructor(&destruct_TRootJet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootJet*)
   {
      return GenerateInitInstanceLocal((::TRootJet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRootJet*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ExRootTreeReader(void *p = 0);
   static void *newArray_ExRootTreeReader(Long_t size, void *p);
   static void delete_ExRootTreeReader(void *p);
   static void deleteArray_ExRootTreeReader(void *p);
   static void destruct_ExRootTreeReader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootTreeReader*)
   {
      ::ExRootTreeReader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExRootTreeReader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ExRootTreeReader", ::ExRootTreeReader::Class_Version(), "ExRootAnalysis/ExRootTreeReader.h", 26,
                  typeid(::ExRootTreeReader), DefineBehavior(ptr, ptr),
                  &::ExRootTreeReader::Dictionary, isa_proxy, 4,
                  sizeof(::ExRootTreeReader) );
      instance.SetNew(&new_ExRootTreeReader);
      instance.SetNewArray(&newArray_ExRootTreeReader);
      instance.SetDelete(&delete_ExRootTreeReader);
      instance.SetDeleteArray(&deleteArray_ExRootTreeReader);
      instance.SetDestructor(&destruct_ExRootTreeReader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootTreeReader*)
   {
      return GenerateInitInstanceLocal((::ExRootTreeReader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExRootTreeReader*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootTreeBranch_Dictionary();
   static void ExRootTreeBranch_TClassManip(TClass*);
   static void delete_ExRootTreeBranch(void *p);
   static void deleteArray_ExRootTreeBranch(void *p);
   static void destruct_ExRootTreeBranch(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootTreeBranch*)
   {
      ::ExRootTreeBranch *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootTreeBranch));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootTreeBranch", "ExRootAnalysis/ExRootTreeBranch.h", 22,
                  typeid(::ExRootTreeBranch), DefineBehavior(ptr, ptr),
                  &ExRootTreeBranch_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootTreeBranch) );
      instance.SetDelete(&delete_ExRootTreeBranch);
      instance.SetDeleteArray(&deleteArray_ExRootTreeBranch);
      instance.SetDestructor(&destruct_ExRootTreeBranch);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootTreeBranch*)
   {
      return GenerateInitInstanceLocal((::ExRootTreeBranch*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExRootTreeBranch*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootTreeBranch_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeBranch*)0x0)->GetClass();
      ExRootTreeBranch_TClassManip(theClass);
   return theClass;
   }

   static void ExRootTreeBranch_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_ExRootTreeWriter(void *p = 0);
   static void *newArray_ExRootTreeWriter(Long_t size, void *p);
   static void delete_ExRootTreeWriter(void *p);
   static void deleteArray_ExRootTreeWriter(void *p);
   static void destruct_ExRootTreeWriter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootTreeWriter*)
   {
      ::ExRootTreeWriter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExRootTreeWriter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ExRootTreeWriter", ::ExRootTreeWriter::Class_Version(), "ExRootAnalysis/ExRootTreeWriter.h", 25,
                  typeid(::ExRootTreeWriter), DefineBehavior(ptr, ptr),
                  &::ExRootTreeWriter::Dictionary, isa_proxy, 4,
                  sizeof(::ExRootTreeWriter) );
      instance.SetNew(&new_ExRootTreeWriter);
      instance.SetNewArray(&newArray_ExRootTreeWriter);
      instance.SetDelete(&delete_ExRootTreeWriter);
      instance.SetDeleteArray(&deleteArray_ExRootTreeWriter);
      instance.SetDestructor(&destruct_ExRootTreeWriter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootTreeWriter*)
   {
      return GenerateInitInstanceLocal((::ExRootTreeWriter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExRootTreeWriter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootResult_Dictionary();
   static void ExRootResult_TClassManip(TClass*);
   static void *new_ExRootResult(void *p = 0);
   static void *newArray_ExRootResult(Long_t size, void *p);
   static void delete_ExRootResult(void *p);
   static void deleteArray_ExRootResult(void *p);
   static void destruct_ExRootResult(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootResult*)
   {
      ::ExRootResult *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootResult));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootResult", "ExRootAnalysis/ExRootResult.h", 21,
                  typeid(::ExRootResult), DefineBehavior(ptr, ptr),
                  &ExRootResult_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootResult) );
      instance.SetNew(&new_ExRootResult);
      instance.SetNewArray(&newArray_ExRootResult);
      instance.SetDelete(&delete_ExRootResult);
      instance.SetDeleteArray(&deleteArray_ExRootResult);
      instance.SetDestructor(&destruct_ExRootResult);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootResult*)
   {
      return GenerateInitInstanceLocal((::ExRootResult*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExRootResult*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootResult_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ExRootResult*)0x0)->GetClass();
      ExRootResult_TClassManip(theClass);
   return theClass;
   }

   static void ExRootResult_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootClassifier_Dictionary();
   static void ExRootClassifier_TClassManip(TClass*);
   static void delete_ExRootClassifier(void *p);
   static void deleteArray_ExRootClassifier(void *p);
   static void destruct_ExRootClassifier(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootClassifier*)
   {
      ::ExRootClassifier *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootClassifier));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootClassifier", "ExRootAnalysis/ExRootClassifier.h", 8,
                  typeid(::ExRootClassifier), DefineBehavior(ptr, ptr),
                  &ExRootClassifier_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootClassifier) );
      instance.SetDelete(&delete_ExRootClassifier);
      instance.SetDeleteArray(&deleteArray_ExRootClassifier);
      instance.SetDestructor(&destruct_ExRootClassifier);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootClassifier*)
   {
      return GenerateInitInstanceLocal((::ExRootClassifier*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExRootClassifier*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootClassifier_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ExRootClassifier*)0x0)->GetClass();
      ExRootClassifier_TClassManip(theClass);
   return theClass;
   }

   static void ExRootClassifier_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ExRootFilter_Dictionary();
   static void ExRootFilter_TClassManip(TClass*);
   static void delete_ExRootFilter(void *p);
   static void deleteArray_ExRootFilter(void *p);
   static void destruct_ExRootFilter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootFilter*)
   {
      ::ExRootFilter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExRootFilter));
      static ::ROOT::TGenericClassInfo 
         instance("ExRootFilter", "ExRootAnalysis/ExRootFilter.h", 13,
                  typeid(::ExRootFilter), DefineBehavior(ptr, ptr),
                  &ExRootFilter_Dictionary, isa_proxy, 4,
                  sizeof(::ExRootFilter) );
      instance.SetDelete(&delete_ExRootFilter);
      instance.SetDeleteArray(&deleteArray_ExRootFilter);
      instance.SetDestructor(&destruct_ExRootFilter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootFilter*)
   {
      return GenerateInitInstanceLocal((::ExRootFilter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExRootFilter*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExRootFilter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ExRootFilter*)0x0)->GetClass();
      ExRootFilter_TClassManip(theClass);
   return theClass;
   }

   static void ExRootFilter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_ExRootFactory(void *p = 0);
   static void *newArray_ExRootFactory(Long_t size, void *p);
   static void delete_ExRootFactory(void *p);
   static void deleteArray_ExRootFactory(void *p);
   static void destruct_ExRootFactory(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExRootFactory*)
   {
      ::ExRootFactory *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExRootFactory >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ExRootFactory", ::ExRootFactory::Class_Version(), "ExRootAnalysis/ExRootFactory.h", 27,
                  typeid(::ExRootFactory), DefineBehavior(ptr, ptr),
                  &::ExRootFactory::Dictionary, isa_proxy, 4,
                  sizeof(::ExRootFactory) );
      instance.SetNew(&new_ExRootFactory);
      instance.SetNewArray(&newArray_ExRootFactory);
      instance.SetDelete(&delete_ExRootFactory);
      instance.SetDeleteArray(&deleteArray_ExRootFactory);
      instance.SetDestructor(&destruct_ExRootFactory);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExRootFactory*)
   {
      return GenerateInitInstanceLocal((::ExRootFactory*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExRootFactory*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TSortableObject::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TSortableObject::Class_Name()
{
   return "TSortableObject";
}

//______________________________________________________________________________
const char *TSortableObject::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSortableObject*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TSortableObject::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSortableObject*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TSortableObject::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSortableObject*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TSortableObject::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSortableObject*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootWeight::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootWeight::Class_Name()
{
   return "TRootWeight";
}

//______________________________________________________________________________
const char *TRootWeight::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootWeight*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootWeight::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootWeight*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootWeight::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootWeight*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootWeight::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootWeight*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootLHEFEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootLHEFEvent::Class_Name()
{
   return "TRootLHEFEvent";
}

//______________________________________________________________________________
const char *TRootLHEFEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootLHEFEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootLHEFEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootLHEFEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootLHEFEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootLHEFEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootLHEFEvent::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootLHEFEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootLHEFParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootLHEFParticle::Class_Name()
{
   return "TRootLHEFParticle";
}

//______________________________________________________________________________
const char *TRootLHEFParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootLHEFParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootLHEFParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootLHEFParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootLHEFParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootLHEFParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootLHEFParticle::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootLHEFParticle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootGenEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootGenEvent::Class_Name()
{
   return "TRootGenEvent";
}

//______________________________________________________________________________
const char *TRootGenEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootGenEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootGenEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootGenEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootGenEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootGenEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootGenEvent::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootGenEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootGenParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootGenParticle::Class_Name()
{
   return "TRootGenParticle";
}

//______________________________________________________________________________
const char *TRootGenParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootGenParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootGenParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootGenParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootGenParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootGenParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootGenParticle::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootGenParticle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootGenJet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootGenJet::Class_Name()
{
   return "TRootGenJet";
}

//______________________________________________________________________________
const char *TRootGenJet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootGenJet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootGenJet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootGenJet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootGenJet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootGenJet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootGenJet::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootGenJet*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootEvent::Class_Name()
{
   return "TRootEvent";
}

//______________________________________________________________________________
const char *TRootEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootEvent::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootMissingET::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootMissingET::Class_Name()
{
   return "TRootMissingET";
}

//______________________________________________________________________________
const char *TRootMissingET::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootMissingET*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootMissingET::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootMissingET*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootMissingET::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootMissingET*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootMissingET::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootMissingET*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootPhoton::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootPhoton::Class_Name()
{
   return "TRootPhoton";
}

//______________________________________________________________________________
const char *TRootPhoton::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootPhoton*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootPhoton::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootPhoton*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootPhoton::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootPhoton*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootPhoton::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootPhoton*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootElectron::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootElectron::Class_Name()
{
   return "TRootElectron";
}

//______________________________________________________________________________
const char *TRootElectron::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootElectron*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootElectron::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootElectron*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootElectron::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootElectron*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootElectron::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootElectron*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootMuon::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootMuon::Class_Name()
{
   return "TRootMuon";
}

//______________________________________________________________________________
const char *TRootMuon::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootMuon*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootMuon::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootMuon*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootMuon::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootMuon*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootMuon::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootMuon*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootTau::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootTau::Class_Name()
{
   return "TRootTau";
}

//______________________________________________________________________________
const char *TRootTau::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootTau*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootTau::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootTau*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootTau::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootTau*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootTau::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootTau*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRootJet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootJet::Class_Name()
{
   return "TRootJet";
}

//______________________________________________________________________________
const char *TRootJet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootJet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRootJet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootJet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootJet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootJet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootJet::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootJet*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ExRootTreeReader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ExRootTreeReader::Class_Name()
{
   return "ExRootTreeReader";
}

//______________________________________________________________________________
const char *ExRootTreeReader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeReader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ExRootTreeReader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeReader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ExRootTreeReader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeReader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ExRootTreeReader::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeReader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ExRootTreeWriter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ExRootTreeWriter::Class_Name()
{
   return "ExRootTreeWriter";
}

//______________________________________________________________________________
const char *ExRootTreeWriter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeWriter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ExRootTreeWriter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeWriter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ExRootTreeWriter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeWriter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ExRootTreeWriter::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootTreeWriter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ExRootFactory::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ExRootFactory::Class_Name()
{
   return "ExRootFactory";
}

//______________________________________________________________________________
const char *ExRootFactory::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootFactory*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ExRootFactory::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExRootFactory*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ExRootFactory::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootFactory*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ExRootFactory::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExRootFactory*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TSortableObject::Streamer(TBuffer &R__b)
{
   // Stream an object of class TSortableObject.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TSortableObject::Class(),this);
   } else {
      R__b.WriteClassBuffer(TSortableObject::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TSortableObject(void *p) {
      delete ((::TSortableObject*)p);
   }
   static void deleteArray_TSortableObject(void *p) {
      delete [] ((::TSortableObject*)p);
   }
   static void destruct_TSortableObject(void *p) {
      typedef ::TSortableObject current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TSortableObject

//______________________________________________________________________________
void TRootWeight::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootWeight.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> Weight;
      R__b.CheckByteCount(R__s, R__c, TRootWeight::IsA());
   } else {
      R__c = R__b.WriteVersion(TRootWeight::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << Weight;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootWeight(void *p) {
      return  p ? new(p) ::TRootWeight : new ::TRootWeight;
   }
   static void *newArray_TRootWeight(Long_t nElements, void *p) {
      return p ? new(p) ::TRootWeight[nElements] : new ::TRootWeight[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootWeight(void *p) {
      delete ((::TRootWeight*)p);
   }
   static void deleteArray_TRootWeight(void *p) {
      delete [] ((::TRootWeight*)p);
   }
   static void destruct_TRootWeight(void *p) {
      typedef ::TRootWeight current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TRootWeight(TBuffer &buf, void *obj) {
      ((::TRootWeight*)obj)->::TRootWeight::Streamer(buf);
   }
} // end of namespace ROOT for class ::TRootWeight

//______________________________________________________________________________
void TRootLHEFEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootLHEFEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootLHEFEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootLHEFEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootLHEFEvent(void *p) {
      return  p ? new(p) ::TRootLHEFEvent : new ::TRootLHEFEvent;
   }
   static void *newArray_TRootLHEFEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TRootLHEFEvent[nElements] : new ::TRootLHEFEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootLHEFEvent(void *p) {
      delete ((::TRootLHEFEvent*)p);
   }
   static void deleteArray_TRootLHEFEvent(void *p) {
      delete [] ((::TRootLHEFEvent*)p);
   }
   static void destruct_TRootLHEFEvent(void *p) {
      typedef ::TRootLHEFEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootLHEFEvent

//______________________________________________________________________________
void TRootLHEFParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootLHEFParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootLHEFParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootLHEFParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootLHEFParticle(void *p) {
      return  p ? new(p) ::TRootLHEFParticle : new ::TRootLHEFParticle;
   }
   static void *newArray_TRootLHEFParticle(Long_t nElements, void *p) {
      return p ? new(p) ::TRootLHEFParticle[nElements] : new ::TRootLHEFParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootLHEFParticle(void *p) {
      delete ((::TRootLHEFParticle*)p);
   }
   static void deleteArray_TRootLHEFParticle(void *p) {
      delete [] ((::TRootLHEFParticle*)p);
   }
   static void destruct_TRootLHEFParticle(void *p) {
      typedef ::TRootLHEFParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootLHEFParticle

//______________________________________________________________________________
void TRootGenEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootGenEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootGenEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootGenEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootGenEvent(void *p) {
      return  p ? new(p) ::TRootGenEvent : new ::TRootGenEvent;
   }
   static void *newArray_TRootGenEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TRootGenEvent[nElements] : new ::TRootGenEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootGenEvent(void *p) {
      delete ((::TRootGenEvent*)p);
   }
   static void deleteArray_TRootGenEvent(void *p) {
      delete [] ((::TRootGenEvent*)p);
   }
   static void destruct_TRootGenEvent(void *p) {
      typedef ::TRootGenEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootGenEvent

//______________________________________________________________________________
void TRootGenParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootGenParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootGenParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootGenParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootGenParticle(void *p) {
      return  p ? new(p) ::TRootGenParticle : new ::TRootGenParticle;
   }
   static void *newArray_TRootGenParticle(Long_t nElements, void *p) {
      return p ? new(p) ::TRootGenParticle[nElements] : new ::TRootGenParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootGenParticle(void *p) {
      delete ((::TRootGenParticle*)p);
   }
   static void deleteArray_TRootGenParticle(void *p) {
      delete [] ((::TRootGenParticle*)p);
   }
   static void destruct_TRootGenParticle(void *p) {
      typedef ::TRootGenParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootGenParticle

//______________________________________________________________________________
void TRootGenJet::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootGenJet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootGenJet::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootGenJet::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootGenJet(void *p) {
      return  p ? new(p) ::TRootGenJet : new ::TRootGenJet;
   }
   static void *newArray_TRootGenJet(Long_t nElements, void *p) {
      return p ? new(p) ::TRootGenJet[nElements] : new ::TRootGenJet[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootGenJet(void *p) {
      delete ((::TRootGenJet*)p);
   }
   static void deleteArray_TRootGenJet(void *p) {
      delete [] ((::TRootGenJet*)p);
   }
   static void destruct_TRootGenJet(void *p) {
      typedef ::TRootGenJet current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootGenJet

//______________________________________________________________________________
void TRootEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootEvent(void *p) {
      return  p ? new(p) ::TRootEvent : new ::TRootEvent;
   }
   static void *newArray_TRootEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TRootEvent[nElements] : new ::TRootEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootEvent(void *p) {
      delete ((::TRootEvent*)p);
   }
   static void deleteArray_TRootEvent(void *p) {
      delete [] ((::TRootEvent*)p);
   }
   static void destruct_TRootEvent(void *p) {
      typedef ::TRootEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootEvent

//______________________________________________________________________________
void TRootMissingET::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootMissingET.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootMissingET::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootMissingET::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootMissingET(void *p) {
      return  p ? new(p) ::TRootMissingET : new ::TRootMissingET;
   }
   static void *newArray_TRootMissingET(Long_t nElements, void *p) {
      return p ? new(p) ::TRootMissingET[nElements] : new ::TRootMissingET[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootMissingET(void *p) {
      delete ((::TRootMissingET*)p);
   }
   static void deleteArray_TRootMissingET(void *p) {
      delete [] ((::TRootMissingET*)p);
   }
   static void destruct_TRootMissingET(void *p) {
      typedef ::TRootMissingET current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootMissingET

//______________________________________________________________________________
void TRootPhoton::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootPhoton.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootPhoton::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootPhoton::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootPhoton(void *p) {
      return  p ? new(p) ::TRootPhoton : new ::TRootPhoton;
   }
   static void *newArray_TRootPhoton(Long_t nElements, void *p) {
      return p ? new(p) ::TRootPhoton[nElements] : new ::TRootPhoton[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootPhoton(void *p) {
      delete ((::TRootPhoton*)p);
   }
   static void deleteArray_TRootPhoton(void *p) {
      delete [] ((::TRootPhoton*)p);
   }
   static void destruct_TRootPhoton(void *p) {
      typedef ::TRootPhoton current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootPhoton

//______________________________________________________________________________
void TRootElectron::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootElectron.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootElectron::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootElectron::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootElectron(void *p) {
      return  p ? new(p) ::TRootElectron : new ::TRootElectron;
   }
   static void *newArray_TRootElectron(Long_t nElements, void *p) {
      return p ? new(p) ::TRootElectron[nElements] : new ::TRootElectron[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootElectron(void *p) {
      delete ((::TRootElectron*)p);
   }
   static void deleteArray_TRootElectron(void *p) {
      delete [] ((::TRootElectron*)p);
   }
   static void destruct_TRootElectron(void *p) {
      typedef ::TRootElectron current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootElectron

//______________________________________________________________________________
void TRootMuon::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootMuon.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootMuon::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootMuon::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootMuon(void *p) {
      return  p ? new(p) ::TRootMuon : new ::TRootMuon;
   }
   static void *newArray_TRootMuon(Long_t nElements, void *p) {
      return p ? new(p) ::TRootMuon[nElements] : new ::TRootMuon[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootMuon(void *p) {
      delete ((::TRootMuon*)p);
   }
   static void deleteArray_TRootMuon(void *p) {
      delete [] ((::TRootMuon*)p);
   }
   static void destruct_TRootMuon(void *p) {
      typedef ::TRootMuon current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootMuon

//______________________________________________________________________________
void TRootTau::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootTau.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootTau::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootTau::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootTau(void *p) {
      return  p ? new(p) ::TRootTau : new ::TRootTau;
   }
   static void *newArray_TRootTau(Long_t nElements, void *p) {
      return p ? new(p) ::TRootTau[nElements] : new ::TRootTau[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootTau(void *p) {
      delete ((::TRootTau*)p);
   }
   static void deleteArray_TRootTau(void *p) {
      delete [] ((::TRootTau*)p);
   }
   static void destruct_TRootTau(void *p) {
      typedef ::TRootTau current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootTau

//______________________________________________________________________________
void TRootJet::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootJet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TRootJet::Class(),this);
   } else {
      R__b.WriteClassBuffer(TRootJet::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootJet(void *p) {
      return  p ? new(p) ::TRootJet : new ::TRootJet;
   }
   static void *newArray_TRootJet(Long_t nElements, void *p) {
      return p ? new(p) ::TRootJet[nElements] : new ::TRootJet[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootJet(void *p) {
      delete ((::TRootJet*)p);
   }
   static void deleteArray_TRootJet(void *p) {
      delete [] ((::TRootJet*)p);
   }
   static void destruct_TRootJet(void *p) {
      typedef ::TRootJet current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRootJet

//______________________________________________________________________________
void ExRootTreeReader::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExRootTreeReader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ExRootTreeReader::Class(),this);
   } else {
      R__b.WriteClassBuffer(ExRootTreeReader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootTreeReader(void *p) {
      return  p ? new(p) ::ExRootTreeReader : new ::ExRootTreeReader;
   }
   static void *newArray_ExRootTreeReader(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootTreeReader[nElements] : new ::ExRootTreeReader[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootTreeReader(void *p) {
      delete ((::ExRootTreeReader*)p);
   }
   static void deleteArray_ExRootTreeReader(void *p) {
      delete [] ((::ExRootTreeReader*)p);
   }
   static void destruct_ExRootTreeReader(void *p) {
      typedef ::ExRootTreeReader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExRootTreeReader

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ExRootTreeBranch(void *p) {
      delete ((::ExRootTreeBranch*)p);
   }
   static void deleteArray_ExRootTreeBranch(void *p) {
      delete [] ((::ExRootTreeBranch*)p);
   }
   static void destruct_ExRootTreeBranch(void *p) {
      typedef ::ExRootTreeBranch current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExRootTreeBranch

//______________________________________________________________________________
void ExRootTreeWriter::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExRootTreeWriter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ExRootTreeWriter::Class(),this);
   } else {
      R__b.WriteClassBuffer(ExRootTreeWriter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootTreeWriter(void *p) {
      return  p ? new(p) ::ExRootTreeWriter : new ::ExRootTreeWriter;
   }
   static void *newArray_ExRootTreeWriter(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootTreeWriter[nElements] : new ::ExRootTreeWriter[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootTreeWriter(void *p) {
      delete ((::ExRootTreeWriter*)p);
   }
   static void deleteArray_ExRootTreeWriter(void *p) {
      delete [] ((::ExRootTreeWriter*)p);
   }
   static void destruct_ExRootTreeWriter(void *p) {
      typedef ::ExRootTreeWriter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExRootTreeWriter

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootResult(void *p) {
      return  p ? new(p) ::ExRootResult : new ::ExRootResult;
   }
   static void *newArray_ExRootResult(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootResult[nElements] : new ::ExRootResult[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootResult(void *p) {
      delete ((::ExRootResult*)p);
   }
   static void deleteArray_ExRootResult(void *p) {
      delete [] ((::ExRootResult*)p);
   }
   static void destruct_ExRootResult(void *p) {
      typedef ::ExRootResult current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExRootResult

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ExRootClassifier(void *p) {
      delete ((::ExRootClassifier*)p);
   }
   static void deleteArray_ExRootClassifier(void *p) {
      delete [] ((::ExRootClassifier*)p);
   }
   static void destruct_ExRootClassifier(void *p) {
      typedef ::ExRootClassifier current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExRootClassifier

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ExRootFilter(void *p) {
      delete ((::ExRootFilter*)p);
   }
   static void deleteArray_ExRootFilter(void *p) {
      delete [] ((::ExRootFilter*)p);
   }
   static void destruct_ExRootFilter(void *p) {
      typedef ::ExRootFilter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExRootFilter

//______________________________________________________________________________
void ExRootFactory::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExRootFactory.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ExRootFactory::Class(),this);
   } else {
      R__b.WriteClassBuffer(ExRootFactory::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExRootFactory(void *p) {
      return  p ? new(p) ::ExRootFactory : new ::ExRootFactory;
   }
   static void *newArray_ExRootFactory(Long_t nElements, void *p) {
      return p ? new(p) ::ExRootFactory[nElements] : new ::ExRootFactory[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExRootFactory(void *p) {
      delete ((::ExRootFactory*)p);
   }
   static void deleteArray_ExRootFactory(void *p) {
      delete [] ((::ExRootFactory*)p);
   }
   static void destruct_ExRootFactory(void *p) {
      typedef ::ExRootFactory current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExRootFactory

namespace ROOT {
   static TClass *setlEExRootTreeBranchmUgR_Dictionary();
   static void setlEExRootTreeBranchmUgR_TClassManip(TClass*);
   static void *new_setlEExRootTreeBranchmUgR(void *p = 0);
   static void *newArray_setlEExRootTreeBranchmUgR(Long_t size, void *p);
   static void delete_setlEExRootTreeBranchmUgR(void *p);
   static void deleteArray_setlEExRootTreeBranchmUgR(void *p);
   static void destruct_setlEExRootTreeBranchmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<ExRootTreeBranch*>*)
   {
      set<ExRootTreeBranch*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<ExRootTreeBranch*>));
      static ::ROOT::TGenericClassInfo 
         instance("set<ExRootTreeBranch*>", -2, "set", 90,
                  typeid(set<ExRootTreeBranch*>), DefineBehavior(ptr, ptr),
                  &setlEExRootTreeBranchmUgR_Dictionary, isa_proxy, 0,
                  sizeof(set<ExRootTreeBranch*>) );
      instance.SetNew(&new_setlEExRootTreeBranchmUgR);
      instance.SetNewArray(&newArray_setlEExRootTreeBranchmUgR);
      instance.SetDelete(&delete_setlEExRootTreeBranchmUgR);
      instance.SetDeleteArray(&deleteArray_setlEExRootTreeBranchmUgR);
      instance.SetDestructor(&destruct_setlEExRootTreeBranchmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<ExRootTreeBranch*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const set<ExRootTreeBranch*>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlEExRootTreeBranchmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<ExRootTreeBranch*>*)0x0)->GetClass();
      setlEExRootTreeBranchmUgR_TClassManip(theClass);
   return theClass;
   }

   static void setlEExRootTreeBranchmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlEExRootTreeBranchmUgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) set<ExRootTreeBranch*> : new set<ExRootTreeBranch*>;
   }
   static void *newArray_setlEExRootTreeBranchmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) set<ExRootTreeBranch*>[nElements] : new set<ExRootTreeBranch*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlEExRootTreeBranchmUgR(void *p) {
      delete ((set<ExRootTreeBranch*>*)p);
   }
   static void deleteArray_setlEExRootTreeBranchmUgR(void *p) {
      delete [] ((set<ExRootTreeBranch*>*)p);
   }
   static void destruct_setlEExRootTreeBranchmUgR(void *p) {
      typedef set<ExRootTreeBranch*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<ExRootTreeBranch*>

namespace ROOT {
   static TClass *maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR_Dictionary();
   static void maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR_TClassManip(TClass*);
   static void *new_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(void *p = 0);
   static void *newArray_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(Long_t size, void *p);
   static void delete_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(void *p);
   static void deleteArray_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(void *p);
   static void destruct_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<TString,pair<TBranch*,TClonesArray*> >*)
   {
      map<TString,pair<TBranch*,TClonesArray*> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<TString,pair<TBranch*,TClonesArray*> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<TString,pair<TBranch*,TClonesArray*> >", -2, "map", 96,
                  typeid(map<TString,pair<TBranch*,TClonesArray*> >), DefineBehavior(ptr, ptr),
                  &maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<TString,pair<TBranch*,TClonesArray*> >) );
      instance.SetNew(&new_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR);
      instance.SetNewArray(&newArray_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR);
      instance.SetDelete(&delete_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR);
      instance.SetDestructor(&destruct_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<TString,pair<TBranch*,TClonesArray*> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<TString,pair<TBranch*,TClonesArray*> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<TString,pair<TBranch*,TClonesArray*> >*)0x0)->GetClass();
      maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) map<TString,pair<TBranch*,TClonesArray*> > : new map<TString,pair<TBranch*,TClonesArray*> >;
   }
   static void *newArray_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) map<TString,pair<TBranch*,TClonesArray*> >[nElements] : new map<TString,pair<TBranch*,TClonesArray*> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(void *p) {
      delete ((map<TString,pair<TBranch*,TClonesArray*> >*)p);
   }
   static void deleteArray_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(void *p) {
      delete [] ((map<TString,pair<TBranch*,TClonesArray*> >*)p);
   }
   static void destruct_maplETStringcOpairlETBranchmUcOTClonesArraymUgRsPgR(void *p) {
      typedef map<TString,pair<TBranch*,TClonesArray*> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<TString,pair<TBranch*,TClonesArray*> >

namespace ROOT {
   static TClass *maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR_Dictionary();
   static void maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR_TClassManip(TClass*);
   static void *new_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(void *p = 0);
   static void *newArray_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(Long_t size, void *p);
   static void delete_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(void *p);
   static void deleteArray_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(void *p);
   static void destruct_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >*)
   {
      map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >));
      static ::ROOT::TGenericClassInfo 
         instance("map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >", -2, "map", 96,
                  typeid(map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >), DefineBehavior(ptr, ptr),
                  &maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >) );
      instance.SetNew(&new_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR);
      instance.SetNewArray(&newArray_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR);
      instance.SetDelete(&delete_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR);
      instance.SetDestructor(&destruct_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >*)0x0)->GetClass();
      maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > > : new map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >;
   }
   static void *newArray_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >[nElements] : new map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(void *p) {
      delete ((map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >*)p);
   }
   static void deleteArray_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(void *p) {
      delete [] ((map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >*)p);
   }
   static void destruct_maplEExRootClassifiermUcOpairlEboolcOmaplEintcOTObjArraymUgRsPgRsPgR(void *p) {
      typedef map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<ExRootClassifier*,pair<bool,map<int,TObjArray*> > >

namespace {
  void TriggerDictionaryInitialization_ExRootAnalysisDict_Impl() {
    static const char* headers[] = {
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.00-odfocd7/include",
"/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/Madgraph/MG5_aMC_v2_3_3/ExRootAnalysis/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TSortableObject;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootWeight;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootLHEFEvent;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootLHEFParticle;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootGenEvent;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootGenParticle;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootGenJet;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootEvent;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootMissingET;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootPhoton;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootElectron;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootMuon;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootTau;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  TRootJet;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  ExRootTreeReader;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  ExRootTreeBranch;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  ExRootTreeWriter;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  ExRootResult;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  ExRootClassifier;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  ExRootFilter;
class __attribute__((annotate("$clingAutoload$src/ExRootAnalysisLinkDef.h")))  ExRootFactory;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"ExRootClassifier", payloadCode, "@",
"ExRootFactory", payloadCode, "@",
"ExRootFilter", payloadCode, "@",
"ExRootResult", payloadCode, "@",
"ExRootTreeBranch", payloadCode, "@",
"ExRootTreeReader", payloadCode, "@",
"ExRootTreeWriter", payloadCode, "@",
"FillChain", payloadCode, "@",
"HistStyle", payloadCode, "@",
"TRootElectron", payloadCode, "@",
"TRootEvent", payloadCode, "@",
"TRootGenEvent", payloadCode, "@",
"TRootGenJet", payloadCode, "@",
"TRootGenParticle", payloadCode, "@",
"TRootJet", payloadCode, "@",
"TRootLHEFEvent", payloadCode, "@",
"TRootLHEFParticle", payloadCode, "@",
"TRootMissingET", payloadCode, "@",
"TRootMuon", payloadCode, "@",
"TRootPhoton", payloadCode, "@",
"TRootTau", payloadCode, "@",
"TRootWeight", payloadCode, "@",
"TSortableObject", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("ExRootAnalysisDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_ExRootAnalysisDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_ExRootAnalysisDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_ExRootAnalysisDict() {
  TriggerDictionaryInitialization_ExRootAnalysisDict_Impl();
}
