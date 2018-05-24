#define private public
#define protected public
/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Delphes.h"

#include "modules/AngularSmearing.h"
#include "modules/PhotonConversions.h"
#include "modules/ParticlePropagator.h"
#include "modules/Efficiency.h"
#include "modules/IdentificationMap.h"
#include "modules/EnergySmearing.h"
#include "modules/MomentumSmearing.h"
#include "modules/ImpactParameterSmearing.h"
#include "modules/TimeSmearing.h"
#include "modules/SimpleCalorimeter.h"
#include "modules/Calorimeter.h"
#include "modules/Isolation.h"
#include "modules/EnergyScale.h"
#include "modules/UniqueObjectFinder.h"
#include "modules/TrackCountingBTagging.h"
#include "modules/BTagging.h"
#include "modules/TauTagging.h"
#include "modules/TreeWriter.h"
#include "modules/Merger.h"
#include "modules/LeptonDressing.h"
#include "modules/PileUpMerger.h"
#include "modules/JetPileUpSubtractor.h"
#include "modules/TrackPileUpSubtractor.h"
#include "modules/TaggingParticlesSkimmer.h"
#include "modules/PileUpJetID.h"
#include "modules/ConstituentFilter.h"
#include "modules/StatusPidFilter.h"
#include "modules/PdgCodeFilter.h"
#include "modules/Cloner.h"
#include "modules/Weighter.h"
#include "modules/Hector.h"
#include "modules/JetFlavorAssociation.h"
#include "modules/JetFakeParticle.h"
#include "modules/ExampleModule.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Delphes+;

#pragma link C++ class AngularSmearing+;
#pragma link C++ class PhotonConversions+;
#pragma link C++ class ParticlePropagator+;
#pragma link C++ class Efficiency+;
#pragma link C++ class IdentificationMap+;
#pragma link C++ class EnergySmearing+;
#pragma link C++ class MomentumSmearing+;
#pragma link C++ class ImpactParameterSmearing+;
#pragma link C++ class TimeSmearing+;
#pragma link C++ class SimpleCalorimeter+;
#pragma link C++ class Calorimeter+;
#pragma link C++ class Isolation+;
#pragma link C++ class EnergyScale+;
#pragma link C++ class UniqueObjectFinder+;
#pragma link C++ class TrackCountingBTagging+;
#pragma link C++ class BTagging+;
#pragma link C++ class TauTagging+;
#pragma link C++ class TreeWriter+;
#pragma link C++ class Merger+;
#pragma link C++ class LeptonDressing+;
#pragma link C++ class PileUpMerger+;
#pragma link C++ class JetPileUpSubtractor+;
#pragma link C++ class TrackPileUpSubtractor+;
#pragma link C++ class TaggingParticlesSkimmer+;
#pragma link C++ class PileUpJetID+;
#pragma link C++ class ConstituentFilter+;
#pragma link C++ class StatusPidFilter+;
#pragma link C++ class PdgCodeFilter+;
#pragma link C++ class Cloner+;
#pragma link C++ class Weighter+;
#pragma link C++ class Hector+;
#pragma link C++ class JetFlavorAssociation+;
#pragma link C++ class JetFakeParticle+;
#pragma link C++ class ExampleModule+;

#endif
// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME tmpdImodulesdIModulesDict

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
   static void *new_Delphes(void *p = 0);
   static void *newArray_Delphes(Long_t size, void *p);
   static void delete_Delphes(void *p);
   static void deleteArray_Delphes(void *p);
   static void destruct_Delphes(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Delphes*)
   {
      ::Delphes *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Delphes >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Delphes", ::Delphes::Class_Version(), "modules/Delphes.h", 40,
                  typeid(::Delphes), DefineBehavior(ptr, ptr),
                  &::Delphes::Dictionary, isa_proxy, 4,
                  sizeof(::Delphes) );
      instance.SetNew(&new_Delphes);
      instance.SetNewArray(&newArray_Delphes);
      instance.SetDelete(&delete_Delphes);
      instance.SetDeleteArray(&deleteArray_Delphes);
      instance.SetDestructor(&destruct_Delphes);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Delphes*)
   {
      return GenerateInitInstanceLocal((::Delphes*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Delphes*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_AngularSmearing(void *p = 0);
   static void *newArray_AngularSmearing(Long_t size, void *p);
   static void delete_AngularSmearing(void *p);
   static void deleteArray_AngularSmearing(void *p);
   static void destruct_AngularSmearing(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AngularSmearing*)
   {
      ::AngularSmearing *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::AngularSmearing >(0);
      static ::ROOT::TGenericClassInfo 
         instance("AngularSmearing", ::AngularSmearing::Class_Version(), "modules/AngularSmearing.h", 36,
                  typeid(::AngularSmearing), DefineBehavior(ptr, ptr),
                  &::AngularSmearing::Dictionary, isa_proxy, 4,
                  sizeof(::AngularSmearing) );
      instance.SetNew(&new_AngularSmearing);
      instance.SetNewArray(&newArray_AngularSmearing);
      instance.SetDelete(&delete_AngularSmearing);
      instance.SetDeleteArray(&deleteArray_AngularSmearing);
      instance.SetDestructor(&destruct_AngularSmearing);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AngularSmearing*)
   {
      return GenerateInitInstanceLocal((::AngularSmearing*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::AngularSmearing*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PhotonConversions(void *p = 0);
   static void *newArray_PhotonConversions(Long_t size, void *p);
   static void delete_PhotonConversions(void *p);
   static void deleteArray_PhotonConversions(void *p);
   static void destruct_PhotonConversions(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PhotonConversions*)
   {
      ::PhotonConversions *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PhotonConversions >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PhotonConversions", ::PhotonConversions::Class_Version(), "modules/PhotonConversions.h", 37,
                  typeid(::PhotonConversions), DefineBehavior(ptr, ptr),
                  &::PhotonConversions::Dictionary, isa_proxy, 4,
                  sizeof(::PhotonConversions) );
      instance.SetNew(&new_PhotonConversions);
      instance.SetNewArray(&newArray_PhotonConversions);
      instance.SetDelete(&delete_PhotonConversions);
      instance.SetDeleteArray(&deleteArray_PhotonConversions);
      instance.SetDestructor(&destruct_PhotonConversions);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PhotonConversions*)
   {
      return GenerateInitInstanceLocal((::PhotonConversions*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::PhotonConversions*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ParticlePropagator(void *p = 0);
   static void *newArray_ParticlePropagator(Long_t size, void *p);
   static void delete_ParticlePropagator(void *p);
   static void deleteArray_ParticlePropagator(void *p);
   static void destruct_ParticlePropagator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ParticlePropagator*)
   {
      ::ParticlePropagator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ParticlePropagator >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ParticlePropagator", ::ParticlePropagator::Class_Version(), "modules/ParticlePropagator.h", 38,
                  typeid(::ParticlePropagator), DefineBehavior(ptr, ptr),
                  &::ParticlePropagator::Dictionary, isa_proxy, 4,
                  sizeof(::ParticlePropagator) );
      instance.SetNew(&new_ParticlePropagator);
      instance.SetNewArray(&newArray_ParticlePropagator);
      instance.SetDelete(&delete_ParticlePropagator);
      instance.SetDeleteArray(&deleteArray_ParticlePropagator);
      instance.SetDestructor(&destruct_ParticlePropagator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ParticlePropagator*)
   {
      return GenerateInitInstanceLocal((::ParticlePropagator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ParticlePropagator*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Efficiency(void *p = 0);
   static void *newArray_Efficiency(Long_t size, void *p);
   static void delete_Efficiency(void *p);
   static void deleteArray_Efficiency(void *p);
   static void destruct_Efficiency(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Efficiency*)
   {
      ::Efficiency *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Efficiency >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Efficiency", ::Efficiency::Class_Version(), "modules/Efficiency.h", 36,
                  typeid(::Efficiency), DefineBehavior(ptr, ptr),
                  &::Efficiency::Dictionary, isa_proxy, 4,
                  sizeof(::Efficiency) );
      instance.SetNew(&new_Efficiency);
      instance.SetNewArray(&newArray_Efficiency);
      instance.SetDelete(&delete_Efficiency);
      instance.SetDeleteArray(&deleteArray_Efficiency);
      instance.SetDestructor(&destruct_Efficiency);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Efficiency*)
   {
      return GenerateInitInstanceLocal((::Efficiency*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Efficiency*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_IdentificationMap(void *p = 0);
   static void *newArray_IdentificationMap(Long_t size, void *p);
   static void delete_IdentificationMap(void *p);
   static void deleteArray_IdentificationMap(void *p);
   static void destruct_IdentificationMap(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::IdentificationMap*)
   {
      ::IdentificationMap *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::IdentificationMap >(0);
      static ::ROOT::TGenericClassInfo 
         instance("IdentificationMap", ::IdentificationMap::Class_Version(), "modules/IdentificationMap.h", 38,
                  typeid(::IdentificationMap), DefineBehavior(ptr, ptr),
                  &::IdentificationMap::Dictionary, isa_proxy, 4,
                  sizeof(::IdentificationMap) );
      instance.SetNew(&new_IdentificationMap);
      instance.SetNewArray(&newArray_IdentificationMap);
      instance.SetDelete(&delete_IdentificationMap);
      instance.SetDeleteArray(&deleteArray_IdentificationMap);
      instance.SetDestructor(&destruct_IdentificationMap);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::IdentificationMap*)
   {
      return GenerateInitInstanceLocal((::IdentificationMap*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::IdentificationMap*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_EnergySmearing(void *p = 0);
   static void *newArray_EnergySmearing(Long_t size, void *p);
   static void delete_EnergySmearing(void *p);
   static void deleteArray_EnergySmearing(void *p);
   static void destruct_EnergySmearing(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EnergySmearing*)
   {
      ::EnergySmearing *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EnergySmearing >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EnergySmearing", ::EnergySmearing::Class_Version(), "modules/EnergySmearing.h", 36,
                  typeid(::EnergySmearing), DefineBehavior(ptr, ptr),
                  &::EnergySmearing::Dictionary, isa_proxy, 4,
                  sizeof(::EnergySmearing) );
      instance.SetNew(&new_EnergySmearing);
      instance.SetNewArray(&newArray_EnergySmearing);
      instance.SetDelete(&delete_EnergySmearing);
      instance.SetDeleteArray(&deleteArray_EnergySmearing);
      instance.SetDestructor(&destruct_EnergySmearing);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EnergySmearing*)
   {
      return GenerateInitInstanceLocal((::EnergySmearing*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::EnergySmearing*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MomentumSmearing(void *p = 0);
   static void *newArray_MomentumSmearing(Long_t size, void *p);
   static void delete_MomentumSmearing(void *p);
   static void deleteArray_MomentumSmearing(void *p);
   static void destruct_MomentumSmearing(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MomentumSmearing*)
   {
      ::MomentumSmearing *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MomentumSmearing >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MomentumSmearing", ::MomentumSmearing::Class_Version(), "modules/MomentumSmearing.h", 36,
                  typeid(::MomentumSmearing), DefineBehavior(ptr, ptr),
                  &::MomentumSmearing::Dictionary, isa_proxy, 4,
                  sizeof(::MomentumSmearing) );
      instance.SetNew(&new_MomentumSmearing);
      instance.SetNewArray(&newArray_MomentumSmearing);
      instance.SetDelete(&delete_MomentumSmearing);
      instance.SetDeleteArray(&deleteArray_MomentumSmearing);
      instance.SetDestructor(&destruct_MomentumSmearing);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MomentumSmearing*)
   {
      return GenerateInitInstanceLocal((::MomentumSmearing*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MomentumSmearing*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ImpactParameterSmearing(void *p = 0);
   static void *newArray_ImpactParameterSmearing(Long_t size, void *p);
   static void delete_ImpactParameterSmearing(void *p);
   static void deleteArray_ImpactParameterSmearing(void *p);
   static void destruct_ImpactParameterSmearing(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ImpactParameterSmearing*)
   {
      ::ImpactParameterSmearing *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ImpactParameterSmearing >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ImpactParameterSmearing", ::ImpactParameterSmearing::Class_Version(), "modules/ImpactParameterSmearing.h", 36,
                  typeid(::ImpactParameterSmearing), DefineBehavior(ptr, ptr),
                  &::ImpactParameterSmearing::Dictionary, isa_proxy, 4,
                  sizeof(::ImpactParameterSmearing) );
      instance.SetNew(&new_ImpactParameterSmearing);
      instance.SetNewArray(&newArray_ImpactParameterSmearing);
      instance.SetDelete(&delete_ImpactParameterSmearing);
      instance.SetDeleteArray(&deleteArray_ImpactParameterSmearing);
      instance.SetDestructor(&destruct_ImpactParameterSmearing);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ImpactParameterSmearing*)
   {
      return GenerateInitInstanceLocal((::ImpactParameterSmearing*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ImpactParameterSmearing*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TimeSmearing(void *p = 0);
   static void *newArray_TimeSmearing(Long_t size, void *p);
   static void delete_TimeSmearing(void *p);
   static void deleteArray_TimeSmearing(void *p);
   static void destruct_TimeSmearing(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TimeSmearing*)
   {
      ::TimeSmearing *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TimeSmearing >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TimeSmearing", ::TimeSmearing::Class_Version(), "modules/TimeSmearing.h", 35,
                  typeid(::TimeSmearing), DefineBehavior(ptr, ptr),
                  &::TimeSmearing::Dictionary, isa_proxy, 4,
                  sizeof(::TimeSmearing) );
      instance.SetNew(&new_TimeSmearing);
      instance.SetNewArray(&newArray_TimeSmearing);
      instance.SetDelete(&delete_TimeSmearing);
      instance.SetDeleteArray(&deleteArray_TimeSmearing);
      instance.SetDestructor(&destruct_TimeSmearing);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TimeSmearing*)
   {
      return GenerateInitInstanceLocal((::TimeSmearing*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TimeSmearing*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SimpleCalorimeter(void *p = 0);
   static void *newArray_SimpleCalorimeter(Long_t size, void *p);
   static void delete_SimpleCalorimeter(void *p);
   static void deleteArray_SimpleCalorimeter(void *p);
   static void destruct_SimpleCalorimeter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimpleCalorimeter*)
   {
      ::SimpleCalorimeter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SimpleCalorimeter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SimpleCalorimeter", ::SimpleCalorimeter::Class_Version(), "modules/SimpleCalorimeter.h", 42,
                  typeid(::SimpleCalorimeter), DefineBehavior(ptr, ptr),
                  &::SimpleCalorimeter::Dictionary, isa_proxy, 4,
                  sizeof(::SimpleCalorimeter) );
      instance.SetNew(&new_SimpleCalorimeter);
      instance.SetNewArray(&newArray_SimpleCalorimeter);
      instance.SetDelete(&delete_SimpleCalorimeter);
      instance.SetDeleteArray(&deleteArray_SimpleCalorimeter);
      instance.SetDestructor(&destruct_SimpleCalorimeter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimpleCalorimeter*)
   {
      return GenerateInitInstanceLocal((::SimpleCalorimeter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::SimpleCalorimeter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Calorimeter(void *p = 0);
   static void *newArray_Calorimeter(Long_t size, void *p);
   static void delete_Calorimeter(void *p);
   static void deleteArray_Calorimeter(void *p);
   static void destruct_Calorimeter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Calorimeter*)
   {
      ::Calorimeter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Calorimeter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Calorimeter", ::Calorimeter::Class_Version(), "modules/Calorimeter.h", 41,
                  typeid(::Calorimeter), DefineBehavior(ptr, ptr),
                  &::Calorimeter::Dictionary, isa_proxy, 4,
                  sizeof(::Calorimeter) );
      instance.SetNew(&new_Calorimeter);
      instance.SetNewArray(&newArray_Calorimeter);
      instance.SetDelete(&delete_Calorimeter);
      instance.SetDeleteArray(&deleteArray_Calorimeter);
      instance.SetDestructor(&destruct_Calorimeter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Calorimeter*)
   {
      return GenerateInitInstanceLocal((::Calorimeter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Calorimeter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Isolation(void *p = 0);
   static void *newArray_Isolation(Long_t size, void *p);
   static void delete_Isolation(void *p);
   static void deleteArray_Isolation(void *p);
   static void destruct_Isolation(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Isolation*)
   {
      ::Isolation *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Isolation >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Isolation", ::Isolation::Class_Version(), "modules/Isolation.h", 40,
                  typeid(::Isolation), DefineBehavior(ptr, ptr),
                  &::Isolation::Dictionary, isa_proxy, 4,
                  sizeof(::Isolation) );
      instance.SetNew(&new_Isolation);
      instance.SetNewArray(&newArray_Isolation);
      instance.SetDelete(&delete_Isolation);
      instance.SetDeleteArray(&deleteArray_Isolation);
      instance.SetDestructor(&destruct_Isolation);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Isolation*)
   {
      return GenerateInitInstanceLocal((::Isolation*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Isolation*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_EnergyScale(void *p = 0);
   static void *newArray_EnergyScale(Long_t size, void *p);
   static void delete_EnergyScale(void *p);
   static void deleteArray_EnergyScale(void *p);
   static void destruct_EnergyScale(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EnergyScale*)
   {
      ::EnergyScale *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EnergyScale >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EnergyScale", ::EnergyScale::Class_Version(), "modules/EnergyScale.h", 36,
                  typeid(::EnergyScale), DefineBehavior(ptr, ptr),
                  &::EnergyScale::Dictionary, isa_proxy, 4,
                  sizeof(::EnergyScale) );
      instance.SetNew(&new_EnergyScale);
      instance.SetNewArray(&newArray_EnergyScale);
      instance.SetDelete(&delete_EnergyScale);
      instance.SetDeleteArray(&deleteArray_EnergyScale);
      instance.SetDestructor(&destruct_EnergyScale);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EnergyScale*)
   {
      return GenerateInitInstanceLocal((::EnergyScale*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::EnergyScale*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_UniqueObjectFinder(void *p = 0);
   static void *newArray_UniqueObjectFinder(Long_t size, void *p);
   static void delete_UniqueObjectFinder(void *p);
   static void deleteArray_UniqueObjectFinder(void *p);
   static void destruct_UniqueObjectFinder(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::UniqueObjectFinder*)
   {
      ::UniqueObjectFinder *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::UniqueObjectFinder >(0);
      static ::ROOT::TGenericClassInfo 
         instance("UniqueObjectFinder", ::UniqueObjectFinder::Class_Version(), "modules/UniqueObjectFinder.h", 39,
                  typeid(::UniqueObjectFinder), DefineBehavior(ptr, ptr),
                  &::UniqueObjectFinder::Dictionary, isa_proxy, 4,
                  sizeof(::UniqueObjectFinder) );
      instance.SetNew(&new_UniqueObjectFinder);
      instance.SetNewArray(&newArray_UniqueObjectFinder);
      instance.SetDelete(&delete_UniqueObjectFinder);
      instance.SetDeleteArray(&deleteArray_UniqueObjectFinder);
      instance.SetDestructor(&destruct_UniqueObjectFinder);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::UniqueObjectFinder*)
   {
      return GenerateInitInstanceLocal((::UniqueObjectFinder*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::UniqueObjectFinder*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TrackCountingBTagging(void *p = 0);
   static void *newArray_TrackCountingBTagging(Long_t size, void *p);
   static void delete_TrackCountingBTagging(void *p);
   static void deleteArray_TrackCountingBTagging(void *p);
   static void destruct_TrackCountingBTagging(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TrackCountingBTagging*)
   {
      ::TrackCountingBTagging *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TrackCountingBTagging >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TrackCountingBTagging", ::TrackCountingBTagging::Class_Version(), "modules/TrackCountingBTagging.h", 36,
                  typeid(::TrackCountingBTagging), DefineBehavior(ptr, ptr),
                  &::TrackCountingBTagging::Dictionary, isa_proxy, 4,
                  sizeof(::TrackCountingBTagging) );
      instance.SetNew(&new_TrackCountingBTagging);
      instance.SetNewArray(&newArray_TrackCountingBTagging);
      instance.SetDelete(&delete_TrackCountingBTagging);
      instance.SetDeleteArray(&deleteArray_TrackCountingBTagging);
      instance.SetDestructor(&destruct_TrackCountingBTagging);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TrackCountingBTagging*)
   {
      return GenerateInitInstanceLocal((::TrackCountingBTagging*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TrackCountingBTagging*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_BTagging(void *p = 0);
   static void *newArray_BTagging(Long_t size, void *p);
   static void delete_BTagging(void *p);
   static void deleteArray_BTagging(void *p);
   static void destruct_BTagging(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BTagging*)
   {
      ::BTagging *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::BTagging >(0);
      static ::ROOT::TGenericClassInfo 
         instance("BTagging", ::BTagging::Class_Version(), "modules/BTagging.h", 39,
                  typeid(::BTagging), DefineBehavior(ptr, ptr),
                  &::BTagging::Dictionary, isa_proxy, 4,
                  sizeof(::BTagging) );
      instance.SetNew(&new_BTagging);
      instance.SetNewArray(&newArray_BTagging);
      instance.SetDelete(&delete_BTagging);
      instance.SetDeleteArray(&deleteArray_BTagging);
      instance.SetDestructor(&destruct_BTagging);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BTagging*)
   {
      return GenerateInitInstanceLocal((::BTagging*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BTagging*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TauTagging(void *p = 0);
   static void *newArray_TauTagging(Long_t size, void *p);
   static void delete_TauTagging(void *p);
   static void deleteArray_TauTagging(void *p);
   static void destruct_TauTagging(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TauTagging*)
   {
      ::TauTagging *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TauTagging >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TauTagging", ::TauTagging::Class_Version(), "modules/TauTagging.h", 45,
                  typeid(::TauTagging), DefineBehavior(ptr, ptr),
                  &::TauTagging::Dictionary, isa_proxy, 4,
                  sizeof(::TauTagging) );
      instance.SetNew(&new_TauTagging);
      instance.SetNewArray(&newArray_TauTagging);
      instance.SetDelete(&delete_TauTagging);
      instance.SetDeleteArray(&deleteArray_TauTagging);
      instance.SetDestructor(&destruct_TauTagging);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TauTagging*)
   {
      return GenerateInitInstanceLocal((::TauTagging*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TauTagging*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TreeWriter(void *p = 0);
   static void *newArray_TreeWriter(Long_t size, void *p);
   static void delete_TreeWriter(void *p);
   static void deleteArray_TreeWriter(void *p);
   static void destruct_TreeWriter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TreeWriter*)
   {
      ::TreeWriter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TreeWriter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TreeWriter", ::TreeWriter::Class_Version(), "modules/TreeWriter.h", 41,
                  typeid(::TreeWriter), DefineBehavior(ptr, ptr),
                  &::TreeWriter::Dictionary, isa_proxy, 4,
                  sizeof(::TreeWriter) );
      instance.SetNew(&new_TreeWriter);
      instance.SetNewArray(&newArray_TreeWriter);
      instance.SetDelete(&delete_TreeWriter);
      instance.SetDeleteArray(&deleteArray_TreeWriter);
      instance.SetDestructor(&destruct_TreeWriter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TreeWriter*)
   {
      return GenerateInitInstanceLocal((::TreeWriter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TreeWriter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Merger(void *p = 0);
   static void *newArray_Merger(Long_t size, void *p);
   static void delete_Merger(void *p);
   static void deleteArray_Merger(void *p);
   static void destruct_Merger(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Merger*)
   {
      ::Merger *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Merger >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Merger", ::Merger::Class_Version(), "modules/Merger.h", 38,
                  typeid(::Merger), DefineBehavior(ptr, ptr),
                  &::Merger::Dictionary, isa_proxy, 4,
                  sizeof(::Merger) );
      instance.SetNew(&new_Merger);
      instance.SetNewArray(&newArray_Merger);
      instance.SetDelete(&delete_Merger);
      instance.SetDeleteArray(&deleteArray_Merger);
      instance.SetDestructor(&destruct_Merger);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Merger*)
   {
      return GenerateInitInstanceLocal((::Merger*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Merger*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_LeptonDressing(void *p = 0);
   static void *newArray_LeptonDressing(Long_t size, void *p);
   static void delete_LeptonDressing(void *p);
   static void deleteArray_LeptonDressing(void *p);
   static void destruct_LeptonDressing(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LeptonDressing*)
   {
      ::LeptonDressing *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LeptonDressing >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LeptonDressing", ::LeptonDressing::Class_Version(), "modules/LeptonDressing.h", 33,
                  typeid(::LeptonDressing), DefineBehavior(ptr, ptr),
                  &::LeptonDressing::Dictionary, isa_proxy, 4,
                  sizeof(::LeptonDressing) );
      instance.SetNew(&new_LeptonDressing);
      instance.SetNewArray(&newArray_LeptonDressing);
      instance.SetDelete(&delete_LeptonDressing);
      instance.SetDeleteArray(&deleteArray_LeptonDressing);
      instance.SetDestructor(&destruct_LeptonDressing);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LeptonDressing*)
   {
      return GenerateInitInstanceLocal((::LeptonDressing*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LeptonDressing*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PileUpMerger(void *p = 0);
   static void *newArray_PileUpMerger(Long_t size, void *p);
   static void delete_PileUpMerger(void *p);
   static void deleteArray_PileUpMerger(void *p);
   static void destruct_PileUpMerger(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PileUpMerger*)
   {
      ::PileUpMerger *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PileUpMerger >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PileUpMerger", ::PileUpMerger::Class_Version(), "modules/PileUpMerger.h", 36,
                  typeid(::PileUpMerger), DefineBehavior(ptr, ptr),
                  &::PileUpMerger::Dictionary, isa_proxy, 4,
                  sizeof(::PileUpMerger) );
      instance.SetNew(&new_PileUpMerger);
      instance.SetNewArray(&newArray_PileUpMerger);
      instance.SetDelete(&delete_PileUpMerger);
      instance.SetDeleteArray(&deleteArray_PileUpMerger);
      instance.SetDestructor(&destruct_PileUpMerger);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PileUpMerger*)
   {
      return GenerateInitInstanceLocal((::PileUpMerger*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::PileUpMerger*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_JetPileUpSubtractor(void *p = 0);
   static void *newArray_JetPileUpSubtractor(Long_t size, void *p);
   static void delete_JetPileUpSubtractor(void *p);
   static void deleteArray_JetPileUpSubtractor(void *p);
   static void destruct_JetPileUpSubtractor(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JetPileUpSubtractor*)
   {
      ::JetPileUpSubtractor *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::JetPileUpSubtractor >(0);
      static ::ROOT::TGenericClassInfo 
         instance("JetPileUpSubtractor", ::JetPileUpSubtractor::Class_Version(), "modules/JetPileUpSubtractor.h", 36,
                  typeid(::JetPileUpSubtractor), DefineBehavior(ptr, ptr),
                  &::JetPileUpSubtractor::Dictionary, isa_proxy, 4,
                  sizeof(::JetPileUpSubtractor) );
      instance.SetNew(&new_JetPileUpSubtractor);
      instance.SetNewArray(&newArray_JetPileUpSubtractor);
      instance.SetDelete(&delete_JetPileUpSubtractor);
      instance.SetDeleteArray(&deleteArray_JetPileUpSubtractor);
      instance.SetDestructor(&destruct_JetPileUpSubtractor);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JetPileUpSubtractor*)
   {
      return GenerateInitInstanceLocal((::JetPileUpSubtractor*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::JetPileUpSubtractor*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TrackPileUpSubtractor(void *p = 0);
   static void *newArray_TrackPileUpSubtractor(Long_t size, void *p);
   static void delete_TrackPileUpSubtractor(void *p);
   static void deleteArray_TrackPileUpSubtractor(void *p);
   static void destruct_TrackPileUpSubtractor(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TrackPileUpSubtractor*)
   {
      ::TrackPileUpSubtractor *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TrackPileUpSubtractor >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TrackPileUpSubtractor", ::TrackPileUpSubtractor::Class_Version(), "modules/TrackPileUpSubtractor.h", 37,
                  typeid(::TrackPileUpSubtractor), DefineBehavior(ptr, ptr),
                  &::TrackPileUpSubtractor::Dictionary, isa_proxy, 4,
                  sizeof(::TrackPileUpSubtractor) );
      instance.SetNew(&new_TrackPileUpSubtractor);
      instance.SetNewArray(&newArray_TrackPileUpSubtractor);
      instance.SetDelete(&delete_TrackPileUpSubtractor);
      instance.SetDeleteArray(&deleteArray_TrackPileUpSubtractor);
      instance.SetDestructor(&destruct_TrackPileUpSubtractor);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TrackPileUpSubtractor*)
   {
      return GenerateInitInstanceLocal((::TrackPileUpSubtractor*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TrackPileUpSubtractor*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TaggingParticlesSkimmer(void *p = 0);
   static void *newArray_TaggingParticlesSkimmer(Long_t size, void *p);
   static void delete_TaggingParticlesSkimmer(void *p);
   static void deleteArray_TaggingParticlesSkimmer(void *p);
   static void destruct_TaggingParticlesSkimmer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TaggingParticlesSkimmer*)
   {
      ::TaggingParticlesSkimmer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TaggingParticlesSkimmer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TaggingParticlesSkimmer", ::TaggingParticlesSkimmer::Class_Version(), "modules/TaggingParticlesSkimmer.h", 41,
                  typeid(::TaggingParticlesSkimmer), DefineBehavior(ptr, ptr),
                  &::TaggingParticlesSkimmer::Dictionary, isa_proxy, 4,
                  sizeof(::TaggingParticlesSkimmer) );
      instance.SetNew(&new_TaggingParticlesSkimmer);
      instance.SetNewArray(&newArray_TaggingParticlesSkimmer);
      instance.SetDelete(&delete_TaggingParticlesSkimmer);
      instance.SetDeleteArray(&deleteArray_TaggingParticlesSkimmer);
      instance.SetDestructor(&destruct_TaggingParticlesSkimmer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TaggingParticlesSkimmer*)
   {
      return GenerateInitInstanceLocal((::TaggingParticlesSkimmer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TaggingParticlesSkimmer*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PileUpJetID(void *p = 0);
   static void *newArray_PileUpJetID(Long_t size, void *p);
   static void delete_PileUpJetID(void *p);
   static void deleteArray_PileUpJetID(void *p);
   static void destruct_PileUpJetID(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PileUpJetID*)
   {
      ::PileUpJetID *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PileUpJetID >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PileUpJetID", ::PileUpJetID::Class_Version(), "modules/PileUpJetID.h", 20,
                  typeid(::PileUpJetID), DefineBehavior(ptr, ptr),
                  &::PileUpJetID::Dictionary, isa_proxy, 4,
                  sizeof(::PileUpJetID) );
      instance.SetNew(&new_PileUpJetID);
      instance.SetNewArray(&newArray_PileUpJetID);
      instance.SetDelete(&delete_PileUpJetID);
      instance.SetDeleteArray(&deleteArray_PileUpJetID);
      instance.SetDestructor(&destruct_PileUpJetID);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PileUpJetID*)
   {
      return GenerateInitInstanceLocal((::PileUpJetID*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::PileUpJetID*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ConstituentFilter(void *p = 0);
   static void *newArray_ConstituentFilter(Long_t size, void *p);
   static void delete_ConstituentFilter(void *p);
   static void deleteArray_ConstituentFilter(void *p);
   static void destruct_ConstituentFilter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ConstituentFilter*)
   {
      ::ConstituentFilter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ConstituentFilter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ConstituentFilter", ::ConstituentFilter::Class_Version(), "modules/ConstituentFilter.h", 38,
                  typeid(::ConstituentFilter), DefineBehavior(ptr, ptr),
                  &::ConstituentFilter::Dictionary, isa_proxy, 4,
                  sizeof(::ConstituentFilter) );
      instance.SetNew(&new_ConstituentFilter);
      instance.SetNewArray(&newArray_ConstituentFilter);
      instance.SetDelete(&delete_ConstituentFilter);
      instance.SetDeleteArray(&deleteArray_ConstituentFilter);
      instance.SetDestructor(&destruct_ConstituentFilter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ConstituentFilter*)
   {
      return GenerateInitInstanceLocal((::ConstituentFilter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ConstituentFilter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_StatusPidFilter(void *p = 0);
   static void *newArray_StatusPidFilter(Long_t size, void *p);
   static void delete_StatusPidFilter(void *p);
   static void deleteArray_StatusPidFilter(void *p);
   static void destruct_StatusPidFilter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StatusPidFilter*)
   {
      ::StatusPidFilter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StatusPidFilter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StatusPidFilter", ::StatusPidFilter::Class_Version(), "modules/StatusPidFilter.h", 38,
                  typeid(::StatusPidFilter), DefineBehavior(ptr, ptr),
                  &::StatusPidFilter::Dictionary, isa_proxy, 4,
                  sizeof(::StatusPidFilter) );
      instance.SetNew(&new_StatusPidFilter);
      instance.SetNewArray(&newArray_StatusPidFilter);
      instance.SetDelete(&delete_StatusPidFilter);
      instance.SetDeleteArray(&deleteArray_StatusPidFilter);
      instance.SetDestructor(&destruct_StatusPidFilter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StatusPidFilter*)
   {
      return GenerateInitInstanceLocal((::StatusPidFilter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::StatusPidFilter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PdgCodeFilter(void *p = 0);
   static void *newArray_PdgCodeFilter(Long_t size, void *p);
   static void delete_PdgCodeFilter(void *p);
   static void deleteArray_PdgCodeFilter(void *p);
   static void destruct_PdgCodeFilter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PdgCodeFilter*)
   {
      ::PdgCodeFilter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PdgCodeFilter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PdgCodeFilter", ::PdgCodeFilter::Class_Version(), "modules/PdgCodeFilter.h", 38,
                  typeid(::PdgCodeFilter), DefineBehavior(ptr, ptr),
                  &::PdgCodeFilter::Dictionary, isa_proxy, 4,
                  sizeof(::PdgCodeFilter) );
      instance.SetNew(&new_PdgCodeFilter);
      instance.SetNewArray(&newArray_PdgCodeFilter);
      instance.SetDelete(&delete_PdgCodeFilter);
      instance.SetDeleteArray(&deleteArray_PdgCodeFilter);
      instance.SetDestructor(&destruct_PdgCodeFilter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PdgCodeFilter*)
   {
      return GenerateInitInstanceLocal((::PdgCodeFilter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::PdgCodeFilter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Cloner(void *p = 0);
   static void *newArray_Cloner(Long_t size, void *p);
   static void delete_Cloner(void *p);
   static void deleteArray_Cloner(void *p);
   static void destruct_Cloner(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Cloner*)
   {
      ::Cloner *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Cloner >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Cloner", ::Cloner::Class_Version(), "modules/Cloner.h", 36,
                  typeid(::Cloner), DefineBehavior(ptr, ptr),
                  &::Cloner::Dictionary, isa_proxy, 4,
                  sizeof(::Cloner) );
      instance.SetNew(&new_Cloner);
      instance.SetNewArray(&newArray_Cloner);
      instance.SetDelete(&delete_Cloner);
      instance.SetDeleteArray(&deleteArray_Cloner);
      instance.SetDestructor(&destruct_Cloner);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Cloner*)
   {
      return GenerateInitInstanceLocal((::Cloner*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Cloner*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Weighter(void *p = 0);
   static void *newArray_Weighter(Long_t size, void *p);
   static void delete_Weighter(void *p);
   static void deleteArray_Weighter(void *p);
   static void destruct_Weighter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Weighter*)
   {
      ::Weighter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Weighter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Weighter", ::Weighter::Class_Version(), "modules/Weighter.h", 37,
                  typeid(::Weighter), DefineBehavior(ptr, ptr),
                  &::Weighter::Dictionary, isa_proxy, 4,
                  sizeof(::Weighter) );
      instance.SetNew(&new_Weighter);
      instance.SetNewArray(&newArray_Weighter);
      instance.SetDelete(&delete_Weighter);
      instance.SetDeleteArray(&deleteArray_Weighter);
      instance.SetDestructor(&destruct_Weighter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Weighter*)
   {
      return GenerateInitInstanceLocal((::Weighter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Weighter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Hector(void *p = 0);
   static void *newArray_Hector(Long_t size, void *p);
   static void delete_Hector(void *p);
   static void deleteArray_Hector(void *p);
   static void destruct_Hector(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Hector*)
   {
      ::Hector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Hector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Hector", ::Hector::Class_Version(), "modules/Hector.h", 36,
                  typeid(::Hector), DefineBehavior(ptr, ptr),
                  &::Hector::Dictionary, isa_proxy, 4,
                  sizeof(::Hector) );
      instance.SetNew(&new_Hector);
      instance.SetNewArray(&newArray_Hector);
      instance.SetDelete(&delete_Hector);
      instance.SetDeleteArray(&deleteArray_Hector);
      instance.SetDestructor(&destruct_Hector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Hector*)
   {
      return GenerateInitInstanceLocal((::Hector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Hector*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_JetFlavorAssociation(void *p = 0);
   static void *newArray_JetFlavorAssociation(Long_t size, void *p);
   static void delete_JetFlavorAssociation(void *p);
   static void deleteArray_JetFlavorAssociation(void *p);
   static void destruct_JetFlavorAssociation(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JetFlavorAssociation*)
   {
      ::JetFlavorAssociation *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::JetFlavorAssociation >(0);
      static ::ROOT::TGenericClassInfo 
         instance("JetFlavorAssociation", ::JetFlavorAssociation::Class_Version(), "modules/JetFlavorAssociation.h", 41,
                  typeid(::JetFlavorAssociation), DefineBehavior(ptr, ptr),
                  &::JetFlavorAssociation::Dictionary, isa_proxy, 4,
                  sizeof(::JetFlavorAssociation) );
      instance.SetNew(&new_JetFlavorAssociation);
      instance.SetNewArray(&newArray_JetFlavorAssociation);
      instance.SetDelete(&delete_JetFlavorAssociation);
      instance.SetDeleteArray(&deleteArray_JetFlavorAssociation);
      instance.SetDestructor(&destruct_JetFlavorAssociation);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JetFlavorAssociation*)
   {
      return GenerateInitInstanceLocal((::JetFlavorAssociation*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::JetFlavorAssociation*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_JetFakeParticle(void *p = 0);
   static void *newArray_JetFakeParticle(Long_t size, void *p);
   static void delete_JetFakeParticle(void *p);
   static void deleteArray_JetFakeParticle(void *p);
   static void destruct_JetFakeParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JetFakeParticle*)
   {
      ::JetFakeParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::JetFakeParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("JetFakeParticle", ::JetFakeParticle::Class_Version(), "modules/JetFakeParticle.h", 38,
                  typeid(::JetFakeParticle), DefineBehavior(ptr, ptr),
                  &::JetFakeParticle::Dictionary, isa_proxy, 4,
                  sizeof(::JetFakeParticle) );
      instance.SetNew(&new_JetFakeParticle);
      instance.SetNewArray(&newArray_JetFakeParticle);
      instance.SetDelete(&delete_JetFakeParticle);
      instance.SetDeleteArray(&deleteArray_JetFakeParticle);
      instance.SetDestructor(&destruct_JetFakeParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JetFakeParticle*)
   {
      return GenerateInitInstanceLocal((::JetFakeParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::JetFakeParticle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ExampleModule(void *p = 0);
   static void *newArray_ExampleModule(Long_t size, void *p);
   static void delete_ExampleModule(void *p);
   static void deleteArray_ExampleModule(void *p);
   static void destruct_ExampleModule(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExampleModule*)
   {
      ::ExampleModule *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExampleModule >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ExampleModule", ::ExampleModule::Class_Version(), "modules/ExampleModule.h", 37,
                  typeid(::ExampleModule), DefineBehavior(ptr, ptr),
                  &::ExampleModule::Dictionary, isa_proxy, 4,
                  sizeof(::ExampleModule) );
      instance.SetNew(&new_ExampleModule);
      instance.SetNewArray(&newArray_ExampleModule);
      instance.SetDelete(&delete_ExampleModule);
      instance.SetDeleteArray(&deleteArray_ExampleModule);
      instance.SetDestructor(&destruct_ExampleModule);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExampleModule*)
   {
      return GenerateInitInstanceLocal((::ExampleModule*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExampleModule*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Delphes::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Delphes::Class_Name()
{
   return "Delphes";
}

//______________________________________________________________________________
const char *Delphes::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Delphes*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Delphes::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Delphes*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Delphes::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Delphes*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Delphes::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Delphes*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr AngularSmearing::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *AngularSmearing::Class_Name()
{
   return "AngularSmearing";
}

//______________________________________________________________________________
const char *AngularSmearing::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AngularSmearing*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int AngularSmearing::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AngularSmearing*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *AngularSmearing::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AngularSmearing*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *AngularSmearing::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AngularSmearing*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PhotonConversions::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PhotonConversions::Class_Name()
{
   return "PhotonConversions";
}

//______________________________________________________________________________
const char *PhotonConversions::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PhotonConversions*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PhotonConversions::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PhotonConversions*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PhotonConversions::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PhotonConversions*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PhotonConversions::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PhotonConversions*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ParticlePropagator::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ParticlePropagator::Class_Name()
{
   return "ParticlePropagator";
}

//______________________________________________________________________________
const char *ParticlePropagator::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ParticlePropagator*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ParticlePropagator::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ParticlePropagator*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ParticlePropagator::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ParticlePropagator*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ParticlePropagator::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ParticlePropagator*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Efficiency::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Efficiency::Class_Name()
{
   return "Efficiency";
}

//______________________________________________________________________________
const char *Efficiency::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Efficiency*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Efficiency::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Efficiency*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Efficiency::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Efficiency*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Efficiency::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Efficiency*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr IdentificationMap::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *IdentificationMap::Class_Name()
{
   return "IdentificationMap";
}

//______________________________________________________________________________
const char *IdentificationMap::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::IdentificationMap*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int IdentificationMap::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::IdentificationMap*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *IdentificationMap::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::IdentificationMap*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *IdentificationMap::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::IdentificationMap*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr EnergySmearing::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EnergySmearing::Class_Name()
{
   return "EnergySmearing";
}

//______________________________________________________________________________
const char *EnergySmearing::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EnergySmearing*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EnergySmearing::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EnergySmearing*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EnergySmearing::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EnergySmearing*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EnergySmearing::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EnergySmearing*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MomentumSmearing::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MomentumSmearing::Class_Name()
{
   return "MomentumSmearing";
}

//______________________________________________________________________________
const char *MomentumSmearing::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MomentumSmearing*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MomentumSmearing::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MomentumSmearing*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MomentumSmearing::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MomentumSmearing*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MomentumSmearing::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MomentumSmearing*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ImpactParameterSmearing::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ImpactParameterSmearing::Class_Name()
{
   return "ImpactParameterSmearing";
}

//______________________________________________________________________________
const char *ImpactParameterSmearing::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ImpactParameterSmearing*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ImpactParameterSmearing::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ImpactParameterSmearing*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ImpactParameterSmearing::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ImpactParameterSmearing*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ImpactParameterSmearing::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ImpactParameterSmearing*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TimeSmearing::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TimeSmearing::Class_Name()
{
   return "TimeSmearing";
}

//______________________________________________________________________________
const char *TimeSmearing::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TimeSmearing*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TimeSmearing::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TimeSmearing*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TimeSmearing::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TimeSmearing*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TimeSmearing::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TimeSmearing*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SimpleCalorimeter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SimpleCalorimeter::Class_Name()
{
   return "SimpleCalorimeter";
}

//______________________________________________________________________________
const char *SimpleCalorimeter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimpleCalorimeter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SimpleCalorimeter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimpleCalorimeter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SimpleCalorimeter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimpleCalorimeter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SimpleCalorimeter::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimpleCalorimeter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Calorimeter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Calorimeter::Class_Name()
{
   return "Calorimeter";
}

//______________________________________________________________________________
const char *Calorimeter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Calorimeter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Calorimeter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Calorimeter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Calorimeter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Calorimeter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Calorimeter::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Calorimeter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Isolation::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Isolation::Class_Name()
{
   return "Isolation";
}

//______________________________________________________________________________
const char *Isolation::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Isolation*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Isolation::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Isolation*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Isolation::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Isolation*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Isolation::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Isolation*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr EnergyScale::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EnergyScale::Class_Name()
{
   return "EnergyScale";
}

//______________________________________________________________________________
const char *EnergyScale::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EnergyScale*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EnergyScale::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EnergyScale*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EnergyScale::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EnergyScale*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EnergyScale::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EnergyScale*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr UniqueObjectFinder::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *UniqueObjectFinder::Class_Name()
{
   return "UniqueObjectFinder";
}

//______________________________________________________________________________
const char *UniqueObjectFinder::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UniqueObjectFinder*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int UniqueObjectFinder::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UniqueObjectFinder*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *UniqueObjectFinder::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UniqueObjectFinder*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *UniqueObjectFinder::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UniqueObjectFinder*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TrackCountingBTagging::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TrackCountingBTagging::Class_Name()
{
   return "TrackCountingBTagging";
}

//______________________________________________________________________________
const char *TrackCountingBTagging::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TrackCountingBTagging*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TrackCountingBTagging::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TrackCountingBTagging*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TrackCountingBTagging::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TrackCountingBTagging*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TrackCountingBTagging::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TrackCountingBTagging*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr BTagging::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *BTagging::Class_Name()
{
   return "BTagging";
}

//______________________________________________________________________________
const char *BTagging::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::BTagging*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int BTagging::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::BTagging*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *BTagging::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::BTagging*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *BTagging::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::BTagging*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TauTagging::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TauTagging::Class_Name()
{
   return "TauTagging";
}

//______________________________________________________________________________
const char *TauTagging::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TauTagging*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TauTagging::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TauTagging*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TauTagging::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TauTagging*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TauTagging::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TauTagging*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TreeWriter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TreeWriter::Class_Name()
{
   return "TreeWriter";
}

//______________________________________________________________________________
const char *TreeWriter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TreeWriter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TreeWriter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TreeWriter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TreeWriter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TreeWriter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TreeWriter::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TreeWriter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Merger::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Merger::Class_Name()
{
   return "Merger";
}

//______________________________________________________________________________
const char *Merger::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Merger*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Merger::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Merger*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Merger::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Merger*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Merger::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Merger*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr LeptonDressing::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LeptonDressing::Class_Name()
{
   return "LeptonDressing";
}

//______________________________________________________________________________
const char *LeptonDressing::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LeptonDressing*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LeptonDressing::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LeptonDressing*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LeptonDressing::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LeptonDressing*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LeptonDressing::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LeptonDressing*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PileUpMerger::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PileUpMerger::Class_Name()
{
   return "PileUpMerger";
}

//______________________________________________________________________________
const char *PileUpMerger::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PileUpMerger*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PileUpMerger::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PileUpMerger*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PileUpMerger::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PileUpMerger*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PileUpMerger::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PileUpMerger*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr JetPileUpSubtractor::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *JetPileUpSubtractor::Class_Name()
{
   return "JetPileUpSubtractor";
}

//______________________________________________________________________________
const char *JetPileUpSubtractor::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JetPileUpSubtractor*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int JetPileUpSubtractor::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JetPileUpSubtractor*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JetPileUpSubtractor::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JetPileUpSubtractor*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JetPileUpSubtractor::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JetPileUpSubtractor*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TrackPileUpSubtractor::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TrackPileUpSubtractor::Class_Name()
{
   return "TrackPileUpSubtractor";
}

//______________________________________________________________________________
const char *TrackPileUpSubtractor::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TrackPileUpSubtractor*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TrackPileUpSubtractor::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TrackPileUpSubtractor*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TrackPileUpSubtractor::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TrackPileUpSubtractor*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TrackPileUpSubtractor::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TrackPileUpSubtractor*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TaggingParticlesSkimmer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TaggingParticlesSkimmer::Class_Name()
{
   return "TaggingParticlesSkimmer";
}

//______________________________________________________________________________
const char *TaggingParticlesSkimmer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TaggingParticlesSkimmer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TaggingParticlesSkimmer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TaggingParticlesSkimmer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TaggingParticlesSkimmer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TaggingParticlesSkimmer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TaggingParticlesSkimmer::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TaggingParticlesSkimmer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PileUpJetID::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PileUpJetID::Class_Name()
{
   return "PileUpJetID";
}

//______________________________________________________________________________
const char *PileUpJetID::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PileUpJetID*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PileUpJetID::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PileUpJetID*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PileUpJetID::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PileUpJetID*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PileUpJetID::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PileUpJetID*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ConstituentFilter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ConstituentFilter::Class_Name()
{
   return "ConstituentFilter";
}

//______________________________________________________________________________
const char *ConstituentFilter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ConstituentFilter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ConstituentFilter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ConstituentFilter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ConstituentFilter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ConstituentFilter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ConstituentFilter::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ConstituentFilter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr StatusPidFilter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StatusPidFilter::Class_Name()
{
   return "StatusPidFilter";
}

//______________________________________________________________________________
const char *StatusPidFilter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StatusPidFilter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StatusPidFilter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StatusPidFilter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StatusPidFilter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StatusPidFilter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StatusPidFilter::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StatusPidFilter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PdgCodeFilter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PdgCodeFilter::Class_Name()
{
   return "PdgCodeFilter";
}

//______________________________________________________________________________
const char *PdgCodeFilter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdgCodeFilter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PdgCodeFilter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdgCodeFilter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PdgCodeFilter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdgCodeFilter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PdgCodeFilter::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdgCodeFilter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Cloner::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Cloner::Class_Name()
{
   return "Cloner";
}

//______________________________________________________________________________
const char *Cloner::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Cloner*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Cloner::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Cloner*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Cloner::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Cloner*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Cloner::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Cloner*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Weighter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Weighter::Class_Name()
{
   return "Weighter";
}

//______________________________________________________________________________
const char *Weighter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Weighter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Weighter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Weighter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Weighter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Weighter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Weighter::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Weighter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Hector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Hector::Class_Name()
{
   return "Hector";
}

//______________________________________________________________________________
const char *Hector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Hector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Hector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Hector::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr JetFlavorAssociation::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *JetFlavorAssociation::Class_Name()
{
   return "JetFlavorAssociation";
}

//______________________________________________________________________________
const char *JetFlavorAssociation::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JetFlavorAssociation*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int JetFlavorAssociation::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JetFlavorAssociation*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JetFlavorAssociation::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JetFlavorAssociation*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JetFlavorAssociation::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JetFlavorAssociation*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr JetFakeParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *JetFakeParticle::Class_Name()
{
   return "JetFakeParticle";
}

//______________________________________________________________________________
const char *JetFakeParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JetFakeParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int JetFakeParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JetFakeParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JetFakeParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JetFakeParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JetFakeParticle::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JetFakeParticle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ExampleModule::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ExampleModule::Class_Name()
{
   return "ExampleModule";
}

//______________________________________________________________________________
const char *ExampleModule::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExampleModule*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ExampleModule::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExampleModule*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ExampleModule::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExampleModule*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ExampleModule::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExampleModule*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Delphes::Streamer(TBuffer &R__b)
{
   // Stream an object of class Delphes.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Delphes::Class(),this);
   } else {
      R__b.WriteClassBuffer(Delphes::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Delphes(void *p) {
      return  p ? new(p) ::Delphes : new ::Delphes;
   }
   static void *newArray_Delphes(Long_t nElements, void *p) {
      return p ? new(p) ::Delphes[nElements] : new ::Delphes[nElements];
   }
   // Wrapper around operator delete
   static void delete_Delphes(void *p) {
      delete ((::Delphes*)p);
   }
   static void deleteArray_Delphes(void *p) {
      delete [] ((::Delphes*)p);
   }
   static void destruct_Delphes(void *p) {
      typedef ::Delphes current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Delphes

//______________________________________________________________________________
void AngularSmearing::Streamer(TBuffer &R__b)
{
   // Stream an object of class AngularSmearing.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(AngularSmearing::Class(),this);
   } else {
      R__b.WriteClassBuffer(AngularSmearing::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_AngularSmearing(void *p) {
      return  p ? new(p) ::AngularSmearing : new ::AngularSmearing;
   }
   static void *newArray_AngularSmearing(Long_t nElements, void *p) {
      return p ? new(p) ::AngularSmearing[nElements] : new ::AngularSmearing[nElements];
   }
   // Wrapper around operator delete
   static void delete_AngularSmearing(void *p) {
      delete ((::AngularSmearing*)p);
   }
   static void deleteArray_AngularSmearing(void *p) {
      delete [] ((::AngularSmearing*)p);
   }
   static void destruct_AngularSmearing(void *p) {
      typedef ::AngularSmearing current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AngularSmearing

//______________________________________________________________________________
void PhotonConversions::Streamer(TBuffer &R__b)
{
   // Stream an object of class PhotonConversions.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PhotonConversions::Class(),this);
   } else {
      R__b.WriteClassBuffer(PhotonConversions::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PhotonConversions(void *p) {
      return  p ? new(p) ::PhotonConversions : new ::PhotonConversions;
   }
   static void *newArray_PhotonConversions(Long_t nElements, void *p) {
      return p ? new(p) ::PhotonConversions[nElements] : new ::PhotonConversions[nElements];
   }
   // Wrapper around operator delete
   static void delete_PhotonConversions(void *p) {
      delete ((::PhotonConversions*)p);
   }
   static void deleteArray_PhotonConversions(void *p) {
      delete [] ((::PhotonConversions*)p);
   }
   static void destruct_PhotonConversions(void *p) {
      typedef ::PhotonConversions current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PhotonConversions

//______________________________________________________________________________
void ParticlePropagator::Streamer(TBuffer &R__b)
{
   // Stream an object of class ParticlePropagator.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ParticlePropagator::Class(),this);
   } else {
      R__b.WriteClassBuffer(ParticlePropagator::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ParticlePropagator(void *p) {
      return  p ? new(p) ::ParticlePropagator : new ::ParticlePropagator;
   }
   static void *newArray_ParticlePropagator(Long_t nElements, void *p) {
      return p ? new(p) ::ParticlePropagator[nElements] : new ::ParticlePropagator[nElements];
   }
   // Wrapper around operator delete
   static void delete_ParticlePropagator(void *p) {
      delete ((::ParticlePropagator*)p);
   }
   static void deleteArray_ParticlePropagator(void *p) {
      delete [] ((::ParticlePropagator*)p);
   }
   static void destruct_ParticlePropagator(void *p) {
      typedef ::ParticlePropagator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ParticlePropagator

//______________________________________________________________________________
void Efficiency::Streamer(TBuffer &R__b)
{
   // Stream an object of class Efficiency.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Efficiency::Class(),this);
   } else {
      R__b.WriteClassBuffer(Efficiency::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Efficiency(void *p) {
      return  p ? new(p) ::Efficiency : new ::Efficiency;
   }
   static void *newArray_Efficiency(Long_t nElements, void *p) {
      return p ? new(p) ::Efficiency[nElements] : new ::Efficiency[nElements];
   }
   // Wrapper around operator delete
   static void delete_Efficiency(void *p) {
      delete ((::Efficiency*)p);
   }
   static void deleteArray_Efficiency(void *p) {
      delete [] ((::Efficiency*)p);
   }
   static void destruct_Efficiency(void *p) {
      typedef ::Efficiency current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Efficiency

//______________________________________________________________________________
void IdentificationMap::Streamer(TBuffer &R__b)
{
   // Stream an object of class IdentificationMap.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(IdentificationMap::Class(),this);
   } else {
      R__b.WriteClassBuffer(IdentificationMap::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_IdentificationMap(void *p) {
      return  p ? new(p) ::IdentificationMap : new ::IdentificationMap;
   }
   static void *newArray_IdentificationMap(Long_t nElements, void *p) {
      return p ? new(p) ::IdentificationMap[nElements] : new ::IdentificationMap[nElements];
   }
   // Wrapper around operator delete
   static void delete_IdentificationMap(void *p) {
      delete ((::IdentificationMap*)p);
   }
   static void deleteArray_IdentificationMap(void *p) {
      delete [] ((::IdentificationMap*)p);
   }
   static void destruct_IdentificationMap(void *p) {
      typedef ::IdentificationMap current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::IdentificationMap

//______________________________________________________________________________
void EnergySmearing::Streamer(TBuffer &R__b)
{
   // Stream an object of class EnergySmearing.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EnergySmearing::Class(),this);
   } else {
      R__b.WriteClassBuffer(EnergySmearing::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EnergySmearing(void *p) {
      return  p ? new(p) ::EnergySmearing : new ::EnergySmearing;
   }
   static void *newArray_EnergySmearing(Long_t nElements, void *p) {
      return p ? new(p) ::EnergySmearing[nElements] : new ::EnergySmearing[nElements];
   }
   // Wrapper around operator delete
   static void delete_EnergySmearing(void *p) {
      delete ((::EnergySmearing*)p);
   }
   static void deleteArray_EnergySmearing(void *p) {
      delete [] ((::EnergySmearing*)p);
   }
   static void destruct_EnergySmearing(void *p) {
      typedef ::EnergySmearing current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EnergySmearing

//______________________________________________________________________________
void MomentumSmearing::Streamer(TBuffer &R__b)
{
   // Stream an object of class MomentumSmearing.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MomentumSmearing::Class(),this);
   } else {
      R__b.WriteClassBuffer(MomentumSmearing::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MomentumSmearing(void *p) {
      return  p ? new(p) ::MomentumSmearing : new ::MomentumSmearing;
   }
   static void *newArray_MomentumSmearing(Long_t nElements, void *p) {
      return p ? new(p) ::MomentumSmearing[nElements] : new ::MomentumSmearing[nElements];
   }
   // Wrapper around operator delete
   static void delete_MomentumSmearing(void *p) {
      delete ((::MomentumSmearing*)p);
   }
   static void deleteArray_MomentumSmearing(void *p) {
      delete [] ((::MomentumSmearing*)p);
   }
   static void destruct_MomentumSmearing(void *p) {
      typedef ::MomentumSmearing current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MomentumSmearing

//______________________________________________________________________________
void ImpactParameterSmearing::Streamer(TBuffer &R__b)
{
   // Stream an object of class ImpactParameterSmearing.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ImpactParameterSmearing::Class(),this);
   } else {
      R__b.WriteClassBuffer(ImpactParameterSmearing::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ImpactParameterSmearing(void *p) {
      return  p ? new(p) ::ImpactParameterSmearing : new ::ImpactParameterSmearing;
   }
   static void *newArray_ImpactParameterSmearing(Long_t nElements, void *p) {
      return p ? new(p) ::ImpactParameterSmearing[nElements] : new ::ImpactParameterSmearing[nElements];
   }
   // Wrapper around operator delete
   static void delete_ImpactParameterSmearing(void *p) {
      delete ((::ImpactParameterSmearing*)p);
   }
   static void deleteArray_ImpactParameterSmearing(void *p) {
      delete [] ((::ImpactParameterSmearing*)p);
   }
   static void destruct_ImpactParameterSmearing(void *p) {
      typedef ::ImpactParameterSmearing current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ImpactParameterSmearing

//______________________________________________________________________________
void TimeSmearing::Streamer(TBuffer &R__b)
{
   // Stream an object of class TimeSmearing.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TimeSmearing::Class(),this);
   } else {
      R__b.WriteClassBuffer(TimeSmearing::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TimeSmearing(void *p) {
      return  p ? new(p) ::TimeSmearing : new ::TimeSmearing;
   }
   static void *newArray_TimeSmearing(Long_t nElements, void *p) {
      return p ? new(p) ::TimeSmearing[nElements] : new ::TimeSmearing[nElements];
   }
   // Wrapper around operator delete
   static void delete_TimeSmearing(void *p) {
      delete ((::TimeSmearing*)p);
   }
   static void deleteArray_TimeSmearing(void *p) {
      delete [] ((::TimeSmearing*)p);
   }
   static void destruct_TimeSmearing(void *p) {
      typedef ::TimeSmearing current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TimeSmearing

//______________________________________________________________________________
void SimpleCalorimeter::Streamer(TBuffer &R__b)
{
   // Stream an object of class SimpleCalorimeter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SimpleCalorimeter::Class(),this);
   } else {
      R__b.WriteClassBuffer(SimpleCalorimeter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SimpleCalorimeter(void *p) {
      return  p ? new(p) ::SimpleCalorimeter : new ::SimpleCalorimeter;
   }
   static void *newArray_SimpleCalorimeter(Long_t nElements, void *p) {
      return p ? new(p) ::SimpleCalorimeter[nElements] : new ::SimpleCalorimeter[nElements];
   }
   // Wrapper around operator delete
   static void delete_SimpleCalorimeter(void *p) {
      delete ((::SimpleCalorimeter*)p);
   }
   static void deleteArray_SimpleCalorimeter(void *p) {
      delete [] ((::SimpleCalorimeter*)p);
   }
   static void destruct_SimpleCalorimeter(void *p) {
      typedef ::SimpleCalorimeter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimpleCalorimeter

//______________________________________________________________________________
void Calorimeter::Streamer(TBuffer &R__b)
{
   // Stream an object of class Calorimeter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Calorimeter::Class(),this);
   } else {
      R__b.WriteClassBuffer(Calorimeter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Calorimeter(void *p) {
      return  p ? new(p) ::Calorimeter : new ::Calorimeter;
   }
   static void *newArray_Calorimeter(Long_t nElements, void *p) {
      return p ? new(p) ::Calorimeter[nElements] : new ::Calorimeter[nElements];
   }
   // Wrapper around operator delete
   static void delete_Calorimeter(void *p) {
      delete ((::Calorimeter*)p);
   }
   static void deleteArray_Calorimeter(void *p) {
      delete [] ((::Calorimeter*)p);
   }
   static void destruct_Calorimeter(void *p) {
      typedef ::Calorimeter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Calorimeter

//______________________________________________________________________________
void Isolation::Streamer(TBuffer &R__b)
{
   // Stream an object of class Isolation.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Isolation::Class(),this);
   } else {
      R__b.WriteClassBuffer(Isolation::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Isolation(void *p) {
      return  p ? new(p) ::Isolation : new ::Isolation;
   }
   static void *newArray_Isolation(Long_t nElements, void *p) {
      return p ? new(p) ::Isolation[nElements] : new ::Isolation[nElements];
   }
   // Wrapper around operator delete
   static void delete_Isolation(void *p) {
      delete ((::Isolation*)p);
   }
   static void deleteArray_Isolation(void *p) {
      delete [] ((::Isolation*)p);
   }
   static void destruct_Isolation(void *p) {
      typedef ::Isolation current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Isolation

//______________________________________________________________________________
void EnergyScale::Streamer(TBuffer &R__b)
{
   // Stream an object of class EnergyScale.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EnergyScale::Class(),this);
   } else {
      R__b.WriteClassBuffer(EnergyScale::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EnergyScale(void *p) {
      return  p ? new(p) ::EnergyScale : new ::EnergyScale;
   }
   static void *newArray_EnergyScale(Long_t nElements, void *p) {
      return p ? new(p) ::EnergyScale[nElements] : new ::EnergyScale[nElements];
   }
   // Wrapper around operator delete
   static void delete_EnergyScale(void *p) {
      delete ((::EnergyScale*)p);
   }
   static void deleteArray_EnergyScale(void *p) {
      delete [] ((::EnergyScale*)p);
   }
   static void destruct_EnergyScale(void *p) {
      typedef ::EnergyScale current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EnergyScale

//______________________________________________________________________________
void UniqueObjectFinder::Streamer(TBuffer &R__b)
{
   // Stream an object of class UniqueObjectFinder.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(UniqueObjectFinder::Class(),this);
   } else {
      R__b.WriteClassBuffer(UniqueObjectFinder::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_UniqueObjectFinder(void *p) {
      return  p ? new(p) ::UniqueObjectFinder : new ::UniqueObjectFinder;
   }
   static void *newArray_UniqueObjectFinder(Long_t nElements, void *p) {
      return p ? new(p) ::UniqueObjectFinder[nElements] : new ::UniqueObjectFinder[nElements];
   }
   // Wrapper around operator delete
   static void delete_UniqueObjectFinder(void *p) {
      delete ((::UniqueObjectFinder*)p);
   }
   static void deleteArray_UniqueObjectFinder(void *p) {
      delete [] ((::UniqueObjectFinder*)p);
   }
   static void destruct_UniqueObjectFinder(void *p) {
      typedef ::UniqueObjectFinder current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::UniqueObjectFinder

//______________________________________________________________________________
void TrackCountingBTagging::Streamer(TBuffer &R__b)
{
   // Stream an object of class TrackCountingBTagging.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TrackCountingBTagging::Class(),this);
   } else {
      R__b.WriteClassBuffer(TrackCountingBTagging::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TrackCountingBTagging(void *p) {
      return  p ? new(p) ::TrackCountingBTagging : new ::TrackCountingBTagging;
   }
   static void *newArray_TrackCountingBTagging(Long_t nElements, void *p) {
      return p ? new(p) ::TrackCountingBTagging[nElements] : new ::TrackCountingBTagging[nElements];
   }
   // Wrapper around operator delete
   static void delete_TrackCountingBTagging(void *p) {
      delete ((::TrackCountingBTagging*)p);
   }
   static void deleteArray_TrackCountingBTagging(void *p) {
      delete [] ((::TrackCountingBTagging*)p);
   }
   static void destruct_TrackCountingBTagging(void *p) {
      typedef ::TrackCountingBTagging current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TrackCountingBTagging

//______________________________________________________________________________
void BTagging::Streamer(TBuffer &R__b)
{
   // Stream an object of class BTagging.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(BTagging::Class(),this);
   } else {
      R__b.WriteClassBuffer(BTagging::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_BTagging(void *p) {
      return  p ? new(p) ::BTagging : new ::BTagging;
   }
   static void *newArray_BTagging(Long_t nElements, void *p) {
      return p ? new(p) ::BTagging[nElements] : new ::BTagging[nElements];
   }
   // Wrapper around operator delete
   static void delete_BTagging(void *p) {
      delete ((::BTagging*)p);
   }
   static void deleteArray_BTagging(void *p) {
      delete [] ((::BTagging*)p);
   }
   static void destruct_BTagging(void *p) {
      typedef ::BTagging current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BTagging

//______________________________________________________________________________
void TauTagging::Streamer(TBuffer &R__b)
{
   // Stream an object of class TauTagging.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TauTagging::Class(),this);
   } else {
      R__b.WriteClassBuffer(TauTagging::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TauTagging(void *p) {
      return  p ? new(p) ::TauTagging : new ::TauTagging;
   }
   static void *newArray_TauTagging(Long_t nElements, void *p) {
      return p ? new(p) ::TauTagging[nElements] : new ::TauTagging[nElements];
   }
   // Wrapper around operator delete
   static void delete_TauTagging(void *p) {
      delete ((::TauTagging*)p);
   }
   static void deleteArray_TauTagging(void *p) {
      delete [] ((::TauTagging*)p);
   }
   static void destruct_TauTagging(void *p) {
      typedef ::TauTagging current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TauTagging

//______________________________________________________________________________
void TreeWriter::Streamer(TBuffer &R__b)
{
   // Stream an object of class TreeWriter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TreeWriter::Class(),this);
   } else {
      R__b.WriteClassBuffer(TreeWriter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TreeWriter(void *p) {
      return  p ? new(p) ::TreeWriter : new ::TreeWriter;
   }
   static void *newArray_TreeWriter(Long_t nElements, void *p) {
      return p ? new(p) ::TreeWriter[nElements] : new ::TreeWriter[nElements];
   }
   // Wrapper around operator delete
   static void delete_TreeWriter(void *p) {
      delete ((::TreeWriter*)p);
   }
   static void deleteArray_TreeWriter(void *p) {
      delete [] ((::TreeWriter*)p);
   }
   static void destruct_TreeWriter(void *p) {
      typedef ::TreeWriter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TreeWriter

//______________________________________________________________________________
void Merger::Streamer(TBuffer &R__b)
{
   // Stream an object of class Merger.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Merger::Class(),this);
   } else {
      R__b.WriteClassBuffer(Merger::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Merger(void *p) {
      return  p ? new(p) ::Merger : new ::Merger;
   }
   static void *newArray_Merger(Long_t nElements, void *p) {
      return p ? new(p) ::Merger[nElements] : new ::Merger[nElements];
   }
   // Wrapper around operator delete
   static void delete_Merger(void *p) {
      delete ((::Merger*)p);
   }
   static void deleteArray_Merger(void *p) {
      delete [] ((::Merger*)p);
   }
   static void destruct_Merger(void *p) {
      typedef ::Merger current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Merger

//______________________________________________________________________________
void LeptonDressing::Streamer(TBuffer &R__b)
{
   // Stream an object of class LeptonDressing.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LeptonDressing::Class(),this);
   } else {
      R__b.WriteClassBuffer(LeptonDressing::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LeptonDressing(void *p) {
      return  p ? new(p) ::LeptonDressing : new ::LeptonDressing;
   }
   static void *newArray_LeptonDressing(Long_t nElements, void *p) {
      return p ? new(p) ::LeptonDressing[nElements] : new ::LeptonDressing[nElements];
   }
   // Wrapper around operator delete
   static void delete_LeptonDressing(void *p) {
      delete ((::LeptonDressing*)p);
   }
   static void deleteArray_LeptonDressing(void *p) {
      delete [] ((::LeptonDressing*)p);
   }
   static void destruct_LeptonDressing(void *p) {
      typedef ::LeptonDressing current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LeptonDressing

//______________________________________________________________________________
void PileUpMerger::Streamer(TBuffer &R__b)
{
   // Stream an object of class PileUpMerger.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PileUpMerger::Class(),this);
   } else {
      R__b.WriteClassBuffer(PileUpMerger::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PileUpMerger(void *p) {
      return  p ? new(p) ::PileUpMerger : new ::PileUpMerger;
   }
   static void *newArray_PileUpMerger(Long_t nElements, void *p) {
      return p ? new(p) ::PileUpMerger[nElements] : new ::PileUpMerger[nElements];
   }
   // Wrapper around operator delete
   static void delete_PileUpMerger(void *p) {
      delete ((::PileUpMerger*)p);
   }
   static void deleteArray_PileUpMerger(void *p) {
      delete [] ((::PileUpMerger*)p);
   }
   static void destruct_PileUpMerger(void *p) {
      typedef ::PileUpMerger current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PileUpMerger

//______________________________________________________________________________
void JetPileUpSubtractor::Streamer(TBuffer &R__b)
{
   // Stream an object of class JetPileUpSubtractor.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(JetPileUpSubtractor::Class(),this);
   } else {
      R__b.WriteClassBuffer(JetPileUpSubtractor::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_JetPileUpSubtractor(void *p) {
      return  p ? new(p) ::JetPileUpSubtractor : new ::JetPileUpSubtractor;
   }
   static void *newArray_JetPileUpSubtractor(Long_t nElements, void *p) {
      return p ? new(p) ::JetPileUpSubtractor[nElements] : new ::JetPileUpSubtractor[nElements];
   }
   // Wrapper around operator delete
   static void delete_JetPileUpSubtractor(void *p) {
      delete ((::JetPileUpSubtractor*)p);
   }
   static void deleteArray_JetPileUpSubtractor(void *p) {
      delete [] ((::JetPileUpSubtractor*)p);
   }
   static void destruct_JetPileUpSubtractor(void *p) {
      typedef ::JetPileUpSubtractor current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::JetPileUpSubtractor

//______________________________________________________________________________
void TrackPileUpSubtractor::Streamer(TBuffer &R__b)
{
   // Stream an object of class TrackPileUpSubtractor.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TrackPileUpSubtractor::Class(),this);
   } else {
      R__b.WriteClassBuffer(TrackPileUpSubtractor::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TrackPileUpSubtractor(void *p) {
      return  p ? new(p) ::TrackPileUpSubtractor : new ::TrackPileUpSubtractor;
   }
   static void *newArray_TrackPileUpSubtractor(Long_t nElements, void *p) {
      return p ? new(p) ::TrackPileUpSubtractor[nElements] : new ::TrackPileUpSubtractor[nElements];
   }
   // Wrapper around operator delete
   static void delete_TrackPileUpSubtractor(void *p) {
      delete ((::TrackPileUpSubtractor*)p);
   }
   static void deleteArray_TrackPileUpSubtractor(void *p) {
      delete [] ((::TrackPileUpSubtractor*)p);
   }
   static void destruct_TrackPileUpSubtractor(void *p) {
      typedef ::TrackPileUpSubtractor current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TrackPileUpSubtractor

//______________________________________________________________________________
void TaggingParticlesSkimmer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TaggingParticlesSkimmer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TaggingParticlesSkimmer::Class(),this);
   } else {
      R__b.WriteClassBuffer(TaggingParticlesSkimmer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TaggingParticlesSkimmer(void *p) {
      return  p ? new(p) ::TaggingParticlesSkimmer : new ::TaggingParticlesSkimmer;
   }
   static void *newArray_TaggingParticlesSkimmer(Long_t nElements, void *p) {
      return p ? new(p) ::TaggingParticlesSkimmer[nElements] : new ::TaggingParticlesSkimmer[nElements];
   }
   // Wrapper around operator delete
   static void delete_TaggingParticlesSkimmer(void *p) {
      delete ((::TaggingParticlesSkimmer*)p);
   }
   static void deleteArray_TaggingParticlesSkimmer(void *p) {
      delete [] ((::TaggingParticlesSkimmer*)p);
   }
   static void destruct_TaggingParticlesSkimmer(void *p) {
      typedef ::TaggingParticlesSkimmer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TaggingParticlesSkimmer

//______________________________________________________________________________
void PileUpJetID::Streamer(TBuffer &R__b)
{
   // Stream an object of class PileUpJetID.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PileUpJetID::Class(),this);
   } else {
      R__b.WriteClassBuffer(PileUpJetID::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PileUpJetID(void *p) {
      return  p ? new(p) ::PileUpJetID : new ::PileUpJetID;
   }
   static void *newArray_PileUpJetID(Long_t nElements, void *p) {
      return p ? new(p) ::PileUpJetID[nElements] : new ::PileUpJetID[nElements];
   }
   // Wrapper around operator delete
   static void delete_PileUpJetID(void *p) {
      delete ((::PileUpJetID*)p);
   }
   static void deleteArray_PileUpJetID(void *p) {
      delete [] ((::PileUpJetID*)p);
   }
   static void destruct_PileUpJetID(void *p) {
      typedef ::PileUpJetID current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PileUpJetID

//______________________________________________________________________________
void ConstituentFilter::Streamer(TBuffer &R__b)
{
   // Stream an object of class ConstituentFilter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ConstituentFilter::Class(),this);
   } else {
      R__b.WriteClassBuffer(ConstituentFilter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ConstituentFilter(void *p) {
      return  p ? new(p) ::ConstituentFilter : new ::ConstituentFilter;
   }
   static void *newArray_ConstituentFilter(Long_t nElements, void *p) {
      return p ? new(p) ::ConstituentFilter[nElements] : new ::ConstituentFilter[nElements];
   }
   // Wrapper around operator delete
   static void delete_ConstituentFilter(void *p) {
      delete ((::ConstituentFilter*)p);
   }
   static void deleteArray_ConstituentFilter(void *p) {
      delete [] ((::ConstituentFilter*)p);
   }
   static void destruct_ConstituentFilter(void *p) {
      typedef ::ConstituentFilter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ConstituentFilter

//______________________________________________________________________________
void StatusPidFilter::Streamer(TBuffer &R__b)
{
   // Stream an object of class StatusPidFilter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StatusPidFilter::Class(),this);
   } else {
      R__b.WriteClassBuffer(StatusPidFilter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StatusPidFilter(void *p) {
      return  p ? new(p) ::StatusPidFilter : new ::StatusPidFilter;
   }
   static void *newArray_StatusPidFilter(Long_t nElements, void *p) {
      return p ? new(p) ::StatusPidFilter[nElements] : new ::StatusPidFilter[nElements];
   }
   // Wrapper around operator delete
   static void delete_StatusPidFilter(void *p) {
      delete ((::StatusPidFilter*)p);
   }
   static void deleteArray_StatusPidFilter(void *p) {
      delete [] ((::StatusPidFilter*)p);
   }
   static void destruct_StatusPidFilter(void *p) {
      typedef ::StatusPidFilter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::StatusPidFilter

//______________________________________________________________________________
void PdgCodeFilter::Streamer(TBuffer &R__b)
{
   // Stream an object of class PdgCodeFilter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PdgCodeFilter::Class(),this);
   } else {
      R__b.WriteClassBuffer(PdgCodeFilter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PdgCodeFilter(void *p) {
      return  p ? new(p) ::PdgCodeFilter : new ::PdgCodeFilter;
   }
   static void *newArray_PdgCodeFilter(Long_t nElements, void *p) {
      return p ? new(p) ::PdgCodeFilter[nElements] : new ::PdgCodeFilter[nElements];
   }
   // Wrapper around operator delete
   static void delete_PdgCodeFilter(void *p) {
      delete ((::PdgCodeFilter*)p);
   }
   static void deleteArray_PdgCodeFilter(void *p) {
      delete [] ((::PdgCodeFilter*)p);
   }
   static void destruct_PdgCodeFilter(void *p) {
      typedef ::PdgCodeFilter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PdgCodeFilter

//______________________________________________________________________________
void Cloner::Streamer(TBuffer &R__b)
{
   // Stream an object of class Cloner.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Cloner::Class(),this);
   } else {
      R__b.WriteClassBuffer(Cloner::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Cloner(void *p) {
      return  p ? new(p) ::Cloner : new ::Cloner;
   }
   static void *newArray_Cloner(Long_t nElements, void *p) {
      return p ? new(p) ::Cloner[nElements] : new ::Cloner[nElements];
   }
   // Wrapper around operator delete
   static void delete_Cloner(void *p) {
      delete ((::Cloner*)p);
   }
   static void deleteArray_Cloner(void *p) {
      delete [] ((::Cloner*)p);
   }
   static void destruct_Cloner(void *p) {
      typedef ::Cloner current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Cloner

//______________________________________________________________________________
void Weighter::Streamer(TBuffer &R__b)
{
   // Stream an object of class Weighter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Weighter::Class(),this);
   } else {
      R__b.WriteClassBuffer(Weighter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Weighter(void *p) {
      return  p ? new(p) ::Weighter : new ::Weighter;
   }
   static void *newArray_Weighter(Long_t nElements, void *p) {
      return p ? new(p) ::Weighter[nElements] : new ::Weighter[nElements];
   }
   // Wrapper around operator delete
   static void delete_Weighter(void *p) {
      delete ((::Weighter*)p);
   }
   static void deleteArray_Weighter(void *p) {
      delete [] ((::Weighter*)p);
   }
   static void destruct_Weighter(void *p) {
      typedef ::Weighter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Weighter

//______________________________________________________________________________
void Hector::Streamer(TBuffer &R__b)
{
   // Stream an object of class Hector.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Hector::Class(),this);
   } else {
      R__b.WriteClassBuffer(Hector::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Hector(void *p) {
      return  p ? new(p) ::Hector : new ::Hector;
   }
   static void *newArray_Hector(Long_t nElements, void *p) {
      return p ? new(p) ::Hector[nElements] : new ::Hector[nElements];
   }
   // Wrapper around operator delete
   static void delete_Hector(void *p) {
      delete ((::Hector*)p);
   }
   static void deleteArray_Hector(void *p) {
      delete [] ((::Hector*)p);
   }
   static void destruct_Hector(void *p) {
      typedef ::Hector current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Hector

//______________________________________________________________________________
void JetFlavorAssociation::Streamer(TBuffer &R__b)
{
   // Stream an object of class JetFlavorAssociation.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(JetFlavorAssociation::Class(),this);
   } else {
      R__b.WriteClassBuffer(JetFlavorAssociation::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_JetFlavorAssociation(void *p) {
      return  p ? new(p) ::JetFlavorAssociation : new ::JetFlavorAssociation;
   }
   static void *newArray_JetFlavorAssociation(Long_t nElements, void *p) {
      return p ? new(p) ::JetFlavorAssociation[nElements] : new ::JetFlavorAssociation[nElements];
   }
   // Wrapper around operator delete
   static void delete_JetFlavorAssociation(void *p) {
      delete ((::JetFlavorAssociation*)p);
   }
   static void deleteArray_JetFlavorAssociation(void *p) {
      delete [] ((::JetFlavorAssociation*)p);
   }
   static void destruct_JetFlavorAssociation(void *p) {
      typedef ::JetFlavorAssociation current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::JetFlavorAssociation

//______________________________________________________________________________
void JetFakeParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class JetFakeParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(JetFakeParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(JetFakeParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_JetFakeParticle(void *p) {
      return  p ? new(p) ::JetFakeParticle : new ::JetFakeParticle;
   }
   static void *newArray_JetFakeParticle(Long_t nElements, void *p) {
      return p ? new(p) ::JetFakeParticle[nElements] : new ::JetFakeParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_JetFakeParticle(void *p) {
      delete ((::JetFakeParticle*)p);
   }
   static void deleteArray_JetFakeParticle(void *p) {
      delete [] ((::JetFakeParticle*)p);
   }
   static void destruct_JetFakeParticle(void *p) {
      typedef ::JetFakeParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::JetFakeParticle

//______________________________________________________________________________
void ExampleModule::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExampleModule.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ExampleModule::Class(),this);
   } else {
      R__b.WriteClassBuffer(ExampleModule::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExampleModule(void *p) {
      return  p ? new(p) ::ExampleModule : new ::ExampleModule;
   }
   static void *newArray_ExampleModule(Long_t nElements, void *p) {
      return p ? new(p) ::ExampleModule[nElements] : new ::ExampleModule[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExampleModule(void *p) {
      delete ((::ExampleModule*)p);
   }
   static void deleteArray_ExampleModule(void *p) {
      delete [] ((::ExampleModule*)p);
   }
   static void destruct_ExampleModule(void *p) {
      typedef ::ExampleModule current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExampleModule

namespace ROOT {
   static TClass *vectorlEvectorlEdoublegRmUgR_Dictionary();
   static void vectorlEvectorlEdoublegRmUgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEdoublegRmUgR(void *p = 0);
   static void *newArray_vectorlEvectorlEdoublegRmUgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEdoublegRmUgR(void *p);
   static void deleteArray_vectorlEvectorlEdoublegRmUgR(void *p);
   static void destruct_vectorlEvectorlEdoublegRmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<double>*>*)
   {
      vector<vector<double>*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<double>*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<double>*>", -2, "vector", 214,
                  typeid(vector<vector<double>*>), DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEdoublegRmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<double>*>) );
      instance.SetNew(&new_vectorlEvectorlEdoublegRmUgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEdoublegRmUgR);
      instance.SetDelete(&delete_vectorlEvectorlEdoublegRmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEdoublegRmUgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEdoublegRmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<double>*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<double>*>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEdoublegRmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<double>*>*)0x0)->GetClass();
      vectorlEvectorlEdoublegRmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEdoublegRmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEdoublegRmUgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<vector<double>*> : new vector<vector<double>*>;
   }
   static void *newArray_vectorlEvectorlEdoublegRmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<vector<double>*>[nElements] : new vector<vector<double>*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEdoublegRmUgR(void *p) {
      delete ((vector<vector<double>*>*)p);
   }
   static void deleteArray_vectorlEvectorlEdoublegRmUgR(void *p) {
      delete [] ((vector<vector<double>*>*)p);
   }
   static void destruct_vectorlEvectorlEdoublegRmUgR(void *p) {
      typedef vector<vector<double>*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<double>*>

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 214,
                  typeid(vector<int>), DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

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
         instance("vector<double>", -2, "vector", 214,
                  typeid(vector<double>), DefineBehavior(ptr, ptr),
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
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

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
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
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

namespace ROOT {
   static TClass *vectorlELong64_tgR_Dictionary();
   static void vectorlELong64_tgR_TClassManip(TClass*);
   static void *new_vectorlELong64_tgR(void *p = 0);
   static void *newArray_vectorlELong64_tgR(Long_t size, void *p);
   static void delete_vectorlELong64_tgR(void *p);
   static void deleteArray_vectorlELong64_tgR(void *p);
   static void destruct_vectorlELong64_tgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Long64_t>*)
   {
      vector<Long64_t> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Long64_t>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Long64_t>", -2, "vector", 214,
                  typeid(vector<Long64_t>), DefineBehavior(ptr, ptr),
                  &vectorlELong64_tgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Long64_t>) );
      instance.SetNew(&new_vectorlELong64_tgR);
      instance.SetNewArray(&newArray_vectorlELong64_tgR);
      instance.SetDelete(&delete_vectorlELong64_tgR);
      instance.SetDeleteArray(&deleteArray_vectorlELong64_tgR);
      instance.SetDestructor(&destruct_vectorlELong64_tgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Long64_t> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<Long64_t>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELong64_tgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Long64_t>*)0x0)->GetClass();
      vectorlELong64_tgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELong64_tgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELong64_tgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<Long64_t> : new vector<Long64_t>;
   }
   static void *newArray_vectorlELong64_tgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<Long64_t>[nElements] : new vector<Long64_t>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELong64_tgR(void *p) {
      delete ((vector<Long64_t>*)p);
   }
   static void deleteArray_vectorlELong64_tgR(void *p) {
      delete [] ((vector<Long64_t>*)p);
   }
   static void destruct_vectorlELong64_tgR(void *p) {
      typedef vector<Long64_t> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Long64_t>

namespace ROOT {
   static TClass *dequelEdoublegR_Dictionary();
   static void dequelEdoublegR_TClassManip(TClass*);
   static void *new_dequelEdoublegR(void *p = 0);
   static void *newArray_dequelEdoublegR(Long_t size, void *p);
   static void delete_dequelEdoublegR(void *p);
   static void deleteArray_dequelEdoublegR(void *p);
   static void destruct_dequelEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const deque<double>*)
   {
      deque<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(deque<double>));
      static ::ROOT::TGenericClassInfo 
         instance("deque<double>", -2, "deque", 735,
                  typeid(deque<double>), DefineBehavior(ptr, ptr),
                  &dequelEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(deque<double>) );
      instance.SetNew(&new_dequelEdoublegR);
      instance.SetNewArray(&newArray_dequelEdoublegR);
      instance.SetDelete(&delete_dequelEdoublegR);
      instance.SetDeleteArray(&deleteArray_dequelEdoublegR);
      instance.SetDestructor(&destruct_dequelEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< deque<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const deque<double>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *dequelEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const deque<double>*)0x0)->GetClass();
      dequelEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void dequelEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_dequelEdoublegR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) deque<double> : new deque<double>;
   }
   static void *newArray_dequelEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) deque<double>[nElements] : new deque<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_dequelEdoublegR(void *p) {
      delete ((deque<double>*)p);
   }
   static void deleteArray_dequelEdoublegR(void *p) {
      delete [] ((deque<double>*)p);
   }
   static void destruct_dequelEdoublegR(void *p) {
      typedef deque<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class deque<double>

namespace {
  void TriggerDictionaryInitialization_ModulesDict_Impl() {
    static const char* headers[] = {
0
    };
    static const char* includePaths[] = {
"external",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.00-odfocd7/include",
"/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/Madgraph/MG5_aMC_v2_3_3/Delphes/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  Delphes;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  AngularSmearing;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  PhotonConversions;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  ParticlePropagator;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  Efficiency;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  IdentificationMap;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  EnergySmearing;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  MomentumSmearing;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  ImpactParameterSmearing;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  TimeSmearing;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  SimpleCalorimeter;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  Calorimeter;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  Isolation;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  EnergyScale;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  UniqueObjectFinder;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  TrackCountingBTagging;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  BTagging;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  TauTagging;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  TreeWriter;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  Merger;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  LeptonDressing;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  PileUpMerger;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  JetPileUpSubtractor;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  TrackPileUpSubtractor;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  TaggingParticlesSkimmer;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  PileUpJetID;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  ConstituentFilter;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  StatusPidFilter;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  PdgCodeFilter;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  Cloner;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  Weighter;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  Hector;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  JetFlavorAssociation;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  JetFakeParticle;
class __attribute__((annotate("$clingAutoload$modules/ModulesLinkDef.h")))  ExampleModule;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"AngularSmearing", payloadCode, "@",
"BTagging", payloadCode, "@",
"Calorimeter", payloadCode, "@",
"Cloner", payloadCode, "@",
"ConstituentFilter", payloadCode, "@",
"Delphes", payloadCode, "@",
"Efficiency", payloadCode, "@",
"EnergyScale", payloadCode, "@",
"EnergySmearing", payloadCode, "@",
"ExampleModule", payloadCode, "@",
"Hector", payloadCode, "@",
"IdentificationMap", payloadCode, "@",
"ImpactParameterSmearing", payloadCode, "@",
"Isolation", payloadCode, "@",
"JetFakeParticle", payloadCode, "@",
"JetFlavorAssociation", payloadCode, "@",
"JetPileUpSubtractor", payloadCode, "@",
"LeptonDressing", payloadCode, "@",
"Merger", payloadCode, "@",
"MomentumSmearing", payloadCode, "@",
"ParticlePropagator", payloadCode, "@",
"PdgCodeFilter", payloadCode, "@",
"PhotonConversions", payloadCode, "@",
"PileUpJetID", payloadCode, "@",
"PileUpMerger", payloadCode, "@",
"SimpleCalorimeter", payloadCode, "@",
"StatusPidFilter", payloadCode, "@",
"TaggingParticlesSkimmer", payloadCode, "@",
"TauTagging", payloadCode, "@",
"TimeSmearing", payloadCode, "@",
"TrackCountingBTagging", payloadCode, "@",
"TrackPileUpSubtractor", payloadCode, "@",
"TreeWriter", payloadCode, "@",
"UniqueObjectFinder", payloadCode, "@",
"Weighter", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("ModulesDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_ModulesDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_ModulesDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_ModulesDict() {
  TriggerDictionaryInitialization_ModulesDict_Impl();
}
