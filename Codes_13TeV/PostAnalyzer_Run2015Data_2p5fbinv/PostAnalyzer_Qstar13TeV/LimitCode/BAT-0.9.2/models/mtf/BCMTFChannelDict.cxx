//
// File generated by /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.10-cms7//bin/rootcint at Sun Feb 14 14:58:36 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME BCMTFChannelDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "BCMTFChannelDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void BCMTFChannel_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void BCMTFChannel_Dictionary();
   static void delete_BCMTFChannel(void *p);
   static void deleteArray_BCMTFChannel(void *p);
   static void destruct_BCMTFChannel(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BCMTFChannel*)
   {
      ::BCMTFChannel *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BCMTFChannel),0);
      static ::ROOT::TGenericClassInfo 
         instance("BCMTFChannel", "./BCMTFChannel.h", 33,
                  typeid(::BCMTFChannel), DefineBehavior(ptr, ptr),
                  0, &BCMTFChannel_Dictionary, isa_proxy, 0,
                  sizeof(::BCMTFChannel) );
      instance.SetDelete(&delete_BCMTFChannel);
      instance.SetDeleteArray(&deleteArray_BCMTFChannel);
      instance.SetDestructor(&destruct_BCMTFChannel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BCMTFChannel*)
   {
      return GenerateInitInstanceLocal((::BCMTFChannel*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BCMTFChannel*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void BCMTFChannel_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::BCMTFChannel*)0x0)->GetClass();
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_BCMTFChannel(void *p) {
      delete ((::BCMTFChannel*)p);
   }
   static void deleteArray_BCMTFChannel(void *p) {
      delete [] ((::BCMTFChannel*)p);
   }
   static void destruct_BCMTFChannel(void *p) {
      typedef ::BCMTFChannel current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BCMTFChannel

/********************************************************
* BCMTFChannelDict.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableBCMTFChannelDict();

extern "C" void G__set_cpp_environmentBCMTFChannelDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("BCMTFChannel.h");
  G__cpp_reset_tagtableBCMTFChannelDict();
}
#include <new>
extern "C" int G__cpp_dllrevBCMTFChannelDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* BCMTFChannel */
static int G__BCMTFChannelDict_218_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCMTFChannel* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new BCMTFChannel((const char*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) BCMTFChannel((const char*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFChannel));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         string* pobj;
         string xobj = ((BCMTFChannel*) G__getstructoffset())->GetName();
         pobj = new string(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFChannel*) G__getstructoffset())->GetData());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFChannel*) G__getstructoffset())->GetTemplate((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFChannel*) G__getstructoffset())->GetSystematicVariation((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 103, (long) ((BCMTFChannel*) G__getstructoffset())->GetFlagChannelActive());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFChannel*) G__getstructoffset())->GetHistUncertaintyBandExpectation());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFChannel*) G__getstructoffset())->GetHistUncertaintyBandPoisson());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((BCMTFChannel*) G__getstructoffset())->GetRangeYMin());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((BCMTFChannel*) G__getstructoffset())->GetRangeYMax());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->SetName((const char*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->SetData((BCMTFTemplate*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->SetHistUncertaintyBandExpectation((TH2D*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->SetHistUncertaintyBandPoisson((TH2D*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->SetFlagChannelActive((bool) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->SetRangeY((double) G__double(libp->para[0]), (double) G__double(libp->para[1]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->AddTemplate((BCMTFTemplate*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->AddSystematicVariation((BCMTFSystematicVariation*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->CalculateHistUncertaintyBandPoisson();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFChannel*) G__getstructoffset())->CalculateUncertaintyBandPoisson((double) G__double(libp->para[0]), (double) G__double(libp->para[1])
, (int) G__int(libp->para[2])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->PrintTemplates((const char*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->PrintTemplate((int) G__int(libp->para[0]), (const char*) G__int(libp->para[1]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->PrintHistUncertaintyBandExpectation((const char*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_24(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->PrintHistUncertaintyBandPoisson((const char*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFChannelDict_218_0_25(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFChannel*) G__getstructoffset())->PrintHistCumulativeUncertaintyBandPoisson((const char*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__BCMTFChannelDict_218_0_26(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   BCMTFChannel* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new BCMTFChannel(*(BCMTFChannel*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFChannel));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef BCMTFChannel G__TBCMTFChannel;
static int G__BCMTFChannelDict_218_0_27(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 0
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (BCMTFChannel*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((BCMTFChannel*) (soff+(sizeof(BCMTFChannel)*i)))->~G__TBCMTFChannel();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (BCMTFChannel*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((BCMTFChannel*) (soff))->~G__TBCMTFChannel();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__BCMTFChannelDict_218_0_28(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCMTFChannel* dest = (BCMTFChannel*) G__getstructoffset();
   *dest = *(BCMTFChannel*) libp->para[0].ref;
   const BCMTFChannel& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* BCMTFChannel */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncBCMTFChannelDict {
 public:
  G__Sizep2memfuncBCMTFChannelDict(): p(&G__Sizep2memfuncBCMTFChannelDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncBCMTFChannelDict::*p)();
};

size_t G__get_sizep2memfuncBCMTFChannelDict()
{
  G__Sizep2memfuncBCMTFChannelDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceBCMTFChannelDict() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableBCMTFChannelDict() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Float_t>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_TVectorTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Double_t>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_TVectorTlEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<BCMTFTemplate*>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<BCMTFSystematicVariation*>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* BCMTFChannel */
static void G__setup_memvarBCMTFChannel(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFChannel));
   { BCMTFChannel *p; p=(BCMTFChannel*)0x1000; if (p) { }
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_string),-1,-1,4,"fName=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFTemplate),-1,-1,4,"fData=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,-1,-1,4,"fRangeYMin=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,-1,-1,4,"fRangeYMax=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgR),G__defined_typename("vector<BCMTFTemplate*>"),-1,4,"fTemplateContainer=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_vectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgR),G__defined_typename("vector<BCMTFSystematicVariation*>"),-1,4,"fSystematicVariationContainer=",0,(char*)NULL);
   G__memvar_setup((void*)0,103,0,0,-1,-1,-1,4,"fFlagChannelActive=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_TH2D),-1,-1,4,"fHistUncertaintyBandExpectation=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BCMTFChannelDictLN_TH2D),-1,-1,4,"fHistUncertaintyBandPoisson=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarBCMTFChannelDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncBCMTFChannel(void) {
   /* BCMTFChannel */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFChannel));
   G__memfunc_setup("BCMTFChannel",1061,G__BCMTFChannelDict_218_0_1, 105, G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFChannel), -1, 0, 1, 1, 1, 0, "C - - 10 - name", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetName",673,G__BCMTFChannelDict_218_0_2, 117, G__get_linked_tagnum(&G__BCMTFChannelDictLN_string), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetData",666,G__BCMTFChannelDict_218_0_3, 85, G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFTemplate), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetTemplate",1116,G__BCMTFChannelDict_218_0_4, 85, G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFTemplate), -1, 0, 1, 1, 1, 0, "i - - 0 - index", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetSystematicVariation",2291,G__BCMTFChannelDict_218_0_5, 85, G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFSystematicVariation), -1, 0, 1, 1, 1, 0, "i - - 0 - index", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetFlagChannelActive",1967,G__BCMTFChannelDict_218_0_6, 103, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetHistUncertaintyBandExpectation",3399,G__BCMTFChannelDict_218_0_7, 85, G__get_linked_tagnum(&G__BCMTFChannelDictLN_TH2D), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetHistUncertaintyBandPoisson",2990,G__BCMTFChannelDict_218_0_8, 85, G__get_linked_tagnum(&G__BCMTFChannelDictLN_TH2D), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetRangeYMin",1162,G__BCMTFChannelDict_218_0_9, 100, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetRangeYMax",1164,G__BCMTFChannelDict_218_0_10, 100, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetName",685,G__BCMTFChannelDict_218_0_11, 121, -1, -1, 0, 1, 1, 1, 0, "C - - 10 - name", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetData",678,G__BCMTFChannelDict_218_0_12, 121, -1, -1, 0, 1, 1, 1, 0, "U 'BCMTFTemplate' - 0 - bctemplate", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetHistUncertaintyBandExpectation",3411,G__BCMTFChannelDict_218_0_13, 121, -1, -1, 0, 1, 1, 1, 0, "U 'TH2D' - 0 - hist", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetHistUncertaintyBandPoisson",3002,G__BCMTFChannelDict_218_0_14, 121, -1, -1, 0, 1, 1, 1, 0, "U 'TH2D' - 0 - hist", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetFlagChannelActive",1979,G__BCMTFChannelDict_218_0_15, 121, -1, -1, 0, 1, 1, 1, 0, "g - - 0 - flag", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetRangeY",882,G__BCMTFChannelDict_218_0_16, 121, -1, -1, 0, 2, 1, 1, 0, 
"d - - 0 - min d - - 0 - max", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("AddTemplate",1093,G__BCMTFChannelDict_218_0_17, 121, -1, -1, 0, 1, 1, 1, 0, "U 'BCMTFTemplate' - 0 - bctemplate", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("AddSystematicVariation",2268,G__BCMTFChannelDict_218_0_18, 121, -1, -1, 0, 1, 1, 1, 0, "U 'BCMTFSystematicVariation' - 0 - variation", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("CalculateHistUncertaintyBandPoisson",3612,G__BCMTFChannelDict_218_0_19, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("CalculateUncertaintyBandPoisson",3204,G__BCMTFChannelDict_218_0_20, 85, G__get_linked_tagnum(&G__BCMTFChannelDictLN_TH1D), -1, 0, 3, 1, 1, 0, 
"d - - 0 - minimum d - - 0 - maximumm "
"i - - 0 - color", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("PrintTemplates",1468,G__BCMTFChannelDict_218_0_21, 121, -1, -1, 0, 1, 1, 1, 0, "C - - 10 - filename", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("PrintTemplate",1353,G__BCMTFChannelDict_218_0_22, 121, -1, -1, 0, 2, 1, 1, 0, 
"i - - 0 - index C - - 10 - filename", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("PrintHistUncertaintyBandExpectation",3636,G__BCMTFChannelDict_218_0_23, 121, -1, -1, 0, 1, 1, 1, 0, "C - - 10 - filename", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("PrintHistUncertaintyBandPoisson",3227,G__BCMTFChannelDict_218_0_24, 121, -1, -1, 0, 1, 1, 1, 0, "C - - 10 - filename", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("PrintHistCumulativeUncertaintyBandPoisson",4282,G__BCMTFChannelDict_218_0_25, 121, -1, -1, 0, 1, 1, 1, 0, "C - - 10 - filename", (char*)NULL, (void*) NULL, 0);
   // automatic copy constructor
   G__memfunc_setup("BCMTFChannel", 1061, G__BCMTFChannelDict_218_0_26, (int) ('i'), G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFChannel), -1, 0, 1, 1, 1, 0, "u 'BCMTFChannel' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~BCMTFChannel", 1187, G__BCMTFChannelDict_218_0_27, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__BCMTFChannelDict_218_0_28, (int) ('u'), G__get_linked_tagnum(&G__BCMTFChannelDictLN_BCMTFChannel), -1, 1, 1, 1, 1, 0, "u 'BCMTFChannel' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncBCMTFChannelDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalBCMTFChannelDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {
}

static void G__cpp_setup_func14() {
}

static void G__cpp_setup_func15() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcBCMTFChannelDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
  G__cpp_setup_func14();
  G__cpp_setup_func15();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__BCMTFChannelDictLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_TVectorTlEfloatgR = { "TVectorT<float>" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_TVectorTlEdoublegR = { "TVectorT<double>" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_TH1D = { "TH1D" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_BCMTFTemplate = { "BCMTFTemplate" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_BCMTFSystematicVariation = { "BCMTFSystematicVariation" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_TH2D = { "TH2D" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_BCMTFChannel = { "BCMTFChannel" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_vectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgR = { "vector<BCMTFTemplate*,allocator<BCMTFTemplate*> >" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<BCMTFTemplate*,allocator<BCMTFTemplate*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_vectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgR = { "vector<BCMTFSystematicVariation*,allocator<BCMTFSystematicVariation*> >" , 99 , -1 };
G__linked_taginfo G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<BCMTFSystematicVariation*,allocator<BCMTFSystematicVariation*> >::iterator>" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableBCMTFChannelDict() {
  G__BCMTFChannelDictLN_string.tagnum = -1 ;
  G__BCMTFChannelDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCMTFChannelDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__BCMTFChannelDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCMTFChannelDictLN_TVectorTlEfloatgR.tagnum = -1 ;
  G__BCMTFChannelDictLN_TVectorTlEdoublegR.tagnum = -1 ;
  G__BCMTFChannelDictLN_TH1D.tagnum = -1 ;
  G__BCMTFChannelDictLN_BCMTFTemplate.tagnum = -1 ;
  G__BCMTFChannelDictLN_BCMTFSystematicVariation.tagnum = -1 ;
  G__BCMTFChannelDictLN_TH2D.tagnum = -1 ;
  G__BCMTFChannelDictLN_BCMTFChannel.tagnum = -1 ;
  G__BCMTFChannelDictLN_vectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgR.tagnum = -1 ;
  G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCMTFChannelDictLN_vectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgR.tagnum = -1 ;
  G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgRcLcLiteratorgR.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableBCMTFChannelDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_string);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_TVectorTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_TVectorTlEdoublegR);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_TH1D);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_BCMTFTemplate);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_BCMTFSystematicVariation);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_TH2D);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_BCMTFChannel),sizeof(BCMTFChannel),-1,33792,(char*)NULL,G__setup_memvarBCMTFChannel,G__setup_memfuncBCMTFChannel);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_vectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFTemplatemUcOallocatorlEBCMTFTemplatemUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_vectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFChannelDictLN_reverse_iteratorlEvectorlEBCMTFSystematicVariationmUcOallocatorlEBCMTFSystematicVariationmUgRsPgRcLcLiteratorgR);
}
extern "C" void G__cpp_setupBCMTFChannelDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupBCMTFChannelDict()");
  G__set_cpp_environmentBCMTFChannelDict();
  G__cpp_setup_tagtableBCMTFChannelDict();

  G__cpp_setup_inheritanceBCMTFChannelDict();

  G__cpp_setup_typetableBCMTFChannelDict();

  G__cpp_setup_memvarBCMTFChannelDict();

  G__cpp_setup_memfuncBCMTFChannelDict();
  G__cpp_setup_globalBCMTFChannelDict();
  G__cpp_setup_funcBCMTFChannelDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncBCMTFChannelDict();
  return;
}
class G__cpp_setup_initBCMTFChannelDict {
  public:
    G__cpp_setup_initBCMTFChannelDict() { G__add_setup_func("BCMTFChannelDict",(G__incsetup)(&G__cpp_setupBCMTFChannelDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initBCMTFChannelDict() { G__remove_setup_func("BCMTFChannelDict"); }
};
G__cpp_setup_initBCMTFChannelDict G__cpp_setup_initializerBCMTFChannelDict;

