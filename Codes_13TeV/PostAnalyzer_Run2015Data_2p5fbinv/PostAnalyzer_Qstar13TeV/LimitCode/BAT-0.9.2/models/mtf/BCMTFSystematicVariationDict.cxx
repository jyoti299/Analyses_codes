//
// File generated by /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.10-cms7//bin/rootcint at Sun Feb 14 14:58:41 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME BCMTFSystematicVariationDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "BCMTFSystematicVariationDict.h"

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
   void BCMTFSystematicVariation_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void BCMTFSystematicVariation_Dictionary();
   static void delete_BCMTFSystematicVariation(void *p);
   static void deleteArray_BCMTFSystematicVariation(void *p);
   static void destruct_BCMTFSystematicVariation(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BCMTFSystematicVariation*)
   {
      ::BCMTFSystematicVariation *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BCMTFSystematicVariation),0);
      static ::ROOT::TGenericClassInfo 
         instance("BCMTFSystematicVariation", "./BCMTFSystematicVariation.h", 31,
                  typeid(::BCMTFSystematicVariation), DefineBehavior(ptr, ptr),
                  0, &BCMTFSystematicVariation_Dictionary, isa_proxy, 0,
                  sizeof(::BCMTFSystematicVariation) );
      instance.SetDelete(&delete_BCMTFSystematicVariation);
      instance.SetDeleteArray(&deleteArray_BCMTFSystematicVariation);
      instance.SetDestructor(&destruct_BCMTFSystematicVariation);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BCMTFSystematicVariation*)
   {
      return GenerateInitInstanceLocal((::BCMTFSystematicVariation*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BCMTFSystematicVariation*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void BCMTFSystematicVariation_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::BCMTFSystematicVariation*)0x0)->GetClass();
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_BCMTFSystematicVariation(void *p) {
      delete ((::BCMTFSystematicVariation*)p);
   }
   static void deleteArray_BCMTFSystematicVariation(void *p) {
      delete [] ((::BCMTFSystematicVariation*)p);
   }
   static void destruct_BCMTFSystematicVariation(void *p) {
      typedef ::BCMTFSystematicVariation current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BCMTFSystematicVariation

/********************************************************
* BCMTFSystematicVariationDict.cxx
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

extern "C" void G__cpp_reset_tagtableBCMTFSystematicVariationDict();

extern "C" void G__set_cpp_environmentBCMTFSystematicVariationDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("BCMTFSystematicVariation.h");
  G__cpp_reset_tagtableBCMTFSystematicVariationDict();
}
#include <new>
extern "C" int G__cpp_dllrevBCMTFSystematicVariationDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* BCMTFSystematicVariation */
static int G__BCMTFSystematicVariationDict_169_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCMTFSystematicVariation* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 3
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new BCMTFSystematicVariation(
(const char*) G__int(libp->para[0]), (const char*) G__int(libp->para[1])
, (int) G__int(libp->para[2]));
   } else {
     p = new((void*) gvp) BCMTFSystematicVariation(
(const char*) G__int(libp->para[0]), (const char*) G__int(libp->para[1])
, (int) G__int(libp->para[2]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFSystematicVariationDict_169_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFSystematicVariation*) G__getstructoffset())->GetHistogramUp((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFSystematicVariationDict_169_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFSystematicVariation*) G__getstructoffset())->GetHistogramDown((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFSystematicVariationDict_169_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFSystematicVariation*) G__getstructoffset())->SetHistogramUp((int) G__int(libp->para[0]), (TH1D*) G__int(libp->para[1]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFSystematicVariationDict_169_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFSystematicVariation*) G__getstructoffset())->SetHistogramDown((int) G__int(libp->para[0]), (TH1D*) G__int(libp->para[1]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFSystematicVariationDict_169_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFSystematicVariation*) G__getstructoffset())->SetHistograms((int) G__int(libp->para[0]), (TH1D*) G__int(libp->para[1])
, (TH1D*) G__int(libp->para[2]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFSystematicVariationDict_169_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFSystematicVariation*) G__getstructoffset())->AddHistogramUp((TH1D*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFSystematicVariationDict_169_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFSystematicVariation*) G__getstructoffset())->AddHistogramDown((TH1D*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFSystematicVariationDict_169_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFSystematicVariation*) G__getstructoffset())->AddHistograms((TH1D*) G__int(libp->para[0]), (TH1D*) G__int(libp->para[1]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__BCMTFSystematicVariationDict_169_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   BCMTFSystematicVariation* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new BCMTFSystematicVariation(*(BCMTFSystematicVariation*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef BCMTFSystematicVariation G__TBCMTFSystematicVariation;
static int G__BCMTFSystematicVariationDict_169_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (BCMTFSystematicVariation*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((BCMTFSystematicVariation*) (soff+(sizeof(BCMTFSystematicVariation)*i)))->~G__TBCMTFSystematicVariation();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (BCMTFSystematicVariation*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((BCMTFSystematicVariation*) (soff))->~G__TBCMTFSystematicVariation();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__BCMTFSystematicVariationDict_169_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCMTFSystematicVariation* dest = (BCMTFSystematicVariation*) G__getstructoffset();
   *dest = *(BCMTFSystematicVariation*) libp->para[0].ref;
   const BCMTFSystematicVariation& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* BCMTFSystematicVariation */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncBCMTFSystematicVariationDict {
 public:
  G__Sizep2memfuncBCMTFSystematicVariationDict(): p(&G__Sizep2memfuncBCMTFSystematicVariationDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncBCMTFSystematicVariationDict::*p)();
};

size_t G__get_sizep2memfuncBCMTFSystematicVariationDict()
{
  G__Sizep2memfuncBCMTFSystematicVariationDict a;
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
extern "C" void G__cpp_setup_inheritanceBCMTFSystematicVariationDict() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableBCMTFSystematicVariationDict() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TH1D*>",117,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlETH1DmUcOallocatorlETH1DmUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETH1DmUcOallocatorlETH1DmUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlETH1DmUcOallocatorlETH1DmUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETH1DmUcOallocatorlETH1DmUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlETH1DmUcOallocatorlETH1DmUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* BCMTFSystematicVariation */
static void G__setup_memvarBCMTFSystematicVariation(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation));
   { BCMTFSystematicVariation *p; p=(BCMTFSystematicVariation*)0x1000; if (p) { }
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlETH1DmUcOallocatorlETH1DmUgRsPgR),G__defined_typename("vector<TH1D*>"),-1,4,"fHistogramUpContainer=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_vectorlETH1DmUcOallocatorlETH1DmUgRsPgR),G__defined_typename("vector<TH1D*>"),-1,4,"fHistogramDownContainer=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_string),-1,-1,4,"fChannelName=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_string),-1,-1,4,"fSystematicName=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarBCMTFSystematicVariationDict() {
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
static void G__setup_memfuncBCMTFSystematicVariation(void) {
   /* BCMTFSystematicVariation */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation));
   G__memfunc_setup("BCMTFSystematicVariation",2367,G__BCMTFSystematicVariationDict_169_0_1, 105, G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation), -1, 0, 3, 1, 1, 0, 
"C - - 10 - channelname C - - 10 - systematicname "
"i - - 0 - nprocesses", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetHistogramUp",1427,G__BCMTFSystematicVariationDict_169_0_2, 85, G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_TH1D), -1, 0, 1, 1, 1, 0, "i - - 0 - index", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetHistogramDown",1638,G__BCMTFSystematicVariationDict_169_0_3, 85, G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_TH1D), -1, 0, 1, 1, 1, 0, "i - - 0 - index", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetHistogramUp",1439,G__BCMTFSystematicVariationDict_169_0_4, 121, -1, -1, 0, 2, 1, 1, 0, 
"i - - 0 - index U 'TH1D' - 0 - hist", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetHistogramDown",1650,G__BCMTFSystematicVariationDict_169_0_5, 121, -1, -1, 0, 2, 1, 1, 0, 
"i - - 0 - index U 'TH1D' - 0 - hist", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetHistograms",1357,G__BCMTFSystematicVariationDict_169_0_6, 121, -1, -1, 0, 3, 1, 1, 0, 
"i - - 0 - index U 'TH1D' - 0 - hist_up "
"U 'TH1D' - 0 - hist_down", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("AddHistogramUp",1404,G__BCMTFSystematicVariationDict_169_0_7, 121, -1, -1, 0, 1, 1, 1, 0, "U 'TH1D' - 0 - hist", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("AddHistogramDown",1615,G__BCMTFSystematicVariationDict_169_0_8, 121, -1, -1, 0, 1, 1, 1, 0, "U 'TH1D' - 0 - hist", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("AddHistograms",1322,G__BCMTFSystematicVariationDict_169_0_9, 121, -1, -1, 0, 2, 1, 1, 0, 
"U 'TH1D' - 0 - hist_up U 'TH1D' - 0 - hist_down", (char*)NULL, (void*) NULL, 0);
   // automatic copy constructor
   G__memfunc_setup("BCMTFSystematicVariation", 2367, G__BCMTFSystematicVariationDict_169_0_10, (int) ('i'), 
G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation), -1, 0, 1, 1, 1, 0, "u 'BCMTFSystematicVariation' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~BCMTFSystematicVariation", 2493, G__BCMTFSystematicVariationDict_169_0_11, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__BCMTFSystematicVariationDict_169_0_12, (int) ('u'), G__get_linked_tagnum(&G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation), -1, 1, 1, 1, 1, 0, "u 'BCMTFSystematicVariation' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncBCMTFSystematicVariationDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalBCMTFSystematicVariationDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
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

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcBCMTFSystematicVariationDict() {
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
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__BCMTFSystematicVariationDictLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__BCMTFSystematicVariationDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCMTFSystematicVariationDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCMTFSystematicVariationDictLN_TH1D = { "TH1D" , 99 , -1 };
G__linked_taginfo G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation = { "BCMTFSystematicVariation" , 99 , -1 };
G__linked_taginfo G__BCMTFSystematicVariationDictLN_vectorlETH1DmUcOallocatorlETH1DmUgRsPgR = { "vector<TH1D*,allocator<TH1D*> >" , 99 , -1 };
G__linked_taginfo G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETH1DmUcOallocatorlETH1DmUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TH1D*,allocator<TH1D*> >::iterator>" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableBCMTFSystematicVariationDict() {
  G__BCMTFSystematicVariationDictLN_string.tagnum = -1 ;
  G__BCMTFSystematicVariationDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCMTFSystematicVariationDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCMTFSystematicVariationDictLN_TH1D.tagnum = -1 ;
  G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation.tagnum = -1 ;
  G__BCMTFSystematicVariationDictLN_vectorlETH1DmUcOallocatorlETH1DmUgRsPgR.tagnum = -1 ;
  G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETH1DmUcOallocatorlETH1DmUgRsPgRcLcLiteratorgR.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableBCMTFSystematicVariationDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__BCMTFSystematicVariationDictLN_string);
   G__get_linked_tagnum_fwd(&G__BCMTFSystematicVariationDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCMTFSystematicVariationDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCMTFSystematicVariationDictLN_TH1D);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__BCMTFSystematicVariationDictLN_BCMTFSystematicVariation),sizeof(BCMTFSystematicVariation),-1,33792,(char*)NULL,G__setup_memvarBCMTFSystematicVariation,G__setup_memfuncBCMTFSystematicVariation);
   G__get_linked_tagnum_fwd(&G__BCMTFSystematicVariationDictLN_vectorlETH1DmUcOallocatorlETH1DmUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFSystematicVariationDictLN_reverse_iteratorlEvectorlETH1DmUcOallocatorlETH1DmUgRsPgRcLcLiteratorgR);
}
extern "C" void G__cpp_setupBCMTFSystematicVariationDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupBCMTFSystematicVariationDict()");
  G__set_cpp_environmentBCMTFSystematicVariationDict();
  G__cpp_setup_tagtableBCMTFSystematicVariationDict();

  G__cpp_setup_inheritanceBCMTFSystematicVariationDict();

  G__cpp_setup_typetableBCMTFSystematicVariationDict();

  G__cpp_setup_memvarBCMTFSystematicVariationDict();

  G__cpp_setup_memfuncBCMTFSystematicVariationDict();
  G__cpp_setup_globalBCMTFSystematicVariationDict();
  G__cpp_setup_funcBCMTFSystematicVariationDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncBCMTFSystematicVariationDict();
  return;
}
class G__cpp_setup_initBCMTFSystematicVariationDict {
  public:
    G__cpp_setup_initBCMTFSystematicVariationDict() { G__add_setup_func("BCMTFSystematicVariationDict",(G__incsetup)(&G__cpp_setupBCMTFSystematicVariationDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initBCMTFSystematicVariationDict() { G__remove_setup_func("BCMTFSystematicVariationDict"); }
};
G__cpp_setup_initBCMTFSystematicVariationDict G__cpp_setup_initializerBCMTFSystematicVariationDict;

