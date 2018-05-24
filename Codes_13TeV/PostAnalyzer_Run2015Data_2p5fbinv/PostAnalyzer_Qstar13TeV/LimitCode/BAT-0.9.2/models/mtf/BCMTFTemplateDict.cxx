//
// File generated by /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.10-cms7//bin/rootcint at Sun Feb 14 14:58:42 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME BCMTFTemplateDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "BCMTFTemplateDict.h"

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
   void BCMTFTemplate_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void BCMTFTemplate_Dictionary();
   static void delete_BCMTFTemplate(void *p);
   static void deleteArray_BCMTFTemplate(void *p);
   static void destruct_BCMTFTemplate(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BCMTFTemplate*)
   {
      ::BCMTFTemplate *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BCMTFTemplate),0);
      static ::ROOT::TGenericClassInfo 
         instance("BCMTFTemplate", "./BCMTFTemplate.h", 30,
                  typeid(::BCMTFTemplate), DefineBehavior(ptr, ptr),
                  0, &BCMTFTemplate_Dictionary, isa_proxy, 0,
                  sizeof(::BCMTFTemplate) );
      instance.SetDelete(&delete_BCMTFTemplate);
      instance.SetDeleteArray(&deleteArray_BCMTFTemplate);
      instance.SetDestructor(&destruct_BCMTFTemplate);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BCMTFTemplate*)
   {
      return GenerateInitInstanceLocal((::BCMTFTemplate*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BCMTFTemplate*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void BCMTFTemplate_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::BCMTFTemplate*)0x0)->GetClass();
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_BCMTFTemplate(void *p) {
      delete ((::BCMTFTemplate*)p);
   }
   static void deleteArray_BCMTFTemplate(void *p) {
      delete [] ((::BCMTFTemplate*)p);
   }
   static void destruct_BCMTFTemplate(void *p) {
      typedef ::BCMTFTemplate current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BCMTFTemplate

/********************************************************
* BCMTFTemplateDict.cxx
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

extern "C" void G__cpp_reset_tagtableBCMTFTemplateDict();

extern "C" void G__set_cpp_environmentBCMTFTemplateDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("BCMTFTemplate.h");
  G__cpp_reset_tagtableBCMTFTemplateDict();
}
#include <new>
extern "C" int G__cpp_dllrevBCMTFTemplateDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* BCMTFTemplate */
static int G__BCMTFTemplateDict_170_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCMTFTemplate* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 2
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new BCMTFTemplate((const char*) G__int(libp->para[0]), (const char*) G__int(libp->para[1]));
   } else {
     p = new((void*) gvp) BCMTFTemplate((const char*) G__int(libp->para[0]), (const char*) G__int(libp->para[1]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_BCMTFTemplate));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFTemplateDict_170_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         string* pobj;
         string xobj = ((BCMTFTemplate*) G__getstructoffset())->GetChannelName();
         pobj = new string(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFTemplateDict_170_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         string* pobj;
         string xobj = ((BCMTFTemplate*) G__getstructoffset())->GetProcessName();
         pobj = new string(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFTemplateDict_170_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((BCMTFTemplate*) G__getstructoffset())->GetEfficiency());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFTemplateDict_170_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFTemplate*) G__getstructoffset())->GetHistogram());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFTemplateDict_170_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCMTFTemplate*) G__getstructoffset())->GetFunctionContainer());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFTemplateDict_170_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((BCMTFTemplate*) G__getstructoffset())->GetNBins());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFTemplateDict_170_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFTemplate*) G__getstructoffset())->SetEfficiency((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFTemplateDict_170_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFTemplate*) G__getstructoffset())->SetHistogram((TH1D*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFTemplateDict_170_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFTemplate*) G__getstructoffset())->SetFunctionContainer((vector<TF1*>*) G__int(libp->para[0]), (int) G__int(libp->para[1]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__BCMTFTemplateDict_170_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   BCMTFTemplate* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new BCMTFTemplate(*(BCMTFTemplate*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_BCMTFTemplate));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef BCMTFTemplate G__TBCMTFTemplate;
static int G__BCMTFTemplateDict_170_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (BCMTFTemplate*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((BCMTFTemplate*) (soff+(sizeof(BCMTFTemplate)*i)))->~G__TBCMTFTemplate();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (BCMTFTemplate*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((BCMTFTemplate*) (soff))->~G__TBCMTFTemplate();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__BCMTFTemplateDict_170_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCMTFTemplate* dest = (BCMTFTemplate*) G__getstructoffset();
   *dest = *(BCMTFTemplate*) libp->para[0].ref;
   const BCMTFTemplate& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* BCMTFTemplate */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncBCMTFTemplateDict {
 public:
  G__Sizep2memfuncBCMTFTemplateDict(): p(&G__Sizep2memfuncBCMTFTemplateDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncBCMTFTemplateDict::*p)();
};

size_t G__get_sizep2memfuncBCMTFTemplateDict()
{
  G__Sizep2memfuncBCMTFTemplateDict a;
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
extern "C" void G__cpp_setup_inheritanceBCMTFTemplateDict() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableBCMTFTemplateDict() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TF1*>",117,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlETF1mUcOallocatorlETF1mUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETF1mUcOallocatorlETF1mUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlETF1mUcOallocatorlETF1mUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETF1mUcOallocatorlETF1mUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlETF1mUcOallocatorlETF1mUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* BCMTFTemplate */
static void G__setup_memvarBCMTFTemplate(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__BCMTFTemplateDictLN_BCMTFTemplate));
   { BCMTFTemplate *p; p=(BCMTFTemplate*)0x1000; if (p) { }
   G__memvar_setup((void*)0,100,0,0,-1,-1,-1,4,"fEfficiency=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_TH1D),-1,-1,4,"fHistogram=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlETF1mUcOallocatorlETF1mUgRsPgR),G__defined_typename("vector<TF1*>"),-1,4,"fFunctionContainer=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"fNBins=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_string),-1,-1,4,"fChannelName=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFTemplateDictLN_string),-1,-1,4,"fProcessName=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarBCMTFTemplateDict() {
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
static void G__setup_memfuncBCMTFTemplate(void) {
   /* BCMTFTemplate */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__BCMTFTemplateDictLN_BCMTFTemplate));
   G__memfunc_setup("BCMTFTemplate",1192,G__BCMTFTemplateDict_170_0_1, 105, G__get_linked_tagnum(&G__BCMTFTemplateDictLN_BCMTFTemplate), -1, 0, 2, 1, 1, 0, 
"C - - 10 - channelname C - - 10 - processname", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetChannelName",1370,G__BCMTFTemplateDict_170_0_2, 117, G__get_linked_tagnum(&G__BCMTFTemplateDictLN_string), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetProcessName",1408,G__BCMTFTemplateDict_170_0_3, 117, G__get_linked_tagnum(&G__BCMTFTemplateDictLN_string), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetEfficiency",1301,G__BCMTFTemplateDict_170_0_4, 100, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetHistogram",1230,G__BCMTFTemplateDict_170_0_5, 85, G__get_linked_tagnum(&G__BCMTFTemplateDictLN_TH1D), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetFunctionContainer",2057,G__BCMTFTemplateDict_170_0_6, 85, G__get_linked_tagnum(&G__BCMTFTemplateDictLN_vectorlETF1mUcOallocatorlETF1mUgRsPgR), G__defined_typename("vector<TF1*>"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetNBins",762,G__BCMTFTemplateDict_170_0_7, 105, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetEfficiency",1313,G__BCMTFTemplateDict_170_0_8, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - eff", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetHistogram",1242,G__BCMTFTemplateDict_170_0_9, 121, -1, -1, 0, 1, 1, 1, 0, "U 'TH1D' - 0 - hist", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetFunctionContainer",2069,G__BCMTFTemplateDict_170_0_10, 121, -1, -1, 0, 2, 1, 1, 0, 
"U 'vector<TF1*,allocator<TF1*> >' 'vector<TF1*>' 0 - funccont i - - 0 - nbins", (char*)NULL, (void*) NULL, 0);
   // automatic copy constructor
   G__memfunc_setup("BCMTFTemplate", 1192, G__BCMTFTemplateDict_170_0_11, (int) ('i'), G__get_linked_tagnum(&G__BCMTFTemplateDictLN_BCMTFTemplate), -1, 0, 1, 1, 1, 0, "u 'BCMTFTemplate' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~BCMTFTemplate", 1318, G__BCMTFTemplateDict_170_0_12, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__BCMTFTemplateDict_170_0_13, (int) ('u'), G__get_linked_tagnum(&G__BCMTFTemplateDictLN_BCMTFTemplate), -1, 1, 1, 1, 1, 0, "u 'BCMTFTemplate' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncBCMTFTemplateDict() {
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
extern "C" void G__cpp_setup_globalBCMTFTemplateDict() {
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

extern "C" void G__cpp_setup_funcBCMTFTemplateDict() {
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
G__linked_taginfo G__BCMTFTemplateDictLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__BCMTFTemplateDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCMTFTemplateDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCMTFTemplateDictLN_TH1D = { "TH1D" , 99 , -1 };
G__linked_taginfo G__BCMTFTemplateDictLN_BCMTFTemplate = { "BCMTFTemplate" , 99 , -1 };
G__linked_taginfo G__BCMTFTemplateDictLN_vectorlETF1mUcOallocatorlETF1mUgRsPgR = { "vector<TF1*,allocator<TF1*> >" , 99 , -1 };
G__linked_taginfo G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETF1mUcOallocatorlETF1mUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TF1*,allocator<TF1*> >::iterator>" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableBCMTFTemplateDict() {
  G__BCMTFTemplateDictLN_string.tagnum = -1 ;
  G__BCMTFTemplateDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCMTFTemplateDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCMTFTemplateDictLN_TH1D.tagnum = -1 ;
  G__BCMTFTemplateDictLN_BCMTFTemplate.tagnum = -1 ;
  G__BCMTFTemplateDictLN_vectorlETF1mUcOallocatorlETF1mUgRsPgR.tagnum = -1 ;
  G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETF1mUcOallocatorlETF1mUgRsPgRcLcLiteratorgR.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableBCMTFTemplateDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__BCMTFTemplateDictLN_string);
   G__get_linked_tagnum_fwd(&G__BCMTFTemplateDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCMTFTemplateDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCMTFTemplateDictLN_TH1D);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__BCMTFTemplateDictLN_BCMTFTemplate),sizeof(BCMTFTemplate),-1,33792,(char*)NULL,G__setup_memvarBCMTFTemplate,G__setup_memfuncBCMTFTemplate);
   G__get_linked_tagnum_fwd(&G__BCMTFTemplateDictLN_vectorlETF1mUcOallocatorlETF1mUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFTemplateDictLN_reverse_iteratorlEvectorlETF1mUcOallocatorlETF1mUgRsPgRcLcLiteratorgR);
}
extern "C" void G__cpp_setupBCMTFTemplateDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupBCMTFTemplateDict()");
  G__set_cpp_environmentBCMTFTemplateDict();
  G__cpp_setup_tagtableBCMTFTemplateDict();

  G__cpp_setup_inheritanceBCMTFTemplateDict();

  G__cpp_setup_typetableBCMTFTemplateDict();

  G__cpp_setup_memvarBCMTFTemplateDict();

  G__cpp_setup_memfuncBCMTFTemplateDict();
  G__cpp_setup_globalBCMTFTemplateDict();
  G__cpp_setup_funcBCMTFTemplateDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncBCMTFTemplateDict();
  return;
}
class G__cpp_setup_initBCMTFTemplateDict {
  public:
    G__cpp_setup_initBCMTFTemplateDict() { G__add_setup_func("BCMTFTemplateDict",(G__incsetup)(&G__cpp_setupBCMTFTemplateDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initBCMTFTemplateDict() { G__remove_setup_func("BCMTFTemplateDict"); }
};
G__cpp_setup_initBCMTFTemplateDict G__cpp_setup_initializerBCMTFTemplateDict;

