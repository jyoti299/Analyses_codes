//
// File generated by /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.10-cms7//bin/rootcint at Sun Feb 14 14:56:16 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME BCH2DDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "BCH2DDict.h"

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
   void BCH2D_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void BCH2D_Dictionary();
   static void *new_BCH2D(void *p = 0);
   static void *newArray_BCH2D(Long_t size, void *p);
   static void delete_BCH2D(void *p);
   static void deleteArray_BCH2D(void *p);
   static void destruct_BCH2D(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BCH2D*)
   {
      ::BCH2D *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BCH2D),0);
      static ::ROOT::TGenericClassInfo 
         instance("BCH2D", "./../BAT/BCH2D.h", 34,
                  typeid(::BCH2D), DefineBehavior(ptr, ptr),
                  0, &BCH2D_Dictionary, isa_proxy, 0,
                  sizeof(::BCH2D) );
      instance.SetNew(&new_BCH2D);
      instance.SetNewArray(&newArray_BCH2D);
      instance.SetDelete(&delete_BCH2D);
      instance.SetDeleteArray(&deleteArray_BCH2D);
      instance.SetDestructor(&destruct_BCH2D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BCH2D*)
   {
      return GenerateInitInstanceLocal((::BCH2D*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BCH2D*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void BCH2D_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::BCH2D*)0x0)->GetClass();
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_BCH2D(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::BCH2D : new ::BCH2D;
   }
   static void *newArray_BCH2D(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::BCH2D[nElements] : new ::BCH2D[nElements];
   }
   // Wrapper around operator delete
   static void delete_BCH2D(void *p) {
      delete ((::BCH2D*)p);
   }
   static void deleteArray_BCH2D(void *p) {
      delete [] ((::BCH2D*)p);
   }
   static void destruct_BCH2D(void *p) {
      typedef ::BCH2D current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BCH2D

/********************************************************
* BCH2DDict.cxx
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

extern "C" void G__cpp_reset_tagtableBCH2DDict();

extern "C" void G__set_cpp_environmentBCH2DDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("../BAT/BCH2D.h");
  G__cpp_reset_tagtableBCH2DDict();
}
#include <new>
extern "C" int G__cpp_dllrevBCH2DDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* BCH2D */
static int G__BCH2DDict_171_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCH2D* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new BCH2D[n];
     } else {
       p = new((void*) gvp) BCH2D[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new BCH2D;
     } else {
       p = new((void*) gvp) BCH2D;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCH2DDictLN_BCH2D));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCH2D* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new BCH2D((TH2D*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) BCH2D((TH2D*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCH2DDictLN_BCH2D));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCH2D*) G__getstructoffset())->GetHistogram());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCH2D*) G__getstructoffset())->SetHistogram((TH2D*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCH2D*) G__getstructoffset())->SetGlobalMode((double*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 4:
      ((BCH2D*) G__getstructoffset())->Print((const char*) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (int) G__int(libp->para[3]));
      G__setnull(result7);
      break;
   case 3:
      ((BCH2D*) G__getstructoffset())->Print((const char*) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2]));
      G__setnull(result7);
      break;
   case 2:
      ((BCH2D*) G__getstructoffset())->Print((const char*) G__int(libp->para[0]), (int) G__int(libp->para[1]));
      G__setnull(result7);
      break;
   case 1:
      ((BCH2D*) G__getstructoffset())->Print((const char*) G__int(libp->para[0]));
      G__setnull(result7);
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 2:
      ((BCH2D*) G__getstructoffset())->Draw((int) G__int(libp->para[0]), (bool) G__int(libp->para[1]));
      G__setnull(result7);
      break;
   case 1:
      ((BCH2D*) G__getstructoffset())->Draw((int) G__int(libp->para[0]));
      G__setnull(result7);
      break;
   case 0:
      ((BCH2D*) G__getstructoffset())->Draw();
      G__setnull(result7);
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCH2D*) G__getstructoffset())->CalculateIntegratedHistogram();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((BCH2D*) G__getstructoffset())->GetLevel((double) G__double(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         vector<int>* pobj;
         vector<int> xobj = ((BCH2D*) G__getstructoffset())->GetNIntervalsY((TH2D*) G__int(libp->para[0]), *(int*) G__Intref(&libp->para[1]));
         pobj = new vector<int>(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCH2D*) G__getstructoffset())->GetLowestBandGraph((TH2D*) G__int(libp->para[0]), *((vector<int>*) G__int(libp->para[1]))));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCH2D*) G__getstructoffset())->GetLowestBandGraph((TH2D*) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         vector<double>* pobj;
         vector<double> xobj = ((BCH2D*) G__getstructoffset())->GetLevelBoundary((double) G__double(libp->para[0]));
         pobj = new vector<double>(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         vector<double>* pobj;
         vector<double> xobj = ((BCH2D*) G__getstructoffset())->GetLevelBoundary((TH2D*) G__int(libp->para[0]), (double) G__double(libp->para[1]));
         pobj = new vector<double>(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCH2D*) G__getstructoffset())->GetBandGraph((double) G__double(libp->para[0]), (double) G__double(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCH2D*) G__getstructoffset())->GetBandGraph((TH2D*) G__int(libp->para[0]), (double) G__double(libp->para[1])
, (double) G__double(libp->para[2])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCH2DDict_171_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((BCH2D*) G__getstructoffset())->GetBandGraphs((TH2D*) G__int(libp->para[0]), *(int*) G__Intref(&libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__BCH2DDict_171_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   BCH2D* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new BCH2D(*(BCH2D*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCH2DDictLN_BCH2D));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef BCH2D G__TBCH2D;
static int G__BCH2DDict_171_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (BCH2D*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((BCH2D*) (soff+(sizeof(BCH2D)*i)))->~G__TBCH2D();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (BCH2D*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((BCH2D*) (soff))->~G__TBCH2D();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__BCH2DDict_171_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCH2D* dest = (BCH2D*) G__getstructoffset();
   *dest = *(BCH2D*) libp->para[0].ref;
   const BCH2D& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* BCH2D */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncBCH2DDict {
 public:
  G__Sizep2memfuncBCH2DDict(): p(&G__Sizep2memfuncBCH2DDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncBCH2DDict::*p)();
};

size_t G__get_sizep2memfuncBCH2DDict()
{
  G__Sizep2memfuncBCH2DDict a;
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
extern "C" void G__cpp_setup_inheritanceBCH2DDict() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableBCH2DDict() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__BCH2DDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCH2DDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCH2DDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCH2DDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCH2DDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__BCH2DDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCH2DDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCH2DDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCH2DDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCH2DDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<int>",117,G__get_linked_tagnum(&G__BCH2DDictLN_vectorlEintcOallocatorlEintgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCH2DDictLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCH2DDictLN_vectorlEintcOallocatorlEintgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCH2DDictLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCH2DDictLN_vectorlEintcOallocatorlEintgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* BCH2D */
static void G__setup_memvarBCH2D(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__BCH2DDictLN_BCH2D));
   { BCH2D *p; p=(BCH2D*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BCH2DDictLN_TH2D),-1,-1,4,"fHistogram=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BCH2DDictLN_TH1D),-1,-1,4,"fIntegratedHistogram=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,-1,-1,4,"fMode[2]=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"fModeFlag=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarBCH2DDict() {
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
static void G__setup_memfuncBCH2D(void) {
   /* BCH2D */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__BCH2DDictLN_BCH2D));
   G__memfunc_setup("BCH2D",323,G__BCH2DDict_171_0_1, 105, G__get_linked_tagnum(&G__BCH2DDictLN_BCH2D), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("BCH2D",323,G__BCH2DDict_171_0_2, 105, G__get_linked_tagnum(&G__BCH2DDictLN_BCH2D), -1, 0, 1, 1, 1, 0, "U 'TH2D' - 0 - h", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetHistogram",1230,G__BCH2DDict_171_0_3, 85, G__get_linked_tagnum(&G__BCH2DDictLN_TH2D), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetHistogram",1242,G__BCH2DDict_171_0_4, 121, -1, -1, 0, 1, 1, 1, 0, "U 'TH2D' - 0 - hist", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetGlobalMode",1282,G__BCH2DDict_171_0_5, 121, -1, -1, 0, 1, 1, 1, 0, "D - - 0 - mode", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Print",525,G__BCH2DDict_171_0_6, 121, -1, -1, 0, 4, 1, 1, 0, 
"C - - 10 - filename i - - 0 '0' options "
"i - - 0 '0' ww i - - 0 '0' wh", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Draw",398,G__BCH2DDict_171_0_7, 121, -1, -1, 0, 2, 1, 1, 0, 
"i - - 0 '0' options g - - 0 'true' drawmode", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("CalculateIntegratedHistogram",2883,G__BCH2DDict_171_0_8, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetLevel",792,G__BCH2DDict_171_0_9, 100, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - p", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetNIntervalsY",1407,G__BCH2DDict_171_0_10, 117, G__get_linked_tagnum(&G__BCH2DDictLN_vectorlEintcOallocatorlEintgRsPgR), G__defined_typename("vector<int>"), 0, 2, 1, 1, 0, 
"U 'TH2D' - 0 - h i - - 1 - nfoundmax", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetLowestBandGraph",1797,G__BCH2DDict_171_0_11, 85, G__get_linked_tagnum(&G__BCH2DDictLN_TGraph), -1, 0, 2, 1, 1, 0, 
"U 'TH2D' - 0 - h u 'vector<int,allocator<int> >' 'vector<int>' 0 - nint", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetLowestBandGraph",1797,G__BCH2DDict_171_0_12, 85, G__get_linked_tagnum(&G__BCH2DDictLN_TGraph), -1, 0, 1, 1, 1, 0, "U 'TH2D' - 0 - h", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetLevelBoundary",1628,G__BCH2DDict_171_0_13, 117, G__get_linked_tagnum(&G__BCH2DDictLN_vectorlEdoublecOallocatorlEdoublegRsPgR), G__defined_typename("vector<double>"), 0, 1, 1, 1, 0, "d - - 0 - level", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetLevelBoundary",1628,G__BCH2DDict_171_0_14, 117, G__get_linked_tagnum(&G__BCH2DDictLN_vectorlEdoublecOallocatorlEdoublegRsPgR), G__defined_typename("vector<double>"), 0, 2, 1, 1, 0, 
"U 'TH2D' - 0 - h d - - 0 - level", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetBandGraph",1159,G__BCH2DDict_171_0_15, 85, G__get_linked_tagnum(&G__BCH2DDictLN_TGraph), -1, 0, 2, 1, 1, 0, 
"d - - 0 - level1 d - - 0 - level2", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetBandGraph",1159,G__BCH2DDict_171_0_16, 85, G__get_linked_tagnum(&G__BCH2DDictLN_TGraph), -1, 0, 3, 1, 1, 0, 
"U 'TH2D' - 0 - h d - - 0 - level1 "
"d - - 0 - level2", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetBandGraphs",1274,G__BCH2DDict_171_0_17, 85, G__get_linked_tagnum(&G__BCH2DDictLN_TGraph), -1, 2, 2, 1, 1, 0, 
"U 'TH2D' - 0 - h i - - 1 - n", (char*)NULL, (void*) NULL, 0);
   // automatic copy constructor
   G__memfunc_setup("BCH2D", 323, G__BCH2DDict_171_0_18, (int) ('i'), G__get_linked_tagnum(&G__BCH2DDictLN_BCH2D), -1, 0, 1, 1, 1, 0, "u 'BCH2D' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~BCH2D", 449, G__BCH2DDict_171_0_19, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__BCH2DDict_171_0_20, (int) ('u'), G__get_linked_tagnum(&G__BCH2DDictLN_BCH2D), -1, 1, 1, 1, 1, 0, "u 'BCH2D' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncBCH2DDict() {
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
extern "C" void G__cpp_setup_globalBCH2DDict() {
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

extern "C" void G__cpp_setup_funcBCH2DDict() {
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
G__linked_taginfo G__BCH2DDictLN_vectorlEdoublecOallocatorlEdoublegRsPgR = { "vector<double,allocator<double> >" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_TH1D = { "TH1D" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_TH2D = { "TH2D" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_TGraph = { "TGraph" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_BCH2D = { "BCH2D" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_vectorlEintcOallocatorlEintgRsPgR = { "vector<int,allocator<int> >" , 99 , -1 };
G__linked_taginfo G__BCH2DDictLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<int,allocator<int> >::iterator>" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableBCH2DDict() {
  G__BCH2DDictLN_vectorlEdoublecOallocatorlEdoublegRsPgR.tagnum = -1 ;
  G__BCH2DDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__BCH2DDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCH2DDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__BCH2DDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCH2DDictLN_TH1D.tagnum = -1 ;
  G__BCH2DDictLN_TH2D.tagnum = -1 ;
  G__BCH2DDictLN_TGraph.tagnum = -1 ;
  G__BCH2DDictLN_BCH2D.tagnum = -1 ;
  G__BCH2DDictLN_vectorlEintcOallocatorlEintgRsPgR.tagnum = -1 ;
  G__BCH2DDictLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableBCH2DDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_vectorlEdoublecOallocatorlEdoublegRsPgR);
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_TH1D);
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_TH2D);
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_TGraph);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__BCH2DDictLN_BCH2D),sizeof(BCH2D),-1,34048,(char*)NULL,G__setup_memvarBCH2D,G__setup_memfuncBCH2D);
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_vectorlEintcOallocatorlEintgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCH2DDictLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR);
}
extern "C" void G__cpp_setupBCH2DDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupBCH2DDict()");
  G__set_cpp_environmentBCH2DDict();
  G__cpp_setup_tagtableBCH2DDict();

  G__cpp_setup_inheritanceBCH2DDict();

  G__cpp_setup_typetableBCH2DDict();

  G__cpp_setup_memvarBCH2DDict();

  G__cpp_setup_memfuncBCH2DDict();
  G__cpp_setup_globalBCH2DDict();
  G__cpp_setup_funcBCH2DDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncBCH2DDict();
  return;
}
class G__cpp_setup_initBCH2DDict {
  public:
    G__cpp_setup_initBCH2DDict() { G__add_setup_func("BCH2DDict",(G__incsetup)(&G__cpp_setupBCH2DDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initBCH2DDict() { G__remove_setup_func("BCH2DDict"); }
};
G__cpp_setup_initBCH2DDict G__cpp_setup_initializerBCH2DDict;

