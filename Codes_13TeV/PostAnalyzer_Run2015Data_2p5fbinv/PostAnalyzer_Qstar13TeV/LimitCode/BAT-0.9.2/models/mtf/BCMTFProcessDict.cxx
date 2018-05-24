//
// File generated by /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.10-cms7//bin/rootcint at Sun Feb 14 14:58:39 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME BCMTFProcessDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "BCMTFProcessDict.h"

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
   void BCMTFProcess_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void BCMTFProcess_Dictionary();
   static void delete_BCMTFProcess(void *p);
   static void deleteArray_BCMTFProcess(void *p);
   static void destruct_BCMTFProcess(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BCMTFProcess*)
   {
      ::BCMTFProcess *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BCMTFProcess),0);
      static ::ROOT::TGenericClassInfo 
         instance("BCMTFProcess", "./BCMTFProcess.h", 27,
                  typeid(::BCMTFProcess), DefineBehavior(ptr, ptr),
                  0, &BCMTFProcess_Dictionary, isa_proxy, 0,
                  sizeof(::BCMTFProcess) );
      instance.SetDelete(&delete_BCMTFProcess);
      instance.SetDeleteArray(&deleteArray_BCMTFProcess);
      instance.SetDestructor(&destruct_BCMTFProcess);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BCMTFProcess*)
   {
      return GenerateInitInstanceLocal((::BCMTFProcess*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BCMTFProcess*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void BCMTFProcess_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::BCMTFProcess*)0x0)->GetClass();
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_BCMTFProcess(void *p) {
      delete ((::BCMTFProcess*)p);
   }
   static void deleteArray_BCMTFProcess(void *p) {
      delete [] ((::BCMTFProcess*)p);
   }
   static void destruct_BCMTFProcess(void *p) {
      typedef ::BCMTFProcess current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BCMTFProcess

/********************************************************
* BCMTFProcessDict.cxx
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

extern "C" void G__cpp_reset_tagtableBCMTFProcessDict();

extern "C" void G__set_cpp_environmentBCMTFProcessDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("BCMTFProcess.h");
  G__cpp_reset_tagtableBCMTFProcessDict();
}
#include <new>
extern "C" int G__cpp_dllrevBCMTFProcessDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* BCMTFProcess */
static int G__BCMTFProcessDict_168_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCMTFProcess* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new BCMTFProcess((const char*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) BCMTFProcess((const char*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCMTFProcessDictLN_BCMTFProcess));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFProcessDict_168_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         string* pobj;
         string xobj = ((BCMTFProcess*) G__getstructoffset())->GetName();
         pobj = new string(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BCMTFProcessDict_168_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((BCMTFProcess*) G__getstructoffset())->SetName((const char*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__BCMTFProcessDict_168_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   BCMTFProcess* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new BCMTFProcess(*(BCMTFProcess*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BCMTFProcessDictLN_BCMTFProcess));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef BCMTFProcess G__TBCMTFProcess;
static int G__BCMTFProcessDict_168_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (BCMTFProcess*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((BCMTFProcess*) (soff+(sizeof(BCMTFProcess)*i)))->~G__TBCMTFProcess();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (BCMTFProcess*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((BCMTFProcess*) (soff))->~G__TBCMTFProcess();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__BCMTFProcessDict_168_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   BCMTFProcess* dest = (BCMTFProcess*) G__getstructoffset();
   *dest = *(BCMTFProcess*) libp->para[0].ref;
   const BCMTFProcess& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* BCMTFProcess */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncBCMTFProcessDict {
 public:
  G__Sizep2memfuncBCMTFProcessDict(): p(&G__Sizep2memfuncBCMTFProcessDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncBCMTFProcessDict::*p)();
};

size_t G__get_sizep2memfuncBCMTFProcessDict()
{
  G__Sizep2memfuncBCMTFProcessDict a;
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
extern "C" void G__cpp_setup_inheritanceBCMTFProcessDict() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableBCMTFProcessDict() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__BCMTFProcessDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFProcessDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFProcessDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFProcessDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFProcessDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__BCMTFProcessDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BCMTFProcessDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFProcessDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BCMTFProcessDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BCMTFProcessDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* BCMTFProcess */
static void G__setup_memvarBCMTFProcess(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__BCMTFProcessDictLN_BCMTFProcess));
   { BCMTFProcess *p; p=(BCMTFProcess*)0x1000; if (p) { }
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__BCMTFProcessDictLN_string),-1,-1,4,"fName=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarBCMTFProcessDict() {
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
static void G__setup_memfuncBCMTFProcess(void) {
   /* BCMTFProcess */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__BCMTFProcessDictLN_BCMTFProcess));
   G__memfunc_setup("BCMTFProcess",1099,G__BCMTFProcessDict_168_0_1, 105, G__get_linked_tagnum(&G__BCMTFProcessDictLN_BCMTFProcess), -1, 0, 1, 1, 1, 0, "C - - 10 - name", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetName",673,G__BCMTFProcessDict_168_0_2, 117, G__get_linked_tagnum(&G__BCMTFProcessDictLN_string), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetName",685,G__BCMTFProcessDict_168_0_3, 121, -1, -1, 0, 1, 1, 1, 0, "C - - 10 - name", (char*)NULL, (void*) NULL, 0);
   // automatic copy constructor
   G__memfunc_setup("BCMTFProcess", 1099, G__BCMTFProcessDict_168_0_4, (int) ('i'), G__get_linked_tagnum(&G__BCMTFProcessDictLN_BCMTFProcess), -1, 0, 1, 1, 1, 0, "u 'BCMTFProcess' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~BCMTFProcess", 1225, G__BCMTFProcessDict_168_0_5, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__BCMTFProcessDict_168_0_6, (int) ('u'), G__get_linked_tagnum(&G__BCMTFProcessDictLN_BCMTFProcess), -1, 1, 1, 1, 1, 0, "u 'BCMTFProcess' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncBCMTFProcessDict() {
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
extern "C" void G__cpp_setup_globalBCMTFProcessDict() {
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

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcBCMTFProcessDict() {
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
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__BCMTFProcessDictLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__BCMTFProcessDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__BCMTFProcessDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCMTFProcessDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__BCMTFProcessDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BCMTFProcessDictLN_BCMTFProcess = { "BCMTFProcess" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableBCMTFProcessDict() {
  G__BCMTFProcessDictLN_string.tagnum = -1 ;
  G__BCMTFProcessDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__BCMTFProcessDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCMTFProcessDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__BCMTFProcessDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BCMTFProcessDictLN_BCMTFProcess.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableBCMTFProcessDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__BCMTFProcessDictLN_string);
   G__get_linked_tagnum_fwd(&G__BCMTFProcessDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFProcessDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BCMTFProcessDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BCMTFProcessDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__BCMTFProcessDictLN_BCMTFProcess),sizeof(BCMTFProcess),-1,33792,(char*)NULL,G__setup_memvarBCMTFProcess,G__setup_memfuncBCMTFProcess);
}
extern "C" void G__cpp_setupBCMTFProcessDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupBCMTFProcessDict()");
  G__set_cpp_environmentBCMTFProcessDict();
  G__cpp_setup_tagtableBCMTFProcessDict();

  G__cpp_setup_inheritanceBCMTFProcessDict();

  G__cpp_setup_typetableBCMTFProcessDict();

  G__cpp_setup_memvarBCMTFProcessDict();

  G__cpp_setup_memfuncBCMTFProcessDict();
  G__cpp_setup_globalBCMTFProcessDict();
  G__cpp_setup_funcBCMTFProcessDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncBCMTFProcessDict();
  return;
}
class G__cpp_setup_initBCMTFProcessDict {
  public:
    G__cpp_setup_initBCMTFProcessDict() { G__add_setup_func("BCMTFProcessDict",(G__incsetup)(&G__cpp_setupBCMTFProcessDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initBCMTFProcessDict() { G__remove_setup_func("BCMTFProcessDict"); }
};
G__cpp_setup_initBCMTFProcessDict G__cpp_setup_initializerBCMTFProcessDict;

