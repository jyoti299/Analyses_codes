//---------------------------------------------------------------------------
#ifndef HiggsAnalysis_CombinedLimit_RooF1_3ParamPdf_h
#define HiggsAnalysis_CombinedLimit_RooF1_3ParamPdf_h
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include <TH1.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

//---------------------------------------------------------------------------
class RooF1_3ParamPdf : public RooAbsPdf
{
public:
   RooF1_3ParamPdf() {} ;
   RooF1_3ParamPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,
		  RooAbsReal& _p2,
		  RooAbsReal& _sqrts);
   RooF1_3ParamPdf(const RooF1_3ParamPdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooF1_3ParamPdf(*this,newname); }
   inline virtual ~RooF1_3ParamPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   RooRealProxy sqrts;        // sqrts
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooF1_3ParamPdf,1) // RazorDijetBinPdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class F1_3ParamFunction: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double pdf = pow(1+x/p[0],p[1])/pow(x/p[0],p[2]);
     return pdf;
   }
   
   double DoEval(double x) const
   {
     double pdf = pow(1+x/pars[0],pars[1])/pow(x/pars[0],pars[2]);
     return pdf;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new F1_3ParamFunction();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 3;
   }
};
