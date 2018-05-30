//---------------------------------------------------------------------------
#ifndef HiggsAnalysis_CombinedLimit_RooF5_5ParamPdf_h
#define HiggsAnalysis_CombinedLimit_RooF5_5ParamPdf_h
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
class RooF5_5ParamPdf : public RooAbsPdf
{
public:
   RooF5_5ParamPdf() {} ;
   RooF5_5ParamPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,
                    RooAbsReal& _p2, RooAbsReal& _p3,
		    RooAbsReal& _p4, RooAbsReal& _sqrts);
   RooF5_5ParamPdf(const RooF5_5ParamPdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooF5_5ParamPdf(*this,newname); }
   inline virtual ~RooF5_5ParamPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   RooRealProxy p3;        // p3
   RooRealProxy p4;        // p4
   RooRealProxy sqrts;        // sqrts
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooF5_5ParamPdf,1) // RazorDijetBinPdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class F5_5ParamFunction: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double pdf = pow(x/p[0],p[1])*exp(1+(p[2]*x/p[0])+p[3]*pow(x/p[0],2)+p[4]*pow(x/p[0],3)); 
     return pdf;
   }
   
   double DoEval(double x) const
   {
     double pdf = pow(x/pars[0],pars[1])*exp(1+(pars[2]*x/pars[0])+pars[3]*pow(x/pars[0],2)+pars[4]*pow(x/pars[0],3));
     return pdf;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new F5_5ParamFunction();
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
      return 5;
   }
};
