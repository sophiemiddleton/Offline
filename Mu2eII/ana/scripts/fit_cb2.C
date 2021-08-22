///////////////////////////////////////////////////////////////////////////////
// cb2: Crystall Ball function : gaussian with power-law tails on both sides
// make it a C++ class to simplify fitting multiple functions
//
// to fit a resolution histogram: 
// ------------------------------
// .L Mu2eII/scripts/fit_cb2.C
// cbx = new cb2("x");
// cbx->fit(file,....)
//
// or
//
// TH1F*hist = ...
// cbx->fit(...)
//
// to find the range corresponding to the resolution peak half-maximum,
// scan a range of momenta with 
// 
// cbx->GetFunc()->Eval(104.5) 
///////////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "TNamed.h"
#include "TF1.h"

//-----------------------------------------------------------------------------
class cb2: public TNamed { 
public:

  TH1F* fHist;    // not owned by the fitter
  TH1F* fHGaus;   // owned by the fitter

  TF1*  fFunc;
  TF1*  fGaus;

  cb2(const char* Name = "x");
  ~cb2();

  static double  gauss(double* X, double* P);
  static double  f_cb2(double* X, double* P);    // double Crystall Ball itself

  void           init_parameters    (TF1* F, 
				     double P0, double P1, double P2, double P3, 
				     double P4, double P5, double P6);

  void           fit_histogram      (TH1F* Hist, TF1* F, double XMin, double XMax);
  void           create_fit_function(TF1*& F, double X0, double XMin, double XMax);
//-----------------------------------------------------------------------------
// two functions doing the same, but with different interfaces
//-----------------------------------------------------------------------------
  void           fit(TH1* Hist, double X0, double XMin, double XMax, double Sigma = -1.);

  void           fit(const char* File, const char* Module, const char* Hist, 
		     double X0, double XMin, double XMax, double Sigma = -1.);

  TF1* GetFunc() { return fFunc; }
  TH1* GetHist() { return (TH1*) fHist; }

};


//-----------------------------------------------------------------------------
cb2::cb2(const char* Name): TNamed(Name,Name) {
  fHist  = nullptr;
  fHGaus = nullptr;
  fFunc  = nullptr;
}

//-----------------------------------------------------------------------------
cb2::~cb2() {
}

//-----------------------------------------------------------------------------
double cb2::gauss(double* X, double* P) {
  double x0, f(0), sig, xdx;

  x0         = P[1];
  sig        = P[2];
  xdx        = (X[0]-x0)/sig;

  f = P[0]*TMath::Exp(-xdx*xdx/2.);
  
  return f;
}

//-----------------------------------------------------------------------------
double cb2::f_cb2(double* X, double* P) {
  double x0, x1, f, sig, alpha1, abs_alpha1, alpha2, abs_alpha2, n1, n2, a, b, dx0, dx1, xdx;

  x0         = P[1];
  sig        = P[2];
  alpha1     = P[3];
  n1         = P[4];
  alpha2     = P[5];
  n2         = P[6];

  dx0        = X[0]-x0;
  xdx        = dx0/sig;

  abs_alpha1 = fabs(alpha1);
  abs_alpha2 = fabs(alpha2);

  if (xdx < -abs_alpha1) {
    a = TMath::Power((n1/abs_alpha1),n1)*TMath::Exp(-(alpha1*alpha1)/2);
    b = n1/abs_alpha1-abs_alpha1;
    f = P[0]*a*TMath::Power((b-xdx),-n1);
  }
  else if (xdx < abs_alpha2) { 
					// gauss
    f = P[0]*TMath::Exp(-xdx*xdx/2.);
  }
  else {
    a = TMath::Power((n2/abs_alpha2),n2)*TMath::Exp(-(alpha2*alpha2)/2);
    b = -n2/abs_alpha2+abs_alpha2;
    f = P[0]*a*TMath::Power((xdx-b),-n2);
  }
  
  return f;
}

//-----------------------------------------------------------------------------
void cb2::init_parameters(TF1* F, 
			  double P0, double P1, double P2, double P3, 
			  double P4, double P5, double P6) {
  F->SetParameter(0,P0);
  F->SetParameter(1,P1);
  F->SetParameter(2,P2);
  F->SetParameter(3,P3);
  F->SetParameter(4,P4);
  F->SetParameter(5,P5);
  F->SetParameter(6,P6);
}

//-----------------------------------------------------------------------------
void cb2::fit_histogram(TH1F* Hist, TF1* F, double XMin, double XMax) {
  Hist->GetXaxis()->SetRangeUser(XMin,XMax);
  Hist->Fit(F,"L","",XMin,XMax);
}

//-----------------------------------------------------------------------------
void cb2::create_fit_function(TF1*& F, double X0, double XMin, double XMax) {

  fFunc = new TF1(Form("%s_f_cb2",GetName()),cb2::f_cb2,XMin,XMax,7);

  F->SetParName(0,"anorm");
  F->SetParName(1,"mean");
  F->SetParName(2,"sigma");
  F->SetParName(3,"alpha1");
  F->SetParName(4,"N1");
  F->SetParName(5,"alpha2");
  F->SetParName(6,"N2");

  F->SetParLimits(2,0.,1.);
  //  F->SetParLimits(5,0.,5.);
  //  F->SetParLimits(6,0.,5.);

  F->SetLineWidth(1);
}

//-----------------------------------------------------------------------------
void cb2::fit(TH1* Hist, double X0, double XMin, double XMax, double Sigma = -1.) {

  TH1F* h      = (TH1F*) Hist;
  double anorm = h->GetEntries()/10;

  printf("anorm = %12.5e\n",anorm);

  create_fit_function(fFunc,X0,XMin,XMax);
  init_parameters    (fFunc,anorm,X0,0.180,1,.4,1.,4);

//   h->Draw();
//   fFunc->Draw("same");

  if (Sigma > 0) {			// effectively, fix Sigma
    fFunc->SetParameter(2,Sigma);
    fFunc->SetParLimits(2,Sigma*(1+1.e-5),Sigma*(1+1.e-5));
  }

  fit_histogram(h,fFunc,XMin,XMax);
					// gaussian component of the fit, to display
  
  fGaus = new TF1(Form("%s_f_gaus",GetName()),cb2::gauss,XMin,XMax,3);

  fGaus->SetParameter(0,fFunc->GetParameter(0));
  fGaus->SetParameter(1,fFunc->GetParameter(1));
  fGaus->SetParameter(2,fFunc->GetParameter(2));

  fHGaus = (TH1F*) h->Clone(Form("%s_h_gaus",GetName()));

  fHGaus->Reset();

  for (int i=1; i<h->GetNbinsX(); i++) {
    float x = h->GetBinCenter(i);
    float y = fGaus->Eval(x);
    fHGaus->SetBinContent(i,y);
    fHGaus->SetBinError  (i,0);
  }
  
  fHGaus->SetLineColor(kBlue+3);
  fHGaus->SetFillColor(kBlue+3);
  fHGaus->SetFillStyle(3003);
  
  fHGaus->Draw("same");
//-----------------------------------------------------------------------------
// calculate the background fraction: high-momentum side , above the gaussian
// use region from alpha2 to XMax
//-----------------------------------------------------------------------------
  double total = Hist->Integral(1,Hist->GetNbinsX());

  double x0     = fFunc->GetParameter("mean"  );
  double sigma  = fFunc->GetParameter("sigma" );
  double alpha2 = fFunc->GetParameter("alpha2");
  double x1     = x0+fabs(sigma*alpha2);

  printf("x0 = %10.4f, sigma = %10.4f alpha2 = %10.4f\n",x0,sigma, alpha2);
 
  double htail = (fFunc->Integral(x1,XMax)-fGaus->Integral(x1,XMax))/Hist->GetBinWidth(1);

  printf("total = %12.5f, HTail = %12.5f, htail/total : %12.5e\n",total,htail,htail/total);
}

//-----------------------------------------------------------------------------
void cb2::fit(const char* File, const char* Module, const char* Hist, 
	      double X0, double XMin, double XMax, double Sigma = -1.) {

  fHist = (TH1F*) gh1(File,Module,Hist)->Clone("h_fit_crystal_ball");

  this->fit(fHist,X0,XMin,XMax,Sigma);
}

