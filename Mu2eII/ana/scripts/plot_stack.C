//Script to make a stack of track histograms, normalized for each physics process

#include "DatasetInfo.C"

//////////////////////////////////////////
//                                      //
// Change data paths to personal files! //
//                                      //
//////////////////////////////////////////

//Configure for 0, 1, or 2 batch mode
int batch_mode_ = 0; //default to no pileup

//Define the path to the histogram directory
const char* data_path          = "/mu2e/data/users/mmackenz/Mu2eII/Mu2eII/histograms/";

/***********
RMC datasets information
k_xydz: k_{max} = xy.z MeV
sxvyz: Spectrum x version yz
Spectrums and versions:
Spectrum 0: Closure approximation
 version 0: Standard closure approximation
 version 1: Signal region only sample, using resampled photon conversions
Spectrum 1: Closure approximation + flat spectrum up to the kinematic endpoint
 version 1: Br(flat > closure end point) / Br(RMC < closure end point) = 1.e-4
Spectrum 2: Closure approximation + exp spectrum up to the kinematic endpoint
 version 1: Br(exp > closure end point) / Br(RMC < closure end point) = 1.e-4, slope = -3
Spectrum 2: Closure approximation + mono-enenergetic photon transition
 version 1: Br(transition) / Br(RMC < closure end point) = 1.e-4, E = 101.85 MeV
Spectrum 3: Mono-energetic photon sample at E = 101.85 MeV
 version 0: Default
***********/
//for normalization and file selection of RMC datasets
double kmax_ = 90.1; 
int rmcSpectrum_ = 0; 
int rmcExternalVersion_ = 1;
int rmcInternalVersion_ = 0;

/////////////////////////////
// Plotting configurations //
/////////////////////////////

bool unitNorm_ = false; //normalize histograms to 1
bool doHists_ = false; //plot individual histograms instead of a stack
bool print_ = false; //whether or not to print the figures to disk

/////////////////////////////
//  Signal configurations  //
/////////////////////////////

double br_conv_ = 1.e-15; //branching ratio of conversion signal
bool doPositron_ = true; //switch between mu- --> e+ and e-

//construct the name of the RMC file from given parameters
TString get_rmc_file(bool isExternal, int spectrum, int version, double kmax) {
  TString name = "rmc";
  if(!isExternal) name += "i1s51b";
  else if(spectrum == 0 && version == 1) name += "e0s61b";
  else if(spectrum == 3) name += "e4s51b";
  else name += "e1s51b";
  // name += batch_mode_;
  name += 0; //FIXME: Include batch mode for RMC when available
  TString k = Form("%.1f", kmax);
  k.ReplaceAll(".","d");
  name += Form(".Mu2eII_rmc_ana.s%iv0%i.k_%s.hist", spectrum, version, k.Data());
  return name;
}

//get file name for the dataset
TString get_file_name(TString dataset) {
  TString name = "b";
  name += batch_mode_;
  if     (dataset == "cele") name = "cele0s51" + name + ".Mu2eII_conv_ana.hist";
  else if(dataset == "cpos") name = "cpos0s51" + name + ".Mu2eII_conv_ana.hist";
  else if(dataset == "cosm") name = "cosm0s91b0.Mu2eII_cosmic_ana.hist"; //use no pileup file
  else if(dataset == "cry3") name = "cry31s91b0.Mu2eII_cosmic_ana.hist"; //use no pileup file
  // else if(dataset == "rpce") name = "rpce0s51" + name + ".Mu2eII_rpc_ana.hist";
  // else if(dataset == "rpci") name = "rpci0s51" + name + ".Mu2eII_rpc_ana.hist";
  else if(dataset == "rpce") name = "rpce0s51b0.Mu2eII_rpc_ana.hist"; //FIXME: swith to pileup
  else if(dataset == "rpci") name = "rpci0s51b0.Mu2eII_rpc_ana.hist"; //FIXME: switch to pileup
  else if(dataset == "pbar") name = "pbar0s11b0.Mu2eII_pbar_ana.hist"; //use no pileup file
  else if(dataset == "dio" ) name = "fele0s51" + name + ".Mu2eII_track_ana.hist";
  else if(dataset == "dio_mdc") name = "e10s721z.Mu2eII_track_ana.hist";
  else if(dataset == "rmce") name = get_rmc_file(true , rmcSpectrum_, rmcExternalVersion_, kmax_);
  else if(dataset == "rmci") name = get_rmc_file(true , rmcSpectrum_, rmcExternalVersion_, kmax_);
  else {
    cout << "Warning in get_file_name! Dataset " << dataset.Data() << " unknown!\n";
    name = "";
  }
  return name;
}

//integral of the RMC closure approximation
double cl_approx_integral(double xmin, double xmax) {
  const double x = xmin;
  const double y = xmax;
  double integral = 1./3.*x*x*(-20.*pow(x,4)+72.*pow(x,3) -105.*x*x + 80.*x - 30.);
  integral       -= 1./3.*y*y*(-20.*pow(y,4)+72.*pow(y,3) -105.*y*y + 80.*y - 30.);
  return integral;
}

//get normalization values for a given dataset
double get_normalization(TString dataset) {
  double scale_conv, scale_rmce, scale_rmci, scale_cosm, scale_rpce, scale_rpci, scale_dio;
  //scale to experiment
  // experiment parameters from:
  // https://mu2e-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=32752&filename=20200416_CEsensitivity.pdf&version=2

  ///////////////////////////
  // Experiment parameters //
  ///////////////////////////
  const double npot_1batch      = 2.8631579e19;
  const double npot_2batch      = 9.0285712e18;
  const double npot             = npot_1batch+npot_2batch; //6.85e19; //3.6e20;
  const double cycle_1batch     = 1.33; //seconds
  const double pot_1batch       = 4.e12; //pot/cycle
  const double cycle_2batch     = 1.4; //seconds
  const double pot_2batch       = 8.e12; //pot/cycle
  const double time_1batch      = npot_1batch/pot_1batch*cycle_1batch;
  const double time_2batch      = npot_2batch/pot_2batch*cycle_2batch;
  const double livetime         = time_1batch+time_2batch; //3.5e6; //6.e7;

  //////////////////////////
  //      Muon stops      //
  //////////////////////////
  const double stop_fraction    = 0.0015;
  const double capture_fraction = 0.61;

  //////////////////////////
  //         RMC          //
  //////////////////////////
  const double br_rmc           = 1.43e-5/cl_approx_integral(57./kmax_, 1.); //br known for > 57 MeV from TRIUMF
  const double br_int_rmc       = 0.0069; //internal RMC only

  //////////////////////////
  //         RPC          //
  //////////////////////////
  const double br_rpc           = 0.0215; //branching ratio of RPC
  const double br_int_rpc       = 0.0069; //internal RPC only
  const double pion_stopping_fraction = 0.00211; //N(pion stops) / N(POT)
  const double pion_survival    = 0.02423; //fraction of pions with T > 450 ns (simulation cut off)

  //////////////////////////
  //   Event window       //
  //////////////////////////
  const double t0Begin          = 700.; //ns
  const double t0End            = 1695.; //ns
  const double t0Event          = 1695.; //time length of an event window

  //////////////////////////
  // Normalization params //
  //////////////////////////
  if(dataset == "cpos") {
    const long ngen_su_conv = 1000000;
    scale_conv = 1./ngen_su_conv;
    scale_conv *= npot*stop_fraction*br_conv_*(capture_fraction);
    return scale_conv;
  }else if(dataset == "cele") {
    const long ngen_su_conv = 1000000;
    scale_conv = 1./ngen_su_conv;
    scale_conv *= npot*stop_fraction*br_conv_*(capture_fraction);
    return scale_conv;
  } else if(dataset == "rmce0") {
    const long   ngen_su_rmce  = 1832902322;
    const double emin_su_rmce  = 85.;
    const double emax_su_rmce  = 93.;
    scale_rmce =  ((emax_su_rmce-emin_su_rmce)/kmax_) / ngen_su_rmce;
    scale_rmce  *= npot*stop_fraction*capture_fraction*br_rmc;
    return scale_rmce;
  } else if(dataset == "rmci1") {
    const long   ngen_su_rmci  = 10000000;
    const double emin_su_rmci  = 57.;
    const double emax_su_rmci  = 102.;
    scale_rmci =  ((emax_su_rmci-emin_su_rmci)/kmax_) / ngen_su_rmci;
    scale_rmci  *= npot*stop_fraction*capture_fraction*br_rmc*br_int_rmc;
    return scale_rmci;
  } else if(dataset == "cosm") {
    const double cosm_livetime_su = 3.83e8;
    scale_cosm = livetime/cosm_livetime_su*(t0End - t0Begin)/t0Event;
    return scale_cosm;
  } else if(dataset == "cry3") {
    const double cosm_livetime_su = 1.28e7;
    scale_cosm = livetime/cosm_livetime_su*(t0End - t0Begin)/t0Event;
    return scale_cosm;
  } else if(dataset == "rpce0") {
    const long ngen_rpce = 100000000;
    scale_rpce = npot*pion_stopping_fraction*pion_survival*br_rpc/ngen_rpce;
    return scale_rpce;
  } else if(dataset == "rpci0") {
    const long ngen_rpci = 1000000;
    scale_rpci = npot*pion_stopping_fraction*pion_survival*br_rpc*br_int_rpc/ngen_rpci;
    return scale_rpci;
  } else if(dataset == "dio") {
    const long ngen_dio = 1000000;
    const double emax_dio = 110.;
    const double emin_dio = 85.;
    scale_dio = npot*stop_fraction*(1.-capture_fraction)/ngen_dio*(emax_dio - emin_dio);
    return scale_dio;
  } else if(dataset == "dio_mdc") {
    const long ngen_dio = 2000000;
    const double emax_dio = 110.;
    const double emin_dio = 75.;
    scale_dio = npot*stop_fraction*(1.-capture_fraction)/ngen_dio*(emax_dio - emin_dio);
    return scale_dio;
  }
  cout << "ERROR! Unknown dataset " << dataset << endl;
  return 1.;
}

//Open and return the file for a given dataset
TFile* get_file(TString dataset) {
  TString path = data_path;
  TString file_name = get_file_name(dataset);
  if(file_name != "")
    return TFile::Open((path + file_name).Data(),"READ");  
  return NULL;
}

//get the stack of background samples properly normalized
THStack* get_background_stack(int set, TString hist, int rebin = 1, TString type = "trk") {
  //open background files
  TFile* frmce = get_file("rmce");
  TFile* frmci = get_file("rmci");
  TFile* fcosm = get_file("cosm");
  TFile* fcosm_hi = get_file("cry3");
  TFile* frpce = get_file("rpce");
  TFile* frpci = get_file("rpci");
  bool useMDC = true; // FIXME: switch to SU2020 DIO when available
  TFile* fdio  = get_file((useMDC) ? "dio_mdc" : "dio"); 
  
  //check they exist
  if(!frmce || !frmci ||!fcosm || !fcosm_hi || !frpce || !frpci || (!doPositron_&&!fdio)) {cout << "Not all histogram files were found!\n"; return NULL;}

  //format path based on module used for histogramming
  TString rmce_path = Form("Ana/Mu2eII_RMCAna/Hist/%s_%i/%s" , type.Data(), set, hist.Data());
  TString rmci_path = Form("Ana/Mu2eII_RMCAna/Hist/%s_%i/%s" , type.Data(), set, hist.Data());
  TString cosm_path = Form("Ana/Mu2eII_CosmicAna/Hist/%s_%i/%s", type.Data(), set, hist.Data());
  TString cry3_path = Form("Ana/Mu2eII_CosmicAna/Hist/%s_%i/%s", type.Data(), set, hist.Data());
  TString rpce_path = Form("Ana/Mu2eII_RPCAna/Hist/%s_%i/%s" , type.Data(), set, hist.Data());
  TString rpci_path = Form("Ana/Mu2eII_RPCAna/Hist/%s_%i/%s" , type.Data(), set, hist.Data());
  TString dio_path  = Form("Ana/Mu2eII_TrackAna/Hist/%s_%i/%s" , type.Data(), set, hist.Data());

  //get histograms from the files, check they exist
  TH1F* hrmce = (TH1F*) frmce->Get(rmce_path.Data());
  if(!hrmce) {cout << "No external RMC histogram!\n"; return NULL;}
  hrmce->SetName("hrmce");  
  TH1F* hrmci = (TH1F*) frmci->Get(rmci_path.Data());
  if(!hrmci) {cout << "No internal RMC histogram!\n"; return NULL;}
  hrmci->SetName("hrmci");
  TH1F* hcosm = (TH1F*) fcosm->Get(cosm_path.Data());
  if(!hcosm) {cout << "No Cosmic lo histogram!\n"; return NULL;}
  hcosm->SetName("hcosm");
  TH1F* hcosm_hi = (TH1F*) fcosm_hi->Get(cry3_path.Data());
  if(!hcosm_hi) {cout << "No Cosmic hi histogram!\n"; return NULL;}
  hcosm_hi->SetName("hcosm_hi");
  TH1F* hrpce = (TH1F*) frpce->Get(rpce_path.Data());
  if(!hrpce) {cout << "No external RPC histogram!\n"; return NULL;}
  hrpce->SetName("hrpce");  
  TH1F* hrpci = (TH1F*) frpci->Get(rpci_path.Data());
  if(!hrpci) {cout << "No internal RPC histogram!\n"; return NULL;}
  hrpci->SetName("hrpci");
  TH1F* hdio = (doPositron_) ? 0 : (TH1F*) fdio->Get(dio_path.Data());
  if(!doPositron_) {
    if(!hdio) {cout << "No internal RPC histogram!\n"; return NULL;}
    hdio->SetName("hdio");
  }
  
  //Rebin the histograms if needed
  if(rebin > 1) {
    hrmce->Rebin(rebin);
    hrmci->Rebin(rebin);
    hcosm->Rebin(rebin);
    hcosm_hi->Rebin(rebin);
    hrpce->Rebin(rebin);
    hrpci->Rebin(rebin);
    if(hdio) hdio->Rebin(rebin);
  }
  
  //Scale the histograms, either to unit of physics based normalizations
  hrmce->Scale(((unitNorm_) ? 1./hrmce->Integral() : get_normalization("rmce0")));
  hrmci->Scale(((unitNorm_) ? 1./hrmci->Integral() : get_normalization("rmci1")));
  hcosm->Scale(((unitNorm_) ? 1./hcosm->Integral() : get_normalization("cosm")));
  hcosm_hi->Scale(((unitNorm_) ? 1./hcosm_hi->Integral() : get_normalization("cry3")));
  hrpce->Scale(((unitNorm_) ? 1./hrpce->Integral() : get_normalization("rpce0")));
  hrpci->Scale(((unitNorm_) ? 1./hrpci->Integral() : get_normalization("rpci0")));
  if(hdio)
    hdio->Scale(((unitNorm_) ? 1./hdio->Integral() : get_normalization((useMDC) ? "dio_mdc" : "dio")));
  
  //Set colors
  hrmce->SetLineWidth(2);
  hrmce->SetLineColor(kRed+3);
  if(!doHists_) hrmce->SetFillColor(kRed+2);
  hrmci->SetLineWidth(2);
  hrmci->SetLineColor(kRed+1);
  if(!doHists_) hrmci->SetFillColor(kRed);
  hcosm->SetLineWidth(2);
  hcosm->SetLineColor(kYellow+2);
  if(!doHists_) hcosm->SetFillColor(kYellow+1);
  hcosm_hi->SetLineWidth(2);
  hcosm_hi->SetLineColor(kOrange+2);
  if(!doHists_) hcosm_hi->SetFillColor(kOrange+1);
  hrpce->SetLineWidth(2);
  hrpce->SetLineColor(kGreen-1);
  if(!doHists_) hrpce->SetFillColor(kGreen-2);
  hrpci->SetLineWidth(2);
  hrpci->SetLineColor(kGreen-5);
  if(!doHists_) hrpci->SetFillColor(kGreen-6);
  if(hdio) {
    hdio->SetLineWidth(2);
    hdio->SetLineColor(kViolet-1);
    if(!doHists_) hdio->SetFillColor(kViolet-2);
  }

  //Set the histogram titles
  hcosm->SetTitle("Cosmics (2025lo)");
  hcosm_hi->SetTitle("Cosmics (2025hi)");
  hrpce->SetTitle("External RPC");
  hrpci->SetTitle("Internal RPC");
  hrmce->SetTitle(Form("External RMC k_{max} = %.1f MeV", kmax_));
  hrmci->SetTitle(Form("Internal RMC k_{max} = %.1f MeV", kmax_));
  if(hdio)
    hdio->SetTitle("DIO");

  //Create the stack
  THStack* hstack = new THStack("hstack", "Mu2e background stack");
  hstack->Add(hcosm);
  hstack->Add(hcosm_hi);
  hstack->Add(hrpce);
  hstack->Add(hrpci);
  hstack->Add(hrmce);
  hstack->Add(hrmci);
  if(hdio) hstack->Add(hdio);
  
  return hstack;
}

//Get the signal histogram
TH1F* get_signal(int set, TString hist, int rebin, TString type) {
  TFile* fconv = (doPositron_) ? get_file("cpos") : get_file("cele");
  if(!fconv) return NULL;
  TString path = Form("Ana/Mu2eII_ConvAna/Hist/%s_%i/%s", type.Data(), set, hist.Data());
  TH1F* hconv = (TH1F*) fconv->Get(path.Data());
  if(!hconv) {cout << "No conversion signal histogram!\n"; return NULL;}
  if(rebin > 1) hconv->Rebin(rebin);
  hconv->Scale((unitNorm_) ? 1./hconv->Integral() : ((doPositron_) ? get_normalization("cpos") : get_normalization("cele")));
  hconv->SetLineWidth(2);
  hconv->SetFillStyle(3003);
  hconv->SetFillColor(kBlue);
  hconv->SetLineColor(kBlue);
  return hconv;
}

TCanvas* plot_stack(TString hist = "p2", int set = 4000, double xmin = 85., double xmax = 93., 
		    double ymin = 1.e-3, double ymax_log = 1.e3, double ymax_lin = 1.,
		    int rebin = 1, TString type = "trk") { 
  THStack* hstack = get_background_stack(set, hist, rebin, type);
  if(!hstack) {cout << "Error! No background stack returned!\n"; return NULL;}
  TH1F* hsignal = get_signal(set, hist, rebin, type);  
  if(!hsignal) {cout << "Error! No signal histogram returned!\n"; return NULL;}
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas(Form("c_%s_%i",hist.Data(), set), Form("c_%s_%i",hist.Data(), set), 1000, 800);  

  TH1F* hax;
  TAxis *xax, *yax;
  if(doHists_) {
    bool first = true;
    for(TObject* h : *(hstack->GetHists())) {
      h->Draw(((first) ? "hist E1" : "hist E1 sames"));
      if(first) {
	first = false;
	hax = (TH1F*) h;
	xax = hax->GetXaxis();
	yax = hax->GetYaxis();
      }
    }
  } else {
    hstack->Draw("hist E1 noclear");
    xax = hstack->GetXaxis();
    yax = hstack->GetYaxis();
  }
  hsignal->Draw("hist E1 same");

  //Configure the axes
  if(xmin < xmax)                 { xax->SetRangeUser(xmin, xmax);}
  if(ymin < ymax_lin && !doHists_){ hstack->SetMinimum(ymin); hstack->SetMaximum(ymax_lin);}
  else if(ymin < ymax_lin)        { yax->SetRangeUser(ymin, ymax_lin);}
  if(doHists_) hax->SetTitle("Mu2e backgrounds");

  //Setup the legend
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9, "", "LP brNDC");
  leg->AddEntry(hsignal, Form("#mu^{-} #rightarrow e^{%s} Br = %.1e", (doPositron_) ? "+" : "-", br_conv_));
  for(TObject* h : *hstack->GetHists())
    leg->AddEntry(h);
  leg->Draw("same");
  c->Modified(); c->Update(); 
  if(print_) {
    gSystem->Exec("[ ! -d figures ] && mkdir figures"); //make the directory if needed
    TString filename = Form("figures/%s_%s_%s%s_%i", type.Data(), hist.Data(), (doHists_) ? "hist" : "stack",
			    (unitNorm_) ? "_norm" : "", set);
    c->SaveAs((filename+".png").Data());
    c->SaveAs((filename+".pdf").Data());
    if(ymin < ymax_log && !doHists_){ hstack->SetMinimum(ymin); hstack->SetMaximum(ymax_log);}
    else if(ymin < ymax_log)        { yax->SetRangeUser(ymin, ymax_log);}
    c->SetLogy();
    c->SaveAs((filename+"_log.png").Data());
    c->SaveAs((filename+"_log.pdf").Data());
  } else if(ymin < ymax_log) {
    if(!doHists_) {hstack->SetMinimum(ymin); hstack->SetMaximum(ymax_log);}
    else          { yax->SetRangeUser(ymin, ymax_log);}
    c->SetLogy(); 
  }
  return c;
}

int print_standard_figures() {
  int status = 0;
  TCanvas* c = 0;
  print_ = true;
  bool prevBatch = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);
  //conversion positron
  c = plot_stack("p2", 3000, 85., 93., 1.e-4, 1.e3, 5);
  status += !c;
  c = plot_stack("p2", 3001, 85., 93., 1.e-4, 1.e3, 5);
  status += !c;
  c = plot_stack("p2", 4000, 85., 93., 1.e-4, 1.e3, 5);
  status += !c;
  c = plot_stack("p2", 4001, 85., 93., 1.e-4, 1.e3, 5);
  status += !c;
  //conversion electron
  doPositron_ = false;
  c = plot_stack("p2", 1000, 100., 107., 1.e-4, 1.e3, 5);
  status += !c;
  c = plot_stack("p2", 1001, 100., 107., 1.e-4, 1.e3, 5);
  status += !c;
  c = plot_stack("p2", 2000, 100., 107., 1.e-4, 1.e3, 5);
  status += !c;
  c = plot_stack("p2", 2001, 100., 107., 1.e-4, 1.e3, 5);
  status += !c;
  gROOT->SetBatch(prevBatch);
  return status;
}
