//Draw stacks using the DataPlotter object
#include "DataPlotter.C"

DataPlotter* dataplotter_ = 0;
TString path_ = "/mu2e/data/users/mmackenz/Mu2eII/Mu2eII/histograms/"; //path to histogram directory

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
//  Signal configurations  //
/////////////////////////////

double br_conv_ = 1.e-15; //branching ratio of conversion signal
bool doPositron_ = true; //switch between mu- --> e+ and e-



//construct the name of the RMC file from given parameters
TString get_rmc_file(bool isExternal, int spectrum, int version, double kmax, int batch_mode) {
  TString name = "rmc";
  if(!isExternal) name += "i1s51b";
  else if(spectrum == 0 && version == 1) name += "e0s61b";
  else if(spectrum == 3) name += "e4s51b";
  else name += "e1s51b";
  name += batch_mode;
  TString k = Form("%.1f", kmax);
  k.ReplaceAll(".","d");
  name += Form(".Mu2eII_rmc_ana.s%iv0%i.k_%s.hist", spectrum, version, k.Data());
  return name;
}

//get file name for the dataset
TString get_file_name(TString dataset, int batch_mode) {
  TString name = "b";
  name += batch_mode;
  if     (dataset == "cele") name = "cele0s61" + name + ".Mu2eII_conv_ana.hist";
  else if(dataset == "cpos") name = "cpos0s51" + name + ".Mu2eII_conv_ana.hist";
  else if(dataset == "cosm") name = "cosm0s91b0.Mu2eII_cosmic_ana.hist"; //use no pileup file
  else if(dataset == "cry3") name = "cry31s91b0.Mu2eII_cosmic_ana.hist"; //use no pileup file
  else if(dataset == "rpce") name = "rpce0s51" + name + ".Mu2eII_rpc_ana.hist";
  else if(dataset == "rpci") name = "rpci0s51" + name + ".Mu2eII_rpc_ana.hist";
  else if(dataset == "pbar") name = "pbar0s11b0.Mu2eII_pbar_ana.hist"; //use no pileup file
  else if(dataset == "dio" ) name = "fele2s51" + name + ".Mu2eII_track_ana.hist";
  else if(dataset == "dio_mdc") name = "e10s721z.Mu2eII_track_ana.hist";
  else if(dataset == "rmce") name = get_rmc_file(true , rmcSpectrum_, rmcExternalVersion_, kmax_, batch_mode);
  else if(dataset == "rmci") name = get_rmc_file(false, rmcSpectrum_, rmcInternalVersion_, kmax_, batch_mode);
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
double get_normalization(TString dataset, int batch_mode, bool scale_set = false) {
  double scale = 1.;

  ////////////////////////////////
  // Batch efficiency reduction //
  ////////////////////////////////
  double batch_1_eff = 226227./231040.; //From trk set 2000 cele0s61b1 / cele0s51b0 
  double batch_2_eff = 215629./231040.; //From trk set 2000 cele0s61b2 / cele0s51b0 

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
  //  Cosmic veto scale   //
  //////////////////////////
  const double cosmHiPosScale =   5./ 9021.; //Scale no veto --> veto for positrons for the hi sample, from trk set 4002 --> 4000
  const double cosmHiEleScale =  12./34429.; //Scale no veto --> veto for electrons for the hi sample, from trk set 2002 --> 2000
  const double cosmLoPosScale = 151./ 1819.; //Scale no veto --> veto for positrons for the lo sample, from trk set 4002 --> 4000
  const double cosmLoEleScale = 212./ 6082.; //Scale no veto --> veto for electrons for the lo sample, from trk set 2002 --> 2000

  //////////////////////////
  // Normalization params //
  //////////////////////////
  if(dataset == "cpos") {
    const long ngen_su_conv = 1000000;
    scale = 1./ngen_su_conv;
    scale *= stop_fraction*br_conv_*(capture_fraction);
    return scale;
  }else if(dataset == "cele") {
    const long ngen_su_conv = 1000000;
    scale = 1./ngen_su_conv;
    scale *= stop_fraction*br_conv_*(capture_fraction);
    return scale;
  } else if(dataset == "rmce0") {
    const long   ngen_su_rmce  = (batch_mode == 0) ? 1832902322 : 1833779978;
    const double emin_su_rmce  = 85.;
    const double emax_su_rmce  = 93.;
    scale =  ((emax_su_rmce-emin_su_rmce)/kmax_) / ngen_su_rmce;
    scale  *= stop_fraction*capture_fraction*br_rmc;
    return scale;
  } else if(dataset == "rmci1") {
    const long   ngen_su_rmci  = 10000000;
    const double emin_su_rmci  = 57.;
    const double emax_su_rmci  = 102.;
    scale  = ((emax_su_rmci-emin_su_rmci)/kmax_) / ngen_su_rmci;
    scale *= stop_fraction*capture_fraction*br_rmc*br_int_rmc;
    return scale;
  } else if(dataset == "cosm") {
    const double cosm_livetime_su = 3.83e8;
    scale = 1./cosm_livetime_su*(t0End - t0Begin)/t0Event;
    if(batch_mode == 1) scale *= batch_1_eff;
    if(batch_mode == 2) scale *= batch_2_eff;
    if(scale_set) scale *= (doPositron_) ? cosmLoPosScale : cosmLoEleScale;
    return scale;
  } else if(dataset == "cry3") {
    const double cosm_livetime_su = 1.28e7;
    scale = 1./cosm_livetime_su*(t0End - t0Begin)/t0Event;
    if(batch_mode == 1) scale *= batch_1_eff;
    if(batch_mode == 2) scale *= batch_2_eff;
    if(scale_set) scale *= (doPositron_) ? cosmHiPosScale : cosmHiEleScale;
    return scale;
  } else if(dataset == "rpce") {
    const long ngen_rpce = 100000000;
    scale = pion_stopping_fraction*pion_survival*br_rpc/ngen_rpce;
    if(scale_set) scale *= 1.171e-8/0.0006745; //reduction from no t0 cut --> t0 > 700 ns
    return scale;
  } else if(dataset == "rpci") {
    const long ngen_rpci = 1000000;
    scale = pion_stopping_fraction*pion_survival*br_rpc*br_int_rpc/ngen_rpci;
    if(scale_set) scale *= 3.31e-8/0.001412; //reduction from no t0 cut --> t0 > 700 ns 3000 --> 3001
    return scale;
  } else if(dataset == "dio") {
    const long ngen_dio = 1000000;
    const double emax_dio = 110.;
    const double emin_dio = 85.;
    scale = stop_fraction*(1.-capture_fraction)/ngen_dio*(emax_dio - emin_dio);
    return scale;
  } else if(dataset == "dio_mdc") {
    const long ngen_dio = 2000000;
    const double emax_dio = 110.;
    const double emin_dio = 75.;
    scale = stop_fraction*(1.-capture_fraction)/ngen_dio*(emax_dio - emin_dio);
    return scale;
  }
  cout << "ERROR! Unknown dataset " << dataset.Data() << endl;
  return 1.;
}

//load data cards into the DataPlotter
int init_dataplotter() {
  std::vector<DataCard_t> cards;
  //constructor:            isOneBatch        fname                       fpath                      label                scale                      isSignal isBeam  color  setOffset
  cards.push_back(DataCard_t(true ,path_+get_file_name("rpce", 1), "Ana/Mu2eII_RPCAna/Hist"   , "RPC(External)", get_normalization("rpce" , 1,  true), false, true , kGreen-2 ,-1));
  cards.push_back(DataCard_t(false,path_+get_file_name("rpce", 1), "Ana/Mu2eII_RPCAna/Hist"   , "RPC(External)", get_normalization("rpce" , 1,  true), false, true , kGreen-2 , 3));
  cards.push_back(DataCard_t(true ,path_+get_file_name("rpci", 1), "Ana/Mu2eII_RPCAna/Hist"   , "RPC(Internal)", get_normalization("rpci" , 1,  true), false, true , kGreen-4 ,-1));
  cards.push_back(DataCard_t(false,path_+get_file_name("rpci", 1), "Ana/Mu2eII_RPCAna/Hist"   , "RPC(Internal)", get_normalization("rpci" , 1,  true), false, true , kGreen-4 , 3));
  cards.push_back(DataCard_t(true ,path_+get_file_name("cosm", 0), "Ana/Mu2eII_CosmicAna/Hist", "Cosmic(lo)"   , get_normalization("cosm" , 1,  true), false, false, kYellow+1, 2));
  cards.push_back(DataCard_t(false,path_+get_file_name("cosm", 0), "Ana/Mu2eII_CosmicAna/Hist", "Cosmic(lo)"   , get_normalization("cosm" , 2,  true), false, false, kYellow+1, 2));
  cards.push_back(DataCard_t(true ,path_+get_file_name("cry3", 0), "Ana/Mu2eII_CosmicAna/Hist", "Cosmic(hi)"   , get_normalization("cry3" , 1,  true), false, false, kOrange+1, 2));
  cards.push_back(DataCard_t(false,path_+get_file_name("cry3", 0), "Ana/Mu2eII_CosmicAna/Hist", "Cosmic(hi)"   , get_normalization("cry3" , 2,  true), false, false, kOrange+1, 2));
  cards.push_back(DataCard_t(true ,path_+get_file_name("rmce", 1), "Ana/Mu2eII_RMCAna/Hist"   , "RMC(External)", get_normalization("rmce0", 1, false), false, true , kRed+2      ));
  cards.push_back(DataCard_t(false,path_+get_file_name("rmce", 1), "Ana/Mu2eII_RMCAna/Hist"   , "RMC(External)", get_normalization("rmce0", 1, false), false, true , kRed+2   , 4));
  cards.push_back(DataCard_t(true ,path_+get_file_name("rmci", 0), "Ana/Mu2eII_RMCAna/Hist"   , "RMC(Internal)", get_normalization("rmci1", 0, false), false, true , kRed        ));
  cards.push_back(DataCard_t(false,path_+get_file_name("rmci", 0), "Ana/Mu2eII_RMCAna/Hist"   , "RMC(Internal)", get_normalization("rmci1", 0, false), false, true , kRed        ));
  cards.push_back(DataCard_t(true ,path_+get_file_name("dio" , 1), "Ana/Mu2eII_TrackAna/Hist" , "DIO"          , get_normalization("dio"  , 1, false), false, true , kViolet-2   ));
  cards.push_back(DataCard_t(false,path_+get_file_name("dio" , 1), "Ana/Mu2eII_TrackAna/Hist" , "DIO"          , get_normalization("dio"  , 1, false), false, true , kViolet-2, 4));
  if(doPositron_) {
    cards.push_back(DataCard_t(true ,path_+get_file_name("cpos", 1), "Ana/Mu2eII_ConvAna/Hist" , "#mu^{-}#rightarrow e^{+}", get_normalization("cpos" , 1), true , true , kBlue));
    cards.push_back(DataCard_t(false,path_+get_file_name("cpos", 1), "Ana/Mu2eII_ConvAna/Hist" , "#mu^{-}#rightarrow e^{+}", get_normalization("cpos" , 1), true , true , kBlue));
  } else {
    cards.push_back(DataCard_t(true ,path_+get_file_name("cele", 1), "Ana/Mu2eII_ConvAna/Hist" , "#mu^{-}#rightarrow e^{-}", get_normalization("cele" , 1), true , true , kBlue));
    cards.push_back(DataCard_t(false,path_+get_file_name("cele", 2), "Ana/Mu2eII_ConvAna/Hist" , "#mu^{-}#rightarrow e^{-}", get_normalization("cele" , 2), true , true , kBlue));
  }
  if(dataplotter_) delete dataplotter_;
  dataplotter_ = new DataPlotter();
  dataplotter_->lumi_[0]     = 2.8631579e19; //one batch
  dataplotter_->lumi_[1]     = 9.0285712e18; //two batch
  const double cycle_1batch  = 1.33; //seconds
  const double pot_1batch    = 4.e12; //pot/cycle
  const double cycle_2batch  = 1.4; //seconds
  const double pot_2batch    = 8.e12; //pot/cycle
  dataplotter_->livetime_[0] = dataplotter_->lumi_[0]/pot_1batch*cycle_1batch;
  dataplotter_->livetime_[1] = dataplotter_->lumi_[1]/pot_2batch*cycle_2batch;
  
  return dataplotter_->AddFiles(cards);
}

//Print standard figures for both conversion channels
int print_standard_plots() {
  int status(0);
  TCanvas* c = 0;

  //mu- --> e- figures
  doPositron_ = false;
  status += init_dataplotter();
  dataplotter_->logy_ = 1;
  c = dataplotter_->PrintStack(PlottingCard_t("p_2", "trk", 1001, 1, 98., 106., 1.e-5, 1.e3, "P (MeV/c)", ""));
  if(c) delete c;
  else ++status;
  c = dataplotter_->PrintStack(PlottingCard_t("p_2", "trk", 2001, 1, 98., 106., 1.e-5, 1.e3, "P (MeV/c)", ""));
  if(c) delete c;
  else ++status;

  //mu- --> e+ figures
  doPositron_ = true;
  status += init_dataplotter();
  dataplotter_->logy_ = 1;
  c = dataplotter_->PrintStack(PlottingCard_t("p_2", "trk", 3001, 1, 85., 93., 1.e-5, 1.e3, "P (MeV/c)", ""));
  if(c) delete c;
  else ++status;
  c = dataplotter_->PrintStack(PlottingCard_t("p_2", "trk", 4001, 1, 85., 93., 1.e-5, 1.e3, "P (MeV/c)", ""));
  if(c) delete c;
  else ++status;

  return status;
}

//use values from CD3 for comparisons
int init_cd3() {
  int status = init_dataplotter();
  if(status) return status;
  dataplotter_->lumi_[0] = 0.; //no one batch
  dataplotter_->lumi_[1] = 3.6e20; //all two batch
  const double cycle_2batch  = 1.4; //seconds
  const double pot_2batch    = 8.e12; //pot/cycle
  dataplotter_->livetime_[0] = 0.;
  dataplotter_->livetime_[1] = dataplotter_->lumi_[1]/pot_2batch*cycle_2batch;
  return 0;
}
