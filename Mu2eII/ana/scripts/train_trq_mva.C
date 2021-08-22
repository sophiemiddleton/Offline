///////////////////////////////////////////////////////////////////////////////
// initially cloned from Dave's code, diverged. 
//
// AlgorithmCode : 0:PAR, 1:DAR
// BkgMode       : 0,1,2,3,4  (+100, if use Z)
//
// location of the input training ntuples: /mu2e/data/projects/Mu2eII/trq_tmva_training
//
// negative tracks:
// ----------------
// fele2s51b1.tmva_training_0000.root
// fele2s51b1.tmva_training_0001.root
//
// positive tracks:
// ----------------
// fpos2s51b1.tmva_training_0000.root
// fpos2s51b1.tmva_training_1000.root
//
// As all tracks come from a single source ntuple, use just one input ntuple, not two
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

//#include "TMVAGui.C"
#include "TMVA/Config.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

enum e_bgr_weight{
  e_unit        = 0, 
  e_linear      = 1, 
  e_exponential = 2,
  e_polynomial  = 3,
  e_exp4        = 4
};

void   init_used_methods (std::map<std::string,int>& Use);
void   book_methods      (TMVA::Factory* factory, TMVA::DataLoader *dataloader, std::map<std::string,int>& Use);
int    trq_tmva          (const char* FileName, const char* TreeName, int AlgorithmCode, float DpPin);
//-----------------------------------------------------------------------------
// main function, called by the user
// Algorithm = 0 : PAR
// Algorithm = 1 : DAR
// DpMin     = 0.7 ...
// input tree doesn't know about the background weighting mode
//-----------------------------------------------------------------------------
int train_trq_mva(const char* DsID = "fele2s51b1", int Algorithm = 0, int Charge=-1, int BkgWeightMode = 0, float DpMin = 0.7) {

  TString dir = "/mu2e/data/projects/Mu2eII/trq_tmva_training";
  TString fn  = dir+Form("/%s.tmva_training_%04i.root",DsID,1000*Algorithm);
  int code    = 1000*Algorithm+BkgWeightMode;
  if (Charge == 1) code = code+100;

  trq_tmva(fn,"tmva_training_tree",code,DpMin);

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/


//-----------------------------------------------------------------------------
// Algorithm     = 0: PAR
// Algorithm     = 1: PAR
// AlgorithmCode = 1000*Algorithm+BkgWeightMode+100(if positive charge)
// Treename      = "tmva_training_tree" (this is the default)
//-----------------------------------------------------------------------------
int trq_tmva(const char* FileName, const char* TreeName, int AlgorithmCode, float DpMin) {

  std::map<std::string,int> Use;   // map defines MVA methods to be trained + tested

  TMVA::Tools::Instance();

  init_used_methods(Use);          // sets everything to zero
  Use["MLP"] = 1;                  // Recommended ANN
  Use["BDT"] = 0; 

  printf("\n==> Start TMVAClassification: input ntuple: %s input tree: %s\n",FileName,TreeName);
//-----------------------------------------------------------------------------
// input file
//-----------------------------------------------------------------------------
  TFile* f      = TFile::Open(FileName);
  TTree* MyTree = (TTree*) f->Get(TreeName);
//-----------------------------------------------------------------------------
// Create a ROOT output file where TMVA will store ntuples, histograms, etc.
//-----------------------------------------------------------------------------
  int code = AlgorithmCode+100*DpMin;
  TString tname       = Form("tmva_training_output_%04i",code);
  TString outfilename = tname+".root";
  TFile*  outputFile  = TFile::Open( outfilename, "RECREATE" );
//-----------------------------------------------------------------------------
// Create factory object. Later you can choose the methods performance of which you'd like to investigate. 
// The factory is the only TMVA object you have to interact with
//
// The first argument is the base of the name of all the weightfiles in the directory weight/
//
// The second argument is the output file for the training results
// All TMVA output can be suppressed by removing the "!" (not) in
// front of the "Silent" argument in the option string
//-----------------------------------------------------------------------------
  TMVA::DataLoader *dataloader = new TMVA::DataLoader(tname);
  
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

  // If you wish to modify default settings
  // (please check "src/Config.h" to see all available global options)
  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

  std::string dirname(tname);
  dirname += "Weights";
  (TMVA::gConfig().GetIONames()).fWeightFileDir = dirname.c_str();

//-----------------------------------------------------------------------------
// Define input variables to be used for training
// note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
// [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
//-----------------------------------------------------------------------------
  dataloader->AddVariable("nactive"     ,  "NActive"          ,"Count"   ,'F');
  dataloader->AddVariable("nafract"     ,  "FActive"          ,"Fraction",'F');
  dataloader->AddVariable("log10(fcons)",  "LogFitCon"        ,"Probability",'F');
  dataloader->AddVariable("momerr"      ,  "FitMomErr"        , "MeV/c"   , 'F');
  dataloader->AddVariable("t0err"       ,  "T0Err"            , "nsec"    , 'F');
  // dataloader->AddVariable("d0"          ,  "D0"               , "mm"      , 'F');
  // dataloader->AddVariable("rmax"        ,  "MaxRadius"        , "mm"      , 'F');
  dataloader->AddVariable("nda_o_na"    ,  "DoubleHitFraction", "Fraction", 'F');
  // dataloader->AddVariable("nda/nd"    ,  "DoubleHitFraction", "Fraction", 'F');
  dataloader->AddVariable("nza_o_na"    ,  "NullHitFraction"  , "Fraction", 'F');
  dataloader->AddVariable("nma_o_nm"    ,  "MatFraction"      , "Fraction", 'F');

  // if (usez) {
  //   dataloader->AddVariable("z1"    ,"Z_first_hit"      , "mm"      , 'F');
  // }

  // You can add so-called "Spectator variables", which are not used in the MVA training,
  // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
  // input variables, the response values of all trained MVAs, and the spectator variables

  Double_t sw = 1.0;
  Double_t bw = 1.0;

  dataloader->AddSignalTree    (MyTree, sw);
  dataloader->AddBackgroundTree(MyTree, bw);

  // If no numbers of events are given, half of the events in the tree are used for training, and the other half for testing:
  //
  //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
  //
  // To also specify the number of testing events, use:
  //    factory->PrepareTrainingAndTestTree( mycut,
  //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
//-----------------------------------------------------------------------------
// different weighting schemes of the background event - they don't change much
//-----------------------------------------------------------------------------
  TString sig_weight("");
  TString bkg_weight("");

  int bkg_weight_mode = AlgorithmCode % 100;

  if      (bkg_weight_mode == e_unit       ) bkg_weight = "1.";
  else if (bkg_weight_mode == e_linear     ) bkg_weight = "max(1.0,5.0*(min(p-pmc,3.0)))";
  else if (bkg_weight_mode == e_exponential) bkg_weight = "max(1.0,exp(2.0*min(p-pmc,3.0)))";
  else if (bkg_weight_mode == e_polynomial ) bkg_weight = "max(1.0,pow(2.0*min(p-pmc,3.0),5.0))";
  else if (bkg_weight_mode == e_exp4       ) bkg_weight = "max(1.0,pow(10.0*min(p-pmc,3.0),5.0))";
  else { 
    printf (" >>> ERROR: unknown background weight mode: %i. EXIT\n",bkg_weight_mode);
    return -1;
  }

  if (sig_weight != "") {
    printf("Signal     weight expression = %s\n",sig_weight.Data());
    dataloader->SetSignalWeightExpression    (sig_weight);
  }
  if (bkg_weight != "") {
    printf("Background weight expression = %s\n",bkg_weight.Data());
    dataloader->SetBackgroundWeightExpression(bkg_weight);
  }
//-----------------------------------------------------------------------------
// signal is defined as the momentum resolution core,
// tail is defined as the high-side tail above 0.70 or 0.60
//-----------------------------------------------------------------------------
//  TCut signal_cuts("(pmc>100)&&(tdip>0.5)&&(tdip<1.0)&&(t0err<5)&&(p-pmc<0.25)&&(p-pmc)>-0.25");
//  TCut bkg_cuts("(pmc>100)&&(tdip>0.5)&&(tdip<1.0)&&(t0err<50)&&(p-pmc)>0.7");

  TCut signal_cuts("(pmc>86)&&(p>85)&&(tdip>0.5)&&(tdip<1.0)&&(t0err<5)&&(fabs(p-pmc)<0.25)");

  TCut bkg_cuts  = Form("(pmc>86)&&(p>85)&&(tdip>0.5)&&(tdip<1.0)&&(t0err<5)&&(p-pmc)>%4.2f",DpMin);

  TString training_opt = "nTrain_Signal=50000:nTrain_Background=20000";
  training_opt        += ":nTest_Signal=50000:nTest_Background=20000";
  //  training_opt        += ":SplitMode=Random:!V:SplitSeed=592309";
  training_opt        += ":SplitMode=Random:!V";

  dataloader->PrepareTrainingAndTestTree(signal_cuts,bkg_cuts,training_opt.Data());

  book_methods (factory,dataloader,Use);
//-----------------------------------------------------------------------------
// For an example of the category classifier usage, see: TMVAClassificationCategory
// optimize settings (configuration) of the MVAs using the set of training events
//-----------------------------------------------------------------------------
  // factory->OptimizeAllMethods("SigEffAt001","Scan");
  // factory->OptimizeAllMethods("ROCIntegral","GA");
// --------------------------------------------------------------------------------------------------
// tell the factory to train, test, and evaluate the MVAs
//-----------------------------------------------------------------------------
//  factory->PrintHelpMessage();
//  factory->EvaluateAllVariables();

  factory->SetVerbose();
  factory->TrainAllMethods();    // Train MVAs using the set of training events
  factory->TestAllMethods();     // Evaluate all MVAs using the set of test events
  factory->EvaluateAllMethods(); // Evaluate and compare performance of all configured MVAs

  outputFile->Close();           // Save output

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  delete factory;
//-----------------------------------------------------------------------------
// Launch the GUI for the root macros
//-----------------------------------------------------------------------------
  if (!gROOT->IsBatch()) {
    TMVA::TMVAGui(outfilename); 
  }

  return 0;
}

//-----------------------------------------------------------------------------
void init_used_methods(std::map<std::string,int>& Use) {

  // --- Cut optimisation
  Use["Cuts"]            = 0;
  Use["CutsD"]           = 0;
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  // 
  // --- 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 0;
  Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  //
  // --- Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 0;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 0;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 0; // k-nearest neighbour method
  //
  // --- Linear Discriminant Analysis
  Use["LD"]              = 0; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
  //
  // --- Function Discriminant analysis
  Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  //
  // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN" ]        = 0; // ROOT's own ANN
  //
  // --- Support Vector Machine 
  Use["SVM"]             = 0;
  // 
  // --- Boosted Decision Trees
  Use["BDT" ]            = 0; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
  // 
  // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 0;
}

//-----------------------------------------------------------------------------
// Book MVA methods
//
// Please lookup the various method configuration options in the corresponding cxx files, eg:
// src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
// it is possible to preset ranges in the option string in which the cut optimisation should be done:
// "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable
//-----------------------------------------------------------------------------
void book_methods(TMVA::Factory* factory, TMVA::DataLoader *dataloader, std::map<std::string,int>& Use) {

  // Cut optimisation
  if (Use["Cuts"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

  if (Use["CutsD"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsD",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

  if (Use["CutsPCA"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsPCA",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

  if (Use["CutsGA"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsGA",
			 "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

  if (Use["CutsSA"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsSA",
			 "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  // Likelihood ("naive Bayes estimator")
  if (Use["Likelihood"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood",
			 "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

  // Decorrelated likelihood
  if (Use["LikelihoodD"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodD",
			 "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

  // PCA-transformed likelihood
  if (Use["LikelihoodPCA"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
			 "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

  // Use a kernel density estimator to approximate the PDFs
  if (Use["LikelihoodKDE"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodKDE",
			 "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

  // Use a variable-dependent mix of splines and kernel density estimator
  if (Use["LikelihoodMIX"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodMIX",
			 "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

  // Test the multi-dimensional probability density estimator
  // here are the options strings for the MinMax and RMS methods, respectively:
  //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
  //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
  if (Use["PDERS"])
    factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERS",
			 "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

  if (Use["PDERSD"])
    factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSD",
			 "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

  if (Use["PDERSPCA"])
    factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSPCA",
			 "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

  // Multi-dimensional likelihood estimator using self-adapting phase-space binning
  if (Use["PDEFoam"])
    factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoam",
			 "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

  if (Use["PDEFoamBoost"])
    factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoamBoost",
			 "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

  // K-Nearest Neighbour classifier (KNN)
  if (Use["KNN"])
    factory->BookMethod(dataloader,  TMVA::Types::kKNN, "KNN",
			"H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

  // H-Matrix (chi2-squared) method
  if (Use["HMatrix"])
    factory->BookMethod( dataloader, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

  // Linear discriminant (same as Fisher discriminant)
  if (Use["LD"])
    factory->BookMethod( dataloader, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  // Fisher discriminant (same as LD)
  if (Use["Fisher"])
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  // Fisher with Gauss-transformed input variables
  if (Use["FisherG"])
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

  // Composite classifier: ensemble (tree) of boosted Fisher classifiers
  if (Use["BoostedFisher"])
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "BoostedFisher", 
			 "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

  // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
  if (Use["FDA_MC"])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MC",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

  if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GA",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

  if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_SA",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  if (Use["FDA_MT"])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

  if (Use["FDA_GAMT"])
    factory->BookMethod(dataloader,  TMVA::Types::kFDA, "FDA_GAMT",
			"H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

  if (Use["FDA_MCMT"])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MCMT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

  // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
  if (Use["MLP"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=N" );

  if (Use["MLPBFGS"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

  if (Use["MLPBNN"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

  // CF(Clermont-Ferrand)ANN
  if (Use["CFMlpANN"])
    factory->BookMethod( dataloader, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

  // Tmlp(Root)ANN
  if (Use["TMlpANN"])
    factory->BookMethod(dataloader,  TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

  // Support Vector Machine
  if (Use["SVM"])
    factory->BookMethod(dataloader,  TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

  // Boosted Decision Trees
  if (Use["BDTG"]) // Gradient Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
			 "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
			 "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );


  if (Use["BDTB"]) // Bagging
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
			 "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

  if (Use["BDTD"]) // Decorrelation + Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
			 "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

  if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTMitFisher",
			 "!H:!V:NTrees=50:nEventsMin=150:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

  // RuleFit -- TMVA implementation of Friedman's method
  if (Use["RuleFit"])
    factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
			 "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
}
