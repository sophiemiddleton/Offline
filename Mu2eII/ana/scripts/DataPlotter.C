// Object to hold information about datasets and create plots
#include "DataCard_t.C"
#include "PlottingCard_t.C"

class DataPlotter {
public:
  DataPlotter() { }
  ~DataPlotter() {
    for(unsigned index = 0; index < files_.size(); ++index) {
      delete files_[index];
    }
  }

  //create a stack of the background processes
  THStack* GetStack(PlottingCard_t &card) {
    unsigned ndatasets = files_.size();
    std::map<TString, TH1F*> hmap; //map label -> hist to combine like labels
    std::vector<TString> labels; //store labels to preserve dataset order
    //loop through datasets
    for(unsigned index = 0; index < ndatasets; ++index) {
      if(isSignal_[index]) continue;
      TH1F* h = (TH1F*) files_[index]->Get(Form("%s_%i/%s", card.type_.Data(), card.set_+setOffsets_[index], card.hist_.Data()));
      if(!h) {
	std::cout << "Histogram " << card.type_.Data() << "_" << card.set_+setOffsets_[index] << "/" << card.hist_.Data() 
		  << " not found in index " << index << " = " << labels_[index].Data() << endl;
	continue;
      }
      //clone to create new object
      h = (TH1F*) h->Clone();
      //apply proper normalization factors
      h->Scale(scales_[index]);

      int era = !isOneBatch_[index]; //scale for the given period
      if(isBeam_[index]) h->Scale(lumi_[era]);
      else               h->Scale(livetime_[era]);
      
      //rebin if needed
      if(card.rebin_ > 1) h->Rebin(card.rebin_); 

      //see if the label is already in the map, add to existing histogram if it is
      auto itr = hmap.find(labels_[index]);
      if(itr != hmap.end()) { //already exists
	itr->second->Add(h);
	h = itr->second;
      } else {hmap[labels_[index]] = h; labels.push_back(labels_[index]);}

      //configure histogram for drawing
      h->SetName(Form("h_%s_%s_%i", labels_[index].Data(), card.hist_.Data(), card.set_));
      h->SetTitle(labels_[index].Data());
      h->SetLineWidth(2);
      h->SetFillColor(colors_[index]);
      h->SetLineColor(colors_[index]+1); //slightly different line color
    }

    //add histograms to a stack
    THStack* hstack = new THStack(Form("s_%s_%i", card.hist_.Data(), card.set_), "");
    for(unsigned index = 0; index < labels.size(); ++index) {
      auto h = hmap[labels[index]];
      if(print_stats_) {
	double res, err;
	int bin1(1), bin2(h->GetNbinsX());
	if(card.xmin_ < card.xmax_) {
	  bin1 = h->FindBin(card.xmin_);
	  bin2 = h->FindBin(card.xmax_);
	}
	res = h->IntegralAndError(bin1, bin2, err);
	std::cout << labels[index].Data() << ": " << res << " +- " << err << std::endl;
      }
      if(h->GetEntries() > 0) //only add if contributes
	hstack->Add(h);
    }

    return hstack;
  }

  //get a list of the signal histograms
  std::vector<TH1F*> GetSignals(PlottingCard_t &card) {
    std::vector<TH1F*> signals;
    unsigned ndatasets = files_.size();
    std::map<TString, TH1F*> hmap; //map label -> hist to combine like labels
    std::vector<TString> labels; //store labels to preserve dataset order

    //loop through datasets
    for(unsigned index = 0; index < ndatasets; ++index) {
      if(!isSignal_[index]) continue; //must be signal
      TH1F* h = (TH1F*) files_[index]->Get(Form("%s_%i/%s", card.type_.Data(), card.set_+setOffsets_[index], card.hist_.Data()));
      if(!h) {
	std::cout << "Histogram " << card.type_.Data() << "_" << card.set_+setOffsets_[index] << "/" << card.hist_.Data() 
		  << " not found in index " << index << " = " << labels_[index].Data() << endl;
	continue;
      }

      h = (TH1F*) h->Clone();
      h->Scale(scales_[index]);
      int era = !isOneBatch_[index]; //scale for the given period
      if(isBeam_[index]) h->Scale(lumi_[era]);
      else               h->Scale(livetime_[era]);

      if(card.rebin_ > 1) h->Rebin(card.rebin_); 

      //see if the label is already in the map, add to existing histogram if it is
      auto itr = hmap.find(labels_[index]);
      if(itr != hmap.end()) { //already exists
	itr->second->Add(h);
	h = itr->second;
      } else {hmap[labels_[index]] = h; labels.push_back(labels_[index]);}

      h->SetName(Form("h_%s_%s_%i", labels_[index].Data(), card.hist_.Data(), card.set_));
      h->SetTitle(labels_[index].Data());
      h->SetLineWidth(3);
      h->SetLineColor(colors_[index]);
    }
    for(unsigned index = 0; index < labels.size(); ++index) {
      auto h = hmap[labels[index]];
      if(print_stats_) {
	double res, err;
	int bin1(1), bin2(h->GetNbinsX());
	if(card.xmin_ < card.xmax_) {
	  bin1 = h->FindBin(card.xmin_);
	  bin2 = h->FindBin(card.xmax_);
	}
	res = h->IntegralAndError(bin1, bin2, err);
	std::cout << labels[index].Data() << ": " << res << " +- " << err << std::endl;
      }
      signals.push_back(h);
    }
    return signals;
  }

  //Add some information to the canvas
  void DrawInfo() {
    TLatex label;
    label.SetNDC();
    label.SetTextFont(72);
    label.SetTextColor(1);

    //add label
    label.SetTextSize(0.06);
    label.SetTextAlign(22);
    label.DrawText(0.29, 0.92, "Mu2e Preliminary");

    //add run info
    label.SetTextSize(0.03);
    label.SetTextAlign(13);
    label.DrawLatex(0.7, 0.975, Form("%.2e POT, %.2e s",lumi_[0]+lumi_[1], livetime_[0]+livetime_[1]));
  }

  //Plot a background stack + signal histograms
  TCanvas* PlotStack(PlottingCard_t card) {
    //get the background stack
    THStack* hstack = GetStack(card);
    if(!hstack) {
      std::cout << "Background stack not found!\n";
      return NULL;
    }
    //get the signal histograms
    std::vector<TH1F*> signals = GetSignals(card);
    TCanvas* c = new TCanvas(Form("c_%s_%i", card.hist_.Data(), card.set_), Form("c_%s_%i", card.hist_.Data(), card.set_), canvas_x_, canvas_y_);
    //draw the stack, preserving drawing info, add error bars
    hstack->Draw("hist E1 noclear");
    //track maximum histogram value while drawing signals
    double mx(hstack->GetMaximum());
    for(auto h : signals) {h->Draw("hist E1 same"); mx = max(mx, h->GetMaximum());}

    //set axis ranges
    if(card.xmin_ < card.xmax_)
      hstack->GetXaxis()->SetRangeUser(card.xmin_, card.xmax_);
    if(card.ymin_ < card.ymax_) {
      hstack->SetMaximum(card.ymax_);
      hstack->SetMinimum(card.ymin_);
    } else
      hstack->SetMaximum(((logy_ > 0) ? 2.*mx : 1.2*mx));

    hstack->GetXaxis()->SetTitle(card.xlabel_.Data());
    TString ylabel = card.ylabel_;
    if(add_bin_width_) ylabel = Form("Entries / %.2f %s", ((TH1F*) hstack->GetStack()->Last())->GetBinWidth(1), ylabel.Data());
    hstack->GetYaxis()->SetTitle(ylabel.Data());

    //add a legend
    TLegend* leg = new TLegend(0.65, 0.68, 0.95, 0.95, "", "LP brNDC");
    for(auto h : signals) leg->AddEntry(h);
    for(auto h : *hstack->GetHists()) leg->AddEntry(h);
    leg->Draw("same");
    
    //Draw info on canvas
    DrawInfo();

    //configure canvas
    if(logy_ > 0) c->SetLogy();
    c->SetGrid();
    c->SetTopMargin(0.05); c->SetRightMargin(0.05);
    c->Modified(); c->Update();
    return c;
  }

  //print a stack canvas
  TCanvas* PrintStack(PlottingCard_t card) {
    TCanvas* c = PlotStack(card);
    if(!c) return c;
    gSystem->Exec("[ ! -d figures ] && mkdir figures"); //make the directory if needed
    TString filename = Form("figures/stack_%s_%s%s_%i.png", card.type_.Data(), card.hist_.Data(), 
			    (logy_ > 0) ? "_log" : "", card.set_);
    c->SaveAs(filename.Data());
    return c;
  }

  //add a data file to plot
  Int_t AddFile(DataCard_t &card) {
    Int_t status(0);
    TFile* f = TFile::Open(card.filename_.Data(), "READ");
    if(!f) {return 1;}
    if(!(f->Get(card.folderpath_.Data()))) {
      std::cout << card.folderpath_.Data() << " not found in " << card.filename_.Data()
		<< " in DataPlotter::" << __func__ << endl;
      return 2;
    }
    files_.push_back((TFile*) f->Get(card.folderpath_.Data()));
    labels_.push_back(card.label_.Data());
    scales_.push_back(card.scale_);
    isBeam_.push_back(card.isBeam_);
    isOneBatch_.push_back(card.isOneBatch_);
    isSignal_.push_back(card.isSignal_);
    colors_.push_back(card.color_);
    setOffsets_.push_back(card.setOffset_);
    return status;
  }

  //add a list of files to plot
  Int_t AddFiles(std::vector<DataCard_t> cards) {
    Int_t status(0);
    for(auto card : cards) status += AddFile(card);
    return status;
  }

  //data members
  std::vector<TFile*>  files_; //list of histogram files
  std::vector<TString> labels_; //label for data file, combine like labels
  std::vector<double>  scales_; //normalization for histograms
  std::vector<bool>    isBeam_; //whether or not scales with POT or livetime
  std::vector<bool>    isOneBatch_; //1-batch mode or 2-batch mode histogram
  std::vector<bool>    isSignal_; //whether or not is signal histogram
  std::vector<int>     colors_; //color to draw with
  std::vector<int>     setOffsets_; //set number offset if using different set for shape
  double lumi_[2]; //0: 1 batch mode lumi 1: 2 batch mode lumi
  double livetime_[2]; //0: 1 batch mode time 1: 2 batch mode time
  double canvas_x_ = 800; //canvas size
  double canvas_y_ = 600;
  int    logy_ = 0; //whether to draw in log Y or not
  bool   add_bin_width_ = true; //add Y axis label with width
  bool   print_stats_ = false; //print statistics estimate in plotted window
};
