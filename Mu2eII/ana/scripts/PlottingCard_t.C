// Object to hold information about datasets

struct PlottingCard_t {
  TString hist_;
  TString type_;
  int     set_;
  int     rebin_;
  double  xmin_;
  double  xmax_;
  double  ymin_;
  double  ymax_;
  TString xlabel_;
  TString ylabel_;

  PlottingCard_t() : hist_(""), type_(""), set_(0), rebin_(1), xmin_(1.), xmax_(-1.), ymin_(1.), ymax_(-1.), xlabel_(""), ylabel_("") {}
  PlottingCard_t(TString hist, TString type) : PlottingCard_t() { hist_=hist; type_=type; }
  PlottingCard_t(TString hist, TString type, int set) : PlottingCard_t(hist,type) { set_=set; }
  PlottingCard_t(TString hist, TString type, int set, int rebin, double xmin, double xmax) : PlottingCard_t(hist,type,set) { rebin_=rebin; xmin_=xmin; xmax_=xmax; }
  PlottingCard_t(TString hist, TString type, int set, int rebin, double xmin, double xmax, double ymin, double ymax) : 
    PlottingCard_t(hist,type,set,rebin,xmin,xmax) { ymin_=ymin; ymax_=ymax; }
  PlottingCard_t(TString hist, TString type, int set, int rebin, double xmin, double xmax, double ymin, double ymax, TString xlabel, TString ylabel) :
    PlottingCard_t(hist,type,set,rebin,xmin,xmax,ymin,ymax) { xlabel_=xlabel; ylabel_=ylabel; }
};

