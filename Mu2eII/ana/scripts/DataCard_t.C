// Object to hold information about datasets

struct DataCard_t {
  float scale_;
  bool  isBeam_; //if beam related, scale by POT, else scale by livetime
  bool  isOneBatch_;
  bool  isSignal_;
  TString filename_;
  TString folderpath_;
  TString label_;
  int     color_;
  int     setOffset_; //in case using a different selection to get shape
  DataCard_t() : scale_(1.), isBeam_(true), isOneBatch_(true), isSignal_(false), filename_(""), folderpath_(""), label_(""), color_(kRed), setOffset_(0) {}
  DataCard_t(float scale) : DataCard_t() { scale_=scale; }
  DataCard_t(TString fname, TString fpath, float scale)  : DataCard_t(scale) { filename_=fname; folderpath_=fpath; }
  DataCard_t(TString fname, TString fpath, TString label, float scale)  : DataCard_t(fname,fpath,scale) { label_=label; }
  DataCard_t(bool isOneBatch, TString fname, TString fpath, TString label, float scale)  : DataCard_t(fname,fpath,label,scale) { isOneBatch_=isOneBatch; }
  DataCard_t(bool isOneBatch, TString fname, TString fpath, TString label, float scale, bool isSignal)  : DataCard_t(isOneBatch,fname,fpath,label,scale) { isSignal_=isSignal; }
  DataCard_t(bool isOneBatch, TString fname, TString fpath, TString label, float scale, bool isSignal, bool isBeam)  : 
    DataCard_t(isOneBatch,fname,fpath,label,scale,isSignal) { isBeam_=isBeam; }
  DataCard_t(bool isOneBatch, TString fname, TString fpath, TString label, float scale, bool isSignal, bool isBeam, int color)  : 
    DataCard_t(isOneBatch,fname,fpath,label,scale,isSignal,isBeam) { color_=color; }
  DataCard_t(bool isOneBatch, TString fname, TString fpath, TString label, float scale, bool isSignal, bool isBeam, int color, int setOffset)  : 
    DataCard_t(isOneBatch,fname,fpath,label,scale,isSignal,isBeam,color) { setOffset_=setOffset; }
};

