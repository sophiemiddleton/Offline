#ifndef __Mu2eII_ana_TPidMva_hh__
#define __Mu2eII_ana_TPidMva_hh__

class TBranch;

namespace Mu2eII {
  struct TPidMvaTrainData_t {
    float    fP;          
    float    fEcl;
    float    fSeedFr;
    float    fNCrystals;
    float    fElePath;
    float    fEleTchDt;
    float    fEleTchDz;
    float    fEleTchDr;
    float    fMuoPath;
    float    fMuoTchDt;
    float    fMuoTchDz;
    float    fMuoTchDr;
  };

  struct TPidMvaTrainBranches_t {
    TBranch*   fP;       
    TBranch*   fEcl;     
    TBranch*   fSeedFr;    
    TBranch*   fNCrystals;
    TBranch*   fElePath; 
    TBranch*   fEleTchDt;
    TBranch*   fEleTchDz;
    TBranch*   fEleTchDr;
    TBranch*   fMuoPath; 
    TBranch*   fMuoTchDt;
    TBranch*   fMuoTchDz;
    TBranch*   fMuoTchDr;
  };
}
#endif
