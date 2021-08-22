#ifndef __Mu2eII_ana_SimPar_t__
#define __Mu2eII_ana_SimPar_t__

namespace Mu2eII {
  struct SimPar_t {
    TGenParticle*  fGenp;		// 
    TSimParticle*  fParticle;		// pointer to the signal particle
    TStepPointMC*  fTFront;		// VD hit at the tracker front
    TStepPointMC*  fTMid;		// VD hit at the tracker center
    TStepPointMC*  fTBack;		// VD hit at the tracker back (+1610 mm)
    int            fStage;              // last stage the particle reached
    double         fWeight;             // weight , either 1 or d2sigma/dp/dcosth
    double         fEleE;               // electron energy
  };
}
#endif
