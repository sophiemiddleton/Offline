#ifndef __Mu2eII_ana_TrackPar_t__
#define __Mu2eII_ana_TrackPar_t__

#include "Stntuple/ana/ParBase_t.hh"

class TStnTrackID;
class TEmuLogLH;

namespace Mu2eII {
class TrackPar_t : public ParBase_t {
public:
  int           fNHPl;
  int           fNEPl;
  int           fNDPl;
  int           fDiskID;               // 0: first disk, 1:second disk
  int           fFitType;              // index for applying corrections
  int           fTrqMvaIndex;          // index of the trk qual MVAL 0=PAR, 1=DAR

  TStnTrack*    fTrack;                // pointer to its own track - convenience
  TStnCluster*  fCluster;              // pointer to the calorimeter cluster
  SimPar_t*     fSimPar;               // pointer to simparticle 
  TStnTrackID*  fTrackID[50];          //
  TEmuLogLH*    fLogLH;

  int           fIDWord[20];

  double        fTrqMvaOut[20];			// outputs of different MVA classifiers; 0: Dave/Andy MVA; 1: CPR

  double        fPidMvaOut[20];			// 
  
  double        fP   ; 			// momentum, corrected to set DPF at 0
  double        fPStOut;                      // MC true mom on exit from ST
  double        fPFront;                      // MC true mom at the tracker front
  float         fDpF ;			// tracker-only resolution (corrected)
  float         fXDpF;                        // delta(Pfront)/sigmaP
  float         fDp0 ;
  float         fDp2 ;
  float         fDpFSt;
//-----------------------------------------------------------------------------
// different weights are used by different modules, so in principle could use 
// the same variable for different weights. Just don't want to do that
//-----------------------------------------------------------------------------
  double        fDioLOWt;			// DIO LO weight
  double        fDioLLWt;			// DIO LL weight (with leading log rad corrections)
  double        fRPCTimeWt;			// RPC timing weight (don't want to use the same memory for different things)
  double        fRMCEnergyWt;		// RMC weight (energy spectrum of the generated photons)
  double        fLumWt;			// luminosity weight
  double        fTotWt;			// total weight
  double        fTotWtRC;			// total weight with rad corrections

  double        fDtZ0;			// delta(T) at z=0
  double        fDtBack;			// delta(T) at z=Z(TT_Back)
  double        fXtZ0;

  double        fEcl;
  int           fNCrystals;
  double        fSeedFr;                // fraction of energy in the seed crystal
  double        fEp;
  double        fDrDzCal;			// dr/dz at the calorimeter intersection point
  double        fDtClZ0;
  double        fEDiskP;
  double        fECalP;
  double        fDx;
  double        fDy;
  double        fDz;
  double        fDt;
  double        fDu;		         	// rotated residuals
  double        fDv;
  double        fChi2Tcm;
  double        fChi2XY;
  double        fChi2T;
  double        fPath;
  
  double        fSinTC;			// angle between the cluster and the track
  double        fDrTC;
  double        fSInt;

  double        fTMean;			    // mean time over the track hits
//-----------------------------------------------------------------------------
// CRV-related parameters
//-----------------------------------------------------------------------------
  int           fCRVSector;
  int           fTwoEndBars;
  float         fXCorrected;                // Corrected X position in top Crv bars

  double        fDtCRVCorrProp;             // CRV time residual using corrections for propagation only
  double        fDtCRVCorrTof;              // CRV time residual using corrections for TOF only
  double        fDtCRVCorr;                 // CRV time residual using corrections for TOF and propagation
  float         fCRVStubSlope;
  float         fCRVStubSlopeMCProduct;
  double        fDtCRV;			    // minimal time residual (CRV)
  double        fZCRV;                      // Z-coordinate of the closest in time CRV stub candidate
  double        fDtCRV2;	            // minimal time residual (CRV)

//upstream track pars
  double        fUeT0;
  double        fUeP0;
  double        fUeDtCRV;
  double        fDtUe;
  double        fDpUe;
  double        fBounceDt;
  double        fBounceDp;
  double        fBounceDchi;

  int           fDNhitsUe;

//-----------------------------------------------------------------------------
// TrkCaloHit parameters
//-----------------------------------------------------------------------------
  double        fTchTime;
  double        fTchPath;
  double        fTchDx;
  double        fTchDy;
  double        fTchDr;
  double        fTchDz;
  double        fTchDt;
};
}
#endif
