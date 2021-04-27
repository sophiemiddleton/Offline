#ifndef __Mu2eII_ana_TrackHist_t_hh
#define __Mu2eII_ana_TrackHist_t_hh

namespace Mu2eII {

  struct TrackHist_t {
    //    TH1F*    fInstLumi[6];              // MC truth lumi info: 0: nominal 1: undo batch weight 2: opposite batch weight
    TH1F*    fP[3];			// total momentum, 3 hists with different binning
    TH1F*    fPdio;                     // 103.85-104.9MeV/c e-
    TH1F*    fPdioMC;                    // same momentum window, generated particle momentum
    TH1F*    fPdioWeight;                // same window, event weight
    TH1F*    fP0;
    TH1F*    fP2;
    TH1F*    fPDio;                      // with the LL DIO weight

    TH1F*    fFitMomErr;
    TH1F*    fPFront;                   // momentum of the corresponding MC particle at the tracker front
    TH1F*    fDpFront;
    TH1F*    fXDpF;                     // DpF/MomErr
    TH1F*    fDpFront0;
    TH1F*    fDpFront2;
    TH1F*    fPStOut;
    TH1F*    fDpFSt;			// P(TT_Hollow) - P(ST_Out)
    TH2F*    fDpFVsZ1;

    TH1F*    fPt;
    TH1F*    fCosTh;
    TH1F*    fChi2;
    TH1F*    fChi2Dof;
    TH2F*    fPVsGenE;                 // total momentum vs generated energy

    TH1F*    fNActive;
    TH1F*    fNaFract;
    TH1F*    fDNa;
    TH1F*    fNWrong;			// MC-only histogram: N(hits) with wrong drift signs
    TH1F*    fNDoublets;
    TH1F*    fNadOverNd;		// fraction of doublets with all hits active
    TH1F*    fNSSD;
    TH1F*    fNOSD;
    TH1F*    fNdOverNa;

    TH1F*    fNssdOverNa;
    TH1F*    fNosdOverNa;
    TH1F*    fNZeroAmb;
    TH1F*    fNzaOverNa;
    TH1F*    fNMatActive;
    TH1F*    fNmaOverNa;
    TH1F*    fNBend;

    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fQ;
    TH1F*    fFitCons[2];		// fit consistency (0 to 1)
    TH1F*    fD0;
    TH1F*    fZ0;
    TH1F*    fTanDip;
    TH1F*    fRMax;
    TH1F*    fDtZ0;			// MC truth: T0-T(MC TMid)
    TH1F*    fXtZ0;                     // pull(dt) at Z=0

    TH1F*    fResid;
    TH1F*    fAlgMask;
					// matching
    TH1F*    fChi2Tcm;
    TH1F*    fChi2XY;
    TH1F*    fChi2T;

    TH1F*    fDt;			// track-cluster residuals
    TH1F*    fDx;
    TH1F*    fDy;
    TH1F*    fDz;
    TH1F*    fDu;
    TH1F*    fDv;
    TH1F*    fPath;

    TH1F*    fECl;
    TH1F*    fSeedFr;
    TH1F*    fNCrystals;
    TH1F*    fEClEKin;
    TH1F*    fEp;
    TH1F*    fDrDzCal;
    TH1F*    fDtClZ0;                   // T(cluster back at Z0)-T_true(Z0)
    TH2F*    fDtClZ0VsECl;              // 
    TH2F*    fDtClZ0VsP;                // 
    TH2F*    fEclVsPath;

    TH2F*    fFConsVsNActive;
    TH1F*    fDaveTrkQual;
    TH1F*    fTrqMvaOut;	        // output of our TRQ MVA
    TH1F*    fDeltaMVA;			// DaveTrkQual-fTrqMvaOut

    TH2F*    fPVsTime;                  // 2D histogram momentum vs time at the center of the tracker (T0) (used for sensitivity calculation)
    TH2F*    fT0vsP;                    // 2D histogram time vs momentum at the center of the tracker (T0) (used for sensitivity calculation)

    TH1F*    fCRVSector;                // CRV sector type
    TH1F*    fDtCRV;                    // delta(T) from CRV stub candidate
    TH1F*    fDtCRV2;                   // delta(T) from CRV stub candidate
    TH2F*    fDtVsZCRV;                 // delta(T) from CRV stub candidate

    TH2F*    fDtVsZCRVCorrProp;         // dT  from CRV stub candidate with bar propagation corrections only
    TH2F*    fDtVsZCRVCorrTof;          // dT  from CRV stub candidate with time of flight corrections only
    TH2F*    fDtVsZCRVCorr;             // dT  from CRV stub candidate with propagation and TOF corrections
    TH2F*    fDt2VsZCRV;                // dT2 vs Z

    TH1F*    fPDGCode;                  // PDG Code of particle which made the most hits
    TH1F*    fDPback;                   // DP at the back of the tracker
    TH1F*    fXCorrected;               // histogram of corrected position from CrvStubPar
    
    TH1F*    fDeUeDt;
    TH1F*    fUeT0;
    TH1F*    fUeP0;
    TH1F*    fDtUe;
    TH1F*    fDpUe;
    TH1F*    fBounceDt;
    TH1F*    fBounceDp;
    TH1F*    fBounceDchi;
    TH1F*    fDNhitsUe;

    TH2F*    fUeDtVsZCRV;

//-----------------------------------------------------------------------------
// TrkCaloHit histograms
//-----------------------------------------------------------------------------
    TH1F*    fTchTime;
    TH1F*    fTchPath;
    TH1F*    fTchDx;
    TH1F*    fTchDy;
    TH1F*    fTchDr;
    TH1F*    fTchDz;
    TH1F*    fTchDt;

    TH2F*    fEpVsPath;
    TH2F*    fEpVsTchDz;
    TH2F*    fTchDtVsDz;
//-----------------------------------------------------------------------------
// MVA-based PID
//-----------------------------------------------------------------------------
    TH1F*    fPidMvaOut;		      // output of the PID MVA
  };
}
#endif
