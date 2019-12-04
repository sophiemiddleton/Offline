//Author: S Middleton
// Purpose: Holds functions for the kalman fitting of Cosmic Tracks 
// Date: Dec 2019

#include "CosmicReco/inc/CosmicKalFit.hh"
#include "CosmicReco/inc/CosmicTrkUtils.hh"

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/Mu2eMaterial.hh"
#include "TrackerConditions/inc/Mu2eDetector.hh"
// utiliites
#include "GeneralUtilities/inc/Angles.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "TrkReco/inc/TrkDef.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
// data
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
// Diagnostic data objects
#include "TrkReco/inc/CosmicKalFitData.hh"

// BaBar
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
//#include "BTrk/TrkBase/TrkHelixUtils.hh" //TODO
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
// Mu2e BaBar
#include "BTrkData/inc/TrkStrawHit.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
// root
#include "TH1F.h"
#include "TTree.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;

namespace mu2e
{
  class CosmicKalFit : public art::EDProducer
  {
 _minnstraws
}

void CosmicKalFit::MakeTrack(StrawResponse::cptr_t srep, 
			 Mu2eDetector::cptr_t detmodel,
			 CosmicKalFitData& kalData){

    if(fitable(*kalData.cosmicKalSeed)){
      
      double flt0 = kalData.cosmicKalSeed->flt0();
      auto kseg = kalData.cosmicKalSeed->nearestSegment(flt0);
      if(kseg->fmin() > kseg->localFlt(flt0) || 
	kseg->fmax() < kseg->localFlt(flt0) ){
	
      }
      
      HepVector pvec(4,0);
      HepSymMatrix pcov(4,0);
      kseg->cosmic().hepVector(pvec);
      kseg->cosmic_cov().symMatrix(pcov);
     
      CosmicLineTraj htraj(pvec,pcov);
      kalData.cosmicTraj = &htraj;
      
      TrkStrawHitVector tshv;
      makeTrkStrawHits(srep,kalData, tshv);
      
      std::vector<DetIntersection> detinter;
      if(_matcorr)makeMaterials(detmodel, tshv,*kalData.cosmicTraj,detinter);
   
      std::vector<TrkHit*> thv(0);
      for(auto ihit = tshv.begin(); ihit != tshv.end(); ++ihit){
        thv.push_back(*ihit);
      }
      
      TrkT0 t0(kalData.cosmicKalSeed->t0()); 
      kalData.krep = new KalRep(htraj, thv, detinter, *this, kalData.cosmicKalSeed->particle(), t0, flt0); //TODO KalRep
      assert(kalData.krep != 0);
    
      kalData.krep->addHistory(TrkErrCode(),"KalFit creation");

      TrkErrCode fitstat = fitTrack(detmodel,kalData);
      kalData.krep->addHistory(fitstat,"KalFit fit");

      if(fitstat.success()){
	fitstat = extendFit(kalData.krep);
	kalData.krep->addHistory(fitstat,"KalFit extension");
      }
    }
  }

CosmicKalFit::fitable(CosmicKalSeed const& kseed){
    return kseed.segments().size() > 0 && kseed.hits().size() >= _minnstraws;
  }

  void CosmicKalFit::makeTrkStrawHits(StrawResponse::cptr_t srep, CosmicKalFitData& kalData, TrkStrawHitVector& tshv ) {

    std::vector<TrkStrawHitSeed>const hseeds = kalData.cosmicseed->trkstrawhits();
    CosmicLineTraj const htraj = *kalData.cosmicTraj;
    
    for(auto ths : hseeds ){
      size_t index = ths.index();
      const ComboHit& strawhit(kalData.chcol->at(index));
      const Straw& straw = _tracker->getStraw(strawhit.strawId());
      TrkStrawHit* trkhit = new TrkStrawHit(srep,strawhit,straw,ths.index(),ths.t0(),ths.trkLen(), _maxpull,_strHitW);
      assert(trkhit != 0);
   
      trkhit->setAmbig(ths.ambig());
     
      TrkErrCode pstat = trkhit->updatePoca(&htraj);
      if(pstat.failure()){
        trkhit->setActivity(false);
      }
      tshv.push_back(trkhit);
    }

    std::sort(tshv.begin(),tshv.end(),fcomp());
  }





}


