//Author: S Middleotn
//Date: Dec 2019
//Purpose: this is a slitghly refactorized version of Helix TrkUtilities, they could probably be combined in future.
// Mu2e
#include "CosmicReco/inc/CosmicTrkUtils.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "GeneralUtilities/inc/Angles.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TrkStraw.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
// BTrk
#include "BTrk/TrkBase/CosmicLineTraj.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalMaterial.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "Mu2eBTrk/inc/DetStrawElem.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "DataProducts/inc/XYZVec.hh"
//C++
#include <math.h>
using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
using namespace std;
namespace mu2e {
  namespace CosmicTrkUtils {

    bool CosmicTrack2Traj (CosmicTrack const& track, HepVector& hpvec, unsigned amsign) { 
      bool retval(false);

      if(hpvec.num_row() == CosmicLineTraj::NHLXPRM) {
	// phi0 is the azimuthal angle of the particle velocity vector at the point
	// of closest approach to the origin.  It's sign also depends on the angular
	// momentum.  To translate from the center, we need to reverse coordinates
	hpvec[CosmicLineTraj::phi0Index] = amsign*atan2(track.GetPOCA().x(),track.GetPOCA().y());
	// d0 describes the distance to the origin at closest approach.
	// It is signed by the particle angular momentum WRT the origin.
	// The Helix fit radial bias is anti-correlated with d0; correct for it here.
	hpvec[CosmicLineTraj::d0Index] = amsign*(sqrt(track.GetPOCA().x()*track.GetPOCA().x()+track.GetPOCA().y()+track.GetPOCA().y())); 
	hpvec[CosmicLineTraj::thetaIndex] = asin(track.GetTrackDirection().y()/sqrt(track.GetTrackDirection().Mag2()));
	hpvec[CosmicLineTraj::z0Index] = track.GetPOCA().z();

	retval = true;
      }
      return retval;
    }



    void fillSegment(CosmicLineTraj const& htraj, BbrVectorErr const& momerr,double dflt, KalSegment& kseg) {
      kseg._fmin = htraj.lowRange();
      kseg._fmax = htraj.hiRange();
      kseg._dflt = dflt;
      kseg._cosmic = htraj.parameters()->parameter(); //added cosmic option to Kseg
      kseg._hcov = htraj.parameters()->covariance();
      kseg._mom = momerr.mag();
      Hep3Vector md = momerr.unit();
      HepVector mdir(3);
      for(size_t icor=0;icor<3;++icor) // CLHEP auto-conversion is broken, FIXME!!
	mdir[icor] = md[icor];

      kseg._momerr = sqrt(momerr.covMatrix().similarity(mdir));
    }

    void fillStraws(const KalRep* krep, std::vector<TrkStraw>& tstraws) {
      tstraws.clear();
      // get material sites from the KalRep
      for(auto isite : krep->siteList()){
	if(isite->kalMaterial() != 0) {
	  const KalMaterial* kmat = isite->kalMaterial();
	  const DetStrawElem* detstraw = dynamic_cast<const DetStrawElem*>(kmat->detElem());
	  if(detstraw != 0){
	    // found a straw: create a TrkStraw object from it
	    // i must recompute POCA since the KalMaterial doesn't cache the hit flight FIXME!
	    TrkPoca poca(krep->traj(),kmat->detIntersection().pathlen,*detstraw->wireTraj(),0);
	    TrkStraw tstraw(detstraw->straw()->id(),
	      kmat->detIntersection().dist, //poca.doca(),
	      kmat->detIntersection().pathlen, // poca.flt1(),
	      poca.flt2(),  // not stored in KalMaterial, FIXME!
	      kmat->detIntersection().pathLength(),
	      detstraw->radiationFraction(kmat->detIntersection()),
	      kmat->momFraction(),
	      isite->isActive() );
	    tstraws.push_back(tstraw);
	  }
	}
      }
    }

    void fillStrawHitSeeds(const KalRep* krep,ComboHitCollection const& chits, std::vector<TrkStrawHitSeed>& hitseeds) {
      // extract the TkrStrawHits from the KalRep
      TrkStrawHitVector tshv;
      convert(krep->hitVector(),tshv);
      // loop over the TrkStrawHits and convert them
      for(auto tsh : tshv ) {
      // find the associated ComboHit
	auto const& chit = chits.at(tsh->index());
	// set the flag according to the status of this hit
	StrawHitFlag hflag = chit.flag();
	if(tsh->isActive())hflag.merge(StrawHitFlag::active);
	if(tsh->poca().status().success())hflag.merge(StrawHitFlag::doca);
	// fill the seed.  I have to protect the class from TrkStrawHit to avoid a circular dependency, FIXME!
	TrkStrawHitSeed seedhit(tsh->index(),
	    tsh->hitT0(), tsh->fltLen(), tsh->hitLen(),
	    tsh->driftRadius(), tsh->signalTime(),
	    tsh->poca().doca(), tsh->ambig(),tsh->driftRadiusErr(), hflag, chit);
	hitseeds.push_back(seedhit);
      }
    }

    double overlap(TimeCluster const& tc1, TimeCluster const& tc2) {
      double hover = overlap(tc1._strawHitIdxs,tc2._strawHitIdxs);
      double norm = std::min(tc1.hits().size(),tc2.hits().size());
      double over = norm*hover;
      // add in CaloCluster; count is as much as all the hits
      if(tc1._caloCluster.isNonnull() && tc2._caloCluster.isNonnull()) {
	if(tc1._caloCluster == tc2._caloCluster)
	  over += norm;
	norm *= 2;
      }
      return over/norm;
    }


    double overlap(CosmicKalSeed const& ks1, CosmicKalSeed const& ks2) {
  // translate hit info into a simple index array.  Only count active hits
      SHIV shiv1, shiv2;
      for(auto tshs : ks1.hits()){
	if(tshs.flag().hasAllProperties(StrawHitFlag::active))
	  shiv1.push_back(tshs.index());
      }
      for(auto tshs : ks2.hits()){
	if(tshs.flag().hasAllProperties(StrawHitFlag::active))
	  shiv2.push_back(tshs.index());
      }
      double hover = overlap(shiv1,shiv2);
      double norm = std::min(shiv1.size(),shiv2.size());
      double over = norm*hover;

      return over/norm;
    }

    double overlap(CosmicKalSeed const& ks, CosmicTrackSeed const& hs){
      SHIV shiv1, shiv2;
      for(auto tshs : ks.hits()){
	if(tshs.flag().hasAllProperties(StrawHitFlag::active))
	  shiv1.push_back(tshs.index());
      }
      for(auto hh : hs.hits()){
      // exclude outliers
	if(!hh.flag().hasAnyProperty(StrawHitFlag::outlier))
	  shiv2.push_back(hh.index());
      }
      double hover = overlap(shiv1,shiv2);
      double norm = std::min(shiv1.size(),shiv2.size());
      double over = norm*hover;

      return over/norm;
    }

    double overlap(CosmicTrackSeed const& hs,TimeCluster const& tc) {
      SHIV shiv;
      for(auto hh : hs.hits()){
      // exclude outliers
	if(!hh.flag().hasAnyProperty(StrawHitFlag::outlier))
	  shiv.push_back(hh.index());
      }
      double hover = overlap(shiv,tc.hits());
      double norm = std::min(shiv.size(),tc.hits().size());
      double over = norm*hover;

      return over/norm;
    }

    void countHits(const std::vector<TrkStrawHitSeed>& hits, unsigned& nhits, unsigned& nactive, unsigned& ndouble, unsigned& ndactive, unsigned& nnullambig) {
      nhits = 0; nactive = 0; ndouble = 0; ndactive = 0; nnullambig = 0;
      static StrawHitFlag active(StrawHitFlag::active);
      for (std::vector<TrkStrawHitSeed>::const_iterator ihit = hits.begin(); ihit != hits.end(); ++ihit) {
	++nhits;
	if (ihit->flag().hasAllProperties(active)) {
	  ++nactive;
	  if (ihit->ambig()==0) {
	    ++nnullambig;
	  }
	}

	const auto& jhit = ihit+1;
	const auto& hhit = ihit-1;
	if( (jhit != hits.end() &&
	     jhit->flag().hasAllProperties(active) &&
	     jhit->strawId().getPlane() == ihit->strawId().getPlane() &&
	     jhit->strawId().getPanel() == ihit->strawId().getPanel() ) ||
	    (hhit >= hits.begin() &&
	     hhit->flag().hasAllProperties(active) &&
	     hhit->strawId().getPlane() == ihit->strawId().getPlane() &&
	     hhit->strawId().getPanel() == ihit->strawId().getPanel() )
	    ) {
	  ++ndouble;
	  if (ihit->flag().hasAllProperties(StrawHitFlag::active)) {
	    ++ndactive;
	  }
	}
      }
    }

    double chisqConsistency(const KalRep* krep) {
      return ChisqConsistency(krep->chisq(),krep->nDof()-1).significanceLevel();
    }
  } // TrkUtilities
}// mu2e
