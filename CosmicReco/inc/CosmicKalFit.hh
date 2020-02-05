//Author: S Middleton
//Date: Dec 2019

#ifndef CosmicReco_CosmicKalFit_HH
#define CosmicReco_CosmicKalFit_HH

#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

// data
#include "RecoDataProducts/inc/ComboHit.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalContext.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BField/BField.hh"
// Mu2e objects
#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "TrkReco/inc/AmbigResolver.hh"
#include "CosmicReco/inc/CosmicKalFitData.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/Mu2eDetector.hh"
#include "TrkReco/inc/TrkPrintUtils.hh"

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// C++
#include <array>
// Framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
namespace mu2e 
{

  class CosmicKalFit : public KalContext
  {
  public:
	
	#ifndef __GCCXML__
    	explicit CosmicKalFit(fhicl::ParameterSet const&);
	#endif/*__GCCXML__*/
	virtual ~CosmicKalFit();

	void MakeTrack(StrawResponse::cptr_t strawResponse, Mu2eDetector::cptr_t detmodel, CosmicKalFitData& kalData);

	virtual const TrkVolume* trkVolume(trkDirection trkdir) const;
	BField const& bField() const;

  	TrkErrCode FitIteration  (Mu2eDetector::cptr_t detmodel, CosmicKalFitData& kalData,int iter); 

	void       setTracker(const Tracker*  Tracker) { _tracker = Tracker; }
	bool       updateT0 (CosmicKalFitData& kalData, int    iter);
	bool       hit_time  (TrkHit* hit, HitT0& hitT0);
	HitT0      krep_hitT0(KalRep* krep, const TrkHit*hit);

//TrkPrintUtils*  printUtils() { return _printUtils; }
  private:

	int _debug;		    // debug level
	float _maxpull;   // maximum pull in TrkHit 
	float _strHitW;
	std::vector<bool> _updatet0; // update t0 ieach iteration?
	unsigned _minnstraws;   // minimum # staws for fit
	std::vector<double> _herr; //hiterror
	unsigned _maxIterations; //max number of iterations
	std::vector<bool> _addmaterial; // look for additional materials along the track
	const mu2e::Tracker*    _tracker;     // straw tracker geometry
	//TrkTimeCalculator _ttcalc;
	int    _annealingStep;
	mutable BField* _bfield;
 	TrkPrintUtils*  _printUtils;
	bool fitable(CosmicKalSeed const& kseed);
	void initT0(CosmicKalFitData& kalData);

	void MakeTrkStrawHits  (StrawResponse::cptr_t strawResponse, CosmicKalFitData& kalData, TrkStrawHitVector& tshv );

	void MakeMaterials( Mu2eDetector::cptr_t detmodel,TrkStrawHitVector const&, CosmicLineTraj const& traj,std::vector<DetIntersection>& dinter);

	TrkErrCode FitTrack    (Mu2eDetector::cptr_t detmodel, CosmicKalFitData& kalData);
	
  };
}
#endif
