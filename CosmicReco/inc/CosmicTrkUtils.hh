//Author: S Middleton
//Date: Dec 2019
//Purpose: Convert offlne cosmic track to BTrk Traj (analogous to TrkUtils for helix)
#ifndef TrkReco_CosmicTrkUtils_HH
#define TrkReco_ComicTrkUtils_HH
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include <vector>

class CosmicLineTraj;
class BbrVectorErr;
class KalRep;
class TrkDifPieceTraj;
namespace mu2e {
  class CosmicTrackSeedHelix;
  class KalSegment;
  class CosmicKalSeed;
  class CosmicTrack;
  class TrkStrawHitSeed;
  class TrkStraw;
  class TimeCluster;
  class ComboHitCollection;
  typedef std::vector<StrawHitIndex> SHIV;
  namespace CosmicTrkUtils{
    bool CosmicTrack2Traj (CosmicTrack const& cosmic, CLHEP::HepVector& hpvec, float amsign);
    void CosmicTrackFromMom(CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& mom, double charge, double Bz, CosmicTrack& cosmic);
    void fillSegment(CosmicLineTraj const& traj, BbrVectorErr const& momerr,double dflt, KalSegment& kseg);
    void fillStrawHitSeeds(const KalRep* krep, ComboHitCollection const& chits, std::vector<TrkStrawHitSeed>& hitseeds); 
    void fillStraws(const KalRep* krep, std::vector<TrkStraw>& straws);
    // compute overlap between 2 time clusters
    double overlap(TimeCluster const& tc1, TimeCluster const& tc2);//TODO - no calo cluster needed for us
    double overlap(CosmicKalSeed const& ks1, CosmicKalSeed const& ks2);
    double overlap(CosmicKalSeed const& ks, CosmicTrackSeed const& hs);
    double overlap(CosmicTrackSeed const& hs,TimeCluster const& tc);
    double overlap(SHIV const& shiv1, SHIV const& shiv2); 
    double chisqConsistency(const KalRep* krep);
  }
}
#endif
