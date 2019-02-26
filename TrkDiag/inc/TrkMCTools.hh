//
// Namespace for collecting tools used in MC truth evaluation
// Original author: Dave Brown (LBNL) 8/10/2016
//
#ifndef TrkDiag_TrkMCTools_hh
#define TrkDiag_TrkMCTools_hh
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "RecoDataProducts/inc/KalSeed.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"
#include "TrkDiag/inc/CaloClusterInfoMC.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "MCDataProducts/inc/KalSeedMC.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "MCDataProducts/inc/PrimaryParticle.hh"

#include <vector>
#include <functional>
namespace mu2e {
  namespace TrkMCTools {
    /////////////////////////////
    // Straw Hit Level Utilities

    // return the earliest StepPointMC associated with the threshold crossing.
    int stepPoint(art::Ptr<StepPointMC>& sp, StrawDigiMC const& mcdigi);
    // determine if a hit was generated by signal (conversion electron)
    bool CEDigi(StrawDigiMC const& mcdigi);
    // count the number of digis produced by the given particle
    unsigned countDigis(art::Ptr<SimParticle> const& pspp, const StrawDigiMCCollection* mcdigis);
    // determine the sim particle which originated most of the hits specified in the specified set of digis
    unsigned primaryParticle(art::Ptr<SimParticle>& pspp, std::vector<StrawHitIndex> const& hits, const StrawDigiMCCollection* mcdigis);
    // determine the sim particle corresponding to a digi
    unsigned simParticle(art::Ptr<SimParticle>& spp, StrawDigiMC const& digimc);

    /////////////////////////////
    // Track Level Utilities
    struct spcount {
      spcount() : _count(0), _acount(0) {}
      spcount(art::Ptr<SimParticle> const& spp,bool active) : _spp(spp), _count(1), _acount(0) {
	if(active)_acount =1; }
      void append(art::Ptr<SimParticle> const& sp,bool active) { if(sp == _spp){
	++_count; if(active)++_acount; } }
      bool operator ==(art::Ptr<SimParticle> const& sp) const { return _spp == sp; }
      art::Ptr<SimParticle> _spp;
      unsigned _count; // counts all hits
      unsigned _acount; // counts active 
    };
// sort by active hits
    struct spcountcomp : public std::binary_function <spcount, spcount, bool> {
      bool operator() (spcount a, spcount b) { return a._acount > b._acount; }
    };

    typedef StepPointMCCollection::const_iterator MCStepItr;
    struct timecomp : public std::binary_function<MCStepItr,MCStepItr, bool> {
      bool operator()(MCStepItr x,MCStepItr y) { return x->time() < y->time(); }
    };    


    // find associated sim particles to a track.  The first returns a hit-weighted vector of
    // all particles, the second just the one with the most hits
    void findMCTrk(const KalSeed& kseed, art::Ptr<SimParticle>& spp, const StrawDigiMCCollection& mcdigis);
    void findMCTrk(const KalSeed& kseed, std::vector<spcount>& sct, const StrawDigiMCCollection& mcdigis);

    // find steps associated with a given SimParticle ID
    void findMCSteps(const StepPointMCCollection& mcsteps, cet::map_vector_key const& trkid, std::vector<int> const& vids, std::vector<MCStepItr>& steps);

    // count types of hits and digis
    //    void countHits(const KalSeed& kseed, const art::Ptr<SimParticle>& spp, const StrawDigiMCCollection& mcdigis, const double& mingood, int& nactive, int& nhits, int& ngood, int& nambig);
    void countDigis(const KalSeedMC& kseedmc, const KalSeed& kseed, int& ndigi, int& digigood, int& ngood);

    // fill various info structs
    void fillTrkInfoMC(const KalSeedMC& kseedmc, const KalSeed& kseed, TrkInfoMC& trkinfomc);
    void fillTrkInfoMCStep(const KalSeedMC& kseedmc, TrkInfoMCStep& trkinfomcstep, const PrimaryParticle& primary);
    void fillTrkInfoMCStep(const KalSeedMC& kseedmc, TrkInfoMCStep& trkinfomcstep,std::vector<int> const& vids);
    void fillHitInfoMCs(const KalSeedMC& kseedmc, std::vector<TrkStrawHitInfoMC>& tshinfomcs);
    void fillHitInfoMC(const KalSeedMC& kseedmc, TrkStrawHitInfoMC& tshinfomc, const TrkStrawHitMC& tshmc);
    void fillCaloClusterInfoMC(CaloClusterMC const& ccmc, CaloClusterInfoMC& ccimc);

  }
}

#endif
