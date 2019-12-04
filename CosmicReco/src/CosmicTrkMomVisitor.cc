//--------------------------------------------------------------------------
// File and Version Information:
//
// Description:
//
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Justin Albert, Steve Schaffner
//
//------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/CosmicTrkMomVisitor.hh"
#include "BTrk/TrkBase/TrkSimpTraj.hh"


//------------------------------------------------------------------------
CosmicTrkMomVisitor::~CosmicTrkMomVisitor() {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
CosmicTrkMomVisitor::CosmicTrkMomVisitor(const TrkSimpTraj& theTraj) {
//------------------------------------------------------------------------
	theTraj.visitAccept(this);
}

//------------------------------------------------------------------------
void
CosmicTrkMomVisitor::trkVisitHelixTraj(const HelixTraj* theTraj) {
//------------------------------------------------------------------------
// set the "array"

  _ht = theTraj;
  _ct = 0;
  _nt = 0;
  _cos = 0;
}

//------------------------------------------------------------------------
void
CosmicTrkMomVisitor::trkVisitCircleTraj(const TrkCircleTraj* theTraj) {
//------------------------------------------------------------------------
// set the "array"

  _ht = 0;
  _ct = theTraj;
  _nt = 0;
  _cos =0;
}

//------------------------------------------------------------------------
void
CosmicTrkMomVisitor::trkVisitNeutTraj(const NeutTraj* theTraj) {
//------------------------------------------------------------------------
// set the "array"

  _ht = 0;
  _ct = 0;
  _nt = theTraj;
  _cos = 0;
}


//------------------------------------------------------------------------
void
CosmicTrkMomVisitor::trkVisitCosmicLineTraj(const CosmicLineTraj* theTraj) {
//------------------------------------------------------------------------
// set the "array"

  _ht = 0;
  _ct = 0;
  _nt = 0;
  _cos = theTraj;
}

