// Data structures used to dump info and then re-sample particles from
// a ROOT tree in some multi-stage job configurations.
//
// Andrei Gaponenko, 2015

#ifndef GeneralUtilities_inc_RSNTIO_hh
#define GeneralUtilities_inc_RSNTIO_hh

namespace mu2e {
  namespace IO {

    //================================================================
    struct StoppedParticleF {
      float x;
      float y;
      float z;
      float t;
      StoppedParticleF() : x(), y(), z(), t() {}

      static const char *branchDescription() {
        return "x/F:y/F:z/F:time/F";
      }

      static unsigned numBranchLeaves() { return 4; }
    };

    //================================================================
    struct StoppedParticleTauNormF {
      float x;
      float y;
      float z;
      float t;
      float px;
      float py;
      float pz;
      float pt;
      float code;
      float id;
      float tauNormalized;
      //float productionP;
      //float productionCosTheta;
      StoppedParticleTauNormF() : x(), y(), z(),  t(), px(), py(), pz(), pt(), code(), id(), tauNormalized(){}//, productionP(), productionCosTheta() {}

      static const char *branchDescription() {
        return "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:pt/F:code/F:id/F:tauNormalized/F";//:productionP/F:productionCosTheta/F";
      }

      static unsigned numBranchLeaves() { return 11;}
    };

    //================================================================
    struct InFlightParticleD {
      double x;
      double y;
      double z;
      double time;
      double px;
      double py;
      double pz;
      int    pdgId;

      InFlightParticleD() : x(), y(), z(), time(), px(), py(), pz(), pdgId() {}

      static const char *branchDescription() {
        return "x/D:y/D:z/D:time/D:px/D:py/D:pz/D:pdgId/I";
      }

      static unsigned numBranchLeaves() { return 8; }
    };

    //================================================================

  } // IO
} // mu2e

#endif/*GeneralUtilities_inc_RSNTIO_hh*/
