//
// SubRun product: the number of muon captures (and, in simulation, stops and POT) in the
// stopping target that one subrun represents. It is the normalisation the CE analysis needs
// (R_mue = N_signal / (N_cap * eps_CE)), delivered per subrun so EventNtuple can sum it over
// whatever subruns an analysis selects.
//
// One product for simulation and data, so downstream code reads both the same way:
//
//   source MCTruth : expected counts from the true beam intensity, filled in the mixing job
//                    by MuonCaptureNormMC:
//                      nStop = sum over microbunches of PBI * stops/POT
//                      nCap  = nStop * captureFraction
//                    pot() is the summed true PBI. The counts are expectation values, so the
//                    statistical uncertainty is zero; the realised Poisson count differs by
//                    about sqrt(nCap), far below any measurement precision.
//
//   source STM / CAPHRI / Combined : N_cap measured from the capture photon lines and read
//                    from the DB. There is no POT measurement in data, so pot() and nStop()
//                    stay at -1 (hasPOT() false). POT is only ever a simulation quantity.
//
// Uncertainties stored per subrun are those that are independent between subruns
// (statistics, and any uncorrelated systematic). Systematics common to every subrun
// (line yields, detector efficiency model, f_cap) are NOT stored here: they belong to the
// method/version and enter the CE fit once, as a single nuisance parameter.
//
// Like GenEventCount and StageNormalization there is no aggregate(): a subrun must be
// produced whole by one job, which is how Mu2e mixing jobs are run.
//
#ifndef DataProducts_MuonCaptureNormalization_hh
#define DataProducts_MuonCaptureNormalization_hh

#include <cmath>
#include <cstdint>
#include <string>

namespace mu2e {

  class MuonCaptureNormalization {
  public:
    // Stored as an int for ROOT I/O; use source() for the typed value.
    enum Source : int { unknown = 0, MCTruth = 1, STM = 2, CAPHRI = 3, Combined = 4 };

    MuonCaptureNormalization() = default;

    MuonCaptureNormalization(Source source, double nCap, double nCapErrStat, double nCapErrSyst,
                             double nStop, double pot, uint64_t nMicrobunches,
                             std::string provenance) :
      source_(source), nCap_(nCap), nCapErrStat_(nCapErrStat), nCapErrSyst_(nCapErrSyst),
      nStop_(nStop), pot_(pot), nMicrobunches_(nMicrobunches), provenance_(std::move(provenance)) {}

    Source   source()        const { return static_cast<Source>(source_); }
    double   nCap()          const { return nCap_; }
    double   nCapErrStat()   const { return nCapErrStat_; }
    double   nCapErrSyst()   const { return nCapErrSyst_; } // uncorrelated between subruns
    double   nStop()         const { return nStop_; }       // -1 if not known (data)
    double   pot()           const { return pot_; }         // -1 if not known (data): MC only
    uint64_t nMicrobunches() const { return nMicrobunches_; }
    // free text: how the numbers were obtained, e.g. the SimEfficiencies2 tags and f_cap for
    // MC truth, or the method and DB version for a measurement
    std::string const& provenance() const { return provenance_; }

    bool hasPOT()   const { return pot_ >= 0.; }
    bool isMCTruth() const { return source() == MCTruth; }
    bool valid()    const { return source() != unknown && nCap_ >= 0.; }

    // Sum two subruns of the same kind. Counts add; the per-subrun uncertainties are
    // independent by construction, so both add in quadrature. Mixing sources gives an
    // invalid (unknown-source) result rather than a silently meaningless number.
    MuonCaptureNormalization& operator+=(MuonCaptureNormalization const& other) {
      if (source_ != other.source_) source_ = unknown;
      nCap_        += other.nCap_;
      nCapErrStat_  = std::hypot(nCapErrStat_, other.nCapErrStat_);
      nCapErrSyst_  = std::hypot(nCapErrSyst_, other.nCapErrSyst_);
      nStop_        = (nStop_ >= 0. && other.nStop_ >= 0.) ? nStop_ + other.nStop_ : -1.;
      pot_          = (hasPOT() && other.hasPOT()) ? pot_ + other.pot_ : -1.;
      nMicrobunches_ += other.nMicrobunches_;
      if (provenance_ != other.provenance_) provenance_ += " | " + other.provenance_;
      return *this;
    }

  private:
    int         source_        = unknown;
    double      nCap_          = -1.;
    double      nCapErrStat_   = 0.;
    double      nCapErrSyst_   = 0.;
    double      nStop_         = -1.;
    double      pot_           = -1.;
    uint64_t    nMicrobunches_ = 0;
    std::string provenance_;
  };

}
#endif
