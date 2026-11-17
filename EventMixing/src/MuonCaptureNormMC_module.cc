//
// MuonCaptureNormMC: MC-truth muon capture normalisation per subrun.
//
// Sums the true proton bunch intensity (PBI) over the subrun and writes the expected numbers
// of stopped and captured muons as a MuonCaptureNormalization SubRun product with source
// MCTruth:
//
//   nStop = sum over microbunches of PBI * stops/POT
//   nCap  = nStop * captureFraction
//
// It belongs in the mixing job, next to MuStopPileupMixer, because that is where every
// microbunch of the subrun exists; triggered or split files further downstream would
// undercount. stops/POT is obtained exactly as MixBackgroundFrames obtains the mean number of
// muon-stop pileup events per POT:
//
//   stops/POT = prod_i SimBookkeeper::getEff(simStageEfficiencyTags_i) * prod_j meanEventsPerPOTFactors_j
//
// with simStageEfficiencyTags = [MuBeamCat, MuminusStopsCat] and
// meanEventsPerPOTFactors = [MuminusTargetStopPrescale], so the truth always matches what was
// mixed (the mixer's own list also contains MuStopPileupCat, the detector-step filter
// efficiency of the pileup pool, which must NOT be included here). The efficiency is looked up
// per event, so a change of SimBookkeeper IOV within a subrun is handled.
//
// captureFraction defaults to 1 - PhysicsParams::getDecayFraction(material), the value the
// Pileup generator uses to choose capture or decay for each stop.
//
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/DataProducts/inc/MuonCaptureNormalization.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/SimulationConditions/inc/SimBookkeeper.hh"

namespace mu2e {

  class MuonCaptureNormMC : public art::EDProducer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> pbiTag { Name("protonBunchIntensityTag"),
        Comment("True proton bunch intensity per microbunch"), art::InputTag("PBISim") };
      fhicl::Sequence<std::string> simStageEfficiencyTags { Name("simStageEfficiencyTags"),
        Comment("SimEfficiencies2 tags whose product gives (prescaled) target stops per POT"),
        std::vector<std::string>{"MuBeamCat", "MuminusStopsCat"} };
      fhicl::Sequence<double> meanEventsPerPOTFactors { Name("meanEventsPerPOTFactors"),
        Comment("Extra factors, as for MuStopPileupMixer: the target-stop prescale") };
      fhicl::Atom<std::string> stoppingTargetMaterial { Name("stoppingTargetMaterial"),
        Comment("Material for the default capture fraction"), "Al" };
      fhicl::OptionalAtom<double> captureFraction { Name("captureFraction"),
        Comment("Override captures per stop; default 1 - PhysicsParams decay fraction") };
      fhicl::Atom<int> debugLevel { Name("debugLevel"), Comment("Debug printout level"), 0 };
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit MuonCaptureNormMC(Parameters const& conf);

    void produce(art::Event& event) override;
    void beginSubRun(art::SubRun& sr) override;
    void endSubRun(art::SubRun& sr) override;

  private:
    art::InputTag            pbiTag_;
    std::vector<std::string> effTags_;
    std::vector<double>      extraFactors_;
    double                   captureFraction_;
    int                      debugLevel_;
    std::string              provenance_;
    ProditionsHandle<SimBookkeeper> simbookkeeperH_;

    // per-subrun accumulators
    double   sumPOT_   = 0.;
    double   sumStops_ = 0.;
    uint64_t nMicrobunches_ = 0;
  };

  MuonCaptureNormMC::MuonCaptureNormMC(Parameters const& conf) :
    art::EDProducer{conf},
    pbiTag_(conf().pbiTag()),
    effTags_(conf().simStageEfficiencyTags()),
    extraFactors_(conf().meanEventsPerPOTFactors()),
    debugLevel_(conf().debugLevel())
  {
    double fcap = 0.;
    if (!conf().captureFraction(fcap)) {
      fcap = 1. - GlobalConstantsHandle<PhysicsParams>()->getDecayFraction(conf().stoppingTargetMaterial());
    }
    if (fcap <= 0. || fcap > 1.) {
      throw cet::exception("MuonCaptureNormMC") << "Unphysical capture fraction " << fcap;
    }
    captureFraction_ = fcap;

    std::ostringstream prov;
    prov << "MCTruth: sum(PBI from " << pbiTag_.encode() << ") x [";
    for (size_t i = 0; i < effTags_.size(); ++i) prov << (i ? " x " : "") << effTags_[i];
    prov << "]";
    for (auto f : extraFactors_) prov << " x " << f;
    prov << "; captureFraction " << captureFraction_;
    provenance_ = prov.str();

    produces<MuonCaptureNormalization, art::InSubRun>();
  }

  void MuonCaptureNormMC::beginSubRun(art::SubRun&) {
    sumPOT_ = 0.;
    sumStops_ = 0.;
    nMicrobunches_ = 0;
  }

  void MuonCaptureNormMC::produce(art::Event& event) {
    auto const& pbi = *event.getValidHandle<ProtonBunchIntensity>(pbiTag_);

    SimBookkeeper const& bookkeeper = simbookkeeperH_.get(event.id());
    double stopsPerPOT = 1.;
    for (auto const& tag : effTags_) {
      double eff = bookkeeper.getEff(tag);
      if (eff < 0.) throw cet::exception("MuonCaptureNormMC") << "Unphysical SimStageEfficiency " << tag << " = " << eff;
      stopsPerPOT *= eff;
    }
    for (auto f : extraFactors_) stopsPerPOT *= f;

    sumPOT_   += pbi.intensity();
    sumStops_ += pbi.intensity() * stopsPerPOT;
    ++nMicrobunches_;

    if (debugLevel_ > 1) {
      mf::LogInfo("MuonCaptureNormMC") << "event " << event.id() << " PBI " << pbi.intensity()
                                       << " stops/POT " << stopsPerPOT;
    }
  }

  void MuonCaptureNormMC::endSubRun(art::SubRun& sr) {
    const double nCap = sumStops_ * captureFraction_;
    // Expectation values: no statistical or per-subrun systematic uncertainty.
    auto norm = std::make_unique<MuonCaptureNormalization>(
      MuonCaptureNormalization::MCTruth, nCap, 0., 0., sumStops_, sumPOT_, nMicrobunches_, provenance_);
    if (debugLevel_ > 0) {
      mf::LogInfo("MuonCaptureNormMC") << "subrun " << sr.id() << ": " << nMicrobunches_
                                       << " microbunches, POT " << sumPOT_ << ", nStop " << sumStops_
                                       << ", nCap " << nCap << " (" << provenance_ << ")";
    }
    sr.put(std::move(norm), art::fullSubRun());
  }

}

DEFINE_ART_MODULE(mu2e::MuonCaptureNormMC)
