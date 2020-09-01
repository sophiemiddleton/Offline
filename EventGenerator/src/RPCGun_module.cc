// Andrei Gaponenko, 2013

// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

// cetlib includes
#include "cetlib_except/exception.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/PionCaptureSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "MCDataProducts/inc/FixedTimeMap.hh"
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "DataProducts/inc/EventWindowMarker.hh"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"

namespace mu2e {

  //================================================================
  class RPCGun : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    BinnedSpectrum spectrum_;

    int                 verbosityLevel_;
    GenId               genId_;
    bool                applySurvivalProbability_;
    double              survivalProbScaling_;
    double              tmin_;

    double              czmin_;
    double              czmax_;
    double              phimin_;
    double              phimax_;

    art::RandomNumberGenerator::base_engine_t& eng_;

    CLHEP::RandGeneral  randSpectrum_;
    CLHEP::RandFlat     randomFlat_;
    RandomUnitSphere    randomUnitSphere_;
    PionCaptureSpectrum pionCaptureSpectrum_;

    RootTreeSampler<IO::StoppedParticleTauNormF> stops_;
    std::unique_ptr<ProtonPulseRandPDF>  protonPulse_;

    bool                doHistograms_;
    fhicl::ParameterSet protonPset_;

    double generateEnergy();

    struct Hist_t {
      TH1F* momentum;
      TH1F* eleMom {nullptr};
      TH1F* posMom {nullptr};
      TH1F* mee;
      TH2F* meeVsE;
      TH1F* meeOverE;                   // M(ee)/E(gamma)
      TH1F* y;                          // splitting function
      TH1F* ptime;                      // proton time
      TH1F* vtime;                      // vertex time
      TH1F* vtime_w;                    // vertex time, survival prob-weighted
    } _hist ;

  public:
    explicit RPCGun(const fhicl::ParameterSet& pset);
    virtual void beginRun(art::Run&   r) override;
    virtual void produce(art::Event& event);
  };

  //================================================================
  RPCGun::RPCGun(const fhicl::ParameterSet& pset)
    : art::EDProducer{pset}
    , psphys_(pset.get<fhicl::ParameterSet>("physics"))
    , spectrum_                  (BinnedSpectrum(psphys_))
    , verbosityLevel_            (pset.get<int>   ("verbosityLevel"))
    , genId_  (GenId::findByName (psphys_.get<std::string>("genId"          )))
    , applySurvivalProbability_  (psphys_.get<bool>("ApplySurvivalProb",false))
    , survivalProbScaling_       (psphys_.get<double>("SurvivalProbScaling",1))
    , tmin_                      (pset.get<double>("tmin"  ))
    , czmin_                     (pset.get<double>("czmin" ))
    , czmax_                     (pset.get<double>("czmax" ))
    , phimin_                    (pset.get<double>("phimin"))
    , phimax_                    (pset.get<double>("phimax"))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randSpectrum_       (eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , randomFlat_         (eng_)
    , randomUnitSphere_   (eng_, czmin_,czmax_,phimin_,phimax_)
    , pionCaptureSpectrum_(&randomFlat_,&randomUnitSphere_)
    //    , randomUnitSphere_(eng_)
    , stops_              (eng_,   pset.get<fhicl::ParameterSet>("pionStops"))
    , doHistograms_       (pset.get<bool>               ("doHistograms"))
    , protonPset_         (pset.get<fhicl::ParameterSet>("randPDFparameters"))
    {
      produces<mu2e::GenParticleCollection>();
      produces<mu2e::EventWeight>();
      produces<mu2e::FixedTimeMap>();

      if (phimax_ > CLHEP::twopi) phimax_ = CLHEP::twopi;
//-----------------------------------------------------------------------------
// make sure the requested genID makes sense
//-----------------------------------------------------------------------------
      if ((genId_ != GenId::InternalRPC) && (genId_ != GenId::ExternalRPC)) {
	throw cet::exception("BADCONFIG") << "RPCGun: wrong process ID: " 
					  << genId_.name() << "BAIL OUT!\n";
      }

      if(verbosityLevel_ > 0) {
        std::cout<<"RPCGun: using = "
                 <<stops_.numRecords()
                 <<" stopped particles"
                 <<std::endl;

        std::cout<<"RPCGun: producing photon " << std::endl;
      }

      if ( doHistograms_ ) {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tfdir = tfs->mkdir( "RPCGun" );

        _hist.momentum = tfdir.make<TH1F>("hmomentum"  , "photon momentum"      , 100,   40.,  140.);
	_hist.ptime    = tfdir.make<TH1F>("ptime"      , "proton time"          , 250, -500., 2000.);
	_hist.vtime    = tfdir.make<TH1F>("vtime"      , "vertex time"          , 200,    0., 2000.);
	_hist.vtime_w  = tfdir.make<TH1F>("vtime_w"    , "vertex time, weighted", 200,    0., 2000.);

        if (genId_ == GenId::InternalRPC) {
          _hist.eleMom   = tfdir.make<TH1F>("hElecMom" , "electron momentum"    , 140,  0. , 140.);
          _hist.posMom   = tfdir.make<TH1F>("hPosiMom" , "positron momentum"    , 140,  0. , 140.);
          _hist.mee      = tfdir.make<TH1F>("hMee"     , "M(e+e-) "             , 200,0.,200.);
          _hist.meeVsE   = tfdir.make<TH2F>("hMeeVsE"  , "M(e+e-) vs E"         , 200,0.,200.,200,0,200);
          _hist.meeOverE = tfdir.make<TH1F>("hMeeOverE", "M(e+e-)/E "           , 200, 0.,1);
          _hist.y        = tfdir.make<TH1F>("hy"       , "y = (ee-ep)/|pe+pp|"  , 200,-1.,1.);
        }
      }
    }

  //================================================================

  void RPCGun::beginRun(art::Run& run) {
    protonPulse_.reset( new ProtonPulseRandPDF( eng_, protonPset_ ) );
  }

  void RPCGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    std::unique_ptr<FixedTimeMap> timemap(new FixedTimeMap);

    ConditionsHandle<AcceleratorParams> accPar("ignored");
    double _mbtime = accPar->deBuncherPeriod;
    double vtime(0), ptime(0);

    IO::StoppedParticleTauNormF stop;
    if (tmin_ > 0) {
      while (true) {
        const auto& tstop = stops_.fire();
	ptime = protonPulse_->fire();
        timemap->SetTime(ptime);
	vtime = ptime+tstop.t;
        if (vtime < 0 || vtime > tmin_) {
          if (applySurvivalProbability_){
            double weight = exp(-tstop.tauNormalized)*survivalProbScaling_;
            if (weight > 1)
              std::cout << "WEIGHT TOO HIGH " << weight << " " << fmod(tstop.t,_mbtime) << std::endl;
            double rand = randomFlat_.fire();
            if (weight > rand){
              stop = tstop;
              break;
            }
          }else{
            stop = tstop;
            break;
          }
        }
      }
    }
    else {
      ptime = protonPulse_->fire();
      timemap->SetTime(ptime);
      if (applySurvivalProbability_){
        while (true){
          const auto& tstop = stops_.fire();
          double weight = exp(-tstop.tauNormalized)*survivalProbScaling_;
          double rand = randomFlat_.fire();
          if (weight > rand){
            stop = tstop;
            break;
          }
        }
      } 
      else {
        stop = stops_.fire();
      }
      vtime = ptime+stop.t;
    }

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    const double energy = generateEnergy();

    if (genId_ == GenId::ExternalRPC) {
//-----------------------------------------------------------------------------
// generate photon - external conversions
//-----------------------------------------------------------------------------
      output->emplace_back( PDGCode::gamma,
                            GenId::ExternalRPC,
                            pos,
                            CLHEP::HepLorentzVector( randomUnitSphere_.fire(energy), energy),
                            stop.t );

      event.put(std::move(output));
      event.put(std::move(timemap));
    } else {
//-----------------------------------------------------------------------------
// generate internal conversions
//-----------------------------------------------------------------------------
      CLHEP::HepLorentzVector mome, momp;
      pionCaptureSpectrum_.getElecPosiVectors(energy,mome,momp);
      // Add particles to list
      auto output = std::make_unique<GenParticleCollection>();
      output->emplace_back(PDGCode::e_minus, GenId::InternalRPC,pos,mome,stop.t);
      output->emplace_back(PDGCode::e_plus , GenId::InternalRPC,pos,momp,stop.t);
      event.put(move(output));
      event.put(std::move(timemap));

      if(doHistograms_){
        _hist.eleMom ->Fill(mome.vect().mag());
        _hist.posMom ->Fill(momp.vect().mag());

        double mee = (mome+momp).m();
        _hist.mee->Fill(mee);
        _hist.meeVsE->Fill(energy,mee);
        _hist.meeOverE->Fill(mee/energy);

        CLHEP::Hep3Vector p = mome.vect()+momp.vect();
        double y = (mome.e()-momp.e())/p.mag();

        _hist.y->Fill(y);
      }
    }

    // Calculate survival probability
    const double weight = exp(-stop.tauNormalized);
    std::unique_ptr<EventWeight> pw(new EventWeight(weight));
    event.put(std::move(pw));

    if ( !doHistograms_ ) return;

    _hist.momentum->Fill(energy);
    _hist.ptime->Fill(ptime);
    _hist.vtime->Fill(vtime);
    _hist.vtime_w->Fill(vtime,weight);

  }

  //================================================================
  double RPCGun::generateEnergy() {
    return spectrum_.sample(randSpectrum_.fire());
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::RPCGun);
