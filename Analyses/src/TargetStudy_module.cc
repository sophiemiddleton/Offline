
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "art/Framework/Principal/Provenance.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectory.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

namespace mu2e {

  class TargetStudy : public art::EDAnalyzer {
  public:

    typedef SimParticleCollection::key_type key_type;

    explicit TargetStudy(fhicl::ParameterSet const& pset) :
      art::EDAnalyzer(pset),
      _nAnalyzed(0),
      _ntpssp(0),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _TrajModuleLabel(pset.get<std::string>("TrajModuleLabel", "g4run")),
      _stepModuleLabel(pset.get<std::string>("stepModuleLabel", "virtualdetector"))
    {
    }

    virtual ~TargetStudy() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);

    void analyze(const art::Event& e);

  private:
    int _nAnalyzed;

    TNtuple* _ntpssp;
    TNtuple* _ntptraj;
    TNtuple* _ntpstep;
    // Module label of the g4 module that produced the particles
    std::string _g4ModuleLabel;
    std::string _TrajModuleLabel;
    std::string _stepModuleLabel;
  };

  void TargetStudy::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

          // Fill the ntuple.
    _ntpssp = tfs->make<TNtuple>( "ntpssp", "SimParticle ntuple",
                                  "run:evt:trk:gid:pdg:pk:ppdg:"
                                  "partGlobalTime:partStartProperTime:partStartGlobalIndex:partStartVolIndex:partStartCreationCode:partStartOriginCC:partStartPosx:partStartPosy:partStartPosz:partStartMomx:partStartMomy:partStartMomz:"
                                  "partEndGlobalTime:partEndProperTime:partEndVolIndex:partEndG4Status:partStoppingCode:partEndPosx:partEndPosy:partEndPosz:partEndMomx:partEndMomy:partEndMomz:"
                                  "partEndKE:partEndnstep:trackLength:startEnergy:trackerEnergy"
                                  );
  _ntptraj = tfs->make<TNtuple>( "ntptraj", "Traj ntuple",
                                  "trajStartPosx:trajStartPosy:trajStartPosz:trajEndPosx:trajEndPosy:trajEndPosz:TrackLength:TotalLength:TrackAngle"
                                  );
                                  
  _ntpstep = tfs->make<TNtuple>( "ntpstep", "StepPoint ntuple",
                                  "EndMomx:EndMomy:EndMomz:EndEnergy"
                                  );
  }
  

  void TargetStudy::beginRun(art::Run const& run){

  }

  void TargetStudy::analyze(const art::Event& event) {
        bool first = true;
        bool inVolume = false;
        ++_nAnalyzed;

        // ntuple buffer.
        float nt[_ntpssp->GetNvar()];
        float nyt[_ntptraj->GetNvar()];
        float nstep[_ntpstep->GetNvar()];
        art::Handle<SimParticleCollection> simPCH;
        event.getByLabel(_g4ModuleLabel, simPCH);
        
        art::Handle<StepPointMCCollection> steps;
        event.getByLabel(_stepModuleLabel, steps);
        
        double tot = 0;
        double angle = 0;
        CLHEP::Hep3Vector StartHitPos, EndHitPos;
        art::Handle<MCTrajectoryCollection> trajcol;
        event.getByLabel(_TrajModuleLabel, trajcol);
        std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
        for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++)
        {
                tot = 0;
                const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();
                string pdgId= to_string(trajectoryIter->first->pdgId());
                for(unsigned int p = 0; p < points.size(); p++){
                        tot += sqrt(points[p].x()*points[p].x() + points[p].y()*points[p].y() + points[p].z()*points[p].z());
                }
                StartHitPos.setX(points[0].x());
                StartHitPos.setY(points[0].y());
                StartHitPos.setZ(points[0].z());
                EndHitPos.setX(points[points.size()-1].x());
                EndHitPos.setY(points[points.size()-1].y());
                EndHitPos.setZ(points[points.size()-1].z());
                angle = atan((points[points.size()-1].x() - points[0].x())/(points[points.size()-1].z() - points[0].z()));
                
        }
        nyt[0] = StartHitPos.x();
        nyt[1] = StartHitPos.y();
        nyt[2] = StartHitPos.z();
        nyt[3] = EndHitPos.x();
        nyt[4] = EndHitPos.y();
        nyt[5] = EndHitPos.z();
        
        double length = sqrt((EndHitPos.z() - StartHitPos.z())*(EndHitPos.z() - StartHitPos.z()) + (EndHitPos.y() - StartHitPos.y())*(EndHitPos.y() - StartHitPos.y()) +(EndHitPos.z() - StartHitPos.z())*(EndHitPos.z() - StartHitPos.z()));
        nyt[6] = length;
        nyt[7] = tot;
       
        nyt[8] = angle;
        _ntptraj->Fill(nyt);
        const SimParticleCollection& simPC = *simPCH;

        for (const auto& simPMVO : simPC) {

                const mu2e::SimParticle& simP = simPMVO.second;

                long unsigned int pkey = 0;
                int ppdg = 0;
                long unsigned int gid(0);

                art::Ptr<SimParticle> const& pptr = simP.parent();
                art::Ptr<GenParticle> const& gptr = simP.genParticle();

                if(pptr) {
                pkey = pptr->id().asUint();
                ppdg = pptr->pdgId();
                }

                if(gptr) {
                gid  = gptr->generatorId().id();
                ppdg = gptr->pdgId();
                }
                double dx = simP.endPosition().x() - simP.startPosition().x();
                double dy = simP.endPosition().y() - simP.startPosition().y();
                double dz = simP.endPosition().z() - simP.startPosition().z();
                double trackLength = sqrt(dx*dx + dy*dy + dz*dz);
                nt[ 0] = event.id().run(); //run
                nt[ 1] = event.id().event(); //event
                nt[ 2] = simP.id().asUint(); //track
                nt[ 3] = gid; //genId
                nt[ 4] = simP.pdgId(); //pdg
                nt[ 5] = pkey; //parent id
                nt[ 6] = ppdg; //partend pdg
                nt[ 7] = simP.startGlobalTime();
                nt[ 8] = simP.startProperTime();
                nt[ 9] = simP.startVolumeIndex();
                nt[10] = simP.startG4Status();
                nt[11] = simP.creationCode().id();
                nt[12] = simP.originParticle().creationCode().id();
                nt[13] = simP.startPosition().x();
                nt[14] = simP.startPosition().y();
                nt[15] = simP.startPosition().z();
                nt[16] = simP.startMomentum().x();
                nt[17] = simP.startMomentum().y();
                nt[18] = simP.startMomentum().z();
                nt[19] = simP.endGlobalTime();
                nt[20] = simP.endProperTime();
                nt[21] = simP.endVolumeIndex();
                nt[22] = simP.endG4Status();
                nt[23] = simP.stoppingCode().id();
                nt[24] = simP.endPosition().x();
                nt[25] = simP.endPosition().y();
                nt[26] = simP.endPosition().z();
                nt[27] = simP.endMomentum().x();
                nt[28] = simP.endMomentum().y();
                nt[29] = simP.endMomentum().z();
                nt[30] = simP.endKineticEnergy();
                nt[31] = simP.nSteps();
                nt[32] = trackLength;
                if(first){
                        nt[33] = sqrt(simP.startMomentum().x()*simP.startMomentum().x() + simP.startMomentum().y()*simP.startMomentum().y() +simP.startMomentum().z()*simP.startMomentum().z() + 0.511*0.511);
                        first = false;
                }
               
                _ntpssp->Fill(nt);

        } // end loop over simparticles
         const StepPointMCCollection& stepPC = *steps;
        
        for (const auto& step : stepPC) {
                if(step.position().z()>8400 and !inVolume){
                        std::cout<<"Found "<<event.id().event()<<std::endl;
                        nstep[0] = step.momentum().x();
                        nstep[1] = step.momentum().y();
                        nstep[2] = step.momentum().z();
                        nstep[3] = sqrt(step.momentum().x()*step.momentum().x() + step.momentum().y()*step.momentum().y() +step.momentum().z()*step.momentum().z() + 0.511*0.511);
                        inVolume = true;
                        _ntpstep->Fill(nstep);
                
                }
                
        }
  }

}  // end namespace mu2e

using mu2e::TargetStudy;
DEFINE_ART_MODULE(TargetStudy);
