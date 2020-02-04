//Author: S Middleton
// Purpose: Holds functions for the kalman fitting of Cosmic Tracks, based on KalFit.cc
// Date: Dec 2019

#include "CosmicReco/inc/CosmicKalFit.hh"
#include "CosmicReco/inc/CosmicTrkUtils.hh"
#include "CosmicReco/inc/CosmicKalFitData.hh"

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/Mu2eMaterial.hh"
#include "TrackerConditions/inc/Mu2eDetector.hh"
// utiliites
#include "GeneralUtilities/inc/Angles.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
// data
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"

// BaBar
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
//#include "BTrk/TrkBase/TrkHelixUtils.hh" //TODO
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
// Mu2e BaBar
#include "BTrkData/inc/TrkStrawHit.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
// root
#include "TH1F.h"
#include "TTree.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;

namespace mu2e
{
	// comparison functor for ordering hits.  This should operate on TrkHit, FIXME!
	struct fcomp : public binary_function<TrkHit*, TrkHit*, bool> {
		bool operator()(TrkHit* x, TrkHit* y) {
		return x->fltLen() < y->fltLen();
		}
	};


	// construct from a parameter set
  	CosmicKalFit::CosmicKalFit(fhicl::ParameterSet const& pset) :
    	_debug(pset.get<int>("debugLevel",0)),
    	_maxpull(pset.get<double>("maxPull",5)),
	_strHitW(pset.get<double>("strawHitT0Weight")),
	_minnstraws(pset.get<double>("minnstraws",2)),
	_herr(pset.get< vector<double> >("hiterr")),
	_maxIterations(pset.get<unsigned>("maxIterations",10)),
	_bfield(0){
		// set KalContext parameters
		_disttol = pset.get<double>("IterationTolerance",0.1);
		_intertol = pset.get<double>("IntersectionTolerance",100.0);
		_maxiter = pset.get<long>("MaxIterations",10);
		_maxinter = pset.get<long>("MaxIntersections",0);
		_matcorr = pset.get<bool>("materialCorrection",true);
		_fieldcorr = pset.get<bool>("fieldCorrection",false);
		_smearfactor = pset.get<double>("SeedSmear",1.0e6);
		_sitethresh = pset.get<double>("SiteMomThreshold",0.2);
		_momthresh = pset.get<double>("MomThreshold",10.0);
		_mingap = pset.get<double>("mingap",1.0);
		_minfltlen = pset.get<double>("MinFltLen",0.1);
		_minmom = pset.get<double>("MinMom",10.0);
		_fltepsilon = pset.get<double>("FltEpsilon",0.001);
		_divergeflt = pset.get<double>("DivergeFlt",1.0e3);
		_mindot = pset.get<double>("MinDot",0.0);
		_maxmomdiff = pset.get<double>("MaxMomDiff",0.5);
		_momfac = pset.get<double>("MomFactor",0.0);
		_maxpardif[0] = _maxpardif[1] = pset.get<double>("MaxParameterDifference",1.0);

		_mindof = pset.get<double>("MinNDOF",10);

		_printUtils = new TrkPrintUtils(pset.get<fhicl::ParameterSet>("printUtils",fhicl::ParameterSet()));
	}

	 CosmicKalFit::~CosmicKalFit(){
	    delete _bfield;
  	}




//-----------------------------------------------------------------------------
// create the track (KalRep) from a cosmic track seed
//-----------------------------------------------------------------------------
	void CosmicKalFit::MakeTrack(StrawResponse::cptr_t srep,Mu2eDetector::cptr_t detmodel, CosmicKalFitData& kalData){

    	if(fitable(*kalData.cosmicKalSeed)){

		double flt0 = kalData.cosmicKalSeed->flt0();
		auto kseg = kalData.cosmicKalSeed->nearestSegment(flt0);
		
		HepVector pvec(4,0);
		HepSymMatrix pcov(4,0);
		kseg->cosmic().hepVector(pvec);
		kseg->cosmic_cov().symMatrix(pcov);

	      	CosmicLineTraj htraj(pvec,pcov);
	      	kalData.cosmicTraj = &htraj;

	      	TrkStrawHitVector tshv;
	      	MakeTrkStrawHits(srep,kalData, tshv);

	      	std::vector<DetIntersection> detinter;
	      	if(_matcorr)MakeMaterials(detmodel, tshv,*kalData.cosmicTraj,detinter);

	      	std::vector<TrkHit*> thv(0);
	      	for(auto ihit = tshv.begin(); ihit != tshv.end(); ++ihit){
			thv.push_back(*ihit);
	      	}

	      	TrkT0 t0(kalData.cosmicKalSeed->t0()); 
	      	kalData.krep = new KalRep(htraj, thv, detinter, *this, kalData.cosmicKalSeed->particle(), t0, flt0); 
	      	assert(kalData.krep != 0);
		if (_debug > 0) {
			char msg[100];
			sprintf(msg,"makeTrack_001 annealing step: %2i",_annealingStep);
			_printUtils->printTrack(kalData.event,kalData.krep,"banner+data+hits",msg);
      		}

	      	kalData.krep->addHistory(TrkErrCode(),"KalFit creation");
		
	      	TrkErrCode fitstat = FitTrack(detmodel,kalData);
	      	kalData.krep->addHistory(fitstat,"KalFit fit");

		}
  	}

	bool CosmicKalFit::fitable(CosmicKalSeed const& kseed){
    		return kseed.segments().size() > 0 && kseed.hits().size() >= _minnstraws;
  	}

  	void CosmicKalFit::MakeTrkStrawHits(StrawResponse::cptr_t srep, CosmicKalFitData& kalData, TrkStrawHitVector& tshv ) {

	    std::vector<TrkStrawHitSeed>const hseeds = kalData.cosmicSeed->trkstrawhits();
	    CosmicLineTraj const htraj = *kalData.cosmicTraj;

	    for(auto ths : hseeds ){
		size_t index = ths.index();
		const ComboHit& strawhit(kalData.chcol->at(index));
		const Straw& straw = _tracker->getStraw(strawhit.strawId());
		TrkStrawHit* trkhit = new TrkStrawHit(srep,strawhit,straw,ths.index(),ths.t0(),ths.trkLen(), _maxpull,_strHitW);
		assert(trkhit != 0);

		trkhit->setAmbig(ths.ambig());

		TrkErrCode pstat = trkhit->updatePoca(&htraj);
	      	if(pstat.failure()){
			trkhit->setActivity(false);
	      	}
	      	tshv.push_back(trkhit);
    	}

    	std::sort(tshv.begin(),tshv.end(),fcomp());
  	}

	void CosmicKalFit::MakeMaterials( Mu2eDetector::cptr_t detmodel,TrkStrawHitVector const& tshv, CosmicLineTraj const& traj,std::vector<DetIntersection>& detinter) {
		// loop over strawhits and extract the straws
		for (auto trkhit : tshv) {
		// find the DetElem associated this straw
		const DetStrawElem* strawelem = detmodel->strawElem(trkhit->straw());
		// create intersection object for this element; it includes all materials
		DetIntersection strawinter;
		strawinter.delem = strawelem;
		strawinter.pathlen = trkhit->fltLen();
		strawinter.thit = trkhit;
		// compute initial intersection: this gets updated each fit iteration
		strawelem->reIntersect(&traj,strawinter);
		detinter.push_back(strawinter);
    		}
  	}

	TrkErrCode CosmicKalFit::FitTrack(Mu2eDetector::cptr_t detmodel, CosmicKalFitData& kalData) {
	    	// loop over external hit errors, ambiguity assignment, t0 toleratnce
	    	TrkErrCode fitstat;
	    	for(size_t iherr=0;iherr < _herr.size(); ++iherr) {
			fitstat = FitIteration(detmodel,kalData,iherr);
			if(_debug > 0) { 
				cout << "Iteration " << iherr 
				     << " NDOF = " << kalData.krep->nDof() 
				     << " Fit Status = " <<  fitstat << endl;

				char msg[200];
				sprintf(msg,"CosmicKalFit::fitTrack Iteration = %2li success = %i",iherr,fitstat.success());
				//_printUtils->PrintTrack(kalData.event,kalData.krep,"banner+data+hits",msg);
			}
			if(!fitstat.success())break;
	    	}
	    	return fitstat;
	}

	TrkErrCode CosmicKalFit::FitIteration(Mu2eDetector::cptr_t detmodel,
				  CosmicKalFitData& kalData, int iter) {

		if (iter == -1) iter =  _herr.size()-1;
		_annealingStep = iter;//used in the printHits routine

		TrkHitVector* thv   = &(kalData.krep->hitVector());
		for (auto itsh=thv->begin();itsh!=thv->end(); ++itsh){
			(*itsh)->setTemperature(_herr[iter]);
		}

		// update t0, and propagate it to the hits
		double oldt0 = kalData.krep->t0()._t0;
		unsigned niter(0);
		bool changed(true);
		TrkErrCode retval = TrkErrCode::succeed;

		KalRep* krep =  kalData.krep;
		//bool    flagMaterialAdded(false);

		while(retval.success() && changed && ++niter < _maxIterations){

		krep->resetFit();
		retval = krep->fit();
		if(! retval.success())break;

		/* unsigned nmat(0);
		if(_addmaterial[iter]){
		nmat = addMaterial(detmodel,krep);
		changed |= nmat>0;
		if (!flagMaterialAdded) flagMaterialAdded=true;
		}*/

		if(_debug > 1) std::cout << "Inner iteration " << niter << " changed = "
		<< changed << " t0 old " << oldt0 << " new " << krep->t0()._t0 <<std::endl;
		oldt0 = krep->t0()._t0;
		}
		if(_debug > 1)
		std::cout << "Outer iteration " << iter << " stopped after "
		<< niter << " iterations" << std::endl;
		// make sure the fit is current
		if(!krep->fitCurrent())retval = krep->fit();
		return retval;
  	}

	BField const& CosmicKalFit::bField() const {
		//GeomHandle<BFieldManager> *bf;
		if(_bfield == 0){
		      GeomHandle<BFieldConfig> bfconf;
			_bfield= 0;//new BFieldFixed(bfconf->getDSUniformValue());
			
		      }
	    	
    		return *_bfield;
  	}

	

}
