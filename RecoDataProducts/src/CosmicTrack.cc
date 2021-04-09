//S. Middleton, Aug 2019
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include <vector>

using namespace std;

TrackParams::TrackParams(){
	A0 = 0.;
	A1 = 0.;
	B0 = 0.;
	B1 = 0.;
	T0 = 0.;
} 

TrackCov::TrackCov(){
  sigA0A1 = 0.
  sigA1A0 = 0.
  sigA0 = 0.;
  sigA1 = 0.;
  sigB0 = 0.;
  sigB0B1 = 0.;
  sigB1B0 = 0.;
  sigB1 = 0.;
} 

TrackAxes::TrackAxes(){
	_XDoublePrime.SetXYZ(0,0,0);
	_YDoublePrime.SetXYZ(0,0,0);
	_ZPrime.SetXYZ(0,0,0);

}
TrackEquation::TrackEquation(){
	Pos.SetXYZ(0,0,0);
	Dir.SetXYZ(0,0,0);
} 

TrackSeedDiag::TrackSeedDiag(){
 	FinalChiX = 0;
 	FinalChiY = 0;
 	FinalChiTot = 0;
 	
 	InitialChiX = 0;
 	InitialChiY = 0;
 	InitialChiTot = 0;
 	
	}

namespace mu2e{

	CosmicTrack::CosmicTrack() {
    		
    InitParams.A0 = 0;
    InitParams.A1 = 0;
    InitParams.B0 = 0;
    InitParams.B1 = 0;
    InitParams.T0 = 0;         	
	 }

	// Destructor
	CosmicTrack::~CosmicTrack() {}
	
	//function to make tuple of POCA info
  std::tuple <XYZVec, double, double> CosmicTrack::GetTrackPOCAInfo() {
    XYZVec const& zpos(0,0,0);
    XYZVec const& zdir(0,0,1);
    XYZVec const& pos0(this->MinuitParams.A0, 0, this->MinuitParams.B0);
    XYZVec const& dir(this->MinuitParams.A1, -1, this->MinuitParams.B1);

    std::tuple <XYZVec, double, double> poca_info;
    TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(pos0, dir, zpos, zdir);
    XYZVec POCA = PCA.pca();
    double DOCA = PCA.dca();
    double AMSIGN = copysign(1.0,PCA.pca().X());
    poca_info = make_tuple(POCA, DOCA, AMSIGN);
    return poca_info;
  }
	    

}
