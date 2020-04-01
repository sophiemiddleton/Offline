// Tests to validate the CosmicBTrk Algorithms
//Usage: ROOT - .L Test.C , TVector3 start(0,0,0), TVector3(1,2,3),
#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
#include "TPolyLine3D.h"
#include "TArrow.h"
#include "TAxis3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include "TVector3.h"
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace ROOT::Math;


TVector3 GetPOCA( TVector3 point,TVector3 start,TVector3 end){

    double tMin = -(start-point).Dot(end-start) /((end-start).Mag2());
	  double POCA_x = start.x() + (end.x()-start.x())*tMin;
	  double POCA_y = start.y() + (end.x()-start.y())*tMin;
	  double POCA_z = start.z() + (end.z()-start.z())*tMin;
	  TVector3 closestPointOnLine;
	  closestPointOnLine.SetXYZ(POCA_x ,POCA_y ,POCA_z);

	 return closestPointOnLine;

}

double GetDOCA(TVector3 closestPointOnLine, TVector3 point){
	return sqrt((closestPointOnLine-point).Mag2());

}

TVector3 direction(double phi0, double theta)
{
	TVector3 Dir;
	Dir.SetXYZ(cos(phi0)*sin(theta),sin(theta)*sin(phi0), cos(theta));
  return Dir;
}

TVector3 position(double f, double phi0, double theta, double z0, double d0, TVector3 Pos0)
{
	double sphi0 = sin(phi0);
	double cphi0 = cos(phi0);
	double x_pos = -1*d0*sphi0+ cos(phi0)*sin(theta)*f + Pos0.x() ;
	double y_pos = d0*cphi0+  sin(theta)*sin(phi0)*f + Pos0.y();
	double z_pos = z0 + cos(theta)*f + Pos0.z();
	TVector3 Pos;
	Pos.SetXYZ(x_pos, y_pos, z_pos);
	return Pos;
}

void TestBTrk( TVector3 Start, TVector3 End) {
	TVector3 Pos0;
	Pos0.SetXYZ(0,0,0); //Reference point is set at 0.
	TVector3 POCA = GetPOCA(Pos0, Start, End);
	unsigned AMSIGN = copysign(1,POCA.X());
	double DOCA = GetDOCA(POCA, Pos0);
	TVector3 TrackDirection;
	TrackDirection = End - Start;
	//TrackDirection.SetXYZ(TrackDirection.X()/TrackDirection.Z(),TrackDirection.Y()/TrackDirection.Z(),1);
	double theta =asin(TrackDirection.Y()/(TrackDirection.Mag()));
	double phi0 = AMSIGN*atan2(POCA.X(),POCA.Y());
	double z0  = POCA.Z();
	double d0 = DOCA;

	cout << "Cosmic Test Parameter:" << endl
	<< "Theta = " <<theta << endl
	<< "phi-0 = " << phi0 << endl
	<< "d0 = " << d0 << endl
	<< "z0 = " << z0 << endl;

	TVector3 Dir = direction(phi0,theta);
	TVector3 TestDirection = TrackDirection - Dir;
	cout<<"Test 1: Direction test : "<<" input direction : "<<TrackDirection.X()<<" "<<endl
	<<TrackDirection.Y()<<" "<<TrackDirection.Z()<<endl
	<<" output direction "<<Dir.X()<<" "<<Dir.Y()<<" "<<Dir.Z()<<endl;

	double f = (End-Start).Mag();
	TVector3 TestPosition = position(1,phi0,theta,z0,d0, Pos0);
	cout<<"Test 2 : given End can I get back to Start :" <<endl;
	TVector3 PosEndTest = TestPosition + f*TestDirection;
	cout<<"start position is: "<<Start.X()<<" "<<Start.Y()<<" "<<Start.Z()<<endl;
	cout<<"test position is:"<<PosEndTest.X()<<" "<<PosEndTest.Y()<<" "<<PosEndTest.Z()<<endl;

}

void derivDeflect(double fltlen, double theta, double d0, int opt)
{

  double thetaIndex;
  double d0Index;
  double phi0Index;
  double z0Index;
  if(opt == 1){

    thetaIndex = 1;
    d0Index = 0;
    phi0Index = 0;
    z0Index = -1*fltlen*(1/sin(theta));
  }
  if(opt == 2){

    thetaIndex = 0;
    d0Index = -1*fltlen;
    phi0Index = 1/sin(theta);
    z0Index = -d0*1/(sin(theta)*tan(theta));

  }
  std::cout<<"Theta "<<thetaIndex<<" d0 "<<d0Index<<endl
  <<" z0 "<<z0Index<<" phi 0 "<<phi0Index<<std::endl;
}


void calcCosmicLineAllCovs(double fltlen, double theta, double phi0, double z0, double d0) {

  double sinphi0 = sin(phi0);
  double cosphi0 = cos(phi0);

  double costheta = cos(theta);
  double sintheta = sin(theta);

  TrkMomVisitor theVisitor(theTraj);
  double pt =0;//TODO

  double d_x_d0 = -1*sinphi0;
  double d_x_phi0 = -d0*cosphi0-fltlen*sintheta*sinphi0;
  double d_x_theta = fltlen*cosphi0*costheta;
  double d_x_z0 = 0;

  double d_y_d0 = cosphi0;
  double d_y_phi0 = -1*d0*sinphi0+fltlen*sintheta*cosphi0;
  double d_y_theta = fltlen*sinphi0*costheta;
  double d_y_z0 = 0;

  double d_z_theta = -1*fltlen*sintheta;
  double d_z_z0 = 1;
  double d_z_d0 = 0;
  double d_z_phi0 = 0;

  double d_px_phi0 = -pt*sinphi0*sintheta;
  double d_py_phi0 = pt*cosphi0*sintheta;
  double d_px_theta = pt*cosphi0*costheta;
  double d_py_theta = pt*sinphi0*costheta;
  double d_pz_phi0 = 0;
  double d_pz_theta = -1*pt*sintheta;
  double d_pz_z0 = 1;

  double m_d0_d0 =  1;
  double m_d0_phi0 = 0;
  double m_d0_theta =  0;
  double m_d0_z0 = 0;

  double m_phi0_phi0 =  1;
  double m_phi0_z0 =  0;
  double m_phi0_theta =  0;

  double m_theta_theta =  1;
  double m_z0_theta =  0;
  double m_z0_z0 =  1;

  double xxCov11 =
    d_x_d0* (  d_x_d0*m_d0_d0 + d_x_phi0*m_d0_phi0 + d_x_theta*m_d0_theta + d_x_z0*m_d0_z0  )   +
    d_x_phi0* (  d_x_d0*m_d0_phi0 + d_x_phi0*m_phi0_phi0 + d_x_theta*m_phi0_theta + d_x_z0*m_phi0_z0  )   +
    d_x_theta* (  d_x_d0*m_d0_theta + d_x_phi0*m_phi0_theta + d_x_theta*m_theta_theta + d_x_z0*m_z0_theta  )   +
    d_x_z0* (  d_x_d0*m_d0_z0 + d_x_phi0*m_phi0_z0 + d_x_theta*m_phi0_theta + d_x_z0*m_z0_z0  )  ;
  double xxCov21 =
    d_y_d0* (  d_x_d0*m_d0_d0 + d_x_phi0*m_d0_phi0 + d_x_theta*m_d0_theta + d_x_z0*m_d0_z0  )   +
    d_y_phi0* (  d_x_d0*m_d0_phi0 + d_x_phi0*m_phi0_phi0 + d_x_theta*m_phi0_theta + d_x_z0*m_phi0_z0  )   +
    d_y_theta* (  d_x_d0*m_d0_theta + d_x_phi0*m_phi0_theta + d_x_theta*m_theta_theta + d_x_z0*m_z0_theta  )   +
    d_y_z0* (  d_x_d0*m_d0_z0 + d_x_phi0*m_phi0_z0 + d_x_theta*m_z0_theta + d_x_z0*m_z0_z0  )   ;
  double xxCov22 =
    d_y_d0* (  d_y_d0*m_d0_d0 + d_y_phi0*m_d0_phi0 + d_y_theta*m_d0_theta + d_y_z0*m_d0_z0  )   +
    d_y_phi0* (  d_y_d0*m_d0_phi0 + d_y_phi0*m_phi0_phi0 + d_y_theta*m_phi0_theta + d_y_z0*m_phi0_z0  )   +
    d_y_theta* (  d_y_d0*m_d0_theta + d_y_phi0*m_phi0_theta + d_y_theta*m_theta_theta + d_y_z0*m_z0_theta  )   +
    d_y_z0* (  d_y_d0*m_d0_z0 + d_y_phi0*m_phi0_z0 + d_y_theta*m_z0_theta + d_y_z0*m_z0_z0  )   ;

  double xxCov31 =
    d_z_d0*(  d_x_d0*m_d0_d0 + d_x_phi0*m_d0_theta + d_x_theta*m_d0_theta + d_x_z0*m_d0_theta )+
    d_z_phi0*(  d_x_d0*m_d0_phi0 + d_x_phi0*m_phi0_phi0 + d_x_theta*m_phi0_theta + d_x_z0*m_phi0_z0 )+d_z_theta*(  d_x_d0*m_d0_theta + d_x_phi0*m_phi0_theta + d_x_theta*m_theta_theta + d_x_z0*m_z0_theta )+
   d_z_z0*(  d_x_d0*m_d0_phi0 + d_x_phi0*m_phi0_z0 + d_x_theta*m_z0_theta + d_x_z0*m_z0_z0 );

  double xxCov32 = d_z_d0*(  d_y_d0*m_d0_d0 + d_y_phi0*m_d0_theta + d_y_theta*m_d0_theta + d_y_z0*m_d0_theta )+
    d_z_phi0*(  d_y_d0*m_d0_phi0 + d_y_phi0*m_phi0_phi0 + d_y_theta*m_phi0_theta + d_y_z0*m_phi0_z0 )+d_z_theta*(  d_y_d0*m_d0_theta + d_y_phi0*m_phi0_theta + d_y_theta*m_theta_theta + d_y_z0*m_z0_theta )+
   d_z_z0*(  d_y_d0*m_d0_phi0 + d_y_phi0*m_phi0_z0 + d_y_theta*m_z0_theta + d_y_z0*m_z0_z0 );

  double xxCov33 = d_z_d0*(  d_z_d0*m_d0_d0 + d_z_phi0*m_d0_theta + d_z_theta*m_d0_theta + d_z_z0*m_d0_theta )+
    d_z_phi0*(  d_z_d0*m_d0_phi0 + d_z_phi0*m_phi0_phi0 + d_z_theta*m_phi0_theta + d_z_z0*m_phi0_z0 )+d_z_theta*(  d_z_d0*m_d0_theta + d_z_phi0*m_phi0_theta + d_z_theta*m_theta_theta + d_z_z0*m_z0_theta )+
   d_z_z0*(  d_z_d0*m_d0_phi0 + d_z_phi0*m_phi0_z0 + d_z_theta*m_z0_theta + d_z_z0*m_z0_z0 );


  double xpCov11 = d_px_theta* (  d_x_d0*m_d0_theta + d_x_phi0*m_phi0_theta + d_x_theta*m_theta_theta + d_x_z0*m_z0_theta  )  +
 d_px_phi0* (  d_x_d0*m_d0_phi0 + d_x_phi0*m_phi0_phi0 + d_x_theta*m_phi0_theta + d_x_z0*m_phi0_z0  ) ;

  double xpCov21 = d_px_theta* (  d_y_d0*m_d0_theta + d_y_phi0*m_phi0_theta + d_y_theta*m_theta_theta + d_y_z0*m_z0_theta  )   +
    d_px_phi0* (  d_y_d0*m_d0_phi0 + d_y_phi0*m_phi0_phi0 + d_y_theta*m_phi0_theta + d_y_z0*m_phi0_z0  )   ;

  double xpCov31 = d_px_theta* (  d_z_d0*m_d0_theta + d_z_phi0*m_phi0_theta + d_z_theta*m_theta_theta + d_z_z0*m_z0_theta  )   +
    d_px_phi0* (  d_z_d0*m_d0_phi0 + d_z_phi0*m_phi0_z0 + d_z_theta*m_phi0_theta + d_z_z0*m_phi0_z0  )   ;

  double xpCov12 =  d_py_theta* (  d_x_d0*m_d0_theta + d_x_phi0*m_phi0_theta + d_x_theta*m_theta_theta + d_x_z0*m_z0_theta  )   +
    d_py_phi0* (  d_x_d0*m_d0_phi0 + d_x_phi0*m_phi0_phi0 + d_x_theta*m_phi0_theta + d_x_z0*m_phi0_z0 )   ;

  double xpCov22 =  d_py_theta* (  d_y_d0*m_d0_theta + d_y_phi0*m_phi0_theta + d_y_theta*m_theta_theta + d_y_z0*m_z0_theta  )   +
    d_py_phi0* (  d_y_d0*m_d0_phi0 + d_y_phi0*m_phi0_phi0 + d_y_theta*m_phi0_theta + d_y_z0*m_phi0_z0 )   ;

  double xpCov32 =  d_py_theta* (  d_z_d0*m_d0_theta + d_z_phi0*m_phi0_theta + d_z_theta*m_theta_theta + d_y_z0*m_z0_theta  )   +
    d_py_phi0* (  d_z_d0*m_d0_phi0 + d_z_phi0*m_phi0_phi0 + d_z_theta*m_phi0_theta + d_z_z0*m_phi0_z0  )   ;

  double xpCov13 =  d_pz_theta* (  d_x_d0*m_d0_theta + d_x_phi0*m_phi0_theta + d_x_theta*m_theta_theta + d_x_z0*m_z0_theta  )   +
    d_pz_phi0* (  d_x_d0*m_d0_phi0 + d_x_phi0*m_phi0_phi0 + d_x_theta*m_phi0_theta + d_x_z0*m_phi0_z0  )   ;

  double xpCov23 =  d_pz_theta* (  d_y_d0*m_d0_theta + d_y_phi0*m_phi0_theta + d_y_theta*m_theta_theta + d_y_z0*m_z0_theta  )   +
    d_pz_phi0* (  d_y_d0*m_d0_phi0 + d_y_phi0*m_phi0_z0 + d_y_theta*m_phi0_theta + d_y_z0*m_phi0_z0 )   ;

  double xpCov33 =  d_pz_theta* (  d_z_d0*m_d0_theta + d_z_phi0*m_phi0_theta + d_z_theta*m_theta_theta + d_z_z0*m_z0_theta  )   +
    d_pz_phi0* (  d_z_d0*m_d0_phi0 + d_z_phi0*m_phi0_phi0 + d_z_theta*m_phi0_theta + d_z_z0*m_phi0_z0)   ;


  double ppCov11 = d_px_theta* ( d_px_theta*m_theta_theta + d_px_phi0*m_phi0_theta  )   ;

  double ppCov21 = d_py_theta* ( d_px_theta*m_theta_theta + d_px_phi0*m_phi0_theta  );

  double ppCov22 = d_py_theta* ( d_py_theta*m_theta_theta + d_py_phi0*m_phi0_theta  );

  double ppCov31 = d_pz_theta* ( d_px_theta*m_theta_theta + d_px_phi0*m_phi0_theta  ) ;

  double ppCov32 = d_pz_theta* ( d_py_theta*m_theta_theta + d_py_phi0*m_phi0_theta  ) ;

  double ppCov33 = d_pz_theta* ( d_pz_theta*m_theta_theta  ) + d_pz_z0*(d_pz_z0*m_z0_z0) ;

}

