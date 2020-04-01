// Tests to validate the CosmicBTrk Algorithms
//
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
	TrackDirection = End-Start;
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
/*
  TCanvas* hcan = new TCanvas("hcan","Cosmics",1000,1000);
//TPolyLine to graph the result
  TPolyLine3D* hel = new TPolyLine3D(nsteps+1);
  TPolyLine3D* invhel = new TPolyLine3D(nsteps+1);

  for(int istep=0;istep<nsteps+1;++istep){
  
    hel->SetPoint(istep, hpos._x, hpos._y, hpos._z);
    HelixPos(invpars,hpos);
    invhel->SetPoint(istep, hpos._x, hpos._y, hpos._z);
  }
  // draw the helix
  if(charge > 0)
    hel->SetLineColor(kBlue);
  else
    hel->SetLineColor(kRed);
  hel->Draw();

  // draw the origin and axes
  TAxis3D* rulers = new TAxis3D();
  rulers->GetXaxis()->SetAxisColor(kBlue);
  rulers->GetXaxis()->SetLabelColor(kBlue);
  rulers->GetYaxis()->SetAxisColor(kCyan);
  rulers->GetYaxis()->SetLabelColor(kCyan);
  rulers->GetZaxis()->SetAxisColor(kOrange);
  rulers->GetZaxis()->SetLabelColor(kOrange);
  rulers->Draw();

  TPolyLine3D* refmom = new TPolyLine3D(2);
  refmom->SetPoint(0,pos._x,pos._y,pos._z);
  double mommag = sqrt(mom._x*mom._x + mom._y*mom._y + mom._z*mom._z);
  double momscale = fabs(pars._pars[_rad])/mommag;
  refmom->SetPoint(1,pos._x + mom._x*momscale,pos._y + mom._y*momscale ,pos._z + mom._z*momscale);
  int momcolor;
  if(charge>0.0)
    momcolor = kRed;
  else
    momcolor = kBlack;
  refmom->SetLineColor(momcolor);
  refmom->Draw();

  TPolyMarker3D* refp = new TPolyMarker3D(1,24);
  refp->SetMarkerColor(momcolor);
  refp->SetPoint(0,pos._x,pos._y,pos._z);
  refp->Draw();
  TPolyMarker3D* refmomp = new TPolyMarker3D(1,22);
  refmomp->SetPoint(0,pos._x + mom._x*momscale,pos._y + mom._y*momscale ,pos._z + mom._z*momscale);
  refmomp->SetMarkerColor(momcolor);
  refmomp->Draw();

  TPolyMarker3D* helixp = new TPolyMarker3D(1,3);
  helixp->SetMarkerColor(kGreen);
  hpos._t = pos._t;
  HelixPos(pars,hpos);
  helixp->SetPoint(0,hpos._x,hpos._y,hpos._z);
  helixp->Draw();

  TPolyMarker3D* testp = new TPolyMarker3D(1,25);
  testp->SetMarkerColor(kOrange);
  testp->SetPoint(0,hpos._x,hpos._y,hpos._z+2*M_PI*pars._pars[_lambda]);
  testp->Draw();

  TPolyMarker3D* startp = new TPolyMarker3D(1,21);
  startp->SetMarkerColor(kBlue);
  hpos._t = tmin;
  HelixPos(pars,hpos);
  startp->SetPoint(0,hpos._x,hpos._y,hpos._z);
  startp->Draw();

  TPolyMarker3D* endp = new TPolyMarker3D(1,22);
  endp->SetMarkerColor(kBlue);
  hpos._t = tmax;
  HelixPos(pars,hpos);
  endp->SetPoint(0,hpos._x,hpos._y,hpos._z);
  endp->Draw();

  TLegend* leg = new TLegend(0.8,0.8,1.0,1.0);
  char title[80];
  snprintf(title,80,"Helix, mass=%3.1g MeV/c^{2}, q=%1.1f",mass,charge);
  leg->AddEntry(hel,title,"L");
  snprintf(title,80,"Initial Momentum =%3.1g MeV/c",momval);
  leg->AddEntry(refmom,title,"L");
  leg->AddEntry(refp,"Initial Position","P");
  snprintf(title,80,"Helix, t=%4.2g ns",pos._t);
  leg->AddEntry(helixp,title,"P");
  snprintf(title,80,"Helix, t=%4.2g ns",pos._t+tmin);
  leg->AddEntry(startp,title,"P");
  snprintf(title,80,"Helix, t=%4.2g ns",pos._t+tmax);
  leg->AddEntry(endp,title,"P");
  leg->Draw();
*/
}
