#include<TEvePointSet.h>

#include "TEveEventDisplay/src/dict_classes/Draw_Interface.h"

using namespace mu2e;
namespace mu2e{

	Draw_Interface::Draw_Interface(){}

	void Draw_Interface::DrawHit(const std::string &pstr, Int_t mColor, Int_t mSize, Int_t n, CLHEP::Hep3Vector HitPos, TEveElementList *list, Geom_Interface *g)
  	{
	
		g->GetTrackerCenter();
		CLHEP::Hep3Vector pointInMu2e = g->PointToTracker(HitPos);
		std::string hstr=" hit %d";
		std::string dstr=" hit# %d\nLayer: %d";
		std::string strlst=pstr+hstr;
		std::string strlab=pstr+dstr;

		TEvePointSet* h = new TEvePointSet(Form(strlst.c_str(),n));
		h->SetTitle(Form(strlab.c_str(),n,hstr));

		std::cout<<"in mu2e : "<<n<<" "<<pointInMu2e.x()/10<<" "<<pointInMu2e.y()/10<<" "<<pointInMu2e.z()/10<<std::endl;
		h->SetMarkerColor(kRed);
		h->SetMarkerSize(mSize);

		h->SetNextPoint(pointInMu2e.x()/10, pointInMu2e.y()/10, pointInMu2e.z()/10); //as GDML is in cm (I think...and hope)
		h->SetMarkerColor(mColor);
		h->SetMarkerSize(mSize);
		list->AddElement(h);
	  }


	void Draw_Interface::DrawCluster(const std::string &pstr, Int_t mColor, Int_t mSize, Int_t n, CLHEP::Hep3Vector ClusterPos, int nDisk, TEveElementList *list, Geom_Interface *g) 
  	{
		
		g->GetCaloCenter(nDisk);
		CLHEP::Hep3Vector pointInMu2e = g->PointToCalo(ClusterPos, nDisk);
		std::string hstr=" cluster %d";
		std::string dstr=" cluster# %d\nLayer: %d";
		std::string strlst=pstr+hstr;
		std::string strlab=pstr+dstr;

		TEvePointSet* h = new TEvePointSet(Form(strlst.c_str(),n));
		h->SetTitle(Form(strlab.c_str(),n,hstr));

		std::cout<<"Cluster in mu2e : "<<n<<" "<<pointInMu2e.x()/10<<" "<<pointInMu2e.y()/10<<" "<<pointInMu2e.z()/10<<std::endl;
		h->SetMarkerColor(kRed);
		h->SetMarkerSize(mSize);

		h->SetNextPoint(pointInMu2e.x()/10, pointInMu2e.y()/10, pointInMu2e.z()/10); 
		h->SetMarkerColor(mColor);
		h->SetMarkerSize(mSize);
		list->AddElement(h);
	  }


   }

