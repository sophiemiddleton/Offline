#ifndef Draw_Interface_h
#define Draw_Interface_h

#include <TObject.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGString.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
#include<TGLViewer.h>
#include <TGMsgBox.h>
// ... libGeom
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TGeoNode.h>
#include <TGeoPhysicalNode.h>
// ... libRGL
#include <TGLViewer.h>
// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "ConfigTools/inc/SimpleConfig.hh"

#include "TEveEventDisplay/src/dict_classes/Geom_Interface.h"

namespace mu2e{
	class Draw_Interface {
            
    public:
		 #ifndef __CINT__
		 explicit Draw_Interface();
		 virtual ~Draw_Interface(){};
		 void DrawHit(const std::string &pstr, Int_t mColor, Int_t mSize, Int_t n, CLHEP::Hep3Vector HitPos, TEveElementList *list, Geom_Interface *g);
		void DrawCluster(const std::string &pstr, Int_t mColor, Int_t mSize, Int_t n, CLHEP::Hep3Vector ClusterPos, int nDisk, TEveElementList *list, Geom_Interface *g);

	   #endif

		 ClassDef(Draw_Interface,0);

	}; //end class def

}//end namespace mu2e

#endif /*Draw_Interface.h*/
