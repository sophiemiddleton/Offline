#ifndef Geom_Interface_h
#define Geom_Interface_h

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


class TGeoManager;
class TGeoVolume;
class TGMainFrame;

namespace mu2e{
	class Geom_Interface {
             
               public:
		  #ifndef __CINT__
		 explicit Geom_Interface();
		 virtual ~Geom_Interface(){};
		 TGeoManager *_geom;
		 TGeoManager* getGeom(TString filename) {
			TGeoManager *geom;
			geom = geom->TGeoManager::Import(filename);
			return geom;
 		}


		void CreateGeomManager();
		void RemoveComponents();
		void toForeground();
		

		CLHEP::Hep3Vector GetTrackerCenter();
		CLHEP::Hep3Vector GetCaloCenter(int nDisk);
		CLHEP::Hep3Vector GetGDMLTrackerCenter();		
		CLHEP::Hep3Vector GetGDMLOffsetFromMu2e();
		CLHEP::Hep3Vector PointToTracker(CLHEP::Hep3Vector point);
		CLHEP::Hep3Vector PointToCalo(CLHEP::Hep3Vector point, int nDsk);
		CLHEP::Hep3Vector PointToGDML(CLHEP::Hep3Vector point);
		CLHEP::Hep3Vector TransformToG4(CLHEP::Hep3Vector vec);
		CLHEP::Hep3Vector TransformToDet(CLHEP::Hep3Vector vec);

		void InsideDS( TGeoNode * node, bool inDSVac );
		void hideTop(TGeoNode* node, int _diagLevel);
		void hideNodesByName(TGeoNode* node, const std::string& str, bool onOff, int _diagLevel) ;
		void hideNodesByMaterial(TGeoNode* node, const std::string& mat, bool onOff);
		void hideBuilding(TGeoNode* node);
        void Heirarchy( TGeoNode * node, std::vector<CLHEP::Hep3Vector> &TransformList );

		art::Event  *_event;
		art::Run    *_run;
		CLHEP::Hep3Vector TrackerG4Origin;
		CLHEP::Hep3Vector CaloG4Origin;
		CLHEP::Hep3Vector TrackMu2eOrigin;
	        #endif
		ClassDef(Geom_Interface,0);

	}; //end class def

}//end namespace mu2e

#endif /*Geom_Interface.h*/
