//Author: S Middleton
//Date: Feb 2020
//Purpose: For carrying out functions related to the GeomServices in TEve event display
#include <TObject.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
//Geom:
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/WorldG4Maker.hh"
#include "GeometryService/inc/TrackerMaker.hh"
#include "GeometryService/inc/Mu2eHallMaker.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
//Mu2e Tracker Geom:
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/Mu2eCoordTransform.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"

#include "TEveEventDisplay/src/dict_classes/Geom_Interface.h"
using namespace mu2e;
namespace mu2e{

	Geom_Interface::Geom_Interface(){}

	CLHEP::Hep3Vector Geom_Interface::GetTrackerCenter(){
                std::string filename("Mu2eG4/geom/geom_common_current.txt");
		SimpleConfig GeomConfig(filename);
		double zCenter  =  GeomConfig.getDouble("mu2e.detectorSystemZ0")*CLHEP::mm;
		double xCenter  = -GeomConfig.getDouble("mu2e.solenoidOffset")*CLHEP::mm;
		CLHEP::Hep3Vector c(xCenter, 0, zCenter);//y=0 as unchanged
		TrackMu2eOrigin = c;
		return c;
	}

	CLHEP::Hep3Vector Geom_Interface::GetCaloCenter(unsigned nDisk){
                std::string filename("Mu2eG4/geom/Calorimeter_CsI.txt");
		SimpleConfig GeomConfig(filename);
		double zCenter;
		if(nDisk==0) zCenter = GeomConfig.getDouble("mu2e.calorimeter.caloMotherZ0")*CLHEP::mm;
		if(nDisk==1) zCenter = GeomConfig.getDouble("mu2e.calorimeter.caloMotherZ1")*CLHEP::mm;
		double xCenter  = -GeomConfig.getDouble("mu2e.solenoidOffset")*CLHEP::mm;
		CLHEP::Hep3Vector c(xCenter, 0, zCenter);//y=0 as unchanged
		TrackMu2eOrigin = c;
		return c;
	}

	CLHEP::Hep3Vector Geom_Interface::GetGDMLOffsetFromMu2e(){
		std::string filename("Mu2eG4/geom/mu2eHall.txt");
		SimpleConfig HallConfig(filename);
      		double zCenter = 0;//-7929;
		double yCenter  =  0;//HallConfig.getDouble("yOfFloorSurface.below.mu2eOrigin")*CLHEP::mm;
		double xCenter  = 1288;//3040;//-HallConfig.getDouble("mu2e.solenoidOffset")*CLHEP::mm;
	        CLHEP::Hep3Vector center(xCenter,yCenter,zCenter);
		TrackerG4Origin = center;
		return center;
	}

	CLHEP::Hep3Vector Geom_Interface::GetGDMLTrackerCenter() {
		//World->HallAir->Ds3->TrackerMother ->(0,0,0)->(0,0,0)->(0,0,1288)->(0,0,0) 
		std::string filename("Mu2eG4/geom/mu2eHall.txt");
		SimpleConfig HallConfig(filename);
		double yCenter  = HallConfig.getDouble("yOfFloorSurface.below.mu2eOrigin")*CLHEP::mm;
		double zCenter  = 1288;
		double xCenter  = 0;
		CLHEP::Hep3Vector c(xCenter, yCenter, zCenter);
	     
		return c;
	}

	CLHEP::Hep3Vector Geom_Interface::PointToTracker(CLHEP::Hep3Vector point){
		CLHEP::Hep3Vector Mu2eTrackerOrigin = Geom_Interface::GetTrackerCenter();
		CLHEP::Hep3Vector PointToTracker(point.x() + Mu2eTrackerOrigin.x(), point.y()+Mu2eTrackerOrigin.y(), point.z() +Mu2eTrackerOrigin.z());
		return PointToTracker;

        }


	CLHEP::Hep3Vector Geom_Interface::PointToGDML(CLHEP::Hep3Vector point){
		//Step 1: Get Tracker origin ~ (-3904,0,10171):
		CLHEP::Hep3Vector Mu2eTrackerOrigin = Geom_Interface::GetTrackerCenter();
		//Step 2: Extract origin from GDML file:
		CLHEP::Hep3Vector GDMLTrackerOrigin = GetGDMLOffsetFromMu2e();
		//Step 3: Transform Tracker origin to the GDML origin. Need to find the transform of tracker ofigin -->GDML origin:
		//CLHEP::Hep3Vector Tracker2GDMLOrigin(GDMLTrackerOrigin + Mu2eTrackerOrigin);
		//Step 4: Transfrom a point first to tracker, then to GDML
		CLHEP::Hep3Vector PointToTracker(point.x() + Mu2eTrackerOrigin.x(), point.y()+Mu2eTrackerOrigin.y(), point.z() +Mu2eTrackerOrigin.z());
		cout<<" in tracker "<<PointToTracker.x()<<" "<<PointToTracker.y()<<" "<<PointToTracker.z()<<endl;
		//CLHEP::Hep3Vector PointToGDML(point.x()+GDMLTrackerOrigin.x(),point.y() +GDMLTrackerOrigin.y(),point.z()+GDMLTrackerOrigin.z());
		CLHEP::Hep3Vector PointToGDML(GDMLTrackerOrigin.x() - point.x(), +GDMLTrackerOrigin.y() - point.y(), GDMLTrackerOrigin.z()-point.z());
		GetGDMLTrackerCenter();

		return PointToGDML;
	}

	CLHEP::Hep3Vector Geom_Interface::TransformToG4(CLHEP::Hep3Vector point){
		/*std::string hallconfigFile("Mu2eG4/geom/mu2eHall.txt");
		SimpleConfig HallConfig(hallconfigFile);
		Mu2eHallMaker hallMaker;
		G4GeometryOptions g4opt(HallConfig);

		std::unique_ptr<Mu2eHall> Hall = hallMaker.makeBuilding(g4opt,HallConfig);
		std::string worldconfigFile("Mu2eG4/geom/geom_common_current.txt");
		SimpleConfig WorldConfig(worldconfigFile);

		WorldG4Maker worldMaker;
                std::unique_ptr<WorldG4> world = worldMaker.make(Hall, WorldConfig);
		CLHEP::Hep3Vector vec2 = world->hallFormalCenterInWorld();
		cout<<vec2.x()<<" "<<vec2.y()<<" "<<vec2.z()<<endl;*/

		return point;

	}

	CLHEP::Hep3Vector Geom_Interface::TransformToDet(CLHEP::Hep3Vector vec1){
		GeomHandle<DetectorSystem> det;
		CLHEP::Hep3Vector vec2 = det->toDetector(vec1);
		return vec2;
    	}


}
