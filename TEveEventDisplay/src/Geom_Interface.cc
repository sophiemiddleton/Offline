//Author: S Middleton
//Date: Feb 2020
//Purpose: For carrying out functions related to the GeomServices and GDML access in TEve event display.


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

    void Geom_Interface::InsideDS( TGeoNode * node, bool inDSVac ){
	    std::string _name = (node->GetVolume()->GetName());
	    if ( node->GetMotherVolume() ) {
		    std::string motherName(node->GetMotherVolume()->GetName());
		    if ( motherName == "DS2Vacuum" || motherName == "DS3Vacuum" ){
		    inDSVac = true;
		    }
	    }
	    if ( inDSVac && _name.find("VirtualDetector_TT_Mid") != 0 ) {
		    node->SetVisibility(kTRUE);
	    } else{
		    node->SetVisibility(kFALSE);
	    }
	    int ndau = node->GetNdaughters();
	    for ( int i=0; i<ndau; ++i ){
		    TGeoNode * dau = node->GetDaughter(i);
		    InsideDS( dau, inDSVac );
      	}

    }

	void Geom_Interface::hideNodesByMaterial(TGeoNode* node, 
					 const std::string& mat, bool onOff) {

        std::string material(node->GetVolume()->GetMaterial()->GetName());
        if ( material.find(mat) != std::string::npos ) node->SetVisibility(onOff);
        int ndau = node->GetNdaughters();
        for ( int i=0; i<ndau; ++i ){
	        TGeoNode * dau = node->GetDaughter(i);
	        hideNodesByMaterial( dau, mat, onOff);
        }

	}

	void Geom_Interface::hideNodesByName(TGeoNode* node, const std::string& str,
				     bool onOff, int _diagLevel) {

		std::string name(node->GetName());
		if ( name.find(str) != std::string::npos ){
			node->SetVisibility(onOff);
			if(_diagLevel > 0) std::cout <<"hiding "<< name << std::endl;
		}
		int ndau = node->GetNdaughters();
		for ( int i=0; i<ndau; ++i ){
			TGeoNode * dau = node->GetDaughter(i);
			hideNodesByName( dau, str, onOff, _diagLevel);
		}

	}


	void Geom_Interface::hideBuilding(TGeoNode* node) {

		static std::vector <std::string> substrings  { "Ceiling",
		"backfill", "dirt", "concrete", "VirtualDetector",
		"pipeType","CRSAluminium","CRV","CRS", "ExtShield", "PSShield"};
		for(auto& i: substrings) hideNodesByName(node,i,kFALSE, 0);

		static std::vector <std::string> materials { "MBOverburden", "CONCRETE"};
		for(auto& i: materials) hideNodesByMaterial(node,i,kFALSE);

	}

	void Geom_Interface::hideTop(TGeoNode* node, int _diagLevel) {
	  	TString name = node->GetName();
	  	if(_diagLevel > 0 and name.Index("Shield")>0) {
			std::cout << name << " " <<  name.Index("mBox_") << std::endl;
	  	}
	  	bool test = false;

		// from gg1
		if(name.Index("mBox_45_")>=0) test = true;
		if(name.Index("mBox_46_")>=0) test = true;
		if(name.Index("mBox_47_")>=0) test = true;
		if(name.Index("mBox_48_")>=0) test = true;
		if(name.Index("mBox_49_")>=0) test = true;
		if(name.Index("mBox_74_")>=0) test = true;

		if(test) {
			std::cout << "turning off " << name << std::endl;
			node->SetVisibility(false);
		}

		// Descend recursively into each daughter TGeoNode.
		int ndau = node->GetNdaughters();
		for ( int i=0; i<ndau; ++i ){
			TGeoNode * dau = node->GetDaughter(i);
			hideTop( dau, _diagLevel );
  		}

	}


    void Geom_Interface::Heirarchy( TGeoNode * node, std::vector<CLHEP::Hep3Vector> &TransformList ){
        std::string _name = (node->GetVolume()->GetName());
        if( _name == "HallAir") {
        cout<<"HallAir Origin IS "<<node->GetMotherVolume()->GetName();
        TGeoVolume *vol = node->GetVolume();
        TGeoBBox *shape = (TGeoBBox*)vol->GetShape();
        Double_t master[3];
        const Double_t *local = shape->GetOrigin();
        if(shape!=NULL){
	        gGeoManager->LocalToMaster(local,master);
	        CLHEP::Hep3Vector hallToworld(master[0], master[1], master[2]);
	        TransformList.push_back(hallToworld);
            }
        }
        if( _name == "DS3Vacuum") {
        cout<<"DS3 Origin IS "<<node->GetMotherVolume()->GetName();
        TGeoVolume *vol = node->GetVolume();
        TGeoBBox *shape = (TGeoBBox*)vol->GetShape();
        Double_t master[3];
        const Double_t *local = shape->GetOrigin();
        if(shape!=NULL){
	        gGeoManager->LocalToMaster(local,master);
	        CLHEP::Hep3Vector DSTohall(master[0], master[1], master[2]);
	        TransformList.push_back(DSTohall);
            }
        }
        if( _name == "TrackerMother") {
        cout<<"Tracker Origin IS "<<node->GetMotherVolume()->GetName();
        TGeoVolume *vol = node->GetVolume();
        TGeoBBox *shape = (TGeoBBox*)vol->GetShape();
        Double_t master[3];
        const Double_t *local = shape->GetOrigin();
        if(shape!=NULL){
	        gGeoManager->LocalToMaster(local,master);
	        CLHEP::Hep3Vector TrackerToDS(master[0], master[1], master[2]);
	        TransformList.push_back(TrackerToDS);
            }

        }

        // Descend into each daughter TGeoNode.
        int ndau = node->GetNdaughters();
        for ( int i=0; i<ndau; ++i ){
        TGeoNode * dau = node->GetDaughter(i);
        Heirarchy(dau, TransformList);
        }
      
    }

	CLHEP::Hep3Vector Geom_Interface::GetTrackerCenter(){
                std::string filename("Mu2eG4/geom/geom_common_current.txt");
		SimpleConfig GeomConfig(filename);
		double zCenter  =  GeomConfig.getDouble("mu2e.detectorSystemZ0")*CLHEP::mm;
		double xCenter  = -GeomConfig.getDouble("mu2e.solenoidOffset")*CLHEP::mm;
		CLHEP::Hep3Vector c(xCenter, 0, zCenter);//y=0 as unchanged
		TrackMu2eOrigin = c;
		return c;
	}

	CLHEP::Hep3Vector Geom_Interface::GetCaloCenter(int nDisk){
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
		
		return point;

	}

	CLHEP::Hep3Vector Geom_Interface::TransformToDet(CLHEP::Hep3Vector vec1){
		GeomHandle<DetectorSystem> det;
		CLHEP::Hep3Vector vec2 = det->toDetector(vec1);
		return vec2;
	}

	CLHEP::Hep3Vector Geom_Interface::PointToCalo(CLHEP::Hep3Vector point, int nDisk){
		CLHEP::Hep3Vector Mu2eCaloOrigin = Geom_Interface::GetCaloCenter(nDisk);
		CLHEP::Hep3Vector PointToCalo(point.x() + Mu2eCaloOrigin.x(), point.y()+Mu2eCaloOrigin.y(), point.z() + Mu2eCaloOrigin.z());
		return PointToCalo;

    }


}
