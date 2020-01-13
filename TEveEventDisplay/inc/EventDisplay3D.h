
// ROOT includes
// ... libCore
#include <TApplication.h>
#include <TString.h>
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
// ... libGeom
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
// ... libEG
#include <TParticle.h>
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
#include <TEvePointSet.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include "art/Framework/Principal/Event.h"
#include <sstream>
#include <memory>

#include "TEveEventDisplay/inc/EvtDisplayUtils.h"



  class EventDisplay3D  {

  public:

    explicit EventDisplay3D();

    void beginJob();
    //void endJob();
    void beginRun();
    void analyze(const art::Event& event);

    class point {
    public:
      float x;
      float y;
      float z;
    };
    
    class traj {
    public:
      int id;
      std::vector<point> tv;
    };

  private:

    // Set by parameter set variables.
    //art::InputTag   gensTag_;
    bool            drawGenTracks_;
    bool            drawHits_;
    Double_t        hitMarkerSize_;
    Double_t        trkMaxR_;
    Double_t        trkMaxZ_;
    Double_t        trkMaxStepSize_;
    Double_t        camRotateCenterH_;
    Double_t        camRotateCenterV_;
    Double_t        camDollyDelta_;

    //art::ServiceHandle<Geometry>          geom_;
    //art::ServiceHandle<PDT>               pdt_;

    EvtDisplayUtils* visutil_;
    TEveGeoShape* fSimpleGeom;

    TEveViewer *fXYView;
    TEveViewer *fRZView;
    TEveProjectionManager *fXYMgr;
    TEveProjectionManager *fRZMgr;
    TEveScene *fDetXYScene;
    TEveScene *fDetRZScene;
    TEveScene *fEvtXYScene;
    TEveScene *fEvtRZScene;

    TGTextEntry      *fTeRun,*fTeEvt;
    TGLabel          *fTlRun,*fTlEvt;

    TEveTrackList *fTrackList;
    TEveElementList *fHitsList;

    void makeNavPanel();

    void hideNodesByName(TGeoNode* node, const std::string& str, bool onOff);
    void hideNodesByMaterial(TGeoNode* node, const std::string& str, bool onOff);
    void hideBuilding(TGeoNode* node);
    void hideTop(TGeoNode* node);
    void hideTop2();
    void InsideDS( TGeoNode * node, bool inDSVac );
    void AddHits(const art::Event& event);
    static const int nn=1000;
    int nev;
    int ntr;
    float x[3][nn+1];
    float c[3][nn+1];
    int   f[4][nn+1];
    
    std::vector<int> vtptr;
    std::vector<traj> vtraj;
    
    bool traj_to_plane(traj& tr, int dir, float w, 
		       float& xp, float& yp);
    
    bool event_traj_to_plane(int iev, int dir, float w, 
			     float& xp, float& yp);
    
    void read_traj_file();
    
    
  };

