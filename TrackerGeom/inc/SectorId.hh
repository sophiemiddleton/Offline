#ifndef TrackerGeom_SectorId_hh
#define TrackerGeom_SectorId_hh

//
// Identifier for a sector.
//

//
// $Id: SectorId.hh,v 1.8 2012/05/14 19:20:45 brownd Exp $
// $Author: brownd $
// $Date: 2012/05/14 19:20:45 $
//
// Original author Rob Kutschke
//

#include <ostream>
#include "TrackerGeom/inc/DeviceId.hh"

namespace mu2e {

  struct SectorId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const SectorId& s );

  class SectorId{

  public:

    SectorId():
      _did(-1),
      _sector(-1){
    }

    SectorId( DeviceId device,
              int sector
              ):
      _did(device),
      _sector(sector){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const DeviceId& getDeviceId() const {
      return _did;
    }

    int getDevice() const {
      return _did;
    }

    int getSector() const {
      return _sector;
    }

    bool operator==(SectorId const& rhs) const{
      return ( _did == rhs._did && _sector == rhs._sector );
    }

    bool operator!=(SectorId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator < (SectorId const& rhs) const {
      return _did < rhs._did || (_did == rhs._did && _sector < rhs._sector);
    }

    friend std::ostream& operator<<(std::ostream& ost,
                                    const SectorId& s ){
      ost << s._did << " " << s._sector;
      return ost;
    }

  private:

    DeviceId _did;
    int      _sector;

  };

}  //namespace mu2e

#endif /* TrackerGeom_SectorId_hh */
