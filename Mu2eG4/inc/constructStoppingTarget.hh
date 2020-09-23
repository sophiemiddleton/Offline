#ifndef Mu2eG4_constructStoppingTarget_hh
#define Mu2eG4_constructStoppingTarget_hh
//
// Free function to construct the stopping targets.
//
//
// Original author Peter Shanahan
//
// Notes:

#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e{

    class SimpleConfig;
    // Initializor
    VolumeInfo constructStoppingTarget( VolumeInfo   const& mother,
                                      SimpleConfig const& config );
    //Foils
    VolumeInfo constructStoppingTarget_foil( VolumeInfo   const& mother,
                                      SimpleConfig const& config );
    //Screen:
    VolumeInfo constructStoppingTarget_screen( VolumeInfo   const& mother,
                                      SimpleConfig const& config );
    //Cylinder:
    VolumeInfo constructStoppingTarget_cylinder( VolumeInfo   const& mother,
                                      SimpleConfig const& config );
    //Hexagon:
    VolumeInfo constructStoppingTarget_hexagon( VolumeInfo   const& mother,
                                      SimpleConfig const& config );
    
}  // end namespace mu2e

#endif /* Mu2eG4_constructStoppingTarget_hh */
