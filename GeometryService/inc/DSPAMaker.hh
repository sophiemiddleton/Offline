#ifndef GeometryService_DSPAMaker_hh
#define GeometryService_DSPAMaker_hh

#include <memory>

namespace mu2e  { class SimpleConfig; }
namespace mu2e  { class DSPA; }
namespace mu2e  { class DetectorSolenoid; }

namespace mu2e {
  class DSPAMaker {
  public:
    static std::unique_ptr<DSPA> make(const SimpleConfig& config, const DetectorSolenoid& ds);
  };
}

#endif/* GeometryService_DSPAMaker_hh */
