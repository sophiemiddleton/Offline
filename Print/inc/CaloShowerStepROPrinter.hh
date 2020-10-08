//
//  Utility class to print CaloShowerStepRO
// 
#ifndef Print_inc_CaloShowerStepROPrinter_hh
#define Print_inc_CaloShowerStepROPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CaloShowerStepROPrinter : public ProductPrinter {
  public:

    typedef std::vector<std::string> vecstr;

    CaloShowerStepROPrinter():_eCut(-1) { }
    CaloShowerStepROPrinter(const ConfigE& conf):ProductPrinter(conf) { 
      _eCut = conf.eCut();
    }

    // do not print if e is below this cut
    void setECut(double e) { _eCut = e; }
    double eCut() const { return _eCut; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<CaloShowerStepROCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CaloShowerStepROCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const CaloShowerStepROCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<CaloShowerStepRO>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::CaloShowerStepRO& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    double _eCut;

  };

}
#endif
