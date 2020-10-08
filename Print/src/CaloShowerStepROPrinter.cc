
#include "Print/inc/CaloShowerStepROPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CaloShowerStepROPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CaloShowerStepROCollection> > vah;
    event.getManyByType(vah);
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CaloShowerStepROCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CaloShowerStepROPrinter::Print(const art::Handle<CaloShowerStepROCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloShowerStepROPrinter::Print(const art::ValidHandle<CaloShowerStepROCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloShowerStepROPrinter::Print(const CaloShowerStepROCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CaloShowerStepROCollection has " << coll.size() << " steps\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CaloShowerStepROPrinter::Print(const art::Ptr<CaloShowerStepRO>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CaloShowerStepROPrinter::Print(const mu2e::CaloShowerStepRO& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(6) << ind;

  if(verbose()==1) {
    os 
      << " " << std::setw( 6) << obj.caloShowerStep()->volumeId() 
      << " " << std::setw( 6) << obj.ROID()
      << " " << std::setw( 8) << std::setprecision(2) << obj.time() 
      << " " << std::setw( 8) << std::setprecision(2) << obj.energy() 
      << std::endl;
  }
}

void 
mu2e::CaloShowerStepROPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CaloShowerStepROPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind    volID   ROID    Time   Energy\n";

}

