//
// swig interface file to wrap c++ code for python
//

%module DataProducts

%include "std_string.i"
%include "math.i"
%include "stdint.i"
%include "std_vector.i"
%include "std_map.i"

%apply const std::string& {std::string* foo};
%ignore operator <<;
%ignore operator int8_t;
//%template(EnumToStringSpars) mu2e::SurfaceIdDetail;
%{
#include<vector>
#include<map>
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"
%}
%include "Offline/DataProducts/inc/StrawEnd.hh"
%include "Offline/DataProducts/inc/StrawId.hh"
%include "Offline/DataProducts/inc/SurfaceId.hh"
%include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"
