//
// swig interface file to wrap c++ code for python
//

%module GeneralUtilities

%include "std_vector.i"
%include "stdint.i"
%include "std_map.i"
%template(vec_double) std::vector<double>;

%{
#include <map>
#include <vector>
#include "Offline/GeneralUtilities/inc/SplineInterpolation.hh"
#include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"
%}

%include "Offline/GeneralUtilities/inc/SplineInterpolation.hh"
%include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"
