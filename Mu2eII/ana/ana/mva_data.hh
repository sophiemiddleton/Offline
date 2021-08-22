//
#ifndef __Mu2eII_ana_mva_data_hh__
#define __Mu2eII_ana_mva_data_hh__

#include "TString.h"
#include "TMVA/Reader.h"

namespace Mu2eII {

class mva_data {
public: 

  struct data_t {
    const  char* fName;
    const  char* fVarNames;            // names of the variables
    const  char* fXmlWeightsFile;
    double       fCutValue;
  };

  int            fTrainingCode;

  data_t         fData;

  float          fVar[100];            // 

  TMVA::Reader*  fMva;

  mva_data();

  mva_data(const char* Dataset, int MVATrainingCode);

  ~mva_data();

  const char* Name             () { return fData.fName;              }
  const char* VarNames         () { return fData.fVarNames;          }
  const char* XmlWeightsFile   () { return fData.fXmlWeightsFile;    }
  double      CutValue         () { return fData.fCutValue;          }
  int         TrainingCode     () { return fTrainingCode;            }

  double      Eval()              { return fMva->EvaluateMVA(fData.fName); }
};

};

#endif
