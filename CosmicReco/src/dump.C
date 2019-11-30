#include <TH1F.h>
#include <TFile.h>
#include <TObject.h>
#include <fstream>
using namespace std;



void dump(){

TFile *f = TFile::Open("test.root", “READ”);
if (!f) return; // just a precaution

TTree *t; f->GetObject("h_cem_Dio_mom__mom", t);

t->SetScanField(0);
t->Scan("*");
t->SaveAs("dump.csv");

delete f; // no longer needed (automatically deletes “t”)
}
