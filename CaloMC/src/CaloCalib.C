// run in ROOT
#include <fstream>
#include <iostream>
#include <algorithm>
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <time.h>
using namespace RooFit;
TString filepath = "/Users/user/ski_examples/CaloCalib/10M.root";
using namespace std;
using namespace TMath;

//================== Select the sample =========================//
TTree* get_data_tree(){
    TFile *f =  new TFile(filepath);
    TTree *t = (TTree*)f->Get("CaloExample/Calo");
    return t;
}

int main(){
  for(unsigned int numVal = 700; numVal < 1348; numVal++){
      std::cout<<" Crystal "<<numVal<<std::endl;
    TString cryNum;
    cryNum = to_string(numVal);

    // create the output file.
    TString ouptName = "mu2e_caloSimu_crySpec" + cryNum + ".root";
    TString ergTitle = "Energy Spectrum of Crystal " + cryNum;
    TString ratTitle = "The Ratio of Energy Deposited in Crystal " + cryNum;
    TFile *ouptFile = new TFile(ouptName, "RECREATE");
    TH1F *h_spec = new TH1F("h_spec", ergTitle, 300, 0.0, 7.5);
    h_spec -> SetXTitle("Energy Deposited [MeV]");
    h_spec -> SetYTitle("#Event / (25 keV)");
    TH1F *h_ratio = new TH1F("h_ratio", ratTitle, 120, 0.0, 1.2);
    h_ratio -> SetXTitle("Ratio of Energy Deposited");
    h_ratio -> SetYTitle("#Event / (0.01)");
    TH1F *h_ntrg = new TH1F("h_ntrg", "Number of Hits on Target Crystal", 10, 0, 10);
    h_ntrg -> SetXTitle("number of target crystal");
    TH1F *h_stim = new TH1F("h_stim", "Time Difference from Same Crystal", 200, 0.0, 100.0);
    h_stim -> SetXTitle("#DeltaT [ns]");
    TH1F *h_time = new TH1F("h_time", "Maxium Time Difference for All Hits", 200, 0.0, 100.0);
    h_time -> SetXTitle("#DeltaT [ns]");
    //TH1F *h_tErg = new TH1F("h_tErg", "Energy Spectrum with Special Time Difference", 300, 0.0, 300);
    TH1F *h_tErg = new TH1F("h_tErg", "Energy Spectrum with Special Time Difference", 300, 0.0, 7.5);
    h_tErg -> SetXTitle("Energy Deposited [MeV]");
    h_tErg -> SetYTitle("#Event / (25 keV)");


    TTree *ergTree = get_data_tree();

    float cryEdep[100];
    int nCry;
    int cryId[100];
    float cryTime[100];
    float cryPosX[100], cryPosY[100], cryPosZ[100];

    ergTree -> SetBranchAddress("nCry", &nCry);
    ergTree -> SetBranchAddress("cryId", cryId);
    ergTree -> SetBranchAddress("cryEdep", cryEdep);
    ergTree -> SetBranchAddress("cryTime", cryTime);
    ergTree -> SetBranchAddress("cryPosX", &cryPosX);
    ergTree -> SetBranchAddress("cryPosY", &cryPosY);
    ergTree -> SetBranchAddress("cryPosZ", &cryPosZ);

    unsigned int nEvt = (int)ergTree -> GetEntries();
    for(unsigned int iEvt=0; iEvt<nEvt; iEvt++)
    {
      ergTree -> GetEntry(iEvt);
      int idExist = 0;
      std::vector<float> sameCryTime;
      std::vector<float> sameCryEdep;
      for(int icry=0; icry<nCry; icry++)
      {
        if(cryId[icry] == numVal)
        {
          idExist += 1;
          sameCryTime.push_back(cryTime[icry]);
          sameCryEdep.push_back(cryEdep[icry]);
        }
      }
      if(idExist == 0) continue;
      sort(sameCryTime.begin(), sameCryTime.end());
      float sameCryDeltaT = sameCryTime.back() - sameCryTime.front();
      h_stim -> Fill(sameCryDeltaT);

      float ergTarget = 0.0;
      float ergRest   = 0.0;
      int nTarget = 0;
      std::vector<float> deltaTime;

      for(int icry=0; icry<nCry; icry++)
      {
        if(cryId[icry] == numVal)
        {
          ergTarget = cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);
          nTarget += 1;
        }
        else
        {
          ergRest += cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);
        }
      }
      sort(deltaTime.begin(), deltaTime.end());
      float difTime = deltaTime.back() - deltaTime.front();
      int numHit = sameCryEdep.size();
      if(difTime > 20)
      {
        for(int iHit = 0; iHit < numHit; iHit++)
        {
          h_tErg -> Fill(sameCryEdep[iHit]);
        }
      }
      if((ergTarget / (ergTarget + ergRest)) >= 0.8 && difTime < 4)
      {
        h_spec -> Fill(ergTarget);
      }
      h_ratio -> Fill(ergTarget / (ergTarget + ergRest));
      h_ntrg -> Fill(nTarget);
      h_time -> Fill(difTime);
    }

  ouptFile -> Write();
  ouptFile -> Close();

  }
    return 0;
}

