#include <iostream>
#include <fstream>
#include <sstream>

#include "TH3D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TBits.h"

#include "Particle.h"
#include "Head.h"

using namespace std;
using namespace xiaohai;

int main(int argc, char *argv[])
{
  ifstream inputList(argv[1]);
  char FileList[256];
  Event<double> *evt = new Event<double>;
  while(inputList.good())
    {
      inputList.getline(FileList, 256);
      if(inputList.fail()) break;
      ifstream inputData(FileList);
      while(inputData.good())
        {
          inputData >> *evt;
          if(inputData.fail()) break;
          cout << *evt << endl;
          // cout << evt->GetMulti() << endl;
          // for (int i = 0; i != evt->GetMulti(); ++i)
          //   {
          //     double PT = evt->GetTrack(i).GetPt();
          //     if(PT > 1e-7)testHist->Fill(PT);
          //   }
        }
      inputData.close();
    }
  TFile *File = new TFile("file.root", "RECREATE");
  testHist->Write();

  delete evt;
  inputList.close();
  return 0;
}
