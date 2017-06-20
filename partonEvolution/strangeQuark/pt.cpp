/********************************************************/
/*                 Created  by  xiaohai                 */
/*                 Telphone : 18501781924               */
/*            E-mail : jinxiaohai@sinap.ac.cn           */
/*            E-mail : xiaohaijin@outlook.com           */
/*   Address : Shanghai Institute of Applied Physics    */
/********************************************************/
//root
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

//c++
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <vector>
#include <fstream>
#include <sstream>
#include "AMPT.h"
#include "head.h"

using namespace std;
const double PI = TMath::Pi();

int main(int argc, char *argv[])
{
  if (argc != 3) return 0;
  char *inputFile = argv[1];
  char *outputFile = argv[2];
  
  TChain *chain = new TChain("AMPT");
  char FileList[1024];
  double dy = 1.0;
  double dpt = 0.2;

  ofstream output_data(outputFile);
  if (!output_data)
    {
      cerr << "parameter error !" << endl;
      return -1;
    }
  /// data list 的处理方式
  if ((static_cast<string>(inputFile)).find(".list") != string::npos)
    {
      ifstream input_data(inputFile);
      if (!input_data)
        {
          cerr << "parameter error !" << endl;
          return -1;
        }
      for (;input_data.good();)
        {
          input_data.getline(FileList,512);
          if  ( input_data.good() )
            {
              /// test the file is ok.
              TFile *ftmp = new TFile(FileList);
              if(!ftmp || !(ftmp->IsOpen()) || !(ftmp->GetNkeys()))
                {
                  cout << "The root file " << FileList << " is broken." << endl;
                }
              else
                {
                  cout << FileList << endl;
                  chain->Add(FileList);
                }
              delete ftmp;
            }
        }
      input_data.close();
    }
  else if((static_cast<string>(inputFile)).find(".root") != string::npos)
    {
      /// test the file is ok.
      TFile *ftmp = new TFile(inputFile);
      if(!ftmp || !(ftmp->IsOpen()) || !(ftmp->GetNkeys()))
        {
          cout << "The root file " << FileList << " is broken." << endl;
        }
      else
        {
          cout << inputFile << endl;
          chain->Add(inputFile);
        }
      delete ftmp;
    }
  else
    {
      cerr << "parameter is error, please check it." << endl;
      return -1;
    }

  TFile *saveFile = new TFile("saveFile.root", "RECREATE");
  AMPT *ampt = new AMPT(chain);
  unsigned entries = chain->GetEntries();
  if(entries >= 1.e+008) return -1;

  unsigned int tracks = 0;

  for (unsigned l = 0; l != entries; ++l)
    {
      /// progress
      static char str[100] = "#";
      unsigned int per_entry = 2 * (entries / 100);
      if (l % per_entry == 0)
        {
          const int per = ceil(((static_cast<double>(l))/entries ) * 100);
          cout << left << "progress:[" << setw(50) << str
               << "] " << per << "%" << flush;
          strcat(str, "#");
          usleep(100000);
          cout << "\r" << flush;
          if (per >= 99)
            {
              cout << endl;
              cout << "Long live Chairman Mao." << endl;
            } 
        }
      /// progress

      chain->GetEntry(l);
      tracks = ampt->Event_multi;
      for (unsigned t = 0; t != tracks; ++t)
        {
	  double Energy = sqrt(  pow(ampt->Mass[t], 2) + pow(ampt->Px[t], 2)
				 + pow(ampt->Py[t], 2) + pow(ampt->Pz[t], 2));
	  TLorentzVector lor(ampt->Px[t], ampt->Py[t], ampt->Pz[t], Energy);
	  if (lor.Rapidity() > -dy && lor.Rapidity() < dy)
	    {
	      if (abs(ampt->ID[t]) == 3)
		{
		  pt->Fill(lor.Pt(), 1./(2.*PI*lor.Pt()*dpt*2*dy));
		}
	      double PHI = atan2(ampt->Py[t], ampt->Px[t]);
	      double v2 = cos(PHI);
	      Strangev2->Fill(lor.Pt(), v2);
	    }
        }
    }

  pt->Scale(1./entries);
  pt->Write();
  Strangev2->Write();
  output_data.close();
  delete chain;
  delete ampt;
  delete saveFile;
  return 0;
}
