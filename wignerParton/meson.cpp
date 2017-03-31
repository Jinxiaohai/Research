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
#include "TLorentzRotation.h"
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
#include "Particle.h"

using namespace std;

int main(int argc, char *argv[])
{
  if (argc != 3) return 0;
  char *inputFile = argv[1];
  char *outputFile = argv[2];
  
  TChain *chain = new TChain("AMPT");
  char FileList[1024];

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
          if (abs(ampt->ID[t]) == 3)
            {
              for (unsigned w = t+1; w != tracks; ++w)
                {
                  if((ampt->ID[w] == 3 && ampt->ID[t] == -3)
                     || (ampt->ID[w] == -3 && ampt->ID[t] == 3))
                    {
                      xiaohai::Particle labTParticle(ampt->ID[t], ampt->Px[t], ampt->Py[t],
                                                     ampt->Pz[t], ampt->Mass[t], ampt->X[t],
                                                     ampt->Y[t], ampt->Z[t], ampt->Time[t]);
                      xiaohai::Particle labWParticle(ampt->ID[w], ampt->Px[w], ampt->Py[w],
                                                     ampt->Pz[w], ampt->Mass[w], ampt->X[w],
                                                     ampt->Y[w], ampt->Z[w], ampt->Time[w]);
              
                      double totalEnergy = labTParticle.GetEnergy() + labWParticle.GetEnergy();
                      double betax = -(labTParticle.GetPx() + labWParticle.GetPx()) / totalEnergy;
                      double betay = -(labTParticle.GetPy() + labWParticle.GetPy()) / totalEnergy;
                      double betaz = -(labTParticle.GetPz() + labWParticle.GetPz()) / totalEnergy;

                      xiaohai::Particle cenTParticle = labTParticle.Boost(betax, betay, betaz);
                      xiaohai::Particle cenWParticle = labWParticle.Boost(betax, betay, betaz);
                      
                    }// in if
                }//for
            }//out if
        }//for
    }//for

  output_data.close();
  delete chain;
  delete ampt;
  delete saveFile;
  return 0;
}
