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
#include "ZPC.h"
#include "Particle.h"

using namespace std;
using namespace xiaohai;

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
  ZPC *zpc = new ZPC(chain);
  unsigned entries = chain->GetEntries();
  if(entries >= 1.e+008) return -1;

  unsigned int tracks = 0;

  for (unsigned l = 0; l != entries; ++l)
    {
      // /// progress
      // static char str[100] = "#";
      // unsigned int per_entry = 2 * (entries / 100);
      // if (l % per_entry == 0)
      //   {
      //     const int per = ceil(((static_cast<double>(l))/entries ) * 100);
      //     cout << left << "progress:[" << setw(50) << str
      //          << "] " << per << "%" << flush;
      //     strcat(str, "#");
      //     usleep(100000);
      //     cout << "\r" << flush;
      //     if (per >= 99)
      //       {
      //         cout << endl;
      //         cout << "Long live Chairman Mao." << endl;
      //       } 
      //   }
      // /// progress

      chain->GetEntry(l);
      tracks = zpc->Event_multi;

      vector<Particle> vecParticle;
      for (unsigned t = 0; t != tracks; ++t)
        {
          Particle particle(zpc->X[t], zpc->Y[t], zpc->Z[t]);
          vecParticle.push_back(particle);
        }

      /// get average x and y
      double x_cm = 0., y_cm = 0., z_cm = 0.;
      for (vector<Particle>::iterator iter = vecParticle.begin();
           iter != vecParticle.end(); ++iter)
        {
          x_cm += (*iter).GetX();
          y_cm += (*iter).GetY();
          z_cm += (*iter).GetZ();
        }
      x_cm = x_cm / vecParticle.size();
      y_cm = y_cm / vecParticle.size();
      z_cm = z_cm / vecParticle.size();
      // cout << x_cm << "   " << y_cm << "   " << z_cm << endl;

      /// 重新设定做标值
      for (vector<Particle>::iterator iter = vecParticle.begin();
           iter != vecParticle.end(); ++iter)
        {
          (*iter).SetX((*iter).GetX()-x_cm);
          (*iter).SetY((*iter).GetY()-y_cm);
          (*iter).SetZ((*iter).GetZ()-z_cm);
        }      

      /// calculate ecc
      double qx_2 = 0., qy_2 = 0., qx_3 = 0., qy_3 = 0., rSquare = 0.;
      double phi = 0., rsq = 0.;
      for (vector<Particle>::iterator iter = vecParticle.begin();
           iter != vecParticle.end(); ++iter)
        {
          phi = (*iter).GetPhi();
          rsq = (*iter).GetPt() * (*iter).GetPt();
          qx_2 += rsq * cos(2*phi);
          qy_2 += rsq * sin(2*phi);
          qx_3 += rsq * cos(3*phi);
          qy_3 += rsq * sin(3*phi);
          rSquare += rsq;
        }
      qx_2 = qx_2 / vecParticle.size();
      qy_2 = qy_2 / vecParticle.size();
      qx_3 = qx_3 / vecParticle.size();
      qy_3 = qy_3 / vecParticle.size();
      rSquare = rSquare / vecParticle.size();
      double epsilon2 = GetEpsilon2(qx_2, qy_2, rSquare);
      double epsilon3 = GetEpsilon3(qx_3, qy_3, rSquare);
      output_data << epsilon2 << "      " << epsilon3 << endl;
    }

  output_data.close();
  delete chain;
  delete zpc;
  delete saveFile;
  return 0;
}
