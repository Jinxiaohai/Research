/********************************************************/
/*                 Created  by  xiaohai                 */
/*                 Telphone : 18501781924               */
/*            E-mail : jinxiaohai@sinap.ac.cn           */
/*            E-mail : xiaohaijin@outlook.com           */
/*   Address : Shanghai Institute of Applied Physics    */
/********************************************************/
//zpc
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
#include "head.h"
using namespace std;

double PI = TMath::Pi(  );
double ANUMBER = 197.0;

int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      std::cout << "need four parameter !!!" << std::endl;
      return -1;
    }
  //parameter
  int track_count;

  int nevent, nrun, multi, NELP, NINP, NELT, NINT;
  double impactpar;
  double id, px, py, pz, am, x, y, z, time;
  double energy;
  std::string others1, others2, others3, others4;
  double rap = 1000.;
  double pt = 1000.;
  TFile *file = new TFile(argv[2], "RECREATE");
  char filelist[512];
  int count = 0;
  int filenumber = 0;
  std::ifstream input(argv[1]);
  if(!input.good())
    {
      std::cerr << "open argv[1] error !!!" << std::endl;
      return -1;
    }
  std::ofstream output(argv[3]);
  if(!output.good())
    {
      std::cerr << "open argv[3] error !!!" << std::endl;
      return -1;
    }
  while (input.good())
    {
      input.getline(filelist,512);
      if (!input.good())break;
      if (input.good())
        {
          std::cout << "read in file " << filelist << std::endl;
          std::ifstream filein(filelist);
          int xiaohai = 0;
          while (filein.good())
            {
              filein >> nevent >> nrun >> multi >> impactpar >> NELP >> NINP >> NELT >> NINT;
              if(!filein.good())break;
              ++xiaohai;
              track_count = 0;
              int num = 0;
              for (int nlines = 0; nlines != multi; ++nlines)
                {
                  filein >> id >> px >> py >> pz >> am >> others1 >> others2 >> others3 >> others4;
                  if(!filein.good())break;
                  ++num;
                  //string::npos 为未找到。不等于未找到及找到了。
                  if (others1.find("*")!=std::string::npos || others2.find("*")!=std::string::npos || \
                      others3.find("*")!=std::string::npos || others4.find("*")!=std::string::npos)
                    {
                      // std::cout << "pipei" << std::endl;
                      x = 10000000;
                      y = 10000000;
                      z = 10000000;
                      time = 10000000;
                    }
                  else
                    {
                      x = atof(others1.c_str());
                      y = atof(others2.c_str());
                      z = atof(others3.c_str());
                      time = atof(others4.c_str());
                    }
                  energy = sqrt(px*px + py*py + pz*pz + am*am);
                  TVector3 R(x, y, z);
                  TVector3 P(px, py, pz);
                  TLorentzVector Mr(x, y, z, time);
                  TLorentzVector Mp(px, py, pz, energy);
                  if(Mp.Pt()<10e-7)continue;
                  rap = Mp.Rapidity();
                  pt = Mp.Pt();
                  
                  /// test
                  if (id == 1)
                    {
                      rap_d->Fill(rap);
                      pt_d->Fill(pt);
                    }
                  if (id == 2)
                    {
                      rap_u->Fill(rap);
                      pt_u->Fill(pt);
                    }
                  if (id == 3)
                    {
                      rap_s->Fill(rap);
                      pt_s->Fill(pt);
                    }
                  if (id == -1)
                    {
                      rap_anti_d->Fill(rap);
                      pt_anti_d->Fill(pt);
                    }
                  if (id == -2)
                    {
                      rap_anti_u->Fill(rap);
                      pt_anti_u->Fill(pt);
                    }
                  if (id == -3)
                    {
                      rap_anti_s->Fill(rap);
                      pt_anti_s->Fill(pt);
                    }




                  
                  ++track_count;
                }//track
              ++count;
              if (num != multi)
                {
                  cout << "error error error !!!!" << endl;
                  cout << count << endl;
                }
            }//event
          if(xiaohai != 10)
            {
              cout << " !!!!!!!!!!!!!!!!!!!!!!!" << endl;
            }
          filein.close();
        }
      ++filenumber;
    }
  input.close();
  std::cout << filenumber << "  files  " << count << "  events !!!" << std::endl;

  /// test
  rap_d->Scale(1./count);
  rap_d->Write();
  rap_u->Scale(1./count);
  rap_u->Write();
  rap_s->Scale(1./count);
  rap_s->Write();
  rap_anti_d->Scale(1./count);
  rap_anti_d->Write();
  rap_anti_u->Scale(1./count);
  rap_anti_u->Write();
  rap_anti_s->Scale(1./count);
  rap_anti_s->Write();

  pt_d->Scale(1./count);
  pt_d->Write();
  pt_u->Scale(1./count);
  pt_u->Write();
  pt_s->Scale(1./count);
  pt_s->Write();
  pt_anti_d->Scale(1./count);
  pt_anti_d->Write();
  pt_anti_u->Scale(1./count);
  pt_anti_u->Write();
  pt_anti_s->Scale(1./count);
  pt_anti_s->Write();

  
  file->Write();
  output.close();
  return 0;
}
