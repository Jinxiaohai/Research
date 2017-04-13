//-*- mode:utf-8 -*-
/********************************************************/
/*                 Created  by  xiaohai                 */
/*                 Telphone : 18501781924               */
/*            E-mail : jinxiaohai@sinap.ac.cn           */
/*            E-mail : xiaohaijin@outlook.com           */
/*   Address : Shanghai Institute of Applied Physics    */
/********************************************************/
//ampt
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
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
double PI = TMath::Pi();
double ANUMBER = 197.0;
int main(int argc, char *argv[])
{
  if (argc != 3)
    {
      cerr << "input three parameter, please !!!" << endl;
      return -1;
    }
  //head of event
  int nevent = 0, nrun = 0, multi = 0, NpartP = 0, NpartT = 0, NELP = 0, NINP = 0, NELT = 0, NINT = 0;
  double impactpar = 0., passhead = 0.;
  //head of track
  double px = 0., py = 0., pz =0., am = 0, x = 0, y = 0., z =0., time = 0.;
  double energy = 0.;
  string others1, others2, others3, others4;
  int id = 0;
  char filelist[512];
  int count = 0;
  int filenum = 0;
  
  ifstream streaminput(argv[1]);
  if(!streaminput.good())
    {
      cerr << "open fileinput error !!!" << endl;
      return -1;
    }
  TFile *file = new TFile(argv[2], "RECREATE");
  
  while (streaminput.good())
    {
      streaminput.getline(filelist,512);
      if(!streaminput.good())break;
      if (streaminput.good())
        {
          cout << filelist << endl;
          ifstream filein(filelist);
          int shijianshu = 0;
          while (filein.good())
            {
              int num = 0;
              filein >> nevent >> nrun >> multi >> impactpar >> NpartP >> NpartT >> NELP >> NINP >> NELT >> NINT >> passhead;
              if (!filein.good())break;
              ++shijianshu;
              for (int nlines = 0; nlines != multi; ++nlines)
                {
                  filein >> id >> px >> py >> pz >> am >> others1 >> others2 >> others3 >> others4;
                  ++num;
                  if(!filein.good())break;
                  //string::npos 为未找到。不等于未找到及找到了。
                  if (others1.find("*")!=string::npos || others2.find("*")!=string::npos || \
                      others3.find("*")!=string::npos || others4.find("*")!=string::npos)
                    {
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
                  TLorentzVector Mp(px, py, pz, energy);
                  TLorentzVector LR(x, y, z, time);
                  if(Mp.Pt()<10e-7)continue;



                }//track
              if (num != multi)
                {
                  cout << "error error error !!!!" << endl;
                }
              ++count;
            }//event
          if(shijianshu != 10)
            {
              cout << "shijianqueshao" << endl;
            }
          filein.close();
          ++filenum;
        }
    }
  cout << count << endl;
  streaminput.close();
  file->Write();
  return 0;
}
