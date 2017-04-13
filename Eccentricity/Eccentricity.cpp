//root
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class PlotFile;
#endif
#ifndef __CINT__
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
#endif

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
#include "Particle.h"

using namespace std;
using namespace xiaohai;

int main(int argc, char *argv[])
{
  if (argc != 3) return 0;
  char *inputFile = argv[1];
  char *outputFile = argv[2];

  char FileList[1024];
  /// event head
  unsigned int eventNum = 0, zeroNum = 0, partP = 0, partT = 0;
  double bimpactor = 0.;
  /// tracks
  double x = 0., y = 0., z = 0.;
  int num = 0, status = 0, currentFlavor = 0, originalFlavor = 0;

  ofstream output_data(outputFile);
  if (!output_data)
    {
      cerr << "parameter error !" << endl;
      return -1;
    }
  ifstream streaminput(inputFile);
  if (!streaminput.good())
    {
      cerr << "open error !!!" << endl;
      return -1;
    }

  while (streaminput)
    {
      streaminput.getline(FileList, 1024);
      ifstream input_data(FileList);
      while (input_data)
        {
          input_data >> eventNum >> zeroNum >> partP >> partT >> bimpactor;
          if(!input_data)
            {
              break;
            }
          vector<Particle> vecParticle;
          for (unsigned i = 0; i != partP + partT; ++i)
            {
              input_data >> x >> y >> num >> status >> z
                         >> currentFlavor >> originalFlavor;
              if (status != 0)
                {
                  Particle particle(x, y, z);
                  vecParticle.push_back(particle);
                }
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
          cout << x_cm << "   " << y_cm << "   " << z_cm << endl;

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
      input_data.close();
    }
  output_data.close();
  streaminput.close();
  return 0;
}
