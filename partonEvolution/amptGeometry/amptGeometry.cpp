/**
 * @file   partoninitial.cpp
 * @author xiaohai <xiaohaijin@outlook.com>
 * @date   Mon Apr 17 08:59:36 2017
 * 
 * @brief  the file was created to get the epsilon and psi
 *         from parton initial data
 */
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
#include "./include/Particle.h"

using namespace std;
using namespace xiaohai;
void Process(vector<Particle>&, vector<Particle>&, double&, double&, double&, double&);

int main(int argc, char *argv[])
{
  char *inputFile = argv[1];
  char *outputFile = argv[2];
  if (argc != 3)
    {
      cout << "needs three parameters ." << endl;
      return -1;
    }
  char fileList[512];
  
  ifstream inputData(inputFile);
  if (!inputData)
    {
      cerr << "open error !!!" << endl;
      return -1;
    }
  ofstream outputData(outputFile);
  if (!outputData)
    {
      cerr << "open error !!!" << endl;
      return -1;
    }
  ofstream psi2_file("./out/psi_2.txt");
  ofstream psi3_file("./out/psi_3.txt");
  ofstream epsilon2_file("./out/epsilon_2.txt");
  ofstream epsilon3_file("./out/epsilon_3.txt");
  //Read input file
  //File format for each event:
  //   <iaevt> <miss> <hadronNum> <bimpact> <npartP> <npartT> <nelp> <ninp> <nelt> <nint>
  //   <ityp> <px> <py> <pz> <xmass> <gx> <gy> <gz> <ft>
  //   ...
  //Note: accroding to the value of pz to choice the forward and backward particle.
  unsigned int iaevt, miss, hadronNum, npartP, npartT, nelp, ninp, nelt, nint;
  double bimpact, passhead;
  int ityp;
  double px, py, pz, xmass, gx, gy, gz, ft;
  
  /// single event
  double e2 = 0., e3 = 0., phi2 = 0., phi3 = 0.;
  /// average
  double e2Ave = 0., e3Ave = 0., psi2Ave = 0., psi3Ave = 0.;
  /// 事件数
  unsigned int eventCount = 0;
  /// average participant nucleons
  unsigned int forwardNum = 0, backwardNum = 0;
  while (inputData)
    {
      inputData.getline(fileList, 512);
      if (!inputData) break;
      ifstream input(fileList);
      while (input)
        {
          input >> iaevt >> miss >> hadronNum >> bimpact >> npartP >> npartT
                >> nelp >> ninp >> nelt >> nint >> passhead;
          if(!input) break;
          vector<Particle> particles_forward;
          vector<Particle> particles_backward;
          ++eventCount;
          for (unsigned i = 0; i != hadronNum; ++i)
            {
              input >> ityp >> px >> py >> pz >> xmass >> gx
                    >> gy >> gz >> ft;
              Particle particle(gx, gy, gz);
              if (pz > 0)
                {
                  particles_forward.push_back(particle);
                }
              else if (pz < 0)
                {
                  particles_backward.push_back(particle);
                }
            }//track
          Process(particles_forward, particles_backward, e2, e3, phi2, phi3);
          forwardNum += particles_forward.size();
          backwardNum += particles_backward.size();
          e2Ave += e2;
          e3Ave += e3;
          psi2Ave += phi2;
          psi3Ave += phi3;
          
          psi2_file << phi2 << endl;
          psi3_file << phi3 << endl;
          epsilon2_file << e2 << endl;
          epsilon3_file << e3 << endl;
        }//event
      input.close();
    }//file
  outputData << "average forward partcles number ==>  "
       << forwardNum / static_cast<double>(eventCount) << endl;
  outputData << "average backward particles number ==>  "
       << backwardNum / static_cast<double>(eventCount) << endl;
  outputData << "average epsilon2 ==>  " << e2Ave / eventCount << endl;
  outputData << "average epsilon3 ==>  " << e3Ave / eventCount << endl;
  outputData << "average psi2 ==>  " << psi2Ave / eventCount << endl;
  outputData << "average psi3 ==>  " << psi3Ave / eventCount << endl;

  inputData.close();
  outputData.close();
  psi2_file.close();
  psi3_file.close();
  epsilon2_file.close();
  epsilon3_file.close();
  return 0;
}

void Process(vector<Particle> &particles_forward, vector<Particle> &particles_backward,
             double &e2, double &e3, double &phi2, double &phi3)
{
  /// change the value of Withproj and Withtarg to get what you want.
  /// the value of (Withproj || withtarg) must be true, otherwise, the
  /// program can not work.
  bool Forward = true;
  bool Backward = true;
  double x_cm = 0.;
  double y_cm = 0.;

  unsigned int Count = 0;

  if(Forward)
    {
      for (vector<Particle>::iterator iter = particles_forward.begin();
           iter != particles_forward.end(); ++iter)
        {
          x_cm += iter->GetX();
          y_cm += iter->GetY();
        }
      Count += particles_forward.size();
    }
  if (Backward)
    {
      for (vector<Particle>::iterator iter = particles_backward.begin();
           iter != particles_backward.end(); ++iter)
        {
          x_cm += iter->GetX();
          y_cm += iter->GetY();
        }
      Count += particles_backward.size();
    }

  x_cm = x_cm / static_cast<double>(Count);
  y_cm = y_cm / static_cast<double>(Count);

  if (Forward)
    {
      for (vector<Particle>::iterator iter = particles_forward.begin();
           iter != particles_forward.end(); ++iter)
        {
          iter->SetX(iter->GetX()-x_cm);
          iter->SetY(iter->GetY()-y_cm);
        }
    }
  if (Backward)
    {
      for (vector<Particle>::iterator iter = particles_backward.begin();
           iter != particles_backward.end(); ++iter)
        {
          iter->SetX(iter->GetX()-x_cm);
          iter->SetY(iter->GetY()-y_cm);
        }
    }

  double qx_2 = 0;
  double qy_2 = 0;
  double qx_3 = 0;
  double qy_3 = 0;
  double phi = 0;
  double rsq = 0;
  double r = 0;

  if (Forward)
    {
      for (vector<Particle>::iterator iter=particles_forward.begin();
           iter != particles_forward.end(); ++iter)
        {
          r = iter->GetPt();
          phi = iter->GetPhi();

          qx_2 += pow(r, 2)*cos(2.*phi);
          qy_2 += pow(r, 2)*sin(2.*phi);
          qx_3 += pow(r, 2)*cos(3.*phi);
          qy_3 += pow(r, 2)*sin(3.*phi);
          rsq += pow(r, 2);
        }
    }
  if (Backward)
    {
      for (vector<Particle>::iterator iter=particles_backward.begin();
           iter != particles_backward.end(); ++iter)
        {
          r = iter->GetPt();
          phi = iter->GetPhi();

          qx_2 += pow(r, 2)*cos(2.*phi);
          qy_2 += pow(r, 2)*sin(2.*phi);
          qx_3 += pow(r, 2)*cos(3.*phi);
          qy_3 += pow(r, 2)*sin(3.*phi);
          rsq += pow(r, 2);
        }
    }

  qx_2 = qx_2 / static_cast<double>(Count);
  qy_2 = qy_2 / static_cast<double>(Count);
  qx_3 = qx_3 / static_cast<double>(Count);
  qy_3 = qy_3 / static_cast<double>(Count);
  rsq = rsq / static_cast<double>(Count);

  e2 = GetEpsilon2(qx_2, qy_2, rsq);
  e3 = GetEpsilon3(qx_3, qy_3, rsq);
  phi2 = GetPhi(qx_2, qy_2, 2);
  phi3 = GetPhi(qx_3, qy_3, 3);
}
