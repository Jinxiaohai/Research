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
#include "Head.h"

using namespace std;

const double PI = TMath::Pi();
TH1D *rapidityhist = new TH1D("rapidityhist", "rapidityhist", 100, -6, 6);
TH1D *pthist = new TH1D("pthist", "pthist", 100, 0, 5);

int main(int argc, char *argv[])
{
  if (argc != 3) return 0;
  char *inputFile = argv[1];
  char *outputFile = argv[2];
  
  TChain *chain = new TChain("AMPT");
  char FileList[1024];
  double RRMS = 6.082;
  double G = 1./54;
  double dy = 1.;
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
      for (unsigned t = 0; t != tracks-2; ++t)
        {
          if (abs(ampt->ID[t]) == 3)
            {
              for (unsigned w = t+1; w != tracks-1; ++w)
                {
                  if((ampt->ID[w] == 3 && ampt->ID[t] == 3)
                     || (ampt->ID[w] == -3 && ampt->ID[t] == -3))
                    {
                      for (unsigned int h = w+1; h != tracks; ++h)
                        {
                          if ((ampt->ID[t] == 3 && ampt->ID[w] == 3 && ampt->ID[h] == 3)
                              || (ampt->ID[t] == -3 && ampt->ID[w] == -3 && ampt->ID[h] == -3))
                            {
                              
                              xiaohai::Particle labTParticle(ampt->ID[t], ampt->Px[t], ampt->Py[t],
                                                             ampt->Pz[t], ampt->Mass[t], ampt->X[t],
                                                             ampt->Y[t], ampt->Z[t], ampt->Time[t]);
                              xiaohai::Particle labWParticle(ampt->ID[w], ampt->Px[w], ampt->Py[w],
                                                             ampt->Pz[w], ampt->Mass[w], ampt->X[w],
                                                             ampt->Y[w], ampt->Z[w], ampt->Time[w]);
                              xiaohai::Particle labHParticle(ampt->ID[h], ampt->Px[h], ampt->Py[h],
                                                             ampt->Pz[h], ampt->Mass[h], ampt->X[h],
                                                             ampt->Y[h], ampt->Z[h], ampt->Time[h]);
              
                              double Rapidity = 1./2. * log((labTParticle.GetEnergy()+labWParticle.GetEnergy()+labHParticle.GetEnergy()
                                                             + labTParticle.GetPz() + labWParticle.GetPz() + labHParticle.GetPz())
                                                            / (labTParticle.GetEnergy()+labWParticle.GetEnergy()+labHParticle.GetEnergy()
                                                               -labTParticle.GetPz() - labWParticle.GetPz() - labHParticle.GetPz()));
                              double transpt = sqrt( pow(labTParticle.GetPx() + labWParticle.GetPx() + labHParticle.GetPx(), 2)
                                                     + pow(labTParticle.GetPy() + labWParticle.GetPy() + labHParticle.GetPy(), 2));
                              if (abs(Rapidity) > dy){
                                continue;
                              }
                              double betax = -(labTParticle.GetPx()+labWParticle.GetPx()+labHParticle.GetPx())
                                /(labTParticle.GetEnergy()+labWParticle.GetEnergy()+labHParticle.GetEnergy());
                              double betay = -(labTParticle.GetPy()+labWParticle.GetPy()+labHParticle.GetPy())
                                /(labTParticle.GetEnergy()+labWParticle.GetEnergy()+labHParticle.GetEnergy());
                              double betaz = -(labTParticle.GetPz()+labWParticle.GetPz()+labHParticle.GetPz())
                                /(labTParticle.GetEnergy()+labWParticle.GetEnergy()+labHParticle.GetEnergy());

                              xiaohai::Particle cenTParticle = labTParticle.Boost(betax, betay, betaz);
                              xiaohai::Particle cenWParticle = labWParticle.Boost(betax, betay, betaz);
                              xiaohai::Particle cenHParticle = labHParticle.Boost(betax, betay, betaz);

                              double TMAX = cenTParticle.GetTime() > (cenWParticle.GetTime() > cenHParticle.GetTime() ?
                                                                      cenWParticle.GetTime() : cenHParticle.GetTime()) ?
                                cenTParticle.GetTime() : (cenWParticle.GetTime() > cenHParticle.GetTime() ?
                                                          cenWParticle.GetTime() : cenHParticle.GetTime());
                              /// equal time                      
                              cenWParticle.SetX(cenWParticle.GetX()
                                                + cenWParticle.GetBetaX() * (TMAX - cenWParticle.GetTime()));
                              cenWParticle.SetY(cenWParticle.GetY()
                                                + cenWParticle.GetBetaY() * (TMAX - cenWParticle.GetTime()));
                              cenWParticle.SetZ(cenWParticle.GetZ()
                                                + cenWParticle.GetBetaZ() * (TMAX - cenWParticle.GetTime()));

                              cenTParticle.SetX(cenTParticle.GetX()
                                                + cenTParticle.GetBetaX() * (TMAX - cenTParticle.GetTime()));
                              cenTParticle.SetY(cenTParticle.GetY()
                                                + cenTParticle.GetBetaY() * (TMAX - cenTParticle.GetTime()));
                              cenTParticle.SetZ(cenTParticle.GetZ()
                                                + cenTParticle.GetBetaZ() * (TMAX - cenTParticle.GetTime()));

                              cenHParticle.SetX(cenHParticle.GetX()
                                                + cenHParticle.GetBetaX() * (TMAX - cenHParticle.GetTime()));
                              cenHParticle.SetY(cenHParticle.GetY()
                                                + cenHParticle.GetBetaY() * (TMAX - cenHParticle.GetTime()));
                              cenHParticle.SetZ(cenHParticle.GetZ()
                                                + cenHParticle.GetBetaZ() * (TMAX - cenHParticle.GetTime()));

                              xiaohai::ThreeDimensionVector<double> r1(cenTParticle.GetX(),
                                                                       cenTParticle.GetY(),
                                                                       cenTParticle.GetZ());
                              // cout << r1.GetX() << "  " << r1.GetY() << "  " << r1.GetZ() << endl;
                              xiaohai::ThreeDimensionVector<double> r2(cenWParticle.GetX(),
                                                                       cenWParticle.GetY(),
                                                                       cenWParticle.GetZ());
                              // cout << r2.GetX() << "  " << r2.GetY() << "  " << r2.GetZ() << endl;
                              xiaohai::ThreeDimensionVector<double> r3(cenHParticle.GetX(),
                                                                       cenHParticle.GetY(),
                                                                       cenHParticle.GetZ());
                              
                              xiaohai::ThreeDimensionVector<double> k1(cenTParticle.GetPx(),
                                                                       cenTParticle.GetPy(),
                                                                       cenTParticle.GetPz());
                              // cout << k1.GetX() << "  " << k1.GetY() << "  " << k1.GetZ() << endl;
                              xiaohai::ThreeDimensionVector<double> k2(cenWParticle.GetPx(),
                                                                       cenWParticle.GetPy(),
                                                                       cenWParticle.GetPz());
                              // cout << k2.GetX() << "  " << k2.GetY() << "  " << k2.GetZ() << endl;
                              xiaohai::ThreeDimensionVector<double> k3(cenHParticle.GetPx(),
                                                                       cenHParticle.GetPy(),
                                                                       cenHParticle.GetPz());
                              xiaohai::ThreeDimensionVector<double> rho = 1./sqrt(2.) * (r1 - r2);
                              xiaohai::ThreeDimensionVector<double> lambda = 1./sqrt(6.) * (r1 + r2 - 2.*r3);
                              xiaohai::ThreeDimensionVector<double> krho = 1./sqrt(2.) * (k1 - k2);
                              xiaohai::ThreeDimensionVector<double> klambda = 1./sqrt(6.) * (k1 + k2 - 2.*k3);
                              double rhoSquare = rho.GetMag2() * 25.6889;
                              double lambdaSquare = lambda.GetMag2() * 25.6889;
                              double krhoSquare = krho.GetMag2();
                              double klambdaSquare = klambda.GetMag2();
                              
                              double omegaSquare = RRMS*RRMS;
                              double FB = G * 8*8 * exp(-(rhoSquare+lambdaSquare)/omegaSquare - (krhoSquare+klambdaSquare)*omegaSquare);
                              // cout << Rho << endl;
                              // fb->Fill(Rho);
                              rapidityhist->Fill(Rapidity, FB);
                              if (transpt > 1e-6){
                                pthist->Fill(transpt, FB);
                              }
                              if (transpt > 1e-5)
                                {
                                  yield->Fill(transpt, FB*1./(2.*PI*transpt*dpt*2*dy));
                                }
                            }//if
                        }// for
                    }//if
                }//for
            }//out if
        }//for
    }//for

  rapidityhist->Write();
  pthist->Write();
  fb->Write();
  yield->Scale(1./entries);
  yield->Write();
  output_data.close();
  delete chain;
  delete ampt;
  delete saveFile;
  return 0;
}
