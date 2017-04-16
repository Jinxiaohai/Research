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
#include "Head.h"

using namespace std;

const double PI = TMath::Pi();
TH1D *rapidityhist = new TH1D("rapidityhist", "rapidityhist", 100, -6, 6);
TH1D *pthist = new TH1D("pthist", "pthist", 100, 0, 5);

int main(int argc, char *argv[])
{
  if (argc != 3) return 0;
  char *inputFile = argv[1];
  char *PsiFile = argv[2];
  
  TChain *chain = new TChain("AMPT");
  char FileList[1024];
  double RRMS = 3.294; /// GeV
  double G = 1./12;
  double dy = 1.;
  double dpt = 0.2;

  ifstream inputPsi(PsiFile);
  if (!inputPsi.good())
    {
      cerr << "open error !!!" << endl;
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
  // cout << entries << endl;

  unsigned int tracks = 0;
  double PHI2 = 0., PHI3 = 0.;

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
      tracks = ampt->Event_multi;

      /**
       * 这里开始构建参与者
       * 平面角.
       */
      inputPsi >> PHI2 >> PHI3;
      double EventNum = 0;

      
      for (unsigned t = 0; t != tracks-1; ++t)
        {
          // cout << ampt->ID[t] << endl;
          if (abs(ampt->ID[t]) == 3)
            {
              for (unsigned w = t+1; w != tracks; ++w)
                {
                  if((ampt->ID[w] == 3 && ampt->ID[t] == -3)
                     || (ampt->ID[w] == -3 && ampt->ID[t] == 3))
                    {
                      // cout << ampt->ID[w] << endl;
                      double tpx = ampt->Px[t], tpy = ampt->Py[t], tpz = ampt->Pz[t], tmass = ampt->Mass[t];
                      double tx = ampt->X[t], ty = ampt->Y[t], tz = ampt->Z[t], ttime = ampt->Time[t];
                      double tenergy = sqrt(tpx*tpx + tpy*tpy + tpz*tpz + tmass*tmass);
                      double wpx = ampt->Px[w], wpy = ampt->Py[w], wpz = ampt->Pz[w], wmass = ampt->Mass[w];
                      double wx = ampt->X[w], wy = ampt->Y[w], wz = ampt->Z[w], wtime = ampt->Time[w];
                      double wenergy = sqrt(wpx*wpx + wpy*wpy + wpz*wpz + wmass*wmass);

                      TLorentzVector tlabmomentum(tpx, tpy, tpz, tenergy);
                      TLorentzVector tlabcoordinate(tx, ty, tz, ttime);
                      TLorentzVector wlabmomentum(wpx, wpy, wpz, wenergy);
                      TLorentzVector wlabcoordinate(wx, wy, wz, wtime);

                      double rapidity = 1./2. * log((tenergy+wenergy+tpz+wpz) / (tenergy+wenergy-tpz-wpz));
                      double transpt = (tlabmomentum + wlabmomentum).Pt();
                      double betax = -(tlabmomentum.Px() + wlabmomentum.Px()) / (tlabmomentum.E() + wlabmomentum.E());
                      double betay = -(tlabmomentum.Py() + wlabmomentum.Py()) / (tlabmomentum.E() + wlabmomentum.E());
                      double betaz = -(tlabmomentum.Pz() + wlabmomentum.Pz()) / (tlabmomentum.E() + wlabmomentum.E());

                      TLorentzRotation ROTA(betax, betay, betaz);

                      TLorentzVector tcenmomentum = ROTA * tlabmomentum;
                      TVector3 tcenmomentum3D = tcenmomentum.Vect();
                      TLorentzVector tcencoordinate = ROTA * tlabcoordinate;
                      TVector3 tcencoordinate3D = tcencoordinate.Vect();
                      
                      TLorentzVector wcenmomentum = ROTA * wlabmomentum;
                      TVector3 wcenmomentum3D = wcenmomentum.Vect();
                      TLorentzVector wcencoordinate = ROTA * wlabcoordinate;
                      TVector3 wcencoordinate3D = wcencoordinate.Vect();

                      double TimeMax = TMath::Max(tcencoordinate.T(), wcencoordinate.T());
                      tcencoordinate3D = tcenmomentum3D * (1./tcenmomentum.E())
                        * (TimeMax-tcencoordinate.T()) + tcencoordinate3D;
                      wcencoordinate3D = wcenmomentum3D * (1./wcenmomentum.E())
                        * (TimeMax-wcencoordinate.T()) + wcencoordinate3D;
                      double rSquare = (tcencoordinate3D-wcencoordinate3D).Mag2();
                      rSquare = 25.6889*rSquare;
                      double kSquare = ((1./2.) * (tcenmomentum3D-wcenmomentum3D)).Mag2();
                      double sigmaSquare = (4./3.) * RRMS * RRMS;

                      double Fb = G * 8 * exp(-rSquare/sigmaSquare - kSquare*sigmaSquare);
                      // cout << Fb << endl;

                      /// calculate flow
                      double pxval = (tlabmomentum+wlabmomentum).Px();
                      double pyval = (tlabmomentum+wlabmomentum).Py();
                      double hadron_phi = atan2(pyval, pxval);
                      double V2 = cos(2.*(hadron_phi-PHI2));
                      double V3 = cos(3.*(hadron_phi-PHI3));
                      
                      
                      rapidityhist->Fill(rapidity, Fb);
                      pthist->Fill(transpt, Fb);
                      if (abs(rapidity) > 1.)
                        {
                          continue;
                        }
                      
                      EventNum += Fb;
                      
                      v2SumHistgram->Fill(transpt, V2*Fb);
                      v3SumHistgram->Fill(transpt, V3*Fb);
                      probablySum->Fill(transpt,   Fb);
                        
                      if (transpt > 1e-7)
                        {
                          yield->Fill(transpt, Fb*1./(2.*PI*transpt*dpt*2*dy));
                        }
                    }// in if
                }//for
            }//out if
        }//for
      EventNumHist->Fill(EventNum);
    }//for

  rapidityhist->Scale(1./entries);
  rapidityhist->Write();
  pthist->Scale(1./entries);
  pthist->Write();
  fb->Write();
  yield->Scale(1./entries);
  yield->Write();

  EventNumHist->Write();
  v2SumHistgram->Write();
  v3SumHistgram->Write();
  probablySum->Write();

  saveFile->Write();
  
  inputPsi.close();
  delete chain;
  delete ampt;
  delete saveFile;
  return 0;
}
