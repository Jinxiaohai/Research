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
  double RRMS = 6.082;
  double G = 1./54;
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
      
      
      for (unsigned t = 0; t != tracks-2; ++t)
        {
          if (abs(ampt->ID[t]) == 3)
            {
              // cout << l << "  ======>  " << ampt->ID[t] << endl;
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
                              double tpx = ampt->Px[t], tpy = ampt->Py[t], tpz = ampt->Pz[t], tmass = ampt->Mass[t];
                              double tx = ampt->X[t], ty = ampt->Y[t], tz = ampt->Z[t], ttime = ampt->Time[t];
                              double tenergy = sqrt(tpx*tpx + tpy*tpy + tpz*tpz + tmass*tmass);
                              double wpx = ampt->Px[w], wpy = ampt->Py[w], wpz = ampt->Pz[w], wmass = ampt->Mass[w];
                              double wx = ampt->X[w], wy = ampt->Y[w], wz = ampt->Z[w], wtime = ampt->Time[w];
                              double wenergy = sqrt(wpx*wpx + wpy*wpy + wpz*wpz + wmass*wmass);
                              double hpx = ampt->Px[h], hpy = ampt->Py[h], hpz = ampt->Pz[h], hmass = ampt->Mass[h];
                              double hx = ampt->X[h], hy = ampt->Y[h], hz = ampt->Z[h], htime = ampt->Time[h];
                              double henergy = sqrt(hpx*hpx + hpy*hpy + hpz*hpz + hmass*hmass);
                              // cout << tmass << "  " << wmass << "  " << hmass << "  " << endl;

                              TLorentzVector tlabmomentum(tpx, tpy, tpz, tenergy);
                              TLorentzVector tlabcoordinate(tx, ty, tz, ttime);
                              TLorentzVector wlabmomentum(wpx, wpy, wpz, wenergy);
                              TLorentzVector wlabcoordinate(wx, wy, wz, wtime);
                              TLorentzVector hlabmomentum(hpx, hpy, hpz, henergy);
                              TLorentzVector hlabcoordinate(hx, hy, hz, htime);

                              double rapidity = 1./2. * log((tenergy+wenergy+henergy + tpz+wpz+hpz) /
                                                            (tenergy+wenergy+henergy -tpz-wpz-hpz));
                              double transpt = (tlabmomentum + wlabmomentum + hlabmomentum).Pt();
                              double betax = -(tlabmomentum.Px() + wlabmomentum.Px() + hlabmomentum.Px())
                                / (tlabmomentum.E() + wlabmomentum.E() + hlabmomentum.E());
                              double betay = -(tlabmomentum.Py() + wlabmomentum.Py() + hlabmomentum.Py())
                                / (tlabmomentum.E() + wlabmomentum.E() + hlabmomentum.E());
                              double betaz = -(tlabmomentum.Pz() + wlabmomentum.Pz() + hlabmomentum.Pz())
                                / (tlabmomentum.E() + wlabmomentum.E() + hlabmomentum.E());

                              TLorentzRotation ROTA(betax, betay, betaz);

                              TLorentzVector tcenmomentum = ROTA * tlabmomentum;
                              TVector3 tcenmomentum3D = tcenmomentum.Vect();
                              TLorentzVector tcencoordinate = ROTA * tlabcoordinate;
                              TVector3 tcencoordinate3D = tcencoordinate.Vect();
                      
                              TLorentzVector wcenmomentum = ROTA * wlabmomentum;
                              TVector3 wcenmomentum3D = wcenmomentum.Vect();
                              TLorentzVector wcencoordinate = ROTA * wlabcoordinate;
                              TVector3 wcencoordinate3D = wcencoordinate.Vect();

                              TLorentzVector hcenmomentum = ROTA * hlabmomentum;
                              TVector3 hcenmomentum3D = hcenmomentum.Vect();
                              TLorentzVector hcencoordinate = ROTA * hlabcoordinate;
                              TVector3 hcencoordinate3D = hcencoordinate.Vect();

                              double TimeMax = TMath::Max(tcencoordinate.T(),
                                                          TMath::Max(wcencoordinate.T(),
                                                                     hcencoordinate.T()));
                              tcencoordinate3D = tcenmomentum3D * (1./tcenmomentum.E())
                                * (TimeMax-tcencoordinate.T()) + tcencoordinate3D;
                              wcencoordinate3D = wcenmomentum3D * (1./wcenmomentum.E())
                                * (TimeMax-wcencoordinate.T()) + wcencoordinate3D;
                              hcencoordinate3D = hcenmomentum3D * (1./hcenmomentum.E())
                                * (TimeMax-hcencoordinate.T()) + hcencoordinate3D;
                              double rhoSquare = (sqrt(1./2.) * (tcencoordinate3D-wcencoordinate3D)).Mag2();
                              rhoSquare = 25.6889 * rhoSquare;//GeV
                              double lambdaSquare = (sqrt(1./6.)
                                                     * (tcencoordinate3D+wcencoordinate3D-2.*hcencoordinate3D)).Mag2();
                              lambdaSquare = 25.6889 * lambdaSquare; //GeV
                              double krhoSquare = (sqrt(1./2.) * (tcenmomentum3D-wcenmomentum3D)).Mag2();
                              double klambdaSquare = (sqrt(1./6.)
                                                      * (tcenmomentum3D+wcenmomentum3D-2.*hcenmomentum3D)).Mag2();
                              double omegaSquare = RRMS * RRMS;

                              cout << rhoSquare/omegaSquare << endl;
                              cout << lambdaSquare/omegaSquare << endl;
                              cout << krhoSquare*omegaSquare << endl;
                              cout << klambdaSquare*omegaSquare << endl;
                              double Fb = G * 8 * 8 * exp(-(rhoSquare+lambdaSquare)/omegaSquare
                                                          - (krhoSquare+klambdaSquare)*omegaSquare);

                              /// calculate flow
                              double pxval = (tlabmomentum+wlabmomentum+hlabmomentum).Px();
                              double pyval = (tlabmomentum+wlabmomentum+hlabmomentum).Py();
                              double hadron_phi = atan2(pyval, pxval);
                              double V2 = cos(2.*(hadron_phi-PHI2));
                              double V3 = cos(3.*(hadron_phi-PHI3));


                              rapidityhist->Fill(rapidity, Fb);
                              pthist->Fill(transpt, Fb);
                              if (abs(rapidity) > dy)
                                {
                                  continue;
                                }
                              EventNum += Fb;
                              
                              v2SumHistgram->Fill(transpt, V2*Fb);
                              v3SumHistgram->Fill(transpt, V3*Fb);
                              probablySum->Fill(transpt,   Fb);

                              
                              if (transpt > 1e-6)
                                {
                                  yield->Fill(transpt, Fb*1./(2.*PI*transpt*dpt*2*dy));
                                }
                            }//if
                        }// for
                    }//if
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
