/**
 * @file   cumulant.cpp
 * @author xiaohai <xiaohaijin@outlook.com>
 * @date   Thu May  4 12:55:25 2017
 * 
 * @brief  Code to test the implementation of 
 *         two and four-particle cumulants on
 *         output from the AMPT generator.
 *         please read the paper : 
 *         FLow analysis with Q-cumulants.
 */
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
#include "Head.h"

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
      for (unsigned t = 0; t != tracks; ++t)
        {
          double energy = sqrt(pow(ampt->Px[t], 2) + pow(ampt->Py[t], 2)
                               + pow(ampt->Pz[t], 2) + pow(ampt->Mass[t], 2));
          TLorentzVector lor(ampt->Px[t], ampt->Py[t], ampt->Pz[t], energy);
          if (lor.Pt() < 0.001) continue;
          if (abs(lor.PseudoRapidity()) > 2) continue;
          int particleId = abs(ampt->ID[t]);
          if (!(particleId == 211 || particleId == 321 || particleId == 2212)) continue;
          xiaohai::Particle p;
          p.id = ampt->ID[t];
          p.px = ampt->Px[t];
          p.py = ampt->Py[t];
          p.pz = ampt->Pz[t];
          p.x = ampt->X[t];
          p.y = ampt->Y[t];
          p.z = ampt->Z[t];
          p.phi = lor.Phi();
          p.rsquare = pow(ampt->X[t], 2) + pow(ampt->Y[t], 2);
          p.pT = lor.Pt();
          p.eta = lor.PseudoRapidity();

          finalparticles.push_back(p);
        }
      if (finalparticles.size() > 50)
        {
          processEvent();
        }
      finalparticles.clear();
    }

  computeFlow();
  output_data.close();
  delete chain;
  delete ampt;
  delete saveFile;
  return 0;
}



void processEvent()
{
  int mult = finalparticles.size();
  int mult_A = 0;
  int mult_B = 0;

  // --- first track loop, q-vectors
  float Q2x = 0;
  float Q2y = 0;
  float Q4x = 0;
  float Q4y = 0;

  float Q2x_A = 0;
  float Q2y_A = 0;
  float Q4x_A = 0;
  float Q4y_A = 0;

  float Q2x_B = 0;
  float Q2y_B = 0;
  float Q4x_B = 0;
  float Q4y_B = 0;

  for (int itrk = 0; itrk < mult; itrk++)
    {
      float phi = finalparticles[itrk].phi;
      float eta = finalparticles[itrk].eta;

      /// 总的粒子的
      Q2x += cos(2 * phi);
      Q2y += sin(2 * phi);
      Q4x += cos(4 * phi);
      Q4y += sin(4 * phi);

      /// 后向快度的
      if (eta < -1 * etaGap)
        {
          Q2x_A += cos(2 * phi);
          Q2y_A += sin(2 * phi);
          Q4x_A += cos(4 * phi);
          Q4y_A += sin(4 * phi);

          mult_A++;
        }
      /// 前向快度的
      else if (eta > etaGap)
        {
          Q2x_B += cos(2 * phi);
          Q2y_B += sin(2 * phi);
          Q4x_B += cos(4 * phi);
          Q4y_B += sin(4 * phi);

          mult_B++;
        }

    } // End of track loop

  /// 公式(29)
  float two = ( Q2x * Q2x + Q2y * Q2y - mult) / (mult * mult - mult);
  float two_etagap = (Q2x_A * Q2x_B + Q2y_A * Q2y_B) / (mult_A * mult_B);
  /// 公式(32)
  float four = calc4_event(Q2x, Q2y, Q4x, Q4y, mult);

  h_db_2->Fill(0.5, two);
  h_db_4->Fill(0.5, four);
  h_db_2_etagap->Fill(0.5, two_etagap);

  /// 二次径迹循环。
  for (int itrk = 0; itrk < mult; itrk++)
    {
      xiaohai::Particle p = finalparticles[itrk];
      float phi = p.phi;
      float pt = p.pT;
      float eta = p.eta;

      float u2x = cos(2 * phi);
      float u2y = sin(2 * phi);

      float twoprime = ( u2x * Q2x + u2y * Q2y - 1) / (mult - 1);
      h_db_2prime->Fill(pt, twoprime);

      float u4x = cos(4 * phi);
      float u4y = sin(4 * phi);

      float fourprime = calc4_track(u2x, u2y, u4x, u4y, Q2x, Q2y, Q4x, Q4y, mult);
      h_db_4prime->Fill(pt, fourprime);

      if (eta < -1 * etaGap)
        {
          float u2x_A = cos(2 * phi);
          float u2y_A = sin(2 * phi);

          float twoprime_etagap = (u2x_A * Q2x_B + u2y_A * Q2y_B) / mult_B;
          h_db_2prime_etagap->Fill(pt, twoprime_etagap);
        }

    } // End of track loop
}

/// formular (32)
float calc4_event(float Xn, float Yn, float X2n, float Y2n, float M)
{

  float Qn2   = Xn * Xn + Yn * Yn;
  float Qn2d  = Xn * Xn - Yn * Yn;

  float one   = Qn2 * Qn2;
  float two   = X2n * X2n + Y2n * Y2n;
  float three = (2 * (X2n * Qn2d + 2 * Y2n * Xn * Yn));
  float four  = 2 * (2 * (M - 2) * Qn2);
  float five  = 2 * M * (M - 3);

  float numerator = one + two - three - four + five;
  float denominator = M * (M - 1) * (M - 2) * (M - 3);

  return numerator / denominator;

}

float calc4_track(float xn, float yn, float x2n, float y2n, float Xn,
                  float Yn, float X2n, float Y2n, float M)
{

  float one   = (xn * Xn + yn * Yn) * (Xn * Xn + Yn * Yn);
  float two   = x2n * Xn * Xn - x2n * Yn * Yn + 2 * y2n * Xn * Yn;
  float three = xn * Xn * X2n + xn * Yn * Y2n - yn * (X2n * Yn - Xn * Y2n);
  float four  = 2 * M * (xn * Xn + yn * Yn);
  float five  = 2 * (Xn * Xn + Yn * Yn);
  float six   = 7 * (xn * Xn + yn * Yn);
  float seven = xn * Xn + yn * Yn;
  float eight = x2n * X2n + y2n * Y2n;
  float nine = 2 * (xn * Xn + yn * Yn);

  float numerator = one - two - three - four - five + six - seven + eight + nine + 2 * M - 6;
  float denominator = (M - 1) * (M - 2) * (M - 3);

  return numerator / denominator;

}

void computeFlow()
{
  float c2_2        = h_db_2->GetBinContent(1);
  float c2_2_etagap = h_db_2_etagap->GetBinContent(1);
  float c2_4        = h_db_4->GetBinContent(1) - 2 * pow(h_db_2->GetBinContent(1), 2);

  cout << "c2_2    " << c2_2 << endl;
  cout << "c2_2 eg " << c2_2_etagap << endl;
  cout << "c2_4    " << c2_4 << endl;

  TProfile *hd2_2 = (TProfile*) h_db_2prime->Clone("hd2_2");
  hd2_2->SetTitle("Two Particle Cumulants");
  hd2_2->Scale(1.0 / sqrt(c2_2));

  TProfile *hd2_2_etagap = (TProfile*) h_db_2prime_etagap->Clone("hd2_2_etagap");
  hd2_2_etagap->SetTitle("Two Particle Cumulanta with Eta Gap");
  hd2_2_etagap->Scale(1.0/sqrt(c2_2_etagap));

  TH1D *h_proj_db_4prime = h_db_4prime->ProjectionX("h_proj_db_4prime", "E");
  TH1D *h_proj_db_2prime = h_db_2prime->ProjectionX("h_proj_db_2prime", "E");

  h_proj_db_2prime->Scale(2 * h_db_2->GetBinContent(1));
  TH1D *hd2_4 = (TH1D*) h_proj_db_4prime->Clone("hd2_4");
  hd2_4->Add(h_proj_db_2prime, -1);

  hd2_4->Scale(-1.0 / pow(-1 * c2_4, 0.75));

  hd2_2_etagap->Draw();
  hd2_2->SetLineColor(kRed);
  hd2_2->Draw("same");
}

void writeHistosToFile(char *outFileName = "")
{
  TFile* fout = new TFile(outFileName, "RECREATE");
  h_db_2prime->Write();
  h_db_4prime->Write();
  h_db_2->Write();
  h_db_4->Write();
  fout->Close();
}
