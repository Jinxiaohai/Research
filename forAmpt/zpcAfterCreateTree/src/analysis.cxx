//--------------------------------
// Load header                              
//--------------------------------

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "math.h"
#include "string.h"

#include <TROOT.h>
#include <TFile.h>
#include "TObject.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMath.h"
#include <TLorentzVector.h>
#include "TVector2.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TComplex.h"

using namespace std;

#endif

#include "AMPT.h"

//--------------------------------
// Declare functions                                
//--------------------------------

int GetCentrality(AMPT* ampt);
int GetCharge(int id);
bool IsRFP(float eta, float pt, int loop);
bool IsPOI(float eta, float pt, int id, int loop);

//--------------------------------
// Main function starts here                               
//--------------------------------

int main(int argc, char** argv)
{
  if (argc!=3) return 1;

  //   char* inputFile;
  //   char* outputFile;
  TString inputFile;
  TString outputFile;
  if (argc==3) {
    inputFile = argv[1];
    outputFile = argv[2];
  }  

  //----------------------------------
  // Open files and add to chain
  //----------------------------------

  int fileNumber = 0;
  char fileList[512];
  TChain* chain = new TChain("AMPT");
  if (inputFile.Contains(".list"))  {
    ifstream* inputStream = new ifstream;
    inputStream->open(inputFile);
    if (!(inputStream)) {
      cout<<"can not open file list"<<endl;
      return 1;
    }
    while (inputStream->good()) {
      inputStream->getline(fileList, 512);
      if (inputStream->eof()) break;
      TFile *fTmp = new TFile(fileList);
      
      if (!fTmp || !(fTmp->IsOpen()) || !(fTmp->GetNkeys())) {
	cout<<"open file list error"<<endl;
	return 1;
      } else {
	cout<<"reading file "<<fileList<<endl;
	chain->Add(fileList);
	fileNumber++;
      }
      delete fTmp;
    }
    cout<<fileNumber<<" files read in"<<endl;
  } else if (inputFile.Contains(".root")) {
    chain->Add(inputFile.Data());
  }

  //--------------------------------
  // define histograms
  //--------------------------------

  TH1F* h1_centrality = new TH1F("centrality","",10,0,10);
  TH1F* h1_eta = new TH1F("eta","",200,-10,10);
  TH1F* h1_phi = new TH1F("phi","",140,-7,7);
  TProfile* p1_subEvent_refCorr = new TProfile("subEventRefCorr","",5,0,5,"s");
  TH1F* h1_subEvent_sumW2 = new TH1F("subEventSumW2","",5,0,5);
  TH3F* h3_mP[2] = {NULL};
  TProfile* p1_refCorrelations[2] = {NULL};
  TH1F* h1_refBinSumW2[2] = {NULL};
  TProfile3D* p3_C14[2] = {NULL};
  TProfile3D* p3_C15[2] = {NULL};
  TProfile3D* p3_T28[2] = {NULL};
  TH3F* h3_DiffBinSumW12[2] = {NULL};
  TH3F* h3_DiffBinSumW22[2] = {NULL};
  TProfile3D* p3_T16abT28[2] = {NULL};  
  for (Int_t i = 0; i < 2; ++i) {
    h3_mP[i] = new TH3F((TString("yield_sub")+=i).Data(), "", 5, 0, 5, 40, 0, 2, 20, -1, 1);
    p1_refCorrelations[i] = new TProfile((TString("refCorrelations_sub")+=i).Data(), "", 8, 0, 8, "s");
    h1_refBinSumW2[i] = new TH1F((TString("refBinSumW2_sub")+=i).Data(), "", 8, 0, 8);
    p3_C14[i] = new TProfile3D((TString("C14_sub")+=i).Data(), "", 5, 0, 5, 40, 0, 2, 20, -1, 1, "s");
    p3_C15[i] = new TProfile3D((TString("C15_sub")+=i).Data(), "", 5, 0, 5, 40, 0, 2, 20, -1, 1, "s");
    p3_T28[i] = new TProfile3D((TString("T28_sub")+=i).Data(), "", 5, 0, 5, 40, 0, 2, 20, -1, 1, "s");
    h3_DiffBinSumW12[i]= new TH3F((TString("diffBinSumW12_sub")+=i).Data(), "", 5, 0, 5, 40, 0, 2, 20, -1, 1);
    h3_DiffBinSumW22[i]= new TH3F((TString("diffBinSumW22_sub")+=i).Data(), "", 5, 0, 5, 40, 0, 2, 20, -1, 1);
    p3_T16abT28[i]= new TProfile3D((TString("T16abT28_sub")+=i).Data(), "", 5, 0, 5, 40, 0, 2, 20, -1, 1);
  }
  
  //--------------------------------
  // loop events
  //--------------------------------

  int nEvents = (int)chain->GetEntries();
  cout<<"Total events : "<<nEvents<<endl;
  AMPT* ampt = new AMPT(chain);
  gRandom = new TRandom3(0);
  Int_t harmonics = 2;
  
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (iEvent%1000==0) cout << "Processing event # " << iEvent << endl;
    ampt->GetEntry(iEvent);

    int nTracks = ampt->Event_multiplicity;
    int cent = GetCentrality(ampt);
    if (cent<0) continue;
    h1_centrality->Fill(cent);

    TH3F* h3_reP = new TH3F("Rep2", "", 5, 0, 5, 40, 0, 2, 20, -1, 1);
    TH3F* h3_imP = new TH3F("Imp2", "", 5, 0, 5, 40, 0, 2, 20, -1, 1);
    TH3F* h3_mP_this_event = new TH3F("mp", "", 5, 0, 5, 40, 0, 2, 20, -1, 1);

    //--------------------------------
    // 1st loop tracks                               
    //--------------------------------

    TComplex q1(0,0), q2(0,0);
    Int_t m1 = 0, m2 = 0;
    for (int iTrack = 0; iTrack < nTracks; ++iTrack) {
      int id = ampt->id[iTrack];
      int charge = GetCharge(id);
      if (charge==0) continue;
      float px = ampt->px[iTrack];
      float py = ampt->py[iTrack];
      float pz = ampt->pz[iTrack];
      TVector3 p(px,py,pz);
      float pt = p.Pt();
      if (p.Mag()<1.e-5 || pt<1.e-5) continue;
      float eta = p.PseudoRapidity();
      h1_eta->Fill(eta);
      float phi = p.Phi();
      h1_phi->Fill(phi);
      
//       if (fabs(eta)>1) continue;
//       if (fabs(eta)<0.3 || fabs(eta)>1) continue;
//       if (pt>2) continue; // 0.15

      TComplex q(TMath::Cos(harmonics*phi), TMath::Sin(harmonics*phi));
      if (eta < 0) {
	q1 += q;
	m1++;
      } else if (eta > 0) {
	q2 += q;
	m2++;
      }
    }

    if (!m1 || !m2) continue;
    p1_subEvent_refCorr->Fill((double)0, (q1*TComplex::Conjugate(q2)).Re() / (m1*m2), m1*m2);
    p1_subEvent_refCorr->Fill(1, q1.Re()/m1, m1);
    p1_subEvent_refCorr->Fill(2, q2.Re()/m2, m2);
    p1_subEvent_refCorr->Fill(3, q1.Im()/m1, m1);
    p1_subEvent_refCorr->Fill(4, q2.Im()/m2, m2);
    h1_subEvent_sumW2->Fill((double)0, m1*m2*m1*m2);
    h1_subEvent_sumW2->Fill(1, m1*m1);
    h1_subEvent_sumW2->Fill(2, m2*m2);
    h1_subEvent_sumW2->Fill(3, m1*m1);
    h1_subEvent_sumW2->Fill(4, m2*m2);
    
    //--------------------------------
    // 2nd loop tracks                               
    //--------------------------------

    for (int iLoop = 0; iLoop < 2; ++iLoop) {
      Double_t sum_cos_nPhi = 0;
      Double_t sum_sin_nPhi = 0;
      Int_t nRFP = 0;

      for (int iTrack = 0; iTrack < nTracks; ++iTrack) {
	int id = ampt->id[iTrack];
	int charge = GetCharge(id);
	if (charge==0) continue;
	float px = ampt->px[iTrack];
	float py = ampt->py[iTrack];
	float pz = ampt->pz[iTrack];
	TVector3 p(px,py,pz);
	float pt = p.Pt();
	if(p.Mag()<1.e-5 || pt<1.e-5) continue;
	float eta = p.PseudoRapidity();
// 	if (eta>1) continue;
	float phi = p.Phi();

// 	if (pt>2) continue; // 0.15

	double cos_nPhi = TMath::Cos(harmonics*phi);
	double sin_nPhi = TMath::Sin(harmonics*phi);

	//--------------------------------
	// 1st loop: η>0.3~RFP, η<0~POI
	// 2nd loop: η<-0.3~RFP, η>0~POI
	//--------------------------------

	bool isRFP = IsRFP(eta, pt, iLoop);
	bool isPOI = IsPOI(eta, pt, id, iLoop);

	if (isRFP) {
	  nRFP++;
          sum_cos_nPhi += cos_nPhi;
	  sum_sin_nPhi += sin_nPhi;
	} else if (isPOI) {
	  h3_reP->Fill(1, pt, eta, cos_nPhi);
	  h3_imP->Fill(1, pt, eta, sin_nPhi);
	  h3_mP_this_event->Fill(1, pt, eta);
	  h3_mP[iLoop]->Fill(1, pt, eta);
	}
      }

      if (nRFP < 2) continue;

      //--------------------------------
      // cumulants
      //--------------------------------

      Double_t M       = nRFP;
      Double_t nComb1  = M;
      Double_t nComb2  = M*(M-1);
      Double_t weight1 = nComb1;
      Double_t weight2 = nComb2;
      // Double_t weight1 = 1;
      // Double_t weight2 = 1;

      TComplex Q2(sum_cos_nPhi, sum_sin_nPhi);
      TComplex Q2Star = TComplex::Conjugate(Q2);
      Double_t Q2Square = Q2.Rho2();
      Double_t coor22   = (Q2Square - M);

      p1_refCorrelations[iLoop]->Fill((double)0,  Q2.Re()/nComb1,    weight1);
      p1_refCorrelations[iLoop]->Fill(1,  Q2.Im()/nComb1,    weight1);
      p1_refCorrelations[iLoop]->Fill(2, coor22/nComb2,     weight2);
      h1_refBinSumW2[iLoop]->Fill((double)0,  weight1*weight1);
      h1_refBinSumW2[iLoop]->Fill(1,  weight1*weight1);
      h1_refBinSumW2[iLoop]->Fill(2, weight2*weight2);
      
      Int_t etaBin[4] = {1,10,11,20}; 
      for (Int_t iPOI = 1; iPOI <= 5; ++iPOI) {
	for (Int_t iPt = 1; iPt <= 40; ++iPt) {
	  for (Int_t iEta = etaBin[iLoop]; iEta <= etaBin[iLoop+2]; ++iEta) {

	    Double_t binPOI = h3_reP->GetXaxis()->GetBinCenter(iPOI);
	    Double_t binPt  = h3_reP->GetYaxis()->GetBinCenter(iPt);
	    Double_t binEta = h3_reP->GetZaxis()->GetBinCenter(iEta);

	    TComplex p2(h3_reP->GetBinContent(iPOI, iPt, iEta),
			h3_imP->GetBinContent(iPOI, iPt, iEta));
	    Double_t mp = h3_mP_this_event->GetBinContent(iPOI, iPt, iEta);

	    Double_t diffNComb1 = mp;
	    Double_t diffNComb2 = mp*M;
	    Double_t diffWeight1 = diffNComb1;
	    Double_t diffWeight2 = diffNComb2;
	    if (diffNComb1 == 0 || diffNComb2 == 0) continue;
	    Double_t diffCoor22  = (p2*Q2Star).Re(); // <2'>

	    p3_C14[iLoop]->Fill(binPOI, binPt, binEta, p2.Re()/diffNComb1,    diffWeight1);
	    p3_C15[iLoop]->Fill(binPOI, binPt, binEta, p2.Im()/diffNComb1,    diffWeight1);
	    p3_T28[iLoop]->Fill(binPOI, binPt, binEta, diffCoor22/diffNComb2, diffWeight2);
	    h3_DiffBinSumW12[iLoop]->Fill(binPOI, binPt, binEta, diffWeight1*diffWeight1);
	    h3_DiffBinSumW22[iLoop]->Fill(binPOI, binPt, binEta, diffWeight2*diffWeight2);
	    p3_T16abT28[iLoop]->Fill(binPOI, binPt, binEta,
				     ((q1*TComplex::Conjugate(q2)).Re()/(m1*m2))*(diffCoor22/diffNComb2),
				     (m1*m2)*diffWeight2);
	  }
	}
      }
    }

    h3_reP->Delete();
    h3_imP->Delete();
    h3_mP_this_event->Delete();
  }

  //--------------------------------
  // Write Hists and exit                               
  //--------------------------------

  TFile* f = new TFile(outputFile, "RECREATE");
  f->cd();

  h1_eta->Write();
  h1_phi->Write();
  h1_centrality->Write();
    
  p1_subEvent_refCorr->Write();
  h1_subEvent_sumW2->Write();
  for (Int_t i = 0; i < 2; ++i) {
    h3_mP[i]->Write();
    p1_refCorrelations[i]->Write();
    h1_refBinSumW2[i]->Write();
    p3_C14[i]->Write();
    p3_C15[i]->Write();
    p3_T28[i]->Write();
    h3_DiffBinSumW12[i]->Write();
    h3_DiffBinSumW22[i]->Write();
    p3_T16abT28[i]->Write();
  }

  f->Close();

  delete chain;
  return 0;
}

//--------------------------------
// Define functions here                               
//--------------------------------

int GetCentrality(AMPT* ampt) // different with StRefMultCorr
{
  float b = ampt->Event_b;

  int cent=-1;
  for (int i=1 ; i<=10 ; ++i) {
    if (b<=sqrt(10*i/100.0)*2.0*pow(197, 1.0/3.0)*1.2) { // ...*1.124
      cent = i-1;
      return cent;
    }
  }

  return cent;
}

int GetCharge(int id)
{
  int charge = 0;

  if (id==11) charge=-1;   // e-
  else if (id==-11) charge=1;   // e+
  else if (id==211) charge=1;   // π+
  else if (id==-211) charge=-1;   // π-
  else if (id==213) charge=1;   // ρ(770)+
  else if (id==321) charge=1;   // Κ+
  else if (id==-321) charge=1;   // Κ-
  else if (id==323) charge=1;   // Κ*(892)+
  else if (id==411) charge=1;   // D+
  else if (id==413) charge=1;   // D*(2010)+
  else if (id==431) charge=1;   // Ds+
  else if (id==433) charge=1;   // Ds*+
  else if (id==521) charge=1;   // B+
  else if (id==523) charge=1;   // B*+
  else if (id==1114) charge=-1; // Δ-
  else if (id==2212) charge=1;  // p+
  else if (id==-2212) charge=1;  // p-
  else if (id==2214) charge=1;  // Δ+
  else if (id==2224) charge=2;  // Δ++
  else if (id==3112) charge=-1; // Σ-
  else if (id==3114) charge=-1; // Σ*-
  else if (id==3222) charge=1;  // Σ+
  else if (id==3224) charge=1;  // Σ*+
  else if (id==3312) charge=-1; // Ξ-
  else if (id==3314) charge=-1; // Ξ*-
  else if (id==3334) charge=-1; // Ω-
  else if (id==4122) charge=1;  // Λc+
  else if (id==4222) charge=2;  // Σc++
  else if (id==4212) charge=1;  // Σc+
  else if (id==4224) charge=2;  // Σc*++
  else if (id==4232) charge=1;  // Ξc+
  else if (id==4322) charge=1;  // Ξ'c+
  else if (id==4324) charge=1;  // Ξc*+

  return charge;
}

bool IsRFP(float eta, float pt, int loop)
{
//   if (loop==0 && eta>0.3 && pt<2) return true;
//   else if (loop==1 && eta<-0.3 && pt<2) return true;
  if (loop==0 && (eta>0.3 || eta<-1.3)) return true;
  else if (loop==1 && (eta<-0.3 || eta>1.3)) return true;
  else return false;
}

bool IsPOI(float eta, float pt, int id, int loop)
{
  if (loop==0 && eta<0 && eta>-1 && id==211) return true;
  else if (loop==1 && eta>0 && eta<1 && id==211) return true;
//   if (loop==0 && eta<0 && pt<2) return true;
//   else if (loop==1 && eta>0 && pt<2) return true;
  else return false;
}
