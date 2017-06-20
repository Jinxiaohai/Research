#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "math.h"
#include "string.h"

#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"

#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TRandom3.h"
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "iomanip"
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
//#include "MYEvent.h"

double PI = TMath::Pi();
double ANUMBER = 197.;
//double pTbin[15]={0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.5,4.0,5.0};

TH1D *proton_pT_dst = new TH1D("proton_pT_dst","proton_pT_dst",50,0,5.0);
TProfile *proton_v2_pT  = new TProfile("proton_v2_pT","proton_v2_pT",30,0,3.0,-1.,1.);
TProfile *proton_v3_pT  = new TProfile("proton_v3_pT","proton_v3_pT",30,0,3.0,-1.,1.);
TProfile *proton_v3_pT_etap  = new TProfile("proton_v3_pT_etap","proton_v3_pT_etap",30,0,3.0,-1.,1.);
TProfile *proton_v1_eta = new TProfile("proton_v1_eta","proton_v1_eta",60,-3.0,3.0,-1.,1.);
TProfile *proton_v1_eta_test = new TProfile("proton_v1_eta_test","proton_v1_eta_test",60,-3.0,3.0,-1.,1.);

TH1D *pion_pT_dst = new TH1D("pion_pT_dst","pion_pT_dst",50,0.,5.0);
TProfile *pion_v2_pT = new TProfile("pion_v2_pT","pion_v2_pT",30,0,3.0,-1.,1.);
TProfile *pion_v3_pT = new TProfile("pion_v3_pT","pion_v3_pT",30,0,3.0,-1.,1.);
TProfile *pion_v3_pT_etap = new TProfile("pion_v3_pT_etap","pion_v3_pT_etap",30,0,3.0,-1.,1.);
TProfile *pion_v1_eta = new TProfile("pion_v1_eta","pion_v1_eta",60,-3.0,3.0,-1.,1.0);

TH1D *ch_pT_dst = new TH1D("ch_pT_dst","ch_pT_dst",50,0.,5.0);
TProfile *ch_v2_pT = new TProfile("ch_v2_pT","ch_v2_pT",30,0.,3.0,-1.,1.);
TProfile *ch_v3_pT = new TProfile("ch_v3_pT","ch_v3_pT",30,0.,3.0,-1.,1.);
TProfile *ch_v3_pT_etap = new TProfile("ch_v3_pT_etap","ch_v3_pT_etap",30,0.,3.0,-1.,1.);
TProfile *ch_v1_eta = new TProfile("ch_v1_eta","ch_v1_eta",30,0.,3.0,-1.,1.);

double EllipticFlow(TVector3 PP);
int Centrality( double impact_parameter);

int main(int argc, char **argv)
{
  char *FileInput=0;
  char *FileOutput=0;

  if(argc!=3 && argc!=1) return 0;

  if(argc==1)
    {
      FileInput = "example.list";
      FileOutput = "example.root";
    }

  if(argc==3)
    {
      FileInput = argv[1];
      FileOutput = argv[2];
    }

  Int_t nevent,nrun,multi,NpartP,NpartT,NELP,NINP,NELT,NINT;
  Float_t impactpar;
  Float_t passhead;
  Float_t px,py,pz,am,x,y,z,time;
  Float_t energy;

  char others1[256],others2[256],others3[256],others4[256];
  Int_t id;

  char outfile[256];
  sprintf(outfile,"%s.root",FileOutput); 
  TFile *file = new TFile(outfile,"RECREATE");


  //read in data
  char FileList[512];
  Int_t count=0;
  Int_t filenum=0;
  ifstream* inputStream = new ifstream;
  inputStream->open(FileInput);
  if (!(inputStream))
    {
      printf("can not open list file\n");
      return 0;
    }
  for(;inputStream->good();){
    inputStream->getline(FileList,512);
    //TFile *ftmp = new TFile(FileList);
    //if(!ftmp)
    if  ( inputStream->good() ){
      printf(" read in file %s\n",FileList);
      //filenum++; 
      ifstream filein;
      filein.open(FileList);

      while(filein.good()){
        filein>>nevent>>nrun>>multi>>impactpar>>NpartP>>NpartT>>NELP>>NINP>>NELT>>NINT>>passhead;
        if(!filein.good()) break;

        for (Int_t nlines=0;nlines<multi;nlines++) {
          //filein >> id >> px >> py >> pz >> am >> x >> y >> z >> time;
          filein >> id >> px >> py >> pz >> am >> others1 >> others2 >> others3 >> others4;
          if(!filein.good()) break;

          if(strstr(others1,"*") || strstr(others2,"*") || strstr(others3,"*") ||strstr(others4,"*"))
            {
              x=1000000.;
              y=1000000.;
              z=1000000.;
              time=1000000.;
            }
          else
            {
              x=atof(others1);
              y=atof(others2);
              z=atof(others3);
              time=atof(others4);
            }

          //TVector3 mP(px,py,pz);
          energy = sqrt(px*px+py*py+pz*pz+am*am);
          TLorentzVector mP(px,py,pz,energy);
          double pseorap=10.;
          if(mP.Pt()<10e-7)continue;
          pseorap=mP.PseudoRapidity();
          double phiangle = mP.Phi();
          //int centrality = Centrality(impactpar);
          //if(centrality>=0 && centrality<80 && id==2212)
          if(fabs(id)==211){//pip, pim
            pion_v1_eta->Fill(pseorap,cos(phiangle));
            if(pseorap>0.0)pion_v3_pT_etap->Fill(mP.Pt(),cos(3.*phiangle));
            if(fabs(pseorap)<1.0){
              pion_pT_dst->Fill(mP.Pt());
              pion_v2_pT->Fill(mP.Pt(),cos(2.*phiangle));
              pion_v3_pT->Fill(mP.Pt(),cos(3.*phiangle));
            }
          }
          if(id==2212)
            {
              proton_v1_eta->Fill(pseorap,cos(phiangle));
              proton_v1_eta_test->Fill(pseorap,px/mP.Pt());
              if(pseorap>0.0)proton_v3_pT_etap->Fill(mP.Pt(),cos(3.*phiangle));

              if(fabs(pseorap)<1.0){
                proton_pT_dst->Fill(mP.Pt());
                proton_v2_pT->Fill(mP.Pt(),cos(2.*phiangle));
                proton_v3_pT->Fill(mP.Pt(),cos(3.*phiangle));
              }
            }
          if(fabs(id)==211||fabs(id)==321||fabs(id)==2212){//pip,pim,kp,pm,p,anti-p
            ch_v1_eta->Fill(pseorap,cos(phiangle));
            if(pseorap>0.0)ch_v3_pT_etap->Fill(mP.Pt(),cos(3*phiangle));
            if(fabs(pseorap)<1.0){
              ch_pT_dst->Fill(mP.Pt());
              ch_v2_pT->Fill(mP.Pt(),cos(2.*phiangle));
              ch_v2_pT->Fill(mP.Pt(),cos(3.*phiangle));
            }
          }
        }
        count++;
      }
      filein.close();
      filenum++;
    }
  }
  inputStream->close();
  cout<< "read in "<<filenum <<" good files w/ "<<count<<" events"<<endl;
  pion_v1_eta->Write();
  pion_pT_dst->Write();
  pion_v2_pT->Write();
  pion_v3_pT->Write();
  pion_v3_pT_etap->Write();

  proton_v1_eta->Write();
  proton_v1_eta_test->Write();
  proton_v3_pT_etap->Write();
  proton_pT_dst->Write();
  proton_v2_pT->Write();
  proton_v3_pT->Write();
	
  ch_v1_eta->Write();
  ch_pT_dst->Write();
  ch_v2_pT->Write();
  ch_v3_pT->Write();
  ch_v3_pT_etap->Write();

  file->Write();

  delete pion_v1_eta;
  delete pion_pT_dst;
  delete pion_v2_pT;
  delete pion_v3_pT;
  delete pion_v3_pT_etap;
  delete proton_v1_eta;
  delete proton_v1_eta_test;
  delete proton_v3_pT_etap;
  delete proton_pT_dst;
  delete proton_v2_pT;
  delete proton_v3_pT;
  delete ch_v1_eta;
  delete ch_v2_pT;
  delete ch_v3_pT;
  delete ch_v3_pT_etap;
  delete ch_pT_dst;
  delete file;

}
//-----------
int Centrality( double impact_parameter)
{
  int central=1000;
  for(int i=1 ; i<=100 ; i++)
    {
      if(impact_parameter<=sqrt(i/100.)*2.*pow(ANUMBER,1./3.)*1.124)
        {
          central = i-1;
          return central;
        }
    }
  return central;
}
//-------elliptic flow-----
double EllipticFlow(TVector3 PP)
{
  if((PP.Px()==0.)&&(PP.Py()!=0.)) {return -1.;}
  else if((PP.Px()!=0.)&&(PP.Py()==0.)) {return 1.;}
  else if((PP.Px()==0.)&&(PP.Py()==0.)) {return 0.;}
  else {return (PP.Px()*PP.Px()-PP.Py()*PP.Py())/PP.Perp2();}
}
