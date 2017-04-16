#include "HIJING_file_parser.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include "ntupler/nt.C"
#include "ntupler/hi.C"


struct newEvent{
public :
  static const int MAXNUCLEON = 1000;   // Declaration of leaf types
   int            event;
   float          b;
   float          beta;
   int            nTrk;
   float          frame_E;
   float          total_E;
   int            Ntrkpos;
   int            Ntrkneg;
   int            boostedNtrkpos;
   int            boostedNtrkneg;
   float          ETpos;
   float          ETneg;
   float          boostedETpos;
   float          boostedETneg;
   float          x;
   float          trkPt[20000]; 
   float          trkPhi[20000]; 
   float          trkEta[20000]; 
   float          trkE[20000];
   float          trkPx[20000];
   float          trkPy[20000];
   float          trkPz[20000];
   float          trkRapidity[20000];
   float          boostedPt[20000]; 
   float          boostedPhi[20000]; 
   float          boostedEta[20000]; 
   float          boostedE[20000];
   float          boostedPx[20000];
   float          boostedPy[20000];
   float          boostedPz[20000];
   float          boostedRapidity[20000];

   float          x1;
   float          x2;
   float          pt1;
   float          pt2;
   float          eta1;
   float          eta2;
   float          phi1;
   float          phi2;
   float          dphi;
   float          etadijet;

   void reset(){
    event=-9;
    frame_E=-9;
    total_E=-9;
    Ntrkpos=-9;
    Ntrkneg=-9;
    boostedNtrkpos=-9;
    boostedNtrkneg=-9;
    ETpos=-9;
    ETneg=-9;
    boostedETpos=-9;
    boostedETneg=-9;
    b=-9;
    beta=-9;
    x=-9;
    nTrk = -9;
    x1=-9;
    x2=-9;
    pt1=-9;
    pt2=-9;
    phi1 = -9;
    phi2 = -9;
    eta1 = -9;
    eta2 = -9;
    etadijet = -9;
    dphi = -9;
   for(int i=0;i<20000;i++){
    trkPt[i] = -900;
    trkPhi[i] = -900;
    trkEta[i] = -900;
    trkE[i] = -900;
    trkPx[i] = -900;
    trkPy[i] = -900;
    trkPz[i] = -900;
    trkRapidity[i] = -900;
    boostedPt[i] = -900;
    boostedPhi[i] = -900;
    boostedEta[i] = -900;
    boostedE[i] = -900;
    boostedPx[i] = -900;
    boostedPy[i] = -900;
    boostedPz[i] = -900;
    boostedRapidity[i] = -900;
   }
   }
   newEvent(){
    reset();
   }
   void setEvent(int fevent, float fb, float fframe_E, float ftotal_E, float fx, float fbeta,float fETpos, float fETneg,float fboostedETpos, float fboostedETneg,int fntrkpos, int fntrkneg,int  fboostedntrkpos, int  fboostedntrkneg,int ntrk, float *trkpt, float *trkphi, float *trketa, float *trkpx, float *trkpy, float *trkpz, float *trke,float *trkrapidity ,float *boostedpt, float *boostedphi, float *boostedeta, float *boostedpx, float *boostedpy, float *boostedpz, float *boostede,float *boostedrapidity,float fpt1,float fpt2,float fphi1, float fphi2, float feta1, float feta2, float fetadijet, float fx1, float fx2){
   // void setEvent(int fevent, float fb, float fframe_E, float ftotal_E, float fx, float fbeta,float fETpos, float fETneg,float fboostedETpos, float fboostedETneg,int fntrkpos, int fntrkneg,int  fboostedntrkpos, int  fboostedntrkneg,int ntrk, float *trkpt, float *trkphi, float *trketa, float *trkpx, float *trkpy, float *trkpz, float *trke,float *trkrapidity ,float *boostedpt, float *boostedphi, float *boostedeta, float *boostedpx, float *boostedpy, float *boostedpz, float *boostede,float *boostedrapidity){
    event = fevent;
    b = fb;
    nTrk = ntrk; 
    x=fx;
    beta=fbeta;
    ETpos=fETpos;
    ETneg=fETneg;
    boostedETpos=fboostedETpos;
    boostedETneg=fboostedETneg;
    Ntrkpos=fntrkpos;
    Ntrkneg=fntrkneg;
    boostedNtrkpos=fboostedntrkpos;
    boostedNtrkneg=fboostedntrkneg;
    pt1=fpt1;
    pt2=fpt2;
    eta1=feta1;
    eta2=feta2;
    phi1=fphi1;
    phi2=fphi2;
    // dphi=fdphi;
    etadijet=fetadijet;
    x1=fx1;
    x2=fx2;
  for(int i=0; i<nTrk;i++){
   trkPt[i] = trkpt[i];
   trkPhi[i] = trkphi[i];
   trkEta[i] = trketa[i];
   trkPx[i] = trkpx[i];
   trkPy[i] = trkpy[i];
   trkPz[i] = trkpz[i];
   trkE[i] = trke[i];
   trkRapidity[i] = trkrapidity[i];
   boostedPt[i] = boostedpt[i];
   boostedPhi[i] = boostedphi[i];
   boostedEta[i] = boostedeta[i];
   boostedPx[i] = boostedpx[i];
   boostedPy[i] = boostedpy[i];
   boostedPz[i] = boostedpz[i];
   boostedE[i] = boostede[i];
   boostedRapidity[i] = boostedrapidity[i];
  }
 }
};

void make_Ntuple(){

 HIJING_file_parser * file = new HIJING_file_parser;
 vector <Event> event_vector = file->get_event_vector();
  
 newEvent evnt; 
 TFile * outfile = new TFile("out.root","recreate");//output file
    
 TTree * fT = new TTree("fT","a tree with HIJING events");
 fT->Branch("event",&evnt.event,"event/I");
 fT->Branch("b",&evnt.b,"b/F");
 fT->Branch("frame_E",&evnt.frame_E,"frame_E/F");
 fT->Branch("total_E",&evnt.total_E,"total_E/F");
 fT->Branch("x",&evnt.x,"x/F");
 fT->Branch("beta",&evnt.beta,"beta/F");
 fT->Branch("boostedETneg",&evnt.boostedETneg,"boostedETneg/F");
 fT->Branch("boostedETpos",&evnt.boostedETpos,"boostedETpos/F");
 fT->Branch("ETneg",&evnt.ETneg,"ETneg/F");
 fT->Branch("ETpos",&evnt.ETpos,"ETpos/F");
 fT->Branch("boostedNtrkneg",&evnt.boostedNtrkneg,"boostedNtrkneg/I");
 fT->Branch("boostedNtrkpos",&evnt.boostedNtrkpos,"boostedNtrkpos/I");
 fT->Branch("Ntrkneg",&evnt.Ntrkneg,"Ntrkneg/I");
 fT->Branch("Ntrkpos",&evnt.Ntrkpos,"Ntrkpos/I");
 fT->Branch("nTrk",&evnt.nTrk,"nTrk/I");
 fT->Branch("trkPt",evnt.trkPt,"trkPt[nTrk]/F");
 fT->Branch("trkPhi",evnt.trkPhi,"trkPhi[nTrk]/F");
 fT->Branch("trkEta",evnt.trkEta,"trkEta[nTrk]/F");
 fT->Branch("trkPx",evnt.trkPx,"trkPx[nTrk]/F");
 fT->Branch("trkPy",evnt.trkPy,"trkPy[nTrk]/F");
 fT->Branch("trkPz",evnt.trkPz,"trkPz[nTrk]/F");
 fT->Branch("trkE",evnt.trkE,"trkE[nTrk]/F");
 fT->Branch("boostedPt",evnt.boostedPt,"boostedPt[nTrk]/F");
 fT->Branch("boostedPhi",evnt.boostedPhi,"boostedPhi[nTrk]/F");
 fT->Branch("boostedEta",evnt.boostedEta,"boostedEta[nTrk]/F");
 fT->Branch("boostedPx",evnt.boostedPx,"boostedPx[nTrk]/F");
 fT->Branch("boostedPy",evnt.boostedPy,"boostedPy[nTrk]/F");
 fT->Branch("boostedPz",evnt.boostedPz,"boostedPz[nTrk]/F");
 fT->Branch("boostedE",evnt.boostedE,"boostedE[nTrk]/F");
 
 
 fT->Branch("eta1",&evnt.eta1,"eta1/F");
 fT->Branch("eta2",&evnt.eta2,"eta2/F");
 fT->Branch("phi1",&evnt.phi1,"phi1/F");
 fT->Branch("phi2",&evnt.phi2,"phi2/F");
 fT->Branch("pt1",&evnt.pt1,"pt1/F");
 fT->Branch("pt2",&evnt.pt2,"pt2/F");
 fT->Branch("x1",&evnt.x1,"x1/F");
 fT->Branch("x2",&evnt.x2,"x2/F");
 // fT->Branch("dphi",&evnt.dphi,"dphi/F");
 fT->Branch("etadijet",&evnt.etadijet,"etadijet/F");
 
 nt * fdijet = new nt("/afs/cern.ch/user/d/dgulhan/workDir/pythiaJets2/merged.root");
 hi * fgen = new hi("/afs/cern.ch/user/d/dgulhan/workDir/pythiaJets2/merged.root");
 
 int            event=0;
 float          b;
 int            nTrk;
 float          frame_E;
 float          total_E;
 float          x;
 float          trkPt[20000]; 
 float          trkPhi[20000]; 
 float          trkEta[20000]; 
 float          trkRapidity[20000]; 
 float          trkE[20000];
 float          trkPx[20000];
 float          trkPy[20000];
 float          trkPz[20000];
 float          boostedPt[20000]; 
 float          boostedPhi[20000]; 
 float          boostedEta[20000]; 
 float          boostedRapidity[20000]; 
 float          boostedE[20000];
 float          boostedPx[20000];
 float          boostedPy[20000];
 float          boostedPz[20000];
 double         beta;
 float          pt1;
 float          pt2;
 float          dphi;
 float          eta1;
 float          eta2;
 float          phi1;
 float          phi2;
 float          etadijet;
 float          x1;
 float          x2;
 
 int ievent=0;
 int nentries = fdijet->GetEntriesFast();  
 bool isfirst=true;
 for(vector <Event>::iterator it = event_vector.begin(); it != event_vector.end(); ++it){
  bool passedDijetCuts=false;
  while(!passedDijetCuts){ 
   fgen->GetEntry(ievent); 
   fdijet->GetEntry(ievent); 
   if(fdijet->pt1>120 && fdijet->pt2>30 && acos(cos(fdijet->phi1-fdijet->phi2))>2.094){
    passedDijetCuts=true;
    pt1=fdijet->pt1;
    pt2=fdijet->pt2;
    eta1=fdijet->eta1;
    eta2=fdijet->eta2;
    phi1=fdijet->phi1;
    phi2=fdijet->phi2;
    dphi=acos(cos(phi1-phi2));
    etadijet=(eta1+eta2)/2;
    x1=fgen->x1;
    x2=fgen->x2;
    if(isfirst){
     isfirst=false;
     passedDijetCuts=false;
    }
   }
   ievent++;
  }
  
  b=it->get_b();
  nTrk=it->get_number_of_particles();
  frame_E=it->get_frame_E();
  total_E=it->get_total_E();
  x=1-pow(it->frame_E/5020,2);
    if(ievent==0) cout<<b<<" "<<nTrk<<" "<<frame_E<<" "<<total_E<<endl;

  float          ETpos=0;
  float          ETneg=0;
  int            Ntrkpos=0;
  int            Ntrkneg=0;
  int            boostedNtrkpos=0;
  int            boostedNtrkneg=0;
  float          boostedETpos=0;
  float          boostedETneg=0;
  beta=sqrt(1-1/(cosh(0.464*log(1-x))*cosh(0.464*log(1-x))));
  // cout<<"beta"<<beta<<endl;
  // cout<<"b= "<<it->get_b()<<" total_E= "<<it->get_total_E()<<" number_of_particles"<<it->get_number_of_particles()<<endl;
  int ntrk=0;
  for(int i=0;i<it->get_number_of_particles();i++){
  // for(int i=0;i<1;i++){
   Particle * current_part=it->get_ieth_particle(i);
   trkPx[ntrk]=current_part->get_px();
   trkPy[ntrk]=current_part->get_py();
   trkPz[ntrk]=current_part->get_pz();
   trkE[ntrk]=current_part->get_E();
   if(trkPx[ntrk]==0) continue;
   TLorentzVector v;
   v.SetPxPyPzE(trkPx[ntrk],trkPy[ntrk],trkPz[ntrk],trkE[ntrk]);
   v.Boost(0,0,beta); 
   // cout<<"beta"<<beta<<endl;
   boostedPx[ntrk]=v.Px();
   boostedPz[ntrk]=v.Py();
   boostedPy[ntrk]=v.Pz();
   boostedPt[ntrk]=v.Pt();
   boostedE[ntrk]=v.E();
   boostedPhi[ntrk]=v.Phi();
   boostedEta[ntrk]=v.Eta();
   boostedRapidity[ntrk]=v.Rapidity();
   trkPt[ntrk]=sqrt(pow(trkPx[ntrk],2)+pow(trkPy[ntrk],2));
   trkEta[ntrk] = 0.5*log((sqrt(trkPt[ntrk]*trkPt[ntrk]+trkPz[ntrk]*trkPz[ntrk])+trkPz[ntrk])/(sqrt(trkPt[ntrk]*trkPt[ntrk]+trkPz[ntrk]*trkPz[ntrk])-trkPz[ntrk]));
   trkPhi[ntrk]=atan2(trkPy[ntrk],trkPx[ntrk]);
   trkRapidity[ntrk]=0.5*log((trkE[ntrk]*trkE[ntrk]+trkPz[ntrk]*trkPz[ntrk])/(trkE[ntrk]*trkE[ntrk]-trkPz[ntrk]*trkPz[ntrk]));
   // cout<<"i="<<i<<" px= "<<current_part->get_px()<<" py= "<<current_part->get_py()<<" pz= "<<current_part->get_pz()<<" E= "<<current_part->get_E()<<endl;
   // cout<<"eta="<<trkEta[i]<<" phi="<<trkPhi[i]<<" pt= "<<trkPt[i]<<" pz "<<trkPz[i]<<endl;
   if(boostedEta[ntrk]>4 && boostedEta[ntrk]<5){
    boostedETpos+=boostedE[ntrk]/cosh(boostedEta[ntrk]);
    boostedNtrkpos++;
   }
   if(boostedEta[ntrk]<-4 && boostedEta[ntrk]>-5){
    boostedETneg+=boostedE[ntrk]/cosh(boostedEta[ntrk]);
    boostedNtrkneg++;
   }
   if(trkEta[ntrk]>4 && trkEta[ntrk]<5){
    ETpos+=trkE[ntrk]/cosh(trkEta[ntrk]);
    Ntrkpos++;
   }
   if(trkEta[ntrk]<-4 && trkEta[ntrk]>-5){
    ETneg+=trkE[ntrk]/cosh(trkEta[ntrk]);
    Ntrkneg++;
   }
   ntrk++;
  }

  evnt.setEvent(event,b,frame_E,total_E,x,beta,ETpos,ETneg,boostedETpos,boostedETneg,Ntrkpos,Ntrkneg,boostedNtrkpos,boostedNtrkneg,ntrk,trkPt,trkPhi,trkEta,trkPx,trkPy,trkPz,trkE,trkRapidity,boostedPt,boostedPhi,boostedEta,boostedPx,boostedPy,boostedPz,boostedE,boostedRapidity,pt1,pt2,phi1,phi2,eta1,eta2,etadijet,x1,x2);
  // evnt.setEvent(event,b,frame_E,total_E,x,beta,ETpos,ETneg,boostedETpos,boostedETneg,Ntrkpos,Ntrkneg,boostedNtrkpos,boostedNtrkneg,ntrk,trkPt,trkPhi,trkEta,trkPx,trkPy,trkPz,trkE,trkRapidity,boostedPt,boostedPhi,boostedEta,boostedPx,boostedPy,boostedPz,boostedE,boostedRapidity);
  event++;
  // cout<<evnt.beta<<endl;
  fT->Fill();
  evnt.reset();
 }
 outfile->cd();
 fT->Write();
 outfile->Close();
}

