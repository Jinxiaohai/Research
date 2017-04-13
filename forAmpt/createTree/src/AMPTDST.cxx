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

const Int_t kMaxtrack = 25000;

typedef struct { 
    Int_t        nevent;
    Int_t        nrun;
    Int_t        multi; // Multiplicity
    Float_t      impactpar;
    Int_t        NpartP;
    Int_t        NpartT;
    Int_t        NELP;
    Int_t        NINP;
    Int_t        NELT;
    Int_t        NINT;

    Int_t        Indx[kMaxtrack];
    Int_t        ID[kMaxtrack];      //KF code of particle
    Float_t      fPx[kMaxtrack];           //X component of the momentum
    Float_t      fPy[kMaxtrack];           //Y component of the momentum
    Float_t      fPz[kMaxtrack];           //Z component of the momentum
    Float_t      fMass[kMaxtrack];        //The mass  of this particle
    Float_t      fX[kMaxtrack];       //X coordinate of the first point
    Float_t      fY[kMaxtrack];       //Y coordinate of the first point
    Float_t      fZ[kMaxtrack];       //Z coordinate of the first point
    Float_t      fT[kMaxtrack];       //T coordinate of the first point

} Cell_t;

double PI = TMath::Pi();
double ANUMBER = 197.;

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

    ifstream filein;
    filein.open(FileInput);
    cout<<FileInput<<endl;

    if(!filein) return(0);

    Int_t nevent,nrun,multi,NpartP,NpartT,NELP,NINP,NELT,NINT;
    Float_t impactpar;
    Float_t passhead;
    Float_t px,py,pz,am,x,y,z,time;
	//Float_t energy;
    
    char others1[256],others2[256],others3[256],others4[256];
    Int_t id;

    char outfile[256];
    sprintf(outfile,"%s.root",FileOutput); 
    TFile *file = new TFile(outfile,"RECREATE");
    TTree *tree = new TTree("AMPT","AMPT DST Tree");
    Cell_t cell;
    tree->Branch("Event",&cell.nevent,"nevent/I:nrun:multi:impactpar/F:NpartP/I:NpartT:NELP:NINP:NELT:NINT");
    tree->Branch("Indx",cell.Indx,"Indx[multi]/I");
    tree->Branch("ID",cell.ID,"ID[multi]/I");
    //  tree->Branch("Momentum",cell.fPx,"fPx[multi]/F:fPy[multi]:fPz[multi]");
    tree->Branch("Px",cell.fPx,"fPx[multi]/F");
    tree->Branch("Py",cell.fPy,"fPy[multi]/F");
    tree->Branch("Pz",cell.fPz,"fPz[multi]/F");  
    tree->Branch("Mass",cell.fMass,"fMass[multi]/F");
    //  tree->Branch("Coodinate",cell.fX,"fX[multi]/F:fY[multi]:fZ[multi]");
    tree->Branch("X",cell.fX,"fX[multi]/F");
    tree->Branch("Y",cell.fY,"fY[multi]/F");
    tree->Branch("Z",cell.fZ,"fZ[multi]/F");  
    tree->Branch("Time",cell.fT,"fT[multi]/F");

    Int_t count=0;
    while(filein.good()){

	if (count%100==0)cout << count << endl;
	filein>>nevent>>nrun>>multi>>impactpar>>NpartP>>NpartT>>NELP>>NINP>>NELT>>NINT>>passhead;
	if(!filein.good()){ cout <<"Nevent----> " << nevent << " Multi---> "<<multi<< endl; break;}
	//if(!filein.good()) break;

	cell.nevent   =nevent;       
	cell.nrun     =nrun;         
	cell.multi    =multi;        
	cell.impactpar=impactpar;    
	cell.NpartP   =NpartP;       
	cell.NpartT   =NpartT;       
	cell.NELP     =NELP;         
	cell.NINP     =NINP;         
	cell.NELT     =NELT;         
	cell.NINT     =NINT;         

	for (Int_t nlines=0;nlines<multi;nlines++) {
	    //filein >> id >> px >> py >> pz >> am >> x >> y >> z >> time;
	    filein >> id >> px >> py >> pz >> am >> others1 >> others2 >> others3 >> others4;
	    //filein >> indx >> id >> px >> py >> pz >> am >> others1 >> others2 >> others3 >> others4;

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
            
		cell.Indx[nlines]    = 1;
	    cell.ID[nlines]    = id;
	    cell.fPx[nlines]   = px;
	    cell.fPy[nlines]   = py;
	    cell.fPz[nlines]   = pz;
	    cell.fMass[nlines] = am;
	    cell.fX[nlines]    = x;
	    cell.fY[nlines]    = y;
	    cell.fZ[nlines]    = z;
	    cell.fT[nlines]    = time;
	}
	tree->Fill();
	count++;
    }
    filein.close();
    file->Write();
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
