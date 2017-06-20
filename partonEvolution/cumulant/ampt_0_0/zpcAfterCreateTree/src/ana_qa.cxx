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
TH2D *parton_xy = new TH2D("parton_xy","parton_xy",2000,-10,10,2000,-10,10);

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
				//filein>>nevent>>nrun>>multi>>impactpar>>NpartP>>NpartT>>NELP>>NINP>>NELT>>NINT>>passhead;//for ampt.dat
				filein>>nevent>>nrun>>multi>>impactpar>>passhead>>passhead>>passhead>>passhead; //for zpc.dat
				if(!filein.good()) break;
				if(impactpar>2.001)continue;//b=0
				if(count==1){break;}//restore the 1st central event
                cout<<"This is the 1st central event, with b::mult "<<impactpar<<" :: "<<multi<<endl;
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
					parton_xy->Fill(x,y);

				}
				count++;
			}
			filein.close();
			filenum++;
		}
	}
	inputStream->close();
	cout<< "read in "<<filenum <<" good files w/ "<<count<<" events"<<endl;

	parton_xy->Write();
	file->Write();

	delete parton_xy;
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
