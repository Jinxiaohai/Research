//***************************************************************************
//Calculting the various order cumulants and the errors based on factoial moments.
//One can also refer the paper: X. Luo, arXiv: 1410.3914.
//**************************************************************************

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif
#ifndef __CINT__
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "math.h"
#include "string.h"

#include "TROOT.h"
#include "TFile.h"

#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TGraphErrors.h"
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
#include <iomanip>
#include "Centrality.h"
#include "EffPara.h"


using namespace std;
#define PI 3.14159 

const int MAXSIZE=800;


double CalCumulants(const char *name, int order, TFile *, int cent);
double Cal_Poisson_base(const char *name, int order, TFile *, int cent);
double Cal_IndProd_base(const char *name, int order, TFile *, int cent);
double Cal_NBD_base(const char *name, int order, TFile *, int cent);


double Moments_Errors(const char *, const char *, TFile *, TFile *,int);
double Moments_IndProd_Errors(const char *, const char *, TFile *, int );
 
int Delta(int a, int b);
double lnchoose(int n, int m);
double C(int, int);
double StirlingS1(int, int);  
double StirlingS2(int, int);
double UnsignS1(int, int);
double Moments4(int m, int n, int s, int t, double ff[][10][10][10]);
double Moments2(int m, int n, double ff[][10]);


double FM(int, int, double ff[][10][10][10]);
 
double CM(int order, double FF[][10], const char *name);
double Cumulant(int n, double FF[][10], const char *name);
double ffff(int i, int j, int k, int h,TFile *f, int cent, int ref);
double DD(int n,double mean, double F[][10],int i, int j, int k, int h);
double D_F2_F4(int r1, int r2, int u, int v, int j, int k);
double D_Moment_F2(int n, int u, int v,double mean, double F[][10]);
double D_CM_F2(int n, int u, int v, double mean,double FF[][10]);
double CovMoment4(int i, int j, int k, int h, int r, int s, int u, int v,double FF[][10][10][10]);
double CovFact4( int i, int j, int k, int h, int r, int s, int u, int v,double FF[][10][10][10]);

TF1 *fun_plow=new TF1("fun_plow","[2]*x*x+[1]*x+[0]",0,800);
TF1 *fun_phigh=new TF1("fun_phigh","[2]*x*x+[1]*x+[0]",0,800);
TF1 *fun_aplow=new TF1("fun_aplow","[2]*x*x+[1]*x+[0]",0,800);
TF1 *fun_aphigh=new TF1("fun_aphigh","[2]*x*x+[1]*x+[0]",0,800);

   
 	

int main(int argc, char *argv[])
{

  for(int i=0;i<=2;i++)
    {
      fun_plow->SetParameter(i,para[0][i]);
      fun_phigh->SetParameter(i,para[1][i]);
      fun_aplow->SetParameter(i,para[2][i]);
      fun_aphigh->SetParameter(i,para[3][i]);
    }


  const  char *name2[10]={"dummy","C1","C2","C3","C4","KV","SD","VM","KDS"};
 	 
 	
  const   char *name_net[20][4]={{"C1","C2","C3","C4"},{"KV","SD","VM","KDS"},{"C1_base_pos"
									       ,"C2_base_pos","C3_base_pos","C4_base_pos"},{"KV_base_pos","SD_base_pos","VM_base_pos","KDS_base_pos"}
				 ,{"C1_base_NBD","C2_base_NBD","C3_base_NBD","C4_base_NBD"},{"KV_base_NBD","SD_base_NBD","VM_base_NBD","KDS_base_NBD"},
				 {"C1_base_Ind","C2_base_Ind","C3_base_Ind","C4_base_Ind"},{"KV_base_Ind","SD_base_Ind","VM_base_Ind","KDS_base_Ind"}
				 ,{"C1_pos_ratio","C2_pos_ratio","C3_pos_ratio","C4_pos_ratio"},{"KV_pos_ratio","SD_pos_ratio","VM_pos_ratio","KDS_pos_ratio"
				 },{"C1_NBD_ratio","C2_NBD_ratio","C3_NBD_ratio","C4_NBD_ratio"},{"KV_NBD_ratio","SD_NBD_ratio","VM_NBD_ratio","KDS_NBD_ratio"
				 },{"C1_Ind_ratio","C2_Ind_ratio","C3_Ind_ratio","C4_Ind_ratio"},{"KV_Ind_ratio","SD_Ind_ratio","VM_Ind_ratio","KDS_Ind_ratio"}};
 	 
 	 

	 
  const char *name="netp";
	 
  TFile *f=new TFile("Fact.root");
  TFile *f4=new TFile("Fact4.root");
 	 

  double Cum[10][10];
  double Poisson_base[10][10];
  double NBD_base[10][10];
  double IndProd_base[10][10];
  
  double Cum_err[10][10];
  
  double IndProd_base_err[10][10];
  
  
    
  double Cum_Poisson_ratio[10][10];
  double Cum_NBD_ratio[10][10];
   
  double Cum_Poisson_ratio_err[10][10];
  double Cum_NBD_ratio_err[10][10];
   
   
   
   
   
  double Cum_IndProd_ratio[10][10];
  double Cum_IndProd_ratio_err[10][10];

  
  double KV[10];
  double SD[10];
  double VM[10];
  double KDS[10];
  
  double KV_err[10];
  double SD_err[10];
  double VM_err[10];
  double KDS_err[10];
  
  
  double KV_Poisson_base[10];
  double SD_Poisson_base[10];
  double VM_Poisson_base[10];
  double KDS_Poisson_base[10];
  
  double KV_NBD_base[10];  
  double SD_NBD_base[10];  
  double VM_NBD_base[10];  
  double KDS_NBD_base[10]; 

  double KV_IndProd_base[10];  
  double SD_IndProd_base[10];  
  double VM_IndProd_base[10];  
  double KDS_IndProd_base[10]; 

  double KV_IndProd_base_err[10];  
  double SD_IndProd_base_err[10];  
  double VM_IndProd_base_err[10];  
  double KDS_IndProd_base_err[10]; 


  double KV_Poisson_ratio[10];
  double SD_Poisson_ratio[10];
  double VM_Poisson_ratio[10];
  double KDS_Poisson_ratio[10];
  
  double KV_Poisson_ratio_err[10];
  double SD_Poisson_ratio_err[10];
  double VM_Poisson_ratio_err[10];
  double KDS_Poisson_ratio_err[10];
  
  
  
  
  double KV_NBD_ratio[10];
  double SD_NBD_ratio[10];
  double VM_NBD_ratio[10];
  double KDS_NBD_ratio[10];
  
  
  double KV_NBD_ratio_err[10];
  double SD_NBD_ratio_err[10];
  double VM_NBD_ratio_err[10];
  double KDS_NBD_ratio_err[10];
  
  double KV_IndProd_ratio[10];
  double SD_IndProd_ratio[10];
  double VM_IndProd_ratio[10];
  double KDS_IndProd_ratio[10];
  
  double KV_IndProd_ratio_err[10];
  double SD_IndProd_ratio_err[10];
  double VM_IndProd_ratio_err[10];
  double KDS_IndProd_ratio_err[10];
 
  TGraphErrors *gr[15];
  
  for(int i=0; i<=9; i++)
    {
	  	
	  
	  	
      KV[i]=0;
      SD[i]=0;
      VM[i]=0;
      KDS[i]=0;
   
      KV_err[i]=0;
      SD_err[i]=0;
      VM_err[i]=0;
      KDS_err[i]=0;
   
      KV_Poisson_base[i]=0;
      SD_Poisson_base[i]=0;
      VM_Poisson_base[i]=0;
      KDS_Poisson_base[i]=0;
  
      KV_NBD_base[i]=0;  
      SD_NBD_base[i]=0;  
      VM_NBD_base[i]=0;  
      KDS_NBD_base[i]=0; 

      KV_IndProd_base[i]=0;  
      SD_IndProd_base[i]=0;  
      VM_IndProd_base[i]=0;  
      KDS_IndProd_base[i]=0; 
 
      KV_IndProd_base_err[i]=0;  
      SD_IndProd_base_err[i]=0;   
      VM_IndProd_base_err[i]=0;  
      KDS_IndProd_base_err[i]=0;
 

      KV_Poisson_ratio[i]=0;
      SD_Poisson_ratio[i]=0;
      VM_Poisson_ratio[i]=0;
      KDS_Poisson_ratio[i]=0;
  
      KV_NBD_ratio[i]=0;
      SD_NBD_ratio[i]=0;
      VM_NBD_ratio[i]=0;
      KDS_NBD_ratio[i]=0;
   
   
      KV_Poisson_ratio_err[i]=0;
      SD_Poisson_ratio_err[i]=0;
      VM_Poisson_ratio_err[i]=0;
      KDS_Poisson_ratio_err[i]=0;
  
      KV_NBD_ratio_err[i]=0;
      SD_NBD_ratio_err[i]=0;
      VM_NBD_ratio_err[i]=0;
      KDS_NBD_ratio_err[i]=0;
   
   
  
      KV_IndProd_ratio[i]=0;
      SD_IndProd_ratio[i]=0;
      VM_IndProd_ratio[i]=0;
      KDS_IndProd_ratio[i]=0;
  
      KV_IndProd_ratio_err[i]=0;
      SD_IndProd_ratio_err[i]=0;
      VM_IndProd_ratio_err[i]=0;
      KDS_IndProd_ratio_err[i]=0;
	  	
      for(int j=0;j<=9;j++)
	{
	  	
	  Cum[i][j]=0;
	  Cum_err[i][j]=0;
	 
	 
	  Poisson_base[i][j]=0;
	  NBD_base[i][j]=0;
   
   
   
	  IndProd_base[i][j]=0;
	  IndProd_base_err[i][j]=0;
	 
	  Cum_Poisson_ratio[i][j]=0;
	  Cum_Poisson_ratio_err[i][j]=0;

	  Cum_NBD_ratio[i][j]=0;
	  Cum_NBD_ratio_err[i][j]=0;
    
	  Cum_IndProd_ratio[i][j]=0;
	  Cum_IndProd_ratio_err[i][j]=0;

	 
	}
	  
    }
	  
	


	 
	
  //	cout<<"******************"<<"Centrality: "<<cent+1<<"***************"<<endl;
  for(int j=1;j<=4;j++)	
    {
	 	   	
	 	  
      cout<<"******************"<<name2[j]<<"*****************"<<endl;
	 	   	
      for(int cent=0;cent<=NumberOfCentrality-1;cent++)	
	{
	  cout<<"******************"<<"Centrality: "<<cent+1<<"***************"<<endl;
 
	  Cum[j][cent]=CalCumulants(name,j,f,cent+1); 
          
 
 	 
 	 
 	 
 	 
	 	  	
	}
    }
	 

  for(int j=1;j<=4;j++)	
    {
      
      if(strcmp(name2[j],argv[1])!=0) continue;
	 	   	
      //		cout<<"******************"<<name2[j]<<"*****************"<<endl;

	 	   	
      for(int cent=0;cent<=NumberOfCentrality-1;cent++)	
	{
	  
	  Cum_err[j][cent]=Moments_Errors(name,name2[j],f,f4,cent+1); 
	  Poisson_base[j][cent]=Cal_Poisson_base(name,j,f,cent+1);          
	  NBD_base[j][cent]=Cal_NBD_base(name,j,f,cent+1); 
 	
	  IndProd_base[j][cent]=Cal_IndProd_base(name,j,f,cent+1);
	  // 	IndProd_base_err[j][cent]=Moments_IndProd_Errors(name,name2[j],f,cent+1);
 	
 	
	  Cum_Poisson_ratio[j][cent]=Cum[j][cent]/Poisson_base[j][cent]; 
	  Cum_Poisson_ratio_err[j][cent]=Cum_err[j][cent]/Poisson_base[j][cent]; 
	  Cum_NBD_ratio[j][cent]=Cum[j][cent]/NBD_base[j][cent];	
	  Cum_NBD_ratio_err[j][cent]=Cum_err[j][cent]/NBD_base[j][cent];  
	  Cum_IndProd_ratio[j][cent]=Cum[j][cent]/IndProd_base[j][cent];   	
	  Cum_IndProd_ratio_err[j][cent]=Cum_err[j][cent]/IndProd_base[j][cent];   	
	}
	 	  
      gr[0]=new TGraphErrors(NumberOfCentrality,Npart,Cum[j],0,Cum_err[j]);
      gr[1]=new TGraphErrors(NumberOfCentrality,Npart,Poisson_base[j],0,0);
      gr[2]=new TGraphErrors(NumberOfCentrality,Npart,NBD_base[j],0,0);
      gr[3]=new TGraphErrors(NumberOfCentrality,Npart,IndProd_base[j],0,IndProd_base_err[j]);
      gr[4]=new TGraphErrors(NumberOfCentrality,Npart,Cum_Poisson_ratio[j],0,Cum_Poisson_ratio_err[j]);
      gr[5]=new TGraphErrors(NumberOfCentrality,Npart,Cum_NBD_ratio[j],0,Cum_NBD_ratio_err[j]);
      gr[6]=new TGraphErrors(NumberOfCentrality,Npart,Cum_IndProd_ratio[j+1],0,Cum_IndProd_ratio_err[j]);
	 	  
	 	
    }
	 
	

 
	   
	 	   	
           
  if(strcmp("KV",argv[1])==0)
    {
 
      for(int cent=0;cent<=NumberOfCentrality-1;cent++)	
	{		
 
	  KV[cent]=Cum[4][cent]/Cum[2][cent]; 		
	  cout<<argv[1]<<": "<<KV[cent]<<endl;
	  KV_err[cent]=Moments_Errors(name,"KV",f,f4,cent+1);  
	  KV_Poisson_base[cent]=Poisson_base[4][cent]/Poisson_base[2][cent];  
	  KV_NBD_base[cent]=NBD_base[4][cent]/NBD_base[2][cent];  
	  KV_IndProd_base[cent]=IndProd_base[4][cent]/IndProd_base[2][cent]; 
	  //KV_IndProd_base_err[cent]=Moments_IndProd_Errors(name,"KV",f,cent+1);
	  KV_Poisson_ratio[cent]=KV[cent]/KV_Poisson_base[cent];
	  KV_Poisson_ratio_err[cent]=KV_err[cent]/KV_Poisson_base[cent];
	  KV_NBD_ratio[cent]=KV[cent]/KV_NBD_base[cent];
	  KV_NBD_ratio_err[cent]=KV_err[cent]/KV_NBD_base[cent];
	  KV_IndProd_ratio[cent]=KV[cent]/KV_IndProd_base[cent]; 
	  KV_IndProd_ratio_err[cent]=KV_err[cent]/KV_IndProd_base[cent]; 
	}



      gr[0]=new TGraphErrors(NumberOfCentrality,Npart,KV,0,KV_err);
      gr[1]=new TGraphErrors(NumberOfCentrality,Npart,KV_Poisson_base,0,0);
      gr[2]=new TGraphErrors(NumberOfCentrality,Npart,KV_NBD_base,0,0);
      gr[3]=new TGraphErrors(NumberOfCentrality,Npart,KV_IndProd_base,0,KV_IndProd_base_err);
      gr[4]=new TGraphErrors(NumberOfCentrality,Npart,KV_Poisson_ratio,0,KV_Poisson_ratio_err);
      gr[5]=new TGraphErrors(NumberOfCentrality,Npart,KV_NBD_ratio,0,KV_NBD_ratio_err);
      gr[6]=new TGraphErrors(NumberOfCentrality,Npart,KV_IndProd_ratio,0,KV_IndProd_ratio_err);

   


    }
 
  if(strcmp("SD",argv[1])==0)
    {
 		
      for(int cent=0;cent<=NumberOfCentrality-1;cent++)	
	{		
	  SD[cent]=Cum[3][cent]/Cum[2][cent]; 		
	  cout<<argv[1]<<": "<<SD[cent]<<endl;
	  SD_err[cent]=Moments_Errors(name,"SD",f,f4,cent+1);  
	  SD_Poisson_base[cent]=Poisson_base[3][cent]/Poisson_base[2][cent];  
	  SD_NBD_base[cent]=NBD_base[3][cent]/NBD_base[2][cent];  
	  SD_IndProd_base[cent]=IndProd_base[3][cent]/IndProd_base[2][cent]; 
	  //SD_IndProd_base_err[cent]=Moments_IndProd_Errors(name,"SD",f,cent+1);
	  SD_Poisson_ratio[cent]=SD[cent]/SD_Poisson_base[cent];
	  SD_Poisson_ratio_err[cent]=SD_err[cent]/SD_Poisson_base[cent];
	  SD_NBD_ratio[cent]=SD[cent]/SD_NBD_base[cent];
	  SD_NBD_ratio_err[cent]=SD_err[cent]/SD_NBD_base[cent];
	  SD_IndProd_ratio[cent]=SD[cent]/SD_IndProd_base[cent]; 
	  SD_IndProd_ratio_err[cent]=SD_err[cent]/SD_IndProd_base[cent]; 
	}
 
 
 
 
      gr[0]=new TGraphErrors(NumberOfCentrality,Npart,SD,0,SD_err);
      gr[1]=new TGraphErrors(NumberOfCentrality,Npart,SD_Poisson_base,0,0);
      gr[2]=new TGraphErrors(NumberOfCentrality,Npart,SD_NBD_base,0,0);
      gr[3]=new TGraphErrors(NumberOfCentrality,Npart,SD_IndProd_base,0,SD_IndProd_base_err);
      gr[4]=new TGraphErrors(NumberOfCentrality,Npart,SD_Poisson_ratio,0,SD_Poisson_ratio_err);
      gr[5]=new TGraphErrors(NumberOfCentrality,Npart,SD_NBD_ratio,0,SD_NBD_ratio_err);
      gr[6]=new TGraphErrors(NumberOfCentrality,Npart,SD_IndProd_ratio,0,SD_IndProd_ratio_err);
 
 
 
 
 
    }
 
  if(strcmp("VM",argv[1])==0)
    {
 		
 		
      for(int cent=0;cent<=NumberOfCentrality-1;cent++)	
	{		
 		
	  VM[cent]=Cum[2][cent]/Cum[1][cent]; 		
	  cout<<argv[1]<<": "<<VM[cent]<<endl;
	  VM_err[cent]=Moments_Errors(name,"VM",f,f4,cent+1);  
	  VM_Poisson_base[cent]=Poisson_base[2][cent]/Poisson_base[1][cent];  
	  VM_NBD_base[cent]=NBD_base[2][cent]/NBD_base[1][cent];  
	  VM_IndProd_base[cent]=IndProd_base[2][cent]/IndProd_base[1][cent]; 
	  //VM_IndProd_base_err[cent]=Moments_IndProd_Errors(name,"VM",f,cent+1);
	  VM_Poisson_ratio[cent]=VM[cent]/VM_Poisson_base[cent];
	  VM_Poisson_ratio_err[cent]=VM_err[cent]/VM_Poisson_base[cent];
	  VM_NBD_ratio[cent]=VM[cent]/VM_NBD_base[cent];
	  VM_NBD_ratio_err[cent]=VM_err[cent]/VM_NBD_base[cent];
	  VM_IndProd_ratio[cent]=VM[cent]/VM_IndProd_base[cent]; 
	  VM_IndProd_ratio_err[cent]=VM_err[cent]/VM_IndProd_base[cent]; 
	}
 
      gr[0]=new TGraphErrors(NumberOfCentrality,Npart,VM,0,VM_err);
      gr[1]=new TGraphErrors(NumberOfCentrality,Npart,VM_Poisson_base,0,0);
      gr[2]=new TGraphErrors(NumberOfCentrality,Npart,VM_NBD_base,0,0);
      gr[3]=new TGraphErrors(NumberOfCentrality,Npart,VM_IndProd_base,0,VM_IndProd_base_err);
      gr[4]=new TGraphErrors(NumberOfCentrality,Npart,VM_Poisson_ratio,0,VM_Poisson_ratio_err);
      gr[5]=new TGraphErrors(NumberOfCentrality,Npart,VM_NBD_ratio,0,VM_NBD_ratio_err);
      gr[6]=new TGraphErrors(NumberOfCentrality,Npart,VM_IndProd_ratio,0,VM_IndProd_ratio_err);

    }
 
  if(strcmp("KDS",argv[1])==0)
    {
 		
      for(int cent=0;cent<=NumberOfCentrality-1;cent++)	
	{		
 		
 		
	  KDS[cent]=Cum[4][cent]/Cum[3][cent]; 		
	  cout<<argv[1]<<": "<<KDS[cent]<<endl;
	  KDS_err[cent]=Moments_Errors(name,"KDS",f,f4,cent+1);  
	  KDS_Poisson_base[cent]=Poisson_base[4][cent]/Poisson_base[3][cent];  
	  KDS_NBD_base[cent]=NBD_base[4][cent]/NBD_base[3][cent];  
	  KDS_IndProd_base[cent]=IndProd_base[4][cent]/IndProd_base[3][cent]; 
	  //KDS_IndProd_base_err[cent]=Moments_IndProd_Errors(name,"KDS",f,cent+1);
	  KDS_Poisson_ratio[cent]=KDS[cent]/KDS_Poisson_base[cent];
	  KDS_Poisson_ratio_err[cent]=KDS_err[cent]/KDS_Poisson_base[cent];
	  KDS_NBD_ratio[cent]=KDS[cent]/KDS_NBD_base[cent];
	  KDS_NBD_ratio_err[cent]=KDS_err[cent]/KDS_NBD_base[cent];
	  KDS_IndProd_ratio[cent]=KDS[cent]/KDS_IndProd_base[cent]; 
	  KDS_IndProd_ratio_err[cent]=KDS_err[cent]/KDS_IndProd_base[cent]; 
	} 
 

          
      gr[0]=new TGraphErrors(NumberOfCentrality,Npart,KDS,0,KDS_err);
      gr[1]=new TGraphErrors(NumberOfCentrality,Npart,KDS_Poisson_base,0,0);
      gr[2]=new TGraphErrors(NumberOfCentrality,Npart,KDS_NBD_base,0,0);
      gr[3]=new TGraphErrors(NumberOfCentrality,Npart,KDS_IndProd_base,0,KDS_IndProd_base_err);
      gr[4]=new TGraphErrors(NumberOfCentrality,Npart,KDS_Poisson_ratio,0,KDS_Poisson_ratio_err);
      gr[5]=new TGraphErrors(NumberOfCentrality,Npart,KDS_NBD_ratio,0,KDS_NBD_ratio_err);
      gr[6]=new TGraphErrors(NumberOfCentrality,Npart,KDS_IndProd_ratio,0,KDS_IndProd_ratio_err);
                       
                       
    }
 



            
           
           
         

    
  TFile *ff=new TFile(Form("Moments_%s.root",argv[1]),"recreate");
    
  ff->cd();
     
  gr[0]->SetName(Form("%s",argv[1]));
  gr[1]->SetName(Form("%s_Poisson_base",argv[1]));  	
  gr[2]->SetName(Form("%s_NBD_base",argv[1]));   	
  gr[3]->SetName(Form("%s_IndProd_base",argv[1]));   	
  gr[4]->SetName(Form("%s_Poisson_ratio",argv[1]));          
  gr[5]->SetName(Form("%s_NBD_ratio",argv[1]));          
  gr[6]->SetName(Form("%s_IndProd_ratio",argv[1]));     
	
	
  for(int i=0;i<=6;i++)
    {
	  	
      gr[i]->Write();
    }
 
  cout<<"Finished !!!"<<endl;
        
  ff->Write();
  ff->Close();  
	  	
	  
  	    
	  
	
	  
       
       
  return 0;
}


double CalCumulants(const char *name, int n, TFile *f, int cent)
{
  //cout<<name<<": "<<n+cent<<endl;

  double F[10][10];
  double Events=0;
  double TotalEvents=0;
  double Cum[10]={0};
	
	
	
  for(int i=0; i<=9; i++)
    {
	  	
      for(int j=0;j<=9;j++)
	{
	  	
	  F[i][j]=0;
	  	
	}
	  
    }
	
	
  for(int ref=1;ref<=1000;ref++)
    {	
		
      Events=((TProfile *)f->Get(Form("F%d%d_Cent%d",1,0,cent)))->GetBinEntries(ref+1);
	

      if(Events<=10) continue;
   
		
      //	cout<<Events<<endl;
	
      for(int i=0; i<=8; i++)
	{
	  	
	  for(int j=0;j<=8;j++)
	    {
	  
	      if((i+j)<=8)
		{  

		  F[i][j]=((TProfile *)f->Get(Form("F%d%d_Cent%d",i,j,cent)))->GetBinContent(ref+1);

		}
	    }
	}
  
     
  
   
      TotalEvents+=Events; 
  
 
      for(int i=1;i<=4;i++)
	{
 			
 			
	  Cum[i]+=Cumulant(i,F,name)*Events;
 			
	}
 		
 

      // cout<<ref<<","<<Events<<","<<(F[1][0]+F[0][1]+F[0][2]-pow(F[1][0]-F[0][1],2)-2*F[1][1]+F[2][0])<<","<<Cumulant(2,F,"netp")<<endl;


    }

 	

  if(n==0) return 1;
  if(n==1) return Cum[1]/TotalEvents;
  if(n==2) return Cum[2]/TotalEvents;	
  if(n==3) return Cum[3]/TotalEvents;
  if(n==4) return Cum[4]/TotalEvents;	
	
	
	
  return 0;	
	
	


}



double Moments_Errors(const char *name, const char *name2, TFile *f, TFile *f4,int cent)
{
  cout<<"Centrality: "<<cent<<", Calculating Errors: "<<name<<","<<name2<<endl;
	
  double F[10][10];
  double Events=0;
  double TotalEvents=0;
  //double NCM[10]={0};
  double Cum_err=0;
  double D[10][10][10][10];
  double ff[10][10][10][10];
	
  for(int i=0; i<=9; i++)
    {
	  	
      for(int j=0;j<=9;j++)
	{
	  	
	  F[i][j]=0;
	  	
	  for(int k=0;k<=9;k++)
	    {
	  		
	      for(int h=0;h<=9;h++)
		{
	  			
		  ff[i][j][k][h]=0;
		  D[i][j][k][h]=0;
	  			
		}
	  		
	    }
	  	
	  	
	  	
	}
	  
    }
	
	 
	
  for(int ref=1;ref<=1000;ref++)
    {	
		
      Events=((TProfile *)f->Get(Form("F%d%d_Cent%d",1,0,cent)))->GetBinEntries(ref+1);

      if(Events<=10) continue;
		
  
      for(int i=0; i<=8; i++)
	{
	  	
	  for(int j=0;j<=8;j++)
	    {
	  
	      if((i+j)>8) continue;
	  	  
	      F[i][j]=((TProfile *)f->Get(Form("F%d%d_Cent%d",i,j,cent)))->GetBinContent(ref+1);
		 
		 
	      for(int k=0;k<=8;k++)
	  	{
	  		
		  for(int h=0;h<=8;h++)
		    {
	  			
		      if((k+h)>8) continue;	
		      if((i+j+k+h)>8) continue;
	  			
	  				
		      ff[i][j][k][h]=((TProfile *)f4->Get(Form("f%d%d%d%d_Cent%d",i,j,k,h,cent)))->GetBinContent(ref+1);
	  				
	  			
	  					 
	  			 
	  			 
		    }
	  			
		}
		 
		  
		 
	    }
	}
   
  
   
      TotalEvents+=Events; 
  
      double meanV=0;
      // double var=0;
      if(strcmp(name,"netp")==0)
	{
	  meanV=F[1][0]-F[0][1];
	}

      if(strcmp(name,"pro")==0)
	{
	  meanV=F[1][0];
	}

      if(strcmp(name,"apro")==0)
	{
	  meanV=F[0][1];
	}

 
      double u2=CM(2, F, name);
      // var=CM(2,F,name);
      if(u2<=0) continue;
	  
 
				
     
    
      if(strcmp(name2,"C1")==0)
	{
	  //	double c1_err=sqrt(var)/sqrt(Events);
		
	  //	double var_c1=var/Events;
	  //	Cum_err+=var_c1*pow(Events,2);
	  if(strcmp(name,"netp")==0)
	    {
	      D[1][0][0][0]=1;
	      D[0][1][0][0]=1;
	      D[0][0][1][0]=-1;	
	      D[0][1][0][1]=-1;	
	    }
			
	  if(strcmp(name,"pro")==0)
	    {
	      D[1][0][0][0]=1;
	      D[0][1][0][0]=1;
			
	    }
		
	  if(strcmp(name,"apro")==0)
	    {
		
	      D[0][0][1][0]=1;	
	      D[0][1][0][1]=1;	
	    }	
				
			
	}
	
      if(strcmp(name2,"C2")==0)
	{
	  //			double c2_err=sqrt(var*(NCM[4]-1))/sqrt(Events)/2;
	  //		double var_c2=var*(NCM[4]-1)/Events/4.;
	  //		Cum_err+=var_c2*pow(Events,2);
				
	
		
	  for(int i=0; i<=8; i++)
	    {
	  	
	      for(int j=0;j<=8;j++)
		{
	
		  if((i+j)>8)continue;
	 	
		  for(int k=0;k<=8;k++)
		    {
		      if((j+k)>8)continue;
	  		 	
		      for(int h=0;h<=8;h++)
	  		{
	  			
			  if((k+h)>8)	continue;
			  if((i+j+k+h)>8) continue;
	  			
	  				
			  D[i][j][k][h]=DD(2,meanV,F,i,j,k,h);
	  			
				  
					
			}
		    }
		}
	    }
			
	}
	
      if(strcmp(name2,"C3")==0)
	{
	  //	double var_c3=pow(var,3)*(NCM[6] -NCM[3]*NCM[3]- 6*(NCM[4]-3) - 9)/Events;

				
	  for(int i=0; i<=8; i++)
	    {
	  	
	      for(int j=0;j<=8;j++)
		{
	
		  if((i+j)>8)continue;
	 	
		  for(int k=0;k<=8;k++)
		    {
		      if((j+k)>8)continue;
	  		 	
		      for(int h=0;h<=8;h++)
	  		{
	  			
			  if((k+h)>8)	continue;
			  if((i+j+k+h)>8) continue;
	  			
	  				
			  D[i][j][k][h]=DD(3,meanV,F,i,j,k,h);
	  			
				  
					
			}
		    }
		}
	    }	
			
	} 
	
      if(strcmp(name2,"C4")==0)
	{
			
	  //	double u2=CM(2, F, name);	
		
	  for(int i=0; i<=8; i++)
	    {
	  	
	      for(int j=0;j<=8;j++)
		{
	
		  if((i+j)>8)continue;
	 	
		  for(int k=0;k<=8;k++)
		    {
		      if((j+k)>8)continue;
	  		 	
		      for(int h=0;h<=8;h++)
	  		{
	  			
			  if((k+h)>8)	continue;
			  if((i+j+k+h)>8) continue;
	  			
	  				
			  D[i][j][k][h]=DD(4,meanV,F,i,j,k,h)-6*u2*DD(2,meanV,F,i,j,k,h);
	  			
				  
					
			}
		    }
		}
	    }



			
	  //double var_c4=pow(var,4)*(99 + 42*(NCM[4]-3) - (NCM[4]-3)*(NCM[4]-3) - 12*NCM[6] + NCM[8] - 8*NCM[5]*NCM[3] + 64*NCM[3]*NCM[3])/Events;
			
    
    
    
	}
	
	
      if(strcmp(name2,"KV")==0)
	{
	
	  //double u2=CM(2, F, name);
	  double u4=CM(4, F, name);
	
	
	
	  for(int i=0; i<=8; i++)
	    {
	  	
	      for(int j=0;j<=8;j++)
		{
	
		  if((i+j)>8)continue;
	 	
		  for(int k=0;k<=8;k++)
		    {
		      if((j+k)>8)continue;
	  		 	
		      for(int h=0;h<=8;h++)
	  		{
	  			
			  if((k+h)>8)	continue;
			  if((i+j+k+h)>8) continue;
	  			
	  				
			  D[i][j][k][h]=DD(4,meanV,F,i,j,k,h)/u2-(u4/pow(u2,2.)+3)*DD(2,meanV,F,i,j,k,h);
	  			
			  //cout<<i<<","<<j<<","<<k<<","<<h<<","<<D[i][j][k][h]<<endl;			  
					
			} 
		    }
		}
	    }

	
	}
	
      if(strcmp(name2,"SD")==0)
	{
	  //		double SD_err=sqrt(var)*sqrt((-9 + NCM[6]-2*NCM[5]*NCM[3] + 9*NCM[3]*NCM[3] + (NCM[4]-3)*(-6 + NCM[3]*NCM[3])))/sqrt(Events);
	  //			double var_SD=var*((-9 + NCM[6]-2*NCM[5]*NCM[3] + 9*NCM[3]*NCM[3] + (NCM[4]-3)*(-6 + NCM[3]*NCM[3])))/Events;
	  //		Cum_err+=var_SD*pow(Events,2);
				
	  // double u2=CM(2, F, name);
	  double u3=CM(3, F, name);
	
	
	
	  for(int i=0; i<=8; i++)
	    {
	  	
	      for(int j=0;j<=8;j++)
		{
	
		  if((i+j)>8)continue;
	 	
		  for(int k=0;k<=8;k++)
		    {
		      if((j+k)>8)continue;
	  		 	
		      for(int h=0;h<=8;h++)
	  		{
	  			
			  if((k+h)>8)	continue;
			  if((i+j+k+h)>8) continue;
	  			
	  				
	  					  				
			  D[i][j][k][h]=DD(3,meanV,F,i,j,k,h)/u2-(u3/pow(u2,2))*DD(2,meanV,F,i,j,k,h);
	  			
				  
					
			}
		    }
		}
	    }
	
	
			
	}
	
	
	
      if(strcmp(name2,"VM")==0)
	{
	  //		double VM_err=sqrt(var*(NCM[4]-1))/sqrt(Events)/2/MeanV;
	  //		double var_VM=var*(NCM[4]-1)/Events/4/pow(Cumulant(1,F,name),2);
	  //		Cum_err+=var_VM*pow(Events,2);
		
	  //double u2=CM(2, F, name);
	
	
	
	  for(int i=0; i<=8; i++)
	    {
	  	
	      for(int j=0;j<=8;j++)
		{
	
		  if((i+j)>8)continue;
	 	
		  for(int k=0;k<=8;k++)
		    {
		      if((j+k)>8)continue;
	  		 	
		      for(int h=0;h<=8;h++)
	  		{
	  			
			  if((k+h)>8)	continue;
			  if((i+j+k+h)>8) continue;
	  			
			  if(meanV==(F[1][0]-F[0][1]))
			    {
	  					  				
			      D[i][j][k][h]=DD(2,meanV,F,i,j,k,h)/meanV-(u2/pow(meanV,2.))*Delta(1,i+j+k+h)*pow(-1.,k+h);
			    }
	  			    
			  if((meanV==F[1][0]))
			    {
	  					  				
			      D[i][j][k][h]=DD(2,meanV,F,i,j,k,h)/meanV-(u2/pow(meanV,2.))*Delta(1,i+j+k+h)*Delta(1,i+j);
			    }
	  			    
			  if((meanV==F[0][1]))
			    {
	  					  				
			      D[i][j][k][h]=DD(2,meanV,F,i,j,k,h)/meanV-(u2/pow(meanV,2.))*Delta(1,i+j+k+h)*Delta(1,k+h);
			    }
	  			     
	  			 
					
			}
		    }
		}
	    }
	
		
		
				
			
	}
	
      if(strcmp(name2,"KDS")==0)
	{
			
	  //	double u2=CM(2, F, name);
	  double u3=CM(3, F, name);
	  double u4=CM(4, F, name);
		
	  // var_KDS=var*((-6*pow(NCM[4]-3,3) + 2*(NCM[4]-3)*(9*NCM[5] - NCM[7])*NCM[3] + (NCM[4]-3)*(NCM[4]-3)*(-9 + 
	  //NCM[6] + 8*NCM[3]*NCM[3])+NCM[3]*NCM[3]*(99 - 12*NCM[6] + NCM[8] - 8*NCM[5]*NCM[3] + 64*NCM[3]*NCM[3])))/pow(NCM[3],2)/pow(NCM[3],4)/Events;

			
	  for(int i=0; i<=8; i++)
	    {
	  	
	      for(int j=0;j<=8;j++)
		{
	
		  if((i+j)>8)continue;
	 	
		  for(int k=0;k<=8;k++)
		    {
		      if((j+k)>8)continue;
	  		 	
		      for(int h=0;h<=8;h++)
	  		{
	  			
			  if((k+h)>8)	continue;
			  if((i+j+k+h)>8) continue;
	  			
	  				
	  					  				
			  D[i][j][k][h]=DD(4,meanV,F,i,j,k,h)/u3-(u4-3*pow(u2,2.))*DD(3,meanV,F,i,j,k,h)/pow(u3,2.)-6*u2*DD(2,meanV,F,i,j,k,h)/u3;
	  			
				  
					
			}
		    }
		}
	    }
				
			
	}
		
		
		
      //***********************************		
	
      double temp_var=0;	
	
      for(int i=0; i<=8; i++)
	{
	  	
	  for(int j=0;j<=8;j++)
	    {
	  	
	      if((i+j)>8)continue;
	
	      for(int k=0;k<=8;k++)
	  	{
	  		
		  if((j+k)>8) continue;
	  		
	  		
		  for(int h=0;h<=8;h++)
		    {
	  			
		      if((k+h)>8) continue;	
		      if((i+j+k+h)>8) continue;
		      if(D[i][j][k][h]==0) continue;	
	  		
	  	
		      for(int r=0; r<=8; r++)
			{


			  if((r+h)>8) continue;
	  	
			  for(int s=0;s<=8;s++)
			    {
	
			      if((r+s)>8) continue;
	
			      for(int u=0;u<=8;u++)
				{
	  		
				  if((u+s)>8) continue;
	  		
				  for(int v=0;v<=8;v++)
				    {
	  			
				      if((u+v)>8) continue;
				      if((r+s+u+v)>8) continue;
				      if(D[r][s][u][v]==0) continue;
				      //double eff=1/(pow(fun_plow->Eval(ref),i+r)*pow(fun_phigh->Eval(ref),j+s)*pow(fun_aplow->Eval(ref),k+u)*pow(fun_aphigh->Eval(ref),h+v));
				      double eff=1/(pow(averEff_low[cent],i+r)*pow(averEff_low[cent],j+s)*pow(averEff_low[cent],k+u)*pow(averEff_low[cent],h+v));
	  			
	  			
				      temp_var+=eff*D[i][j][k][h]*D[r][s][u][v]*CovFact4(i,j,k,h,r,s,u,v,ff)/Events;
	  		
				    }
				}
			    }
	  		}
		    }
		}
	    }
	}

	

      //	cout<<ref<<","<<Events<<","<<CM(2, F, name)<<","<<temp_var<<endl;
      if(temp_var>1.e5) continue;
      Cum_err+=temp_var*pow(Events,2);
	
    }
	
  //	cout<<Cum_err<<endl;
  //p E	cout<<sqrt(Cum_err)/TotalEvents<<endl;
  cout<<"TotalEvents"<<":"<<TotalEvents<<endl;
  cout<<"The errors of "<<name<<":"<<sqrt(Cum_err)/TotalEvents<<endl;
  return sqrt(Cum_err)/TotalEvents;
	
	
	
}


double Cal_IndProd_base(const char *name, int n, TFile *f, int cent)
{
  double F[10][10];
  double Events=0;
  double TotalEvents=0;
  double Cum[10]={0};
	
	
	
  for(int i=0; i<=9; i++)
    {
	  	
      for(int j=0;j<=9;j++)
	{
	  	
	  F[i][j]=0;
	  	
	}
	  
    }
	
	
  for(int ref=1;ref<=1000;ref++)
    {	
		
      Events=((TProfile *)f->Get(Form("F%d%d_Cent%d",1,0,cent)))->GetBinEntries(ref+1);

      if(Events<=5) continue;
		
  
      for(int i=0; i<=8; i++)
	{
	  	
	  for(int j=0;j<=8;j++)
	    {
	  
	      if((i+j)<=8)
	  	{  
		  F[i][j]=(((TProfile *)f->Get(Form("F%d%d_Cent%d",i,0,cent)))->GetBinContent(ref))*(((TProfile *)f->Get(Form("F%d%d_Cent%d",0,j,cent)))->GetBinContent(ref));
		}
	    }
	}
   
      TotalEvents+=Events; 
  
    
 
		
 	
 
      for(int i=1;i<=4;i++)
	{
 			
 			
	  Cum[i]+=Cumulant(i,F,name)*Events;
 			
	}
 		
 		
 		   
   
   
   
   
   
	
    }

  if(n==0) return 1;
  if(n==1) return Cum[1]/TotalEvents;
  if(n==2) return Cum[2]/TotalEvents;	
  if(n==3) return Cum[3]/TotalEvents;
  if(n==4) return Cum[4]/TotalEvents;	
	
	
  return 0;	
	
	
	
}

double Moments_IndProd_Errors(const char *name, const char *name2, TFile *f, int cent)
{
	
  //	cout<<"Centrality: "<<cent<<", Invoking IndProd errors: "<<name<<","<<name2<<endl;

	
  double F[10][10];
  double Events=0;
  double TotalEvents=0;
  double NCM[10]={0};
  double Cum_err=0;
	
	
	
  for(int i=0; i<=9; i++)
    {
	  	
      for(int j=0;j<=9;j++)
	{
	  	
	  F[i][j]=0;
	  	
	}
	  
    }
	
	
  for(int ref=1;ref<=1000;ref++)
    {	
		
      Events=((TProfile *)f->Get(Form("F%d%d_Cent%d",1,0,cent)))->GetBinEntries(ref+1);

      if(Events<=5) continue;
		
  
      for(int i=0; i<=8; i++)
	{
	  	
	  for(int j=0;j<=8;j++)
	    {
	  
	      if((i+j)<=8)
	  	{  
		  F[i][j]=(((TProfile *)f->Get(Form("F%d%d_Cent%d",i,0,cent)))->GetBinContent(ref))*(((TProfile *)f->Get(Form("F%d%d_Cent%d",0,j,cent)))->GetBinContent(ref));
		}
	    }
	}
   
   
      TotalEvents+=Events; 
  
   
 
 
      double var=0;
 
  
 
      var=CM(2,F,name);
      if(var<=0) continue;
	  
 
      for(int i=3;i<=8;i++)
	{
 			
 			
	  NCM[i]=CM(i,F,name)/pow(sqrt(var),i);
 						
 			
	}
 		 		
 		 
    
      if(strcmp(name2,"C1")==0)
	{
	  //	double c1_err=sqrt(var)/sqrt(Events);
	  double var_c1=var/Events;
	  Cum_err+=var_c1*pow(Events,2);
				
			
	}
	
      if(strcmp(name2,"C2")==0)
	{
	  //		double c2_err=sqrt(var*(NCM[4]-1))/sqrt(Events)/2;
	  double var_c2=var*(NCM[4]-1)/Events/4.;
	  Cum_err+=var_c2*pow(Events,2);
				
			
	}
	
      if(strcmp(name2,"C3")==0)
	{
	  //		double c3_err=pow(var,3/2.)*sqrt(NCM[6] -NCM[3]*NCM[3]- 6*(NCM[4]-3) - 9)/sqrt(Events);
	  double var_c3=pow(var,3)*(NCM[6] -NCM[3]*NCM[3]- 6*(NCM[4]-3) - 9)/Events;
	  Cum_err+=var_c3*pow(Events,2);
				
			
	}
	
      if(strcmp(name2,"C4")==0)
	{
			
	  //	if(var<=0||(99 + 42*(NCM[4]-3) - (NCM[4]-3)*(NCM[4]-3) - 12*NCM[6] + NCM[8] - 8*NCM[5]*NCM[3] + 64*NCM[3]*NCM[3])<0) continue;
	  //		double c4_err=pow(var,2)*sqrt(99 + 42*(NCM[4]-3) - (NCM[4]-3)*(NCM[4]-3) - 12*NCM[6] + NCM[8] - 8*NCM[5]*NCM[3] + 64*NCM[3]*NCM[3])/sqrt(Events);
	  //	cout<<Events<<","<<99 + 42*(NCM[4]-3) - (NCM[4]-3)*(NCM[4]-3) - 12*NCM[6] + NCM[8] - 8*NCM[5]*NCM[3] + 64*NCM[3]*NCM[3]<<endl;
	  double var_c4=pow(var,4)*(99 + 42*(NCM[4]-3) - (NCM[4]-3)*(NCM[4]-3) - 12*NCM[6] + NCM[8] - 8*NCM[5]*NCM[3] + 64*NCM[3]*NCM[3])/Events;
	  Cum_err+=var_c4*pow(Events,2);
				
			
	}
	
	
      if(strcmp(name2,"KV")==0)
	{
	  //		double KV_err=var*sqrt(NCM[8]-2*(NCM[4]+3)*NCM[6]-8*NCM[3]*NCM[5]+8*pow(NCM[3],2)*(NCM[4]+5)+(pow(NCM[4]-3,3)+15*pow(NCM[4]-3,2)+72*(NCM[4]-3)+99))/sqrt(Events);
	  double var_KV=pow(var,2)*(NCM[8]-2*(NCM[4]+3)*NCM[6]-8*NCM[3]*NCM[5]+8*pow(NCM[3],2)*(NCM[4]+5)+(pow(NCM[4]-3,3)+15*pow(NCM[4]-3,2)+72*(NCM[4]-3)+99))/Events;
	  Cum_err+=var_KV*pow(Events,2);
				
	  //	cout<<ref<<","<<NCM[8]-2*(NCM[4]+3)*NCM[6]-8*NCM[3]*NCM[5]+8*pow(NCM[3],2)*(NCM[4]+5)+(pow(NCM[4]-3,3)+15*pow(NCM[4]-3,2)+72*(NCM[4]-3)+99)<<endl;
	  //	 cout<<Events<<","<<NCM[8]<<","<<NCM[7]<<","<<NCM[6]<<","<<NCM[5]<<","<<NCM[4]<<","<<NCM[3]<<","<<var<<endl;
	
	  // cout<<Events<<","<<var<<endl;
	
	}
	
      if(strcmp(name2,"SD")==0)
	{
	  //		double SD_err=sqrt(var)*sqrt((-9 + NCM[6]-2*NCM[5]*NCM[3] + 9*NCM[3]*NCM[3] + (NCM[4]-3)*(-6 + NCM[3]*NCM[3])))/sqrt(Events);
	  double var_SD=var*((-9 + NCM[6]-2*NCM[5]*NCM[3] + 9*NCM[3]*NCM[3] + (NCM[4]-3)*(-6 + NCM[3]*NCM[3])))/Events;
	  Cum_err+=var_SD*pow(Events,2);
				
			
	}
	
	
	
      if(strcmp(name2,"VM")==0)
	{
	  //			double VM_err=sqrt(var*(NCM[4]-1))/sqrt(Events)/2/MeanV;
	  double var_VM=var*(NCM[4]-1)/Events/4/pow(Cumulant(1,F,name),2);
	  Cum_err+=var_VM*pow(Events,2);
				
			
	}
	
      if(strcmp(name2,"KDS")==0)
	{
			
			
	  //			double KDS_err=sqrt(var)*sqrt((-6*pow(NCM[4]-3,3) + 2*(NCM[4]-3)*(9*NCM[5] - NCM[7])*NCM[3] + (NCM[4]-3)*(NCM[4]-3)*(-9 + 
	  //			NCM[6] + 8*NCM[3]*NCM[3])+NCM[3]*NCM[3]*(99 - 12*NCM[6] + NCM[8] - 8*NCM[5]*NCM[3] + 64*NCM[3]*NCM[3])))/pow(NCM[3],2)/sqrt(Events);
			
	  double var_KDS=var*((-6*pow(NCM[4]-3,3) + 2*(NCM[4]-3)*(9*NCM[5] - NCM[7])*NCM[3] + (NCM[4]-3)*(NCM[4]-3)*(-9 + 
														     NCM[6] + 8*NCM[3]*NCM[3])+NCM[3]*NCM[3]*(99 - 12*NCM[6] + NCM[8] - 8*NCM[5]*NCM[3] + 64*NCM[3]*NCM[3])))/pow(NCM[3],2)/pow(NCM[3],4)/Events;

	  Cum_err+=var_KDS*pow(Events,2);
				
			
	}
		
	
	
	
    }
	
  //	cout<<Cum_err<<endl;
  //p E	cout<<sqrt(Cum_err)/TotalEvents<<endl;
  return sqrt(Cum_err)/TotalEvents;
	
	
}



double Cal_Poisson_base(const char *name, int n, TFile *f, int cent)
{
	
  double F[10][10];
  double Events=0;
  double TotalEvents=0;
  double Cum[10]={0};

	
	
  for(int i=0; i<=9; i++)
    {
	  	
      for(int j=0;j<=9;j++)
	{
	  	
	  F[i][j]=0;
	  	
	}
	  
    }
	
	
  for(int ref=1;ref<=1000;ref++)
    {	
		
      Events=((TProfile *)f->Get(Form("F%d%d_Cent%d",1,0,cent)))->GetBinEntries(ref+1);

      if(Events<=5) continue;
		
  
      for(int i=0; i<=8; i++)
	{
	  	
	  for(int j=0;j<=8;j++)
	    {
	  
	      if((i+j)<=8)
	  	{  
		  F[i][j]=((TProfile *)f->Get(Form("F%d%d_Cent%d",i,j,cent)))->GetBinContent(ref+1);
		}
	    }
	}
   
 
   
      TotalEvents+=Events; 
  
      double MeanV=F[1][0]-F[0][1];
      double TotalV=F[1][0]+F[0][1];
   
      if(strcmp("pro",name)==0)
 	{
 		
	  Cum[1]+=F[1][0]*Events;
	  Cum[2]+=F[1][0]*Events;
	  Cum[3]+=F[1][0]*Events;
 	  Cum[4]+=F[1][0]*Events;
 		
 		
 		
 	}
   
      if(strcmp("apro",name)==0)
 	{
 		
	  Cum[1]+=F[0][1]*Events;
	  Cum[2]+=F[0][1]*Events;
	  Cum[3]+=F[0][1]*Events;
 	  Cum[4]+=F[0][1]*Events;
 		
 		
 		
 	}
     
      if(strcmp("netp",name)==0)
 	{
 		 


	  Cum[1]+=MeanV*Events;
 		
	  Cum[2]+=TotalV*Events;
 		
	  Cum[3]+=MeanV*Events;
 	  
 	  Cum[4]+=TotalV*Events;
 		
 		
 	}
     
      if(strcmp("totp",name)==0)
 	{
 	
		
	  Cum[1]+=TotalV*Events;
 		
	  Cum[2]+=TotalV*Events;
 		
	  Cum[3]+=TotalV*Events;
 	  
 	  Cum[4]+=TotalV*Events;
 		
 		
 	}
       
   
   
   
   
   
	
    }

  if(n==0) return 1;
  if(n==1) return Cum[1]/TotalEvents;
  if(n==2) return Cum[2]/TotalEvents;	
  if(n==3) return Cum[3]/TotalEvents;
  if(n==4) return Cum[4]/TotalEvents;	
  return 0;	
	
}



double Cal_NBD_base(const char *name, int n, TFile *f, int cent)
{
  double F[10][10];
  double Events=0;
  double TotalEvents=0;
  double Cum[10]={0};
	
	

	
	
  for(int i=0; i<=9; i++)
    {
	  	
      for(int j=0;j<=9;j++)
	{
	  	
	  F[i][j]=0;
	  	
	}
	  
    }
	
	
  for(int ref=1;ref<=1000;ref++)
    {	
		
      Events=((TProfile *)f->Get(Form("F%d%d_Cent%d",1,0,cent)))->GetBinEntries(ref+1);

      if(Events<=5) continue;
		
  
      for(int i=0; i<=8; i++)
	{
	  	
	  for(int j=0;j<=8;j++)
	    {
	  
	      if((i+j)<=8)
	  	{  
		  F[i][j]=((TProfile *)f->Get(Form("F%d%d_Cent%d",i,j,cent)))->GetBinContent(ref+1);
		}
	    }
	}
   
  
      double eps_p=(F[2][0]+F[1][0]-F[1][0]*F[1][0])/F[1][0];
      double eps_pbar=(F[0][2]+F[0][1]-F[0][1]*F[0][1])/F[0][1];
      double MeanV=F[1][0]-F[0][1];
      double TotalV=F[1][0]+F[0][1];
   
      if(F[1][0]==0||F[0][1]==0) continue;
   
      TotalEvents+=Events; 
   
   
      if(strcmp("pro",name)==0)
 	{
 		
 		
 		
	  Cum[1]+=F[1][0]*Events;
	  Cum[2]+=(F[2][0]+F[1][0]-F[1][0]*F[1][0])*Events;
 		
 		
	  Cum[3]+=eps_p*F[1][0]*(2*eps_p-1)*Events;
 	  Cum[4]+=eps_p*F[1][0]*(6*eps_p*eps_p-6*eps_p+1)*Events;
 		
 		
 		
 	}
   
      if(strcmp("apro",name)==0)
 	{
 		
	  Cum[1]+=F[0][1]*Events;
	  Cum[2]+=(F[0][2]+F[0][1]-F[0][1]*F[0][1])*Events;
 		
 		
	  Cum[3]+=eps_pbar*F[0][1]*(2*eps_pbar-1)*Events;
 	  Cum[4]+=eps_pbar*F[0][1]*(6*eps_pbar*eps_pbar-6*eps_pbar+1)*Events;
 		
 		
 		
 	}
     
      if(strcmp("netp",name)==0)
 	{
 		 
 		 


	  Cum[1]+=MeanV*Events;
 		
	  Cum[2]+=(eps_p*F[1][0]+eps_pbar*F[0][1])*Events;
 		
	  Cum[3]+=(eps_p*F[1][0]*(2*eps_p-1)-eps_pbar*F[0][1]*(2*eps_pbar-1))*Events;
 	  
 	  Cum[4]+=(eps_p*F[1][0]*(6*eps_p*eps_p-6*eps_p+1)+eps_pbar*F[0][1]*(6*eps_pbar*eps_pbar-6*eps_pbar+1))*Events;
 		
 		
 	}
     
      if(strcmp("totp",name)==0)
 	{
 	
		
	  Cum[1]+=TotalV*Events;
 		
	  Cum[2]+=(eps_p*F[1][0]+eps_pbar*F[0][1])*Events;
 		
	  Cum[3]+=(eps_p*F[1][0]*(2*eps_p-1)+eps_pbar*F[0][1]*(2*eps_pbar-1))*Events;
 	  
 	  Cum[4]+=(eps_p*F[1][0]*(6*eps_p*eps_p-6*eps_p+1)+eps_pbar*F[0][1]*(6*eps_pbar*eps_pbar-6*eps_pbar+1))*Events;
 		
 		
 	}
       
   
   
   
   
   
	
    }

  if(n==0) return 1;
  if(n==1) return Cum[1]/TotalEvents;
  if(n==2) return Cum[2]/TotalEvents;	
  if(n==3) return Cum[3]/TotalEvents;
  if(n==4) return Cum[4]/TotalEvents;	
	
	
  return 0;	
	
}






//*************************Useful function**********************************

double StirlingS2(int n,int k)
{
	
  if(k>n) { return 0;}
  if(k==n)  {return 1;}
  if(k<n)
    {	
  
      if(k==0) return 0;
      if(k==1) return 1;
   	
      return  (StirlingS2(n-1,k-1)+k*StirlingS2(n-1,k));	
  
    }
  	
  return 0; 	
 
}




double UnsignS1(int n,int k)
{
	 
  if(k>n){ return 0;}
  if(k==n) {return 1;}
  if(k<n)
    {	
      if(k==0) {return 0;}
  	
      return (UnsignS1(n-1,k-1)+(n-1)*UnsignS1(n-1,k));	
   
   
   
    }
 
  return 0;
    
}




double StirlingS1(int n,int k)
{
	 
	    	
  return ((n-k)%2==0 ? 1 : -1 )*UnsignS1(n,k);
   
   
   
 
    
}



double lnchoose(int n, int m)                 
{                                         
  if (m > n)                            
    {                                     
      return 0;                         
    }                                     
  if (m < n/2.0)                        
    {                                     
      m = n-m;                          
    }                                     
  double s1 = 0;                        
  for (int i=m+1; i<=n; i++)            
    {                                     
      s1 += log((double)i);             
    }                                     
  double s2 = 0;                        
  int ub = n-m;                         
  for (int i=2; i<=ub; i++)             
    {                                     
      s2 += log((double)i);             
    }                                     
  return s1-s2;                         
}                                         
                                               
 
double C(int n, int m)               
{                                         
  if (m > n)                            
    {                                     
      return 0;                         
    }                                     
  return exp(lnchoose(n, m));           
}                                         





double Moments2(int m, int n, double ff[][10])
{
	
  double counts=0;
	
  for(int i=0; i<=m;i++)
    {
      for(int j=0; j<=n;j++)
	{
		
	  counts+=StirlingS2(m,i)*StirlingS2(n,j)*ff[i][j];
		
		
		
	}
	
	
    }


  return counts;

}


double Moments4(int m, int n, int s, int t, double ff[][10][10][10])
{
	
  double counts=0;
	
  for(int i=0; i<=m;i++)
    {
      for(int j=0; j<=n;j++)
	{
		
	  for(int k=0; k<=s;k++)
	    {
		
	      for(int h=0; h<=t;h++)
		{
		
		
		  counts+=StirlingS2(m,i)*StirlingS2(n,j)*StirlingS2(s,k)*StirlingS2(t,h)*ff[i][j][k][h];
		
		
		
		}
	
	
	    }
	}

    }	
  return counts;
	
	
	
}


double FM(int m, int n, double ff[][10][10][10])
{
  //	cout<<"Invoking FM..."<<endl;
  double counts=0;
	
		
  for(int i=0; i<=m;i++)
    {
			
			
      for(int j=0; j<=n;j++)
	{
	
		
	  for(int r=0;r<=i;r++)
	    {
	   	
	 
	   	
	      for(int s=0;s<=j;s++)
	   	{
	   		
	   		
	   		
	   		
	   		
		  counts+=StirlingS1(m,i)*StirlingS1(n,j)*C(i,r)*C(j,s)*Moments4(i-r,r,j-s,s,ff);
	   		
	   		
	   		
	   	}
	   	
	   	
	    }
		
		
	}
	
	
    }
	
  return counts;
	
	
	
}

double CM(int order, double FF[][10], const char *name)
{
	
  double counts=0;
  double MeanV=0;
  if(strcmp("pro",name)==0) 	MeanV=FF[1][0];
  if(strcmp("apro",name)==0) MeanV=FF[0][1];
  if(strcmp("netp",name)==0)	MeanV=FF[1][0]-FF[0][1]; 
  if(strcmp("totp",name)==0) MeanV=FF[1][0]+FF[0][1];
		
	
	
  if(strcmp("netp",name)==0)
    {
	
      for(int i=0;i<=order;i++)
	{
		
	  for(int r=0; r<=order-i; r++)
	    {
			
			
			
		
	      counts+=C(order,i)*pow(-1.,i)*pow(MeanV,i)*C(order-i,r)*pow(-1.,r)*Moments2(order-i-r,r,FF);
		
			
			
			
	    }
			
		
		
	}
	
    }



  else if(strcmp("totp",name)==0)
    {
	
	
      for(int i=0;i<=order;i++)
	{
		
	  for(int r=0; r<=order-i; r++)
	    {
			
			
			
		
	      counts+=C(order,i)*pow(-1.,i)*pow(MeanV,i)*C(order-i,r)*Moments2(order-i-r,r,FF);
		
			
			
			
	    }
			
		
		
	}
	
	
	
	
    }



  else if(strcmp("pro",name)==0)
    {
	
	
      for(int i=0;i<=order;i++)
	{
		
			
	  counts+=C(order,i)*pow(-1.,i)*pow(MeanV,i)*Moments2(order-i,0,FF);
		
	
	}
	
	
	
	
    }


  else 
    {
	
	
      for(int i=0;i<=order;i++)
	{
		
			
	  counts+=C(order,i)*pow(-1.,i)*pow(MeanV,i)*Moments2(0,order-i,FF);
		
	
	}
	
	
	
	
    }

	
  return counts;
	
	
	
	
}


double Cumulant(int n, double FF[][10], const char *name)
{
	
  double m1=0;
  double m2=0;
	
	
  if(strcmp("netp",name)==0)
    {
	
      for(int i=0;i<=n;i++)
	{
		
			
			
		
	  m1+=C(n,i)*pow(-1.,i)*Moments2(n-i,i,FF);
		
			
			
			
	}
		
		
		
		
      for(int m=1;m<=n-1;m++)
	{
		
			
	  for(int i=0;i<=n-m;i++)	
	    {
		
	      m2+=C(n-1,m-1)*Cumulant(m,FF,name)*C(n-m,i)*pow(-1.,i)*Moments2(n-m-i,i,FF);
	
	    }
			
			
			
	}
			
			
		

	
    }



  else if(strcmp("totp",name)==0)
    {
	
      for(int i=0;i<=n;i++)
	{
		
			
			
		
	  m1+=C(n,i)*Moments2(n-i,i,FF);
		
			
			
			
	}
		
		
		
		
      for(int m=1;m<=n-1;m++)
	{
		
			
	  for(int i=0;i<=n-m;i++)	
	    {
		
	      m2+=C(n-1,m-1)*Cumulant(m,FF,name)*C(n-m,i)*Moments2(n-m-i,i,FF);
	
	    }
			
			
			
	}
			
			
	
	
    }



  else if(strcmp("pro",name)==0)
    {
	

			
		
      m1=Moments2(n,0,FF);
	
		
      for(int m=1;m<=n-1;m++)
	{
		
			
		
		
	  m2+=C(n-1,m-1)*Cumulant(m,FF,name)*Moments2(n-m,0,FF);
	
	 
			
			
	}
			
			
	
	
	
	
    }


  else 
    {
		
      m1=Moments2(0,n,FF);
		
			
		
      for(int m=1;m<=n-1;m++)
	{
		
			
		
		
	  m2+=C(n-1,m-1)*Cumulant(m,FF,name)*Moments2(0,n-m,FF);
	
	 			
			
	}
			
	
	
	
	
    }

	
  return m1-m2;
	
	
	
	
}


double CovFact4( int i, int j, int k, int h, int r, int s, int u, int v,double FF[][10][10][10])
{


  double count=0;


  for(int x1=0;x1<=i;x1++)
    {
      		
      for(int x2=0;x2<=j;x2++)	
	{  	
      		
	  for(int x3=0;x3<=k;x3++)
	    {
	      for(int x4=0;x4<=h;x4++)
		{
      		
      		
		  for(int y1=0;y1<=r;y1++)
		    {
      				 
		      for(int y2=0;y2<=s;y2++)
			{ 
      				 	 
			  for(int y3=0;y3<=u;y3++)
      				 	 	 
			    {	 
			      for(int y4=0;y4<=v;y4++)
				{
      					   
      		 
      		
				  count+=StirlingS1(i,x1)*StirlingS1(j,x2)*StirlingS1(k,x3)*StirlingS1(h,x4)*StirlingS1(r,y1)*StirlingS1(s,y2)*StirlingS1(u,y3)*StirlingS1(v,y4)*CovMoment4(x1,x2,x3,x4,y1,y2,y3,y4,FF);
      		
      		
      		
				}
			    }
			}
		    }
		}
	    }
	}
    }


  return count;


}





double CovMoment4(int i, int j, int k, int h, int r, int s, int u, int v,double FF[][10][10][10])
{




  return Moments4(i+r,j+s,k+u,h+v,FF)-Moments4(i,j,k,h,FF)*Moments4(r,s,u,v,FF);






}





double D_CM_F2(int n, int u, int v, double mean,double F[][10])
{
	
  double count1=0;
	
	
  for(int r=0;r<=n;r++)	
    {
		
		
		
      count1+=C(n,r)*pow(-1.,r)*D_Moment_F2(n-r,u,v,mean,F)*pow(mean,r);
		
		
    }
	
	
  double count2=0;
	
  if((u+v)==1)
    {


      if(mean==(F[1][0]-F[0][1]))
	{	
	  for(int r=0;r<=n;r++)	
	    {
		
	      for(int k=0;k<=n-r;k++)
		{
		
		
		  count2+=pow(-1.,v)*C(n,r)*pow(-1.,r)*C(n-r,k)*pow(-1.,k)*Moments2(n-r-k,k,F)*r*pow(mean,r-1);
			
		}
			
	    }
	

	}




      if((mean==F[1][0]))
	{
	  for(int r=0;r<=n;r++)
	    {



	      count2+=C(n,r)*pow(-1.,r)*Moments2(n-r,0,F)*r*pow(mean,r-1);


	    }


	}


      if((mean==F[0][1]))
	{
	  for(int r=0;r<=n;r++)
	    {



	      count2+=C(n,r)*pow(-1.,r)*Moments2(0,n-r,F)*r*pow(mean,r-1);


	    }


	}






    }




  return count1+count2;



}




double D_Moment_F2(int n, int u, int v,double mean, double F[][10])
{
	
  double count=0;
  if(mean==(F[1][0]-F[0][1]))
    {
      for(int r=0;r<=n;r++)	
	{
		
		
	  count+=C(n,r)*pow(-1.,r)*StirlingS2(n-r,u)*StirlingS2(r,v);
		
		
	}
	
	
	
    }


  if((mean==F[1][0]))
    {


      count=StirlingS2(n,u);



    }

  if((mean==F[0][1]))
    {


      count=StirlingS2(n,v);



    }




  return count;


	
}


double D_F2_F4(int r1, int r2, int u, int v, int j, int k)
{
	
  double count=0;
	

  for(int i1=0;i1<=r1;i1++)
    {
		
      for(int i2=0;i2<=r2;i2++)
	{
			
	  for(int s=0;s<=i1;s++)
	    {
	
	      for(int t=0;t<=i2;t++)
		{
	
	
	
		  count+=StirlingS1(r1,i1)*StirlingS1(r2,i2)*C(i1,s)*C(i2,t)*StirlingS2(i1-s,u)*StirlingS2(s,v)*StirlingS2(i2-t,j)*StirlingS2(t,k);
	
	
		}


	    }
	}
 
    }

	
  return count;	
	
}

double DD(int n,double mean, double F[][10],int i, int j, int k, int h)
{

  double count=0;
	
	
  if(mean==(F[1][0]-F[0][1]))
    {
			
      for(int u=0;u<=n;u++)
	{
	  for(int v=0;v<=n;v++)
	    {
	 	
	      if((u+v)>n) continue;
	 		
	      count+=D_F2_F4(u,v,i,j,k,h)*D_CM_F2(n,u,v,mean,F);
		
	    }
	}
    } 

  if(mean==F[1][0])
    {
			
      for(int u=0;u<=n;u++)
	{
	  for(int v=0;v<=0;v++)
	    {
	 	
	      if((u+v)>n) continue;
	 		
	      count+=D_F2_F4(u,v,i,j,k,h)*D_CM_F2(n,u,v,mean,F);
		
	    }
	}
    }

  if((mean==F[0][1]))
    {
			
      for(int u=0;u<=0;u++)
	{
	  for(int v=0;v<=n;v++)
	    {
	 	
	      if((u+v)>n) continue;
	      //	cout<<"test"<<endl;
	 		
	      count+=D_F2_F4(u,v,i,j,k,h)*D_CM_F2(n,u,v,mean,F);
		
	    }
	}
    }

  return count;



}

int Delta(int a, int b)
{

  return (a==b) ? 1 : 0;


}


double ffff(int i, int j, int k, int h,TFile *f, int cent, int ref)
{
	
  return ((TProfile *)f->Get(Form("F%d%d%d%d_Cent%d",i,j,k,h,cent)))->GetBinContent(ref+1);
	
	
	
	
}



