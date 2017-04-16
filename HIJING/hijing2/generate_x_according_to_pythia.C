// #include <algorithm> 
// #include <functional> 
// #include <cctype>
// #include <locale>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
// #include <vector>
// #include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// #include <time.h>
#include "ntupler/nt.C"
#include "ntupler/hi.C"

void generate_x_according_to_pythia()
{
 int nEvt=10000;
 double bMin=0;
 double bMax=5;
 string f="InputFilesWPYTHIA/input";
 ofstream myfile; 
 // srand(time(NULL)); 
 // int r = rand();
 nt * fdijet = new nt("/afs/cern.ch/user/d/dgulhan/workDir/pythiaJets2/merged.root");
 hi * fgen = new hi("/afs/cern.ch/user/d/dgulhan/workDir/pythiaJets2/merged.root");
 int nentries = fdijet->GetEntriesFast();  
 // int nentries = fgen->GetEntriesFast(); 
 int ievent=0;
 int i=0;
 // int nentries=100;
 while(ievent<nEvt && i<nentries){
 fgen->GetEntry(i); 
 fdijet->GetEntry(i); 
 // for(int i=0;i<nEvt;i++){
   if(fdijet->pt1>120 && fdijet->pt2>30 && acos(cos(fdijet->phi1-fdijet->phi2))>2.094){
    std::stringstream fout;
    fout << f.data() << "_p" << ievent << ".txt";  
    myfile.open(fout.str().c_str());
    cout << "file " << fout.str() << " opened" << endl;
    // double RAND_MAX=1;
    double x=fgen->x1;
    // double x=0;
    // while (x>0.99){
     // x=((double) rand() / (RAND_MAX));
    // }
    double EFRM=5020*sqrt(1-x);
    myfile << "32,"<<i<<"\n";
    myfile << "32,CMS,"<<EFRM <<"\n";
    myfile << "32,P,A\n";
    myfile << "32,1,0,208,82\n";
    myfile << "32,1\n";
    myfile << "32,1\n";
    myfile << "32,0\n";
    myfile << "32,"<<bMin<<","<<bMax<<"\n";
    myfile << "x="<<x<<"\n";
    myfile.close();
    ievent++;
   }
   i++;
 }
 cout << "Done!" << endl;
}