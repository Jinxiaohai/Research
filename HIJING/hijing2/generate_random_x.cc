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
#include <time.h>

using namespace std;
// int main()
int main(void)
{
 int nEvt=10000;
 double bMin=0;
 double bMax=5;
 string f="InputFiles/input";
 ofstream myfile; 
 srand(time(NULL)); 
 // int r = rand();
 for(int i=0;i<nEvt;i++){
   std::stringstream fout;
   fout << f.data() << "_p" << i << ".txt";  
   myfile.open(fout.str().c_str());
   cout << "file " << fout.str() << " opened" << endl;
   // double RAND_MAX=1;
   double x=1;
   // while (x>0.99){
    x=((double) rand() / (RAND_MAX));
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
 }
 cout << "Done!" << endl;
 return 0;
}