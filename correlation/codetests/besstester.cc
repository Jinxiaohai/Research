#include <iostream>
#include <cstdlib>
#include <cmath>

#include "sf.cc"
#include "misc.cc"

const double PI=3.14159265358979323844;
using namespace std;

int main(){
  int n;
  double z,ans;
  
 TRYNEW:
  cout << "What is z?\n";
  cin >> z;
  ans=Bessel::j0(z);
  cout << "j0(" << z << ")= " << ans << "=" << sin(z)/z << endl;
  ans=Bessel::j1(z);
  cout << "j1(" << z << ")= " << ans << "=" << -cos(z)/z+sin(z)/(z*z) << endl;


  cout << "What is n? (Choose n>0)";
  cin >> n;
  ans=Bessel::jn(n,z);
  cout << "j(" << n << "," << z << ")= " << ans << endl;

  cout << "Enter 0 to quit, other to continue\n";
  cin >> n;
  if(n!=0) goto TRYNEW;

  return 0;
  
}

