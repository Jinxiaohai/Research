/********************************************************/
/*                 Created  by  xiaohai                 */
/*                 Telphone : 18501781924               */
/*            E-mail : jinxiaohai@sinap.ac.cn           */
/*            E-mail : xiaohaijin@outlook.com           */
/*   Address : Shanghai Institute of Applied Physics    */
/********************************************************/
/********************************************************/
/*
 * regex     : regular expression
 * iomanip   : manipulator(操作器)
 */
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <complex>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <regex>
#include <random>
#include <ctime>
#include "Particle.h"

using namespace std;
using namespace xiaohai;

int main(int argc, char *argv[])
{
  ThreeDimensionVector<double> vect1(1, 2, 3);
  cout << vect1.GetX() << "  " << vect1.GetY() << "  " << vect1.GetZ() << endl;
  cout << vect1.GetMag() << "  " << vect1.GetMag2() << "  " << vect1.GetPt() << endl;
  ThreeDimensionVector<double> vect2(3, 2, 1), vect3;
  vect3 = vect1 + vect2;
  vect3 = vect3 * 2.;
  cout << vect3.GetX() << "  " << vect3.GetY() << "  " << vect3.GetZ() << endl;
  cout << vect3.GetMag() << "  " << vect3.GetMag2() << "  " << vect3.GetPt() << endl;
  vect3 = -vect3;
  cout << vect3.GetX() << "  " << vect3.GetY() << "  " << vect3.GetZ() << endl;
  return 0;
}
