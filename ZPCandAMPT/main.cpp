#include <fstream>
#include <iostream>
#include <exception>

#include "ZPCandAMPT.h"

using namespace std;

int main(int argc, char *argv[])
{
  ifstream input("./ampt/ana/zpc.dat");
  xiaohai::ZPCandAMPT zpc;
  input >> zpc;
  for (int i = 0; i != zpc.GetZPCandAMPTMulti(); ++i)
    {
      cout << zpc.GetZPCandAMPTId(i) << endl;
    }
  return 0;
}
