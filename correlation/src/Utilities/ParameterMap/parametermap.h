#ifndef PARAMETERMAP_H
#define PARAMETERMAP_H

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

//---------------------------------------------------------------
//This code helps one parse the map that contains configuration 
// parameters.
//
//Example:
// 
// MyCoolClass(parametermap m){  //m has the parameters and is passed in
// int importantParameter = parameter::getI(m,"nameOfParameter",-1);
//  ...
//
//  If "nameOfParameter is not a key in the map, -1 is returned since that
// is the 3rd argment of the getI function.
//
//MH 22 jun04
//---------------------------------------------------------------

//This code only works with a map of the type below.  The type def is
//to make it easy to remember.
typedef  map<string,string> parametermap;

//These functions are all in the namespace parameter.
namespace parameter {
  bool   getB(parametermap ,string ,bool);
  int    getI(parametermap ,string ,int);
  string getS(parametermap ,string ,string);
  double getD(parametermap ,string ,double);
  vector< double > getV(parametermap, string, double);
  vector< string > getVS(parametermap, string, string);
  vector< vector< double > > getM(parametermap, string, double);
  void set(parametermap&, string, double);
  void set(parametermap&, string, int);
  void set(parametermap&, string, bool);
  void set(parametermap&, string, string);
  void set(parametermap&, string, char*);
  void set(parametermap&, string, vector< double >);
  void set(parametermap&, string, vector< string >);
  void set(parametermap&, string, vector< vector< double > >);
  void ReadParsFromFile(parametermap&, char *filename);
  void PrintPars(parametermap&);
};

#endif

