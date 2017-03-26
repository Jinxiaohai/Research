#ifndef PARAMETERMAP_CC
#define PARAMETERMAP_CC

#include "parametermap.h" 
#include <iostream>

//Returns an integer from the map.
int parameter::getI(parametermap m,string key,int def)
{
  int param;
  map<string,string>::iterator itr; 
  itr = m.find(key); 
  if(itr!=m.end()){
    stringstream ss(itr->second); 
    ss>>param; 
  }else{
    param = def; 
  }
  return param;
}

//Returns an bool from the map.
bool parameter::getB(parametermap m,string key,bool def)
{
  bool param;
  map<string,string>::iterator itr; 
  itr = m.find(key); 
  if(itr!=m.end()){
    stringstream ss(itr->second); 
    ss>>param; 
  }else{
    param = def; 
  }
  return param;
}

//Returns a string from the map.
string parameter::getS(parametermap m,string key,string def)
{
  string param;
  map<string,string>::iterator itr; 
  itr = m.find(key); 
  if(itr!=m.end()){
    param = itr->second;
  }else{
    param = def; 
  }
  return param;
}

//Returns a double from the map.
double parameter::getD(parametermap m,string key,double def)
{
  double param;
  map<string,string>::iterator itr; 
  itr = m.find(key); 
  if(itr!=m.end()){
    stringstream ss(itr->second); 
    ss>>param; 
  }else{
    param = def; 
  }
  return param;
}

//Returns a STL Vector from the map
vector< double > parameter::getV(parametermap m, string key, double def){
    vector< double > vec(0);
    double tmp;
    map<string,string>::iterator itr; 
    itr = m.find(key); 
    if(itr!=m.end()){
        stringstream ss(itr->second);
        while (ss>>tmp){vec.push_back(tmp);}
    }else{
        vec.push_back(def); 
    }
    return vec;

}

//Returns a STL Vector from the map
vector< string > parameter::getVS(parametermap m, string key, string def){
    vector< string > vec(0);
    string tmp;
    map<string,string>::iterator itr; 
    itr = m.find(key); 
    if(itr!=m.end()){
        stringstream ss(itr->second);
        while (ss>>tmp){vec.push_back(tmp);}
    }else{
        vec.push_back(def); 
    }
    return vec;

}

//Returns a matrix from the map
vector< vector< double > > parameter::getM(parametermap m, string key, double def){
    vector< vector< double > > mtx;
    double tmp;
    map<string,string>::iterator itr; 
    itr = m.find(key); 
    if(itr!=m.end()){
        stringstream ss(itr->second);
        string line("");
        while (ss.good()){
            vector< double > vec;
            while ((line=="")&&ss.good()) {
                getline(ss,line); 
                stringstream buf(line);
                while (buf>>tmp) {vec.push_back(tmp);}
            }
            line = "";
            if (vec.size()!=0) mtx.push_back(vec);
        }
    }else{
        vector< double > vec(1,0.); 
        mtx.push_back(vec); 
    }
    return mtx;
}

//Adds a double to the map.
void parameter::set(parametermap& m,string key,double val)
{
    string sval;
    stringstream ss; 
    ss<<val;ss>>sval; 
    set(m,key,sval);
}
//Adds an int to the map.
void parameter::set(parametermap& m,string key,int val)
{
    string sval;
    stringstream ss; 
    ss<<val;ss>>sval; 
    set(m,key,sval);
}
//Adds a bool to the map.
void parameter::set(parametermap& m,string key,bool val)
{
    string sval;
    stringstream ss; 
    ss<<val;ss>>sval; 
    set(m,key,sval);
}
//Adds a char* to the map.
void parameter::set(parametermap& m,string key,char* val)
{
    string sval(val);
    set(m,key,sval);
}
//Adds a string to the map.
void parameter::set(parametermap& m,string key,string val)
{
    map<string,string>::iterator itr; 
    itr = m.find(key); 
    if(itr!=m.end()){
        m[key]=val;
    } else {
        m.insert(make_pair(key,val));
    }
}

//Adds a vector to the map.
void parameter::set(parametermap& m, string key, vector< double > val){
    stringstream ss;
    for (int i=0;i<static_cast<int>(val.size());++i) {ss<<"    "<<val[i]<<"\n";}
    set(m,key,ss.str());
}

//Adds a vector to the map.
void parameter::set(parametermap& m, string key, vector< string > val){
    stringstream ss;
    for (int i=0;i<static_cast<int>(val.size());++i) {ss<<"    "<<val[i]<<"\n";}
    set(m,key,ss.str());
}

//Adds a matrix to the map.
void parameter::set(parametermap& m, string key, vector< vector< double > > val){
    stringstream ss;
    for (int i=0;i<static_cast<int>(val.size());++i) {
        ss << "    ";
        for (int j=0;j<static_cast<int>(val[i].size());++j) {ss<<val[i][j]<<" ";}
        ss << "\n";
    }
    set(m,key,ss.str());
}

// Read parameters from file in format "type  key value", e.g.,
// int nqmax 50
// A # sign denotes a comment line
void parameter::ReadParsFromFile(parametermap& m,char *filename){
  ifstream parsfile;
  char line[120];
  string type,key,value;
  stringstream ss;
  parsfile.open(filename);
  if(! parsfile){
    printf("attempting to read non-existent parameter file %s\n",filename);
    exit(1);
  }

  while(parsfile.getline(line,100)){
    if(line[0]!='#' && strlen(line)>4){
      ss << line;
      ss >> type >> key >> value;
      //cout << "type=" << type << ", key=" 
      //<< key << ", value=" << value << endl;
      if(type=="int"){
	set(m,key,atoi(value.c_str()));
      }
      if(type=="double"){
	set(m,key,atof(value.c_str()));
      }
      if(type=="bool"){
	set(m,key,bool(atoi(value.c_str())));
      }
      if(type=="string"){
	set(m,key,value);
      }
      ss.clear();
    }
  }
  parsfile.close();
}

void parameter::PrintPars(parametermap& pars){
  map<string,string>::iterator itr;
  for(itr=pars.begin(); itr !=pars.end(); ++itr){
    cout << itr->first << " = " << itr->second << endl;
  }
}

#endif
