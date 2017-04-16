#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

struct Particle
{
 private:
  int ind;
  int mom_ind;
  int pdg_id;
  int status;
  double px;
  double py;
  double pz;
  double E;
 public:
  Particle(int ind, int mom_ind, int pdg_id, int status, double px, double py, double pz, double E)
  {
   this->ind = ind;
   this->mom_ind = mom_ind;
   this->pdg_id = pdg_id;
   this->status = status;
   this->px = px;
   this->py = py;
   this->pz = pz;
   this->E = E;
  }
  int get_ind(){return ind;}
  int get_mom_ind(){return mom_ind;}
  int get_pdg_id(){return pdg_id;}
  int get_status(){return status;}
  double get_px(){return px;}
  double get_py(){return py;}
  double get_pz(){return pz;}
  double get_E(){return E;}
};

struct Event
{
 private:
  Particle **particles;//!SET BY READING FROM THE DATA FILE
  int number_of_particles;
  double b;
  double frame_E;
  double total_E;
  int n_stable;

 public:
  void reset()
  {
   particles = NULL;
   number_of_particles=0;
   frame_E=0;
   b=0;
   total_E=0;
   n_stable=0;
  }
  
  Event()
  {
   reset();
  }
  
  ~Event()
  {
   
  }
  void set_frame_E(double frame_E){this->frame_E=frame_E;}
  void set_b(double b){this->b=b;}
  void set_total_E(double total_E){this->total_E=total_E;}
  void set_n_stable(int n_stable){this->n_stable=n_stable;}
  double get_frame_E(){return frame_E;}
  double get_b(){return b;}
  double get_total_E(){return total_E;}
  int get_n_stable(){return n_stable;}
  int get_number_of_particles(){return number_of_particles;}
  void add_particle(Particle *particle)
  {
   if(particles == NULL)
   {
    particles = (Particle **)calloc(2,sizeof(Particle *));
    particles[0] = particle;
	  particles[1] = NULL;
   }
   else
   {
    particles = (Particle **)realloc(particles, (number_of_particles+2)*sizeof(Particle *));
	  particles[number_of_particles] = particle;
  	particles[number_of_particles+1] = NULL;
   }
   number_of_particles++;
  }
  Particle * get_ieth_particle(int i){return particles[i];}
};

// trim from start
static inline string &ltrim(string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(string &s) {
        s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline string &trim(string &s) {
        return ltrim(rtrim(s));
}

class HIJING_file_parser
{
 private:
 vector <Event> events;

 public:
 vector <Event> get_event_vector()
 {
  return events;
 } 

 // DataFile_Parser(const char *file_path)
 HIJING_file_parser()
 { 
 for(int ifile=1;ifile<10000;ifile++)
 {
  int ievent=0;
  bool hasjet=false;
  // cout<<"file: "<<file_path<<endl;
  // ifstream input_file(file_path);
  std::stringstream file_path;
  string f="OutputFiles/output_";
  file_path << f.data() << ifile << ".txt";  
  ifstream input_file(file_path.str().c_str());
  int state=0;
  Event current;
  string line;
  while(getline(input_file,line))
  {
   // cout<<line<<endl;
   trim(line);
   // if(line.length()==0) continue;
   // if(ievent==200000) break;
   if(state==0 and string::npos != line.find("BEGINNINGOFEVENT"))
   {
    state=1;
    continue;
   }
   if(state==1 and string::npos==line.find("PARTICLES"))
   {
    cout<<"Beginning of event"<<endl;
    int n_stable;
    double total_E;
    double frame_E;
    double b;
    int k;
    // getline(input_file,line);
    // trim(line);
    stringstream line_s;
    cout<<line<<endl;
    line_s << line;
    line_s >> k >> n_stable >> total_E;
    getline(input_file,line);
  	trim(line);
    stringstream line_ss;
    cout<<line<<endl;
    line_ss << line;
    line_ss >> frame_E >> b;
    cout << frame_E<<" "<<b<<endl;
    // state==1;
    current.set_n_stable(n_stable);
    current.set_total_E(total_E);
    current.set_frame_E(frame_E);
    current.set_b(b);
   }
   if(state==1 and string::npos != line.find("PARTICLES"))
   {
    state=2;
    cout<<"Beginning of particles"<<endl;
    continue;
   }
   if(state==2 and string::npos==line.find("END"))
   {
    int ind;
    int mom_ind;
    int pdg_id;
    int status;
    double px;
    double py;
    double pz;
    double E;    
    stringstream line_s;
    line_s << line;
    line_s >> pdg_id >> mom_ind >> status;
    getline(input_file,line);
  	trim(line);
    stringstream line_ss;
    line_ss << line;
    line_ss >> ind >> px >> py;
    getline(input_file,line);
  	trim(line);
    stringstream line_sss;
    line_sss << line;
    line_sss >> ind >> pz >> E;
    Particle *particle = new Particle(ind, mom_ind, pdg_id, status, px, py, pz, E);
    current.add_particle(particle);
   }
   // else if(state==2 and string::npos != line.find("END"))
   if(state==2 and string::npos != line.find("END"))
   { 
    cout<<"found the end"<<endl;
    state=0;
  	events.push_back(current);
    current.reset();
    ievent++;
    cout<<ievent<<endl;
   }
   // else
   // {
     // cout <<"Unexpected else. Exiting..." << endl;
     // exit(1);
   // }
  }
  input_file.close();
  cout<<"ifile="<<ifile<<endl;
 }
 }
};