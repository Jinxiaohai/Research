#ifndef HEAD_H
#define HEAD_H
#include <vector>
using namespace std;

void processEvent();
float calc4_event(float Xn, float Yn, float X2n, float Y2n, float M);
float calc4_track(float xn, float yn, float x2n, float y2n, float Xn,
                  float Yn, float X2n, float Y2n, float M);
void computeFlow();

const unsigned int NOB = 100;
const unsigned int NBIN = 50;
namespace xiaohai
{
  struct Particle
  {
    int id;
    float px;
    float py;
    float pz;
    float x;
    float y;
    float z;
    float phi;
    float rsquare;
    float pT;
    float eta;
  };
}
vector<xiaohai::Particle> finalparticles;

//事件的特征变量
int npart          = 0;
int npartsum       = 0;
int nspectator     = 0;
int numevent       = 0;
int process_number = 0;

//Event-wise reference flow
//Must be cleared after every event
/// reference流
float q2x  = 0;
float q2y  = 0;
float q4x  = 0;
float q4y  = 0;
float sb_2 = 0;
float sb_4 = 0;

int evtnumber = 0;

float etaGap = 0.75;

TProfile *h_db_2      = new TProfile("h_db_2", "h_db_2", 1, 0, 1, -500, 500);
TProfile *h_db_4      = new TProfile("h_db_4", "h_db_4", 1, 0, 1, -500, 500);
TProfile *h_db_2prime = new TProfile("h_db_2prime", "h_db_2prime", NBIN, 0, 5, -500, 500);
TProfile *h_db_4prime = new TProfile("h_db_4prime", "h_db_4prime", NBIN, 0, 5, -500, 500);

TProfile *h_db_2_etagap = new TProfile("h_db_2_etagap", "h_db_2_etagap", 1, 0, 1, -500, 500);
TProfile *h_db_2prime_etagap = new TProfile("h_db_2prime_etagap", "h_db_2prime_etagap", NBIN, 0, 5, -500, 500);

TH1F *h_sb_2prime = new TH1F("h_sb_2prime", "h_sb_2prime", NBIN, 0, 5);
TH1F *h_sb_4prime = new TH1F("h_sb_4prime", "h_sb_4prime", NBIN, 0, 5);
TH1F *hn          = new TH1F("hn", "hn", NBIN, 0, 5);
TH1F *hp2x        = new TH1F("hp2x", "hp2x", NBIN, 0, 5);
TH1F *hp2y        = new TH1F("hp2y", "hp2y", NBIN, 0, 5);
TH1F *hp4x        = new TH1F("hp4x", "hp4x", NBIN, 0, 5);
TH1F *hp4y        = new TH1F("hp4y", "hp4y", NBIN, 0, 5);

TH1F *hEta         = new TH1F("hEta", "hEta;#eta;Counts", NOB, -6, 6);
TH1F *hpT          = new TH1F("hpT", "hpT;p_{T};Counts", NOB, 0, 10);
TH1F *hpT_EtaCut   = new TH1F("hpT_EtaCut", "hpT_EtaCut;p_{T};Counts", NOB, 0, 20);
TH1F *hPhi         = new TH1F("hPhi", "hPhi;#phi;Counts", NOB, -TMath::Pi(), TMath::Pi());


#endif /* HEAD_H */
