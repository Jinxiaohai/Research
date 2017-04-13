#ifndef HEAD_H
#define HEAD_H

/// rapidity
TH1D *rap_d = new TH1D("rap_d", "rap_d", 100, -10, 10);
TH1D *rap_u = new TH1D("rap_u", "rap_u", 100, -10, 10);
TH1D *rap_s = new TH1D("rap_s", "rap_s", 100, -10, 10);
TH1D *rap_anti_d = new TH1D("rap_anti_d", "rap_anti_d", 100, -10, 10);
TH1D *rap_anti_u = new TH1D("rap_anti_u", "rap_anti_u", 100, -10, 10);
TH1D *rap_anti_s = new TH1D("rap_anti_s", "rap_anti_s", 100, -10, 10);

/// PT
TH1D *pt_d = new TH1D("pt_d", "pt_d", 100, 0, 10);
TH1D *pt_u = new TH1D("pt_u", "pt_u", 100, 0, 10);
TH1D *pt_s = new TH1D("pt_s", "pt_s", 100, 0, 10);
TH1D *pt_anti_d = new TH1D("pt_anti_d", "pt_anti_d", 100, 0, 10);
TH1D *pt_anti_u = new TH1D("pt_anti_u", "pt_anti_u", 100, 0, 10);
TH1D *pt_anti_s = new TH1D("pt_anti_s", "pt_anti_s", 100, 0, 10);

#endif /* HEAD_H */
