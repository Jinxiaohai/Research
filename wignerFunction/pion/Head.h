#ifndef HEAD_H
#define HEAD_H

TH1D *fb = new TH1D("fb", "fb", 100, 0, 0.3);
TH1D *yield = new TH1D("yield", "yield", 25, 0, 5);
TH1D *EventNumHist = new TH1D("EventNumHist", "EventNumHist", 500, 0, 20000);


TH1D *pt_yield      = new TH1D("pt_yield",      "pt_yield",      25, 0, 5);
TH1D *v2SumHistgram = new TH1D("v2SumHistgram", "v2SumHistgram", 25, 0, 5);
TH1D *v3SumHistgram = new TH1D("v3SumHistgram", "v3SumHistgram", 25, 0, 5);

TH1D *probablySum = new TH1D("probablySum", "ProbablySum", 25, 0, 5);

#endif /* HEAD_H */
