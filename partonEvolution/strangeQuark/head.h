#ifndef HEAD_H
#define HEAD_H

TH1D *pt = new TH1D("pt", "pt", 25, 0, 5);
TProfile *Strangev2 = new TProfile("Strangev2", "Strangev2", 25, 0, 5, 0, 1);

#endif /* HEAD_H */
