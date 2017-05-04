/**
 * @file   ampt.cpp
 * @author xiaohai <xiaohaijin@outlook.com>
 * @date   Fri Apr 28 08:16:41 2017
 * 
 * @brief  The file was created to draw
 *         a streamplot vs heatmap.
 */
#include <cmath>

int backward()
{
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 0, 0, 900, 900);
  TPad *pad = new TPad("pad", "pad", 0.05, 0.05, 0.95, 0.95);
  pad->Draw();
  pad->cd();
  ifstream input("../../toydata/zpc.dat");
  const int CellNum = 14;
  double range  = 10.;/// [-10, 10)
  TH2D *scatter = new TH2D("scatter", "scatter", 80, -10, 10, 80, -10, 10);
  unsigned int nevent = 0, nrun = 0, multi = 0,
    nelp = 0, ninp = 0, nelt = 0, nint = 0;
  double impactpar = 0.;
  int id = 0;
  double px = 0., py = 0., pz = 0., xmass = 0., x = 0, y = 0., z = 0., Time = 0.;
  double cell[CellNum][CellNum][4] = {0};
  /// cell[][][0]  vx
  /// cell[][][1]  vy
  /// cell[][][2]  theta
  /// cell[][][3]  count
  double delta = (2*range)/CellNum;
  input >> nevent >> nrun >> multi >> impactpar
        >> nelp >> ninp >> nelt >> nint;
  for (int i = 0; i != multi; ++i)
    {
      input >> id >> px >> py >> pz >> xmass >> x >> y >> z >> Time;
      if (!input) break;
      scatter->Fill(x, y);
      if (pz < 0)
        {
          for (int j = 0; j != CellNum; ++j)
            {
              if (x > j*delta-range && x < (j+1)*delta-range)
                {
                  for (int k = 0; k != CellNum; ++k)
                    {
                      if (y > k*delta-range && y < (k+1)*delta-range)
                        {
                          cell[j][k][0] += px;
                          cell[j][k][1] += py;
                          cell[j][k][3] += 1;
                        }
                    }
                }
            }
        }
    }
  scatter->GetXaxis()->SetTitle("x [fm]");
  scatter->GetXaxis()->SetTitleSize(0.035);
  scatter->GetXaxis()->CenterTitle(true);
  scatter->GetYaxis()->SetTitle("y [fm]");
  scatter->GetYaxis()->SetTitleSize(0.035);
  scatter->GetYaxis()->CenterTitle(true);

  scatter->Draw("cont0");
  scatter->SetTitle("");

  /// TPaveLabel
  TPaveLabel *par = new TPaveLabel(-8, 15.5, 8, 17.8, "streamplot vs heatmap");
  par->SetFillColor(42);
  par->Draw();

  /// streamplot
  for (int i = 0; i != CellNum; ++i)
    {
      for (int j = 0; j != CellNum; ++j)
        {
          cell[i][j][2] = atan2(cell[i][j][1], cell[i][j][0]);
        }
    }
  /// fixed length
  // double length = sqrt(delta*delta + delta*delta);
  for (int i = 0; i != CellNum; ++i)
    {
      for (int j = 0; j != CellNum; ++j)
        {
          if (cell[i][j][3] > 8)
            {
              double starx = i*delta + delta/2 - range;
              double stary = j*delta + delta/2 - range;
              // double length = sqrt(cell[i][j][0]*cell[i][j][0]
              //                      + cell[i][j][1]*cell[i][j][1])/cell[i][j][3];
              double length = sqrt(cell[i][j][0]*cell[i][j][0]
                                   + cell[i][j][1]*cell[i][j][1]);
              /// the parameter is to adjust the length of the arrow.
              double endx = starx + length*cos(cell[i][j][2])/35;
              double endy = stary + length*sin(cell[i][j][2])/35;
              TArrow *arr = new TArrow(starx, stary, endx, endy, 0.008, "|>");
              arr->SetAngle(45);
              arr->SetLineWidth(3);
              arr->SetLineColor(1);
              arr->Draw();
              pad->Update();
            }
        }
    }
  return 0;
}
