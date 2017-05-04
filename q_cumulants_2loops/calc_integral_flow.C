//--------------------------------
// Calc Integral Flow
//--------------------------------

void calc_integral_flow(TString inputFileName="output.root", double ptLow=0, double ptUp=2.0, double etaLow=-1, double etaUp=1)
{
  TFile* inputFile = new TFile(inputFileName.Data(), "OPEN");
  if (!inputFile->IsOpen()) cout<<"File Open Error"<<endl;

  enum POI
  {
    POI_CHARGED, POI_PIPLUS, POI_PIMINUS, POI_KPLUS, POI_KMINUS,
    POI_NUM=5
  };
  //   enum POI {
  //     POI_CHARGED, POI_PIPLUS, POI_PIMINUS, POI_KPLUS, POI_KMINUS, POI_PPLUS, POI_PMINUS,
  //     POI_NUM=7
  //   };
  POI poi = POI_PIPLUS;
  //   double ptBin[21] = {0.0, 0.2, 0.4, 0.6, 0.8,
  // 		      1.0, 1.2, 1.4, 1.6, 1.8,
  // 		      2.0, 2.2, 2.4, 2.6, 2.8,
  // 		      3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
  //   double ptBin[11] = {0.0, 0.2, 0.4, 0.6, 0.8,
  // 		      1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
  double ptBin[21] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
		      1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};

  //   TH3D* h3_vn2 = {NULL};
  //   TH3D* h3_yield[9] = {NULL};
  //   TH1D* h1_vn2POIPt[9] = {NULL};
  //   TH1D* h1_yldPOIPt[9] = {NULL};
  //   TH1D* h1_vn2POIPtRebin[9] = {NULL};
  //   TH1D* h1_yldPOIPtRebin[9] = {NULL};


  
  // /// 前后向快度
  // TH3D* h3_vn2 = (TH3D*)inputFile->Get("v22Hist");
  // TH3D* h3_yield = (TH3D*)inputFile->Get("yldHist");

  /// 后向快度
  TH3D* h3_vn2 = (TH3D*)inputFile->Get("vn2Hist_sub0");
  TH3D* h3_yield = (TH3D*)inputFile->Get("yldHist_sub0");

  // /// 前向快度
  // TH3D* h3_vn2 = (TH3D*)inputFile->Get("vn2Hist_sub1");
  // TH3D* h3_yield = (TH3D*)inputFile->Get("yldHist_sub1");


  
  int ptBinLow = h3_vn2->GetYaxis()->FindBin(ptLow + 1.0e-6);
  int ptBinUp = h3_vn2->GetYaxis()->FindBin(ptUp - 1.0e-6);
  int etaBinLow = h3_vn2->GetZaxis()->FindBin(etaLow + 1.0e-6);
  int etaBinUp = h3_vn2->GetZaxis()->FindBin(etaUp - 1.0e-6);
  TH1D* h1_vn2POIPt = new TH1D("vn2_pt", "", 40, 0, 2);
  TH1D* h1_yldPOIPt = new TH1D("yield_pt", "", 40, 0, 2);

  //--------------------------------
  // Integrate eta
  //--------------------------------

  for (int iPtBin=ptBinLow; iPtBin<=ptBinUp; ++iPtBin)
    {
      double vSum      = 0.;
      double yieldSum  = 0.;
      double error2sum = 0.;
      double vSqrSum   = 0.;
      for (int iEtaBin=etaBinLow; iEtaBin<=etaBinUp; ++iEtaBin)
        {
          double v     = h3_vn2->GetBinContent(poi+1, iPtBin, iEtaBin);
          double verr  = h3_vn2->GetBinError(poi+1, iPtBin, iEtaBin);
          double yield = h3_yield->GetBinContent(poi+1, iPtBin, iEtaBin);
          if (TMath::IsNaN(v)) continue;
          if (TMath::IsNaN(verr)) continue;
          if (TMath::IsNaN(yield)) continue;
          if (yield<5) continue;
          if (v != 0)
            {
              vSum      += yield * v;
              yieldSum  += yield;
              error2sum += pow(yield*verr, 2.);
              vSqrSum   += yield * (yield*pow(verr, 2.) + pow(v, 2.));
            }
        }

      if (yieldSum)
        {
          h1_vn2POIPt->SetBinContent(iPtBin, vSum / yieldSum);
          h1_vn2POIPt->SetBinError(iPtBin, sqrt(vSqrSum/yieldSum - pow(vSum/yieldSum, 2.))
                                   / sqrt(yieldSum));
          h1_yldPOIPt->SetBinContent(iPtBin, yieldSum);
        }
    }

  //--------------------------------
  // Integrate pT
  //--------------------------------
    
  h1_vn2POIPtRebin = new TH1D("vn2POIPtRebin", "", 20, ptBin);
  h1_yldPOIPtRebin = new TH1D("yldPOIPtRebin", "", 20, ptBin);
  for (int i=0; i<20; ++i)
    {
      int singlePtBinLow = h1_vn2POIPt->GetXaxis()->FindBin(ptBin[i] + 1.0e-6);
      int singlePtBinUp = h1_vn2POIPt->GetXaxis()->FindBin(ptBin[i+1] - 1.0e-6);
      double vSum      = 0.;
      double yieldSum  = 0.;
      double error2sum = 0.;
      double vSqrSum   = 0.;
      for (int k = singlePtBinLow; k <= singlePtBinUp; ++k)
        {
          double v     = h1_vn2POIPt->GetBinContent(k);
          double verr  = h1_vn2POIPt->GetBinError(k);
          double yield = h1_yldPOIPt->GetBinContent(k);
          if (TMath::IsNaN(v)) continue;
          if (TMath::IsNaN(verr)) continue;
          if (TMath::IsNaN(yield)) continue;
          if (yield<5) continue;
          if (v != 0)
            {
              vSum      += yield * v;
              yieldSum  += yield;
              error2sum += pow(yield*verr, 2.);
              vSqrSum   += yield * (yield*pow(verr, 2.) + pow(v, 2.));
            }
        }
      //            cout<<iCent<<" "<<vSum<<" "<<yieldSum<<endl;
      if (yieldSum)
        {
          h1_vn2POIPtRebin->SetBinContent(i+1, vSum / yieldSum);
          h1_vn2POIPtRebin->SetBinError(i+1, sqrt(vSqrSum/yieldSum - pow(vSum/yieldSum, 2.))
                                        / sqrt(yieldSum));
          h1_yldPOIPtRebin->SetBinContent(i+1, yieldSum);
        }
    }
  //--------------------------------
  // Draw
  //--------------------------------
    
  TCanvas* allCentIn1 = new TCanvas("allCentIn1", "", 700, 600);
  //   allCentIn1->Divide(3, 3);
  //   TH2D* dummy2 = new TH2D("dummy2", ";p_{T};v_{n}", 20, 0, ptBin[20], 10, -0.05, 0.2);
  TH2D* dummy2 = new TH2D("dummy2", ";p_{T};v_{n}", 10, 0, ptBin[20], 10, -0.05, 0.2);
  dummy2->Draw();
  gPad->SetLeftMargin(0.2);

  h1_vn2POIPtRebin->SetMarkerStyle(kOpenCircle);
  h1_vn2POIPtRebin->SetMarkerSize(1.3);
  h1_vn2POIPtRebin->SetMarkerColor(kBlack);
  h1_vn2POIPtRebin->Draw("same p");
    
  allCentIn1->Update();
  allCentIn1->Print("200_melting_b8_v2.pdf");

  TFile qa("200_melting_b8_v2.root","RECREATE");
  h1_vn2POIPtRebin->Write();
  qa.Close();


  //   TH1D* h1_copy[9] = {NULL};
  //   for (int i=0; i<9; ++i) {

  //     allCentIn1->cd(i+1);
  //     dummy2->Draw();
  //     gPad->SetLeftMargin(0.2);

  //     h1_vn2POIPtRebin->SetMarkerStyle(kOpenCircle);
  //     h1_vn2POIPtRebin->SetMarkerSize(1.3);
  //     h1_vn2POIPtRebin->SetMarkerColor(kBlack);
  //     h1_vn2POIPtRebin->Draw("same p");
  //     h1_copy[i] = (TH1D*)h1_vn2POIPtRebin->Clone("h1_copy");
    
  //     //     h1_vn2POIPt[i]->SetMarkerStyle(kOpenCircle);
  //     //     h1_vn2POIPt[i]->SetMarkerSize(1.1);
  //     //     h1_vn2POIPt[i]->SetMarkerColor(kBlack);
  //     //     h1_vn2POIPt[i]->Draw("same p");
  //     //     h1_copy[i] = (TH1D*)h1_vn2POIPt[i]->Clone("h1_copy");
      
  //     h1_copy[i]->SetMarkerStyle(kFullCircle);
  //     h1_copy[i]->SetMarkerSize(1.2);
  //     h1_copy[i]->SetMarkerColor(kRed);
  //     h1_copy[i]->Draw("same p");
      
  //     allCentIn1->Update();
  //     //     allCentIn1->Print("test.pdf");
  //   }   
}
