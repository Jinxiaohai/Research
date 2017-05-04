//--------------------------------
// Calc v_{n}{2}
//--------------------------------

void calc_flow(TString inputFileName="test.root")
{
  TFile* inFile = new TFile(inputFileName.Data(), "OPEN");
  if (!inFile->IsOpen())
    {
     cout<<"File Open Error"<<endl; 
    }
  TFile* outFile = new TFile("output.root", "RECREATE");

  /// bin number and bin.
  const int	NUM_POI	    = 5;
  const int     NUM_ETABIN  = 20;
  const double  FLOW_ETALOW = -1.0;
  const double  FLOW_ETAUP  = 1.0;
  const int     NUM_PTBIN   = 40;
  const double  FLOW_PTLOW  = 0;
  const double  FLOW_PTUP   = 2.0;

  //--------------------------------
  //
  //--------------------------------

  TProfile* pSubEventRefCorr = (TProfile*)inFile->Get("subEventRefCorr");
  TH1D* hSubEventSumW2 = (TH1D*)inFile->Get("subEventSumW2");
  TProfile* pRefCorrelations[2] = {NULL};
  TH1D* hRefBinSumW2[2] = {NULL};
  TProfile3D* pC14[2] = {NULL};
  TProfile3D* pC15[2] = {NULL};
  TProfile3D* pT28[2] = {NULL};
  TProfile3D* pT16abT28[2] = {NULL};
  TH3D* hDiffBinSumW12[2] = {NULL};
  TH3D* hDiffBinSumW22[2] = {NULL};
  TH3D* hYield[2] = {NULL};
  
  TH3D* subV22Hist[2] = {NULL};
  /// get 出存取的数据
  for (Int_t i_loop = 0; i_loop < 2; ++i_loop)
    {
      TString subPostfix("_sub");
      subPostfix += i_loop;
      pRefCorrelations[i_loop] = (TProfile*)inFile->
        Get(TString("refCorrelations").Append(subPostfix).Data());
      hRefBinSumW2[i_loop] = (TH1D*)inFile->Get(TString("refBinSumW2").Append(subPostfix).Data());
      pC14[i_loop] = (TProfile3D*)inFile->Get(TString("C14").Append(subPostfix).Data());
      pC15[i_loop] = (TProfile3D*)inFile->Get(TString("C15").Append(subPostfix).Data());
      pT28[i_loop] = (TProfile3D*)inFile->Get(TString("T28").Append(subPostfix).Data());
      hDiffBinSumW12[i_loop] = (TH3D*)inFile->Get(TString("diffBinSumW12").Append(subPostfix).Data());
      hDiffBinSumW22[i_loop] = (TH3D*)inFile->Get(TString("diffBinSumW22").Append(subPostfix).Data());
      pT16abT28[i_loop] = (TProfile3D*)inFile->Get(TString("T16abT28").Append(subPostfix).Data());
      hYield[i_loop]  = (TH3D*)inFile->Get(TString("yield").Append(subPostfix).Data());
      TString yldName(hYield[i_loop]->GetName());
      yldName.ReplaceAll("yield", "yldHist");
      hYield[i_loop]->SetName(yldName.Data());
    
      subV22Hist[i_loop] = new TH3D(TString("vn2Hist").Append(subPostfix).Data(),
                                    "v22Hist;# POI;p_{T} (GeV/c);#eta",
                                    NUM_POI, 0, NUM_POI, NUM_PTBIN, FLOW_PTLOW, FLOW_PTUP,
                                    NUM_ETABIN, FLOW_ETALOW, FLOW_ETAUP);
      Double_t subEvtMean[5];
      Double_t subEvtMeanErr2[5];
      Double_t subEvtEffEntries[5];
      /// 取出数值，赋给上面的三个数组。
      for (Int_t i = 0; i < 5; ++i)
        {
          Int_t iBin          = i+1;
          subEvtMean[i]       = pSubEventRefCorr->GetBinContent(iBin);
          subEvtEffEntries[i] = TMath::Power(pSubEventRefCorr->GetBinEntries(iBin), 2)
            / hSubEventSumW2->GetBinContent(iBin);
          subEvtMeanErr2[i]   = TMath::Power(pSubEventRefCorr->GetBinError(iBin), 2)
            / (subEvtEffEntries[i]- 1);
        }
    
      Double_t c22ab = subEvtMean[0] - subEvtMean[1]*subEvtMean[2] - subEvtMean[3]*subEvtMean[4]; 

      //--------------------------------
      //
      //--------------------------------
    
      Double_t mean_ref[3] = {0}; // C2 C3 T16
      Double_t meanErr2_ref[3] = {0};
      Double_t effEntries_ref[3] = {0};
      for ( Int_t i = 0; i < 3; ++i )
        {
          Int_t iBin        = i + 1;
          mean_ref[i]       = pRefCorrelations[i_loop]->GetBinContent(iBin);
          effEntries_ref[i] = TMath::Power(pRefCorrelations[i_loop]->GetBinEntries(iBin), 2)
            / hRefBinSumW2[i_loop]->GetBinContent(iBin); // \sum{w_i}^2 / \sum{w_i^2}
          meanErr2_ref[i]   = TMath::Power(pRefCorrelations[i_loop]->GetBinError(iBin), 2)
            / (effEntries_ref[i] - 1);
        }

      Double_t mean_diff[3] = {0}; // C14 C15 T28
      Double_t meanErr2_diff[3] = {0};
      Double_t effEntries_diff[3] = {0};
      for ( Int_t iPOI = 0; iPOI < 5; ++iPOI )
        {
          for ( Int_t iPt = 0; iPt < 40; ++iPt )
            {
              for ( Int_t iEta = 0; iEta < 20; ++iEta )
                {
                  Int_t globalBinNum = pC14[i_loop]->GetBin(iPOI+1, iPt+1, iEta+1);
                  mean_diff[0]       = pC14[i_loop]->GetBinContent(globalBinNum);
                  effEntries_diff[0] = TMath::Power(pC14[i_loop]->GetBinEntries(globalBinNum), 2)
                    / hDiffBinSumW12[i_loop]->GetBinContent(globalBinNum); // \sum{w_i}^2 / \sum{w_i^2}
                  meanErr2_diff[0]   = TMath::Power(pC14[i_loop]->GetBinError(globalBinNum), 2)
                    / (effEntries_diff[0] - 1);
                  mean_diff[1]       = pC15[i_loop]->GetBinContent(globalBinNum);
                  effEntries_diff[1] = TMath::Power(pC15[i_loop]->GetBinEntries(globalBinNum), 2)
                    / hDiffBinSumW12[i_loop]->GetBinContent(globalBinNum); // \sum{w_i}^2 / \sum{w_i^2}
                  meanErr2_diff[1]   = TMath::Power(pC15[i_loop]->GetBinError(globalBinNum), 2)
                    / (effEntries_diff[1] - 1);
	  
                  mean_diff[2]       = pT28[i_loop]->GetBinContent(globalBinNum);
                  effEntries_diff[2] = TMath::Power(pT28[i_loop]->GetBinEntries(globalBinNum), 2)
                    / hDiffBinSumW22[i_loop]->GetBinContent(globalBinNum); // \sum{w_i}^2 / \sum{w_i^2}
                  meanErr2_diff[2]   = TMath::Power(pT28[i_loop]->GetBinError(globalBinNum), 2)
                    / (effEntries_diff[2] - 1);
                  Double_t d22 = mean_diff[2] - mean_diff[0]*mean_ref[0] - mean_diff[1]*mean_ref[1];
                  Double_t v22 = d22 / TMath::Sqrt(c22ab);
                  subV22Hist[i_loop]->SetBinContent(iPOI+1, iPt+1, iEta+1, v22);

                  // stat. error
                  Double_t v22Err = 0;

                  v22Err += subEvtMeanErr2[0] * TMath::Power(-(mean_diff[0]*mean_ref[0]) - mean_diff[1]*mean_ref[1] + mean_diff[2],2)/(4.*TMath::Power(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0],3));
                  v22Err += subEvtMeanErr2[1] * (TMath::Power(subEvtMean[2],2)*TMath::Power(-(mean_diff[0]*mean_ref[0]) - mean_diff[1]*mean_ref[1] + mean_diff[2],2))/(4.*TMath::Power(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0],3));
                  v22Err += subEvtMeanErr2[2] * (TMath::Power(subEvtMean[1],2)*TMath::Power(-(mean_diff[0]*mean_ref[0]) - mean_diff[1]*mean_ref[1] + mean_diff[2],2))/(4.*TMath::Power(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0],3));
                  v22Err += subEvtMeanErr2[3] * (TMath::Power(subEvtMean[4],2)*TMath::Power(-(mean_diff[0]*mean_ref[0]) - mean_diff[1]*mean_ref[1] + mean_diff[2],2))/(4.*TMath::Power(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0],3));
                  v22Err += subEvtMeanErr2[4] * (TMath::Power(subEvtMean[3],2)*TMath::Power(-(mean_diff[0]*mean_ref[0]) - mean_diff[1]*mean_ref[1] + mean_diff[2],2))/(4.*TMath::Power(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0],3));
	  
                  v22Err += meanErr2_ref[0]  * TMath::Power(mean_diff[0],2)/(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0]);
                  v22Err += meanErr2_ref[1]  * TMath::Power(mean_diff[1],2)/(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0]);
                  v22Err += meanErr2_diff[0] * TMath::Power(mean_ref[0],2)/(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0]);
                  v22Err += meanErr2_diff[1] * TMath::Power(mean_ref[1],2)/(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0]);
                  v22Err += meanErr2_diff[2] * 1/(-(subEvtMean[1]*subEvtMean[2]) - subEvtMean[3]*subEvtMean[4] + subEvtMean[0]);
 
                  // Cov(<2^{ab}>, <2'>)
                  Double_t covEffEntries = pT16abT28[i_loop]->GetBinEntries(globalBinNum) /
                    (pSubEventRefCorr->GetBinEntries(1) * pT28[i_loop]->GetBinEntries(globalBinNum));
                  Double_t covT16abT28 = (pT16abT28[i_loop]->GetBinContent(globalBinNum) -
                                          pSubEventRefCorr->GetBinContent(1) * pT28[i_loop]->GetBinContent(globalBinNum)) / (1 - covEffEntries);
	  
                  v22Err += 2 * ((mean_diff[0]*mean_ref[0] + mean_diff[1]*mean_ref[1] - mean_diff[2])/(2.*TMath::Power(subEvtMean[1]*subEvtMean[2] + subEvtMean[3]*subEvtMean[4] - subEvtMean[0],2)))
                    * covEffEntries * covT16abT28;
	  
                  v22Err = TMath::IsNaN(v22Err) ? 0 : TMath::Sqrt(v22Err); 
                  subV22Hist[i_loop]->SetBinError(iPOI+1, iPt+1, iEta+1, v22Err);                
                }/// eta
            }/// pt
        }/// poi
    }/// i_loop

  //--------------------------------
  //
  //--------------------------------

  if (!subV22Hist[0] || !subV22Hist[1])
    {
      cerr << "Error: sub-event flow not available" << endl;
      return;
    }

  TH3D* v22Hist = new TH3D("v22Hist", "v22Hist;# POI;p_{T} (GeV/c);#eta",
			   5, 0, 5, 40, 0, 2, 20, -1, 1);
  TH3D* yldHist = new TH3D("yldHist", "yldHist;# POI;p_{T} (GeV/c);#eta",
			   5, 0, 5, 40, 0, 2, 20, -1, 1);

  for (Int_t i_loop = 0; i_loop < 2; ++i_loop)
    {
      v22Hist->Add(subV22Hist[i_loop]);
      yldHist->Add(hYield[i_loop]);
    }

  /// 前后向的快度
  outFile->cd();
  if (v22Hist)
    {
      v22Hist->Write();
      yldHist->Write();
      cout << "Hist: " << v22Hist->GetName() << " is saved!" << endl;
    }

  /// 后向快度  
  if (subV22Hist[0])
    {
      subV22Hist[0]->Write();
      hYield[0]->Write();
      cout << "Hist: " << subV22Hist[0]->GetName() << " is saved!" << endl;
    }
  /// 前向快度
  if (subV22Hist[1])
    {
      subV22Hist[1]->Write();
      hYield[1]->Write();
      cout << "Hist: " << subV22Hist[1]->GetName() << " is saved!" << endl;
    }

  inFile->Close();
  outFile->Close();
}
