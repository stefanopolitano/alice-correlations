#include <TAxis.h>
Float_t gpTMin = 1.;
Float_t gpTMax = 49.99;
Float_t gZVtxRange = 7;

void SetupRanges(CorrelationContainer *obj)
{
  obj->setEtaRange(-1.5, 1.5);
  obj->setPtRange(gpTMin, gpTMax);
  if (gZVtxRange > 0)
    obj->setZVtxRange(-gZVtxRange + 0.01, gZVtxRange - 0.01);
}
void extractEfficiency(const char *fileNameEfficiency = "AnalysisResults.root", const char *outputFile = "Efficiency.root", const char *folder = "correlation-task", Bool_t condenseMultiplicity = kFALSE, Bool_t extrapolatePt = kTRUE)
{
  gStyle->SetOptStat(1111111);
  CorrelationContainer::EfficiencyStep step1 = CorrelationContainer::MC;
  CorrelationContainer::EfficiencyStep step2 = CorrelationContainer::RecoAll;

  auto file = TFile::Open(outputFile, "RECREATE");
  file->Close();

  const Double_t assocPtArr[10] = {0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0};

  auto inputFile = TFile::Open(fileNameEfficiency);
  auto h = (CorrelationContainer *)inputFile->Get(Form("%s/sameEvent", folder));
  auto hMixed = (CorrelationContainer *)inputFile->Get(Form("%s/mixedEvent", folder));

  THn *efficiency4D, *efficiency4D_new;

  SetupRanges(h);

  efficiency4D = (THn *)h->getTrackEfficiencyND(step1, step2);
  cout <<efficiency4D->GetAxis(1)->GetBinWidth(1)<<endl;

  Int_t nBins[] = {efficiency4D->GetAxis(0)->GetNbins(), efficiency4D->GetAxis(1)->GetNbins(), 1, efficiency4D->GetAxis(3)->GetNbins()};
  Double_t xMin[] = {efficiency4D->GetAxis(0)->GetXmin(), efficiency4D->GetAxis(1)->GetXmin(), 0, efficiency4D->GetAxis(3)->GetXmin()};
  Double_t xMax[] = {efficiency4D->GetAxis(0)->GetXmax(), efficiency4D->GetAxis(1)->GetXmax(), 101, efficiency4D->GetAxis(3)->GetXmax()};
  efficiency4D_new = new THnF("efficiency4D_new", "", 4, nBins, xMin, xMax);
  efficiency4D_new->SetBinEdges(1,assocPtArr);

  if (condenseMultiplicity) // condensing centrality
  {
    for (int i = 0; i < 4; i++)
    efficiency4D_new->GetAxis(i)->SetTitle(efficiency4D->GetAxis(i)->GetTitle());
    efficiency4D_new->RebinnedAdd(efficiency4D);
    efficiency4D->Reset();
    efficiency4D = (THnF *)efficiency4D_new->Clone("efficiency4D");
  }

  double maxEffValue = 5.;
  for (int bin0 = 1; bin0 <= efficiency4D->GetAxis(0)->GetNbins(); bin0++)
    for (int bin1 = 1; bin1 <= efficiency4D->GetAxis(1)->GetNbins(); bin1++)
      for (int bin2 = 1; bin2 <= efficiency4D->GetAxis(2)->GetNbins(); bin2++)
        for (int bin3 = 1; bin3 <= efficiency4D->GetAxis(3)->GetNbins(); bin3++)
        {
          nBins[0] = bin0;
          nBins[1] = bin1;
          nBins[2] = bin2;
          nBins[3] = bin3;
          if (efficiency4D->GetBinContent(nBins) > maxEffValue)
          {
            LOGF(info, "Nulling %d %d %d %d %.2f %.2f %.2f %.2f which was %f", bin0, bin1, bin2, bin3, efficiency4D->GetAxis(0)->GetBinCenter(bin0), efficiency4D->GetAxis(1)->GetBinCenter(bin1), efficiency4D->GetAxis(2)->GetBinCenter(bin2), efficiency4D->GetAxis(3)->GetBinCenter(bin3), efficiency4D->GetBinContent(nBins));
            efficiency4D->SetBinContent(nBins, 0);
          }
        }

  const float fitRangeBegin = 5.01;
  const float fitRangeEnd = 14.99;
  const float extendRangeBegin = 8.01;
  Bool_t verbose = kTRUE;
  if (extrapolatePt) // extrapolating pT
  {
    Printf("Extrapolating high pT...");
    for (int bin0 = 1; bin0 <= efficiency4D->GetAxis(0)->GetNbins(); bin0++)
      for (int bin2 = 1; bin2 <= efficiency4D->GetAxis(2)->GetNbins(); bin2++)
        for (int bin3 = 1; bin3 <= efficiency4D->GetAxis(3)->GetNbins(); bin3++)
        {
          efficiency4D->GetAxis(0)->SetRange(bin0, bin0);
          efficiency4D->GetAxis(2)->SetRange(bin2, bin2);
          efficiency4D->GetAxis(3)->SetRange(bin3, bin3);

          if (gRandom->Uniform() < 0.01)
            verbose = kTRUE;
          TH1 *proj = efficiency4D->Projection(1);
          proj->SetName(Form("Proj_%i_%i_%i",bin0,bin2,bin3));
          if (proj->Integral() <= 0)
            continue;
          if (proj->Integral(proj->FindBin(fitRangeBegin), proj->FindBin(fitRangeEnd)) <= 0)
            continue;

          LOGF(info, "%d %d %d %f", bin0, bin2, bin3, proj->Integral(proj->FindBin(fitRangeBegin), proj->FindBin(fitRangeEnd)));
          if (verbose)
          {
            new TCanvas;
            proj->Draw();
          }
          proj->Fit("pol0", (verbose) ? "+" : "Q0+", "SAME", fitRangeBegin, fitRangeEnd);

          if (!proj->GetFunction("pol0"))
            continue;

          float trackingEff = proj->GetFunction("pol0")->GetParameter(0);

          for (int bin1 = 1; bin1 <= efficiency4D->GetAxis(1)->GetNbins(); bin1++)
          {
            if (efficiency4D->GetAxis(1)->GetBinCenter(bin1) < extendRangeBegin)
              continue;
            nBins[0] = bin0;
            nBins[1] = bin1;
            nBins[2] = bin2;
            nBins[3] = bin3;
            LOGF(info, "Setting %d %d %d %d %.2f %.2f %.2f %.2f to %f which was %f", bin0, bin1, bin2, bin3, efficiency4D->GetAxis(0)->GetBinCenter(bin0), efficiency4D->GetAxis(1)->GetBinCenter(bin1), efficiency4D->GetAxis(2)->GetBinCenter(bin2), efficiency4D->GetAxis(3)->GetBinCenter(bin3), trackingEff, efficiency4D->GetBinContent(nBins));
            efficiency4D->SetBinContent(nBins, trackingEff);
            if (verbose)
              verbose = kFALSE;
          }
          efficiency4D->GetAxis(0)->UnZoom();
          efficiency4D->GetAxis(1)->UnZoom();
          efficiency4D->GetAxis(2)->UnZoom();
          efficiency4D->GetAxis(3)->UnZoom();
        }
  }

  file = TFile::Open(outputFile, "UPDATE");
  efficiency4D->SetName("ccdb_object");
  efficiency4D->Write();
  file->Close();

  delete h;
  delete hMixed;
}
