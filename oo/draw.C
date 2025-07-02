void draw() {
    auto file = TFile::Open("flow.root");

    auto v22NoMass = (TH3F*) file->Get("v22NoMass");

    for (int icent = 2; icent <= v22NoMass->GetNbinsZ(); icent++) {
        auto v2 = new TGraphErrors;
        for (int ipt = 1; ipt <= v22NoMass->GetNbinsX(); ipt++) {
            float v22 = v22NoMass->GetBinContent(ipt, ipt, icent);
            float v22err = v22NoMass->GetBinError(ipt, ipt, icent);
            if (v22 > 0) {
                v2->SetPoint(v2->GetN(), v22NoMass->GetXaxis()->GetBinCenter(ipt), std::sqrt(v22));
                v2->SetPointError(v2->GetN()-1, 0, 0.5 / std::sqrt(v22) * v22err);
            }
        }
        v2->SetLineColor(icent);
        v2->Draw((icent == 2) ? "AL*" : "LSAME*");
    }
}
