void draw() {
    auto file = TFile::Open("flow.root");

    int colors[] = { 1, 2, 3, 4, 6 };

    const char* names[] = { "v22NoMassTemplate", "v22NoMassLMSubtraction"};
    for (std::size_t pass = 0; pass < std::size(names); ++pass) {

        auto v22NoMass = (TH3F*) file->Get(names[pass]);

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
            v2->SetLineColor(colors[icent-2]);
            v2->SetMarkerColor(colors[icent-2]);
            v2->SetLineStyle(pass+1);
            v2->Draw((icent == 2 && pass == 0) ? "AL*" : "LSAME*");
            v2->GetYaxis()->SetRangeUser(0, 0.25);
            printf("\nPass = %zu, Cent = %d:\n", pass, icent);
            v2->Print();
        }
    }
}
