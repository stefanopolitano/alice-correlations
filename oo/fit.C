static TH1D *gpsimMass;
static TGraphErrors *gpvsb;

static const double massLimits[] = {1.7,2.0};

template<bool massFit, bool vnFit>
double SimFit(const double *par){
	//invariant mass fit
	double chi2_1 = 0.0, chi2_2 = 0.0;
	if(massFit){
		int a = gpsimMass->GetXaxis()->FindBin(massLimits[0]);
		int b = gpsimMass->GetXaxis()->FindBin(massLimits[1]);
		for(int i = a; i < b; ++i){
			double x = gpsimMass->GetBinCenter(i);
			double delta = gpsimMass->GetBinContent(i)
				-(par[0]+par[1]*x+par[2]*x*x+par[3]*TMath::Gaus(x,par[4],par[5]));
			double sigErr = gpsimMass->GetBinError(i);
			chi2_1 += delta*delta/(sigErr*sigErr);
		}
	}

	//vsb(m) = vs * (Ns/(Nb+Ns))(m) + vb(m) * (Nb/(Nb+Ns))(m), vb(m) = a*m+b
	if(vnFit){
		for(uint i = 0; i < gpvsb->GetN(); ++i){
			double x = gpvsb->GetX()[i];
			double Nbkg = par[0]+par[1]*x+par[2]*x*x;
			double Nsig = par[3]*TMath::Gaus(x,par[4],par[5]);
			double delta = gpvsb->GetY()[i]
				-(par[6]*Nsig/(Nsig+Nbkg)+(par[7]*x+par[8])*Nbkg/(Nsig+Nbkg)); //vn^B modeled with linear dependency, according to CMS
			double vsbErr = gpvsb->GetEY()[i];
			chi2_2 += delta*delta/(vsbErr*vsbErr);
		}
	}
	//chi2_2 /= double(gpvsb->GetN());
	return chi2_1+chi2_2;
}

TF1* extract(TH1D *phf1, const char *pnamepf, bool write = true){
	const TString fitName = Form("fit%s",pnamepf);
	
	auto fit = new TF1(fitName,"[0]+[1]*(1 + 2*[2]*TMath::Cos(x) + 2*[3]*TMath::Cos(2*x) + 2*[4]*TMath::Cos(3*x) + 2*[5]*TMath::Cos(4*x) + 2*[6]*TMath::Cos(5*x))");
	fit->SetParNames("czyam","c","v1","v2","v3","v4","v5");

	// fit->FixParameter(5,0.0); // v4 := 0
	// fit->FixParameter(6,0.0); // v5 := 0
	fit->FixParameter(0,0.0); // offset := 0

	// new TCanvas;
	const char* option = "Q0SE"; //"0QSE";
	TFitResultPtr r = phf1->Fit(fit,option,"",-TMath::Pi()/2.0,3.0/2.0*TMath::Pi());

	// fit->Print("V");

	if(write){
		fit->SetRange(-TMath::Pi()/2.0,3.0/2.0*TMath::Pi());
		fit->Write();
	}

	return fit;
}

TF1* extractSubtraction(TH1D *phf1, const char *pnamepf, TF1* subtract, bool write = true){
	const TString fitName = Form("fit%s",pnamepf);

	auto fit = new TF1(fitName,"[0]+[1]*(1 + 2*[2]*TMath::Cos(x)  + 2*[3] *TMath::Cos(2*x) + 2*[4] *TMath::Cos(3*x) + 2*[5] *TMath::Cos(4*x) + 2*[6]*TMath::Cos(5*x)) + "
		                  "[7]*([8]+[9]*(1 + 2*[10]*TMath::Cos(x) + 2*[11]*TMath::Cos(2*x) + 2*[12]*TMath::Cos(3*x) + 2*[13]*TMath::Cos(4*x) + 2*[14]*TMath::Cos(5*x)))"	
	);
	const char* parNames[] = {"czyam","c","v1","v2","v3","v4","v5","scaling",
					          "czyam_lm","c_lm","v1_lm","v2_lm","v3_lm","v4_lm","v5_lm" };
	for (int i=0; i<15; i++)
		fit->SetParName(i, parNames[i]);

	fit->FixParameter(2,0.0); // v1 := 0
	fit->FixParameter(5,0.0); // v4 := 0
	fit->FixParameter(6,0.0); // v5 := 0
	fit->FixParameter(0,0.0); // offset := 0

	fit->SetParameter(7, 1.6);
	fit->SetParLimits(7, 1.0, 3);

	for (int i=0; i<7; i++)
		fit->FixParameter(i+8, subtract->GetParameter(i));

	// new TCanvas;
	TString option = "0SE"; //"0SE";
	printf("Pass 1 (ignore the warnings):\n");
	phf1->Fit(fit,option + "Q","",-TMath::Pi()/2.0,3.0/2.0*TMath::Pi());
	printf("Pass 2 (check the output):\n");
	TFitResultPtr r = phf1->Fit(fit,option,"",-TMath::Pi()/2.0,3.0/2.0*TMath::Pi());

	// fit->Print("V");

	if(write){
		fit->SetRange(-TMath::Pi()/2.0,3.0/2.0*TMath::Pi());
		fit->Write();
	}

	return fit;
}

bool kSubtraction = kTRUE;

void fit(const char *pinFileName = "dphi_corr.root", const char *poutFileName = "flow.root", double absDeltaEtaMin = 1.0, double absDeltaEtaMax = 1.6, double absDeltaPhiMax = 1.3){
	TFile *pf = new TFile(pinFileName,"read");

	TTree *paxes = (TTree*)pf->Get("axes");
	UInt_t Ncent, Nmass, Ntrig, Nassoc;
	paxes->SetBranchAddress("Ncentrality",&Ncent);
	paxes->SetBranchAddress("Nmass",&Nmass);
	paxes->SetBranchAddress("NleadingPt",&Ntrig);
	paxes->SetBranchAddress("NassocPt",&Nassoc);
	paxes->GetEntry(0);

	float *pcentBins = new float[Ncent+1];
	float *pmassBins = new float[Nmass+1];
	float *pmeanCents = new float[Ncent+1];
	float *passocPt = new float[Nassoc+1];
	float *ptrigPt = new float[Ntrig+1];
	paxes->SetBranchAddress("centrality",pcentBins);
	paxes->SetBranchAddress("mass",pmassBins);
	paxes->SetBranchAddress("meanCent",pmeanCents);
	paxes->SetBranchAddress("assocPt",passocPt);
	paxes->SetBranchAddress("leadingPt",ptrigPt);
	paxes->GetEntry(0);

	printf("Assoc bins: ----- %u\n",Nassoc);
	printf("Trig bins: ----- %u\n",Ntrig);
	printf("Cent bins: ----- %u\n",Ncent);
	printf("Mass bins: ----- %u\n",Nmass);
	for(uint i = 0; i < Ncent; ++i)
		printf("bin %u = %lf\n",i,pmeanCents[i]);

	auto absDeltaEtaMinOrig = absDeltaEtaMin;
	auto absDeltaEtaMaxOrig = absDeltaEtaMax;

	absDeltaEtaMin += 0.001;
	absDeltaEtaMax -= 0.001;

	ROOT::Math::Minimizer *pmin = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
	pmin->SetMaxFunctionCalls(10000000); //for Minuit/Minuit2
	pmin->SetMaxIterations(10000); //for GSL
	pmin->SetTolerance(0.001);
	pmin->SetPrintLevel(1);

	ROOT::Math::Functor f(&SimFit<true,false>,6); //mass fit only
	ROOT::Math::Functor fStage2(&SimFit<false,true>,9); //vn fit only with fixed mass params
	ROOT::Math::Functor fStage3(&SimFit<true,true>,9); //refine both together dependently
	//ROOT::Math::Functor f(&SimFit,6);

	TH1 *pmult = pf->Get<TH1>("multiplicity");

	auto v22NoMassTemplate = new TH3F("v22NoMassTemplate", ";pT,trig;pT,assoc;multiplicity,v22", int(Ntrig), ptrigPt, int(Nassoc), passocPt, int(Ncent), pcentBins);
	auto v22NoMassLMSubtraction = (TH3*) v22NoMassTemplate->Clone("v22NoMassLMSubtraction");

	TFile *pfout = new TFile(poutFileName,"recreate");
	pfout->cd();

	for(uint itrig = 0; itrig < Ntrig; ++itrig){
		for(uint iassoc = itrig; iassoc <= itrig; ++iassoc){
		// for(uint iassoc = 0; iassoc <= 0; ++iassoc){
			TGraphErrors grv22(Ncent); //v22, no non-flow subtraction
			TGraphErrors grv22nfsub(Ncent-1); //v22
			TGraph grnfNAssoc(Ncent);
			TGraph grnfYJet(Ncent);
			double TempSub, TempSubErr;

			TF1* lowestFit = nullptr;
			TH1* lowestHist = nullptr;

			for(uint ib = 1; ib < Ncent; ++ib){ // skip lowest
			// for(int ib = Ncent.q-1; ib >= 0; --ib){ //for centrality, start from highest bin in order to handle the LM template first
				gpsimMass = (TH1D*)pf->Get(Form("mass_%u_%u",itrig,ib));

				TGraphErrors *pvsbMass = new TGraphErrors(Nmass); //mass dependent V22 Fourier coefficient
				TH1D *pjetYieldRaw;
				for(int imass = -1; imass < int(Nmass); ++imass){
					const TString masspf = imass < 0?"":Form("_mass%d",imass);
					const TString namepf = Form("_%u_%u_%u%s",itrig,iassoc,ib,masspf.Data());
					//
					TH2D *ph = (TH2D*)pf->Get(Form("dphi%s",namepf.Data())); //2d histogram to be normalized to DeltaEta
					if(!ph){
						printf("No histograms corresponding mult bin %u. (itrig=%u, iassoc=%u)\n",ib,itrig,iassoc);
						continue;
					}

					int a = ph->GetYaxis()->FindBin(absDeltaEtaMin);
					int b = ph->GetYaxis()->FindBin(absDeltaEtaMax);
					TH1D *pp = ph->ProjectionX(Form("proj_dphi_P%s",namepf.Data()),a,b,"e"); //positive side DeltaPhi long-range

					a = ph->GetYaxis()->FindBin(-absDeltaEtaMax);
					b = ph->GetYaxis()->FindBin(-absDeltaEtaMin);
					TH1D *pn = ph->ProjectionX(Form("proj_dphi_N%s",namepf.Data()),a,b,"e"); //negative side DeltaPhi long-range

					TH1D *phf1 = (TH1D*)pp->Clone(Form("proj_dphi%s",namepf.Data())); //combined DeltaPhi distribution (no CZYAM)
					phf1->Add(pp,pn,0.5,0.5);
					phf1->Rebin(2);
					phf1->Scale(0.5);

					double normPhi = 2.0*(absDeltaEtaMaxOrig-absDeltaEtaMinOrig);
					phf1->Scale(1.0/normPhi); //normalize to long-range width (small delta-factor)

					phf1->Write();

					a = ph->GetXaxis()->FindBin(-absDeltaPhiMax);
					b = ph->GetXaxis()->FindBin(absDeltaPhiMax);
					TH1D *pd = ph->ProjectionY(Form("proj_deta%s",namepf.Data())); //near side DeltaEta projection

					double normEta = 2.0*(absDeltaPhiMax);
					pd->Scale(1.0/normEta);

					if(imass < 0)
						pjetYieldRaw = (TH1D*)pd->Clone(); //yield needed later for the non-flow subtraction
					pd->Write();

					float v2template, v2templateerr;
					float scaling, scalingerr;
					float v2lmsub, v2lmsuberr;
					
					if (!lowestFit) {
						auto fit = extract(phf1, namepf.Data());
						v2template = v2lmsub = fit->GetParameter(3);
						v2templateerr = v2lmsuberr = fit->GetParError(3);
						scaling = scalingerr = 0;

						lowestFit = fit;
						lowestHist = ph;
					} else {
						auto fit = extractSubtraction(phf1, namepf.Data(), lowestFit);
						v2template = fit->GetParameter(3);
						v2templateerr = fit->GetParError(3);
						scaling = fit->GetParameter(7);
						scalingerr = fit->GetParError(7);

						// convert to LM subtraction by rescaling
						float factor = fit->GetParameter(1) / (fit->GetParameter(1) + scaling * fit->GetParameter(9));
						v2lmsub = v2template * factor;
						v2lmsuberr = v2templateerr * factor; // TODO uncertainty propagation of param 1 and 9

						// draw a subtracted 2d histogram with the scale factor (looks horrible in pp, let's see in OO)
						if (false) {
							new TCanvas;
							auto clone = (TH2*) ph->Clone("clone");
							clone->Add(lowestHist, -1.0 * scaling);
							clone->Draw("SURF1");
							return;
						}
					}
					printf("Fitting %d %d %d --> v22(template) = %.4f +- %.4f | v22(LM subtraction) = %.4f +- %.4f | scaling: %.2f +- %.2f\n", itrig, iassoc, ib, v2template, v2templateerr, v2lmsub, v2lmsuberr, scaling, scalingerr);

					if(imass >= 0){
						gpsimMass->GetXaxis()->SetRangeUser(pmassBins[imass],pmassBins[imass+1]);
						double meanMass = gpsimMass->GetMean();
						pvsbMass->SetPoint(imass,meanMass,v2template); //TODO: the point probably has to be set correctly according to the distribution of D0 (mean of the bins)
						pvsbMass->SetPointError(imass,0.0f,v2templateerr);
					} else {
						// TODO here we still lack the extraction of the reference flow over the full pT range
						
						v22NoMassTemplate->SetBinContent(itrig+1, iassoc+1, ib+1, v2template);
						v22NoMassTemplate->SetBinError(itrig+1, iassoc+1, ib+1, v2templateerr);

						v22NoMassLMSubtraction->SetBinContent(itrig+1, iassoc+1, ib+1, v2lmsub);
						v22NoMassLMSubtraction->SetBinError(itrig+1, iassoc+1, ib+1, v2lmsuberr);
					}
				}
				pvsbMass->Write(Form("v2sb_%u_%u_%u",itrig,iassoc,ib));

				if (Nmass == 0)
					continue;

				// -- MASS fit -- Verbatim code from JP, NOT YET CHECKED

				//simultaneous mass + vn sig+bkg fit
				gpvsb = pvsbMass;

				//stage 1: mass fit only
				pmin->Clear();
				pmin->SetFunction(f);
				pmin->SetVariable(0,"bgC",500.0,10.0);
				pmin->SetVariable(1,"bgB",500.0,10.0);
				pmin->SetVariable(2,"bgA",500.0,10.0);
				pmin->SetVariable(3,"gausA",500.0,0.1);
				pmin->SetVariable(4,"mean",1.86,0.01);
				pmin->SetVariableLimits(4,1.85,1.87);
				pmin->SetVariable(5,"sigma",0.01,0.01);
				pmin->SetVariableLimits(5,0.0001,0.1);//0.001,0.1);
				/*pmin->SetVariable(6,"vnsig",0.1,0.01);
				pmin->SetVariableLowerLimit(6,0.0);
				pmin->SetVariable(7,"vnbkgA",0.1,0.01);
				pmin->SetVariable(8,"vnbkgB",0.1,0.01);*/
				//pmin->SetPrintLevel(0);
				printf("------------------------------------------\n");
				pmin->Minimize();

				//stage 2: vn sig+bkg with fixed mass
				pmin->SetFunction(fStage2);
				for(uint vi = 0; vi < 6; ++vi)
					pmin->FixVariable(vi);
				pmin->SetVariable(6,"vnsig",0.01,0.01);
				pmin->SetVariable(7,"vnbkgA",0.01,0.01);
				pmin->SetVariable(8,"vnbkgB",0.01,0.01);
				pmin->Minimize();

				//stage 3: simultaneous refinement
				pmin->SetFunction(fStage3);
				for(uint vi = 0; vi < 6; ++vi)
					pmin->ReleaseVariable(vi);
				pmin->Minimize();
				printf("------------------------------------------\n");

				//separate the mass fit components to signal and background
				TF1 fitMass(Form("fitMassSim_%u_%u_%u",itrig,iassoc,ib),"[0]+[1]*x+[2]*x*x+gaus(3)");
				fitMass.SetLineColor(kBlue);
				fitMass.SetNpx(10000);
				for(uint i = 0; i < 6; ++i)
					fitMass.SetParameter(i,pmin->X()[i]);
				fitMass.Write();

				TF1 fitBkg(Form("fitBkgSim_%u_%u_%u",itrig,iassoc,ib),"[0]+[1]*x+[2]*x*x");
				fitBkg.SetLineColor(kRed);
				fitBkg.SetNpx(10000);
				for(uint i = 0; i < 3; ++i)
					fitBkg.SetParameter(i,pmin->X()[i]);
				fitBkg.Write();

				//Create the V2signal+V2bkg decomposition
				TF1 fitVnTot(Form("fitVnTot_%u_%u_%u",itrig,iassoc,ib),"([6]*gaus(3)/(gaus(3)+([0]+[1]*x+[2]*x*x))+([7]*x+[8])*([0]+[1]*x+[2]*x*x)/(gaus(3)+([0]+[1]*x+[2]*x*x)))",massLimits[0],massLimits[1]+0.01);
				fitVnTot.SetLineColor(kRed);
				fitVnTot.SetNpx(10000);
				for(uint i = 0; i < 9; ++i){
					fitVnTot.SetParameter(i,pmin->X()[i]);
					fitVnTot.SetParError(i,pmin->Errors()[i]);
				}
				fitVnTot.Write();
				//printf("mean---%f\nsigma--%f\n",pmin->X()[4],pmin->X()[5]);

				printf("^^^Magnitude at D0 mass: %lf, V2 = %lf pm %lf\n",fitVnTot.Eval(1.8),fitVnTot.GetParameter(6),fitVnTot.GetParError(6));

				//calculate the jet yield needed for non-flow subtraction
				double DeltaEtaFitRange = 1.7;
				TF1 fitJet(Form("fitJet_%u_%u_%u",itrig,iassoc,ib),
					"[3]+[2]*[0]/(2*[1]*TMath::Gamma(1/[0]))*TMath::Exp(-TMath::Power(TMath::Abs(x)/[1],[0]))",-DeltaEtaFitRange,DeltaEtaFitRange);
				fitJet.SetParameter(0,0.1);
				fitJet.SetParameter(1,0.1);
				fitJet.SetParameter(2,0.01);
				if(itrig >= 1){
					fitJet.SetParLimits(0,1e-5,100.0);
					fitJet.SetParLimits(2,1e-5,100.0);
				}
				fitJet.SetParameter(3,pjetYieldRaw->GetMinimum(1e-5));
				pjetYieldRaw->Fit(&fitJet,"0QSE","",-DeltaEtaFitRange,DeltaEtaFitRange);
				double dzyam = fitJet.GetParameter(3);
				for(uint i = 1; i <= pjetYieldRaw->GetXaxis()->GetNbins(); ++i){
					double y = pjetYieldRaw->GetBinContent(i);
					pjetYieldRaw->SetBinContent(i,y-dzyam);
				}
				pjetYieldRaw->Write(Form("proj_deta_allMass_Dzyam_%u_%u_%u",itrig,iassoc,ib));
				printf("**** DZYAM %u, %u, %u = %lf ||| %lf, %lf, %lf, %lf\n",itrig,iassoc,ib,dzyam,fitJet.GetParameter(0),fitJet.GetParameter(1),fitJet.GetParameter(2),fitJet.GetParameter(3));
				fitJet.SetParameter(3,0.0);
				fitJet.Write(Form("fitJet_deta_allMass_%u_%u_%u",itrig,iassoc,ib));
				//grv22
				grv22.SetPoint(ib,pmeanCents[ib],fitVnTot.GetParameter(6));
				grv22.SetPointError(ib,0.0,fitVnTot.GetParError(6));

				//ASSOC needs to be event normalized
				TH1 *passocYield = pf->Get<TH1>(Form("assocYield_%u",ib));
				double Nassoc = passocYield->Integral(passocYield->FindBin(passocPt[iassoc]+0.01),passocYield->FindBin(passocPt[iassoc+1]-0.01));
				Nassoc /= pmult->Integral(pmult->FindBin(pcentBins[ib]),pmult->FindBin(pcentBins[ib+1]));
				double Yjet = pjetYieldRaw->Integral(pjetYieldRaw->FindBin(-absDeltaEtaMax),pjetYieldRaw->FindBin(absDeltaEtaMax));

				grnfNAssoc.SetPoint(ib,pmeanCents[ib],Nassoc);
				grnfYJet.SetPoint(ib,pmeanCents[ib],Yjet);

				if(ib < Ncent-1){
					grv22nfsub.SetPoint(ib,pmeanCents[ib],fitVnTot.GetParameter(6)-TempSub*Yjet/Nassoc);
					double e = fitVnTot.GetParError(6);
					double t = TempSubErr*Yjet/Nassoc;
					grv22nfsub.SetPointError(ib,0.0,TMath::Sqrt(e*e+t*t));
				}else{
					TempSub = fitVnTot.GetParameter(6)*Nassoc/Yjet;
					TempSubErr = fitVnTot.GetParError(6)*Nassoc/Yjet;
				}

				delete pvsbMass;
			}
			grv22.Write(Form("V2_%u_%u",itrig,iassoc)); //v22 from the fits
			grv22nfsub.Write(Form("V2NFsub_%u_%u",itrig,iassoc)); //v22 with the CMS non-flow subtraction
			grnfNAssoc.Write(Form("nfNassoc_%u_%u",itrig,iassoc)); //Assoc yield
			grnfYJet.Write(Form("nfYJet_%u_%u",itrig,iassoc)); //Jet yield
		}
	}

	// delete []pcentBins;
	// delete []pmassBins;
	// delete []pmeanCents;
	// delete []passocPt;

	v22NoMassTemplate->Write();
	v22NoMassLMSubtraction->Write();

	pfout->Close();
	delete pfout;

	pf->Close();
	delete pf;
}

