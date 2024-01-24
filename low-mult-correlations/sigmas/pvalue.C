#include <stdio.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>

//first row mult, second row 95% CL
// CR1 and CR2 values
// const static double ALEPHdata[][4] = {
// 	{8.818258164,0.001201959,8.915567601,0.000114976},
// 	{15.80360608,9.68782e-06,15.81986143,9.57e-07},
// 	{23.38995667,5.203e-05,23.40621203,5.1374e-06},
// 	{32.57823882,0.002266584,32.5942715,0.000231013},
// 	{37.17649941,0.008319748,37.19253208,0.000847959}
// };
// same as above, adjusted for third row
// const static double ALEPHdata[][4] = {
// 	{8.818258164,0.001201959,	0.95, 0},
// 	{15.80360608,9.68782e-06,	0.95, 0},
// 	{23.38995667,5.203e-05,		0.95, 0},
// 	{32.57823882,0.002266584,	0.95, 0},
// 	{37.17649941,0.008319748,	0.95, 0}
// };

// Hepdata, thrust axis: https://www.hepdata.net/record/ins1737859
// const static double ALEPHdata[][4] = {
// 	{8.9,  0.001195229,  0.95, 0},
// 	{15.8, 1e-5,         0.96, 0},
// 	{23.4, 5.221903e-05, 0.95, 0},
// 	{32.6, 0.002284824,  0.95, 0},
// 	{37.2, 0.008600333,  0.95, 0}
// };

// Hepdata, lab axis: https://www.hepdata.net/record/ins1737859
// const static double ALEPHdata[][4] = {
// 	{8.9,  1e-5,         0.963, 0},
// 	{15.8, 1e-5,         0.9997, 0},
// 	{23.4, 4.449485e-05, 0.95, 0},
// 	{32.6, 0.0003234487, 0.95, 0},
// 	{37.2, 0.07806029,   0.95, 0}
// };

// By hand read of https://arxiv.org/pdf/2312.05084.pdf
const static double ALEPHdata[][4] = {
	{15.8, 1.0e-7, 0.99,  0},
	{23.7, 1.0e-7, 0.984, 0},
	{32.9, 1.0e-7, 0.99,  0},
	{43.3, 4.57e-3, 0.95,  0},
	{53.8, 3.8372e-2+3.7494123e-2, 0.95,  0} //2sigma
};

std::tuple<double, double, double> Interpolate(double ALICEmult, const TGraphErrors *pstat, const TGraphErrors *psyst){
	double y = pstat->Eval(ALICEmult);
	TGraph yerrgr(pstat->GetN(),pstat->GetX(),pstat->GetEY());
	double yerr = yerrgr.Eval(ALICEmult);
	TGraph ysysgr(psyst->GetN(),psyst->GetX(),psyst->GetEY());
	double ysys = ysysgr.Eval(ALICEmult);
	return std::make_tuple(y,yerr,ysys);
}

/*
Parallelization might be needed for the new ALEPH result for which passing the delta conditions is extremely rare.
This can be done as below by summing the counts for individual jobs.

seq 1001 1032 | xargs -P32 -I{} root -l -q 'pvalue.C(0,5,{},1000000000)' | grep -Po '[0-9]+/[0-9]+' | awk -F'/' 'BEGIN{sumc=0;sumt=0;}{sumc+=$1;sumt+=$2;print "Done: " $0}END{print "Total " sumc "/" sumt; system("root -l -q -e '"'"'printf(\"%lf\",TMath::ErfInverse(1.0-(double)" sumc "/(double)" sumt ")*sqrt(2.0));'"'"'");}'
*/
void pvalue(int from, int to, uint seed=1000, uint Niter = 1000000000){
	gRandom->SetSeed(seed);
	TFile *pf = new TFile("/tmp/bootstrapGraphsMerged.root","read");
	TGraphErrors *pstat = (TGraphErrors*)pf->Get("gr_boot_stat_max");
	TGraphErrors *psyst = (TGraphErrors*)pf->Get("gr_boot_together_max_reduced");
	TGraphErrors *ptot = (TGraphErrors*)pf->Get("gr_boot_together_max");

	//ALEPH -> ALICE multiplicity conversion factors
	double acceptanceCorr_pp = 0.576209; //pythia pp
	//double acceptanceCorr_ee = 0.786990; //pythia ee (91 GeV)
	double acceptanceCorr_ee = 0.70; //pythia ee (209 GeV)
	
	//double acceptanceCorr = acceptanceCorr_ee;
	double acceptanceCorr = acceptanceCorr_pp;

	const UInt_t n = sizeof(ALEPHdata)/(4*8); //number of ALEPH (interp.) points
	double y[n], yerr[n], ysys[n], alephSigma68[n];
	for(uint i = 0; i < n; ++i){
		alephSigma68[i] = ALEPHdata[i][1]/(1.4142*TMath::ErfInverse(ALEPHdata[i][2])); //Half-gaussian sigma
		double x1 = acceptanceCorr*ALEPHdata[i][0];
		std::tie(y[i],yerr[i],ysys[i]) = Interpolate(x1,pstat,psyst);
		if (false) { // consistency check 
			y[i] = alephSigma68[i];
			yerr[i] = 0;
			ysys[i] = 0;
		}
		//auto [y1,yerr1,ysys1] = Interpolate(x1,pstat,psyst);
	}

	uint c3 = 0;
	for(UInt_t i = 0; i < Niter; ++i){
		double s = gRandom->Gaus(0.0,1.0);
		bool pass = true;
		for(int j = from; j < to; ++j){
			if (true) { // JF
				double Y = y[j] + gRandom->Gaus(0.0,yerr[j])+s*ysys[j];
				double a = gRandom->Gaus(0.0,alephSigma68[j]);
				if (a < Y){ //all deltas must be > y
					pass = false;
					break;
				}
				// if (i % 10000 == 0) {
				// 	printf("%d | %f \t %f | %d\n", j, Y, a, pass);
				// }
			} else { // Jasper
				double Y = gRandom->Gaus(0.0,yerr[j])+s*ysys[j];
				double a = gRandom->Gaus(0.0,alephSigma68[j]);
				double delta = Y-a; //TMath::Abs(Y-a);
				if(delta <= y[j]){ //all deltas must be > y
					pass = false;
					break;
				}
			}
		}
    
    	if(pass)
    		c3++;
	}

	printf("NSigma=%lf, (p-combination) (counts = %u/%u)\n",
		TMath::ErfInverse(1.0-(double)c3/(double)Niter)*sqrt(2),c3,Niter);

	pf->Close();
	delete pf;
}

/*int pvalue() {
	//for (int i=0; i<5; i++)
	//	pvalue5(i, i+1);

	printf("All points combined:\n");
	pvalue5(0, 4);
	// pvalue5(1, 3);

	return 0;
}*/

