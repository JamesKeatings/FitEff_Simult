#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMatrixTSym.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "Fit/FitResult.h"
#include "TLegend.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

const int neffpars = 5;
const int nnormpars = 2;
const int npars = neffpars + nnormpars;
const int E0 = 325.;
const int nbin1 = 16;
const int nbin2 = 7;
const int nbin3 = 3;
const int Estart = 1;
const int Eend = 4400;

int iparEff[neffpars] = { 0, 1, 2, 3, 4 };
int iparEffBa[neffpars+1] = { 0, 1, 2, 3, 4, 5 };
int iparEffBeam[neffpars+1] = { 0, 1, 2, 3, 4, 6 };

struct GlobalChi2 { 
	GlobalChi2(  ROOT::Math::IMultiGenFunction & f1_Eu,
		ROOT::Math::IMultiGenFunction & f1_Ba,
		ROOT::Math::IMultiGenFunction & f1_Beam,
		ROOT::Math::IMultiGenFunction & f1_BaNorm,
		ROOT::Math::IMultiGenFunction & f1_BeamNorm) :
	fChi2_1(&f1_Eu), fChi2_2(&f1_Ba), fChi2_3(&f1_Beam), fChi2_4(&f1_BaNorm) , fChi2_5(&f1_BeamNorm) {}

	double operator() (const double *p) const {
		double p1[neffpars];
		for (int i = 0; i < neffpars; ++i) p1[i] = p[iparEff[i]];
		double p2[neffpars+1];
		for (int i = 0; i < neffpars+1; ++i) p2[i] = p[iparEffBa[i]];
		double p3[neffpars+1];
		for (int i = 0; i < neffpars+1; ++i) p3[i] = p[iparEffBeam[i]];
		
		double p4[1], p5[1];
		p4[0] = p[iparEffBa[5]];
		p5[0] = p[iparEffBeam[5]];

      		return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3) + (*fChi2_4)(p4) + (*fChi2_5)(p5);
   	}

	const  ROOT::Math::IMultiGenFunction * fChi2_1;
	const  ROOT::Math::IMultiGenFunction * fChi2_2;
	const  ROOT::Math::IMultiGenFunction * fChi2_3;
	const  ROOT::Math::IMultiGenFunction * fChi2_4;
	const  ROOT::Math::IMultiGenFunction * fChi2_5;
};


double ExpFit( double *x, double *par ) {
	Double_t f = par[0];
	for( int i = 1; i < neffpars; i++ )
		f += par[i] * pow(TMath::Log(x[0]/E0),i);
	f = TMath::Exp(f);
	
	return f;
}

double ExpFitBa( double *x, double *par ) {
	Double_t f = par[0];
	for( int i = 1; i < neffpars; i++ )
		f += par[i] * pow(TMath::Log(x[0]/E0),i);
	f = TMath::Exp(f);
	
	return f*par[iparEffBa[5]];
}

double ExpFitBeam( double *x, double *par ) {
	Double_t f = par[0];
	for( int i = 1; i < neffpars; i++ )
		f += par[i] * pow(TMath::Log(x[0]/E0),i);
	f = TMath::Exp(f);
	
	return f*par[iparEffBeam[5]];
}

double NormFuncBa( double *x, double *par ) {
	
	return par[iparEffBa[5]];
}

double NormFuncBeam( double *x, double *par ) {
	
	return par[iparEffBeam[5]];
}

double ExpFitErr( double *x, double *par ) {

	double effpar[neffpars];
	for( int m = 0; m < neffpars; m++ )
		effpar[m] = par[npars*npars+m];

	Double_t f = 0, g = 0;
	Double_t h = pow(ExpFit( x, effpar ),2);
	for( int m = 0; m < neffpars; m++ ) {
		for( int n = 0; n < neffpars; n++ ) {
			g = h*pow(TMath::Log(x[0]/E0),m);
			g *= pow(TMath::Log(x[0]/E0),n);
			g *= par[m*npars+n];
			f += g;
		}
	}
	
	return TMath::Sqrt(f);
}

void FitEff_Simult(){
	ifstream ifile1;
	ifile1.open("EffValues_Eu.dat");
	if(!ifile1.is_open()){cerr << "ifile1 not open" << endl;}
	ifstream ifile2;
	ifile2.open("EffValues_Ba.dat");
	if(!ifile2.is_open()){cerr << "ifile2 not open" << endl;}
	ifstream ifile3;
	ifile3.open("EffValues_InBeam.dat");
	if(!ifile3.is_open()){cerr << "ifile3 not open" << endl;}
	double y1[nbin1], yerr1[nbin1],y2[nbin2],yerr2[nbin2],y3[nbin3],yerr3[nbin3];
	double x1[nbin1],x2[nbin2],x3[nbin3];
	double xerr1[nbin1]={0}; double xerr2[nbin2]={0}; double xerr3[nbin3]={0};
	double a,b,c;
	int i = 0;
	ifile1 >> a >> b >> c;
	while(!ifile1.eof()){
		x1[i] = a;
		y1[i] = b;
		yerr1[i] = c;
		i++;
		ifile1 >> a >> b >> c;
	}
	i = 0;
	ifile2 >> a >> b >> c;
	while(!ifile2.eof()){
		x2[i] = a;
		y2[i] = b;
		yerr2[i] = c;
		i++;
		ifile2 >> a >> b >> c;
	}
	i = 0;
	ifile3 >> a >> b >> c;
	while(!ifile3.eof()){
		x3[i] = a;
		y3[i] = b;
		yerr3[i] = c;
		i++;
		ifile3 >> a >> b >> c;
	}

	int nbin = nbin1 + nbin2 + nbin3;

	TCanvas *c1 = new TCanvas();
	TGraphErrors *gr_Eu = new TGraphErrors(nbin1,x1,y1,xerr1,yerr1);
	TGraphErrors *gr_Ba = new TGraphErrors(nbin2,x2,y2,xerr2,yerr2);
	TGraphErrors *gr_Beam = new TGraphErrors(nbin3,x3,y3,xerr3,yerr3);
	TGraphErrors *gr_BaNorm = new TGraphErrors();
	TGraphErrors *gr_BeamNorm = new TGraphErrors();
	gr_BaNorm->SetPoint( 0, 0, 1.);
	gr_BaNorm->SetPointError( 0, 0, 0.01);
	gr_BeamNorm->SetPoint( 0, 0, 1.);
	gr_BeamNorm->SetPointError( 0, 0, 0.01);
	
	TF1 *f1_Eu = new TF1( "ExpFit_Eu", ExpFit, Estart, Eend, neffpars );
	TF1 *f1_Ba = new TF1( "ExpFit_Ba", ExpFitBa, Estart, Eend, neffpars+1 );
	TF1 *f1_Beam = new TF1( "ExpFit_Bean", ExpFitBeam, Estart, Eend, neffpars+1 );
	TF1 *f1_BaNorm = new TF1( "ExpFit_BaNorm", NormFuncBa, Estart, Eend, 1 );
	TF1 *f1_BeamNorm = new TF1( "ExpFit_BeamNorm", NormFuncBeam, Estart, Eend, 1 );
	TF1 *f2 = new TF1("ExpFitErr", ExpFitErr, Estart, Eend, npars*(npars+1));

	ROOT::Math::WrappedMultiTF1 wfEffEu(*f1_Eu,1);
 	ROOT::Math::WrappedMultiTF1 wfEffBa(*f1_Ba,1);
	ROOT::Math::WrappedMultiTF1 wfEffBeam(*f1_Beam,1);
 	ROOT::Math::WrappedMultiTF1 wfEffBaNorm(*f1_BaNorm,1);
 	ROOT::Math::WrappedMultiTF1 wfEffBeamNorm(*f1_BeamNorm,1);

	ROOT::Fit::DataOptions opt_Eu;
	ROOT::Fit::DataOptions opt_Ba;
	ROOT::Fit::DataOptions opt_Beam;
	ROOT::Fit::DataOptions opt_BaNorm;
	ROOT::Fit::DataOptions opt_BeamNorm; 
	ROOT::Fit::DataRange range; 
	range.SetRange(Estart,Eend);

	ROOT::Fit::BinData data_Eu(opt_Eu,range);
	ROOT::Fit::BinData data_Ba(opt_Ba,range);
	ROOT::Fit::BinData data_Beam(opt_Beam,range);
	ROOT::Fit::BinData data_BaNorm(opt_BaNorm,range);
	ROOT::Fit::BinData data_BeamNorm(opt_BeamNorm,range);

	ROOT::Fit::FillData(data_Eu, gr_Eu);
	ROOT::Fit::FillData(data_Ba, gr_Ba);
	ROOT::Fit::FillData(data_Beam, gr_Beam);
	ROOT::Fit::FillData(data_BaNorm, gr_BaNorm); 
	ROOT::Fit::FillData(data_BeamNorm, gr_BeamNorm); 

	ROOT::Fit::Chi2Function chi2_Eu(data_Eu, wfEffEu);
	ROOT::Fit::Chi2Function chi2_Ba(data_Ba, wfEffBa);
	ROOT::Fit::Chi2Function chi2_Beam(data_Beam, wfEffBeam);
	ROOT::Fit::Chi2Function chi2_BaNorm(data_BaNorm, wfEffBaNorm);
	ROOT::Fit::Chi2Function chi2_BeamNorm(data_BeamNorm, wfEffBeamNorm);

	GlobalChi2 globalChi2(chi2_Eu, chi2_Ba, chi2_Beam, chi2_BaNorm, chi2_BeamNorm);

	ROOT::Fit::Fitter fitter;

	double par0[npars] = {2.,-0.8,-0.04,0.1,0.03,1,1};
	fitter.Config().SetParamsSettings(npars,par0);

	fitter.Config().SetMinimizer("Minuit2","Migrad");
	fitter.FitFCN(npars,globalChi2,0,data_Eu.Size()+data_Ba.Size()+data_Beam.Size()+data_BaNorm.Size()+data_BeamNorm.Size(),true);
	ROOT::Fit::FitResult fitres = fitter.Result();

	fitres.NormalizeErrors();
	double chisq = fitres.Chi2() / fitres.Ndf();

	double parArray[npars*npars+2*npars];
	for( int i = 0; i < npars; i++ ) {
		for( int j = 0; j < npars; j++){
			parArray[i*npars+j] = fitres.CovMatrix(i,j);		
			//cout << fitres.CovMatrix(i,j) << endl;;
		}
	}
	for( int i = 0; i < npars; i++ ) {
		par0[i] = fitres.Value(i);
		parArray[npars*npars+i] = fitres.Value(i);
	}
	f2->SetParameters( parArray );

	TGraph *gFinal = new TGraph(Eend-Estart);
	TGraph *gLow = new TGraph(Eend-Estart);
	TGraph *gUpp = new TGraph(Eend-Estart);
	TGraphAsymmErrors *gShade = new TGraphAsymmErrors(2*Eend-Estart);

	double xval[1];
	double eff, err, eff2, err2;

	for( int i = Estart; i < Eend; i++ ) {
		xval[0] = (double)i;
		eff = ExpFit(xval,par0);
		gFinal->SetPoint(i-Estart, i , eff);
		err = f2->Eval(i);
		if( i == 1905) cout << "efficiency @ 1000 keV = " << eff << " +/- " << err << endl;
		gLow->SetPoint( i-Estart, i, eff-err );
		gUpp->SetPoint( i-Estart, i, eff+err );
	}
	
	double end[1];

	for( int i = Estart; i < Eend; i++ ) {
		xval[0] = (double)i;
		end[0] = (double)Eend - (double)i;
		gShade->SetPoint( Eend+i-Estart, i, (ExpFit(xval,par0))+(f2->Eval(i)) );
		gShade->SetPoint( i-Estart, Eend-i, (ExpFit(end,par0))-(f2->Eval(Eend-i)) );
	}

	double  xvalue, yvalue;
	for( int i = 0; i < gr_Ba->GetN(); i++ ) {
		gr_Ba->GetPoint( i, xvalue, yvalue );
		gr_Ba->SetPoint( i, xvalue, yvalue/par0[iparEffBa[5]] );
		yvalue = gr_Ba->GetErrorY( i );
		gr_Ba->SetPointError( i, 0, yvalue/par0[iparEffBa[5]] );
	}
	for( int i = 0; i < gr_Beam->GetN(); i++ ) {
		gr_Beam->GetPoint( i, xvalue, yvalue );
		gr_Beam->SetPoint( i, xvalue, yvalue/par0[iparEffBeam[5]] );
		yvalue = gr_Beam->GetErrorY( i );
		gr_Beam->SetPointError( i, 0, yvalue/sqrt(par0[iparEffBeam[5]]) );
	}

	gShade->SetFillStyle(1001);
   	//gShade->SetFillColor(46);
	gShade->SetFillColorAlpha(46,0.5);

	gFinal->SetLineColor(2);
	gFinal->SetLineWidth(3);
	gLow->SetLineStyle(10);
	gLow->SetLineColorAlpha(1,0.3);
	gLow->SetLineWidth(2);
	gUpp->SetLineStyle(10);
	gUpp->SetLineColorAlpha(1,0.3);
	gUpp->SetLineWidth(2);
	gr_Eu->SetMarkerColor(8);
	gr_Eu->SetMarkerStyle(24);
	gr_Eu->SetMarkerSize(2);
	gr_Eu->SetLineWidth(2);
	gr_Ba->SetMarkerColor(4);
	gr_Ba->SetMarkerStyle(25);
	gr_Ba->SetMarkerSize(2);
	gr_Ba->SetLineWidth(2);
	gr_Beam->SetMarkerColor(6);
	gr_Beam->SetMarkerStyle(26);
	gr_Beam->SetMarkerSize(2);
	gr_Beam->SetLineWidth(2);

	TMultiGraph *mg = new TMultiGraph();

	TLegend* leg = new TLegend(0.5,0.5, .9, .9);
	leg->AddEntry(gr_Eu, "^{152}Eu source peaks", "lep");
	leg->AddEntry(gr_Ba, "^{133}Ba source peaks", "lep");
	leg->AddEntry(gr_Beam, "^{132}Sn in-beam peaks", "lep");
	leg->AddEntry(gFinal, "Fit results", "l");
	leg->AddEntry(gShade, "Error region", "f");
	leg->AddEntry(gLow, "Upper/lower boundaries", "l");

	mg->Add(gShade,"f");
	mg->Add(gLow,"Cf");
	mg->Add(gUpp,"Cf");
	mg->Add(gFinal,"C");
	mg->Add(gr_Eu,"P");
	mg->Add(gr_Ba,"P");
	mg->Add(gr_Beam,"P");
	mg->Draw("A");
	mg->GetXaxis()->SetRangeUser(0,4500);
	mg->GetYaxis()->SetRangeUser(0,20);
	c1->SetGridy();
	c1->SetGridx();

	leg->Draw("same");

	mg->GetXaxis()->SetTitle("Energy (keV)");
	mg->GetYaxis()->SetTitle("Efficiency (%)");
	mg->GetXaxis()->SetTickLength(0.015);
	mg->GetYaxis()->SetTickLength(0.015);
	mg->GetXaxis()->SetTitleSize(0.045);
	mg->GetYaxis()->SetTitleSize(0.045);
	mg->GetYaxis()->SetTitleOffset(0.75);
	c1->Update();

	//fitres.Print(std::cout,false);

	for( int i = 0; i < npars; i++ ) {
			cout << fitres.GetParameterName(i) << " = " << fitres.Value(i);
			cout << " (+" << fitres.UpperError(i) << ";-" << fitres.LowerError(i);
			cout << ")" << endl;
	}
	cout << "Chi2/NDF = " << chisq << endl;
}
