#include <fstream>
#include <string>
#include <sstream>

#include "TFile.h"
#include "TH2D.h"
#include "TPRegexp.h"
#include "TObjArray.h"

void
xsec()
{
  using namespace std;

  TFile* output = new TFile("xsec.root", "recreate");
  TH2D* gsq_B_LO = new TH2D("gsq_B_LO", "LO;squark mass (GeV);gluino mass (GeV);xsec (pb)", 17, 350., 2050., 17, 350., 2050.);
  TH2D* gsq_W_LO = new TH2D("gsq_W_LO", "LO;squark mass (GeV);gluino mass (GeV);xsec (pb)", 17, 350., 2050., 17, 350., 2050.);
  TH2D* gsq_B_NLO = new TH2D("gsq_B_NLO", "NLO;squark mass (GeV);gluino mass (GeV);xsec (pb)", 17, 350., 2050., 17, 350., 2050.);
  TH2D* gsq_W_NLO = new TH2D("gsq_W_NLO", "NLO;squark mass (GeV);gluino mass (GeV);xsec (pb)", 17, 350., 2050., 17, 350., 2050.);

  TPRegexp pat("^[ ]*[0-9]+[ ]+([0-9]+)[ ]+([0-9]+)[ ]+[0-9]+[ ]+[0-9]+[ ]LO:[ ]([0-9e.+-]+)[ ][+][ ]([0-9e.+-]+)[ ]-[ ]([0-9e.+-]+)[ ]NLO:[ ]([0-9e.+-]+)[ ]([0-9e.+-]+)[ ]([0-9e.+-]+)");
  std::string line;

  std::ifstream Binput("Spectra_gsq_B_8TeV.xsec");
  while(std::getline(Binput, line), Binput.good()){
    TObjArray* arr = pat.MatchS(line.c_str());
    if(arr->GetEntries() != 9) continue;
    double msq = TString(arr->At(1)->GetName()).Atof();
    double mgl = TString(arr->At(2)->GetName()).Atof();
    double x = TString(arr->At(3)->GetName()).Atof();
    double eu = TString(arr->At(4)->GetName()).Atof();
    double el = TString(arr->At(5)->GetName()).Atof();
    int bin = gsq_B_LO->FindBin(msq, mgl);
    gsq_B_LO->SetBinContent(bin, x);
    gsq_B_LO->SetBinError(bin, (eu + el) / 2.);

    x = TString(arr->At(6)->GetName()).Atof();
    eu = TString(arr->At(7)->GetName()).Atof();
    el = TString(arr->At(8)->GetName()).Atof();
    gsq_B_NLO->SetBinContent(bin, x);
    gsq_B_NLO->SetBinError(bin, (eu + el) / 2.);
  }
  Binput.close();

  std::ifstream Winput("Spectra_gsq_W_8TeV.xsec");
  while(std::getline(Winput, line), Winput.good()){
    TObjArray* arr = pat.MatchS(line.c_str());
    if(arr->GetEntries() != 9) continue;
    double msq = TString(arr->At(1)->GetName()).Atof();
    double mgl = TString(arr->At(2)->GetName()).Atof();
    double x = TString(arr->At(3)->GetName()).Atof();
    double eu = TString(arr->At(4)->GetName()).Atof();
    double el = TString(arr->At(5)->GetName()).Atof();
    int bin = gsq_W_LO->FindBin(msq, mgl);
    gsq_W_LO->SetBinContent(bin, x);
    gsq_W_LO->SetBinError(bin, (eu + el) / 2.);

    x = TString(arr->At(6)->GetName()).Atof();
    eu = TString(arr->At(7)->GetName()).Atof();
    el = TString(arr->At(8)->GetName()).Atof();
    gsq_W_NLO->SetBinContent(bin, x);
    gsq_W_NLO->SetBinError(bin, (eu + el) / 2.);
  }
  Winput.close();

  output->Write();

  delete output;
}
