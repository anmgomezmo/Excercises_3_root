#include <iostream>
#include <string>
#include <stdio.h>
#include <vector>

void TTbarSel_8TeV_MC_Data(){ //void function with the same name as script or file code

TFile *file = TFile::Open("../Data_8TeV.root");

TTree *tree = (TTree*) file->Get("mini");
//tree->Print();

// activate variables to use

//scalar variables
Bool_t e_trig;
Bool_t mu_trig;
Bool_t good_vtx;
UInt_t lep_n;
UInt_t jet_n;
Float_t MET;
Float_t MET_phi;
Float_t evtw;

// vectorial or array variables
Float_t lep_pt[10];
Float_t lep_eta[10];
Float_t lep_phi[10];
Float_t lep_E[10];
Float_t lep_ptcone30[10];
Float_t lep_etcone20[10];
Int_t lep_type[10];

Float_t jet_pt[10];
Float_t jet_eta[10];
Float_t jet_jvf[10];
Float_t jet_MV1[10];

// assign address of data to variables
tree->SetBranchAddress("trigE", &e_trig);
tree->SetBranchAddress("trigM", &mu_trig);
tree->SetBranchAddress("hasGoodVertex", &good_vtx);
tree->SetBranchAddress("lep_n", &lep_n);
tree->SetBranchAddress("jet_n", &jet_n);
tree->SetBranchAddress("met_et", &MET);
tree->SetBranchAddress("met_phi", &MET_phi);
tree->SetBranchAddress("mcWeight", &evtw);

tree->SetBranchAddress("lep_pt", &lep_pt);
tree->SetBranchAddress("lep_eta", &lep_eta);
tree->SetBranchAddress("lep_phi", &lep_phi);
tree->SetBranchAddress("lep_E", &lep_E);
tree->SetBranchAddress("lep_type", &lep_type);
tree->SetBranchAddress("lep_ptcone30", &lep_ptcone30);
tree->SetBranchAddress("lep_etcone20", &lep_etcone20);

tree->SetBranchAddress("jet_pt", &jet_pt);
tree->SetBranchAddress("jet_eta", &jet_eta);
tree->SetBranchAddress("jet_jvf", &jet_jvf);
tree->SetBranchAddress("jet_MV1", &jet_MV1);


//create canvas
TCanvas *canvas = new TCanvas("Canvas","",800,600);

//create leading jet pt and all jets histograms
TH1F *cutflow = new TH1F("Cutflow","Cutflow; Cut; Events",10,0,10);
TFile f("hist_njets_data.root","new");
TH1F *hist_njets_data = new TH1F("Number of jets","n-jets; Jet multiplicity; Events",10,0,10);
//TH1F *aux = new TH1F("Number of good leptons","n-leptons; n-leptons; Events",10,0,10);



//loop and fill the histograms
int nentries, nbytes, i;
nentries = (Int_t)tree->GetEntries();

int cut1 = 0;
int cut2 = 0;
int cut3 = 0;
int cut4 = 0;
int cut5 = 0;
int cut6 = 0;
int cut7 = 0;
int cut8 = 0;



for (i = 0; i < nentries; i++) {

  nbytes = tree->GetEntry(i);
  vector <int> index;

  //First cut: Good vertex
  if (!good_vtx) continue;
  cut1++;
  cutflow->Fill(1);

  //Second cut: Trigger
  if (!e_trig && !mu_trig) continue;
  cut2++;
  cutflow->Fill(2);

  //Preselection of good leptons
  int n_mu = 0; //number of good muons
  int n_el = 0; //number of good electrons
  int n_lep = 0; //number of good leptons


  //Loop over leptons
  for (unsigned i = 0; i < lep_n; i++) {
    if (lep_pt[i] < 25000.0) continue;
    if (lep_ptcone30[i]/lep_pt[i] > 0.15) continue;
    if (lep_etcone20[i]/lep_pt[i] > 0.15) continue;
    if (lep_type[i] == 13 && TMath::Abs(lep_eta[i]) < 2.5) {
      n_mu++;
    }
    //To complete: Add electrons and extract the index for the good lepton
    float center = (1.52+1.37)/2;
    float width = (1.52-1.37)/2;
    if (lep_type[i] == 11 && TMath::Abs(lep_eta[i]) < 2.47 ) {
      if (TMath::Abs(TMath::Abs(lep_eta[i]) - center) > width) {
        n_el++;
      }
    }
    index.push_back(i);
    n_lep++;
  }

  /*std::cout << index.size() << std::endl;
  aux->Fill(index.size());*/

  //Select events with only 1 good lepton and fill the cutflow histogram
  //Third cut (one good lepton):
  if (n_lep != 1) continue;
  cutflow->Fill(3);
  cut3++;


  int n_jets = 0;
  int n_bjets = 0;

  //Number of jets distribution
  //hist_njets->Fill(jet_n);

  //Fourth cut: At least 4 jets
  if (jet_n < 4) continue;
  cutflow->Fill(4);
  cut4++;

  for (unsigned j = 0; j < jet_n; j++) {
    // To complete: apply jet cuts to find the good jets
    if (jet_pt[j] < 25000.) continue;
    if (TMath::Abs(jet_eta[j]) > 2.5) continue; //Eta cut
    //JVF cut
    if (jet_pt[j] < 50000. && TMath::Abs(jet_eta[j]) < 2.4 && TMath::Abs(jet_jvf[j]) > 0.5) {
      n_jets++;
    }
    if (jet_pt[j] > 50000. && TMath::Abs(jet_jvf[j]) > 0.5) {
      n_jets++;
    }
    if (jet_MV1[j] < 0.7892) continue; //cut on 0.7892 MV1 and count the number of b-jets
      n_bjets++;
  }

  //Fifth cut: At least 4 good jets
  if (n_jets < 4) continue;
  cutflow->Fill(5);
  cut5++;

  //Sixth cut: at least two b-jet
  if (n_bjets < 2) continue;
  cutflow->Fill(6);
  cut6++;

  //Seventh cut: MET > 30 GeV
  if (MET < 30000.) continue;
  cutflow->Fill(7);
  cut7++;

  //TLorentzVector definitions
  TLorentzVector Lepton = TLorentzVector();
  TLorentzVector MeT = TLorentzVector();

  //To complete: Lorentz vectors for the lepton and MET. Use SetPtEtaPhiE().
  Lepton.SetPtEtaPhiE(lep_pt[index[0]],lep_eta[index[0]],lep_phi[index[0]],lep_E[index[0]]);
  MeT.SetPtEtaPhiE(0,0,MET_phi,MET);

  //Calculation of the mTW using TLorenzt vectors
  float mTW = sqrt(2*Lepton.Pt()*MeT.E()*(1-cos(Lepton.DeltaPhi(MeT))));

  //Eight cut: mTW > 30 GeV
  if (mTW < 30000.) continue;
  cutflow->Fill(8);
  cut8++;



  index.clear();

  hist_njets_data->Fill(jet_n);

}

hist_njets_data->Write();

std::cout << "Done!" << std::endl;
std::cout << "All events:" << "\t" << nentries << std::endl;
std::cout << "Cut1:" << "\t" << cut1 << std::endl;
std::cout << "Cut2:" << "\t" << cut2 << std::endl;
std::cout << "Cut3:" << "\t" << cut3 << std::endl;
std::cout << "Cut4:" << "\t" << cut4 << std::endl;
std::cout << "Cut5:" << "\t" << cut5 << std::endl;
std::cout << "Cut6:" << "\t" << cut6 << std::endl;
std::cout << "Cut7:" << "\t" << cut7 << std::endl;
std::cout << "Cut8:" << "\t" << cut8 << std::endl;

/*aux->Draw();
canvas->Draw();*/


//Draw histograms
/*cutflow->Draw();
canvas->Draw();*/

/*hist_njets_data->Draw();
canvas->Draw();*/

//////////////////////////////////////////////////////////////////////////////////////


TFile *file1 = TFile::Open("../ttbar_8TeV.root");
TTree *tree1 = (TTree*) file1->Get("mini");
//tree1->Print();

// activate variables to use

//scalar variables
Bool_t e_trig1;
Bool_t mu_trig1;
Bool_t good_vtx1;
UInt_t lep_n1;
UInt_t jet_n1;
Float_t MET1;
Float_t MET_phi1;
Float_t evtw1;

// vectorial or array variables
Float_t lep_pt1[10];
Float_t lep_eta1[10];
Float_t lep_phi1[10];
Float_t lep_E1[10];
Float_t lep_ptcone301[10];
Float_t lep_etcone201[10];
Int_t lep_type1[10];

Float_t jet_pt1[10];
Float_t jet_eta1[10];
Float_t jet_jvf1[10];
Float_t jet_MV11[10];

Float_t scaleFactor_PILEUP;

// assign address of data to variables
tree1->SetBranchAddress("trigE", &e_trig1);
tree1->SetBranchAddress("trigM", &mu_trig1);
tree1->SetBranchAddress("hasGoodVertex", &good_vtx1);
tree1->SetBranchAddress("lep_n", &lep_n1);
tree1->SetBranchAddress("jet_n", &jet_n1);
tree1->SetBranchAddress("met_et", &MET1);
tree1->SetBranchAddress("met_phi", &MET_phi1);
tree1->SetBranchAddress("mcWeight", &evtw1);

tree1->SetBranchAddress("lep_pt", &lep_pt1);
tree1->SetBranchAddress("lep_eta", &lep_eta1);
tree1->SetBranchAddress("lep_phi", &lep_phi1);
tree1->SetBranchAddress("lep_E", &lep_E1);
tree1->SetBranchAddress("lep_type", &lep_type1);
tree1->SetBranchAddress("lep_ptcone30", &lep_ptcone301);
tree1->SetBranchAddress("lep_etcone20", &lep_etcone201);

tree1->SetBranchAddress("jet_pt", &jet_pt1);
tree1->SetBranchAddress("jet_eta", &jet_eta1);
tree1->SetBranchAddress("jet_jvf", &jet_jvf1);
tree1->SetBranchAddress("jet_MV1", &jet_MV11);
tree1->SetBranchAddress("scaleFactor_PILEUP", &scaleFactor_PILEUP);

//create leading jet pt and all jets histograms
TH1F *cutflow1 = new TH1F("Cutflow","Cutflow; Cut; Events",10,0,10);
TH1F *hist_njets_mc = new TH1F("Number of jets","n-jets; Jet multiplicity; Events",10,0,10);
//TH1F *aux = new TH1F("Number of good leptons","n-leptons; n-leptons; Events",10,0,10);



//loop and fill the histograms
int nentries1, nbytes1, i1;
nentries1 = (Int_t)tree1->GetEntries();

int cut11 = 0;
int cut12 = 0;
int cut13 = 0;
int cut14 = 0;
int cut15 = 0;
int cut16 = 0;
int cut17 = 0;
int cut18 = 0;


for (i1 = 0; i1 < nentries1; i1++) {

  nbytes1 = tree1->GetEntry(i1);
  vector <int> index1;

  //First cut: Good vertex
  if (!good_vtx1) continue;
  cut11++;
  cutflow1->Fill(1);

  //Second cut: Trigger
  if (!e_trig1 && !mu_trig1) continue;
  cut12++;
  cutflow1->Fill(2);

  //Preselection of good leptons
  int n_mu1 = 0; //number of good muons
  int n_el1 = 0; //number of good electrons
  int n_lep1 = 0; //number of good leptons

  //Loop over leptons
  for (unsigned i = 0; i < lep_n1; i++) {
    if (lep_pt1[i] < 25000.0) continue;
    if (lep_ptcone301[i]/lep_pt1[i] > 0.15) continue;
    if (lep_etcone201[i]/lep_pt1[i] > 0.15) continue;
    if (lep_type1[i] == 13 && TMath::Abs(lep_eta1[i]) < 2.5) {
      n_mu1++;
    }
    //To complete: Add electrons and extract the index for the good lepton
    float center = (1.52+1.37)/2;
    float width = (1.52-1.37)/2;
    if (lep_type1[i] == 11 && TMath::Abs(lep_eta1[i]) < 2.47 ) {
      if (TMath::Abs(TMath::Abs(lep_eta1[i]) - center) > width) {
        n_el1++;
      }
    }
    index1.push_back(i);
    n_lep1++;
  }

  /*std::cout << index.size() << std::endl;
  aux->Fill(index.size());*/

  //Select events with only 1 good lepton and fill the cutflow histogram
  //Third cut (one good lepton):
  if (n_lep1 != 1) continue;
  cutflow1->Fill(3);
  cut13++;


  int n_jets1 = 0;
  int n_bjets1 = 0;

  //Number of jets distribution
  //hist_njets->Fill(jet_n);

  //Fourth cut: At least 4 jets
  if (jet_n1 < 4) continue;
  cutflow1->Fill(4);
  cut14++;

  for (unsigned j = 0; j < jet_n1; j++) {
    // To complete: apply jet cuts to find the good jets
    if (jet_pt1[j] < 25000.) continue;
    if (TMath::Abs(jet_eta1[j]) > 2.5) continue; //Eta cut
    //JVF cut
    if (jet_pt1[j] < 50000. && TMath::Abs(jet_eta1[j]) < 2.4 && TMath::Abs(jet_jvf1[j]) > 0.5) {
      n_jets1++;
    }
    if (jet_pt1[j] > 50000. && TMath::Abs(jet_jvf1[j]) > 0.5) {
      n_jets1++;
    }
    if (jet_MV11[j] < 0.7892) continue; //cut on 0.7892 MV1 and count the number of b-jets
      n_bjets1++;
  }

  //Fifth cut: At least 4 good jets
  if (n_jets1 < 4) continue;
  cutflow1->Fill(5);
  cut15++;

  //Sixth cut: at least two b-jet
  if (n_bjets1 < 2) continue;
  cutflow1->Fill(6);
  cut16++;

  //Seventh cut: MET > 30 GeV
  if (MET1 < 30000.) continue;
  cutflow1->Fill(7);
  cut17++;

  //TLorentzVector definitions
  TLorentzVector Lepton1 = TLorentzVector();
  TLorentzVector MeT1 = TLorentzVector();

  //To complete: Lorentz vectors for the lepton and MET. Use SetPtEtaPhiE().
  Lepton1.SetPtEtaPhiE(lep_pt1[index1[0]],lep_eta1[index1[0]],lep_phi1[index1[0]],lep_E1[index1[0]]);
  MeT1.SetPtEtaPhiE(0,0,MET_phi1,MET1);

  //Calculation of the mTW using TLorenzt vectors
  float mTW1 = sqrt(2*Lepton1.Pt()*MeT1.E()*(1-cos(Lepton1.DeltaPhi(MeT1))));

  //Eight cut: mTW > 30 GeV
  if (mTW1 < 30000.) continue;
  cutflow1->Fill(8);
  cut18++;



  index1.clear();

  Float_t evtw = scaleFactor_PILEUP*(1000.*137.29749)/(49761200.21*0.072212854);
  hist_njets_mc->Fill(jet_n1,evtw);

  }

std::cout << "Done!" << std::endl;
std::cout << "All events:" << "\t" << nentries1 << std::endl;
std::cout << "Cut1:" << "\t" << cut11 << std::endl;
std::cout << "Cut2:" << "\t" << cut12 << std::endl;
std::cout << "Cut3:" << "\t" << cut13 << std::endl;
std::cout << "Cut4:" << "\t" << cut14 << std::endl;
std::cout << "Cut5:" << "\t" << cut15 << std::endl;
std::cout << "Cut6:" << "\t" << cut16 << std::endl;
std::cout << "Cut7:" << "\t" << cut17 << std::endl;
std::cout << "Cut8:" << "\t" << cut18 << std::endl;

/*aux->Draw();
canvas->Draw();*/

//Draw histograms
/*cutflow->Draw();
canvas->Draw();*/

/*hist_njets_mc->Draw();
canvas->Draw();*/

hist_njets_data->SetMarkerStyle(20);
hist_njets_data->SetMarkerColor(kRed);
hist_njets_data->Draw("p");
//hist_njets_mc->SetMarkerStyle(21);
//hist_njets_mc->SetMarkerColor(kBlue);
hist_njets_mc->Draw("same");
canvas->Draw();

}
