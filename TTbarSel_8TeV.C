#include <iostream>
#include <string>
#include <stdio.h>
#include <vector>

void TTbarSel_8TeV(){ //void function with the same name as script or file code

  TFile *file = TFile::Open("../Data_8TeV.root");

  TTree *tree = (TTree*) file->Get("mini");
  tree->Print();

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

  TFile f("hist_data.root","new");
  TH1F *hist_nlep = new TH1F("Number of leptons","n-leptons; Lepton multiplicity; Events",4,0,4);
  TH1F *hist_lep_pt = new TH1F("Leptons pt","Lepton pt; pt (GeV); Events",20,0,200);
  TH1F *hist_lep_trackisolation = new TH1F("Leptons track isolation","Lepton track; ptcone/pt; Events",10,0,0.2);
  TH1F *hist_lep_calorisolation = new TH1F("Leptons calorimeter isolation","Lepton calorimeter; etcone/pt; Events",10,0,0.2);
  TH1F *hist_lep_eta = new TH1F("Leptons eta","Lepton eta; eta; Events",20,-10,10);
  TH1F *hist_lep_e = new TH1F("Leptons E","Lepton E; E (GeV); Events",20,0,200);
  TH1F *hist_nele = new TH1F("Number of electrons","Number of electrons; Electron multiplicity; Events",5,0,5);
  TH1F *hist_nmuon = new TH1F("Number of muons","Number of muons; Muon multiplicity; Events",5,0,5);
  TH1F *hist_jetn = new TH1F("Number of jets","n-jets; Jet multiplicity; Events",10,0,10);
  TH1F *hist_njets = new TH1F("Number of good jets","n-jets; Jet multiplicity; Events",6,4,10);
  TH1F *hist_jet_pt = new TH1F("Jet pt","Jet pt; pt (GeV); Events",20,0,200);
  TH1F *hist_jet_eta = new TH1F("Jet eta","Jet eta; eta; Events",20,-10,10);
  TH1F *hist_jet_jvf = new TH1F("Jet jvf","Jet jvf; jvf; Events",20,-1,2);
  TH1F *hist_jet_mv1 = new TH1F("MV1","MV1; mv1; Events",10,0,1);
  TH1F *hist_nbjets = new TH1F("Number of b jets","n-bjets; Jet multiplicity; Events",4,2,6);
  TH1F *hist_met = new TH1F("MET value","MET; MET (GeV); Events",20,30,200);
  TH1F *hist_mtw = new TH1F("mTW value","mTW; mTW (GeV); Events",20,30,200);

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
    //hist_jetn->Fill(jet_n);

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
    MeT.SetPtEtaPhiE(MET,0,MET_phi,MET);

    //Calculation of the mTW using TLorenzt vectors
    float mTW = sqrt(2*Lepton.Pt()*MeT.E()*(1-cos(Lepton.DeltaPhi(MeT))));

    //Eight cut: mTW > 30 GeV
    if (mTW < 30000.) continue;
    cutflow->Fill(8);
    cut8++;


    // Fill histograms to compare with MC

    hist_nlep->Fill(n_lep);
    hist_lep_pt->Fill(lep_pt[index[0]]/1000.);
    hist_lep_trackisolation->Fill(lep_ptcone30[index[0]]/lep_pt[index[0]]);
    hist_lep_calorisolation->Fill(lep_etcone20[index[0]]/lep_pt[index[0]]);
    hist_lep_eta->Fill(lep_phi[index[0]]);
    hist_lep_e->Fill(lep_E[index[0]]/1000.);
    hist_nele->Fill(n_el);
    hist_nmuon->Fill(n_mu);
    hist_jetn->Fill(jet_n);
    hist_njets->Fill(n_jets);
    for (unsigned i = 0; i < n_jets; i++) {
      hist_jet_pt->Fill(jet_pt[i]/1000.);
      hist_jet_eta->Fill(jet_eta[i]);
      hist_jet_jvf->Fill(jet_jvf[i]);
      hist_jet_mv1->Fill(jet_MV1[i]);
    }
    hist_nbjets->Fill(n_bjets);
    hist_met->Fill(MET/1000.);
    hist_mtw->Fill(mTW/1000.);


    index.clear();

  }

  hist_nlep->Write();
  hist_lep_pt->Write();
  hist_lep_trackisolation->Write();
  hist_lep_calorisolation->Write();
  hist_lep_eta->Write();
  hist_lep_e->Write();
  hist_nele->Write();
  hist_nmuon->Write();
  hist_jetn->Write();
  hist_njets->Write();
  hist_jet_pt->Write();
  hist_jet_eta->Write();
  hist_jet_jvf->Write();
  hist_jet_mv1->Write();
  hist_nbjets->Write();
  hist_met->Write();
  hist_mtw->Write();



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
  cutflow->Draw();
  canvas->Draw();

  /*hist_njets_data->Draw();
  canvas->Draw();*/

}
