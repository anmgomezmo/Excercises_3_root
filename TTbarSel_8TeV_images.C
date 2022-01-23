#include <iostream>
#include <string>
#include <stdio.h>
#include <vector>

void TTbarSel_8TeV_images(){ //void function with the same name as script or file code

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

  TH1F *hist_lep_pt = new TH1F("Leptons pt","Lepton pt; pt (GeV); Events",100,0,150);
  TH1F *hist_lep_trackisolation = new TH1F("Leptons track isolation","Lepton track; ptcone/pt; Events",30,0,1);
  TH1F *hist_lep_calorisolation = new TH1F("Leptons calorimeter isolation","Lepton calorimeter; etcone/pt; Events",30,0,1);
  TH1F *hist_lep_eta = new TH1F("Leptons eta","Lepton eta; eta; Events",40,-10,10);
  TH1F *hist_jet_number = new TH1F("Number of jets","Number of jets; Jet multiplicity; Events",10,0,10);
  TH1F *hist_jet_pt = new TH1F("Jet pt","Jet pt; pt (Gev); Events",100,0,150);
  TH1F *hist_jet_eta = new TH1F("Jet eta","Jet eta; eta; Events",40,-10,10);
  TH1F *hist_jet_jvf = new TH1F("Jet jvf","Jet jvf; jvf; Events",20,-2,3);
  TH1F *hist_jet_mv1 = new TH1F("MV1","MV1; mv1; Events",20,0,1);
  TH1F *hist_bjet_number = new TH1F("Number of bjets","Number of bjets; Jet multiplicity; Events",6,0,6);
  TH1F *hist_met = new TH1F("MET","MET; met (GeV); Events",100,0,150);
  TH1F *hist_mtw = new TH1F("mTW","mTW; mTW (GeV); Events",100,0,200);


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

    //First cut: Good vertex
    if (!good_vtx) continue;
    cut1++;
    cutflow->Fill(1);

    //Second cut: Trigger
    if (!e_trig && !mu_trig) continue;
    cut2++;
    cutflow->Fill(2);


    //Loop over leptons
    for (unsigned i = 0; i < lep_n; i++) {
      hist_lep_pt->Fill(lep_pt[i]/1000.);
      hist_lep_trackisolation->Fill(lep_ptcone30[i]/lep_pt[i]);
      hist_lep_calorisolation->Fill(lep_etcone20[i]/lep_pt[i]);
      hist_lep_eta->Fill(lep_eta[i]);


      //TLorentzVector definitions
      TLorentzVector Lepton = TLorentzVector();
      TLorentzVector MeT = TLorentzVector();

      //To complete: Lorentz vectors for the lepton and MET. Use SetPtEtaPhiE().
      Lepton.SetPtEtaPhiE(lep_pt[i],lep_eta[i],lep_phi[i],lep_E[i]);
      MeT.SetPtEtaPhiE(MET,0,MET_phi,MET);

      //Calculation of the mTW using TLorenzt vectors
      float mTW = sqrt(2*Lepton.Pt()*MeT.E()*(1-cos(Lepton.DeltaPhi(MeT))));
      hist_mtw->Fill(mTW/1000.);
    }

    int n_bjet = 0;
    hist_jet_number->Fill(jet_n);
    for (unsigned j = 0; j < jet_n; j++) {
        hist_jet_pt->Fill(jet_pt[j]/1000.);
        hist_jet_eta->Fill(jet_eta[j]);
        hist_jet_jvf->Fill(jet_jvf[j]);
        hist_jet_mv1->Fill(jet_MV1[j]);
        if (jet_MV1[j] > 0.7892){
            n_bjet++;
        }
    }
    hist_bjet_number->Fill(n_bjet);
    hist_met->Fill(MET/1000.);



  }



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

  hist_lep_pt->SetMarkerStyle(20);
  hist_lep_pt->SetMarkerColor(kRed);
  hist_lep_pt->SetStats(0);
  hist_lep_pt->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/1_lep_pt.png");
  canvas->Print("Data/Data.pdf(","Title:Leptons pt");

  hist_lep_trackisolation->SetMarkerStyle(20);
  hist_lep_trackisolation->SetMarkerColor(kRed);
  hist_lep_trackisolation->SetStats(0);
  hist_lep_trackisolation->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/2_lep_track.png");
  canvas->Print("Data/Data.pdf","Title:Leptons track");

  hist_lep_calorisolation->SetMarkerStyle(20);
  hist_lep_calorisolation->SetMarkerColor(kRed);
  hist_lep_calorisolation->SetStats(0);
  hist_lep_calorisolation->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/3_lep_calo.png");
  canvas->Print("Data/Data.pdf","Title:Leptons calorimeter");

  hist_lep_eta->SetMarkerStyle(20);
  hist_lep_eta->SetMarkerColor(kRed);
  hist_lep_eta->SetStats(0);
  hist_lep_eta->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/4_lep_eta.png");
  canvas->Print("Data/Data.pdf","Title:Leptons eta");

  hist_jet_number->SetMarkerStyle(20);
  hist_jet_number->SetMarkerColor(kRed);
  hist_jet_number->SetStats(0);
  hist_jet_number->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/5_jet_number.png");
  canvas->Print("Data/Data.pdf","Title:Number of jets");

  hist_jet_pt->SetMarkerStyle(20);
  hist_jet_pt->SetMarkerColor(kRed);
  hist_jet_pt->SetStats(0);
  hist_jet_pt->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/6_jet_pt.png");
  canvas->Print("Data/Data.pdf","Title:Jets pt");

  hist_jet_eta->SetMarkerStyle(20);
  hist_jet_eta->SetMarkerColor(kRed);
  hist_jet_eta->SetStats(0);
  hist_jet_eta->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/7_jet_eta.png");
  canvas->Print("Data/Data.pdf","Title:Jets eta");

  hist_jet_jvf->SetMarkerStyle(20);
  hist_jet_jvf->SetMarkerColor(kRed);
  hist_jet_jvf->SetStats(0);
  hist_jet_jvf->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/8_jet_jvf.png");
  canvas->Print("Data/Data.pdf","Title:Jets jvf");

  hist_jet_mv1->SetMarkerStyle(20);
  hist_jet_mv1->SetMarkerColor(kRed);
  hist_jet_mv1->SetStats(0);
  hist_jet_mv1->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/9_jet_mv1.png");
  canvas->Print("Data/Data.pdf","Title:Jets mv1");

  hist_bjet_number->SetMarkerStyle(20);
  hist_bjet_number->SetMarkerColor(kRed);
  hist_bjet_number->SetStats(0);
  hist_bjet_number->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/10_bjet_number.png");
  canvas->Print("Data/Data.pdf","Title:Number of b-Jets");

  hist_met->SetMarkerStyle(20);
  hist_met->SetMarkerColor(kRed);
  hist_met->SetStats(0);
  hist_met->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/11_met.png");
  canvas->Print("Data/Data.pdf","Title:MET");

  hist_mtw->SetMarkerStyle(20);
  hist_mtw->SetMarkerColor(kRed);
  hist_mtw->SetStats(0);
  hist_mtw->Draw("p");
  canvas->Draw();
  canvas->SaveAs("Data/12_mtw.png");
  canvas->Print("Data/Data.pdf)","Title:mTW");


}
