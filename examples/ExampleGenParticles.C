/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
# include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void FillHistos(std::string sample)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  TString tsample(sample);
  //  if (sample =="XXQQQQ_m50_ctau10000mm")  chain.Add(TString::Format("/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/Delphes/XXQQQQ_m50_ctau10000mm.root")); //livia
  if (sample =="XXQQQQ_m50_ctau10000mm")  chain.Add(TString::Format("/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/Delphes/XXQQQQ_m50_ctau10000mm_190520.root")); //matthew
  if (sample =="XXQQQQ_m50_ctau100mm")  chain.Add(TString::Format("/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/Delphes/XXQQQQ_m50_ctau100mm.root")); //matthew
  if (sample =="XXQQQQ_m50_ctau1mm")  chain.Add(TString::Format("/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/Delphes/XXQQQQ_m50_ctau1mm.root")); //matthew
  if(sample=="test")chain.Add("/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/Delphes/XXQQQQ_m50_ctau100mm/ggF-H-S1S2-Sdecay-qqnunu_UFO_ctau-10p0_412439_37093606_GENSIM_chunk0.root");
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  
  TClonesArray *branchJet = treeReader->UseBranch("JetPUPPI");
  TClonesArray *branchTracks = treeReader->UseBranch("EFlowTrack"); 

  TClonesArray *branchVtx = treeReader->UseBranch("Vertex");



  //output tree and histos
  const Float_t nsConv = 1E9; 
  const Float_t lightspeed = 2.998E2;
  
  TFile* fout;
  if (sample =="XXQQQQ_m50_ctau10000mm")fout= new TFile(TString::Format("fout_XXQQQQ_m50_ctau10000mm.root"), "RECREATE");
  if (sample =="XXQQQQ_m50_ctau100mm")fout= new TFile(TString::Format("fout_XXQQQQ_m50_ctau100mm.root"), "RECREATE");
  if (sample =="XXQQQQ_m50_ctau1mm")fout= new TFile(TString::Format("fout_XXQQQQ_m50_ctau1mm.root"), "RECREATE");
  if(sample=="test")fout= new TFile("test.root", "RECREATE");

  TTree t1("t1","a simple Tree with simple variables");

  Int_t nJetsGen_;
  t1.Branch("nJetsGen_",&nJetsGen_,"nJetsGen_/I");
  int nDisplJetsGen_;
  t1.Branch("nDisplJetsGen_",&nDisplJetsGen_,"nDisplJetsGen_/I");
  int nDisplJetsGen0_;
  t1.Branch("nDisplJetsGen0_",&nDisplJetsGen0_,"nDisplJetsGen0_/I");
  const Int_t nJetsGenMax=20;
  Float_t jetGenPT_[nJetsGenMax];
  t1.Branch("jetGenPT_", jetGenPT_, "jetGenPT_[nJetsGen_]/F");
  Float_t jetGenEta_[nJetsGenMax];
  t1.Branch("jetGenEta_", jetGenEta_, "jetGenEta_[nJetsGen_]/F");
  Float_t jetGenPhi_[nJetsGenMax];
  t1.Branch("jetGenPhi_", jetGenPhi_, "jetGenPhi_[nJetsGen_]/F");
  Float_t jetGenMass_[nJetsGenMax];
  t1.Branch("jetGenMass_", jetGenMass_, "jetGenMass_[nJetsGen_]/F");
  Float_t jetGenT_[nJetsGenMax];
  t1.Branch("jetGenT_", jetGenT_, "jetGenT_[nJetsGen_]/F");
  Float_t jetGenT0_[nJetsGenMax];
  t1.Branch("jetGenT0_", jetGenT0_, "jetGenT0_[nJetsGen_]/F");
  Float_t jetGenTuncorr_[nJetsGenMax];
  t1.Branch("jetGenTuncorr_", jetGenTuncorr_, "jetGenTuncorr_[nJetsGen_]/F");
  int nTracksInJetGen_[nJetsGenMax];
  t1.Branch("nTracksInJetGen_",&nTracksInJetGen_,"nTracksInJetGen_[nJetsGen_]/I");
  Int_t nTrk_;
  t1.Branch("nTrk_",&nTrk_,"nTrk_/I");
  Float_t trkT_[2000];
  t1.Branch("trkT_", trkT_, "trkT_[nTrk_]/F");
  Float_t trkEta_[2000];
  t1.Branch("trkEta_", trkEta_, "trkEta_[nTrk_]/F");
  Float_t trkPhi_[2000];
  t1.Branch("trkPhi_", trkPhi_, "trkPhi_[nTrk_]/F");
  Float_t trkT0_[2000];
  t1.Branch("trkT0_", trkT0_, "trkT0_[nTrk_]/F");
  Float_t trkDMTDToPV_[2000];
  t1.Branch("trkDMTDToPV_", trkDMTDToPV_, "trkDMTDToPV_[nTrk_]/F");
  Float_t trkTuncorr_[2000];
  t1.Branch("trkTuncorr_", trkTuncorr_, "trkTuncorr_[nTrk_]/F");
  



  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    std::cout<<"-------------------------------------------------------"<<std::endl;    
    if(entry % 100 == 0)std::cout<<entry<< " of " <<numberOfEntries<<std::endl;





    //save privary vertex
    Vertex  *vtx = (Vertex*) branchVtx->At(0);
    double    PVX = vtx->X;
    double    PVY = vtx->Y;
    double    PVZ = vtx->Z;


    //consider gen particles
   //use gen jet
    int nJetsGen=0;
    int nDisplJetsGen=0;
    int nDisplJetsGen0=0;
    Float_t jetGenPT[20];
    Float_t jetGenEta[20];
    Float_t jetGenPhi[20];
    Float_t jetGenMass[20];
    Float_t jetGenT[20];
    Float_t jetGenT0[20];
    Float_t jetGenTuncorr[20];
    Int_t nTracksInJetGen[20];
    Float_t trkDMTDToPV[2000];
    Float_t trkT[2000];
    Float_t trkEta[2000];
    Float_t trkPhi[2000];
    Float_t trkT0[2000];
    Float_t trkTuncorr[2000];
     int nTrk=0;
  


   GenParticle *particle, *mother;
    for(int i = 0; i < branchGenParticle->GetEntriesFast(); ++i)
      {
	particle = (GenParticle*) branchGenParticle->At(i);
	int motherPID=-999;
	float motherEta=999.;
	float motherPhi = 999.;
	float motherX=999.;
	float motherY=999.;
	float motherZ=999.;
	if(particle->M1>0){
	  mother = (GenParticle *) branchGenParticle->At(particle->M1);
	  motherPID= mother->PID;
	  motherEta=mother->Eta;
	  motherPhi=mother->Phi;
	  motherX=mother->X;
	  motherY=mother->Y;
	  motherZ=mother->Z;

      }
	int pdgid = particle->PID;
	int status =particle->Status;
	int charge = particle->Charge;
	float eta = particle->Eta;
	float phi = particle->Phi;
	float pt = particle->PT;
	TVector3* gp_vec = new TVector3();
	gp_vec->SetPtEtaPhi(pt,eta,phi);

	if(motherPID>900000&&status==23){
	  //std::cout<<" motherPID: "<<motherPID<< " X: "<<motherX<<" Y: "<<motherY<<" Z: "<<motherZ<<" eta: "<<motherEta<< " phi: "<<motherPhi<< std::endl;  
	  //	  std::cout<<"pdgid: "<<pdgid<<" X: "<<particle->X<<" Y: "<<particle->Y<<" Z: "<<particle->Z<<" eta: "<<particle->Eta<< " phi: "<<particle->Phi<< std::endl;
	  std::cout<<" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "<< std::endl;
	  double jetGenTime =0;                                                                                                                                                                     
	  double jetGenTimeCorrected =0;                                                                                                                                                            
	  double jetGenTimeCorrected0 =0;                                                                                                                                                            
	  int nTracksInGenJet=0;       

	  for(Int_t t = 0; t< branchTracks->GetEntriesFast(); t++){
	    
	    Track *tt = (Track *) branchTracks->At(t);
	    
	    if(tt->PT<1) continue;
	    
	    //track distance to MTD (assuming origin at 0,0,0)                                                                                                                                              
	    float distanceMTDToPV= TMath::Sqrt(pow((tt->XMTD-PVX),2)+pow((tt->YMTD-PVY),2)+pow((tt->ZMTD-PVZ),2));
	    float distanceMTDTo0= TMath::Sqrt(pow((tt->XMTD-0),2)+pow((tt->YMTD-0),2)+pow((tt->ZMTD-0),2));
	    float displacement = distanceMTDToPV;
	    float timeCorrection = distanceMTDToPV/lightspeed;
	    float timeCorrection0 = distanceMTDTo0/lightspeed;
	
	    TVector3* tt_vec = new TVector3();
	    tt_vec->SetPtEtaPhi(tt->PT,tt->Eta,tt->Phi);

	    double deltaR = gp_vec->DeltaR(*tt_vec);
	    if(deltaR>0.2)continue;
	    nTracksInGenJet++;

	    if(abs(tt->Eta)<1.4442)std::cout<<" tmtd: "<<tt->TMTD<< " x: "<<tt->XMTD<< " y: "<<tt->YMTD<<" z: "<<tt->ZMTD<<std::endl;
	    
	    trkDMTDToPV[nTrk] =distanceMTDToPV;
	    trkEta[nTrk]= tt->Eta;
	    trkPhi[nTrk]= tt->Phi;
	    trkT[nTrk]= tt->TMTD-timeCorrection;
	    trkT0[nTrk]= tt->TMTD-timeCorrection0;
	    trkTuncorr[nTrk]= tt->TMTD;
	    //	    std::cout<<trkTuncorr[nTrk]<<" "<<" MTDTOPV: "<<trkDMTDToPV[nTrk]<< " timecorr: "<<timeCorrection<<std::endl;
	    nTrk++;
	    //std::cout<<tt->T<<" "<<tt->TMTD<<" "<<tt->TECAL<<std::endl;
	    jetGenTime+=tt->TMTD;
	    jetGenTimeCorrected+=(tt->TMTD-timeCorrection);
	    jetGenTimeCorrected0+=(tt->TMTD-timeCorrection0);
	    
      }
	  jetGenTime/=nTracksInGenJet;                                                                                                                                                              
	  jetGenTimeCorrected/=nTracksInGenJet;                                                                                                                                                     
	  jetGenTimeCorrected0/=nTracksInGenJet;                                                                                                                                                     
	  jetGenTuncorr[nJetsGen]=jetGenTime;                                                                                                                                                       
	  jetGenT[nJetsGen]=jetGenTimeCorrected;                                                                                                                                                    
	  jetGenT0[nJetsGen]=jetGenTimeCorrected0;                                                                                                                                                    
	  nTracksInJetGen[nJetsGen]=nTracksInGenJet;                                                                                                                                                	  //      std::cout<<"----> "<<jetGenTime<<" "<<jetGenTimeCorrected<<std::endl;                                                                                                             
	  nJetsGen++;                                                                                                                                                                               
	  if(jetGenTimeCorrected>0.2) nDisplJetsGen++;     
	  if(jetGenTimeCorrected0>0.2) nDisplJetsGen0++;     
	}
      }


    /*    for(Int_t j = 0; j< branchGenJet->GetEntriesFast() ; j++) {

      Jet *gj = (Jet *) branchGenJet->At(j);
      //      std::cout<<jj<<" "<<jj->PT<<" "<<jj->X<<" "<<jj->XOuter<<" "<<jj->Y<<" "<<jj->YOuter<<" "<<jj->Z<<" "<<jj->ZOuter<<std::endl;                                                           
      if(gj->PT<10 || abs(gj->Eta)>3.0) continue;


      jetGenPT[nJetsGen]=gj->PT;
      jetGenEta[nJetsGen]=gj->Eta;
      jetGenPhi[nJetsGen]=gj->Phi;
      jetGenMass[nJetsGen]=gj->Mass;

      TVector3* gj_vec = new TVector3();
      gj_vec->SetPtEtaPhi(gj->PT,gj->Eta,gj->Phi);


      double jetGenTime =0;
      double jetGenTimeCorrected =0;
      int nTracksInGenJet=0;
      

      std::cout<<"========================"<<std::endl;                                                                                                                                         
      for(Int_t t = 0; t< branchTracks->GetEntriesFast(); t++){

        Track *tt = (Track *) branchTracks->At(t);

        if(tt->PT<1) continue;

        //track distance to MTD (assuming origin at 0,0,0)                                                                                                                                              
        float distanceToMTD = TMath::Sqrt(tt->XMTD*tt->XMTD+tt->YMTD*tt->YMTD+tt->ZMTD*tt->ZMTD);
	float distanceTo000= TMath::Sqrt(tt->X*tt->X+tt->Y*tt->Y+tt->Z*tt->Z);
	float displacement = distanceToMTD-distanceTo000;
        float timeCorrection = distanceToMTD/lightspeed;



        h_trkTime->Fill(tt->TMTD*nsConv-timeCorrection);
        TVector3* tt_vec = new TVector3();
        tt_vec->SetPtEtaPhi(tt->PT,tt->Eta,tt->Phi);

        double deltaR = gj_vec->DeltaR(*tt_vec);
        if(deltaR>0.5)continue;
        nTracksInGenJet++;
        //      std::cout<<"trk time: "<<tt->TOuter<<std::endl;                                                                                                                                         
	trkDToMTD[nTrk] = distanceToMTD;
	trkDTo000[nTrk] = distanceTo000;
	trkT[nTrk]= tt->TMTD*nsConv-timeCorrection;
	trkTuncorr[nTrk]= tt->TMTD*nsConv;
	//	std::cout<<trkT[nTrk]<<" "<<trkTuncorr[nTrk]<<std::endl;
	nTrk++;
	//std::cout<<tt->T<<" "<<tt->TMTD<<" "<<tt->TECAL<<std::endl;
	jetGenTime+=tt->TMTD*nsConv;
        jetGenTimeCorrected+=(tt->TMTD*nsConv-timeCorrection);

      }

      jetGenTime/=nTracksInGenJet;
      jetGenTimeCorrected/=nTracksInGenJet;
      jetGenTuncorr[nJetsGen]=jetGenTime;
      jetGenT[nJetsGen]=jetGenTimeCorrected;
      nTracksInJetGen[nJetsGen]=nTracksInGenJet;
        
      //      std::cout<<"----> "<<jetGenTime<<" "<<jetGenTimeCorrected<<std::endl;
      nJetsGen++;
      if(jetGenTimeCorrected>0.2) nDisplJetsGen++;
    }
    */

    //fill tree                                                                                                                                                                                         
    nJetsGen_=nJetsGen;
    nDisplJetsGen_=nDisplJetsGen;
    nDisplJetsGen0_=nDisplJetsGen0;
    for(int i =0;i<nJetsGen;i++){
      /*  jetGenPT_[i]=jetGenPT[i];
      jetGenEta_[i]=jetGenEta[i];
      jetGenPhi_[i]=jetGenPhi[i];
      jetGenMass_[i]=jetGenMass[i];
      */
      jetGenT_[i]=jetGenT[i];
      jetGenT0_[i]=jetGenT0[i];
      jetGenTuncorr_[i]=jetGenTuncorr[i];
      nTracksInJetGen_[i]=nTracksInJetGen[i];
    }
    nTrk_=nTrk;
    for(int i =0;i<nTrk;i++){
      trkDMTDToPV_[i]=trkDMTDToPV[i];
      trkEta_[i]=trkEta[i];
      trkPhi_[i]=trkPhi[i];
      trkT_[i]=trkT[i];
      trkT0_[i]=trkT0[i];
      trkTuncorr_[i]=trkTuncorr[i];
    }

    t1.Fill();

  }
  fout->cd();
  
  t1.Write();
  fout->Write();
  fout->Close();
}




void FillAllHistos(){
  //FillHistos("XXQQQQ_m50_ctau10000mm");
  //  FillHistos("XXQQQQ_m50_ctau100mm");
  //FillHistos("XXQQQQ_m50_ctau1mm");
  FillHistos("test");
}



void MakeAllPlots(){

  TFile* f1 = TFile::Open("fout_XXQQQQ_m50_ctau100mm.root", "READ");
  TTree* t1 = (TTree*) f1->Get("t1");
  TFile* f2 = TFile::Open("fout_XXQQQQ_m50_ctau100mm.root", "READ");
  TTree* t2 = (TTree*) f2->Get("t1");
  TFile* f3 = TFile::Open("fout_XXQQQQ_m50_ctau100mm.root", "READ");
  TTree* t3 = (TTree*) f3->Get("t1");


  TH1F* h1_t = new TH1F("h1_t","", 160, -0.5, 5); ;
  t1->Draw("jetGenT_>>h1_t");
  h1_t->SetLineColor(kAzure+1);
  h1_t->GetXaxis()->SetTitle("Jet Time [ns]");
  TH1F* h2_t = new TH1F("h2_t","", 160, -0.5, 5); ;
  t2->Draw("jetGenT_>>h2_t");
  h2_t->SetLineColor(kSpring);
  TH1F* h3_t = new TH1F("h3_t","", 160, -0.5, 5); ;
  t3->Draw("jetGenT_>>h3_t");
  h3_t->SetLineColor(kOrange+8);


  TH1F* h1_tuncorr = new TH1F("h1_tuncorr","", 100, 0., 25.); ;
  t1->Draw("jetGenTuncorr_>>h1_tuncorr");
  h1_tuncorr->SetLineColor(kAzure+1);
  h1_tuncorr->GetXaxis()->SetTitle("Jet Time [ns]");
  TH1F* h2_tuncorr = new TH1F("h2_tuncorr","", 100, 0., 25.); ;
  t2->Draw("jetGenTuncorr_>>h2_tuncorr");
  h2_tuncorr->SetLineColor(kSpring);
  TH1F* h3_tuncorr = new TH1F("h3_tuncorr","", 100, 0., 25.); ;
  t3->Draw("jetGenTuncorr_>>h3_tuncorr");
  h3_tuncorr->SetLineColor(kOrange+8);



  TH1F* h1_tt = new TH1F("h1_tt","", 160, -0.5, 5); ;
  t1->Draw("trkT_>>h1_tt");
  h1_tt->SetLineColor(kAzure+1);
  h1_tt->GetXaxis()->SetTitle("Tracks Time [ns]");
  TH1F* h2_tt = new TH1F("h2_tt","", 160, -0.5, 5); ;
  t2->Draw("trkT_>>h2_tt");
  h2_tt->SetLineColor(kSpring);
  TH1F* h3_tt = new TH1F("h3_tt","", 160, -0.5, 5); ;
  t3->Draw("trkT_>>h3_tt");
  h3_tt->SetLineColor(kOrange+8);

  TH1F* h1_ttunc = new TH1F("h1_ttunc","", 160, 0., 25.); ;
  t1->Draw("trkTuncorr_>>h1_ttunc");
  h1_ttunc->SetLineColor(kAzure+1);
  h1_ttunc->GetXaxis()->SetTitle("Tracks Time [ns]");
  TH1F* h2_ttunc = new TH1F("h2_ttunc","", 160, 0., 25.); ;
  t2->Draw("trkTuncorr_>>h2_ttunc");
  h2_ttunc->SetLineColor(kSpring);
  TH1F* h3_ttunc = new TH1F("h3_ttunc","", 160, 0., 25.); ;
  t3->Draw("trkTuncorr_>>h3_ttunc");
  h3_ttunc->SetLineColor(kOrange+8);



  TH1F* h1_ntrk = new TH1F("h1_ntrk","", 50, 0, 50); ;
  t1->Draw("nTracksInJetGen_>>h1_ntrk");
  h1_ntrk->SetLineColor(kAzure+1);
  h1_ntrk->GetXaxis()->SetTitle("# tracks in Jet");
  TH1F* h2_ntrk = new TH1F("h2_ntrk","", 50, 0, 50); ;
  t2->Draw("nTracksInJetGen_>>h2_ntrk");
  h2_ntrk->SetLineColor(kSpring);
  TH1F* h3_ntrk = new TH1F("h3_ntrk","", 50, 0, 50); ;
  t3->Draw("nTracksInJetGen_>>h3_ntrk");
  h3_ntrk->SetLineColor(kOrange+8);




  TH1F* h1_njet = new TH1F("h1_njet","", 10, 0, 10); ;
  t1->Draw("nJetsGen_>>h1_njet");
  h1_njet->SetLineColor(kAzure+1);
  h1_njet->GetXaxis()->SetTitle("# Gen Jets");
  TH1F* h2_njet = new TH1F("h2_njet","", 10, 0, 10); ;
  t2->Draw("nJetsGen_>>h2_njet");
  h2_njet->SetLineColor(kSpring);
  TH1F* h3_njet = new TH1F("h3_njet","", 10, 0, 10); ;
  t3->Draw("nJetsGen_>>h3_njet");
  h3_njet->SetLineColor(kOrange+8);

  TH1F* h1_ndispljet = new TH1F("h1_ndispljet","", 10, 0, 10); ;
  t1->Draw("nDisplJetsGen_>>h1_ndispljet");
  h1_ndispljet->SetLineColor(kAzure+1);
  h1_ndispljet->GetXaxis()->SetTitle("# Displaced Gen Jets");
  TH1F* h2_ndispljet = new TH1F("h2_ndispljet","", 10, 0, 10); ;
  t2->Draw("nDisplJetsGen_>>h2_ndispljet");
  h2_ndispljet->SetLineColor(kSpring);
  TH1F* h3_ndispljet = new TH1F("h3_ndispljet","", 10, 0, 10); ;
  t3->Draw("nDisplJetsGen_>>h3_ndispljet");
  h3_ndispljet->SetLineColor(kOrange+8);






  TLegend* leg = new TLegend(0.6, 0.75, 0.8, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  //  leg->AddEntry(h1_t, "c#tau = 10000 mm", "L");
  //leg->AddEntry(h2_t, "c#tau = 100 mm", "L");
  leg->AddEntry(h3_t, "c#tau = 100 mm", "L");


  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  h1_t->DrawNormalized();
  leg->Draw("same");
  h2_t->DrawNormalized("same");
  h3_t->DrawNormalized("same");
  c->SaveAs("~/www/DisplacedJets/time.png");
  c->SaveAs("~/www/DisplacedJets/time.pdf");
  c->SetLogy();
  c->SaveAs("~/www/DisplacedJets/time_log.png");
  c->SaveAs("~/www/DisplacedJets/time_log.pdf");


  h1_tuncorr->DrawNormalized();
  leg->Draw("same");
  h2_tuncorr->DrawNormalized("same");
  h3_tuncorr->DrawNormalized("same");
  c->SaveAs("~/www/DisplacedJets/timeuncorr.png");
  c->SaveAs("~/www/DisplacedJets/timeuncorr.pdf");
  c->SetLogy();
  c->SaveAs("~/www/DisplacedJets/timeuncorr_log.png");
  c->SaveAs("~/www/DisplacedJets/timeuncorr_log.pdf");


  h1_tt->DrawNormalized();
  leg->Draw("same");
  h2_tt->DrawNormalized("same");
  h3_tt->DrawNormalized("same");
  c->SaveAs("~/www/DisplacedJets/time_tt.png");
  c->SaveAs("~/www/DisplacedJets/time_tt.pdf");
  c->SetLogy();
  c->SaveAs("~/www/DisplacedJets/time_tt_log.png");
  c->SaveAs("~/www/DisplacedJets/time_tt_log.pdf");

 
  h1_ttunc->DrawNormalized();
  leg->Draw("same");
  h2_ttunc->DrawNormalized("same");
  h3_ttunc->DrawNormalized("same");
  c->SaveAs("~/www/DisplacedJets/time_ttunc.png");
  c->SaveAs("~/www/DisplacedJets/time_ttunc.pdf");
  c->SetLogy();
  c->SaveAs("~/www/DisplacedJets/time_ttunc_log.png");
  c->SaveAs("~/www/DisplacedJets/time_ttunc_log.pdf");



  h1_ntrk->DrawNormalized();
  leg->Draw("same");
  h2_ntrk->DrawNormalized("same");
  h3_ntrk->DrawNormalized("same");
  c->SaveAs("~/www/DisplacedJets/ntrk.png");
  c->SaveAs("~/www/DisplacedJets/ntrk.pdf");



  h1_njet->DrawNormalized();
  leg->Draw("same");
  h2_njet->DrawNormalized("same");
  h3_njet->DrawNormalized("same");
  c->SaveAs("~/www/DisplacedJets/njet.png");
  c->SaveAs("~/www/DisplacedJets/njet.pdf");


  h1_ndispljet->DrawNormalized();
  leg->Draw("same");
  h2_ndispljet->DrawNormalized("same");
  h3_ndispljet->DrawNormalized("same");
  c->SaveAs("~/www/DisplacedJets/ndispljet.png");
  c->SaveAs("~/www/DisplacedJets/ndispljet.pdf");





}
