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
  if (sample =="XXQQQQ_m50_ctau10000mm")  chain.Add(TString::Format("/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/Delphes/XXQQQQ_m50_ctau10000mm_NewDelphes_160519.root")); //matthew
  if (sample =="XXQQQQ_m50_ctau100mm")  chain.Add(TString::Format("/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/Delphes/XXQQQQ_m50_ctau100mm.root")); //matthew
  if (sample =="XXQQQQ_m50_ctau1mm")  chain.Add(TString::Format("/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/Delphes/XXQQQQ_m50_ctau1mm.root")); //matthew

  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
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
  TTree t1("t1","a simple Tree with simple variables");

  Int_t nJetsGen_;
  t1.Branch("nJetsGen_",&nJetsGen_,"nJetsGen_/I");
  int nDisplJetsGen_;
  t1.Branch("nDisplJetsGen_",&nDisplJetsGen_,"nDisplJetsGen_/I");
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
  Float_t jetGenTuncorr_[nJetsGenMax];
  t1.Branch("jetGenTuncorr_", jetGenTuncorr_, "jetGenTuncorr_[nJetsGen_]/F");
  int nTracksInJetGen_[nJetsGenMax];
  t1.Branch("nTracksInJetGen_",&nTracksInJetGen_,"nTracksInJetGen_[nJetsGen_]/I");
  Int_t nTrk_;
  t1.Branch("nTrk_",&nTrk_,"nTrk_/I");
  Float_t trkT_[2000];
  t1.Branch("trkT_", trkT_, "trkT_[nTrk_]/F");
  Float_t trkTuncorr_[2000];
  t1.Branch("trkTuncorr_", trkTuncorr_, "trkTuncorr_[nTrk_]/F");
  


  Int_t nJets_;
  t1.Branch("nJets_",&nJets_,"nJets_/I");
  int nGenJets_;
  t1.Branch("nGenJets_",&nGenJets_,"nGenJets_/I");
  int nDisplJets_;
  t1.Branch("nDisplJets_",&nDisplJets_,"nDisplJets_/I");
  const Int_t nJetsMax=20;
  int nTracksInJet_[nJetsMax];
  t1.Branch("nTracksInJet_",&nTracksInJet_,"nTracksInJet_[nJets_]/I");
  Float_t jetPT_[nJetsMax];
  t1.Branch("jetPT_", jetPT_, "jetPT_[nJets_]/F");
  Float_t jetEta_[nJetsMax];
  t1.Branch("jetEta_", jetEta_, "jetEta_[nJets_]/F");
  Float_t jetPhi_[nJetsMax];
  t1.Branch("jetPhi_", jetPhi_, "jetPhi_[nJets_]/F");
  Float_t jetMass_[nJetsMax];
  t1.Branch("jetMass_", jetMass_, "jetMass_[nJets_]/F");
  
  Int_t jetIsGenMatched_[nJetsMax];
  t1.Branch("jetIsGenMatched_", jetIsGenMatched_, "jetIsGenMatched_[nJets_]/I");
  Int_t jetGenMatchedIndex_[nJetsMax];
  t1.Branch("jetGenMatchedIndex_", jetGenMatchedIndex_, "jetGenMatchedIndex_[nJets_]/I");
  Float_t jetDRMatching_[nJetsMax];
  t1.Branch("jetDRMatching_", jetDRMatching_, "jetDRMatching_[nJets_]/F");

  
  Float_t jetT_[nJetsMax];
  t1.Branch("jetT_", jetT_, "jetT_[nJets_]/F");
  Float_t jetAlphaMax_[nJetsMax];
  t1.Branch("jetAlphaMax_", jetAlphaMax_, "jetAlphaMax_[nJets_]/F");
  Float_t jetIP2D_[nJetsMax];
  t1.Branch("jetIP2D_", jetIP2D_, "jetIP2D_[nJets_]/F");
  Float_t jetTheta2D_[nJetsMax];
  t1.Branch("jetTheta2D_", jetTheta2D_, "jetTheta2D_[nJets_]/F");

  TH1F* h_nGenJets = new TH1F("h_nGenJets", "", 10, 0, 10);
  h_nGenJets->Sumw2();

  TH1F* h_nJets = new TH1F("h_nJets", "", 10, 0, 10);
  h_nJets->Sumw2();

  TH1F* h_nTracksInJet = new TH1F("h_nTracksInJet", "", 30, 0, 30);
  h_nTracksInJet->Sumw2();

  TH1F* h_jetTime = new TH1F("h_jetTime", "", 1000, 0, 25);
  h_jetTime->Sumw2();

  TH1F* h_jetTimeCorrected = new TH1F("h_jetTimeCorrected", "", 200, -2, 2);
  h_jetTimeCorrected->Sumw2();

  TH1F* h_trkTime = new TH1F("h_trkTime", "", 1000, 0, 25);
  h_trkTime->Sumw2();

  TH1F* h_IP2D = new TH1F("h_IP2D", "", 50, -10., 5.);
  h_IP2D->Sumw2();

  TH1F* h_alpha = new TH1F("h_alpha", "", 20, 0.,1.);
  h_alpha->Sumw2();

  TH1F* h_theta2D = new TH1F("h_theta2D", "", 40,-3.5, 2.5);
  h_theta2D->Sumw2();


  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    std::cout<<"-------------------------------------------------------"<<std::endl;    
    if(entry % 100 == 0)std::cout<<entry<< " of " <<numberOfEntries<<std::endl;

    //use gen jet
   
    int nJetsGen=0;
    int nDisplJetsGen=0;
    Float_t jetGenPT[20];
    Float_t jetGenEta[20];
    Float_t jetGenPhi[20];
    Float_t jetGenMass[20];
    Float_t jetGenT[20];
    Float_t jetGenTuncorr[20];
    Int_t nTracksInJetGen[20];
    Float_t trkT[2000];
    Float_t trkTuncorr[2000];
    int nTrk=0;
    for(Int_t j = 0; j< branchGenJet->GetEntriesFast() ; j++) {

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
        float timeCorrection = distanceToMTD/lightspeed;



        h_trkTime->Fill(tt->TMTD*nsConv-timeCorrection);
        TVector3* tt_vec = new TVector3();
        tt_vec->SetPtEtaPhi(tt->PT,tt->Eta,tt->Phi);

        double deltaR = gj_vec->DeltaR(*tt_vec);
        if(deltaR>0.5)continue;
        nTracksInGenJet++;
        //      std::cout<<"trk time: "<<tt->TOuter<<std::endl;                                                                                                                                         
	trkT[nTrk]= tt->TMTD*nsConv-timeCorrection;
	trkTuncorr[nTrk]= tt->TMTD*nsConv;
	std::cout<<trkT[nTrk]<<" "<<trkTuncorr[nTrk]<<std::endl;
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
        
      std::cout<<"----> "<<jetGenTime<<" "<<jetGenTimeCorrected<<std::endl;
      nJetsGen++;
      if(jetGenTimeCorrected>0.2) nDisplJetsGen++;
    }


    //fill tree                                                                                                                                                                                         
    nJetsGen_=nJetsGen;
    nDisplJetsGen_=nDisplJetsGen;
    for(int i =0;i<nJetsGen;i++){
      jetGenPT_[i]=jetGenPT[i];
      jetGenEta_[i]=jetGenEta[i];
      jetGenPhi_[i]=jetGenPhi[i];
      jetGenMass_[i]=jetGenMass[i];
      jetGenT_[i]=jetGenT[i];
      jetGenTuncorr_[i]=jetGenTuncorr[i];
      nTracksInJetGen_[i]=nTracksInJetGen[i];
    }
    nTrk_=nTrk;
    for(int i =0;i<nTrk;i++){
      trkT_[i]=trkT[i];
      trkTuncorr_[i]=trkTuncorr[i];
    }


    int nJets=0;
    Float_t jetPT[20];
    Float_t jetEta[20];
    Float_t jetPhi[20];
    Float_t jetMass[20];
    Float_t jetT[20];
    Float_t jetAlphaMax[20];
    Float_t jetIP2D[20];
    Float_t jetTheta2D[20];
    Int_t jetIsGenMatched[20];
    Int_t jetGenMatchedIndex[20];
    Float_t jetDRMatching[20];
    Int_t nTracksInJet[20];
    int nGenJets=0;
    int nDisplJets=0;






    //count # jets pt >10 GeV in good acceptance
    for(Int_t j = 0; j< branchJet->GetEntriesFast() ; j++) {

      Jet *jj = (Jet *) branchJet->At(j);
      
      if(jj->PT<10 || abs(jj->Eta)>3.0) continue;
      
      jetPT[nJets]=jj->PT;
      jetEta[nJets]=jj->Eta;
      jetPhi[nJets]=jj->Phi;
      jetMass[nJets]=jj->Mass;

      TVector3* jj_vec = new TVector3();
      jj_vec->SetPtEtaPhi(jj->PT,jj->Eta,jj->Phi);


      //gen matching
     
      double deltaRMin=999;
      int genJetIndex=999;
      for(Int_t j = 0; j< branchGenJet->GetEntriesFast() ; j++) {
	
	Jet *gj = (Jet *) branchGenJet->At(j);
	//      std::cout<<jj<<" "<<jj->PT<<" "<<jj->X<<" "<<jj->XOuter<<" "<<jj->Y<<" "<<jj->YOuter<<" "<<jj->Z<<" "<<jj->ZOuter<<std::endl;
	if(gj->PT<10 || abs(gj->Eta)>3.0) continue;
	nGenJets++;
	TVector3* gj_vec = new TVector3();
        gj_vec->SetPtEtaPhi(gj->PT,gj->Eta,gj->Phi);

        double deltaR = jj_vec->DeltaR(*gj_vec);
        if(deltaR>deltaRMin)continue;
	
	deltaRMin=deltaR;
	genJetIndex=j;
      }
      //      std::cout<<deltaRMin<<" "<<genJetIndex<<std::endl;
      jetIsGenMatched[nJets]=0;
      jetGenMatchedIndex[nJets]=999;
      jetDRMatching[nJets]=deltaRMin;

      if(deltaRMin<0.4){
      jetIsGenMatched[nJets]=1;
      jetGenMatchedIndex[nJets]=genJetIndex;
      jetDRMatching[nJets]=deltaRMin;
      }
      //      std::cout<<jetIsGenMatched[nJets]<<" "<<jetGenMatchedIndex[nJets]<<std::endl;

      //count # tracks in jet
      double jetTime =0;
      double jetTimeCorrected =0;
      int nTracksInJ=0;
      double ptAllTracks=0;
      double IP2D=0;
      double theta2D=0;
      std::vector<double> theta2Ds;
      std::vector<double> IP2Ds;

      //      std::cout<<"========================"<<std::endl;
      for(Int_t t = 0; t< branchTracks->GetEntriesFast(); t++){
	
	Track *tt = (Track *) branchTracks->At(t);
	
	if(tt->PT<1) continue;
	
	//track distance to MTD (assuming origin at 0,0,0)
	float distanceToMTD = TMath::Sqrt(tt->XMTD*tt->XMTD+tt->YMTD*tt->YMTD+tt->ZMTD*tt->ZMTD);
	float timeCorrection = distanceToMTD/lightspeed;
	
	h_trkTime->Fill(tt->TMTD*nsConv-timeCorrection);
	TVector3* tt_vec = new TVector3();
	tt_vec->SetPtEtaPhi(tt->PT,tt->Eta,tt->Phi);

	double deltaR = jj_vec->DeltaR(*tt_vec);
	if(deltaR>0.5)continue;
	nTracksInJ++;
	//	std::cout<<"trk time: "<<tt->TOuter<<std::endl;
	jetTime+=tt->TMTD*nsConv;
	jetTimeCorrected+=(tt->TMTD*nsConv-timeCorrection);
	ptAllTracks+=tt->PT;
	//	std::cout<<tt->ErrorT<<" "<<tt->T<<std::endl;
	//compute theta2D
	TVector3 tt_innerpos;
	tt_innerpos.SetXYZ(tt->X, tt->Y, tt->Z);
	Vertex *vtx = (Vertex *) branchVtx->At(0);

	TVector3 primvtx_pos;
	primvtx_pos.SetXYZ(vtx->X, vtx->Y, vtx->Z);
	
	TVector3 deltaPos;
	deltaPos= tt_innerpos-primvtx_pos;

	double mag2DeltaPos = TMath::Sqrt((deltaPos.X()*deltaPos.X()) + (deltaPos.Y()*deltaPos.Y()));
	double mag2Mom = TMath::Sqrt((tt_vec->X()*tt_vec->X()) + (tt_vec->Y()*tt_vec->Y()));
	double theta2D = TMath::ACos((deltaPos.X()*tt_vec->X()+deltaPos.Y()*tt_vec->Y())/(mag2Mom*mag2DeltaPos));

	theta2Ds.push_back(theta2D);

	//compute IP2D
	TVector3 primvtx_posErr;
        primvtx_posErr.SetXYZ(vtx->ErrorX, vtx->ErrorY, vtx->ErrorZ);
	//	std::cout<<vtx->ErrorX<<" "<<vtx->ErrorY<<" "<<vtx->ErrorZ<<std::endl;
	TVector3 tt_innerposErr;
        tt_innerposErr.SetXYZ(tt->ErrorD0/2, tt->ErrorD0/2, tt->ErrorDZ);
	//	std::cout<<tt->ErrorD0/2<<" "<<tt->ErrorDZ<<std::endl;
	TVector3 deltaPosErr;
        deltaPosErr= tt_innerposErr+primvtx_posErr;
	double mag2DeltaPosErr = TMath::Sqrt((deltaPosErr.X()*deltaPosErr.X()) + (deltaPosErr.Y()*deltaPosErr.Y()));
	IP2D=mag2DeltaPos;
	IP2Ds.push_back(IP2D);

	//	std::cout<<IP2D<<" "<<mag2DeltaPos<<" "<<mag2DeltaPosErr<<std::endl;
      }


      std::sort(theta2Ds.begin(),theta2Ds.end());
      double medianTheta2D = 0;
      double log10medianTheta2D = 0;
      if (theta2Ds.size() > 0){
	medianTheta2D = theta2Ds[theta2Ds.size()/2];
	log10medianTheta2D = log10(medianTheta2D);
      }


      std::sort(IP2Ds.begin(),IP2Ds.end());
      double medianIP2D = 0;
      double log10medianIP2D = 0;
      if (IP2Ds.size() > 0){
	medianIP2D = IP2Ds[IP2Ds.size()/2];
	log10medianIP2D =log10(medianIP2D);
      }



      double ptAllPVTracks=0;
      double nTracksPVTemp=0;
      double ptPVTracksMax =0;
      int nTracksPV=0;
      double alphaMax=0;
      
      //compute alpha: need to sum pt of all trks from PV inside the jet cone
      for(Int_t v = 0; v < branchVtx->GetEntriesFast(); ++v)
	{
	  Vertex *vtx = (Vertex*) branchVtx->At(v);
	  double ptPVTracks = 0.;
	  int nTracksPVTemp = 0;
	  //	  std::cout<<"vtx: "<<v<<std::endl;
	  for(int c = 0; c < vtx->Constituents.GetEntriesFast(); ++c)
	    {
	     TObject *object = vtx->Constituents.At(c);
	     //	     std::cout<<"const: "<<c<<std::endl;
	      // Check if the constituent is accessible
	      if(object == 0) continue;

	      if(object->IsA() == Jet::Class())
		{
		  //		  std::cout<<"jet"<<std::endl;
		  Track *track = (Track*) object;
		  if (track->PT < 1)continue; 
		  TVector3 pvTrackVecTemp;
		  pvTrackVecTemp.SetPtEtaPhi(track->PT,track->Eta,track->Phi);
		  //If pv track associated with jet add pt to ptPVTracks
		  if (pvTrackVecTemp.DeltaR(*jj_vec) > 0.4)continue;
		      ptPVTracks += track->PT;
		      ptAllPVTracks += track->PT;
		      nTracksPVTemp ++;
		}//if is track
	    }//loop costituents
	  
	  if (ptPVTracks > ptPVTracksMax) {
	    ptPVTracksMax = ptPVTracks;
	    nTracksPV = nTracksPVTemp;
	  }//check max
       	}//loop vtx
      
      alphaMax = ptPVTracksMax/ptAllTracks;
      jetTime/=nTracksInJ;
      jetTimeCorrected/=nTracksInJ;


      h_jetTime->Fill(jetTime);
      h_jetTimeCorrected->Fill(jetTimeCorrected);
      h_nTracksInJet->Fill(nTracksInJ);
      h_theta2D->Fill(log10medianTheta2D);
      h_IP2D->Fill(log10medianIP2D);
      h_alpha->Fill(alphaMax);
    
      jetT[nJets]=jetTimeCorrected;
      //      jetErrorT[nJets]=jetTimeErr;
      jetAlphaMax[nJets]=alphaMax;
      jetIP2D[nJets]=log10medianIP2D;
      jetTheta2D[nJets]=log10medianTheta2D;
      nTracksInJet[nJets]=nTracksInJ;
      nJets++;

      if(jetTimeCorrected>0.2) nDisplJets++;
    }

    h_nGenJets->Fill(nGenJets);
    h_nJets->Fill(nJets);
    
    //fill tree
    nJets_=nJets;
    nGenJets_=nGenJets;
    nDisplJets_=nDisplJets;
    for(int i =0;i<nJets;i++){
      jetPT_[i]=jetPT[i];
      jetEta_[i]=jetEta[i];
      jetPhi_[i]=jetPhi[i];
      jetMass_[i]=jetMass[i];
      jetT_[i]=jetT[i];
      //      jetErrorT_[i]=jetErrorT[i];
      jetAlphaMax_[i]=jetAlphaMax[i];
      jetIP2D_[i]=jetIP2D[i];
      jetTheta2D_[i]=jetTheta2D[i];
      jetIsGenMatched_[i]=jetIsGenMatched[i];
      jetGenMatchedIndex_[i]=jetGenMatchedIndex[i];
      jetDRMatching_[i]=jetDRMatching[i];
    
      nTracksInJet_[i]=nTracksInJet[i];
}
    t1.Fill();

  }
  fout->cd();
  // Show resulting histograms
  h_nJets->GetXaxis()->SetTitle("# jets");
  h_jetTime->GetXaxis()->SetTitle("Jet Time [s]");
  h_jetTimeCorrected->GetXaxis()->SetTitle("Jet Time Corrected [s]");
  h_nTracksInJet->GetXaxis()->SetTitle("# tracks in jet");
  h_theta2D->GetXaxis()->SetTitle("log_{10}(#hat{#Theta}_{2D})");
  h_IP2D->GetXaxis()->SetTitle("log_{10}(#hat{IP}_{2D})");
  h_alpha->GetXaxis()->SetTitle("Jet #alpha_{max}");
  
  t1.Write();
  h_nGenJets->Write();
  h_nJets->Write();
  h_trkTime->Write();
  h_jetTime->Write();
  h_jetTimeCorrected->Write();
  h_nTracksInJet->Write();
  h_theta2D->Write();
  h_IP2D->Write();
  h_alpha->Write();
  fout->Write();
  fout->Close();
}




void FillAllHistos(){
  FillHistos("XXQQQQ_m50_ctau10000mm");
  FillHistos("XXQQQQ_m50_ctau100mm");
  FillHistos("XXQQQQ_m50_ctau1mm");

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
