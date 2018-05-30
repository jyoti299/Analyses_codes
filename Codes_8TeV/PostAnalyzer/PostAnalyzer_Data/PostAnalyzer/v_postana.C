#define PostAnalyzerData_cxx
#include "PostAnalyzerData.h"

using namespace std;
using namespace ROOT;

int main(){
  //PostAnalyzerData *a = new PostAnalyzerData();
  PostAnalyzerData a;
  a.Loop();
  return 0;
}


void PostAnalyzerData::Loop()
{
  // Give values to all the parameters at the beginning
  //Luminosity
  Lumi           =  19740.0 ; //11361.7 // (pb^{-1})

  //Vertex Selection
  Cvertex_z      =  24.0 ; //(cm)
  Cvertex_ndof   =   4.  ;
  Cvertex_rho    =   2.0 ; //(cm)

  //Fiducial Cuts
  photon_pt_cut  = 170.0 ; //(GeV)
  photon_eta_cut = 1.4442;

  jet_pt_cut     = 170.0 ; //(GeV)
  jet_eta_cut    = 3.0   ;

  mass_cut       = 560.0 ; //(GeV)

  PhotJetDPhi_cut    = 1.5   ;
  PhotJetDEta_cut    = 2.0   ;

  Bool_t passHLT;
  Bool_t nonScraping, primaryVtx;
  Bool_t tightJetID;

  //Output file
  f1 = new TFile("Postanalyzerdata.root", "recreate");

  const int nbins = 13;
  TString  CutFlowLabel[nbins] ={ "Total", "HLT", "Scraping", "PrimaryVtx", "PhotonID", "PhotonPt", "PhotonEta", "JetID", "JetPt ", "JetEta ", "Dphi\
1.5", "DEta1.0", "MassCut"};
  Double_t CutFlowNumber[nbins]={  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } ;
  Double_t CutFlowBins[nbins+1]={  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
  h_CutFlowTable    =new TH1F("h_CutFlowTable"," cut flow of selection",nbins,CutFlowBins);


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  cout<<" Will Analyze = "<<nentries<<" events"<<endl;

  BookHistos();

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    PC = -1;
    JC = -1;
    goodVertex = 0;

    passHLT     = PassHLT(passedHLT150, passed90HLTall);
    nonScraping = NonScraping();
    primaryVtx  = PrimaryVertex(goodVertex);

    // Selecting photon candidate using tight photon ID.
    foundPhoton.clear();
    for(int ipho=0; ipho<Photon_n; ++ipho){
      if( NoSpike(ipho)   &&  TightPhotonPFIso(ipho) ){
	foundPhoton.push_back(ipho);
      }
    }
    if(foundPhoton.size() != 0) PC = foundPhoton[0];

    // Selecting Jet candidate using tight Jet ID.
    foundJet.clear();
    if(PC >= 0){
      for(int ijet=0; ijet<pfJet_n; ++ijet){
	Float_t drp = -1.0;
	drp = getDR(pfJet_eta[ijet],Photon_eta[PC],pfJet_phi[ijet], Photon_phi[PC]);
	if( drp > 0.5  &&  TightJetID(ijet) ){
	  foundJet.push_back(ijet);
	}
      }
    }
    if(foundJet.size() != 0) JC = foundJet[0];

    //Kinematical Cuts  ----------------
    Bool_t PhoPt    = ( Photon_pt[PC] > photon_pt_cut );
    Bool_t PhoEtaEB = ( fabs(Photon_sc_eta[PC]) < photon_eta_cut);
    //      Bool_t PhoEtaEB = ( fabs(Photon_sc_eta[PC]) < photon_eta_cut && Photon_isEB[PC]);

    Bool_t JetPt    = ( pfJet_pt[JC] > jet_pt_cut );
    Bool_t JetEtaEB = ( fabs(pfJet_eta[JC]) < jet_eta_cut);

    Bool_t GJDeltaPhi = getDPhi(Photon_phi[PC], pfJet_phi[JC])    > PhotJetDPhi_cut;
    Bool_t GJDeltaEta = getDEta(Photon_sc_eta[PC], pfJet_eta[JC]) < PhotJetDEta_cut;

    Bool_t MassCut  = ( getMass(PC,JC) > mass_cut );


    CutFlowNumber[0]++;

    if(passHLT){
      CutFlowNumber[1]++;

      if(nonScraping){
	CutFlowNumber[2]++;

	if(primaryVtx){
	  CutFlowNumber[3]++;

	  if(PC > -1){
	    CutFlowNumber[4]++;

	    if(PhoPt){
	      CutFlowNumber[5]++;
    
	      if(PhoEtaEB){
		CutFlowNumber[6]++;

		if(JC > -1){
		  CutFlowNumber[7]++;

		  if(JetPt){
		    CutFlowNumber[8]++;

		    if(JetEtaEB){
		      CutFlowNumber[9]++;
		      cout << "for entry = " << jentry << "and Photon no. = " << PC << "and jet no. = " << JC <<  ", delta phi is " << getDPhi(Photon_phi[PC], pfJet_phi[JC]) << endl;

		      if(GJDeltaPhi){
			CutFlowNumber[10]++;

			if(GJDeltaEta){
			  CutFlowNumber[11]++;

			  if(MassCut){
			    CutFlowNumber[12]++;
			    //			cout << getDPhi(Photon_phi[PC], pfJet_phi[JC]) << endl;  
			    h_nPhoton->Fill(Photon_n);
			    h_nJet->Fill(pfJet_n);
			  
			    h_PC->Fill(PC);
			    h_JC->Fill(JC);
			  
			    h_ptPhoton->Fill(Photon_pt[PC]);
			    h_etaPhoton->Fill(Photon_sc_eta[PC]);			  
                            h_ptJet->Fill(pfJet_pt[JC]);
			    h_etaJet->Fill(pfJet_eta[JC]);
			    h_mass_bin25->Fill(getMass(PC,JC));
			    h_Photon_SigmaIetaIeta->Fill(Photon_SigmaIetaIeta[PC]); 
        		    h_DR_PhotonJet->Fill(getDR(Photon_sc_eta[PC], pfJet_eta[JC],Photon_phi[PC], pfJet_phi[JC]));
			    h_dEta->Fill(getDEta(Photon_sc_eta[PC],pfJet_eta[JC]));
			    h_dphi->Fill(getDPhi(Photon_phi[PC],pfJet_phi[JC]));
			  		 
						 						 
			  } //MassCut
			}//GJDeltaEta
		      }//GJDeltaPhi
		    }//JetEtaEB
		  }//JetPt
		}//JC
	      }//PhoEtaEB
	    }//PhoPt
	  }//PC
	}//PV
      }//nonScraping
    }//passHLT

  }
    cout << " Total no. of events = " << CutFlowNumber[0] << endl;
    cout << "No. of events passing HLT = " << CutFlowNumber[1] << endl;
    cout << "No. of events passing Scraping = " << CutFlowNumber[2] << endl;
    cout << "No. of events passing PrimaryVtx = " << CutFlowNumber[3] << endl;
    cout << "No. of events passing PhotonID = " << CutFlowNumber[4] << endl;
    cout << "No. of events passing PhotonPt = " << CutFlowNumber[5] << endl;
    cout << "No. of events passing PhotonEta = " << CutFlowNumber[6] << endl;
    cout << "No. of events passing JetID = " << CutFlowNumber[7] << endl;
    cout << "No. of events passing JetPt = " << CutFlowNumber[8] << endl;
    cout << "No. of events passing JetEta = " << CutFlowNumber[9] << endl;
    cout << "No. of events passing Dphi1.5 = " << CutFlowNumber[10] << endl;
    cout << "No. of events passing DEta1.0 = " << CutFlowNumber[11] << endl;
    cout << "No. of events passing MassCut = " << CutFlowNumber[12] << endl;
}
