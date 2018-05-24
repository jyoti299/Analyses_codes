SignalEfficiencies(){

    const Int_t nSigFiles = 11;
    Int_t NoCut, Final;

    TString InSignalPath = "./InputFiles/ABCDEF/IntSignal/";
    TString OutPath = "./OutputFiles/ABCDEF/";
    
    TString FileSignals[nSigFiles] = { 
	"QstarToGJ_M_700.root",
	"QstarToGJ_M_1000.root",
	"QstarToGJ_M_1200.root",
	"QstarToGJ_M_1500.root",
	"QstarToGJ_M_1700.root",
	"QstarToGJ_M_2000.root",
	"QstarToGJ_M_2500.root",
	"QstarToGJ_M_3000.root",
	"QstarToGJ_M_3500.root",
	"QstarToGJ_M_4000.root",
	"QstarToGJ_M_4500.root"
    };


    const Int_t nvar = 2;
    TString var[nvar] = {"h_CutFlowTable", "h_CutExpFlowTable"};
    TString PlotsName[nvar] = {"CutFlowTable", "CutExpFlowTable"};
    TH1F *hSig[nSigFiles];

    // Script for plotting and saving 
    for(Int_t i = 1 ; i < nvar ; ++i){
//    for(Int_t i = 0 ; i < 1 ; ++i){
	for(Int_t j = 0 ; j < nSigFiles ; ++j){
	    TFile* fileIn = new TFile(InSignalPath+FileSignals[j],"READ");
	    hSig[j]  = (TH1F*)fileIn->Get(var[i]);   //Signal
	}
    }

    Int_t nbins = hSig[0]->GetNbinsX();

    std::vector<Double_t> ExpEntries_Sig700, ExpEntries_Sig1000, ExpEntries_Sig1200, ExpEntries_Sig1500, ExpEntries_Sig1700, ExpEntries_Sig2000, ExpEntries_Sig2500, ExpEntries_Sig3000, ExpEntries_Sig3500, ExpEntries_Sig4000, ExpEntries_Sig4500;
    std::vector<TString> Labels;

    ofstream fout(OutPath+"SignalEfficiencies.tex");
    fout<<"Cuts     &  700  & 1000 & 1200 & 1500 & 1700 & 2000 & 2500 & 3000 & 3500 & 4000 & 4500 \\\\"<<endl<<"\\hline \\hline"<<endl;

    for(Int_t k = 1; k <= nbins ; ++k){ 
	ExpEntries_Sig700.push_back(hSig[0]->GetBinContent(k));
	ExpEntries_Sig1000.push_back(hSig[1]->GetBinContent(k));
	ExpEntries_Sig1200.push_back(hSig[2]->GetBinContent(k));
	ExpEntries_Sig1500.push_back(hSig[3]->GetBinContent(k));
	ExpEntries_Sig1700.push_back(hSig[4]->GetBinContent(k));
	ExpEntries_Sig2000.push_back(hSig[5]->GetBinContent(k));
	ExpEntries_Sig2500.push_back(hSig[6]->GetBinContent(k));
	ExpEntries_Sig3000.push_back(hSig[7]->GetBinContent(k));
	ExpEntries_Sig3500.push_back(hSig[8]->GetBinContent(k));
	ExpEntries_Sig4000.push_back(hSig[9]->GetBinContent(k));
	ExpEntries_Sig4500.push_back(hSig[10]->GetBinContent(k));

	Labels.push_back(hSig[0]->GetXaxis()->GetBinLabel(k));

	fout<<Labels[k-1]<<" & "<<ExpEntries_Sig700[k-1]<<" &  "<<ExpEntries_Sig1000[k-1]<<" &  "<<ExpEntries_Sig1200[k-1]<<" &  "<<ExpEntries_Sig1500[k-1]<<" &  "<<ExpEntries_Sig1700[k-1]<<" &  "<<ExpEntries_Sig2000[k-1]<<" &  "<<ExpEntries_Sig2500[k-1]<<" &  "<<ExpEntries_Sig3000[k-1]<<" &  "<<ExpEntries_Sig3500[k-1]<<" &  "<<ExpEntries_Sig4000[k-1]<<" &  "<<ExpEntries_Sig4500[k-1]<<" \\"<<"\\"<<endl;


	if(Labels[k-1] == "Total") NoCut = k-1;
	if(Labels[k-1] == "MassCut") Final = k-1;
    }

   //std::cout<<NoCut<<"  "<<Final<<std::endl;
    fout<<"\\hline"<<std::endl;
    fout<<"Efficiency" <<" & "<<ExpEntries_Sig700[Final]/ExpEntries_Sig700[NoCut]<<" &  "<<ExpEntries_Sig1000[Final]/ExpEntries_Sig1000[NoCut]<<" &  "<<ExpEntries_Sig1200[Final]/ExpEntries_Sig1200[NoCut]<<" &  "<<ExpEntries_Sig1500[Final]/ExpEntries_Sig1500[NoCut]<<" &  "<<ExpEntries_Sig1700[Final]/ExpEntries_Sig1700[NoCut]<<" &  "<<ExpEntries_Sig2000[Final]/ExpEntries_Sig2000[NoCut]<<" &  "<<ExpEntries_Sig2500[Final]/ExpEntries_Sig2500[NoCut]<<" &  "<<ExpEntries_Sig3000[Final]/ExpEntries_Sig3000[NoCut]<<" &  "<<ExpEntries_Sig3500[Final]/ExpEntries_Sig3500[NoCut]<<" &  "<<ExpEntries_Sig4000[Final]/ExpEntries_Sig4000[NoCut]<<" &  "<<ExpEntries_Sig4500[Final]/ExpEntries_Sig4500[NoCut]<<" \\"<<"\\"<<endl;
    fout<<"\\hline"<<std::endl;
 
    fout.close();
}
