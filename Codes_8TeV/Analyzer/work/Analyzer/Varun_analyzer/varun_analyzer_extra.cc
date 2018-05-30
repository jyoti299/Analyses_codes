
    ////.............................................muon collection info.................................................
    if(runmuons_){
	edm::Handle<edm::View<pat::Muon> > muonHandle;
	iEvent.getByLabel(muoLabel_,muonHandle);
	vector <pat::Muon> mymuon_container;

	const edm::View<pat::Muon> & muons = *muonHandle;   // const ... &, we don't make a copy of it!
	for(edm::View<pat::Muon>::const_iterator muon = muons.begin(); muon!=muons.end(); ++muon){
	    //if(muon->isGlobalMuon())
	    mymuon_container.push_back(*muon);
	}
	Muon_n = 0;
	for(unsigned int x=0;x < min(mymuon_container.size(),MaxN); x++){
	    muon_pt[x]  = mymuon_container[x].pt();
	    muon_energy[x]  = mymuon_container[x].energy();
	    muon_px[x]  = mymuon_container[x].px();
	    muon_py[x]  = mymuon_container[x].py();
	    muon_pz[x]  = mymuon_container[x].pz();
	    muon_phi[x] = correct_phi(mymuon_container[x].phi());
	    muon_eta[x] = mymuon_container[x].eta();
	    muon_charge[x] = mymuon_container[x].charge();
	    muon_vx[x]  = mymuon_container[x].vx();
	    muon_vy[x]  = mymuon_container[x].vy();
	    muon_vz[x]  = mymuon_container[x].vz();

	    muon_trackIso[x] = mymuon_container[x].trackIso();
	    muon_ecalIso[x] = mymuon_container[x].ecalIso();
	    muon_hcalIso[x] = mymuon_container[x].hcalIso();
	    muon_relIso[x] =(muon_trackIso[x] + muon_ecalIso[x] + muon_hcalIso[x])/muon_pt[x];

	    muon_normChi2[x]= -99;
	    muon_validHits[x]= -99;

	    if(mymuon_container[x].globalTrack().isNonnull() ){
		muon_normChi2[x] = mymuon_container[x].normChi2();
		muon_validHits[x] = mymuon_container[x].globalTrack()->hitPattern().numberOfValidMuonHits();
	    }

	    muon_numberOfMatches[x] = mymuon_container[x].numberOfMatches();

	    muon_tkHits[x] =-99;
	    muon_pixHits[x]=-99;

	    if(mymuon_container[x].track().isNonnull() ){
		muon_tkHits[x] = mymuon_container[x].track()->numberOfValidHits();
		muon_pixHits[x] = mymuon_container[x].track()->hitPattern().numberOfValidPixelHits();
	    }


	    //tia's stuff
	    muon_isGlobalMuon[x] =     mymuon_container[x].isGlobalMuon();
	    muon_isTrackerMuon[x] =    mymuon_container[x].isTrackerMuon();
	    muon_isStandAloneMuon[x] = mymuon_container[x].isStandAloneMuon();

	    muon_OuterTrack_InnerPoint_x[x] = 0.;
	    muon_OuterTrack_InnerPoint_y[x] = 0.;
	    muon_OuterTrack_InnerPoint_z[x] = 0.;
	    muon_OuterTrack_InnerPoint_px[x] = 0.;
	    muon_OuterTrack_InnerPoint_py[x] = 0.;
	    muon_OuterTrack_InnerPoint_pz[x] = 0.;
	    muon_OuterTrack_OuterPoint_x[x] = 0.;
	    muon_OuterTrack_OuterPoint_y[x] = 0.;
	    muon_OuterTrack_OuterPoint_z[x] = 0.;
	    muon_OuterTrack_OuterPoint_px[x] = 0.;
	    muon_OuterTrack_OuterPoint_py[x] = 0.;
	    muon_OuterTrack_OuterPoint_pz[x] = 0.;
	    muon_InnerTrack_InnerPoint_x[x] = 0.;
	    muon_InnerTrack_InnerPoint_y[x] = 0.;
	    muon_InnerTrack_InnerPoint_z[x] = 0.;
	    muon_InnerTrack_InnerPoint_px[x] = 0.;
	    muon_InnerTrack_InnerPoint_py[x] = 0.;
	    muon_InnerTrack_InnerPoint_pz[x] = 0.;
	    muon_InnerTrack_OuterPoint_x[x] = 0.;
	    muon_InnerTrack_OuterPoint_y[x] = 0.;
	    muon_InnerTrack_OuterPoint_z[x] = 0.;
	    muon_InnerTrack_OuterPoint_px[x] = 0.;
	    muon_InnerTrack_OuterPoint_py[x] = 0.;
	    muon_InnerTrack_OuterPoint_pz[x] = 0.;

	    if(isAOD_){
		//FIXME: Seems needed  top and bottom referene point, but not sure how to do that !!!
		reco::TrackRef moTrkref;
		if((muon_isGlobalMuon[x]) || (muon_isTrackerMuon[x])){
		    moTrkref = mymuon_container[x].innerTrack();
		    muon_InnerTrack_isNonnull[x] = mymuon_container[x].innerTrack().isNonnull();
		}
		else{   
		    moTrkref = mymuon_container[x].outerTrack();
		    muon_OuterTrack_isNonnull[x] =   mymuon_container[x].outerTrack().isNonnull();
		}
		//For StandAlone
		muon_OuterPoint_x[x]=0.0;
		muon_OuterPoint_y[x]=0.0;
		muon_OuterPoint_z[x]=0.0;
		//For Global,Tracker
		muon_InnerPoint_x[x]=0.0;
		muon_InnerPoint_y[x]=0.0;
		muon_InnerPoint_z[x]=0.0;

		if((muon_OuterTrack_isNonnull[x])){//stand-alone
		    muon_OuterPoint_x[x]= moTrkref->referencePoint().x();
		    muon_OuterPoint_y[x]= moTrkref->referencePoint().y();
		    muon_OuterPoint_z[x]= moTrkref->referencePoint().z();
		}
		if(muon_InnerTrack_isNonnull[x]){//global,tracker
		    muon_InnerPoint_x[x]= moTrkref->referencePoint().x();
		    muon_InnerPoint_y[x]= moTrkref->referencePoint().y();
		    muon_InnerPoint_z[x]= moTrkref->referencePoint().z();
		}
	    }//if(AOD_)


	    if(!isAOD_){
		muon_InnerTrack_isNonnull[x] =   mymuon_container[x].innerTrack().isNonnull();
		muon_OuterTrack_isNonnull[x] =   mymuon_container[x].outerTrack().isNonnull();


		if(mymuon_container[x].innerTrack().isNonnull()){
		    muon_InnerTrack_InnerPoint_x[x] = mymuon_container[x].innerTrack()->innerPosition().x();
		    muon_InnerTrack_InnerPoint_y[x] = mymuon_container[x].innerTrack()->innerPosition().y();
		    muon_InnerTrack_InnerPoint_z[x] = mymuon_container[x].innerTrack()->innerPosition().z();
		    muon_InnerTrack_InnerPoint_px[x] = mymuon_container[x].innerTrack()->innerMomentum().x();
		    muon_InnerTrack_InnerPoint_py[x] = mymuon_container[x].innerTrack()->innerMomentum().y();
		    muon_InnerTrack_InnerPoint_pz[x] = mymuon_container[x].innerTrack()->innerMomentum().z();
		    muon_InnerTrack_OuterPoint_x[x] = mymuon_container[x].innerTrack()->outerPosition().x();
		    muon_InnerTrack_OuterPoint_y[x] = mymuon_container[x].innerTrack()->outerPosition().y();
		    muon_InnerTrack_OuterPoint_z[x] = mymuon_container[x].innerTrack()->outerPosition().z();
		    muon_InnerTrack_OuterPoint_px[x] = mymuon_container[x].innerTrack()->outerMomentum().x();
		    muon_InnerTrack_OuterPoint_py[x] = mymuon_container[x].innerTrack()->outerMomentum().y();
		    muon_InnerTrack_OuterPoint_pz[x] = mymuon_container[x].innerTrack()->outerMomentum().z();
		}
		if(mymuon_container[x].outerTrack().isNonnull()){
		    muon_OuterTrack_InnerPoint_x[x] = mymuon_container[x].outerTrack()->innerPosition().x();
		    muon_OuterTrack_InnerPoint_y[x] = mymuon_container[x].outerTrack()->innerPosition().y();
		    muon_OuterTrack_InnerPoint_z[x] = mymuon_container[x].outerTrack()->innerPosition().z();
		    muon_OuterTrack_InnerPoint_px[x] = mymuon_container[x].outerTrack()->innerMomentum().x();
		    muon_OuterTrack_InnerPoint_py[x] = mymuon_container[x].outerTrack()->innerMomentum().y();
		    muon_OuterTrack_InnerPoint_pz[x] = mymuon_container[x].outerTrack()->innerMomentum().z();
		    muon_OuterTrack_OuterPoint_x[x] = mymuon_container[x].outerTrack()->outerPosition().x();
		    muon_OuterTrack_OuterPoint_y[x] = mymuon_container[x].outerTrack()->outerPosition().y();
		    muon_OuterTrack_OuterPoint_z[x] = mymuon_container[x].outerTrack()->outerPosition().z();
		    muon_OuterTrack_OuterPoint_px[x] = mymuon_container[x].outerTrack()->outerMomentum().x();
		    muon_OuterTrack_OuterPoint_py[x] = mymuon_container[x].outerTrack()->outerMomentum().y();
		    muon_OuterTrack_OuterPoint_pz[x] = mymuon_container[x].outerTrack()->outerMomentum().z();
		}
	    }//if(!isAOD_)
	    Muon_n++;
	}//end of for loop
    }// if runmuons_
    /////....................muon info upto here.......................................................................
    ///xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    ///..............................cosmic muons info.....................................................................................
    if(runcosmicmuons_){
	//cosmic muon
	edm::Handle<reco::MuonCollection>cosmicMuonHandle;
	iEvent.getByLabel(cosMuoLabel_,cosmicMuonHandle);
	const reco::MuonCollection & cosmicmuons = *cosmicMuonHandle;
	vector <reco::Muon> mycosmicmuon_container;

	for(reco::MuonCollection::const_iterator cosmuon = cosmicmuons.begin(); cosmuon!=cosmicmuons.end(); ++cosmuon){
	    mycosmicmuon_container.push_back(*cosmuon);
	}

	CosmicMuon_n = 0;
	for(unsigned int x=0;x < min(mycosmicmuon_container.size(),MaxN);x++){

	    cosmicmuon_pt[x]  = mycosmicmuon_container[x].pt();
	    cosmicmuon_energy[x]  = mycosmicmuon_container[x].energy();
	    cosmicmuon_px[x]  = mycosmicmuon_container[x].px();
	    cosmicmuon_py[x]  = mycosmicmuon_container[x].py();
	    cosmicmuon_pz[x]  = mycosmicmuon_container[x].pz();
	    cosmicmuon_phi[x] = correct_phi(mycosmicmuon_container[x].phi());
	    cosmicmuon_eta[x] = mycosmicmuon_container[x].eta();
	    cosmicmuon_charge[x] = mycosmicmuon_container[x].charge();

	    //tia's stuff
	    cosmicmuon_isGlobalMuon[x] =     mycosmicmuon_container[x].isGlobalMuon();
	    cosmicmuon_isTrackerMuon[x] =    mycosmicmuon_container[x].isTrackerMuon();
	    cosmicmuon_isStandAloneMuon[x] = mycosmicmuon_container[x].isStandAloneMuon();

	    cosmicmuon_OuterTrack_InnerPoint_x[x] = 0.;
	    cosmicmuon_OuterTrack_InnerPoint_y[x] = 0.;
	    cosmicmuon_OuterTrack_InnerPoint_z[x] = 0.;
	    cosmicmuon_OuterTrack_InnerPoint_px[x] = 0.;
	    cosmicmuon_OuterTrack_InnerPoint_py[x] = 0.;
	    cosmicmuon_OuterTrack_InnerPoint_pz[x] = 0.;
	    cosmicmuon_OuterTrack_OuterPoint_x[x] = 0.;
	    cosmicmuon_OuterTrack_OuterPoint_y[x] = 0.;
	    cosmicmuon_OuterTrack_OuterPoint_z[x] = 0.;
	    cosmicmuon_OuterTrack_OuterPoint_px[x] = 0.;
	    cosmicmuon_OuterTrack_OuterPoint_py[x] = 0.;
	    cosmicmuon_OuterTrack_OuterPoint_pz[x] = 0.;
	    cosmicmuon_InnerTrack_InnerPoint_x[x] = 0.;
	    cosmicmuon_InnerTrack_InnerPoint_y[x] = 0.;
	    cosmicmuon_InnerTrack_InnerPoint_z[x] = 0.;
	    cosmicmuon_InnerTrack_InnerPoint_px[x] = 0.;
	    cosmicmuon_InnerTrack_InnerPoint_py[x] = 0.;
	    cosmicmuon_InnerTrack_InnerPoint_pz[x] = 0.;
	    cosmicmuon_InnerTrack_OuterPoint_x[x] = 0.;
	    cosmicmuon_InnerTrack_OuterPoint_y[x] = 0.;
	    cosmicmuon_InnerTrack_OuterPoint_z[x] = 0.;
	    cosmicmuon_InnerTrack_OuterPoint_px[x] = 0.;
	    cosmicmuon_InnerTrack_OuterPoint_py[x] = 0.;
	    cosmicmuon_InnerTrack_OuterPoint_pz[x] = 0.;


	    if(isAOD_){
		reco::TrackRef cmoTrkref = mycosmicmuon_container[x].outerTrack();
		cosmicmuon_OuterTrack_isNonnull[x] = mycosmicmuon_container[x].outerTrack().isNonnull();    

		cosmicmuon_OuterPoint_x[x]=0.0;
		cosmicmuon_OuterPoint_y[x]=0.0;
		cosmicmuon_OuterPoint_z[x]=0.0;
		// standalone muon variables
		if(cosmicmuon_OuterTrack_isNonnull[x]){
		    cosmicmuon_OuterPoint_x[x]= cmoTrkref->referencePoint().x();
		    cosmicmuon_OuterPoint_y[x]= cmoTrkref->referencePoint().y();
		    cosmicmuon_OuterPoint_z[x]= cmoTrkref->referencePoint().z();}
	    }//if(AOD_)

	    if(!isAOD_){ 
		cosmicmuon_InnerTrack_isNonnull[x] =   mycosmicmuon_container[x].innerTrack().isNonnull();
		cosmicmuon_OuterTrack_isNonnull[x] =   mycosmicmuon_container[x].outerTrack().isNonnull();

		if(mycosmicmuon_container[x].innerTrack().isNonnull()){
		    cosmicmuon_InnerTrack_InnerPoint_x[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().x();
		    cosmicmuon_InnerTrack_InnerPoint_y[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().y();
		    cosmicmuon_InnerTrack_InnerPoint_z[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().z();
		    cosmicmuon_InnerTrack_InnerPoint_px[x] = mycosmicmuon_container[x].innerTrack()->innerMomentum().x();
		    cosmicmuon_InnerTrack_InnerPoint_py[x] = mycosmicmuon_container[x].innerTrack()->innerMomentum().y();
		    cosmicmuon_InnerTrack_InnerPoint_pz[x] = mycosmicmuon_container[x].innerTrack()->innerMomentum().z();
		    cosmicmuon_InnerTrack_OuterPoint_x[x] = mycosmicmuon_container[x].innerTrack()->outerPosition().x();
		    cosmicmuon_InnerTrack_OuterPoint_y[x] = mycosmicmuon_container[x].innerTrack()->outerPosition().y();
		    cosmicmuon_InnerTrack_OuterPoint_z[x] = mycosmicmuon_container[x].innerTrack()->outerPosition().z();
		    cosmicmuon_InnerTrack_OuterPoint_px[x] = mycosmicmuon_container[x].innerTrack()->outerMomentum().x();
		    cosmicmuon_InnerTrack_OuterPoint_py[x] = mycosmicmuon_container[x].innerTrack()->outerMomentum().y();
		    cosmicmuon_InnerTrack_OuterPoint_pz[x] = mycosmicmuon_container[x].innerTrack()->outerMomentum().z();
		}
		if(mycosmicmuon_container[x].outerTrack().isNonnull()){
		    cosmicmuon_OuterTrack_InnerPoint_x[x] = mycosmicmuon_container[x].outerTrack()->innerPosition().x();
		    cosmicmuon_OuterTrack_InnerPoint_y[x] = mycosmicmuon_container[x].outerTrack()->innerPosition().y();
		    cosmicmuon_OuterTrack_InnerPoint_z[x] = mycosmicmuon_container[x].outerTrack()->innerPosition().z();
		    cosmicmuon_OuterTrack_InnerPoint_px[x] = mycosmicmuon_container[x].outerTrack()->innerMomentum().x();
		    cosmicmuon_OuterTrack_InnerPoint_py[x] = mycosmicmuon_container[x].outerTrack()->innerMomentum().y();
		    cosmicmuon_OuterTrack_InnerPoint_pz[x] = mycosmicmuon_container[x].outerTrack()->innerMomentum().z();
		    cosmicmuon_OuterTrack_OuterPoint_x[x] = mycosmicmuon_container[x].outerTrack()->outerPosition().x();
		    cosmicmuon_OuterTrack_OuterPoint_y[x] = mycosmicmuon_container[x].outerTrack()->outerPosition().y();
		    cosmicmuon_OuterTrack_OuterPoint_z[x] = mycosmicmuon_container[x].outerTrack()->outerPosition().z();
		    cosmicmuon_OuterTrack_OuterPoint_px[x] = mycosmicmuon_container[x].outerTrack()->outerMomentum().x();
		    cosmicmuon_OuterTrack_OuterPoint_py[x] = mycosmicmuon_container[x].outerTrack()->outerMomentum().y();
		    cosmicmuon_OuterTrack_OuterPoint_pz[x] = mycosmicmuon_container[x].outerTrack()->outerMomentum().z();
		}
	    }//if(!isAOD_)
	    CosmicMuon_n++;
	}//end of for loop
    }//if runcosmicmuons_
    ////////.......................................cosmic muon info upto here.......................................................
    //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    if(runCSCseg_){ //Store CSC segments
	// Get the CSC Geometry :
	ESHandle<CSCGeometry> cscGeom;
	iSetup.get<MuonGeometryRecord>().get(cscGeom);
	// get CSC segment collection
	Handle<CSCSegmentCollection> cscSegments;
	iEvent.getByLabel(cscLabel_, cscSegments);
	CSCseg_n = 0;

	for(CSCSegmentCollection::const_iterator it=cscSegments->begin(); CSCseg_n<10000 && it != cscSegments->end(); it++){
	    CSCDetId id  = (CSCDetId)it->cscDetId();
	    LocalVector segDir = it->localDirection();
	    LocalPoint localPos = it->localPosition();

	    // global transformation
	    const CSCChamber* cscchamber = cscGeom->chamber(id);
	    if (cscchamber) {

		GlobalPoint globalPosition = cscchamber->toGlobal(localPos);
		GlobalVector globalDirection = cscchamber->toGlobal(segDir);

		int nhits = 0;
		float timeSum = 0;
		// Get the CSC recHits that contribute to this segment.
		std::vector<CSCRecHit2D> theseRecHits = it->specificRecHits();
		for( vector<CSCRecHit2D>::const_iterator iRH =theseRecHits.begin(); iRH != theseRecHits.end(); iRH++){
		    if( !(iRH->isValid()) ) continue;  // only interested in valid hits
		    nhits++;
		    float rhTime = iRH->tpeak();
		    timeSum += rhTime;
		}//end rechit loop

		float segmentTime = -999;
		if (nhits>0)
		    segmentTime = timeSum/nhits;

		CSCseg_time[CSCseg_n] = segmentTime;
		CSCseg_x[CSCseg_n]    = globalPosition.x();
		CSCseg_y[CSCseg_n]    = globalPosition.y();
		CSCseg_z[CSCseg_n]    = globalPosition.z();
		CSCseg_phi[CSCseg_n]  = globalPosition.phi();

		CSCseg_DirectionX[CSCseg_n] = globalDirection.x();
		CSCseg_DirectionY[CSCseg_n] = globalDirection.y();
		CSCseg_DirectionZ[CSCseg_n] = globalDirection.z();

		CSCseg_n++;
	    }
	}// loop over segs
    }// runCSCseg_


    if(runRPChit_){
	/*
	// Get the RPC Geometry  
	edm::ESHandle<RPCGeometry> rpcGeo;
	iSetup.get<MuonGeometryRecord>().get(rpcGeo);

	//get RPC hit collection
	edm::Handle<RPCRecHitCollection> rpcHits;
	iEvent.getByLabel("cscSegments",rpcHits);

	for(RPCRecHitCollection::const_iterator RecHitsIt = rpcRecHitRange->begin(); RecHitsIt!=rpcRecHitRange->end(); RecHitsIt++){


	}//loop over rpc hits
	 */

	//get geometry of tracking stuff
	//edm::ESHandle<GlobalTrackingGeometry> geometry;
	edm::ESHandle<RPCGeometry> geometry;
	//iSetup.get<GlobalTrackingGeometryRecord>().get(geometry);
	iSetup.get<MuonGeometryRecord>().get(geometry);
	//get RPC hit collection
	edm::Handle<RPCRecHitCollection> rpcHits;
	iEvent.getByLabel(rpcLabel_,rpcHits);

	RPChit_n=0;
	for(RPCRecHitCollection::const_iterator RecHitsIt = rpcHits->begin(); RPChit_n<10000 && RecHitsIt!=rpcHits->end(); RecHitsIt++){
	    const GeomDet *det = geometry->idToDet(RecHitsIt->geographicalId());
	    GlobalPoint pos = det->toGlobal(RecHitsIt->localPosition());
	    RPChit_x[RPChit_n] = pos.x();
	    RPChit_y[RPChit_n] = pos.y();
	    RPChit_z[RPChit_n] = pos.z();
	    //BunchX is in units of 25ns. BunchX=0 is the expected value for collision muons.
	    RPChit_BunchX[RPChit_n] = RecHitsIt->BunchX();
	    RPChit_n++;
	}//loop over rpc hits

    }//runRPChit_

    //calomet variables  
    if(runmet_){
	edm::Handle<edm::View<pat::MET> > metHandle;
	iEvent.getByLabel(metLabel_,metHandle);
	if ( metHandle.isValid() ){
	    const edm::View<pat::MET> & met = *metHandle;
	    //edm::View<pat::MET>::const_iterator met;
	    CaloMetSig                            = met[0].mEtSig();
	    CaloEtFractionHadronic                = met[0].etFractionHadronic();
	    CaloEmEtFraction                      = met[0].emEtFraction();
	    CaloHadEtInHB                         = met[0].hadEtInHB();
	    CaloHadEtInHO                         = met[0].hadEtInHO();
	    CaloHadEtInHE                         = met[0].hadEtInHE();
	    CaloHadEtInHF                         = met[0].hadEtInHF();
	    CaloEmEtInEB                          = met[0].emEtInEB();
	    CaloEmEtInEE                          = met[0].emEtInEE();
	    CaloEmEtInHF                          = met[0].emEtInHF();
	    CaloMetEz                             = met[0].e_longitudinal();
	    CaloMaxEtInEmTowers                   = met[0].maxEtInEmTowers();
	    CaloMaxEtInHadTowers                  = met[0].maxEtInHadTowers();
	    // 0=full corr1=uncorrNone  2=uncorrALL 3=uncorrJES  4=uncorrMUON  5=TAU

	    CaloMetPt[0]                          = met[0].pt();
	    CaloMetPx[0]                          = met[0].px();
	    CaloMetPy[0]                          = met[0].py();
	    CaloMetPhi[0]                         = correct_phi(met[0].phi());
	    CaloMetSumEt[0]                       = met[0].sumEt();

	    CaloMetPt[1]                          = met[0].uncorrectedPt(pat::MET::uncorrNONE); 
	    CaloMetPhi[1]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrNONE)); 
	    CaloMetPx[1]                          = met[0].corEx(pat::MET::uncorrNONE); 
	    CaloMetPy[1]                          = met[0].corEy(pat::MET::uncorrNONE); 
	    CaloMetSumEt[1]                       = met[0].corSumEt(pat::MET::uncorrNONE); 

	    CaloMetPt[2]                          = met[0].uncorrectedPt(pat::MET::uncorrALL); 
	    CaloMetPhi[2]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrALL)); 
	    CaloMetPx[2]                          = met[0].corEx(pat::MET::uncorrALL); 
	    CaloMetPy[2]                          = met[0].corEy(pat::MET::uncorrALL); 
	    CaloMetSumEt[2]                       = met[0].corSumEt(pat::MET::uncorrALL); 

	    CaloMetPt[3]                          = met[0].uncorrectedPt(pat::MET::uncorrJES); 
	    CaloMetPhi[3]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrJES)); 
	    CaloMetPx[3]                          = met[0].corEx(pat::MET::uncorrJES); 
	    CaloMetPy[3]                          = met[0].corEy(pat::MET::uncorrJES); 
	    CaloMetSumEt[3]                       = met[0].corSumEt(pat::MET::uncorrJES); 

	    CaloMetPt[4]                          = met[0].uncorrectedPt(pat::MET::uncorrMUON); 
	    CaloMetPhi[4]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrMUON)); 
	    CaloMetPx[4]                          = met[0].corEx(pat::MET::uncorrMUON); 
	    CaloMetPy[4]                          = met[0].corEy(pat::MET::uncorrMUON); 
	    CaloMetSumEt[4]                       = met[0].corSumEt(pat::MET::uncorrMUON); 

	    CaloMetPt[5]                          = met[0].uncorrectedPt(pat::MET::uncorrTAU); 
	    CaloMetPhi[5]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrTAU)); 
	    CaloMetPx[5]                          = met[0].corEx(pat::MET::uncorrTAU); 
	    CaloMetPy[5]                          = met[0].corEy(pat::MET::uncorrTAU); 
	    CaloMetSumEt[5]                       = met[0].corSumEt(pat::MET::uncorrTAU); 
	    if(runphotons_==1)
		if (myphoton_container.size()!=0)
		    Delta_phi                       = fabs(reco::deltaPhi(correct_phi(met[0].phi()),correct_phi(myphoton_container[0].phi())));
	    if(rungenmet_){
		const reco::GenMET *genMet = met[0].genMET();
		genMetPt     = genMet->et();
		genMetPhi    = correct_phi(genMet->phi());
		genMetSumEt  = genMet->sumEt();
		genMetPx     = genMet->px();
		genMetPy     = genMet->py();
		if(runphotons_==1)
		    if (myphoton_container.size()!=0)
			Delta_phiGEN                        = fabs(reco::deltaPhi(correct_phi(genMet->phi()),correct_phi(myphoton_container[0].phi())));
	    }
	}
    }
    if(runPFmet_){
	edm::Handle<edm::View<pat::MET> > metPFHandle;
	iEvent.getByLabel(PFmetLabel_,metPFHandle);
	const edm::View<pat::MET> & metsPF = *metPFHandle;
	if ( metPFHandle.isValid() ){
	    // 0=full corr1=uncorrNone  2=uncorrALL 3=uncorrJES  4=uncorrMUON  5=TAU
	    PFMetPt[0]                          = metsPF[0].pt();
	    PFMetPx[0]                          = metsPF[0].px();
	    PFMetPy[0]                          = metsPF[0].py();
	    PFMetPhi[0]                         = correct_phi(metsPF[0].phi());
	    PFMetSumEt[0]                       = metsPF[0].sumEt();

	    PFMetPt[1]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrNONE); 
	    PFMetPhi[1]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrNONE)); 
	    PFMetPx[1]                          = metsPF[0].corEx(pat::MET::uncorrNONE); 
	    PFMetPy[1]                          = metsPF[0].corEy(pat::MET::uncorrNONE); 
	    PFMetSumEt[1]                       = metsPF[0].corSumEt(pat::MET::uncorrNONE); 

	    PFMetPt[2]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrALL); 
	    PFMetPhi[2]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrALL)); 
	    PFMetPx[2]                          = metsPF[0].corEx(pat::MET::uncorrALL); 
	    PFMetPy[2]                          = metsPF[0].corEy(pat::MET::uncorrALL); 
	    PFMetSumEt[2]                       = metsPF[0].corSumEt(pat::MET::uncorrALL); 

	    PFMetPt[3]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrJES); 
	    PFMetPhi[3]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrJES)); 
	    PFMetPx[3]                          = metsPF[0].corEx(pat::MET::uncorrJES); 
	    PFMetPy[3]                          = metsPF[0].corEy(pat::MET::uncorrJES); 
	    PFMetSumEt[3]                       = metsPF[0].corSumEt(pat::MET::uncorrJES); 

	    PFMetPt[4]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrMUON); 
	    PFMetPhi[4]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrMUON)); 
	    PFMetPx[4]                          = metsPF[0].corEx(pat::MET::uncorrMUON); 
	    PFMetPy[4]                          = metsPF[0].corEy(pat::MET::uncorrMUON); 
	    PFMetSumEt[4]                       = metsPF[0].corSumEt(pat::MET::uncorrMUON); 

	    PFMetPt[5]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrTAU); 
	    PFMetPhi[5]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrTAU)); 
	    PFMetPx[5]                          = metsPF[0].corEx(pat::MET::uncorrTAU); 
	    PFMetPy[5]                          = metsPF[0].corEy(pat::MET::uncorrTAU); 
	    PFMetSumEt[5]                       = metsPF[0].corSumEt(pat::MET::uncorrTAU); 

	    if(runphotons_==1)
		if (myphoton_container.size()!=0)
		    Delta_phiPF  = fabs(reco::deltaPhi(PFMetPhi[0],correct_phi(myphoton_container[0].phi())));
	}
	else{
	    LogWarning("METEventSelector") << "No Met results for InputTag " ;
	    return;
	}
    }
    if(runTCmet_){
	edm::Handle<edm::View<pat::MET> > metTCHandle;
	iEvent.getByLabel(TCmetLabel_,metTCHandle);
	const edm::View<pat::MET> & metsTC = *metTCHandle;
	if ( metTCHandle.isValid() ){
	    // 0=full corr1=uncorrNone  2=uncorrALL 3=uncorrJES  4=uncorrMUON  5=TAU
	    TCMetPt[0]                          = metsTC[0].pt();
	    TCMetPx[0]                          = metsTC[0].px();
	    TCMetPy[0]                          = metsTC[0].py();
	    TCMetPhi[0]                         = correct_phi(metsTC[0].phi());
	    TCMetSumEt[0]                       = metsTC[0].sumEt();

	    TCMetPt[1]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrNONE); 
	    TCMetPhi[1]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrNONE)); 
	    TCMetPx[1]                          = metsTC[0].corEx(pat::MET::uncorrNONE); 
	    TCMetPy[1]                          = metsTC[0].corEy(pat::MET::uncorrNONE); 
	    TCMetSumEt[1]                       = metsTC[0].corSumEt(pat::MET::uncorrNONE); 

	    TCMetPt[2]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrALL); 
	    TCMetPhi[2]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrALL)); 
	    TCMetPx[2]                          = metsTC[0].corEx(pat::MET::uncorrALL); 
	    TCMetPy[2]                          = metsTC[0].corEy(pat::MET::uncorrALL); 
	    TCMetSumEt[2]                       = metsTC[0].corSumEt(pat::MET::uncorrALL); 

	    TCMetPt[3]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrJES); 
	    TCMetPhi[3]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrJES)); 
	    TCMetPx[3]                          = metsTC[0].corEx(pat::MET::uncorrJES); 
	    TCMetPy[3]                          = metsTC[0].corEy(pat::MET::uncorrJES); 
	    TCMetSumEt[3]                       = metsTC[0].corSumEt(pat::MET::uncorrJES); 

	    TCMetPt[4]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrMUON); 
	    TCMetPhi[4]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrMUON)); 
	    TCMetPx[4]                          = metsTC[0].corEx(pat::MET::uncorrMUON); 
	    TCMetPy[4]                          = metsTC[0].corEy(pat::MET::uncorrMUON); 
	    TCMetSumEt[4]                       = metsTC[0].corSumEt(pat::MET::uncorrMUON); 

	    TCMetPt[5]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrTAU); 
	    TCMetPhi[5]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrTAU)); 
	    TCMetPx[5]                          = metsTC[0].corEx(pat::MET::uncorrTAU); 
	    TCMetPy[5]                          = metsTC[0].corEy(pat::MET::uncorrTAU); 
	    TCMetSumEt[5]                       = metsTC[0].corSumEt(pat::MET::uncorrTAU); 

	    if(runphotons_==1)
		if (myphoton_container.size()!=0)
		    Delta_phiTC  = fabs(reco::deltaPhi(TCMetPhi[0],correct_phi(myphoton_container[0].phi()))); 
	}
	else{
	    LogWarning("METEventSelector") << "No Met results for InputTag " ;
	    return;
	} 

    }//end of if(runTCmet_)



    if(runtaus_){
	edm::Handle<edm::View<pat::Tau> > tauHandle;
	iEvent.getByLabel(tauLabel_,tauHandle);
	vector <pat::Tau> mytau_container;


	const edm::View<pat::Tau> & taus = *tauHandle;   // const ... &, we don't make a copy of it!
	for(edm::View<pat::Tau>::const_iterator tau = taus.begin(); tau!=taus.end(); ++tau){
	    mytau_container.push_back(*tau);
	}
	Tau_n = 0;

	if(mytau_container.size() >1){std::sort(mytau_container.begin(),mytau_container.end(),PtSortCriteriumtau());}

	for(unsigned int x=0;x < min(mytau_container.size(), (MaxN/2));x++){
	    tau_pt[x]  = mytau_container[x].pt();
	    tau_energy[x]  = mytau_container[x].energy();
	    tau_px[x]  = mytau_container[x].px();
	    tau_py[x]  = mytau_container[x].py();
	    tau_pz[x]  = mytau_container[x].pz();
	    tau_vx[x]  = mytau_container[x].vx();
	    tau_vy[x]  = mytau_container[x].vy();
	    tau_vz[x]  = mytau_container[x].vz();
	    tau_phi[x] = correct_phi(mytau_container[x].phi());
	    tau_charge[x] = mytau_container[x].charge();



	    if(runDetailTauInfo_){//-------Bhawna need this infor

		nPions[x]       =0;                  
		nPi0[x]         =0;                  
		nPhotons[x]     =0;                  
		oneProng0Pi0[x] =0;                  
		oneProng1Pi0[x] =0;                  
		oneProng2Pi0[x] =0;                  
		threeProng0Pi0[x]=0;                 
		threeProng1Pi0[x]=0;                 
		tauelectron[x]   =0;                 
		taumuon[x]       =0;      

		if(mytau_container[x].genJet())
		{      
		    genTauDecayMode = JetMCTagUtils::genTauDecayMode(*mytau_container[x].genJet());
		    genTauDecayMode1.push_back(JetMCTagUtils::genTauDecayMode(*mytau_container[x].genJet()));

		    if (genTauDecayMode == "oneProng0Pi0")
		    { oneProng0Pi0[x]=1;
			oneProng0Pi0Pt[x] = (mytau_container[x]).genJet()->pt() ;
			oneProng0Pi0Eta[x] = (mytau_container[x]).genJet()->eta();
			oneProng0Pi0Phi[x] = (mytau_container[x]).genJet()->phi();
		    }  

		    if (genTauDecayMode == "oneProng1Pi0")oneProng1Pi0[x]=1;
		    if(genTauDecayMode == "oneProng2Pi0")oneProng2Pi0[x]=1;
		    if (genTauDecayMode == "threeProng0Pi0"){ threeProng0Pi0[x]=1;
			nthreeProng0Pi0++;}
			if (genTauDecayMode == "threeProng1Pi0"){ threeProng1Pi0[x]=1;
			    nthreeProng1Pi0++;}
			    if (genTauDecayMode == "electron"){tauelectron[x]=1;
				ntauelectron++;}
				if (genTauDecayMode == "muon"){ taumuon[x]=1;
				    ntaumuon++;}

				    if (  (genTauDecayMode == "oneProng0Pi0"   ||   genTauDecayMode == "oneProng1Pi0"   ||
						genTauDecayMode == "oneProng2Pi0"   ||   genTauDecayMode == "threeProng0Pi0" ||
						genTauDecayMode == "threeProng1Pi0" ||   genTauDecayMode == "electron"       ||
						genTauDecayMode == "muon"
					  )  )   
				    {     
					genHadTauPt[x] = (mytau_container[x]).genJet()->pt() ;
					genHadTauEta[x] = (mytau_container[x]).genJet()->eta();
					genHadTauPhi[x] = (mytau_container[x]).genJet()->phi();

					if( (mytau_container[x].genJet()->getGenConstituents()).size() !=0 )
					{  
					    genParticleList = mytau_container[x].genJet()->getGenConstituents();
					    std::vector <const GenParticle*>  genParticleList = mytau_container[x].genJet()->getGenConstituents();

					    for ( int i = 0; i < 5; i++)
					    { PionPdgId[x][i]   = -999;
						PionPt[x][i]      = -999;
						PionEta[x][i]     = -999;
						PionPhi[x][i]     = -999;
						PhotonPdgId[x][i] = -999;
						PhotonPt[x][i]    = -999;
						PhotonEta[x][i]   = -999;
						PionPhi[x][i]     = -999;
						Pi0PdgId[x][i]    = -999;
						Pi0Pt[x][i]       = -999;
						Pi0Eta[x][i]      = -999;
						Pi0Phi[x][i]      = -999;
					    }                 

					    for (int iList = 0; iList <int(genParticleList.size()); iList++)
					    {

						if(abs(genParticleList[iList]->pdgId()) == 211)
						{ if(nPions[x]>5) continue;
						    PionPt[x][nPions[x]] = genParticleList[iList]->pt();
						    PionEta[x][nPions[x]] = genParticleList[iList]->eta();
						    PionPhi[x][nPions[x]] = genParticleList[iList]->phi();
						    PionPdgId[x][nPions[x]] = genParticleList[iList]->pdgId();
						    nPions[x]++;
						}
						if(abs(genParticleList[iList]->pdgId()) == 22)
						{ if(nPhotons[x]> 5)continue;
						    PhotonPt[x][nPhotons[x]] = genParticleList[iList]->pt();
						    PhotonEta[x][nPhotons[x]] = genParticleList[iList]->eta();
						    PhotonPhi[x][nPhotons[x]] = genParticleList[iList]->phi();
						    PhotonPdgId[x][nPhotons[x]]  = genParticleList[iList]->pdgId();
						    nPhotons[x]++;
						}
						if(abs(genParticleList[iList]->pdgId()) == 111)
						{ if(nPi0[x]> 5 ) continue;
						    Pi0Pt[x][nPi0[x]] = genParticleList[iList]->pt();
						    if(debug_)cout<<"Pi0Pt  = "<<Pi0Pt[x][nPi0[x]]<<endl;
						    Pi0Eta[x][nPi0[x]] = genParticleList[iList]->eta();
						    Pi0Phi[x][nPi0[x]] = genParticleList[iList]->phi();
						    Pi0PdgId[x][nPi0[x]]  = genParticleList[iList]->pdgId();
						    nPi0[x]++;
						}
					    }//end of  for (int iList = 0; iList <int(genParticleList.size()); iList++)


					}//end of if( (mytau_container[x].genJet()->getGenConstituents()).size() !=0 )
				    }//TauDeca mode
		}//genJet is true
	    }//if(runDetailTauInfor_)//---------Bhawana need


	    Tau_n++;
	}//end of for loop
    }//Run taus_



    //Below is the uncleaned photon colleciton
    std::vector<pat::Photon> ucphoton_container;
    ucphoton_container.clear();

    if(runucphotons_){

	edm::Handle<reco::PhotonCollection> UCphotonH;
	bool found = iEvent.getByLabel(inputTagUCPhotons_,UCphotonH);

	unsigned nTypes=3;
	IsoDepositVals photonIsoValPFId(nTypes);
	for (size_t j = 0; j<inputTagIsoValUCPhotonsPFId_.size(); ++j) {
	    iEvent.getByLabel(inputTagIsoValUCPhotonsPFId_[j], photonIsoValPFId[j]);
	}

	// Photons - from reco
	const IsoDepositVals * photonIsoVals = &photonIsoValPFId;
	unsigned nrecopho=UCphotonH->size();

	edm::Handle<edm::View<pat::Photon> > ucphoHandle;
	iEvent.getByLabel(ucphoLabel_,ucphoHandle);
	edm::View<pat::Photon>::const_iterator ucphoton;

	set<DetId> ucHERecHitSet;
	ucHERecHit_subset_n = 0;
	ucnpho=0;
	for(ucphoton = ucphoHandle->begin();ucphoton!=ucphoHandle->end();++ucphoton){

	    for(unsigned ipho=0; ipho<nrecopho;++ipho) {
		reco::PhotonRef myPhotonRef(UCphotonH,ipho);

		if (myPhotonRef->et() != ucphoton->et()) continue;

		ucphoElectronveto[npho] = !ConversionTools::hasMatchedPromptElectron(myPhotonRef->superCluster(), hElectrons, hConversions, beamspot.position());

		uccharged03 =  (*(*photonIsoVals)[0])[myPhotonRef];
		ucphoton03 = (*(*photonIsoVals)[1])[myPhotonRef];
		ucneutral03 = (*(*photonIsoVals)[2])[myPhotonRef];

		ucpho_PFisocharged03[ucnpho] =  ((*(*photonIsoVals)[0])[myPhotonRef]);
		ucpho_PFisophoton03[ucnpho]  = ((*(*photonIsoVals)[1])[myPhotonRef]);
		ucpho_PFisoneutral03[ucnpho] = ((*(*photonIsoVals)[2])[myPhotonRef]);
		ucpho_PFphotonssum03[ucnpho] = (uccharged03+ucphoton03+ucneutral03);
		ucnpho++;
		//std::cout<<"inside npho loop" << std::endl;
		//std::cout<<"UC Photon pt/eta/phi:" << ucphoton->et()<<"\t"<< ucphoton->eta() << "\t"<< ucphoton->phi()<<  std::endl;
	    }
	    ucphoton_container.push_back(*ucphoton) ;
	}
	ucPhoton_n = 0;
	if(ucphoton_container.size()!=0){
	    //std::cout<<"inside x loop" << std::endl;
	    for(unsigned int x_uc=0; x_uc < min(ucphoton_container.size(), MaxN);x_uc++){
		std::cout<<"UCPhoton pt/eta/phi:" << ucphoton_container[x_uc].et()<<"\t"<<ucphoton_container[x_uc].eta() << "\t"<< ucphoton_container[x_uc].phi()<<  std::endl;
		ucpho_E[x_uc]                     =  ucphoton_container[x_uc].energy();
		ucpho_pt[x_uc]                    =  ucphoton_container[x_uc].pt();
		//cout<<"Uncleaned Photon pT  = "<<ucpho_pt[x_uc]<<endl;
		ucpho_px[x_uc]                    =  ucphoton_container[x_uc].px();
		ucpho_py[x_uc]                    =  ucphoton_container[x_uc].py();
		ucpho_pz[x_uc]                    =  ucphoton_container[x_uc].pz();
		ucpho_vx[x_uc]                    =  ucphoton_container[x_uc].vx();
		ucpho_vy[x_uc]                    =  ucphoton_container[x_uc].vy();
		ucpho_vz[x_uc]                    =  ucphoton_container[x_uc].vz();
		ucpho_et[x_uc]                    =  ucphoton_container[x_uc].et();
		ucpho_eta[x_uc]                   =  ucphoton_container[x_uc].eta();
		ucpho_phi[x_uc]                   =  correct_phi(ucphoton_container[x_uc].phi());
		ucpho_theta[x_uc]                 =  ucphoton_container[x_uc].theta();
		ucpho_r9[x_uc]                    =  ucphoton_container[x_uc].r9();
		ucpho_e1x5[x_uc]                  =  ucphoton_container[x_uc].e1x5();
		ucpho_e2x5[x_uc]                  =  ucphoton_container[x_uc].e2x5();
		ucpho_e3x3[x_uc]                  =  ucphoton_container[x_uc].e3x3();
		ucpho_e5x5[x_uc]                  =  ucphoton_container[x_uc].e5x5();
		ucpho_maxEnergyXtal[x_uc]         =  ucphoton_container[x_uc].maxEnergyXtal();
		ucpho_SigmaEtaEta[x_uc]           =  ucphoton_container[x_uc].sigmaEtaEta();        
		ucpho_SigmaIetaIeta[x_uc]         =  ucphoton_container[x_uc].sigmaIetaIeta();        
		ucpho_r1x5[x_uc]                  =  ucphoton_container[x_uc].r1x5();
		ucpho_r2x5[x_uc]                  =  ucphoton_container[x_uc].r2x5();
		ucpho_size[x_uc]                  =  ucphoton_container[x_uc].superCluster()->clustersSize();
		ucpho_sc_energy[x_uc]             =  ucphoton_container[x_uc].superCluster()->energy();
		ucpho_sc_eta[x_uc]                =  ucphoton_container[x_uc].superCluster()->eta();
		ucpho_sc_phi[x_uc]                =  correct_phi(ucphoton_container[x_uc].superCluster()->phi());
		ucpho_sc_x[x_uc]                  =  ucphoton_container[x_uc].superCluster()->x();
		ucpho_sc_y[x_uc]                  =  ucphoton_container[x_uc].superCluster()->y();
		ucpho_sc_z[x_uc]                  =  ucphoton_container[x_uc].superCluster()->z();
		ucpho_sc_etaWidth[x_uc]           =  ucphoton_container[x_uc].superCluster()->etaWidth();
		ucpho_sc_phiWidth[x_uc]           =  ucphoton_container[x_uc].superCluster()->phiWidth();
		ucpho_HoE[x_uc]                   =  ucphoton_container[x_uc].hadronicOverEm();              
		ucpho_HoEnew[x_uc]                   =  ucphoton_container[x_uc].hadTowOverEm();
		ucpho_ecalRecHitSumEtConeDR03[x_uc]      =  ucphoton_container[x_uc].ecalRecHitSumEtConeDR03();
		ucpho_hcalTowerSumEtConeDR03[x_uc]       =  ucphoton_container[x_uc].hcalTowerSumEtConeDR03();
		ucpho_trkSumPtHollowConeDR03[x_uc]       =  ucphoton_container[x_uc].trkSumPtHollowConeDR03();
		ucpho_trkSumPtSolidConeDR03[x_uc]        =  ucphoton_container[x_uc].trkSumPtSolidConeDR03();
		ucpho_nTrkSolidConeDR03[x_uc]            = ucphoton_container[x_uc].nTrkSolidConeDR03();
		ucpho_nTrkHollowConeDR03[x_uc]           = ucphoton_container[x_uc].nTrkHollowConeDR03();
		ucpho_hcalDepth1TowerSumEtConeDR03[x_uc] = ucphoton_container[x_uc].hcalDepth1TowerSumEtConeDR03();
		ucpho_hcalDepth2TowerSumEtConeDR03[x_uc] = ucphoton_container[x_uc].hcalDepth2TowerSumEtConeDR03();
		ucpho_ecalRecHitSumEtConeDR04[x_uc]      =  ucphoton_container[x_uc].ecalRecHitSumEtConeDR04();
		ucpho_hcalTowerSumEtConeDR04[x_uc]       =  ucphoton_container[x_uc].hcalTowerSumEtConeDR04();
		ucpho_trkSumPtHollowConeDR04[x_uc]       =  ucphoton_container[x_uc].trkSumPtHollowConeDR04();
		ucpho_trkSumPtSolidConeDR04[x_uc]        =  ucphoton_container[x_uc].trkSumPtSolidConeDR04();
		ucpho_nTrkSolidConeDR04[x_uc]            = ucphoton_container[x_uc].nTrkSolidConeDR04();
		ucpho_nTrkHollowConeDR04[x_uc]           = ucphoton_container[x_uc].nTrkHollowConeDR04();
		ucpho_hcalDepth1TowerSumEtConeDR04[x_uc] = ucphoton_container[x_uc].hcalDepth1TowerSumEtConeDR04();
		ucpho_hcalDepth2TowerSumEtConeDR04[x_uc] = ucphoton_container[x_uc].hcalDepth2TowerSumEtConeDR04();
		ucpho_hasPixelSeed[x_uc]                 = ucphoton_container[x_uc].hasPixelSeed(); 
		ucpho_isEB[x_uc]                         = ucphoton_container[x_uc].isEB(); 
		ucpho_isEE[x_uc]                         = ucphoton_container[x_uc].isEE();
		ucpho_isEBGap[x_uc]                      = ucphoton_container[x_uc].isEBGap(); 
		ucpho_isEEGap[x_uc]                      = ucphoton_container[x_uc].isEEGap(); 
		ucpho_isEBEEGap[x_uc]                    = ucphoton_container[x_uc].isEBEEGap(); 
		ucpho_hasConvTrk[x_uc]                   = ucphoton_container[x_uc].hasConversionTracks();

		//Add MIP Variable for each unclenaed photon
		ucpho_mipChi2[x_uc]                      = ucphoton_container[x_uc].mipChi2();                                                                                      
		ucpho_mipTotEnergy[x_uc]                 = ucphoton_container[x_uc].mipTotEnergy();
		ucpho_mipSlope[x_uc]                     = ucphoton_container[x_uc].mipSlope();
		ucpho_mipIntercept[x_uc]                 = ucphoton_container[x_uc].mipIntercept();
		ucpho_mipNhitCone[x_uc]                  = ucphoton_container[x_uc].mipNhitCone();
		ucpho_mipIsHalo[x_uc]                    = ucphoton_container[x_uc].mipIsHalo();


		if(ucphoton_container[x_uc].genParticleRef().isNonnull()){
		    matchucpho_E[x_uc]                =  ucphoton_container[x_uc].genPhoton()->energy();
		    matchucpho_pt[x_uc]               =  ucphoton_container[x_uc].genPhoton()->pt();
		    matchucpho_eta[x_uc]              =  ucphoton_container[x_uc].genPhoton()->eta();
		    matchucpho_phi[x_uc]              =  correct_phi(ucphoton_container[x_uc].genPhoton()->phi());
		    matchucpho_px[x_uc]               =  ucphoton_container[x_uc].genPhoton()->px();
		    matchucpho_py[x_uc]               =  ucphoton_container[x_uc].genPhoton()->py();
		    matchucpho_pz[x_uc]               =  ucphoton_container[x_uc].genPhoton()->pz();
		}
		else{
		    matchucpho_E[x_uc]                = -99.;
		    matchucpho_pt[x_uc]               = -99.;
		    matchucpho_eta[x_uc]              = -99.;
		    matchucpho_phi[x_uc]              = -99.;
		    matchucpho_px[x_uc]               = -99.;
		    matchucpho_py[x_uc]               = -99.;
		    matchucpho_pz[x_uc]               = -99.;
		}
		ismatcheducpho[x_uc]                =  ucphoton_container[x_uc].genParticleRef().isNonnull();

		ucpho_nTracks[x_uc]                    = 9999;
		ucpho_isConverted[x_uc]                = false;
		ucpho_pairInvariantMass[x_uc]          = -99.;
		ucpho_pairCotThetaSeparation[x_uc]     = -99.;
		ucpho_pairMomentum_x[x_uc]             = -99.;
		ucpho_pairMomentum_y[x_uc]             = -99.;
		ucpho_pairMomentum_z[x_uc]             = -99.;
		ucpho_conv_vx[x_uc]                    = -99.;
		ucpho_conv_vy[x_uc]                    = -99.;
		ucpho_conv_vz[x_uc]                    = -99.;
		ucpho_EoverP[x_uc]                     = -99.;
		ucpho_zOfPrimaryVertex[x_uc]           = -99.;
		ucpho_distOfMinimumApproach[x_uc]      = -99.;
		ucpho_dPhiTracksAtVtx[x_uc]            = -99.;
		ucpho_dPhiTracksAtEcal[x_uc]           = -99.;
		ucpho_dEtaTracksAtEcal[x_uc]           = -99.;

		reco::ConversionRefVector conversions_uc   = ucphoton_container[x_uc].conversions();
		for (unsigned int iConv_uc=0; iConv_uc<conversions_uc.size(); iConv_uc++) {
		    reco::ConversionRef aConv_uc=conversions_uc[iConv_uc];
		    if ( aConv_uc->nTracks() <2 ) continue; 
		    if ( aConv_uc->conversionVertex().isValid() ){
			ucpho_nTracks[x_uc]                    = aConv_uc->nTracks();
			ucpho_isConverted[x_uc]                = aConv_uc->isConverted();
			ucpho_pairInvariantMass[x_uc]          = aConv_uc->pairInvariantMass();
			ucpho_pairCotThetaSeparation[x_uc]     = aConv_uc->pairCotThetaSeparation();
			ucpho_pairMomentum_x[x_uc]             = aConv_uc->pairMomentum().x();
			ucpho_pairMomentum_y[x_uc]             = aConv_uc->pairMomentum().y();
			ucpho_pairMomentum_z[x_uc]             = aConv_uc->pairMomentum().z();
			ucpho_conv_vx[x_uc]                    = aConv_uc->conversionVertex().x();
			ucpho_conv_vy[x_uc]                    = aConv_uc->conversionVertex().y();
			ucpho_conv_vz[x_uc]                    = aConv_uc->conversionVertex().z();
			ucpho_EoverP[x_uc]                     = aConv_uc->EoverP();
			ucpho_zOfPrimaryVertex[x_uc]           = aConv_uc->zOfPrimaryVertexFromTracks();
			ucpho_distOfMinimumApproach[x_uc]      = aConv_uc->distOfMinimumApproach();
			ucpho_dPhiTracksAtVtx[x_uc]            = aConv_uc->dPhiTracksAtVtx();
			ucpho_dPhiTracksAtEcal[x_uc]           = aConv_uc->dPhiTracksAtEcal();
			ucpho_dEtaTracksAtEcal[x_uc]           = aConv_uc->dEtaTracksAtEcal();
		    }//end of if ( aConv_uc->conversionVertex().isValid() )
		}//end of for (unsigned int iConv_uc=0; iConv_uc<conversions.size(); iConv_uc++)


		if(runHErechit_ && ucpho_isEB[x_uc]){

		    //Store HE rechits
		    //edm::Handle<HBHERecHitCollection> hcalRecHitHandle;
		    //iEvent.getByLabel(hcalrechitLabel_, hcalRecHitHandle);
		    //const HBHERecHitCollection *hbhe =  hcalRecHitHandle.product();

		    for(HBHERecHitCollection::const_iterator hh_uc = hbhe->begin(); hh_uc != hbhe->end() && ucHERecHit_subset_n<10000; hh_uc++){
			HcalDetId id(hh_uc->detid());
			if (id.subdet()==2){

			    const CaloCellGeometry *hbhe_cell = caloGeom->getGeometry(hh_uc->id());
			    Global3DPoint hbhe_position = hbhe_cell->getPosition();

			    if(fabs(deltaPhi(ucpho_sc_phi[x_uc],correct_phi(hbhe_position.phi())) ) < 0.5
				    && hh_uc->energy()>1.){
				//find the detid in the set
				set<DetId>::const_iterator ucHERecHitChecker = ucHERecHitSet.find(hh_uc->detid());
				//if detid is not found in the set,(we reached the end), save info!
				if(ucHERecHitChecker == ucHERecHitSet.end()){
				    ucHERecHitSet.insert(hh_uc->detid());
				    ucHERecHit_subset_detid[ucHERecHit_subset_n]  = hh_uc->detid();
				    ucHERecHit_subset_energy[ucHERecHit_subset_n] = hh_uc->energy();
				    ucHERecHit_subset_time[ucHERecHit_subset_n]   = hh_uc->time();
				    ucHERecHit_subset_depth[ucHERecHit_subset_n]  = id.depth();
				    ucHERecHit_subset_phi[ucHERecHit_subset_n]    = correct_phi(hbhe_position.phi());
				    ucHERecHit_subset_eta[ucHERecHit_subset_n]	  = hbhe_position.eta();
				    ucHERecHit_subset_x[ucHERecHit_subset_n]      = hbhe_position.x();
				    ucHERecHit_subset_y[ucHERecHit_subset_n]      = hbhe_position.y();
				    ucHERecHit_subset_z[ucHERecHit_subset_n]      = hbhe_position.z();
				    ucHERecHit_subset_n++;
				}//check to see if hit is already saved
			    }//if delta dphi from photon is small and E>1 try to save
			}
		    }
		}// runHErechit_ && ucpho_isEB


		ucPhoton_n++;
	    }//end of for loop over x_uc
	}//if(ucphoton_container.size!=0) 


	//to get the photon hit information from every crystal of SC
	if(runrechit_){ 
	    Handle<EcalRecHitCollection> Brechit;  //barrel
	    Handle<EcalRecHitCollection> Erechit;  //endcap
	    iEvent.getByLabel(rechitBLabel_,Brechit);
	    iEvent.getByLabel(rechitELabel_,Erechit);

	    edm::ESHandle<EcalSeverityLevelAlgo> sevlv_uc;
	    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv_uc);
	    const EcalSeverityLevelAlgo* sevLevel_uc = sevlv_uc.product();


	    //this will be needed later for swiss corss
	    EcalClusterLazyTools lazyTool(iEvent, iSetup,rechitBLabel_, rechitELabel_ );

	    const EcalRecHitCollection* barrelRecHits_uc= Brechit.product();
	    const EcalRecHitCollection* endcapRecHits_uc= Erechit.product();


	    edm::ESHandle<CaloTopology> pTopology;
	    iSetup.get<CaloTopologyRecord>().get(theCaloTopo_);
	    const CaloTopology *topology = theCaloTopo_.product();

	    if(ucphoton_container.size()!=0){

		for(unsigned int x_uc=0; x_uc < ucphoton_container.size();x_uc++){   

		    std::vector< std::pair<DetId, float> >  PhotonHit_DetIds  = ucphoton_container[x_uc].superCluster()->hitsAndFractions();
		    std::vector<CrystalInfo> uccrystalinfo_container;
		    uccrystalinfo_container.clear();
		    CrystalInfo crystal_uc;
		    float uc_timing_avg =0.0;
		    int uc_ncrys   = 0;
		    uc_ncrysPhoton[x_uc]= 0;
		    vector< std::pair<DetId, float> >::const_iterator detitr_uc;
		    for(detitr_uc = PhotonHit_DetIds.begin(); detitr_uc != PhotonHit_DetIds.end(); ++detitr_uc){
			if (((*detitr_uc).first).det() == DetId::Ecal && ((*detitr_uc).first).subdetId() == EcalBarrel) {
			    EcalRecHitCollection::const_iterator j_uc= Brechit->find(((*detitr_uc).first));
			    EcalRecHitCollection::const_iterator thishit;
			    if ( j_uc!= Brechit->end())  thishit = j_uc;
			    if ( j_uc== Brechit->end()){
				continue;
			    }
			    EBDetId detId  = (EBDetId)((*detitr_uc).first);
			    crystal_uc.rawId  = thishit->id().rawId();
			    crystal_uc.energy = thishit->energy();
			    crystal_uc.time   = thishit->time();
			    crystal_uc.timeErr= thishit->timeError();
			    crystal_uc.recoFlag = thishit->recoFlag(); 
			    crystal_uc.ieta   = detId.ieta();
			    crystal_uc.iphi   = detId.iphi();
			    if(crystal_uc.energy > 0.1){
				uc_timing_avg  = uc_timing_avg + crystal_uc.time;
				uc_ncrys++;
			    }  
			}//end of if ((*detitr_uc).det() == DetId::Ecal && (*detitr_uc).subdetId() == EcalBarrel)
			else if (((*detitr_uc).first).det() == DetId::Ecal && ((*detitr_uc).first).subdetId() == EcalEndcap){
			    EcalRecHitCollection::const_iterator j_uc= Erechit->find(((*detitr_uc).first));
			    EcalRecHitCollection::const_iterator thishit;
			    if ( j_uc!= Erechit->end())  thishit = j_uc;
			    if ( j_uc== Erechit->end()){
				continue;
			    }
			    EEDetId detId  = (EEDetId)((*detitr_uc).first);
			    crystal_uc.energy  = thishit->energy();
			    crystal_uc.time    = thishit->time();
			    crystal_uc.timeErr = thishit->timeError();
			    crystal_uc.recoFlag = thishit->recoFlag(); 
			    crystal_uc.rawId  = 999;
			    crystal_uc.ieta   = -99;
			    crystal_uc.iphi   = -99;
			    if(crystal_uc.energy > 0.1){
				uc_timing_avg  = uc_timing_avg + crystal_uc.time;
				uc_ncrys++;
			    } 
			}//end of if EcalEndcap
			uccrystalinfo_container.push_back(crystal_uc);  
		    }//End loop over detids

		    std::sort(uccrystalinfo_container.begin(),uccrystalinfo_container.end(),EnergySortCriterium());
		    //Without taking into account uncertainty, this time makes no sense.
		    if (uc_ncrys !=0) uc_timing_avg = uc_timing_avg/(float)uc_ncrys;
		    else uc_timing_avg = -99.;
		    uc_ncrysPhoton[x_uc] = uccrystalinfo_container.size(); 
		    ucpho_timingavg_xtal[x_uc]      = uc_timing_avg;
		    for (unsigned int y_uc =0; y_uc < 100.;y_uc++){
			ucpho_timing_xtal[x_uc][y_uc]         = -99.;
			ucpho_energy_xtal[x_uc][y_uc]         = -99.;
			ucpho_ieta_xtalEB[x_uc][y_uc]         = -99;
			ucpho_iphi_xtalEB[x_uc][y_uc]         = -99;
			ucpho_recoFlag_xtalEB[x_uc][y_uc]     = -99;
			ucpho_timeError_xtal[x_uc][y_uc]      = -99.;
		    }//end of for (unsigned int y_uc =0; y_uc < uccrystalinfo_container.size();y_uc++)
		    for (unsigned int y_uc =0; y_uc < uccrystalinfo_container.size() && y_uc < 100;y_uc++){ 
			ucpho_timing_xtal[x_uc][y_uc]         = uccrystalinfo_container[y_uc].time;
			ucpho_timeError_xtal[x_uc][y_uc]      = uccrystalinfo_container[y_uc].timeErr;
			ucpho_energy_xtal[x_uc][y_uc]         = uccrystalinfo_container[y_uc].energy;
			ucpho_ieta_xtalEB[x_uc][y_uc]           = uccrystalinfo_container[y_uc].ieta;
			ucpho_iphi_xtalEB[x_uc][y_uc]           = uccrystalinfo_container[y_uc].iphi;
			ucpho_recoFlag_xtalEB[x_uc][y_uc]           = uccrystalinfo_container[y_uc].recoFlag;
		    }//end of for (unsigned int y_uc =0; y_uc < uccrystalinfo_container.size();y_uc++
		    const reco::BasicCluster& seedClus = *(ucphoton_container[x_uc].superCluster()->seed());




		    if(ucphoton_container[x_uc].isEB()){
			std::vector<float> showershapes_barrel = EcalClusterTools::roundnessBarrelSuperClusters(*(ucphoton_container[x_uc].superCluster()),*barrelRecHits_uc,0);
			ucpho_roundness[x_uc]    = (float)showershapes_barrel[0];
			ucpho_angle[x_uc]        = (float)showershapes_barrel[1];
			ucpho_s9[x_uc]           = ucpho_energy_xtal[x_uc][0]/ucpho_e3x3[x_uc];

			//-New way to get swiss cross
			const reco::CaloClusterPtr  seed = ucphoton_container[x_uc].superCluster()->seed();
			DetId id = lazyTool.getMaximum(*seed).first;
			float uc_swissCross=-99.;

			const EcalRecHitCollection & rechits_uc = ( ucphoton_container[x_uc].isEB() ? *Brechit : *Brechit);
			EcalRecHitCollection::const_iterator it = rechits_uc.find( id );

			if( it != rechits_uc.end() ){uc_swissCross = EcalTools::swissCross( id, rechits_uc, 0.08, true);}

			ucpho_swissCross[x_uc]=uc_swissCross;
			ucpho_e2e9[x_uc]      = -99.;
			ucpho_e6e2[x_uc]      = -99.;    
			ucpho_e4e1[x_uc]      = -99.;        
			ucpho_e2e9[x_uc]      = GetE2OverE9(id,rechits_uc);
			ucpho_e6e2[x_uc]      = Gete6e2( id, rechits_uc);
			ucpho_e4e1[x_uc]      = e4e1(id, rechits_uc);

			if(debug_ && 1-ucpho_swissCross[x_uc]/ucpho_maxEnergyXtal[x_uc] > 0.95) 
			    cout<<"This photon candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl;

			vector<float> stdCov_uc = EcalClusterTools::covariances(seedClus,&(*barrelRecHits_uc),&(*topology),&(*caloGeom));
			vector<float> crysCov_uc = EcalClusterTools::localCovariances(seedClus,&(*barrelRecHits_uc),&(*topology),
				flagExcluded_,severitieExcluded_,
				sevLevel_uc, 4.7);

			ucpho_SigmaEtaPhi[x_uc]   = sqrt(stdCov_uc[1]);
			ucpho_SigmaIetaIphi[x_uc] = sqrt(crysCov_uc[1]);    
			ucpho_SigmaPhiPhi[x_uc]   = sqrt(stdCov_uc[2]);
			ucpho_SigmaIphiIphi[x_uc] = sqrt(crysCov_uc[2]);

		    }//end of if(ucphoton_container[x_uc].isEB())
		    else{ 
			ucpho_roundness[x_uc]   = -99.;
			ucpho_angle[x_uc]       = -99.;
			ucpho_s9[x_uc]          = ucpho_energy_xtal[x_uc][0]/ucpho_e3x3[x_uc];

			//----New way to get the swiss cross
			const reco::CaloClusterPtr  seed = ucphoton_container[x_uc].superCluster()->seed();
			DetId id = lazyTool.getMaximum(*seed).first;
			float uc_swissCross=-99.;
			const EcalRecHitCollection & rechits_uc = ( ucphoton_container[x_uc].isEB() ? *Erechit : *Erechit);
			EcalRecHitCollection::const_iterator it = rechits_uc.find( id );

			if( it != rechits_uc.end() ) {uc_swissCross = EcalTools::swissCross( id, rechits_uc, 0.08, true);
			}                                    
			ucpho_swissCross[x_uc]= uc_swissCross;
			ucpho_e2e9[x_uc]      = -99.;
			ucpho_e6e2[x_uc]      = -99.;
			ucpho_e4e1[x_uc]      = -99.;

			ucpho_e2e9[x_uc]      = GetE2OverE9(id,rechits_uc); 
			ucpho_e6e2[x_uc]      = Gete6e2( id, rechits_uc);
			ucpho_e4e1[x_uc]      = e4e1(id, rechits_uc);

			if(debug_ && 1-ucpho_swissCross[x_uc]/ucpho_maxEnergyXtal[x_uc] > 0.95) {
			    cout<<"This photon candidate is an ECAL spike identified by Swiss Cross algorithm." << endl;
			    cout<<"This would be weird since there aren't spikes in the endcap of ECAL"<<endl; 
			}

			vector<float> stdCov_uc = EcalClusterTools::covariances(seedClus,&(*endcapRecHits_uc),&(*topology),&(*caloGeom));
			vector<float> crysCov_uc = EcalClusterTools::localCovariances(seedClus,&(*endcapRecHits_uc),&(*topology),
				flagExcluded_,severitieExcluded_,
				sevLevel_uc, 4.7);
			ucpho_SigmaEtaPhi[x_uc]   = sqrt(stdCov_uc[1]);
			ucpho_SigmaIetaIphi[x_uc] = sqrt(crysCov_uc[1]);   
			ucpho_SigmaPhiPhi[x_uc]   = sqrt(stdCov_uc[2]);
			ucpho_SigmaIphiIphi[x_uc] = sqrt(crysCov_uc[2]);

		    }//end of else (if !EB)


		}//end of for loop over x_uc
	    }//if(ucphoton_container.size!=0)


	}//if(runrechit_)
    }//if(runphotons_)

