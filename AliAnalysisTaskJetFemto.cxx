#include "AliAnalysisTaskJetFemto.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliAODTrackSelection.h"
#include "AliESDtrackCuts.h"
#include "TLorentzVector.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "THnSparse.h"
#include "AliVVZERO.h"
#include "TObjArray.h"
// #include "AliAODv0.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>


using namespace std;
ClassImp(AliAnalysisTaskJetFemto)

//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::AliAnalysisTaskJetFemto():
AliAnalysisTaskSE()
// JET RECONSTRUCTION
,fAODEvent(nullptr)
,fPIDResponse(nullptr)
,fAODTrackCuts(nullptr)
,fESDTrackCuts(nullptr)
,fJetOutput(nullptr)
,fQAList(nullptr)
,fRunData(kFALSE)
,hNumberOfEvents(nullptr)
,hNumberOfTracks(nullptr)
,hLeadingIDs(nullptr)
,hLeadingEta(nullptr)
,hTPCnsigma(nullptr)
,hTOFnsigma(nullptr)
,hDCAxy(nullptr)
,hDCAz(nullptr)
,hDCAzJet(nullptr)
,hFullPt(nullptr)
,hJetPt(nullptr)
,hJetProtonPt(nullptr)
,hLeadingPt(nullptr)
,fMaximumPt(1.0) // Adjust for serious analysis
,fJetRadius(0.3)
,fTrack()
// FEMTOSCOPY
,fFemtoOutput()
,fIsMC(false)
,fTrigger(AliVEvent::kINT7)
,fEventCuts()
,fTrackCutsProton()
,fTrackCutsAntiproton()
,fConfig()
,fPairCleaner()
,fPartColl()
,fGTI()
,fTrackBufferSize()
{}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::AliAnalysisTaskJetFemto(const char *name, bool isMC):
AliAnalysisTaskSE(name)
// JET RECONSTRUCTION
,fAODEvent(nullptr)
,fPIDResponse(nullptr)
,fAODTrackCuts(nullptr)
,fESDTrackCuts(nullptr)
,fJetOutput(nullptr)
,fQAList(nullptr)
,fRunData(kFALSE)
,hNumberOfEvents(nullptr)
,hNumberOfTracks(nullptr)
,hLeadingIDs(nullptr)
,hLeadingEta(nullptr)
,hTPCnsigma(nullptr)
,hTOFnsigma(nullptr)
,hDCAxy(nullptr)
,hDCAz(nullptr)
,hDCAzJet(nullptr)
,hFullPt(nullptr)
,hJetPt(nullptr)
,hJetProtonPt(nullptr)
,hLeadingPt(nullptr)
,fMaximumPt(1.0)
,fJetRadius(0.3)
// FEMTOSCOPY
,fFemtoOutput()
,fIsMC(isMC)
,fEvent()
,fTrack()
,fTrigger(AliVEvent::kINT7)
,fEventCuts()
,fTrackCutsAntiproton()
,fTrackCutsProton()
,fConfig()
,fPairCleaner()
,fPartColl()
,fGTI()
,fTrackBufferSize(2000)
{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::~AliAnalysisTaskJetFemto()  {
    // JET RECONSTRUCTION
    fJetOutput->Clear();
    delete fAODEvent;
    delete fPIDResponse;
    delete fAODTrackCuts;
    delete fESDTrackCuts;
    delete fJetOutput;
    delete fQAList;
    delete hNumberOfEvents;
    delete hNumberOfTracks;
    delete hLeadingIDs;
    delete hLeadingEta;
    delete hTPCnsigma;
    delete hTOFnsigma;
    delete hDCAxy;
    delete hDCAz;
    delete hDCAzJet;
    delete hFullPt;
    delete hJetPt;
    delete hJetProtonPt;
    delete hLeadingPt;
    // FEMTOSCOPY
    delete fFemtoOutput;
    delete fEvent;
    delete fTrack;
    delete fEventCuts;
    delete fTrackCutsProton;
    delete fTrackCutsAntiproton;
    delete fConfig;
    delete fPairCleaner;
    delete fPartColl;
    delete fGTI;
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserCreateOutputObjects()  {
    // JET RECONSTRUCTION
    //Create Output List
    fJetOutput  = new TList();
    fQAList     = new TList();
    fJetOutput  -> SetOwner();
    fQAList     -> SetOwner();

    // Histograms ------------------------------------------------------------------------
    // Counters
    hNumberOfEvents = new TH1F ("hNumberOfEvents","hNumberOfEvents", 20,0,20);
    hNumberOfEvents -> Sumw2();
    fJetOutput      -> Add(hNumberOfEvents);
    hNumberOfTracks = new TH1F ("hNumberOfTracks","hNumberOfTracks", 20,0,20);
    hNumberOfTracks -> Sumw2();
    fJetOutput      -> Add(hNumberOfTracks);

    // Leading Track
    hLeadingIDs     = new TH2F ("hLeadingIDs","hLeadingIDs", 66001,-33000,33000, 100,0,10);
    hLeadingIDs     -> Sumw2();
    fJetOutput      -> Add(hLeadingIDs);
    hLeadingEta     = new TH1F ("hLeadingEta","hLeadingEta", 201,-10,10);
    hLeadingEta     -> Sumw2();
    fJetOutput      -> Add(hLeadingEta);

    // eta
    hDeltaEta       = new TH1F ("hDeltaEta", "hDeltaEta", 100,0,2);
    hDeltaEta       -> Sumw2();
    fJetOutput      -> Add(hDeltaEta);
    hFullEta        = new TH1F ("hFullEta", "hFullEta", 201,-10,10);
    hFullEta        -> Sumw2();
    fJetOutput      -> Add(hFullEta);

    // pT
    hFullPt         = new TH1D ("hFullPt", "hFullPt", 00,0,10);
    hFullPt         -> Sumw2();
    fJetOutput      -> Add(hFullPt);
    hJetPt          = new TH1D ("hJetPt", "hJetPt", 300,0,10);
    hJetPt          -> Sumw2();
    fJetOutput      -> Add(hJetPt);
    hJetProtonPt    = new TH1D ("hJetProtonPt", "hJetProtonPt", 300,0,10);
    hJetProtonPt    -> Sumw2();
    fJetOutput      -> Add(hJetProtonPt);
    hLeadingPt      = new TH1D ("hLeadingPt", "hLeadingPt", 300,0,10);
    hLeadingPt      -> Sumw2();
    fJetOutput      -> Add(hLeadingPt);

    // nSigma
    hTPCnsigma      = new TH2F ("hTPCnsigma", "hTPCnsigma", 100,0,10, 500,0,50);
    hTPCnsigma      -> Sumw2();
    fJetOutput      -> Add(hTPCnsigma);
    hTOFnsigma      = new TH2F ("hTOFnsigma", "hTOFnsigma", 100,0,10, 500,0,50);
    hTOFnsigma      -> Sumw2();
    fJetOutput      -> Add(hTOFnsigma);

    // DCA
    hDCAxy          = new TH1F ("hDCAxy", "hDCAxy", 300,0,1);
    hDCAxy          -> Sumw2();
    fJetOutput      -> Add(hDCAxy);
    hDCAz           = new TH1F ("hDCAz", "hDCAz", 300,0,1);
    hDCAz           -> Sumw2();
    fJetOutput      -> Add(hDCAz);
    hDCAzJet        = new TH1F ("hDCAzJet", "hDCAzJet", 300,0,1);
    hDCAzJet        -> Sumw2();
    fJetOutput      -> Add(hDCAzJet);

    hDeltaTheta     = new TH1F ("hDeltaTheta", "hDeltaTheta", 200,0,2);
    hDeltaTheta     -> Sumw2();
    fJetOutput      -> Add(hDeltaTheta);

    //------------------------------------------------------------------------------------
    //Track Selection
    fESDTrackCuts = new AliESDtrackCuts("fESDTrackCuts");
    fESDTrackCuts -> SetAcceptKinkDaughters(kFALSE);
    fESDTrackCuts -> SetRequireTPCRefit(kTRUE);
    fESDTrackCuts -> SetRequireITSRefit(kTRUE);
    fESDTrackCuts -> SetMinNCrossedRowsTPC(70);
    fESDTrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fESDTrackCuts -> SetMinNClustersITS(2);
    fESDTrackCuts -> SetMaxChi2PerClusterITS(4);
    fESDTrackCuts -> SetMaxChi2PerClusterTPC(4);

    fAODTrackCuts = new AliAODTrackSelection(fESDTrackCuts, 1);

    // FEMTOSCOPY
    /* fFemtoOutput = new TList();
    fFemtoOutput->SetName("Output");
    fFemtoOutput->SetOwner();

    fTrack=new AliFemtoDreamTrack();
    fTrack->SetUseMCInfo(fIsMC);
    if (!fEventCuts) {
      AliFatal("Event Cuts not set!");
    }
    fEventCuts->InitQA();
    fFemtoOutput->Add(fEventCuts->GetHistList());

    float DummyFloat = fEventCuts->GetlowerPtBoundSpherCalc();
    fEvent=new AliFemtoDreamEvent(false,true,fTrigger,true,DummyFloat);
    fFemtoOutput->Add(fEvent->GetEvtCutList());

    if (!fTrackCutsProton) {
        AliFatal("Track Cuts for Proton not set!");
    }
    fTrackCutsProton->Init();
    fTrackCutsProton->SetName("Proton");
    fFemtoOutput->Add(fTrackCutsProton->GetQAHists());
    if (fTrackCutsProton->GetIsMonteCarlo()) {
        fTrackCutsProton->SetMCName("MCProton");
        fFemtoOutput->Add(fTrackCutsProton->GetMCQAHists());
    }
    if (!fTrackCutsAntiproton) {
        AliFatal("Track Cuts for Antiproton not set!");
    }
    fTrackCutsAntiproton->Init();
    fTrackCutsAntiproton->SetName("Antiproton");
    fFemtoOutput->Add(fTrackCutsAntiproton->GetQAHists());
    if (fTrackCutsAntiproton->GetIsMonteCarlo()) {
        fTrackCutsAntiproton->SetMCName("MCAntiproton");
        fFemtoOutput->Add(fTrackCutsAntiproton->GetMCQAHists());
    }

    fPairCleaner=new AliFemtoDreamPairCleaner(0,0,fConfig->GetMinimalBookingME());
    fFemtoOutput->Add(fPairCleaner->GetHistList());

    fPartColl=new AliFemtoDreamPartCollection(fConfig,fConfig->GetMinimalBookingME());
    fFemtoOutput->Add(fPartColl->GetHistList());
    fFemtoOutput->Add(fPartColl->GetQAList()); */

    //Post Data
    PostData(1, fJetOutput);
    PostData(2, fFemtoOutput);
    PostData(3, fQAList);
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserExec(Option_t *)  {

    //Get Input Event
    if (!GetEvent()) return;

    if (fRunData) RunData ();

}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::RunData()  {
    // TRACK SELECTION
    AliAODEvent *fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAODEvent) Fatal("AliAnalysisTaskJetFemto::RunData", "No AOD event found, check the event handler.");

    //Initialisation
    Double_t pTmax(0);
    Int_t leading_ID;
    Double_t leadingEta; // |eta| to be precise
    AliAODTrack leading_particle;
    vector<Int_t> particle_ID;
    vector<AliAODTrack> particles;

    //Preliminary Particle Selection
    for (Int_t i=0; i<fAODEvent->GetNumberOfTracks(); i++) {
        AliAODTrack *track = static_cast<AliAODTrack*>(fAODEvent->GetTrack(i));
        hNumberOfTracks -> Fill(0.5);
        hFullPt         -> Fill(track->Pt());

        Double_t nSigmaTPC_prot = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kProton);
        Double_t nSigmaTOF_prot = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kProton);
        Double_t nSigmaITS_prot = fPIDResponse -> NumberOfSigmasITS(track, AliPID::kProton);
        hTPCnsigma -> Fill(track->Pt(),nSigmaTPC_prot);
        hTOFnsigma -> Fill(track->Pt(),nSigmaTOF_prot);
        Double_t eta = track->Eta();
        hFullEta -> Fill(eta);

        if(TMath::Abs(eta)>0.8) continue; // check eta
        hNumberOfTracks -> Fill(1.5);

        if(track->Pt()>pTmax) { // update leading track
            pTmax = track->Pt();
            leading_particle = *track;
            leading_ID = i;
            leadingEta = TMath::Abs(track->Eta());
            hNumberOfTracks -> Fill(2.5);
        }

        if(!PassedTrackSelection(track)) continue; // check track cuts
        hNumberOfTracks -> Fill(3.5);
        hDCAz -> Fill(TMath::Abs(GetDCAtoPrimaryVertex(track,1)));

        particle_ID.emplace_back(i);
        particles.emplace_back(*track);
    }

    hLeadingIDs->Fill(leading_ID, pTmax);

    if ((Int_t)particle_ID.size()<2) return;
    hLeadingEta -> Fill(leading_particle.Eta());
    if (leadingEta>0.4397) return; // check eta and leading eta
    if (pTmax<fMaximumPt) return; // check max pT
    hNumberOfEvents -> Fill(1.5); // how many events have a jet
    hLeadingPt -> Fill(leading_particle.Pt());

    // JET RECONSTRUCTION
    //Initialisation for Jet Reconstruction
    vector<Int_t> jet_particle_ID;
    vector<AliAODTrack> jet;
    vector<AliAODTrack> jetprotons;
    jet_particle_ID.emplace_back(leading_ID); // start by adding leading track to particle list

    Int_t exit(0);
    Int_t nPartAssociated(0);

    //4-Momentum of Leading Particle
    TLorentzVector P_jet (0.,0.,0.,0.);
    P_jet.SetPxPyPzE(leading_particle.Px(),leading_particle.Py(),leading_particle.Pz(),leading_particle.E());
    Double_t pt_leading = leading_particle.Pt();

    //Jet Finder (anti-kT) ----------------------------------------------------------------------------------
    do {
        //Initialisation
        Double_t distance_jet_min(1e+08);
        Double_t distance_bkg_min(1e+08);
        Int_t label_jet_particle(0);
        AliAODTrack jet_candidate;
        Int_t i_jet_particle(0);

        for (Int_t i=0 ; i<(Int_t)particle_ID.size() ; i++)  { //Loop over all selected Particles
            hNumberOfTracks->Fill(4.5);
            // Skip Leading Particle & Elements already associated to the Jet
            if (particle_ID[i]==leading_ID) continue;
            if (particle_ID[i]==-1) continue;
            hNumberOfTracks->Fill(5.5);

            //Get Particle 4-Momentum
            TLorentzVector P_particle(0,0,0,0);
            AliAODTrack particle = (AliAODTrack) particles.at(i);
            P_particle.SetPxPyPzE(particle.Px(),particle.Py(),particle.Pz(),particle.E());

            //Variables
            Double_t one_over_pt2_part = 1.0/(P_particle.Pt()*P_particle.Pt());
            Double_t one_over_pt2_lead = 1.0/(P_jet.Pt()*P_jet.Pt());
            Double_t deltaY   = P_particle.Rapidity()-P_jet.Rapidity();
            Double_t deltaPhi = TVector2::Phi_0_2pi(P_particle.Phi()-P_jet.Phi());
            Double_t Delta2   = deltaY*deltaY + deltaPhi*deltaPhi;

            //Distances
            Double_t distance_jet = one_over_pt2_lead*Delta2/(fJetRadius*fJetRadius);
            Double_t distance_bkg = one_over_pt2_part;

            //Find Minimum Distance Jet
            if (distance_jet<distance_jet_min)  {// particle's dist to jet is smaller than min dist
                distance_jet_min=distance_jet;      // update jet distance
                label_jet_particle=particle_ID[i];  // remember ID
                jet_candidate = particle;           // find closest candidate in for-loop
                i_jet_particle = i;                 // remember i outside the loop
                hNumberOfTracks->Fill(6.5);
            }

            //Find Minimum Distance Bkg
            if (distance_bkg<distance_bkg_min)  { // particle's dist to bkg is smaller than min bkg dist?
                distance_bkg_min=distance_bkg;// update min bkg distance
                hNumberOfTracks->Fill(7.5);
            }
        } //for (Int_t i=0 ; i<(Int_t)particle_ID.size() ; i++)

        //Looped over selected Particles and received the next Jet Candidate----------------------------------------

        if (distance_jet_min<=distance_bkg_min)  {// if distance track-jet smaller than track-bkg, add to jet

            //Add Particle to Jet
            jet_particle_ID.emplace_back(label_jet_particle);
            jet.emplace_back(jet_candidate);

            //Update 4-Momentum of Jet
            TLorentzVector P_i(0,0,0,0);
            Double_t px_i = jet_candidate.Px();
            Double_t py_i = jet_candidate.Py();
            Double_t pz_i = jet_candidate.Pz();
            Double_t E_i  = jet_candidate.E();
            P_i.SetPxPyPzE(px_i,py_i,pz_i,E_i);
            P_jet = P_jet + P_i;

            //Remove Element
            particle_ID[i_jet_particle] = -1;
            nPartAssociated++;
        } //if (distance_jet_min<=distance_bkg_min)

        if (nPartAssociated>=((Int_t)particle_ID.size()-1)) exit=1; // gone through all the particles
        if (distance_jet_min>distance_bkg_min) exit=2; // no jet candidates remaining, dist bkg > dist jet

    } while (exit==0);

    //Jet Proton Selection
    // ResetGlobalTrackReference();
    for (Int_t i=0; i<(Int_t)jet.size();i++) {
        hJetPt->Fill(jet.at(i).Pt());
        Double_t deltaEta = TMath::Abs(jet.at(i).Eta() - leading_particle.Eta());
        Double_t deltaTheta = TMath::Abs(jet.at(i).Theta() - leading_particle.Theta());

        // if(!IsProton(jet.at(i))) continue; // collect jet protons
        // PROTON CUT INTEGRATION------------------------------------
        hJetProtonPt->Fill(jet.at(i).Pt());
        hDeltaEta->Fill(deltaEta);
        hDeltaTheta->Fill(deltaTheta);
        hDCAzJet->Fill(TMath::Abs(GetDCAtoPrimaryVertex(&jet.at(i),1)));
        AliAODTrack jetProton = jet.at(i);
        jetprotons.emplace_back(jetProton);
        // StoreGlobalTrackReference(&jetProton);
    } //for (Int_t i=0; i<(Int_t)jet.size();i++)
    
    if ((Int_t)jetprotons.size()<2) return;
    hNumberOfEvents -> Fill(2.5); // how many events have more than 1 proton in a jet

    // FEMTOSCOPY
    /* fTrack->SetGlobalTrackInfo(fGTI,fTrackBufferSize);
    static std::vector<AliFemtoDreamBasePart> Protons;
    Protons.clear();
    static std::vector<AliFemtoDreamBasePart> Antiprotons;
    Antiprotons.clear();

    for (int iTrack = 0;iTrack<jetprotons.size();++iTrack) {
    AliAODTrack *track=static_cast<AliAODTrack*>(&jetprotons.at(iTrack));
    if (!track) {
        AliFatal("No Standard AOD");
        return;
    }
    fTrack->SetTrack(track);

    if (fTrackCutsProton->isSelected(fTrack)) {
        Protons.push_back(*fTrack);
    }
    if (fTrackCutsAntiproton->isSelected(fTrack)) {
        Antiprotons.push_back(*fTrack);
    }
    }

    fPairCleaner->ResetArray();
    fPairCleaner->StoreParticle(Protons);
    fPairCleaner->StoreParticle(Antiprotons);

    fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),fEvent->GetZVertex(),
                         fEvent->GetRefMult08(),fEvent->GetV0MCentrality()); */

    //Check
    if (fAODEvent->GetNumberOfTracks()==0) hNumberOfEvents->Fill(18.5); // how many events have no tracks

    //Post Output Data
    PostData(1, fJetOutput);
    PostData(2, fFemtoOutput);
    PostData(3, fQAList);
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFemto::GetEvent ()  {
    //Get Input Event
    fAODEvent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODEvent) return kFALSE;
    hNumberOfEvents -> Fill(0.5); // how many events are there

    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    return kTRUE;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFemto::IsProton(AliAODTrack track) {
    //Initialisation
    Bool_t isProton=(kFALSE);

    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (&track,AliPID::kProton);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (&track,AliPID::kProton);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (&track,AliPID::kProton);
    Double_t pt = track.Pt();

    //Selection
    if (!(pt<0.7 && TMath::Abs(nsigmaITS)<3.0 && TMath::Abs(nsigmaTPC)<3.0)) return isProton;
    if (!(pt>0.7 && TMath::Abs(nsigmaTPC)<2.0 && TMath::Abs(nsigmaTOF)<2.0)) return isProton;

    isProton=kTRUE;
    return isProton;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFemto::PassedTrackSelection (AliAODTrack *track)  {

    //Initialisation
    Bool_t passedTrkSelection=(kFALSE);

    //Basic Track Selection
    if (!fAODTrackCuts->IsTrackAccepted(track)) return passedTrkSelection;

    //Fixed Cuts
    if (!track->HasPointOnITSLayer(0))                  return passedTrkSelection;
    if (!track->HasPointOnITSLayer(1))                  return passedTrkSelection;

    //Track Variables
    Int_t nTPCcr      = track->GetTPCCrossedRows();
    Int_t nTPCfind     = track->GetTPCNclsF();
    Int_t nITScls     = track->GetITSNcls();
    Int_t nTPCcls     = track->GetTPCNcls();
    Int_t nTPCclsdEdx = track->GetTPCsignalN();
    Double_t chi2ITS  = track->GetITSchi2();
    Double_t chi2TPC  = track->GetTPCchi2();
    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(track,1));
    Double_t cr_over_findable = (Double_t)nTPCcr/((Double_t)nTPCfind);
    Double_t chi2TPC_NDF      = chi2TPC/((Double_t)nTPCcls) ;

    //Set of Cuts
    Int_t nTPCcr_min = 70;
    Int_t nITScls_min = 2;
    Int_t nTPCclsdEdx_min = 50;
    Double_t chi2ITS_max = 4.;
    Double_t chi2TPC_NDF_max = 4.;
    Double_t dcaxy_max = 0.1;
    Double_t dcaz_max = 1.;
    Double_t cr_over_findable_min = 0.7;

    //Selections
    if (nTPCcr<nTPCcr_min)                     return passedTrkSelection;
    if (nITScls<nITScls_min)                   return passedTrkSelection;
    if (nTPCclsdEdx<nTPCclsdEdx_min)           return passedTrkSelection;
    if (cr_over_findable<cr_over_findable_min)   return passedTrkSelection;
    if (chi2TPC_NDF>chi2TPC_NDF_max)           return passedTrkSelection;
    if (chi2ITS>chi2ITS_max)                   return passedTrkSelection;
    if (DCAxy>dcaxy_max)                       return passedTrkSelection;
    if (DCAz>dcaz_max)                         return passedTrkSelection;

    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetFemto::GetDCAtoPrimaryVertex (AliAODTrack *track, Int_t index)  {

    Double_t dca[2]; // 0 = DCAxy, 1 = DCAz
    Double_t covMatrix[3];
    if (!track->PropagateToDCA (fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField(),10000,dca,covMatrix)) return -999;
    hDCAxy->Fill(dca[0]);

    return dca[index];
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetFemto::GetRapidity (AliAODTrack *track, Double_t mass)  {

    //Initialisation
    Double_t y(-999);

    //Rapidity Calculation
    Double_t p  = track->P();
    Double_t pz = track->Pz();
    Double_t E = TMath::Sqrt(mass*mass + p*p);
    if (E == TMath::Abs(pz)) return -999;
    y = 0.5*TMath::Log((E+pz)/(E-pz));

    return y;
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::ResetGlobalTrackReference(){
    for(UShort_t i=0;i<fTrackBufferSize;i++)
    {
        fGTI[i]=0;
    }
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::StoreGlobalTrackReference(AliAODTrack *track){    
    // Check that ID is positive
    const int trackID = track->GetID();
    if(trackID<0){
        return;
    }
    
    // Check ID is not too big for buffer
    if(trackID>=fTrackBufferSize){
        printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
               ,trackID,fTrackBufferSize);
        return;
    }
    
    // Warn if we overwrite a track
    if(fGTI[trackID])
    {
        // Seems like there are FilterMap 0 tracks
        // that have zero TPCNcls, don't store these!
        if( (!track->GetFilterMap()) && (!track->GetTPCNcls()) ){
            return;
        }
        // Imagine the other way around, the zero map zero clusters track
        // is stored and the good one wants to be added. We ommit the warning
        // and just overwrite the 'bad' track
        if( fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()  ){
            // If we come here, there's a problem
            printf("Warning! global track info already there!");
            printf("         TPCNcls track1 %u track2 %u",
                   (fGTI[trackID])->GetTPCNcls(),track->GetTPCNcls());
            printf("         FilterMap track1 %u track2 %u\n",
                   (fGTI[trackID])->GetFilterMap(),track->GetFilterMap());
        }
    } // Two tracks same id
    
    // // There are tracks with filter bit 0,
    // // do they have TPCNcls stored?
    // if(!track->GetFilterMap()){
    //   printf("Filter map is zero, TPCNcls: %u\n"
    //     ,track->GetTPCNcls());
    // }
    
    // Assign the pointer
    (fGTI[trackID]) = track;
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::Terminate(Option_t *)  {

    fJetOutput = dynamic_cast<TList*> (GetOutputData(1));
    if (!fJetOutput) return;
}
//__________________________________________________________________________________________________________________________________
