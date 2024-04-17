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
#include "fastjet/ClusterSequence.hh"
// #include "FastJet3.h"


using namespace std;
// using namespace fastjet;
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
,hEventProtocol(nullptr)
,hTrackProtocol(nullptr)
,hLeadingEta(nullptr)
,hTPCnsigma(nullptr)
,hTOFnsigma(nullptr)
,hITSnsigma(nullptr)
,hTPCnsigmaBins{}
,hTOFnsigmaBins{}
,hDCAxy(nullptr)
,hDCAz(nullptr)
,hDCAzJet(nullptr)
,hFullPt(nullptr)
,hFullConePt(nullptr)
,hFullConePtSlow(nullptr)
,hFastJetPt(nullptr)
,hJetRap(nullptr)
,hTotalJetPt(nullptr)
,hConeBkgPt(nullptr)
,hJetProtonPt(nullptr)
,hJetPionPt(nullptr)
,hNumPartInJet(nullptr)
,hNumPartInFastJet(nullptr)
,hNumJets(nullptr)
,hLeadingPt(nullptr)
,hConeRadius(nullptr)
,hConeRadiusAngular(nullptr)
,hConeRadiusFast(nullptr)
,hDeltaTheta(nullptr)
,hFullEta(nullptr)
,hDeltaY(nullptr)
,fMaximumPt(5.0) // Adjust for serious analysis
,fJetRadius(0.4)
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
,hEventProtocol(nullptr)
,hTrackProtocol(nullptr)
,hLeadingEta(nullptr)
,hTPCnsigma(nullptr)
,hTOFnsigma(nullptr)
,hITSnsigma(nullptr)
,hTPCnsigmaBins{}
,hTOFnsigmaBins{}
,hDCAxy(nullptr)
,hDCAz(nullptr)
,hDCAzJet(nullptr)
,hFullPt(nullptr)
,hFullConePt(nullptr)
,hFullConePtSlow(nullptr)
,hFastJetPt(nullptr)
,hJetRap(nullptr)
,hTotalJetPt(nullptr)
,hConeBkgPt(nullptr)
,hJetProtonPt(nullptr)
,hJetPionPt(nullptr)
,hNumPartInJet(nullptr)
,hNumPartInFastJet(nullptr)
,hNumJets(nullptr)
,hLeadingPt(nullptr)
,hConeRadius(nullptr)
,hConeRadiusAngular(nullptr)
,hConeRadiusFast(nullptr)
,hDeltaTheta(nullptr)
,hFullEta(nullptr)
,hDeltaY(nullptr)
,fMaximumPt(5.0)
,fJetRadius(0.4)
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
    delete hEventProtocol;
    delete hTrackProtocol;
    delete hLeadingEta;
    delete hTPCnsigma;
    delete hTOFnsigma;
    delete hITSnsigma;
    delete hDCAxy;
    delete hDCAz;
    delete hDCAzJet;
    delete hFullPt;
    delete hFullConePt;
    delete hFullConePtSlow;
    delete hFastJetPt;
    delete hJetRap;
    delete hTotalJetPt;
    delete hConeBkgPt;
    delete hJetProtonPt;
    delete hJetPionPt;
    delete hNumPartInJet;
    delete hNumPartInFastJet;
    delete hNumJets;
    delete hLeadingPt;
    delete hConeRadius;
    delete hConeRadiusAngular;
    delete hConeRadiusFast;
    delete hDeltaTheta;
    delete hFullEta;
    delete hDeltaY;
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
    
    for (int i=0; i<33; i++) {
        delete hTPCnsigmaBins[i];
    }
    for (int i=0; i<23; i++) {
        delete hTOFnsigmaBins[i];
    }
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
    hNumberOfEvents = new TH1F ("hNumberOfEvents","hNumberOfEvents", 1,0,1);
    fJetOutput      -> Add(hNumberOfEvents);
    hEventProtocol  = new TH1F ("hEventProtocol","hEventProtocol", 20,0,20);
    hEventProtocol  -> Sumw2();
    fJetOutput      -> Add(hEventProtocol);
    hTrackProtocol  = new TH1F ("hTrackProtocol","hTrackProtocol", 20,0,20);
    hTrackProtocol  -> Sumw2();
    fJetOutput      -> Add(hTrackProtocol);

    // (Pseudo)Rapidity
    hDeltaY       = new TH1F ("hDeltaY", "hDeltaY", 100,0,2); // difference between etas of jet and jet candidate
    hDeltaY       -> Sumw2();
    fJetOutput      -> Add(hDeltaY);
    hFullEta        = new TH1F ("hFullEta", "hFullEta", 201,-10,10); // eta of all particles
    hFullEta        -> Sumw2();
    fJetOutput      -> Add(hFullEta);
    hLeadingEta     = new TH1F ("hLeadingEta","hLeadingEta", 201,-10,10);
    hLeadingEta     -> Sumw2();
    fJetOutput      -> Add(hLeadingEta);
    hJetRap         = new TH1D ("hJetRap", "hJetRap", 201,-2,2);
    hJetRap         -> Sumw2();
    fJetOutput      -> Add(hJetRap);

    // pT
    hFullPt         = new TH1D ("hFullPt", "hFullPt", 1000,0,50); // pT of all particles
    hFullPt         -> Sumw2();
    fJetOutput      -> Add(hFullPt);
    hFullConePt     = new TH1D ("hFullConePt", "hFullConePt", 1000,0,50); // pT of jet + background
    hFullConePt     -> Sumw2();
    fJetOutput      -> Add(hFullConePt);
    hFullConePtSlow = new TH1D ("hFullConePtSlow", "hFullConePtSlow", 1000,0,50); // pT of jet + background
    hFullConePtSlow -> Sumw2();
    fJetOutput      -> Add(hFullConePtSlow);
    hFastJetPt      = new TH1D ("hFastJetPt", "hFastJetPt", 1000,0,50); // pT of jet using FastJet
    hFastJetPt      -> Sumw2();
    fJetOutput      -> Add(hFastJetPt);
    hJetPt          = new TH1D ("hJetPt", "hJetPt", 1000,0,50); // pT of jet using seeded anti-kT
    hJetPt          -> Sumw2();
    fJetOutput      -> Add(hJetPt);
    hConeBkgPt      = new TH1D ("hConeBkgPt", "hConeBkgPt", 1000,0,50); // pT of background
    hConeBkgPt      -> Sumw2();
    fJetOutput      -> Add(hConeBkgPt);
    hJetProtonPt    = new TH1D ("hJetProtonPt", "hJetProtonPt", 1000,0,50); // pT of protons in jet
    hJetProtonPt    -> Sumw2();
    fJetOutput      -> Add(hJetProtonPt);
    hJetPionPt      = new TH1D ("hJetPionPt", "hJetPionPt", 1000,0,50); // pT of pions in jet
    hJetPionPt      -> Sumw2();
    fJetOutput      -> Add(hJetPionPt);
    hLeadingPt      = new TH1D ("hLeadingPt", "hLeadingPt", 1000,0,50); // pT of leading tracks
    hLeadingPt      -> Sumw2();
    fJetOutput      -> Add(hLeadingPt);
    hTotalJetPt     = new TH1D ("hTotalJetPt", "hTotalJetPt", 10000,0,500); // pT of total jet
    hTotalJetPt     -> Sumw2();
    fJetOutput      -> Add(hTotalJetPt);

    // nSigma
    hTPCnsigma      = new TH2F ("hTPCnsigma", "hTPCnsigma", 100,-5,5, 100,-5,5);
    hTPCnsigma      -> Sumw2();
    fJetOutput      -> Add(hTPCnsigma);
    hTOFnsigma      = new TH2F ("hTOFnsigma", "hTOFnsigma", 100,-5,5, 100,-5,5);
    hTOFnsigma      -> Sumw2();
    fJetOutput      -> Add(hTOFnsigma);
    hITSnsigma      = new TH2F ("hITSnsigma", "hITSnsigma", 100,-5,5, 100,-5,5);
    hITSnsigma      -> Sumw2();
    fJetOutput      -> Add(hITSnsigma);
    for (int i=0; i<33; i++) {
        hTPCnsigmaBins[i] = new TH1F (Form("hTPCnsigmaBins[%d]",i),Form("hTPCnsigmaBins[%d]",i),100,-5,5);
        hTPCnsigmaBins[i] -> Sumw2();
        fJetOutput        -> Add(hTPCnsigmaBins[i]);
    }
    for (int i=0; i<23; i++) {
        hTOFnsigmaBins[i] = new TH1F (Form("hTOFnsigmaBins[%d]",i),Form("hTOFnsigmaBins[%d]",i),100,-5,5);
        hTOFnsigmaBins[i] -> Sumw2();
        fJetOutput        -> Add(hTOFnsigmaBins[i]);
    }

    // DCA
    hDCAxy          = new TH1F ("hDCAxy", "hDCAxy", 200,-1,1);
    hDCAxy          -> Sumw2();
    fJetOutput      -> Add(hDCAxy);
    hDCAz           = new TH1F ("hDCAz", "hDCAz", 200,-1,1);
    hDCAz           -> Sumw2();
    fJetOutput      -> Add(hDCAz);
    hDCAzJet        = new TH1F ("hDCAzJet", "hDCAzJet", 200,-1,1);
    hDCAzJet        -> Sumw2();
    fJetOutput      -> Add(hDCAzJet);

    hDeltaTheta     = new TH1F ("hDeltaTheta", "hDeltaTheta", 200,0,2);
    hDeltaTheta     -> Sumw2();
    fJetOutput      -> Add(hDeltaTheta);
    hConeRadius     = new TH1D ("hConeRadius","hConeRadius", 200,0,2); // stable jet axis, i.e. after jet reconstruction
    hConeRadius     -> Sumw2();
    fJetOutput          -> Add(hConeRadius);
    hConeRadiusAngular  = new TH1D ("hConeRadiusAngularStable","hConeRadiusAngularStable", 200,0,2);
    hConeRadiusAngular  -> Sumw2();
    fJetOutput          -> Add(hConeRadiusAngular);
    hConeRadiusFast     = new TH1D ("hConeRadiusFast", "hConeRadiusFast", 200,0,2);
    hConeRadiusFast     -> Sumw2();
    fJetOutput          -> Add(hConeRadiusFast);
    hNumPartInJet       = new TH1I ("hNumPartInJet","hNumPartInJet", 200,0,200); // number of particles per jet
    hNumPartInJet       -> Sumw2();
    fJetOutput          -> Add(hNumPartInJet);
    hNumPartInFastJet   = new TH1I ("hNumPartInFastJet","hNumPartInFastJet", 200,0,200); // number of particles per jet (FastJet)
    hNumPartInFastJet   -> Sumw2();
    fJetOutput          -> Add(hNumPartInFastJet);
    hNumJets            = new TH1I ("hNumJets", "hNumJets", 10, 0, 10); // how many jets FastJet collects per event
    hNumJets            -> Sumw2();
    fJetOutput          -> Add(hNumJets);

    //------------------------------------------------------------------------------------
    //Track Selection
    fESDTrackCuts = new AliESDtrackCuts("fESDTrackCuts");
    fESDTrackCuts -> SetAcceptKinkDaughters(kFALSE);
    fESDTrackCuts -> SetRequireTPCRefit(kTRUE);
    fESDTrackCuts -> SetRequireITSRefit(kTRUE);
    fESDTrackCuts -> SetMinNCrossedRowsTPC(70);
    fESDTrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fESDTrackCuts -> SetMinNClustersITS(2);
    fESDTrackCuts -> SetMaxChi2PerClusterITS(36);
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

    if (fRunData) RunData();

}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::RunData()  {
    // TRACK SELECTION
    AliAODEvent *fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAODEvent) Fatal("AliAnalysisTaskJetFemto::RunData", "No AOD event found, check the event handler.");

    //Initialisation
    // Double_t pTmax(0);
    // Int_t leading_ID;
    // Double_t leadingEta; // |eta|
    // AliAODTrack leading_particle;
    vector<Int_t> particle_ID;
    vector<AliAODTrack> particles;
    vector<fastjet::PseudoJet> jetInput;
    Double_t pTbinningTPC[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
    Double_t pTbinningTOF[] = {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0};

    //Preliminary Particle Selection
    for (Int_t i=0; i<fAODEvent->GetNumberOfTracks(); i++) {
        AliAODTrack *track = static_cast<AliAODTrack*>(fAODEvent->GetTrack(i));
        hTrackProtocol -> Fill(0.5);
        hFullPt         -> Fill(track->Pt());

        Double_t eta = track->Eta();
        hFullEta -> Fill(eta);

        if(TMath::Abs(eta)>0.8) continue;
        hTrackProtocol -> Fill(1.5);

        /* if(track->Pt()>pTmax) { // update leading track
            pTmax = track->Pt();
            leading_particle = *track;
            leading_ID = i;
            leadingEta = TMath::Abs(track->Eta());
            hTrackProtocol -> Fill(2.5);
        } */

        if(!PassedTrackSelection(track)) continue; // check track cuts

        Double_t nSigmaTPC_prot = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kProton);
        Double_t nSigmaTOF_prot = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kProton);
        Double_t nSigmaITS_prot = fPIDResponse -> NumberOfSigmasITS(track, AliPID::kProton);
        for (Int_t i=0; i<33; i++) {
            // if (track->Charge()>0) continue;
            if (track->Pt()<0.1) {
                hTPCnsigmaBins[0] -> Fill(nSigmaTPC_prot);
            }
            else if (track->Pt()<pTbinningTPC[i] && track->Pt()>=pTbinningTPC[i-1]) {
                hTPCnsigmaBins[i] -> Fill(nSigmaTPC_prot);
            }
        }
        for (Int_t i=0; i<23; i++) { //------------------------------------------------------------------------------------------------------
            // if (track->Charge()>0) continue;
            if (track->Pt()<0.5) {
                hTOFnsigmaBins[0] -> Fill(nSigmaTOF_prot);
            } 
            else if (track->Pt()<pTbinningTOF[i] && track->Pt()>=pTbinningTOF[i-1]) {
                hTOFnsigmaBins[i] -> Fill(nSigmaTOF_prot);
            }
        }
        hTPCnsigma -> Fill(track->Pt(),nSigmaTPC_prot);
        hTOFnsigma -> Fill(track->Pt(),nSigmaTOF_prot);
        hITSnsigma -> Fill(track->Pt(),nSigmaITS_prot);

        hTrackProtocol -> Fill(3.5);
        hDCAz          -> Fill(GetDCAtoPrimaryVertex(track,1));
        hDCAxy         -> Fill(GetDCAtoPrimaryVertex(track,0));

        particle_ID.emplace_back(i);
        particles.emplace_back(*track);
        jetInput.emplace_back(fastjet::PseudoJet(track->Px(),track->Py(),track->Pz(),track->E()));
    }

    if ((Int_t)particle_ID.size()<2) return;
    hEventProtocol -> Fill(1.5); // how many events have more than 1 particle
    // if (leadingEta>0.43972699396) return; // keep entire cone within |eta|<0.8
    // if (pTmax<fMaximumPt) return; // check max pT
    // hLeadingEta    -> Fill(leading_particle.Eta());
    // hLeadingPt     -> Fill(leading_particle.Pt());

    // JET RECONSTRUCTION
    TLorentzVector jetAxis;
    Double_t jetPt(0);
    Int_t mainJetID;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, fJetRadius);
    fastjet::ClusterSequence cs(jetInput, jet_def);
    vector<fastjet::PseudoJet> jet = sorted_by_pt(cs.inclusive_jets());

    if (jet.size()<1) return;
    hEventProtocol -> Fill(2.5); // how many events have a jet
    hNumJets -> Fill(jet.size());

    for (Int_t i=0; i<jet.size(); i++) {
        if (jet[i].pt()>jetPt) { // select highest pT jet
            Double_t pT = TMath::Sqrt(jet[i].px()*jet[i].px()+jet[i].py()*jet[i].py());
            jetAxis.SetPtEtaPhiE(pT,jet[i].eta(),jet[i].phi(),jet[i].E());
            mainJetID   = i;
        }
        hJetRap     -> Fill(jet[mainJetID].rap());
        hTotalJetPt -> Fill(jet[mainJetID].pt());
    } 

    if (jet[mainJetID].pt()<10) return;
    hEventProtocol -> Fill(3.5); // how many events have a jet with pT > 10 GeV
    vector<fastjet::PseudoJet> constituents = jet[mainJetID].constituents();

    for (Int_t j=0; j<constituents.size(); j++) { // fill jet radius histogram
        hFullConePt       -> Fill(constituents[j].pt());
        Double_t DeltaPhi = TVector2::Phi_0_2pi(constituents[j].phi()-jet[mainJetID].phi());
        Double_t DeltaEta = constituents[j].eta() - jet[mainJetID].eta();
        Double_t Delta    = TMath::Sqrt(DeltaPhi*DeltaPhi-DeltaEta*DeltaEta);
        hConeRadiusFast   -> Fill(Delta);
    }
    
    //Collect Background Particles
    for (Int_t i=0; i<particles.size(); i++) {
        Double_t newPhi   = jetAxis.Phi() + TMath::Pi()/2;
        TLorentzVector UE_Axis;
        UE_Axis.SetPtEtaPhiE(jetAxis.Pt(),jetAxis.Eta(),newPhi,jetAxis.E());
        Double_t DeltaPhi = TVector2::Phi_0_2pi(particles.at(i).Phi() - UE_Axis.Phi());
        Double_t DeltaY   = particles.at(i).Y(); // Delta between phi-pT plane and particle
        Double_t Radius   = TMath::Sqrt(DeltaPhi*DeltaPhi + DeltaY*DeltaY);

        if(Radius>fJetRadius) continue;
        hConeBkgPt->Fill(particles.at(i).Pt());
    }

    //Initialisation for Jet Reconstruction
    /* vector<Int_t> jet_particle_ID;
    vector<AliAODTrack> slowJet;
    jet_particle_ID.emplace_back(leading_ID); // start by adding leading track to particle list

    Int_t exit(0);
    Int_t nPartAssociated(0);

    //4-Momentum of Leading Particle
    TLorentzVector P_leading (0.,0.,0.,0.);
    P_leading.SetPxPyPzE(leading_particle.Px(),leading_particle.Py(),leading_particle.Pz(),leading_particle.E());
    TLorentzVector P_jet = P_leading; */

    //Jet Finder (anti-kT) ----------------------------------------------------------------------------------
    /* do {
        //Initialisation
        Double_t distance_jet_min(1e+08);
        Double_t distance_bkg_min(1e+08);
        AliAODTrack jet_candidate;
        Int_t i_jet_particle(0);

        for (Int_t i=0 ; i<(Int_t)particle_ID.size() ; i++)  { //Loop over all selected Particles
            hTrackProtocol->Fill(4.5);
            // Skip Leading Particle & Elements already associated to the Jet
            if (particle_ID[i]==leading_ID) continue;
            if (particle_ID[i]==-1) continue; // POTENTIAL ERROR --- does deleting particles work -------------------------------------------
            hTrackProtocol->Fill(5.5);

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
                distance_jet_min   = distance_jet;      // update jet distance
                jet_candidate      = particle;           // find closest candidate in for-loop
                i_jet_particle     = i;                 // remember i outside the loop
                hTrackProtocol     -> Fill(6.5);
            }

            //Find Minimum Distance Bkg
            if (distance_bkg<distance_bkg_min)  { // particle's dist to bkg is smaller than min bkg dist?
                distance_bkg_min=distance_bkg;// update min bkg distance
                hTrackProtocol->Fill(7.5);
            }
        } //for (Int_t i=0 ; i<(Int_t)particle_ID.size() ; i++)
        //Looped over selected Particles and received the next Jet Candidate

        if (distance_jet_min<=distance_bkg_min)  {// if distance track-jet smaller than track-bkg, add to jet

            //Add Particle to Jet
            slowJet.emplace_back(jet_candidate);

            //Update 4-Momentum of Jet
            TLorentzVector P_i(0,0,0,0);
            Double_t px_i = jet_candidate.Px();
            Double_t py_i = jet_candidate.Py();
            Double_t pz_i = jet_candidate.Pz();
            Double_t E_i  = jet_candidate.E();
            P_i.SetPxPyPzE(px_i,py_i,pz_i,E_i);
            P_jet = P_jet + P_i;

            //Record Jet Radius
            Double_t DeltaY        = P_jet.Rapidity()-P_i.Rapidity();
            Double_t DeltaPhi      = TVector2::Phi_0_2pi(P_i.Phi()-P_jet.Phi());
            Double_t coneRadiusAng = TMath::Sqrt(DeltaY*DeltaY + DeltaPhi*DeltaPhi);
            Double_t coneRadius    = 1/P_jet.Pt()*coneRadiusAng/fJetRadius;
            hConeRadius         -> Fill(coneRadius);
            hConeRadiusAngular  -> Fill(coneRadiusAng);
            hDeltaY             -> Fill(DeltaY);


            //Remove Element
            particle_ID[i_jet_particle] = -1;
            nPartAssociated++;
        } //if (distance_jet_min<=distance_bkg_min)

        if (nPartAssociated>=((Int_t)particle_ID.size()-1)) exit=1; // gone through all the particles --- in reality this condition should never trigger
        if (distance_jet_min>distance_bkg_min) exit=2; // no jet candidates remaining, dist bkg > dist jet

    } while (exit==0); */

    hFastJetPt -> Add(hFullConePt,1);
    hFastJetPt -> Add(hConeBkgPt,-1);
    // hJetPt -> Add(hFullConePtSlow,1);
    // hJetPt -> Add(hConeBkgPt,-1);
    // hNumPartInJet -> Fill(slowJet.size());

    //Jet Proton Selection
    // ResetGlobalTrackReference();

    /* vector<AliAODTrack> jetprotons;
    vector<AliAODTrack> jetpions;
    for (Int_t i=0; i<(Int_t)slowJet.size();i++) {
        hFullConePtSlow -> Fill(slowJet.at(i).Pt());
        Double_t deltaTheta = TMath::Abs(slowJet.at(i).Theta() - leading_particle.Theta());

        if(IsProton(slowJet.at(i)) && ) {// collect jet protons
            // PROTON CUT INTEGRATION------------------------------------
            hJetProtonPt -> Fill(slowJet.at(i).Pt());
            hDeltaTheta  -> Fill(deltaTheta);
            hDCAzJet     -> Fill(GetDCAtoPrimaryVertex(&slowJet.at(i),1));
            AliAODTrack jetProton = slowJet.at(i);
            jetprotons.emplace_back(jetProton);
            // StoreGlobalTrackReference(&jetProton);
        }
        if(IsPion(slowJet.at(i))) {// collect jet pions
            // PION CUT INTEGRATION--------------------------------------
            hJetPionPt -> Fill(slowJet.at(i).Pt());
            // hDCAzJet     -> Fill(GetDCAtoPrimaryVertex(&slowJet.at(i),1));
            AliAODTrack jetPion = slowJet.at(i);
            jetpions.emplace_back(jetPion);
            // StoreGlobalTrackReference(&jetPion);
        }
    } //for (Int_t i=0; i<(Int_t)slowJet.size();i++)
    
    if ((Int_t)jetprotons.size()<2) return; */
    hEventProtocol -> Fill(4.5); // how many events have more than 1 proton in a jet

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
    if (fAODEvent->GetNumberOfTracks()==0) hEventProtocol->Fill(18.5); // how many events have no tracks

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
    hEventProtocol -> Fill(0.5);

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
    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(&track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(&track,1));

    //Set of Cuts
    Double_t dcaxy_max = 0.1;
    Double_t dcaz_max = 0.1; // go lower to reduce # of secondaries in jets? or is pT cut on jet tracks already enough?

    if (DCAxy>dcaxy_max)                        return isProton;
    if (DCAz>dcaz_max)                          return isProton;

    //Selection
    if (!(pt<0.7 && TMath::Abs(nsigmaITS)<3.0 && TMath::Abs(nsigmaTPC)<5.0)) return isProton; // tighten probably?
    if (!(pt>0.7 && TMath::Abs(nsigmaTPC)<4.0 && TMath::Abs(nsigmaTOF)<10.0)) return isProton;

    isProton=kTRUE;
    return isProton;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFemto::IsPion(AliAODTrack track) {
    //Initialisation
    Bool_t isPion=(kFALSE);

    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (&track,AliPID::kPion);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (&track,AliPID::kPion);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (&track,AliPID::kPion);
    Double_t pt = track.Pt();

    //Selection
    if (!(pt<0.7 && TMath::Abs(nsigmaITS)<3.0 && TMath::Abs(nsigmaTPC)<3.0)) return isPion; // adjust for pion
    if (!(pt>0.7 && TMath::Abs(nsigmaTPC)<2.0 && TMath::Abs(nsigmaTOF)<2.0)) return isPion;

    isPion=kTRUE;
    return isPion;
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
    Int_t nTPCfind    = track->GetTPCNclsF();
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
    Double_t dcaxy_max = 0.3;
    Double_t dcaz_max = 0.5; // go lower to reduce # of secondaries in jets? or is pT cut on jet tracks already enough?

    if (DCAxy>dcaxy_max)                        return passedTrkSelection;
    if (DCAz>dcaz_max)                          return passedTrkSelection;

    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetFemto::GetDCAtoPrimaryVertex (AliAODTrack *track, Int_t index)  {

    Double_t dca[2]; // 0 = DCAxy, 1 = DCAz
    Double_t covMatrix[3];
    if (!track->PropagateToDCA (fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField(),10000,dca,covMatrix)) return -999;

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
