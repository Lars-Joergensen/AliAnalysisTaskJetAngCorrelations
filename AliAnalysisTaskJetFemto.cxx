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
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"

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
,hNumberOfJets(nullptr)
,hNumberOfHighPtJets(nullptr)
,hEventProtocol(nullptr)
,hTrackProtocol(nullptr)
,hIDsafe(nullptr)
,hIDfail(nullptr)
,hJetParticleID(nullptr)
,hParticleID(nullptr)
,hUserIndex(nullptr)
,hTPCnsigma(nullptr)
,hTPCnsigmaProton(nullptr)
,hTOFnsigma(nullptr)
,hTOFnsigmaProton(nullptr)
,hITSnsigma(nullptr)
,hTPCnsigmaBins{}
,hTOFnsigmaBins{}
,hDCAxyFullEvent(nullptr)
,hDCAzFullEvent(nullptr)
,hDCAzJetProton(nullptr)
,hDCAzJetProtons{}
,hPtJetProtonDCAz{}
,hTPCnsigmaJetProtonDCAz{}
,hTOFnsigmaJetProtonDCAz{}
,hPtFullEvent(nullptr)
,hPtJetParticle(nullptr)
,hJetRapidity(nullptr)
,hPtTotalJet(nullptr)
,hPtSubtractedJet(nullptr)
,hPtJetProton(nullptr)
,hPtJetPion(nullptr)
,hPtDiff(nullptr)
,hNumberOfParticlesInJet(nullptr)
,hNumberOfJetsInEvent(nullptr)
,hNumberOfTracksBeforeCuts(nullptr)
,hJetConeRadius(nullptr)
,hEtaFullEvent(nullptr)
// ,hDeltaY(nullptr)
,fJetRadius(0.4)
,fMinJetPt(10.0)
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
,hNumberOfJets(nullptr)
,hNumberOfHighPtJets(nullptr)
,hEventProtocol(nullptr)
,hTrackProtocol(nullptr)
,hIDsafe(nullptr)
,hIDfail(nullptr)
,hJetParticleID(nullptr)
,hParticleID(nullptr)
,hUserIndex(nullptr)
,hTPCnsigma(nullptr)
,hTPCnsigmaProton(nullptr)
,hTOFnsigma(nullptr)
,hTOFnsigmaProton(nullptr)
,hITSnsigma(nullptr)
,hTPCnsigmaBins{}
,hTOFnsigmaBins{}
,hDCAxyFullEvent(nullptr)
,hDCAzFullEvent(nullptr)
,hDCAzJetProton(nullptr)
,hDCAzJetProtons{}
,hPtJetProtonDCAz{}
,hTPCnsigmaJetProtonDCAz{}
,hTOFnsigmaJetProtonDCAz{}
,hPtFullEvent(nullptr)
,hPtJetParticle(nullptr)
,hJetRapidity(nullptr)
,hPtTotalJet(nullptr)
,hPtSubtractedJet(nullptr)
,hPtJetProton(nullptr)
,hPtJetPion(nullptr)
,hPtDiff(nullptr)
,hNumberOfParticlesInJet(nullptr)
,hNumberOfJetsInEvent(nullptr)
,hNumberOfTracksBeforeCuts(nullptr)
,hJetConeRadius(nullptr)
,hEtaFullEvent(nullptr)
// ,hDeltaY(nullptr)
,fJetRadius(0.4)
,fMinJetPt(10.0)
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
    delete hNumberOfJets;
    delete hNumberOfHighPtJets;
    delete hEventProtocol;
    delete hTrackProtocol;
    delete hIDsafe;
    delete hIDfail;
    delete hJetParticleID;
    delete hParticleID;
    delete hUserIndex;
    delete hTPCnsigma;
    delete hTPCnsigmaProton;
    delete hTOFnsigma;
    delete hTOFnsigmaProton;
    delete hITSnsigma;
    delete hDCAxyFullEvent;
    delete hDCAzFullEvent;
    delete hDCAzJetProton;
    delete hPtFullEvent;
    delete hPtJetParticle;
    delete hJetRapidity;
    delete hPtTotalJet;
    delete hPtSubtractedJet;
    delete hPtJetProton;
    delete hPtJetPion;
    delete hPtDiff;
    delete hNumberOfParticlesInJet;
    delete hNumberOfJetsInEvent;
    delete hNumberOfTracksBeforeCuts;
    delete hJetConeRadius;
    delete hEtaFullEvent;
    // delete hDeltaY;
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
    for (int i=0; i<10; i++) {
        delete hDCAzJetProtons[i];
        delete hPtJetProtonDCAz[i];
        delete hTPCnsigmaJetProtonDCAz[i];
        delete hTOFnsigmaJetProtonDCAz[i];
    }
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserCreateOutputObjects()  {
    // JET RECONSTRUCTION
    //Create Output List
    fJetOutput   = new TList();
    fFemtoOutput = new TList();
    fQAList      = new TList();
    fJetOutput   -> SetOwner();
    fFemtoOutput -> SetOwner();
    fQAList      -> SetOwner();


    // Histograms ------------------------------------------------------------------------
    // Counters
    hNumberOfEvents = new TH1F ("hNumberOfEvents","hNumberOfEvents", 1,0,1);
    fJetOutput      -> Add(hNumberOfEvents);
    hNumberOfJets   = new TH1F ("hNumberOfJets","hNumberOfJets", 1,0,1);
    fJetOutput      -> Add(hNumberOfJets);
    hNumberOfHighPtJets = new TH1F ("hNumberOfHighPtJets","hNumberOfHighPtJets", 1,0,1);
    fJetOutput      -> Add(hNumberOfHighPtJets);
    hEventProtocol  = new TH1F ("hEventProtocol","hEventProtocol", 20,0,20);
    fJetOutput      -> Add(hEventProtocol);
    hTrackProtocol  = new TH1F ("hTrackProtocol","hTrackProtocol", 20,0,20);
    fJetOutput      -> Add(hTrackProtocol);
    hJetParticleID  = new TH1I ("hJetParticleID","hJetParticleID",400,-200,200);
    fJetOutput      -> Add(hJetParticleID);
    hParticleID     = new TH1I ("hParticleID","hParticleID",500,0,500);
    fJetOutput      -> Add(hParticleID);
    hUserIndex      = new TH1I ("hUserIndex","hUserIndex",500,0,500);
    fJetOutput      -> Add(hUserIndex);
    hIDsafe         = new TH1F ("hIDsafe","hIDsafe", 1,0,1);
    fJetOutput      -> Add(hIDsafe);
    hIDfail         = new TH1F ("hIDfail","hIDfail", 2,0,2);
    fJetOutput      -> Add(hIDfail);


    // (Pseudo)Rapidity
    // hDeltaY       = new TH1F ("hDeltaY", "hDeltaY", 100,0,2); // difference between etas of jet and jet candidate
    // hDeltaY       -> Sumw2();
    // fJetOutput      -> Add(hDeltaY);
    hEtaFullEvent   = new TH1F ("hEtaFullEvent", "hEtaFullEvent", 200,-10,10); // eta of all particles
    hEtaFullEvent   -> Sumw2();
    fJetOutput      -> Add(hEtaFullEvent);
    hJetRapidity    = new TH1D ("hJetRapidity", "hJetRapidity", 200,-2,2);
    hJetRapidity    -> Sumw2();
    fJetOutput      -> Add(hJetRapidity);

    // pT
    hPtFullEvent    = new TH1D ("hPtFullEvent", "hPtFullEvent", 2000,0,100); // pT of all particles
    hPtFullEvent    -> Sumw2();
    fJetOutput      -> Add(hPtFullEvent);
    hPtJetParticle  = new TH1D ("hPtJetParticle", "hPtJetParticle", 2000,0,100); // pT of jet using FastJet
    hPtJetParticle  -> Sumw2();
    fJetOutput      -> Add(hPtJetParticle);
    hPtSubtractedJet = new TH1D ("hPtSubtractedJet","hPtSubtractedJet", 2000,0,100); // pT of background subtracted jet
    hPtSubtractedJet -> Sumw2();
    fJetOutput       -> Add(hPtSubtractedJet);
    hPtJetProton    = new TH1D ("hPtJetProton", "hPtJetProton", 2000,0,100); // pT of protons in jet
    hPtJetProton    -> Sumw2();
    fJetOutput      -> Add(hPtJetProton);
    hPtJetPion      = new TH1D ("hPtJetPion", "hPtJetPion", 2000,0,100); // pT of pions in jet
    hPtJetPion      -> Sumw2();
    fJetOutput      -> Add(hPtJetPion);
    hPtTotalJet     = new TH1D ("hPtTotalJet", "hPtTotalJet", 10000,0,500); // pT of total jet
    hPtTotalJet     -> Sumw2();
    fJetOutput      -> Add(hPtTotalJet);
    hPtDiff         = new TH1D ("hPtDiff","hPtDiff",1000000,-0.0005,0.0005);
    hPtDiff         -> Sumw2();
    fJetOutput      -> Add(hPtDiff);

    // nSigma
    hTPCnsigma       = new TH2F ("hTPCnsigma", "hTPCnsigma", 10000,-50,50, 100,-5,5);
    hTPCnsigma       -> Sumw2();
    fJetOutput       -> Add(hTPCnsigma);
    hTPCnsigmaProton = new TH2F ("hTPCnsigmaProton", "hTPCnsigmaProton", 10000,-50,50, 100,-5,5);
    hTPCnsigmaProton -> Sumw2();
    fJetOutput       -> Add(hTPCnsigmaProton);
    hTOFnsigma       = new TH2F ("hTOFnsigma", "hTOFnsigma", 10000,-50,50, 100,-10,10);
    hTOFnsigma       -> Sumw2();
    fJetOutput       -> Add(hTOFnsigma);
    hTOFnsigmaProton = new TH2F ("hTOFnsigmaProton", "hTOFnsigmaProton", 10000,-50,50, 100,-10,10);
    hTOFnsigmaProton -> Sumw2();
    fJetOutput       -> Add(hTOFnsigmaProton);
    hITSnsigma       = new TH2F ("hITSnsigma", "hITSnsigma", 100,-5,5, 100,-5,5);
    hITSnsigma       -> Sumw2();
    fJetOutput       -> Add(hITSnsigma);
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
    hDCAxyFullEvent          = new TH1F ("hDCAxyFullEvent", "hDCAxyFullEvent", 200,-1,1);
    hDCAxyFullEvent          -> Sumw2();
    fJetOutput      -> Add(hDCAxyFullEvent);
    hDCAzFullEvent           = new TH1F ("hDCAzFullEvent", "hDCAzFullEvent", 200,-2.4,2.4);
    hDCAzFullEvent           -> Sumw2();
    fJetOutput      -> Add(hDCAzFullEvent);
    hDCAzJetProton        = new TH1F ("hDCAzJetProton", "hDCAzJetProton", 200,-2.4,2.4);
    hDCAzJetProton        -> Sumw2();
    fJetOutput      -> Add(hDCAzJetProton);
    for (int i=0; i<10; i++) {
        hDCAzJetProtons[i] = new TH1F (Form("hDCAzJetProtons[%d]",i), Form("hDCAzJetProtons[%d]",i), 200,-0.5,0.5);
        hDCAzJetProtons[i] -> Sumw2();
        fJetOutput      -> Add(hDCAzJetProtons[i]);
        hPtJetProtonDCAz[i] = new TH1D (Form("hPtJetProtonDCAz[%d]",i),Form("hPtJetProtonDCAz[%d]",i), 2000,0,100);
        hPtJetProtonDCAz[i] -> Sumw2();
        fJetOutput          -> Add(hPtJetProtonDCAz[i]);
        hTPCnsigmaJetProtonDCAz[i] = new TH2F (Form("hTPCnsigmaJetProtonDCAz[%d]",i),Form("hTPCnsigmaJetProtonDCAz[%d]",i), 10000,-50,50,100,-5,5);
        hTPCnsigmaJetProtonDCAz[i] -> Sumw2();
        fJetOutput                 -> Add(hTPCnsigmaJetProtonDCAz[i]);
        hTOFnsigmaJetProtonDCAz[i] = new TH2F (Form("hTOFnsigmaJetProtonDCAz[%d]",i),Form("hTOFnsigmaJetProtonDCAz[%d]",i), 10000,-50,50,100,-10,10);
        hTOFnsigmaJetProtonDCAz[i] -> Sumw2();
        fJetOutput                 -> Add(hTOFnsigmaJetProtonDCAz[i]);
    }

    hJetConeRadius            = new TH1D ("hJetConeRadius", "hJetConeRadius", 200,0,2);
    hJetConeRadius            -> Sumw2();
    fJetOutput                -> Add(hJetConeRadius);
    hNumberOfParticlesInJet   = new TH1I ("hNumberOfParticlesInJet","hNumberOfParticlesInJet", 200,0,200); // number of particles per jet
    hNumberOfParticlesInJet   -> Sumw2();
    fJetOutput                -> Add(hNumberOfParticlesInJet);
    hNumberOfJetsInEvent      = new TH1I ("hNumberOfJetsInEvent", "hNumberOfJetsInEvent", 10, 0, 10); // how many jets FastJet collects per event
    hNumberOfJetsInEvent      -> Sumw2();
    fJetOutput                -> Add(hNumberOfJetsInEvent);
    hNumberOfTracksBeforeCuts = new TH1I ("hNumberOfTracksBeforeCuts", "hNumberOfTracksBeforeCuts", 2000,0,2000);
    hNumberOfTracksBeforeCuts -> Sumw2();
    fJetOutput                -> Add(hNumberOfTracksBeforeCuts);

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
    hNumberOfTracksBeforeCuts -> Fill(fAODEvent->GetNumberOfTracks());
    //Initialisation
    vector<Int_t> particle_ID;
    map<int,AliAODTrack> particles;
    vector<fastjet::PseudoJet> jetInput;
    vector<fastjet::PseudoJet> jets;
    vector<fastjet::PseudoJet> constituents;
    int index = 0;
    fastjet::PseudoJet hardestJet(0.,0.,0.,0.);
    fastjet::PseudoJet subtractedJet(0.,0.,0.,0.);
    particles.clear();
    jetInput.clear();
    jets.clear();
    constituents.clear();

    Double_t pTbinningTPC[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
    Double_t pTbinningTOF[] = {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0};

    //Preliminary Particle Selection
    for (Int_t i=0; i<fAODEvent->GetNumberOfTracks(); i++) {
        AliAODTrack *trackPointer = static_cast<AliAODTrack*>(fAODEvent->GetTrack(i));
        if(!trackPointer) continue; // shouldn't be needed but better safe than sorry
        AliAODTrack track = *trackPointer;
        hTrackProtocol -> Fill(0.5); // how many tracks are there before cuts
        Double_t eta = track.Eta();

        if(TMath::Abs(eta)>0.8) continue;
        hTrackProtocol -> Fill(1.5); // how many tracks with |eta| < 0.8 are there

        if(!PassedTrackSelection(track)) continue; // check track cuts
        hTrackProtocol -> Fill(2.5); // how many tracks are there after basic cuts
        hEtaFullEvent -> Fill(eta);
        hPtFullEvent        -> Fill(track.Pt());

        //Full Event nSigma
        Double_t nSigmaTPC_proton = fPIDResponse -> NumberOfSigmasTPC(&track, AliPID::kProton);
        Double_t nSigmaTOF_proton = fPIDResponse -> NumberOfSigmasTOF(&track, AliPID::kProton);
        Double_t nSigmaITS_proton = fPIDResponse -> NumberOfSigmasITS(&track, AliPID::kProton);
        for (Int_t i=0; i<33; i++) {
            if (track.Charge()>0) continue;
            if (track.Pt()<0.1) {
                hTPCnsigmaBins[0] -> Fill(nSigmaTPC_proton);
            }
            else if (track.Pt()<pTbinningTPC[i] && track.Pt()>=pTbinningTPC[i-1]) {
                hTPCnsigmaBins[i] -> Fill(nSigmaTPC_proton);
            }
        }
        for (Int_t i=0; i<23; i++) { //------------------------------------------------------------------------------------------------------
            if (track.Charge()>0) continue;
            if (track.Pt()<0.5) {
                hTOFnsigmaBins[0] -> Fill(nSigmaTOF_proton);
            }
            else if (track.Pt()<pTbinningTOF[i] && track.Pt()>=pTbinningTOF[i-1]) {
                hTOFnsigmaBins[i] -> Fill(nSigmaTOF_proton);
            }
        }
        hTPCnsigma -> Fill(track.Pt(),(nSigmaTPC_proton*track.Charge()));
        hTOFnsigma -> Fill(track.Pt(),(nSigmaTOF_proton*track.Charge()));
        hITSnsigma -> Fill(track.Pt(),(nSigmaITS_proton*track.Charge()));

        hTrackProtocol -> Fill(3.5);
        hDCAzFullEvent          -> Fill(GetDCAtoPrimaryVertex(track,1));
        hDCAxyFullEvent         -> Fill(GetDCAtoPrimaryVertex(track,0));

        hParticleID -> Fill(i);
        particle_ID.emplace_back(index);
        particles[index] = track;
        fastjet::PseudoJet inputPseudoJet(track.Px(),track.Py(),track.Pz(),track.E());
        inputPseudoJet.set_user_index(index);
        jetInput.emplace_back(inputPseudoJet);
        index++;
    }

    if ((Int_t)particles.size()<2) return;
    hEventProtocol -> Fill(1.5); // how many events have more than 1 particle

    // JET RECONSTRUCTION -------------------------------------------------------------------------------------------------------------------
    double ghost_maxrap = 1.0; // change? -> not necessarily
    int ghost_repeat    = 1;
    double ghost_area   = 0.005;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, fJetRadius);
    fastjet::JetDefinition jetDefBkg(fastjet::kt_algorithm, 0.5); // maybe increase to R=0.6
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(ghost_maxrap,ghost_repeat,ghost_area));
    fastjet::AreaDefinition areaDefBkg(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghost_maxrap));
    fastjet::ClusterSequenceArea clusterSeq(jetInput, jetDef, areaDef);
    jets = sorted_by_pt(clusterSeq.inclusive_jets());
    if (jets.size()==0) return;
    hEventProtocol -> Fill(2.5); // how many events have a jet
    hNumberOfJetsInEvent -> Fill(jets.size());
    for (int i=0; i<jets.size(); i++) {
        hNumberOfJets -> Fill(0.5);
    }

    hardestJet = jets[0]; // jets sorted by pT from highest to lowest
    hJetRapidity     -> Fill(hardestJet.rap());
    hPtTotalJet -> Fill(hardestJet.pt());

    if (hardestJet.constituents().size()<2) return; // no point if the jet only has 1 particle
    if (hardestJet.pt()<fMinJetPt) return;
    hNumberOfHighPtJets -> Fill(0.5);
    hEventProtocol -> Fill(3.5); // how many events have a jet with pT > 10 GeV
    constituents = hardestJet.constituents();

    for (Int_t j=0; j<constituents.size(); j++) { // fill jet radius histogram
        hPtJetParticle            -> Fill(constituents[j].pt());
        Double_t DeltaPhi = TVector2::Phi_0_2pi(constituents[j].phi()-hardestJet.phi());
        Double_t DeltaEta = constituents[j].eta() - hardestJet.eta();
        Double_t Delta    = TMath::Sqrt(DeltaPhi*DeltaPhi-DeltaEta*DeltaEta);
        hJetConeRadius       -> Fill(Delta);
    }

    //Background Subtraction
    fastjet::Selector selector = fastjet::SelectorAbsEtaMax(1.0) * (!fastjet::SelectorNHardest(2));
    fastjet::JetMedianBackgroundEstimator bkgEst(selector, jetDefBkg, areaDefBkg);
    fastjet::Subtractor subtractor(&bkgEst);
    bkgEst.set_particles(jetInput);

    /* for (auto& jet : jets) { // loops over jets, makes ref variable jet
        if (!clusterSeq.is_pure_ghost(jet)) {
        totalJetAreaPhys += jet.area();
        }
        totalAreaCovered += jet.area();
    }
    double occupancyFactor = totalAreaCovered > 0 ? totalJetAreaPhys / totalAreaCovered : 1.;
        rho *= occupancyFactor;
        rhoM *= occupancyFactor; */

    subtractedJet = subtractor(hardestJet);
    if (!subtractedJet.has_constituents()) return;
    for (int i=0; i<subtractedJet.constituents().size(); i++) {
        hPtSubtractedJet -> Fill(subtractedJet.constituents()[i].pt());
    }

    //Jet Proton Selection
    // ResetGlobalTrackReference();

    vector<AliAODTrack> jetProtons;
    // vector<AliAODTrack> jetPions;
    for (Int_t i=0; i<(Int_t)constituents.size();i++) {
        fastjet::PseudoJet pseudoParticle = constituents[i];
        int id = pseudoParticle.user_index();
        hUserIndex -> Fill(pseudoParticle.user_index());

        hTrackProtocol -> Fill(4.5);
        if (id < 1) continue;
        hTrackProtocol -> Fill(5.5);
        if (id % 1 != 0) continue;
        hTrackProtocol -> Fill(6.5);
        hJetParticleID -> Fill(id);
        double jetParticlePt = 0;
        // try {
            AliAODTrack jetParticle = particles[id];
            hIDsafe -> Fill(0.5);
            jetParticlePt = (double)jetParticle.Pt();
        // }
        /* catch (const std::out_of_range& e) {
            hIDfail -> Fill(0.5);
            if (particles.size()-id < 0) hIDfail -> Fill(1.5);
        } */
        double diff = (double)pseudoParticle.pt() - jetParticlePt;
        hPtDiff -> Fill(diff);
        if(IsProton(jetParticle)) {// collect protons in jet
            // PROTON CUT INTEGRATION------------------------------------
            // if (jetParticle.Pt()<2) continue;
            Double_t nSigmaTPC_proton = fPIDResponse -> NumberOfSigmasTPC(&jetParticle, AliPID::kProton);
            Double_t nSigmaTOF_proton = fPIDResponse -> NumberOfSigmasTOF(&jetParticle, AliPID::kProton);
            hTPCnsigma     -> Fill(jetParticle.Pt(),(nSigmaTPC_proton*jetParticle.Charge()));
            hTOFnsigma     -> Fill(jetParticle.Pt(),(nSigmaTOF_proton*jetParticle.Charge()));
            hPtJetProton   -> Fill(jetParticle.Pt());
            hDCAzJetProton -> Fill(GetDCAtoPrimaryVertex(jetParticle,1));
            jetProtons.emplace_back(jetParticle);
            // StoreGlobalTrackReference(&jetParticle);
        }
        // if(IsPion(jetParticle)) {// collect pions in jet
            // PION CUT INTEGRATION--------------------------------------
            // hPtJetPion -> Fill(jetParticlePt);
            // hDCAzJetPion -> Fill(GetDCAtoPrimaryVertex(jetParticle,1));
            // jetPions.emplace_back(jetParticle);
        // }
    } //for (Int_t i=0; i<(Int_t)constituents.size();i++)

    if ((Int_t)jetProtons.size()<2) return;
    hEventProtocol -> Fill(4.5); // how many events have more than 1 proton in a jet

    /* for (Int_t i=0; i<jetProtons.size(); i++) {
        for (Int_t j=0; j<10; j++) {
            Double_t DCAz = GetDCAtoPrimaryVertex(jetProtons[i],1);
            if (DCAz > (0.5 - j*0.05)) continue;
            Double_t nSigmaTPC_proton_DCAz = fPIDResponse -> NumberOfSigmasTPC(&jetProtons[i], AliPID::kProton);
            Double_t nSigmaTOF_proton_DCAz = fPIDResponse -> NumberOfSigmasTOF(&jetProtons[i], AliPID::kProton);
            hDCAzJetProtons[j] -> Fill(DCAz);
            hPtJetProtonDCAz[j] -> Fill(jetProtons[i].Pt());
            hTPCnsigmaJetProtonDCAz[j] -> Fill(jetProtons[i].Pt(),nSigmaTPC_proton_DCAz);
            hTOFnsigmaJetProtonDCAz[j] -> Fill(jetProtons[i].Pt(),nSigmaTOF_proton_DCAz);
        }
    } */

    // FEMTOSCOPY
    // fTrack->SetGlobalTrackInfo(fGTI,fTrackBufferSize);
    /* static std::vector<AliFemtoDreamBasePart> Protons;
    Protons.clear();
    static std::vector<AliFemtoDreamBasePart> Antiprotons;
    Antiprotons.clear();

    for (int iTrack = 0;iTrack<jetProtons.size();++iTrack) {
    AliAODTrack *track=static_cast<AliAODTrack*>(&jetProtons.at(iTrack));
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
    PostData(3, fQAList); // do these still get executed even if return was triggered earlier? logically no, but how then do we post data?
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
    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(track,1));

    //Set of Cuts
    Double_t dcaxy_max = 0.1;
    Double_t dcaz_max = 2.4; // go lower to reduce # of secondaries in jets? or is pT cut on jet tracks already enough?

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
Bool_t AliAnalysisTaskJetFemto::PassedTrackSelection (AliAODTrack track)  {

    //Initialisation
    Bool_t passedTrkSelection=(kFALSE);

    //Basic Track Selection
    if (!fAODTrackCuts->IsTrackAccepted(&track)) return passedTrkSelection;

    //Fixed Cuts
    if (!track.HasPointOnITSLayer(0))                  return passedTrkSelection;
    if (!track.HasPointOnITSLayer(1))                  return passedTrkSelection;

    //Track Variables
    // Int_t nTPCcr      = track->GetTPCCrossedRows();
    // Int_t nTPCfind    = track->GetTPCNclsF();
    // Int_t nITScls     = track->GetITSNcls();
    // Int_t nTPCcls     = track->GetTPCNcls();
    // Int_t nTPCclsdEdx = track->GetTPCsignalN();
    // Double_t chi2ITS  = track->GetITSchi2();
    // Double_t chi2TPC  = track->GetTPCchi2();
    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(track,1));
    // Double_t cr_over_findable = (Double_t)nTPCcr/((Double_t)nTPCfind);
    // Double_t chi2TPC_NDF      = chi2TPC/((Double_t)nTPCcls) ;

    //Set of Cuts
    Double_t dcaxy_max = 0.3;
    Double_t dcaz_max = 2.4; // go lower to reduce # of secondaries in jets? or is pT cut on jet tracks already enough?

    if (DCAxy>dcaxy_max)                        return passedTrkSelection;
    if (DCAz>dcaz_max)                          return passedTrkSelection;

    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetFemto::GetDCAtoPrimaryVertex (AliAODTrack track, Int_t index)  {

    Double_t dca[2]; // 0 = DCAxy, 1 = DCAz
    Double_t covMatrix[3];
    if (!track.PropagateToDCA (fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField(),10000,dca,covMatrix)) return -999;

    return dca[index];
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetFemto::GetRapidity (AliAODTrack track, Double_t mass)  {

    //Initialisation
    Double_t y(-999);

    //Rapidity Calculation
    Double_t p  = track.P();
    Double_t pz = track.Pz();
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
    if (trackID<0) return;

    // Check ID is not too big for buffer
    if (trackID>=fTrackBufferSize) {
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
        // is stored and the good one wants to be added. We omit the warning
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
