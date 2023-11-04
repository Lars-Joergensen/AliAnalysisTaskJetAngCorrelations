#include "AliAnalysisTaskJetFemto.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "TLorentzVector.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "THnSparse.h"
#include "AliVVZERO.h"
#include "TObjArray.h"
#include "AliESDv0.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
ClassImp(AliAnalysisTaskJetFemto)

//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::AliAnalysisTaskJetFemto():
AliAnalysisTaskSE(),
fESDEvent(nullptr),
fPIDResponse(nullptr),
fESDtrackCuts{},
fESDtrackCuts_V0daugh(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fIsITSrecalib(kFALSE),
fRunData(kFALSE),
fMatchingEff(kFALSE),
fESDEventSelection(),
hMean(nullptr),
hWidth(nullptr),
hNumberOfEvents(nullptr),
hMultiplicity(nullptr),
hMultDistribution(nullptr),
hTPCnsigma(nullptr),
hTOFnsigma(nullptr),
hTPCnsigma_vs_rap(nullptr),
hTOFnsigma_vs_rap(nullptr),
hITSnsigma(nullptr),
hAntiprotonsTPC(nullptr),
// hAntiprotonsTOF(nullptr),
hnSigmaProtons_vs_Pt(nullptr)
{}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::AliAnalysisTaskJetFemto(const char *name):
AliAnalysisTaskSE(name),
fESDEvent(nullptr),
fPIDResponse(nullptr),
fESDtrackCuts{},
fESDtrackCuts_V0daugh(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fIsITSrecalib(kFALSE),
fRunData(kFALSE),
fMatchingEff(kFALSE),
fESDEventSelection(),
hMean(nullptr),
hWidth(nullptr),
hNumberOfEvents(nullptr),
hMultiplicity(nullptr),
hMultDistribution(nullptr),
hTPCnsigma(nullptr),
hTOFnsigma(nullptr),
hTPCnsigma_vs_rap(nullptr),
hTOFnsigma_vs_rap(nullptr),
hITSnsigma(nullptr),
hAntiprotonsTPC(nullptr),
// hAntiprotonsTOF(nullptr),
hnSigmaProtons_vs_Pt(nullptr)
{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::~AliAnalysisTaskJetFemto()  {

    fOutputList->Clear();
    delete fESDEvent;
    delete fPIDResponse;
    for (Int_t isyst=0 ; isyst<50 ; isyst++) { delete fESDtrackCuts[isyst]; }

    delete fESDtrackCuts_V0daugh;
    delete fOutputList;
    delete fQAList;
    delete hNumberOfEvents;
    delete hMultiplicity;
    delete hMultDistribution;
    delete hTPCnsigma;
    delete hTOFnsigma;
    delete hTPCnsigma_vs_rap;
    delete hTOFnsigma_vs_rap;
    delete hITSnsigma;
    delete hAntiprotonsTPC;
    // delete hAntiprotonsTOF;
    delete hnSigmaProtons_vs_Pt;
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserCreateOutputObjects()  {

    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();

    //Event Selection
    fESDEventSelection.AddQAplotsToList(fQAList);
    fESDEventSelection.SetManualMode();
    fESDEventSelection.fRequireTrackVertex = true;
    fESDEventSelection.fMinVtz = -10.f;
    fESDEventSelection.fMaxVtz = 10.f;
    fESDEventSelection.fMaxDeltaSpdTrackAbsolute = 0.5f;
    fESDEventSelection.fMaxResolutionSPDvertex = 0.25f;
    fESDEventSelection.fTriggerMask = (AliVEvent::kINT7);
    fESDEventSelection.fRejectDAQincomplete = true;
    fESDEventSelection.fSPDpileupMinContributors = 3;
    fESDEventSelection.fSPDpileupMinZdist = 0.8;
    fESDEventSelection.fSPDpileupNsigmaZdist = 3.;
    fESDEventSelection.fSPDpileupNsigmaDiamXY = 2.;
    fESDEventSelection.fSPDpileupNsigmaDiamZ = 5.;
    fESDEventSelection.fTrackletBGcut = true;


    //Event Counter and Centrality Distribution
    hNumberOfEvents = new TH1F ("hNumberOfEvents","",20,0,20);
    hMultiplicity   = new TH1F ("hMultiplicity","",100,0,100);
    hNumberOfEvents -> Sumw2();
    hMultiplicity   -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    fOutputList -> Add(hMultiplicity);

    hMultDistribution         = new TH1I ("hMultDistribution","",200,0,200);
    hMultDistribution -> Sumw2();
    fOutputList -> Add(hMultDistribution);

    //Binning nsigma_{TPC}
    Int_t    bins_nsigmaTPC[4]  = { 50,   7,  20,   100 };
    Double_t xmin_nsigmaTPC[4]  = {  0, 0.0, 0.0, -10.0 };
    Double_t xmax_nsigmaTPC[4]  = { 50, 0.7, 1.0,  10.0 };

    //Binning nsigma_{TOF}
    Int_t    bins_nsigmaTOF[4]  = { 50,   7,  30,   200 };
    Double_t xmin_nsigmaTOF[4]  = {  0, 0.0, 0.5, -20.0 };
    Double_t xmax_nsigmaTOF[4]  = { 50, 0.7, 2.0,  20.0 };

    //nsigma_{TPC}
    hTPCnsigma = new THnSparseF ("hTPCnsigma","",4, bins_nsigmaTPC, xmin_nsigmaTPC, xmax_nsigmaTPC);
    hTPCnsigma -> Sumw2();
    fOutputList -> Add (hTPCnsigma);

    //nsigma_{TOF}
    hTOFnsigma = new THnSparseF ("hTOFnsigma","",4, bins_nsigmaTOF, xmin_nsigmaTOF, xmax_nsigmaTOF);
    hTOFnsigma -> Sumw2();
    fOutputList -> Add (hTOFnsigma);


    //**************************************** Comparison y>0 vs. y<0 ****************************************

    //Binning nsigma_{TPC}
    Int_t    bins_nsigmaTPC_vs_y[3]  = {   14,  20,   100 };
    Double_t xmin_nsigmaTPC_vs_y[3]  = { -0.7, 0.0, -10.0 };
    Double_t xmax_nsigmaTPC_vs_y[3]  = {  0.7, 1.0,  10.0 };

    //Binning nsigma_{TOF}
    Int_t    bins_nsigmaTOF_vs_y[3]  = {   14,  30,   200 };
    Double_t xmin_nsigmaTOF_vs_y[3]  = { -0.7, 0.5, -20.0 };
    Double_t xmax_nsigmaTOF_vs_y[3]  = {  0.7, 2.0,  20.0 };

    //nsigma_{TPC}
    hTPCnsigma_vs_rap = new THnSparseF ("hTPCnsigma_vs_rap","",3, bins_nsigmaTPC_vs_y, xmin_nsigmaTPC_vs_y, xmax_nsigmaTPC_vs_y);
    hTPCnsigma_vs_rap -> Sumw2();
    fOutputList -> Add (hTPCnsigma_vs_rap);

    //nsigma_{TOF}
    hTOFnsigma_vs_rap = new THnSparseF ("hTOFnsigma_vs_rap","",3, bins_nsigmaTOF_vs_y,xmin_nsigmaTOF_vs_y, xmax_nsigmaTOF_vs_y);
    hTOFnsigma_vs_rap -> Sumw2();
    fOutputList -> Add (hTOFnsigma_vs_rap);


    //2D ITS Re-calibration Map
    if (fIsITSrecalib)  {
        hITSnsigma = new TH3F ("hITSnsigma","",20,-1.0,1.0,20,0.0,1.0,200,-5,5);//eta,p,nsigma_{ITS}
        hITSnsigma  -> Sumw2();
        fOutputList -> Add (hITSnsigma);
    }

    //Track Selection
    for (Int_t isyst=0 ; isyst<50 ; isyst++)  {
        fESDtrackCuts[isyst] = new AliESDtrackCuts(Form("fESDtrackCuts[%d]",isyst));
        fESDtrackCuts[isyst] -> SetAcceptKinkDaughters(kFALSE);
        fESDtrackCuts[isyst] -> SetRequireTPCRefit(kTRUE);
        fESDtrackCuts[isyst] -> SetRequireITSRefit(kTRUE);
        fESDtrackCuts[isyst] -> SetMinNCrossedRowsTPC(50);
        fESDtrackCuts[isyst] -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
        fESDtrackCuts[isyst] -> SetMinNClustersITS(1);
        fESDtrackCuts[isyst] -> SetMaxChi2PerClusterITS(36);
        fESDtrackCuts[isyst] -> SetMaxChi2PerClusterTPC(10);
    }

    fESDtrackCuts_V0daugh = new AliESDtrackCuts("fESDtrackCuts_V0daugh");
    fESDtrackCuts_V0daugh -> SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts_V0daugh -> SetRequireTPCRefit(kTRUE);
    fESDtrackCuts_V0daugh -> SetRequireITSRefit(kTRUE);
    fESDtrackCuts_V0daugh -> SetMinNCrossedRowsTPC(80);
    fESDtrackCuts_V0daugh -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fESDtrackCuts_V0daugh -> SetMaxChi2PerClusterTPC(4.0);
    fESDtrackCuts_V0daugh -> SetMaxChi2PerClusterITS(36);
    fESDtrackCuts_V0daugh -> SetMinNClustersITS(4);


    //Intervals
    Double_t pt[] = {0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
    Double_t y[]  = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7};
    const Int_t nBins_pt = sizeof(pt)/sizeof(Double_t)-1;
    const Int_t nBins_y  = sizeof(y)/sizeof(Double_t)-1;

    //Matching Efficiency
    hAntiprotonsTPC = new TH2F ("hAntiprotonsTPC","",nBins_pt,pt,nBins_y,y);
    // hAntiprotonsTOF = new TH2F ("hAntiprotonsTOF","",nBins_pt,pt,nBins_y,y);
    hAntiprotonsTPC -> Sumw2();
    // hAntiprotonsTOF -> Sumw2();
    fOutputList -> Add (hAntiprotonsTPC);
    // fOutputList -> Add (hAntiprotonsTOF);

    //QA
    hnSigmaProtons_vs_Pt = new TH2F ("hnSigmaProtons_vs_Pt","",200,0,2,200,-10,10);
    hnSigmaProtons_vs_Pt -> Sumw2();
    fOutputList -> Add (hnSigmaProtons_vs_Pt);


    //Post Data
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserExec(Option_t *)  {

    //Get Input Event
    if (!GetEvent()) return;

    if (fRunData) RunData ();
    if (fMatchingEff) MatchingEff ();
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::RunData()  {

    //TPC Pre-selection
    Double_t nsigmaTPCmax[50] = {3.0,2.62945,2.69333,3.1122,3.28161,2.70606,2.89633,2.70347,3.31368,2.65839,2.93157,3.10837,3.43406,2.51501,3.49459,3.27995,3.05988,2.89139,2.60827,3.16441,2.69537,2.71876,2.67132,2.80894,3.42022,3.23984,2.96282,2.84138,2.78734,3.102,3.15273,2.95713,2.90697,3.05029,2.86195,3.16312,3.43193,2.57391,3.38719,3.177,3.38128,3.40739,2.83285,2.7041,3.05584,2.58265,2.99982,3.40591,2.64627,3.45544};

    //ITS Pre-selection
    Double_t nsigmaITSmax[50] = {3.0,2.62945,2.69333,3.1122,3.28161,2.70606,2.89633,2.70347,3.31368,2.65839,2.93157,3.10837,3.43406,2.51501,3.49459,3.27995,3.05988,2.89139,2.60827,3.16441,2.69537,2.71876,2.67132,2.80894,3.42022,3.23984,2.96282,2.84138,2.78734,3.102,3.15273,2.95713,2.90697,3.05029,2.86195,3.16312,3.43193,2.57391,3.38719,3.177,3.38128,3.40739,2.83285,2.7041,3.05584,2.58265,2.99982,3.40591,2.64627,3.45544};

    //Proton Mass
    const Double_t mass = AliPID::ParticleMass(AliPID::kProton);

    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fESDEvent->GetNumberOfTracks() ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDEvent->GetTrack(i);
        if (!track) continue;
        if (track->Pt()>2.0)   continue;
        if (track->Charge()>0) continue;

        //Variables
        Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
        Double_t eta       = track->Eta();
        Double_t pt        = track->Pt();
        Double_t p         = track->P();
        Double_t y         = GetRapidity(track,mass);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();

        for (Int_t isyst=0 ; isyst<50 ; isyst++)  {

            //Variables
            Double_t var_tpc[4] = {static_cast<double>(isyst)+0.5,TMath::Abs(y),pt,nsigmaTPC};
            Double_t var_tof[4] = {static_cast<double>(isyst)+0.5,TMath::Abs(y),pt,nsigmaTOF};
            Double_t tpc_rap[3] = {y,pt,nsigmaTPC};
            Double_t tof_rap[3] = {y,pt,nsigmaTOF};

            //TOF Analysis
            if (pt<0.5)       continue;
            if (!hasTOFhit)   continue;
            if (length<350.0) continue;
            if (TMath::Abs(nsigmaTPC)>nsigmaTPCmax[isyst]) continue;

            //Fill TOF Spectra
            hTOFnsigma -> Fill (var_tof);
            if (isyst==0) hTOFnsigma_vs_rap -> Fill (tof_rap);
        }
    }

    //Check
    if (fESDEvent->GetNumberOfTracks()==0) hNumberOfEvents->Fill(18.5);

    //Post Output Data
    PostData(1, fOutputList);
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::MatchingEff()  {

    //Select Clean Antiproton Candidates
    vector<Int_t> prot_ID;

    //Proton Mass
    const Double_t mass = AliPID::ParticleMass(AliPID::kProton);

    //Loop Over Reconstructed V0s
    for (Int_t i=0 ; i<fESDEvent->GetNumberOfV0s() ; i++)  {

        //Get V0 Candidate
        AliESDv0 *V0 = (AliESDv0*)fESDEvent->GetV0(i);
        if (!V0) continue;
        if (V0->GetOnFlyStatus()) continue;//Offline V0s Selection

        //Get V0 Daughters
        AliESDtrack *posTrack = (AliESDtrack*) fESDEvent->GetTrack(V0->GetPindex());
        AliESDtrack *negTrack = (AliESDtrack*) fESDEvent->GetTrack(V0->GetNindex());
        if (!posTrack) continue;
        if (!negTrack) continue;
        if (posTrack->Charge() == negTrack->Charge()) continue;
        if (posTrack->GetID()  == negTrack->GetID() ) continue;

        //Store Candidate IDs
        prot_ID.push_back(V0->GetNindex());
    }

    //Fill Histogram for AntiProtons
    Int_t nProtons = (Int_t)prot_ID.size();
    for (Int_t i=0 ; i<nProtons ; i++)  {

        //Get Track
        AliESDtrack *track = static_cast<AliESDtrack*>(fESDEvent->GetTrack(prot_ID[i]));
        if (track->Pt()>2.0)   continue;
        if (track->Charge()>0) continue;

        //Variables
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
        // Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
        Double_t pt        = track->Pt();
        Double_t y         = GetRapidity(track,mass);
        // Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        // Double_t length    = track->GetIntegratedLength();

        //Fill QA
        hnSigmaProtons_vs_Pt -> Fill (pt,nsigmaTPC);

        //Fill TPC
        if (TMath::Abs(nsigmaTPC)>3.0) continue;
        hAntiprotonsTPC -> Fill (pt,TMath::Abs(y));

        //TOF Analysis
        // if (!hasTOFhit)     continue;
        // if (length<350.0)   continue;
        // if (nsigmaTOF<-3.0) continue;
        // if (nsigmaTOF>+3.5) continue;

        //Fill TOF
        // hAntiprotonsTOF -> Fill (pt,TMath::Abs(y));
    }


    //Post Output Data
    PostData(1, fOutputList);
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFemto::GetEvent ()  {

    //Get Input Event
    fESDEvent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDEvent) return kFALSE;
    hNumberOfEvents -> Fill(0.5);

    //Standard Event Cuts
    if (!fESDEventSelection.AcceptEvent(fESDEvent)) {
        PostData(2, fQAList);
        return kFALSE;
    }
    hNumberOfEvents -> Fill(1.5);

    //Reject Events with Incomplete DAQ
    if (fESDEvent->IsIncompleteDAQ()) return kFALSE;
    hNumberOfEvents -> Fill(2.5);

    //V0 Timing Decision
    AliVVZERO *vzeroData = fESDEvent->GetVZEROData();
    if (!(vzeroData->GetV0ADecision()) || !(vzeroData->GetV0CDecision())) return kFALSE;
    hNumberOfEvents -> Fill(3.5);

    //Pileup Rejection
    Int_t nClustersLayer0 = fESDEvent->GetNumberOfITSClusters(0);
    Int_t nClustersLayer1 = fESDEvent->GetNumberOfITSClusters(1);
    Int_t nTracklets      = fESDEvent->GetMultiplicity()->GetNumberOfTracklets();
    if ((nClustersLayer0 + nClustersLayer1) > 65.0 + (Double_t)nTracklets*4.0) return kFALSE;
    hNumberOfEvents -> Fill(4.5);

    //Primary Vertex Tracks
    AliESDVertex *vertex_tracks = (AliESDVertex*) fESDEvent->GetPrimaryVertexTracks();
    if (!vertex_tracks) return kFALSE;
    hNumberOfEvents -> Fill(5.5);

    //Vertex Contributors Tracks
    if ( vertex_tracks->GetNContributors() < 1 ) return kFALSE;
    hNumberOfEvents -> Fill(6.5);

    //Primary Vertex SPD
    AliESDVertex *vertex_SPD = (AliESDVertex*) fESDEvent->GetPrimaryVertexSPD();
    if (!vertex_SPD) return kFALSE;
    hNumberOfEvents -> Fill(7.5);

    //Vertex Contributors SPD
    if ( vertex_SPD->GetNContributors() < 1 ) return kFALSE;
    hNumberOfEvents -> Fill(8.5);

    //SPD Pile-up in Mult Bins
    if (fESDEvent->IsPileupFromSPDInMultBins()) return kFALSE;
    hNumberOfEvents -> Fill(9.5);

    //Cut on Z-Vertex Resolution
    if (TMath::Abs(vertex_SPD->GetZ() - vertex_tracks->GetZ()) > 0.3) return kFALSE;
    hNumberOfEvents -> Fill(10.5);

    //Primary Vertex Selection
    if ( vertex_tracks->GetZ() < -10.0 ) return kFALSE;
    if ( vertex_tracks->GetZ() > +10.0 ) return kFALSE;
    hNumberOfEvents -> Fill(11.5);

    //Multiplicity
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDEvent->FindListObject("MultSelection");
    if( !multiplicitySelection) return kFALSE;
    hNumberOfEvents -> Fill(12.5);

    //Multiplicity Distribution
    Double_t mult_percentile = multiplicitySelection->GetMultiplicityPercentile("V0M");
    hMultiplicity -> Fill(mult_percentile);

    Int_t mult(0);
    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fESDEvent->GetNumberOfTracks() ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDEvent->GetTrack(i);
        if (!track) continue;

        mult++;
    }

    hMultDistribution -> Fill (mult);

    //Selection of Multiplicity Range
    if (mult_percentile <   0) return kFALSE;
    if (mult_percentile > 100) return kFALSE;
    hNumberOfEvents -> Fill(13.5);

    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    return kTRUE;
}
//__________________________________________________________________________________________________________________________________
/* Double_t AliAnalysisTaskJetFemto::GetDecayLengthV0 (AliESDv0 *V0)  {

    //Initialization
    Double_t decayLengthV0 = 0.0;

    //Secondary Vertex Position
    Double_t secVertex[3] = { 0.0, 0.0, 0.0 };
    V0->GetXYZ(secVertex[0],secVertex[1],secVertex[2]);

    //Primary Vertex Position
    AliESDVertex *vertex = (AliESDVertex*) fESDEvent->GetPrimaryVertex();
    Double_t primVertex[3] = { 0.0, 0.0, 0.0 };
    vertex->GetXYZ(primVertex);

    //Decay Length
    Double_t Dx = primVertex[0]-secVertex[0];
    Double_t Dy = primVertex[1]-secVertex[1];
    Double_t Dz = primVertex[2]-secVertex[2];
    decayLengthV0 = TMath::Sqrt(Dx*Dx + Dy*Dy + Dz*Dz);

    return decayLengthV0;
} */
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetFemto::GetRapidity (AliESDtrack *track, Double_t mass)  {

    //Initialization
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
void AliAnalysisTaskJetFemto::Terminate(Option_t *)  {

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//__________________________________________________________________________________________________________________________________
