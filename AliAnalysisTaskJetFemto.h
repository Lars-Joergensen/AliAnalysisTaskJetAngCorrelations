#ifndef AliAnalysisTaskJetFemto_cxx
#define AliAnalysisTaskJetFemto_cxx

//======================== Jet Femtoscopy ========================//
//                                                                //
//    Reconstructing jets and performing a femtoscopic analysis   //
//             on them to find the source size in jets             //
//                                                                //
//================================================================//

#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliAODTrackSelection.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrack.h"
#include "fastjet/PseudoJet.hh"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TObjString.h"
#include "THnSparse.h"
// #include "RTypes.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <vector>

//___________________________________________________________________________________________________________________________________________
class AliAnalysisTaskJetFemto : public AliAnalysisTaskSE {

public:
    AliAnalysisTaskJetFemto();
    AliAnalysisTaskJetFemto(const char *name, bool isMC);
    virtual ~AliAnalysisTaskJetFemto();

    //Running Mode
    void SetRunningMode (Bool_t runData) {
        fRunData      = runData;
    }

    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    void RunData();
    void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts=evtCuts;};
    void SetTrackCutsProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProton=trkCuts;};
    void SetTrackCutsAntiproton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiproton=trkCuts;};
    void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig=config;};
    void SetTrigger(UInt_t trigger) { fTrigger = trigger;};
    void SetIsMC(bool isMC) { fIsMC = isMC;};

    //User Functions
    Bool_t   GetEvent                   ();
    Bool_t   IsProton                   (AliAODTrack  track);
    Bool_t   IsPion                     (AliAODTrack  track);
    Bool_t   PassedTrackSelection       (AliAODTrack  track);
    Double_t GetDCAtoPrimaryVertex      (AliAODTrack  track, Int_t index);
    Double_t GetRapidity                (AliAODTrack  track, Double_t mass);

private:
    void ResetGlobalTrackReference();
    void StoreGlobalTrackReference(AliAODTrack *track);
    bool fIsMC;                     //
    AliAODEvent     *fAODEvent;     //!
    AliPIDResponse  *fPIDResponse;  //!
    AliAODTrackSelection *fAODTrackCuts;         //!
    AliESDtrackCuts *fESDTrackCuts; //!
    TList           *fJetOutput;    //!
    TList           *fQAList;       //!
    TList           *fFemtoOutput;  //!
    Bool_t           fRunData;      //
    UInt_t           fTrigger;      //
    AliFemtoDreamEvent *fEvent;                  //!
    AliFemtoDreamTrack *fTrack;                  //!
    AliFemtoDreamEventCuts *fEventCuts;          //
    AliFemtoDreamTrackCuts *fTrackCutsProton;    //
    AliFemtoDreamTrackCuts *fTrackCutsAntiproton;//
    AliFemtoDreamCollConfig *fConfig;            //
    AliFemtoDreamPairCleaner *fPairCleaner;      //!
    AliFemtoDreamPartCollection *fPartColl;      //!
    AliAODTrack** fGTI;             //!

    // Histograms
    TH1F *hNumberOfEvents;          //!
    TH1F *hNumberOfJets;            //!
    TH1F *hNumberOfHighPtJets;      //!
    TH1F *hEventProtocol;           //!
    TH1F *hTrackProtocol;           //!
    TH1I *hJetParticleID;           //!
    TH1I *hParticleID;              //!
    TH1I *hUserIndex;               //!
    TH1F *hIDsafe;                  //!
    TH1F *hIDfail;                  //!

    // TH1F *hDeltaY;                  //!
    TH1F *hEtaFullEvent;            //!

    TH1D *hPtFullEvent;             //!
    TH1D *hFullConePt;              //!
    TH1D *hPtJetParticle;           //!
    TH1D *hJetRapidity;             //!
    TH1D *hPtTotalJet;              //!
    TH1D *hConeBkgPt;               //!
    TH1D *hPtSubtractedJet;         //!
    TH1D *hPtJetProton;             //!
    TH1D *hPtJetPion;               //!
    TH1D *hPtDiff;                  //!

    TH1I *hNumberOfParticlesInJet;  //!
    TH1I *hNumberOfJetsInEvent;     //!
    TH1I *hNumberOfTracksBeforeCuts;//!

    TH1D *hJetConeRadius;           //!

    TH2F *hTPCnsigma;               //!
    TH2F *hTPCnsigmaProton;         //!
    TH2F *hTOFnsigma;               //!
    TH2F *hTOFnsigmaProton;         //!
    TH2F *hITSnsigma;               //!
    TH1F *hTPCnsigmaBins[33];       //!
    TH1F *hTOFnsigmaBins[23];       //!

    TH1F *hDCAxyFullEvent;          //!
    TH1F *hDCAzFullEvent;           //!
    TH1F *hDCAzJetProton;                 //!
    TH1F *hDCAzJetProtons[10];      //!
    TH1D *hPtJetProtonDCAz[10];     //!
    TH2F *hTPCnsigmaJetProtonDCAz[10]; //!
    TH2F *hTOFnsigmaJetProtonDCAz[10]; //!

    Double_t fJetRadius;
    Double_t fMinJetPt;
    int fTrackBufferSize;           //

    AliAnalysisTaskJetFemto(const AliAnalysisTaskJetFemto&);
    AliAnalysisTaskJetFemto& operator=(const AliAnalysisTaskJetFemto&);

    ClassDef(AliAnalysisTaskJetFemto, 1);

};
//___________________________________________________________________________________________________________________________________________
#endif


//___________________________________________________________________________________________________________________________________________