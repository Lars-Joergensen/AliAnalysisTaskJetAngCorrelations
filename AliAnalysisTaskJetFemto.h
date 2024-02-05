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
    Bool_t   PassedTrackSelection       (AliAODTrack *track);
    Double_t GetDCAtoPrimaryVertex      (AliAODTrack *track, Int_t index);
    Double_t GetRapidity                (AliAODTrack *track, Double_t mass);

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
    TH1F *hNumberOfTracks;          //!
    TH2F *hLeadingIDs;              //!
    TH1F *hLeadingEta;              //!
    TH1F *hDeltaEta;                //!
    TH1F *hFullEta;                 //!
    TH1D *hFullPt;                  //!
    TH1D *hJetPt;                   //!
    TH1D *hJetProtonPt;             //!
    TH1D *hLeadingPt;               //!
    // v not yet used v
    TH1F *hProtonYield;             //!
    TH2F *hTPCnsigma;               //!
    TH2F *hTOFnsigma;               //!
    TH1F *hDCAxy;                   //!
    TH1F *hDCAz;                    //!
    TH1F *hDCAzJet;                 //!
    TH1F *hDeltaTheta;              //!

    Double_t fJetRadius;
    Double_t fMaximumPt;
    int fTrackBufferSize;           //

    AliAnalysisTaskJetFemto(const AliAnalysisTaskJetFemto&);
    AliAnalysisTaskJetFemto& operator=(const AliAnalysisTaskJetFemto&);

    ClassDef(AliAnalysisTaskJetFemto, 1);

};
//___________________________________________________________________________________________________________________________________________

#endif
