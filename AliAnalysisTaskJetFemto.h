#ifndef AliAnalysisTaskJetFemto_cxx
#define AliAnalysisTaskJetFemto_cxx

//======================== Jet Femtoscopy ========================//
//                                                                //
//    Reconstructing jets and performing a femtoscopic analysis   //
//             on them to find the source size in jets             //
//                                                                //
//================================================================//

#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
// #include "AliAODtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <vector>

//___________________________________________________________________________________________________________________________________________
class AliAnalysisTaskJetFemto : public AliAnalysisTaskSE {

public:
    AliAnalysisTaskJetFemto();
    AliAnalysisTaskJetFemto(const char *name);
    virtual ~AliAnalysisTaskJetFemto();

    //Running Mode
    void SetRunningMode (Bool_t runData) {
        fRunData      = runData;
    }

    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    void  RunData();

    //User Functions
    Bool_t  GetEvent                ();
    Bool_t  IsHighPurityProton      (AliAODTrack *track);
    Bool_t  IsProtonCandidate       (AliAODTrack *track);
    Double_t GetRapidity            (AliAODTrack *track, Double_t mass);
    Double_t Minimum                (Double_t x1, Double_t x2);
    Double_t GetDCAtoPrimaryVertex  (AliAODTrack *track, Int_t index);

private:
    AliAODEvent     *fAODEvent;     //!
    AliPIDResponse  *fPIDResponse;  //!
    TList           *fOutputList;   //!
    TList           *fQAList;       //!
    Bool_t           fRunData;      //

    TH1F *hNumberOfEvents;          //!
    TH1F *hNumberOfTracks;          //!
    TH1I *hParticleIDs;             //!
    TH1I *hLeadingIDs;              //!
    TH1D *hFullPt;                  //!
    TH1D *hJetPt;                   //!
    TH1F *hProtonYield;             //!
    TH1F *hTPCnsigma;               //!
    TH1F *hTOFnsigma;               //!
    TH1F *hDCAxy;                   //!
    // ^ actually THnSparse (THnSparseF in cxx)

    // vector<AliAODTrack*> leadingParticles; //!
    Double_t fJetRadius;

    AliAnalysisTaskJetFemto(const AliAnalysisTaskJetFemto&);
    AliAnalysisTaskJetFemto& operator=(const AliAnalysisTaskJetFemto&);

    ClassDef(AliAnalysisTaskJetFemto, 1);

};
//___________________________________________________________________________________________________________________________________________

#endif
