#ifndef AliAnalysisTaskJetFemto_cxx
#define AliAnalysisTaskJetFemto_cxx

//======================== Antiprotons vs. Rapidity ========================//
//                                                                          //
//    Antiproton production vs. transverse momentum and rapidity            //
//    for the calculation of the coalescence parameters B_{2} and B_{3}.    //
//                                                                          //
//==========================================================================//

#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliESDVertex.h"
#include "AliEventCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

//____________________________________________________________________________________________________________________________________________________
class AliAnalysisTaskJetFemto : public AliAnalysisTaskSE {

public:
    AliAnalysisTaskJetFemto();
    AliAnalysisTaskJetFemto(const char *name);
    virtual ~AliAnalysisTaskJetFemto();

    //Running Mode
    void SetRunningMode (Bool_t isITSrecalib, Bool_t runData, Bool_t matchingEff) {

        fIsITSrecalib = isITSrecalib;
        fRunData      = runData;
        fMatchingEff  = matchingEff;
    }

    //Set ITS Recalibration maps
    void SetITSRecalibrationMaps (TH2F *hITSnsigma_Mean, TH2F *hITSnsigma_Width)  {

        hMean  = hITSnsigma_Mean;
        hWidth = hITSnsigma_Width;
    }

    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    void  RunData();
    void  MatchingEff();

    //User Functions
    Bool_t   GetEvent ();
    Bool_t   PassedTrackSelection        (AliESDtrack *track, Int_t isyst);
    Bool_t   IsHighPurityProton          (AliESDtrack *track);
    Bool_t   IsProtonCandidate           (AliESDtrack *track);
    Double_t GetDCAtoPrimaryVertex       (AliESDtrack *track, Int_t index);
    Double_t GetRapidity                 (AliESDtrack *track, Double_t mass);
    Double_t GetRecalibratedITSnsigma    (Double_t nsigma, Double_t eta, Double_t p);
    Bool_t   IsPionCandidate             (AliESDtrack *track);
    //Bool_t   PassedV0Selection           (AliESDv0 *V0);
    //Bool_t   PassedTrackSelectionV0daugh (AliESDtrack *track);
    Bool_t   PassedAntiLambdaSelection   (AliESDv0 *V0);
    Double_t GetDecayLengthV0            (AliESDv0 *V0);
    Double_t MassLambda                  (TVector3 Ppion, TVector3 Pprot);

    //Standard Event Selection
    AliEventCuts  fESDEventSelection;//

private:
    AliESDEvent     *fESDEvent;//!
    AliPIDResponse  *fPIDResponse;//!
    AliESDtrackCuts *fESDtrackCuts[50];//!
    AliESDtrackCuts *fESDtrackCuts_V0daugh;//!
    TList           *fOutputList;//!
    TList           *fQAList;//!
    Bool_t           fIsITSrecalib;//
    Bool_t           fRunData;//
    Bool_t           fMatchingEff;//

    TH2F *hMean;//
    TH2F *hWidth;//

    //Event Counter and Centrality Distribution
    TH1F *hNumberOfEvents;//!
    TH1F *hMultiplicity;//!
    TH1I *hMultDistribution;//!

    //n-Dimensional Histograms
    THnSparse *hTPCnsigma;//!
    THnSparse *hTOFnsigma;//!
    THnSparse *hDCAxy;//!

    //2D ITS Recalibration Map
    TH3F *hITSnsigma;//!

    //n-Dimensional Histograms (y>0 vs. y<0)
    THnSparse *hTPCnsigma_vs_rap;//!
    THnSparse *hTOFnsigma_vs_rap;//!

    //Matching Efficiency
    TH2F *hAntiprotonsTPC;//!
    TH2F *hAntiprotonsTOF;//!

    //QA
    TH2F *hnSigmaProtons_vs_Pt;//!

    AliAnalysisTaskJetFemto(const AliAnalysisTaskJetFemto&);
    AliAnalysisTaskJetFemto& operator=(const AliAnalysisTaskJetFemto&);

    ClassDef(AliAnalysisTaskJetFemto, 1);

};
//____________________________________________________________________________________________________________________________________________________

#endif
