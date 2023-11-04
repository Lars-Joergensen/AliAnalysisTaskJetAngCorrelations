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
    void SetRunningMode (Bool_t runData) {
        fRunData      = runData;
    }

    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    void  RunData();

    //User Functions
    Bool_t   GetEvent ();
    // Double_t GetRapidity                 (AliESDtrack *track, Double_t mass);

private:
    AliESDEvent     *fESDEvent;//!
    AliPIDResponse  *fPIDResponse;//!
    TList           *fOutputList;//!
    TList           *fQAList;//!
    Bool_t           fRunData;//
    TH1F *hNumberOfEvents;//!

    AliAnalysisTaskJetFemto(const AliAnalysisTaskJetFemto&);
    AliAnalysisTaskJetFemto& operator=(const AliAnalysisTaskJetFemto&);

    ClassDef(AliAnalysisTaskJetFemto, 1);

};
//____________________________________________________________________________________________________________________________________________________

#endif
