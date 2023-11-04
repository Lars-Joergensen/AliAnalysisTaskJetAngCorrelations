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
fOutputList(nullptr),
fQAList(nullptr),
fRunData(kFALSE),
hNumberOfEvents(nullptr)
{}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::AliAnalysisTaskJetFemto(const char *name):
AliAnalysisTaskSE(name),
fESDEvent(nullptr),
fPIDResponse(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fRunData(kFALSE),
hNumberOfEvents(nullptr)
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
    delete fOutputList;
    delete fQAList;
    delete hNumberOfEvents;
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserCreateOutputObjects()  {

    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();

    //Event Counter and Centrality Distribution
    hNumberOfEvents = new TH1F ("hNumberOfEvents","",20,0,20);
    hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);

    //Post Data
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserExec(Option_t *)  {

    //Get Input Event
    if (!GetEvent()) return;

    if (fRunData) RunData ();
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::RunData()  {
    //Check
    if (fESDEvent->GetNumberOfTracks()==0) hNumberOfEvents->Fill(18.5);

    //Post Output Data
    PostData(1, fOutputList);
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFemto::GetEvent ()  {

    //Get Input Event
    fESDEvent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDEvent) return kFALSE;
    hNumberOfEvents -> Fill(0.5);

    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    return kTRUE;
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::Terminate(Option_t *)  {

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//__________________________________________________________________________________________________________________________________
