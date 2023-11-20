#ifdef __CLING__

R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/train/AddAODHandler.C>
#include <ANALYSIS/macros/AddTaskPIDResponse.C>

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>

R__LOAD_LIBRARY(AliAnalysisTaskJetFemto_cxx.so)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskJetFemto.h"
#endif
//===========================================================================================================================================

//Definitions
#define ALIPHYSICS_VER  "vAN-20220531_ROOT6-1"
#define GRIDWorkingDir  "JET_FEMTOSCOPY_WORKING_DIR"
#define AnalysisMacro   "AnalysisJetFemtoscopy"
#define AnalysisTask    "AliAnalysisTaskJetFemto"
#define TTL             20000
#define nRunsPerMaster  1

//Functions
AliAnalysisGrid *CreateAlienHandler (const char *mode, Bool_t merge );
void LoadAnalysisTask (AliAnalysisManager *mgr);
void EventHandler     (AliAnalysisManager *mgr);
void LoadPIDResponse();

//Test
Bool_t test = (kFALSE);

//___________________________________________________________________________________________________________________________________________
void LoadAnalysisTask (AliAnalysisManager *mgr)  {

    //Run Mode
    Bool_t isITSrecalib = (kFALSE);
    Bool_t runData      = (kTRUE);

    //Load Analysis Task
    gROOT->LoadMacro("AliAnalysisTaskJetFemto.cxx+g");

    //Input container
    AliAnalysisDataContainer *input = mgr -> GetCommonInputContainer();

    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();

    //Analysis Task
    AliAnalysisTaskJetFemto *task = new AliAnalysisTaskJetFemto ("task_JetReconstruction_Femtoscopy");
    // task -> SelectCollisionCandidates (AliVEvent::kINT7);
    task -> AliAnalysisTaskJetFemto::SetRunningMode(runData);
    mgr -> AddTask(task);
    mgr -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput(task,1,mgr->CreateContainer("Jets", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr -> ConnectOutput(task,2,mgr->CreateContainer("QAPlots", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
}
//___________________________________________________________________________________________________________________________________________
void runAnalysisJetFemto (const char *mode="full", Bool_t merge=kTRUE)  {

    //Grid Connection
    TGrid::Connect("alien://");

    //Alien Handler
    AliAnalysisGrid *alienGridHandler = CreateAlienHandler (mode,merge);
    if (!alienGridHandler) return;

    //Analysis Manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager");
    mgr->SetGridHandler(alienGridHandler);

    //Event Handler, PID
    EventHandler(mgr);
    LoadPIDResponse();

    //Analysis Task
    LoadAnalysisTask(mgr);

    //Start analysis
    if (!mgr->InitAnalysis())  return;
    mgr->PrintStatus();
    mgr->StartAnalysis("grid");
};
//___________________________________________________________________________________________________________________________________________
AliAnalysisGrid *CreateAlienHandler (const char *mode, Bool_t merge )  {

    //Alien Handler
    AliAnalysisAlien *alien = new AliAnalysisAlien();
    alien->SetOverwriteMode();
    alien->SetCheckCopy(kFALSE);
    alien->SetRunMode(mode);
    alien->SetNtestFiles(10);
    alien->SetAPIVersion("V1.1x");
    alien->SetAliPhysicsVersion(ALIPHYSICS_VER);
    alien->AddIncludePath("$ALICE_PHYSICS/include");
    alien->AddIncludePath("$ALICE_ROOT/include");
    alien->SetGridDataDir("/alice/data/2011/LHC11h_2");
    alien->SetRunPrefix("000");
    alien->SetDataPattern("*ESDs/pass2/AOD145/*AOD.root");
    alien->AddRunNumber(167813);
    alien->SetNrunsPerMaster(nRunsPerMaster);
    alien->SetGridWorkingDir(Form("%s",GRIDWorkingDir));
    alien->SetGridOutputDir("OUTPUT");
    alien->SetAnalysisSource(Form("%s.cxx",AnalysisTask));
    alien->SetAdditionalLibs(Form("%s.cxx  %s.h",AnalysisTask,AnalysisTask));
    alien->SetMergeViaJDL(merge);
    alien->SetMaxMergeStages(1);
    alien->SetAnalysisMacro(Form("%s.C",AnalysisMacro));
    alien->SetSplitMaxInputFileNumber(50);
    alien->SetMasterResubmitThreshold(90);
    alien->SetTTL(TTL);
    alien->SetExecutable(Form("%s.sh",AnalysisMacro));
    alien->SetInputFormat("xml-single");
    alien->SetJDLName(Form("%s.jdl",AnalysisMacro));
    alien->SetMergeExcludes("EventStat_temp.root");
    alien->SetPrice(1);
    alien->SetSplitMode("se");
    return alien;
}
//___________________________________________________________________________________________________________________________________________
void EventHandler (AliAnalysisManager *mgr)  {

    AliAODInputHandler *inputEventHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(inputEventHandler);
}
//___________________________________________________________________________________________________________________________________________
void LoadPIDResponse ()  {

    Bool_t isMC = kFALSE;
    AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(isMC);
}
//___________________________________________________________________________________________________________________________________________
