#ifdef __CLING__

R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/train/AddESDHandler.C>
#include <ANALYSIS/macros/AddTaskPIDResponse.C>

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>

R__LOAD_LIBRARY(AliAnalysisTaskJetFemto_cxx.so)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskJetFemto.h"
#endif
//================================================================================================================

//Definitions
#define ALIPHYSICS_VER  "vAN-20220531_ROOT6-1"
#define GRIDWorkingDir  "JET_FEMTOSCOPY_WORKING_DIR"
#define AnalysisMacro   "AnalysisJetFemtoscopy"
#define AnalysisTask    "AliAnalysisTaskJetFemto"
#define TTL             20000
#define nRunsPerMaster  1

//Functions
AliAnalysisGrid *CreateAlienHandler (Int_t iChunk, const char *mode, Bool_t merge );
void LoadAnalysisTask (Int_t iChunk, AliAnalysisManager *mgr);
void EventHandler     (AliAnalysisManager *mgr);
void LoadPIDResponse();
// void LoadCentrality();
void SetInputRuns (Int_t iChunk, AliAnalysisAlien *alien, const char *mode);

//Test
Bool_t test = (kFALSE);

//______________________________________________________________________________________________________________________________________________________
void LoadAnalysisTask (Int_t iChunk, AliAnalysisManager *mgr)  {

    //Run Mode
    Bool_t isITSrecalib = (kFALSE);
    Bool_t runData      = (kFALSE);

    //Load Analysis Task
    gROOT->LoadMacro("AliAnalysisTaskJetFemto.cxx+g");

    //Input container
    AliAnalysisDataContainer *input = mgr -> GetCommonInputContainer();

    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();

    //Analysis Task
    AliAnalysisTaskJetFemto *task = new AliAnalysisTaskJetFemto ("task_AntiProtons_vs_y");
    task -> SelectCollisionCandidates (AliVEvent::kINT7);
    task -> AliAnalysisTaskJetFemto::SetRunningMode(runData);
    mgr -> AddTask(task);
    mgr -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput(task,1,mgr->CreateContainer("Antiprotons", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr -> ConnectOutput(task,2,mgr->CreateContainer("QAPlots", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
}
//______________________________________________________________________________________________________________________________________________________
void runAnalysisJetFemto (Int_t iChunk=0, const char *mode="full", Bool_t merge=kTRUE)  {

    //Grid Connection
    TGrid::Connect("alien://");

    //Alien Handler
    AliAnalysisGrid *alienHandler = CreateAlienHandler (iChunk,mode,merge);
    if (!alienHandler) return;

    //Analysis Manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager");
    mgr->SetGridHandler(alienHandler);

    //Event Handler, PID & Centrality
    EventHandler(mgr);
    // LoadCentrality ();
    LoadPIDResponse();

    //Analysis Task
    LoadAnalysisTask(iChunk,mgr);

    //Start analysis
    if (!mgr->InitAnalysis())  return;
    mgr->PrintStatus();
    mgr->StartAnalysis("grid");
};
//______________________________________________________________________________________________________________________________________________________
AliAnalysisGrid *CreateAlienHandler (Int_t iChunk, const char *mode, Bool_t merge )  {

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
    alien->SetGridDataDir("/alice/data/2018/LHC18m");
    alien->SetDataPattern("/pass2/*/AliESDs.root");
    alien->SetRunPrefix("000");
    SetInputRuns (iChunk,alien,mode);
    alien->SetNrunsPerMaster(nRunsPerMaster);
    alien->SetGridWorkingDir(Form("%s_%d",GRIDWorkingDir,iChunk));
    alien->SetGridOutputDir("OUTPUT");
    alien->SetAnalysisSource(Form("%s.cxx",AnalysisTask));
    alien->SetAdditionalLibs(Form("%s.cxx  %s.h",AnalysisTask,AnalysisTask));
    alien->SetMergeViaJDL(merge);
    alien->SetMaxMergeStages(2);
    alien->SetAnalysisMacro(Form("%s_%d.C",AnalysisMacro,iChunk));
    alien->SetSplitMaxInputFileNumber(50);
    alien->SetMasterResubmitThreshold(90);
    alien->SetTTL(TTL);
    alien->SetExecutable(Form("%s_%d.sh",AnalysisMacro,iChunk));
    alien->SetInputFormat("xml-single");
    alien->SetJDLName(Form("%s_%d.jdl",AnalysisMacro,iChunk));
    alien->SetMergeExcludes("EventStat_temp.root");
    alien->SetPrice(1);
    alien->SetSplitMode("se");
    return alien;
}
//______________________________________________________________________________________________________________________________________________________
void EventHandler (AliAnalysisManager *mgr)  {

    AliESDInputHandler *inputH = new AliESDInputHandler();
    mgr->SetInputEventHandler(inputH);
}
//______________________________________________________________________________________________________________________________________________________
void LoadPIDResponse ()  {

    Bool_t isMC = kFALSE;
    AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(isMC);
}
//______________________________________________________________________________________________________________________________________________________
/* void LoadCentrality()  {

    Bool_t lCalibration = kFALSE;
    AliMultSelectionTask *MultSelTask = AddTaskMultSelection(lCalibration);
} */
//______________________________________________________________________________________________________________________________________________________
void SetInputRuns (Int_t iChunk, AliAnalysisAlien *alien, const char *mode)  {

    //Run List: RunList_LHC18m_pass2_CentralBarrelTracking_electronPID.txt
    Int_t run0[] = { 292839, 292836, 292834, 292832, 292831, 292811, 292810, 292809, 292804, 292803, 292752, 292750, 292748, 292747, 292744, 292739, 292737, 292704, 292701, 292698, 292696, 292695, 292693, 292586, 292584, 292563, 292560, 292559, 292557, 292554, 292553, 292526, 292524, 292523, 292521, 292500, 292497, 292496, 292495, 292457, 292456, 292434, 292432, 292430, 292429, 292428, 292406, 292405, 292398, 292397, 292298, 292273, 292265, 292242, 292241, 292240, 292192, 292168, 292167, 292166, 292164, 292163, 292162, 292161, 292160, 292140, 292115, 292114
    };

    Int_t run1[] = { 292109 };

    Int_t run2[] = { 291116 };

    //Number of Runs
    Int_t nRuns(0);
    if (iChunk==0) nRuns = sizeof(run0)/sizeof(Int_t);
    if (iChunk==1) nRuns = sizeof(run1)/sizeof(Int_t);
    if (iChunk==2) nRuns = sizeof(run2)/sizeof(Int_t);

    //Loop Over Runs
    nRuns=1;
    //if (test) nRuns=1;
    for ( Int_t iRun=0 ; iRun<nRuns ; iRun++ )  {

        if (iChunk==0) alien->AddRunNumber(run0[iRun]);
        if (iChunk==1) alien->AddRunNumber(run1[iRun]);
        if (iChunk==2) alien->AddRunNumber(run2[iRun]);
    }
}
//______________________________________________________________________________________________________________________________________________________
