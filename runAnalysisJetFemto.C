#ifdef __CLING__

R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/train/AddESDHandler.C>
#include <ANALYSIS/macros/AddTaskPIDResponse.C>

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <OADB/macros/AddTaskPhysicsSelection.C>
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
#define GRIDWorkingDir  "ANTIPROTONS_MATCHING_EFFICIENCY_DATA"
#define AnalysisMacro   "AnalysisAntiProtons"
#define AnalysisTask    "AliAnalysisTaskJetFemto"
#define TTL             20000
#define nRunsPerMaster  1

//Functions
AliAnalysisGrid *CreateAlienHandler (Int_t iChunk, const char *mode, Bool_t merge );
void LoadAnalysisTask (Int_t iChunk, AliAnalysisManager *mgr);
void EventHandler     (AliAnalysisManager *mgr);
void LoadPhysicsSelection();
void LoadPIDResponse();
void LoadCentrality();
void SetInputRuns (Int_t iChunk, AliAnalysisAlien *alien, const char *mode);

//Test
Bool_t test = (kFALSE);

//______________________________________________________________________________________________________________________________________________________
void LoadAnalysisTask (Int_t iChunk, AliAnalysisManager *mgr)  {

    //Run Mode
    Bool_t isITSrecalib = (kFALSE);
    Bool_t runData      = (kFALSE);
    Bool_t matchingEff  = (kTRUE);

    //ITS recalibration map
    TFile *inputfile = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/nSigmaITS_Recalib/ITSrecalibrationMaps_data.root");
    TH2F *hITSnsigma_Mean  = (TH2F*)inputfile->Get ("hITSnsigma_Mean");
    TH2F *hITSnsigma_Width = (TH2F*)inputfile->Get ("hITSnsigma_Width");

    //Load Analysis Task
    gROOT->LoadMacro("AliAnalysisTaskJetFemto.cxx+g");

    //Input container
    AliAnalysisDataContainer *input = mgr -> GetCommonInputContainer();

    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();

    //Analysis Task
    AliAnalysisTaskJetFemto *task = new AliAnalysisTaskJetFemto ("task_AntiProtons_vs_y");
    task -> SelectCollisionCandidates (AliVEvent::kINT7);
    task -> AliAnalysisTaskJetFemto::SetITSRecalibrationMaps (hITSnsigma_Mean,hITSnsigma_Width);
    task -> AliAnalysisTaskJetFemto::SetRunningMode(isITSrecalib,runData,matchingEff);
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
    LoadPhysicsSelection ();
    LoadCentrality ();
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
void LoadPhysicsSelection()  {

    Bool_t isMC = kFALSE;
    AliPhysicsSelectionTask *PhySel = AddTaskPhysicsSelection(isMC);
}
//______________________________________________________________________________________________________________________________________________________
void LoadPIDResponse ()  {

    Bool_t isMC = kFALSE;
    AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(isMC);
}
//______________________________________________________________________________________________________________________________________________________
void LoadCentrality()  {

    Bool_t lCalibration = kFALSE;
    AliMultSelectionTask *MultSelTask = AddTaskMultSelection(lCalibration);
}
//______________________________________________________________________________________________________________________________________________________
void SetInputRuns (Int_t iChunk, AliAnalysisAlien *alien, const char *mode)  {

    //Run List: RunList_LHC18m_pass2_CentralBarrelTracking_electronPID.txt
    Int_t run0[] = { 292839, 292836, 292834, 292832, 292831, 292811, 292810, 292809, 292804, 292803, 292752, 292750, 292748, 292747, 292744, 292739, 292737, 292704, 292701, 292698, 292696, 292695, 292693, 292586, 292584, 292563, 292560, 292559, 292557, 292554, 292553, 292526, 292524, 292523, 292521, 292500, 292497, 292496, 292495, 292457, 292456, 292434, 292432, 292430, 292429, 292428, 292406, 292405, 292398, 292397, 292298, 292273, 292265, 292242, 292241, 292240, 292192, 292168, 292167, 292166, 292164, 292163, 292162, 292161, 292160, 292140, 292115, 292114
    };

    Int_t run1[] = { 292109, 292108, 292107, 292106, 292081, 292080, 292077, 292075, 292067, 292062, 292061, 292060, 292040, 292012, 291982, 291977, 291976, 291953, 291948, 291946, 291945, 291944, 291943, 291942, 291803, 291796, 291795, 291769, 291768, 291766, 291762, 291760, 291756, 291755, 291729, 291706, 291698, 291697, 291690, 291665, 291661, 291657, 291626, 291624, 291622, 291618, 291615, 291614, 291590, 291485, 291484, 291482, 291481, 291457, 291456, 291453, 291451, 291447, 291424, 291420, 291417, 291416, 291402, 291400, 291399, 291397, 291377, 291375, 291363, 291362, 291361, 291360, 291286, 291285, 291284, 291282, 291266, 291265, 291263, 291257, 291240, 291209, 291188, 291143 };

    Int_t run2[] = { 291116, 291111, 291110, 291101, 291100, 291093, 291069, 291066, 291065, 291041, 291037, 291035, 291006, 291005, 291004, 291003, 291002, 290980, 290979, 290976, 290975, 290974, 290948, 290944, 290943, 290941, 290935, 290932, 290895, 290894, 290888, 290887, 290886, 290862, 290860, 290853, 290848, 290846, 290843, 290841, 290790, 290787, 290766, 290689, 290687, 290665, 290660, 290645, 290632, 290627, 290615, 290614, 290613, 290612, 290590, 290588, 290553, 290550, 290549, 290544, 290540, 290539, 290538, 290501, 290500, 290499, 290469, 290467, 290459, 290458, 290456, 290427, 290426, 290425, 290423, 290412, 290411, 290404, 290401, 290399, 290376, 290375, 290374, 290350, 290327, 290323, 290300, 290298, 290297, 290294, 290293, 290254, 290253, 290223, 290222 };

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
