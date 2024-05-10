// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnalysisTaskJetSourceSize.h"
#include "AliAnalysisAlien.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "TChain.h"
#include </home/lars/alice/sw/ubuntu2004_x86-64/AliRoot/f00f49055a7b3a21cb11fcc41ddf25a617fc782a_O2-local1/ANALYSIS/macros/train/AddAODHandler.C>
#include </home/lars/alice/sw/ubuntu2004_x86-64/AliRoot/f00f49055a7b3a21cb11fcc41ddf25a617fc782a_O2-local1/ANALYSIS/macros/AddTaskPIDResponse.C>
#include </home/lars/alice/sw/ubuntu2004_x86-64/AliPhysics/0_O2-local1/OADB/macros/AddTaskPhysicsSelection.C>
#include </home/lars/alice/sw/ubuntu2004_x86-64/AliPhysics/0_O2-local1/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>

#define ALIPHYSICS_VERSION "vAN-20220531_ROOT6-1"
#define GRIDWORKINGDIR "JET_SOURCE_SIZE_WORKING_DIR"
#define NUMFOLDERS 1 // 1 -- 74

void SetInputRuns(AliAnalysisAlien *alien, const char *mode, int effort);

void runAnalysisJetSourceSize(TString mode="full", Bool_t merge=kFALSE)// full terminate test;  kFALSE kTRUE
{
    Bool_t local    = kFALSE; // run analysis locally (kTRUE) or on grid (kFALSE)
    Bool_t isMC     = kFALSE;
    Bool_t lCalib   = kFALSE;
    int effort      = 1;

    // since we will compile a class, tell root where to look for headers
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskJetSourceSize");
    AliAODInputHandler *inputHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(inputHandler);

    AddTaskPhysicsSelection(isMC);
    AddTaskMultSelection(lCalib);
    AddTaskPIDResponse(isMC);

    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisTaskJetSourceSize.cxx++g");
    AliAnalysisTaskJetSourceSize *task = reinterpret_cast<AliAnalysisTaskJetSourceSize*>(gInterpreter->ExecuteMacro("AddTaskJetSourceSize.C"));
#else
    gROOT->LoadMacro("AliAnalysisTaskJetSourceSize.cxx++g");
    gROOT->LoadMacro("AddTaskJetSourceSize.C");
    AliAnalysisTaskJetSourceSize *task = AddTaskJetSourceSize();
#endif

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(0);
    mgr->PrintStatus();
    // mgr->SetUseProgressBar(1, 25);

    if(local) {
        TChain* chain = new TChain("aodTree");
        for (int i=0; i<NUMFOLDERS; i++) {
            if (i<9) chain->Add(Form("/home/lars/alice/test-data/00%d/AliAOD.root",i+1));
            if (i>=9) chain->Add(Form("/home/lars/alice/test-data/0%d/AliAOD.root",i+1));
        }
        mgr->StartAnalysis("local", chain);
    } else {
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        alienHandler->SetAdditionalLibs("AliAnalysisTaskJetSourceSize.cxx AliAnalysisTaskJetSourceSize.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskJetSourceSize.cxx");
        alienHandler->SetAliPhysicsVersion(ALIPHYSICS_VERSION); // update! (when is the code included in daily tag)
        alienHandler->SetAPIVersion("V1.1x");
        alienHandler->SetGridDataDir("/alice/data/2018/LHC18b");
        alienHandler->SetDataPattern("*pass2/AOD/*AOD.root");
        alienHandler->SetRunPrefix("000");
        SetInputRuns(alienHandler, mode, effort);
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("JetSourceSize.sh");
        alienHandler->SetTTL(20000);
        alienHandler->SetJDLName("JetSourceSize.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(merge);

        alienHandler->SetGridWorkingDir(GRIDWORKINGDIR);
        alienHandler->SetGridOutputDir("OUTPUT");
        alienHandler->SetRunMode(mode);
        alienHandler->SetNtestFiles(1);

        mgr->SetGridHandler(alienHandler);
        mgr->StartAnalysis("grid");
    }
}

void SetInputRuns (AliAnalysisAlien *alien, const char *mode, int effort) {
    Int_t run0[] = { 285447, 285396, 285365, 285364, 285347, 285328, 285327, 285224, 285222, 285203, 285202, 285200, 285165, 285127, 285125, 285108, 285106, 285066, 285065, 285064, 285015, 285014, 285013, 285012, 285011, 285010, 285009, 285008 }; // [28]
    Int_t nRuns(28);
    nRuns = sizeof(run0)/sizeof(Int_t);

    if(effort==0) {
        alien->AddRunNumber(285365);
    } else if(effort==1) {
        for(Int_t iRun=0; iRun<nRuns; iRun++) {
            alien->AddRunNumber(run0[iRun]);
        }
    } else {
        printf("#####         No valid option: effort!         #####\n");
    }
}
