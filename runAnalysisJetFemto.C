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
#include "AliFemtoDreamTrackCuts.h"
#endif
//===========================================================================================================================================

//Definitions
#define ALIPHYSICS_VER  "vAN-20220531_ROOT6-1"
#define GRIDWorkingDir  "JET_FEMTOSCOPY_WORKING_DIR"
#define AnalysisMacro   "AnalysisJetFemtoscopy"
#define AnalysisTask    "AliAnalysisTaskJetFemto"
#define TTL             5000
#define nRunsPerMaster  1

//Functions
AliAnalysisGrid *CreateAlienHandler (const char *mode, Bool_t merge, int effort );
void LoadAnalysisTask (AliAnalysisManager *mgr, TString CentEst, bool isRun1, AliFemtoDreamEventCuts *evtCuts, AliFemtoDreamCollConfig *config, AliFemtoDreamTrackCuts *fTrackCutsProton, AliFemtoDreamTrackCuts *fTrackCutsAntiproton, int effort);
void EventHandler     (AliAnalysisManager *mgr);
void LoadPIDResponse();
void SetInputRuns(AliAnalysisAlien *alien, const char *mode, int effort);

//Test
Bool_t test = (kFALSE);

//___________________________________________________________________________________________________________________________________________
void LoadAnalysisTask (AliAnalysisManager *mgr, TString CentEst, bool isRun1, AliFemtoDreamEventCuts *evtCuts, AliFemtoDreamCollConfig *config, AliFemtoDreamTrackCuts *fTrackCutsProton, AliFemtoDreamTrackCuts *fTrackCutsAntiproton, int effort)  {

    //Run Mode
    Bool_t runData      = (kTRUE);
    bool isMC           = (false);

    //Load Analysis Task
    gROOT->LoadMacro("AliAnalysisTaskJetFemto.cxx+g");

    //Input container
    AliAnalysisDataContainer *input = mgr -> GetCommonInputContainer();

    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();

    //Analysis Task
    AliAnalysisTaskJetFemto *task = new AliAnalysisTaskJetFemto ("task_JetReconstruction_Femtoscopy", isMC);

    if(CentEst == "kInt7"){// INT7 = BIT(1)
        if(isRun1||(effort==0)) {
            task->SetTrigger(AliVEvent::kMB);
            task->SelectCollisionCandidates(AliVEvent::kMB);
            std::cout << "Added kMB Trigger \n";
        } else {
            task->SetTrigger(AliVEvent::kINT7);
            task->SelectCollisionCandidates(AliVEvent::kINT7);
            std::cout << "Added kINT7 Trigger \n";
        }
    } else if (CentEst == "kHM") {
        task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
        std::cout << "Added kHighMultV0 Trigger \n";
    }else{
        std::cout << "=====================================================================" << std::endl;
        std::cout << "=====================================================================" << std::endl;
        std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
        std::cout << "=====================================================================" << std::endl;
        std::cout << "=====================================================================" << std::endl;
    }

    /* task->SetEventCuts(evtCuts);
    task->SetTrackCutsProton(fTrackCutsProton);
    task->SetTrackCutsAntiproton(fTrackCutsAntiproton);
    task->SetCollectionConfig(config);
    task->SetIsMC(isMC); */

    task -> AliAnalysisTaskJetFemto::SetRunningMode(runData);
    mgr -> AddTask(task);
    mgr -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput(task,1,mgr->CreateContainer("Jets", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr -> ConnectOutput(task,2,mgr->CreateContainer("Femtoscopy", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr -> ConnectOutput(task,3,mgr->CreateContainer("QAPlots", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
}


//###########################################################################################################################################
//###########################################################################################################################################
//                                                            Analysis Task
//###########################################################################################################################################
//###########################################################################################################################################


void runAnalysisJetFemto (const char *mode="full", Bool_t merge=kTRUE )  {
    bool isMC               = false;
    bool isRun1             = false;
    bool oldreject          = false;
    bool MCtemplatefit      = false;
    bool doSharedCut        = false;
    float fSpherDown        = 0.7;
    float fdPhidEta         = 0.01;
    TString CentEst         = "kINT7";
    const char *cutVar      = "0";
    int effort              = 1;
    


    //Grid Connection
    TGrid::Connect("alien://");

    //Alien Handler
    AliAnalysisGrid *alienGridHandler = CreateAlienHandler (mode,merge,effort);
    if (!alienGridHandler) return;

    //Analysis Manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager");
    mgr->SetGridHandler(alienGridHandler);

    //===============================================================================================================
    // FEMTOSCOPY
    AliFemtoDreamEventCuts *evtCuts = new AliFemtoDreamEventCuts;
    evtCuts->CleanUpMult(false,false,false,true);
    if(isRun1 && oldreject) {
        evtCuts->SetCutMinContrib(3);
        evtCuts->SetZVtxPosition(-8., 8.);
    } else {
        evtCuts->SetCutMinContrib(2);
        evtCuts->SetZVtxPosition(-10., 10.);
    }
    evtCuts->SetSphericityCuts(fSpherDown, 1.0, 0.5);

    AliFemtoDreamTrackCuts *fTrackCutsProton = new AliFemtoDreamTrackCuts();
    fTrackCutsProton->SetIsMonteCarlo(isMC);
    fTrackCutsProton->SetCutCharge(1);
    fTrackCutsProton->SetPtRange(0.14, 4.0);
    fTrackCutsProton->SetEtaRange(-0.8, 0.8);
    fTrackCutsProton->SetNClsTPC(80);
    fTrackCutsProton->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
    if ( !MCtemplatefit ) {
        fTrackCutsProton->SetFilterBit(96); // Filterbit 5+6
        fTrackCutsProton->SetDCAVtxZ(0.3);
        fTrackCutsProton->SetDCAVtxXY(0.3);
    } else {
        fTrackCutsProton->SetFilterBit(128); // Filterbit 7
    }
    if ( doSharedCut ) { fTrackCutsProton->SetCutSharedCls(true);}
    fTrackCutsProton->SetNClsTPC(80); // In Indico + additional Chi²/NDF <4
    fTrackCutsProton->SetPID(AliPID::kProton, 0.5);
    fTrackCutsProton->SetRejLowPtPionsTOF(true);
    fTrackCutsProton->SetMinimalBooking(false);
    fTrackCutsProton->SetPlotDCADist(true);

    if ( isMC && MCtemplatefit ) {
        //fTrackCutsProton->SetPlotContrib(true);
        fTrackCutsProton->CheckParticleMothers(true);
        fTrackCutsProton->SetPlotDCADist(true);
        //fTrackCutsProton->SetOriginMultiplicityHists(true);
        fTrackCutsProton->SetFillQALater(false); //Be careful about this flag! When the MinimalBooking is set
    }

    AliFemtoDreamTrackCuts *fTrackCutsAntiproton = new AliFemtoDreamTrackCuts();
    fTrackCutsAntiproton->SetIsMonteCarlo(isMC);
    fTrackCutsAntiproton->SetCutCharge(-1);
    fTrackCutsAntiproton->SetPtRange(0.14, 4.0);
    fTrackCutsAntiproton->SetEtaRange(-0.8, 0.8);
    fTrackCutsAntiproton->SetNClsTPC(80);
    fTrackCutsAntiproton->SetDCAReCalculation(true);
    if ( !MCtemplatefit ) {
        fTrackCutsAntiproton->SetFilterBit(96);
        fTrackCutsAntiproton->SetDCAVtxZ(0.3);
        fTrackCutsAntiproton->SetDCAVtxXY(0.3);
    } else {
        fTrackCutsAntiproton->SetFilterBit(128);
    }
    if ( doSharedCut ) { fTrackCutsAntiproton->SetCutSharedCls(true);}
    fTrackCutsAntiproton->SetNClsTPC(80);
    fTrackCutsAntiproton->SetPID(AliPID::kProton, 0.5);
    fTrackCutsAntiproton->SetRejLowPtPionsTOF(true); // check if func exists
    fTrackCutsAntiproton->SetMinimalBooking(false);
    fTrackCutsAntiproton->SetPlotDCADist(true);

    if ( isMC && MCtemplatefit ) {
        //fTrackCutsAntiproton->SetPlotContrib(true);
        fTrackCutsAntiproton->CheckParticleMothers(true);
        fTrackCutsAntiproton->SetPlotDCADist(true);
        //fTrackCutsAntiproton->SetOriginMultiplicityHists(true);
        fTrackCutsAntiproton->SetFillQALater(false); //Be careful about this flag! When the MinimalBooking is set
    }
    /* std::vector<int> PDGParticles;
    PDGParticles.push_back(2212); // p+
    PDGParticles.push_back(-2212); // p-

    std::vector<float> ZVtxBins;
    if(isRun1) ZVtxBins.push_back(-10);
    ZVtxBins.push_back(-8);
    ZVtxBins.push_back(-6);
    ZVtxBins.push_back(-4);
    ZVtxBins.push_back(-2);
    ZVtxBins.push_back(0);
    ZVtxBins.push_back(2);
    ZVtxBins.push_back(4);
    ZVtxBins.push_back(6);
    ZVtxBins.push_back(8);
    if(isRun1) ZVtxBins.push_back(10);
    //The Multiplicity bins are set here
    std::vector<int> MultBins;
    MultBins.push_back(0);
    MultBins.push_back(18);
    MultBins.push_back(30);

    std::vector<int> NBins;
    NBins.push_back(750);
    NBins.push_back(750);
    NBins.push_back(750);
    std::vector<float> kMin;
    //minimum k* value
    kMin.push_back(0.);
    kMin.push_back(0.);
    kMin.push_back(0.);
    //maximum k* value
    std::vector<float> kMax;
    kMax.push_back(3.);
    kMax.push_back(3.);
    kMax.push_back(3.);
    //pair rejection ---> CHECK AND ADJUST
    std::vector<bool> closeRejection;
    closeRejection.push_back(true); // pi+ pi+
    closeRejection.push_back(false); // pi+ pi-
    closeRejection.push_back(true); // pi- pi-

    // if (suffix == "5") {
    //     //Deactivate the ClosePairRejection
    //     fdPhidEta=0.;
    //     closeRejection.clear();
    //     closeRejection.push_back(false); // pi+ pi+
    //     closeRejection.push_back(false); // pi+ pi-
    //     closeRejection.push_back(false); // pi- pi-
    // }

    //QA plots for tracks
    std::vector<int> pairQA;
    pairQA.push_back(11); // pi+ pi+
    pairQA.push_back(11); // pi+ pi-
    pairQA.push_back(11); // pi- pi- */

    AliFemtoDreamCollConfig *config=new AliFemtoDreamCollConfig("Femto","Femto");
    /* config->SetZBins(ZVtxBins);
    config->SetMultBins(MultBins);
    config->SetMultBinning(true);
    config->SetPDGCodes(PDGParticles);
    config->SetNBinsHist(NBins);
    config->SetMinKRel(kMin);
    config->SetMaxKRel(kMax);
    config->SetClosePairRejection(closeRejection);
    config->SetDeltaEtaMax(fdPhidEta); // https://alice-notes.web.cern.ch/system/files/notes/analysis/616/2018-08-10-NotepPb.pdf
    config->SetDeltaPhiMax(fdPhidEta);
    config->SetExtendedQAPairs(pairQA);
    //Here we set the mixing depth.
    config->SetMixingDepth(10); // AN
    config->SetkTBinning(true);
    config->SetmTBinning(true);
    config->SetMinimalBookingME(false);
    config->SetdPhidEtaPlots(true);
    config->SetkTandMultBinning(true);
    config->SetkTandMultPtBinning(true);
    if(isMC) config->SetkTandMultMCTrueBinning(true);
    config->SetdPhidEtaPlotsSmallK(true);
    config->SetPhiEtaBinnign(true);
    
    if (isMC) {
        config->SetMomentumResolution(true);
        } else {
        std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
        }
        if (isMC) {
        config->SetPhiEtaBinnign(true);
        } else {
        std::cout << "You are trying to request the Eta Phi Plots without MC Info; fix it wont work! \n";
    } */

    //===============================================================================================================

    //Event Handler, PID & Centrality
    
    EventHandler(mgr);
    LoadPIDResponse();

    //Analysis Task
    LoadAnalysisTask(mgr, CentEst, isRun1, evtCuts, config, fTrackCutsProton, fTrackCutsAntiproton, effort);

    //Start analysis
    if (!mgr->InitAnalysis())  return;
    mgr->PrintStatus();
    mgr->StartAnalysis("grid");
};
//___________________________________________________________________________________________________________________________________________
AliAnalysisGrid *CreateAlienHandler (const char *mode, Bool_t merge, int effort )  {

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
    alien->SetGridDataDir("/alice/data/2018/LHC18b");
    alien->SetRunPrefix("000");
    alien->SetDataPattern("*pass2/AOD/*AliAOD.root");
    SetInputRuns(alien, mode, effort);
    alien->SetNrunsPerMaster(nRunsPerMaster);
    alien->SetGridWorkingDir(Form("%s",GRIDWorkingDir));
    alien->SetGridOutputDir("OUTPUT");
    alien->SetAnalysisSource(Form("%s.cxx",AnalysisTask));
    alien->SetAdditionalLibs(Form("%s.cxx  %s.h",AnalysisTask,AnalysisTask));
    alien->SetMergeViaJDL(merge);
    alien->SetMaxMergeStages(2);
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
void SetInputRuns (AliAnalysisAlien *alien, const char *mode, int effort) {
    Int_t run0[] = { 285447, 285396, 285365, 285364, 285347, 285328, 285327, 285224, 285222, 285203, 285202, 285200, 285165, 285127, 285125, 285108, 285106, 285066, 285065, 285064, 285015, 285014, 285013, 285012, 285011, 285010, 285009, 285008 };
    Int_t nRuns(0);
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
