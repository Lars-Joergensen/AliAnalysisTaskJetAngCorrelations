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
    Bool_t runData      = (kTRUE);

    //Load Analysis Task
    gROOT->LoadMacro("AliAnalysisTaskJetFemto.cxx+g");

    //Input container
    AliAnalysisDataContainer *input = mgr -> GetCommonInputContainer();

    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();

    //Analysis Task
    AliAnalysisTaskJetFemto *task = new AliAnalysisTaskJetFemto ("task_JetReconstruction_Femtoscopy");
    task -> SelectCollisionCandidates (AliVEvent::kMB); // removed in working base but throws execution errors if removed here; kMB = minimum bias, kINT7 = BIT(1)
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

    /* AliFemtoDreamEventCuts *evtCuts = new AliFemtoDreamEventCuts;
    evtCuts->CleanUpMult(false,false,false,true);
    if(isRun1 && oldreject) { // add bool
        evtCuts->SetCutMinContrib(3);
        evtCuts->SetZVtxPosition(-8., 8.);
    } else {
        evtCuts->SetCutMinContrib(2);
        evtCuts->SetZVtxPosition(-10., 10.);
    }
    evtCuts->SetSphericityCuts(fSpherDown, 1.0, 0.5);

    AliFemtoDreamTrackCuts *fTrackCutsPosPion=new AliFemtoDreamTrackCuts();
    fTrackCutsPosPion->SetIsMonteCarlo(isMC);
    fTrackCutsPosPion->SetCutCharge(1);
    fTrackCutsPosPion->SetPtRange(0.14, 4.0);
    fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
    fTrackCutsPosPion->SetNClsTPC(80);

    fTrackCutsPosPion->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
    if ( !MCtemplatefit ) { // add bool
        fTrackCutsPosPion->SetFilterBit(96); // Filterbit 5+6
        fTrackCutsPosPion->SetDCAVtxZ(0.3);
        fTrackCutsPosPion->SetDCAVtxXY(0.3);
    } else {
        fTrackCutsPosPion->SetFilterBit(128); // Filterbit 7
    }
    if ( doSharedCut ) { fTrackCutsPosPion->SetCutSharedCls(true);} // add bool
    fTrackCutsPosPion->SetNClsTPC(80); // In Indico + additional Chi²/NDF <4
    fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
    fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
    fTrackCutsPosPion->SetMinimalBooking(false);
    fTrackCutsPosPion->SetPlotDCADist(true);

    if ( isMC && MCtemplatefit ) {
        //fTrackCutsPosPion->SetPlotContrib(true);
        fTrackCutsPosPion->CheckParticleMothers(true);
        fTrackCutsPosPion->SetPlotDCADist(true);
        //fTrackCutsPosPion->SetOriginMultiplicityHists(true);
        fTrackCutsPosPion->SetFillQALater(false); //Be careful about this flag! When the MinimalBooking is set
    }

    AliFemtoDreamTrackCuts *fTrackCutsNegPion=new AliFemtoDreamTrackCuts();
    fTrackCutsNegPion->SetIsMonteCarlo(isMC);
    fTrackCutsNegPion->SetCutCharge(-1);
    fTrackCutsNegPion->SetPtRange(0.14, 4.0);
    fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
    fTrackCutsNegPion->SetNClsTPC(80);
    fTrackCutsNegPion->SetDCAReCalculation(true);
    if ( !MCtemplatefit ) {
        fTrackCutsNegPion->SetFilterBit(96);
        fTrackCutsNegPion->SetDCAVtxZ(0.3);
        fTrackCutsNegPion->SetDCAVtxXY(0.3);
    } else {
        fTrackCutsNegPion->SetFilterBit(128);
    }
    if ( doSharedCut ) { fTrackCutsNegPion->SetCutSharedCls(true);}
    fTrackCutsNegPion->SetNClsTPC(80);
    fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
    fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
    fTrackCutsNegPion->SetMinimalBooking(false);
    fTrackCutsNegPion->SetPlotDCADist(true);

    if ( isMC && MCtemplatefit ) {
        //fTrackCutsNegPion->SetPlotContrib(true);
        fTrackCutsNegPion->CheckParticleMothers(true);
        fTrackCutsNegPion->SetPlotDCADist(true);
        //fTrackCutsNegPion->SetOriginMultiplicityHists(true);
        fTrackCutsNegPion->SetFillQALater(false); //Be careful about this flag! When the MinimalBooking is set
    }
    std::vector<int> PDGParticles;
    PDGParticles.push_back(211); // pi+
    PDGParticles.push_back(-211); // pi-

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
    //pair rejection
    std::vector<bool> closeRejection;
    closeRejection.push_back(true); // pi+ pi+
    closeRejection.push_back(false); // pi+ pi-
    closeRejection.push_back(true); // pi- pi-

    if (suffix == "5") {
        //Deactivate the ClosePairRejection
        fdPhidEta=0.;
        closeRejection.clear();
        closeRejection.push_back(false); // pi+ pi+
        closeRejection.push_back(false); // pi+ pi-
        closeRejection.push_back(false); // pi- pi-
    }

    //QA plots for tracks
    std::vector<int> pairQA;
    pairQA.push_back(11); // pi+ pi+
    pairQA.push_back(11); // pi+ pi-
    pairQA.push_back(11); // pi- pi-

    //QA plots for tracks
    std::vector<int> pairQA;
    pairQA.push_back(11); // pi+ pi+
    pairQA.push_back(11); // pi+ pi-
    pairQA.push_back(11); // pi- pi-

    AliFemtoDreamCollConfig *config=new AliFemtoDreamCollConfig("Femto","Femto");
    config->SetZBins(ZVtxBins);
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
    config->SetPhiEtaBinnign(true); */

//===============================================================================================================

    //Event Handler, PID & Centrality
    EventHandler(mgr);
    LoadPIDResponse();

    /* if(CentEst == "kInt7"){
        if(isRun1) {
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

    task->SetEventCuts(evtCuts);
    task->SetTrackCutsPosPion(fTrackCutsPosPion);
    task->SetTrackCutsNegPion(fTrackCutsNegPion);
    task->SetCollectionConfig(config);
    task->SetIsMC(isMC); */

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