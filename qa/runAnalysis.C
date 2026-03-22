#include <TChain.h>
#include <TString.h>
#include <TSystemDirectory.h>

#include <AliAnalysisAlien.h>
#include <AliAnalysisManager.h>
#include <AliDataFile.h>
#include <AliESDInputHandler.h>
#include <AliMCEventHandler.h>

#include <AliAnalysisTaskPIDResponse.h>
#include <AliMultSelectionTask.h>
#include <AliPhysicsSelectionTask.h>

#include "AliAnalysisTaskSexaReactionQA.h"

void runAnalysis(const TString &InputPath,           // what comes before the RN
                 int RunNumber,                      // single run number
                 int NDirs = 1,                      // number of subdirs per run
                 long long Local_LimitToNEvents = 0  // 0 means all events
) {

    // # Start # //

    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");

    auto *mgr = new AliAnalysisManager("SexaReactionQA");

    std::cout << "INFO  !! runAnalysis.C !! Created AliAnalysisManager" << '\n';

    // # Input Handlers # //

    auto *esdH = new AliESDInputHandler();
    esdH->SetNeedField();  // necessary to get GoldenChi2
    mgr->SetInputEventHandler(esdH);

    auto *mcH = new AliMCEventHandler();
    mcH->SetReadTR(false);
    mgr->SetMCtruthEventHandler(mcH);

    std::cout << "INFO  !! runAnalysis.C !! Passed creation of input handlers" << '\n';

    // # Add Helper Tasks # //

    // Reference:
    // https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsPhysSel
    bool applyPileupCuts = false;
    TString TaskPhysicsSelection_Options = TString::Format("(%i, %i)", (int)true, (int)applyPileupCuts);
    auto *TaskPhysicsSelection = reinterpret_cast<AliPhysicsSelectionTask *>(
        gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C" + TaskPhysicsSelection_Options));
    if (TaskPhysicsSelection == nullptr) return;

    TString TaskCentrality_Options = "";  // nothing
    auto *TaskCentrality = reinterpret_cast<AliMultSelectionTask *>(
        gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C" + TaskCentrality_Options));
    if (TaskCentrality == nullptr) return;
    TString ProductionName = "LHC15o";
    if (RunNumber >= 295581 && RunNumber <= 296689) ProductionName = "LHC18q";
    if (RunNumber >= 296690 && RunNumber <= 300000) ProductionName = "LHC18r";
    TaskCentrality->SetAlternateOADBforEstimators(ProductionName + "-DefaultMC-HIJING_V0fix");

    // References:
    // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PIDInAnalysis
    // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/TPCSplines
    // https://alisw.github.io/git-advanced/#how-to-use-large-data-files-for-analysis
    bool pid_is_mc = true;
    bool pid_auto_mc_esd = true;
    bool pid_tune_on_data = true;
    TString pid_reco_pass = "";
    bool pid_cache_pid = false;
    TString pid_path_oadb =
        AliDataFile::GetFileNameOADB("COMMON/PID/data/TPCPIDResponseOADB_pileupCorr.root");  // NOTE: large storage file not in common repos
    TString pid_path_eta_maps = "$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCetaMaps_pileupCorr.root";
    TString pid_det_response = TString::Format("TPC-OADB:%s;TPC-Maps:%s", pid_path_oadb.Data(), pid_path_eta_maps.Data());
    TString TaskPIDResponse_Options = TString::Format("(%i, %i, %i, \"%s\", %i, \"%s\")", (int)pid_is_mc, (int)pid_auto_mc_esd, (int)pid_tune_on_data,
                                                      pid_reco_pass.Data(), pid_cache_pid, pid_det_response.Data());
    auto *TaskPIDResponse = reinterpret_cast<AliAnalysisTaskPIDResponse *>(
        gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C" + TaskPIDResponse_Options));
    if (TaskPIDResponse == nullptr) return;

    std::cout << "INFO  !! runAnalysis.C !! Passed addition of helper tasks" << '\n';

    // # Add Main Task # //

    gInterpreter->LoadMacro("AliAnalysisTaskSexaReactionQA.cxx++g");

    auto *TaskSexaReactionQA = reinterpret_cast<AliAnalysisTaskSexaReactionQA *>(gInterpreter->ExecuteMacro("AddTaskSexaReactionQA.C"));
    if (TaskSexaReactionQA == nullptr) return;

    std::cout << "INFO  !! runAnalysis.C !! Passed addition of main task" << '\n';

    // # Init Analysis Manager # //

    mgr->SetDebugLevel(0);
    if (!mgr->InitAnalysis()) return;

    std::cout << "INFO  !! runAnalysis.C !! Passed InitAnalysis" << '\n';

    // # Start Analysis # //

    TChain *chain = new TChain("esdTree");
    TString FilePath = "";
    for (int DN = 1; DN <= NDirs; ++DN) {
        FilePath = TString::Format("%s/%i/%03i/AliESDs.root", InputPath.Data(), RunNumber, DN);
        std::cout << "INFO  !! runAnalysis.C !! Adding file " << FilePath << '\n';
        chain->AddFile(FilePath);
    }

    if (Local_LimitToNEvents == 0) {
        mgr->StartAnalysis("local", chain);  // read all events
    } else {
        mgr->StartAnalysis("local", chain, Local_LimitToNEvents);  // read first NEvents
    }

    std::cout << "INFO  !! runAnalysis.C !! Passed StartAnalysis" << '\n';
}
