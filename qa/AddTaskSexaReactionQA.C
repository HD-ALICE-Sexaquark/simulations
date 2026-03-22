#include <AliAnalysisManager.h>

#include "AliAnalysisTaskSexaReactionQA.h"

AliAnalysisTaskSexaReactionQA *AddTaskSexaReactionQA() {

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (mgr == nullptr) return nullptr;

    auto *task = new AliAnalysisTaskSexaReactionQA("AliAnalysisTaskSexaReactionQA");
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer("Hists", TList::Class(), mgr->kOutputContainer, mgr->GetCommonFileName()));
    mgr->ConnectOutput(task, 2, mgr->CreateContainer("Events", TTree::Class(), mgr->kOutputContainer, mgr->GetCommonFileName()));

    return task;
}
