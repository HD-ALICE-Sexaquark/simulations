#include <TChain.h>

#include <Math/Vector4D.h>
namespace RMath = ROOT::Math;

#include <AliAnalysisManager.h>
#include <AliAnalysisUtils.h>
#include <AliESDInputHandler.h>
#include <AliESDVertex.h>
#include <AliESDtrack.h>
#include <AliExternalTrackParam.h>
#include <AliInputEventHandler.h>
#include <AliLog.h>
#include <AliMultSelection.h>
#include <AliTPCPIDResponse.h>
#include <AliVEvent.h>

#include "AliAnalysisTaskSexaReactionQA.h"

ClassImp(AliAnalysisTaskSexaReactionQA);

// Constructor, called locally.
AliAnalysisTaskSexaReactionQA::AliAnalysisTaskSexaReactionQA(const char* name)
    : AliAnalysisTaskSE{name},
      //
      fMC{nullptr},
      fMC_PrimaryVertex{nullptr},
      fESD{nullptr},
      fPrimaryVertex{nullptr},
      fPIDResponse{nullptr},
      //
      fRunNumber{0},
      fEventNumberInFile{0},
      fCentrality{0.},
      fMagneticField{0.},
      //
      fOutputHists{nullptr},
      fHist_Resolution_PV_X{nullptr},
      fHist_Resolution_PV_Y{nullptr},
      fHist_Resolution_PV_Z{nullptr},
      fHist_TruePions_NSigmaPion{nullptr},
      fHist_TrueKaons_NSigmaKaon{nullptr},
      fHist_TrueProtons_NSigmaProton{nullptr},
      fHist_MC_AntiSexaquark_Mass{nullptr},
      fHist_MC_AntiSexaquark_Pt{nullptr},
      //
      fOutputTree{nullptr},
      //
      tEvent_PV_Xv{0.},
      tEvent_PV_Yv{0.},
      tEvent_PV_Zv{0.},
      tEvent_PV_CovMatrix{},
      tEvent_NTracks{0},
      tEvent_NTPCClusters{0},
      //
      tEvent_MC_PV_Xv{0.},
      tEvent_MC_PV_Yv{0.},
      tEvent_MC_PV_Zv{0.},
      //
      tMC_McIndex{},
      tMC_PdgCode{},
      tMC_Mother_McIndex{},
      tMC_X{},
      tMC_Y{},
      tMC_Z{},
      tMC_Px{},
      tMC_Py{},
      tMC_Pz{},
      tMC_E{},
      tMC_Status{},
      tMC_Generator{},
      tMC_IsPrimary{},
      tMC_IsSecFromMat{},
      tMC_IsSecFromWeak{},
      //
      tTrack_EsdIndex{},
      tTrack_McIndex{},
      tTrack_X{},
      tTrack_Y{},
      tTrack_Z{},
      tTrack_Px{},
      tTrack_Py{},
      tTrack_Pz{},
      tTrack_Charge{},
      tTrack_DCAxy{},
      tTrack_DCAz{},
      tTrack_NSigmaPion{},
      tTrack_NSigmaKaon{},
      tTrack_NSigmaProton{},
      tTrack_SigmaX2{},
      tTrack_SigmaXY{},
      tTrack_SigmaY2{},
      tTrack_SigmaXZ{},
      tTrack_SigmaYZ{},
      tTrack_SigmaZ2{},
      tTrack_SigmaXPx{},
      tTrack_SigmaYPx{},
      tTrack_SigmaZPx{},
      tTrack_SigmaPx2{},
      tTrack_SigmaXPy{},
      tTrack_SigmaYPy{},
      tTrack_SigmaZPy{},
      tTrack_SigmaPxPy{},
      tTrack_SigmaPy2{},
      tTrack_SigmaXPz{},
      tTrack_SigmaYPz{},
      tTrack_SigmaZPz{},
      tTrack_SigmaPxPz{},
      tTrack_SigmaPyPz{},
      tTrack_SigmaPz2{},
      tTrack_TPC_DCAxy{},
      tTrack_TPC_DCAz{},
      tTrack_TPC_NCrossedRows{},
      tTrack_TPC_NClusters{},
      tTrack_TPC_NClustersLC{},
      tTrack_TPC_NClustersFound{},
      tTrack_TPC_NClustersShared{},
      tTrack_TPC_Chi2{},
      tTrack_TPC_Chi2Constrained{},
      tTrack_TPC_Signal{},
      tTrack_TPC_SignalTunedOnData{},
      tTrack_TPC_SignalSigma{},
      tTrack_TPC_SignalCorrected{},
      tTrack_TPC_ESignalPion{},
      tTrack_TPC_ESigmaPion{},
      tTrack_TPC_ESignalKaon{},
      tTrack_TPC_ESigmaKaon{},
      tTrack_TPC_ESignalProton{},
      tTrack_TPC_ESigmaProton{},
      tTrack_TPC_SignalN{},
      tTrack_TPC_PointsFirst{},
      tTrack_TPC_PointsIndexMax{},
      tTrack_TPC_PointsLast{},
      tTrack_TPC_PointsMaxDens{} {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());  // fOutputHists
    DefineOutput(2, TTree::Class());  // fOutputTree
}

// Empty I/O constructor. Non-persistent members are initialized to their default values from here.
AliAnalysisTaskSexaReactionQA::AliAnalysisTaskSexaReactionQA() : AliAnalysisTaskSexaReactionQA{""} {}

// Destructor.
// NOTE: if `TList::SetOwner(true)` was called, the TList destructor should delete all objects added to it.
AliAnalysisTaskSexaReactionQA::~AliAnalysisTaskSexaReactionQA() {
    delete fOutputHists;
    delete fOutputTree;
}

// # Executed at Runtime # //

// Create output objects, called once at RUNTIME ~ execution on Grid.
void AliAnalysisTaskSexaReactionQA::UserCreateOutputObjects() {

    auto* man = AliAnalysisManager::GetAnalysisManager();
    if (man == nullptr) AliFatal("AliAnalysisManager couldn't be found.");

    auto* inputHandler = dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler());
    if (inputHandler == nullptr) AliFatal("AliESDInputHandler couldn't be found.");

    fPIDResponse = inputHandler->GetPIDResponse();

    // Prepare output list and histograms //

    fOutputHists = new TList();
    fOutputHists->SetOwner(true);

    fHist_Resolution_PV_X = new TH1D("Resolution_PV_X", ";Resolution_PV_X;Counts", 100, -10., 10.);
    fOutputHists->Add(fHist_Resolution_PV_X);
    fHist_Resolution_PV_Y = new TH1D("Resolution_PV_Y", ";Resolution_PV_Y;Counts", 100, -10., 10.);
    fOutputHists->Add(fHist_Resolution_PV_Y);
    fHist_Resolution_PV_Z = new TH1D("Resolution_PV_Z", ";Resolution_PV_Z;Counts", 100, -10., 10.);
    fOutputHists->Add(fHist_Resolution_PV_Z);
    fHist_TruePions_NSigmaPion = new TH1D("TruePions_NSigmaPion", ";TruePions_NSigmaPion;Counts", 100, -10., 10.);
    fOutputHists->Add(fHist_TruePions_NSigmaPion);
    fHist_TrueKaons_NSigmaKaon = new TH1D("TrueKaons_NSigmaKaon", ";TrueKaons_NSigmaKaon;Counts", 100, -10., 10.);
    fOutputHists->Add(fHist_TrueKaons_NSigmaKaon);
    fHist_TrueProtons_NSigmaProton = new TH1D("TrueProtons_NSigmaProton", ";TrueProtons_NSigmaProton;Counts", 100, -10., 10.);
    fOutputHists->Add(fHist_TrueProtons_NSigmaProton);
    fHist_MC_AntiSexaquark_Mass = new TH1D("MC_AntiSexaquark_Mass", ";MC_AntiSexaquark_Mass;Counts", 32, 1.59, 2.15);
    fOutputHists->Add(fHist_MC_AntiSexaquark_Mass);
    fHist_MC_AntiSexaquark_Pt = new TH1D("MC_AntiSexaquark_Pt", ";MC_AntiSexaquark_Pt;Counts", 100, 0., 10.);
    fOutputHists->Add(fHist_MC_AntiSexaquark_Pt);

    // Prepare output tree //

    fOutputTree = new TTree("Events", "Events");
    CreateEventsBranches();
    CreateMCBranches();
    CreateTracksBranches();

    // Post data //

    PostData(1, fOutputHists);
    PostData(2, fOutputTree);
}

// Main function, called per each event at RUNTIME ~ execution on Grid.
void AliAnalysisTaskSexaReactionQA::UserExec(Option_t* option) {

    // events //

    ProcessEvent();

    // mc particles //

    ProcessMCParticles();

    // tracks //

    ProcessTracks();

    // end of event //

    fOutputTree->Fill();
    ClearBranches_MC();
    ClearBranches_Tracks();

    // Post data //

    PostData(1, fOutputHists);
    PostData(2, fOutputTree);
}

// # Events # //

void AliAnalysisTaskSexaReactionQA::ProcessEvent() {

    // Load MC Event //

    fMC = MCEvent();
    if (fMC == nullptr) AliFatal("AliMCEvent couldn't be found.");
    fPIDResponse->SetCurrentMCEvent(fMC);  // enable PID tuned on data

    // Load Reconstructed Event, PV and Magnetic Field //

    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (fESD == nullptr) AliFatal("AliESDEvent couldn't be found.");

    // Assign metadata branches (1) //

    fRunNumber = fESD->GetRunNumber();
    fEventNumberInFile = fESD->GetEventNumberInFile();
    fPrimaryVertex = fESD->GetPrimaryVertex();

    // Assign more branches (2) //

    auto* MultSelection = dynamic_cast<AliMultSelection*>(fESD->FindListObject("MultSelection"));
    if (MultSelection == nullptr) AliFatal("AliMultSelection couldn't be found.");
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    fMagneticField = static_cast<float>(fESD->GetMagneticField());

    // Assign rest of branches (3) //

    tEvent_PV_Xv = static_cast<float>(fPrimaryVertex->GetX());
    tEvent_PV_Yv = static_cast<float>(fPrimaryVertex->GetY());
    tEvent_PV_Zv = static_cast<float>(fPrimaryVertex->GetZ());

    double PV_CovMatrix[6]{};
    fPrimaryVertex->GetCovarianceMatrix(PV_CovMatrix);
    for (int i = 0; i < 6; ++i) tEvent_PV_CovMatrix[i] = static_cast<float>(PV_CovMatrix[i]);

    tEvent_NTracks = fESD->GetNumberOfTracks();
    tEvent_NTPCClusters = fESD->GetNumberOfTPCClusters();

    // Assign MC branches (4) //

    fMC_PrimaryVertex = fMC->GetPrimaryVertex();
    tEvent_MC_PV_Xv = static_cast<float>(fMC_PrimaryVertex->GetX());
    tEvent_MC_PV_Yv = static_cast<float>(fMC_PrimaryVertex->GetY());
    tEvent_MC_PV_Zv = static_cast<float>(fMC_PrimaryVertex->GetZ());

    // Fill hists //

    fHist_Resolution_PV_X->Fill(tEvent_PV_Xv - tEvent_MC_PV_Xv);
    fHist_Resolution_PV_Y->Fill(tEvent_PV_Yv - tEvent_MC_PV_Yv);
    fHist_Resolution_PV_Z->Fill(tEvent_PV_Zv - tEvent_MC_PV_Zv);
}

// # Trees # //

// Add branches to `fOutputTree`.
void AliAnalysisTaskSexaReactionQA::CreateEventsBranches() {
    fOutputTree->Branch("RunNumber", &fRunNumber);
    fOutputTree->Branch("EventNumber", &fEventNumberInFile);
    fOutputTree->Branch("Centrality", &fCentrality);
    fOutputTree->Branch("MagneticField", &fMagneticField);
    fOutputTree->Branch("MC_PV_Xv", &tEvent_MC_PV_Xv);
    fOutputTree->Branch("MC_PV_Yv", &tEvent_MC_PV_Yv);
    fOutputTree->Branch("MC_PV_Zv", &tEvent_MC_PV_Zv);
    fOutputTree->Branch("PV_Xv", &tEvent_PV_Xv);
    fOutputTree->Branch("PV_Yv", &tEvent_PV_Yv);
    fOutputTree->Branch("PV_Zv", &tEvent_PV_Zv);
    fOutputTree->Branch("PV_CovMatrix", &tEvent_PV_CovMatrix, "PV_CovMatrix[6]/F");
    fOutputTree->Branch("NTracks", &tEvent_NTracks);
    fOutputTree->Branch("NTPCClusters", &tEvent_NTPCClusters);
}

// Add branches to the MC tree.
void AliAnalysisTaskSexaReactionQA::CreateMCBranches() {
    fOutputTree->Branch("MC_McIndex", &tMC_McIndex);
    fOutputTree->Branch("MC_PdgCode", &tMC_PdgCode);
    fOutputTree->Branch("MC_Mother_McIndex", &tMC_Mother_McIndex);
    fOutputTree->Branch("MC_X", &tMC_X);
    fOutputTree->Branch("MC_Y", &tMC_Y);
    fOutputTree->Branch("MC_Z", &tMC_Z);
    fOutputTree->Branch("MC_Px", &tMC_Px);
    fOutputTree->Branch("MC_Py", &tMC_Py);
    fOutputTree->Branch("MC_Pz", &tMC_Pz);
    fOutputTree->Branch("MC_E", &tMC_E);
    fOutputTree->Branch("MC_Status", &tMC_Status);
    fOutputTree->Branch("MC_Generator", &tMC_Generator);
    fOutputTree->Branch("MC_IsPrimary", &tMC_IsPrimary);
    fOutputTree->Branch("MC_IsSecFromMat", &tMC_IsSecFromMat);
    fOutputTree->Branch("MC_IsSecFromWeak", &tMC_IsSecFromWeak);
}

// Add branches to the tracks tree.
void AliAnalysisTaskSexaReactionQA::CreateTracksBranches() {
    fOutputTree->Branch("Track_EsdIndex", &tTrack_EsdIndex);
    fOutputTree->Branch("Track_McIndex", &tTrack_McIndex);
    fOutputTree->Branch("Track_X", &tTrack_X);
    fOutputTree->Branch("Track_Y", &tTrack_Y);
    fOutputTree->Branch("Track_Z", &tTrack_Z);
    fOutputTree->Branch("Track_Px", &tTrack_Px);
    fOutputTree->Branch("Track_Py", &tTrack_Py);
    fOutputTree->Branch("Track_Pz", &tTrack_Pz);
    fOutputTree->Branch("Track_Charge", &tTrack_Charge);
    fOutputTree->Branch("Track_DCAxy", &tTrack_DCAxy);
    fOutputTree->Branch("Track_DCAz", &tTrack_DCAz);
    fOutputTree->Branch("Track_NSigmaPion", &tTrack_NSigmaPion);
    fOutputTree->Branch("Track_NSigmaKaon", &tTrack_NSigmaKaon);
    fOutputTree->Branch("Track_NSigmaProton", &tTrack_NSigmaProton);
    fOutputTree->Branch("Track_SigmaX2", &tTrack_SigmaX2);
    fOutputTree->Branch("Track_SigmaXY", &tTrack_SigmaXY);
    fOutputTree->Branch("Track_SigmaY2", &tTrack_SigmaY2);
    fOutputTree->Branch("Track_SigmaXZ", &tTrack_SigmaXZ);
    fOutputTree->Branch("Track_SigmaYZ", &tTrack_SigmaYZ);
    fOutputTree->Branch("Track_SigmaZ2", &tTrack_SigmaZ2);
    fOutputTree->Branch("Track_SigmaXPx", &tTrack_SigmaXPx);
    fOutputTree->Branch("Track_SigmaYPx", &tTrack_SigmaYPx);
    fOutputTree->Branch("Track_SigmaZPx", &tTrack_SigmaZPx);
    fOutputTree->Branch("Track_SigmaPx2", &tTrack_SigmaPx2);
    fOutputTree->Branch("Track_SigmaXPy", &tTrack_SigmaXPy);
    fOutputTree->Branch("Track_SigmaYPy", &tTrack_SigmaYPy);
    fOutputTree->Branch("Track_SigmaZPy", &tTrack_SigmaZPy);
    fOutputTree->Branch("Track_SigmaPxPy", &tTrack_SigmaPxPy);
    fOutputTree->Branch("Track_SigmaPy2", &tTrack_SigmaPy2);
    fOutputTree->Branch("Track_SigmaXPz", &tTrack_SigmaXPz);
    fOutputTree->Branch("Track_SigmaYPz", &tTrack_SigmaYPz);
    fOutputTree->Branch("Track_SigmaZPz", &tTrack_SigmaZPz);
    fOutputTree->Branch("Track_SigmaPxPz", &tTrack_SigmaPxPz);
    fOutputTree->Branch("Track_SigmaPyPz", &tTrack_SigmaPyPz);
    fOutputTree->Branch("Track_SigmaPz2", &tTrack_SigmaPz2);
    fOutputTree->Branch("Track_TPC_DCAxy", &tTrack_TPC_DCAxy);
    fOutputTree->Branch("Track_TPC_DCAz", &tTrack_TPC_DCAz);
    fOutputTree->Branch("Track_TPC_NCrossedRows", &tTrack_TPC_NCrossedRows);
    fOutputTree->Branch("Track_TPC_NClusters", &tTrack_TPC_NClusters);
    fOutputTree->Branch("Track_TPC_NClustersLC", &tTrack_TPC_NClustersLC);
    fOutputTree->Branch("Track_TPC_NClustersFound", &tTrack_TPC_NClustersFound);
    fOutputTree->Branch("Track_TPC_NClustersShared", &tTrack_TPC_NClustersShared);
    fOutputTree->Branch("Track_TPC_Chi2", &tTrack_TPC_Chi2);
    fOutputTree->Branch("Track_TPC_Chi2Constrained", &tTrack_TPC_Chi2Constrained);
    fOutputTree->Branch("Track_TPC_Signal", &tTrack_TPC_Signal);
    fOutputTree->Branch("Track_TPC_SignalTunedOnData", &tTrack_TPC_SignalTunedOnData);
    fOutputTree->Branch("Track_TPC_SignalSigma", &tTrack_TPC_SignalSigma);
    fOutputTree->Branch("Track_TPC_SignalCorrected", &tTrack_TPC_SignalCorrected);
    fOutputTree->Branch("Track_TPC_ESignalPion", &tTrack_TPC_ESignalPion);
    fOutputTree->Branch("Track_TPC_ESigmaPion", &tTrack_TPC_ESigmaPion);
    fOutputTree->Branch("Track_TPC_ESignalKaon", &tTrack_TPC_ESignalKaon);
    fOutputTree->Branch("Track_TPC_ESigmaKaon", &tTrack_TPC_ESigmaKaon);
    fOutputTree->Branch("Track_TPC_ESignalProton", &tTrack_TPC_ESignalProton);
    fOutputTree->Branch("Track_TPC_ESigmaProton", &tTrack_TPC_ESigmaProton);
    fOutputTree->Branch("Track_TPC_SignalN", &tTrack_TPC_SignalN);
    fOutputTree->Branch("Track_TPC_PointsFirst", &tTrack_TPC_PointsFirst);
    fOutputTree->Branch("Track_TPC_PointsIndexMax", &tTrack_TPC_PointsIndexMax);
    fOutputTree->Branch("Track_TPC_PointsLast", &tTrack_TPC_PointsLast);
    fOutputTree->Branch("Track_TPC_PointsMaxDens", &tTrack_TPC_PointsMaxDens);
}

// # MC Generated # //

// Loop over MC particles in a single event.
void AliAnalysisTaskSexaReactionQA::ProcessMCParticles() {
    const int n_mc = fMC->GetNumberOfTracks();
    ReserveBranches_MC(n_mc);
    // read mc particles //
    for (int mc_idx = 0; mc_idx < n_mc; ++mc_idx) {
        auto* mcPart = dynamic_cast<AliMCParticle*>(fMC->GetTrack(mc_idx));
        if (mcPart == nullptr) continue;
        // -- fill branches
        tMC_McIndex.emplace_back(mc_idx);
        tMC_PdgCode.emplace_back(mcPart->PdgCode());
        tMC_Mother_McIndex.emplace_back(mcPart->GetMother());
        tMC_X.emplace_back(static_cast<float>(mcPart->Xv()));
        tMC_Y.emplace_back(static_cast<float>(mcPart->Yv()));
        tMC_Z.emplace_back(static_cast<float>(mcPart->Zv()));
        tMC_Px.emplace_back(static_cast<float>(mcPart->Px()));
        tMC_Py.emplace_back(static_cast<float>(mcPart->Py()));
        tMC_Pz.emplace_back(static_cast<float>(mcPart->Pz()));
        tMC_E.emplace_back(static_cast<float>(mcPart->E()));
        tMC_Status.emplace_back(static_cast<int>(mcPart->MCStatusCode()));
        tMC_Generator.emplace_back(static_cast<char>(mcPart->GetGeneratorIndex()));
        tMC_IsPrimary.emplace_back(static_cast<char>(mcPart->IsPhysicalPrimary()));
        tMC_IsSecFromMat.emplace_back(static_cast<char>(mcPart->IsSecondaryFromMaterial()));
        tMC_IsSecFromWeak.emplace_back(static_cast<char>(mcPart->IsSecondaryFromWeakDecay()));

        // Fill Histograms //

        if (mcPart->PdgCode() == -2112 && mcPart->MCStatusCode() >= 600 && mcPart->MCStatusCode() < 700 && mcPart->GetMother() == -1) {
            RMath::PxPyPzEVector anti_sexaquark(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->E());
            fHist_MC_AntiSexaquark_Mass->Fill(anti_sexaquark.M());
            fHist_MC_AntiSexaquark_Pt->Fill(anti_sexaquark.Pt());
        }
    }  // end of loop over MC particles
}

void AliAnalysisTaskSexaReactionQA::ReserveBranches_MC(size_t size) {
    tMC_McIndex.reserve(size);
    tMC_PdgCode.reserve(size);
    tMC_Mother_McIndex.reserve(size);
    tMC_X.reserve(size);
    tMC_Y.reserve(size);
    tMC_Z.reserve(size);
    tMC_Px.reserve(size);
    tMC_Py.reserve(size);
    tMC_Pz.reserve(size);
    tMC_E.reserve(size);
    tMC_Status.reserve(size);
    tMC_Generator.reserve(size);
    tMC_IsPrimary.reserve(size);
    tMC_IsSecFromMat.reserve(size);
    tMC_IsSecFromWeak.reserve(size);
}

// Clear MC branches.
void AliAnalysisTaskSexaReactionQA::ClearBranches_MC() {
    tMC_McIndex.clear();
    tMC_PdgCode.clear();
    tMC_Mother_McIndex.clear();
    tMC_X.clear();
    tMC_Y.clear();
    tMC_Z.clear();
    tMC_Px.clear();
    tMC_Py.clear();
    tMC_Pz.clear();
    tMC_E.clear();
    tMC_Status.clear();
    tMC_Generator.clear();
    tMC_IsPrimary.clear();
    tMC_IsSecFromMat.clear();
    tMC_IsSecFromWeak.clear();
}

// # Reconstructed # //

// Loop over the reconstructed tracks in a single event.
void AliAnalysisTaskSexaReactionQA::ProcessTracks() {
    const int n_tracks = fESD->GetNumberOfTracks();
    ReserveBranches_Tracks(n_tracks);
    for (int esd_idx = 0; esd_idx < n_tracks; ++esd_idx) {
        // get track //
        auto* track = fESD->GetTrack(esd_idx);
        // get info //
        const auto* inner_param = track->GetInnerParam();  // NOTE: already protected in `PassesTrackSelection`
        if (inner_param == nullptr) continue;
        double position[3]{};
        inner_param->GetXYZ(position);
        double momentum[3]{};
        inner_param->GetPxPyPz(momentum);
        float dca[2]{};
        float dca_cov[3]{};
        track->GetImpactParameters(dca, dca_cov);
        double cov_xyz_pxpypz[21]{};
        track->GetCovarianceXYZPxPyPz(cov_xyz_pxpypz);
        // fill branches //
        tTrack_EsdIndex.emplace_back(esd_idx);
        int mc_idx = std::abs(track->GetLabel());
        tTrack_McIndex.emplace_back(mc_idx);
        tTrack_X.emplace_back(static_cast<float>(position[0]));
        tTrack_Y.emplace_back(static_cast<float>(position[1]));
        tTrack_Z.emplace_back(static_cast<float>(position[2]));
        tTrack_Px.emplace_back(static_cast<float>(momentum[0]));
        tTrack_Py.emplace_back(static_cast<float>(momentum[1]));
        tTrack_Pz.emplace_back(static_cast<float>(momentum[2]));
        tTrack_Charge.emplace_back(inner_param->Charge());
        // -- dca //
        tTrack_DCAxy.emplace_back(dca[0]);
        tTrack_DCAz.emplace_back(dca[1]);
        // -- pid //
        float n_sigmas_pion = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kPion);      // NOTE: includes corrections
        float n_sigmas_kaon = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kKaon);      // NOTE: includes corrections
        float n_sigmas_proton = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kProton);  // NOTE: includes corrections
        tTrack_NSigmaPion.emplace_back(n_sigmas_pion);
        tTrack_NSigmaKaon.emplace_back(n_sigmas_kaon);
        tTrack_NSigmaProton.emplace_back(n_sigmas_proton);
        // -- cov. matrix //
        tTrack_SigmaX2.emplace_back(static_cast<float>(cov_xyz_pxpypz[0]));
        tTrack_SigmaXY.emplace_back(static_cast<float>(cov_xyz_pxpypz[1]));
        tTrack_SigmaY2.emplace_back(static_cast<float>(cov_xyz_pxpypz[2]));
        tTrack_SigmaXZ.emplace_back(static_cast<float>(cov_xyz_pxpypz[3]));
        tTrack_SigmaYZ.emplace_back(static_cast<float>(cov_xyz_pxpypz[4]));
        tTrack_SigmaZ2.emplace_back(static_cast<float>(cov_xyz_pxpypz[5]));
        tTrack_SigmaXPx.emplace_back(static_cast<float>(cov_xyz_pxpypz[6]));
        tTrack_SigmaYPx.emplace_back(static_cast<float>(cov_xyz_pxpypz[7]));
        tTrack_SigmaZPx.emplace_back(static_cast<float>(cov_xyz_pxpypz[8]));
        tTrack_SigmaPx2.emplace_back(static_cast<float>(cov_xyz_pxpypz[9]));
        tTrack_SigmaXPy.emplace_back(static_cast<float>(cov_xyz_pxpypz[10]));
        tTrack_SigmaYPy.emplace_back(static_cast<float>(cov_xyz_pxpypz[11]));
        tTrack_SigmaZPy.emplace_back(static_cast<float>(cov_xyz_pxpypz[12]));
        tTrack_SigmaPxPy.emplace_back(static_cast<float>(cov_xyz_pxpypz[13]));
        tTrack_SigmaPy2.emplace_back(static_cast<float>(cov_xyz_pxpypz[14]));
        tTrack_SigmaXPz.emplace_back(static_cast<float>(cov_xyz_pxpypz[15]));
        tTrack_SigmaYPz.emplace_back(static_cast<float>(cov_xyz_pxpypz[16]));
        tTrack_SigmaZPz.emplace_back(static_cast<float>(cov_xyz_pxpypz[17]));
        tTrack_SigmaPxPz.emplace_back(static_cast<float>(cov_xyz_pxpypz[18]));
        tTrack_SigmaPyPz.emplace_back(static_cast<float>(cov_xyz_pxpypz[19]));
        tTrack_SigmaPz2.emplace_back(static_cast<float>(cov_xyz_pxpypz[20]));
        // -- get tpc dca
        float tpc_dca_xy = 0.;
        float tpc_dca_z = 0.;
        track->GetImpactParametersTPC(tpc_dca_xy, tpc_dca_z);
        // -- fill tpc branches (1)
        tTrack_TPC_DCAxy.emplace_back(tpc_dca_xy);
        tTrack_TPC_DCAz.emplace_back(tpc_dca_z);
        tTrack_TPC_NCrossedRows.emplace_back(track->GetTPCCrossedRows());
        tTrack_TPC_NClusters.emplace_back(track->GetTPCNcls());
        tTrack_TPC_NClustersLC.emplace_back(track->GetTPCncls());
        tTrack_TPC_NClustersFound.emplace_back(track->GetTPCNclsF());
        tTrack_TPC_NClustersShared.emplace_back(track->GetTPCnclsS());
        tTrack_TPC_Chi2.emplace_back(static_cast<float>(track->GetTPCchi2()));
        tTrack_TPC_Chi2Constrained.emplace_back(static_cast<float>(track->GetConstrainedChi2TPC()));
        tTrack_TPC_Signal.emplace_back(track->GetTPCsignal());  // NOTE: uncorrected
        tTrack_TPC_SignalTunedOnData.emplace_back(static_cast<float>(track->GetTPCsignalTunedOnData()));
        tTrack_TPC_SignalSigma.emplace_back(track->GetTPCsignalSigma());
        // -- get corrected tpc pid info
        double bethe_pion = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC, track, AliPID::kPion);
        double e_sigma_pion = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC, track, AliPID::kPion);
        double bethe_kaon = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC, track, AliPID::kKaon);
        double e_sigma_kaon = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC, track, AliPID::kKaon);
        double bethe_proton = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC, track, AliPID::kProton);
        double e_sigma_proton = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC, track, AliPID::kProton);
        //    -- pion
        tTrack_TPC_SignalCorrected.emplace_back(e_sigma_pion * n_sigmas_pion + bethe_pion);
        tTrack_TPC_ESignalPion.emplace_back(bethe_pion);
        tTrack_TPC_ESigmaPion.emplace_back(e_sigma_pion);
        //    -- kaon
        tTrack_TPC_ESignalKaon.emplace_back(bethe_kaon);
        tTrack_TPC_ESigmaKaon.emplace_back(e_sigma_kaon);
        //    -- proton
        tTrack_TPC_ESignalProton.emplace_back(bethe_proton);
        tTrack_TPC_ESigmaProton.emplace_back(e_sigma_proton);
        // -- fill tpc branches (2)
        tTrack_TPC_SignalN.emplace_back(track->GetTPCsignalN());
        tTrack_TPC_PointsFirst.emplace_back(track->GetTPCPoints(0));
        tTrack_TPC_PointsIndexMax.emplace_back(track->GetTPCPoints(1));
        tTrack_TPC_PointsLast.emplace_back(track->GetTPCPoints(2));
        tTrack_TPC_PointsMaxDens.emplace_back(track->GetTPCPoints(3));

        // Fill Histograms //

        auto* mcPart = dynamic_cast<AliMCParticle*>(fMC->GetTrack(mc_idx));

        if (std::abs(mcPart->PdgCode()) == 211) fHist_TruePions_NSigmaPion->Fill(n_sigmas_pion);
        if (std::abs(mcPart->PdgCode()) == 321) fHist_TrueKaons_NSigmaKaon->Fill(n_sigmas_kaon);
        if (std::abs(mcPart->PdgCode()) == 2212) fHist_TrueProtons_NSigmaProton->Fill(n_sigmas_proton);
    }  // end of loop over tracks
}

// Reserve the branches of the reconstructed tracks.
void AliAnalysisTaskSexaReactionQA::ReserveBranches_Tracks(size_t size) {
    tTrack_EsdIndex.reserve(size);
    tTrack_McIndex.reserve(size);
    tTrack_X.reserve(size);
    tTrack_Y.reserve(size);
    tTrack_Z.reserve(size);
    tTrack_Px.reserve(size);
    tTrack_Py.reserve(size);
    tTrack_Pz.reserve(size);
    tTrack_Charge.reserve(size);
    tTrack_DCAxy.reserve(size);
    tTrack_DCAz.reserve(size);
    tTrack_NSigmaPion.reserve(size);
    tTrack_NSigmaKaon.reserve(size);
    tTrack_NSigmaProton.reserve(size);
    tTrack_SigmaX2.reserve(size);
    tTrack_SigmaXY.reserve(size);
    tTrack_SigmaY2.reserve(size);
    tTrack_SigmaXZ.reserve(size);
    tTrack_SigmaYZ.reserve(size);
    tTrack_SigmaZ2.reserve(size);
    tTrack_SigmaXPx.reserve(size);
    tTrack_SigmaYPx.reserve(size);
    tTrack_SigmaZPx.reserve(size);
    tTrack_SigmaPx2.reserve(size);
    tTrack_SigmaXPy.reserve(size);
    tTrack_SigmaYPy.reserve(size);
    tTrack_SigmaZPy.reserve(size);
    tTrack_SigmaPxPy.reserve(size);
    tTrack_SigmaPy2.reserve(size);
    tTrack_SigmaXPz.reserve(size);
    tTrack_SigmaYPz.reserve(size);
    tTrack_SigmaZPz.reserve(size);
    tTrack_SigmaPxPz.reserve(size);
    tTrack_SigmaPyPz.reserve(size);
    tTrack_SigmaPz2.reserve(size);
    tTrack_TPC_DCAxy.reserve(size);
    tTrack_TPC_DCAz.reserve(size);
    tTrack_TPC_NCrossedRows.reserve(size);
    tTrack_TPC_NClusters.reserve(size);
    tTrack_TPC_NClustersLC.reserve(size);
    tTrack_TPC_NClustersFound.reserve(size);
    tTrack_TPC_NClustersShared.reserve(size);
    tTrack_TPC_Chi2.reserve(size);
    tTrack_TPC_Chi2Constrained.reserve(size);
    tTrack_TPC_Signal.reserve(size);
    tTrack_TPC_SignalTunedOnData.reserve(size);
    tTrack_TPC_SignalSigma.reserve(size);
    tTrack_TPC_SignalCorrected.reserve(size);
    tTrack_TPC_ESignalPion.reserve(size);
    tTrack_TPC_ESigmaPion.reserve(size);
    tTrack_TPC_ESignalKaon.reserve(size);
    tTrack_TPC_ESigmaKaon.reserve(size);
    tTrack_TPC_ESignalProton.reserve(size);
    tTrack_TPC_ESigmaProton.reserve(size);
    tTrack_TPC_SignalN.reserve(size);
    tTrack_TPC_PointsFirst.reserve(size);
    tTrack_TPC_PointsIndexMax.reserve(size);
    tTrack_TPC_PointsLast.reserve(size);
    tTrack_TPC_PointsMaxDens.reserve(size);
}

// Clear the branches of the reconstructed tracks.
void AliAnalysisTaskSexaReactionQA::ClearBranches_Tracks() {
    tTrack_EsdIndex.clear();
    tTrack_McIndex.clear();
    tTrack_X.clear();
    tTrack_Y.clear();
    tTrack_Z.clear();
    tTrack_Px.clear();
    tTrack_Py.clear();
    tTrack_Pz.clear();
    tTrack_Charge.clear();
    tTrack_DCAxy.clear();
    tTrack_DCAz.clear();
    tTrack_NSigmaPion.clear();
    tTrack_NSigmaKaon.clear();
    tTrack_NSigmaProton.clear();
    tTrack_SigmaX2.clear();
    tTrack_SigmaXY.clear();
    tTrack_SigmaY2.clear();
    tTrack_SigmaXZ.clear();
    tTrack_SigmaYZ.clear();
    tTrack_SigmaZ2.clear();
    tTrack_SigmaXPx.clear();
    tTrack_SigmaYPx.clear();
    tTrack_SigmaZPx.clear();
    tTrack_SigmaPx2.clear();
    tTrack_SigmaXPy.clear();
    tTrack_SigmaYPy.clear();
    tTrack_SigmaZPy.clear();
    tTrack_SigmaPxPy.clear();
    tTrack_SigmaPy2.clear();
    tTrack_SigmaXPz.clear();
    tTrack_SigmaYPz.clear();
    tTrack_SigmaZPz.clear();
    tTrack_SigmaPxPz.clear();
    tTrack_SigmaPyPz.clear();
    tTrack_SigmaPz2.clear();
    tTrack_TPC_DCAxy.clear();
    tTrack_TPC_DCAz.clear();
    tTrack_TPC_NCrossedRows.clear();
    tTrack_TPC_NClusters.clear();
    tTrack_TPC_NClustersLC.clear();
    tTrack_TPC_NClustersFound.clear();
    tTrack_TPC_NClustersShared.clear();
    tTrack_TPC_Chi2.clear();
    tTrack_TPC_Chi2Constrained.clear();
    tTrack_TPC_Signal.clear();
    tTrack_TPC_SignalTunedOnData.clear();
    tTrack_TPC_SignalSigma.clear();
    tTrack_TPC_SignalCorrected.clear();
    tTrack_TPC_ESignalPion.clear();
    tTrack_TPC_ESigmaPion.clear();
    tTrack_TPC_ESignalKaon.clear();
    tTrack_TPC_ESigmaKaon.clear();
    tTrack_TPC_ESignalProton.clear();
    tTrack_TPC_ESigmaProton.clear();
    tTrack_TPC_SignalN.clear();
    tTrack_TPC_PointsFirst.clear();
    tTrack_TPC_PointsIndexMax.clear();
    tTrack_TPC_PointsLast.clear();
    tTrack_TPC_PointsMaxDens.clear();
}
