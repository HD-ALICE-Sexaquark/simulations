#pragma once

#include <array>
#include <cstdlib>
#include <vector>

#include <TH1.h>
#include <TList.h>
#include <TTree.h>

#include <AliAnalysisTaskSE.h>
#include <AliVVertex.h>

#include <AliESDEvent.h>

#include <AliMCEvent.h>
#include <AliMCEventHandler.h>
#include <AliMCParticle.h>

#include <AliPIDResponse.h>

class AliPIDResponse;

class AliAnalysisTaskSexaReactionQA final : public AliAnalysisTaskSE {
   public:
    AliAnalysisTaskSexaReactionQA(const char* name);
    AliAnalysisTaskSexaReactionQA();
    ~AliAnalysisTaskSexaReactionQA();

    AliAnalysisTaskSexaReactionQA(const AliAnalysisTaskSexaReactionQA&) = delete;
    AliAnalysisTaskSexaReactionQA& operator=(const AliAnalysisTaskSexaReactionQA&) = delete;
    AliAnalysisTaskSexaReactionQA(AliAnalysisTaskSexaReactionQA&&) = delete;
    AliAnalysisTaskSexaReactionQA& operator=(AliAnalysisTaskSexaReactionQA&&) = delete;

    // Main ~ executed at runtime //
    void UserCreateOutputObjects();
    void UserExec(Option_t* option);
    void Terminate(Option_t* option) {}

    // Tree //
    void CreateEventsBranches();
    void CreateMCBranches();
    void CreateTracksBranches();

    // Events //
    void ProcessEvent();

    // MC Particles //
    void ProcessMCParticles();
    void ReserveBranches_MC(size_t size);
    void ClearBranches_MC();

    // Tracks //
    void ProcessTracks();
    void ReserveBranches_Tracks(size_t size);
    void ClearBranches_Tracks();

   private:
    // AliRoot Objects //
    AliMCEvent* fMC;                      //! MC event
    const AliVVertex* fMC_PrimaryVertex;  //! MC gen. (or true) primary vertex
    AliESDEvent* fESD;                    //! reconstructed event
    const AliESDVertex* fPrimaryVertex;   //! primary vertex
    AliPIDResponse* fPIDResponse;         //! pid response object

    unsigned int fRunNumber;          //! run number
    unsigned int fEventNumberInFile;  //! event number
    float fCentrality;                //! centrality percentile
    float fMagneticField;             //! magnetic field

    // # Output //

    // QA Histograms //
    TList* fOutputHists;                   //!
    TH1D* fHist_Resolution_PV_X;           //!
    TH1D* fHist_Resolution_PV_Y;           //!
    TH1D* fHist_Resolution_PV_Z;           //!
    TH1D* fHist_TruePions_NSigmaPion;      //!
    TH1D* fHist_TrueKaons_NSigmaKaon;      //!
    TH1D* fHist_TrueProtons_NSigmaProton;  //!
    TH1D* fHist_MC_AntiSexaquark_Mass;     //!
    TH1D* fHist_MC_AntiSexaquark_Pt;       //!

    // Trees //
    TTree* fOutputTree;  //!
    // -- Event properties //
    float tEvent_PV_Xv;                        //!
    float tEvent_PV_Yv;                        //!
    float tEvent_PV_Zv;                        //!
    std::array<float, 6> tEvent_PV_CovMatrix;  //!
    int tEvent_NTracks;                        //!
    int tEvent_NTPCClusters;                   //!
    // -- MC Event properties //
    float tEvent_MC_PV_Xv;  //!
    float tEvent_MC_PV_Yv;  //!
    float tEvent_MC_PV_Zv;  //!
    // -- MC particles properties //
    std::vector<int> tMC_McIndex;         //!
    std::vector<int> tMC_PdgCode;         //!
    std::vector<int> tMC_Mother_McIndex;  //!
    std::vector<float> tMC_X;             //! origin x-vertex
    std::vector<float> tMC_Y;             //! origin y-vertex
    std::vector<float> tMC_Z;             //! origin z-vertex
    std::vector<float> tMC_Px;            //!
    std::vector<float> tMC_Py;            //!
    std::vector<float> tMC_Pz;            //!
    std::vector<float> tMC_E;             //!
    std::vector<int> tMC_Status;          //!
    std::vector<char> tMC_Generator;      //! 0: HIJING, 1: anti-neutron injector, 2: anti-sexaquark reaction
    std::vector<char> tMC_IsPrimary;      //!
    std::vector<char> tMC_IsSecFromMat;   //!
    std::vector<char> tMC_IsSecFromWeak;  //!
    // -- Tracks properties //
    std::vector<int> tTrack_EsdIndex;        //!
    std::vector<int> tTrack_McIndex;         //!
    std::vector<float> tTrack_X;             //! inner parametrization
    std::vector<float> tTrack_Y;             //! inner parametrization
    std::vector<float> tTrack_Z;             //! inner parametrization
    std::vector<float> tTrack_Px;            //! inner parametrization
    std::vector<float> tTrack_Py;            //! inner parametrization
    std::vector<float> tTrack_Pz;            //! inner parametrization
    std::vector<int> tTrack_Charge;          //!
    std::vector<float> tTrack_DCAxy;         //! pre-calculated DCA wrt PV
    std::vector<float> tTrack_DCAz;          //! pre-calculated DCA wrt PV
    std::vector<float> tTrack_NSigmaPion;    //!
    std::vector<float> tTrack_NSigmaKaon;    //!
    std::vector<float> tTrack_NSigmaProton;  //!
    std::vector<float> tTrack_SigmaX2;       //!
    std::vector<float> tTrack_SigmaXY;       //!
    std::vector<float> tTrack_SigmaY2;       //!
    std::vector<float> tTrack_SigmaXZ;       //!
    std::vector<float> tTrack_SigmaYZ;       //!
    std::vector<float> tTrack_SigmaZ2;       //!
    std::vector<float> tTrack_SigmaXPx;      //!
    std::vector<float> tTrack_SigmaYPx;      //!
    std::vector<float> tTrack_SigmaZPx;      //!
    std::vector<float> tTrack_SigmaPx2;      //!
    std::vector<float> tTrack_SigmaXPy;      //!
    std::vector<float> tTrack_SigmaYPy;      //!
    std::vector<float> tTrack_SigmaZPy;      //!
    std::vector<float> tTrack_SigmaPxPy;     //!
    std::vector<float> tTrack_SigmaPy2;      //!
    std::vector<float> tTrack_SigmaXPz;      //!
    std::vector<float> tTrack_SigmaYPz;      //!
    std::vector<float> tTrack_SigmaZPz;      //!
    std::vector<float> tTrack_SigmaPxPz;     //!
    std::vector<float> tTrack_SigmaPyPz;     //!
    std::vector<float> tTrack_SigmaPz2;      //!
    // extra tpc info //
    std::vector<float> tTrack_TPC_DCAxy;              //!
    std::vector<float> tTrack_TPC_DCAz;               //!
    std::vector<float> tTrack_TPC_NCrossedRows;       //!
    std::vector<float> tTrack_TPC_NClusters;          //!
    std::vector<float> tTrack_TPC_NClustersLC;        //! `LC` stands for the lowercase method `GetTPCncls`
    std::vector<float> tTrack_TPC_NClustersFound;     //!
    std::vector<float> tTrack_TPC_NClustersShared;    //!
    std::vector<float> tTrack_TPC_Chi2;               //!
    std::vector<float> tTrack_TPC_Chi2Constrained;    //! TPC-only chi2 at the primary vertex
    std::vector<float> tTrack_TPC_Signal;             //!
    std::vector<float> tTrack_TPC_SignalTunedOnData;  //! only available for MC
    std::vector<float> tTrack_TPC_SignalSigma;        //!
    std::vector<float> tTrack_TPC_SignalCorrected;    //!
    std::vector<float> tTrack_TPC_ESignalPion;        //! `E` stands for `Expected`
    std::vector<float> tTrack_TPC_ESigmaPion;         //!
    std::vector<float> tTrack_TPC_ESignalKaon;        //!
    std::vector<float> tTrack_TPC_ESigmaKaon;         //!
    std::vector<float> tTrack_TPC_ESignalProton;      //!
    std::vector<float> tTrack_TPC_ESigmaProton;       //!
    std::vector<int> tTrack_TPC_SignalN;              //!
    std::vector<int> tTrack_TPC_PointsFirst;          //!
    std::vector<int> tTrack_TPC_PointsIndexMax;       //!
    std::vector<int> tTrack_TPC_PointsLast;           //!
    std::vector<int> tTrack_TPC_PointsMaxDens;        //!

    ClassDef(AliAnalysisTaskSexaReactionQA, 1);
};
