// Generator of Sexaquark-nucleon interactions

#include <set>
#include <vector>

#include "TDatabasePDG.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPDGCode.h"

#include "AliConst.h"
#include "AliGenEventHeader.h"
#include "AliGenSexaquarkReaction.h"
#include "AliLog.h"
#include "AliPDG.h"
#include "AliRun.h"

#define VERBOSE_MODE 1
#define PDG_SEXAQUARK -900000020
#define PDG_GEANTINO 0
#define MASS_SEXAQUARK 1.8

// macro to integrate class for interactive use in ROOT
ClassImp(AliGenSexaquarkReaction);

//_____________________________________________________________________________
AliGenSexaquarkReaction::AliGenSexaquarkReaction()
    : AliGenerator(), fNucleonPDG(0), fDecayChannelPDG(0), fRadiusMin(0), fRadiusMax(0), fPtDistribution(0) {
    //
    // Default constructor:
    // - Set all kinematic ranges
    // - As no specific nucleon or decay channel was provided, all reactions will be produced randomly
    //
    fName = "AliGenSexaquarkReaction";
    fTitle = "Generator of AntiSexaquark-Nucleon Reactions";
    InitReactionsInfo();
}

//_____________________________________________________________________________
AliGenSexaquarkReaction::~AliGenSexaquarkReaction() {
    //
    // Default destructor
    //
    delete fPtDistribution;
}

//_____________________________________________________________________________
void AliGenSexaquarkReaction::Init() {
    //
    // Initialisation, check consistency of selected ranges
    //
    if (TestBit(kPtRange) && TestBit(kMomentumRange)) {
        Fatal("Init", "You should not set the momentum range and the pt range!\n");
    }
    if ((!TestBit(kPtRange)) && (!TestBit(kMomentumRange))) {
        Fatal("Init", "You should set either the momentum or the pt range!\n");
    }
    if ((TestBit(kYRange) && TestBit(kThetaRange)) || (TestBit(kYRange) && TestBit(kEtaRange)) ||
        (TestBit(kEtaRange) && TestBit(kThetaRange))) {
        Fatal("Init", "You should only set the range of one of these variables: y, eta or theta\n");
    }
    if ((!TestBit(kYRange)) && (!TestBit(kEtaRange)) && (!TestBit(kThetaRange))) {
        Fatal("Init", "You should set the range of one of these variables: y, eta or theta\n");
    }
#if VERBOSE_MODE
    // print information
    PrintTitle();
    PrintParameters();
#endif
}

//_____________________________________________________________________________
void AliGenSexaquarkReaction::InitReactionsInfo() {
    //
    // Prepare the allowed decay channels to the set of vectors
    //
    kSexaquarkNeutronProducts = {{-3122, 310},             // AntiLambda, K0
                                 {-3122, 321, -211, 111},  // AntiLambda, K+, Pi-, Pi0
                                 {-3122, 310, -211, 211},  // AntiLambda, K0, Pi-, Pi+
                                 {-3122, 310, 111, 111},   // AntiLambda, K0, Pi-, Pi+
                                 {-2212, 310, 310, 211},   // AntiProton, K0, K0, Pi+
                                 {-2212, 310, 321, 111},   // AntiProton, K0, K+, Pi0
                                 {-3312, -211}};           // XiBaryon+, Pi-
    kSexaquarkProtonProducts = {{-3122, 321},              // AntiLambda, K+
                                {-3122, 321, -211, 211},   // AntiLambda, K+, Pi-, Pi+
                                {-3122, 321, 111, 111},    // AntiLambda, K+, Pi0, Pi0
                                {-3122, 310, 211, 111},    // AntiLambda, K0, Pi+, Pi0
                                {-2212, 321, 321, 111},    // AntiProton, K+, K+, Pi0
                                {-2212, 321, 310, 211}};   // AntiProton, K+, K0, Pi+
}

//_____________________________________________________________________________
void AliGenSexaquarkReaction::SetReaction(Int_t NucleonPDG, std::vector<Int_t> DecayChannelPDG) {
    //
    //
    //
    // 1) Prepare nucleon
    if (NucleonPDG == 2112 || NucleonPDG == 2212) {
        fNucleonPDG = NucleonPDG;
    } else {
        // print error and exit
        AliFatal("ERROR: input nucleon PDG should be 2112 (neutron) or 2212 (proton).");
    }
    // 2) Check if the input decay channel is in the sets
    if (fNucleonPDG == 2112 && (Int_t)kSexaquarkNeutronProducts.count(DecayChannelPDG) > 0) {
        fDecayChannelPDG = DecayChannelPDG;
    } else if (fNucleonPDG == 2212 && (Int_t)kSexaquarkProtonProducts.count(DecayChannelPDG) > 0) {
        fDecayChannelPDG = DecayChannelPDG;
    } else {
        // print error and exit
        AliFatal("ERROR: input decay channel is not allowed by QCD.");
    }
}

//_____________________________________________________________________________
void AliGenSexaquarkReaction::Generate() {
    //
    // Generate one Sexaquark-nucleon interaction, and push event to MC Stack
    //

    // 1) Prepare random numbers
    GetRandom()->SetSeed(0);
    Float_t fRandom[4];
    Rndm(fRandom, 4);

    // 3) Set Sexaquark momentum
    std::vector<Float_t> myIntegral;
    std::vector<Float_t> myAlpha;
    std::vector<Float_t> myBeta;
    std::vector<Float_t> myGamma;
    Float_t Pt_S;
    if (fPtDistribution) {
        Pt_S = GetRandomFromFunction(myIntegral, myAlpha, myBeta, myGamma, fRandom[0]);
    } else {
        Pt_S = fPtMin + fRandom[0] * (fPtMax - fPtMin);
    }
    // Float_t M_S = TDatabasePDG::Instance()->GetParticle(PDG_SEXAQUARK)->Mass();
    Float_t M_S = MASS_SEXAQUARK;
    Float_t Mt_S = TMath::Sqrt(M_S * M_S + Pt_S * Pt_S);                     // transverse mass (derived from inv. mass and Pt)
    Float_t Y_S = fYMin + fRandom[1] * (fYMax - fYMin);                      // rapidity (uniform distribution)
    Float_t Phi_S = fPhiMin + fRandom[2] * (fPhiMax - fPhiMin);              // azimuthal angle (uniform distribution) (in radians)
    Float_t Radius_S = fRadiusMin + fRandom[3] * (fRadiusMax - fRadiusMin);  // radius (uniform distribution) (in cm)
    Float_t Px_S = Pt_S * TMath::Cos(Phi_S);                                 // Px
    Float_t Py_S = Pt_S * TMath::Sin(Phi_S);                                 // Py
    Float_t Pl_S = Mt_S * TMath::SinH(Y_S);                                  // longitudinal momentum = Pz (derived from trans. mass and Y)
    Float_t E_S = Mt_S * TMath::CosH(Y_S);                                   // energy (derived from trans. mass and Y)
    TLorentzVector Sexaquark(Px_S, Py_S, Pl_S, E_S);

    AliInfo(" > Generated AntiSexaquark:");
    AliInfoF("   P        = %.4f", Sexaquark.P());
    AliInfoF("   Px,Py,Pz = (%.4f, %.4f, %.4f)", Sexaquark.Px(), Sexaquark.Py(), Sexaquark.Pz());
    AliInfoF("   E        = %.4f", Sexaquark.E());
    AliInfoF("   M        = %.4f", Sexaquark.M());
    AliInfoF("   Pt       = %.4f", Pt_S);
    AliInfoF("   Mt       = %.4f", Mt_S);
    AliInfoF("   Y        = %.4f", Y_S);
    AliInfoF("   Theta    = %.4f", Sexaquark.Theta());
    AliInfoF("   Phi      = %.4f", Phi_S);
    AliInfoF("   Radius   = %.4f", Radius_S);

    // 4) Set nucleon momentum, at rest (no nuclear Fermi motion)
    TLorentzVector Nucleon(0., 0., 0., TDatabasePDG::Instance()->GetParticle(fNucleonPDG)->Mass());

    // 5) Set primary vertex
    TVector3 PrimaryVertex(0., 0., 0.);
    Float_t PrimaryTime = 0.;

    // 6) Set secondary vertex
    // let's start the transformation, we already have
    // - radius_S: 2D Radius -- radius in cylindrical coordinates
    // - phi_S: phi -- the same in both spherical and cylindrical coordinates
    Float_t Theta_S = Sexaquark.Theta();  // polar angle (in radians)
    Float_t Vz_S = Radius_S / TMath::Tan(Theta_S); // z, in both cartesian and cylindrical coordinates
    Float_t Radius3D_S = TMath::Sqrt(TMath::Power(Radius_S, 2) + TMath::Power(Vz_S, 2)); // 3D Radius -- radius in spherical coordinates
    Float_t Vx_S = Radius3D_S * TMath::Sin(Theta_S) * TMath::Cos(Phi_S); // x, in cartesian coordinates
    Float_t Vy_S = Radius3D_S * TMath::Sin(Theta_S) * TMath::Sin(Phi_S); // y, in cartesian coordinates
    TVector3 SecondaryVertex(Vx_S, Vy_S, Vz_S);

    Float_t Beta_S = Sexaquark.P() / Sexaquark.E();
    Float_t SecondaryTime = SecondaryVertex.Mag() / (Beta_S * TMath::Ccgs());

    // 7) Create array of masses (must be a double)
    Double_t DecayChannelMasses[(Int_t)fDecayChannelPDG.size()];
    for (Int_t i = 0; i < (Int_t)fDecayChannelPDG.size(); i++) {
        DecayChannelMasses[i] = TDatabasePDG::Instance()->GetParticle(fDecayChannelPDG[i])->Mass();
    }

    // ROOT class to generate n-body event with constant cross-section
    TGenPhaseSpace Generator;

    // Randomly determine the momenta of the reaction products
    TLorentzVector W = Sexaquark + Nucleon;
    Generator.SetDecay(W, (Int_t)fDecayChannelPDG.size(), &DecayChannelMasses[0]);

    Float_t weight = (Float_t)Generator.Generate();

    Int_t TrackCounter;

#if VERBOSE_MODE
    // verbose
    AliInfo("   Vertex            R        Theta          Phi");
    AliInfoF("  Primary %12.7f %12.7f %12.7f", PrimaryVertex.Perp(), PrimaryVertex.Theta(), PrimaryVertex.Phi());
    AliInfoF("Secondary %12.7f %12.7f %12.7f", SecondaryVertex.Perp(), SecondaryVertex.Theta(), SecondaryVertex.Phi());
    AliInfo("");
    AliInfo("  I           Id           Px           Py           Pz            E            M");
    // AliInfoF("%3i %12i %12.7f %12.7f %12.7f %12.7f %12.7f", -1, PDG_SEXAQUARK, Sexaquark.Px(), Sexaquark.Py(), Sexaquark.Pz(),
    //  Sexaquark.E(), Sexaquark.M());
    // AliInfoF("%3i %12i %12.7f %12.7f %12.7f %12.7f %12.7f", -1, fNucleonPDG, Nucleon.Px(), Nucleon.Py(), Nucleon.Pz(), Nucleon.E(),
    //  Nucleon.M());
#endif

    // add sexaquark to the stack
    /*
    PushTrack(fTrackIt,           // done: track final state particles
              -1,                 // parent: -1 = primary
              PDG_GEANTINO,       //
              Sexaquark.Px(),     // px: Px of the final-state particle (in GeV/c)
              Sexaquark.Py(),     // py: Py of the final-state particle (in GeV/c)
              Sexaquark.Pz(),     // pz: Pz of the final-state particle (in GeV/c)
              Sexaquark.E(),      // E: energy of the final-state particle
              PrimaryVertex.X(),  // vx: x-component of position of interaction vertex
              PrimaryVertex.Y(),  // vy: y-component of position of interaction vertex
              PrimaryVertex.Z(),  // vz: z-component of position of interaction vertex
              PrimaryTime,        // vt: time of interaction vertex (in seconds)
              0.,                 // polx: x-component of polarisation vector, 0 by default
              0.,                 // poly: y-component of polarisation vector, 0 by default
              0.,                 // polz: z-component of polarisation vector, 0 by default
              kPNoProcess,        // mech: production mechanism
              TrackCounter,       // counter of tracks (incremented internally by AliStack)
              1);                 // weight of the event
    */

    // get 4-momentum vectors of the reaction products and put them on the Geant stack
    for (Int_t i = 0; i < (Int_t)fDecayChannelPDG.size(); i++) {

        // add particle to the stack
        PushTrack(fTrackIt,                     // done: track final state particles
                  -1,                           // parent: the idea is that is the sexaquark
                                                // (for the moment, it's -1 = primary)
                  fDecayChannelPDG[i],          // pdg of the final-state particle
                  Generator.GetDecay(i)->Px(),  // px: Px of the final-state particle (in GeV/c)
                  Generator.GetDecay(i)->Py(),  // py: Py of the final-state particle (in GeV/c)
                  Generator.GetDecay(i)->Pz(),  // pz: Pz of the final-state particle (in GeV/c)
                  Generator.GetDecay(i)->E(),   // E: energy of the final-state particle
                  SecondaryVertex.X(),          // vx: x-component of position of interaction vertex
                  SecondaryVertex.Y(),          // vy: y-component of position of interaction vertex
                  SecondaryVertex.Z(),          // vz: z-component of position of interaction vertex
                  SecondaryTime - PrimaryTime,  // vt: time of interaction vertex (in seconds)
                  0.,                           // polx: x-component of polarisation vector, 0 by default
                  0.,                           // poly: y-component of polarisation vector, 0 by default
                  0.,                           // polz: z-component of polarisation vector, 0 by default
                  kPHInhelastic,                // mech: production mechanism (Hadronic Inelastic Scattering)
                  TrackCounter,                 // counter of tracks (incremented internally by AliStack)
                  weight,                       // weight of the event
                  6);                           // is: generation status code
                                                // (0: unstable, 1: stable, 6: from sexaquark)
        //     1 * (TMath::Abs(fDecayChannelPDG[i]) == 211 || TMath::Abs(fDecayChannelPDG[i]) == 2212 ||
        //            TMath::Abs(fDecayChannelPDG[i]) == 321)
#if VERBOSE_MODE
        // verbose
        AliInfoF("%3i %12i %12.7f %12.7f %12.7f %12.7f %12.7f", TrackCounter, fDecayChannelPDG[i], Generator.GetDecay(i)->Px(),
                 Generator.GetDecay(i)->Py(), Generator.GetDecay(i)->Pz(), Generator.GetDecay(i)->E(), Generator.GetDecay(i)->M());
#endif
    }

    // unknown things with AliGenEventHeader and Containers and gAlice ???
    TArrayF auxSecondaryVertex(3);
    auxSecondaryVertex[0] = SecondaryVertex.X();
    auxSecondaryVertex[1] = SecondaryVertex.Y();
    auxSecondaryVertex[2] = SecondaryVertex.Z();

    AliGenEventHeader *header = new AliGenEventHeader("BOX");
    header->SetPrimaryVertex(auxSecondaryVertex);          // the real vertex of the interaction!
    header->SetNProduced((Int_t)fDecayChannelPDG.size());  // amount of particles produced?
    header->SetInteractionTime(SecondaryTime - PrimaryTime);

    // passes header either to the container or to gAlice
    if (fContainer) {
        header->SetName(fName);
        fContainer->AddHeader(header);
    } else {
        gAlice->SetGenEventHeader(header);
    }

#if VERBOSE_MODE
    AliInfo("");
#endif
}

//_____________________________________________________________________________
void AliGenSexaquarkReaction::SetPtDistribution(TF1 &PtDistribution) {
    //
    // Clone the pT distribution into a new address that's owned by the class
    //
    TF1 *tmp = new TF1(PtDistribution);
    fPtDistribution = tmp;
    if (fPtDistribution) {
        SetPtRange(fPtDistribution->GetXmin(), fPtDistribution->GetXmax());
    }
}

//_____________________________________________________________________________
void AliGenSexaquarkReaction::GenerateN(Int_t N) {
    //
    // Generate N sexaquark-nucleon reactions
    //
    for (Int_t i = 0; i < N; i++) {
#if VERBOSE_MODE
        AliInfo("");
        AliInfo("\tEvent listing (summary)");
        AliInfo("");
#endif
        Generate();
    }
}

//_____________________________________________________________________________
void AliGenSexaquarkReaction::PrintParameters() {
    //
    // Print the parameters of the current instance
    //
    AliInfo("> Starting new AntiSexaquark-Nucleon Reaction with the following parameters:");
    AliInfo("");
    AliInfoF("  Chosen nucleon : %i", fNucleonPDG);
    TString DecayChannel_str = "{";
    for (Int_t p = 0; p < (Int_t)fDecayChannelPDG.size(); p++) {
        DecayChannel_str += Form("%i", fDecayChannelPDG[p]);
        if (p + 1 != (Int_t)fDecayChannelPDG.size()) DecayChannel_str += ", ";
    }
    DecayChannel_str += "}";
    AliInfoF("  Chosen decay   : %s", DecayChannel_str.Data());
    AliInfo("");
    AliInfo("> Kinematical ranges:");
    AliInfo("");
    if (fPtDistribution) {
        AliInfoF("  Pt : %s", fPtDistribution->GetExpFormula().Data());
    }
    AliInfoF("  %f < Pt < %f", fPtMin, fPtMax);
    AliInfoF("  %f < Phi < %f", fPhiMin, fPhiMax);
    AliInfoF("  %f < Y < %f", fYMin, fYMax);
    AliInfoF("  %f < Radius < %f", fRadiusMin, fRadiusMax);
    AliInfo("");
}

//_____________________________________________________________________________
void AliGenSexaquarkReaction::PrintTitle() {
    //
    // A small ASCII divertiment
    //
    std::cout << "     ___    ___ ______          _____                                        __   ____                  __  _            "
              << std::endl;
    std::cout
        << "    /   |  / (_) ____/__  ____ / ___/___  _  ______ _____ ___  ______ ______/ /__/ __ \\___  ____ ______/ /_(_)___  ____  "
        << std::endl;
    std::cout << "   / /| | / / / / __/ _ \\/ __ \\\\__ \\/ _ \\| |/_/ __ `/ __ `/ / / / __ `/ ___/ //_/ /_/ / _ \\/ __ `/ ___/ __/ / __ "
                 "\\/ __ \\ "
              << std::endl;
    std::cout << "  / ___ |/ / / /_/ /  __/ / / /__/ /  __/>  </ /_/ / /_/ / /_/ / /_/ / /  / ,< / _, _/  __/ /_/ / /__/ /_/ / /_/ / / / / "
              << std::endl;
    std::cout << " /_/  |_/_/_/\\____/\\___/_/ /_/____/\\___/_/|_|\\__,_/\\__, /\\__,_/\\__,_/_/  /_/|_/_/ "
                 "|_|\\___/\\__,_/\\___/\\__/_/\\____/_/ /_/  "
              << std::endl;
    std::cout << "                                                     /_/                                                                 "
              << std::endl;
    std::cout << std::endl;
}

//_____________________________________________________________________________
Bool_t AliGenSexaquarkReaction::ComputeCdfTable(std::vector<Float_t> &Integral_vec, std::vector<Float_t> &Alpha_vec,
                                                std::vector<Float_t> &Beta_vec, std::vector<Float_t> &Gamma_vec) {
    //
    //
    //
    if (!fPtDistribution) {
        return kFALSE;
    }

    Int_t Npx = fPtDistribution->GetNpx();

    Integral_vec.resize(Npx + 1);
    Alpha_vec.resize(Npx + 1);
    Beta_vec.resize(Npx);
    Gamma_vec.resize(Npx);
    Integral_vec[0] = 0;
    Alpha_vec[Npx] = 0;

    Float_t integ;
    Int_t intNegative = 0;
    Int_t i;
    Bool_t logbin = kFALSE;
    Float_t dx;
    Float_t xmin = fPtDistribution->GetXmin();
    Float_t xmax = fPtDistribution->GetXmax();
    dx = (xmax - xmin) / Npx;

    std::vector<Float_t> xx(Npx + 1);
    for (i = 0; i < Npx; i++) {
        xx[i] = xmin + i * dx;
    }
    xx[Npx] = xmax;
    for (i = 0; i < Npx; i++) {
        if (logbin) {
            integ = fPtDistribution->Integral(TMath::Power(10, xx[i]), TMath::Power(10, xx[i + 1]), 0.0);
        } else {
            integ = fPtDistribution->Integral(xx[i], xx[i + 1], 0.0);
        }
        if (integ < 0) {
            intNegative++;
            integ = -integ;
        }
        Integral_vec[i + 1] = Integral_vec[i] + integ;
    }
    if (intNegative > 0) {
        Warning("GetRandom", "function:%s has %d negative values: abs assumed", fPtDistribution->GetName(), intNegative);
    }
    if (Integral_vec[Npx] == 0) {
        Error("GetRandom", "Integral of function is zero");
        return kFALSE;
    }
    Float_t total = Integral_vec[Npx];
    for (i = 1; i <= Npx; i++) {  // normalize integral to 1
        Integral_vec[i] /= total;
    }
    // the integral r for each bin is approximated by a parabola
    //  x = alpha + beta*r +gamma*r**2
    // compute the coefficients alpha, beta, gamma for each bin
    Float_t x0, r1, r2, r3;
    for (i = 0; i < Npx; i++) {
        x0 = xx[i];
        r2 = Integral_vec[i + 1] - Integral_vec[i];
        if (logbin) {
            r1 = fPtDistribution->Integral(TMath::Power(10, x0), TMath::Power(10, x0 + 0.5 * dx), 0.0) / total;
        } else {
            r1 = fPtDistribution->Integral(x0, x0 + 0.5 * dx, 0.0) / total;
        }
        r3 = 2 * r2 - 4 * r1;
        if (TMath::Abs(r3) > 1e-8) {
            Gamma_vec[i] = r3 / (dx * dx);
        } else {
            Gamma_vec[i] = 0;
        }
        Beta_vec[i] = r2 / dx - Gamma_vec[i] * dx;
        Alpha_vec[i] = x0;
        Gamma_vec[i] *= 2;
    }
    return kTRUE;
}

//_____________________________________________________________________________
Float_t AliGenSexaquarkReaction::GetRandomFromFunction(std::vector<Float_t> &Integral_vec, std::vector<Float_t> &Alpha_vec,
                                                       std::vector<Float_t> &Beta_vec, std::vector<Float_t> &Gamma_vec, Float_t r) {
    //
    // r: random number
    //

    if (!fPtDistribution) {
        return 0;
    }

    Int_t Npx = fPtDistribution->GetNpx();

    //  Check if integral array must be build
    if (Integral_vec.size() == 0) {
        Bool_t ret = ComputeCdfTable(Integral_vec, Alpha_vec, Beta_vec, Gamma_vec);
        if (!ret) {
            return TMath::QuietNaN();
        }
    }

    // return random number
    Int_t bin = TMath::BinarySearch(Npx, Integral_vec.data(), r);
    Float_t rr = r - Integral_vec[bin];

    Float_t yy;
    if (Gamma_vec[bin] != 0) {
        yy = (-Beta_vec[bin] + TMath::Sqrt(Beta_vec[bin] * Beta_vec[bin] + 2 * Gamma_vec[bin] * rr)) / Gamma_vec[bin];
    } else {
        yy = rr / Beta_vec[bin];
    }
    Float_t x = Alpha_vec[bin] + yy;
    if (Alpha_vec[Npx] > 0) {
        return TMath::Power(10, x);
    }
    return x;
}
