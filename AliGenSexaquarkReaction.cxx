#include <cmath>

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include <Math/Functor.h>
#include <Math/Integrator.h>
namespace RMath = ROOT::Math;

#include "AliGenEventHeader.h"
#include "AliGenSexaquarkReaction.h"
#include "AliLog.h"
#include "AliRun.h"

#define VERBOSE_MODE 0

ClassImp(AliGenSexaquarkReaction);  // integrate class for interactive use in ROOT

// Set default kinematic ranges for injected anti-sexaquark.
// Set default radius range for secondary vertex.
// Struck nucleon momentum assigned with simple gaussian.
AliGenSexaquarkReaction::AliGenSexaquarkReaction(int n_reactions, double sexaquark_mass, Sexaquark::EReactionChannel reaction_channel)
    : AliGenerator{},
      fRndm{0},
      fBlastwaveFunctions{},
      fReactionProducts_Mass{},     // set in `CacheReactionInfo`
      fReactionProducts_PdgCode{},  // set in `CacheReactionInfo`
      fStruckNucleon_Mass{},        // set in `CacheReactionInfo`
      fRadiusMin{},
      fRadiusMax{},
      fFermiMomentum{},
      fFermiMomentumSigma{},
      fSexaquarkMass{sexaquark_mass},
      fStruckNucleon_PdgCode{},  // set in `CacheReactionInfo`
      fN_ReactionProducts{},     // set in `CacheReactionInfo`
      fReactionsPerTrigger{n_reactions},
      fReactionChannel{reaction_channel} {

    SetDefaultNameTitle();
    SetDefaultRanges();

    InitFermiMomentumInfo();
    CacheReactionInfo();
    InitBlastwaveFunctions();

    SetNumberParticles(fReactionsPerTrigger * (2 + fN_ReactionProducts));
    AliConfig::Instance()->Add(this);
}

// Default constructor.
AliGenSexaquarkReaction::AliGenSexaquarkReaction()  //
    : AliGenSexaquarkReaction{1, 1.8, Sexaquark::EReactionChannel::kA} {}

// Set default name and title.
void AliGenSexaquarkReaction::SetDefaultNameTitle() {
    fName = "AliGenSexaquarkReaction";
    fTitle = "Generator of AntiSexaquark-Nucleon Reactions";
}

// Initialization of the generator.
void AliGenSexaquarkReaction::Init() {

    // check consistency of selected ranges
    if (TestBit(kPtRange) && TestBit(kMomentumRange)) {
        Fatal("Init", "You should not set the momentum range and the pt range!\n");
    }
    if ((!TestBit(kPtRange)) && (!TestBit(kMomentumRange))) {
        Fatal("Init", "You should set either the momentum or the pt range!\n");
    }
    if ((TestBit(kYRange) && TestBit(kThetaRange)) || (TestBit(kYRange) && TestBit(kEtaRange)) || (TestBit(kEtaRange) && TestBit(kThetaRange))) {
        Fatal("Init", "You should only set the range of one of these variables: y, eta or theta\n");
    }
    if ((!TestBit(kYRange)) && (!TestBit(kEtaRange)) && (!TestBit(kThetaRange))) {
        Fatal("Init", "You should set the range of one of these variables: y, eta or theta\n");
    }

    PrintParameters();
}

// Derive information from input `fReactionChannel`.
void AliGenSexaquarkReaction::CacheReactionInfo() {

    // -- struck nucleon
    fStruckNucleon_Mass = Sexaquark::Particle_Mass[Sexaquark::RChannel_NucleonPID[fReactionChannel]];
    fStruckNucleon_PdgCode = Sexaquark::Particle_PdgCode[Sexaquark::RChannel_NucleonPID[fReactionChannel]];

    // -- reaction products
    for (auto p : Sexaquark::RChannel_ProductsPID[fReactionChannel]) {
        fReactionProducts_PdgCode.push_back(Sexaquark::Particle_PdgCode[p]);
        fReactionProducts_Mass.push_back(Sexaquark::Particle_Mass[p]);
    }
    fN_ReactionProducts = static_cast<int>(fReactionProducts_PdgCode.size());
}

// Set default kinematic ranges for the anti-sexaquark, and radius range for the interaction vertex.
void AliGenSexaquarkReaction::SetDefaultRanges() {

    // -- anti-sexaquark
    SetPtRange(0., 10.);    // (GeV/c)
    SetPhiRange(0., 360.);  // (degrees) // NOTE: get converted into radians here)
    SetYRange(-0.8, 0.8);

    // -- interaction vertex
    SetRadiusRange(5., 180.);  // cm
}

// Initialise the Fermi momentum information.
void AliGenSexaquarkReaction::InitFermiMomentumInfo() {
    fFermiMomentum = 0.250;  // (GeV/c)
    fFermiMomentumSigma = fFermiMomentum / std::sqrt(5.);
}

void AliGenSexaquarkReaction::InitBlastwaveFunctions() {
    fBlastwaveFunctions.reserve(Sexaquark::kNCentralityClasses);
    for (size_t cc = 0; cc < Sexaquark::kNCentralityClasses; ++cc) {
        auto fcn = std::make_unique<TF1>(TString::Format("Fcn_BlastWave_%.2f_%s", fSexaquarkMass, Sexaquark::Centrality_Str[cc].data()),  //
                                         BlastWaveFcn, fPtMin, fPtMax, 4);
        fcn->SetParameter(0, fSexaquarkMass);
        fcn->SetParameter(1, Sexaquark::KineticFOTemp[cc]);
        fcn->SetParameter(2, Sexaquark::AvgBetaTemp[cc]);
        fcn->SetParameter(3, Sexaquark::RadialFlowPower[cc]);
        fBlastwaveFunctions.push_back(std::move(fcn));
    }
}

// Generate n_triggers.
void AliGenSexaquarkReaction::GenerateN(int n_triggers) {

    TGenPhaseSpace ThisPhaseSpace;
    int mult = n_triggers * fReactionsPerTrigger;

    for (int i = 0; i < mult; ++i) {

#if VERBOSE_MODE
        AliInfoF("\tSummary of Reaction #%i", i);
        AliInfo("");
#endif

        // 1) Prepare anti-sexaquark
        int ReactionID = Sexaquark::ReactionID_Offset + i;
        // -- assign random centrality, hoping the distribution averages out
        //    NOTE: width-weighted draw
        double fake_centrality = fRndm.Uniform(90.);
        int centrality_bin = 0;
        if (fake_centrality < 5.) {
            centrality_bin = 0;
        } else if (fake_centrality < 10.) {
            centrality_bin = 1;
        } else {
            centrality_bin = 2 + static_cast<int>((fake_centrality - 10.) / 10.);
        }
        double Pt_S = fBlastwaveFunctions[centrality_bin]->GetRandom(&fRndm);
        double Mt_S = std::sqrt(fSexaquarkMass * fSexaquarkMass + Pt_S * Pt_S);  // transverse mass (derived from inv. mass and Pt)
        double Y_S = fYMin + fRndm.Rndm() * (fYMax - fYMin);                     // rapidity (uniform distribution)
        double Phi_S = fPhiMin + fRndm.Rndm() * (fPhiMax - fPhiMin);             // azimuthal angle (uniform distribution) (in radians)
        double Px_S = Pt_S * std::cos(Phi_S);
        double Py_S = Pt_S * std::sin(Phi_S);
        double Pl_S = Mt_S * std::sinh(Y_S);  // longitudinal momentum = Pz (derived from trans. mass and Y)
        double E_S = Mt_S * std::cosh(Y_S);   // energy (derived from trans. mass and Y)
        TLorentzVector Sexaquark(Px_S, Py_S, Pl_S, E_S);

#if VERBOSE_MODE
        AliInfo("Generated AntiSexaquark:");
        AliInfoF(">> Px,Py,Pz = (%.4f, %.4f, %.4f)", Sexaquark.Px(), Sexaquark.Py(), Sexaquark.Pz());
        AliInfoF(">> M        = %.4f", Sexaquark.M());
        AliInfoF(">> Y        = %.4f", Y_S);
        AliInfoF(">> Theta    = %.4f", Sexaquark.Theta());
        AliInfoF(">> Phi      = %.4f", Phi_S);
#endif

        // 2) Prepare struck nucleon
        double Px_N = 0.;
        double Py_N = 0.;
        double Pz_N = 0.;
        GetFermiMomentum(Px_N, Py_N, Pz_N);
        double E_N = std::sqrt(Px_N * Px_N + Py_N * Py_N + Pz_N * Pz_N + fStruckNucleon_Mass * fStruckNucleon_Mass);
        TLorentzVector Nucleon(Px_N, Py_N, Pz_N, E_N);

        // 3) Generate reaction
        // determine momenta of reaction products
        TLorentzVector W = Sexaquark + Nucleon;
        ThisPhaseSpace.SetDecay(W, fN_ReactionProducts, &fReactionProducts_Mass[0]);
        double weight = (double)ThisPhaseSpace.Generate();

        // 4) Set secondary vertex
        double Radius_S = std::sqrt(fRadiusMin * fRadiusMin + fRndm.Rndm() * (fRadiusMax * fRadiusMax - fRadiusMin * fRadiusMin));
#if VERBOSE_MODE
        AliInfoF(">> Radius   = %.4f", Radius_S);
#endif
        double Theta_S = Sexaquark.Theta();          // polar angle (in radians)
        double Vx_S = Radius_S * std::cos(Phi_S);    // x, in cartesian coordinates
        double Vy_S = Radius_S * std::sin(Phi_S);    // y, in cartesian coordinates
        double Vz_S = Radius_S / std::tan(Theta_S);  // z, in both cartesian and cylindrical coordinates (unprotected on purpose)

        TVector3 SecondaryVertex(Vx_S, Vy_S, Vz_S);

        double Beta_S = Sexaquark.P() / Sexaquark.E();
        double SecondaryTime = (SecondaryVertex.Mag() * 1E-7) / (Beta_S * TMath::Ccgs());  // (in ns)

#if VERBOSE_MODE
        // -- print anti-sexaquark and nucleon information
        AliInfoF("%i,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e", ReactionID, Sexaquark.Px(), Sexaquark.Py(), Sexaquark.Pz(), Nucleon.Px(), Nucleon.Py(),
                 Nucleon.Pz());
#endif

#if VERBOSE_MODE
        AliInfo("   Vertex            R        Theta          Phi");
        AliInfoF("Secondary %12.7f %12.7f %12.7f", SecondaryVertex.Perp(), SecondaryVertex.Theta(), SecondaryVertex.Phi());
        AliInfo("");
        AliInfo("  I          PDG           Px           Py           Pz            E            M");
#endif

        // -- define track counter
        int nt = 0;

        // -- push anti-sexaquark
        PushTrack(0,                                                                // don't transport
                  -1,                                                               // no parent
                  Sexaquark::Particle_PdgCode[Sexaquark::EParticle::kAntiNeutron],  // NOTE: pretend it's an anti-neutron
                  Sexaquark.Px(), Sexaquark.Py(), Sexaquark.Pz(), Sexaquark.E(),    //
                  fVertex[0], fVertex[1], fVertex[2], fTime,                        //
                  0., 0., 0.,                                                       //
                  TMCProcess::kPUserDefined,                                        //
                  nt,                                                               //
                  weight,                                                           //
                  ReactionID);                                                      //

        // -- push struck nucleon
        PushTrack(0,                                                                     // don't transport
                  -1,                                                                    // no parent
                  fStruckNucleon_PdgCode,                                                //
                  Nucleon.Px(), Nucleon.Py(), Nucleon.Pz(), Nucleon.E(),                 //
                  SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fTime,  //
                  0., 0., 0.,                                                            //
                  TMCProcess::kPUserDefined,                                             //
                  nt,                                                                    //
                  weight,                                                                //
                  ReactionID);                                                           //

        // -- push reaction products
        for (int j = 0; j < fN_ReactionProducts; ++j) {
            PushTrack(1,                                 // transport
                      -1,                                // NOTE: pretend they are primaries
                      fReactionProducts_PdgCode[j],      // PDG code of the product particle
                      ThisPhaseSpace.GetDecay(j)->Px(),  // Px of the product particle (in GeV/c)
                      ThisPhaseSpace.GetDecay(j)->Py(),  // Py of the product particle (in GeV/c)
                      ThisPhaseSpace.GetDecay(j)->Pz(),  // Pz of the product particle (in GeV/c)
                      ThisPhaseSpace.GetDecay(j)->E(),   // energy of the product particle (in GeV)
                      SecondaryVertex.X(),               // vx: x-component of position of interaction vertex (in cm)
                      SecondaryVertex.Y(),               // vy: y-component of position of interaction vertex (in cm)
                      SecondaryVertex.Z(),               // vz: z-component of position of interaction vertex (in cm)
                      SecondaryTime,                     // vt: time of interaction vertex (in ns)
                      0.,                                // polx: x-component of polarisation vector, 0 by default
                      0.,                                // poly: y-component of polarisation vector, 0 by default
                      0.,                                // polz: z-component of polarisation vector, 0 by default
                      TMCProcess::kPHInhelastic,         // mech: production mechanism (Hadronic Inelastic Scattering)
                      nt,                                // output particle index
                      weight,                            // weight of the event
                      ReactionID);                       // generation status code
#if VERBOSE_MODE
            AliInfoF("%3i %12i %12.7f %12.7f %12.7f %12.7f %12.7f", nt, fReactionProducts_PdgCode[j], ThisPhaseSpace.GetDecay(j)->Px(),
                     ThisPhaseSpace.GetDecay(j)->Py(), ThisPhaseSpace.GetDecay(j)->Pz(), ThisPhaseSpace.GetDecay(j)->E(),
                     ThisPhaseSpace.GetDecay(j)->M());
#endif
        }

#if VERBOSE_MODE
        AliInfo("");
#endif
    }  // end of loop

    // Copied from AliGenBox //

    AliGenEventHeader* header = new AliGenEventHeader("SEXAQUARK");
    header->SetPrimaryVertex(fVertex);
    header->SetNProduced(n_triggers * fReactionsPerTrigger * (2 + fN_ReactionProducts));
    header->SetInteractionTime(fTime);

    // passes header either to the container or to gAlice (ownership transferred)
    if (fContainer) {
        header->SetName(fName);
        fContainer->AddHeader(header);
    } else {
        gAlice->SetGenEventHeader(header);
    }
    // no leak, because ownership is transferred in both branches
}

// Generate a single trigger.
void AliGenSexaquarkReaction::Generate() { GenerateN(1); }

// Print the parameters of the current instance.
void AliGenSexaquarkReaction::PrintParameters() {
    AliInfo("> Initializing AliGenSexaquarkReaction");
    AliInfo("");
    AliInfoF("  Number of Reactions = %i", fReactionsPerTrigger);
    AliInfoF("  Chosen Reaction Channel : %c", Sexaquark::RChannel_Char[fReactionChannel]);
    AliInfoF("  * Nucleon : %i", fStruckNucleon_PdgCode);
    TString ReactionProducts_str = "{";
    for (int p = 0; p < fN_ReactionProducts; ++p) {
        ReactionProducts_str += Form("%i", fReactionProducts_PdgCode[p]);
        if (p + 1 != fN_ReactionProducts) ReactionProducts_str += ", ";
    }
    ReactionProducts_str += "}";
    AliInfoF("  * Reaction Products : %s", ReactionProducts_str.Data());
    AliInfo("");
    AliInfo("  Injected AntiSexaquark:");
    AliInfoF("  * M = %.2f (GeV/c^2)", fSexaquarkMass);
    AliInfoF("  * %f < Pt < %f (GeV/c)", fPtMin, fPtMax);
    AliInfoF("  * %f < Phi < %f (rad)", fPhiMin, fPhiMax);
    AliInfoF("  * %f < Y < %f", fYMin, fYMax);
    AliInfo("");
    AliInfo("  Secondary Vertex:");
    AliInfoF("  * %f < Radius < %f (cm)", fRadiusMin, fRadiusMax);
    AliInfo("");
}

// Based on the current model for the Fermi momentum, get a random value.
// - uses: fFermiMomentumSigma
// - output: Px, Py, Pz
void AliGenSexaquarkReaction::GetFermiMomentum(double& Px, double& Py, double& Pz) {
    Px = fRndm.Gaus(0., fFermiMomentumSigma);
    Py = fRndm.Gaus(0., fFermiMomentumSigma);
    Pz = fRndm.Gaus(0., fFermiMomentumSigma);
}

double AliGenSexaquarkReaction::BlastWaveFcn(double* x, double* par) {

    double pt = x[0];

    double m0 = par[0];
    double t_kin = par[1];
    double avg_beta_temp = par[2];
    double np = par[3];

    double mt = std::sqrt(pt * pt + m0 * m0);
    double beta_s = (np + 2) * avg_beta_temp / 2.;

    static RMath::Integrator integrator;

    auto BlastWaveIntegrand = [&](double r) -> double {
        double rho = std::atanh(beta_s * std::pow(r, np));
        return r * TMath::BesselI0(pt * std::sinh(rho) / t_kin) * TMath::BesselK1(mt * std::cosh(rho) / t_kin);
    };

    RMath::Functor1D functor(BlastWaveIntegrand);
    integrator.SetFunction(functor);

    return pt * mt * integrator.Integral(0., 1.);
}
