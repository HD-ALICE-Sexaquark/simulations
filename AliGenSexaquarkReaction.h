#ifndef ALIGENSEXAQUARKREACTION_H
#define ALIGENSEXAQUARKREACTION_H

/* Copyright(c) 1998-2023, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <array>
#include <memory>
#include <string_view>
#include <vector>

#include <TF1.h>
#include <TRandomGen.h>

#include "AliGenerator.h"

namespace Sexaquark {

// Constants //

inline constexpr int ReactionID_Offset = 600;

// Particles DB //

enum EParticle {
    kPiMinus,
    kPiPlus,
    kPiZero,
    kNegKaon,
    kPosKaon,
    kAntiProton,
    kProton,
    kAntiNeutron,
    kNeutron,
    kKaonZeroShort,
    kAntiLambda,
    kLambda,
    kXiPlus,
    kNParticles
};
inline constexpr std::array<double, kNParticles> Particle_Mass{0.13957040, 0.13957040, 0.13497680, 0.49367700, 0.49367700,  //
                                                               0.93827210, 0.93827210, 0.93956540, 0.93956540, 0.49761100,
                                                               1.1156830,  1.1156830,  1.3217100};
inline constexpr std::array<int, kNParticles> Particle_PdgCode{-211,  211,  111,   -321, 321,  //
                                                               -2212, 2212, -2112, 2112, 310,  //
                                                               -3122, 3122, -3312};
inline constexpr std::array<std::string_view, kNParticles> Particle_Name{"PiMinus",    "PiPlus", "PiZero",      "NegKaon", "PosKaon",
                                                                         "AntiProton", "Proton", "AntiNeutron", "Neutron", "KaonZeroShort",
                                                                         "AntiLambda", "Lambda", "XiPlus"};

// Reaction Channels //

enum EReactionChannel {
    kA,  // AntiS + Neutron -> AntiLambda, K0
    kB,  // AntiS + Neutron -> AntiLambda, K0, Pi-, Pi+
    kC,  // AntiS + Neutron -> AntiProton, K0, K0, Pi+
    kD,  // AntiS + Proton -> AntiLambda, K+
    kE,  // AntiS + Proton -> AntiLambda, K+, Pi-, Pi+
    kF,  // AntiS + Proton -> AntiProton, K+, K0, Pi+
    kG,  // AntiS + Neutron -> Xi+, Pi-
    kH,  // AntiS + Proton -> AntiProton, K+, K+, Pi0
    kNChannels
};
inline constexpr std::array<char, kNChannels> RChannel_Char{'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
inline constexpr std::array<EParticle, kNChannels> RChannel_NucleonPID{EParticle::kNeutron, EParticle::kNeutron, EParticle::kNeutron,
                                                                       EParticle::kProton,  EParticle::kProton,  EParticle::kProton,
                                                                       EParticle::kNeutron, EParticle::kProton};
inline const std::array<std::vector<EParticle>, kNChannels> RChannel_ProductsPID{
    {{EParticle::kAntiLambda, EParticle::kKaonZeroShort},                                                 //
     {EParticle::kAntiLambda, EParticle::kKaonZeroShort, EParticle::kPiMinus, EParticle::kPiPlus},        //
     {EParticle::kAntiProton, EParticle::kKaonZeroShort, EParticle::kKaonZeroShort, EParticle::kPiPlus},  //
     {EParticle::kAntiLambda, EParticle::kPosKaon},                                                       //
     {EParticle::kAntiLambda, EParticle::kPosKaon, EParticle::kPiMinus, EParticle::kPiPlus},              //
     {EParticle::kAntiProton, EParticle::kPosKaon, EParticle::kKaonZeroShort, EParticle::kPiPlus},        //
     {EParticle::kXiPlus, EParticle::kPiMinus},                                                           //
     {EParticle::kAntiProton, EParticle::kPosKaon, EParticle::kPosKaon, EParticle::kPiZero}}};

constexpr EReactionChannel ReactionChannelFromChar(char c) {
    for (size_t i = 0; i < kNChannels; ++i) {
        if (RChannel_Char[i] == c) return static_cast<EReactionChannel>(i);
    }
    return EReactionChannel::kA;  // fallback
}

// Blast-wave Function //

// -- Parameters taken from:
//    "Production of charged pions, kaons, and (anti-)protons in Pb-Pb and inelastic pp collisions at √sNN = 5.02 TeV"
//    ALICE Collab. PRC 101, 044907 (2020)
inline constexpr size_t kNCentralityClasses = 10;
inline constexpr std::array<std::string_view, kNCentralityClasses> Centrality_Str{"0-5",   "5-10",  "10-20", "20-30", "30-40",  //
                                                                                  "40-50", "50-60", "60-70", "70-80", "80-90"};

// kinetic freeze-out temperature (in GeV)
inline constexpr std::array<float, kNCentralityClasses> KineticFOTemp{0.090, 0.091, 0.094, 0.097, 0.101, 0.108, 0.115, 0.129, 0.147, 0.161};

// transverse velocity profile
inline constexpr std::array<float, kNCentralityClasses> AvgBetaTemp{0.663, 0.660, 0.655, 0.643, 0.622, 0.595, 0.557, 0.506, 0.435, 0.355};

// radial flow power
inline constexpr std::array<float, kNCentralityClasses> RadialFlowPower{0.735, 0.736, 0.739, 0.771, 0.828, 0.908, 1.052, 1.262, 1.678, 2.423};

}  // namespace Sexaquark

// Injector of the products of AntiSexaquark-Nucleon interactions
// - anti-sexaquark kinematics under a flat distribution
// - radius of secondary vertex chosen from a flat distribution
// - struck nucleon momentum assigned with simple gaussian
// -- based on /EVGEN/AliGenBox[.cxx/.h]
// -- author: aborquez@cern.ch
class AliGenSexaquarkReaction : public AliGenerator {

   public:
    AliGenSexaquarkReaction(int n_reactions, double sexaquark_mass, Sexaquark::EReactionChannel reaction_channel);
    AliGenSexaquarkReaction();  // NOTE: needed by the framework
    virtual ~AliGenSexaquarkReaction() = default;
    virtual void GenerateN(int n_triggers);
    virtual void Generate();
    virtual void Init();

    // Set minimum and maximum radius for uniform distribution to locate interaction (or secondary) vertex.
    void SetRadiusRange(double radiusMin, double radiusMax) {
        fRadiusMin = radiusMin;
        fRadiusMax = radiusMax;
    }

   private:
    void SetDefaultRanges();
    void SetDefaultNameTitle();
    void CacheReactionInfo();
    void PrintParameters();
    void InitFermiMomentumInfo();
    void GetFermiMomentum(double &px, double &py, double &pz);
    void InitBlastwaveFunctions();
    static double BlastWaveFcn(double *x, double *par);

    TRandomMT64 fRndm;                                      //
    std::vector<std::unique_ptr<TF1>> fBlastwaveFunctions;  //
    std::vector<double> fReactionProducts_Mass;             // inv. masses of reaction products
    std::vector<int> fReactionProducts_PdgCode;             // PDG codes of reaction products
    double fStruckNucleon_Mass;                             // invariant mass of struck nucleon
    double fRadiusMin;                                      // minimum radius
    double fRadiusMax;                                      // maximum radius
    double fFermiMomentum;                                  // central value of Fermi Momentum
    double fFermiMomentumSigma;                             // sigma of Fermi Momentum distribution
    double fSexaquarkMass;                                  // (input) sexaquark mass
    int fStruckNucleon_PdgCode;                             // PDG code of struck nucleon
    int fN_ReactionProducts;                                // number of product particles per reaction
    int fReactionsPerTrigger;                               // (input) number of reactions per trigger
    Sexaquark::EReactionChannel fReactionChannel;           // (input) reaction channel

    ClassDef(AliGenSexaquarkReaction, 3);  // integrate class for interactive use in ROOT
};

#endif
