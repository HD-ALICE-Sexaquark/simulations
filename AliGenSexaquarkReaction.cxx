#include <fstream>
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

ClassImp(AliGenSexaquarkReaction);  // integrate class for interactive use in ROOT

/*
 Default constructor:
 - Set all kinematic ranges
 - As no specific nucleon or decay channel was provided, all reactions will be produced randomly
 */
AliGenSexaquarkReaction::AliGenSexaquarkReaction()
    : AliGenerator(), fStruckNucleonPDG(0), fSexaquarkMass(0.), fReactionProductsPDG(0), fRadiusMin(0.), fRadiusMax(0.) {
    fName = "AliGenSexaquarkReaction";
    fTitle = "Generator of AntiSexaquark-Nucleon Reactions";
    InitReactionsInfo();
    InitFermiMomentumInfo();
    SetDefaultRanges();
    fCurrentReactionID = 600;
}

/*
 Default constructor:
 - Set all kinematic ranges
 - As no specific nucleon or decay channel was provided, all reactions will be produced randomly
 */
AliGenSexaquarkReaction::AliGenSexaquarkReaction(Int_t NReactions, Float_t SexaquarkMass, Char_t ReactionChannel)
    : AliGenerator(NReactions),
      fStruckNucleonPDG(0),
      fSexaquarkMass(SexaquarkMass),
      fReactionChannel(ReactionChannel),
      fReactionProductsPDG(0),
      fRadiusMin(0.),
      fRadiusMax(0.) {
    fName = "AliGenSexaquarkReaction";
    fTitle = "Generator of AntiSexaquark-Nucleon Reactions";
    InitReactionsInfo();
    InitFermiMomentumInfo();
    SetDefaultRanges();
    fCurrentReactionID = 600;
}

/*
 Default destructor
*/
AliGenSexaquarkReaction::~AliGenSexaquarkReaction() { delete fFermiMomentumModel; }

/*
 Initialisation of the generator
*/
void AliGenSexaquarkReaction::Init() {

    // check consistency of selected ranges
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

    // assign current nucleon and reaction channel
    fStruckNucleonPDG = fReactionChannelsMap[fReactionChannel][0][0];
    fReactionProductsPDG = fReactionChannelsMap[fReactionChannel][1];

#if VERBOSE_MODE
    PrintParameters();
#endif
}

/*
 Prepare the allowed decay channels to the set of vectors
 */
void AliGenSexaquarkReaction::SetDefaultRanges() {
    SetRadiusRange(5., 180.);  // cm
    SetPtRange(0., 5.);        // GeV/c
    SetPhiRange(0., 360.);     // degrees
    SetYRange(-0.8, 0.8);      // no units
}

/*
 Prepare the allowed decay channels to the set of vectors
 */
void AliGenSexaquarkReaction::InitReactionsInfo() {
    fReactionChannelsMap['A'] = {{2112}, {-3122, 310}};             // AntiS + Neutron -> AntiLambda, K0
    fReactionChannelsMap['B'] = {{2112}, {-3122, 310, -211, 211}};  // AntiS + Neutron -> AntiLambda, K0, Pi-, Pi+
    fReactionChannelsMap['C'] = {{2112}, {-2212, 310, 310, 211}};   // AntiS + Neutron -> AntiProton, K0, K0, Pi+
    fReactionChannelsMap['D'] = {{2212}, {-3122, 321}};             // AntiS + Proton -> AntiLambda, K+
    fReactionChannelsMap['E'] = {{2212}, {-3122, 321, -211, 211}};  // AntiS + Proton -> AntiLambda, K+, Pi-, Pi+
    fReactionChannelsMap['F'] = {{2212}, {-2212, 321, 310, 211}};   // AntiS + Proton -> AntiProton, K+, K0, Pi+
    fReactionChannelsMap['G'] = {{2112}, {-3312, -211}};            // AntiS + Neutron -> Xi+, Pi-
    fReactionChannelsMap['H'] = {{2212}, {-2212, 321, 321, 111}};   // AntiS + Proton -> AntiProton, K+, K+, Pi0
}

/*
 Initialise the Fermi momentum information
 */
void AliGenSexaquarkReaction::InitFermiMomentumInfo() {
    // central value and error obtained from:
    // Povh, Rith, Scholz, Zetsche, Rodejohann. "Particles and Nuclei: An Introduction to the Physical Concepts". 7th edition (Springer,
    // 2015)
    fFermiMomentum = 0.250;       // GeV
    fFermiMomentumError = 0.005;  // GeV
    fFermiMomentumModel = new TF1("Fermi Momentum Model", "[0]*exp(-0.5*((x-[1])/[2])**2) ", fFermiMomentum - 5 * fFermiMomentumError,
                                  fFermiMomentum + 5 * fFermiMomentumError);
    fFermiMomentumModel->SetParameter(0, 1.);
    fFermiMomentumModel->SetParameter(1, fFermiMomentum);
    fFermiMomentumModel->SetParameter(2, fFermiMomentumError);
}

/*
 Generate one anti-sexaquark into nucleon interaction, and push event to MC Stack
 */
void AliGenSexaquarkReaction::Generate() {
    // (pending) add an identifier for every reaction, and write the injected sexaquark and struck nucleon's momentum into a csv file

    // prepare random numbers
    GetRandom()->SetSeed(0);
    Float_t fRandom[4];
    Rndm(fRandom, 4);

    // 1) Prepare AntiSexaquark
    Float_t Pt_S = fPtMin + fRandom[0] * (fPtMax - fPtMin);
    Float_t Mt_S = TMath::Sqrt(fSexaquarkMass * fSexaquarkMass + Pt_S * Pt_S);  // transverse mass (derived from inv. mass and Pt)
    Float_t Y_S = fYMin + fRandom[1] * (fYMax - fYMin);                         // rapidity (uniform distribution)
    Float_t Phi_S = fPhiMin + fRandom[2] * (fPhiMax - fPhiMin);                 // azimuthal angle (uniform distribution) (in radians)
    Float_t Px_S = Pt_S * TMath::Cos(Phi_S);                                    // Px
    Float_t Py_S = Pt_S * TMath::Sin(Phi_S);                                    // Py
    Float_t Pl_S = Mt_S * TMath::SinH(Y_S);  // longitudinal momentum = Pz (derived from trans. mass and Y)
    Float_t E_S = Mt_S * TMath::CosH(Y_S);   // energy (derived from trans. mass and Y)
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
    Float_t Px_N = 0.;
    Float_t Py_N = 0.;
    Float_t Pz_N = 0.;
    GetFermiMomentum(Px_N, Py_N, Pz_N);
    Float_t M_N = TDatabasePDG::Instance()->GetParticle(fStruckNucleonPDG)->Mass();
    Float_t E_N = TMath::Sqrt(TMath::Power(Px_N, 2) + TMath::Power(Py_N, 2) + TMath::Power(Pz_N, 2) + TMath::Power(M_N, 2));
    TLorentzVector Nucleon;
    Nucleon.SetPxPyPzE(Px_N, Py_N, Pz_N, E_N);

    // 3) Prepare reaction
    // create array of masses of reaction products (must be a double)
    Double_t ReactionProductsMasses[(Int_t)fReactionProductsPDG.size()];
    for (Int_t i = 0; i < (Int_t)fReactionProductsPDG.size(); i++) {
        ReactionProductsMasses[i] = TDatabasePDG::Instance()->GetParticle(fReactionProductsPDG[i])->Mass();
    }
    // ROOT class to generate n-body event with constant cross-section
    TGenPhaseSpace ThisPhaseSpace;
    // determine momenta of reaction products
    TLorentzVector W = Sexaquark + Nucleon;
    ThisPhaseSpace.SetDecay(W, (Int_t)fReactionProductsPDG.size(), &ReactionProductsMasses[0]);
    Float_t weight = (Float_t)ThisPhaseSpace.Generate();

    // 5) Set primary vertex (taken from AliGenCocktail)
    TVector3 PrimaryVertex(fVertex[0], fVertex[1], fVertex[2]);  // (in cm)
    Float_t PrimaryTime = fTime;                                 // (in ns)

    // 6) Set secondary vertex
    // let's start the transformation, we already have
    // - radius_S: 2D Radius -- radius in cylindrical coordinates
    // - phi_S: phi -- the same in both spherical and cylindrical coordinates
    Float_t Radius_S = fRadiusMin + fRandom[3] * (fRadiusMax - fRadiusMin);  // radius (uniform distribution) (in cm)
#if VERBOSE_MODE
    AliInfoF(">> Radius   = %.4f", Radius_S);
#endif
    Float_t Theta_S = Sexaquark.Theta();            // polar angle (in radians)
    Float_t Vz_S = Radius_S / TMath::Tan(Theta_S);  // z, in both cartesian and cylindrical coordinates
    Float_t Radius3D_S = TMath::Sqrt(TMath::Power(Radius_S, 2) + TMath::Power(Vz_S, 2));  // 3D Radius -- radius in spherical coordinates
    Float_t Vx_S = Radius3D_S * TMath::Sin(Theta_S) * TMath::Cos(Phi_S);                  // x, in cartesian coordinates
    Float_t Vy_S = Radius3D_S * TMath::Sin(Theta_S) * TMath::Sin(Phi_S);                  // y, in cartesian coordinates
    TVector3 SecondaryVertex(Vx_S, Vy_S, Vz_S);

    // add displacement from PV
    SecondaryVertex = SecondaryVertex + PrimaryVertex;

    Float_t Beta_S = Sexaquark.P() / Sexaquark.E();
    Float_t SecondaryTime = (SecondaryVertex.Mag() * TMath::Power(10, -7)) / (Beta_S * TMath::Ccgs());  // (in ns)

    /* Output anti-sexaquark and nucleon information into an external CSV file */

    std::ofstream fOutputFile;
    TString fOutputFilename = "this_event.csv";
    fOutputFile.open(fOutputFilename, std::ofstream::app);

    TString auxStr = Form("%i,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",      //
                          fCurrentReactionID,                                       //
                          PrimaryVertex.X(), PrimaryVertex.Y(), PrimaryVertex.Z(),  //
                          Sexaquark.Px(), Sexaquark.Py(), Sexaquark.Pz(),           //
                          Nucleon.Px(), Nucleon.Py(), Nucleon.Pz());
    fOutputFile << auxStr;
    fOutputFile.close();

    AliInfoF(">> %s has been generated.", fOutputFilename.Data());

#if VERBOSE_MODE
    // verbose
    AliInfo("   Vertex            R        Theta          Phi");
    AliInfoF("  Primary %12.7f %12.7f %12.7f", PrimaryVertex.Perp(), PrimaryVertex.Theta(), PrimaryVertex.Phi());
    AliInfoF("Secondary %12.7f %12.7f %12.7f", SecondaryVertex.Perp(), SecondaryVertex.Theta(), SecondaryVertex.Phi());
    AliInfo("");
    AliInfo("  I          PDG           Px           Py           Pz            E            M");
#endif

    // get 4-momentum vectors of the reaction products and put them on the Geant stack
    Int_t TrackCounter = 0;
    for (Int_t i = 0; i < (Int_t)fReactionProductsPDG.size(); i++) {

        // add particle to the stack
        PushTrack(fTrackIt,                          // done: track final state particles
                  -1,                                // parent: the idea is that is the sexaquark
                                                     // (for the moment, it's -1 = primary)
                  fReactionProductsPDG[i],           // pdg of the product particle
                  ThisPhaseSpace.GetDecay(i)->Px(),  // px: Px of the product particle (in GeV/c)
                  ThisPhaseSpace.GetDecay(i)->Py(),  // py: Py of the product particle (in GeV/c)
                  ThisPhaseSpace.GetDecay(i)->Pz(),  // pz: Pz of the product particle (in GeV/c)
                  ThisPhaseSpace.GetDecay(i)->E(),   // E: energy of the product particle (in GeV)
                  SecondaryVertex.X(),               // vx: x-component of position of interaction vertex (in cm)
                  SecondaryVertex.Y(),               // vy: y-component of position of interaction vertex (in cm)
                  SecondaryVertex.Z(),               // vz: z-component of position of interaction vertex (in cm)
                  SecondaryTime - PrimaryTime,       // vt: time of interaction vertex (in seconds)
                  0.,                                // polx: x-component of polarisation vector, 0 by default
                  0.,                                // poly: y-component of polarisation vector, 0 by default
                  0.,                                // polz: z-component of polarisation vector, 0 by default
                  kPHInhelastic,                     // mech: production mechanism (Hadronic Inelastic Scattering)
                  TrackCounter,                      // counter of tracks (incremented internally by AliStack)
                  weight,                            // weight of the event
                  fCurrentReactionID);               // is: generation status code
                                                     // (0: unstable, 1: stable, 6: from sexaquark)
        //     1 * (TMath::Abs(fReactionProductsPDG[i]) == 211 || TMath::Abs(fReactionProductsPDG[i]) == 2212 ||
        //            TMath::Abs(fReactionProductsPDG[i]) == 321)
#if VERBOSE_MODE
        // verbose
        AliInfoF("%3i %12i %12.7f %12.7f %12.7f %12.7f %12.7f", TrackCounter, fReactionProductsPDG[i], ThisPhaseSpace.GetDecay(i)->Px(),
                 ThisPhaseSpace.GetDecay(i)->Py(), ThisPhaseSpace.GetDecay(i)->Pz(), ThisPhaseSpace.GetDecay(i)->E(),
                 ThisPhaseSpace.GetDecay(i)->M());
#endif
    }

    // unknown things with AliGenEventHeader and Containers and gAlice ???
    TArrayF auxSecondaryVertex(3);
    auxSecondaryVertex[0] = SecondaryVertex.X();
    auxSecondaryVertex[1] = SecondaryVertex.Y();
    auxSecondaryVertex[2] = SecondaryVertex.Z();

    AliGenEventHeader* header = new AliGenEventHeader("BOX");
    header->SetPrimaryVertex(auxSecondaryVertex);              // the real vertex of the interaction!
    header->SetNProduced((Int_t)fReactionProductsPDG.size());  // amount of particles produced?
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

/*
 Generate N sexaquark-nucleon reactions
 */
void AliGenSexaquarkReaction::GenerateN(Int_t N) {
    for (Int_t i = 0; i < N; i++) {
#if VERBOSE_MODE
        AliInfo("");
        AliInfo("\tEvent listing (summary)");
        AliInfo("");
#endif
        fCurrentReactionID = fCurrentReactionID + i;
        Generate();
    }
}

/*
 Print the parameters of the current instance
*/
void AliGenSexaquarkReaction::PrintParameters() {
    AliInfo("> Starting new AntiSexaquark-Nucleon Reaction with the following parameters:");
    AliInfo("");
    AliInfoF("  Reaction Channel : %c", fReactionChannel);
    AliInfoF("  >> Nucleon : %i", fStruckNucleonPDG);
    TString ReactionProducts_str = "{";
    for (Int_t p = 0; p < (Int_t)fReactionProductsPDG.size(); p++) {
        ReactionProducts_str += Form("%i", fReactionProductsPDG[p]);
        if (p + 1 != (Int_t)fReactionProductsPDG.size()) ReactionProducts_str += ", ";
    }
    ReactionProducts_str += "}";
    AliInfoF("  >> Reaction Products : %s", ReactionProducts_str.Data());
    AliInfo("");
    AliInfo("Injected AntiSexaquark Kinematics:");
    AliInfoF("  M = %f", fSexaquarkMass);
    AliInfoF("  %f < Pt < %f", fPtMin, fPtMax);
    AliInfoF("  %f < Phi < %f", fPhiMin, fPhiMax);
    AliInfoF("  %f < Y < %f", fYMin, fYMax);
    AliInfoF("  %f < Radius < %f", fRadiusMin, fRadiusMax);
    AliInfo("");
}

/*
 Based on the current model for the Fermi momentum, get a random value
 - uses: fFermiMomentum, fFermiMomentumError, fFermiMomentumModel
 - output: Px, Py, Pz
 */
void AliGenSexaquarkReaction::GetFermiMomentum(Float_t& Px, Float_t& Py, Float_t& Pz) {

    Float_t P_N = fFermiMomentumModel->GetRandom();
    // given the total momentum, set phi,theta as uniform and uncorrelated variables
    Float_t random_N[2];
    Rndm(random_N, 2);
    Float_t Phi_N = 2 * TMath::Pi() * random_N[0];  // azimuthal angle (uniform distribution) (in radians)
    Float_t Theta_N = TMath::Pi() * random_N[1];    // polar angle (uniform distribution) (in radians)
    // then, assign px,py,pz
    Px = P_N * TMath::Cos(Phi_N) * TMath::Sin(Theta_N);
    Py = P_N * TMath::Sin(Phi_N) * TMath::Sin(Theta_N);
    Pz = P_N * TMath::Cos(Theta_N);
}
