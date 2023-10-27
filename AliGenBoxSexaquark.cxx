#include "TDatabasePDG.h"
#include "TF1.h"
#include "TMath.h"
#include "TPDGCode.h"

#include "AliConst.h"
#include "AliGenBoxSexaquark.h"
#include "AliGenEventHeader.h"
#include "AliLog.h"
#include "AliPDG.h"
#include "AliRun.h"

#define VERBOSE_MODE 1
#define PDG_SEXAQUARK -900000020

ClassImp(AliGenBoxSexaquark);

//_____________________________________________________________________________
AliGenBoxSexaquark::AliGenBoxSexaquark() : AliGenerator(), fIpart(PDG_SEXAQUARK), fEtaMin(0), fEtaMax(0), fPtDistribution(0) {
    //
    // Default constructor
    //
}

//_____________________________________________________________________________
AliGenBoxSexaquark::AliGenBoxSexaquark(Int_t npart)
    : AliGenerator(npart),
      fIpart(PDG_SEXAQUARK),
      fEtaMin(0),
      fEtaMax(0),
      fRandomOffset(kFALSE),
      fRMinOffset(0),
      fRMaxOffset(0),
      fZMinOffset(0),
      fZMaxOffset(0),
      fPtDistribution(0) {
    //
    // Standard constructor
    //
    fName = "AliGenBoxSexaquark";
    fTitle = "Generator of AntiSexaquarks";
}

//_____________________________________________________________________________
AliGenBoxSexaquark::~AliGenBoxSexaquark() {
    //
    // Destructor
    //
    delete fPtDistribution;
}

//_____________________________________________________________________________
void AliGenBoxSexaquark::Generate() {
    //
    // Generate one trigger (fNpart particles)
    //
    GenerateN(1);
}

//_____________________________________________________________________________
void AliGenBoxSexaquark::GenerateN(Int_t ntimes) {
    //
    // Generate ntimes triggers
    //   total ntimes*fNpart particles
    //

    Float_t polar[3] = {0, 0, 0};

    Float_t origin[3];
    Float_t time;
    Float_t p[3];
    Int_t i, j, nt;
    Float_t pmom, theta, phi, pt;
    Float_t y, mt;
    Float_t energy;

    // Prepare random numbers:
    // if the seed is zero, then the seed is set to a random value
    // which in case of TRandom depends on the lowest 4 bytes of UUID
    GetRandom()->SetSeed(0);
    Float_t random[6];

    std::vector<Float_t> myIntegral;
    std::vector<Float_t> myAlpha;
    std::vector<Float_t> myBeta;
    std::vector<Float_t> myGamma;

    for (j = 0; j < 3; j++) origin[j] = fOrigin[j];
    time = fTimeOrigin;
    if (fVertexSmear == kPerEvent && !fRandomOffset) {
        Vertex();
        for (j = 0; j < 3; j++) {
            origin[j] = fVertex[j];
        }
        time = fTime;
    }

    Float_t m = TDatabasePDG::Instance()->GetParticle(fIpart)->Mass();

    Int_t mult = fNpart * ntimes;
    for (i = 0; i < mult; i++) {
        Rndm(random, 3);

        if (TestBit(kYRange)) {
            y = fYMin + random[0] * (fYMax - fYMin);

            if (fPtDistribution) {
                pt = GetRandomFromFunction(myIntegral, myAlpha, myBeta, myGamma, random[1]);
                mt = TMath::Sqrt(pt * pt + m * m);
            } else if (TestBit(kMomentumRange)) {
                pmom = fPMin + random[1] * (fPMax - fPMin);
                mt = TMath::Sqrt(pmom * pmom + m * m) / TMath::CosH(y);
                pt = TMath::Sqrt(mt * mt - m * m);
            } else {
                pt = fPtMin + random[1] * (fPtMax - fPtMin);
                mt = TMath::Sqrt(pt * pt + m * m);
            }

            phi = fPhiMin + random[2] * (fPhiMax - fPhiMin);
            p[0] = pt * TMath::Cos(phi);
            p[1] = pt * TMath::Sin(phi);
            p[2] = mt * TMath::SinH(y);
            pmom = TMath::Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        } else {
            if (TestBit(kThetaRange)) {
                theta = fThetaMin + random[0] * (fThetaMax - fThetaMin);
            } else {
                Float_t eta = fEtaMin + random[0] * (fEtaMax - fEtaMin);
                theta = 2. * TMath::ATan(TMath::Exp(-eta));
            }

            if (fPtDistribution) {
                pt = GetRandomFromFunction(myIntegral, myAlpha, myBeta, myGamma, random[1]);
                pmom = pt / TMath::Sin(theta);
            } else if (TestBit(kMomentumRange)) {
                pmom = fPMin + random[1] * (fPMax - fPMin);
                pt = pmom * TMath::Sin(theta);
            } else {
                pt = fPtMin + random[1] * (fPtMax - fPtMin);
                pmom = pt / TMath::Sin(theta);
            }

            phi = fPhiMin + random[2] * (fPhiMax - fPhiMin);
            p[0] = pt * TMath::Cos(phi);
            p[1] = pt * TMath::Sin(phi);
            p[2] = pmom * TMath::Cos(theta);
        }

        if (fVertexSmear == kPerTrack || fRandomOffset) {
            if (fRandomOffset) {
                Rndm(random, 2);
                float roffs = random[0] * (fRMaxOffset - fRMinOffset) + fRMinOffset;
                origin[0] = fOrigin[0] + roffs * TMath::Cos(phi);
                origin[1] = fOrigin[1] + roffs * TMath::Sin(phi);
                origin[2] = fOrigin[2] + random[1] * (fZMaxOffset - fZMinOffset) + fZMinOffset;
            } else {
                Rndm(random, 6);
                for (j = 0; j < 3; j++) {
                    origin[j] = fOrigin[j] +
                                fOsigma[j] * TMath::Cos(2 * random[2 * j] * TMath::Pi()) * TMath::Sqrt(-2 * TMath::Log(random[2 * j + 1]));
                }
            }

            Rndm(random, 2);
            time = fTimeOrigin +
                   fOsigma[2] / TMath::Ccgs() * TMath::Cos(2 * random[0] * TMath::Pi()) * TMath::Sqrt(-2 * TMath::Log(random[1]));
        }

        PushTrack(fTrackIt, -1, fIpart, p, origin, polar, time, kPNoProcess, nt, 1., 1);  // previously it was kPPrimary

        energy = TMath::Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + m * m);

#if VERBOSE_MODE
        AliInfo("          Id           Px           Py           Pz           Pt            E            M");
        AliInfoF("%12i %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f", fIpart, p[0], p[1], p[2], pt, energy, m);
#endif
    }

    AliGenEventHeader *header = new AliGenEventHeader("BOX");
    header->SetPrimaryVertex(fVertex);
    header->SetNProduced(mult);
    header->SetInteractionTime(fTime);

    // Passes header either to the container or to gAlice
    if (fContainer) {
        header->SetName(fName);
        fContainer->AddHeader(header);
    } else {
        gAlice->SetGenEventHeader(header);
    }
}

//_____________________________________________________________________________
void AliGenBoxSexaquark::Init() {
    // Initialisation, check consistency of selected ranges
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
    AliInfo("");
    AliInfo("> Generating Anti Sexaquarks with the following kinematics:");
    AliInfo("");
    AliInfoF("  PDG Code        = %i", fIpart);
    if (fPtDistribution) {
        AliInfoF("  Pt Distribution = %s", fPtDistribution->GetExpFormula().Data());
    }
    AliInfoF("  %f < Pt < %f", fPtMin, fPtMax);
    AliInfoF("  %f < Phi < %f", fPhiMin, fPhiMax);
    AliInfoF("  %f < Y < %f", fYMin, fYMax);
    AliInfo("");
#endif
}

//_____________________________________________________________________________
void AliGenBoxSexaquark::SetRandomOffset(float rmin, float rmax, float zmin, float zmax) {
    if (rmin < 0 || rmin > rmax) {
        AliFatalF("Wrong random R offset request: RMin: %.3f RMax: %.3f", rmin, rmax);
    }
    if (zmin > zmax) {
        AliFatalF("Wrong random Z offset request: ZMin: %.3f ZMax: %.3f", zmin, zmax);
    }
    fRMaxOffset = rmax;
    fRMinOffset = rmin;
    fZMaxOffset = zmax;
    fZMinOffset = zmin;
    fRandomOffset = kTRUE;
    AliInfoF("Origin of every track will be randomly shifted by %.2f<R<%.2f and %.2f<Z<%.2f", fRMinOffset, fRMaxOffset, fZMinOffset,
             fZMaxOffset);
}

//_____________________________________________________________________________
void AliGenBoxSexaquark::SetPtDistribution(TF1 &PtDistribution) {
    TF1 *tmp = new TF1(PtDistribution);
    fPtDistribution = tmp;
    if (fPtDistribution) {
        SetPtRange(fPtDistribution->GetXmin(), fPtDistribution->GetXmax());
    }
}

//_____________________________________________________________________________
Bool_t AliGenBoxSexaquark::ComputeCdfTable(std::vector<Float_t> &Integral_vec, std::vector<Float_t> &Alpha_vec,
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
Float_t AliGenBoxSexaquark::GetRandomFromFunction(std::vector<Float_t> &Integral_vec, std::vector<Float_t> &Alpha_vec,
                                                   std::vector<Float_t> &Beta_vec, std::vector<Float_t> &Gamma_vec, Float_t r) {

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
    // Float_t r = rng->Rndm();
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
