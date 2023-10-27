#ifndef ALIGENBOXSEXAQUARK_H
#define ALIGENBOXSEXAQUARK_H

// Generator of sexaquarks in a preset kinematic range (flat distribution)

#include "AliGenerator.h"

class AliGenBoxSexaquark : public AliGenerator {
   public:
    AliGenBoxSexaquark();
    AliGenBoxSexaquark(Int_t npart);
    virtual ~AliGenBoxSexaquark();
    virtual void GenerateN(Int_t ntimes);
    virtual void Generate();
    virtual void Init();
    virtual void SetEtaRange(Float_t etamin, Float_t etamax) {
        SetBit(kEtaRange);
        fEtaMin = etamin;
        fEtaMax = etamax;
    }
    virtual void SetPart(Int_t part) { fIpart = part; }
    virtual void SetParticleType(Int_t part) { SetPart(part); }
    virtual void SetSeed(UInt_t /*seed*/) { ; }
    virtual void SetPtDistribution(TF1 &PtDistribution);

    Bool_t ComputeCdfTable(std::vector<Float_t> &Integral_vec, std::vector<Float_t> &Alpha_vec, std::vector<Float_t> &Beta_vec,
                           std::vector<Float_t> &Gamma_vec);
    Float_t GetRandomFromFunction(std::vector<Float_t> &Integral_vec, std::vector<Float_t> &Alpha_vec, std::vector<Float_t> &Beta_vec,
                                   std::vector<Float_t> &Gamma_vec, Float_t r);

    void SetRandomOffset(float rmin, float rmax, float zmin, float zmax);

   protected:
    Int_t fIpart;          // Particle type
    Float_t fEtaMin;       // Minimum eta
    Float_t fEtaMax;       // Maximum eta
    Bool_t fRandomOffset;  // add random offset to origin per track
    Float_t fRMinOffset;   // min R of the offset
    Float_t fRMaxOffset;   // max R of the offset
    Float_t fZMinOffset;   // max Z of the offset
    Float_t fZMaxOffset;   // max Z of the offset

    TF1 *fPtDistribution;  // model to choose random Pt from

    ClassDef(AliGenBoxSexaquark, 3)
};

#endif
