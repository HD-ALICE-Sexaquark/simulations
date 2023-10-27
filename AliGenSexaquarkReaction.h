#ifndef ALIGENSEXAQUARKREACTION_H
#define ALIGENSEXAQUARKREACTION_H

#include "AliGenerator.h"

class AliGenSexaquarkReaction : public AliGenerator {
   public:
    AliGenSexaquarkReaction();
    virtual ~AliGenSexaquarkReaction();
    virtual void GenerateN(Int_t N);
    virtual void Generate();
    virtual void Init();
    void SetRadiusRange(Float_t RadiusMin, Float_t RadiusMax) {
        fRadiusMin = RadiusMin;
        fRadiusMax = RadiusMax;
    };
    virtual void SetSeed(UInt_t /*seed*/) { ; }
    virtual void SetPtDistribution(TF1 &PtDistribution);

    void InitReactionsInfo();
    void SetReaction(Int_t NucleonPDG, std::vector<Int_t> DecayChannelPDG);

    // Pt-related
    Bool_t ComputeCdfTable(std::vector<Float_t> &Integral_vec, std::vector<Float_t> &Alpha_vec, std::vector<Float_t> &Beta_vec,
                           std::vector<Float_t> &Gamma_vec);
    Float_t GetRandomFromFunction(std::vector<Float_t> &Integral_vec, std::vector<Float_t> &Alpha_vec, std::vector<Float_t> &Beta_vec,
                                   std::vector<Float_t> &Gamma_vec, Float_t r);

    void PrintParameters();
    void PrintTitle();

   protected:
    Int_t fNucleonPDG;                    // pdg code of target/nucleon
    std::vector<Int_t> fDecayChannelPDG;  // vector that contains the pdg codes of the decay channel
    Float_t fRadiusMin;                  // minimum radius
    Float_t fRadiusMax;                  // maximum radius

    // Sets of vectors to describe possible decay channels
    // chosen to rely on the efficiency of std::set.count()
    std::set<std::vector<Int_t>> kSexaquarkNeutronProducts;
    std::set<std::vector<Int_t>> kSexaquarkProtonProducts;

    TF1 *fPtDistribution;  // model to choose random Pt from

    // integrate class for interactive use in ROOT
    ClassDef(AliGenSexaquarkReaction, 3);
};

#endif
