#ifndef ALIGENSEXAQUARKREACTION_H
#define ALIGENSEXAQUARKREACTION_H

/* Copyright(c) 1998-2023, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenerator.h"

/*
 Injector of anti-sexaquark into nucleon interactions
 * with an injected anti-sexaquark under a flat kinematic range
 * interaction vertex at under range
 * struck nucleon momentum assigned with simple gaussian

 -- based on /EVGEN/AliGenBox[.cxx/.h]

 author: aborquez@cern.ch
 */
class AliGenSexaquarkReaction : public AliGenerator {

    /* Copied from AliGenBox */

   public:
    AliGenSexaquarkReaction();
    AliGenSexaquarkReaction(Int_t NReactions, Float_t SexaquarkMass, Char_t ReactionChannel);
    virtual ~AliGenSexaquarkReaction();
    virtual void GenerateN(Int_t N);
    virtual void Generate();
    virtual void Init();
    virtual void SetSeed(UInt_t /*seed*/) { ; }

    /* Sexaquark */

   public:
    /*
     Set minimum and maximum radius for uniform distribution to locate interaction (or secondary) vertex
     */
    void SetDefaultRanges();
    void SetRadiusRange(Float_t radiusMin, Float_t radiusMax) {
        fRadiusMin = radiusMin;
        fRadiusMax = radiusMax;
    };
    void InitReactionsInfo();
    void PrintParameters();
    Int_t fCurrentReactionID;  // ID of current reaction
   private:
    std::map<Char_t, std::vector<std::vector<Int_t>>> fReactionChannelsMap;  // PDG codes of the products of the current reaction
    Float_t fRadiusMin;                                                      // minimum radius
    Float_t fRadiusMax;                                                      // maximum radius
    Int_t fStruckNucleonPDG;                                                 // PDG code of struck nucleon
    std::vector<Int_t> fReactionProductsPDG;                                 // vector of PDG codes of reaction products

    /* Fermi Momentum */

   public:
    void InitFermiMomentumInfo();
    void GetFermiMomentum(Float_t &px, Float_t &py, Float_t &pz);

   private:
    Float_t fFermiMomentum;       // central value of Fermi Momentum
    Float_t fFermiMomentumError;  // sigma value of Fermi Momentum
    TF1 *fFermiMomentumModel;

    /* Input Options */

   private:
    Int_t fNReactions;
    Float_t fSexaquarkMass;
    Char_t fReactionChannel;

    ClassDef(AliGenSexaquarkReaction, 3);  // integrate class for interactive use in ROOT
};

#endif
