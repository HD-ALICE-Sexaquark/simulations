AliGenerator *GeneratorCustom() {

    // cocktail
    AliGenCocktail *ctl = (AliGenCocktail *)GeneratorCocktail("Hijing+Sexaquark");

    // bkg
    AliGenHijing *hij = (AliGenHijing *)GeneratorHijing();
    hij->SetPtRange(0., 5.);
    hij->SetYRange(-0.8, 0.8);
    ctl->AddGenerator(hij, "Hijing", 1.);

    // signal
    AliGenSexaquarkReaction *sr = new AliGenSexaquarkReaction();
    sr->SetReaction();
    sr->SetPtRange(0., 5.);
    sr->SetPhiRange(0., 360.);
    sr->SetYRange(-0.8, 0.8);
    sr->SetRadiusRange(5., 180.);
    ctl->AddGenerator(sr, "Sexaquark", 1.);

    return ctl;
}
