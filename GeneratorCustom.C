AliGenerator *GeneratorCustom(Char_t reaction_channel = 'A', Double_t mass_sexaquark = 1.8)
{

    AliGenCocktail *ctl = (AliGenCocktail *)GeneratorCocktail("Hijing+AntiNeutron+AntiSexaquark");

    AliGenHijing *hij = (AliGenHijing *)GeneratorHijing();
    ctl->AddGenerator(hij, "Hijing", 1.);

    AliGenerator *ani = GeneratorInjector(40, -2112, 0., 6., -0.8, 0.8);
    ctl->AddGenerator(ani, "Injector (Anti-Neutron)", 1.);

    AliGenSexaquarkReaction *sr = new AliGenSexaquarkReaction(20, mass_sexaquark, reaction_channel);
    sr->SetPtRange(0., 5.);    // GeV/c
    sr->SetPhiRange(0., 360.); // degrees
    sr->SetYRange(-0.8, 0.8);
    sr->SetRadiusRange(5., 180.); // cm

    ctl->AddGenerator(sr, Form("Anti-Sexaquark Interaction (Channel %c)", reaction_channel), 1.);

    return ctl;
}
