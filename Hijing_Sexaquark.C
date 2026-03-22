AliGenerator *GeneratorCustom(TString opt = "") {

    // `opt` should be formatted as "<reaction_channel>:<mass_sexaquark>"
    // where:
    // - <reaction_channel> should be replaced by a single-digit char from 'A' to 'H'
    //   - 'A' : AntiSexaquark + Neutron -> AntiLambda, K0
    //   - 'D' : AntiSexaquark + Proton  -> AntiLambda, K+
    //   - 'H' : AntiSexaquark + Proton  -> AntiProton, K+, K+, Pi0
    // - <mass_sexaquark> corresponds to the mass of the injected anti-sexaquark
    //   - we're interested in: 1.73, 1.8, 1.87, 1.94, 2.01
    // examples of valid `opt`: "A:1.8", "D:1.73", "H:2.01", etc.

    // Parse option string //

    // -- default values
    char char_r_channel = 'A';
    double m_sexaquark = 1.8;  // (GeV/c^2)

    // -- retrieve options from separated by ':'
    if (opt != "") {
        TObjArray *tokens = opt.Tokenize(":");
        if (tokens->GetEntries() == 2) {
            char_r_channel = ((TObjString *)tokens->At(0))->GetString()[0];
            m_sexaquark = ((TObjString *)tokens->At(1))->GetString().Atof();
        }
        delete tokens;
    }
    Sexaquark::EReactionChannel r_channel = Sexaquark::ReactionChannelFromChar(char_r_channel);

    // Initialize generator //

    auto *sr = new AliGenSexaquarkReaction(100, m_sexaquark, r_channel);
    sr->SetPtRange(0., 10.);    // GeV/c
    sr->SetPhiRange(0., 360.);  // degrees
    sr->SetYRange(-0.8, 0.8);
    sr->SetRadiusRange(5., 180.);  // cm

    return sr;
}
