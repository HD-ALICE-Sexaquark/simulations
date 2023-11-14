#include "TDirectoryFile.h"
#include "TFile.h"
#include "TParticle.h"
#include "TTree.h"

/*
 Small macro to read Kinematics.root files
 */
void parse_kinematics(TString input_filename = "Kinematics.root",  //
                      TString dirName = "Event0",                  //
                      TString treeName = "TreeK", ) {

    TFile* file = TFile::Open(input_filename);
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Cannot open file " << input_filename << std::endl;
        return;
    }

    TDirectoryFile* dir = dynamic_cast<TDirectoryFile*>(file->Get(dirName));
    if (!dir) {
        std::cerr << "ERROR: Cannot find directory " << dirName << " in file " << input_filename << std::endl;
        file->Close();
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(dir->Get(treeName));
    if (!tree) {
        std::cerr << "ERROR: Cannot find tree " << treeName << " in directory " << dirName << std::endl;
        file->Close();
        return;
    }

    TParticle* particle = new TParticle();  // assuming TParticle is the class used in the TTree

    tree->SetBranchAddress("Particles", &particle);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        particle->Print();
    }

    file->Close();
}
