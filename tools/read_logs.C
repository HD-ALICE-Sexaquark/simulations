void read_logs(TString server = "S10") {  // "S10" or "S12" or "S14"

    TString pathToA = "../samples/";
    TString line1 = "SIMULATION TIME:";
    TString line2 = "RECONSTRUCTION TIME:";
    TString line3 = "QA TRAIN TIME:";

    TSystemDirectory dirA("dirA", pathToA);
    TList *listOfDirA = dirA.GetListOfFiles();
    if (!listOfDirA) return;

    TIter nextDirA(listOfDirA);
    TObject *objDirA = nullptr;

    Double_t n_valid_files = 0.;
    Double_t seconds_valid_files = 0.;

    while ((objDirA = nextDirA())) {

        TString dirName = objDirA->GetName();
        if (dirName(dirName.Length() - 3, dirName.Length() - 1) != server) continue;

        TString pathToB = pathToA + objDirA->GetName() + "/";

        TSystemDirectory dirB("dirB", pathToB);
        TList *listOfDirB = dirB.GetListOfFiles();
        if (!listOfDirB) return;

        TIter nextDirB(listOfDirB);
        TObject *objDirB = nullptr;

        while ((objDirB = nextDirB())) {

            TString logFilePath = pathToB + objDirB->GetName() + "/stdout.log";
            // std::cout << "File: " << logFilePath << std::endl; // debug

            std::ifstream logFile(logFilePath);
            if (!logFile.is_open()) {
                // std::cerr << "Unable to open file: " << logFilePath << std::endl; // debug
                continue;
            }

            std::string str_line;
            Bool_t foundLine1 = kFALSE;
            Bool_t foundLine2 = kFALSE;
            Bool_t foundLine3 = kFALSE;

            Double_t seconds_sim = 0.;
            Double_t seconds_reco = 0.;
            Double_t seconds_qa = 0.;
            Double_t seconds_total = 0.;

            while (std::getline(logFile, str_line)) {

                TString tstr_line = str_line;
                if (tstr_line.Contains(line1)) {
                    foundLine1 = kTRUE;
                    std::istringstream iss(str_line);
                    std::string token;
                    std::getline(iss, token, ':');
                    std::getline(iss, token, ':');
                    seconds_sim = std::stod(token);
                    // std::cout << "-> token = " << std::stod(token) << std::endl;
                } else if (tstr_line.Contains(line2)) {
                    foundLine2 = kTRUE;
                    std::istringstream iss(str_line);
                    std::string token;
                    std::getline(iss, token, ':');
                    std::getline(iss, token, ':');
                    seconds_reco = std::stod(token);
                    // std::cout << "-> token = " << std::stod(token) << std::endl;
                } else if (tstr_line.Contains(line3)) {
                    foundLine3 = kTRUE;
                    std::istringstream iss(str_line);
                    std::string token;
                    std::getline(iss, token, ':');
                    std::getline(iss, token, ':');
                    seconds_qa = std::stod(token);
                    // std::cout << "-> token = " << std::stod(token) << std::endl;
                }

                // std::cout << "-> count = " << std::stoi(token) << std::endl;
            }

            seconds_total = seconds_sim + seconds_reco + seconds_qa;

            logFile.close();

            if (foundLine1 && foundLine2 && foundLine3) {
                // std::cout << "The three lines have been found in file: " << logFilePath << std::endl;
                n_valid_files += 1.;
                seconds_valid_files += seconds_total;
            } else {
                std::cout << logFilePath << std::endl;
            }

        }  // end of loop over dirs B
    }      // end of loop over dirs A

    std::cout << "n_valid_files = " << n_valid_files << std::endl;
    std::cout << "seconds_valid_files = " << seconds_valid_files << std::endl;
    std::cout << "seconds_per_file = " << seconds_valid_files / n_valid_files << std::endl;
    std::cout << "seconds_per_event = " << seconds_valid_files / n_valid_files / 10. << std::endl;
    std::cout << "n_desired_events * seconds_per_event = " << 7.4E6 * seconds_valid_files / n_valid_files / 10. << std::endl;
    std::cout << "n_desired_events * days_per_event = " << 7.4E6 * seconds_valid_files / n_valid_files / 10. / 60. / 60. / 24. << std::endl;
    std::cout << "CPU time = " << 7.4E6 * seconds_valid_files / n_valid_files / 10. / 60. / 60. / 24. / 10000. << " days per 10k CPUs"
              << std::endl;
}
