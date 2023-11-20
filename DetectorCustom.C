void DetectorCustom() {

  Int_t year = atoi(gSystem->Getenv("CONFIG_YEAR"));
  UInt_t mask = strtol(gSystem->Getenv("CONFIG_DETECTORMASK"), 0, 16);

  iFRAME = 1;
  iABSO = 1;
  iDIPO = 1;
  iHALL = 1;
  iMAG = 1;
  iPIPE = 1;
  iSHIL = 1;
  iACORDE = 1;
  iAD = 0;
  iEMCAL = 1;
  iFMD = 1;
  iMUON = 1;
  iPHOS = 1;
  iPMD = 1;
  iZDC = 1;
  iTPC = 1;
  iTOF = 1;
  iITS = 1;
  iTRD = 1;
  iT0 = 1;
  iVZERO = 1;
  iHMPID = 1;
}
