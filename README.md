simulations
===========

## Environment

- modern installation and environment of **AliRoot**, **AliPhysics** and **AliDPG**
- `ALICE_DATA` -- path to ALICE large data files ([reference](https://alisw.github.io/git-advanced/#how-to-use-large-data-files-for-analysis))
- `FSSR_ROOT_DIR` -- main repository dir
- `FSSR_TASK_DIR` -- `${FSSR_ROOT_DIR}/qa`
- `FSSR_SLURM_DIR` -- `${FSSR_ROOT_DIR}/slurm`
- `FSSR_OUTPUT_DIR` -- `${FSSR_ROOT_DIR}/output`

## Files

`AliGenSexaquarkReaction.cxx`, `AliGenSexaquarkReaction.h` -- Must be copied into the `EVGEN/` directory of [**AliRoot**](https://github.com/alisw/AliRoot)

`Hijing_Sexaquark.C` -- Must be copied into the `MC/CustomGenerators/PWGLF/` directory of [**AliDPG**](https://github.com/alisw/AliDPG).

`DetectorCustom.C` -- Custom detector configuration. It turns off most detectors.

`scripts/download_ocdb.sh` -- Must be executed before `quick_test.sh` and the `farm-pi/` scripts. It will create and fill a `${FSSR_ROOT_DIR}/ocdb` directory.

`scripts/quick_test.sh` -- (requires `download_ocdb.sh`) Interactive quick test of MC production.

`scripts/farm-pi/` -- (requires `download_ocdb.sh`) Scripts to send batch simulations via Slurm.

`qa/` -- ALICE analysis task to do quick QA of simulations.

`doc/` -- List of run numbers.

`scripts/releaser.sh` -- Utility script to have working in a screen/tmux session to release stuck jobs in Slurm.

`scripts/patch_alisw.sh` -- Utility script that copies `AliGenSexaquarkReaction.cxx`+`AliGenSexaquarkReaction.h` and `Hijing_Sexaquark.C` into the `AliRoot` and `AliDPG` directories, respectively (before being built)
