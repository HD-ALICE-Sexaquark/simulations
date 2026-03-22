simulations
===========

## Environment

- modern installation and environment of [**AliRoot**](https://github.com/alisw/AliRoot), [**AliPhysics**](https://github.com/alisw/AliPhysics) and [**AliDPG**](https://github.com/alisw/AliDPG)
- `ALICE_DATA` -- path to ALICE large data files ([reference](https://alisw.github.io/git-advanced/#how-to-use-large-data-files-for-analysis))
- `FSSR_ROOT_DIR` -- main repository dir
- `FSSR_SLURM_DIR` -- `${FSSR_ROOT_DIR}/slurm`
- `FSSR_OUTPUT_DIR` -- `${FSSR_ROOT_DIR}/output`
- `FSSR_TASK_DIR` -- `${FSSR_ROOT_DIR}/qa`, required by `qa/exec.sh`
- `LOCAL_OCDB_DIR` -- directory that contains OCDB files, filled by `scripts/download_ocdb.sh`

## Files

`AliGenSexaquarkReaction.cxx`, `AliGenSexaquarkReaction.h` -- Must be copied into the `EVGEN/` directory of [**AliRoot**](https://github.com/alisw/AliRoot)

`Hijing_Sexaquark.C` -- Must be copied into the `MC/CustomGenerators/PWGLF/` directory of [**AliDPG**](https://github.com/alisw/AliDPG).

`DetectorCustom.C` -- Custom detector configuration. It turns off most detectors.

`scripts/download_ocdb.sh` -- Must be executed before `quick_test.sh` and the `farm-pi/` scripts. It will create and fill a `${LOCAL_OCDB_DIR}` directory.

`scripts/quick_test.sh` -- (requires `download_ocdb.sh`) Interactive quick test of MC production.

`scripts/farm-pi/`, `scripts/farm-gsi/` -- (require `download_ocdb.sh`) Scripts to send batch simulations via Slurm in local or GSI cluster, respectively

`scripts/patch_alisw.sh` -- Utility script that copies `AliGenSexaquarkReaction.cxx`+`AliGenSexaquarkReaction.h` and `Hijing_Sexaquark.C` into the `AliRoot` and `AliDPG` directories, respectively (before being built)

`scripts/releaser.sh` -- Utility script to have working in a screen/tmux session to release stuck jobs in Slurm.

`qa/` -- ALICE analysis task to do quick QA of simulations.

`doc/` -- Lists of run numbers.

