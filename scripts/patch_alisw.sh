#!/bin/bash

# `patch_alisw.sh`
# ===============
# Quick utility to patch AliRoot + AliDPG with current changes, before they are built.

set -euo pipefail

if [[ -z ${FSSR_ROOT_DIR:-} ]]; then echo "error: missing env. var. FSSR_ROOT_DIR"; exit 1; fi

if [[ $# -ne 2 ]]; then echo "usage: $0 <path aliroot> <path alidpg>"; exit 1; fi
aliroot_path=$1
alidpg_path=$2

cp -v "${FSSR_ROOT_DIR}/AliGenSexaquarkReaction.cxx" "${aliroot_path}/EVGEN/"
cp -v "${FSSR_ROOT_DIR}/AliGenSexaquarkReaction.h" "${aliroot_path}/EVGEN/"
cp -v "${FSSR_ROOT_DIR}/Hijing_Sexaquark.C" "${alidpg_path}/MC/CustomGenerators/PWGLF/"
