#!/bin/bash

# `quick_test.sh`
# ===============
# Quick utility to test the simulations generation and reconstruction.

set -euo pipefail

if [[ -z ${ALIDPG_ROOT:-} ]]; then echo "error: missing env. var. ALIDPG_ROOT"; exit 1; fi
if [[ -z ${FSSR_ROOT_DIR:-} ]]; then echo "error: missing env. var. FSSR_ROOT_DIR"; exit 1; fi
if [[ -z ${LOCAL_OCDB_DIR:-} ]]; then echo "error: missing env. var. LOCAL_OCDB_DIR"; exit 1; fi

if [[ $# -ne 3 ]]; then echo "usage: $0 <run number> <r_channel> <s_mass>"; exit 1; fi
run_number=$1
r_channel=$2
s_mass=$3

# resolve metadata
pass_number="3"
data_prod="LHC18q"
data_prod_year="2018"
bkg_mc_prod="LHC20e3"
if grep -qx "${run_number}" "${FSSR_ROOT_DIR}/doc/LHC18r_pass3_rn.txt" 2>/dev/null; then
    data_prod="LHC18r"
elif grep -qx "${run_number}" "${FSSR_ROOT_DIR}/doc/LHC15o_pass2_rn.txt" 2>/dev/null; then
    pass_number="2"
    data_prod="LHC15o"
    data_prod_year="2015"
    bkg_mc_prod="LHC20j6"
fi

output_dir="test/${r_channel}${s_mass}/${run_number}/001"
mkdir -p "${output_dir}"

cd "${output_dir}"
ln -sf "${FSSR_ROOT_DIR}/DetectorCustom.C" .
ln -sf "${LOCAL_OCDB_DIR}/${bkg_mc_prod}/${run_number}/OCDBsim.root" .
ln -sf "${LOCAL_OCDB_DIR}/${bkg_mc_prod}/${run_number}/OCDBrec.root" .

# set env vars for metadata #

export ALIEN_JDL_LPMANCHORPASSNAME="pass${pass_number}"
export ALIEN_JDL_LPMRUNNUMBER="${run_number}"
export ALIEN_JDL_LPMPRODUCTIONTYPE="MC"
export ALIEN_JDL_LPMINTERACTIONTYPE="PbPb"
export ALIEN_JDL_LPMPRODUCTIONTAG="LHC26t"
export ALIEN_JDL_LPMANCHORRUN="${run_number}"
export ALIEN_JDL_LPMANCHORPRODUCTION="${data_prod}"
export ALIEN_JDL_LPMANCHORYEAR="${data_prod_year}"

export PRODUCTION_METADATA=""
PRODUCTION_METADATA+="LPMRawPass=${pass_number};"
PRODUCTION_METADATA+="LPMProductionType=${ALIEN_JDL_LPMPRODUCTIONTYPE};"
PRODUCTION_METADATA+="LPMProductionTag=${ALIEN_JDL_LPMPRODUCTIONTAG};"
PRODUCTION_METADATA+="LPMRunNumber=${ALIEN_JDL_LPMRUNNUMBER};"
PRODUCTION_METADATA+="LPMAnchorRun=${ALIEN_JDL_LPMRUNNUMBER};"
PRODUCTION_METADATA+="LPMAnchorProduction=${ALIEN_JDL_LPMANCHORPRODUCTION};"
PRODUCTION_METADATA+="LPMAnchorPassName=${ALIEN_JDL_LPMANCHORPASSNAME};"
PRODUCTION_METADATA+="LPMInteractionType=${ALIEN_JDL_LPMINTERACTIONTYPE};" # not needed by AliProdInfo
PRODUCTION_METADATA+="LPMAnchorYear=${ALIEN_JDL_LPMANCHORYEAR};" # not needed by AliProdInfo

export CONFIG_HLT=""

# main command #

${ALIDPG_ROOT}/bin/aliroot_dpgsim.sh --run ${run_number} \
                                     --mode sim,rec \
                                     --uid 1 \
                                     --nevents 10 \
                                     --generator PWGLF:Hijing_Sexaquark:${r_channel}${s_mass} \
                                     --simulation SimulationDefaultIonTail \
                                     --system Pb-Pb \
                                     --detector Custom
