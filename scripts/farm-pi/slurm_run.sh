#!/bin/bash

# `slurm_run.sh`
# ==============
# NOTE: don't execute this script directly, it is meant to be used by `farm-pi/slurm_submit.sh`

#SBATCH --partition=main
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

set -euo pipefail

# check environment
if [[ -z ${ALIDPG_ROOT:-} ]]; then echo "error: missing env var ALIDPG_ROOT"; exit 1; fi
if [[ -z ${FSSR_ROOT_DIR:-} ]]; then echo "error: missing env var FSSR_ROOT_DIR"; exit 1; fi
if [[ -z ${FSSR_OUTPUT_DIR:-} ]]; then echo "error: missing env var FSSR_OUTPUT_DIR"; exit 1; fi
if [[ -z ${FSSR_SLURM_DIR:-} ]]; then echo "error: missing env var FSSR_SLURM_DIR"; exit 1; fi
if [[ -z ${LOCAL_OCDB_DIR:-} ]]; then echo "error: missing env var LOCAL_OCDB_DIR"; exit 1; fi
# -- per run number options
if [[ -z ${RUN_NUMBERS_STR:-} ]]; then echo "error: missing env var RUN_NUMBERS_STR"; exit 1; fi
if [[ -z ${DIR_NUMBERS_STR:-} ]]; then echo "error: missing env var DIR_NUMBERS_STR"; exit 1; fi

if [[ $# -ne 2 ]]; then echo "error: needs 2 arguments <r_channel> <s_mass>"; exit 1; fi
r_channel=$1
s_mass=$2

# convert strings into arrays (NOTE: because Slurm)
read -ra RUN_NUMBERS_ARR <<< "${RUN_NUMBERS_STR}"
read -ra DIR_NUMBERS_ARR <<< "${DIR_NUMBERS_STR}"

# pick specific rn+dn pair
n_runs=${#RUN_NUMBERS_ARR[@]}
run_idx=$(( SLURM_ARRAY_TASK_ID % n_runs ))
dir_idx=$(( SLURM_ARRAY_TASK_ID / n_runs ))
run_number=${RUN_NUMBERS_ARR[${run_idx}]}
dir_number=${DIR_NUMBERS_ARR[${dir_idx}]}

# resolve gen. purp. and new mc production names
pass_number="3"
data_prod="LHC18q"
data_prod_year="2018"
bkg_mc_prod="LHC20e3"
new_mc_prod="LHC26a"
if grep -qx "${run_number}" "${FSSR_ROOT_DIR}/doc/LHC18r_pass3_rn.txt" 2>/dev/null; then
    data_prod="LHC18r"
elif grep -qx "${run_number}" "${FSSR_ROOT_DIR}/doc/LHC15o_pass2_rn.txt" 2>/dev/null; then
    pass_number="2"
    data_prod="LHC15o"
    data_prod_year="2015"
    bkg_mc_prod="LHC20j6"
    new_mc_prod="LHC26b"
fi

output_dir=${FSSR_OUTPUT_DIR}/${new_mc_prod}/${r_channel}${s_mass}/${run_number}/${dir_number}
mkdir -p "${output_dir}"

cd "${output_dir}"

ln -sf "${FSSR_ROOT_DIR}/DetectorCustom.C" .
ln -sf "${LOCAL_OCDB_DIR}/${bkg_mc_prod}/${run_number}/OCDBsim.root" .
ln -sf "${LOCAL_OCDB_DIR}/${bkg_mc_prod}/${run_number}/OCDBrec.root" .

# log hack (https://unix.stackexchange.com/a/585453)
tmp_logfile=${FSSR_SLURM_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log # = ${FSSR_SLURM_DIR}/tmp/%A_%a.log
slurm_subdir=${FSSR_SLURM_DIR}/${new_mc_prod}_${r_channel}${s_mass}
mkdir -p "${slurm_subdir}"
ln -f "${tmp_logfile}" "${slurm_subdir}/${run_number}_${dir_number}.log"

# set extra env vars #

export ALIEN_JDL_LPMANCHORPASSNAME="pass${pass_number}"
export ALIEN_JDL_LPMRUNNUMBER="${run_number}"
export ALIEN_JDL_LPMPRODUCTIONTYPE="MC"
export ALIEN_JDL_LPMINTERACTIONTYPE="PbPb"
export ALIEN_JDL_LPMPRODUCTIONTAG="${new_mc_prod}"
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

${ALIDPG_ROOT}/bin/aliroot_dpgsim.sh --run "${run_number}" \
                                     --mode sim,rec \
                                     --uid "${SLURM_ARRAY_JOB_ID}${SLURM_ARRAY_TASK_ID}" \
                                     --nevents 10 \
                                     --generator "PWGLF:Hijing_Sexaquark:${r_channel}${s_mass}" \
                                     --simulation SimulationDefaultIonTail \
                                     --system Pb-Pb \
                                     --detector Custom

# remove symbolic links
rm -v DetectorCustom.C OCDBsim.root OCDBrec.root

# remove intermediate files
rm -rfv GRP/ AliESDfriends.root geometry.root ./*.dat grpdump.sh Merged.QA.Data.root QA* Run*.root TOFQA.root Trigger.root rec.log simwatch.log recwatch.log

# end of log hack
rm "${tmp_logfile}"
