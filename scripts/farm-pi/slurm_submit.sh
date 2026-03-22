#!/bin/bash

# `slurm_submit.sh`
# =================
# Bash script to submit simulations jobs in batches to the Slurm farm.
# Output files will be:
# - ${FSSR_OUTPUT_DIR}/<new_mc_prod_name>/<r_channel><s_mass>/<run_number>/<dir_number>/AliESDs.root
# - ${FSSR_OUTPUT_DIR}/<new_mc_prod_name>/<r_channel><s_mass>/<run_number>/<dir_number>/galice.root
# - ${FSSR_OUTPUT_DIR}/<new_mc_prod_name>/<r_channel><s_mass>/<run_number>/<dir_number>/Kinematics.root
# - ${FSSR_OUTPUT_DIR}/<new_mc_prod_name>/<r_channel><s_mass>/<run_number>/<dir_number>/sim.log
# Log files will be:
# - ${FSSR_SLURM_DIR}/<new_mc_prod_name>_<r_channel><s_mass>/<run_number>_<dir_number>.log

set -euo pipefail

# hardcoded options #

MAX_PARALLEL_JOBS=96
N_DIR_NUMBERS=20

# functions #

print_usage() {
    echo "USAGE: "
    echo "  ./slurm_submit.sh --channel <reaction_channel> --mass <sexaquark_mass> --rn <rn1,rn2,...>"
    echo "  ./slurm_submit.sh --channel <reaction_channel> --mass <sexaquark_mass> --rn <rn_file>"
    echo "where:"
    echo "  <reaction_channel> = reaction to produce, it can be:"
    echo "                       * A : AntiS + N -> AntiL + K0"
    echo "                       * D : AntiS + P -> AntiL + Kp"
    echo "                       * E : AntiS + P -> AntiL + Kp + pim + pip"
    echo "                       * H : AntiS + N -> AntiP + Kp + Kp + pi0"
    echo "  <sexaquark_mass>   = mass of anti-sexaquark to generate, it can be:"
    echo "                       * 1.73"
    echo "                       * 1.8"
    echo "                       * 1.87"
    echo "                       * 1.94"
    echo "                       * 2.01"
    echo "  <rn1,rn2,...>      = run numbers to anchor, separated by comma, e.g."
    echo "                       --rn 297595,296623,246994"
    echo "  <rn_file>          = file that contains the run numbers to anchor, separated by new-lines, e.g."
    echo "                       --rn ../../doc/LHC15o_pass2.txt"
    echo "                       --rn ../../doc/LHC18qr_pass3.txt"
    echo "EXAMPLES:"
    echo "  ./slurm_submit.sh --channel A --mass 1.8 --rn 244917,244918,295585,295586,296690,296691"
    echo "  ./slurm_submit.sh --channel A --mass 1.8 --rn ../../doc/rn.txt"
}

process_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --channel) reaction_channel="$2"; shift 2 ;;
            --mass)    sexaquark_mass="$2";   shift 2 ;;
            --rn)      input_rn="$2";         shift 2 ;;
            *)         echo "error: unrecognized argument: $1"; print_usage; exit 1 ;;
        esac
    done
}

# check environment
if [[ -z ${ALIDPG_ROOT:-} ]]; then echo "error: missing env var ALIDPG_ROOT"; exit 1; fi
if [[ -z ${FSSR_ROOT_DIR:-} ]]; then echo "error: missing env var FSSR_ROOT_DIR"; exit 1; fi
if [[ -z ${FSSR_SLURM_DIR:-} ]]; then echo "error: missing env var FSSR_SLURM_DIR"; exit 1; fi
if [[ -z ${FSSR_OUTPUT_DIR:-} ]]; then echo "error: missing env var FSSR_OUTPUT_DIR"; exit 1; fi
if [[ -z ${LOCAL_OCDB_DIR:-} ]]; then echo "error: missing env var LOCAL_OCDB_DIR"; exit 1; fi

# command-line arguments
if [[ $# -ne 6 ]]; then print_usage; exit 1; fi
process_args "$@"

# validate input args
[[ " A D E H " == *" ${reaction_channel} "* ]] || { echo "error: invalid reaction channel"; exit 1; }
[[ " 1.73 1.8 1.87 1.94 2.01 " == *" ${sexaquark_mass} "* ]] || { echo "error: invalid sexaquark mass"; exit 1; }

# define strings (NOTE: not arrays, because Slurm)
export RUN_NUMBERS_STR
export DIR_NUMBERS_STR
DIR_NUMBERS_STR=$(seq -s " " -f "%03g" 1 "${N_DIR_NUMBERS}")

if [[ -e ${input_rn} ]]; then
    # rn case: file
    while read -r run_number; do
        RUN_NUMBERS_STR+="${run_number} "
    done < "${input_rn}"
else
    # rn case: comma-separated list
    RUN_NUMBERS_STR=${input_rn//,/ } # into a space-separated list
fi
n_total_rn=$(echo "${RUN_NUMBERS_STR}" | wc -w)
array_max=$((N_DIR_NUMBERS * n_total_rn - 1))

# prepare tmp dir for log hack
tmp_slurm_dir="${FSSR_SLURM_DIR}/tmp"
mkdir -p "${tmp_slurm_dir}"

sbatch \
    --output="${tmp_slurm_dir}/%A_%a.log" \
    --array="0-${array_max}%${MAX_PARALLEL_JOBS}" \
    -- "${FSSR_ROOT_DIR}/scripts/farm-pi/slurm_run.sh" "${reaction_channel}" "${sexaquark_mass}"

echo "$0 @ ${HOSTNAME} :: a total of $(( N_DIR_NUMBERS * n_total_rn )) jobs have been submitted"
