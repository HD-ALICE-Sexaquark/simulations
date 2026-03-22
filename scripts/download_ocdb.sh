#!/bin/bash

# `download_ocdb.sh`
# ==================
# Quick utility bash script to download the OCDB files from the grid.

set -euo pipefail

if [[ -z ${ALIPHYSICS_VERSION:-} ]]; then echo "error: missing env var ALIPHYSICS_VERSION"; exit 1; fi
if [[ -z ${ALIROOT_VERSION:-} ]]; then echo "error: missing env var ALIROOT_VERSION"; exit 1; fi
if [[ -z ${ALIDPG_VERSION:-} ]]; then echo "error: missing env var ALIDPG_VERSION"; exit 1; fi

if [[ $# -ne 2 ]]; then echo "usage: $0 <run numbers file> <target dir>"; exit 1; fi
run_numbers_file=$1
target_ocdb_dir=$2

prod_equiv="LHC20e3"
prod_year=2020
if [[ $(basename "${run_numbers_file}") == "LHC15o_pass2.txt" ]]; then
    prod_equiv="LHC20j6"
fi

mkdir -p "${target_ocdb_dir}"

remote_ocdb_path="/alice/sim/${prod_year}/${prod_equiv}/OCDB"

# loop over run numbers
while read -r run_number; do
    mkdir -p "${target_ocdb_dir}/${prod_equiv}/${run_number}"

    for type in "sim" "rec"; do

        remote_file=${remote_ocdb_path}/${run_number}/OCDB${type}.root
        local_file=${target_ocdb_dir}/${prod_equiv}/${run_number}/OCDB${type}.root

        if ! alien.py ls "${remote_file}" &>/dev/null; then
            echo "warning: remote file not found: ${remote_file}, skipping"
            continue
        fi

        if [[ ! -e ${local_file} ]]; then
            alien.py cp alien://"${remote_file}" file://"${local_file}"
        else
            read -r md5_remote_file _ < <(alien.py md5sum "${remote_file}")
            if echo "${md5_remote_file}  ${local_file}" | md5sum -c --status; then
                echo "MD5(${remote_file}) == MD5(${local_file})"
            else
                alien.py cp alien://"${remote_file}" file://"${local_file}"
            fi
        fi
    done
done < "${run_numbers_file}"
