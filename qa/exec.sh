#!/bin/bash

if [[ -z ${FSSR_TASK_DIR} ]]; then echo "missing FSSR_TASK_DIR" ; exit 1; fi

# hardcoded
export ATTEMPT_NAME="from_prod"
export SIMULATION_SET="A1.8"
export INPUT_PATH="${HOME}/some/alice_newsim/simulations/output/LHC26a/A1.8"
# export INPUT_PATH="${HOME}/some/alice_newsim/simulations/scripts"
# export INPUT_PATH="${HOME}/some/alice_newsim/simulations/scripts/test"
export RUN_NUMBER=295585
export LOCAL_N_DIRS=10
# export LOCAL_N_DIRS=1
export LOCAL_LIMIT_N_EVENTS=0

echo "exec.sh :: INPUT_PATH           = \"${INPUT_PATH}\""
echo "exec.sh :: RUN_NUMBER           = ${RUN_NUMBER}"
echo "exec.sh :: LOCAL_N_DIRS         = ${LOCAL_N_DIRS}"
echo "exec.sh :: LOCAL_LIMIT_N_EVENTS = ${LOCAL_LIMIT_N_EVENTS}"

attempt_dir=${FSSR_TASK_DIR}/attempts/${ATTEMPT_NAME}
mkdir -p "${attempt_dir}"

cp "${FSSR_TASK_DIR}/AliAnalysisTaskSexaReactionQA.cxx" ${attempt_dir}/
cp "${FSSR_TASK_DIR}/AliAnalysisTaskSexaReactionQA.h" ${attempt_dir}/
cp "${FSSR_TASK_DIR}/AddTaskSexaReactionQA.C" ${attempt_dir}/
cp "${FSSR_TASK_DIR}/runAnalysis.C" ${attempt_dir}/

cd "${attempt_dir}" || exit

analysis_options="("
analysis_options+="\"${INPUT_PATH}\","
analysis_options+="${RUN_NUMBER},"
analysis_options+="${LOCAL_N_DIRS},"
analysis_options+="${LOCAL_LIMIT_N_EVENTS}"
analysis_options+=")"

aliroot_command="aliroot -l -b -q runAnalysis.C${analysis_options}"
echo "${aliroot_command}"
${aliroot_command} 2>&1 | tee analysis.log

cd "${FSSR_TASK_DIR}" || exit
