#!/bin/bash

###########################################
##                                       ##
##  Script to send sexaquark simulation  ##
##      jobs to the ht_condor farm       ##
##                                       ##
###########################################

# 21.Nov.2023
## A. Borquez

#############
# Functions #
#############

function print_help() {
    echo "SCRIPT: send_production.sh"
    echo "=========================="
    echo "./send_production.sh --N <N> --channel <channel> --mass <sexaquark_mass> --rn <rn1,rn2,...> --serv <serv>"
    echo "where:"
    echo "  <N>              = amount of events to generate in a single job"
    echo "  <channel>        = reaction to produce, it can be:"
    echo "                     * A : AntiS + N -> AntiL + K0"
    echo "                     * D : AntiS + P -> AntiL + Kp"
    echo "                     * E : AntiS + P -> AntiL + Kp + pim + pip"
    echo "                     * H : AntiS + N -> AntiP + Kp + Kp + pi0"
    echo "  <sexaquark_mass> = mass of anti-sexaquark to generate, it can be:"
    echo "                     * 1.73"
    echo "                     * 1.8"
    echo "                     * 1.87"
    echo "                     * 1.94"
    echo "                     * 2.01"
    echo "  <rn>             = run numbers to anchor, separated by comma, for example:"
    echo "                     --rn 297595,296623,246994"
    echo "                     (IMPORTANT: make sure you store the respective <RUN_NUMBER>_sim.root and <RUN_NUMBER>_rec.root files inside tools/)"
    echo "  <serv>           = select which machine to use: alice-serv<serv>, it can be:"
    echo "                     * 10"
    echo "                     * 12"
    echo "                     * 13"
    echo "                     * 14"
}

function process_args() {
    arr=("$@")
    ic=0
    while [[ $ic -le $((${#arr[@]}-1)) ]]; do
        if [[ "${arr[$ic]}" == "--N" ]]; then
            N_EVENTS=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--channel" ]]; then
            CHANNEL=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--mass" ]]; then
            SEXA_MASS=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--rn" ]]; then
            RUN_NUMBER_LIST=${arr[$((ic+1))]}
            RUN_NUMBER_ARR=(${RUN_NUMBER_LIST//,/ }) # convert comma-separated list into array of numbers
        elif [[ "${arr[$ic]}" == "--serv" ]]; then
            SERV=${arr[$((ic+1))]}
        else
            echo "ERROR: unrecognized argument: ${arr[$((ic))]}."
            print_help
            exit 1
        fi
        ((ic+=2))
    done
}

function prepare_sim_script() {

    RUN_NUMBER=${1}

    # copy env. scripts
    cp ${SIMDIR}/tools/alice_env.sh .

    # copy OCDB files
    cp ${SIMDIR}/tools/${RUN_NUMBER}_OCDBrec.root OCDBrec.root
    cp ${SIMDIR}/tools/${RUN_NUMBER}_OCDBsim.root OCDBsim.root

    # create run/command file (within OUTDIR)
    run_file="run.sh"
    echo "#!/bin/bash"                                     > ${run_file}
    echo ""                                               >> ${run_file}
    echo "export ALIBUILD_WORK_DIR=${ALIBUILD_WORK_DIR}"  >> ${run_file}
    echo 'eval `/usr/local/bin/alienv shell-helper`'      >> ${run_file}
    echo "export PATH=${HOME}/.local/bin:${PATH}"         >> ${run_file}
    echo ""                                               >> ${run_file}
    echo "source alice_env.sh"                            >> ${run_file}
    echo ""                                               >> ${run_file}

    # prepare final command to run simulations
    echo '${ALIDPG_ROOT}/bin/aliroot_dpgsim.sh --run '${RUN_NUMBER}' --mode sim,rec,qa --uid 1 --nevents '${N_EVENTS}' --generator PWGLF:Hijing_Sexaquark:'${CHANNEL}${SEXA_MASS}' --simulation SimulationDefaultIonTail --detector NoAD --system Pb-Pb' >> ${run_file}
    echo ""                                               >> ${run_file}

    # compress GRP/ into a tar (condor being condor...)
    # and delete it, to save some space
    echo "tar czf GRP.tar.gz GRP/"                        >> ${run_file}
    echo "rm -rfv GRP/"                                   >> ${run_file}
    chmod a+x ${run_file}

    # print info about output directory and copied files (debug purposes)
    echo "${PWD}"
    for file in *.*; do echo "  ${file}"; done
}

function prepare_job() {

    # create job file (within OUTDIR)
    job_file="${1}"
    echo "executable            = run.sh"                                   > ${job_file}
    echo "output                = stdout.log"                              >> ${job_file}
    echo "error                 = stderr.log"                              >> ${job_file}
    echo "log                   = log.condor"                              >> ${job_file}
    echo "should_transfer_files = YES"                                     >> ${job_file}
    echo "transfer_input_files  = alice_env.sh,OCDBrec.root,OCDBsim.root"  >> ${job_file}
    echo 'requirements          = (Machine == "alice-serv'${SERV}'")'      >> ${job_file}
    echo "queue" >> ${job_file}
}

########
# Main #
########

# define main sim. directory
if [[ ${PWD##*/} != "simulations" ]]; then
    echo "ERROR: you must be located at the simulations/ dir."
fi
SIMDIR=$(readlink -f ${PWD})

# check environment
if [[ -z ${ALIBUILD_WORK_DIR} ]]; then
    echo "ERROR: please set your ALIBUILD_WORK_DIR env."
fi

# create env. script
if [[ ! -e tools/alice_env.sh ]]; then
    alienv printenv AliDPG/latest,AliRoot/latest,AliPhysics/latest > tools/alice_env.sh
    chmod a+x tools/alice_env.sh
fi

# process input
argArray=("$@")
process_args "${argArray[@]}"

# start loop, each run corresponds to a different job
for CURRENT_RN in ${RUN_NUMBER_ARR[@]}; do

    # create and enter output dir
    # (important to do so, because condor requires input files to be in the same dir)
    OUTDIR=${SIMDIR}/samples/${CHANNEL}${SEXA_MASS}_S${SERV}/${CURRENT_RN}
    mkdir -p ${OUTDIR}

    cd ${OUTDIR}

    CURRENT_JOB_NAME="${CHANNEL}${SEXA_MASS}_${CURRENT_RN}"

    prepare_sim_script ${CURRENT_RN}
    prepare_job ${CURRENT_JOB_NAME}.condor

    # submit job to the farm
    condor_submit ${CURRENT_JOB_NAME}.condor
done

cd ${SIMDIR}
