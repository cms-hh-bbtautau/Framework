#!/usr/bin/env bash

function run_cmd {
    "$@"
    RESULT=$?
    if (( $RESULT != 0 )); then
        echo "Error while running '$@'"
        kill -INT $$
    fi
}

do_install_cmssw() {
    local this_file="$( [ ! -z "$ZSH_VERSION" ] && echo "${(%):-%x}" || echo "${BASH_SOURCE[0]}" )"
    local this_dir="$( cd "$( dirname "$this_file" )" && pwd )"

    export SCRAM_ARCH=$1
    local CMSSW_VER=$2
    local mode=$3
    if ! [ -f "$this_dir/soft/$CMSSW_VER/.installed" ]; then
        run_cmd mkdir -p "$this_dir/soft"
        run_cmd cd "$this_dir/soft"
        run_cmd source /cvmfs/cms.cern.ch/cmsset_default.sh
        if [ -d $CMSSW_VER ]; then
            echo "Removing incomplete $CMSSW_VER installation..."
            run_cmd rm -rf $CMSSW_VER
        fi
        echo "Creating $CMSSW_VER area for in $PWD ..."
        run_cmd scramv1 project CMSSW $CMSSW_VER
        run_cmd cd $CMSSW_VER/src
        run_cmd eval `scramv1 runtime -sh`
        if [ "$mode" == "ana" ]; then
            run_cmd mkdir -p Framework/NanoProd/
            run_cmd ln -s "$this_dir/NanoProd" Framework/NanoProd/python
            run_cmd mkdir -p HHTools
            run_cmd ln -s "$this_dir/HHbtag" HHTools/HHbtag
            run_cmd mkdir -p TauAnalysis
            run_cmd ln -s "$this_dir/ClassicSVfit" TauAnalysis/ClassicSVfit
            run_cmd ln -s "$this_dir/SVfitTF" TauAnalysis/SVfitTF
            run_cmd mkdir -p HHKinFit2
            run_cmd ln -s "$this_dir/HHKinFit2" HHKinFit2/HHKinFit2
        fi
        if [ "$mode" == "comb" ]; then
            run_cmd git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
            run_cmd cd HiggsAnalysis/CombinedLimit
            run_cmd git checkout v9.1.0
            run_cmd cd ../..
            run_cmd bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/main/CombineTools/scripts/sparse-checkout-ssh.sh)
        fi
        run_cmd scram b -j8
        run_cmd touch "$this_dir/soft/$CMSSW_VER/.installed"
    fi
}

do_install_inference() {
    local this_file="$( [ ! -z "$ZSH_VERSION" ] && echo "${(%):-%x}" || echo "${BASH_SOURCE[0]}" )"
    local this_dir="$( cd "$( dirname "$this_file" )" && pwd )"

    if ! [ -f $this_dir/soft/git-lfs/.installed ]; then
        run_cmd mkdir -p "$this_dir/soft"
        run_cmd cd "$this_dir/soft"
        run_cmd wget https://github.com/git-lfs/git-lfs/releases/download/v3.3.0/git-lfs-linux-amd64-v3.3.0.tar.gz
        run_cmd tar xzf git-lfs-linux-amd64-v3.3.0.tar.gz
        run_cmd mkdir -p git-lfs
        run_cmd mv git-lfs-3.3.0/git-lfs git-lfs/
        run_cmd rm -r git-lfs-3.3.0 git-lfs-linux-amd64-v3.3.0.tar.gz
        run_cmd touch $this_dir/soft/git-lfs/.installed
    fi
    export PATH=$this_dir/soft/git-lfs:$PATH
    if ! [ -f "$this_dir/soft/inference/.installed" ]; then
        run_cmd mkdir -p "$this_dir/soft"
        run_cmd cd "$this_dir/soft"
        run_cmd source /cvmfs/cms.cern.ch/cmsset_default.sh
        if [ -d inference ]; then
            echo "Removing incomplete inference installation..."
            run_cmd rm -rf inference
        fi

        run_cmd git clone -b resonant --recursive ssh://git@gitlab.cern.ch:7999/acarvalh/inference.git
        run_cmd cd inference
        run_cmd source setup.sh
        run_cmd touch "$this_dir/soft/inference/.installed"
    fi
}

install_cmssw() {
    local this_file="$( [ ! -z "$ZSH_VERSION" ] && echo "${(%):-%x}" || echo "${BASH_SOURCE[0]}" )"
    local this_dir="$( cd "$( dirname "$this_file" )" && pwd )"
    local scram_arch=$1
    local cmssw_version=$2
    local os_version=$3
    local mode=$4
    if [[ $os_version < 8 ]] ; then
        local env_cmd=cmssw-cc$os_version
    else
        local env_cmd=cmssw-el$os_version
    fi
    if ! [ -f "$this_dir/soft/$CMSSW_VER/.installed" ]; then
        run_cmd $env_cmd --command-to-run /usr/bin/env -i HOME=$HOME bash "$this_file" install_cmssw $scram_arch $cmssw_version $mode
    fi
}

install_inference() {
    local this_file="$( [ ! -z "$ZSH_VERSION" ] && echo "${(%):-%x}" || echo "${BASH_SOURCE[0]}" )"
    local this_dir="$( cd "$( dirname "$this_file" )" && pwd )"
    local os_version=$1
    local mode=$4
    if [[ $os_version < 8 ]] ; then
        local env_cmd=cmssw-cc$os_version
    else
        local env_cmd=cmssw-el$os_version
    fi
    if ! [ -f "$this_dir/soft/inference/.installed" ]; then
        run_cmd $env_cmd --command-to-run /usr/bin/env -i HOME=$HOME bash "$this_file" install_inference
    fi
}

action() {
    # determine the directory of this file
    local this_file="$( [ ! -z "$ZSH_VERSION" ] && echo "${(%):-%x}" || echo "${BASH_SOURCE[0]}" )"
    local this_dir="$( cd "$( dirname "$this_file" )" && pwd )"

    export PYTHONPATH="$this_dir:$PYTHONPATH"
    export LAW_HOME="$this_dir/.law"
    export LAW_CONFIG_FILE="$this_dir/config/law.cfg"

    export ANALYSIS_PATH="$this_dir"
    export ANALYSIS_DATA_PATH="$ANALYSIS_PATH/data"
    export X509_USER_PROXY="$ANALYSIS_DATA_PATH/voms.proxy"
    export CENTRAL_STORAGE="/eos/home-k/kandroso/cms-hh-bbtautau"
    export VDAMANTE_STORAGE="/eos/home-v/vdamante/HH_bbtautau_resonant_Run2"
    export ANALYSIS_BIG_DATA_PATH="$CENTRAL_STORAGE/tmp/$(whoami)/data"

    run_cmd mkdir -p "$ANALYSIS_DATA_PATH"

    local os_version=$(cat /etc/os-release | grep VERSION_ID | sed -E 's/VERSION_ID="([0-9]+).*"/\1/')
    local default_cmssw_ver=CMSSW_13_0_15
    local cmb_cmssw_ver=CMSSW_11_3_4
    export DEFAULT_CMSSW_BASE="$ANALYSIS_PATH/soft/$default_cmssw_ver"
    export CMB_CMSSW_BASE="$ANALYSIS_PATH/soft/$cmb_cmssw_ver"

    run_cmd install_cmssw el9_amd64_gcc11 $default_cmssw_ver 9 ana
    run_cmd install_cmssw slc7_amd64_gcc10 $cmb_cmssw_ver 7 comb
    run_cmd install_inference 7

    if [ ! -z $ZSH_VERSION ]; then
        autoload bashcompinit
        bashcompinit
    fi
    source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.02/x86_64-almalinux9.3-gcc114-opt/bin/thisroot.sh
    export MAMBA_ROOT_PREFIX=/afs/cern.ch/work/k/kandroso/micromamba
    eval "$($MAMBA_ROOT_PREFIX/micromamba shell hook -s posix)"
    micromamba activate hh
    run_cmd source /afs/cern.ch/user/m/mrieger/public/law_sw/setup.sh
    source "$( law completion )"

    alias cmsEnv="env -i HOME=$HOME ANALYSIS_PATH=$ANALYSIS_PATH ANALYSIS_DATA_PATH=$ANALYSIS_DATA_PATH X509_USER_PROXY=$X509_USER_PROXY CENTRAL_STORAGE=$CENTRAL_STORAGE ANALYSIS_BIG_DATA_PATH=$ANALYSIS_BIG_DATA_PATH DEFAULT_CMSSW_BASE=$DEFAULT_CMSSW_BASE $ANALYSIS_PATH/RunKit/cmsEnv.sh"
    alias cmbEnv="env -i HOME=$HOME ANALYSIS_PATH=$ANALYSIS_PATH ANALYSIS_DATA_PATH=$ANALYSIS_DATA_PATH X509_USER_PROXY=$X509_USER_PROXY CENTRAL_STORAGE=$CENTRAL_STORAGE ANALYSIS_BIG_DATA_PATH=$ANALYSIS_BIG_DATA_PATH DEFAULT_CMSSW_BASE=$CMB_CMSSW_BASE /cvmfs/cms.cern.ch/common/cmssw-cc7 -- $ANALYSIS_PATH/RunKit/cmsEnv.sh"
    alias infEnv="env -i HOME=$HOME ANALYSIS_PATH=$ANALYSIS_PATH ANALYSIS_DATA_PATH=$ANALYSIS_DATA_PATH X509_USER_PROXY=$X509_USER_PROXY CENTRAL_STORAGE=$CENTRAL_STORAGE ANALYSIS_BIG_DATA_PATH=$ANALYSIS_BIG_DATA_PATH INF_PATH=$ANALYSIS_PATH/soft/inference /cvmfs/cms.cern.ch/common/cmssw-cc7 -- $ANALYSIS_PATH/inf_env.sh"
}

if [ "X$1" = "Xinstall_cmssw" ]; then
    do_install_cmssw "${@:2}"
elif [ "X$1" = "Xinstall_inference" ]; then
    do_install_inference "${@:2}"
else
    action "$@"
fi