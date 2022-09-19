#!/usr/bin/env bash

function run_cmd {
    "$@"
    RESULT=$?
    if (( $RESULT != 0 )); then
        echo "Error while running '$@'"
        kill -INT $$
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
    export CENTRAL_STORAGE="/eos/home-v/vdamante/HH_bbtautau_resonant_Run2"
    export ANALYSIS_BIG_DATA_PATH="$CENTRAL_STORAGE/tmp/$(whoami)/data"
    export PATH=$PATH:$HOME/.local/bin:$ANALYSIS_PATH/scripts

    # source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_101 x86_64-centos7-gcc8-opt
    # conda activate bbtt

    local os_version=$(cat /etc/os-release | grep VERSION_ID | sed -E 's/VERSION_ID="([0-9]+)"/\1/')
    if [ $os_version = "7" ]; then
        export SCRAM_ARCH=slc7_amd64_gcc10
    else
        export SCRAM_ARCH=el8_amd64_gcc10
    fi
    CMSSW_VER=CMSSW_12_4_8
    if ! [ -f "$this_dir/soft/CentOS$os_version/$CMSSW_VER/.installed" ]; then
        run_cmd mkdir -p "$this_dir/soft/CentOS$os_version"
        run_cmd cd "$this_dir/soft/CentOS$os_version"
        if [ -d $CMSSW_VER ]; then
            echo "Removing incomplete $CMSSW_VER installation..."
            run_cmd rm -rf $CMSSW_VER
            run_cmd rm -f "$this_dir/soft/CentOS$os_version/bin/python"
        fi
        echo "Creating $CMSSW_VER area for CentOS$os_version in $PWD ..."
        run_cmd scramv1 project CMSSW $CMSSW_VER
        run_cmd cd $CMSSW_VER/src
        run_cmd eval `scramv1 runtime -sh`
        run_cmd mkdir -p Framework/NanoProd/
        run_cmd ln -s "$this_dir/NanoProd" Framework/NanoProd/python
        run_cmd scram b -j8
        run_cmd cd "$this_dir"
        run_cmd mkdir -p "$this_dir/soft/CentOS$os_version/bin"
        run_cmd ln -s $(which python3) "$this_dir/soft/CentOS$os_version/bin/python"
        touch "$this_dir/soft/CentOS$os_version/$CMSSW_VER/.installed"
    else
        run_cmd cd "$this_dir/soft/CentOS$os_version/$CMSSW_VER/src"
        run_cmd eval `scramv1 runtime -sh`
        run_cmd cd "$this_dir"
    fi
    mkdir -p "$ANALYSIS_DATA_PATH"
    export PATH="$this_dir/soft/CentOS$os_version/bin:$PATH"
    source /afs/cern.ch/user/m/mrieger/public/law_sw/setup.sh
    source "$( law completion )"
}
action
