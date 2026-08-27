#!/bin/bash

#SBATCH -D /lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/log
#SBATCH -p mephi
##SBATCH -t 24:00:00
#SBATCH -J QnTools
#SBATCH -a 1-1
#SBATCH --mem-per-cpu=64G
#SBATCH -o /lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/log/%A_%a.log
#SBATCH --exclude=n02p002,n02p069,n02p064
date
hostname
output_dir=/lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/lambda/out
cd $output_dir
qa=/lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/lambda/out/qa.root
source /cvmfs/nica.jinr.ru/sw/os/login.sh
source /cvmfs/nica.jinr.ru/centos7/bmn/env.sh
source /cvmfs/nica.jinr.ru/centos7/bmnroot/23.08.0/bmnroot_config.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/user/v/vtroshin/bmn_hyperons/PFSimple/install/lib/:/lustre/home/user/v/vtroshin/bmn_hyperons/PFSimple/install/external/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/user/v/vtroshin/bmn_hyperons/QnTools/install/lib:/lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/build/


# PLAIN
time hadd -n 2 /lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/lambda/out/qa.root qa_*

echo "The End."
