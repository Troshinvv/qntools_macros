#!/bin/bash

#SBATCH -D /lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/log
#SBATCH -p cascade
##SBATCH -t 24:00:00
#SBATCH -J QnTools
#SBATCH -a 1-155
#SBATCH --mem-per-cpu=4G
#SBATCH -o /lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/log/%A_%a.log
#SBATCH --exclude=n02p007,n02p006,n02p010,n02p039,n02p012,n02p021,n02p025,n02p003
list_dir=/lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/lambda/list_dir/
output_dir=/lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/lambda/out
id=$SLURM_ARRAY_TASK_ID
input_list=/lustre/home/user/v/vtroshin/bmn_hyperons/candidates_k_short_2.list
split -l 4 -d -a 4 --additional-suffix=.txt $input_list $list_dir
file_list=$( ls $list_dir | head -n $id | tail -n 1 )

mkdir -p $output_dir
cd $output_dir
mkdir $id
cd $id
qa=/lustre/home/user/v/vtroshin/bmn_hyperons/qn_bmn/lambda/out/qa.root
source /cvmfs/nica.jinr.ru/sw/os/login.sh
source /cvmfs/nica.jinr.ru/centos7/bmn/env.sh
source /cvmfs/nica.jinr.ru/centos7/bmnroot/dev/bmnroot_config.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/user/v/vtroshin/bmn_hyperons/PFSimple/install/lib/:/lustre/home/user/v/vtroshin/bmn_hyperons/PFSimple/install/external/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/user/v/vtroshin/bmn_hyperons/QnTools/install/lib:/lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/build/


# PLAIN
#time /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/build/correct /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/macro/k_short_correct.cc $list_dir/$file_list ${qa}
#time mv $output_dir/$id/qa.root $output_dir/qa_${id}.root

# RECENTERING
#time /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/build/correct /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/macro/lambda_correct.cc $list_dir/$file_list
# TWIST AND RESCALING
#time /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/build/correct /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/macro/lambda_correct.cc $list_dir/$file_list

time /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/build/correlate /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/macro/k_short_correlate.cc correction_out.root
time mv $output_dir/$id/corr.root $output_dir/corr_${id}.root

echo "The End."
