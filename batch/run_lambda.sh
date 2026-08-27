#!/bin/bash

#SBATCH -D /scratch2/troshin/qn/log
#SBATCH -p nica
##SBATCH -t 24:00:00
#SBATCH -J QnTools
#SBATCH -a 1-300
#SBATCH --mem-per-cpu=4G
#SBATCH -o /scratch2/troshin/qn/log/%A_%a.log
#SBATCH --exclude=ncx112,ncx114,ncx115,ncx121,ncx127,ncx132,ncx145,ncx146,ncx148,ncx151,ncx153,ncx154,ncx155,ncx156,ncx157,ncx158,ncx159,ncx160,ncx161,ncx163,ncx164,ncx165,ncx166,ncx168,ncx171,ncx172,ncx175,ncx181,ncx184,ncx206,ncx212,ncx214,ncx216,ncx222,ncx223,ncx225,ncx227,ncx228
date
hostname
list_dir=/scratch2/troshin/qn/list_dir/
output_dir=/scratch2/troshin/qn/out
id=$SLURM_ARRAY_TASK_ID
#NP_lambda_candidates_1040.list
input_list=/lhep/users/vtroshin/candidates.list
split -l 1 -d -a 4 --additional-suffix=.txt $input_list $list_dir
file_list=$( ls $list_dir | head -n $id | tail -n 1 )
##107_6
mkdir -p $output_dir
cd $output_dir
mkdir $id
cd $id
qa=/scratch2/troshin/qn/out/qa.root
eff_file=/lhep/users/vtroshin/JAM_bmn_eff_map.root
source /cvmfs/nica.jinr.ru/sw/os/login.sh legacy
module add mpddev
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lhep/users/vtroshin/PFSimple/install/lib/:/lhep/users/vtroshin/PFSimple/install/external/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lhep/users/vtroshin/QnTools/install/lib:/lhep/users/vtroshin/qntools_macros/build/


# PLAIN
time /lhep/users/vtroshin/qntools_macros/build/correct /lhep/users/vtroshin/qntools_macros/macro/lambda_correct.cc $list_dir/$file_list ${eff_file} ${qa}
time mv $output_dir/$id/qa.root $output_dir/qa_${id}.root
# RECENTERING
#time /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/build/correct /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/macro/lambda_correct.cc $list_dir/$file_list
# TWIST AND RESCALING
#time /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/build/correct /lustre/home/user/v/vtroshin/bmn_hyperons/qntools_macros/macro/lambda_correct.cc $list_dir/$file_list

time /lhep/users/vtroshin/qntools_macros/build/correlate /lhep/users/vtroshin/qntools_macros/macro/lambda_correlate.cc correction_out.root
time mv $output_dir/$id/corr.root $output_dir/corr_${id}.root

echo "The End."
