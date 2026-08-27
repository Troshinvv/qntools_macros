#!/bin/bash

#SBATCH -D /nica/mpd1/vtroshin/bmn_protons/preliminary/log
#SBATCH -p nica
#SBATCH -t 24:10:00
#SBATCH -J QnTools
#SBATCH -a 1-753
#SBATCH --mem-per-cpu=8G
#SBATCH -o /nica/mpd1/vtroshin/bmn_protons/preliminary/log/%A_%a.log
#SBATCH --exclude=ncx111,ncx112,ncx113,ncx115,ncx117,ncx121,ncx127,ncx169,ncx171,ncx172,ncx181,ncx185,ncx203,ncx207,ncx208,ncx211,ncx213,ncx215,ncx216,ncx217,ncx222,ncx223,ncx224,ncx225,ncx226,ncx227,ncx177
date
hostname
list_dir=/nica/mpd1/vtroshin/bmn_protons/preliminary/list_dir/
output_dir=/nica/mpd1/vtroshin/bmn_protons/preliminary/out_old_fix
id=$SLURM_ARRAY_TASK_ID
input_list=/lhep/users/vtroshin/run8_vf_24.04_all.list
split -l 36 -d -a 4 --additional-suffix=.txt $input_list $list_dir
#65_400
file_list=$( ls $list_dir | head -n $id | tail -n 1 )

efficiency_dir=/lhep/users/vtroshin/preliminary_qn/qntools_macros/efficiency

mkdir -p $output_dir
cd $output_dir
mkdir $id
cd $id
cp $efficiency_dir/efficiency.2024.04.03.root .
cp $output_dir/first_qa.root qa.root
sleep 10
source /cvmfs/nica.jinr.ru/sw/os/login.sh
source /cvmfs/bmn.jinr.ru/config/x86_64-centos7/cluster_config.sh
source /cvmfs/bmn.jinr.ru/bmnroot/25.09.0/x86_64-centos7/bmnroot_config.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lhep/users/vtroshin/QnTools/install/lib:/lhep/users/vtroshin/qntools_macros/build/
sleep 10
echo "Plain started"
time /lhep/users/vtroshin/qntools_macros/build/correct /lhep/users/vtroshin/preliminary_qn/qntools_macros/macro/run8_proton_correct_clean.cc $list_dir/$file_list efficiency.2024.04.03.root 
#echo "Correlate started"
time mv $output_dir/$id/qa.root $output_dir/qa_${id}.root
#time /lhep/users/vtroshin/qntools_macros/build/correlate /lhep/users/vtroshin/preliminary_qn/qntools_macros/macro/run8_proton_correlate_clean.cc correction_out.root
#time mv $output_dir/$id/corr.root $output_dir/corr_${id}.root
#time mv $output_dir/$id/out_qa.root $output_dir/out_qa_${id}.root
echo "The End."




