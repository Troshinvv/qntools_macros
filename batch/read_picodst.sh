#!/bin/bash

#SBATCH -D /nica/mpd1/vtroshin/bmn_protons/log
#SBATCH -p nica
#SBATCH -t 24:00:00
#SBATCH -J qa
#SBATCH -a 1-812
#SBATCH --mem-per-cpu=8G
#SBATCH -o /nica/mpd1/vtroshin/bmn_protons/log/%A_%a.log
#SBATCH --exclude=ncx111,ncx112,ncx113,ncx115,ncx117,ncx121,ncx127,ncx169,ncx171,ncx172,ncx181,ncx185,ncx203,ncx207,ncx208,ncx211,ncx213,ncx215,ncx216,ncx217,ncx222,ncx223,ncx224,ncx225,ncx226,ncx227
date
hostname
list_dir=/nica/mpd1/vtroshin/bmn_protons/list_dir/
output_dir=/nica/mpd1/vtroshin/bmn_protons/qa
id=$SLURM_ARRAY_TASK_ID
input_list=/lhep/users/vtroshin/run8_vf_25.09_140826.list
split -l 32 -d -a 4 --additional-suffix=.txt $input_list $list_dir
#65_400
file_list=$( ls $list_dir | head -n $id | tail -n 1 )

efficiency_dir=/lhep/users/vtroshin/qntools_macros/efficiency
calib_dir=/lhep/users/vtroshin/qntools_macros/centrality

mkdir -p $output_dir
cd $output_dir
mkdir $id
cd $id
cp $efficiency_dir/efficiency.2022.01.25.root .
cp $calib_dir/run8_25.09_corrections.root .
eff_file=efficiency.2022.01.25.root
cent_file=run8_25.09_corrections.root
sleep 10
source /cvmfs/nica.jinr.ru/sw/os/login.sh legacy
source /cvmfs/bmn.jinr.ru/config/x86_64-centos7/cluster_config.sh
source /cvmfs/bmn.jinr.ru/bmnroot/25.09.0/x86_64-centos7/bmnroot_config.sh
sleep 10
time root -l -b -q /lhep/users/vtroshin/qntools_macros/macro/final_push/run8_proton_qa.cc'("'$list_dir/$file_list'","'${eff_file}'","'${cent_file}'")'
time mv $output_dir/$id/rec_lambda_sim_qa.root $output_dir/rec_lambda_sim_qa_${id}.root 
echo "The End."




