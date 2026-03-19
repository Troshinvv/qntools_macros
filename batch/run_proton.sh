#!/bin/bash

#SBATCH -D /scratch2/troshin/bmn_protons/qn_bmn/log
#SBATCH -p nica
#SBATCH -t 24:00:00
#SBATCH -J QnTools
#SBATCH -a 1-400
#SBATCH --mem-per-cpu=32G
#SBATCH -o /scratch2/troshin/bmn_protons/qn_bmn/log/%A_%a.log
#SBATCH --exclude=ncx111,ncx112,ncx115,ncx117,ncx121,ncx127,ncx158,ncx159,ncx171,ncx181,ncx207,ncx214,ncx216,ncx222,ncx223,ncx224,ncx225,ncx227,ncx153,ncx123,ncx113,ncx138,ncx213,ncx142,ncx163,ncx168
date
hostname
list_dir=/scratch2/troshin/bmn_protons/qn_bmn/list_dir/
output_dir=/scratch2/troshin/bmn_protons/qn_bmn/out/runid_cent
id=$SLURM_ARRAY_TASK_ID
input_list=/scratch2/troshin/bmn_protons/run8.2025.04.0.list
split -l 73 -d -a 4 --additional-suffix=.txt $input_list $list_dir
#28394_142
#26807_134
#27520_138
#29209_146
file_list=$( ls $list_dir | head -n $id | tail -n 1 )

efficiency_dir=/scratch2/troshin/bmn_protons/qntools_macros/efficiency
calib_dir=/scratch2/troshin/bmn_protons/qntools_macros/centrality

mkdir -p $output_dir
cd $output_dir
mkdir $id
cd $id
qa=$output_dir/second_qa.root
source /cvmfs/nica.jinr.ru/sw/os/login.sh
source /cvmfs/bmn.jinr.ru/config/x86_64-centos7/cluster_config.sh
source /cvmfs/bmn.jinr.ru/bmnroot/25.09.0/x86_64-centos7/bmnroot_config.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch2/troshin/bmn_protons/QnTools/install/lib:/scratch2/troshin/bmn_protons/qntools_macros/build/

echo "Plain started"
time /scratch2/troshin/bmn_protons/qntools_macros/build/correct /scratch2/troshin/bmn_protons/qntools_macros/macro/run8_proton_correct_clean_runid.cc $list_dir/$file_list ${efficiency_dir}/efficiency.2022.01.25.root $calib_dir/CorrRunId.root ${qa}
echo "Recentered started"
time /scratch2/troshin/bmn_protons/qntools_macros/build/correct /scratch2/troshin/bmn_protons/qntools_macros/macro/run8_proton_correct_clean_runid.cc $list_dir/$file_list ${efficiency_dir}/efficiency.2022.01.25.root $calib_dir/CorrRunId.root ${qa}
echo "Rescaled started"
time /scratch2/troshin/bmn_protons/qntools_macros/build/correct /scratch2/troshin/bmn_protons/qntools_macros/macro/run8_proton_correct_clean_runid.cc $list_dir/$file_list ${efficiency_dir}/efficiency.2022.01.25.root $calib_dir/CorrRunId.root ${qa}
echo "Correlate started"
#time mv $output_dir/$id/qa.root $output_dir/qa_${id}.root
time /scratch2/troshin/bmn_protons/qntools_macros/build/correlate /scratch2/troshin/bmn_protons/qntools_macros/macro/run8_proton_correlate_clean.cc correction_out.root
time mv $output_dir/$id/corr.root $output_dir/corr_${id}.root

echo "The End."




