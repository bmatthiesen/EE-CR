#!/bin/sh

# Copyright (C) 2014-2016 Bho Matthiesen, Alessio Zappone
# 
# This program is used in the article:
# 
# Alessio Zappone, Bho Matthiesen, and Eduard Jorswieck, "Energy Efficiency in
# MIMO Underlay and Overlay Device-to-Device Communications and Cognitive Radio
# Systems," IEEE Transactions on Signal Processing, vol. 65, no. 4, pp.
# 1026-1041 Feb. 2017, https://doi.org/10.1109/TSP.2016.2626249
# 
# 
# License:
# This program is licensed under the GPLv2 license. If you in any way use this
# code for research that results in publications, please cite our original
# article listed above.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


#SBATCH --array=1-27
#SBATCH --time=47:59:59
#SBATCH --mem-per-cpu 7875
#SBATCH --cpus-per-task=1
#SBATCH --mail-user 'bho.matthiesen@tu-dresden.de'
#SBATCH --mail-type ALL
#SBATCH -J "underlay_EE"

export UNDERLAY_HPC_SAVEDIR="/scratch/p_mimo/bho/underlay_OOBIF-2"
export UNDERLAY_WP_SUFFIX="_PIF1=-10dBW_d1=1000m_PIF2=-30dBW_d2=600m"
#export UNDERLAY_WP_SUFFIX="_PIF1=-10dBW_d1=1000m_PIF2=-20dBW_d2=600m"
#export UNDERLAY_WP_SUFFIX="_PIF1=0dBW_d1=1000m_PIF2=-30dBW_d2=600m"
#export UNDERLAY_WP_SUFFIX="_PIF1=0dBW_d1=1000m_PIF2=-20dBW_d2=600m"

module load matlab
srun matlab -nodesktop -nodisplay -nosplash -r underlay_EE
