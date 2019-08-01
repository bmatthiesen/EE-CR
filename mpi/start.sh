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

echo "Doing git push in 5 seconds"
sleep 5
git push

echo "\n\ngit pull on host1"
ssh host1 'cd sources && git pull'

echo "\n\ngit pull on host2"
ssh host2 'cd data/work/src && git pull'

echo "\n\nstarting MPI on host2"
ssh host2 'cd MPI && screen -L -m bash -c "source ompi1.2.7/env.sh && mpirun -n 16 --hostfile host_file ./mpitb_matlab mpi/manager.m"'
