Energy Efficiency in MIMO Underlay and Overlay Device-to-Device Communications and Cognitive Radio Systems
==================

This code package is related to the following scientific article:

Alessio Zappone, Bho Matthiesen, and Eduard Jorswieck, "[Energy Efficiency in MIMO Underlay and Overlay Device-to-Device Communications and Cognitive Radio Systems](https://doi.org/10.1109/TSP.2016.2626249)," *IEEE Transactions on Signal Processing*, vol. 65, no. 4, pp. 1026-1041 Feb. 2017.


## Abstract of Article

This paper addresses the problem of resource allocation for systems in which a primary and a secondary link share the available spectrum by an underlay or overlay approach. After observing that such a scenario models both cognitive radio and D2D communications, we formulate the problem as the maximization of the secondary energy efficiency subject to a minimum rate requirement for the primary user. This leads to challenging nonconvex, fractional problems. In the underlay scenario, we obtain the global solution by means of a suitable reformulation. In the overlay scenario, two algorithms are proposed. The first one yields a resource allocation fulfilling the first-order optimality conditions of the resource allocation problem, by solving a sequence of easier fractional programs. The second one enjoys a weaker optimality claim, but an even lower computational complexity. Numerical results demonstrate the merits of the proposed algorithms both in terms of energy-efficient performance and complexity, also showing that the two proposed algorithms for the overlay scenario perform very similarly, despite the different complexity.

## Remarks on the Contents of this Code Package

We have decided to relase this code although it is in no shape to do so. This is due to several requests we have received since the first publication of this article. It is provided as is, without any documentation or support. 

The overlay code should be located in the root folder and `mpi/`. If I recall correctly, we have used the [bcMPI](www.bluecollarcomputing.org/applications/bcMPI/index.shtml) toolbox and [Open MPI](https://www.open-mpi.org) to run it. The start script is `mpi/manager.m`. The underlay code should be easier to run and is located in `underlay/`. In addition, the [CVX](http://cvxr.com/cvx) toolbox is required, possibly with a professional license that is free of charge for [academic users](http://cvxr.com/cvx/academic/). Good luck!


## Acknowledgements

The work of A. Zappone was supported
by the German Research Foundation, Deutsche Forschungsgemeinschaft (DFG), under grant ZA 747/1-3.

The work of E. Jorswieck and B. Matthiesen was supported in part by
the DFG in the [Collaborative Research Center 912 "Highly Adaptive Energy-Efficient Computing"](https://tu-dresden.de/ing/forschung/sfb912).


## License and Referencing

This program is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

