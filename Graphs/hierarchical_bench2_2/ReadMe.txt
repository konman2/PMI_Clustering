
This program is an implementation of the algorithm described in the paper "Directed, weighted and overlapping benchmark graphs for community detection algorithms", written by Andrea Lancichinetti and Santo Fortunato. In particular, this program is to produce binary networks with overlapping nodes and hierarchies.
Each feedback is very welcome. If you have found a bug or have problems, or want to give advises, please contact us:

andrea.lancichinetti@isi.it
fortunato@isi.it


Turin, 9 February 2011
---------------------------------------------------------------





-------------------- How to compile -----------------------------------
In order to compile, type:

make

-------------------- How to run the program ---------------------------


To run the program, type:

./benchmark [FLAG] [P]


[FLAG]		[P]

1.	-N              [number of nodes]
2.	-k              [average degree]
3.	-maxk           [maximum degree]
4.	-t1             [minus exponent for the degree sequence]


5.	-minc           [minimum for the micro community sizes]
6.	-maxc           [maximum for the micro community sizes]
7.	-on             [number of overlapping nodes (micro communities only)]
8.	-om             [number of memberships of the overlapping nodes (micro communities only)]

9.	-t2             [minus exponent for the community size distribution]

10.	-minC           [minimum for the macro community size]
11.	-maxC           [maximum for the macro community size]


12.	-mu1            [mixing parameter for the macro communities (see below)]
13.	-mu2            [mixing parameter for the micro communities (see below)]
----------------------



The flags from 5. to 8. concern the micro-communities. -t2 refers to micro and macro community distribution.

-mu1 is the fraction of links between nodes belonging to different macro-communities.
-mu2 is the fraction of links between nodes belonging to the same macro but not micro community.

More specifically, each node i has k(i) stubs: mu1 * k(i) of these stubs link to nodes of different macro-communities, and mu2 * k(i) link to nodes belonging to the same macro but not micro community.
So, (1 -mu2 -mu1) *k(i) stubs link to nodes belonging to the same micro and macro community.



-------------------- Examples ---------------------------

Example1:
./hbenchmark -N 10000 -k 20 -maxk 50 -mu2 0.3 -minc 20 -maxc 50 -minC 100 -maxC 1000 -mu1 0.1
Example2:
./hbenchmark -f flags.dat

Recapture [Hier-1] and [Hier-3] paper:
./hbenchmark -N 512 -k 48 -maxk 48 -minc 32 -maxc 32 -minC 128 -maxC 128 -mu1 0.33333333333 -mu2 0.33333333333

More realistic small graph for intuition:
hier7 (small 2*4*16)
./hbenchmark -N 128 -k 12 -maxk 12 -minc 16 -maxc 16 -minC 32 -maxC 32 -mu1 0.1 -mu2 0.1

hier8 (large 4*4*32)
./hbenchmark -N 512 -k 48 -maxk 48 -minc 32 -maxc 32 -minC 128 -maxC 128 -mu1 0.2 -mu2 0.2

hier10 (middle 4*4*16)
./hbenchmark -N 256 -k 24 -maxk 24 -minc 16 -maxc 16 -minC 64 -maxC 64 -mu1 0.2 -mu2 0.2

hier18,hier721 (small 3*3*12)
./hbenchmark -N 108 -k 12 -maxk 12 -minc 12 -maxc 12 -minC 36 -maxC 36 -mu1 0.1 -mu2 0.16

hier18_tailed (small 3*3*12)
./hbenchmark -N 108 -k 12 -maxk 18 -minc 12 -maxc 12 -minC 36 -maxC 36 -mu1 0.2 -mu2 0.2

hier18_varied (200 nodes)
./hbenchmark -N 200 -k 12 -maxk 18 -minc 4 -maxc 16 -minC 20 -maxC 100 -mu1 0.3 -mu2 0.3

hier_perfect, hier433 (4*4*10)
./hbenchmark -N 160 -k 10 -maxk 10 -minc 10 -maxc 10 -minC 40 -maxC 40 -mu1 0.3 -mu2 0.3

hier532d10 (4*4*10)
./hbenchmark -N 160 -k 10 -maxk 10 -minc 10 -maxc 10 -minC 40 -maxC 40 -mu1 0.2 -mu2 0.3

hier523d10 (4*4*10)
./hbenchmark -N 160 -k 10 -maxk 10 -minc 10 -maxc 10 -minC 40 -maxC 40 -mu1 0.3 -mu2 0.2

hier424d10 (4*4*10)
./hbenchmark -N 160 -k 10 -maxk 10 -minc 10 -maxc 10 -minC 40 -maxC 40 -mu1 0.4 -mu2 0.2

hier622 (3*3*12)
./hbenchmark -N 108 -k 12 -maxk 12 -minc 12 -maxc 12 -minC 36 -maxC 36 -mu1 0.2 -mu2 0.2

hier514 (3*3*12)
./hbenchmark -N 108 -k 12 -maxk 12 -minc 12 -maxc 12 -minC 36 -maxC 36 -mu1 0.4 -mu2 0.1

hier424 (3*3*12)
./hbenchmark -N 108 -k 12 -maxk 12 -minc 12 -maxc 12 -minC 36 -maxC 36 -mu1 0.4 -mu2 0.2

hier7310 (3*3*12)
./hbenchmark -N 108 -k 12 -maxk 12 -minc 12 -maxc 12 -minC 36 -maxC 36 -mu1 0.5 -mu2 0.15

hier1235 (3*3*12)
./hbenchmark -N 108 -k 12 -maxk 12 -minc 12 -maxc 12 -minC 36 -maxC 36 -mu1 0.25 -mu2 0.15

hier721l (3*3*24)
./hbenchmark -N 216 -k 24 -maxk 24 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.1 -mu2 0.16

hier721l2 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.1 -mu2 0.16

hier721l3 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.13 -mu2 0.17

hier721l4 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.18 -mu2 0.18

hier721l5 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.15 -mu2 0.15

hier721l6 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.1 -mu2 0.2

hier721l7 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.15 -mu2 0.24

hier721l8 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.2 -mu2 0.16

hier721l9 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.24 -mu2 0.16

hier721l10 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.3 -mu2 0.16

hier721l11 (3*3*24)
./hbenchmark -N 216 -k 12 -maxk 12 -minc 24 -maxc 24 -minC 72 -maxC 72 -mu1 0.3 -mu2 0.1
-------------------- Output ---------------------------

Please note that the community size distribution can be modified by the program to satisfy several constraints (a warning will be displayed).

The program will produce three files:

1) network.dat contains the list of edges (nodes are labelled from 1 to the number of nodes; the edges are ordered and repeated twice, i.e. source-target and target-source).

2) community_first_level.dat contains a list of the nodes and their membership for the micro-communities

3) community_second_level.dat is the same thing for the macro-communities




-------------------- Seed for the random number generator ---------------------------


-In the file time_seed.dat you can edit the seed which generates the random numbers. After reading, the program will increase this number by 1 (this is done to generate different networks running the program again and again). If the file is erased, it will be produced by the program again.
