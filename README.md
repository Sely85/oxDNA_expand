oxDNA2dump
==========

Welcome to oxDNA_exapand code!  
This code converts oxDNA[1] dump
(obtained with LAMMPS[2]) configurations into a dump file that
contains explicitly base, stacking and backbone sites[3]. These file
can be read by OVITO visualization software[4].

[1] https://dna.physics.ox.ac.uk
[2] https://lammps.sandia.gov
[3] https://dna.physics.ox.ac.uk/index.php/Documentation 
[4] www.ovito.org

BUILD (Linux)
-------------
g++ -lm oxDNA_exapand.cpp -o oxDNA_exapand


USAGE
-----
./oxDNA_expand  <oxDNAdump>


EXAMPLE 
------- 

As an example, a file called "example.dump" has been attached. It
containes two double strand polyC placed at increasing distance
(101 frames). Using the following command

./oxDNA_example example.dump

will write a file "oxDNAconv_expand.dump" in which there are three
times the beads of the original file. For each "original" bead, the
following site have been created:
- backbone site (type 7), 
- stacking site (type 9),
- base site (adenine type 1, cytosine type 2, guanine type 3, thymine
  type 4).

A tentative shape (for OVITO rendering purpose) has been created for
each site: the base has an ellipsoidal shape, containing the stacking
as a small sphere, while the backbone is spherical.