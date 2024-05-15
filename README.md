# iPPA

Inverse Primitive Path Analysis package for LAMMPS.

# Paper

The iPPA package was published in [Computer Physics Communications]{10.1016/j.cpc.2024.109209} for BibTeX entry see below.

# Compiling

Assuming your LAMMPS source code is located in ~user/lammps/, then copy the IPPA folder into  ~user/lammps/src/ 
Note you do not want to copy the content of the folder into source, but the whole folder.

Add IPPA folder as a package known by cmake edit ~user/lammps/cmake/CMakeLists.txt . There is a long statement
listing all packages which corresponds one-to-one with folders in the src directory, and add IPPA to the end as shown below.

```
set(STANDARD_PACKAGES
  ADIOS
  AMOEBA
:
  YAFF
  *IPPA*)
```

To add the IPPA package when you compile with cmake, add *-DPKG_IPPA=yes* to your cmake command.

# LAMMPS simulations with iPPA

To freeze a KG configuration use this LAMMPS snippet:
(For the complete code see the example folder)

```
variable tfreeze     equal 0.001        #temperature
variabl  tdampfreeze equal ${mass}/50   #friction
variabe  tstepfreeze equal 0.005        #timestep
variale  runfreeze   equal 1000         #number of steps to run

#setup force field;
pair_style      lj/cut ${cutoff}
pair_coeff      * * 1.0 1.0
pair_modify     shift yes

#modified FENE potential which halts if bond distance > 1.40
bond_style      fenehalt
bond_coeff      *   ${kspring} 1.5 0.0 1.0 1.40

angle_style     cosine
angle_coeff     * ${kappa}
special_bonds lj 1 1 1

#run Langevin simulation to freeze configuration
variable unity equal 1.0
thermo_style custom step temp epair ebond eangle c_bmin  c_bavg  c_bmax  c_pmin  v_unity
thermo 1

fix dynamics   mobile langevin ${tfreeze} ${tfreeze}   ${tdampfreeze}   325439524
timestep ${tstepfreeze}
run ${runfreeze}
unfix dynamics
```

To convert frozen KG system to its PPA mesh using iPPA force field transformation with topology violation checks:

```
variable atypebead  string "1*2"        #normal beads
variable atypetopo  string "3"          #dummy beads for topology check.
variable nwindow    equal  10           #number of bonds
variable u0         equal  200          #ucap 100
variable power      equal  log(${nwindow})/log(2)      # half of the time for final bond.
variable geometry   string ""           # ""=linear chain  or "circular" for circular chains.

variable ppatemp     equal  0.001       #temperature
variable tdampppa    equal  ${mass}/50  #friction
variable tstepppa    equal  0.005       #timestep
variable runippa     equal  10000       #number of steps to run

comm_modify cutoff 6

#define switching force field with topology monitoring
pair_style      hybrid  wca/ppa   window ${nwindow} alpha ${power}  lambda  1.0  ${geometry}  u0 ${u0}                  &
                        topo      window ${nwindow} alpha ${power}  lambda  1.0  ${geometry}  cutoff 2.0  bondtype *  debug
pair_coeff      ${atypebead}  ${atypebead}    wca/ppa  1.0 1.0
pair_coeff      ${atypebead}  ${atypetopo}    none
pair_coeff      ${atypetopo}  ${atypetopo}    topo
special_bonds lj 1 1 1

bond_style      poly06
bond_coeff      *  0     0  0  50 0 0  0   100
angle_style     none

#setup force field transformation from KG to PPA
variable ramp_down equal "ramp(1.00,0.00)"
thermo 1
thermo_style custom step temp epair ebond eangle c_bmin  c_bavg  c_bmax  c_pmin f_topo[2] v_ramp_down  

timestep ${tstepppa}

fix iPPA      all    adapt  1  pair wca/ppa lambda ${atypebead} ${atypebead} v_ramp_down  &
                               pair topo    lambda ${atypetopo} ${atypetopo} v_ramp_down
fix dynamics   mobile langevin ${ppatemp} ${ppatemp}   ${tdampppa}   325439524

run ${runippa}
unfix iPPA
unfix dynamics

```

To convert a mesh into a KG melt using the iPPA force field transformation, the code is nearly
identical to the code above, except lambda starts at 0.0 and the ramp goes 0 to 1:

```
variable runppa     equal  10000       #number of steps to run

comm_modify cutoff 6

#setup force field
pair_style      hybrid  wca/ppa   window ${nwindow} alpha ${power}  lambda  0.0  ${geometry}  u0 ${u0}                  &
                        topo      window ${nwindow} alpha ${power}  lambda  0.0  ${geometry}  cutoff 2.0  bondtype *  debug
pair_coeff      ${atypebead}  ${atypebead}    wca/ppa  1.0 1.0
pair_coeff      ${atypebead}  ${atypetopo}    none
pair_coeff      ${atypetopo}  ${atypetopo}    topo
special_bonds lj 1 1 1

bond_style      poly06
bond_coeff      *  0     0  0  50 0 0  0   100
angle_style     none

#run force field transformation from PPA to KG
variable ramp_up   equal "ramp(0.00,1.00)"
thermo     1
thermo_style custom step temp epair ebond eangle c_bmin  c_bavg  c_bmax  c_pmin   f_topo[2] v_ramp_up

fix iPPA      all    adapt  1  pair wca/ppa lambda ${atypebead} ${atypebead} v_ramp_up  &
                               pair topo    lambda ${atypetopo} ${atypetopo} v_ramp_up
fix dynamics  mobile langevin ${ppatemp} ${ppatemp}   ${tdampppa}   325439524

timestep ${tstepppa}
run     ${runppa}

unfix dynamics
unfix iPPA

```

Finally to heat up a the KG configuration resulting from the iPPA pushoff:

```
variable tdampwarmup    equal  ${mass}/200    #very high friction
variable tstepwarmup    equal 0.0001          #tiny time step
variable runwarmup      equal 1000            #short simulation

variable unity equal 1.0
thermo_style custom step temp epair ebond eangle c_bmin  c_bavg  c_bmax  c_pmin  v_unity
thermo 1

comm_modify cutoff 1

pair_style      lj/cut ${cutoff}
pair_coeff      * * 1.0 1.0
pair_modify     shift yes

bond_style      fenehalt
bond_coeff      * 30.0 1.5 0.0 1.0 1.49

angle_style     cosine
angle_coeff     * ${kappa}
special_bonds lj 1 1 1
#this is a matter of taste

fix dynamics   mobile langevin 1.0 1.0   ${tdampwarmup}   325439524
timestep ${tstepwarmup}
run     ${runwarmup}
unfix dynamics
```

Note that "special_bonds lj 1 1 1" is required during the force field transformations, since the WCA interaction between bonded beads is switched off
as part of the transformation of the pair energy. To keep a consistent definition of bond and pair energy, I also use this special bonds setting
for the KG simulations by disabling the WCA contribution to the FENE potential.

When pair topo is used it internally creates a fix topo that can be used for diagnostics via f_topo[1] to f_topo[4]. The first integers exported are
1) instanteous (intra+intermolecular) topology violations at the current time step, 2) accumulated (intra+intermolecular) topology violations so far,
3) instanteous intermolecular topology violations at the current time step, 4) accumulated intermolecular topology violations so far. Pair topo is
informed about the switching process, no topology checks is performed between bonds that within the window where interactions are switched off. Using
the difference between 2 and 4 gives an idea of how many topology violations are bonds cutting through other bonds on the same chain.

# LAMMPS simulations with KG with topology check

The topology checking code can be used e.g. during fast deformations to ensure that no bond
cuts through another during the deformation. (For the complete code see the example folder)

```
variable    atypebead  string "1*2"
variable    atypetopo  string "3"

:
comm_modify cutoff 4
variable cutoff  equal  2^(1/6)

pair_style      hybrid  lj/cut ${cutoff}                 &
                        topo      window 1    lambda  1.0   cutoff 2.0  bondtype *  debug
pair_coeff      ${atypebead}  ${atypebead}    lj/cut  1.0 1.0
pair_coeff      ${atypebead}  ${atypetopo}    none
pair_coeff      ${atypetopo}  ${atypetopo}    topo

bond_style      fene
bond_coeff      * 30.0 1.5 0.0 1.0

angle_style     cosine
angle_coeff     * ${kappa}
special_bonds lj 1 1 1

thermo_style custom step temp epair ebond eangle c_bmin  c_bavg  c_bmax  c_pmin  f_topo[2] f_topo[4]
:
```

The code is very similar to the force field switches above, except the lambda parameter is constant 1.0
signifying all beads interact as in the usual force field.  Only change is that instead of using
wca/ppa pair potential, we use standard lj/cut potential. For consistency of energy output, I have chosen
to switch off the WCA contribution to the FENE potential using special_bonds 1 1 1, but the standard
choice could also be used.


# SRP code

The iPPA topology violation counting code is based on a modified version of the [LAMMPS SRP package]{https://docs.lammps.org/pair_srp.html} of Tim Sirk et al.
This code works by inserting dummy beads at the center of each bond, and then using the neighbor lists of these dummy bonds to identify spatially neiboring bonds. The definition of pair_style topo automaticlly brings a fix topo to life. Fix topo augments particle information so they also carry their position at the previous time step, and various topology violation counters that are locally incremented when pair style topo identifies topology violations.

See T. W. Sirk, Y.R. Slizoberg, J.K. Brennan, M. Lisal, and J.W. Andzelm. (2012).
["An enhanced entangled polymer model for dissipative particle dynamics."]{https://doi.org/10.1063/1.3698476} The Journal of Chemical Physics, 136, 134903.


# iPPA paper

The reference for the iPPA package is

```
@article{Svaneborg2024iPPA,
  title={Inverse Primitive Path Analysis},
  author={Svaneborg, Carsten}
  journal={Computer Physics Communications},
  volume={XX},
  pages={XXXXX},
  year={2024},
  doi={10.1016/j.cpc.2024.109209},
  url={https://doi.org/10.1016/j.cpc.2024.109209}
}


```
