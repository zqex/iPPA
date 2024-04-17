# iPPA
Inverse Primitive Path Analysis package for LAMMPS

# Compiling

Assuming your LAMMPS source code is located in ~user/lammps/, then copy the IPPA folder into   ~user/lammps/src/ note do not copy the content of the folder into source, but the whole folder.

Add IPPA folder as a package known by cmake edit ~user/lammps/cmake/CMakeLists.txt .  There is a long statement listing all packages which corresponds one-to-one with folders in the src directory, and add IPPA to the end as shown below.

set(STANDARD_PACKAGES
  ADIOS
  AMOEBA
:
  VORONOI
  VTK
  YAFF
  IPPA)

To add the IPPA package when you compile with cmake, add -DPKG_IPPA=yes to your cmake command.

