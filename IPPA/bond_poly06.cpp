/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "bond_poly06.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondPoly06::BondPoly06(LAMMPS *lmp) : Bond(lmp) 
{
}

/* ---------------------------------------------------------------------- */

BondPoly06::~BondPoly06()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(r0);
    memory->destroy(k1);
    memory->destroy(k2);
    memory->destroy(k3);
    memory->destroy(k4);
    memory->destroy(k5);
    memory->destroy(k6);
  }
}

/* ---------------------------------------------------------------------- */

void BondPoly06::compute(int eflag, int vflag)
{
  int i1, i2, n, type;
  double delx, dely, delz, ebond, fbond;
  double rsq, r, dr, dr2, dr3, dr4, dr5, dr6,  de_bond;

  ebond = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);
    dr = r - r0[type];
    dr2 = dr*dr;
    dr3 = dr2*dr;
    dr4 = dr2*dr2;
    dr5 = dr2*dr3;
    dr6 = dr3*dr3;

    // force & energy

    if (r > 0.0)
      fbond = -(k1[type]+2.0*k2[type]*dr + 3.0*k3[type]*dr2 + 4.0*k4[type]*dr3 + 5.0*k5[type]*dr4 + 6.0*k6[type]*dr5)/r;
    else
      fbond = 0.0;

    if (eflag) ebond = k0[type]+k1[type]*dr+k2[type]*dr2 + k3[type]*dr3 + k4[type]*dr4 + k5[type]*dr5 + k6[type]*dr6;

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx * fbond;
      f[i1][1] += dely * fbond;
      f[i1][2] += delz * fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx * fbond;
      f[i2][1] -= dely * fbond;
      f[i2][2] -= delz * fbond;
    }

    if (evflag) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondPoly06::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(k0, np1, "bond:k0");
  memory->create(k1, np1, "bond:k1");
  memory->create(k2, np1, "bond:k2");
  memory->create(k3, np1, "bond:k3");
  memory->create(k4, np1, "bond:k4");
  memory->create(k5, np1, "bond:k4");
  memory->create(k6, np1, "bond:k5");
  memory->create(r0, np1, "bond:r0");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondPoly06::coeff(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double r0_one = utils::numeric(FLERR,arg[1],false,lmp);
  double k0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double k1_one = utils::numeric(FLERR,arg[3],false,lmp);
  double k2_one = utils::numeric(FLERR,arg[4],false,lmp);
  double k3_one = utils::numeric(FLERR,arg[5],false,lmp);
  double k4_one = utils::numeric(FLERR,arg[6],false,lmp);
  double k5_one = utils::numeric(FLERR,arg[7],false,lmp);
  double k6_one = utils::numeric(FLERR,arg[8],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k0[i] = k0_one;
    k1[i] = k1_one;
    k2[i] = k2_one;
    k3[i] = k3_one;
    k4[i] = k4_one;
    k5[i] = k5_one;
    k6[i] = k6_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondPoly06::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondPoly06::write_restart(FILE *fp)
{
  fwrite(&r0[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&k0[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&k1[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&k2[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&k3[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&k4[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&k5[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&k6[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondPoly06::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &r0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &k0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &k1[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &k2[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &k3[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &k4[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &k5[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &k6[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&r0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k1[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k2[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k3[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k4[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k5[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k6[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondPoly06::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
 fprintf(fp,"%d %g %g %g %g %g %g %g %g\n", i, r0[i], k0[i], k1[i], k2[i], k3[i], k4[i], k5[i], k6[i] );
}

/* ---------------------------------------------------------------------- */

double BondPoly06::single(int type, double rsq, int i, int j, double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double dr2 = dr*dr;
  double dr3 = dr2*dr;
  double dr4 = dr2*dr2;
  double dr5 = dr2*dr3;
  double dr6 = dr3*dr3;

  fforce=0;
  if (r > 0.0) fforce = -(k1[type]+2.0*k2[type]*dr + 3.0*k3[type]*dr2 + 4.0*k4[type]*dr3+ 5.0*k5[type]*dr4+ 6.0*k6[type]*dr5)/r;
  return (k0[type]+k1[type]*dr+k2[type]*dr2 + k3[type]*dr3 + k4[type]*dr4 + k5[type]*dr5 + k6[type]*dr6);
}
