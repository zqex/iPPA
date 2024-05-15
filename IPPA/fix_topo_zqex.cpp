// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing authors: Timothy Sirk (ARL), Pieter in't Veld (BASF)
------------------------------------------------------------------------- */

/*
    Modified from MISC/fix_srp.cpp
    
    Increased the per atom array from 2 to 3.
    
    SRP dummy beads introduced at the middle of each bond carry the two tags
    for identifying which end beads define the bond.
    
    Real beads carry the position of the previous timestep.

    topologyviolation counts the global number of topology violations this timestep,
    topologyviolationaccumulator counts the total topology violations so far.
*/

#include "fix_topo_zqex.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixTopo::FixTopo(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
 if (comm->me == 0)
    utils::logmesg(lmp,"FixTopo: Creating FixTopo");

  // settings
  nevery=1;
  peratom_freq = 1;
  time_integrate = 0;
  create_attribute = 0;
  comm_border = 4;

  // restart settings
  restart_global = 1;
  restart_peratom = 1;
  restart_pbc = 1;

  // ZQEX: Changed 2 to 3.
  // Dummy particles for bonds carry two ends beads.
  // Non-dummy particles carry their past x,y,z, coordinate.
  // per-atom array width 3
  peratom_flag = 1;
  size_peratom_cols = 3;

  // Counters for topology violations
  vector_flag=1;
  size_vector=4; 
  extvector=0;
  global_freq =1;
  debug=false;  

  // All topology violations  (intra+inter)
  topologyviolations=0;
  topologyviolationsaccumulator=0;  

  // intermolecular violations.
  topologyviolationsinter=0;
  topologyviolationsinteraccumulator=0;  

  // initial allocation of atom-based array
  // register with Atom class
  array = nullptr;
  flip = nullptr;
  FixTopo::grow_arrays(atom->nmax);

  // extends pack_exchange()
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART); // restart
  atom->add_callback(Atom::BORDER);

  // initialize to illegal values so we capture
  btype = -1;
  bptype = -1;

  pair_name = "topo";

  // zero
  for (int i = 0; i < atom->nmax; i++)
    for (int m = 0; m < 3; m++)
      array[i][m] = 0.0;

  for (int i = 0; i < atom->nmax; i++) flip[i]=false;
}

/* ---------------------------------------------------------------------- */

FixTopo::~FixTopo()
{
 if (comm->me == 0)
    utils::logmesg(lmp,"FixTopo: Destroying topo\n");

  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);
  atom->delete_callback(id,Atom::BORDER);
  memory->destroy(array);
  memory->destroy(flip);
}

/* ---------------------------------------------------------------------- */

int FixTopo::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_EXCHANGE;
  mask |= POST_RUN;
  mask |= POST_FORCE;

  mask |= MIN_PRE_FORCE;
  mask |= MIN_PRE_EXCHANGE;
  mask |= MIN_POST_FORCE;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTopo::init()
{
 if (comm->me == 0)
    utils::logmesg(lmp,"FixTopo: Reached init, atoms = {}\n",atom->natoms);

 // std::cout << "me=" << comm->me << "  reached FixTopo::init\n" << std::flush;
  if (!force->pair_match("^hybrid",0))
    error->all(FLERR,"FixTopo: Cannot use pair {} without pair_style hybrid", pair_name);

  if (modify->get_fix_by_style("^rigid").size() > 0)
    error->all(FLERR,"FixTopo: Pair {} is not compatible with rigid fixes.", pair_name);

  if ((bptype < 1) || (bptype > atom->ntypes))
    error->all(FLERR,"FixTopo: Illegal bond particle type");

  // this fix must come before any fix which migrates atoms in its pre_exchange()
  // because this fix's pre_exchange() creates per-atom data structure
  // that data must be current for atom migration to carry it along

  for (auto &ifix : modify->get_fix_list()) {
    if (ifix == this) break;
    if (ifix->pre_exchange_migrate)
      error->all(FLERR,"Fix {} comes after a fix which migrates atoms in pre_exchange", style);
  }

  // setup neigh exclusions for diff atom types
  // bond particles do not interact with other types
  // type bptype only interacts with itself

  for (int z = 1; z < atom->ntypes; z++) {
    if (z == bptype)
      continue;
    neighbor->modify_params(fmt::format("exclude type {} {}",z,bptype));
  }
  // std::cout << "me=" << comm->me << "  done FixTopo::init\n"  << std::flush;
}

/* ----------------------------------------------------------------------
   insert bond particles
------------------------------------------------------------------------- */

void FixTopo::setup_pre_force(int /*zz*/)
{
  double **x = atom->x;
  double **xold;
  tagint *tag = atom->tag;
  tagint *tagold;
  int *type = atom->type;
  int* dlist;
  AtomVec *avec = atom->avec;
  int **bondlist = neighbor->bondlist;

  int nlocal, nlocal_old;
  nlocal = nlocal_old = atom->nlocal;
  bigint nall = atom->nlocal + atom->nghost;
  int nbondlist = neighbor->nbondlist;
  int i,j,n;

  // make a copy of all coordinates and tags
  // that is consistent with the bond list as
  // atom->x will be affected by creating/deleting atoms.
  // also compile list of local atoms to be deleted.

  memory->create(xold,nall,3,"fix_srp:xold");
  memory->create(tagold,nall,"fix_srp:tagold");
  memory->create(dlist,nall,"fix_srp:dlist");

  for (i = 0; i < nall; i++) {
    xold[i][0] = x[i][0];
    xold[i][1] = x[i][1];
    xold[i][2] = x[i][2];
    tagold[i]=tag[i];
    dlist[i] = (type[i] == bptype) ? 1 : 0;
    
    if (type[i] == bptype)
        for (n = 0; n < 2; n++)  array[i][n] = 0.0;
      else
        {
          array[i][0]=x[i][0];
          array[i][1]=x[i][1];
          array[i][2]=x[i][2];
        }
  }

  // delete local atoms flagged in dlist
  i = 0;
  int ndel = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      avec->copy(nlocal-1,i,1);
      dlist[i] = dlist[nlocal-1];
      nlocal--;
      ndel++;
    } else i++;
  }

  atom->nlocal = nlocal;
  memory->destroy(dlist);

  int nadd = 0;
  double rsqold = 0.0;
  double delx, dely, delz, rmax, rsq, rsqmax;
  double xone[3];

  for (n = 0; n < nbondlist; n++) {

    // consider only the user defined bond type
    // btype of zero considers all bonds
    if (btype > 0 && bondlist[n][2] != btype)
      continue;

    i = bondlist[n][0];
    j = bondlist[n][1];

    // position of bond i
    xone[0] = (xold[i][0] + xold[j][0])*0.5;
    xone[1] = (xold[i][1] + xold[j][1])*0.5;
    xone[2] = (xold[i][2] + xold[j][2])*0.5;

    // record longest bond
    // this used to set ghost cutoff
    delx = xold[j][0] - xold[i][0];
    dely = xold[j][1] - xold[i][1];
    delz = xold[j][2] - xold[i][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > rsqold) rsqold = rsq;

    // make one particle for each bond
    // i is local
    // if newton bond, always make particle
    // if j is local, always make particle
    // if j is ghost, decide from tag

    if ((force->newton_bond) || (j < nlocal_old) || (tagold[i] > tagold[j])) {
      atom->natoms++;
      avec->create_atom(bptype,xone);
      // pack tag i/j into buffer for comm
      array[atom->nlocal-1][0] = static_cast<double>(tagold[i]);
      array[atom->nlocal-1][1] = static_cast<double>(tagold[j]);
      nadd++;
    }
  }

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  // free temporary storage
  memory->destroy(xold);
  memory->destroy(tagold);

  int nadd_all = 0, ndel_all = 0;
  MPI_Allreduce(&ndel,&ndel_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nadd,&nadd_all,1,MPI_INT,MPI_SUM,world);
  if (comm->me == 0)
    error->message(FLERR,"FixTopo: Removed/inserted {}/{} bond particles.", ndel_all,nadd_all);

  // check ghost comm distances
  // warn and change if shorter from estimate
  // ghost atoms must be present for bonds on edge of neighbor cutoff
  // extend cutghost slightly more than half of the longest bond
  MPI_Allreduce(&rsqold,&rsqmax,1,MPI_DOUBLE,MPI_MAX,world);
  rmax = sqrt(rsqmax);
  double cutneighmax_srp = neighbor->cutneighmax + 0.51*rmax;

  double length0,length1,length2;
  if (domain->triclinic) {
    double *h_inv = domain->h_inv;
    length0 = sqrt(h_inv[0]*h_inv[0] + h_inv[5]*h_inv[5] + h_inv[4]*h_inv[4]);
    length1 = sqrt(h_inv[1]*h_inv[1] + h_inv[3]*h_inv[3]);
    length2 = h_inv[2];
  } else length0 = length1 = length2 = 1.0;

  // find smallest cutghost.
  // comm->cutghost is stored in fractional coordinates for triclinic
  double cutghostmin = comm->cutghost[0]/length0;
  if (cutghostmin > comm->cutghost[1]/length1)
    cutghostmin = comm->cutghost[1]/length1;
  if (cutghostmin > comm->cutghost[2]/length2)
    cutghostmin = comm->cutghost[2]/length2;

  // stop if cutghost is insufficient
  if (cutneighmax_srp > cutghostmin)
    error->all(FLERR,"FixTopo: Communication cutoff too small for fix {}. "
               "Need {:.8}, current {:.8}", style, cutneighmax_srp, cutghostmin);

  // assign tags for new atoms, update map
  atom->tag_extend();
  if (atom->map_style != Atom::MAP_NONE) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // put new particles in the box before exchange
  // move owned to new procs
  // get ghosts
  // build neigh lists again

  // if triclinic, lambda coords needed for pbc, exchange, borders
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  // back to box coords
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();
  neighbor->build(1);
  neighbor->ncalls = 0;

  // new atom counts

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  // zero all forces

  for (i = 0; i < nall; i++)
    atom->f[i][0] = atom->f[i][1] = atom->f[i][2] = 0.0;

  // do not include bond particles in thermo output
  // remove them from all groups. set their velocity to zero.

  for (i=0; i< nlocal; i++)
    if (atom->type[i] == bptype) {
      atom->mask[i] = 0;
      atom->v[i][0] = atom->v[i][1] = atom->v[i][2] = 0.0;
    }
   // std::cout << "me=" << comm->me << "  end FixTopo::setup_pre-force\n" << std::flush;

  // Reset all flips.
  for (i = 0; i < nall; i++) flip[i]=false;
}

/* ----------------------------------------------------------------------
   set position of bond particles
------------------------------------------------------------------------- */

void FixTopo::pre_exchange()
/*
  pre_exchange  updates the coordinates of the dummy bond-beads. Hence we do not need to worry about updating array,
  since this ís only used for storage of real beads.
*/

{
/*
  if (comm->me == 0)
    utils::logmesg(lmp,"Reached pre_exchange.\n");
*/
  // update ghosts
  comm->forward_comm();

  // reassign bond particle coordinates to midpoint of bonds
  // only need to do this before neigh rebuild
  double **x=atom->x;
  int i,j;
  int nlocal = atom->nlocal;

  for (int ii = 0; ii < nlocal; ii++) {
    if (atom->type[ii] != bptype) continue;
    
    i = atom->map(static_cast<tagint>(array[ii][0]));
    if (i < 0) error->all(FLERR,"FixTopo: failed to map atom");
    i = domain->closest_image(ii,i);

    j = atom->map(static_cast<tagint>(array[ii][1]));
    if (j < 0) error->all(FLERR,"FixTopo: failed to map atom");
    j = domain->closest_image(ii,j);

    // position of bond particle ii
    x[ii][0] = (x[i][0] + x[j][0])*0.5;
    x[ii][1] = (x[i][1] + x[j][1])*0.5;
    x[ii][2] = (x[i][2] + x[j][2])*0.5;

/*
    utils::logmesg(lmp,fmt::format("FixTopo: pre_exchange updating  ii={} tag={} at {} {} {}   {}\n",ii,atom->tag[ii],x[ii][0],x[ii][1],x[ii][2], atom->type[ii] ))/ 
    utils::logmesg(lmp,fmt::format("                                i={} tag={} at {} {} {}   {}\n",i,atom->tag[i],x[i][0],x[i][1],x[i][2], atom->type[i] ));
    utils::logmesg(lmp,fmt::format("                                j={} tag={} at {} {} {}   {}\n",j,atom->tag[j],x[j][0],x[j][1],x[j][2], atom->type[j] ));
*/        

// Occationally I observe that 0,0,0 positions appear in the list. Increasing the skin size, or reducing the dmax removes this bug,
// which seems to be related to neighbor communication.

    if (update->ntimestep>0 && x[i][0]==0 && x[i][1]==0 && x[i][2]==0 )
       utils::logmesg(lmp,fmt::format("fix_topo Warning  n={} atom i={} tag={} typ={} pos= {} {} {} pos=zero can indicate too large steps has been applied e.g. during minimization reduce dmax\n",update->ntimestep, i,atom->tag[i],atom->tag[i],x[i][0],x[i][1],x[i][2] ));

    if (update->ntimestep>0 && x[j][0]==0 && x[j][1]==0 && x[j][2]==0 )
       utils::logmesg(lmp,fmt::format("fix_topo Warning  n={} atom j={} tag={} typ={} pos= {} {} {} pos=zero can indicate too large steps has been applied e.g. during minimization reduce dmax\n",update->ntimestep, j,atom->tag[j],atom->tag[j],x[j][0],x[j][1],x[j][2] ));

  }     
   // std::cout << "me=" << comm->me << "  end FixTopo::pre_exchange\n" << std::flush;
}

void FixTopo::min_pre_exchange()
{
  pre_exchange();
}

void FixTopo::post_force(int /* what */)
// Added by zqex: Store current coordinate in array for reach real particle for detecting topology violations.
// This occurs just before integrating the forces to move the particles, hence in the next time step
// array contains the previous coordinates.
{
/*
  if (comm->me == 0)
    utils::logmesg(lmp,"FixTopo:: Storing previous coordinates\n");
*/
  // WE need to retrieve flip state for ghosts and OR them with the local flip state.
  if (doFlip) comm->reverse_comm(this,1);

  // std::cout << "Now post_force was called! with " << atom->nmax << "\n" << std::flush;
  double **x=atom->x;
  for (int ii = 0; ii < atom->nlocal; ii++) {
    if (atom->type[ii] == bptype) continue;

/*
    if (atom->type[ii] == 2)
      {
        utils::logmesg(lmp,fmt::format("FixTopo: post_force  {} {} at {} {} {}\n",ii,atom->tag[ii],x[ii][0],x[ii][1],x[ii][2]));
      }
*/

    if (!doFlip)
         {
           array[ii][0] = x[ii][0];
           array[ii][1] = x[ii][1];
           array[ii][2] = x[ii][2];
         }
     else if (!flip[ii])
         {
           array[ii][0] = x[ii][0];
           array[ii][1] = x[ii][1];
           array[ii][2] = x[ii][2];
         }
      else
         {  // doFlip && flip[ii]
           if (debug)
              {
                std::string fnam = fmt::format("DEBUG_topo.me{}.txt",comm->me);
                std::ofstream fo(fnam.c_str(), std::ios_base::app);
                fo << fmt::format("{} FLIPPED {}\n", update->ntimestep, atom->tag[ii]);
                fo.close();
              }
 
           x[ii][0] = array[ii][0];
           x[ii][1] = array[ii][1];
           x[ii][2] = array[ii][2];
                        
           atom->v[ii][0]*=-1.0;
           atom->v[ii][1]*=-1.0;
           atom->v[ii][2]*=-1.0;
         }

  }

  // Reset flip data.
  bigint nall = atom->nlocal + atom->nghost;
  if (doFlip)
      for (int ii = 0; ii < nall; ii++) flip[ii]=false;

   // std::cout << "me=" << comm->me << "  end FixTopo::pre_exchange\n" << std::flush;
}

void FixTopo::min_post_force( int z )
{
   post_force(z);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixTopo::memory_usage()
{
  double bytes = (double)atom->nmax*4 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixTopo::grow_arrays(int nmax)
{
  memory->grow(array,nmax,3,"fix_srp:array");
  memory->grow(flip,nmax,"Fix_topo:flip");
  array_atom = array;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
   called when move to new proc
------------------------------------------------------------------------- */

void FixTopo::copy_arrays(int i, int j, int /*delflag*/)
{
  for (int m = 0; m < 3; m++)
    array[j][m] = array[i][m];

  flip[j]=flip[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values
   called when atom is created
------------------------------------------------------------------------- */

void FixTopo::set_arrays(int i)
{
  array[i][0] = -1;
  array[i][1] = -1;
  array[i][2] = -1;

  flip[i]=false;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixTopo::pack_exchange(int i, double *buf)
{
  int m = 0;  
  buf[m++] = array[i][0];
  buf[m++] = array[i][1];
  buf[m++] = array[i][2];
  buf[m++] = flip[i];

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixTopo::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  array[nlocal][0] = buf[m++];
  array[nlocal][1] = buf[m++];
  array[nlocal][2] = buf[m++];
  flip[nlocal]=buf[m++];
  return m;
}
/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixTopo::pack_border(int n, int *list, double *buf)
{
  // pack buf for border com
  int i,j;
  int m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = array[j][0];
      buf[m++] = array[j][1];
      buf[m++] = array[j][2];
      buf[m++] = flip[j];
    }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixTopo::unpack_border(int n, int first, double *buf)
{
  // unpack buf into array
  int i,last;
  int m = 0;
  last = first + n;

  for (i = first; i < last; i++) {
    array[i][0] = buf[m++];
    array[i][1] = buf[m++];
    array[i][2] = buf[m++];
    flip[i]=buf[m++];
  }
  return m;
}

// Transfer flip information from ghosts back to the core that owns the particle
int FixTopo::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

    for (i = first; i < last; i++)
      buf[m++] = ubuf(flip[i]).d;
  return m;
}

/* ---------------------------------------------------------------------- */

// Receive flip information from ghosts.
// IF any ghosts are flipped, then their owner atom becomes flipped.
void FixTopo::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
      j = list[i];
      flip[j] = flip[j] || (bool)ubuf(buf[m++]).i;
    }
}


/* ----------------------------------------------------------------------
   remove particles after run
------------------------------------------------------------------------- */

void FixTopo::post_run()
{
  // all bond particles are removed after each run
  // useful for write_data and write_restart commands
  // since those commands occur between runs

  bigint natoms_previous = atom->natoms;
  int nlocal = atom->nlocal;
  int* dlist;
  memory->create(dlist,nlocal,"fix_srp:dlist");

  for (int i = 0; i < nlocal; i++) {
    if (atom->type[i] == bptype)
      dlist[i] = 1;
    else
      dlist[i] = 0;
  }

  // delete local atoms flagged in dlist
  // reset nlocal

  AtomVec *avec = atom->avec;

  int i = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      avec->copy(nlocal-1,i,1);
      dlist[i] = dlist[nlocal-1];
      nlocal--;
    } else i++;
  }

  atom->nlocal = nlocal;
  memory->destroy(dlist);

  // reset atom->natoms
  // reset atom->map if it exists
  // set nghost to 0 so old ghosts won't be mapped

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (atom->map_style != Atom::MAP_NONE) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // print before and after atom count

  bigint ndelete = natoms_previous - atom->natoms;

  if (comm->me == 0)
    utils::logmesg(lmp,"Deleted {} atoms, new total = {}\n",ndelete,atom->natoms);

  // verlet calls box_too_small_check() in post_run
  // this check maps all bond partners
  // therefore need ghosts

  // need to convert to lambda coords before apply pbc
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->setup();
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  // change back to box coordinates
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTopo::pack_restart(int i, double *buf)
{
  int m = 1;
  // pack buf[0] this way because other fixes unpack it
  buf[m++] = array[i][0];
  buf[m++] = array[i][1];
  buf[m++] = array[i][2];
  buf[m++] = flip[i];
  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixTopo::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) {
    m += static_cast<int> (extra[nlocal][m]);
  }

  m++;
  array[nlocal][0] = extra[nlocal][m++];
  array[nlocal][1] = extra[nlocal][m++];
  array[nlocal][2] = extra[nlocal][m++];
  flip[nlocal] = extra[nlocal][m++];
}
/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTopo::maxsize_restart()
{
  return 5;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTopo::size_restart(int /*nlocal*/)
{
  return 5;
}

/* ----------------------------------------------------------------------
   pack global state of Fix
------------------------------------------------------------------------- */

void FixTopo::write_restart(FILE *fp)
{
  int n = 0;
  double list[7];
  list[n++] = comm->cutghostuser;
  list[n++] = btype;
  list[n++] = bptype;
  list[n++] = topologyviolations;
  list[n++] = topologyviolationsaccumulator;
  list[n++] = topologyviolationsinter;
  list[n++] = topologyviolationsinteraccumulator;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }

 if (comm->me == 0 && debug)
    utils::logmesg(lmp,"Fix restart writing data to restart file\n");
}

/* ----------------------------------------------------------------------
   use info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTopo::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;

  comm->cutghostuser = static_cast<double> (list[n++]);
  btype = static_cast<int> (list[n++]);
  bptype = static_cast<int> (list[n++]);
  topologyviolations = static_cast<int> (list[n++]);
  topologyviolationsaccumulator = static_cast<int> (list[n++]);
  topologyviolationsinter = static_cast<int> (list[n++]);
  topologyviolationsinteraccumulator = static_cast<int> (list[n++]);

 if (comm->me == 0 && debug)
    utils::logmesg(lmp,"Fix restart reading data from restart file\n");
}

/* ----------------------------------------------------------------------
   interface with pair class
   pair srp sets the bond type in this fix
------------------------------------------------------------------------- */

int FixTopo::modify_param(int /*narg*/, char **arg)
{
  if (strcmp(arg[0],"btype") == 0) {
    btype = utils::inumeric(FLERR,arg[1],false,lmp);
    return 2;
  }
  if (strcmp(arg[0],"bptype") == 0) {
    bptype = utils::inumeric(FLERR,arg[1],false,lmp);
    return 2;
  }
  return 0;
}

double FixTopo::compute_vector(int n)
{
   // std::cout << "me=" << comm->me << "  start FixTopo::compute_Vector\n" << std::flush;

  int all=0;
  if (n==0)       MPI_Allreduce(&topologyviolations               ,&all,1,MPI_INT,MPI_SUM,world);
  else if (n==1)  MPI_Allreduce(&topologyviolationsaccumulator    ,&all,1,MPI_INT,MPI_SUM,world);    
  else if (n==2)  MPI_Allreduce(&topologyviolationsinter            ,&all,1,MPI_INT,MPI_SUM,world);    
  else if (n==3)  MPI_Allreduce(&topologyviolationsinteraccumulator ,&all,1,MPI_INT,MPI_SUM,world);    

  return (double)all;
}

