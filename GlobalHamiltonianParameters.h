#ifndef GHP_H
#define GHP_H

const int d = 2,                              // size of one-site Hilbert space
          nCouplingConstants = 3,               // number of coupling constants
          nCouplingOperators = 3,               // number of coupling operators
          nIndepCouplingOperators = 2,
              // number of that are independent - must be <= nCouplingOperators
          nSiteTypes = 3,
              // number of distinct kinds of sites (e.g. size of Bravais basis)
          farthestInterbasisBond = 2,
                          // number of lattice bases connected by farthest bond
          farthestNeighborCoupling = nSiteTypes * farthestInterbasisBond;
        // max coupling size after the system has been straighted out to a line

// only one of these next three lines should be uncommented, depending on
// nSiteTypes:
// #define basisSize2
#define basisSize3
// #define basisSizeGreaterThan3

// only one of these next two lines should be uncommented, depending on whether
// the Hamiltonian has real or complex elements:
#define realHamiltonian
// #define complexHamiltonian

// only one of these next two lines should be uncommented, depending on whether
// the observable matrices have real or complex elements:
#define realObservables
// #define complexObservables

#endif
