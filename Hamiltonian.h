#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "main.h"

#define kp kroneckerProduct
#define Id(size) MatrixXd::Identity(size, size)
#define Id_d Matrix<double, d, d>::Identity()   // one-site identity matrix

typedef std::vector<MatrixD_t, Eigen::aligned_allocator<MatrixD_t>> vecMatD_t;

class Hamiltonian
{
    public:
        Hamiltonian();
        void setParams(const std::vector<double>& couplingConstants,
                       int targetQNumIn, int lSysIn);
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    private:
        Eigen::Matrix<double, nSiteTypes, nSiteTypes> intrabasisCouplings;
        // the row gives the type of site within the lattice basis
        // the column gives the length of the coupling in the stretched-out chain
        Eigen::Matrix<double, nSiteTypes, farthestInterbasisBond + 1>
            interbasisCouplings;
        // the row gives the type of site within the lattice basis
        // the column gives the number of lattice bases the bond crosses
        std::vector<int> oneSiteQNums;              // one-site quantum numbers
        int targetQNum,                              // targeted quantum number
            lSys;                                      // current system length
        vecMatD_t siteBasisH2;                 // site-basis coupling operators
        
        MatrixX_t blockAdjacentSiteJoin(int siteType, int jType,
                                        const std::vector<MatrixX_t>&
                                            offIRhoBasisH2,
                                        bool intrabasisBond) const,
                            // jType corresponds to the straightened-out chain,
                            // not the real-space coupling constants
                  lBlockrSiteJoin(int siteType, int jType,
                                  const std::vector<MatrixX_t>& offIRhoBasisH2,
                                  int compm, bool intrabasisBond) const,
                  lSiterBlockJoin(int siteType, int jType, int m,
                                  const std::vector<MatrixX_t>&
                                      compOffIRhoBasisH2,
                                  bool intrabasisBond) const,
                  siteSiteJoin(int siteType, int m, int compm) const,
                                           // joins the two free sites together
                  blockBlockJoin(int siteType, int l, int comp_l,
                                 const std::vector<std::vector<MatrixX_t>>&
                                     rhoBasisH2,
                                 const std::vector<std::vector<MatrixX_t>>&
                                     compRhoBasisH2) const;
                                               // joins the two blocks together
    
    friend class TheBlock;
};

#endif
