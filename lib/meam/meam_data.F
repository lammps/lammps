
      module meam_data

      integer, parameter :: maxelt = 5

c cutforce = force cutoff
c cutforcesq = force cutoff squared

      real*8 cutforce,cutforcesq

c Ec_meam = cohesive energy
c re_meam = nearest-neighbor distance
c Omega_meam = atomic volume
c B_meam = bulk modulus
c Z_meam = number of first neighbors for reference structure
c ielt_meam = atomic number of element
c A_meam = adjustable parameter
c alpha_meam = sqrt(9*Omega*B/Ec)
c rho0_meam = density scaling parameter
c delta_meam = heat of formation for alloys
c beta[0-3]_meam = electron density constants 
c t[0-3]_meam = coefficients on densities in Gamma computation
c rho_ref_meam = background density for reference structure
c ibar_meam(i) = selection parameter for Gamma function for elt i,
c lattce_meam(i,j) = lattce configuration for elt i or alloy (i,j)
c neltypes = maximum number of element type defined
c eltind = index number of pair (similar to Voigt notation; ij = ji)
c phir = pair potential function array
c phirar[1-6] = spline coeffs
c attrac_meam = attraction parameter in Rose energy
c repuls_meam = repulsion parameter in Rose energy
c nn2_meam = 1 if second nearest neighbors are to be computed, else 0
c zbl_meam = 1 if zbl potential for small r to be use, else 0
c emb_lin_neg = 1 if linear embedding function for rhob to be used, else 0
c bkgd_dyn = 1 if reference densities follows Dynamo, else 0
c Cmin_meam, Cmax_meam = min and max values in screening cutoff
c rc_meam = cutoff distance for meam
c delr_meam = cutoff region for meam
c ebound_meam = factor giving maximum boundary of sceen fcn ellipse
c augt1 = flag for whether t1 coefficient should be augmented
c ialloy = flag for newer alloy formulation (as in dynamo code)
c mix_ref_t = flag to recover "old" way of computing t in reference config
c erose_form = selection parameter for form of E_rose function
c gsmooth_factor = factor determining length of G smoothing region
c vind[23]D = Voight notation index maps for 2 and 3D
c v2D,v3D = array of factors to apply for Voight notation

c nr,dr = pair function discretization parameters
c nrar,rdrar = spline coeff array parameters

      real*8 Ec_meam(maxelt,maxelt),re_meam(maxelt,maxelt)
      real*8 Omega_meam(maxelt),Z_meam(maxelt)
      real*8 A_meam(maxelt),alpha_meam(maxelt,maxelt),rho0_meam(maxelt)
      real*8 delta_meam(maxelt,maxelt)
      real*8 beta0_meam(maxelt),beta1_meam(maxelt)
      real*8 beta2_meam(maxelt),beta3_meam(maxelt)
      real*8 t0_meam(maxelt),t1_meam(maxelt)
      real*8 t2_meam(maxelt),t3_meam(maxelt)
      real*8 rho_ref_meam(maxelt)
      integer ibar_meam(maxelt),ielt_meam(maxelt)
      character*3 lattce_meam(maxelt,maxelt)
      integer nn2_meam(maxelt,maxelt)
      integer zbl_meam(maxelt,maxelt)
      integer eltind(maxelt,maxelt)
      integer neltypes

      real*8, allocatable :: phir(:,:)

      real*8, allocatable :: phirar(:,:),phirar1(:,:),phirar2(:,:),
     $     phirar3(:,:),phirar4(:,:),phirar5(:,:),phirar6(:,:)

      real*8 attrac_meam(maxelt,maxelt),repuls_meam(maxelt,maxelt)
 
      real*8 Cmin_meam(maxelt,maxelt,maxelt)
      real*8 Cmax_meam(maxelt,maxelt,maxelt)
      real*8 rc_meam,delr_meam,ebound_meam(maxelt,maxelt)
      integer augt1, ialloy, mix_ref_t, erose_form
      integer emb_lin_neg, bkgd_dyn
      real*8  gsmooth_factor
      integer vind2D(3,3),vind3D(3,3,3)
      integer v2D(6),v3D(10)

      integer nr,nrar
      real*8 dr,rdrar

      end module
