      parameter (natmax=10000000,nelmax=12)
      common /lat/ natoms,ntypes,rv(6,natmax),itype(natmax),
     *   perlb(3),perub(3),perlen(3),xy,xz,yz,ilatseed,ntag(natmax),
     *   nntype(nelmax),amass(nelmax),ielement(nelmax),
     *   nhitcards,nhittag(natmax),nw_del,natoms0,numneigh(natmax),
     *   neigh(24,natmax)
