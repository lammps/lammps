#include "meam_kokkos.h"

template<class DeviceType>
void MEAMKokkos<DeviceType>::meam_setup_done(void)
{
  MemoryKokkos *memoryKK = (MemoryKokkos *)memory;

  memoryKK->destroy_kokkos(k_phirar6,phirar6);
  memoryKK->destroy_kokkos(k_phirar5,phirar5);
  memoryKK->destroy_kokkos(k_phirar4,phirar4);
  memoryKK->destroy_kokkos(k_phirar3,phirar3);
  memoryKK->destroy_kokkos(k_phirar2,phirar2);
  memoryKK->destroy_kokkos(k_phirar1,phirar1);
  memoryKK->destroy_kokkos(k_phirar,phirar);
  memoryKK->destroy_kokkos(k_phir,phir);

  memoryKK->create_kokkos(k_phir, phir, (neltypes * (neltypes + 1)) / 2, nr, "pair:phir");
  memoryKK->create_kokkos(k_phirar, phirar, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar");
  memoryKK->create_kokkos(k_phirar1, phirar1, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar1");
  memoryKK->create_kokkos(k_phirar2, phirar2, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar2");
  memoryKK->create_kokkos(k_phirar3, phirar3, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar3");
  memoryKK->create_kokkos(k_phirar4, phirar4, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar4");
  memoryKK->create_kokkos(k_phirar5, phirar5, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar5");
  memoryKK->create_kokkos(k_phirar6, phirar6, (neltypes * (neltypes + 1)) / 2, nr, "pair:phirar6");

  h_phir = k_phir.h_view;
  h_phirar = k_phirar.h_view;
  h_phirar1 = k_phirar1.h_view;
  h_phirar2 = k_phirar2.h_view;
  h_phirar3 = k_phirar3.h_view;
  h_phirar4 = k_phirar4.h_view;
  h_phirar5 = k_phirar5.h_view;
  h_phirar6 = k_phirar6.h_view;

  for (int i = 0; i <(neltypes * (neltypes + 1)) / 2; i++)
    for(int j = 0; j < nr; j++)
    {
        h_phir(i,j) = phir[i][j];
        h_phirar(i,j) = phirar[i][j];
        h_phirar1(i,j) = phirar1[i][j];
        h_phirar2(i,j) = phirar2[i][j];
        h_phirar3(i,j) = phirar3[i][j];
        h_phirar4(i,j) = phirar4[i][j];
        h_phirar5(i,j) = phirar5[i][j];
        h_phirar6(i,j) = phirar6[i][j];
    }
  k_phir.template modify<LMPHostType>();
  k_phir.template sync<DeviceType>();  
  d_phir = k_phir.template view<DeviceType>();
  k_phirar.template modify<LMPHostType>();
  k_phirar.template sync<DeviceType>();  
  d_phirar = k_phirar.template view<DeviceType>();
  k_phirar1.template modify<LMPHostType>();
  k_phirar1.template sync<DeviceType>();  
  d_phirar1 = k_phirar1.template view<DeviceType>();
  k_phirar2.template modify<LMPHostType>();
  k_phirar2.template sync<DeviceType>();  
  d_phirar2 = k_phirar2.template view<DeviceType>();
  k_phirar3.template modify<LMPHostType>();
  k_phirar3.template sync<DeviceType>();  
  d_phirar3 = k_phirar3.template view<DeviceType>();
  k_phirar4.template modify<LMPHostType>();
  k_phirar4.template sync<DeviceType>();  
  d_phirar4 = k_phirar4.template view<DeviceType>();
  k_phirar5.template modify<LMPHostType>();
  k_phirar5.template sync<DeviceType>();  
  d_phirar5 = k_phirar5.template view<DeviceType>();
  k_phirar6.template modify<LMPHostType>();
  k_phirar6.template sync<DeviceType>();  
  d_phirar6 = k_phirar6.template view<DeviceType>();
}
