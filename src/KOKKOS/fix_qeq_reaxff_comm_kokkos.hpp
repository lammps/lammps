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
   Contributing authors: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

/*
   EXPLANATION OF WHY THESE FUNCTIONS ARE NEEDED:

   In parallel sparse matrix-vector multiplication with half neighbor lists:

   1. Each processor computes H*s for its local atoms
   2. Due to half neighbor lists, some interactions are "owned" by processor A
      but affect atoms owned by processor B
   3. Processor A computes these contributions but can't directly update
      processor B's atoms
   4. REVERSE COMMUNICATION solves this:
      - pack_reverse_comm: Processor A packs contributions for atoms owned by B
      - These get sent to processor B
      - unpack_reverse_comm: Processor B accumulates these contributions
        into its local result

   The key insight: Forward communication distributes data, reverse communication
   accumulates results.

   Example with 2 processors:
   - Proc 0 owns atoms 1-1000, has ghosts 1001-1010
   - Proc 1 owns atoms 1001-2000, has ghosts 991-1000
   
   When Proc 0 computes H*s:
   - It computes contributions for atoms 1-1000 (local) ✓
   - It also computes contributions for atoms 1001-1010 (ghosts)
   - These ghost contributions must be sent back to Proc 1
   - Proc 1 then adds these to its local computation for atoms 1001-1010

  UNDERSTANDING THE MATHEMATICAL FLOW:

   Think of reverse communication as solving this distributed computing problem:
   "I computed some results that belong to atoms owned by other processors -
   how do I get those results to where they belong and combine them properly?"

   Here's the step-by-step flow during matrix-vector multiplication H*s:

   1. COMPUTATION PHASE: Each processor computes H*s for its local atoms
      - Some interactions involve ghost atoms (owned by other processors)
      - These computations produce partial contributions for those ghost atoms

   2. PACK PHASE (pack_reverse_comm_kokkos):
      - Each processor packs up the partial contributions it computed
      - These are contributions for atoms it doesn't own but computed interactions with
      - Think: "Here are the partial results I computed for your atoms"

   3. COMMUNICATION PHASE (handled by LAMMPS communication infrastructure):
      - MPI sends these packed contributions to the processors that own those atoms
      - Multiple processors might send contributions for the same atom

   4. UNPACK PHASE (unpack_reverse_comm_kokkos):
      - Each processor receives contributions for atoms it owns
      - These contributions get ADDED together (not replaced!)
      - Think: "I'm collecting all the partial results computed for my atoms"

   The atomic_add operation is crucial because multiple processors might have
   computed contributions for the same atom, and we need to sum them all up
   to get the correct total contribution.

   WHY LAMBDAS ARE BETTER HERE:

   1. READABILITY: The logic is right there in the function, easy to understand
   2. MAINTAINABILITY: No need to hunt through the code for operator() definitions
   3. FLEXIBILITY: Each lambda can capture exactly what it needs
   4. DEBUGGING: Easier to set breakpoints and inspect values
   5. MODERN C++: Follows current best practices for Kokkos development

   The named kernels ("QEq::pack_reverse_comm") help with profiling and debugging -
   when you look at performance traces, you'll see exactly which kernels are
   taking time.

*/

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                                      int /*pbc_flag*/, int * /*pbc*/)
{
  if (pack_flag == 1) {
    k_s.template sync<LMPHostType>();
    for (int m = 0; m < n; m++)
      buf[m] = k_s.h_view(list[m]);
  } else if (pack_flag == 2) {
    for (int m = 0; m < n; m++)
      buf[m] = atom->q[list[m]];
  }
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  if (pack_flag == 1) {
    for (int m = 0, i = first; m < n; m++, i++)
      k_s.h_view(i) = buf[m];
    k_s.template modify<LMPHostType>();
  } else if (pack_flag == 2) {
    for (int m = 0, i = first; m < n; m++, i++)
      atom->q[i] = buf[m];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_forward_comm_kokkos(
    int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_xfloat_1d &k_buf,
    int pbc_flag, int *pbc)
{
  // Current pack flag determines what data to send
  int current_pack_flag = pack_flag;
  
  // Sync dual views
  k_sendlist.template sync<DeviceType>();
  k_buf.template sync<DeviceType>();
  
  // Create device views
  auto d_sendlist = k_sendlist.template view<DeviceType>();
  auto d_buf = k_buf.template view<DeviceType>();
  
  // Sync appropriate data based on pack_flag
  if (pack_flag == 1) {
    k_s.template sync<DeviceType>();
    auto d_s = k_s.template view<DeviceType>();
    
    // Pack s values
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int& i) {
        int j = d_sendlist(i);
        d_buf(i) = d_s(j);
      });
  } else if (pack_flag == 2) {
    // For charges, we might need to sync appropriately
    atomKK->sync(Device,Q_MASK);

    auto d_q = atomKK->k_q.view<DeviceType>();

    // Pack q values
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int& i) {
        int j = d_sendlist(i);
        d_buf(i) = d_q(j);
      });
  }
  
  // Mark buffer as modified
  k_buf.template modify<DeviceType>();
  
  // Return number of values packed
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_forward_comm_kokkos(
    int n, int first_in, DAT::tdual_xfloat_1d &k_buf)
{
  // Current pack flag determines what data to unpack
  int current_pack_flag = pack_flag;
  
  // Sync buffer
  k_buf.template sync<DeviceType>();
  auto d_buf = k_buf.template view<DeviceType>();
  
  // Store first index (offset for ghost atoms)
  int first = first_in;
  
  if (current_pack_flag == 1) {
    // Sync s values
    k_s.template sync<DeviceType>();
    auto d_s = k_s.template view<DeviceType>();
    
    // Unpack s values
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int& i) { d_s(i + first) = d_buf(i); });

    // Mark as modified
    k_s.template modify<DeviceType>();
  } else if (current_pack_flag == 2) {
    // Sync charge values
    atomKK->sync(Device,Q_MASK);

    auto d_q = atomKK->k_q.view<DeviceType>();

    // Unpack q values
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int& i) { d_q(i + first) = d_buf(i); });
    
    // Mark as modified
    atomKK->modified(Device,Q_MASK);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  // Host version: pack contributions from ghost atoms to send back to owners
  // During matrix-vector multiply, we computed partial contributions for ghost atoms
  // These need to be sent back to the processors that actually own those atoms
  
  k_o.template sync<LMPHostType>();
  
  // Simple loop to pack the data - starting from 'first' ghost atom
  int m = 0;
  for (int i = first; i < first + n; i++) {
    buf[m++] = k_o.h_view(i);
  }
  
  return m;  // Return number of values packed
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  // Host version: unpack and ACCUMULATE contributions from other processors
  // This is the crucial difference from forward comm - we add, not replace
  
  k_o.template sync<LMPHostType>();
  
  // Add each contribution to the appropriate owned atom
  for (int m = 0; m < n; m++) {
    k_o.h_view(list[m]) += buf[m];  // Key: += for accumulation
  }
  
  k_o.template modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_reverse_comm_kokkos(int n, int first_in,
                                                             DAT::tdual_xfloat_1d &buf)
{
  // Kokkos version with proper signature and lambda
  // This packs contributions from ghost atoms (starting at first_in)
  // that need to be sent back to their owning processors
  
  // Set up device views - this is the standard pattern
  buf.template sync<DeviceType>();
  k_o.template sync<DeviceType>();
  
  auto d_buf = buf.template view<DeviceType>();
  auto d_o = k_o.template view<DeviceType>();
  
  // The lambda captures the offset (first_in) and packs sequential data
  Kokkos::parallel_for("QEq::pack_reverse_comm", 
    Kokkos::RangePolicy<DeviceType>(0, n),
    KOKKOS_LAMBDA(const int i) {
      // Pack the contribution from ghost atom (first_in + i)
      // This represents partial matrix-vector product contribution
      // that was computed locally but belongs to an atom owned elsewhere
      d_buf(i) = d_o(first_in + i);
    });
  
  // Mark buffer as modified on device
  buf.template modify<DeviceType>();
  
  return n;  // Return number of values packed
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_reverse_comm_kokkos(int n,
                                                                DAT::tdual_int_1d k_sendlist,
                                                                DAT::tdual_xfloat_1d &buf)
{
  // This accumulates contributions sent from other processors
  // k_sendlist tells us which atoms to update with the received data
  
  // Set up device views
  k_sendlist.template sync<DeviceType>();
  buf.template sync<DeviceType>();
  k_o.template sync<DeviceType>();
  
  auto d_sendlist = k_sendlist.template view<DeviceType>();
  auto d_buf = buf.template view<DeviceType>();
  auto d_o = k_o.template view<DeviceType>();
  
  // Use lambda with atomic operations for thread-safe accumulation
  // This is where the mathematical magic happens - we're summing up
  // partial contributions from multiple processors
  Kokkos::parallel_for("QEq::unpack_reverse_comm",
    Kokkos::RangePolicy<DeviceType>(0, n), 
    KOKKOS_LAMBDA(const int i) {
      // Get the local atom index that should receive this contribution
      int atom_idx = d_sendlist(i);
      
      // Atomically add the contribution to avoid race conditions
      // Multiple threads might be updating the same atom simultaneously
      // This ensures we get the correct sum of all contributions
      Kokkos::atomic_add(&d_o(atom_idx), d_buf(i));
    });
  
  // Mark the output array as modified on device
  k_o.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_xfloat_2d &k_buf,
   DAT::tdual_int_1d k_exchange_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace space)
{
  // Sync all dual views to ensure device has latest data
  k_buf.template sync<DeviceType>();
  k_copylist.template sync<DeviceType>();
  k_exchange_sendlist.template sync<DeviceType>();
  k_theta.template sync<DeviceType>();
  k_theta_dot.template sync<DeviceType>();
  
  // Create device views from dual views
  auto d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  auto d_copylist = k_copylist.template view<DeviceType>();
  auto d_exchange_sendlist = k_exchange_sendlist.template view<DeviceType>();
  
  // Get device access to extended Lagrangian variables
  auto d_theta = k_theta.template view<DeviceType>();
  auto d_theta_dot = k_theta_dot.template view<DeviceType>();
  
  copymode = 1;
  
  // Optimization: Use a more efficient packing approach
  // Buffer structure is known: 2 values per atom (theta, theta_dot)
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nsend),
    KOKKOS_LAMBDA(const int& isend) {
      // Get the actual atom index from our send list
      const int i = d_exchange_sendlist(isend);
      const int buf_offset = isend * 2;
      
      // Pack theta and theta_dot into the buffer
      d_buf(buf_offset) = d_theta(i);
      d_buf(buf_offset + 1) = d_theta_dot(i);
      
      // Handle any copy operations if needed
      const int j = d_copylist(isend);
      if (j >= 0) {
        // Copy from i to j
        d_theta(j) = d_theta(i);
        d_theta_dot(j) = d_theta_dot(i);
      }
    });
  
  copymode = 0;
  
  // Mark views as modified on device
  k_theta.template modify<DeviceType>();
  k_theta_dot.template modify<DeviceType>();
  
  // Sync buffer if needed based on execution space
  k_buf.template modify<DeviceType>();
  if (space == Host) k_buf.template sync<LMPHostType>();
  else k_buf.template sync<LMPDeviceType>();
  
  // Return the number of values packed (2 per atom)
  return nsend*2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_exchange_kokkos(
   DAT::tdual_xfloat_2d &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
   int nrecv1, int nextrarecv1, ExecutionSpace space)
{
  // Sync input dual views to device
  k_buf.template sync<DeviceType>();
  k_indices.template sync<DeviceType>();
  
  // Sync extended Lagrangian variables
  k_theta.template sync<DeviceType>();
  k_theta_dot.template sync<DeviceType>();
  
  // Create device views
  auto d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  auto d_indices = k_indices.template view<DeviceType>();
  
  // Get device access to extended Lagrangian variables
  auto d_theta = k_theta.template view<DeviceType>();
  auto d_theta_dot = k_theta_dot.template view<DeviceType>();
  
  copymode = 1;
  
  // Optimization: Improved unpacking with better memory access pattern
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nrecv),
    KOKKOS_LAMBDA(const int& i) {
      // Get the atom index this data belongs to
      int index = d_indices(i);
      
      // Only unpack if this is a valid atom (index > -1)
      if (index > -1) {
        const int buf_offset = i * 2;
        // Unpack in a memory-coalesced manner
        d_theta(index) = d_buf(buf_offset);
        d_theta_dot(index) = d_buf(buf_offset + 1);
      }
    });
  
  copymode = 0;
  
  // Mark views as modified on device
  k_theta.template modify<DeviceType>();
  k_theta_dot.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_theta.template sync<LMPHostType>();
  k_theta_dot.template sync<LMPHostType>();

  k_theta.resize(nmax);
  k_theta_dot.resize(nmax);

  k_theta.template modify<LMPHostType>();
  k_theta_dot.template modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::copy_arrays(int i, int j, int /*delflag*/)
{
  k_theta.template sync<LMPHostType>();
  k_theta_dot.template sync<LMPHostType>();

  k_theta.h_view(j) = k_theta.h_view(i);
  k_theta_dot.h_view(j) = k_theta_dot.h_view(i);

  k_theta.template modify<LMPHostType>();
  k_theta_dot.template modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  k_theta.sync_device();
  k_theta_dot.sync_device();

  Sorter.sort(LMPDeviceType(), k_theta.d_view);
  Sorter.sort(LMPDeviceType(), k_theta_dot.d_view);

  k_theta.modify_device();
  k_theta_dot.modify_device();
}
