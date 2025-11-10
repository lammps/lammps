#ifndef COLVARATOMS_SOA_H
#define COLVARATOMS_SOA_H

#include "colvarmodule.h"
#include "colvardeps.h"
#include "colvar_rotation_derivative.h"
#include <unordered_map>
#include <mutex>

/**
 * @brief Store the information of a group of atoms in a structure-of-arrays
 *        (SoA) style
 *
 * This class keeps the information of a group of atoms. These information
 * includes the internal indices in the colvarproxy class, positions,
 * charges, velocities, masses, gradients, total forces and atom IDs of
 * all atoms in the atom group. The memory layout of this class is different
 * from the traditional \link cvm::atom_group \endlink class as this class
 * uses the SoA layout to group each atom properties in their own array
 * instead of using a single array of \link cvm::atom \endlink
 * (array-of-structures or AoS). According to various micro-benchmarks,
 * the use of SoA can significantly improve the
 * performance. There are two disadvantages of adopting SoA:
 *
 * 1. The insertion or delete of atoms in an object of \link cvm::atom_group
 *    \endlink is not as handy as the AoS class \link cvm::atom_group \endlink.
 *    To solve this issue, this nested class \p cvm::atom_group::atom_modifier
 *    has been implemented as a temporary AoS class that shares the same
 *    interfaces like \p cvm::atom_group::atom_modifier::add_atom_numbers
 *    as the traditional AoS class \link cvm::atom_group \endlink. To
 *    add or remove atoms in the SoA class, you can do the following:
 *    @code
 *    {
 *      auto modify_atoms = atoms->get_atom_modifier();
 *      modify_atoms->add_atom_numbers(conf);
 *    }
 *    @endcode
 *    Specifically, when \p get_atom_modifier is called, the temporary AoS atom
 *    group will be created based on this SoA atom group, and after the
 *    modification of the temporary AoS atom group is completed, the
 *    destructor \p cvm::atom_group::atom_modifier::~atom_modifier
 *    will synchronize the changes back to this SoA atom group.
 *
 * 2. It would be difficult to reorder the atoms inside an SoA group.
 *    In MD engines that need to move atoms between cells/patches (like NAMD),
 *    it is better to transform the atoms back to AoS before moving (migration).
 *    Fortunately Colvars does not need this as the order of atoms inside
 *    a group never changes.
 *
 */
class cvm::atom_group: public colvarparse, public colvardeps {
public:
  /**
   *  @brief A helper function to re-arrange the a vector of
   *         \link cvm::atom_pos \endlink (in AoS style xyz...xyz) to
   *         an SoA vector (x...xy...yz...z)
   */
  static std::vector<cvm::real> pos_aos_to_soa(const std::vector<cvm::atom_pos>& aos_in);
  /**
   * @brief Default constructor
   *
   * @attention Will call \p init_dependencies() to initialize the
   *            dependency tree.
   */
  atom_group();
  /**
   * @brief Create a group object, assign a name to it
   */
  atom_group(char const *key_in);
  /**
   * @brief Destructor
   *
   * @attention Will call \p clear_soa() to de-reference all atoms (if not scalable)
   *            from the proxy class.
   */
  ~atom_group() override;
  /**
   * @brief Set default values for common flags
   */
  int init();
  /**
   * @brief Update data required to calculate cvc's
   *
   * The data updated includes atomic charges and masses, and also the total
   * charge and total mass of this atom group.
   */
  int setup();
  /**
   * @brief Initialize dependency tree
   */
  int init_dependencies() override;
  /// \brief Implementation of the feature list accessor for atom group
  const std::vector<feature *> &features() const override { return ag_features; }
  std::vector<feature *> &modify_features() override { return ag_features; }
  static void delete_features()
  {
    for (size_t i = 0; i < ag_features.size(); i++) {
      delete ag_features[i];
    }
    ag_features.clear();
  }
  /**
   * @brief A simplified class of \p cvm::atom that can be used with
   *        \p cvm::atom_group::atom_modifier
   */
  struct simple_atom {
    int proxy_index;
    int id;
    cvm::real mass;
    cvm::real charge;
    cvm::atom_pos pos;
    cvm::rvector vel;
    cvm::rvector total_force;
    cvm::rvector grad;
  };
  /**
   * @brief Initialize an instance of \p cvm::atom_group::simple_atom
   *        from \p colvarproxy
   * @param p The pointer to the \p colvarproxy instance.
   * @param residue The residue number of the atom
   * @param atom_name The name of the atom in the residue
   * @param segment_id For PSF topologies, the segment identifier; for other
   *                   type of topologies, may not be required
   */
  static simple_atom init_atom_from_proxy(
    colvarproxy* const p,
    cvm::residue_id const &residue,
    std::string const     &atom_name,
    std::string const     &segment_id);
  /**
   * @brief Initialize an instance of \p cvm::atom_group::simple_atom
   *        from \p colvarproxy
   * @param p The pointer to the \p colvarproxy instance.
   * @param atom_number Atom index in the system topology (1-based)
   */
  static simple_atom init_atom_from_proxy(
    colvarproxy* const p,
    int atom_number);
  /**
   * @brief Initialize an instance of \p cvm::atom_group::simple_atom
   *        from \p colvarproxy
   * @param p The pointer to the \p colvarproxy instance.
   * @param simple_atom The existing object of \p cvm::atom_group::simple_atom
   *
   * @attention The \p cvm::atom_group::simple_atom::proxy_index of
   *            input atom must be known. The reference count in the proxy
   *            will be increased.
   */
  static simple_atom init_atom_from_proxy(
    colvarproxy* const p,
    const simple_atom& atom);
  /**
   * @brief The temporary AoS interface to modify \p cvm::atom_group
   */
  class atom_modifier {
  private:
    /**
     * @brief Pointer to the SoA atom group to be modified
     */
    cvm::atom_group* m_ag;
    /**
     * @brief Internal atom IDs (populated during initialization)
     */
    std::vector<int> m_atoms_ids;
    /**
     * @brief AoS layout of the atoms
     */
    std::vector<simple_atom> m_atoms;
    /**
     * @brief A map to avoid doubly counting atoms
     */
    std::unordered_map<int, int> m_atoms_ids_count;
    /**
     * @brief Total mass of the atom group
     */
    cvm::real m_total_mass;
    /**
     * @brief Total charge of the atom group
     */
    cvm::real m_total_charge;
    /**
     * @brief Populate this temporary AoS object with the atoms from \p m_ag
     */
    void update_from_soa();
    /**
     * @brief Synchronize the atoms from this temporary AoS object to \p m_ag
     */
    void sync_to_soa() const;
  public:
    using atom_iter = decltype(m_atoms)::iterator;
    using const_atom_iter = decltype(m_atoms)::const_iterator;
    /**
     * @brief Construct from \p cvm::atom_group
     *
     * @attention Will call \p cvm::atom_group::atom_modifier::update_from_soa()
     */
    atom_modifier(cvm::atom_group* ag);
    /**
     * @brief Destructor
     *
     * @attention Will call \p cvm::atom_group::atom_modifier::sync_to_soa()
     */
    ~atom_modifier();
    /**
     * @name STL-style iterators and accessors
     */
    ///@{
    inline simple_atom & operator [] (size_t const i){return m_atoms[i];}
    inline simple_atom const & operator [] (size_t const i) const {return m_atoms[i];}
    inline atom_iter begin(){return m_atoms.begin();}
    inline const_atom_iter begin() const {return m_atoms.begin();}
    inline atom_iter end() {return m_atoms.end();}
    inline const_atom_iter end() const {return m_atoms.end();}
    inline size_t size() const {return m_atoms.size();}
    ///@}
    /**
     * @name Add atom(s) to the atom group
     */
    ///@{
    int add_atom(simple_atom const &a);
    int add_atom_numbers(std::string const &numbers_conf);
    int add_atoms_of_group(const atom_group *ag);
    int add_index_group(std::string const &index_group_name, bool silent = false);
    int add_atom_numbers_range(std::string const &range_conf);
    int add_atom_name_residue_range(std::string const &psf_segid,
                                    std::string const &range_conf);
    int add_atom_id(int aid);
    ///@}
    /**
     * @brief Remove an atom object from this group
     */
    int remove_atom(atom_iter ai);
  };
  /**
   * @brief Get the atom modifier object associated with this SoA atom group.
   *
   * @attention You can get only one instance of \p cvm::atom_group::atom_modifier
   *            in a scope. To prevent acquiring multiple instances of \p atom_modifier
   *            at the same time which may corrupt the SoA layout. The associated
   *            \p cvm::atom_group::atom_modifier will try locking \p mutex_lock.
   */
  atom_modifier get_atom_modifier() {
    return atom_modifier(this);
  }
  /**
   * @brief Remove all atoms from this SoA group and the proxy
   *
   * @attention \p atoms_ids , \p sorted_atoms_ids , \p sorted_atoms_ids_map
   *            total mass, total charge, reference positions and
   *            the fitting group are kept.
   */
  void clear_soa();
  /**
   * @brief Initialize the group by looking up its configuration
   */
  int parse(std::string const &conf);
  /**
   * @brief Initialize the fitting group by looking up its configuration
   */
  int parse_fitting_options(std::string const &group_conf);
  /**
   * @brief Set this group as a dummy group (no actual atoms)
   */
  int set_dummy();
  /**
   * @brief If this group is dummy, set the corresponding position
   */
  int set_dummy_pos(cvm::atom_pos const &pos);
  /**
   * @brief Update the total mass of the atom group
   */
  void update_total_mass();
  /**
   * @brief Update the total mass of the group
   */
  void update_total_charge();
  /**
   * @brief Print the updated the total mass and charge of a group.
   *
   * This is needed in case the hosting MD code has an option to
   * change atom masses after their initialization.
   */
  void print_properties(std::string const &colvar_name, int i, int j);
  /**
   * @brief Internal atom IDs (populated during initialization)
   */
  inline std::vector<int> const &ids() const
  {
    return atoms_ids;
  }
  std::string const print_atom_ids() const;
  /**
   * @brief Allocates and populates sorted_ids and sorted_ids_map
   */
  int create_sorted_ids();
  /**
   * @brief Sorted internal atom IDs (populated on-demand by create_sorted_ids);
   *
   * Used to read coordinate files
   */
  inline std::vector<int> const &sorted_ids() const
  {
    return sorted_atoms_ids;
  }
  /// Map entries of sorted_atoms_ids onto the original positions in the group
  inline std::vector<int> const &sorted_ids_map() const
  {
    return sorted_atoms_ids_map;
  }
  /**
   * @brief Detect whether two groups share atoms
   *
   * If yes, returns 1-based number of a common atom; else, returns 0
   */
  static int overlap(const atom_group &g1, const atom_group &g2);
  /**
   * @brief Get the current positions
   */
  void read_positions();
  /**
   * @brief (Re)calculate the optimal roto-translation
   */
  void calc_apply_roto_translation();
  /**
   * @brief Setup the \p rot_deriv object for the calculation of derivative of
   *        the optimal roto-translation
   */
  void setup_rotation_derivative();
  /**
   * @brief Save aside the center of geometry of the reference positions,
   *       then subtract it from them
   *
   * In this way it will be possible to use ref_pos also for the
   * rotational fit.
   * This is called either by atom_group::parse or by CVCs that assign
   * reference positions (eg. RMSD, eigenvector).
   */
  void center_ref_pos();
  /**
   * @brief Rotate all atoms by a rotation matrix
   */
  void rotate(const cvm::rmatrix& rot_mat);
  /**
   * @brief Move all positions
   */
  void apply_translation(cvm::rvector const &t);
  /**
   * @brief Get the current velocities
   *
   * this must be called always *after* read_positions();
   * if f_ag_rotate is defined, the same rotation applied
   * to the coordinates will be used.
   */
  void read_velocities();
  /**
   * @brief Get the current total_forces
   *
   * this must be called always *after* read_positions();
   * if f_ag_rotate is defined, the same rotation applied
   * to the coordinates will be used.
   */
  void read_total_forces();
  /**
   * @brief Recompute all mutable quantities that are required to compute CVCs
   */
  int calc_required_properties();
  /**
   * @brief Calculate the center of geometry of the atomic positions, assuming
   *        that they are already pbc-wrapped
   */
  int calc_center_of_geometry();
  /**
   * @brief Return a copy of the current atom positions
   */
  std::vector<cvm::real> positions() const;
  /**
   * @brief Return the center of geometry of the atomic positions
   */
  inline cvm::atom_pos center_of_geometry() const
  {
    return cog;
  }
  /**
   * @brief Calculate the center of mass of the atomic positions, assuming that
   * they are already pbc-wrapped
   */
  int calc_center_of_mass();
  /**
   * @brief Return the center of mass (COM) of the atomic positions
   */
  inline cvm::atom_pos center_of_mass() const
  {
    return com;
  }
  /**
   * @brief Return previously gradient of scalar variable with respect to the
   * COM
   */
  inline cvm::rvector center_of_mass_scalar_gradient() const
  {
    return scalar_com_gradient;
  }
  /**
   * @brief Return a copy of the current atom positions, shifted by a constant vector
   */
  std::vector<cvm::real> positions_shifted(cvm::rvector const &shift) const;
  /**
   * @brief Return a copy of the current atom velocities
   */
  std::vector<cvm::real> velocities() const;
  /**
   * @brief Calculate the dipole of the atom group around the specified center
   */
  int calc_dipole(cvm::atom_pos const &dipole_center);
  /**
   * @brief Return the (previously calculated) dipole of the atom group
   */
  inline cvm::rvector dipole() const
  {
    return dip;
  }
  /**
   * @brief Return a copy of the total forces
   */
  std::vector<cvm::real> total_forces() const;
  /**
   * @brief Return a copy of the aggregated total force on the group
   */
  cvm::rvector total_force() const;
  /**
   * @brief Shorthand: save the specified gradient on each atom,
   * weighting with the atom mass (mostly used in combination with
   * \p cvm::atom_group::center_of_mass() )
   */
  void set_weighted_gradient(cvm::rvector const &grad);
  /**
   * @brief Calculate the derivatives of the fitting transformation
   */
  void calc_fit_gradients();
  /*! @brief  Actual implementation of `calc_fit_gradients` and
   *          `calc_fit_forces`. The template is
   *          used to avoid branching inside the loops in case that the CPU
   *          branch prediction is broken (or further migration to GPU code).
   *  @tparam B_ag_center Centered the reference to origin? This should follow
   *          the value of `is_enabled(f_ag_center)`.
   *  @tparam B_ag_rotate Calculate the optimal rotation? This should follow
   *          the value of `is_enabled(f_ag_rotate)`.
   *  @tparam main_force_accessor_T The type of accessor of the main
   *          group forces or gradients acting on the rotated frame.
   *  @tparam fitting_force_accessor_T The type of accessor of the fitting group
   *          forces or gradients.
   *  @param accessor_main The accessor of the main group forces or gradients.
   *         accessor_main(i) should return the i-th force or gradient of the
   *         rotated main group.
   *  @param accessor_fitting The accessor of the fitting group forces or gradients.
   *         accessor_fitting(j, v) should store/apply the j-th atom gradient or
   *         force in the fitting group.
   *
   *  This function is used to (i) project the gradients of CV with respect to
   *  rotated main group atoms to fitting group atoms, or (ii) project the forces
   *  on rotated main group atoms to fitting group atoms, by the following two steps
   *  (using the goal (ii) for example):
   *  (1) Loop over the positions of main group atoms and call cvm::quaternion::position_derivative_inner
   *      to project the forces on rotated main group atoms to the forces on quaternion.
   *  (2) Loop over the positions of fitting group atoms, compute the gradients of
   *      \f$\mathbf{q}\f$ with respect to the position of each atom, and then multiply
   *      that with the force on \f$\mathbf{q}\f$ (chain rule).
   */
  template <bool B_ag_center, bool B_ag_rotate,
            typename main_force_accessor_T,
            typename fitting_force_accessor_T>
  void calc_fit_forces_impl(
    main_force_accessor_T accessor_main,
    fitting_force_accessor_T accessor_fitting) const;

/*! @brief  Calculate or apply the fitting group forces from the main group forces.
 *  @tparam main_force_accessor_T The type of accessor of the main
 *          group forces or gradients.
 *  @tparam fitting_force_accessor_T The type of accessor of the fitting group
 *          forces or gradients.
 *  @param accessor_main The accessor of the main group forces or gradients.
 *         accessor_main(i) should return the i-th force or gradient of the
 *         main group.
 *  @param accessor_fitting The accessor of the fitting group forces or gradients.
 *         accessor_fitting(j, v) should store/apply the j-th atom gradient or
 *         force in the fitting group.
 *
 *  This function just dispatches the parameters to calc_fit_forces_impl that really
 *  performs the calculations.
 */
  template <typename main_force_accessor_T, typename fitting_force_accessor_T>
  void calc_fit_forces(
    main_force_accessor_T accessor_main,
    fitting_force_accessor_T accessor_fitting) const;
  /**
   * @brief Used by a (scalar) colvar to apply its force on its \link
   * atom_group \endlink members
   *
   * The (scalar) force is multiplied by the colvar gradient for each
   * atom; this should be used when a colvar with scalar \link
   * colvarvalue \endlink type is used (this is the most frequent
   * case: for colvars with a non-scalar type, the most convenient
   * solution is to sum together the Cartesian forces from all the
   * colvar components, and use apply_force() or apply_forces()).  If
   * the group is being rotated to a reference frame (e.g. to express
   * the colvar independently from the solute rotation), the
   * gradients are temporarily rotated to the original frame.
   */
  void apply_colvar_force(cvm::real const &force);
  /**
   * @brief Apply a force "to the center of mass", i.e. the force is
   * distributed on each atom according to its mass
   *
   * If the group is being rotated to a reference frame (e.g. to
   * express the colvar independently from the solute rotation), the
   * force is rotated back to the original frame.  Colvar gradients
   * are not used, either because they were not defined (e.g because
   * the colvar has not a scalar value) or the biases require to
   * micromanage the force.
   * This function will be phased out eventually, in favor of
   * apply_colvar_force() once that is implemented for non-scalar values
   */
  void apply_force(cvm::rvector const &force);
  /**
   * Implements possible actions to be carried out
   * when a given feature is enabled
   * This overloads the base function in colvardeps
   */
  void do_feature_side_effects(int id) override;
  /**
   * @brief Return the number of atoms of this group
   */
  size_t size() const {return num_atoms;}
  /**
   * @brief Clear the atom positions, velocities, gradients and total forces
   */
  void reset_atoms_data();
  /**
   * @brief Setup reference positions from AOS
   */
  void set_ref_pos_from_aos(const std::vector<cvm::atom_pos>& pos_aos);
  /**
   * @name Accessors to atom IDs
   */
  ///@{
  inline int& id(size_t i) {return atoms_ids[i];}
  inline const int& id(size_t i) const {return atoms_ids[i];}
  ///@}
  /**
   * @name Accessors to positions
   */
  ///@{
  inline cvm::real& pos_x(size_t i) {return atoms_pos[i];}
  inline cvm::real& pos_y(size_t i) {return atoms_pos[i + num_atoms];}
  inline cvm::real& pos_z(size_t i) {return atoms_pos[i + 2 * num_atoms];}
  inline const cvm::real& pos_x(size_t i) const {return atoms_pos[i];}
  inline const cvm::real& pos_y(size_t i) const {return atoms_pos[i + num_atoms];}
  inline const cvm::real& pos_z(size_t i) const {return atoms_pos[i + 2 * num_atoms];}
  ///@}
  /**
   * @name Accessors to velocities
   */
  ///@{
  inline cvm::real& vel_x(size_t i) {return atoms_vel[i];}
  inline cvm::real& vel_y(size_t i) {return atoms_vel[i + num_atoms];}
  inline cvm::real& vel_z(size_t i) {return atoms_vel[i + 2 * num_atoms];}
  inline const cvm::real& vel_x(size_t i) const {return atoms_vel[i];}
  inline const cvm::real& vel_y(size_t i) const {return atoms_vel[i + num_atoms];}
  inline const cvm::real& vel_z(size_t i) const {return atoms_vel[i + 2 * num_atoms];}
  ///@}
  /**
   * @name Accessors to gradients
   */
  ///@{
  inline cvm::real& grad_x(size_t i) {return atoms_grad[i];}
  inline cvm::real& grad_y(size_t i) {return atoms_grad[i + num_atoms];}
  inline cvm::real& grad_z(size_t i) {return atoms_grad[i + 2 * num_atoms];}
  inline const cvm::real& grad_x(size_t i) const {return atoms_grad[i];}
  inline const cvm::real& grad_y(size_t i) const {return atoms_grad[i + num_atoms];}
  inline const cvm::real& grad_z(size_t i) const {return atoms_grad[i + 2 * num_atoms];}
  ///@}
  /**
   * @name Accessors to total forces
   */
  ///@{
  inline cvm::real& total_force_x(size_t i) {return atoms_total_force[i];}
  inline cvm::real& total_force_y(size_t i) {return atoms_total_force[i + num_atoms];}
  inline cvm::real& total_force_z(size_t i) {return atoms_total_force[i + 2 * num_atoms];}
  inline const cvm::real& total_force_x(size_t i) const {return atoms_total_force[i];}
  inline const cvm::real& total_force_y(size_t i) const {return atoms_total_force[i + num_atoms];}
  inline const cvm::real& total_force_z(size_t i) const {return atoms_total_force[i + 2 * num_atoms];}
  ///@}
  /**
   * @name Accessors to reference positions
   */
  ///@{
  inline cvm::real& ref_pos_x(size_t i) {return ref_pos[i];}
  inline cvm::real& ref_pos_y(size_t i) {return ref_pos[i + num_ref_pos];}
  inline cvm::real& ref_pos_z(size_t i) {return ref_pos[i + 2 * num_ref_pos];}
  inline const cvm::real& ref_pos_x(size_t i) const {return ref_pos[i];}
  inline const cvm::real& ref_pos_y(size_t i) const {return ref_pos[i + num_ref_pos];}
  inline const cvm::real& ref_pos_z(size_t i) const {return ref_pos[i + 2 * num_ref_pos];}
  ///@}
  /**
   * @name Accessors to positions before rotation
   */
  ///@{
  inline cvm::real& pos_unrotated_x(size_t i) {return atoms_pos_unrotated[i];}
  inline cvm::real& pos_unrotated_y(size_t i) {return atoms_pos_unrotated[i + num_atoms];}
  inline cvm::real& pos_unrotated_z(size_t i) {return atoms_pos_unrotated[i + 2 * num_atoms];}
  inline const cvm::real& pos_unrotated_x(size_t i) const {return atoms_pos_unrotated[i];}
  inline const cvm::real& pos_unrotated_y(size_t i) const {return atoms_pos_unrotated[i + num_atoms];}
  inline const cvm::real& pos_unrotated_z(size_t i) const {return atoms_pos_unrotated[i + 2 * num_atoms];}
  ///@}
  /**
   * @name Accessors to masses
   */
  ///@{
  inline cvm::real& mass(size_t i) {return atoms_mass[i];}
  inline const cvm::real& mass(size_t i) const {return atoms_mass[i];}
  ///@}
  /**
   * @name Accessors to weights
   */
  ///@{
  inline cvm::real& weight(size_t i) {return atoms_weight[i];}
  inline const cvm::real& weight(size_t i) const {return atoms_weight[i];}
  ///@}
  /**
   * @name Accessors to charges
   */
  ///@{
  inline cvm::real& charge(size_t i) {return atoms_charge[i];}
  inline const cvm::real& charge(size_t i) const {return atoms_charge[i];}
  ///@}
  /**
   * @name Accessors to fit gradients
   */
  ///@{
  inline cvm::real& fit_gradients_x(size_t i) {
    atom_group *group_for_fit = fitting_group ? fitting_group : this;
    return group_for_fit->fit_gradients[i];}
  inline cvm::real& fit_gradients_y(size_t i) {
    atom_group *group_for_fit = fitting_group ? fitting_group : this;
    return group_for_fit->fit_gradients[i + group_for_fit->size()];}
  inline cvm::real& fit_gradients_z(size_t i) {
    atom_group *group_for_fit = fitting_group ? fitting_group : this;
    return group_for_fit->fit_gradients[i + 2 * group_for_fit->size()];}
  inline const cvm::real& fit_gradients_x(size_t i) const {
    const atom_group *group_for_fit = fitting_group ? fitting_group : this;
    return group_for_fit->fit_gradients[i];}
  inline const cvm::real& fit_gradients_y(size_t i) const {
    const atom_group *group_for_fit = fitting_group ? fitting_group : this;
    return group_for_fit->fit_gradients[i + group_for_fit->size()];}
  inline const cvm::real& fit_gradients_z(size_t i) const {
    const atom_group *group_for_fit = fitting_group ? fitting_group : this;
    return group_for_fit->fit_gradients[i + 2 * group_for_fit->size()];}
  ///@}
  /**
   * @name Accessors to group forces
   */
  ///@{
  inline cvm::real& group_forces_x(size_t i) {return group_forces[i];}
  inline cvm::real& group_forces_y(size_t i) {return group_forces[i + num_atoms];}
  inline cvm::real& group_forces_z(size_t i) {return group_forces[i + 2 * num_atoms];}
  inline const cvm::real& group_forces_x(size_t i) const {return group_forces[i];}
  inline const cvm::real& group_forces_y(size_t i) const {return group_forces[i + num_atoms];}
  inline const cvm::real& group_forces_z(size_t i) const {return group_forces[i + 2 * num_atoms];}
  ///@}
  /**
   * @brief Read-only operator[]
   */
  inline simple_atom operator[](size_t i) const {
    return simple_atom{
      /*.proxy_index = */atoms_index[i],
      /*.id = */atoms_ids[i],
      /*.mass = */atoms_mass[i],
      /*.charge = */atoms_charge[i],
      /*.pos = */{pos_x(i), pos_y(i), pos_z(i)},
      /*.vel = */{vel_x(i), vel_y(i), vel_z(i)},
      /*.total_force = */{total_force_x(i), total_force_y(i), total_force_z(i)},
      /*.grad = */{grad_x(i), grad_y(i), grad_z(i)}
    };
  }

  /*! @class group_force_object
   *  @brief A helper class for applying forces on an atom group in a way that
   *         is aware of the fitting group. NOTE: you are encouraged to use
   *         get_group_force_object() to get an instance of group_force_object
   *         instead of constructing directly.
   */
  class group_force_object {
  public:
    /*! @brief Constructor of group_force_object
     *  @param ag The pointer to the atom group that forces will be applied on.
     */
    group_force_object(cvm::atom_group* ag);
    /*! @brief Destructor of group_force_object
     */
    ~group_force_object();
    /*! @brief Apply force to atom i
     *  @param i The i-th of atom in the atom group.
     *  @param force The force being added to atom i.
     *
     * The function can be used as follows,
     * @code
     *       // In your colvar::cvc::apply_force() loop of a component:
     *       auto ag_force = atoms->get_group_force_object();
     *       for (ia = 0; ia < atoms->size(); ia++) {
     *         const cvm::rvector f = compute_force_on_atom_ia();
     *         ag_force.add_atom_force(ia, f);
     *       }
     * @endcode
     * There are actually two scenarios under the hood:
     * (i) If the atom group does not have a fitting group, then the force is
     *     added to atom i directly;
     * (ii) If the atom group has a fitting group, the force on atom i will just
     *      be temporary stashed into ag->group_forces. At the end of the loop
     *      of apply_force(), the destructor ~group_force_object() will be called,
     *      which then call apply_force_with_fitting_group(). The forces on the
     *      main group will be rotated back by multiplying ag->group_forces with
     *      the inverse rotation. The forces on the fitting group (if
     *      enableFitGradients is on) will be calculated by calling
     *      calc_fit_forces.
     */
    void add_atom_force(size_t i, const cvm::rvector& force);
  private:
    cvm::atom_group* m_ag;
    cvm::atom_group* m_group_for_fit;
    bool m_has_fitting_force;
    void apply_force_with_fitting_group();
  };
  group_force_object get_group_force_object();
public:
  /// \brief Optional name to reuse properties of this in other groups
  std::string name;
  /// \brief If this option is on, this group merely acts as a wrapper
  /// for a fixed position; any calls to atoms within or to
  /// functions that return disaggregated data will fail
  bool b_dummy;
  /// \brief Keyword used to define the group
  // TODO Make this field part of the data structures that link a group to a CVC
  std::string key;
  /// \brief If f_ag_center or f_ag_rotate is true, use this group to
  /// define the transformation (default: this group itself)
  cvm::atom_group *fitting_group;
  /// The rotation calculated automatically if f_ag_rotate is defined
  cvm::rotation rot;
  /// \brief Don't apply any force on this group (use its coordinates
  /// only to calculate a colvar)
  bool noforce;
  /// \brief Indicates that the user has explicitly set centerToReference or
  /// rotateReference, and the corresponding reference:
  /// cvc's (eg rmsd, eigenvector) will not override the user's choice
  bool b_user_defined_fit;
  /// \brief Derivatives of the fitting transformation
  std::vector<cvm::real> fit_gradients;
  /// Total mass of the atom group
  cvm::real total_mass;
  /// Total charge of the atom group
  cvm::real total_charge;
  /// Rotation derivative;
  rotation_derivative* rot_deriv;
  /// \brief Implementation of the feature list for atom group
  static std::vector<feature *> ag_features;
private:
  /// \brief Number of atoms
  size_t num_atoms;
  /// \brief SOA atom indices (size: num_atoms)
  std::vector<int> atoms_index;
  /// \brief SOA atom positions (size: 3 * num_atoms)
  std::vector<cvm::real> atoms_pos;
  /// \brief SOA atom charges (size: num_atoms)
  std::vector<cvm::real> atoms_charge;
  /// \brief SOA atom velocities (size: 3 * num_atoms)
  std::vector<cvm::real> atoms_vel;
  /// \brief SOA atom mass (size: num_atoms)
  std::vector<cvm::real> atoms_mass;
  /// \brief SOA atom gradients (size: 3 * num_atoms)
  std::vector<cvm::real> atoms_grad;
  /// \brief SOA atom total forces (size: 3 * num_atoms)
  std::vector<cvm::real> atoms_total_force;
  /// \brief Atom masses divided by total mass (size: num_atoms)
  std::vector<cvm::real> atoms_weight;
  /// \brief Internal atom IDs for host code
  std::vector<int> atoms_ids;
  /// \brief Sorted list of internal atom IDs (populated on-demand by
  /// create_sorted_ids); used to read coordinate files
  std::vector<int> sorted_atoms_ids;
  /// \brief Map entries of sorted_atoms_ids onto the original positions in the group
  std::vector<int> sorted_atoms_ids_map;
  /// \brief Dummy atom position
  cvm::atom_pos dummy_atom_pos;
  /// \brief Index in the colvarproxy arrays (if the group is scalable)
  int index;
  /// \brief The temporary forces acting on the main group atoms.
  ///        Currently this is only used for calculating the fitting group forces for
  ///        non-scalar components.
  std::vector<cvm::real> group_forces;
  /// \brief use reference coordinates for f_ag_center or f_ag_rotate
  std::vector<cvm::real> ref_pos;
  size_t num_ref_pos; // TODO: Do I really need this?
  /// \brief Center of geometry of the reference coordinates; regardless
  /// of whether f_ag_center is true, ref_pos is centered to zero at
  /// initialization, and ref_pos_cog serves to center the positions
  cvm::atom_pos ref_pos_cog;
  /// \brief Center of geometry
  cvm::atom_pos cog;
  /// \brief Center of geometry before any fitting
  cvm::atom_pos cog_orig;
  /// \brief Unrotated atom positions for fit gradients
  std::vector<cvm::real> atoms_pos_unrotated;
  /// \brief Center of mass
  cvm::atom_pos com;
  /// \brief The derivative of a scalar variable with respect to the COM
  // TODO for scalable calculations of more complex variables (e.g. rotation),
  // use a colvarvalue of vectors to hold the entire derivative
  cvm::rvector scalar_com_gradient;
  /// Dipole moment of the atom group
  cvm::rvector dip;
  /// \brief Lock for modifier
  std::mutex modify_lock;
};

#endif // COLVARATOMS_SOA_H
