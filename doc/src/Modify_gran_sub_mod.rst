Granular Sub-Model styles
===============================

In granular models, particles are spheres with a finite radius and rotational
degrees of freedom as further described in the
:doc:`Howto granular page <Howto_granular>`. Interactions between pair of
particles or particles and walls may therefore depend on many different modes
of motion as described in :doc:`pair granular <pair_granular>` and
:doc:`fix wall/gran <fix_wall_gran>`. In both cases, the exchange of forces,
torques, and heat flow between two types of bodies is defined using a
GranularModel class. The GranularModel class organizes the details of an
interaction using a series of granular sub-models each of which describe a
particular interaction mode (e.g. normal forces or rolling friction). From a
parent GranSubMod class, several types of sub-model classes are derived:

* GranSubModNormal: normal force sub-model
* GranSubModDamping: normal damping sub-model
* GranSubModTangential: tangential forces and sliding friction sub-model
* GranSubModRolling: rolling friction sub-model
* GranSubModTwisting: twisting friction sub-model
* GranSubModHeat: heat conduction sub-model

For each type of sub-model, more classes are further derived, each describing a
specific implementation. For instance, from the GranSubModNormal class the
GranSubModNormalHooke, GranSubModNormalHertz, and GranSubModNormalJKR classes
are derived which calculate Hookean, Hertzian, or JKR normal forces,
respectively.  This modular structure simplifies the addition of new granular
contact models as one only needs to create a new GranSubMod class without
having to modify the more complex PairGranular, FixGranWall, and GranularModel
classes. Most GranSubMod methods are also already defined by the parent classes
so new contact models typically only require edits to a few relevant methods
(e.g. methods that define coefficients and calculate forces).

Each GranSubMod class has a pointer to both the LAMMPS class and the GranularModel
class which owns it, ``lmp`` and ``gm``, respectively. The GranularModel class
includes several public variables that describe the geometry/dynamics of the
contact such as

.. list-table::

   * - ``xi`` and ``xj``
     - Positions of the two contacting bodies
   * - ``vi`` and ``vj``
     - Velocities of the two contacting bodies
   * - ``omegai`` and ``omegaj``
     - Angular velocities of the two contacting bodies
   * - ``dx`` and ``nx``
     - The displacement and normalized displacement vectors
   * - ``r``, ``rsq``, and ``rinv``
     - The distance, distance squared, and inverse distance
   * - ``radsum``
     - The sum of particle radii
   * - ``vr``, ``vn``, and ``vt``
     - The relative velocity vector and its normal and tangential components
   * - ``wr``
     - The relative rotational velocity

These quantities, among others, are calculated in the ``GranularModel->check_contact()``
and ``GranularModel->calculate_forces()`` methods which can be referred to for more
details.

To create a new GranSubMod class, it is recommended that one first looks at similar
GranSubMod classes. All GranSubMod classes share several general methods which may
need to be defined

.. list-table::

   * - ``GranSubMod->mix_coeff()``
     - Optional method to define how coefficients are mixed for different atom types. By default, coefficients are mixed using a geometric mean.
   * - ``GranSubMod->coeffs_to_local()``
     - Parses coefficients to define local variables. Run once at model construction.
   * - ``GranSubMod->init()``
     - Optional method to define local variables after other GranSubMod types were created. For instance, this method may be used by a tangential model that derives parameters from the normal model.

There are also several type-specific methods

.. list-table::

   * - ``GranSubModNormal->touch()``
     - Optional method to test when particles are in contact. By default, this is when particles overlap.
   * - ``GranSubModNormal->pulloff_distance()``
     - Optional method to return the distance at which particles stop interacting. By default, this is when particles no longer overlap.
   * - ``GranSubModNormal->calculate_radius()``
     - Optional method to return the radius of the contact. By default, this returns the radius of the geometric cross section.
   * - ``GranSubModNormal->set_fncrit()``
     - Optional method that defines the critical force to break the contact used by some tangential, rolling, and twisting sub-models. By default, this is the current total normal force including damping.
   * - ``GranSubModNormal->calculate_forces()``
     - Required method that returns the normal contact force
   * - ``GranSubModDamping->calculate_forces()``
     - Required method that returns the normal damping force
   * - ``GranSubModTangential->calculate_forces()``
     - Required method that calculates tangential forces/torques
   * - ``GranSubModTwisting->calculate_forces()``
     - Required method that calculates twisting friction forces/torques
   * - ``GranSubModRolling->calculate_forces()``
     - Required method that calculates rolling friction forces/torques
   * - ``GranSubModHeat->calculate_heat()``
     - Required method that returns the rate of heat flow

As an example, say one wanted to create a new normal force option that consisted
of a Hookean force with a piecewise stiffness. This could be done by adding a new
set of files ``gran_sub_mod_custom.h``:

.. code-block:: c++

   #ifdef GranSubMod_CLASS
   // clang-format off
   GranSubModStyle(hooke/piecewise,GranSubModNormalHookePiecewise,NORMAL);
   // clang-format on
   #else

   #ifndef GRAN_SUB_MOD_CUSTOM_H_
   #define GRAN_SUB_MOD_CUSTOM_H_

   #include "gran_sub_mod.h"
   #include "gran_sub_mod_normal.h"

   namespace LAMMPS_NS {
   namespace Granular_NS {
     class GranSubModNormalHookePiecewise : public GranSubModNormal {
      public:
       GranSubModNormalHookePiecewise(class GranularModel *, class LAMMPS *);
       void coeffs_to_local() override;
       double calculate_forces() override;
      protected:
       double k1, k2, delta_switch;
     };

   }    // namespace Granular_NS
   }    // namespace LAMMPS_NS

   #endif /*GRAN_SUB_MOD_CUSTOM_H_ */
   #endif /*GRAN_SUB_MOD_CLASS_H_ */


and ``gran_sub_mod_custom.cpp``

.. code-block:: c++

   #include "gran_sub_mod_custom.h"
   #include "gran_sub_mod_normal.h"
   #include "granular_model.h"

   using namespace LAMMPS_NS;
   using namespace Granular_NS;

   GranSubModNormalHookePiecewise::GranSubModNormalHookePiecewise(GranularModel *gm, LAMMPS *lmp) :
       GranSubModNormal(gm, lmp)
   {
     num_coeffs = 4;
   }

   /* ---------------------------------------------------------------------- */

   void GranSubModNormalHookePiecewise::coeffs_to_local()
   {
     k1 = coeffs[0];
     k2 = coeffs[1];
     damp = coeffs[2];
     delta_switch = coeffs[3];
   }

   /* ---------------------------------------------------------------------- */

   double GranSubModNormalHookePiecewise::calculate_forces()
   {
     double Fne;
     if (gm->delta >= delta_switch) {
       Fne = k1 * delta_switch + k2 * (gm->delta - delta_switch);
     } else {
       Fne = k1 * gm->delta;
     }
     return Fne;
   }

