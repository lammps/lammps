Granular Sub-model (GSM) styles
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
parent GSM class, several types of sub-model classes are derived:

* GSMNormal: normal force sub-model
* GSMDamping: normal damping sub-model
* GSMTangential: tangential forces and sliding friction sub-model
* GSMRolling: rolling friction sub-model
* GSMTwisting: twisting friction sub-model
* GSMHeat: heat conduction sub-model

For each type of sub-model, more classes are further derived, each describing
a specific implementation. For instance, from the GSMNormal class the
GSMNormalHooke, GSMNormalHertz, and GSMNormalJKR classes are derived which
calculate Hookean, Hertzian, or JKR normal forces, respectively. This modular
structure simplifies the addition of new granular contact models as as one only
needs to create a new GSM class without having to modify the more complex
PairGranular, FixGranWall, and GranularModel classes. Most GSM methods are also
already defined by the parent classes so new contact models typically only require
edits to a few relevant methods (e.g. methods that define coefficients and
calculate forces).

Each GSM class has a pointer to both the LAMMPS class and the GranularModel
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

To create a new GSM class, it is recommended that one first looks at similar GSM
classes. All GSM classes share several general methods which may need to be defined

.. list-table::

   * - ``GSM->mix_coeff()``
     - Optional method to define how coefficients are mixed for different atom types. By default, coefficients are mixed using a geometric mean.
   * - ``GSM->coeffs_to_local()``
     - Parses coefficients to define local variables. Run once at model construction.
   * - ``GSM->init()``
     - Optional method to define local variables after other GSM types were created. For instance, this method may be used by a tangential model that derives parameters from the normal model.

There are also several type-specific methods

.. list-table::

   * - ``GSMNormal->touch()``
     - Optional method to test when particles are in contact. By default, this is when particles overlap.
   * - ``GSMNormal->pulloff_distance()``
     - Optional method to return the distance at which particles stop interacting. By default, this is when particles no longer overlap.
   * - ``GSMNormal->calculate_area()``
     - Optional method to return the surface area of the contact. By default, this returns the geometric cross section.
   * - ``GSMNormal->set_fncrit()``
     - Optional method that defines the critical force to break the contact used by some tangential, rolling, and twisting sub-models. By default, this is the current total normal force including damping.
   * - ``GSMNormal->calculate_forces()``
     - Required method that returns the normal contact force
   * - ``GSMDamping->calculate_forces()``
     - Required method that returns the normal damping force
   * - ``GSMTangential->calculate_forces()``
     - Required method that calculates tangential forces/torques
   * - ``GSMTwisting->calculate_forces()``
     - Required method that calculates twisting friction forces/torques
   * - ``GSMRolling->calculate_forces()``
     - Required method that calculates rolling friction forces/torques
   * - ``GSMHeat->calculate_heat()``
     - Required method that returns the rate of heat flow

As an example, say one wanted to create a new normal force option that consisted
of a Hookean force with a piecewise stiffness. This could be done by adding a new
set of files ``gsm_custom.h``:

.. code-block:: c++

   #ifdef GSM_CLASS
   // clang-format off
   GSMStyle(hooke/piecewise,
            GSMNormalHookePiecewise,
            NORMAL);
   // clang-format on
   #else

   #ifndef GSM_CUSTOM_H_
   #define GSM_CUSTOM_H_

   #include "gsm.h"
   #include "gsm_normal.h"

   namespace LAMMPS_NS {
   namespace Granular_NS {

   class GSMNormalHookePiecewise : public GSMNormal {
    public:
     GSMNormalHookePiecewise(class GranularModel *, class LAMMPS *);
     void coeffs_to_local() override;
     double calculate_forces();
    protected:
     double k1, k2, delta_switch;
   };

   }    // namespace Granular_NS
   }    // namespace LAMMPS_NS

   #endif /*GSM_CUSTOM_H_ */
   #endif /*GSM_CLASS_H_ */


and ``gsm_custom.cpp``

.. code-block:: c++

   #include "gsm_custom.h"
   #include "gsm_normal.h"
   #include "granular_model.h"

   using namespace LAMMPS_NS;
   using namespace Granular_NS;

   GSMNormalHookePiecewise::GSMNormalHookePiecewise(GranularModel *gm, LAMMPS *lmp) :  GSMNormal(gm, lmp)
   {
     num_coeffs = 4;
   }

   /* ---------------------------------------------------------------------- */

   void GSMNormalHookePiecewise::coeffs_to_local()
   {
     k1 = coeffs[0];
     k2 = coeffs[1];
     damp = coeffs[2];
     delta_switch = coeffs[3];
   }

   /* ---------------------------------------------------------------------- */

   double GSMNormalHookePiecewise::calculate_forces()
   {
     double Fne;
     if (gm->delta >= delta_switch) {
       Fne = k1 * delta_switch + k2 * (gm->delta - delta_switch);
     } else {
       Fne = k1 * gm->delta;
     }
     return Fne;
   }

