#ifndef CBPOTENTIAL_H
#define CBPOTENTIAL_H

namespace ATC
{
    //! Enumerated list of interaction types
    enum Interaction{PAIRWISE=1, EAM=2, THREE_BDY=4, ANGLE_BND=8};
    //! Booleans that enable types of terms the potential uses.
    struct Interactions {
      //! Enables up to 3 interaction types.  (order independent)
      Interactions(int a=0, int b=0, int c=0);
      bool pairwise;      //!< Pairwise interaction terms exist.
      bool embedding;     //!< Embedding interaction terms (EAM) exist.
      bool three_body;    //!< Three-body interaction terms exist.
      bool angle_bending; //!< Angle bending interaction terms exist.
    };

  //! Compute bond angle jik from the squared length of vectors ij,ik,kj.
  double calculate_theta(double ij2, double ik2, double jk2);

  /**
   *  @class  CbPotential
   *  @brief  Base class for computing Cauchy-Born quantities for an interatomic potential material
   *          (assume all terms return 0)
   */

  class CbPotential
  {
  protected:
    //! CbPotential base constructor:
    //! Initializes which terms are included in energy computation.
    //@param potential_terms Switches for atomic interaction terms. 
    CbPotential(Interactions interaction_terms) : terms(interaction_terms) {}
  public:
    virtual ~CbPotential() {}   
    const Interactions terms;  //!< switches for types of potential terms.

    //! Returns the minimum distance that all interactions get neglected.
    virtual double cutoff_radius() const=0;

    //! @name Pairwise interaction term and derivatives.
    //@{
    virtual double phi    (const double & /* r */) const { return 0.0; }
    virtual double phi_r  (const double &r) const;
    virtual double phi_rr (const double &r) const;
    virtual double phi_rrr(const double &r) const;
    //@}

    //! @name Embedding terms. Electron cloud density and embedding functions
    //@{
    virtual double rho   (const double & /* r */) const { return 0.0; }
    virtual double rho_r (const double &r) const;
    virtual double rho_rr(const double &r) const;
    virtual double rho_rrr(const double &r) const;
    virtual double F   (const double & /* p */) const { return 0.0; }
    virtual double F_p (const double &p) const;
    virtual double F_pp(const double &p) const;
    virtual double F_ppp(const double &p) const;
    //@}

    //! @name Three-body terms and derivatives
    //@{
    virtual double phi3   (const double & /* q */) const {return 0.0; }
    virtual double phi3_q (const double &q) const;
    virtual double phi3_qq(const double &q) const;
    //@}

  };
}
#endif
