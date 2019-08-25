#ifndef QUADRATURE_H
#define QUADRATURE_H

namespace ATC {
/** 
  *  @class Quadrature 
  *  @brief create quadrature lists  
*/
class Quadrature {
  public:
  /** Static instance of this class */
  static Quadrature * instance();

  /** Destroy */
  static void Destroy();

  /** domain of integration is -1 to 1 */
  void set_line_quadrature(const int ng, double* xg, double* wg);

  protected:
    Quadrature();
  private:
    static Quadrature * myInstance_;
};
}
#endif
