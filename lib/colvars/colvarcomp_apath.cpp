#if (__cplusplus >= 201103L)

#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"

colvar::aspathCV::aspathCV(std::string const &conf): CVBasedPath(conf) {
    function_type = "aspathCV";
    cvm::log(std::string("Total number of frames: ") + cvm::to_str(total_reference_frames) + std::string("\n"));
    std::vector<cvm::real> p_weights(cv.size(), 1.0);
    get_keyval(conf, "weights", p_weights, std::vector<cvm::real>(cv.size(), 1.0));
    x.type(colvarvalue::type_scalar);
    use_explicit_gradients = true;
    cvm::real p_lambda;
    get_keyval(conf, "lambda", p_lambda, -1.0);
    ArithmeticPathCV::ArithmeticPathBase<colvarvalue, cvm::real, ArithmeticPathCV::path_sz::S>::initialize(cv.size(), total_reference_frames, p_lambda, ref_cv[0], p_weights);
    cvm::log(std::string("Lambda is ") + cvm::to_str(lambda) + std::string("\n"));
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
        cvm::log(std::string("The weight of CV ") + cvm::to_str(i_cv) + std::string(" is ") + cvm::to_str(weights[i_cv]) + std::string("\n"));
    }
}

void colvar::aspathCV::updateDistanceToReferenceFrames() {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_value();
    }
    for (size_t i_frame = 0; i_frame < ref_cv.size(); ++i_frame) {
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            colvarvalue ref_cv_value(ref_cv[i_frame][i_cv]);
            colvarvalue current_cv_value(cv[i_cv]->value());
            if (current_cv_value.type() == colvarvalue::type_scalar) {
                frame_element_distances[i_frame][i_cv] = 0.5 * cv[i_cv]->dist2_lgrad(cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np)), ref_cv_value.real_value);
            } else {
                frame_element_distances[i_frame][i_cv] = 0.5 * cv[i_cv]->dist2_lgrad(cv[i_cv]->sup_coeff * current_cv_value, ref_cv_value);
            }
        }
    }
}

void colvar::aspathCV::calc_value() {
    if (lambda < 0) {
        // this implies that the user may not set a valid lambda value
        // so recompute it by the suggested value in Parrinello's paper
        cvm::log("A non-positive value of lambda is detected, which implies that it may not set in the configuration.\n");
        cvm::log("This component (aspathCV) will recompute a value for lambda following the suggestion in the origin paper.\n");
        std::vector<cvm::real> rmsd_between_refs(total_reference_frames - 1, 0.0);
        computeDistanceBetweenReferenceFrames(rmsd_between_refs);
        reComputeLambda(rmsd_between_refs);
        cvm::log("Ok, the value of lambda is updated to " + cvm::to_str(lambda));
    }
    computeValue();
    x = s;
}

void colvar::aspathCV::calc_gradients() {
    computeDerivatives();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_gradients();
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)) {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = dsdx[i_cv][j_elem] * factor_polynomial * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                    }
                }
            }
        }
    }
}

void colvar::aspathCV::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)
        ) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            colvarvalue cv_force = dsdx[i_cv] * force.real_value * factor_polynomial;
            cv[i_cv]->apply_force(cv_force);
        }
    }
}

colvar::aspathCV::~aspathCV() {}

colvar::azpathCV::azpathCV(std::string const &conf): CVBasedPath(conf) {
    function_type = "azpathCV";
    cvm::log(std::string("Total number of frames: ") + cvm::to_str(total_reference_frames) + std::string("\n"));
    std::vector<cvm::real> p_weights(cv.size(), 1.0);
    get_keyval(conf, "weights", p_weights, std::vector<cvm::real>(cv.size(), 1.0));
    x.type(colvarvalue::type_scalar);
    use_explicit_gradients = true;
    cvm::real p_lambda;
    get_keyval(conf, "lambda", p_lambda, -1.0);
    ArithmeticPathCV::ArithmeticPathBase<colvarvalue, cvm::real, ArithmeticPathCV::path_sz::Z>::initialize(cv.size(), total_reference_frames, p_lambda, ref_cv[0], p_weights);
    cvm::log(std::string("Lambda is ") + cvm::to_str(lambda) + std::string("\n"));
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
        cvm::log(std::string("The weight of CV ") + cvm::to_str(i_cv) + std::string(" is ") + cvm::to_str(weights[i_cv]) + std::string("\n"));
    }
}

void colvar::azpathCV::updateDistanceToReferenceFrames() {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_value();
    }
    for (size_t i_frame = 0; i_frame < ref_cv.size(); ++i_frame) {
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            colvarvalue ref_cv_value(ref_cv[i_frame][i_cv]);
            colvarvalue current_cv_value(cv[i_cv]->value());
            if (current_cv_value.type() == colvarvalue::type_scalar) {
                frame_element_distances[i_frame][i_cv] = 0.5 * cv[i_cv]->dist2_lgrad(cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np)), ref_cv_value.real_value);
            } else {
                frame_element_distances[i_frame][i_cv] = 0.5 * cv[i_cv]->dist2_lgrad(cv[i_cv]->sup_coeff * current_cv_value, ref_cv_value);
            }
        }
    }
}

void colvar::azpathCV::calc_value() {
    if (lambda < 0) {
        // this implies that the user may not set a valid lambda value
        // so recompute it by the suggested value in Parrinello's paper
        cvm::log("A non-positive value of lambda is detected, which implies that it may not set in the configuration.\n");
        cvm::log("This component (azpathCV) will recompute a value for lambda following the suggestion in the origin paper.\n");
        std::vector<cvm::real> rmsd_between_refs(total_reference_frames - 1, 0.0);
        computeDistanceBetweenReferenceFrames(rmsd_between_refs);
        reComputeLambda(rmsd_between_refs);
        cvm::log("Ok, the value of lambda is updated to " + cvm::to_str(lambda));
    }
    computeValue();
    x = z;
}

void colvar::azpathCV::calc_gradients() {
    computeDerivatives();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_gradients();
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)) {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = dzdx[i_cv][j_elem] * factor_polynomial * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                    }
                }
            }

        }
    }
}

void colvar::azpathCV::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)
        ) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            const colvarvalue cv_force = dzdx[i_cv] * force.real_value * factor_polynomial;
            cv[i_cv]->apply_force(cv_force);
        }
    }
}

colvar::azpathCV::~azpathCV() {}

#endif
