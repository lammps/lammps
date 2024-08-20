// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvar_arithmeticpath.h"


struct ArithmeticPathImpl: public ArithmeticPathCV::ArithmeticPathBase<cvm::real> {
    std::vector<std::vector<colvarvalue>> frame_element_distances;
    std::vector<std::vector<colvarvalue>> dsdx;
    std::vector<std::vector<colvarvalue>> dzdx;
    template <typename T>
    void updateCartesianDistanceToReferenceFrames(T* obj) {
        for (size_t i_frame = 0; i_frame < obj->reference_frames.size(); ++i_frame) {
            for (size_t i_atom = 0; i_atom < obj->atoms->size(); ++i_atom) {
                frame_element_distances[i_frame][i_atom] = (*(obj->comp_atoms[i_frame]))[i_atom].pos - obj->reference_frames[i_frame][i_atom];
            }
        }
    }
    template <typename T>
    void updateCVDistanceToReferenceFrames(T* obj) {
        for (size_t i_cv = 0; i_cv < obj->cv.size(); ++i_cv) {
            obj->cv[i_cv]->calc_value();
        }
        for (size_t i_frame = 0; i_frame < obj->ref_cv.size(); ++i_frame) {
            for (size_t i_cv = 0; i_cv < obj->cv.size(); ++i_cv) {
                colvarvalue ref_cv_value(obj->ref_cv[i_frame][i_cv]);
                colvarvalue current_cv_value(obj->cv[i_cv]->value());
                if (current_cv_value.type() == colvarvalue::type_scalar) {
                    frame_element_distances[i_frame][i_cv] = 0.5 * obj->cv[i_cv]->dist2_lgrad(obj->cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, obj->cv[i_cv]->sup_np)), ref_cv_value.real_value);
                } else {
                    frame_element_distances[i_frame][i_cv] = 0.5 * obj->cv[i_cv]->dist2_lgrad(obj->cv[i_cv]->sup_coeff * current_cv_value, ref_cv_value);
                }
            }
        }
    }
    ArithmeticPathImpl(size_t p_num_elements, size_t p_total_frames, cvm::real p_lambda, const std::vector<cvm::real>& p_weights) {
        ArithmeticPathCV::ArithmeticPathBase<cvm::real>::initialize(p_num_elements, p_total_frames, p_lambda, p_weights);
        frame_element_distances.resize(p_total_frames, std::vector<colvarvalue>(p_num_elements, colvarvalue(colvarvalue::Type::type_notset)));
        dsdx.resize(p_total_frames, std::vector<colvarvalue>(p_num_elements, colvarvalue(colvarvalue::Type::type_notset)));
        dzdx.resize(p_total_frames, std::vector<colvarvalue>(p_num_elements, colvarvalue(colvarvalue::Type::type_notset)));
    }
    cvm::real get_lambda() const {return lambda;}
    cvm::real compute_s() {
        cvm::real s;
        computeValue(frame_element_distances, &s, nullptr);
        return s;
    }
    cvm::real compute_z() {
        cvm::real z;
        computeValue(frame_element_distances, nullptr, &z);
        return z;
    }
    void compute_s_derivatives() {
        computeDerivatives<colvarvalue>(frame_element_distances, &dsdx, nullptr);
    }
    void compute_z_derivatives() {
        computeDerivatives<colvarvalue>(frame_element_distances, nullptr, &dzdx);
    }
    // for debug gradients of implicit sub CVs
    template <typename T>
    colvarvalue compute_s_analytical_derivative_ij(size_t i, size_t j, cvm::real eps, T* obj) const {
        ArithmeticPathImpl tmp_left(*this), tmp_right(*this);
        const size_t value_size = frame_element_distances[i][j].size();
        colvarvalue result(frame_element_distances[i][j].type());
        colvarvalue ref_cv_value(obj->ref_cv[i][j]);
        for (size_t k = 0; k < value_size; ++k) {
            // get the current CV value
            colvarvalue current_cv_value(obj->cv[j]->value());
            // save the old values in frame element distance matrices
            const auto saved_left = tmp_left.frame_element_distances[i][j][k];
            const auto saved_right = tmp_right.frame_element_distances[i][j][k];
            // update frame element distance matrices
            if (current_cv_value.type() == colvarvalue::type_scalar) {
                tmp_left.frame_element_distances[i][j] = 0.5 * obj->cv[j]->dist2_lgrad(obj->cv[j]->sup_coeff * (cvm::pow(current_cv_value.real_value - eps, obj->cv[j]->sup_np)), ref_cv_value.real_value);
                tmp_right.frame_element_distances[i][j] = 0.5 * obj->cv[j]->dist2_lgrad(obj->cv[j]->sup_coeff * (cvm::pow(current_cv_value.real_value + eps, obj->cv[j]->sup_np)), ref_cv_value.real_value);
            } else {
                current_cv_value[k] -= eps;
                tmp_left.frame_element_distances[i][j] = 0.5 * obj->cv[j]->dist2_lgrad(obj->cv[j]->sup_coeff * current_cv_value, ref_cv_value);
                current_cv_value[k] += eps + eps;
                tmp_right.frame_element_distances[i][j] = 0.5 * obj->cv[j]->dist2_lgrad(obj->cv[j]->sup_coeff * current_cv_value, ref_cv_value);
            }
            const cvm::real s_left = tmp_left.compute_s();
            const cvm::real s_right = tmp_right.compute_s();
            result[k] = (s_right - s_left) / (2.0 * eps);
            tmp_left.frame_element_distances[i][j][k] = saved_left;
            tmp_right.frame_element_distances[i][j][k] = saved_right;
        }
        return result;
    }
    template <typename T>
    colvarvalue compute_z_analytical_derivative_ij(size_t i, size_t j, cvm::real eps, T* obj) const {
        ArithmeticPathImpl tmp_left(*this), tmp_right(*this);
        const size_t value_size = frame_element_distances[i][j].size();
        colvarvalue result(frame_element_distances[i][j].type());
        colvarvalue ref_cv_value(obj->ref_cv[i][j]);
        for (size_t k = 0; k < value_size; ++k) {
            // get the current CV value
            colvarvalue current_cv_value(obj->cv[j]->value());
            // save the old values in frame element distance matrices
            const auto saved_left = tmp_left.frame_element_distances[i][j][k];
            const auto saved_right = tmp_right.frame_element_distances[i][j][k];
            // update frame element distance matrices
            if (current_cv_value.type() == colvarvalue::type_scalar) {
                tmp_left.frame_element_distances[i][j] = 0.5 * obj->cv[j]->dist2_lgrad(obj->cv[j]->sup_coeff * (cvm::pow(current_cv_value.real_value - eps, obj->cv[j]->sup_np)), ref_cv_value.real_value);
                tmp_right.frame_element_distances[i][j] = 0.5 * obj->cv[j]->dist2_lgrad(obj->cv[j]->sup_coeff * (cvm::pow(current_cv_value.real_value + eps, obj->cv[j]->sup_np)), ref_cv_value.real_value);
            } else {
                current_cv_value[k] -= eps;
                tmp_left.frame_element_distances[i][j] = 0.5 * obj->cv[j]->dist2_lgrad(obj->cv[j]->sup_coeff * current_cv_value, ref_cv_value);
                current_cv_value[k] += eps + eps;
                tmp_right.frame_element_distances[i][j] = 0.5 * obj->cv[j]->dist2_lgrad(obj->cv[j]->sup_coeff * current_cv_value, ref_cv_value);
            }
            const cvm::real z_left = tmp_left.compute_z();
            const cvm::real z_right = tmp_right.compute_z();
            result[k] = (z_right - z_left) / (2.0 * eps);
            tmp_left.frame_element_distances[i][j][k] = saved_left;
            tmp_right.frame_element_distances[i][j][k] = saved_right;
        }
        return result;
    }
};

colvar::aspath::aspath()
{
    set_function_type("aspath");
    x.type(colvarvalue::type_scalar);
}


int colvar::aspath::init(std::string const &conf)
{
    int error_code = CartesianBasedPath::init(conf);
    if (error_code != COLVARS_OK) return error_code;
    cvm::log(std::string("Total number of frames: ") + cvm::to_str(total_reference_frames) + std::string("\n"));
    cvm::real p_lambda;
    get_keyval(conf, "lambda", p_lambda, -1.0);
    const size_t num_atoms = atoms->size();
    std::vector<cvm::real> p_weights(num_atoms, std::sqrt(1.0 / num_atoms));
    // ArithmeticPathCV::ArithmeticPathBase<cvm::atom_pos, cvm::real, ArithmeticPathCV::path_sz::S>::initialize(num_atoms, total_reference_frames, p_lambda, reference_frames[0], p_weights);
    if (impl_) impl_.reset();
    impl_ = std::unique_ptr<ArithmeticPathImpl>(new ArithmeticPathImpl(num_atoms, total_reference_frames, p_lambda, p_weights));
    cvm::log(std::string("Lambda is ") + cvm::to_str(impl_->get_lambda()) + std::string("\n"));
    return error_code;
}

colvar::aspath::~aspath() {}

void colvar::aspath::calc_value() {
    if (impl_->get_lambda() < 0) {
        // this implies that the user may not set a valid lambda value
        // so recompute it by the suggested value in Parrinello's paper
        cvm::log("A non-positive value of lambda is detected, which implies that it may not set in the configuration.\n");
        cvm::log("This component (aspath) will recompute a value for lambda following the suggestion in the origin paper.\n");
        std::vector<cvm::real> rmsd_between_refs(total_reference_frames - 1, 0.0);
        computeDistanceBetweenReferenceFrames(rmsd_between_refs);
        impl_->reComputeLambda(rmsd_between_refs);
        cvm::log("Ok, the value of lambda is updated to " + cvm::to_str(impl_->get_lambda()));
    }
    impl_->updateCartesianDistanceToReferenceFrames(this);
    x = impl_->compute_s();
}

void colvar::aspath::calc_gradients() {
    impl_->compute_s_derivatives();
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            (*(comp_atoms[i_frame]))[i_atom].grad += impl_->dsdx[i_frame][i_atom];
        }
    }
}

void colvar::aspath::apply_force(colvarvalue const &force) {
    cvm::real const &F = force.real_value;
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        (*(comp_atoms[i_frame])).apply_colvar_force(F);
    }
}

colvar::azpath::azpath()
{
    set_function_type("azpath");
    x.type(colvarvalue::type_scalar);
}

int colvar::azpath::init(std::string const &conf)
{
    int error_code = CartesianBasedPath::init(conf);
    if (error_code != COLVARS_OK) return error_code;
    cvm::log(std::string("Total number of frames: ") + cvm::to_str(total_reference_frames) + std::string("\n"));
    x.type(colvarvalue::type_scalar);
    cvm::real p_lambda;
    get_keyval(conf, "lambda", p_lambda, -1.0);
    const size_t num_atoms = atoms->size();
    std::vector<cvm::real> p_weights(num_atoms, std::sqrt(1.0 / num_atoms));
    if (impl_) impl_.reset();
    impl_ = std::unique_ptr<ArithmeticPathImpl>(new ArithmeticPathImpl(num_atoms, total_reference_frames, p_lambda, p_weights));
    cvm::log(std::string("Lambda is ") + cvm::to_str(impl_->get_lambda()) + std::string("\n"));
    return error_code;
}

colvar::azpath::~azpath() {}

void colvar::azpath::calc_value() {
    if (impl_->get_lambda() < 0) {
        // this implies that the user may not set a valid lambda value
        // so recompute it by the suggested value in Parrinello's paper
        cvm::log("A non-positive value of lambda is detected, which implies that it may not set in the configuration.\n");
        cvm::log("This component (azpath) will recompute a value for lambda following the suggestion in the origin paper.\n");
        std::vector<cvm::real> rmsd_between_refs(total_reference_frames - 1, 0.0);
        computeDistanceBetweenReferenceFrames(rmsd_between_refs);
        impl_->reComputeLambda(rmsd_between_refs);
        cvm::log("Ok, the value of lambda is updated to " + cvm::to_str(impl_->get_lambda()));
    }
    impl_->updateCartesianDistanceToReferenceFrames(this);
    x = impl_->compute_z();
}

void colvar::azpath::calc_gradients() {
    impl_->compute_z_derivatives();
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            (*(comp_atoms[i_frame]))[i_atom].grad += impl_->dzdx[i_frame][i_atom];
        }
    }
}

void colvar::azpath::apply_force(colvarvalue const &force) {
    cvm::real const &F = force.real_value;
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        (*(comp_atoms[i_frame])).apply_colvar_force(F);
    }
}

colvar::aspathCV::aspathCV()
{
    set_function_type("aspathCV");
    x.type(colvarvalue::type_scalar);
}

int colvar::aspathCV::init(std::string const &conf)
{
    int error_code = CVBasedPath::init(conf);
    if (error_code != COLVARS_OK) return error_code;
    cvm::log(std::string("Total number of frames: ") + cvm::to_str(total_reference_frames) + std::string("\n"));
    std::vector<cvm::real> p_weights(cv.size(), 1.0);
    get_keyval(conf, "weights", p_weights, std::vector<cvm::real>(cv.size(), 1.0));
    use_explicit_gradients = true;
    cvm::real p_lambda;
    get_keyval(conf, "lambda", p_lambda, -1.0);
    if (impl_) impl_.reset();
    impl_ = std::unique_ptr<ArithmeticPathImpl>(new ArithmeticPathImpl(cv.size(), total_reference_frames, p_lambda, p_weights));
    cvm::log(std::string("Lambda is ") + cvm::to_str(impl_->get_lambda()) + std::string("\n"));
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
        cvm::log(std::string("The weight of CV ") + cvm::to_str(i_cv) + std::string(" is ") + cvm::to_str(p_weights[i_cv]) + std::string("\n"));
    }
    return error_code;
}

colvar::aspathCV::~aspathCV() {}

void colvar::aspathCV::calc_value() {
    if (impl_->get_lambda() < 0) {
        // this implies that the user may not set a valid lambda value
        // so recompute it by the suggested value in Parrinello's paper
        cvm::log("A non-positive value of lambda is detected, which implies that it may not set in the configuration.\n");
        cvm::log("This component (aspathCV) will recompute a value for lambda following the suggestion in the origin paper.\n");
        std::vector<cvm::real> rmsd_between_refs(total_reference_frames - 1, 0.0);
        computeDistanceBetweenReferenceFrames(rmsd_between_refs);
        impl_->reComputeLambda(rmsd_between_refs);
        cvm::log("Ok, the value of lambda is updated to " + cvm::to_str(impl_->get_lambda()));
    }
    impl_->updateCVDistanceToReferenceFrames(this);
    x = impl_->compute_s();
}

void colvar::aspathCV::calc_gradients() {
    impl_->compute_s_derivatives();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_gradients();
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            // compute the gradient (grad) with respect to the i-th CV
            colvarvalue grad(cv[i_cv]->value().type());
            // sum up derivatives with respect to all frames
            for (size_t m_frame = 0; m_frame < impl_->dsdx.size(); ++m_frame) {
                // dsdx is the derivative of s with respect to the m-th frame
                grad += impl_->dsdx[m_frame][i_cv];
            }
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = grad[j_elem] * factor_polynomial * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                    }
                }
            }
        }
    }
}

void colvar::aspathCV::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            // compute the gradient (grad) with respect to the i-th CV
            colvarvalue grad(cv[i_cv]->value().type());
            for (size_t m_frame = 0; m_frame < impl_->dsdx.size(); ++m_frame) {
                // dsdx is the derivative of s with respect to the m-th frame
                grad += impl_->dsdx[m_frame][i_cv];
            }
            grad *= factor_polynomial;
            cv[i_cv]->apply_force(force.real_value * grad);
            // try my best to debug gradients even if the sub-CVs do not have explicit gradients
            if (is_enabled(f_cvc_debug_gradient)) {
                cvm::log("Debugging gradients for " + description +
                         " with respect to sub-CV " + cv[i_cv]->description +
                         ", which has no explicit gradient with respect to its own input(s)");
                colvarvalue analytical_grad(cv[i_cv]->value().type());
                for (size_t m_frame = 0; m_frame < impl_->dsdx.size(); ++m_frame) {
                    analytical_grad += impl_->compute_s_analytical_derivative_ij(
                        m_frame, i_cv, cvm::debug_gradients_step_size, this);
                }
                cvm::log("dx(actual) = "+cvm::to_str(analytical_grad, 21, 14)+"\n");
                cvm::log("dx(interp) = "+cvm::to_str(grad, 21, 14)+"\n");
                cvm::log("|dx(actual) - dx(interp)|/|dx(actual)| = "+
                  cvm::to_str((analytical_grad - grad).norm() /
                              (analytical_grad).norm(), 12, 5)+"\n");
            }
        }
    }
}

colvar::azpathCV::azpathCV()
{
    set_function_type("azpathCV");
    x.type(colvarvalue::type_scalar);
}

int colvar::azpathCV::init(std::string const &conf)
{
    int error_code = CVBasedPath::init(conf);
    if (error_code != COLVARS_OK) return error_code;
    cvm::log(std::string("Total number of frames: ") + cvm::to_str(total_reference_frames) + std::string("\n"));
    std::vector<cvm::real> p_weights(cv.size(), 1.0);
    get_keyval(conf, "weights", p_weights, std::vector<cvm::real>(cv.size(), 1.0));
    use_explicit_gradients = true;
    cvm::real p_lambda;
    get_keyval(conf, "lambda", p_lambda, -1.0);
    if (impl_) impl_.reset();
    impl_ = std::unique_ptr<ArithmeticPathImpl>(new ArithmeticPathImpl(cv.size(), total_reference_frames, p_lambda, p_weights));
    cvm::log(std::string("Lambda is ") + cvm::to_str(impl_->get_lambda()) + std::string("\n"));
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
        cvm::log(std::string("The weight of CV ") + cvm::to_str(i_cv) + std::string(" is ") + cvm::to_str(p_weights[i_cv]) + std::string("\n"));
    }
    return error_code;
}

void colvar::azpathCV::calc_value() {
    if (impl_->get_lambda() < 0) {
        // this implies that the user may not set a valid lambda value
        // so recompute it by the suggested value in Parrinello's paper
        cvm::log("A non-positive value of lambda is detected, which implies that it may not set in the configuration.\n");
        cvm::log("This component (azpathCV) will recompute a value for lambda following the suggestion in the origin paper.\n");
        std::vector<cvm::real> rmsd_between_refs(total_reference_frames - 1, 0.0);
        computeDistanceBetweenReferenceFrames(rmsd_between_refs);
        impl_->reComputeLambda(rmsd_between_refs);
        cvm::log("Ok, the value of lambda is updated to " + cvm::to_str(impl_->get_lambda()));
    }
    impl_->updateCVDistanceToReferenceFrames(this);
    x = impl_->compute_z();
}

void colvar::azpathCV::calc_gradients() {
    impl_->compute_z_derivatives();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_gradients();
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            // compute the gradient (grad) with respect to the i-th CV
            colvarvalue grad(cv[i_cv]->value().type());
            // sum up derivatives with respect to all frames
            for (size_t m_frame = 0; m_frame < impl_->dzdx.size(); ++m_frame) {
                // dzdx is the derivative of z with respect to the m-th frame
                grad += impl_->dzdx[m_frame][i_cv];
            }
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = grad[j_elem] * factor_polynomial * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                    }
                }
            }
        }
    }
}

void colvar::azpathCV::apply_force(colvarvalue const &force) {
    // the PCV component itself is a scalar, so force should be scalar
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            // compute the gradient (grad) with respect to the i-th CV
            colvarvalue grad(cv[i_cv]->value().type());
            for (size_t m_frame = 0; m_frame < impl_->dzdx.size(); ++m_frame) {
                // dzdx is the derivative of z with respect to the m-th frame
                grad += impl_->dzdx[m_frame][i_cv];
            }
            grad *= factor_polynomial;
            cv[i_cv]->apply_force(force.real_value * grad);
            // try my best to debug gradients even if the sub-CVs do not have explicit gradients
            if (is_enabled(f_cvc_debug_gradient)) {
                cvm::log("Debugging gradients for " + description +
                         " with respect to sub-CV " + cv[i_cv]->description +
                         ", which has no explicit gradient with respect to its own input(s)");
                colvarvalue analytical_grad(cv[i_cv]->value().type());
                for (size_t m_frame = 0; m_frame < impl_->dzdx.size(); ++m_frame) {
                    analytical_grad += impl_->compute_z_analytical_derivative_ij(
                        m_frame, i_cv, cvm::debug_gradients_step_size, this);
                }
                cvm::log("dx(actual) = "+cvm::to_str(analytical_grad, 21, 14)+"\n");
                cvm::log("dx(interp) = "+cvm::to_str(grad, 21, 14)+"\n");
                cvm::log("|dx(actual) - dx(interp)|/|dx(actual)| = "+
                  cvm::to_str((analytical_grad - grad).norm() /
                              (analytical_grad).norm(), 12, 5)+"\n");
            }
        }
    }
}

colvar::azpathCV::~azpathCV() {}

