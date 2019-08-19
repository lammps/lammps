#if (__cplusplus >= 201103L)

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

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

namespace GeometricPathCV {
void init_string_cv_map(std::map<std::string, std::function<colvar::cvc* (const std::string& conf)>>& string_cv_map);
}

bool compareColvarComponent(colvar::cvc *i, colvar::cvc *j)
{
    return i->name < j->name;
}

colvar::CartesianBasedPath::CartesianBasedPath(std::string const &conf): cvc(conf), atoms(nullptr), reference_frames(0) {
    // Parse selected atoms
    atoms = parse_group(conf, "atoms");
    has_user_defined_fitting = false;
    std::string fitting_conf;
    if (key_lookup(conf, "fittingAtoms", &fitting_conf)) {
        has_user_defined_fitting = true;
    }
    // Lookup reference column of PDB
    // Copied from the RMSD class
    std::string reference_column;
    double reference_column_value;
    if (get_keyval(conf, "refPositionsCol", reference_column, std::string(""))) {
        bool found = get_keyval(conf, "refPositionsColValue", reference_column_value, 0.0);
        if (found && reference_column_value == 0.0) {
          cvm::error("Error: refPositionsColValue, "
                     "if provided, must be non-zero.\n");
          return;
        }
    }
    // Lookup all reference frames
    bool has_frames = true;
    total_reference_frames = 0;
    while (has_frames) {
        std::string reference_position_file_lookup = "refPositionsFile" + cvm::to_str(total_reference_frames + 1);
        if (key_lookup(conf, reference_position_file_lookup.c_str())) {
            std::string reference_position_filename;
            get_keyval(conf, reference_position_file_lookup.c_str(), reference_position_filename, std::string(""));
            std::vector<cvm::atom_pos> reference_position(atoms->size());
            cvm::load_coords(reference_position_filename.c_str(), &reference_position, atoms, reference_column, reference_column_value);
            reference_frames.push_back(reference_position);
            ++total_reference_frames;
        } else {
            has_frames = false;
        }
    }
    // Setup alignment to compute RMSD with respect to reference frames
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        cvm::atom_group* tmp_atoms = parse_group(conf, "atoms");
        if (!has_user_defined_fitting) {
            // Swipe from the rmsd class
            tmp_atoms->b_center = true;
            tmp_atoms->b_rotate = true;
            tmp_atoms->ref_pos = reference_frames[i_frame];
            tmp_atoms->center_ref_pos();
            tmp_atoms->enable(f_ag_fit_gradients);
            tmp_atoms->rot.request_group1_gradients(tmp_atoms->size());
            tmp_atoms->rot.request_group2_gradients(tmp_atoms->size());
            comp_atoms.push_back(tmp_atoms);
        } else {
            // parse a group of atoms for fitting
            std::string fitting_group_name = std::string("fittingAtoms") + cvm::to_str(i_frame);
            cvm::atom_group* tmp_fitting_atoms = new cvm::atom_group(fitting_group_name.c_str());
            tmp_fitting_atoms->parse(fitting_conf);
            tmp_fitting_atoms->disable(f_ag_scalable);
            tmp_fitting_atoms->disable(f_ag_scalable_com);
            tmp_fitting_atoms->fit_gradients.assign(tmp_fitting_atoms->size(), cvm::atom_pos(0.0, 0.0, 0.0));
            std::string reference_position_file_lookup = "refPositionsFile" + cvm::to_str(i_frame + 1);
            std::string reference_position_filename;
            get_keyval(conf, reference_position_file_lookup.c_str(), reference_position_filename, std::string(""));
            std::vector<cvm::atom_pos> reference_fitting_position(tmp_fitting_atoms->size());
            cvm::load_coords(reference_position_filename.c_str(), &reference_fitting_position, tmp_fitting_atoms, reference_column, reference_column_value);
            // setup the atom group for calculating
            tmp_atoms->b_center = true;
            tmp_atoms->b_rotate = true;
            tmp_atoms->b_user_defined_fit = true;
            tmp_atoms->disable(f_ag_scalable);
            tmp_atoms->disable(f_ag_scalable_com);
            tmp_atoms->ref_pos = reference_fitting_position;
            tmp_atoms->center_ref_pos();
            tmp_atoms->enable(f_ag_fit_gradients);
            tmp_atoms->enable(f_ag_fitting_group);
            tmp_atoms->fitting_group = tmp_fitting_atoms;
            tmp_atoms->rot.request_group1_gradients(tmp_fitting_atoms->size());
            tmp_atoms->rot.request_group2_gradients(tmp_fitting_atoms->size());
            reference_fitting_frames.push_back(reference_fitting_position);
            comp_atoms.push_back(tmp_atoms);
        }
    }
    x.type(colvarvalue::type_scalar);
    // Don't use implicit gradient
    enable(f_cvc_explicit_gradient);
}

colvar::CartesianBasedPath::~CartesianBasedPath() {
    for (auto it_comp_atoms = comp_atoms.begin(); it_comp_atoms != comp_atoms.end(); ++it_comp_atoms) {
        if (*it_comp_atoms != nullptr) {
            delete (*it_comp_atoms);
            (*it_comp_atoms) = nullptr;
        }
    }
}

void colvar::CartesianBasedPath::computeReferenceDistance(std::vector<cvm::real>& result) {
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        cvm::real frame_rmsd = 0.0;
        for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            frame_rmsd += ((*(comp_atoms[i_frame]))[i_atom].pos - reference_frames[i_frame][i_atom]).norm2();
        }
        frame_rmsd /= cvm::real(atoms->size());
        frame_rmsd = cvm::sqrt(frame_rmsd);
        result[i_frame] = frame_rmsd;
    }
}

colvar::gspath::gspath(std::string const &conf): CartesianBasedPath(conf) {
    function_type = "gspath";
    get_keyval(conf, "useSecondClosestFrame", use_second_closest_frame, true);
    if (use_second_closest_frame == true) {
        cvm::log(std::string("Geometric path s(σ) will use the second closest frame to compute s_(m-1)\n"));
    } else {
        cvm::log(std::string("Geometric path s(σ) will use the neighbouring frame to compute s_(m-1)\n"));
    }
    get_keyval(conf, "useThirdClosestFrame", use_third_closest_frame, false);
    if (use_third_closest_frame == true) {
        cvm::log(std::string("Geometric path s(σ) will use the third closest frame to compute s_(m+1)\n"));
    } else {
        cvm::log(std::string("Geometric path s(σ) will use the neighbouring frame to compute s_(m+1)\n"));
    }
    GeometricPathCV::GeometricPathBase<cvm::atom_pos, cvm::real, GeometricPathCV::path_sz::S>::initialize(atoms->size(), cvm::atom_pos(), total_reference_frames, use_second_closest_frame, use_third_closest_frame);
    cvm::log(std::string("Geometric pathCV(s) is initialized.\n"));
    cvm::log(std::string("Geometric pathCV(s) loaded ") + cvm::to_str(reference_frames.size()) + std::string(" frames.\n"));
}

void colvar::gspath::updateReferenceDistances() {
    computeReferenceDistance(frame_distances);
}

void colvar::gspath::prepareVectors() {
    size_t i_atom;
    for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        // v1 = s_m - z
        v1[i_atom] = reference_frames[min_frame_index_1][i_atom] - (*(comp_atoms[min_frame_index_1]))[i_atom].pos;
        // v2 = z - s_(m-1)
        v2[i_atom] = (*(comp_atoms[min_frame_index_2]))[i_atom].pos - reference_frames[min_frame_index_2][i_atom];
    }
    if (min_frame_index_3 < 0 || min_frame_index_3 > M) {
        cvm::atom_pos reference_cog_1, reference_cog_2;
        for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            reference_cog_1 += reference_frames[min_frame_index_1][i_atom];
            reference_cog_2 += reference_frames[min_frame_index_2][i_atom];
        }
        reference_cog_1 /= reference_frames[min_frame_index_1].size();
        reference_cog_2 /= reference_frames[min_frame_index_2].size();
        std::vector<cvm::atom_pos> tmp_reference_frame_1(reference_frames[min_frame_index_1].size());
        std::vector<cvm::atom_pos> tmp_reference_frame_2(reference_frames[min_frame_index_2].size());
        for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            tmp_reference_frame_1[i_atom] = reference_frames[min_frame_index_1][i_atom] - reference_cog_1;
            tmp_reference_frame_2[i_atom] = reference_frames[min_frame_index_2][i_atom] - reference_cog_2;
        }
        if (has_user_defined_fitting) {
            cvm::atom_pos reference_fitting_cog_1, reference_fitting_cog_2;
            for (i_atom = 0; i_atom < reference_fitting_frames[min_frame_index_1].size(); ++i_atom) {
                reference_fitting_cog_1 += reference_fitting_frames[min_frame_index_1][i_atom];
                reference_fitting_cog_2 += reference_fitting_frames[min_frame_index_2][i_atom];
            }
            reference_fitting_cog_1 /= reference_fitting_frames[min_frame_index_1].size();
            reference_fitting_cog_2 /= reference_fitting_frames[min_frame_index_2].size();
            std::vector<cvm::atom_pos> tmp_reference_fitting_frame_1(reference_fitting_frames[min_frame_index_1].size());
            std::vector<cvm::atom_pos> tmp_reference_fitting_frame_2(reference_fitting_frames[min_frame_index_2].size());
            for (i_atom = 0; i_atom < reference_fitting_frames[min_frame_index_1].size(); ++i_atom) {
                tmp_reference_fitting_frame_1[i_atom] = reference_fitting_frames[min_frame_index_1][i_atom] - reference_fitting_cog_1;
                tmp_reference_fitting_frame_2[i_atom] = reference_fitting_frames[min_frame_index_2][i_atom] - reference_fitting_cog_2;
            }
            rot_v3.calc_optimal_rotation(tmp_reference_fitting_frame_1, tmp_reference_fitting_frame_2);
        } else {
            rot_v3.calc_optimal_rotation(tmp_reference_frame_1, tmp_reference_frame_2);
        }
        for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            v3[i_atom] = rot_v3.q.rotate(tmp_reference_frame_1[i_atom]) - tmp_reference_frame_2[i_atom];
        }
    } else {
        cvm::atom_pos reference_cog_1, reference_cog_3;
        for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            reference_cog_1 += reference_frames[min_frame_index_1][i_atom];
            reference_cog_3 += reference_frames[min_frame_index_3][i_atom];
        }
        reference_cog_1 /= reference_frames[min_frame_index_1].size();
        reference_cog_3 /= reference_frames[min_frame_index_3].size();
        std::vector<cvm::atom_pos> tmp_reference_frame_1(reference_frames[min_frame_index_1].size());
        std::vector<cvm::atom_pos> tmp_reference_frame_3(reference_frames[min_frame_index_3].size());
        for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            tmp_reference_frame_1[i_atom] = reference_frames[min_frame_index_1][i_atom] - reference_cog_1;
            tmp_reference_frame_3[i_atom] = reference_frames[min_frame_index_3][i_atom] - reference_cog_3;
        }
        if (has_user_defined_fitting) {
            cvm::atom_pos reference_fitting_cog_1, reference_fitting_cog_3;
            for (i_atom = 0; i_atom < reference_fitting_frames[min_frame_index_1].size(); ++i_atom) {
                reference_fitting_cog_1 += reference_fitting_frames[min_frame_index_1][i_atom];
                reference_fitting_cog_3 += reference_fitting_frames[min_frame_index_3][i_atom];
            }
            reference_fitting_cog_1 /= reference_fitting_frames[min_frame_index_1].size();
            reference_fitting_cog_3 /= reference_fitting_frames[min_frame_index_3].size();
            std::vector<cvm::atom_pos> tmp_reference_fitting_frame_1(reference_fitting_frames[min_frame_index_1].size());
            std::vector<cvm::atom_pos> tmp_reference_fitting_frame_3(reference_fitting_frames[min_frame_index_3].size());
            for (i_atom = 0; i_atom < reference_fitting_frames[min_frame_index_1].size(); ++i_atom) {
                tmp_reference_fitting_frame_1[i_atom] = reference_fitting_frames[min_frame_index_1][i_atom] - reference_fitting_cog_1;
                tmp_reference_fitting_frame_3[i_atom] = reference_fitting_frames[min_frame_index_3][i_atom] - reference_fitting_cog_3;
            }
            rot_v3.calc_optimal_rotation(tmp_reference_fitting_frame_1, tmp_reference_fitting_frame_3);
        } else {
            rot_v3.calc_optimal_rotation(tmp_reference_frame_1, tmp_reference_frame_3);
        }
        for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            // v3 = s_(m+1) - s_m
            v3[i_atom] = tmp_reference_frame_3[i_atom] - rot_v3.q.rotate(tmp_reference_frame_1[i_atom]);
        }
    }
}

void colvar::gspath::calc_value() {
    computeValue();
    x = s;
}

void colvar::gspath::calc_gradients() {
    computeDerivatives();
    cvm::rvector tmp_atom_grad_v1, tmp_atom_grad_v2;
    // dS(v1, v2(r), v3) / dr = ∂S/∂v1 * dv1/dr + ∂S/∂v2 * dv2/dr
    // dv1/dr = [fitting matrix 1][-1, ..., -1]
    // dv2/dr = [fitting matrix 2][1, ..., 1]
    // ∂S/∂v1 = ± (∂f/∂v1) / (2M)
    // ∂S/∂v2 = ± (∂f/∂v2) / (2M)
    // dS(v1, v2(r), v3) / dr = -1.0 * ± (∂f/∂v1) / (2M) + ± (∂f/∂v2) / (2M)
    for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        tmp_atom_grad_v1[0] = -1.0 * sign * 0.5 * dfdv1[i_atom][0] / M;
        tmp_atom_grad_v1[1] = -1.0 * sign * 0.5 * dfdv1[i_atom][1] / M;
        tmp_atom_grad_v1[2] = -1.0 * sign * 0.5 * dfdv1[i_atom][2] / M;
        tmp_atom_grad_v2[0] = sign * 0.5 * dfdv2[i_atom][0] / M;
        tmp_atom_grad_v2[1] = sign * 0.5 * dfdv2[i_atom][1] / M;
        tmp_atom_grad_v2[2] = sign * 0.5 * dfdv2[i_atom][2] / M;
        (*(comp_atoms[min_frame_index_1]))[i_atom].grad += tmp_atom_grad_v1;
        (*(comp_atoms[min_frame_index_2]))[i_atom].grad += tmp_atom_grad_v2;
    }
}

void colvar::gspath::apply_force(colvarvalue const &force) {
    // The force applied to this CV is scalar type
    cvm::real const &F = force.real_value;
    (*(comp_atoms[min_frame_index_1])).apply_colvar_force(F);
    (*(comp_atoms[min_frame_index_2])).apply_colvar_force(F);
}

colvar::gzpath::gzpath(std::string const &conf): CartesianBasedPath(conf) {
    function_type = "gzpath";
    get_keyval(conf, "useSecondClosestFrame", use_second_closest_frame, true);
    if (use_second_closest_frame == true) {
        cvm::log(std::string("Geometric path z(σ) will use the second closest frame to compute s_(m-1)\n"));
    } else {
        cvm::log(std::string("Geometric path z(σ) will use the neighbouring frame to compute s_(m-1)\n"));
    }
    get_keyval(conf, "useThirdClosestFrame", use_third_closest_frame, false);
    if (use_third_closest_frame == true) {
        cvm::log(std::string("Geometric path z(σ) will use the third closest frame to compute s_(m+1)\n"));
    } else {
        cvm::log(std::string("Geometric path z(σ) will use the neighbouring frame to compute s_(m+1)\n"));
    }
    bool b_use_z_square = false;
    get_keyval(conf, "useZsquare", b_use_z_square, false);
    if (b_use_z_square == true) {
        cvm::log(std::string("Geometric path z(σ) will use the square of distance from current frame to path compute z\n"));
    }
    GeometricPathCV::GeometricPathBase<cvm::atom_pos, cvm::real, GeometricPathCV::path_sz::Z>::initialize(atoms->size(), cvm::atom_pos(), total_reference_frames, use_second_closest_frame, use_third_closest_frame, b_use_z_square);
    // Logging
    cvm::log(std::string("Geometric pathCV(z) is initialized.\n"));
    cvm::log(std::string("Geometric pathCV(z) loaded ") + cvm::to_str(reference_frames.size()) + std::string(" frames.\n"));
}

void colvar::gzpath::updateReferenceDistances() {
    computeReferenceDistance(frame_distances);
}

void colvar::gzpath::prepareVectors() {
    cvm::atom_pos reference_cog_1, reference_cog_2;
    size_t i_atom;
    for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        reference_cog_1 += reference_frames[min_frame_index_1][i_atom];
        reference_cog_2 += reference_frames[min_frame_index_2][i_atom];
    }
    reference_cog_1 /= reference_frames[min_frame_index_1].size();
    reference_cog_2 /= reference_frames[min_frame_index_2].size();
    std::vector<cvm::atom_pos> tmp_reference_frame_1(reference_frames[min_frame_index_1].size());
    std::vector<cvm::atom_pos> tmp_reference_frame_2(reference_frames[min_frame_index_2].size());
    for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        tmp_reference_frame_1[i_atom] = reference_frames[min_frame_index_1][i_atom] - reference_cog_1;
        tmp_reference_frame_2[i_atom] = reference_frames[min_frame_index_2][i_atom] - reference_cog_2;
    }
    std::vector<cvm::atom_pos> tmp_reference_fitting_frame_1;
    std::vector<cvm::atom_pos> tmp_reference_fitting_frame_2;
    if (has_user_defined_fitting) {
        cvm::atom_pos reference_fitting_cog_1, reference_fitting_cog_2;
        for (i_atom = 0; i_atom < reference_fitting_frames[min_frame_index_1].size(); ++i_atom) {
            reference_fitting_cog_1 += reference_fitting_frames[min_frame_index_1][i_atom];
            reference_fitting_cog_2 += reference_fitting_frames[min_frame_index_2][i_atom];
        }
        reference_fitting_cog_1 /= reference_fitting_frames[min_frame_index_1].size();
        reference_fitting_cog_2 /= reference_fitting_frames[min_frame_index_2].size();
        tmp_reference_fitting_frame_1.resize(reference_fitting_frames[min_frame_index_1].size());
        tmp_reference_fitting_frame_2.resize(reference_fitting_frames[min_frame_index_2].size());
        for (i_atom = 0; i_atom < reference_fitting_frames[min_frame_index_1].size(); ++i_atom) {
            tmp_reference_fitting_frame_1[i_atom] = reference_fitting_frames[min_frame_index_1][i_atom] - reference_fitting_cog_1;
            tmp_reference_fitting_frame_2[i_atom] = reference_fitting_frames[min_frame_index_2][i_atom] - reference_fitting_cog_2;
        }
        rot_v4.calc_optimal_rotation(tmp_reference_fitting_frame_1, tmp_reference_fitting_frame_2);
    } else {
        rot_v4.calc_optimal_rotation(tmp_reference_frame_1, tmp_reference_frame_2);
    }
    for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        v1[i_atom] = reference_frames[min_frame_index_1][i_atom] - (*(comp_atoms[min_frame_index_1]))[i_atom].pos;
        v2[i_atom] = (*(comp_atoms[min_frame_index_2]))[i_atom].pos - reference_frames[min_frame_index_2][i_atom];
        // v4 only computes in gzpath
        // v4 = s_m - s_(m-1)
        v4[i_atom] = rot_v4.q.rotate(tmp_reference_frame_1[i_atom]) - tmp_reference_frame_2[i_atom];
    }
    if (min_frame_index_3 < 0 || min_frame_index_3 > M) {
        v3 = v4;
    } else {
        cvm::atom_pos reference_cog_3;
        for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            reference_cog_3 += reference_frames[min_frame_index_3][i_atom];
        }
        reference_cog_3 /= reference_frames[min_frame_index_3].size();
        std::vector<cvm::atom_pos> tmp_reference_frame_3(reference_frames[min_frame_index_3].size());
        for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            tmp_reference_frame_3[i_atom] = reference_frames[min_frame_index_3][i_atom] - reference_cog_3;
        }
        if (has_user_defined_fitting) {
            cvm::atom_pos reference_fitting_cog_3;
            for (i_atom = 0; i_atom < reference_fitting_frames[min_frame_index_3].size(); ++i_atom) {
                reference_fitting_cog_3 += reference_fitting_frames[min_frame_index_3][i_atom];
            }
            reference_fitting_cog_3 /= reference_fitting_frames[min_frame_index_3].size();
            std::vector<cvm::atom_pos> tmp_reference_fitting_frame_3(reference_fitting_frames[min_frame_index_3].size());
            for (i_atom = 0; i_atom < reference_fitting_frames[min_frame_index_3].size(); ++i_atom) {
                tmp_reference_fitting_frame_3[i_atom] =  reference_fitting_frames[min_frame_index_3][i_atom] - reference_fitting_cog_3;
            }
            rot_v3.calc_optimal_rotation(tmp_reference_fitting_frame_1, tmp_reference_fitting_frame_3);
        } else {
            rot_v3.calc_optimal_rotation(tmp_reference_frame_1, tmp_reference_frame_3);
        }
        for (i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            // v3 = s_(m+1) - s_m
            v3[i_atom] = tmp_reference_frame_3[i_atom] - rot_v3.q.rotate(tmp_reference_frame_1[i_atom]);
        }
    }
}

void colvar::gzpath::calc_value() {
    computeValue();
    x = z;
}

void colvar::gzpath::calc_gradients() {
    computeDerivatives();
    cvm::rvector tmp_atom_grad_v1, tmp_atom_grad_v2;
    for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        tmp_atom_grad_v1 = -1.0 * dzdv1[i_atom];
        tmp_atom_grad_v2 = dzdv2[i_atom];
        (*(comp_atoms[min_frame_index_1]))[i_atom].grad += tmp_atom_grad_v1;
        (*(comp_atoms[min_frame_index_2]))[i_atom].grad += tmp_atom_grad_v2;
    }
}

void colvar::gzpath::apply_force(colvarvalue const &force) {
    // The force applied to this CV is scalar type
    cvm::real const &F = force.real_value;
    (*(comp_atoms[min_frame_index_1])).apply_colvar_force(F);
    (*(comp_atoms[min_frame_index_2])).apply_colvar_force(F);
}

colvar::linearCombination::linearCombination(std::string const &conf): cvc(conf) {
    GeometricPathCV::init_string_cv_map(string_cv_map);
    // Lookup all available sub-cvcs
    for (auto it_cv_map = string_cv_map.begin(); it_cv_map != string_cv_map.end(); ++it_cv_map) {
        if (key_lookup(conf, it_cv_map->first.c_str())) {
            std::vector<std::string> sub_cvc_confs;
            get_key_string_multi_value(conf, it_cv_map->first.c_str(), sub_cvc_confs);
            for (auto it_sub_cvc_conf = sub_cvc_confs.begin(); it_sub_cvc_conf != sub_cvc_confs.end(); ++it_sub_cvc_conf) {
                cv.push_back((it_cv_map->second)(*(it_sub_cvc_conf)));
            }
        }
    }
    // Sort all sub CVs by their names
    std::sort(cv.begin(), cv.end(), compareColvarComponent);
    for (auto it_sub_cv = cv.begin(); it_sub_cv != cv.end(); ++it_sub_cv) {
        for (auto it_atom_group = (*it_sub_cv)->atom_groups.begin(); it_atom_group != (*it_sub_cv)->atom_groups.end(); ++it_atom_group) {
            register_atom_group(*it_atom_group);
        }
    }
    x.type(cv[0]->value());
    x.reset();
    use_explicit_gradients = true;
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
    }
    if (!use_explicit_gradients) {
        disable(f_cvc_explicit_gradient);
    }
}

cvm::real colvar::linearCombination::getPolynomialFactorOfCVGradient(size_t i_cv) const {
    cvm::real factor_polynomial = 1.0;
    if (cv[i_cv]->value().type() == colvarvalue::type_scalar) {
        factor_polynomial = cv[i_cv]->sup_coeff * cv[i_cv]->sup_np * cvm::pow(cv[i_cv]->value().real_value, cv[i_cv]->sup_np - 1);
    } else {
        factor_polynomial = cv[i_cv]->sup_coeff;
    }
    return factor_polynomial;
}

colvar::linearCombination::~linearCombination() {
    for (auto it = cv.begin(); it != cv.end(); ++it) {
        delete (*it);
    }
}

void colvar::linearCombination::calc_value() {
    x.reset();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_value();
        colvarvalue current_cv_value(cv[i_cv]->value());
        // polynomial combination allowed
        if (current_cv_value.type() == colvarvalue::type_scalar) {
            x += cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
        } else {
            x += cv[i_cv]->sup_coeff * current_cv_value;
        }
    }
}

void colvar::linearCombination::calc_gradients() {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_gradients();
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)) {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = factor_polynomial * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                    }
                }
            }
        }
    }
}

void colvar::linearCombination::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // If this CV us explicit gradients, then atomic gradients is already calculated
        // We can apply the force to atom groups directly
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)
        ) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            // Compute factors for polynomial combinations
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            colvarvalue cv_force = force.real_value * factor_polynomial;
            cv[i_cv]->apply_force(cv_force);
        }
    }
}

colvar::CVBasedPath::CVBasedPath(std::string const &conf): cvc(conf) {
    GeometricPathCV::init_string_cv_map(string_cv_map);
    // Lookup all available sub-cvcs
    for (auto it_cv_map = string_cv_map.begin(); it_cv_map != string_cv_map.end(); ++it_cv_map) {
        if (key_lookup(conf, it_cv_map->first.c_str())) {
            std::vector<std::string> sub_cvc_confs;
            get_key_string_multi_value(conf, it_cv_map->first.c_str(), sub_cvc_confs);
            for (auto it_sub_cvc_conf = sub_cvc_confs.begin(); it_sub_cvc_conf != sub_cvc_confs.end(); ++it_sub_cvc_conf) {
                cv.push_back((it_cv_map->second)(*(it_sub_cvc_conf)));
            }
        }
    }
    // Sort all sub CVs by their names
    std::sort(cv.begin(), cv.end(), compareColvarComponent);
    // Register atom groups and determine the colvar type for reference
    std::vector<colvarvalue> tmp_cv;
    for (auto it_sub_cv = cv.begin(); it_sub_cv != cv.end(); ++it_sub_cv) {
        for (auto it_atom_group = (*it_sub_cv)->atom_groups.begin(); it_atom_group != (*it_sub_cv)->atom_groups.end(); ++it_atom_group) {
            register_atom_group(*it_atom_group);
        }
        colvarvalue tmp_i_cv((*it_sub_cv)->value());
        tmp_i_cv.reset();
        tmp_cv.push_back(tmp_i_cv);
    }
    // Read path file
    // Lookup all reference CV values
    std::string path_filename;
    get_keyval(conf, "pathFile", path_filename);
    cvm::log(std::string("Reading path file: ") + path_filename + std::string("\n"));
    std::ifstream ifs_path(path_filename);
    if (!ifs_path.is_open()) {
        cvm::error("Error: failed to open path file.\n");
    }
    std::string line;
    const std::string token(" ");
    total_reference_frames = 0;
    while (std::getline(ifs_path, line)) {
        std::vector<std::string> fields;
        split_string(line, token, fields);
        size_t num_value_required = 0;
        for (size_t i_cv = 0; i_cv < tmp_cv.size(); ++i_cv) {
            const size_t value_size = tmp_cv[i_cv].size();
            num_value_required += value_size;
            cvm::log(std::string("Reading CV ") + cv[i_cv]->name + std::string(" with ") + cvm::to_str(value_size) + std::string(" value(s)\n"));
            if (num_value_required <= fields.size()) {
                size_t start_index = num_value_required - value_size;
                for (size_t i = start_index; i < num_value_required; ++i) {
                    tmp_cv[i_cv][i] = std::atof(fields[i].c_str());
                    cvm::log(fields[i] + std::string(" "));
                }
                cvm::log(std::string("\n"));
            } else {
                cvm::error("Error: incorrect format of path file.\n");
            }
        }
        if (!fields.empty()) {
            ref_cv.push_back(tmp_cv);
            ++total_reference_frames;
        }
    }
    x.type(colvarvalue::type_scalar);
    use_explicit_gradients = true;
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
    }
    if (!use_explicit_gradients) {
        disable(f_cvc_explicit_gradient);
    }
}

void colvar::CVBasedPath::computeReferenceDistance(std::vector<cvm::real>& result) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_value();
    }
    for (size_t i_frame = 0; i_frame < ref_cv.size(); ++i_frame) {
        cvm::real rmsd_i = 0.0;
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            colvarvalue ref_cv_value(ref_cv[i_frame][i_cv]);
            colvarvalue current_cv_value(cv[i_cv]->value());
            // polynomial combination allowed
            if (current_cv_value.type() == colvarvalue::type_scalar) {
                // wrapping is already in dist2
                rmsd_i += cv[i_cv]->dist2(cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np)), ref_cv_value.real_value);
            } else {
                rmsd_i += cv[i_cv]->dist2(cv[i_cv]->sup_coeff * current_cv_value, ref_cv_value);
            }
        }
        rmsd_i /= cvm::real(cv.size());
        rmsd_i = cvm::sqrt(rmsd_i);
        result[i_frame] = rmsd_i;
    }
}

cvm::real colvar::CVBasedPath::getPolynomialFactorOfCVGradient(size_t i_cv) const {
    cvm::real factor_polynomial = 1.0;
    if (cv[i_cv]->value().type() == colvarvalue::type_scalar) {
        factor_polynomial = cv[i_cv]->sup_coeff * cv[i_cv]->sup_np * cvm::pow(cv[i_cv]->value().real_value, cv[i_cv]->sup_np - 1);
    } else {
        factor_polynomial = cv[i_cv]->sup_coeff;
    }
    return factor_polynomial;
}

colvar::CVBasedPath::~CVBasedPath() {
    for (auto it = cv.begin(); it != cv.end(); ++it) {
        delete (*it);
    }
}

colvar::gspathCV::gspathCV(std::string const &conf): CVBasedPath(conf) {
    function_type = "gspathCV";
    cvm::log(std::string("Total number of frames: ") + cvm::to_str(total_reference_frames) + std::string("\n"));
    // Initialize variables for future calculation
    get_keyval(conf, "useSecondClosestFrame", use_second_closest_frame, true);
    if (use_second_closest_frame == true) {
        cvm::log(std::string("Geometric path s(σ) will use the second closest frame to compute s_(m-1)\n"));
    } else {
        cvm::log(std::string("Geometric path s(σ) will use the neighbouring frame to compute s_(m-1)\n"));
    }
    get_keyval(conf, "useThirdClosestFrame", use_third_closest_frame, false);
    if (use_third_closest_frame == true) {
        cvm::log(std::string("Geometric path s(σ) will use the third closest frame to compute s_(m+1)\n"));
    } else {
        cvm::log(std::string("Geometric path s(σ) will use the neighbouring frame to compute s_(m+1)\n"));
    }
    GeometricPathCV::GeometricPathBase<colvarvalue, cvm::real, GeometricPathCV::path_sz::S>::initialize(cv.size(), ref_cv[0], total_reference_frames, use_second_closest_frame, use_third_closest_frame);
    x.type(colvarvalue::type_scalar);
    use_explicit_gradients = true;
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
    }
    if (!use_explicit_gradients) {
        cvm::log("Geometric path s(σ) will use implicit gradients.\n");
        disable(f_cvc_explicit_gradient);
    }
}

colvar::gspathCV::~gspathCV() {}

void colvar::gspathCV::updateReferenceDistances() {
    computeReferenceDistance(frame_distances);
}

void colvar::gspathCV::prepareVectors() {
    // Compute v1, v2 and v3
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // values of sub-cvc are computed in update_distances
        // cv[i_cv]->calc_value();
        colvarvalue f1_ref_cv_i_value(ref_cv[min_frame_index_1][i_cv]);
        colvarvalue f2_ref_cv_i_value(ref_cv[min_frame_index_2][i_cv]);
        colvarvalue current_cv_value(cv[i_cv]->value());
        // polynomial combination allowed
        if (current_cv_value.type() == colvarvalue::type_scalar) {
            v1[i_cv] = f1_ref_cv_i_value.real_value - cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
            v2[i_cv] = cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np)) - f2_ref_cv_i_value.real_value;
        } else {
            v1[i_cv] = f1_ref_cv_i_value - cv[i_cv]->sup_coeff * current_cv_value;
            v2[i_cv] = cv[i_cv]->sup_coeff * current_cv_value - f2_ref_cv_i_value;
        }
        cv[i_cv]->wrap(v1[i_cv]);
        cv[i_cv]->wrap(v2[i_cv]);
    }
    if (min_frame_index_3 < 0 || min_frame_index_3 > M) {
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            v3[i_cv] = ref_cv[min_frame_index_1][i_cv] - ref_cv[min_frame_index_2][i_cv];
            cv[i_cv]->wrap(v3[i_cv]);
        }
    } else {
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            v3[i_cv] = ref_cv[min_frame_index_3][i_cv] - ref_cv[min_frame_index_1][i_cv];
            cv[i_cv]->wrap(v3[i_cv]);
        }
    }
}

void colvar::gspathCV::calc_value() {
    computeValue();
    x = s;
}

void colvar::gspathCV::calc_gradients() {
    computeDerivatives();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // No matter whether the i-th cv uses implicit gradient, compute it first.
        cv[i_cv]->calc_gradients();
        // If the gradient is not implicit, then add the gradients to its atom groups
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)) {
            // Temporary variables storing gradients
            colvarvalue tmp_cv_grad_v1(cv[i_cv]->value());
            colvarvalue tmp_cv_grad_v2(cv[i_cv]->value());
            // Compute factors for polynomial combinations
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            // Loop over all elements of the corresponding colvar value
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                // ds/dz, z = vector of CVs
                tmp_cv_grad_v1[j_elem] = -1.0 * sign * 0.5 * dfdv1[i_cv][j_elem] / M;
                tmp_cv_grad_v2[j_elem] = sign * 0.5 * dfdv2[i_cv][j_elem] / M;
                // Apply the gradients to the atom groups in i-th cv
                // Loop over all atom groups
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    // Loop over all atoms in the k-th atom group
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        // Chain rule
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = factor_polynomial * ((*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad * tmp_cv_grad_v1[j_elem] + (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad * tmp_cv_grad_v2[j_elem]);
                    }
                }
            }
        }
    }
}

void colvar::gspathCV::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // If this CV us explicit gradients, then atomic gradients is already calculated
        // We can apply the force to atom groups directly
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)
        ) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            // Temporary variables storing gradients
            colvarvalue tmp_cv_grad_v1(cv[i_cv]->value());
            colvarvalue tmp_cv_grad_v2(cv[i_cv]->value());
            // Compute factors for polynomial combinations
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                // ds/dz, z = vector of CVs
                tmp_cv_grad_v1[j_elem] = -1.0 * sign * 0.5 * dfdv1[i_cv][j_elem] / M;
                tmp_cv_grad_v2[j_elem] = sign * 0.5 * dfdv2[i_cv][j_elem] / M;
            }
            colvarvalue cv_force = force.real_value * factor_polynomial * (tmp_cv_grad_v1 + tmp_cv_grad_v2);
            cv[i_cv]->apply_force(cv_force);
        }
    }
}

colvar::gzpathCV::gzpathCV(std::string const &conf): CVBasedPath(conf) {
    function_type = "gzpathCV";
    cvm::log(std::string("Total number of frames: ") + cvm::to_str(total_reference_frames) + std::string("\n"));
    // Initialize variables for future calculation
    M = cvm::real(total_reference_frames - 1);
    m = 1.0;
    get_keyval(conf, "useSecondClosestFrame", use_second_closest_frame, true);
    if (use_second_closest_frame == true) {
        cvm::log(std::string("Geometric path z(σ) will use the second closest frame to compute s_(m-1)\n"));
    } else {
        cvm::log(std::string("Geometric path z(σ) will use the neighbouring frame to compute s_(m-1)\n"));
    }
    get_keyval(conf, "useThirdClosestFrame", use_third_closest_frame, false);
    if (use_third_closest_frame == true) {
        cvm::log(std::string("Geometric path z(σ) will use the third closest frame to compute s_(m+1)\n"));
    } else {
        cvm::log(std::string("Geometric path z(σ) will use the neighbouring frame to compute s_(m+1)\n"));
    }
    bool b_use_z_square = false;
    get_keyval(conf, "useZsquare", b_use_z_square, false);
    if (b_use_z_square == true) {
        cvm::log(std::string("Geometric path z(σ) will use the square of distance from current frame to path compute z\n"));
    }
    GeometricPathCV::GeometricPathBase<colvarvalue, cvm::real, GeometricPathCV::path_sz::Z>::initialize(cv.size(), ref_cv[0], total_reference_frames, use_second_closest_frame, use_third_closest_frame, b_use_z_square);
    x.type(colvarvalue::type_scalar);
    use_explicit_gradients = true;
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
    }
    if (!use_explicit_gradients) {
        cvm::log("Geometric path z(σ) will use implicit gradients.\n");
        disable(f_cvc_explicit_gradient);
    }
}

colvar::gzpathCV::~gzpathCV() {
}

void colvar::gzpathCV::updateReferenceDistances() {
    computeReferenceDistance(frame_distances);
}

void colvar::gzpathCV::prepareVectors() {
    // Compute v1, v2 and v3
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // values of sub-cvc are computed in update_distances
        // cv[i_cv]->calc_value();
        colvarvalue f1_ref_cv_i_value(ref_cv[min_frame_index_1][i_cv]);
        colvarvalue f2_ref_cv_i_value(ref_cv[min_frame_index_2][i_cv]);
        colvarvalue current_cv_value(cv[i_cv]->value());
        // polynomial combination allowed
        if (current_cv_value.type() == colvarvalue::type_scalar) {
            v1[i_cv] = f1_ref_cv_i_value.real_value - cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
            v2[i_cv] = cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np)) - f2_ref_cv_i_value.real_value;
        } else {
            v1[i_cv] = f1_ref_cv_i_value - cv[i_cv]->sup_coeff * current_cv_value;
            v2[i_cv] = cv[i_cv]->sup_coeff * current_cv_value - f2_ref_cv_i_value;
        }
        v4[i_cv] = f1_ref_cv_i_value - f2_ref_cv_i_value;
        cv[i_cv]->wrap(v1[i_cv]);
        cv[i_cv]->wrap(v2[i_cv]);
        cv[i_cv]->wrap(v4[i_cv]);
    }
    if (min_frame_index_3 < 0 || min_frame_index_3 > M) {
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            v3[i_cv] = ref_cv[min_frame_index_1][i_cv] - ref_cv[min_frame_index_2][i_cv];
            cv[i_cv]->wrap(v3[i_cv]);
        }
    } else {
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            v3[i_cv] = ref_cv[min_frame_index_3][i_cv] - ref_cv[min_frame_index_1][i_cv];
            cv[i_cv]->wrap(v3[i_cv]);
        }
    }
}

void colvar::gzpathCV::calc_value() {
    computeValue();
    x = z;
}

void colvar::gzpathCV::calc_gradients() {
    computeDerivatives();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // No matter whether the i-th cv uses implicit gradient, compute it first.
        cv[i_cv]->calc_gradients();
        // If the gradient is not implicit, then add the gradients to its atom groups
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)) {
            // Temporary variables storing gradients
            colvarvalue tmp_cv_grad_v1 = -1.0 * dzdv1[i_cv];
            colvarvalue tmp_cv_grad_v2 =  1.0 * dzdv2[i_cv];
            // Compute factors for polynomial combinations
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                // Apply the gradients to the atom groups in i-th cv
                // Loop over all atom groups
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    // Loop over all atoms in the k-th atom group
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        // Chain rule
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = factor_polynomial * ((*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad * tmp_cv_grad_v1[j_elem] + (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad * tmp_cv_grad_v2[j_elem]);
                    }
                }
            }
        }
    }
}

void colvar::gzpathCV::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // If this CV us explicit gradients, then atomic gradients is already calculated
        // We can apply the force to atom groups directly
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        }
        else {
            colvarvalue tmp_cv_grad_v1 = -1.0 * dzdv1[i_cv];
            colvarvalue tmp_cv_grad_v2 =  1.0 * dzdv2[i_cv];
            // Temporary variables storing gradients
            // Compute factors for polynomial combinations
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            colvarvalue cv_force = force.real_value * factor_polynomial * (tmp_cv_grad_v1 + tmp_cv_grad_v2);
            cv[i_cv]->apply_force(cv_force);
        }
    }
}

void GeometricPathCV::init_string_cv_map(std::map<std::string, std::function<colvar::cvc* (const std::string& subcv_conf)>>& string_cv_map) {
    string_cv_map["distance"]              = [](const std::string& conf){return new colvar::distance(conf);};
    string_cv_map["dihedral"]              = [](const std::string& conf){return new colvar::dihedral(conf);};
    string_cv_map["angle"]                 = [](const std::string& conf){return new colvar::angle(conf);};
    string_cv_map["rmsd"]                  = [](const std::string& conf){return new colvar::rmsd(conf);};
    string_cv_map["gyration"]              = [](const std::string& conf){return new colvar::gyration(conf);};
    string_cv_map["inertia"]               = [](const std::string& conf){return new colvar::inertia(conf);};
    string_cv_map["inertiaZ"]              = [](const std::string& conf){return new colvar::inertia_z(conf);};
    string_cv_map["tilt"]                  = [](const std::string& conf){return new colvar::tilt(conf);};
    string_cv_map["distanceZ"]             = [](const std::string& conf){return new colvar::distance_z(conf);};
    string_cv_map["distanceXY"]            = [](const std::string& conf){return new colvar::distance_xy(conf);};
    string_cv_map["polarTheta"]            = [](const std::string& conf){return new colvar::polar_theta(conf);};
    string_cv_map["polarPhi"]              = [](const std::string& conf){return new colvar::polar_phi(conf);};
    string_cv_map["distanceVec"]           = [](const std::string& conf){return new colvar::distance_vec(conf);};
    string_cv_map["orientationAngle"]      = [](const std::string& conf){return new colvar::orientation_angle(conf);};
    string_cv_map["distancePairs"]         = [](const std::string& conf){return new colvar::distance_pairs(conf);};
    string_cv_map["dipoleMagnitude"]       = [](const std::string& conf){return new colvar::dipole_magnitude(conf);};
    string_cv_map["coordNum"]              = [](const std::string& conf){return new colvar::coordnum(conf);};
    string_cv_map["selfCoordNum"]          = [](const std::string& conf){return new colvar::selfcoordnum(conf);};
    string_cv_map["dipoleAngle"]           = [](const std::string& conf){return new colvar::dipole_angle(conf);};
    string_cv_map["orientation"]           = [](const std::string& conf){return new colvar::orientation(conf);};
    string_cv_map["orientationProj"]       = [](const std::string& conf){return new colvar::orientation_proj(conf);};
    string_cv_map["eigenvector"]           = [](const std::string& conf){return new colvar::eigenvector(conf);};
    string_cv_map["cartesian"]             = [](const std::string& conf){return new colvar::cartesian(conf);};
    string_cv_map["alpha"]                 = [](const std::string& conf){return new colvar::alpha_angles(conf);};
    string_cv_map["dihedralPC"]            = [](const std::string& conf){return new colvar::dihedPC(conf);};
    string_cv_map["linearCombination"]     = [](const std::string& conf){return new colvar::linearCombination(conf);};
}

#endif
