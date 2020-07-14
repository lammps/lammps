#ifndef GEOMETRICPATHCV_H
#define GEOMETRICPATHCV_H
// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include "colvarmodule.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <map>

namespace GeometricPathCV {

enum path_sz {S, Z};

template <typename element_type, typename scalar_type, path_sz path_type>
class GeometricPathBase {
private:
    struct doCompareFrameDistance {
        doCompareFrameDistance(const GeometricPathBase& obj): m_obj(obj) {}
        const GeometricPathBase& m_obj;
        bool operator()(const size_t& i1, const size_t& i2) {
            return m_obj.frame_distances[i1] < m_obj.frame_distances[i2];
        }
    };
protected:
    scalar_type v1v1;
    scalar_type v2v2;
    scalar_type v3v3;
    scalar_type v4v4;
    scalar_type v1v3;
    scalar_type v1v4;
    scalar_type f;
    scalar_type dx;
    scalar_type s;
    scalar_type z;
    scalar_type zz;
    std::vector<element_type> v1;
    std::vector<element_type> v2;
    std::vector<element_type> v3;
    std::vector<element_type> v4;
    std::vector<element_type> dfdv1;
    std::vector<element_type> dfdv2;
    std::vector<element_type> dzdv1;
    std::vector<element_type> dzdv2;
    std::vector<scalar_type> frame_distances;
    std::vector<size_t> frame_index;
    bool use_second_closest_frame;
    bool use_third_closest_frame;
    bool use_z_square;
    long min_frame_index_1;
    long min_frame_index_2;
    long min_frame_index_3;
    long sign;
    double M;
    double m;
public:
    GeometricPathBase(size_t vector_size, const element_type& element = element_type(), size_t total_frames = 1, bool p_use_second_closest_frame = true, bool p_use_third_closest_frame = false, bool p_use_z_square = false);
    GeometricPathBase(size_t vector_size, const std::vector<element_type>& elements, size_t total_frames = 1, bool p_use_second_closest_frame = true, bool p_use_third_closest_frame = false, bool p_use_z_square = false);
    GeometricPathBase() {}
    virtual ~GeometricPathBase() {}
    virtual void initialize(size_t vector_size, const element_type& element = element_type(), size_t total_frames = 1, bool p_use_second_closest_frame = true, bool p_use_third_closest_frame = false, bool p_use_z_square = false);
    virtual void initialize(size_t vector_size, const std::vector<element_type>& elements, size_t total_frames = 1, bool p_use_second_closest_frame = true, bool p_use_third_closest_frame = false, bool p_use_z_square = false);
    virtual void prepareVectors() = 0;
    virtual void updateDistanceToReferenceFrames() = 0;
    virtual void compute();
    virtual void determineClosestFrames();
    virtual void computeValue();
    virtual void computeDerivatives();
};

template <typename element_type, typename scalar_type, path_sz path_type>
GeometricPathBase<element_type, scalar_type, path_type>::GeometricPathBase(size_t vector_size, const element_type& element, size_t total_frames, bool p_use_second_closest_frame, bool p_use_third_closest_frame, bool p_use_z_square) {
    initialize(vector_size, element, total_frames, p_use_second_closest_frame, p_use_third_closest_frame, p_use_z_square);
}

template <typename element_type, typename scalar_type, path_sz path_type>
GeometricPathBase<element_type, scalar_type, path_type>::GeometricPathBase(size_t vector_size, const std::vector<element_type>& elements, size_t total_frames, bool p_use_second_closest_frame, bool p_use_third_closest_frame, bool p_use_z_square) {
    initialize(vector_size, elements, total_frames, p_use_second_closest_frame, p_use_third_closest_frame, p_use_z_square);
}

template <typename element_type, typename scalar_type, path_sz path_type>
void GeometricPathBase<element_type, scalar_type, path_type>::initialize(size_t vector_size, const element_type& element, size_t total_frames, bool p_use_second_closest_frame, bool p_use_third_closest_frame, bool p_use_z_square) {
    v1v1 = scalar_type();
    v2v2 = scalar_type();
    v3v3 = scalar_type();
    v4v4 = scalar_type();
    v1v3 = scalar_type();
    v1v4 = scalar_type();
    f = scalar_type();
    dx = scalar_type();
    z = scalar_type();
    zz = scalar_type();
    sign = 0;
    v1.resize(vector_size, element);
    v2.resize(vector_size, element);
    v3.resize(vector_size, element);
    v4.resize(vector_size, element);
    dfdv1.resize(vector_size, element);
    dfdv2.resize(vector_size, element);
    dzdv1.resize(vector_size, element);
    dzdv2.resize(vector_size, element);
    frame_distances.resize(total_frames);
    frame_index.resize(total_frames);
    for (size_t i_frame = 0; i_frame < frame_index.size(); ++i_frame) {
        frame_index[i_frame] = i_frame;
    }
    use_second_closest_frame = p_use_second_closest_frame;
    use_third_closest_frame = p_use_third_closest_frame;
    use_z_square = p_use_z_square;
    M = static_cast<scalar_type>(total_frames - 1);
    m = static_cast<scalar_type>(1.0);
}

template <typename element_type, typename scalar_type, path_sz path_type>
void GeometricPathBase<element_type, scalar_type, path_type>::initialize(size_t /* vector_size */, const std::vector<element_type>& elements, size_t total_frames, bool p_use_second_closest_frame, bool p_use_third_closest_frame, bool p_use_z_square) {
    v1v1 = scalar_type();
    v2v2 = scalar_type();
    v3v3 = scalar_type();
    v4v4 = scalar_type();
    v1v3 = scalar_type();
    v1v4 = scalar_type();
    f = scalar_type();
    dx = scalar_type();
    z = scalar_type();
    zz = scalar_type();
    sign = 0;
    v1 = elements;
    v2 = elements;
    v3 = elements;
    v4 = elements;
    dfdv1 = elements;
    dfdv2 = elements;
    dzdv1 = elements;
    dzdv2 = elements;
    frame_distances.resize(total_frames);
    frame_index.resize(total_frames);
    for (size_t i_frame = 0; i_frame < frame_index.size(); ++i_frame) {
        frame_index[i_frame] = i_frame;
    }
    use_second_closest_frame = p_use_second_closest_frame;
    use_third_closest_frame = p_use_third_closest_frame;
    use_z_square = p_use_z_square;
    M = static_cast<scalar_type>(total_frames - 1);
    m = static_cast<scalar_type>(1.0);
}

template <typename element_type, typename scalar_type, path_sz path_type>
void GeometricPathBase<element_type, scalar_type, path_type>::compute() {
    computeValue();
    computeDerivatives();
}

template <typename element_type, typename scalar_type, path_sz path_type>
void GeometricPathBase<element_type, scalar_type, path_type>::determineClosestFrames() {
    // Find the closest and the second closest frames
    std::sort(frame_index.begin(), frame_index.end(), doCompareFrameDistance(*this));
    // Determine the sign
    sign = static_cast<long>(frame_index[0]) - static_cast<long>(frame_index[1]);
    if (sign > 1) {
        // sigma(z) is on the left side of the closest frame
        sign = 1;
    } else if (sign < -1) {
        // sigma(z) is on the right side of the closest frame
        sign = -1;
    }
    if (cvm::fabs(static_cast<long>(frame_index[0]) - static_cast<long>(frame_index[1])) > 1) {
        std::cout << "Warning: Geometrical pathCV relies on the assumption that the second closest frame is the neighbouring frame\n";
        std::cout << "         Please check your configuration or increase restraint on z(σ)\n";
        for (size_t i_frame = 0; i_frame < frame_index.size(); ++i_frame) {
            std::cout << "Frame index: " << frame_index[i_frame] << " ; optimal RMSD = " << frame_distances[frame_index[i_frame]] << "\n";
        }
    }
    min_frame_index_1 = frame_index[0];                                                         // s_m
    min_frame_index_2 = use_second_closest_frame ? frame_index[1] : min_frame_index_1 - sign;   // s_(m-1)
    min_frame_index_3 = use_third_closest_frame ? frame_index[2] : min_frame_index_1 + sign;    // s_(m+1)
    m = static_cast<double>(frame_index[0]);
}

template <typename element_type, typename scalar_type, path_sz path_type>
void GeometricPathBase<element_type, scalar_type, path_type>::computeValue() {
    updateDistanceToReferenceFrames();
    determineClosestFrames();
    prepareVectors();
    v1v1 = scalar_type();
    v2v2 = scalar_type();
    v3v3 = scalar_type();
    v1v3 = scalar_type();
    if (path_type == Z) {
        v1v4 = scalar_type();
        v4v4 = scalar_type();
    }
    for (size_t i_elem = 0; i_elem < v1.size(); ++i_elem) {
        v1v1 += v1[i_elem] * v1[i_elem];
        v2v2 += v2[i_elem] * v2[i_elem];
        v3v3 += v3[i_elem] * v3[i_elem];
        v1v3 += v1[i_elem] * v3[i_elem];
        if (path_type == Z) {
            v1v4 += v1[i_elem] * v4[i_elem];
            v4v4 += v4[i_elem] * v4[i_elem];
        }
    }
    f = (cvm::sqrt(v1v3 * v1v3 - v3v3 * (v1v1 - v2v2)) - v1v3) / v3v3;
    if (path_type == Z) {
        dx = 0.5 * (f - 1);
        zz = v1v1 + 2 * dx * v1v4 + dx * dx * v4v4;
        if (use_z_square) {
            z = zz;
        } else {
            z = cvm::sqrt(cvm::fabs(zz));
        }
    }
    if (path_type == S) {
        s = m/M + static_cast<double>(sign) * ((f - 1) / (2 * M));
    }
}

template <typename element_type, typename scalar_type, path_sz path_type>
void GeometricPathBase<element_type, scalar_type, path_type>::computeDerivatives() {
    const scalar_type factor1 = 1.0 / (2.0 * v3v3 * cvm::sqrt(v1v3 * v1v3 - v3v3 * (v1v1 - v2v2)));
    const scalar_type factor2 = 1.0 / v3v3;
    for (size_t i_elem = 0; i_elem < v1.size(); ++i_elem) {
        // Compute the derivative of f with vector v1
        dfdv1[i_elem] = factor1 * (2.0 * v1v3 * v3[i_elem] - 2.0 * v3v3 * v1[i_elem]) - factor2 * v3[i_elem];
        // Compute the derivative of f with respect to vector v2
        dfdv2[i_elem] = factor1 * (2.0 * v3v3 * v2[i_elem]);
        // dZ(v1(r), v2(r), v3) / dr = ∂Z/∂v1 * dv1/dr + ∂Z/∂v2 * dv2/dr
        // dv1/dr = [fitting matrix 1][-1, ..., -1]
        // dv2/dr = [fitting matrix 2][1, ..., 1]
        // ∂Z/∂v1 = 1/(2*z) * (2v1 + (f-1)v4 + (v1⋅v4)∂f/∂v1 + v4^2 * 1/4 * 2(f-1) * ∂f/∂v1)
        // ∂Z/∂v2 = 1/(2*z) * ((v1⋅v4)∂f/∂v2 + v4^2 * 1/4 * 2(f-1) * ∂f/∂v2)
        if (path_type == Z) {
            if (use_z_square) {
                dzdv1[i_elem] = 2.0 * v1[i_elem] + (f-1) * v4[i_elem] + v1v4 * dfdv1[i_elem] + v4v4 * 0.25 * 2.0 * (f-1) * dfdv1[i_elem];
                dzdv2[i_elem] = v1v4 * dfdv2[i_elem] + v4v4 * 0.25 * 2.0 * (f-1) * dfdv2[i_elem];
            } else {
                if (z > static_cast<scalar_type>(0)) {
                    dzdv1[i_elem] = (1.0 / (2.0 * z)) * (2.0 * v1[i_elem] + (f-1) * v4[i_elem] + v1v4 * dfdv1[i_elem] + v4v4 * 0.25 * 2.0 * (f-1) * dfdv1[i_elem]);
                    dzdv2[i_elem] = (1.0 / (2.0 * z)) * (v1v4 * dfdv2[i_elem] + v4v4 * 0.25 * 2.0 * (f-1) * dfdv2[i_elem]);
                } else {
                    // workaround at z = 0
                    dzdv1[i_elem] = 0;
                    dzdv2[i_elem] = 0;
                }
            }
        }
    }
}

}

#endif // GEOMETRICPATHCV_H
