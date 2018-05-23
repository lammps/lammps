// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVAR_UIESTIMATOR_H
#define COLVAR_UIESTIMATOR_H

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include <typeinfo>

// only for colvar module!
// when integrated into other code, just remove this line and "...cvm::backup_file(...)"
#include "colvarmodule.h"

namespace UIestimator {
    const int Y_SIZE = 21;            // defines the range of extended CV with respect to a given CV
                                      // For example, CV=10, width=1, Y_SIZE=21, then eCV=[0-20], having a size of 21
    const int HALF_Y_SIZE = 10;
    const int EXTENDED_X_SIZE = HALF_Y_SIZE;
    const double EPSILON = 0.000001;   // for comparison of float numbers

    class n_matrix {   // Stores the distribution matrix of n(x,y)

    public:
        n_matrix() {}
        n_matrix(const std::vector<double> & lowerboundary,   // lowerboundary of x
            const std::vector<double> & upperboundary,   // upperboundary of
            const std::vector<double> & width,           // width of x
            const int y_size) {          // size of y, for example, ysize=7, then when x=1, the distribution of y in [-2,4] is considered

            int i;

            this->lowerboundary = lowerboundary;
            this->upperboundary = upperboundary;
            this->width = width;
            this->dimension = lowerboundary.size();
            this->y_size = y_size;     // keep in mind the internal (spare) matrix is stored in diagonal form
            this->y_total_size = int(std::pow(double(y_size), double(dimension)) + EPSILON);

            // the range of the matrix is [lowerboundary, upperboundary]
            x_total_size = 1;
            for (i = 0; i < dimension; i++) {
                x_size.push_back(int((upperboundary[i] - lowerboundary[i]) / width[i] + EPSILON));
                x_total_size *= x_size[i];
            }

            // initialize the internal matrix
            matrix.reserve(x_total_size);
            for (i = 0; i < x_total_size; i++) {
                matrix.push_back(std::vector<int>(y_total_size, 0));
            }

            temp.resize(dimension);
        }

        int inline get_value(const std::vector<double> & x, const std::vector<double> & y) {
            return matrix[convert_x(x)][convert_y(x, y)];
        }

        void inline set_value(const std::vector<double> & x, const std::vector<double> & y, const int value) {
            matrix[convert_x(x)][convert_y(x,y)] = value;
        }

        void inline increase_value(const std::vector<double> & x, const std::vector<double> & y, const int value) {
            matrix[convert_x(x)][convert_y(x,y)] += value;
        }

    private:
        std::vector<double> lowerboundary;
        std::vector<double> upperboundary;
        std::vector<double> width;
        int dimension;
        std::vector<int> x_size;       // the size of x in each dimension
        int x_total_size;              // the size of x of the internal matrix
        int y_size;                    // the size of y in each dimension
        int y_total_size;              // the size of y of the internal matrix

        std::vector<std::vector<int> > matrix;  // the internal matrix

        std::vector<int> temp;         // this vector is used in convert_x and convert_y to save computational resource

        int i, j;

        int convert_x(const std::vector<double> & x) {       // convert real x value to its interal index
            for (i = 0; i < dimension; i++) {
                temp[i] = int((x[i] - lowerboundary[i]) / width[i] + EPSILON);
            }

            int index = 0;
            for (i = 0; i < dimension; i++) {
                if (i + 1 < dimension) {
                    int x_temp = 1;
                    for (j = i + 1; j < dimension; j++)
                        x_temp *= x_size[j];
                    index += temp[i] * x_temp;
                }
                else
                    index += temp[i];
            }
            return index;
        }

        int convert_y(const std::vector<double> & x, const std::vector<double> & y) {       // convert real y value to its interal index

            int i;

            for (i = 0; i < dimension; i++) {
                temp[i] = round((round(y[i] / width[i] + EPSILON) - round(x[i] / width[i] + EPSILON)) + (y_size - 1) / 2 + EPSILON);
            }

            int index = 0;
            for (i = 0; i < dimension; i++) {
                if (i + 1 < dimension)
                    index += temp[i] * int(std::pow(double(y_size), double(dimension - i - 1)) + EPSILON);
                else
                    index += temp[i];
            }
            return index;
        }

        double round(double r) {
            return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
        }
    };

    // vector, store the sum_x, sum_x_square, count_y
    template <typename T>
    class n_vector {

    public:
        n_vector() {}
        n_vector(const std::vector<double> & lowerboundary,   // lowerboundary of x
            const std::vector<double> & upperboundary,   // upperboundary of
            const std::vector<double> & width,                // width of x
            const int y_size,           // size of y, for example, ysize=7, then when x=1, the distribution of y in [-2,4] is considered
            const T & default_value) {         //   the default value of T

            this->width = width;
            this->dimension = lowerboundary.size();

            x_total_size = 1;
            for (int i = 0; i < dimension; i++) {
                this->lowerboundary.push_back(lowerboundary[i] - (y_size - 1) / 2 * width[i] - EPSILON);
                this->upperboundary.push_back(upperboundary[i] + (y_size - 1) / 2 * width[i] + EPSILON);

                x_size.push_back(int((this->upperboundary[i] - this->lowerboundary[i]) / this->width[i] + EPSILON));
                x_total_size *= x_size[i];
            }

            // initialize the internal vector
            vector.resize(x_total_size, default_value);

            temp.resize(dimension);
        }

        const T inline get_value(const std::vector<double> & x) {
            return vector[convert_x(x)];
        }

        void inline set_value(const std::vector<double> & x, const T value) {
            vector[convert_x(x)] = value;
        }

        void inline increase_value(const std::vector<double> & x, const T value) {
            vector[convert_x(x)] += value;
        }
    private:
        std::vector<double> lowerboundary;
        std::vector<double> upperboundary;
        std::vector<double> width;
        int dimension;
        std::vector<int> x_size;       // the size of x in each dimension
        int x_total_size;              // the size of x of the internal matrix

        std::vector<T> vector;  // the internal vector

        std::vector<int> temp;         // this vector is used in convert_x and convert_y to save computational resource

        int convert_x(const std::vector<double> & x) {       // convert real x value to its interal index

            int i, j;

            for (i = 0; i < dimension; i++) {
                temp[i] = int((x[i] - lowerboundary[i]) / width[i] + EPSILON);
            }

            int index = 0;
            for (i = 0; i < dimension; i++) {
                if (i + 1 < dimension) {
                    int x_temp = 1;
                    for (j = i + 1; j < dimension; j++)
                        x_temp *= x_size[j];
                    index += temp[i] * x_temp;
                }
                else
                    index += temp[i];
            }
            return index;
        }
    };

    class UIestimator {     // the implemension of UI estimator

    public:
        UIestimator() {}

        //called when (re)start an eabf simulation
        UIestimator(const std::vector<double> & lowerboundary,
            const std::vector<double> & upperboundary,
            const std::vector<double> & width,
            const std::vector<double> & krestr,                // force constant in eABF
            const std::string & output_filename,              // the prefix of output files
            const int output_freq,
            const bool restart,                              // whether restart from a .count and a .grad file
            const std::vector<std::string> & input_filename,   // the prefixes of input files
            const double temperature) {

            // initialize variables
            this->lowerboundary = lowerboundary;
            this->upperboundary = upperboundary;
            this->width = width;
            this->krestr = krestr;
            this->output_filename = output_filename;
            this->output_freq = output_freq;
            this->restart = restart;
            this->input_filename = input_filename;
            this->temperature = temperature;

            int i, j;

            dimension = lowerboundary.size();

            for (i = 0; i < dimension; i++) {
                sum_x.push_back(n_vector<double>(lowerboundary, upperboundary, width, Y_SIZE, 0.0));
                sum_x_square.push_back(n_vector<double>(lowerboundary, upperboundary, width, Y_SIZE, 0.0));

                x_av.push_back(n_vector<double>(lowerboundary, upperboundary, width, Y_SIZE, 0.0));
                sigma_square.push_back(n_vector<double>(lowerboundary, upperboundary, width, Y_SIZE, 0.0));
            }

            count_y = n_vector<int>(lowerboundary, upperboundary, width, Y_SIZE, 0);
            distribution_x_y = n_matrix(lowerboundary, upperboundary, width, Y_SIZE);

            grad = n_vector<std::vector<double> >(lowerboundary, upperboundary, width, 1, std::vector<double>(dimension, 0.0));
            count = n_vector<int>(lowerboundary, upperboundary, width, 1, 0);

            written = false;
            written_1D = false;

            if (dimension == 1) {
                std::vector<double> upperboundary_temp = upperboundary;
                upperboundary_temp[0] = upperboundary[0] + width[0];
                oneD_pmf = n_vector<double>(lowerboundary, upperboundary_temp, width, 1, 0.0);
            }

            if (restart == true) {
                input_grad = n_vector<std::vector<double> >(lowerboundary, upperboundary, width, 1, std::vector<double>(dimension, 0.0));
                input_count = n_vector<int>(lowerboundary, upperboundary, width, 1, 0);

                // initialize input_Grad and input_count
                // the loop_flag is a n-dimensional vector, increae from lowerboundary to upperboundary when looping
                std::vector<double> loop_flag(dimension, 0);
                for (i = 0; i < dimension; i++) {
                    loop_flag[i] = lowerboundary[i];
                }

                i = 0;
                while (i >= 0) {
                    for (j = 0; j < dimension; j++) {
                        input_grad.set_value(loop_flag, std::vector<double>(dimension,0));
                    }
                    input_count.set_value(loop_flag, 0);

                    // iterate over any dimensions
                    i = dimension - 1;
                    while (i >= 0) {
                        loop_flag[i] += width[i];
                        if (loop_flag[i] > upperboundary[i] - width[i] + EPSILON) {
                            loop_flag[i] = lowerboundary[i];
                            i--;
                        }
                        else
                            break;
                    }
                }
                read_inputfiles(input_filename);
            }
        }

        ~UIestimator() {}

        // called from MD engine every step
        bool update(const int step, std::vector<double> x, std::vector<double> y) {

            int i;

            if (step % output_freq == 0) {
                calc_pmf();
                write_files();
                //write_interal_data();
            }

            for (i = 0; i < dimension; i++) {
                // for dihedral RC, it is possible that x = 179 and y = -179, should correct it
                // may have problem, need to fix
                if (x[i] > 150 && y[i] < -150) {
                    y[i] += 360;
                }
                if (x[i] < -150 && y[i] > 150) {
                    y[i] -= 360;
                }

                if (x[i] < lowerboundary[i] - EXTENDED_X_SIZE * width[i] + EPSILON || x[i] > upperboundary[i] + EXTENDED_X_SIZE * width[i] - EPSILON \
                    || y[i] - x[i] < -HALF_Y_SIZE * width[i] + EPSILON || y[i] - x[i] > HALF_Y_SIZE * width[i] - EPSILON \
                    || y[i] - lowerboundary[i] < -HALF_Y_SIZE * width[i] + EPSILON || y[i] - upperboundary[i] > HALF_Y_SIZE * width[i] - EPSILON)
                    return false;
            }

            for (i = 0; i < dimension; i++) {
                sum_x[i].increase_value(y, x[i]);
                sum_x_square[i].increase_value(y, x[i] * x[i]);
            }
            count_y.increase_value(y, 1);

            for (i = 0; i < dimension; i++) {
                // adapt colvars precision
                if (x[i] < lowerboundary[i] + EPSILON || x[i] > upperboundary[i] - EPSILON)
                    return false;
            }
            distribution_x_y.increase_value(x, y, 1);

            return true;
        }

        // update the output_filename
        void update_output_filename(const std::string& filename) {
            output_filename = filename;
        }

    private:
        std::vector<n_vector<double> > sum_x;                        // the sum of x in each y bin
        std::vector<n_vector<double> > sum_x_square;                 // the sum of x in each y bin
        n_vector<int> count_y;                              // the distribution of y
        n_matrix distribution_x_y;   // the distribution of <x, y> pair

        int dimension;

        std::vector<double> lowerboundary;
        std::vector<double> upperboundary;
        std::vector<double> width;
        std::vector<double> krestr;
        std::string output_filename;
        int output_freq;
        bool restart;
        std::vector<std::string> input_filename;
        double temperature;

        n_vector<std::vector<double> > grad;
        n_vector<int> count;

        n_vector<double> oneD_pmf;

        n_vector<std::vector<double> > input_grad;
        n_vector<int> input_count;

        // used in double integration
        std::vector<n_vector<double> > x_av;
        std::vector<n_vector<double> > sigma_square;

        bool written;
        bool written_1D;

        // calculate gradients from the internal variables
        void calc_pmf() {
            int norm;
            int i, j, k;

            std::vector<double> loop_flag(dimension, 0);
            for (i = 0; i < dimension; i++) {
                loop_flag[i] = lowerboundary[i] - HALF_Y_SIZE * width[i];
            }

            i = 0;
            while (i >= 0) {
                norm = count_y.get_value(loop_flag) > 0 ? count_y.get_value(loop_flag) : 1;
                for (j = 0; j < dimension; j++) {
                    x_av[j].set_value(loop_flag, sum_x[j].get_value(loop_flag) / norm);
                    sigma_square[j].set_value(loop_flag, sum_x_square[j].get_value(loop_flag) / norm - x_av[j].get_value(loop_flag) * x_av[j].get_value(loop_flag));
                }

                // iterate over any dimensions
                i = dimension - 1;
                while (i >= 0) {
                    loop_flag[i] += width[i];
                    if (loop_flag[i] > upperboundary[i] + HALF_Y_SIZE * width[i] - width[i] + EPSILON) {
                        loop_flag[i] = lowerboundary[i] - HALF_Y_SIZE * width[i];
                        i--;
                    }
                    else
                        break;
                }
            }

            // double integration
            std::vector<double> av(dimension, 0);
            std::vector<double> diff_av(dimension, 0);

            std::vector<double> loop_flag_x(dimension, 0);
            std::vector<double> loop_flag_y(dimension, 0);
            for (i = 0; i < dimension; i++) {
                loop_flag_x[i] = lowerboundary[i];
                loop_flag_y[i] = loop_flag_x[i] - HALF_Y_SIZE * width[i];
            }

            i = 0;
            while (i >= 0) {
                norm = 0;
                for (k = 0; k < dimension; k++) {
                    av[k] = 0;
                    diff_av[k] = 0;
                    loop_flag_y[k] = loop_flag_x[k] - HALF_Y_SIZE * width[k];
                }

                int j = 0;
                while (j >= 0) {
                    norm += distribution_x_y.get_value(loop_flag_x, loop_flag_y);
                    for (k = 0; k < dimension; k++) {
                        if (sigma_square[k].get_value(loop_flag_y) > EPSILON || sigma_square[k].get_value(loop_flag_y) < -EPSILON)
                            av[k] += distribution_x_y.get_value(loop_flag_x, loop_flag_y) * ( (loop_flag_x[k] + 0.5 * width[k]) - x_av[k].get_value(loop_flag_y)) / sigma_square[k].get_value(loop_flag_y);

                        diff_av[k] += distribution_x_y.get_value(loop_flag_x, loop_flag_y) * (loop_flag_x[k] - loop_flag_y[k]);
                    }

                    // iterate over any dimensions
                    j = dimension - 1;
                    while (j >= 0) {
                        loop_flag_y[j] += width[j];
                        if (loop_flag_y[j] > loop_flag_x[j] + HALF_Y_SIZE * width[j] - width[j] + EPSILON) {
                            loop_flag_y[j] = loop_flag_x[j] - HALF_Y_SIZE * width[j];
                            j--;
                        }
                        else
                            break;
                    }
                }

                std::vector<double> grad_temp(dimension, 0);
                for (k = 0; k < dimension; k++) {
                    diff_av[k] /= (norm > 0 ? norm : 1);
                    av[k] = cvm::boltzmann() * temperature * av[k] / (norm > 0 ? norm : 1);
                    grad_temp[k] = av[k] - krestr[k] * diff_av[k];
                }
                grad.set_value(loop_flag_x, grad_temp);
                count.set_value(loop_flag_x, norm);

                // iterate over any dimensions
                i = dimension - 1;
                while (i >= 0) {
                    loop_flag_x[i] += width[i];
                    if (loop_flag_x[i] > upperboundary[i] - width[i] + EPSILON) {
                        loop_flag_x[i] = lowerboundary[i];
                        i--;
                    }
                    else
                        break;
                }
            }
        }


        // calculate 1D pmf
        void calc_1D_pmf()
        {
            std::vector<double> last_position(1, 0);
            std::vector<double> position(1, 0);

            double min = 0;
            double dG = 0;
            double i;

            oneD_pmf.set_value(lowerboundary, 0);
            last_position = lowerboundary;
            for (i = lowerboundary[0] + width[0]; i < upperboundary[0] + EPSILON; i += width[0]) {
                position[0] = i + EPSILON;
                if (restart == false || input_count.get_value(last_position) == 0) {
                    dG = oneD_pmf.get_value(last_position) + grad.get_value(last_position)[0] * width[0];
                }
                else {
                    dG = oneD_pmf.get_value(last_position) + ((grad.get_value(last_position)[0] * count.get_value(last_position) + input_grad.get_value(last_position)[0] * input_count.get_value(last_position)) / (count.get_value(last_position) + input_count.get_value(last_position))) * width[0];
                }
                if (dG < min)
                    min = dG;
                oneD_pmf.set_value(position, dG);
                last_position[0] = i + EPSILON;
            }

            for (i = lowerboundary[0]; i < upperboundary[0] + EPSILON; i += width[0]) {
                position[0] = i + EPSILON;
                oneD_pmf.set_value(position, oneD_pmf.get_value(position) - min);
            }
        }

        // write 1D pmf
        void write_1D_pmf() {
            std::string pmf_filename = output_filename + ".UI.pmf";

            // only for colvars module!
            if (written_1D) cvm::backup_file(pmf_filename.c_str());

            std::ostream* ofile_pmf = cvm::proxy->output_stream(pmf_filename.c_str());

            std::vector<double> position(1, 0);
            for (double i = lowerboundary[0]; i < upperboundary[0] + EPSILON; i += width[0]) {
                *ofile_pmf << i << " ";
                position[0] = i + EPSILON;
                *ofile_pmf << oneD_pmf.get_value(position) << std::endl;
            }
            cvm::proxy->close_output_stream(pmf_filename.c_str());

            written_1D = true;
        }

        // write heads of the output files
        void writehead(std::ostream& os) const {
            os << "# " << dimension << std::endl;
            for (int i = 0; i < dimension; i++) {
                os << "# " << lowerboundary[i] << " " << width[i] << " " << int((upperboundary[i] - lowerboundary[i]) / width[i] + EPSILON) << " " << 0 << std::endl;
            }
            os << std::endl;
        }

        // write interal data, used for testing
        void write_interal_data() {
            std::string internal_filename = output_filename + ".UI.internal";

            std::ostream* ofile_internal = cvm::proxy->output_stream(internal_filename.c_str());

            std::vector<double> loop_flag(dimension, 0);
            for (int i = 0; i < dimension; i++) {
                loop_flag[i] = lowerboundary[i];
            }

            int n = 0;
            while (n >= 0) {
                for (int j = 0; j < dimension; j++) {
                    *ofile_internal << loop_flag[j] + 0.5 * width[j] << " ";
                }

                for (int k = 0; k < dimension; k++) {
                    *ofile_internal << grad.get_value(loop_flag)[k] << " ";
                }

                std::vector<double> ii(dimension,0);
                for (double i = loop_flag[0] - 10; i < loop_flag[0] + 10 + EPSILON; i+= width[0]) {
                    for (double j = loop_flag[1] - 10; j< loop_flag[1] + 10 + EPSILON; j+=width[1]) {
                        ii[0] = i;
                        ii[1] = j;
                        *ofile_internal << i <<" "<<j<<" "<< distribution_x_y.get_value(loop_flag,ii)<< " ";
                    }
                }
                *ofile_internal << std::endl;

                // iterate over any dimensions
                n = dimension - 1;
                while (n >= 0) {
                    loop_flag[n] += width[n];
                    if (loop_flag[n] > upperboundary[n] - width[n] + EPSILON) {
                        loop_flag[n] = lowerboundary[n];
                        n--;
                    }
                    else
                        break;
                }
            }
            cvm::proxy->close_output_stream(internal_filename.c_str());
        }

        // write output files
        void write_files() {
            std::string grad_filename = output_filename + ".UI.grad";
            std::string hist_filename = output_filename + ".UI.hist.grad";
            std::string count_filename = output_filename + ".UI.count";

            int i, j;
//
            // only for colvars module!
            if (written) cvm::backup_file(grad_filename.c_str());
            //if (written) cvm::backup_file(hist_filename.c_str());
            if (written) cvm::backup_file(count_filename.c_str());

            std::ostream* ofile = cvm::proxy->output_stream(grad_filename.c_str());
            std::ostream* ofile_hist = cvm::proxy->output_stream(hist_filename.c_str(), std::ios::app);
            std::ostream* ofile_count = cvm::proxy->output_stream(count_filename.c_str());

            writehead(*ofile);
            writehead(*ofile_hist);
            writehead(*ofile_count);

            if (dimension == 1) {
                calc_1D_pmf();
                write_1D_pmf();
            }

            std::vector<double> loop_flag(dimension, 0);
            for (i = 0; i < dimension; i++) {
                loop_flag[i] = lowerboundary[i];
            }

            i = 0;
            while (i >= 0) {
                for (j = 0; j < dimension; j++) {
                    *ofile << loop_flag[j] + 0.5 * width[j] << " ";
                    *ofile_hist << loop_flag[j] + 0.5 * width[j] << " ";
                    *ofile_count << loop_flag[j] + 0.5 * width[j] << " ";
                }

                if (restart == false) {
                    for (j = 0; j < dimension; j++) {
                        *ofile << grad.get_value(loop_flag)[j] << " ";
                        *ofile_hist << grad.get_value(loop_flag)[j] << " ";
                    }
                    *ofile << std::endl;
                    *ofile_hist << std::endl;
                    *ofile_count << count.get_value(loop_flag) << " " <<std::endl;
                }
                else {
                    double final_grad = 0;
                    for (j = 0; j < dimension; j++) {
                        int total_count_temp = (count.get_value(loop_flag) + input_count.get_value(loop_flag));
                        if (input_count.get_value(loop_flag) == 0)
                            final_grad = grad.get_value(loop_flag)[j];
                        else
                            final_grad = ((grad.get_value(loop_flag)[j] * count.get_value(loop_flag) + input_grad.get_value(loop_flag)[j] * input_count.get_value(loop_flag)) / total_count_temp);
                        *ofile << final_grad << " ";
                        *ofile_hist << final_grad << " ";
                    }
                    *ofile << std::endl;
                    *ofile_hist << std::endl;
                    *ofile_count << (count.get_value(loop_flag) + input_count.get_value(loop_flag)) << " " <<std::endl;
                }

                // iterate over any dimensions
                i = dimension - 1;
                while (i >= 0) {
                    loop_flag[i] += width[i];
                    if (loop_flag[i] > upperboundary[i] - width[i] + EPSILON) {
                        loop_flag[i] = lowerboundary[i];
                        i--;
                        *ofile << std::endl;
                        *ofile_hist << std::endl;
                        *ofile_count << std::endl;
                    }
                    else
                        break;
                }
            }
            cvm::proxy->close_output_stream(grad_filename.c_str());
            cvm::proxy->close_output_stream(hist_filename.c_str());
            cvm::proxy->close_output_stream(count_filename.c_str());

            written = true;
        }

        // read input files
        void read_inputfiles(const std::vector<std::string> input_filename)
        {
            char sharp;
            double nothing;
            int dimension_temp;
            int i, j, k, l, m;

            std::vector<double> loop_bin_size(dimension, 0);
            std::vector<double> position_temp(dimension, 0);
            std::vector<double> grad_temp(dimension, 0);
            int count_temp = 0;
            for (i = 0; i < int(input_filename.size()); i++) {
                int size = 1 , size_temp = 0;

                std::string count_filename = input_filename[i] + ".UI.count";
                std::string grad_filename = input_filename[i] + ".UI.grad";

                std::ifstream count_file(count_filename.c_str(), std::ios::in);
                std::ifstream grad_file(grad_filename.c_str(), std::ios::in);

                count_file >> sharp >> dimension_temp;
                grad_file >> sharp >> dimension_temp;

                for (j = 0; j < dimension; j++) {
                    count_file >> sharp >> nothing >> nothing >> size_temp >> nothing;
                    grad_file >> sharp >> nothing >> nothing >> nothing >> nothing;
                    size *= size_temp;
                }

                for (j = 0; j < size; j++) {
                    do {
                        for (k = 0; k < dimension; k++) {
                            count_file >> position_temp[k];
                            grad_file >> nothing;
                        }

                        for (l = 0; l < dimension; l++) {
                            grad_file >> grad_temp[l];
                        }
                        count_file >> count_temp;
                    }
                    while (position_temp[i] < lowerboundary[i] - EPSILON || position_temp[i] > upperboundary[i] + EPSILON);

                    if (count_temp == 0) {
                        continue;
                    }

                    for (m = 0; m < dimension; m++) {
                        grad_temp[m] = (grad_temp[m] * count_temp + input_grad.get_value(position_temp)[m] * input_count.get_value(position_temp)) / (count_temp + input_count.get_value(position_temp));
                    }
                    input_grad.set_value(position_temp, grad_temp);
                    input_count.increase_value(position_temp, count_temp);
                }

                count_file.close();
                grad_file.close();
            }
        }
    };
}

#endif
