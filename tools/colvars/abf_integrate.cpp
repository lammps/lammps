/****************************************************************
 * abf_integrate                                                *
 * Integrate n-dimensional PMF from discrete gradient grid      *
 * Jerome Henin <jerome.henin@ibpc.fr>                          *
 ****************************************************************/

#include "abf_data.h"
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>

char *parse_cl(int argc, char *argv[], unsigned int *nsteps, double *temp,
               bool * meta, double *hill, double *hill_fact);
double compute_deviation(ABFdata * data, bool meta, double kT);

int main(int argc, char *argv[])
{
    char *data_file;
    char *out_file;
    unsigned int step, nsteps, total, out_freq;
    int *pos, *dpos, *newpos;
    unsigned int offset, newoffset;
    double dA;
    double temp;
    double mbeta;
    bool meta;
    double hill, hill_fact, hill_min;
    double rmsd, rmsd_old, rmsd_rel_change, convergence_limit;
    bool converged;
    unsigned int scale_hill_step;

    // Setting default values
    nsteps = 0;
    temp = 500;
    meta = true;
    hill = 0.01;
    hill_fact = 0.5;
    hill_min = 0.0005;

    convergence_limit = -0.001;

    if (!(data_file = parse_cl(argc, argv, &nsteps, &temp, &meta, &hill, &hill_fact))) {
        std::cerr << "\nabf_integrate: MC-based integration of multidimensional free energy gradient\n";
        std::cerr << "Version 20160420\n\n";
        std::cerr << "Syntax: " << argv[0] <<
            " <filename> [-n <nsteps>] [-t <temp>] [-m [0|1] (metadynamics)]"
            " [-h <hill_height>] [-f <variable_hill_factor>]\n\n";
        exit(1);
    }

    if (meta) {
        std::cout << "\nUsing metadynamics-style sampling with hill height: " << hill << "\n";
        if (hill_fact) {
            std::cout << "Varying hill height by factor " << hill_fact << "\n";
        }
    } else {
        std::cout << "\nUsing unbiased MC sampling\n";
    }

    if (nsteps) {
        std::cout << "Sampling " << nsteps << " steps at temperature " << temp << "\n\n";
        out_freq = nsteps / 10;
        scale_hill_step = nsteps / 2;
        converged = true;
    } else {
        std::cout << "Sampling until convergence at temperature " << temp << "\n\n";
        out_freq = 1000000;
        converged = false;
    }

    // Inverse temperature in (kcal/mol)-1
    mbeta = -1 / (0.001987 * temp);

    ABFdata data(data_file);

    if (!nsteps) {
        scale_hill_step = 2000 * data.scalar_dim;
        nsteps = 2 * scale_hill_step;
        std::cout << "Setting minimum number of steps to " << nsteps << "\n";
    }

    srand(time(NULL));

    pos = new int[data.Nvars];
    dpos = new int[data.Nvars];
    newpos = new int[data.Nvars];

    // TODO: this will be an infinite loop if there is no sampling
    // it would be more robust to build a list of non-empty bins, and
    // pick from there (or just treat the special case of no sampling
    // and output a null PMF)
    do {
      for (int i = 0; i < data.Nvars; i++) {
          pos[i] = rand() % data.sizes[i];
      }
      offset = data.offset(pos);
    } while ( !data.allowed (offset) );

    rmsd = compute_deviation(&data, meta, 0.001987 * temp);
    std::cout << "\nInitial gradient RMS is " << rmsd << "\n";

    total = 0;
    for (step = 1; (step <= nsteps || !converged); step++) {

        if ( step % out_freq == 0) {
            rmsd_old = rmsd;
            rmsd = compute_deviation(&data, meta, 0.001987 * temp);
            rmsd_rel_change = (rmsd - rmsd_old) / (rmsd_old * double (out_freq)) * 1000000.0;
            std::cout << "Step " << step << " ; gradient RMSD is " << rmsd
                      << " ; relative change per 1M steps " << rmsd_rel_change;
            if ( rmsd_rel_change > convergence_limit && step >= nsteps ) {
                converged = true;
            }

            if (meta && hill_fact && step > scale_hill_step && hill > hill_min ) {
                hill *= hill_fact;
                std::cout << " - changing hill height to " << hill << "\n";
            } else {
                std::cout << "\n";
            }
        }

        offset = data.offset(pos);
        data.histogram[offset]++;
        if (meta) {
            data.bias[offset] += hill;
        }

        const double *grad = data.gradients + offset * data.Nvars;

        int not_accepted = 1;
        while (not_accepted) {
            dA = 0.0;
            total++;
            for (int i = 0; i < data.Nvars; i++) {
                dpos[i] = rand() % 3 - 1;
                newpos[i] = pos[i] + dpos[i];
                data.wrap(newpos[i], i);
                if (newpos[i] == pos[i])
                    dpos[i] = 0;

                if (dpos[i]) {
                    dA += grad[i] * dpos[i] * data.widths[i];
                    // usefulness of the interpolation below depends on
                    // where the grid points are for the histogram wrt to the gradients
                    // If done, it has to be done in all directions
                    // the line below is useless
                    //dA += 0.5 * (newgrad[i] + grad[i]) * dpos[i] * data.widths[i];
                }
            }

            newoffset = data.offset(newpos);
            if (meta) {
                dA += data.bias[newoffset] - data.bias[offset];
            }

            if ( data.allowed (newoffset) && (((float) rand()) / RAND_MAX < exp(mbeta * dA)) )  {
                // Accept move
                for (int i = 0; i < data.Nvars; i++) {
                    pos[i] = newpos[i];
                    not_accepted = 0;
                }
            }
        }
    }
    std::cout << "Run " << total << " total iterations; acceptance ratio is "
        << double (step) / double (total)
        << " ; final gradient RMSD is " << compute_deviation(&data, meta, 0.001987 * temp) << "\n";

    out_file = new char[strlen(data_file) + 8];

    if (meta) {
        sprintf(out_file, "%s.pmf", data_file);
        std::cout << "Writing PMF to file " << out_file << "\n";
        data.write_bias(out_file);
    }

    // TODO write a PMF for unbiased MC, too...
    sprintf(out_file, "%s.histo", data_file);
    std::cout << "Writing sampling histogram to file " << out_file << "\n";
    data.write_histogram(out_file);

    sprintf(out_file, "%s.est", data_file);
    std::cout << "Writing estimated FE gradient to file " << out_file << "\n";
    data.write_field(data.estimate, out_file);

    sprintf(out_file, "%s.dev", data_file);
    std::cout << "Writing FE gradient deviation to file " << out_file << "\n\n";
    data.write_field(data.deviation, out_file);

    delete [] pos;
    delete [] dpos;
    delete [] newpos;
    delete [] out_file;
    exit(0);
}


double compute_deviation(ABFdata * data, bool meta, double kT)
{
    // Computing deviation between gradients differentiated from pmf
    // and input data
    // NOTE: this is mostly for output, hence NOT performance-critical
    double        *dev = data->deviation;
    double        *est = data->estimate;
    const double  *grad = data->gradients;
    int           *pos, *newpos;
    double        rmsd = 0.0;
    unsigned int  offset, newoffset;
    double        sum;
    int           c;
    bool          moved;
    unsigned int  norm = 0; // number of data points summmed

    pos = new int[data->Nvars];
    newpos = new int[data->Nvars];

    for (int i = 0; i < data->Nvars; i++)
        pos[i] = 0;

    for (offset = 0; offset < data->scalar_dim; offset++) {
        for (int i = data->Nvars - 1; i > 0; i--) {
            if (pos[i] == data->sizes[i]) {
                pos[i] = 0;
                pos[i - 1]++;
            }
        }

        if (data->allowed (offset)) {
          for (int i = 0; i < data->Nvars; i++)
              newpos[i] = pos[i];

          for (int i = 0; i < data->Nvars; i++) {
              est[i] = 0.0;
              sum = 0.0;          // sum of finite differences on two sides (if not on edge of the grid)
              c = 0;              // count of summed values

              newpos[i] = pos[i] - 1;
              moved = data->wrap(newpos[i], i);
              newoffset = data->offset(newpos);
              if ( moved && data->allowed (newoffset) ) {
                  if (meta) {
                      sum = (data->bias[newoffset] - data->bias[offset]) / data->widths[i];
                      c++;
                  } else {
                      if (data->histogram[offset] && data->histogram[newoffset]) {
                          sum = kT * log(double (data->histogram[newoffset]) /
                                            double (data->histogram[offset])) / data->widths[i];
                          c++;
                      }
                  }
              }

              newpos[i] = pos[i] + 1;
              moved = data->wrap(newpos[i], i);
              newoffset = data->offset(newpos);
              if ( moved && data->allowed (newoffset) ) {
                  if (meta) {
                      sum += (data->bias[offset] - data->bias[newoffset]) / data->widths[i];
                      c++;
                  } else {
                      if (data->histogram[offset] && data->histogram[newoffset]) {
                          sum += kT * log(double (data->histogram[offset]) /
                                             double (data->histogram[newoffset])) / data->widths[i];
                          c++;
                      }
                  }
              }

              newpos[i] = pos[i]; // Go back to initial position for next dimension

              est[i] = (c ? sum/double(c) : 0.0);
              dev[i] = grad[i] - est[i];
              rmsd += dev[i] * dev[i];
              norm++;
          }
        }

        pos[data->Nvars - 1]++; // move on to next point
        est += data->Nvars;
        dev += data->Nvars;
        grad += data->Nvars;
    }

    delete [] pos;
    delete [] newpos;

    return sqrt(rmsd / norm);
}


char *parse_cl(int argc, char *argv[], unsigned int *nsteps, double *temp,
               bool * meta, double *hill, double *hill_fact)
{
    int meta_int;

    // getting default value for the integer
    meta_int = (*meta ? 1 : 0);

    // "Syntax: " << argv[0] << " <filename> [-n <nsteps>] [-t <temp>] [-m [0|1] (metadynamics)] [-h <hill_height>]\n";
    if (argc < 2) {
        return NULL;
    }

    for (int i = 2; i + 1 < argc; i += 2) {
        if (argv[i][0] != '-') {
            return NULL;
        }
        switch (argv[i][1]) {
        case 'n':
            if (sscanf(argv[i + 1], "%u", nsteps) != 1)
                return NULL;
            break;
        case 't':
            if (sscanf(argv[i + 1], "%lf", temp) != 1)
                return NULL;
            break;
        case 'm':
            if (sscanf(argv[i + 1], "%u", &meta_int) != 1)
                return NULL;
            break;
        case 'h':
            if (sscanf(argv[i + 1], "%lf", hill) != 1)
                return NULL;
            break;
        case 'f':
            if (sscanf(argv[i + 1], "%lf", hill_fact) != 1)
                return NULL;
            break;
        default:
            return NULL;
        }
    }

    *meta = (meta_int != 0);
    return argv[1];
}
