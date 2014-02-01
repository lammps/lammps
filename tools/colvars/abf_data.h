/// \file integrate.h General headers for ABF_integrate

#include <iostream>
#include <vector>

#define MIN_SAMPLES 1

/// Free energy gradients class
class ABFdata {

  protected:
    /// Sizes of (i-1) dimension blocks
    /// computed as Prod_(j<i) sizes[j]
    int *blocksizes;
    /// Minimum values of each variable
    double *mins;

  public:
    int Nvars;
    /// Free energy gradients (vector field)
    double *gradients;
    /// Sampling from the ABF calculation
    unsigned int *count;
    /// Bin widths
    double *widths;

    unsigned int scalar_dim;
    unsigned int vec_dim;
    unsigned int *histogram;

    /// History-dependent bias
    double *bias;

    /// Estimate of the FE gradient computed
    /// from MtD bias or histogram in standard MC
    double *estimate;

    /// Deviation between starting free energy gradient and 
    /// estimated one
    double *deviation;

    void write_histogram(const char *fileName);
    void write_bias(const char *fileName);
    void write_field(double *field, const char *fileName);

    /// Grid sizes
    int *sizes;

    /// Flag stating if each variable is periodic
    int *PBC;

    /// Constructor: reads from a file
     ABFdata(const char *gradFileName);
    ~ABFdata();

    /// \brief Returns an offset for scalar fields based on a n-index.
    /// multiply by Nvars to get an offset in a Nvars-vector field
    unsigned int offset(const int *);

    inline bool wrap(int &pos, int i);

    /// Decides if an offset is outside the allowed region based on the ABF sampling
    inline bool allowed(unsigned int offset);
};


inline bool ABFdata::wrap(int &pos, int i)
{
    if (PBC[i]) {
        if (pos == -1) {
            pos = sizes[i] - 1;
            return true;
        }
        if (pos == sizes[i]) {
            pos = 0;
            return true;
        }
    } else {
        // No PBC
        if (pos == -1) {
            pos = 0;
            return false;
        }
        if (pos == sizes[i]) {
            pos = sizes[i] - 1;
            return false;
        }
    }
    return true;
}

inline bool ABFdata::allowed(unsigned int offset) {
    return count[offset] > MIN_SAMPLES;
}
