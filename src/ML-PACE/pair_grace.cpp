//
// Created by Yury Lysogorskiy on 01.12.23.
//
#ifndef NO_GRACE_TF

#include "pair_grace.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"

#include <cstring>
#include <exception>
#include <numeric>
#include "yaml-cpp/yaml.h"

#include "ace-evaluator/ace_arraynd.h"

#include "utils_pace.h"

#include <cppflow/tensor.h>

#include <tensorflow/c/c_api.h>

const std::string DEFAULT_INPUT_PREFIX = "serving_default_";


namespace LAMMPS_NS {

    struct ACETPImpl {
        ACETPImpl() : model(nullptr) {}

        ~ACETPImpl() {
            delete model;
        }

        cppflow::model *model;
    };

    bool check_tf_graph_input_presented(cppflow::model *model, const std::string &op_name) {
        try {
            auto op_shape = model->get_operation_shape(op_name);
            return true;
        } catch (std::runtime_error &exc) {
            return false;
        }
    }


}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */
PairGRACE::PairGRACE(LAMMPS *lmp) : Pair(lmp) {
    single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;

    aceimpl = new ACETPImpl;

    scale = nullptr;

    chunksize = 4096;

    data_timer.init();
    tp_timer.init();

    no_virial_fdotr_compute = 1;
}


/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */
PairGRACE::~PairGRACE() {
    if (copymode) return;

    delete aceimpl;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(scale);
    }
    auto data_t = (double) data_timer.as_microseconds();
    auto tp_t = (double) tp_timer.as_microseconds();

    utils::logmesg(lmp,
                   "[GRACE:debug, proc #{:d}]: Data preparation timer: {:g} mcs, graph execution time: {:g} mcs, data preparation time fraction: {:.2f} %\n",
                   comm->me, data_t, tp_t, (data_t / (data_t + tp_t) * 1e2));
}

/* ---------------------------------------------------------------------- */
void PairGRACE::allocate() {
    allocated = 1;
    int n = atom->ntypes + 1;

    memory->create(setflag, n, n, "pair:setflag");
    memory->create(cutsq, n, n, "pair:cutsq");
    memory->create(scale, n, n, "pair:scale");
    map = new int[n];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
void PairGRACE::settings(int narg, char **arg) {
    if (narg > 3) utils::missing_cmd_args(FLERR, "pair_style grace", error);

    // ACE potentials are parameterized in metal units
    if (strcmp("metal", update->unit_style) != 0)
        error->all(FLERR, "GRACE potentials require 'metal' units");

    auto tf_version = TF_Version();
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] TF version: {}\n", tf_version);

    int iarg = 0;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "chunksize") == 0) {
            chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
        } else if (strcmp(arg[iarg], "padding") == 0) {
            neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;

        } else if (strcmp(arg[iarg], "pad_verbose") == 0) {
            pad_verbose = true;
            iarg += 1;
        } else if (strcmp(arg[iarg], "pair_forces") == 0) {
            pair_forces = true;
            iarg += 1;
            if (comm->me == 0)
                utils::logmesg(lmp, "[GRACE] Pair forces are ON \n");
        } else if (strcmp(arg[iarg], "max_number_of_reduction") == 0) {
            max_number_of_reduction = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
            if (comm->me == 0)
                utils::logmesg(lmp, "[GRACE] Maximum number of recompilation during padding reduction: {}\n",
                               max_number_of_reduction);
        } else if (strcmp(arg[iarg], "reduce_padding") == 0) {
            reducing_neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
            if (comm->me == 0)
                utils::logmesg(lmp, "[GRACE] Reducing padding fraction: {}\n", reducing_neigh_padding_fraction);
        } else
            error->all(FLERR, "[GRACE] Unknown pair_style grace keyword: {}", arg[iarg]);
    }

    do_padding = (neigh_padding_fraction > 0);
    if(do_padding)
        if (comm->me == 0)
        utils::logmesg(lmp, "[GRACE] Neighbour padding is ON, padding fraction: {}, max padding fraction before reduction: {}, max number of reduction(s): {}\n",
                       neigh_padding_fraction, reducing_neigh_padding_fraction, max_number_of_reduction);

    if (!pair_forces && comm->nprocs > 1) {
        pair_forces = true;
        if (comm->me == 0)
            utils::logmesg(lmp,
                           "[GRACE] ENFORCE pair-force mode to ON, because number of processes {} is more than one.\n",
                           comm->nprocs);
    }

    if(pair_forces) {
        // mark as available centroid stress flag
        centroidstressflag = CENTROID_AVAIL;
    } else {
        centroidstressflag = CENTROID_NOTAVAIL;
    }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGRACE::coeff(int narg, char **arg) {

    if (!allocated) allocate();

    map_element2type(narg - 3, arg + 3);
    auto potential_path = std::string(arg[2]);

    //load potential file
    delete aceimpl->model;
    //load potential file
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Loading {}\n", potential_path);
    // load cppflow model
    aceimpl->model = new cppflow::model(potential_path);

    // read elements from metadata.yaml
    YAML_PACE::Node metadata_yaml = YAML_PACE::LoadFile(potential_path + "/metadata.yaml");
    auto elements_yaml = metadata_yaml["chemical_symbols"];
    elements_name = elements_yaml.as<std::vector<std::string>>();
    nelements = (int) elements_name.size();
    for (int mu = 0; mu < nelements; mu++) {
        elements_to_index_map[elements_name.at(mu)] = mu;
    }
    cutoff = metadata_yaml["cutoff"].as<double>();

    if (metadata_yaml["cutoff_matrix"]) {
        cutoff_matrix = metadata_yaml["cutoff_matrix"].as<vector<vector<double>>>();
        // assert square size of matrix
        if (cutoff_matrix.size() != nelements)
            error->all(FLERR,
                       "[GRACE] cutoff_matrix is provided, but it's size ({}) is not equal to number of elements ({})\n",
                       cutoff_matrix.size(), nelements);

        for (const auto &v: cutoff_matrix)
            if (v.size() != nelements)
                error->all(FLERR,
                           "[GRACE] cutoff_matrix is provided, but it's row size ({}) is not equal to number of elements ({})\n",
                           v.size(), nelements);

        if (comm->me == 0)
            utils::logmesg(lmp, "[GRACE] Custom cutoff matrix is loaded\n");
        is_custom_cutoffs = true;
    }
    if (comm->me == 0) {
        utils::logmesg(lmp, "[GRACE] Model loaded\n");
    }


    // read args that map atom types to PACE elements
    // map[i] = which element the Ith atom type is, -1 if not mapped
    // map[0] is not used

    const int ntypes = atom->ntypes;
    element_type_mapping.resize(ntypes + 1);
    // elements to species-type map
    for (int i = 1; i <= ntypes; i++) {
        char *elemname = arg[2 + i];
        if (strcmp(elemname, "NULL") == 0) {
            // species_type=-1 value will not reach ACE Evaluator::compute_atom,
            // but if it will ,then error will be thrown there
            element_type_mapping[i] = -1;
            map[i] = -1;
            if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Skipping LAMMPS atom type #{}(NULL)\n", i);
        } else {
            int atomic_number = PACE::AtomicNumberByName(elemname);
            if (atomic_number == -1) error->all(FLERR, "[GRACE] '{}' is not a valid element\n", elemname);
            int mu = elements_to_index_map.at(elemname);
            if (mu != -1) {
                if (comm->me == 0)
                    utils::logmesg(lmp, "[GRACE] Mapping LAMMPS atom type #{}({}) -> ACE species type #{}\n", i,
                                   elemname, mu);
                map[i] = mu;
                // set up LAMMPS atom type to ACE species  mapping for ace evaluator
                element_type_mapping[i] = mu;
            } else {
                error->all(FLERR, "[GRACE] Element {} is not supported by ACE-potential from file {}", elemname,
                           potential_path);
            }
        }
    }

    // initialize scale factor
    for (int i = 1; i <= ntypes; i++) {
        for (int j = i; j <= ntypes; j++) scale[i][j] = 1.0;
    }

    if (is_custom_cutoffs) {
        // matrix of size [ntypes+1][ntypes+1]
        cutoff_matrix_per_lammps_type.resize(ntypes + 1, vector<double>(ntypes + 1));
        double min_cutoff = 1e99, max_cutoff = 0;
        for (int i = 1; i <= ntypes; i++) {
            for (int j = 1; j <= ntypes; j++) {
                auto val = cutoff_matrix[element_type_mapping[i]][element_type_mapping[j]];
                cutoff_matrix_per_lammps_type[i][j] = val;
                if (val < min_cutoff) min_cutoff = val;
                if (val > max_cutoff) max_cutoff = val;
            }
        }

        if (comm->me == 0)
            utils::logmesg(lmp, "[GRACE] Custom cutoffs: min={}, max={}\n", min_cutoff, max_cutoff);

    }

//    auto operations_vec = aceimpl->model->get_operations();
//    std::cout<<"List of operations "<<std::endl;
//    for(const auto& op_name: operations_vec){
//        std::cout<<"Operation `"<<op_name<<"`"<<endl;
//    }

    //
    has_map_atoms_to_structure_op = check_tf_graph_input_presented(aceimpl->model,
                                                                   "serving_default_map_atoms_to_structure");

    has_nstruct_total_op = check_tf_graph_input_presented(aceimpl->model,
                                                          "serving_default_n_struct_total");
    has_mu_i_op = check_tf_graph_input_presented(aceimpl->model,
                                                 "serving_default_mu_i");

}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGRACE::init_style() {
    if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace requires atom IDs");
    if (force->newton_pair == 0) error->all(FLERR, "Pair style grace requires newton pair on");

    // request a full neighbor list
    neighbor->add_request(this, NeighConst::REQ_FULL);

    // request atom map (maybe?)
    if (atom->map_style == Atom::MAP_NONE) {
        atom->map_init();
        atom->map_set();
    }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGRACE::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
    //cutoff from the basis set's radial functions settings
    scale[j][i] = scale[i][j];
    if (is_custom_cutoffs) {
        return cutoff_matrix_per_lammps_type[i][j];
    } else
        return cutoff;
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairGRACE::extract(const char *str, int &dim) {
    dim = 2;
    if (strcmp(str, "scale") == 0) return (void *) scale;
    return nullptr;
}


/* ---------------------------------------------------------------------- */
/**
 * signature_def['serving_default']:
  The given SavedModel SignatureDef contains the following input(s):
    +inputs['atomic_mu_i'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_atomic_mu_i:0
    +inputs['batch_tot_nat'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_batch_tot_nat:0
    +inputs['batch_tot_nat_real'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_batch_tot_nat_real:0
    +inputs['bond_vector'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: serving_default_bond_vector:0
    +inputs['ind_i'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_ind_i:0
    +inputs['ind_j'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_ind_j:0
    +inputs['map_atoms_to_structure'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_map_atoms_to_structure:0
    +inputs['mu_j'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_mu_j:0
    +inputs['n_struct_total'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_n_struct_total:0


  The given SavedModel SignatureDef contains the following output(s):
    outputs['atomic_energy'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 1)
        name: StatefulPartitionedCall:0
    outputs['total_energy'] tensor_info:
        dtype: DT_DOUBLE
        shape: (1, 1)
        name: StatefulPartitionedCall:1
    outputs['total_f'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: StatefulPartitionedCall:2
    outputs['virial'] tensor_info:
        dtype: DT_DOUBLE
        shape: (6)
        name: StatefulPartitionedCall:3
    outputs['z_pair_f'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: StatefulPartitionedCall:4

  Method name is: tensorflow/serving/predict
 */
void PairGRACE::compute(int eflag, int vflag) {
    int i, j, ii, jj, inum, jnum;
    double delx, dely, delz, evdwl;
    double fij[3];
    int *ilist, *jlist, *numneigh, **firstneigh;

    ev_init(eflag, vflag);

    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
//    tagint *tag = atom->tag;

    // number of atoms in cell
    int nlocal = atom->nlocal;
    int n_fake_atoms; // no fake atoms needed. fake bonds will be 1e6
    int n_real_neighbours;
    int n_fake_neighbours;

    int newton_pair = force->newton_pair;

    // inum: length of the neighborlists list
    inum = list->inum;

    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;

    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;

    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;

    data_timer.start();
    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;

    if (do_padding) {
        if (nlocal > tot_atoms) {
            n_fake_atoms = static_cast<int>(std::round(nlocal * neigh_padding_fraction));
            n_fake_atoms = std::max(n_fake_atoms, 1);
            tot_atoms = nlocal + n_fake_atoms; // add fake neighbours
            if (pad_verbose)
                utils::logmesg(lmp,
                               "[GRACE] Atoms padding: new num. of atoms = {} (incl. {:.3f}% fake atoms)\n",
                               tot_atoms, 100. * (double) n_fake_atoms / tot_atoms);
        } else {
            // TODO: select tot_atoms based on previous padding history
        }
    } else {
        tot_atoms = nlocal;
    }
    // atomic_mu_i: per-atom species type + padding with type[0]
    std::vector<int32_t> atomic_mu_i_vector(tot_atoms, element_type_mapping[type[0]]);
    for (i = 0; i < nlocal; ++i)
        atomic_mu_i_vector[i] = element_type_mapping[type[i]];

    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "atomic_mu_i" + ":0",
                        cppflow::tensor(atomic_mu_i_vector, {tot_atoms}));

    //    map_atoms_to_structure
    if (has_map_atoms_to_structure_op) {
        inputs.emplace_back(DEFAULT_INPUT_PREFIX + "map_atoms_to_structure" + ":0",
                            cppflow::tensor(std::vector<int32_t>(tot_atoms, 0), {tot_atoms}));
    }

    // batch_nat = number of extened atoms + padding
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_tot_nat" + ":0",
                        cppflow::tensor(std::vector<int32_t>{tot_atoms}, {}));

    // batch_nreal_atoms_per_structure: number of extened atoms (w/o padding)
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_tot_nat_real" + ":0",
                        cppflow::tensor(std::vector<int32_t>{nlocal}, {}));

    // ind_i, ind_j: bonds
    //determine the maximum number of neighbours (within cutoff)
    std::vector<int> actual_jnum(inum, 0);
    double cutoff_sq = cutoff * cutoff;
    int type_i, type_j;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        type_i = type[i];
        double xtmp = x[i][0];
        double ytmp = x[i][1];
        double ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        int cur_actual_jnum = 0;
        for (jj = 0; jj < jnum; ++jj) {
            j = jlist[jj];
            j &= NEIGHMASK;
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            if (is_custom_cutoffs) {
                type_j = type[j];
                double cur_cutoff = cutoff_matrix_per_lammps_type[type_i][type_j];
                cutoff_sq = cur_cutoff * cur_cutoff;
            }
            const double rsq = delx * delx + dely * dely + delz * delz;
            if (rsq < cutoff_sq) {
                cur_actual_jnum += 1;
            }
        }
        actual_jnum[ii] = cur_actual_jnum;
    }
    n_real_neighbours = std::accumulate(actual_jnum.begin(), actual_jnum.end(), 0);

    std::vector<int> actual_jnum_shift(actual_jnum.size(), 0);
    std::partial_sum(actual_jnum.begin(), actual_jnum.end() - 1, actual_jnum_shift.begin() + 1);

    if (do_padding) {
        // try to find value that strictly greater than n_real_neigh
        auto upper_bound_tot_neighbours = tot_neighbours_set.upper_bound(n_real_neighbours);
        if (upper_bound_tot_neighbours == tot_neighbours_set.end()) { // not found
            // if no previous larger element - create new
            n_fake_neighbours = static_cast<int>(std::round(n_real_neighbours * neigh_padding_fraction));
            n_fake_neighbours = std::max(n_fake_neighbours, 1);
            tot_neighbours = n_real_neighbours + n_fake_neighbours; // add fake neighbours
            tot_neighbours_set.insert(tot_neighbours);
            if (pad_verbose)
                utils::logmesg(lmp,
                               "[GRACE] Neighbours padding: extending new num. of neighbours = {} (incl. {:.3f}% fake neighbours)\n",
                               tot_neighbours, 100. * (double) n_fake_neighbours / tot_neighbours);
        } else {
            // upper bound found
            tot_neighbours = *upper_bound_tot_neighbours;
            n_fake_neighbours = tot_neighbours - n_real_neighbours;

            // if upper bound found is FIRST bound then check for too much fake neighbours
            if (upper_bound_tot_neighbours == tot_neighbours_set.begin()) {
                // no smaller tot_neighbours, check if too many n_fake neighbours
                // and limit of recompilation is not reached
                if (n_fake_neighbours > std::round(n_real_neighbours * reducing_neigh_padding_fraction) &&
                    (num_of_reductions < max_number_of_reduction || max_number_of_reduction == -1)) {
                    // if too many fake neighbours  - reduce
                    tot_neighbours = n_real_neighbours;
                    tot_neighbours_set.insert(tot_neighbours);
                    num_of_reductions++;
                    if (pad_verbose)
                        utils::logmesg(lmp,
                                       "[GRACE] Neighbours padding: reducing new num. of neighbours = {}\n",
                                       tot_neighbours);
                }
            }
        }

    } else {
        // no padding
        tot_neighbours = n_real_neighbours;
    }

//    utils::logmesg(lmp,"[GRACE-DEBUG, #{}] tot_atoms={}, tot_neighbours={} \n", comm->me, tot_atoms,  tot_neighbours);

    std::vector<int32_t> ind_i_vector(tot_neighbours);
    std::vector<int32_t> ind_j_vector(tot_neighbours);
    std::vector<int32_t> mu_i_vector(tot_neighbours);
    std::vector<int32_t> mu_j_vector(tot_neighbours);
    std::vector<double> bond_vector(3 * tot_neighbours, 1e6);

    for (ii = 0; ii < inum; ++ii) {
        i = ilist[ii];
        type_i = type[i];
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        int tot_ind = actual_jnum_shift[ii];
        for (jj = 0; jj < jnum; ++jj) {
            j = jlist[jj];
            j &= NEIGHMASK;
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            if (is_custom_cutoffs) {
                type_j = type[j];
                double cur_cutoff = cutoff_matrix_per_lammps_type[type_i][type_j];
                cutoff_sq = cur_cutoff * cur_cutoff;
            }
            const double rsq = delx * delx + dely * dely + delz * delz;
            if (rsq < cutoff_sq) {
                ind_i_vector[tot_ind] = i;
                // remap j to j_local
                int j_local = atom->map(atom->tag[j]);
                ind_j_vector[tot_ind] = j_local;
                mu_i_vector[tot_ind] = element_type_mapping[type[i]];
                mu_j_vector[tot_ind] = element_type_mapping[type[j]];
                double bondx = atom->x[j][0] - atom->x[i][0];
                double bondy = atom->x[j][1] - atom->x[i][1];
                double bondz = atom->x[j][2] - atom->x[i][2];


                bond_vector[3 * tot_ind + 0] = bondx;
                bond_vector[3 * tot_ind + 1] = bondy;
                bond_vector[3 * tot_ind + 2] = bondz;
                ++tot_ind;
            }
        }
    }

    // add fake bonds
    int fake_atom_ind = tot_atoms - 1;
    for (int tot_ind = n_real_neighbours; tot_ind < tot_neighbours; tot_ind++) {
        ind_i_vector[tot_ind] = fake_atom_ind; // fake atom ind
        ind_j_vector[tot_ind] = fake_atom_ind; // fake atom ind
        mu_i_vector[tot_ind] = 0;
        mu_j_vector[tot_ind] = 0;
    }

    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "ind_i" + ":0",
                        cppflow::tensor(ind_i_vector, {tot_neighbours}));
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "ind_j" + ":0",
                        cppflow::tensor(ind_j_vector, {tot_neighbours}));


    // mu_i, mu_j: bonds
    if (has_mu_i_op) {
        inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_i" + ":0",
                            cppflow::tensor(mu_i_vector, {tot_neighbours}));
    }
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_j" + ":0",
                        cppflow::tensor(mu_j_vector, {tot_neighbours}));


    // num_struc: 1
    if (has_nstruct_total_op) {
        inputs.emplace_back(DEFAULT_INPUT_PREFIX + "n_struct_total" + ":0",
                            cppflow::tensor(std::vector<int32_t>{1}, {}));
    }

    // vector_offsets: 1
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "bond_vector" + ":0",
                        cppflow::tensor(bond_vector, {tot_neighbours, 3}));


    data_timer.stop();
    tp_timer.start();
    vector<string> output_names = {
            "StatefulPartitionedCall:0", // atomic_energy [nat,1]
            "StatefulPartitionedCall:1", // total_energy [-1, 1]
            "StatefulPartitionedCall:2", // total_f [n_at, 3]
            "StatefulPartitionedCall:3", // virial [6]
    };
    // add it optionally
    if (pair_forces)
        output_names.emplace_back("StatefulPartitionedCall:4");// pair_f [n_bonds, 3]

    //CALL MODEL
    std::vector<cppflow::tensor> output = aceimpl->model->operator()(
            inputs,
            output_names
    );
    tp_timer.stop();
//    std::cout << "Ave.timing: " << (double) tp_timer.as_microseconds() / nlocal << " mcs/at" << std::endl;

    data_timer.start();
    auto &e_out = output[0]; // atomic_energy
    auto e_tens = e_out.get_tensor();
    const double *e_data = static_cast<double *>(TF_TensorData(e_tens.get()));

//    auto &te_out = output[1]; // total_energy

    auto &total_f_out = output[2]; // total_f
    auto total_f_tens = total_f_out.get_tensor();
    const double *total_f_data = static_cast<double *>(TF_TensorData(total_f_tens.get()));


    if (!pair_forces) {
        for (ii = 0; ii < inum; ii++) {
            i = ilist[ii];

            const int itype = type[i];
            double fx = total_f_data[ii * 3 + 0];
            double fy = total_f_data[ii * 3 + 1];
            double fz = total_f_data[ii * 3 + 2];


            f[i][0] += scale[itype][itype] * fx;
            f[i][1] += scale[itype][itype] * fy;
            f[i][2] += scale[itype][itype] * fz;


            // tally energy contribution
            if (eflag_either) {
                // evdwl = energy of atom I
                evdwl = scale[itype][itype] * e_data[i];
                ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
            }
        } // end for(ii)

        // virial order, seems to be OK
        if (vflag_global) {
            auto &v_out = output[3]; // virial
            auto v_tens = v_out.get_tensor();
            const double *v_data = static_cast<double *>(TF_TensorData(v_tens.get()));

//            ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], -delx, -dely, -delz);
            virial[0] += v_data[0];
            virial[1] += v_data[1];
            virial[2] += v_data[2];
            virial[3] += v_data[3];
            virial[4] += v_data[4];
            virial[5] += v_data[5];
        }
    } else {
        // pair forces
        auto &f_out = output[4];

        auto f_tens = f_out.get_tensor();
        const double *f_data = static_cast<double *>(TF_TensorData(f_tens.get()));
        int tot_ind = 0;
        for (ii = 0; ii < inum; ++ii) {
            i = ilist[ii];
            type_i = type[i];
            const double xtmp = x[i][0];
            const double ytmp = x[i][1];
            const double ztmp = x[i][2];
            jlist = firstneigh[i];
            jnum = numneigh[i];

            for (jj = 0; jj < jnum; ++jj) {
                j = jlist[jj];
                j &= NEIGHMASK;
                delx = xtmp - x[j][0];
                dely = ytmp - x[j][1];
                delz = ztmp - x[j][2];
                if (is_custom_cutoffs) {
                    type_j = type[j];
                    double cur_cutoff = cutoff_matrix_per_lammps_type[type_i][type_j];
                    cutoff_sq = cur_cutoff * cur_cutoff;
                }
                const double rsq = delx * delx + dely * dely + delz * delz;
                if (rsq < cutoff_sq) {
                    // why "-" ?
                    fij[0] = -scale[type_i][type_i] * f_data[tot_ind];
                    fij[1] = -scale[type_i][type_i] * f_data[tot_ind + 1];
                    fij[2] = -scale[type_i][type_i] * f_data[tot_ind + 2];
                    tot_ind += 3;

                    f[i][0] += fij[0];
                    f[i][1] += fij[1];
                    f[i][2] += fij[2];
                    f[j][0] -= fij[0];
                    f[j][1] -= fij[1];
                    f[j][2] -= fij[2];

                    // tally per-atom virial contribution, OpenMP critical !
                    if (vflag_either) {
                        ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], delx, dely, delz);


                        // Centroid Stress
                        if (cvflag_atom) {
                            const double fx = fij[0];
                            const double fy = fij[1];
                            const double fz = fij[2];

                            cvatom[i][0] += 0.5 * delx * fx; // xx
                            cvatom[i][1] += 0.5 * dely * fy; // yy
                            cvatom[i][2] += 0.5 * delz * fz; // zz
                            cvatom[i][3] += 0.5 * delx * fy; // xy
                            cvatom[i][4] += 0.5 * delx * fz; // xz
                            cvatom[i][5] += 0.5 * dely * fz; // yz
                            cvatom[i][6] += 0.5 * dely * fx; // yx
                            cvatom[i][7] += 0.5 * delz * fx; // zx
                            cvatom[i][8] += 0.5 * delz * fy; // zy


                            cvatom[j][0] += 0.5 * delx * fx; // xx
                            cvatom[j][1] += 0.5 * dely * fy; // yy
                            cvatom[j][2] += 0.5 * delz * fz; // zz
                            cvatom[j][3] += 0.5 * delx * fy; // xy
                            cvatom[j][4] += 0.5 * delx * fz; // xz
                            cvatom[j][5] += 0.5 * dely * fz; // yz
                            cvatom[j][6] += 0.5 * dely * fx; // yx
                            cvatom[j][7] += 0.5 * delz * fx; // zx
                            cvatom[j][8] += 0.5 * delz * fy; // zy
                        }
                    }
                }
            } // loop over neighbours

            // tally energy contribution
            if (eflag_either) {
                // evdwl = energy of atom I
                evdwl = scale[type_i][type_i] * e_data[i];
                ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
            }

        } // loop over atoms -i

        if (vflag_fdotr) virial_fdotr_compute();
    }


    data_timer.stop();
    // end modifications YL
}

#endif //#ifndef NO_GRACE_TF
