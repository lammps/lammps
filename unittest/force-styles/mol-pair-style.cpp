// unit tests for pair styles intended for molecular systems

#include "lammps.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "info.h"
#include "input.h"
#include "universe.h"
#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <iostream>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "yaml.h"

using ::testing::StartsWith;
using ::testing::HasSubstr;

typedef struct {
    double x,y,z;
} coord_t;

typedef struct {
    double xx,yy,zz,xy,xz,yz;
} stress_t;

enum state_value {
    START,
    ACCEPT_KEY,
    ACCEPT_VALUE,
    STOP,
    ERROR,
};

struct parser_state {
    enum state_value state;
    int accepted;
    std::string key;
};

class TestConfig {
public:
    std::string lammps_version;
    std::string date_generated;
    std::vector<std::string> pre_commands;
    std::vector<std::string> post_commands;
    std::string input_file;
    std::string pair_style;
    std::vector<std::string> pair_coeff;
    int natoms;
    double init_vdwl;
    double run_vdwl;
    double init_coul;
    double run_coul;
    stress_t init_stress;
    stress_t run_stress;
    std::vector<coord_t> init_forces;
    std::vector<coord_t> run_forces;
    std::vector<coord_t> run_energy;
    TestConfig() : lammps_version(""),
                   date_generated(""),
                   input_file(""),
                   pair_style("zero"),
                   natoms(0),
                   init_vdwl(0),
                   run_vdwl(0),
                   init_coul(0),
                   run_coul(0),
                   init_stress({0,0,0,0,0,0}),
                   run_stress({0,0,0,0,0,0}) {
        pre_commands.clear();
        post_commands.clear();
        pair_coeff.clear();
        init_forces.clear();
        run_forces.clear();
        run_energy.clear();
    }
} test_config;

// default floating point error margin
double float_epsilon = 1.0e-14;

#define EXPECT_FP_EQ_WITH_EPS(val1,val2,eps)                \
    do {                                                    \
        const double diff = fabs(val1-val2);                \
        const double div = std::max(fabs(val1),fabs(val2)); \
        const double err = (div == 0.0) ? diff : diff/div;  \
        EXPECT_PRED_FORMAT2(::testing::DoubleLE, err, eps); \
    } while (0);

LAMMPS_NS::LAMMPS *init_lammps(int argc, char **argv, const TestConfig &cfg)
{
    LAMMPS_NS::LAMMPS *lmp;

    lmp = new LAMMPS_NS::LAMMPS(argc, argv, MPI_COMM_WORLD);
    for (std::size_t i=0; i < cfg.pre_commands.size(); ++i)
        lmp->input->one(cfg.pre_commands[i].c_str());
    lmp->input->file(cfg.input_file.c_str());

    // determine if pair style is available while applying suffix if active
    LAMMPS_NS::Info *info = new LAMMPS_NS::Info(lmp);
    std::size_t found = cfg.pair_style.find(" ");
    std::string style;
    if ((found > 0) && (found != std::string::npos)) {
        style = cfg.pair_style.substr(0,found);
    } else {
        style = cfg.pair_style;
    }
    if (lmp->suffix_enable) {
        style += "/";
        style += lmp->suffix;
    }
    if (!info->has_style("pair", style)) {
        test_config.pair_style = style;  // for error message
        delete info;
        delete lmp;
        return NULL;
    }

    std::string cmd("pair_style ");
    cmd += cfg.pair_style;
    lmp->input->one(cmd.c_str());
    for (std::size_t i=0; i < cfg.pair_coeff.size(); ++i) {
        cmd = "pair_coeff " + cfg.pair_coeff[i];
        lmp->input->one(cmd.c_str());
    }
    for (std::size_t i=0; i < cfg.post_commands.size(); ++i)
        lmp->input->one(cfg.post_commands[i].c_str());
    lmp->input->one("run 0 post no");
    return lmp;
}

void run_lammps(LAMMPS_NS::LAMMPS *lmp)
{
    lmp->input->one("fix 1 all nve");
    lmp->input->one("compute pe all pe/atom");
    lmp->input->one("compute sum all reduce sum c_pe");
    lmp->input->one("thermo_style custom step temp pe press c_sum");
    lmp->input->one("thermo 2");
    lmp->input->one("run 4 post no");
}

int consume_event(struct parser_state *s, yaml_event_t *event)
{
    s->accepted = 0;
    switch (s->state) {
      case START:
          switch (event->type) {
            case YAML_MAPPING_START_EVENT:
                s->state = ACCEPT_KEY;
                break;
            case YAML_SCALAR_EVENT:
            case YAML_SEQUENCE_START_EVENT:
                s->state = ERROR;
                break;
            case YAML_STREAM_END_EVENT:
                s->state = STOP;
                break;
            case YAML_STREAM_START_EVENT:
            case YAML_DOCUMENT_START_EVENT:
            case YAML_DOCUMENT_END_EVENT:
                // ignore
                break;
            default:
                std::cerr << "UNHANDLED YAML EVENT:  " << event->type << std::endl;
                s->state = ERROR;
                break;
          }
          break;

      case ACCEPT_KEY:
          switch (event->type) {
            case YAML_SCALAR_EVENT:
                s->key =(char *) event->data.scalar.value;
                s->state = ACCEPT_VALUE;
                break;
            case YAML_MAPPING_END_EVENT:
                s->state = STOP;
                break;
            default:
                std::cerr << "UNHANDLED YAML EVENT (key): " << event->type
                          << "\nVALUE: " << event->data.scalar.value
                          << std::endl;
                s->state = ERROR;
                break;
          }
          break;

      case ACCEPT_VALUE:
          switch (event->type) {
            case YAML_SCALAR_EVENT:
                s->accepted = 1;
                s->state = ACCEPT_KEY;
                break;
            default:
                std::cerr << "UNHANDLED YAML EVENT (value): " << event->type
                          << "\nVALUE: " << event->data.scalar.value
                          << std::endl;
                s->state = ERROR;
                break;
          }
          break;

      case ERROR:
      case STOP:
          break;
    }
    return (s->state == ERROR) ? 0 : 1;
}

int parse(const char *infile)
{
    FILE *fp;
    yaml_parser_t parser;
    yaml_event_t event;
    struct parser_state state;

    fp = fopen(infile,"r");
    if (!fp) {
        std::cerr << "Cannot open yaml file '" << infile
                  << "': " << strerror(errno) << std::endl;
        return 1;
    }
    yaml_parser_initialize(&parser);
    yaml_parser_set_input_file(&parser, fp);

    int retval = 0;
    state.state = START;
    do {
        if (!yaml_parser_parse(&parser, &event)) {
            retval = 1;
            state.state = STOP;
        }

        if (!consume_event(&state, &event)) {
            retval = 1;
            state.state = STOP;
        }

        if (state.accepted) {
            if (state.key == "pre_commands") {
                test_config.pre_commands.clear();
                std::string data  = (char *)event.data.scalar.value;
                std::size_t first = 0;
                std::size_t found = data.find("\n");
                while (found != std::string::npos) {
                    std::string line(data.substr(first,found));
                    test_config.pre_commands.push_back(line);
                    data = data.substr(found+1);
                    found = data.find("\n");
                }
            } else if (state.key == "post_commands") {
                test_config.post_commands.clear();
                std::string data  = (char *)event.data.scalar.value;
                std::size_t first = 0;
                std::size_t found = data.find("\n");
                while (found != std::string::npos) {
                    std::string line(data.substr(first,found));
                    test_config.post_commands.push_back(line);
                    data = data.substr(found+1);
                    found = data.find("\n");
                }
            } else if (state.key == "lammps_version") {
                test_config.lammps_version = (char *)event.data.scalar.value;
            } else if (state.key == "date_generated") {
                test_config.date_generated = (char *)event.data.scalar.value;
            } else if (state.key == "input_file") {
                test_config.input_file = (char *)event.data.scalar.value;
            } else if (state.key == "pair_style") {
                test_config.pair_style = (char *)event.data.scalar.value;
            } else if (state.key == "pair_coeff") {
                test_config.pair_coeff.clear();
                std::string data  = (char *)event.data.scalar.value;
                std::size_t first = 0;
                std::size_t found = data.find("\n");
                while (found != std::string::npos) {
                    std::string line(data.substr(first,found));
                    test_config.pair_coeff.push_back(line);
                    data = data.substr(found+1);
                    found = data.find("\n");
                }
            } else if (state.key == "natoms") {
                test_config.natoms = atoi((char *)event.data.scalar.value);
            } else if (state.key == "init_vdwl") {
                test_config.init_vdwl = atof((char *)event.data.scalar.value);
            } else if (state.key == "run_vdwl") {
                test_config.run_vdwl = atof((char *)event.data.scalar.value);
            } else if (state.key == "init_coul") {
                test_config.init_coul = atof((char *)event.data.scalar.value);
            } else if (state.key == "run_coul") {
                test_config.run_coul = atof((char *)event.data.scalar.value);
            } else if (state.key == "init_stress") {
                stress_t stress;
                sscanf((char *)event.data.scalar.value,
                       "%lg %lg %lg %lg %lg %lg",
                       &stress.xx, &stress.yy, &stress.zz,
                       &stress.xy, &stress.xz, &stress.yz);
                test_config.init_stress = stress;
            } else if (state.key == "run_stress") {
                stress_t stress;
                sscanf((char *)event.data.scalar.value,
                       "%lg %lg %lg %lg %lg %lg",
                       &stress.xx, &stress.yy, &stress.zz,
                       &stress.xy, &stress.xz, &stress.yz);
                test_config.run_stress = stress;
            } else if (state.key == "init_forces") {
                test_config.init_forces.clear();
                test_config.init_forces.resize(test_config.natoms+1);
                std::string data  = (char *)event.data.scalar.value;
                std::size_t first = 0;
                std::size_t found = data.find("\n");
                while (found != std::string::npos) {
                    std::string line(data.substr(first,found));
                    coord_t xyz;
                    int tag;
                    sscanf(line.c_str(), "%d %lg %lg %lg", &tag, &xyz.x, &xyz.y, &xyz.z);
                    test_config.init_forces[tag] = xyz;
                    data = data.substr(found+1);
                    found = data.find("\n");
                }
            } else if (state.key == "run_forces") {
                test_config.run_forces.clear();
                test_config.run_forces.resize(test_config.natoms+1);
                std::string data  = (char *)event.data.scalar.value;
                std::size_t first = 0;
                std::size_t found = data.find("\n");
                while (found != std::string::npos) {
                    std::string line(data.substr(first,found));
                    coord_t xyz;
                    int tag;
                    sscanf(line.c_str(), "%d %lg %lg %lg", &tag, &xyz.x, &xyz.y, &xyz.z);
                    test_config.run_forces[tag] = xyz;
                    data = data.substr(found+1);
                    found = data.find("\n");
                }
            } else std::cerr << "Ignoring unknown key/value pair: " << state.key
                             << " = " << event.data.scalar.value << std::endl;
        }
        yaml_event_delete(&event);
    } while (state.state != STOP);

    yaml_parser_delete(&parser);
    fclose(fp);
    return retval;
}


class YamlWriter {
    FILE *fp;
    yaml_emitter_t emitter;
    yaml_event_t   event;
public:
    YamlWriter(const char * outfile) {
        yaml_emitter_initialize(&emitter);
        fp = fopen(outfile, "w");
        if (!fp) {
            perror(__FILE__);
            return;
        }

        yaml_emitter_set_output_file(&emitter, fp);

        yaml_stream_start_event_initialize(&event, YAML_UTF8_ENCODING);
        yaml_emitter_emit(&emitter, &event);
        yaml_document_start_event_initialize(&event, NULL, NULL, NULL, 0);
        yaml_emitter_emit(&emitter, &event);
        yaml_mapping_start_event_initialize(&event, NULL,
                                            (yaml_char_t *)YAML_MAP_TAG,
                                            1, YAML_ANY_MAPPING_STYLE);
        yaml_emitter_emit(&emitter, &event);
    }

    ~YamlWriter() {
        yaml_mapping_end_event_initialize(&event);
        yaml_emitter_emit(&emitter, &event);
        yaml_document_end_event_initialize(&event, 0);
        yaml_emitter_emit(&emitter, &event);
        yaml_stream_end_event_initialize(&event);
        yaml_emitter_emit(&emitter, &event);
        yaml_emitter_delete(&emitter);
        fclose(fp);
    }

    void emit(const std::string &key, const double value){
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *) key.c_str(),
                                    key.size(), 1, 0,
                                    YAML_PLAIN_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
        char buf[256];
        snprintf(buf,256,"%.15g",value);
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *)buf,
                                    strlen(buf), 1, 0,
                                    YAML_PLAIN_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
    }

    void emit(const std::string &key, const long value){
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *) key.c_str(),
                                    key.size(), 1, 0,
                                    YAML_PLAIN_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
        char buf[256];
        snprintf(buf,256,"%ld",value);
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *)buf,
                                    strlen(buf), 1, 0,
                                    YAML_PLAIN_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
    }

    void emit(const std::string &key, const int value){
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *) key.c_str(),
                                    key.size(), 1, 0,
                                    YAML_PLAIN_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
        char buf[256];
        snprintf(buf,256,"%d",value);
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *)buf,
                                    strlen(buf), 1, 0,
                                    YAML_PLAIN_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
    }

    void emit(const std::string &key, const std::string &value)
    {
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *) key.c_str(),
                                    key.size(), 1, 0,
                                    YAML_PLAIN_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *)value.c_str(),
                                    value.size(), 1, 0,
                                    YAML_PLAIN_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
    }

    void emit_block(const std::string &key, const std::string &value){
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *) key.c_str(),
                                    key.size(), 1, 0,
                                    YAML_PLAIN_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
        yaml_scalar_event_initialize(&event, NULL,
                                    (yaml_char_t *)YAML_STR_TAG,
                                    (yaml_char_t *)value.c_str(),
                                    value.size(), 1, 0,
                                    YAML_LITERAL_SCALAR_STYLE);
        yaml_emitter_emit(&emitter, &event);
    }
};

void generate(const char *outfile) {

    // initialize molecular system geometry
    const char *args[] = {"MolPairStyle", "-log", "none", "-echo", "screen", "-nocite" };
    char **argv = (char **)args;
    int argc = sizeof(args)/sizeof(char *);
    LAMMPS_NS::LAMMPS *lmp = init_lammps(argc,argv,test_config);
    if (!lmp) {
        std::cerr << "Pair style: " << test_config.pair_style << " is not available."
            "in this LAMMPS configuration\n";
        fclose(fp);
        return;
    }

    const int natoms = lmp->atom->natoms;
    const int bufsize = 256;
    char buf[bufsize];
    std::string block("");

    YamlWriter writer(outfile);

    // lammps_version
    writer.emit("lammps_version", lmp->universe->version);

    // date_generated
    std::time_t now = time(NULL);
    block = ctime(&now);
    block = block.substr(0,block.find("\n")-1);
    writer.emit("date_generated", block);

    // pre_commands
    block.clear();
    for (auto command :  test_config.pre_commands) {
        block += command + "\n";
    }
    writer.emit_block("pre_commands", block);

    // post_commands
    block.clear();
    for (auto command : test_config.post_commands) {
        block += command + "\n";
    }
    writer.emit_block("post_commands", block);

    // input_file
    writer.emit("input_file", test_config.input_file);

    // pair_style
    writer.emit("pair_style", test_config.pair_style);

    // pair_coeff
    block.clear();
    for (std::size_t i=0; i < test_config.pair_coeff.size(); ++i) {
        block += test_config.pair_coeff[i] + "\n";
    }
    writer.emit_block("pair_coeff", block);

    // natoms
    writer.emit("natoms", natoms);

    // init_vdwl
    writer.emit("init_vdwl", lmp->force->pair->eng_vdwl);

    // init_coul
    writer.emit("init_coul", lmp->force->pair->eng_coul);

    // init_stress
    double *stress = lmp->force->pair->virial;
    snprintf(buf,bufsize,"%.15g %.15g %.15g %.15g %.15g %.15g",
             stress[0],stress[1],stress[2],stress[3],stress[4],stress[5]);
    writer.emit("init_stress", buf);

    // init_forces
    block.clear();
    double **f = lmp->atom->f;
    LAMMPS_NS::tagint *tag = lmp->atom->tag;
    for (int i=0; i < natoms; ++i) {
        snprintf(buf,bufsize,"%d %.15g %.15g %.15g\n",
                 (int)tag[i], f[i][0], f[i][1], f[i][2]);
        block += buf;
    }
    writer.emit_block("init_forces", block);

    // do a few steps of MD
    run_lammps(lmp);

    // run_vdwl
    writer.emit("run_vdwl", lmp->force->pair->eng_vdwl);

    // run_could
    writer.emit("run_coul", lmp->force->pair->eng_coul);

    // run_stress
    stress = lmp->force->pair->virial;
    snprintf(buf,bufsize,"%.15g %.15g %.15g %.15g %.15g %.15g",
             stress[0],stress[1],stress[2],stress[3],stress[4],stress[5]);
    writer.emit("run_stress", buf);

    block.clear();
    f = lmp->atom->f;
    tag = lmp->atom->tag;
    for (int i=0; i < natoms; ++i) {
        snprintf(buf,bufsize,"%d %.15g %.15g %.15g\n",
                 (int)tag[i], f[i][0], f[i][1], f[i][2]);
        block += buf;
    }
    writer.emit("run_forces", block);

    delete lmp;
    return;
}

TEST(MolPairStyle, plain) {

    const char *args[] = {"MolPairStyle", "-log", "none", "-echo", "screen", "-nocite" };
    char **argv = (char **)args;
    int argc = sizeof(args)/sizeof(char *);
    ::testing::internal::CaptureStdout();
    LAMMPS_NS::LAMMPS *lmp = init_lammps(argc,argv,test_config);
    std::string output = ::testing::internal::GetCapturedStdout();
    if (!lmp) GTEST_SKIP();
    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms,nlocal);

    double **f=lmp->atom->f;
    LAMMPS_NS::tagint *tag=lmp->atom->tag;
    const std::vector<coord_t> &f_ref = test_config.init_forces;
    for (int i=0; i < nlocal; ++i) {
        EXPECT_FP_EQ_WITH_EPS(f[i][0], f_ref[tag[i]].x, float_epsilon);
        EXPECT_FP_EQ_WITH_EPS(f[i][1], f_ref[tag[i]].y, float_epsilon);
        EXPECT_FP_EQ_WITH_EPS(f[i][2], f_ref[tag[i]].z, float_epsilon);
    }

    LAMMPS_NS::Pair *pair = lmp->force->pair;
    double *stress = pair->virial;
    EXPECT_FP_EQ_WITH_EPS(stress[0], test_config.init_stress.xx, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[1], test_config.init_stress.yy, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[2], test_config.init_stress.zz, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[3], test_config.init_stress.xy, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[4], test_config.init_stress.xz, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[5], test_config.init_stress.yz, float_epsilon);
    
    EXPECT_FP_EQ_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(pair->eng_coul, test_config.init_coul, float_epsilon);

    ::testing::internal::CaptureStdout();
    run_lammps(lmp);
    ::testing::internal::GetCapturedStdout();

    f = lmp->atom->f;
    stress = pair->virial;
    const std::vector<coord_t> &f_run = test_config.run_forces;
    ASSERT_EQ(nlocal+1,f_run.size());
    for (int i=0; i < nlocal; ++i) {
        EXPECT_FP_EQ_WITH_EPS(f[i][0], f_run[tag[i]].x, float_epsilon*100);
        EXPECT_FP_EQ_WITH_EPS(f[i][1], f_run[tag[i]].y, float_epsilon*100);
        EXPECT_FP_EQ_WITH_EPS(f[i][2], f_run[tag[i]].z, float_epsilon*100);
    }

    stress = pair->virial;
    EXPECT_FP_EQ_WITH_EPS(stress[0], test_config.run_stress.xx, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[1], test_config.run_stress.yy, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[2], test_config.run_stress.zz, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[3], test_config.run_stress.xy, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[4], test_config.run_stress.xz, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[5], test_config.run_stress.yz, float_epsilon*100);
    
    EXPECT_FP_EQ_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(pair->eng_coul, test_config.run_coul, float_epsilon*100);

    ::testing::internal::CaptureStdout();
    delete lmp;
    ::testing::internal::GetCapturedStdout();
};

TEST(MolPairStyle, omp) {
    if (!LAMMPS_NS::LAMMPS::is_installed_pkg("USER-OMP")) GTEST_SKIP();
    const char *args[] = {"MolPairStyle", "-log", "none", "-echo", "screen",
                          "-nocite", "-pk", "omp", "4", "-sf", "omp"};
    char **argv = (char **)args;
    int argc = sizeof(args)/sizeof(char *);
    ::testing::internal::CaptureStdout();
    LAMMPS_NS::LAMMPS *lmp = init_lammps(argc,argv,test_config);
    std::string output = ::testing::internal::GetCapturedStdout();
    if (!lmp) GTEST_SKIP();
    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms,nlocal);

    double **f=lmp->atom->f;
    LAMMPS_NS::tagint *tag=lmp->atom->tag;
    const std::vector<coord_t> &f_ref = test_config.init_forces;
    for (int i=0; i < nlocal; ++i) {
        EXPECT_FP_EQ_WITH_EPS(f[i][0], f_ref[tag[i]].x, float_epsilon);
        EXPECT_FP_EQ_WITH_EPS(f[i][1], f_ref[tag[i]].y, float_epsilon);
        EXPECT_FP_EQ_WITH_EPS(f[i][2], f_ref[tag[i]].z, float_epsilon);
    }

    LAMMPS_NS::Pair *pair = lmp->force->pair;
    double *stress = pair->virial;
    EXPECT_FP_EQ_WITH_EPS(stress[0], test_config.init_stress.xx, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[1], test_config.init_stress.yy, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[2], test_config.init_stress.zz, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[3], test_config.init_stress.xy, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[4], test_config.init_stress.xz, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(stress[5], test_config.init_stress.yz, float_epsilon);

    EXPECT_FP_EQ_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, float_epsilon);
    EXPECT_FP_EQ_WITH_EPS(pair->eng_coul, test_config.init_coul, float_epsilon);

    ::testing::internal::CaptureStdout();
    run_lammps(lmp);
    ::testing::internal::GetCapturedStdout();

    f = lmp->atom->f;
    stress = pair->virial;
    const std::vector<coord_t> &f_run = test_config.run_forces;
    ASSERT_EQ(nlocal+1,f_run.size());
    for (int i=0; i < nlocal; ++i) {
        EXPECT_FP_EQ_WITH_EPS(f[i][0], f_run[tag[i]].x, float_epsilon*100);
        EXPECT_FP_EQ_WITH_EPS(f[i][1], f_run[tag[i]].y, float_epsilon*100);
        EXPECT_FP_EQ_WITH_EPS(f[i][2], f_run[tag[i]].z, float_epsilon*100);
    }

    stress = pair->virial;
    EXPECT_FP_EQ_WITH_EPS(stress[0], test_config.run_stress.xx, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[1], test_config.run_stress.yy, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[2], test_config.run_stress.zz, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[3], test_config.run_stress.xy, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[4], test_config.run_stress.xz, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(stress[5], test_config.run_stress.yz, float_epsilon*100);

    EXPECT_FP_EQ_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, float_epsilon*100);
    EXPECT_FP_EQ_WITH_EPS(pair->eng_coul, test_config.run_coul, float_epsilon*100);

    ::testing::internal::CaptureStdout();
    delete lmp;
    ::testing::internal::GetCapturedStdout();
};

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);

    if ((argc != 2) && (argc != 4)) {
        std::cerr << "usage: " << argv[0] << " <testfile.yaml> "
            "[--gen <newfile.yaml> | --eps <floating-point epsilon>]" << std::endl;
        return 1;
    }
    if (parse(argv[1])) {
        std::cerr << "Error parsing yaml file: " << argv[1] << std::endl;
        return 2;
    }

    if (argc == 4) {
        if (strcmp(argv[2],"--gen") == 0) {
            generate(argv[3]);
            return 0;
        } else if (strcmp(argv[2],"--eps") == 0) {
            float_epsilon = atof(argv[3]);
        } else {
            std::cerr << "usage: " << argv[0] << " <testfile.yaml> "
                "[--gen <newfile.yaml> | --eps <floating-point epsilon>]" << std::endl;
            return 1;
        }
    }
    return RUN_ALL_TESTS();
    return 0;
}
