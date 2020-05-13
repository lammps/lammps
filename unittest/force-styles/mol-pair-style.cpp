// unit tests for pair styles intended for molecular systems

#include "lammps.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "input.h"
#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cerrno>
#include <string>
#include <vector>
#include <iostream>

#include "gtest/gtest.h"
#include "yaml.h"

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
    bool pair_newton_on;
    bool pair_newton_off;
    bool has_global_energy;
    bool has_global_stress;
    bool has_atom_energy;
    bool has_atom_stress;
    std::vector<std::string> pre_commands;
    std::vector<std::string> post_commands;
    std::string input_file;
    std::string pair_style;
    std::vector<std::string> pair_coeff;
    int natoms;
    std::vector<coord_t> init_forces;
    double init_vdwl;
    double init_coul;
    stress_t init_stress;
    TestConfig() : pair_newton_on(false),
                   pair_newton_off(false),
                   has_global_energy(false),
                   has_global_stress(false),
                   has_atom_energy(false),
                   has_atom_stress(false),
                   input_file(""),
                   pair_style("zero"),
                   natoms(0),
                   init_vdwl(0),
                   init_coul(0),
                   init_stress({0,0,0,0,0,0}) {
        pre_commands.clear();
        post_commands.clear();
        pair_coeff.clear();
        init_forces.clear();
    }
} test_config;


LAMMPS_NS::LAMMPS *init_lammps(int argc, char **argv, const TestConfig &cfg)
{
    LAMMPS_NS::LAMMPS *lmp;

    int flag;
    MPI_Initialized(&flag);
    if (!flag) MPI_Init(&argc,&argv);

    lmp = new LAMMPS_NS::LAMMPS(argc, argv, MPI_COMM_WORLD);
    for (std::size_t i=0; i < cfg.pre_commands.size(); ++i)
        lmp->input->one(cfg.pre_commands[i].c_str());
    lmp->input->file(cfg.input_file.c_str());
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
                std::cerr << "UNKNOWN YAML EVENT:  " << event->type << std::endl;
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
                std::cerr << "UNKNOWN YAML EVENT (key): " << event->type
                          << "\nVALUE: " << event->data.scalar.value
                          << std::endl;
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
                std::cerr << "UNKNOWN YAML EVENT (value): " << event->type
                          << "\nVALUE: " << event->data.scalar.value
                          << std::endl;
                break;
          }
          break;

      case ERROR:
      case STOP:
          break;
    }
    return (s->state == ERROR) ? 0 : 1;
}

bool get_bool_from_data(yaml_char_t *data, int len)
{
    bool rv = false;
    char *value = new char[len+1];

    for (int i=0; i < len; ++i) {
        value[i] = tolower(data[i]);
    }
    value[len] = '\0';

    if ((strcmp(value,"yes") == 0) ||
        (strcmp(value,"true") == 0) ||
        (strcmp(value,"on") == 0) ||
        (atoi(value) != 0)) {
        rv = true;
    }
    delete[] value;
    return rv;
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
    try {
        state.state = START;
        do {
            if (!yaml_parser_parse(&parser, &event))
                throw;

            if (!consume_event(&state, &event))
                throw;

            if (state.accepted) {
                if (state.key == "pair_newton_on") {
                    test_config.pair_newton_on =
                        get_bool_from_data(event.data.scalar.value,
                                           event.data.scalar.length);
                } else if (state.key == "pair_newton_off") {
                    test_config.pair_newton_off =
                        get_bool_from_data(event.data.scalar.value,
                                           event.data.scalar.length);
                } else if (state.key == "global_stress") {
                    test_config.has_global_stress =
                        get_bool_from_data(event.data.scalar.value,
                                           event.data.scalar.length);
                } else if (state.key == "global_energy") {
                    test_config.has_global_energy =
                        get_bool_from_data(event.data.scalar.value,
                                           event.data.scalar.length);
                } else if (state.key == "atom_stress") {
                    test_config.has_atom_stress =
                        get_bool_from_data(event.data.scalar.value,
                                           event.data.scalar.length);
                } else if (state.key == "atom_energy") {
                    test_config.has_atom_energy =
                        get_bool_from_data(event.data.scalar.value,
                                           event.data.scalar.length);
                } else if (state.key == "pre_commands") {
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
                } else if (state.key == "init_coul") {
                    test_config.init_coul = atof((char *)event.data.scalar.value);
                } else if (state.key == "init_stress") {
                    stress_t stress;
                    sscanf((char *)event.data.scalar.value,
                           "%lg %lg %lg %lg %lg %lg",
                           &stress.xx, &stress.yy, &stress.zz,
                           &stress.xy, &stress.xz, &stress.yz);
                    test_config.init_stress = stress;
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
                } else std::cerr << "Ignoring unknown key/value pair: " << state.key
                                 << " = " << event.data.scalar.value << std::endl;
            }
            yaml_event_delete(&event);
        } while (state.state != STOP);
    } catch (...) {
        retval = 1;
    }

    yaml_parser_delete(&parser);
    fclose(fp);
    return retval;
}

void emit_boolean(yaml_emitter_t *emitter, yaml_event_t *event,
                  const std::string &key, bool value)
{
    yaml_scalar_event_initialize(event, NULL,
                                 (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *) key.c_str(),
                                 key.size(), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(emitter, event);
    std::string out = value ? "true" : "false";
    yaml_scalar_event_initialize(event, NULL,
                                 (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *) out.c_str(),
                                 out.size(), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(emitter, event);
}

void emit_double(yaml_emitter_t *emitter, yaml_event_t *event,
                 const std::string &key, const double value)
{
    yaml_scalar_event_initialize(event, NULL,
                                 (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *) key.c_str(),
                                 key.size(), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(emitter, event);
    char buf[256];
    snprintf(buf,256,"%.15g",value);
    yaml_scalar_event_initialize(event, NULL,
                                 (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *)buf,
                                 strlen(buf), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(emitter, event);
}

void emit_string(yaml_emitter_t *emitter, yaml_event_t *event,
                 const std::string &key, const std::string &value)
{
    yaml_scalar_event_initialize(event, NULL,
                                 (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *) key.c_str(),
                                 key.size(), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(emitter, event);
    yaml_scalar_event_initialize(event, NULL,
                                 (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *)value.c_str(),
                                 value.size(), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(emitter, event);
}

void emit_block(yaml_emitter_t *emitter, yaml_event_t *event,
                const std::string &key, const std::string &value)
{
    yaml_scalar_event_initialize(event, NULL,
                                 (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *) key.c_str(),
                                 key.size(), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(emitter, event);
    yaml_scalar_event_initialize(event, NULL,
                                 (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *)value.c_str(),
                                 value.size(), 1, 0,
                                 YAML_LITERAL_SCALAR_STYLE);
    yaml_emitter_emit(emitter, event);
}

void generate(const char *outfile)
{
    FILE *fp;
    yaml_emitter_t emitter;
    yaml_event_t   event;

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

    emit_boolean(&emitter, &event,
                 "pair_newton_on", test_config.pair_newton_on);
    emit_boolean(&emitter, &event,
                 "pair_newton_off", test_config.pair_newton_off);
    emit_boolean(&emitter, &event,
                 "global_stress", test_config.has_global_stress);
    emit_boolean(&emitter, &event,
                 "global_energy", test_config.has_global_energy);
    emit_boolean(&emitter, &event,
                 "atom_stress", test_config.has_atom_stress);
    emit_boolean(&emitter, &event,
                 "atom_energy", test_config.has_atom_energy);

    std::string block("");
    for (std::size_t i=0; i < test_config.pre_commands.size(); ++i) {
        block += test_config.pre_commands[i] + "\n";
    }
    emit_block(&emitter, &event, "pre_commands", block);
    block.clear();
    for (std::size_t i=0; i < test_config.post_commands.size(); ++i) {
        block += test_config.post_commands[i] + "\n";
    }
    emit_block(&emitter, &event, "post_commands", block);
    emit_string(&emitter, &event, "input_file", test_config.input_file);

    emit_string(&emitter, &event, "pair_style", test_config.pair_style);
    block.clear();
    for (std::size_t i=0; i < test_config.pair_coeff.size(); ++i) {
        block += test_config.pair_coeff[i] + "\n";
    }
    emit_block(&emitter, &event, "pair_coeff", block);
    
    const char *args[] = {"MolPairStyle", "-log", "none", "-echo", "screen", "-nocite" };
    char **argv = (char **)args;
    int argc = sizeof(args)/sizeof(char *);
    LAMMPS_NS::LAMMPS *lmp = init_lammps(argc,argv,test_config);

    const int natoms = lmp->atom->natoms;
    const int bufsize = 256;
    char buf[bufsize];
    snprintf(buf,bufsize,"%d", natoms);
    emit_string(&emitter, &event, "natoms", buf);
    emit_double(&emitter, &event, "init_vdwl", lmp->force->pair->eng_vdwl);
    emit_double(&emitter, &event, "init_coul", lmp->force->pair->eng_coul);

    double *stress = lmp->force->pair->virial;
    snprintf(buf,bufsize,"%.15g %.15g %.15g %.15g %.15g %.15g",
             stress[0],stress[1],stress[2],stress[3],stress[4],stress[5]);
    emit_string(&emitter, &event, "init_stress", buf);

    block = "";
    double **f = lmp->atom->f;
    LAMMPS_NS::tagint *tag = lmp->atom->tag;
    for (int i=0; i < natoms; ++i) {
        snprintf(buf,bufsize,"%d %.15g %.15g %.15g\n",
                 (int)tag[i], f[i][0], f[i][1], f[i][2]);
        block += std::string(buf);
    }
    emit_block(&emitter, &event, "init_forces", block);

    delete lmp;
    yaml_mapping_end_event_initialize(&event);
    yaml_emitter_emit(&emitter, &event);
    yaml_document_end_event_initialize(&event, 0);
    yaml_emitter_emit(&emitter, &event);
    yaml_stream_end_event_initialize(&event);
    yaml_emitter_emit(&emitter, &event);
    yaml_emitter_delete(&emitter);
    fclose(fp);
    return;
}

#if 0
class MolPairStyle : public ::testing::Test
{
protected:
    LAMMPS_NS::LAMMPS *lmp;
    MolPairStyle() {
        const char *args[] = {"LAMMPS_test",
                              "-log", "none",
                              "-echo", "screen",
                              "-nocite" };
        char **argv = (char **)args;
        int argc = sizeof(args)/sizeof(char *);

        int flag;
        MPI_Initialized(&flag);
        if (!flag) MPI_Init(&argc,&argv);

        ::testing::internal::CaptureStdout();
        lmp = new LAMMPS_NS::LAMMPS(argc, argv, MPI_COMM_WORLD);
        ::testing::internal::GetCapturedStdout();
    }
    ~MolPairStyle() override {}

    void SetUp() override {
        lmp->input->one("clear");
    }
    void TearDown() override {
    }
};

TEST_F(MolPairStyle, initial) {
    double **f=lmp->atom->f;
    double **v=lmp->atom->v;
    const int nlocal = lmp->atom->nlocal;

    // abort if running in parallel
    ASSERT_EQ(lmp->atom->natoms,nlocal);

    lmp->input->one("velocity all set 0.0 0.0 0.0");
    lmp->input->one("run 0 post no");
    for (int i=0; i < nlocal; ++i) {
        EXPECT_DOUBLE_EQ(f[i][0],0.0);
        EXPECT_DOUBLE_EQ(f[i][1],0.0);
        EXPECT_DOUBLE_EQ(f[i][2],0.0);
        EXPECT_DOUBLE_EQ(v[i][0],0.0);
        EXPECT_DOUBLE_EQ(v[i][1],0.0);
        EXPECT_DOUBLE_EQ(v[i][2],0.0);
    }
};

#endif
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    if ((argc != 2) && (argc != 4)) {
        std::cerr << "usage: " << argv[0] << " <testfile.yaml> "
            "[--gen <newfile.yaml>]" << std::endl;
        return 1;
    }
    if (parse(argv[1]))
        return 2;

    if (argc == 4) {
        if (strcmp(argv[2],"--gen") == 0) {
            generate(argv[3]);
            return 0;
        } else {
            std::cerr << "usage: " << argv[0] << " <testfile.yaml> "
                "[--gen <newfile.yaml>]" << std::endl;
            return 1;
        }
    }
//    return RUN_ALL_TESTS();
    return 0;
}
