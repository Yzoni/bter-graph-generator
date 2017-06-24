#include <iostream>
#include <assert.h>
#include <python3.5m/Python.h>
#include <vector>
#include <fstream>
#include <getopt.h>

#include "spdlog/spdlog.h"
#include "BterPhasesGpu.h"
#include "BterPhasesSeq.h"

#define VALUES

using namespace bter;
namespace spd = spdlog;

void setupLogger() {
    auto console = spd::stdout_color_mt("logger");
    console->info("Logging started");
}

void setupEnvironment() {
    char *dirname = get_current_dir_name();
    std::string dir = std::string(dirname);
    const char *path = dir.c_str();

    setenv("PYTHONPATH", path, 1);
    spd::get("logger")->info("Python path set to {}", path);
}

struct Parameters {
    int number_of_vertices;
    int max_degree_bound;
    int average_degree_target;
    float max_clustering_coefficient_target;
    float global_clustering_coefficient_target;

    Parameters(int number_of_vertices,
               int max_degree_bound,
               int average_degree_target,
               float max_clustering_coefficient_target,
               float global_clustering_coefficient_target) :
            number_of_vertices(number_of_vertices),
            max_degree_bound(max_degree_bound),
            average_degree_target(average_degree_target),
            max_clustering_coefficient_target(max_clustering_coefficient_target),
            global_clustering_coefficient_target(global_clustering_coefficient_target) {}
};

// PyObject -> Vector
std::vector<double> pyListToVector(PyObject *incoming) {
    std::vector<double> data;
    data.reserve(100000);
    if (PyList_Check(incoming)) {
        for (Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
            PyObject *value = PyList_GetItem(incoming, i);
            data.push_back(PyFloat_AsDouble(value));
        }
    } else {
        std::cerr << "Passed PyObject pointer was not a list or tuple!" << std::endl;
    }

    return data;
}

void parameterInitialize(Parameters *parameters, std::vector<double> *nd_vector, std::vector<double> *ccd_vector) {

    PyObject *module = PyImport_ImportModule("parameters.search");
    assert(module != NULL);

    PyObject *klass = PyObject_GetAttrString(module, "ParameterSearch");
    assert(klass != NULL);

    PyObject *instance = PyObject_CallFunction(klass, "dddffd",
                                               parameters->number_of_vertices,
                                               parameters->max_degree_bound,
                                               parameters->average_degree_target,
                                               parameters->max_clustering_coefficient_target,
                                               parameters->global_clustering_coefficient_target, 1);
    assert(instance != NULL);

    PyObject *result_nd = PyObject_CallMethod(instance, "run_nd", "(iiiffi)",
                                              parameters->number_of_vertices,
                                              parameters->max_degree_bound,
                                              parameters->average_degree_target,
                                              parameters->max_clustering_coefficient_target,
                                              parameters->global_clustering_coefficient_target, 1);
    assert(result_nd != NULL);

    PyObject *result_ccd = PyObject_CallMethod(instance, "run_ccd", "(iiiffi)",
                                               parameters->number_of_vertices,
                                               parameters->max_degree_bound,
                                               parameters->average_degree_target,
                                               parameters->max_clustering_coefficient_target,
                                               parameters->global_clustering_coefficient_target, 1);
    assert(result_ccd != NULL);

    spd::get("logger")->info("Python finished executing, start copying list to vector");
    *nd_vector = pyListToVector(result_nd);
    spd::get("logger")->info("ND vector copied");
    *ccd_vector = pyListToVector(result_ccd);
    spd::get("logger")->info("CD vector copied");
}

void singleBenchmarkGpu(std::ofstream &outfile_timing, std::ofstream &outfile_edges, Parameters *parameters) {

    std::vector<double> nd_vector;
    std::vector<double> ccd_vector;

    parameterInitialize(parameters, &nd_vector, &ccd_vector);

    // Get pointer to start of vector array
    double *nd = &nd_vector[0];
    double *cd = &ccd_vector[0];

    std::cout << "nd_vector: \n";
    for (std::vector<double>::const_iterator i = nd_vector.begin(); i != nd_vector.end(); ++i)
        std::cout << *i << ' ';


    std::cout << "\nccd_vector: \n";
    for (std::vector<double>::const_iterator i = ccd_vector.begin(); i != ccd_vector.end(); ++i)
        std::cout << *i << ' ';
    std::cout << "\n";

    spd::get("logger")->info("Allocate new arrays");

    /*
     * Start C++ library
     */
    double beta = 1;
    size_t dmax = nd_vector.size();

    int *id = new int[dmax]{};
    double *wd = new double[dmax]{};
    double *rdfill = new double[dmax]{};
    double *ndfill = new double[dmax]{};
    double *wg = new double[dmax]{};
    double *ig = new double[dmax]{};
    double *bg = new double[dmax]{};
    double *ng = new double[dmax]{};
    int *ndprime = new int[dmax]{};

    spd::get("logger")->info("Start BTER setup");
    BTERSetupResult bterSetupResult{
            id, ndprime,
            wd, rdfill, ndfill, wg, ig, bg, ng
    };

    BTERSetup bterSetup(nd, cd, &beta, dmax, &bterSetupResult);
    bterSetup.run();
    spd::get("logger")->info("Finished BTER setup");

    std::cout << "BG: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.bg[i] << " ";
    std::cout << "\n";
    std::cout << "ig: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.ig[i] << " ";
    std::cout << "\n";
    std::cout << "ng: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.ng[i] << " ";
    std::cout << "\n";
    std::cout << "wg: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.wg[i] << " ";
    std::cout << "\n";
    std::cout << "wd: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.wd[i] << " ";
    std::cout << "\n";
    std::cout << "rdfill: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.rdfill[i] << " ";
    std::cout << "\n";
    std::cout << "ndfill: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.ndfill[i] << " ";
    std::cout << "\n";


    // Setup timer
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;

    // COMPUTE PHASES
    spd::get("logger")->info("Start computing phases");
    start = std::chrono::system_clock::now();
    BterPhasesGpu bterPhasesGpu(&bterSetupResult, dmax, nd, cd);
    bterPhasesGpu.computeSamples();
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    spd::get("logger")->info("Finished computing phase, took {} seconds", elapsed_seconds.count());

    // PHASE ONE
    spd::get("logger")->info("Start phase one");
    start = std::chrono::system_clock::now();
    int *phase_one_i = new int[bterPhasesGpu.bterSamples.s1];
    int *phase_one_j = new int[bterPhasesGpu.bterSamples.s1];
    bterPhasesGpu.phaseOne(phase_one_i, phase_one_j);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    outfile_timing << parameters->number_of_vertices << "," << elapsed_seconds.count() << ",";
    spd::get("logger")->info("Finished phase one, took {} seconds", elapsed_seconds.count());

    // PHASE TWO
    spd::get("logger")->info("Start phase two");
    start = std::chrono::system_clock::now();
    int *phase_two_i = new int[bterPhasesGpu.bterSamples.s2];
    int *phase_two_j = new int[bterPhasesGpu.bterSamples.s2];
    bterPhasesGpu.phaseTwo(phase_two_i, phase_two_j);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    outfile_timing << elapsed_seconds.count() << std::endl;
    spd::get("logger")->info("Finished phase two, took {} seconds", elapsed_seconds.count());

#ifdef VALUES
    spd::get("logger")->info("Phase one write vertices to file");
    for (int i = 0; i < bterPhasesGpu.bterSamples.s1; ++i)
        outfile_edges << phase_one_i[i] << ";" << phase_one_j[i] << "\n";

    std::cout << std::endl;

    spd::get("logger")->info("Phase two write vertices to file");
    for (int i = 0; i < bterPhasesGpu.bterSamples.s2; ++i)
        outfile_edges << phase_two_i[i] << ";" << phase_two_j[i] << " \n";
#endif

    spd::get("logger")->info("Freeing memory");
    delete[] id;
    delete[] wd;
    delete[] rdfill;
    delete[] ndfill;
    delete[] wg;
    delete[] ig;
    delete[] bg;
    delete[] ng;
    delete[] ndprime;

    delete[] phase_one_i;
    delete[] phase_one_j;

    delete[] phase_two_i;
    delete[] phase_two_j;
}

void singleBenchmarkSeq(std::ofstream &outfile_timing, std::ofstream &outfile_edges, Parameters *parameters) {

    std::vector<double> nd_vector;
    std::vector<double> ccd_vector;

    parameterInitialize(parameters, &nd_vector, &ccd_vector);

    // Get pointer to start of vector array
    double *nd = &nd_vector[0];
    double *cd = &ccd_vector[0];

    std::cout << "nd_vector: \n";
    for (std::vector<double>::const_iterator i = nd_vector.begin(); i != nd_vector.end(); ++i)
        std::cout << *i << ' ';


    std::cout << "\nccd_vector: \n";
    for (std::vector<double>::const_iterator i = ccd_vector.begin(); i != ccd_vector.end(); ++i)
        std::cout << *i << ' ';
    std::cout << "\n";

    spd::get("logger")->info("Allocate new arrays");

    /*
     * Start C++ library
     */
    double beta = 1;
    size_t dmax = nd_vector.size();

    int *id = new int[dmax]{};
    double *wd = new double[dmax]{};
    double *rdfill = new double[dmax]{};
    double *ndfill = new double[dmax]{};
    double *wg = new double[dmax]{};
    double *ig = new double[dmax]{};
    double *bg = new double[dmax]{};
    double *ng = new double[dmax]{};
    int *ndprime = new int[dmax]{};

    spd::get("logger")->info("Start BTER setup");
    BTERSetupResult bterSetupResult{
            id, ndprime,
            wd, rdfill, ndfill, wg, ig, bg, ng
    };

    BTERSetup bterSetup(nd, cd, &beta, dmax, &bterSetupResult);
    bterSetup.run();

    std::cout << "BG: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.bg[i] << " ";
    std::cout << "\n";
    std::cout << "ig: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.ig[i] << " ";
    std::cout << "\n";
    std::cout << "ng: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.ng[i] << " ";
    std::cout << "\n";
    std::cout << "wg: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.wg[i] << " ";
    std::cout << "\n";
    std::cout << "wd: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.wd[i] << " ";
    std::cout << "\n";
    std::cout << "rdfill: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.rdfill[i] << " ";
    std::cout << "\n";
    std::cout << "ndfill: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.ndfill[i] << " ";
    std::cout << "\n";
    std::cout << "ndprime: \n";
    for (int i = 0; i < dmax; ++i) std::cout << bterSetupResult.ndprime[i] << " ";
    std::cout << "\n";
    spd::get("logger")->info("Finished BTER setup");
    double nd_sum = std::accumulate(nd, &nd[dmax], 0.0, std::plus<double>());
    spd::get("logger")->info("Desired number nodes: {}", nd_sum);
    spd::get("logger")->info("Dmax: {}", dmax);
    double bg_sum = std::accumulate(bg, &bg[dmax], 0.0, std::plus<double>());
    spd::get("logger")->info("# Blocks: {}", bg_sum);
    double wg_sum = std::accumulate(wg, &wg[dmax], 0.0, std::plus<double>());
    spd::get("logger")->info("Phase 1 weigth: {}", wg_sum);
    double wd_sum = std::accumulate(wd, &wd[dmax], 0.0, std::plus<double>());
    spd::get("logger")->info("Phase 2 weigth: {}", wd_sum);

    // Setup timer
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;

    // COMPUTE PHASES
    spd::get("logger")->info("Start computing phases");
    start = std::chrono::system_clock::now();
    BterPhasesSeq bterPhasesSeq(&bterSetupResult, dmax, nd, cd);
    bterPhasesSeq.computeSamples();
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    spd::get("logger")->info("Finished computing phase, took {} seconds", elapsed_seconds.count());

    // PHASE ONE
    spd::get("logger")->info("Start phase one");
    start = std::chrono::system_clock::now();
    int *phase_one_i = new int[bterPhasesSeq.bterSamples.s1];
    int *phase_one_j = new int[bterPhasesSeq.bterSamples.s1];
    bterPhasesSeq.phaseOne(phase_one_i, phase_one_j);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    outfile_timing << parameters->number_of_vertices << "," << elapsed_seconds.count() << ",";
    spd::get("logger")->info("Finished phase one, took {} seconds", elapsed_seconds.count());

    // PHASE TWO
    spd::get("logger")->info("Start phase two");
    start = std::chrono::system_clock::now();
    int *phase_two_i = new int[bterPhasesSeq.bterSamples.s2];
    int *phase_two_j = new int[bterPhasesSeq.bterSamples.s2];
    bterPhasesSeq.phaseTwo(phase_two_i, phase_two_j);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    outfile_timing << elapsed_seconds.count() << std::endl;
    spd::get("logger")->info("Finished phase two, took {} seconds", elapsed_seconds.count());

#ifdef VALUES
    spd::get("logger")->info("Phase one write vertices to file");
    for (int i = 0; i < bterPhasesSeq.bterSamples.s1; ++i)
        outfile_edges << phase_one_i[i] << ";" << phase_one_j[i] << "\n";

    std::cout << std::endl;

    spd::get("logger")->info("Phase two write vertices to file");
    for (int i = 0; i < bterPhasesSeq.bterSamples.s2; ++i)
        outfile_edges << phase_two_i[i] << ";" << phase_two_j[i] << " \n";
#endif

    spd::get("logger")->info("Freeing memory");
    delete[] id;
    delete[] wd;
    delete[] rdfill;
    delete[] ndfill;
    delete[] wg;
    delete[] ig;
    delete[] bg;
    delete[] ng;
    delete[] ndprime;

    delete[] phase_one_i;
    delete[] phase_one_j;

    delete[] phase_two_i;
    delete[] phase_two_j;
}

int main(int argc, char *argv[]) {

    int gpu_flag = 0;

    int number_of_nodes = -1;
    int max_degree_bound = -1;
    int average_degree_target = -1;
    float max_clustering_coefficient_target = -1;
    float global_clustering_coefficient_target = -1;

    int c;
    int option_count = 7;
    while (1) {

        int option_index = 0;
        static struct option long_options[] = {
                {"number_of_nodes",                      required_argument, 0, 'n'},
                {"max_degree_bound",                     required_argument, 0, 'b'},
                {"average_degree_target",                required_argument, 0, 'd'},
                {"max_clustering_coefficient_target",    required_argument, 0, 'm'},
                {"global_clustering_coefficient_target", required_argument, 0, 'c'},
                {"gpu", no_argument, 0, 'g'},
                {"help", no_argument, 0, 'h'}
        };
        c = getopt_long(argc, argv, "n:b:d:m:c:gh", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 'n':
                number_of_nodes = atoi(optarg);
                break;
            case 'b':
                max_degree_bound = atoi(optarg);
                break;
            case 'd':
                average_degree_target = atoi(optarg);
                break;
            case 'm':
                max_clustering_coefficient_target = std::stof(optarg);
                break;
            case 'c':
                global_clustering_coefficient_target = std::stof(optarg);
                break;
            case 'g':
                gpu_flag = 1;
                break;
            case 'h':
                fprintf(stdout,"Valid options are: \n");
                for (int i = 0; i < option_count; ++i) {
                    fprintf(stdout, "-%c, --%s \n", long_options[i].val, long_options[i].name);
                }
            case '?':
               if (isprint(optopt)) {
                   fprintf(stderr, "Unknown option `-%c'.\n", optopt);
               }

                return 1;
            default:
                abort();
        }
    }

    if (number_of_nodes == -1 || max_degree_bound == -1 || average_degree_target == -1 ||
        max_clustering_coefficient_target == -1 || global_clustering_coefficient_target == -1) {
        printf("Parameters are mandatory");
        exit(1);
    }

    setupLogger();

    setupEnvironment();

    Parameters *parameters = new Parameters(number_of_nodes,
                                            max_degree_bound,
                                            average_degree_target,
                                            max_clustering_coefficient_target,
                                            global_clustering_coefficient_target);


    spd::get("logger")->info("Environment setup");
    spd::get("logger")->info("Starting Python");
    Py_Initialize();

    std::ofstream outfile_timing;
    outfile_timing.open("gpu_bench.csv");
    outfile_timing << "nnodes;phase_one_seconds;phase_two_seconds \n";

    std::ofstream outfile_edgelist;
    outfile_edgelist.open("edge_list.csv");

    if (gpu_flag == 1) {
        singleBenchmarkGpu(outfile_timing, outfile_edgelist, parameters);
    } else {
        singleBenchmarkSeq(outfile_timing, outfile_edgelist, parameters);
    }

    outfile_timing.close();
    outfile_edgelist.close();

    Py_Finalize();

    return 0;
}

