#include <iostream>
#include <assert.h>
#include <python3.4m/Python.h>
#include <vector>
#include <fstream>

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

void parameterInitialize(Parameters parameters, std::vector<double> *nd_vector, std::vector<double> *ccd_vector) {

    PyObject *module = PyImport_ImportModule("parameters.search");
    assert(module != NULL);

    PyObject *klass = PyObject_GetAttrString(module, "ParameterSearch");
    assert(klass != NULL);

    PyObject *instance = PyObject_CallFunction(klass, "dddffd",
                                               parameters.number_of_vertices,
                                               parameters.max_degree_bound,
                                               parameters.average_degree_target,
                                               parameters.max_clustering_coefficient_target,
                                               parameters.global_clustering_coefficient_target, 1);
    assert(instance != NULL);

    PyObject *result_nd = PyObject_CallMethod(instance, "run_nd", "(iiiffi)",
                                              parameters.number_of_vertices,
                                              parameters.max_degree_bound,
                                              parameters.average_degree_target,
                                              parameters.max_clustering_coefficient_target,
                                              parameters.global_clustering_coefficient_target, 1);
    assert(result_nd != NULL);

    PyObject *result_ccd = PyObject_CallMethod(instance, "run_ccd", "(iiiffi)",
                                               parameters.number_of_vertices,
                                               parameters.max_degree_bound,
                                               parameters.average_degree_target,
                                               parameters.max_clustering_coefficient_target,
                                               parameters.global_clustering_coefficient_target, 1);
    assert(result_ccd != NULL);

    spd::get("logger")->info("Python finished executing, start copying list to vector");
    *nd_vector = pyListToVector(result_nd);
    spd::get("logger")->info("ND vector copied");
    *ccd_vector = pyListToVector(result_ccd);
    spd::get("logger")->info("CD vector copied");
}

void singleBenchmarkGpu(std::ofstream &outfile_timing, std::ofstream &outfile_edges, Parameters parameters) {

    std::vector<double> nd_vector;
    std::vector<double> ccd_vector;

    parameterInitialize(parameters, &nd_vector, &ccd_vector);

    // Get pointer to start of vector array
    double *nd = &nd_vector[0];
    double *cd = &ccd_vector[0];

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

    // Setup timer
    std::chrono::time_point <std::chrono::system_clock> start, end;
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
    outfile_timing << parameters.number_of_vertices << "," << elapsed_seconds.count() << ",";
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

void singleBenchmarkSeq(std::ofstream &outfile_timing, std::ofstream &outfile_edges, Parameters parameters) {

    std::vector<double> nd_vector;
    std::vector<double> ccd_vector;

    parameterInitialize(parameters, &nd_vector, &ccd_vector);

    // Get pointer to start of vector array
    double *nd = &nd_vector[0];
    double *cd = &ccd_vector[0];

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

    // Setup timer
    std::chrono::time_point <std::chrono::system_clock> start, end;
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
    outfile_timing << parameters.number_of_vertices << "," << elapsed_seconds.count() << ",";
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

void benchmarkSeq(int max_nnodes, int run_count, Parameters parameters) {

    std::ofstream outfile_timing;
    outfile_timing.open("cpu_bench.csv");
    outfile_timing << "nnodes,phase_one_seconds,phase_two_seconds" << std::endl;

    std::ofstream outfile_edgelist;
    outfile_edgelist.open("edge_list.csv");

    int batch_size;
    for (int i = 1; i < run_count; ++i) {
        batch_size = (max_nnodes / run_count) * i;
        singleBenchmarkSeq(outfile_timing, outfile_edgelist, parameters);
    }

    outfile_timing.close();
}

void benchmarkGpu(int max_nnodes, int run_count, Parameters parameters) {

    std::ofstream outfile_timing;
    outfile_timing.open("gpu_bench.csv");
    outfile_timing << "nnodes,phase_one_seconds,phase_two_seconds" << std::endl;

    std::ofstream outfile_edgelist;
    outfile_edgelist.open("edge_list.csv");

    int batch_size;
    for (int i = 1; i < run_count; ++i) {
        batch_size = (max_nnodes / run_count) * i;
        singleBenchmarkGpu(outfile_timing, outfile_edgelist, parameters);
    }

    outfile_timing.close();
}

int main() {

    setupLogger();

    setupEnvironment();
    spd::get("logger")->info("Environment setup");
    spd::get("logger")->info("Starting Python");
    Py_Initialize();

    std::ofstream outfile_timing;
    outfile_timing.open("gpu_bench.csv");
    outfile_timing << "nnodes;phase_one_seconds;phase_two_seconds \n";

    std::ofstream outfile_edgelist;
    outfile_edgelist.open("edge_list.csv");

    Parameters parameters = {
            1000000,
            100000,
            32,
            0.95,
            0.15
    };
    singleBenchmarkGpu(outfile_timing, outfile_edgelist, parameters);

    outfile_timing.close();
    outfile_edgelist.close();

    Py_Finalize();

    return 0;
}

