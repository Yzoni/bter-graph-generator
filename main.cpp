#include <iostream>
#include <assert.h>
#include <python3.4m/Python.h>
#include <vector>

#include "spdlog/spdlog.h"
#include "BterPhases.h"

//#define VALUES

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


int main() {

    setupLogger();

    /*
     * Get parameters using python
     */

    spd::get("logger")->info("Setting up environment");
    setupEnvironment();
    spd::get("logger")->info("Environment setup");

    spd::get("logger")->info("Starting Python");
    Py_Initialize();

    PyObject *module = PyImport_ImportModule("parameters.search");
    assert(module != NULL);

    PyObject *klass = PyObject_GetAttrString(module, "ParameterSearch");
    assert(klass != NULL);

    PyObject *instance = PyObject_CallFunction(klass, "dddffd", 100000, 40000, 32, 0.95, 0.15, 1);
    assert(instance != NULL);

    PyObject *result_nd = PyObject_CallMethod(instance, "run_nd", "(iiiffi)", 100000, 40000, 32, 0.95, 0.15, 1);
    assert(result_nd != NULL);

    PyObject *result_ccd = PyObject_CallMethod(instance, "run_ccd", "(iiiffi)", 100000, 40000, 32, 0.95, 0.15, 1);
    assert(result_ccd != NULL);

    spd::get("logger")->info("Python finished executing, start copying list to vector");
    std::vector<double> nd_vector = pyListToVector(result_nd);
    spd::get("logger")->info("ND vector copied");
    std::vector<double> ccd_vector = pyListToVector(result_ccd);
    spd::get("logger")->info("CD vector copied");

    Py_Finalize();

    // Get pointer to start of vector array
    double *nd = &nd_vector[0];
    double *cd = &ccd_vector[0];

    spd::get("logger")->info("Allocate new arrays");

    /*
     * Start C++ library
     */
    double beta = 1;
    int dmax = nd_vector.size();

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
    BterPhases bterPhases(&bterSetupResult, dmax, nd, cd);
    bterPhases.computeSamples();
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    spd::get("logger")->info("Finished computing phase, took {} seconds", elapsed_seconds.count());

    // PHASE ONE
    spd::get("logger")->info("Start phase one");
    start = std::chrono::system_clock::now();
    int *phase_one_i = new int[bterPhases.bterSamples.s1];
    int *phase_one_j = new int[bterPhases.bterSamples.s1];
    bterPhases.phaseOneGpu(phase_one_i, phase_one_j);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    spd::get("logger")->info("Finished phase one, took {} seconds", elapsed_seconds.count());

    // PHASE TWO
    spd::get("logger")->info("Start phase two");
    start = std::chrono::system_clock::now();
    int *phase_two_i = new int[bterPhases.bterSamples.s2];
    int *phase_two_j = new int[bterPhases.bterSamples.s2];
    bterPhases.phaseTwoGpu(phase_two_i, phase_two_j);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    spd::get("logger")->info("Finished phase two, took {} seconds", elapsed_seconds.count());

#ifdef VALUES
    std::cout << "\nPHASE ONE: \n";
    for (int i = 0; i < bterPhases.bterSamples.s1; ++i)
        std::cout << "[" << phase_one_i[i] << " - " << phase_one_j[i] << "] ";

    std::cout << std::endl;

    std::cout << "\nPHASE TWO: \n";
    for (int i = 0; i < bterPhases.bterSamples.s2; ++i)
        std::cout << "[" << phase_two_i[i] << " - " << phase_two_j[i] << "] ";
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

    return 0;
}

