#include <iostream>
#include <assert.h>
#include <python3.4m/Python.h>
#include <vector>

#include "BterPhases.h"

using namespace bter;

void setupLogger() {

}

void setupEnvironment() {
    char *dirname = get_current_dir_name();
    std::string dir = std::string(dirname);
    const char *path = dir.c_str();

    setenv("PYTHONPATH", path, 1);
    std::cout << "Python path set to: " << std::endl;
    std::cout << path << std::endl;
}

// PyObject -> Vector
std::vector<double> pyListToVector(PyObject *incoming) {
    std::vector<double> data;
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
    setupEnvironment();

    Py_Initialize();

    PyObject *module = PyImport_ImportModule("parameters.search");
    assert(module != NULL);

    PyObject *klass = PyObject_GetAttrString(module, "ParameterSearch");
    assert(klass != NULL);

    PyObject *instance = PyObject_CallFunction(klass, "dddffd", 20, 10, 3, 0.95, 0.15, 1);
    assert(instance != NULL);

    PyObject *result_nd = PyObject_CallMethod(instance, "run_nd", "(iiiffi)", 20, 10, 3, 0.95, 0.15, 1);
    assert(result_nd != NULL);

    PyObject *result_ccd = PyObject_CallMethod(instance, "run_ccd", "(iiiffi)", 20, 10, 3, 0.95, 0.15, 1);
    assert(result_ccd != NULL);

    std::vector<double> nd_vector = pyListToVector(result_nd);
    std::vector<double> ccd_vector = pyListToVector(result_ccd);

    Py_Finalize();

    // Get pointer to start of vector array
    double *nd = &nd_vector[0];
    double *cd = &ccd_vector[0];

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

    BTERSetupResult bterSetupResult{
            id, ndprime,
            wd, rdfill, ndfill, wg, ig, bg, ng
    };

    BTERSetup bterSetup(nd, cd, &beta, dmax, &bterSetupResult);
    bterSetup.run();

    // Phases
    BterPhases bterPhases(&bterSetupResult, dmax, nd, cd);
    bterPhases.computeSamples();

    int *phase_one_i = new int[bterPhases.bterSamples.s1];
    int *phase_one_j = new int[bterPhases.bterSamples.s1];
    bterPhases.phaseOneGpu(phase_one_i, phase_one_j);

    int *phase_two_i = new int[bterPhases.bterSamples.s2];
    int *phase_two_j = new int[bterPhases.bterSamples.s2];
    bterPhases.phaseTwoGpu(phase_two_i, phase_two_j);

    std::cout << "\nPHASE ONE: \n";
    for (int i = 0; i < bterPhases.bterSamples.s1; ++i)
        std::cout << "[" << phase_one_i[i] << " - " << phase_one_j[i] << "] ";

    std::cout << std::endl;

    std::cout << "\nPHASE TWO: \n";
    for (int i = 0; i < bterPhases.bterSamples.s2; ++i)
        std::cout << "[" << phase_two_i[i] << " - " << phase_two_j[i] << "] ";

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

//    std::chrono::time_point<std::chrono::system_clock> start, end;
//    start = std::chrono::system_clock::now();
//    std::cout << "f(42) = " << fibonacci(42) << '\n';
//    end = std::chrono::system_clock::now();
//
//    std::chrono::duration<double> elapsed_seconds = end-start;
//    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
//
//    std::cout << "finished computation at " << std::ctime(&end_time)
//              << "elapsed time: " << elapsed_seconds. << "s\n";
    return 0;
}

