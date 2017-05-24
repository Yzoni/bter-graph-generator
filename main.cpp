#include <iostream>
#include <assert.h>
#include <python3.5/Python.h>
#include <vector>

#include "BterPhases.h"

using namespace bter;

void setupLogger() {

}

void setupEnvironment() {
    char *dirname = get_current_dir_name();
    std::string dir = std::string(dirname);
    std::string full_dir = (dir + std::string("/untitled"));
    const char *path = full_dir.c_str();

    setenv("PYTHONPATH", path, 1);
}

// PyObject -> Vector
std::vector<double> listTupleToVector(PyObject *incoming) {
    std::vector<double> data;
    if (PyTuple_Check(incoming)) {
        for (Py_ssize_t i = 0; i < PyTuple_Size(incoming); i++) {
            PyObject *value = PyTuple_GetItem(incoming, i);
            data.push_back(PyFloat_AsDouble(value));
        }
    } else {
        if (PyList_Check(incoming)) {
            for (Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
                PyObject *value = PyList_GetItem(incoming, i);
                data.push_back(PyFloat_AsDouble(value));
            }
        } else {
            std::cerr << "Passed PyObject pointer was not a list or tuple!" << std::endl;
        }
    }
    return data;
}


int main() {

    setupLogger();

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

    std::vector<double> nd_vector = listTupleToVector(result_nd);
    std::vector<double> ccd_vector = listTupleToVector(result_ccd);

    Py_Finalize();

    double* nd = &nd_vector[0];
    double* ccd = &ccd_vector[0];

    double beta = 1;
    int dmax = nd.size();

    int id[dmax] = new int[dmax];
    double wd[DMAX] = new double[dmax];
    double rdfill[DMAX] = new double[dmax];
    double ndfill[DMAX] = new double[dmax];
    double wg[DMAX] = new double[dmax];
    double ig[DMAX] = new double[dmax];
    double bg[DMAX] = new double[dmax];
    double ng[DMAX] = new double[dmax];
    int ndprime[DMAX] = new double[dmax];

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
    bterPhases.phaseOneSeq(phase_one_i, phase_one_j);

    int *phase_two_i = new int[bterPhases.bterSamples.s2];
    int *phase_two_j = new int[bterPhases.bterSamples.s2];
    bterPhases.phaseTwoSeq(phase_two_i, phase_two_j);

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

