// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>
#include <cstring>
#include <vector>
#include <cmath>
#include <numeric>


double t = 0;
double tFinal = 0;
double tPlot = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double **x;

/**
 * Equivalent to x storing the velocities.
 */
double **v;

/**
 * One mass entry per molecule/particle.
 */
double *mass;

/**
 * Global time step size used.
 */
double timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double maxV;

/**
 * Minimum distance between two elements.
 */
double minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char **argv) {
    NumberOfBodies = (argc - 4) / 7;

    x = new double *[NumberOfBodies];
    v = new double *[NumberOfBodies];
    mass = new double[NumberOfBodies];

    int readArgument = 1;

    tPlotDelta = std::stof(argv[readArgument]);
    readArgument++;
    tFinal = std::stof(argv[readArgument]);
    readArgument++;
    timeStepSize = std::stof(argv[readArgument]);
    readArgument++;

    for (int i = 0; i < NumberOfBodies; i++) {
        x[i] = new double[3];
        v[i] = new double[3];

        x[i][0] = std::stof(argv[readArgument]);
        readArgument++;
        x[i][1] = std::stof(argv[readArgument]);
        readArgument++;
        x[i][2] = std::stof(argv[readArgument]);
        readArgument++;

        v[i][0] = std::stof(argv[readArgument]);
        readArgument++;
        v[i][1] = std::stof(argv[readArgument]);
        readArgument++;
        v[i][2] = std::stof(argv[readArgument]);
        readArgument++;

        mass[i] = std::stof(argv[readArgument]);
        readArgument++;

        if (mass[i] <= 0.0) {
            std::cerr << "invalid mass for body " << i << std::endl;
            exit(-2);
        }
    }

    std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

    if (tPlotDelta <= 0.0) {
        std::cout << "plotting switched off" << std::endl;
        tPlot = tFinal + 1.0;
    } else {
        std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
        tPlot = 0.0;
    }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
    videoFile.open("result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"
              << std::endl
              << "<Collection>";
}


/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
    videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
    static int counter = -1;
    counter++;
    std::stringstream filename;
    filename << "result-" << counter << ".vtp";
    std::ofstream out(filename.str().c_str());
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

    for (int i = 0; i < NumberOfBodies; i++) {
        out << x[i][0]
            << " "
            << x[i][1]
            << " "
            << x[i][2]
            << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>" << std::endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>"
              << std::endl;
}


/**
 * This is the only operation you are allowed to change in the assignment.
 */

#define MAKE_BUFFER(name) static auto* name = new double[NumberOfBodies]();
#define ZERO_BUFFER(name) std::memset(name, 0, sizeof(double)*NumberOfBodies);

void updateBody() {
    MAKE_BUFFER(forceX); MAKE_BUFFER(forceY); MAKE_BUFFER(forceZ);

    static std::vector<std::vector<int>> buckets(10);

    if (t == 0) {
        // If first timestep initalise the buckets.
        for (int ii = 0; ii < 10; ++ii) {
            buckets[ii] = std::vector<int>();
            buckets[ii].reserve(NumberOfBodies);
        }

        // On first iteration first bucket contains all the particles.
        for (int ii = 0; ii < NumberOfBodies; ++ii) {
            buckets[0].push_back(ii);
        }
    }

    maxV = 0;
    minDx = std::numeric_limits<double>::max();

    for (auto b = 0; b < buckets.size(); ++b) {
        auto& bucket = buckets[b];

        if (bucket.empty()) { continue; }

        auto bckNum = b+1;
        auto dt = timeStepSize/bckNum;

        for (auto step = 0; step < bckNum; ++step) {
            // Compute force felt on particles in the bucket.
            for (int p : bucket) {
                for (int ii = 0; ii < NumberOfBodies; ++ii) {
                    if (ii == p) { continue; }

                    const double dx = x[p][0] - x[ii][0], dy = x[p][1] - x[ii][1], dz = x[p][2] - x[ii][2];
                    const double distSqrd = dx * dx + dy * dy + dz * dz, distance = sqrt(distSqrd);

                    const double k = mass[ii] * mass[p] / (distSqrd * distance);
                    const double Fx = k * dx, Fy = k * dy, Fz = k * dz;

                    minDx = std::min(minDx, distance);

                    forceX[ii] += Fx;
                    forceY[ii] += Fy;
                    forceZ[ii] += Fz;
                }
            }

            // Update Position for each particle
            for (int p : bucket) {
                const auto aX = forceX[p]/mass[p], aY = forceY[p]/mass[p], aZ = forceZ[p]/mass[p];

                v[p][0] += dt*aX;
                v[p][1] += dt*aY;
                v[p][2] += dt*aZ;

                x[p][0] += dt*v[p][0];
                x[p][1] += dt*v[p][1];
                x[p][2] += dt*v[p][2];
            }
        }

        bucket.clear();
    }

    ZERO_BUFFER(forceX); ZERO_BUFFER(forceY); ZERO_BUFFER(forceZ);
    t += timeStepSize;

    // --------- See if any particles have collided.
    const double collisionDistanceThreshold = 0.01 * 0.01;
    const double vBucket = 10;

    // Check if an particles are inside of each other.
    for (auto ii = 0; ii < NumberOfBodies; ++ii) {
        for (auto j = ii + 1; j < NumberOfBodies; ++j) {
            const double dx = x[j][0] - x[ii][0], dy = x[j][1] - x[ii][1], dz = x[j][2] - x[ii][2];
            const double distSqrd = dx * dx + dy * dy + dz * dz;

            if (distSqrd > collisionDistanceThreshold) {
                continue;
            }

            // Particles ii and j have collided.
            // We merge particles ii and j into the slot ii in x, v, mass
            const double M = mass[ii] + mass[j];

            for (int k = 0; k < 3; ++k) {
                v[ii][k] = (mass[ii] * v[ii][k] + mass[j] * v[j][k]) / M;
                x[ii][k] = 0.5 * (x[ii][k] + x[j][k]);
            }

            mass[ii] = M;

            if (j != NumberOfBodies - 1) {
                // We then swap the now dead information at j with the info at NumberOfBodies-1 in x, v, mass
                std::swap(mass[j], mass[NumberOfBodies - 1]);
                std::swap(v[j], v[NumberOfBodies - 1]);
                std::swap(x[j], x[NumberOfBodies - 1]);
            }

            NumberOfBodies -= 1;
        }

        auto velocity = sqrt(v[ii][0] * v[ii][0] + v[ii][1] * v[ii][1] + v[ii][2] * v[ii][2]);
        maxV = std::max(maxV, velocity);
        buckets[std::min((int)std::floor(velocity/vBucket), 9)].push_back(ii);
    }

    if (NumberOfBodies == 1) {
        std::cout << x[0][0] << "," << x[0][1] << "," << x[0][2] << std::endl;
    }

}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char **argv) {
    if (argc == 1) {
        std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
                  << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting"
                  << std::endl
                  << "  final-time      simulated time (greater 0)" << std::endl
                  << "  dt              time step size (greater 0)" << std::endl
                  << std::endl
                  << "Examples:" << std::endl
                  << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1"
                  << std::endl
                  << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one"
                  << std::endl
                  << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture"
                  << std::endl
                  << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup"
                  << std::endl
                  << std::endl
                  << "In this naive code, only the first body moves" << std::endl;

        return -1;
    } else if ((argc - 4) % 7 != 0) {
        std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)"
                  << std::endl;
        std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
        std::cerr << "run without arguments for usage instruction" << std::endl;
        return -2;
    }

    std::cout << std::setprecision(15);

    setUp(argc, argv);

    openParaviewVideoFile();

    int snapshotCounter = 0;
    if (t > tPlot) {
        printParaviewSnapshot();
        std::cout << "plotted initial setup" << std::endl;
        tPlot = tPlotDelta;
    }

    int timeStepCounter = 0;
    while (t <= tFinal && NumberOfBodies > 1) {
        updateBody();
        timeStepCounter++;
        if (t >= tPlot) {
            printParaviewSnapshot();
            std::cout << "plot next snapshot"
                      << ",\t time step=" << timeStepCounter
                      << ",\t t=" << t
                      << ",\t dt=" << timeStepSize
                      << ",\t v_max=" << maxV
                      << ",\t dx_min=" << minDx
                      << std::endl;

            tPlot += tPlotDelta;
        }
    }

    printParaviewSnapshot();

    closeParaviewVideoFile();

    return 0;
}
