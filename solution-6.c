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
 * Number of buckets
 */
 const int NumberOfBuckets = 10;

/*
 * buckets[i] is the bucket particle i belongs to
 */
int* buckets;

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
    buckets = new int[NumberOfBodies]();

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
void updateBody() {
    auto *aX = new double[NumberOfBodies]();
    auto *aY = new double[NumberOfBodies]();
    auto *aZ = new double[NumberOfBodies]();

    const double vBucket = 130.0 / NumberOfBuckets;

    #pragma omp parallel for
    for (auto ii = 0; ii < NumberOfBodies; ++ii) {
        buckets[ii] = 0;
        const auto vi = v[ii][0]*v[ii][0] + v[ii][1]*v[ii][1] + v[ii][2]*v[ii][2];

        for (int j = NumberOfBuckets - 1; j >= 1; --j) {
            if (vi >= j * vBucket * j * vBucket) {
                buckets[ii] = j;
                break;
            }
        }
    }

    minDx = std::numeric_limits<double>::max();
    maxV = 0;

    #pragma omp parallel for reduction(max:maxV)
    for (auto bucketNum = 0; bucketNum < NumberOfBuckets; ++bucketNum) {
        const int steps = 1 << bucketNum;
        const auto dt = timeStepSize / steps;
        
        for (auto step = 0; step < steps; ++step) {
            // Compute the force felt on particles in bucket
            for (auto ii = 0; ii < NumberOfBodies; ++ii) {
                if (buckets[ii] != bucketNum) { continue; }

                const double m = mass[ii], px = x[ii][0], py = x[ii][1], pz = x[ii][2];        
                double Fx = 0, Fy = 0, Fz = 0;                                                       
                                                 
                for (int j = 0; j < NumberOfBodies; ++j) {                                           
                    if (ii == j) { continue; }                                                       
                                                         
                    const double dx = x[j][0] - px, dy = x[j][1] - py, dz = x[j][2] - pz;      
                    const double distSqrd = dx * dx + dy * dy + dz * dz, distance = sqrt(distSqrd);  
                                                             
                    const double k = m * mass[j] / (distSqrd * distance);                            
                                                             
                    Fx += k*dx;                                                                      
                    Fy += k*dy;                                                                      
                    Fz += k*dz;                                                                      
                                                             
                    minDx = std::min(minDx, distance);                                         
                }
                
                aX[ii] = Fx/m;
                aY[ii] = Fy/m;
                aZ[ii] = Fz/m;
            }

            // Update the positions & velocities
            for (auto ii = 0; ii < NumberOfBodies; ++ii) {
                if (buckets[ii] != bucketNum) { continue; }

                x[ii][0] += dt*v[ii][0];
                x[ii][1] += dt*v[ii][1];
                x[ii][2] += dt*v[ii][2];

                v[ii][0] += dt*aX[ii];
                v[ii][1] += dt*aY[ii];
                v[ii][2] += dt*aZ[ii];

                maxV = std::max(maxV, v[ii][0]*v[ii][0] + v[ii][1]*v[ii][1] + v[ii][2]*v[ii][2]);
            }

            if (minDx <= 0.01) {
                // Check for collisions
                for (auto ii = 0; ii < NumberOfBodies; ++ii) {
                    if (buckets[ii] != bucketNum) { continue; }

                    const auto px = x[ii][0], py = x[ii][1], pz = x[ii][2], m = mass[ii];

                    for (int j = ii+1; j < NumberOfBodies; ++j) {
                        if (ii == j) { continue; }

                        const auto dx = px-x[j][0], dy = py-x[j][1], dz = pz-x[j][2];
                        const auto distSqrd = dx*dx + dy*dy + dz*dz;

                        if (distSqrd > 0.01*0.01) {
                            continue;
                        }

                        // Particles ii and j have collided.
                        // We merge particles ii and j into the slot ii in x, v, mass
                        const auto M = m + mass[j];
                        const auto ki = m/M, kj = mass[j]/M;

                        v[ii][0] = v[ii][0]*ki + v[j][0]*kj;
                        v[ii][1] = v[ii][1]*ki + v[j][1]*kj;
                        v[ii][2] = v[ii][2]*ki + v[j][2]*kj;

                        x[ii][0] = x[ii][0]*ki + x[j][0]*kj;
                        x[ii][1] = x[ii][1]*ki + x[j][1]*kj;
                        x[ii][2] = x[ii][2]*ki + x[j][2]*kj;

                        mass[ii] = M;

                        mass[j] = mass[NumberOfBodies-1];
                        v[j] = v[NumberOfBodies-1];
                        x[j] = x[NumberOfBodies-1];
                        buckets[j] = buckets[NumberOfBodies-1];

                        --NumberOfBodies;
                    }
                }
            }
        }
    }

    delete[] aX;
    delete[] aY;
    delete[] aZ;

    maxV = std::sqrt(maxV);

    t += timeStepSize;

    if (NumberOfBodies == 1) {
        std::cout << x[0][0] << "," << x[0][1] << "," << x[0][2] << std::endl;
        t = tFinal+1;
        tPlot = t+1;
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
    while (t <= tFinal) {
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

    return 0;
}
