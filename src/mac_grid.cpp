#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#include <cstdlib>

#undef max
#undef min

#include <fstream>

// THREADS
#include <thread>

#define MULTITHREADING 0 // USE ONLY IF theDim[MACGrid::Z]>=4

// Globals
MACGrid target;


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;//true

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++)

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--)

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++)


MACGrid::MACGrid() {
    initialize();
}

MACGrid::MACGrid(const MACGrid &orig) {
    mU = orig.mU;
    mV = orig.mV;
    mW = orig.mW;
    mP = orig.mP;
    mD = orig.mD;
    mT = orig.mT;
}

MACGrid &MACGrid::operator=(const MACGrid &orig) {
    if (&orig == this) {
        return *this;
    }
    mU = orig.mU;
    mV = orig.mV;
    mW = orig.mW;
    mP = orig.mP;
    mD = orig.mD;
    mT = orig.mT;

    return *this;
}

MACGrid::~MACGrid() {
}

void MACGrid::reset() {
    mU.initialize();
    mV.initialize();
    mW.initialize();
    mP.initialize();
    mD.initialize();
    mT.initialize(0.0);

    solidCells.initialize(0.0);

    calculateAMatrix();
    calculatePreconditioner(AMatrix);
}

void MACGrid::initialize() {
    reset();
}

void MACGrid::updateSources() {
    // Set initial values for density, temperature, velocity
    int mul = 1;
    int minx = 0*mul, miny = 0*mul, minz = 0*mul;
    int maxx = 3*mul, maxy = 3*mul, maxz = 0*mul;

    for (int i = minx; i <= maxx; i++) {
        for (int j = miny; j <= maxy; j++) {
            for (int k = minz; k <= maxz; k++) {
                mV(i, j, k) = 3.0;
                //mU(i, j, k) = 2.0;
                mD(i, j, k) = 0.8;
                mT(i, j, k) = 3.0;
            }
        }
    }
    currentStep++;

    initializeSolids();

    // Refresh particles in source.
    for (int i = minx; i <= maxx; i++) {
        for (int j = miny; j <= maxy; j++) {
            for (int k = minz; k <= maxz; k++) {
                vec3 cell_center(theCellSize * (i + 0.5), theCellSize * (j + 0.5), theCellSize * (k + 0.5));
                for (int p = 0; p < 10; p++) {
                    double a = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    double b = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    double c = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    vec3 shift(a, b, c);
                    vec3 xp = cell_center + shift;
                    rendering_particles.push_back(xp);
                    rendering_particles_colIdx.push_back(0);
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if 0 // TURN ON OR OFF FOR MULTIPLE SOURCES
    // Set initial values for density, temperature, velocity
    minx = 29*mul, miny = 0*mul, minz = 0*mul;
    maxx = 34*mul, maxy = 5*mul, maxz = 0*mul;

    for (int i = minx; i <= maxx; i++) {
        for (int j = miny; j <= maxy; j++) {
            for (int k = minz; k <= maxz; k++) {
                mV(i, j, k) = 2.0;
                mU(i, j, k) = 2.0;
                mD(i, j, k) = 0.9;
                mT(i, j, k) = 2.0;
            }
        }
    }

    // Refresh particles in source.
    for (int i = minx; i <= maxx; i++) {
        for (int j = miny; j <= maxy; j++) {
            for (int k = minz; k <= maxz; k++) {
                vec3 cell_center(theCellSize * (i + 0.5), theCellSize * (j + 0.5), theCellSize * (k + 0.5));
                for (int p = 0; p < 10; p++) {
                    double a = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    double b = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    double c = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    vec3 shift(a, b, c);
                    vec3 xp = cell_center + shift;
                    rendering_particles.push_back(xp);
                    rendering_particles_colIdx.push_back(1);
                }
            }
        }
    }
#endif

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if 0
    // Set initial values for density, temperature, velocity
    minx = 6*mul, miny = 0*mul, minz = 6*mul;
    maxx = 7*mul, maxy = 1*mul, maxz = 7*mul;

    for (int i = minx; i <= maxx; i++) {
        for (int j = miny; j <= maxy; j++) {
            for (int k = minz; k <= maxz; k++) {
                mU(i, j, k) = -2.0;
                mV(i, j, k) = 2.0;
                mW(i, j, k) = -1.0;
                mD(i, j, k) = 1.0;
                mT(i, j, k) = 1.0;
            }
        }
    }

    // Refresh particles in source.
    for (int i = minx; i <= maxx; i++) {
        for (int j = miny; j <= maxy; j++) {
            for (int k = minz; k <= maxz; k++) {
                vec3 cell_center(theCellSize * (i + 0.5), theCellSize * (j + 0.5), theCellSize * (k + 0.5));
                for (int p = 0; p < 10; p++) {
                    double a = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    double b = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    double c = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    vec3 shift(a, b, c);
                    vec3 xp = cell_center + shift;
                    rendering_particles.push_back(xp);
                    rendering_particles_colIdx.push_back(2);
                }
            }
        }
    }
#endif

}

void MACGrid::initializeSolids() {

// CUBE
#if 0 // TURN ON OR OFF FOR OBSTACLES
    {
        float mul = 1;
        int minx = 0*mul, miny = 8*mul, minz = 0*mul;
        int maxx = 3*mul, maxy = 10*mul, maxz = 0*mul;
        for (int i = minx; i <= maxx; i++) {
            for (int j = miny; j <= maxy; j++) {
                for (int k = minz; k <= maxz; k++) {
                    solidCells(i, j, k) = 1.0;
                }
            }
        }
    }
#endif

// CUBE
#if 0
    {
        float mul = 2;
        int minx = 13*mul, miny = 12*mul, minz = 13*mul;
        int maxx = 18*mul, maxy = 15*mul, maxz = 18*mul;
        for (int i = minx; i <= maxx; i+=4) {
            for (int k = minz; k <= maxz; k+=4) {
                for (int j = miny; j <= maxy; j++) {
                    solidCells(i, j, k) = 1.0;
                    solidCells(i+1, j, k) = 1.0;
                    solidCells(i, j, k+1) = 1.0;
                    solidCells(i+1, j, k+1) = 1.0;
                }
            }
        }
    }
#endif
}

void MACGrid::advectVelocityThreadX(int tid, double dt) {
    //std::cout << "Thread: " << tid << std::endl;

    FOR_EACH_FACE {
        // X
        if (isValidFace(MACGrid::X, i, j, k)) {
            vec3 currentPt = getFacePosition(MACGrid::X, i, j, k);
            vec3 oldPt = getRewoundPosition(currentPt, dt);
            vec3 newVel = getVelocity(oldPt);
            target.mU(i, j, k) = newVel[0];
        }
    }
    mU = target.mU;
}

void MACGrid::advectVelocityThreadY(int tid, double dt) {
    //std::cout << "Thread: " << tid << std::endl;

    FOR_EACH_FACE {
        // Y
        if (isValidFace(MACGrid::Y, i, j, k)) {
            vec3 currentPt = getFacePosition(MACGrid::Y, i, j, k);
            vec3 oldPt = getRewoundPosition(currentPt, dt);
            vec3 newVel = getVelocity(oldPt);
            target.mV(i, j, k) = newVel[1];
        }
    }
    mV = target.mV;
}

void MACGrid::advectVelocityThreadZ(int tid, double dt) {
    //std::cout << "Thread: " << tid << std::endl;

    FOR_EACH_FACE {
        // Z
        if (isValidFace(MACGrid::Z, i, j, k)) {
            vec3 currentPt = getFacePosition(MACGrid::Z, i, j, k);
            vec3 oldPt = getRewoundPosition(currentPt, dt);
            vec3 newVel = getVelocity(oldPt);
            target.mW(i, j, k) = newVel[2];
        }
    }
    mW = target.mW;
}

void MACGrid::advectVelocity(double dt) {
#if MULTITHREADING
    std::thread t[3];
    t[0] = std::thread(&MACGrid::advectVelocityThreadX, this, 0, dt);
    t[1] = std::thread(&MACGrid::advectVelocityThreadY, this, 1, dt);
    t[2] = std::thread(&MACGrid::advectVelocityThreadZ, this, 2, dt);

    for (int i = 0; i < 3; i++) {
        t[i].join();
    }
#else
    FOR_EACH_FACE {
        // X
        if (isValidFace(MACGrid::X, i, j, k)) {
            vec3 currentPt = getFacePosition(MACGrid::X, i, j, k);
            vec3 oldPt = getRewoundPosition(currentPt, dt);
            vec3 newVel = getVelocity(oldPt);
            target.mU(i, j, k) = newVel[0];
        }
        // Y
        if (isValidFace(MACGrid::Y, i, j, k)) {
            vec3 currentPt = getFacePosition(MACGrid::Y, i, j, k);
            vec3 oldPt = getRewoundPosition(currentPt, dt);
            vec3 newVel = getVelocity(oldPt);
            target.mV(i, j, k) = newVel[1];
        }
        // Z
        if (isValidFace(MACGrid::Z, i, j, k)) {
            vec3 currentPt = getFacePosition(MACGrid::Z, i, j, k);
            vec3 oldPt = getRewoundPosition(currentPt, dt);
            vec3 newVel = getVelocity(oldPt);
            target.mW(i, j, k) = newVel[2];
        }
    }
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
#endif
}



void MACGrid::advectTemperatureThread(int tid, double dt, int kmin, int kmax) {
    for(int k = kmin; k < kmax; k++)
    for(int j = 0; j < theDim[MACGrid::Y]; j++)
    for(int i = 0; i < theDim[MACGrid::X]; i++) {
        vec3 currentPt = getCenter(i, j, k);
        vec3 oldPt = getRewoundPosition(currentPt, dt);
        double newTemp = getTemperature(oldPt);
        target.mT(i, j, k) = newTemp;
    }
}

void MACGrid::advectTemperature(double dt) {
#if MULTITHREADING
    std::thread t[4];
    int block = theDim[MACGrid::Z] / 4;
    for (int i = 0; i < 4; i++) {
        t[i] = std::thread(&MACGrid::advectTemperatureThread, this, i, dt, i * block , (i+1) * block);
    }

    for (int i = 0; i < 4; i++) {
        t[i].join();
    }
#else
    FOR_EACH_CELL {
        vec3 currentPt = getCenter(i, j, k);
        vec3 oldPt = getRewoundPosition(currentPt, dt);
        double newTemp = getTemperature(oldPt);
        target.mT(i, j, k) = newTemp;
    }
#endif
    mT = target.mT;
}

void MACGrid::advectRenderingParticlesThread(int tid, double dt, int kmin, int kmax) {
    for (size_t p = kmin; p < kmax; p++) {
        vec3 currentPosition = rendering_particles[p];
        vec3 currentVelocity = getVelocity(currentPosition);
        vec3 nextPosition = currentPosition + currentVelocity * dt;
        vec3 clippedNextPosition = clipToGrid(nextPosition, currentPosition);
        // Keep going...
        vec3 nextVelocity = getVelocity(clippedNextPosition);
        vec3 averageVelocity = (currentVelocity + nextVelocity) / 2.0;
        vec3 betterNextPosition = currentPosition + averageVelocity * dt;
        vec3 clippedBetterNextPosition = clipToGrid(betterNextPosition, currentPosition);
        rendering_particles[p] = clippedBetterNextPosition;
        rendering_particles_vel[p] = averageVelocity;
    }
}

void MACGrid::advectRenderingParticles(double dt) {
    rendering_particles_vel.resize(rendering_particles.size());
#if MULTITHREADING
    std::thread t[4];
    int block = rendering_particles.size() / 4;
    for (int i = 0; i < 4; i++) {
        t[i] = std::thread(&MACGrid::advectRenderingParticlesThread, this, i, dt, i * block, (i+1) * block);
    }

    for (int i = 0; i < 4; i++) {
        t[i].join();
    }
#else
    for (size_t p = 0; p < rendering_particles.size(); p++) {
        vec3 currentPosition = rendering_particles[p];
        vec3 currentVelocity = getVelocity(currentPosition);
        vec3 nextPosition = currentPosition + currentVelocity * dt;
        vec3 clippedNextPosition = clipToGrid(nextPosition, currentPosition);
        // Keep going...
        vec3 nextVelocity = getVelocity(clippedNextPosition);
        vec3 averageVelocity = (currentVelocity + nextVelocity) / 2.0;
        vec3 betterNextPosition = currentPosition + averageVelocity * dt;
        vec3 clippedBetterNextPosition = clipToGrid(betterNextPosition, currentPosition);
        rendering_particles[p] = clippedBetterNextPosition;
        rendering_particles_vel[p] = averageVelocity;
    }
#endif
}

void MACGrid::advectDensityThread(int tid, double dt, int kmin, int kmax) {
    for(int k = kmin; k < kmax; k++)
    for(int j = 0; j < theDim[MACGrid::Y]; j++)
    for(int i = 0; i < theDim[MACGrid::X]; i++) {
        vec3 currentPt = getCenter(i, j, k);
        vec3 oldPt = getRewoundPosition(currentPt, dt);
        double newDensity = getDensity(oldPt);
        target.mD(i, j, k) = newDensity;
    }
}

void MACGrid::advectDensity(double dt) {
#if MULTITHREADING
    std::thread t[4];
    int block = theDim[MACGrid::Z] / 4;
    for (int i = 0; i < 4; i++) {
        t[i] = std::thread(&MACGrid::advectDensityThread, this, i, dt, i * block, (i+1) * block);
    }

    for (int i = 0; i < 4; i++) {
        t[i].join();
    }
#else
    FOR_EACH_CELL {
        vec3 currentPt = getCenter(i, j, k);
        vec3 oldPt = getRewoundPosition(currentPt, dt);
        double newDensity = getDensity(oldPt);
        target.mD(i, j, k) = newDensity;
    }
#endif
    mD = target.mD;
}

void MACGrid::computeBouyancyThread(int tid, double dt, int kmin, int kmax) {
    for(int k = kmin; k < kmax; k++)
    for(int j = 0; j < theDim[MACGrid::Y]+1; j++)
    for(int i = 0; i < theDim[MACGrid::X]+1; i++) {
        if(isValidFace(MACGrid::Y, i, j, k)) {
            vec3 facePos = getFacePosition(MACGrid::Y, i, j, k);
            double buoyancy = - theBuoyancyAlpha * getDensity(facePos)
                           + theBuoyancyBeta * (getTemperature(facePos) - theBuoyancyAmbientTemperature);
            target.mV(i, j, k) += dt * buoyancy;
        }
    }
}

void MACGrid::computeBuoyancy(double dt) {
    target.mV = mV;
#if MULTITHREADING
    std::thread t[4];
    int block = theDim[MACGrid::Z] / 4;
    for (int i = 0; i < 3; i++) {
        t[i] = std::thread(&MACGrid::computeBouyancyThread, this, i, dt,
                           i * block, (i+1) * block);
    }
    t[3] = std::thread(&MACGrid::computeBouyancyThread, this, 3, dt,
                       3 * block, (3+1) * block + 1);

    for (int i = 0; i < 4; i++) {
        t[i].join();
    }
#else
    FOR_EACH_FACE{
        if(isValidFace(MACGrid::Y, i, j, k)) {
            vec3 facePos = getFacePosition(MACGrid::Y, i, j, k);
            double buoyancy = - theBuoyancyAlpha * getDensity(facePos)
                              + theBuoyancyBeta * (getTemperature(facePos) - theBuoyancyAmbientTemperature);
            target.mV(i, j, k) += dt * buoyancy;
        }
    }
#endif
    mV = target.mV;
}

void MACGrid::computeVorticityConfinementThread_CellCenterVel(int tid, double dt, int kmin, int kmax,
                                                              GridData *u, GridData *v, GridData *w) {
    for(int k = kmin; k < kmax; k++)
    for(int j = 0; j < theDim[MACGrid::Y]; j++)
    for(int i = 0; i < theDim[MACGrid::X]; i++) {
        (*u)(i, j, k) = (mU(i + 1, j, k) + mU(i, j, k)) / 2;
        (*v)(i, j, k) = (mV(i, j + 1, k) + mV(i, j, k)) / 2;
        (*w)(i, j, k) = (mW(i, j, k + 1) + mW(i, j, k)) / 2;
    }
}

void MACGrid::computeVorticityConfinementThread_CellCenterVort(int tid, double dt, int kmin, int kmax,
                                                               GridData *u, GridData *v, GridData *w,
                                                               GridData *x, GridData *y, GridData *z, GridData *len) {
    double inv2h = 1.0 / (2.0 * theCellSize);

    for(int k = kmin; k < kmax; k++)
    for(int j = 0; j < theDim[MACGrid::Y]; j++)
    for(int i = 0; i < theDim[MACGrid::X]; i++) {
        (*x)(i, j, k) = inv2h * ((*w)(i, j + 1, k) - (*w)(i, j - 1, k) - (*v)(i, j, k + 1) + (*v)(i, j, k - 1));
        (*y)(i, j, k) = inv2h * ((*u)(i, j, k + 1) - (*u)(i, j, k - 1) - (*w)(i + 1, j, k) + (*w)(i - 1, j, k));
        (*z)(i, j, k) = inv2h * ((*v)(i + 1, j, k) - (*v)(i - 1, j, k) - (*u)(i, j + 1, k) + (*u)(i, j - 1, k));
        (*len)(i, j, k) = vec3((*x)(i, j, k), (*y)(i, j, k), (*z)(i, j, k)).Length();
    }
}

void MACGrid::computeVorticityConfinementThread_Update(int tid, double dt, int kmin, int kmax,
                                              GridData *x, GridData *y, GridData *z, GridData *len) {
    double inv2h = 1.0 / (2.0 * theCellSize);

    for(int k = kmin; k < kmax; k++)
    for(int j = 0; j < theDim[MACGrid::Y]; j++)
    for(int i = 0; i < theDim[MACGrid::X]; i++) {
        vec3 vortGrad((*len)(i+1, j, k) - (*len)(i-1, j, k),
                      (*len)(i, j+1, k) - (*len)(i, j-1, k),
                      (*len)(i, j, k+1) - (*len)(i, j, k-1));
        vortGrad *= inv2h;

        vec3 N = vortGrad.Normalize();

        vec3 fconf = theVorticityEpsilon * theCellSize *
                     N.Cross(vec3((*x)(i, j, k), (*y)(i, j, k), (*z)(i, j, k)));

        // Apply vortcity confinement to faces..
        // X
        if(isValidFace(MACGrid::X, i, j, k)) {
            target.mU(i, j, k) += dt * fconf[0]/2;
        }
        // Y
        if(isValidFace(MACGrid::Y, i, j, k)) {
            target.mV(i, j, k) += dt * fconf[1]/2;
        }
        // Z
        if(isValidFace(MACGrid::Z, i, j, k)) {
            target.mW(i, j, k) += dt * fconf[2]/2;
        }
        // X+1
        if(isValidFace(MACGrid::X, i+1, j, k)) {
            target.mU(i+1, j, k) += dt * fconf[0]/2;
        }
        // Y+1
        if(isValidFace(MACGrid::Y, i, j+1, k)) {
            target.mV(i, j+1, k) += dt * fconf[1]/2;
        }
        // Z+1
        if(isValidFace(MACGrid::Z, i, j, k+1)) {
            target.mW(i, j, k+1) += dt * fconf[2]/2;
        }
    }
}

void MACGrid::computeVorticityConfinement(double dt) {
    GridData u; u.initialize(); // vel at cell center
    GridData v; v.initialize();
    GridData w; w.initialize();

    GridData vorticityX; vorticityX.initialize(); // vorticity at cell center
    GridData vorticityY; vorticityY.initialize();
    GridData vorticityZ; vorticityZ.initialize();
    GridData vorticityLen; vorticityLen.initialize();

    // Compute the vel at cell centers
#if MULTITHREADING
    std::thread t[4];
    int block = theDim[MACGrid::Z] / 4;
    for (int i = 0; i < 4; i++) {
        t[i] = std::thread(&MACGrid::computeVorticityConfinementThread_CellCenterVel,
                           this, i, dt, i * block, (i+1) * block, &u, &v, &w);
    }
    for (int i = 0; i < 4; i++) {
        t[i].join();
    }
#else
    FOR_EACH_CELL {
        u(i, j, k) = (mU(i + 1, j, k) + mU(i, j, k)) / 2;
        v(i, j, k) = (mV(i, j + 1, k) + mV(i, j, k)) / 2;
        w(i, j, k) = (mW(i, j, k + 1) + mW(i, j, k)) / 2;
    }
#endif

    // Compute the vorticity at cell centers.. eq 5.6
#if MULTITHREADING
    for (int i = 0; i < 4; i++) {
        t[i] = std::thread(&MACGrid::computeVorticityConfinementThread_CellCenterVort,
                           this, i, dt, i * block, (i+1) * block, &u, &v, &w,
                           &vorticityX, &vorticityY, &vorticityZ, &vorticityLen);
    }
    for (int i = 0; i < 4; i++) {
        t[i].join();
    }
#else
    double inv2h = 1.0 / (2.0 * theCellSize);
    FOR_EACH_CELL {
        vorticityX(i, j, k) = inv2h * (w(i, j + 1, k) - w(i, j - 1, k) - v(i, j, k + 1) + v(i, j, k - 1));
        vorticityY(i, j, k) = inv2h * (u(i, j, k + 1) - u(i, j, k - 1) - w(i + 1, j, k) + w(i - 1, j, k));
        vorticityZ(i, j, k) = inv2h * (v(i + 1, j, k) - v(i - 1, j, k) - u(i, j + 1, k) + u(i, j - 1, k));
        vorticityLen(i, j, k) = vec3(vorticityX(i, j, k), vorticityY(i, j, k), vorticityZ(i, j, k)).Length();
    }
#endif

    // Compute the vorticity at faces and add to vels
    target.mU = mU;
    target.mV = mV;
    target.mW = mW;

#if MULTITHREADING
    for (int i = 0; i < 4; i++) {
        t[i] = std::thread(&MACGrid::computeVorticityConfinementThread_Update,
                           this, i, dt, i * block, (i+1) * block,
                           &vorticityX, &vorticityY, &vorticityZ, &vorticityLen);
    }
    for (int i = 0; i < 4; i++) {
        t[i].join();
    }
#else
    //double inv2h = 1.0 / (2.0 * theCellSize);
    FOR_EACH_CELL {
        // eq 5.7
        vec3 vortGrad(vorticityLen(i+1, j, k) - vorticityLen(i-1, j, k),
                      vorticityLen(i, j+1, k) - vorticityLen(i, j-1, k),
                      vorticityLen(i, j, k+1) - vorticityLen(i, j, k-1));
        vortGrad *= inv2h;

        vec3 N = vortGrad.Normalize();

        vec3 fconf = theVorticityEpsilon * theCellSize *
                     N.Cross(vec3(vorticityX(i, j, k), vorticityY(i, j, k), vorticityZ(i, j, k)));

        // Apply vortcity confinement to faces..
        // X
        if(isValidFace(MACGrid::X, i, j, k)) {
            target.mU(i, j, k) += dt * fconf[0]/2;
        }
        // Y
        if(isValidFace(MACGrid::Y, i, j, k)) {
            target.mV(i, j, k) += dt * fconf[1]/2;
        }
        // Z
        if(isValidFace(MACGrid::Z, i, j, k)) {
            target.mW(i, j, k) += dt * fconf[2]/2;
        }
        // X+1
        if(isValidFace(MACGrid::X, i+1, j, k)) {
            target.mU(i+1, j, k) += dt * fconf[0]/2;
        }
        // Y+1
        if(isValidFace(MACGrid::Y, i, j+1, k)) {
            target.mV(i, j+1, k) += dt * fconf[1]/2;
        }
        // Z+1
        if(isValidFace(MACGrid::Z, i, j, k+1)) {
            target.mW(i, j, k+1) += dt * fconf[2]/2;
        }
    }
#endif

    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::addExternalForces(double dt) {
    computeBuoyancy(dt);
    computeVorticityConfinement(dt);
}


void MACGrid::projectThread_Divergence(int tid, double dt, int kmin, int kmax,
                              GridData *d) {
    for(int k = kmin; k < kmax; k++)
    for(int j = 0; j < theDim[MACGrid::Y]; j++)
    for(int i = 0; i < theDim[MACGrid::X]; i++) {
        double velLowX = (i > 0) ? mU(i, j, k) : 0.0;
        double velHighX = (i + 1 < theDim[MACGrid::X]) ? mU(i + 1, j, k) : 0.0;
        double velLowY = (j > 0) ? mV(i, j, k) : 0.0;
        double velHighY = (j + 1 < theDim[MACGrid::Y]) ? mV(i, j + 1, k) : 0.0;
        double velLowZ = (k > 0) ? mW(i, j, k) : 0.0;
        double velHighZ = (k + 1 < theDim[MACGrid::Z]) ? mW(i, j, k + 1) : 0.0;
        (*d)(i, j, k) = -((velHighX - velLowX) + (velHighY - velLowY) + (velHighZ - velLowZ)) / theCellSize;
    }
}

void MACGrid::projectThread_UpdateX(int tid, double dt, GridData *p){
    target.mU = mU;
    double constMultiplier1 = dt / (theAirDensity * theCellSize);
    FOR_EACH_FACE {
        // X
        if (isValidFace(MACGrid::X, i, j, k)) {
            if (isValidCell(i - 1, j, k) && isValidCell(i, j, k)) {
                double deltaPX = (*p)(i, j, k) - (*p)(i - 1, j, k);
                target.mU(i, j, k) -= constMultiplier1 * deltaPX;
            }
            else{
                target.mU(i, j, k) = 0.0;
            }
        }
    }
    mU = target.mU;
}

void MACGrid::projectThread_UpdateY(int tid, double dt, GridData *p){
    target.mV = mV;
    double constMultiplier1 = dt / (theAirDensity * theCellSize);
    FOR_EACH_FACE {
        // Y
        if(isValidFace(MACGrid::Y, i, j, k)) {
            if (isValidCell(i, j - 1, k) && isValidCell(i, j, k)) {
                double deltaPY = (*p)(i, j, k) - (*p)(i, j - 1, k);
                target.mV(i, j, k) -= constMultiplier1 * deltaPY;
            }
            else{
                target.mV(i, j, k) = 0.0;
            }
        }
    }

    mV = target.mV;
}

void MACGrid::projectThread_UpdateZ(int tid, double dt, GridData *p){
    target.mW = mW;
    double constMultiplier1 = dt / (theAirDensity * theCellSize);
    FOR_EACH_FACE {
        // Z
        if(isValidFace(MACGrid::Z, i, j, k)) {
            if (isValidCell(i, j, k - 1) && isValidCell(i, j, k)) {
                double deltaPZ = (*p)(i, j, k) - (*p)(i, j, k - 1);
                target.mW(i, j, k) -= constMultiplier1 * deltaPZ;
            }
            else{
                target.mW(i, j, k) = 0.0;
            }
        }
    }

    mW = target.mW;
}

void MACGrid::project(double dt) {
    // AP = D

    double constMultiplier = theAirDensity * theCellSize * theCellSize / dt; // sig notes: page 31..

    GridData p = GridData(); p.initialize();
    GridData d = GridData(); d.initialize();

    // Compute divergence d1
#if MULTITHREADING
    std::thread t[4];
    int block = theDim[MACGrid::Z] / 4;
    for (int i = 0; i < 4; i++) {
        t[i] = std::thread(&MACGrid::projectThread_Divergence,
                           this, i, dt, i * block, (i+1) * block, &d);
    }
    for (int i = 0; i < 4; i++) {
        t[i].join();
    }
#else
    FOR_EACH_CELL {
        double velLowX = (i > 0 && solidCells(i, j, k) != 1.0) ? mU(i, j, k) : 0.0;
        double velHighX = (i + 1 < theDim[MACGrid::X] && solidCells(i+1, j, k) != 1.0) ? mU(i + 1, j, k) : 0.0;
        double velLowY = (j > 0 && solidCells(i, j, k) != 1.0) ? mV(i, j, k) : 0.0;
        double velHighY = (j + 1 < theDim[MACGrid::Y] && solidCells(i, j+1, k) != 1.0) ? mV(i, j + 1, k) : 0.0;
        double velLowZ = (k > 0 && solidCells(i, j, k) != 1.0) ? mW(i, j, k) : 0.0;
        double velHighZ = (k + 1 < theDim[MACGrid::Z] && solidCells(i, j, k+1) != 1.0) ? mW(i, j, k + 1) : 0.0;
        d(i, j, k) = -((velHighX - velLowX) + (velHighY - velLowY) + (velHighZ - velLowZ)) / theCellSize;
    }
#endif

    // Compute pressure p
    preconditionedConjugateGradient(AMatrix, p, d, 200, 0.01);
    FOR_EACH_CELL {
        p(i, j, k) *= constMultiplier;
    }
    target.mP = mP;

    // Update velocities using p.. refer class notes for equations
#if MULTITHREADING
    t[0] = std::thread(&MACGrid::projectThread_UpdateX, this, 0, dt, &p);
    t[1] = std::thread(&MACGrid::projectThread_UpdateY, this, 1, dt, &p);
    t[2] = std::thread(&MACGrid::projectThread_UpdateZ, this, 2, dt, &p);

    for (int i = 0; i < 3; i++) {
        t[i].join();
    }

#else
    target.mU = mU;
    target.mV = mV;
    target.mW = mW;
    double constMultiplier1 = dt / (theAirDensity * theCellSize);
    FOR_EACH_FACE {
        // X
        if (isValidFace(MACGrid::X, i, j, k)) {
            if (isValidCell(i - 1, j, k) && isValidCell(i, j, k)) {
                double deltaPX = p(i, j, k) - p(i - 1, j, k);
                target.mU(i, j, k) -= constMultiplier1 * deltaPX;
            }
            else{
                target.mU(i, j, k) = 0.0;
            }
        }

        // Y
        if(isValidFace(MACGrid::Y, i, j, k)) {
            if (isValidCell(i, j - 1, k) && isValidCell(i, j, k)) {
                double deltaPY = p(i, j, k) - p(i, j - 1, k);
                target.mV(i, j, k) -= constMultiplier1 * deltaPY;
            }
            else{
                target.mV(i, j, k) = 0.0;
            }
        }

        // Z
        if(isValidFace(MACGrid::Z, i, j, k)) {
            if (isValidCell(i, j, k - 1) && isValidCell(i, j, k)) {
                double deltaPZ = p(i, j, k) - p(i, j, k - 1);
                target.mW(i, j, k) -= constMultiplier1 * deltaPZ;
            }
            else{
                target.mW(i, j, k) = 0.0;
            }
        }
    }

    mU = target.mU;
    mV = target.mV;
    mW = target.mW;

#endif

#ifdef _DEBUG
    // Check border velocities:
    FOR_EACH_FACE {
        if (isValidFace(MACGrid::X, i, j, k)) {

            if (i == 0) {
                if (abs(target.mU(i,j,k)) > 0.0000001) {
                    PRINT_LINE( "LOW X:  " << target.mU(i,j,k) );
                    //target.mU(i,j,k) = 0;
                }
            }

            if (i == theDim[MACGrid::X]) {
                if (abs(target.mU(i,j,k)) > 0.0000001) {
                    PRINT_LINE( "HIGH X: " << target.mU(i,j,k) );
                    //target.mU(i,j,k) = 0;
                }
            }

        }
        if (isValidFace(MACGrid::Y, i, j, k)) {


            if (j == 0) {
                if (abs(target.mV(i,j,k)) > 0.0000001) {
                    PRINT_LINE( "LOW Y:  " << target.mV(i,j,k) );
                    //target.mV(i,j,k) = 0;
                }
            }

            if (j == theDim[MACGrid::Y]) {
                if (abs(target.mV(i,j,k)) > 0.0000001) {
                    PRINT_LINE( "HIGH Y: " << target.mV(i,j,k) );
                    //target.mV(i,j,k) = 0;
                }
            }

        }
        if (isValidFace(MACGrid::Z, i, j, k)) {

            if (k == 0) {
                if (abs(target.mW(i,j,k)) > 0.0000001) {
                    PRINT_LINE( "LOW Z:  " << target.mW(i,j,k) );
                    //target.mW(i,j,k) = 0;
                }
            }

            if (k == theDim[MACGrid::Z]) {
                if (abs(target.mW(i,j,k)) > 0.0000001) {
                    PRINT_LINE( "HIGH Z: " << target.mW(i,j,k) );
                    //target.mW(i,j,k) = 0;
                }
            }
        }
    }
#endif

    // Then save the result to our object
    mP = target.mP;


#ifdef _DEBUG
    // IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
    // TODO: Fix duplicate code:
    FOR_EACH_CELL {
        // Construct the vector of divergences d:
         double velLowX = mU(i,j,k);
         double velHighX = mU(i+1,j,k);
         double velLowY = mV(i,j,k);
         double velHighY = mV(i,j+1,k);
         double velLowZ = mW(i,j,k);
         double velHighZ = mW(i,j,k+1);
         double divergence = ((velHighX - velLowX) + (velHighY - velLowY) + (velHighZ - velLowZ)) / theCellSize;
         if (abs(divergence) > 0.02 ) {
             PRINT_LINE("WARNING: Divergent! ");
             PRINT_LINE("Divergence: " << divergence);
             PRINT_LINE("Cell: " << i << ", " << j << ", " << k);
         }
    }
#endif

}

vec3 MACGrid::getVelocity(const vec3 &pt) {
    vec3 vel;
    vel[0] = getVelocityX(pt);
    vel[1] = getVelocityY(pt);
    vel[2] = getVelocityZ(pt);
    return vel;
}

double MACGrid::getVelocityX(const vec3 &pt) {
    return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3 &pt) {
    return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3 &pt) {
    return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3 &pt) {
    return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3 &pt) {
    return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k) {
    double xstart = theCellSize / 2.0;
    double ystart = theCellSize / 2.0;
    double zstart = theCellSize / 2.0;

    double x = xstart + i * theCellSize;
    double y = ystart + j * theCellSize;
    double z = zstart + k * theCellSize;
    return vec3(x, y, z);
}


vec3 MACGrid::getRewoundPosition(const vec3 &currentPosition, const double dt) {

    /*
    // EULER (RK1):
    vec3 currentVelocity = getVelocity(currentPosition);
    vec3 rewoundPosition = currentPosition - currentVelocity * dt;
    vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
    return clippedRewoundPosition;
    */

    // HEUN / MODIFIED EULER (RK2):
    vec3 currentVelocity = getVelocity(currentPosition);
    vec3 rewoundPosition = currentPosition - currentVelocity * dt;
    vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
    // Keep going...
    vec3 rewoundVelocity = getVelocity(clippedRewoundPosition);
    vec3 averageVelocity = (currentVelocity + rewoundVelocity) / 2.0;
    vec3 betterRewoundPosition = currentPosition - averageVelocity * dt;
    vec3 clippedBetterRewoundPosition = clipToGrid(betterRewoundPosition, currentPosition);
    return clippedBetterRewoundPosition;

}


vec3 MACGrid::clipToGrid(const vec3 &outsidePoint, const vec3 &insidePoint) {
    /*
    // OLD:
    vec3 rewindPosition = outsidePoint;
    if (rewindPosition[0] < 0) rewindPosition[0] = 0; // TEMP!
    if (rewindPosition[1] < 0) rewindPosition[1] = 0; // TEMP!
    if (rewindPosition[2] < 0) rewindPosition[2] = 0; // TEMP!
    if (rewindPosition[0] > theDim[MACGrid::X]) rewindPosition[0] = theDim[MACGrid::X]; // TEMP!
    if (rewindPosition[1] > theDim[MACGrid::Y]) rewindPosition[1] = theDim[MACGrid::Y]; // TEMP!
    if (rewindPosition[2] > theDim[MACGrid::Z]) rewindPosition[2] = theDim[MACGrid::Z]; // TEMP!
    return rewindPosition;
    */

    vec3 clippedPoint = outsidePoint;

    for (int i = 0; i < 3; i++) {
        if (clippedPoint[i] < 0) {
            vec3 distance = clippedPoint - insidePoint;
            double newDistanceI = 0 - insidePoint[i];
            double ratio = newDistanceI / distance[i];
            clippedPoint = insidePoint + distance * ratio;
        }
        if (clippedPoint[i] > getSize(i)) {
            vec3 distance = clippedPoint - insidePoint;
            double newDistanceI = getSize(i) - insidePoint[i];
            double ratio = newDistanceI / distance[i];
            clippedPoint = insidePoint + distance * ratio;
        }
    }

#ifdef _DEBUG
    // Make sure the point is now in the grid:
    if (clippedPoint[0] < 0 || clippedPoint[1] < 0 || clippedPoint[2] < 0 || clippedPoint[0] > getSize(0) || clippedPoint[1] > getSize(1) || clippedPoint[2] > getSize(2)) {
        PRINT_LINE("WARNING: Clipped point is outside grid!");
    }
#endif

    return clippedPoint;

}


double MACGrid::getSize(int dimension) {
    return theDim[dimension] * theCellSize;
}


int MACGrid::getCellIndex(int i, int j, int k) {
    return i + j * theDim[MACGrid::X] + k * theDim[MACGrid::Y] * theDim[MACGrid::X];
}


int MACGrid::getNumberOfCells() {
    return theDim[MACGrid::X] * theDim[MACGrid::Y] * theDim[MACGrid::Z];
}


bool MACGrid::isValidCell(int i, int j, int k) {
    if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
        return false;
    }

    if (i < 0 || j < 0 || k < 0) {
        return false;
    }

    if(solidCells(i, j, k) == 1.0) {
        return false;
    }

    return true;
}


bool MACGrid::isValidFace(int dimension, int i, int j, int k) {
    if (dimension == 0) {
        if (i > theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
            return false;
        }
        if(solidCells(i-1, j, k) == 1.0 && solidCells(i, j, k) == 1.0) {
            return false;
        }
    } else if (dimension == 1) {
        if (i >= theDim[MACGrid::X] || j > theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
            return false;
        }
        if(solidCells(i, j-1, k) == 1.0 && solidCells(i, j, k) == 1.0) {
            return false;
        }
    } else if (dimension == 2) {
        if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k > theDim[MACGrid::Z]) {
            return false;
        }
        if(solidCells(i, j, k-1) == 1.0 && solidCells(i, j, k) == 1.0) {
            return false;
        }
    }

    if (i < 0 || j < 0 || k < 0) {
        return false;
    }

    return true;
}


vec3 MACGrid::getFacePosition(int dimension, int i, int j, int k) {
    if (dimension == 0) {
        return vec3(i * theCellSize, (j + 0.5) * theCellSize, (k + 0.5) * theCellSize);
    } else if (dimension == 1) {
        return vec3((i + 0.5) * theCellSize, j * theCellSize, (k + 0.5) * theCellSize);
    } else if (dimension == 2) {
        return vec3((i + 0.5) * theCellSize, (j + 0.5) * theCellSize, k * theCellSize);
    }

    return vec3(0, 0, 0); //???

}

void MACGrid::calculateAMatrix() {

    FOR_EACH_CELL {

        int numFluidNeighbors = 0;
        if (i - 1 >= 0) {
            AMatrix.plusI(i - 1, j, k) = -1;
            numFluidNeighbors++;
        }
        if (i + 1 < theDim[MACGrid::X]) {
            AMatrix.plusI(i, j, k) = -1;
            numFluidNeighbors++;
        }
        if (j - 1 >= 0) {
            AMatrix.plusJ(i, j - 1, k) = -1;
            numFluidNeighbors++;
        }
        if (j + 1 < theDim[MACGrid::Y]) {
            AMatrix.plusJ(i, j, k) = -1;
            numFluidNeighbors++;
        }
        if (k - 1 >= 0) {
            AMatrix.plusK(i, j, k - 1) = -1;
            numFluidNeighbors++;
        }
        if (k + 1 < theDim[MACGrid::Z]) {
            AMatrix.plusK(i, j, k) = -1;
            numFluidNeighbors++;
        }
        // Set the diagonal:
        AMatrix.diag(i, j, k) = numFluidNeighbors;
    }
}

bool MACGrid::preconditionedConjugateGradient(const GridDataMatrix &A, GridData &p, const GridData &d,
                                         int maxIterations, double tolerance) {
    // Solves Ap = d for p.

    FOR_EACH_CELL {
        p(i, j, k) = 0.0; // Initial guess p = 0.
    }

    GridData r = d; // Residual vector.

    /*
    PRINT_LINE("r: ");
    FOR_EACH_CELL {
        PRINT_LINE(r(i,j,k));
    }
    */
    GridData z;
    z.initialize();
    applyPreconditioner(r, A, z); // Auxillary vector.
    /*
    PRINT_LINE("z: ");
    FOR_EACH_CELL {
        PRINT_LINE(z(i,j,k));
    }
    */

    GridData s = z; // Search vector;

    double sigma = dotProduct(z, r);

    for (int iteration = 0; iteration < maxIterations; iteration++) {

        double rho = sigma; // According to TA. Here???

        apply(A, s, z); // z = applyA(s);

        double alpha = rho / dotProduct(z, s);

        GridData alphaTimesS;
        alphaTimesS.initialize();
        multiply(alpha, s, alphaTimesS);
        add(p, alphaTimesS, p);
        //p += alpha * s;

        GridData alphaTimesZ;
        alphaTimesZ.initialize();
        multiply(alpha, z, alphaTimesZ);
        subtract(r, alphaTimesZ, r);
        //r -= alpha * z;

        if (maxMagnitude(r) <= tolerance) {
            //PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
            return true; //return p;
        }

        applyPreconditioner(r, A, z); // z = applyPreconditioner(r);

        double sigmaNew = dotProduct(z, r);

        double beta = sigmaNew / rho;

        GridData betaTimesS;
        betaTimesS.initialize();
        multiply(beta, s, betaTimesS);
        add(z, betaTimesS, s);
        //s = z + beta * s;

        sigma = sigmaNew;
    }

    PRINT_LINE("PCG didn't converge!");
    return false;

}

void MACGrid::calculatePreconditionerThread(int tid, GridDataMatrix *A, int kmin, int kmax) {
    double tuning = 0.97;
    for(int k = kmin; k < kmax; k++)
    for(int j = 0; j < theDim[MACGrid::Y]; j++)
    for(int i = 0; i < theDim[MACGrid::X]; i++) {
        double a = A->plusI(i - 1, j, k) * precon(i - 1, j, k);
        double b = A->plusJ(i, j - 1, k) * precon(i, j - 1, k);
        double c = A->plusK(i, j, k - 1) * precon(i, j, k - 1);
        double e = A->diag(i, j, k) - a * a - b * b - c * c
                   - tuning * (a * precon(i - 1, j, k) * (A->plusJ(i - 1, j, k) + A->plusK(i - 1, j, k))
                               + b * precon(i, j - 1, k) * (A->plusI(i, j - 1, k) + A->plusK(i, j - 1, k))
                               + c * precon(i, j, k - 1) * (A->plusI(i, j, k - 1) + A->plusJ(i, j, k - 1)));

        precon(i, j, k) = 1.0 / sqrt(e + EPSILON);
    }
}

void MACGrid::calculatePreconditioner(GridDataMatrix &A) {

    precon.initialize();

    // TODO: Build the modified incomplete Cholesky preconditioner following Fig 4.2 on page 36 of Bridson's 2007 SIGGRAPH fluid course notes.
    //       This corresponds to filling in precon(i,j,k) for all cells

#if MULTITHREADING
    std::thread t[4];
    int block = theDim[MACGrid::Z] / 4;
    for (int i = 0; i < 4; i++) {
        t[i] = std::thread(&MACGrid::calculatePreconditionerThread, this, i, &A, i * block, (i+1) * block);
    }

    for (int i = 0; i < 4; i++) {
        t[i].join();
    }
#else
    double tuning = 0.97;
    FOR_EACH_CELL {
        double a = A.plusI(i - 1, j, k) * precon(i - 1, j, k);
        double b = A.plusJ(i, j - 1, k) * precon(i, j - 1, k);
        double c = A.plusK(i, j, k - 1) * precon(i, j, k - 1);
        double e = A.diag(i, j, k) - a * a - b * b - c * c
                   - tuning * (a * precon(i - 1, j, k) * (A.plusJ(i - 1, j, k) + A.plusK(i - 1, j, k))
                               + b * precon(i, j - 1, k) * (A.plusI(i, j - 1, k) + A.plusK(i, j - 1, k))
                               + c * precon(i, j, k - 1) * (A.plusI(i, j, k - 1) + A.plusJ(i, j, k - 1)));

        precon(i, j, k) = 1.0 / sqrt(e + EPSILON);
    }
#endif
}


void MACGrid::applyPreconditioner(const GridData &r, const GridDataMatrix &A, GridData &z) {

    // TODO: change if(0) to if(1) after you implement calculatePreconditoner function.

    if (1) {

        // APPLY THE PRECONDITIONER:
        // Solve Lq = r for q:
        GridData q;
        q.initialize();
        FOR_EACH_CELL {
            //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
            double t = r(i, j, k) - A.plusI(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k)
                       - A.plusJ(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k)
                       - A.plusK(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
            q(i, j, k) = t * precon(i, j, k);
            //}
        }
        // Solve L^Tz = q for z:
        FOR_EACH_CELL_REVERSE {
            //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
            double t = q(i, j, k) - A.plusI(i, j, k) * precon(i, j, k) * z(i + 1, j, k)
                       - A.plusJ(i, j, k) * precon(i, j, k) * z(i, j + 1, k)
                       - A.plusK(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
            z(i, j, k) = t * precon(i, j, k);
            //}
        }
    } else {
        // Unpreconditioned CG: Bypass preconditioner:
        z = r;
        return;
    }

}


double MACGrid::dotProduct(const GridData &vector1, const GridData &vector2) {

    double result = 0.0;

    FOR_EACH_CELL {
                result += vector1(i, j, k) * vector2(i, j, k);
            }

    return result;
}


void MACGrid::add(const GridData &vector1, const GridData &vector2, GridData &result) {

    FOR_EACH_CELL {
                result(i, j, k) = vector1(i, j, k) + vector2(i, j, k);
            }

}


void MACGrid::subtract(const GridData &vector1, const GridData &vector2, GridData &result) {

    FOR_EACH_CELL {
                result(i, j, k) = vector1(i, j, k) - vector2(i, j, k);
            }

}


void MACGrid::multiply(const double scalar, const GridData &vector, GridData &result) {

    FOR_EACH_CELL {
                result(i, j, k) = scalar * vector(i, j, k);
            }

}


double MACGrid::maxMagnitude(const GridData &vector) {

    double result = 0.0;

    FOR_EACH_CELL {
                if (abs(vector(i, j, k)) > result) result = abs(vector(i, j, k));
            }

    return result;
}


void MACGrid::apply(const GridDataMatrix &matrix, const GridData &vector, GridData &result) {

    FOR_EACH_CELL { // For each row of the matrix.

        double diag = 0;
        double plusI = 0;
        double plusJ = 0;
        double plusK = 0;
        double minusI = 0;
        double minusJ = 0;
        double minusK = 0;

        diag = matrix.diag(i, j, k) * vector(i, j, k);
        if (isValidCell(i + 1, j, k)) plusI = matrix.plusI(i, j, k) * vector(i + 1, j, k);
        if (isValidCell(i, j + 1, k)) plusJ = matrix.plusJ(i, j, k) * vector(i, j + 1, k);
        if (isValidCell(i, j, k + 1)) plusK = matrix.plusK(i, j, k) * vector(i, j, k + 1);
        if (isValidCell(i - 1, j, k)) minusI = matrix.plusI(i - 1, j, k) * vector(i - 1, j, k);
        if (isValidCell(i, j - 1, k)) minusJ = matrix.plusJ(i, j - 1, k) * vector(i, j - 1, k);
        if (isValidCell(i, j, k - 1)) minusK = matrix.plusK(i, j, k - 1) * vector(i, j, k - 1);

        result(i, j, k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
    }
}

void MACGrid::saveSmoke(const char *fileName) {
    std::ofstream fileOut(fileName);
    if (fileOut.is_open()) {
        FOR_EACH_CELL {
                    fileOut << mD(i, j, k) << std::endl;
                }
        fileOut.close();
    }
}

void MACGrid::saveParticle(std::string filename) {
    Partio::ParticlesDataMutable *parts = Partio::create();
    Partio::ParticleAttribute posH, vH, cH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);
    cH = parts->addAttribute("c", Partio::VECTOR, 3);

    for (unsigned int i = 0; i < rendering_particles.size(); i++) {
        int idx = parts->addParticle();
        float *p = parts->dataWrite<float>(posH, idx);
        float *v = parts->dataWrite<float>(vH, idx);
        float *c = parts->dataWrite<float>(cH, idx);
        for (int k = 0; k < 3; k++) {
            p[k] = rendering_particles[i][k];
            v[k] = rendering_particles_vel[i][k];
        }

        if (rendering_particles_colIdx[i]==0) { // cream
            c[0] = 0.957; c[1] = 0.945; c[2] = 0.859;
        }
        else if (rendering_particles_colIdx[i]==1) { // orange
            c[0] = 1.0; c[1] = 0.341; c[2] = 0.161;
        }
        else if (rendering_particles_colIdx[i]==2) { // cyan??
            c[0] = 0.0; c[1] = 0.580; c[2] = 0.580;
        }
    }

    Partio::write(filename.c_str(), *parts);
    parts->release();
}

void MACGrid::saveDensity(std::string filename) {
    Partio::ParticlesDataMutable *density_field = Partio::create();
    Partio::ParticleAttribute posH, rhoH;
    posH = density_field->addAttribute("position", Partio::VECTOR, 3);
    rhoH = density_field->addAttribute("density", Partio::VECTOR, 1);
    FOR_EACH_CELL {
                int idx = density_field->addParticle();
                float *p = density_field->dataWrite<float>(posH, idx);
                float *rho = density_field->dataWrite<float>(rhoH, idx);
                vec3 cellCenter = getCenter(i, j, k);
                for (int l = 0; l < 3; l++) {
                    p[l] = cellCenter[l];
                }
                rho[0] = getDensity(cellCenter);
            }
    Partio::write(filename.c_str(), *density_field);
    density_field->release();
}

void MACGrid::draw(const Camera &c) {
    drawWireGrid();
    if (theDisplayVel) drawVelocities();
    if (theRenderMode == CUBES) drawSmokeCubes(c);
    else drawSmoke(c);
}

void MACGrid::drawVelocities() {
    // draw line at each center
    //glColor4f(0.0, 1.0, 0.0, 1.0);
    glBegin(GL_LINES);
    FOR_EACH_CELL {
                vec3 pos = getCenter(i, j, k);
                vec3 vel = getVelocity(pos);
                if (vel.Length() > 0.0001) {
                    //vel.Normalize();
                    vel *= theCellSize / 2.0;
                    vel += pos;
                    glColor4f(1.0, 1.0, 0.0, 1.0);
                    glVertex3dv(pos.n);
                    glColor4f(0.0, 1.0, 0.0, 1.0);
                    glVertex3dv(vel.n);
                }
            }
    glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k) {

    double value = mD(i, j, k);
    vec4 coldColor(0.5, 0.5, 1.0, value);
    vec4 hotColor(1.0, 0.5, 0.5, value);
    return LERP(coldColor, hotColor, mT(i, j, k));


    /*
    // OLD:
    double value = mD(i, j, k);
    return vec4(1.0, 0.9, 1.0, value);
    */
}

vec4 MACGrid::getRenderColor(const vec3 &pt) {
    double value = getDensity(pt);
    vec4 coldColor(0.5, 0.5, 1.0, value);
    vec4 hotColor(1.0, 0.5, 0.5, value);
    return LERP(coldColor, hotColor, getTemperature(pt));

    /*
    // OLD:
    double value = getDensity(pt);
    return vec4(1.0, 1.0, 1.0, value);
    */
}

void MACGrid::drawZSheets(bool backToFront) {
    // Draw K Sheets from back to front
    double back = (theDim[2]) * theCellSize;
    double top = (theDim[1]) * theCellSize;
    double right = (theDim[0]) * theCellSize;

    double stepsize = theCellSize * 0.25;

    double startk = back - stepsize;
    double endk = 0;
    double stepk = -theCellSize;

    if (!backToFront) {
        startk = 0;
        endk = back;
        stepk = theCellSize;
    }

    for (double k = startk; backToFront ? k > endk : k < endk; k += stepk) {
        for (double j = 0.0; j < top;) {
            glBegin(GL_QUAD_STRIP);
            for (double i = 0.0; i <= right; i += stepsize) {
                vec3 pos1 = vec3(i, j, k);
                vec3 pos2 = vec3(i, j + stepsize, k);

                vec4 color1 = getRenderColor(pos1);
                vec4 color2 = getRenderColor(pos2);

                glColor4dv(color1.n);
                glVertex3dv(pos1.n);

                glColor4dv(color2.n);
                glVertex3dv(pos2.n);
            }
            glEnd();
            j += stepsize;

            glBegin(GL_QUAD_STRIP);
            for (double i = right; i >= 0.0; i -= stepsize) {
                vec3 pos1 = vec3(i, j, k);
                vec3 pos2 = vec3(i, j + stepsize, k);

                vec4 color1 = getRenderColor(pos1);
                vec4 color2 = getRenderColor(pos2);

                glColor4dv(color1.n);
                glVertex3dv(pos1.n);

                glColor4dv(color2.n);
                glVertex3dv(pos2.n);
            }
            glEnd();
            j += stepsize;
        }
    }
}

void MACGrid::drawXSheets(bool backToFront) {
    // Draw K Sheets from back to front
    double back = (theDim[2]) * theCellSize;
    double top = (theDim[1]) * theCellSize;
    double right = (theDim[0]) * theCellSize;

    double stepsize = theCellSize * 0.25;

    double starti = right - stepsize;
    double endi = 0;
    double stepi = -theCellSize;

    if (!backToFront) {
        starti = 0;
        endi = right;
        stepi = theCellSize;
    }

    for (double i = starti; backToFront ? i > endi : i < endi; i += stepi) {
        for (double j = 0.0; j < top;) {
            glBegin(GL_QUAD_STRIP);
            for (double k = 0.0; k <= back; k += stepsize) {
                vec3 pos1 = vec3(i, j, k);
                vec3 pos2 = vec3(i, j + stepsize, k);

                vec4 color1 = getRenderColor(pos1);
                vec4 color2 = getRenderColor(pos2);

                glColor4dv(color1.n);
                glVertex3dv(pos1.n);

                glColor4dv(color2.n);
                glVertex3dv(pos2.n);
            }
            glEnd();
            j += stepsize;

            glBegin(GL_QUAD_STRIP);
            for (double k = back; k >= 0.0; k -= stepsize) {
                vec3 pos1 = vec3(i, j, k);
                vec3 pos2 = vec3(i, j + stepsize, k);

                vec4 color1 = getRenderColor(pos1);
                vec4 color2 = getRenderColor(pos2);

                glColor4dv(color1.n);
                glVertex3dv(pos1.n);

                glColor4dv(color2.n);
                glVertex3dv(pos2.n);
            }
            glEnd();
            j += stepsize;
        }
    }
}


void MACGrid::drawSmoke(const Camera &c) {
    vec3 eyeDir = c.getBackward();
    double zresult = fabs(Dot(eyeDir, vec3(1, 0, 0)));
    double xresult = fabs(Dot(eyeDir, vec3(0, 0, 1)));
    //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

    if (zresult < xresult) {
        drawZSheets(c.getPosition()[2] < 0);
    } else {
        drawXSheets(c.getPosition()[0] < 0);
    }
}

void MACGrid::drawSmokeCubes(const Camera &c) {
    std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
    FOR_EACH_CELL {
                MACGrid::Cube cube;
                cube.color = getRenderColor(i, j, k);
                cube.pos = getCenter(i, j, k);
                cube.dist = DistanceSqr(cube.pos, c.getPosition());
                cubes.insert(make_pair(cube.dist, cube));
            }

    // Draw cubes from back to front
    std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
    for (it = cubes.begin(); it != cubes.end(); ++it) {
        drawCube(it->second);
    }
}

void MACGrid::drawWireGrid() {
    // Display grid in light grey, draw top & bottom

    double xstart = 0.0;
    double ystart = 0.0;
    double zstart = 0.0;
    double xend = theDim[0] * theCellSize;
    double yend = theDim[1] * theCellSize;
    double zend = theDim[2] * theCellSize;

    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glColor3f(0.25, 0.25, 0.25);

    glBegin(GL_LINES);
    for (int i = 0; i <= theDim[0]; i++) {
        double x = xstart + i * theCellSize;
        glVertex3d(x, ystart, zstart);
        glVertex3d(x, ystart, zend);

        glVertex3d(x, yend, zstart);
        glVertex3d(x, yend, zend);
    }

    for (int i = 0; i <= theDim[2]; i++) {
        double z = zstart + i * theCellSize;
        glVertex3d(xstart, ystart, z);
        glVertex3d(xend, ystart, z);

        glVertex3d(xstart, yend, z);
        glVertex3d(xend, yend, z);
    }

    glVertex3d(xstart, ystart, zstart);
    glVertex3d(xstart, yend, zstart);

    glVertex3d(xend, ystart, zstart);
    glVertex3d(xend, yend, zstart);

    glVertex3d(xstart, ystart, zend);
    glVertex3d(xstart, yend, zend);

    glVertex3d(xend, ystart, zend);
    glVertex3d(xend, yend, zend);
    glEnd();
    glPopAttrib();

    glEnd();
}

#define LEN 0.5

void MACGrid::drawFace(const MACGrid::Cube &cube) {
    glColor4dv(cube.color.n);
    glPushMatrix();
    glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);
    glScaled(theCellSize, theCellSize, theCellSize);
    glBegin(GL_QUADS);
    glNormal3d(0.0, 0.0, 1.0);
    glVertex3d(-LEN, -LEN, LEN);
    glVertex3d(-LEN, LEN, LEN);
    glVertex3d(LEN, LEN, LEN);
    glVertex3d(LEN, -LEN, LEN);
    glEnd();
    glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube &cube) {
    glColor4dv(cube.color.n);
    glPushMatrix();
    glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);
    glScaled(theCellSize, theCellSize, theCellSize);
    glBegin(GL_QUADS);
    glNormal3d(0.0, -1.0, 0.0);
    glVertex3d(-LEN, -LEN, -LEN);
    glVertex3d(-LEN, -LEN, LEN);
    glVertex3d(LEN, -LEN, LEN);
    glVertex3d(LEN, -LEN, -LEN);

    glNormal3d(0.0, 0.0, -0.0);
    glVertex3d(-LEN, -LEN, -LEN);
    glVertex3d(-LEN, LEN, -LEN);
    glVertex3d(LEN, LEN, -LEN);
    glVertex3d(LEN, -LEN, -LEN);

    glNormal3d(-1.0, 0.0, 0.0);
    glVertex3d(-LEN, -LEN, -LEN);
    glVertex3d(-LEN, -LEN, LEN);
    glVertex3d(-LEN, LEN, LEN);
    glVertex3d(-LEN, LEN, -LEN);

    glNormal3d(0.0, 1.0, 0.0);
    glVertex3d(-LEN, LEN, -LEN);
    glVertex3d(-LEN, LEN, LEN);
    glVertex3d(LEN, LEN, LEN);
    glVertex3d(LEN, LEN, -LEN);

    glNormal3d(0.0, 0.0, 1.0);
    glVertex3d(-LEN, -LEN, LEN);
    glVertex3d(-LEN, LEN, LEN);
    glVertex3d(LEN, LEN, LEN);
    glVertex3d(LEN, -LEN, LEN);

    glNormal3d(1.0, 0.0, 0.0);
    glVertex3d(LEN, -LEN, -LEN);
    glVertex3d(LEN, -LEN, LEN);
    glVertex3d(LEN, LEN, LEN);
    glVertex3d(LEN, LEN, -LEN);
    glEnd();
    glPopMatrix();
}
