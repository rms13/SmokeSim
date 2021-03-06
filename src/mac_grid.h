#ifndef MACGrid_H_
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)

#include "open_gl_headers.h" 
#include "vec.h"
#include "grid_data.h"
#include "grid_data_matrix.h" 
#include <Partio.h>
#include <vector>
class Camera;

class MACGrid
{

public:
	MACGrid();
	~MACGrid();
	MACGrid(const MACGrid& orig);
	MACGrid& operator=(const MACGrid& orig);

	void reset();

	void draw(const Camera& c);
	void updateSources();
	void advectVelocity(double dt);
	void addExternalForces(double dt);
	void project(double dt);
	void advectTemperature(double dt);
	void advectDensity(double dt);
	void advectRenderingParticles(double dt);

	void initializeSolids();

    // MULTITHREADING: VELOCITY
    void advectVelocityThreadX(int tid, double dt);
    void advectVelocityThreadY(int tid, double dt);
    void advectVelocityThreadZ(int tid, double dt);
    // MULTITHREADING: TEMPERATURE
    void advectTemperatureThread(int tid, double dt, int kmin, int kmax);
    // MULTITHREADING: DENSITY
    void advectDensityThread(int tid, double dt, int kmin, int kmax);
    // MULTITHREADING: BUOYANCY
    void computeBouyancyThread(int tid, double dt, int kmin, int kmax);
    // MULTITHREADING: VORTICITY
    void computeVorticityConfinementThread_CellCenterVel(int tid, double dt, int kmin, int kmax,
                                                         GridData *u, GridData *v, GridData *w);
    void computeVorticityConfinementThread_CellCenterVort(int tid, double dt, int kmin, int kmax,
                                                          GridData *u, GridData *v, GridData *w,
                                                          GridData *x, GridData *y, GridData *z, GridData *len);
    void computeVorticityConfinementThread_Update(int tid, double dt, int kmin, int kmax,
                                                  GridData *x, GridData *y, GridData *z, GridData *len);
    // MULTITHREADING: PROJECTION
    void projectThread_Divergence(int tid, double dt, int kmin, int kmax,
                                  GridData *d);
    void projectThread_UpdateX(int tid, double dt, GridData *p);
    void projectThread_UpdateY(int tid, double dt, GridData *p);
    void projectThread_UpdateZ(int tid, double dt, GridData *p);

    // MULTITHREADING: calculatePreconditionerThread
    void calculatePreconditionerThread(int tid, GridDataMatrix *A, int kmin, int kmax);

    // MULTITHREADING: advectRenderingParticles
    void advectRenderingParticlesThread(int tid, double dt, int kmin, int kmax);

public:

	// rendering particles
	std::vector<vec3> rendering_particles;
	std::vector<vec3> rendering_particles_vel;
	std::vector<int> rendering_particles_colIdx;

	enum RenderMode { CUBES, SHEETS };
	static RenderMode theRenderMode;
	static bool theDisplayVel;

	void saveSmoke(const char* fileName);
	void saveParticle(std::string filename);
	void saveDensity(std::string filename);

protected:

	// Setup
	void initialize();

	// Simulation
	void computeBuoyancy(double dt);
	void computeVorticityConfinement(double dt);

	// Rendering
	struct Cube { vec3 pos; vec4 color; double dist; };
	void drawWireGrid();
	void drawSmokeCubes(const Camera& c);
	void drawSmoke(const Camera& c);
	void drawCube(const MACGrid::Cube& c);
	void drawFace(const MACGrid::Cube& c);
	void drawVelocities();
	vec4 getRenderColor(int i, int j, int k);
	vec4 getRenderColor(const vec3& pt);
	void drawZSheets(bool backToFront);
	void drawXSheets(bool backToFront);

	// GridData accessors
	enum Direction { X, Y, Z };
	vec3 getVelocity(const vec3& pt);
	double getVelocityX(const vec3& pt);
	double getVelocityY(const vec3& pt);
	double getVelocityZ(const vec3& pt);
	double getTemperature(const vec3& pt);
	double getDensity(const vec3& pt);
	vec3 getCenter(int i, int j, int k);

	
	vec3 getRewoundPosition(const vec3 & currentPosition, const double dt);
	vec3 clipToGrid(const vec3& outsidePoint, const vec3& insidePoint);
	double getSize(int dimension);
	int getCellIndex(int i, int j, int k);
	int getNumberOfCells();
	bool isValidCell(int i, int j, int k);
	bool isValidFace(int dimension, int i, int j, int k);
	vec3 getFacePosition(int dimension, int i, int j, int k);
	void calculateAMatrix();
	bool preconditionedConjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
	void calculatePreconditioner(GridDataMatrix & A);
	void applyPreconditioner(const GridData & r, const GridDataMatrix & A, GridData & z);
	double dotProduct(const GridData & vector1, const GridData & vector2);
	void add(const GridData & vector1, const GridData & vector2, GridData & result);
	void subtract(const GridData & vector1, const GridData & vector2, GridData & result);
	void multiply(const double scalar, const GridData & vector, GridData & result);
	double maxMagnitude(const GridData & vector);
	void apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result);


	GridDataX mU; // X component of velocity, stored on X faces, size is (dimX+1)*dimY*dimZ
	GridDataY mV; // Y component of velocity, stored on Y faces, size is dimX*(dimY+1)*dimZ
	GridDataZ mW; // Z component of velocity, stored on Z faces, size is dimX*dimY*(dimZ+1)
	GridData mP;  // Pressure, stored at grid centers, size is dimX*dimY*dimZ
	GridData mD;  // Density, stored at grid centers, size is dimX*dimY*dimZ
	GridData mT;  // Temperature, stored at grid centers, size is dimX*dimY*dimZ

    GridData solidCells;
	
	GridDataMatrix AMatrix;
	GridData precon;

	int currentStep;
};

#endif