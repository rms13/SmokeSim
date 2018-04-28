

#ifndef CONSTANTS_H
#define CONSTANTS_H



#define LERP(a,b,t) (1-t)*a + t*b

extern const int numFrames;

// Don't modify the values of these here.
// Modify the values of these in Constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const double theCellSize;
extern const double theAirDensity;
extern const double theBuoyancyAlpha;
extern const double theBuoyancyBeta;	
extern const double theBuoyancyAmbientTemperature;
extern const double theVorticityEpsilon;





#endif