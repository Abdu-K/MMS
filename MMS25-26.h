#ifndef MMS25_26
#define MMS25_26
#define BLOCK_SIZE 500

typedef struct {
    int numberOfBins;//-1
    int *bins;//0
    double minimum;//0
    double maximum;//0
    double binWidth;//0
} Histogram;

typedef struct {
    int numberOfMinimumPositions;
    int *minimumPositionArray;
    int numberOfMaximumPositions;
    int *maximumPositionArray;
} LocalExtrema;

typedef struct {
    int numberOfSamples;
    double *samples;
    double area;
    double mean;
    LocalExtrema *localExtrema;
} MMSignal;

//AUFGABE 1
double interpolateLine(double x1, double y1, double x2, double y2, double xb);
double *scaleValuesInArray(int numberOfValues, double *values, double minimum, double scalingFactor);
double *createSineArray(int totalSamples, int samplesPerPeriod, double amplitude);
int writeArrayFile(char *filePath, double *array, int arrayLength);
int readArrayFile(char *fileName, double *values);
MMSignal *createSignal_array(int numberOfValues, double *values);
MMSignal *createSignal_file(char *fileName);
void deleteMMSignal(MMSignal *In);
void writeSignal(MMSignal *In, char * fileName);
MMSignal *createSineSignal(int totalSamples, int samplesPerPeriod, double amplitude);

//AUFGABE 2
int *getHistogram(int numberOfValues, double *values, int numberOfBins);
Histogram *createHistogram_empty();
Histogram *createHistogram_bins(int numberOfBins);
Histogram *createHistogram_array(int numberOfValues, double *values, int numberOfBins);
void deleteHistogram(Histogram *In);
double computeArea(MMSignal *In);
double computeMean(MMSignal *In);
double computeStandardDeviation(MMSignal *In);
double computeMedian(MMSignal *In);
LocalExtrema *computeExtrema(MMSignal *In);
double computeEntropy(Histogram *histogramIn);

//Aufgabe 3
MMSignal *convoluteSignals(MMSignal *In1, MMSignal *In2);
MMSignal *approximateGaussianBellCurve(int pascalLineNumber); 

//Aufgabe 4
//Direction = 1 Forward FourierTransformation, Direction=-1 Backward FourierTransformation
void dft(int numberOfValues, double* realIn, double* imaginaryIn, double* realOut, double* imaginaryOut,int Direction);
void getCartesianToPolar(int numberOfValues, double* realIn, double* imaginaryIn, double* amplitudesOut, double* angelsOut);
void getPolarToCartesian(int numberOfValues, double* amplitudesIn, double* angelsIn, double* realOut, double* imaginaryOut);
#endif


