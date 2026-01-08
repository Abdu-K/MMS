#define _USE_MATH_DEFINES
#include "MMS25-26.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

// ---------------------------
// Helper: whitespace trim
// ---------------------------
static char *trim(char *s)
{
    if (!s) return s;
    while (*s == ' ' || *s == '\t' || *s == '\n' || *s == '\r') s++;

    char *end = s + strlen(s);
    while (end > s && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\n' || end[-1] == '\r'))
        end--;
    *end = '\0';
    return s;
}

static int maxInt(int a, int b)
{
    return (a > b) ? a : b;
}

static int minInt(int a, int b)
{
    return (a < b) ? a : b;
}

// ====================================================
// AUFGABE 1
// ====================================================

// Funktion 1: lineare Interpolation
double interpolateLine(double x1, double y1, double x2, double y2, double xb)
{
    if (x2 == x1) return y1; // Schutz gegen Division durch 0
    return y1 + ((y2 - y1) / (x2 - x1)) * (xb - x1);
}

// Funktion 2: Werte umskalieren (neues Minimum + Faktor)
double *scaleValuesInArray(int numberOfValues, double *values, double minimum, double scalingFactor)
{
    if (values == NULL || numberOfValues <= 0) return NULL;

    double *out = malloc(sizeof(double) * numberOfValues);
    if (out == NULL) return NULL;

    for (int i = 0; i < numberOfValues; i++)
        out[i] = minimum + values[i] * scalingFactor;

    return out; // caller: free(out)
}

// Funktion 3: Sinusfolge erzeugen
double *createSineArray(int totalSamples, int samplesPerPeriod, double amplitude)
{
    if (totalSamples <= 0 || samplesPerPeriod <= 0) return NULL;

    double *arr = (double *)malloc(sizeof(double) * totalSamples);
    if (arr == NULL) return NULL;

    for (int i = 0; i < totalSamples; i++)
    {
        double phase = 2.0 * M_PI * (double)i / (double)samplesPerPeriod;
        arr[i] = amplitude * sin(phase);
    }

    return arr; // caller: free(arr)
    // return the array of the
}

// Funktion 4: Array als ASCII in Datei (eine Zahl pro Zeile)
int writeArrayFile(char *filePath, double *array, int arrayLength)
{
    if (filePath == NULL || array == NULL || arrayLength <= 0) return -1;

    FILE *f = fopen(filePath, "w");
    if (!f) return -1;

    for (int i = 0; i < arrayLength; i++)
    {
        // viele signifikante Stellen => robust beim Einlesen
        if (fprintf(f, "%.17g\n", array[i]) < 0)
        {
            fclose(f);
            return -1;
        }
    }

    fclose(f);
    return arrayLength; // geschrieben
}

// Funktion 5: Array aus Datei einlesen (gemäß Funktion 4)
// Liest bis BLOCK_SIZE oder EOF. values muss BLOCK_SIZE groß sein.
int readArrayFile(char *fileName, double *values)
{
    if (fileName == NULL || values == NULL) return -1;

    FILE *f = fopen(fileName, "r");
    if (!f) return -1;

    int count = 0;
    char line[256];

    while (count < BLOCK_SIZE && fgets(line, (int)sizeof(line), f))
    {
        char *t = trim(line);
        if (*t == '\0') continue; // leere Zeile ignorieren

        errno = 0;
        char *endptr = NULL;
        double v = strtod(t, &endptr);

        if (endptr == t || errno != 0)
            continue; // ungültige Zeile ignorieren

        values[count++] = v;
    }

    fclose(f);
    return count; // wie viele gelesen
}

// ----------------------------
// Signal-Objekt Funktionen
// ----------------------------

// Signal aus Array erzeugen (kopiert Werte)
MMSignal *createSignal_array(int numberOfValues, double *values)
{
    if (values == NULL || numberOfValues <= 0) return NULL;

    MMSignal *sig = (MMSignal *)malloc(sizeof(MMSignal));
    if (!sig) return NULL;

    sig->samples = (double *)malloc(sizeof(double) * numberOfValues);
    if (!sig->samples)
    {
        free(sig);
        return NULL;
    }

    for (int i = 0; i < numberOfValues; i++)
        sig->samples[i] = values[i];

    sig->numberOfSamples = numberOfValues;

    // Merkmale werden erst in Aufgabe 2 berechnet -> Initialwerte
    sig->area = 0.0;
    sig->mean = 0.0;
    sig->localExtrema = NULL;

    return sig;
}

// Signal aus Datei erzeugen
// Unterstützt zwei Formate:
// A) Nur Werte (eine Zahl/Zeile)  -> alle Zahlen sind Samples
// B) Erste Zeile = Integer N      -> danach N Zahlen als Samples
MMSignal *createSignal_file(char *fileName)
{
    if (fileName == NULL) return NULL;

    FILE *f = fopen(fileName, "r");
    if (!f) return NULL;

    char firstLine[256];
    if (!fgets(firstLine, (int)sizeof(firstLine), f))
    {
        fclose(f);
        return NULL;
    }

    char *firstTrim = trim(firstLine);
    if (*firstTrim == '\0')
    {
        fclose(f);
        return NULL;
    }

    // Kandidat: Header = N
    char *endInt = NULL;
    long nHeader = strtol(firstTrim, &endInt, 10);
    int firstIsPureInt = (endInt != firstTrim && trim(endInt)[0] == '\0' && nHeader >= 0);

    // Restliche Werte dynamisch einlesen
    int cap = 256;
    int n = 0;
    double *tmp = (double *)malloc(sizeof(double) * cap);
    if (!tmp)
    {
        fclose(f);
        return NULL;
    }

    char line[256];
    while (fgets(line, (int)sizeof(line), f))
    {
        char *t = trim(line);
        if (*t == '\0') continue;

        errno = 0;
        char *endptr = NULL;
        double v = strtod(t, &endptr);
        if (endptr == t || errno != 0) continue;

        if (n >= cap)
        {
            cap *= 2;
            double *grown = (double *)realloc(tmp, sizeof(double) * cap);
            if (!grown)
            {
                free(tmp);
                fclose(f);
                return NULL;
            }
            tmp = grown;
        }
        tmp[n++] = v;
    }

    fclose(f);

    // Entscheiden: Header-Format oder reines Array?
    if (firstIsPureInt && (long)n == nHeader)
    {
        // Format B: Header war N, tmp enthält genau N Samples
        MMSignal *sig = createSignal_array(n, tmp);
        free(tmp);
        return sig;
    }
    else
    {
        // Format A: erste Zeile ist Sample + tmp sind weitere Samples
        errno = 0;
        char *endD = NULL;
        double firstSample = strtod(firstTrim, &endD);
        if (endD == firstTrim || errno != 0)
        {
            free(tmp);
            return NULL;
        }

        double *all = (double *)malloc(sizeof(double) * (n + 1));
        if (!all)
        {
            free(tmp);
            return NULL;
        }

        all[0] = firstSample;
        for (int i = 0; i < n; i++)
            all[i + 1] = tmp[i];

        MMSignal *sig = createSignal_array(n + 1, all);

        free(all);
        free(tmp);
        return sig;
    }
}

// Signal in Datei schreiben (mit Header: numberOfSamples in erster Zeile)
void writeSignal(MMSignal *In, char *fileName)
{
    if (In == NULL || fileName == NULL || In->samples == NULL || In->numberOfSamples <= 0) return;

    FILE *f = fopen(fileName, "w");
    if (!f) return;

    fprintf(f, "%d\n", In->numberOfSamples);
    for (int i = 0; i < In->numberOfSamples; i++)
        fprintf(f, "%.17g\n", In->samples[i]);

    fclose(f);
}

// Speicher freigeben
void deleteMMSignal(MMSignal *In)
{
    if (!In) return;

    if (In->samples) free(In->samples);

    if (In->localExtrema)
    {
        if (In->localExtrema->minimumPositionArray) free(In->localExtrema->minimumPositionArray);
        if (In->localExtrema->maximumPositionArray) free(In->localExtrema->maximumPositionArray);
        free(In->localExtrema);
    }

    free(In);
}

// Sinus als MMSignal erzeugen
MMSignal *createSineSignal(int totalSamples, int samplesPerPeriod, double amplitude)
{
    double *arr = createSineArray(totalSamples, samplesPerPeriod, amplitude);
    if (!arr) return NULL;

    MMSignal *sig = createSignal_array(totalSamples, arr);
    free(arr);
    return sig;
}

// ====================================================
// AUFGABE 2
// ====================================================

int *getHistogram(int numberOfValues, double *values, int numberOfBins)
{
    if (values == NULL || numberOfValues <= 0 || numberOfBins <= 0) return NULL;

    double min = values[0];
    double max = values[0];
    for (int i = 1; i < numberOfValues; i++)
    {
        if (values[i] < min) min = values[i];
        if (values[i] > max) max = values[i];
    }

    int *bins = (int *)calloc((size_t)numberOfBins, sizeof(int));
    if (!bins) return NULL;

    if (max == min)
    {
        bins[0] = numberOfValues;
        return bins;
    }

    double binWidth = (max - min) / (double)numberOfBins;
    for (int i = 0; i < numberOfValues; i++)
    {
        int idx = (int)((values[i] - min) / binWidth);
        if (idx < 0) idx = 0;
        if (idx >= numberOfBins) idx = numberOfBins - 1;
        bins[idx]++;
    }

    return bins;
}

Histogram *createHistogram_empty()
{
    Histogram *hist = (Histogram *)malloc(sizeof(Histogram));
    if (!hist) return NULL;

    hist->numberOfBins = -1;
    hist->bins = NULL;
    hist->minimum = 0.0;
    hist->maximum = 0.0;
    hist->binWidth = 0.0;
    return hist;
}

Histogram *createHistogram_bins(int numberOfBins)
{
    if (numberOfBins <= 0) return NULL;

    Histogram *hist = (Histogram *)malloc(sizeof(Histogram));
    if (!hist) return NULL;

    hist->numberOfBins = numberOfBins;
    hist->minimum = 0.0;
    hist->maximum = 0.0;
    hist->binWidth = 0.0;
    hist->bins = (int *)calloc((size_t)numberOfBins, sizeof(int));
    if (!hist->bins)
    {
        free(hist);
        return NULL;
    }
    return hist;
}

Histogram *createHistogram_array(int numberOfValues, double *values, int numberOfBins)
{
    if (values == NULL || numberOfValues <= 0 || numberOfBins <= 0) return NULL;

    double min = values[0];
    double max = values[0];
    for (int i = 1; i < numberOfValues; i++)
    {
        if (values[i] < min) min = values[i];
        if (values[i] > max) max = values[i];
    }

    Histogram *hist = createHistogram_bins(numberOfBins);
    if (!hist) return NULL;

    hist->minimum = min;
    hist->maximum = max;
    hist->binWidth = (max == min) ? 0.0 : (max - min) / (double)numberOfBins;

    if (max == min)
    {
        hist->bins[0] = numberOfValues;
        return hist;
    }

    for (int i = 0; i < numberOfValues; i++)
    {
        int idx = (int)((values[i] - min) / hist->binWidth);
        if (idx < 0) idx = 0;
        if (idx >= numberOfBins) idx = numberOfBins - 1;
        hist->bins[idx]++;
    }

    return hist;
}

void deleteHistogram(Histogram *In)
{
    if (!In) return;
    if (In->bins) free(In->bins);
    free(In);
}

double computeArea(MMSignal *In)
{
    if (In == NULL || In->samples == NULL || In->numberOfSamples <= 0) return 0.0;

    double area = 0.0;
    for (int i = 0; i < In->numberOfSamples; i++)
        area += In->samples[i];

    In->area = area;
    return area;
}

double computeMean(MMSignal *In)
{
    if (In == NULL || In->samples == NULL || In->numberOfSamples <= 0) return 0.0;

    double area = computeArea(In);
    double mean = area / (double)In->numberOfSamples;
    In->mean = mean;
    return mean;
}

double computeStandardDeviation(MMSignal *In)
{
    if (In == NULL || In->samples == NULL || In->numberOfSamples <= 0) return 0.0;

    double mean = computeMean(In);
    double sumSq = 0.0;
    for (int i = 0; i < In->numberOfSamples; i++)
    {
        double diff = In->samples[i] - mean;
        sumSq += diff * diff;
    }

    return sqrt(sumSq / (double)In->numberOfSamples);
}

double computeMedian(MMSignal *In)
{
    if (In == NULL || In->samples == NULL || In->numberOfSamples <= 0) return 0.0;

    int n = In->numberOfSamples;
    double *sorted = (double *)malloc(sizeof(double) * n);
    if (!sorted) return 0.0;

    for (int i = 0; i < n; i++)
        sorted[i] = In->samples[i];

    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (sorted[j] < sorted[i])
            {
                double tmp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = tmp;
            }
        }
    }

    double median = 0.0;
    if (n % 2 == 0)
        median = (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
    else
        median = sorted[n / 2];

    free(sorted);
    return median;
}

LocalExtrema *computeExtrema(MMSignal *In)
{
    if (In == NULL || In->samples == NULL || In->numberOfSamples < 3) return NULL;

    int n = In->numberOfSamples;
    int minCount = 0;
    int maxCount = 0;

    for (int i = 1; i < n - 1; i++)
    {
        if (In->samples[i - 1] > In->samples[i] && In->samples[i] < In->samples[i + 1])
            minCount++;
        if (In->samples[i - 1] < In->samples[i] && In->samples[i] > In->samples[i + 1])
            maxCount++;
    }

    LocalExtrema *ext = (LocalExtrema *)malloc(sizeof(LocalExtrema));
    if (!ext) return NULL;

    ext->numberOfMinimumPositions = minCount;
    ext->numberOfMaximumPositions = maxCount;
    ext->minimumPositionArray = minCount > 0 ? (int *)malloc(sizeof(int) * minCount) : NULL;
    ext->maximumPositionArray = maxCount > 0 ? (int *)malloc(sizeof(int) * maxCount) : NULL;

    if ((minCount > 0 && !ext->minimumPositionArray) || (maxCount > 0 && !ext->maximumPositionArray))
    {
        if (ext->minimumPositionArray) free(ext->minimumPositionArray);
        if (ext->maximumPositionArray) free(ext->maximumPositionArray);
        free(ext);
        return NULL;
    }

    int minIdx = 0;
    int maxIdx = 0;
    for (int i = 1; i < n - 1; i++)
    {
        if (In->samples[i - 1] > In->samples[i] && In->samples[i] < In->samples[i + 1])
            ext->minimumPositionArray[minIdx++] = i;
        if (In->samples[i - 1] < In->samples[i] && In->samples[i] > In->samples[i + 1])
            ext->maximumPositionArray[maxIdx++] = i;
    }

    if (In->localExtrema)
    {
        if (In->localExtrema->minimumPositionArray) free(In->localExtrema->minimumPositionArray);
        if (In->localExtrema->maximumPositionArray) free(In->localExtrema->maximumPositionArray);
        free(In->localExtrema);
    }
    In->localExtrema = ext;

    return ext;
}

double computeEntropy(Histogram *histogramIn)
{
    if (!histogramIn || !histogramIn->bins || histogramIn->numberOfBins <= 0) return 0.0;

    int total = 0;
    for (int i = 0; i < histogramIn->numberOfBins; i++)
        total += histogramIn->bins[i];

    if (total == 0) return 0.0;

    double entropy = 0.0;
    for (int i = 0; i < histogramIn->numberOfBins; i++)
    {
        if (histogramIn->bins[i] == 0) continue;
        double p = (double)histogramIn->bins[i] / (double)total;
        entropy -= p * (log(p) / log(2.0));
    }

    return entropy;
}

// ====================================================
// AUFGABE 3
// ====================================================

MMSignal *convoluteSignals(MMSignal *In1, MMSignal *In2)
{
    if (!In1 || !In2) return NULL;
    if (!In1->samples || !In2->samples) return NULL;
    if (In1->numberOfSamples <= 0 || In2->numberOfSamples <= 0) return NULL;

    int n1 = In1->numberOfSamples;
    int n2 = In2->numberOfSamples;
    int nOut = n1 + n2 - 1;

    MMSignal *out = (MMSignal *)malloc(sizeof(MMSignal));
    if (!out) return NULL;

    out->samples = (double *)calloc((size_t)nOut, sizeof(double));
    if (!out->samples)
    {
        free(out);
        return NULL;
    }

    out->numberOfSamples = nOut;
    out->area = 0.0;
    out->mean = 0.0;
    out->localExtrema = NULL;

    for (int n = 0; n < nOut; n++)
    {
        int kMin = maxInt(0, n - (n2 - 1));
        int kMax = minInt(n, n1 - 1);
        double sum = 0.0;

        for (int k = kMin; k <= kMax; k++)
            sum += In1->samples[k] * In2->samples[n - k];

        out->samples[n] = sum;
    }

    return out;
}

MMSignal *approximateGaussianBellCurve(int pascalLineNumber)
{
    if (pascalLineNumber < 0) return NULL;

    int n = pascalLineNumber;
    int len = n + 1;

    MMSignal *sig = (MMSignal *)malloc(sizeof(MMSignal));
    if (!sig) return NULL;

    sig->samples = (double *)malloc(sizeof(double) * (size_t)len);
    if (!sig->samples)
    {
        free(sig);
        return NULL;
    }

    sig->numberOfSamples = len;
    sig->area = 0.0;
    sig->mean = 0.0;
    sig->localExtrema = NULL;

    long double c = 1.0L;
    long double sum = 0.0L;

    sig->samples[0] = 1.0;
    sum += c;

    for (int k = 1; k <= n; k++)
    {
        c = c * (long double)(n - (k - 1)) / (long double)k;
        sig->samples[k] = (double)c;
        sum += c;
    }

    if (sum > 0.0L)
    {
        for (int k = 0; k < len; k++)
            sig->samples[k] = (double)((long double)sig->samples[k] / sum);
    }

    return sig;
}

// ====================================================
// AUFGABE 4
// ====================================================

void dft(int numberOfValues,
         double *realIn, double *imaginaryIn,
         double *realOut, double *imaginaryOut,
         int Direction)
{
    if (numberOfValues <= 0) return;
    if (!realIn || !realOut || !imaginaryOut) return;

    int dir = (Direction == -1) ? -1 : 1;

    for (int l = 0; l < numberOfValues; l++)
    {
        double sumRe = 0.0;
        double sumIm = 0.0;

        for (int k = 0; k < numberOfValues; k++)
        {
            double a = realIn[k];
            double b = (imaginaryIn != NULL) ? imaginaryIn[k] : 0.0;

            double angle = (2.0 * M_PI * (double)l * (double)k) / (double)numberOfValues;
            angle *= (double)dir;

            double c = cos(angle);
            double s = sin(angle);

            sumRe += a * c - b * s;
            sumIm += a * s + b * c;
        }

        if (dir == -1)
        {
            sumRe /= (double)numberOfValues;
            sumIm /= (double)numberOfValues;
        }

        realOut[l] = sumRe;
        imaginaryOut[l] = sumIm;
    }
}

void getCartesianToPolar(int numberOfValues,
                         double *realIn, double *imaginaryIn,
                         double *amplitudesOut, double *angelsOut)
{
    if (numberOfValues <= 0) return;
    if (!realIn || !imaginaryIn || !amplitudesOut || !angelsOut) return;

    for (int i = 0; i < numberOfValues; i++)
    {
        double re = realIn[i];
        double im = imaginaryIn[i];
        amplitudesOut[i] = sqrt(re * re + im * im);
        angelsOut[i] = atan2(im, re);
    }
}

void getPolarToCartesian(int numberOfValues,
                         double *amplitudesIn, double *angelsIn,
                         double *realOut, double *imaginaryOut)
{
    if (numberOfValues <= 0) return;
    if (!amplitudesIn || !angelsIn || !realOut || !imaginaryOut) return;

    for (int i = 0; i < numberOfValues; i++)
    {
        double amplitude = amplitudesIn[i];
        double angle = angelsIn[i];
        realOut[i] = amplitude * cos(angle);
        imaginaryOut[i] = amplitude * sin(angle);
    }
}
