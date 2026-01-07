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

    double *out = (double *)malloc(sizeof(double) * numberOfValues);
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
