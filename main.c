#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MMS25-26.h"

static void print_double_array(const char *label, const double *arr, int len)
{
    printf("%s", label);
    for (int i = 0; i < len; i++)
        printf(" %.5f", arr[i]);
    printf("\n");
}

static void print_int_array(const char *label, const int *arr, int len)
{
    printf("%s", label);
    for (int i = 0; i < len; i++)
        printf(" %d", arr[i]);
    printf("\n");
}

int main(void)
{
    printf("--- Aufgabe 1 ---\n");
    // F1: Interpolation
    printf("interpolateLine: %.2f\n", interpolateLine(0, 0, 10, 10, 5));

    // F2: Skalierung
    double raw[5] = {0, 1, 2, 3, 4};
    double *scaled = scaleValuesInArray(5, raw, -1.0, 0.5);
    print_double_array("scaled:", scaled, 5);
    free(scaled);

    // F3: Sinusarray
    double *sin_arr = createSineArray(20, 10, 1.0);
    print_double_array("sine array:", sin_arr, 10);

    // F4: Array als Datei schreiben
    int written = writeArrayFile("sin_array.txt", sin_arr, 20);
    printf("writeArrayFile wrote: %d\n", written);

    // F5: Array aus Datei lesen
    double buf[BLOCK_SIZE];
    int read_count = readArrayFile("sin_array.txt", buf);
    printf("readArrayFile read: %d\n", read_count);

    // F6: Signal aus Array erzeugen
    MMSignal *sig_from_array = createSignal_array(5, raw);
    if (sig_from_array)
        printf("createSignal_array samples=%d\n", sig_from_array->numberOfSamples);

    // F7: Signal aus Datei erzeugen (Headerformat)
    MMSignal *sine_signal = createSineSignal(16, 8, 2.0);
    writeSignal(sine_signal, "sine_signal.txt");
    MMSignal *sig_from_file = createSignal_file("sine_signal.txt");
    if (sig_from_file)
        printf("createSignal_file samples=%d first=%.4f\n", sig_from_file->numberOfSamples, sig_from_file->samples[0]);

    // F8: Signal freigeben
    deleteMMSignal(sig_from_array);
    deleteMMSignal(sig_from_file);
    deleteMMSignal(sine_signal);
    free(sin_arr);

    printf("\n--- Aufgabe 2 ---\n");
    // F11: Histogramm als Array
    double values[8] = {1.2, 1.5, 2.0, 2.0, 2.4, 3.1, 3.2, 4.0};
    int *bins = getHistogram(8, values, 4);
    print_int_array("getHistogram bins:", bins, 4);
    free(bins);

    // F12/F13/F14/F15: Histogramm Struct
    Histogram *hist_empty = createHistogram_empty();
    if (hist_empty)
        printf("createHistogram_empty bins=%d\n", hist_empty->numberOfBins);
    deleteHistogram(hist_empty);

    Histogram *hist_bins = createHistogram_bins(5);
    if (hist_bins)
        printf("createHistogram_bins bins=%d\n", hist_bins->numberOfBins);
    deleteHistogram(hist_bins);

    Histogram *hist_array = createHistogram_array(8, values, 4);
    if (hist_array)
    {
        print_int_array("createHistogram_array bins:", hist_array->bins, hist_array->numberOfBins);
        printf("computeEntropy: %.5f\n", computeEntropy(hist_array));
    }
    deleteHistogram(hist_array);

    // F16-F20: Signal Analyse
    MMSignal *analysis_sig = createSignal_array(8, values);
    printf("computeArea: %.4f\n", computeArea(analysis_sig));
    printf("computeMean: %.4f\n", computeMean(analysis_sig));
    printf("computeStandardDeviation: %.4f\n", computeStandardDeviation(analysis_sig));
    printf("computeMedian: %.4f\n", computeMedian(analysis_sig));

    LocalExtrema *ext = computeExtrema(analysis_sig);
    if (ext)
    {
        print_int_array("minima:", ext->minimumPositionArray, ext->numberOfMinimumPositions);
        print_int_array("maxima:", ext->maximumPositionArray, ext->numberOfMaximumPositions);
    }
    deleteMMSignal(analysis_sig);

    printf("\n--- Aufgabe 3 ---\n");
    // F22: Faltung
    double a_vals[3] = {1, 2, 1};
    double b_vals[2] = {1, -1};
    MMSignal *a_sig = createSignal_array(3, a_vals);
    MMSignal *b_sig = createSignal_array(2, b_vals);
    MMSignal *conv = convoluteSignals(a_sig, b_sig);
    if (conv)
    {
        print_double_array("convolution:", conv->samples, conv->numberOfSamples);
    }
    deleteMMSignal(a_sig);
    deleteMMSignal(b_sig);
    deleteMMSignal(conv);

    // F23: Gauss-Approximation
    MMSignal *gauss = approximateGaussianBellCurve(6);
    if (gauss)
        print_double_array("gaussian:", gauss->samples, gauss->numberOfSamples);
    deleteMMSignal(gauss);

    printf("\n--- Aufgabe 4 ---\n");
    // F24/F25/F26: DFT und Polar/Kartesisch
    int n = 8;
    double real_in[8] = {1, 0, -1, 0, 1, 0, -1, 0};
    double imag_in[8] = {0};
    double real_out[8] = {0};
    double imag_out[8] = {0};

    dft(n, real_in, imag_in, real_out, imag_out, 1);
    print_double_array("DFT real:", real_out, n);
    print_double_array("DFT imag:", imag_out, n);

    double amp[8] = {0};
    double ang[8] = {0};
    getCartesianToPolar(n, real_out, imag_out, amp, ang);
    print_double_array("amplitudes:", amp, n);
    print_double_array("angles:", ang, n);

    double real_back[8] = {0};
    double imag_back[8] = {0};
    getPolarToCartesian(n, amp, ang, real_back, imag_back);
    print_double_array("reconstructed real:", real_back, n);
    print_double_array("reconstructed imag:", imag_back, n);

    return 0;
}
