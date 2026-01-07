#include <stdio.h>
#include <stdlib.h>
#include "MMS25-26.h"

int main(void)
{
    // 1) Interpolation
    printf("interpolateLine: %.2f\n", interpolateLine(0, 0, 10, 10, 5)); // 5.00

    // 2) Skalierung
    double raw[5] = {0, 1, 2, 3, 4};
    double *scaled = scaleValuesInArray(5, raw, 0.0, 0.5);
    for (int i = 0; i < 5; i++)
        printf("scaled[%d]=%.2f\n", i, scaled[i]);
    free(scaled);

    // 3) Sinusarray + 4) Export + 5) Import
    double *sinArr = createSineArray(50, 10, 1.0);
    int w = writeArrayFile("sin_array.txt", sinArr, 50);
    printf("writeArrayFile wrote: %d\n", w);

    double buf[BLOCK_SIZE];
    int r = readArrayFile("sin_array.txt", buf);
    printf("readArrayFile read: %d\n", r);

    free(sinArr);

    // Signal schreiben/lesen
    MMSignal *s1 = createSineSignal(40, 20, 2.0);
    writeSignal(s1, "sin_signal.txt");
    deleteMMSignal(s1);

    MMSignal *s2 = createSignal_file("sin_signal.txt");
    if (s2)
    {
        printf("createSignal_file: samples=%d first=%.3f\n", s2->numberOfSamples, s2->samples[0]);
        deleteMMSignal(s2);
    }

    return 0;
}
