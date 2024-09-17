#include "prostt5.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    ProstT5 *model = prostt5_load("model/", false, true, false);
    if (model != NULL) {
        const char *result = prostt5_predict(model, "PIAQIHILEGRSDEQKETLIREVS");
        if (result != NULL) {
            printf("Prediction: %s\n", result);
            free((void*)result); // Free the CString after use
        }
        prostt5_free(model);
    } else {
        printf("Failed to load model.\n");
    }
    return 0;
}
