#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define SIZE 32
#define ENCRYPTED_SIZE 16
#define MAX_OUTPUTS 8000
typedef unsigned char u8;
typedef unsigned int u32;

/* --- Confusion table and diffusion values --- */
u8 confusion[512] = {
    0xac, 0xd1, 0x25, 0x94, 0x1f, 0xb3, 0x33, 0x28, 0x7c, 0x2b, 0x17, 0xbc, 0xf6, 0xb0, 0x55, 0x5d,
    0x8f, 0xd2, 0x48, 0xd4, 0xd3, 0x78, 0x62, 0x1a, 0x02, 0xf2, 0x01, 0xc9, 0xaa, 0xf0, 0x83, 0x71,
    0x72, 0x4b, 0x6a, 0xe8, 0xe9, 0x42, 0xc0, 0x53, 0x63, 0x66, 0x13, 0x4a, 0xc1, 0x85, 0xcf, 0x0c,
    0x24, 0x76, 0xa5, 0x6e, 0xd7, 0xa1, 0xec, 0xc6, 0x04, 0xc2, 0xa2, 0x5c, 0x81, 0x92, 0x6c, 0xda,
    0xc6, 0x86, 0xba, 0x4d, 0x39, 0xa0, 0x0e, 0x8c, 0x8a, 0xd0, 0xfe, 0x59, 0x96, 0x49, 0xe6, 0xea,
    0x69, 0x30, 0x52, 0x1c, 0xe0, 0xb2, 0x05, 0x9b, 0x10, 0x03, 0xa8, 0x64, 0x51, 0x97, 0x02, 0x09,
    0x8e, 0xad, 0xf7, 0x36, 0x47, 0xab, 0xce, 0x7f, 0x56, 0xca, 0x00, 0xe3, 0xed, 0xf1, 0x38, 0xd8,
    0x26, 0x1c, 0xdc, 0x35, 0x91, 0x43, 0x2c, 0x74, 0xb4, 0x61, 0x9d, 0x5e, 0xe9, 0x4c, 0xbf, 0x77,
    0x16, 0x1e, 0x21, 0x1d, 0x2d, 0xa9, 0x95, 0xb8, 0xc3, 0x8d, 0xf8, 0xdb, 0x34, 0xe1, 0x84, 0xd6,
    0x0b, 0x23, 0x4e, 0xff, 0x3c, 0x54, 0xa7, 0x78, 0xa4, 0x89, 0x33, 0x6d, 0xfb, 0x79, 0x27, 0xc4,
    0xf9, 0x40, 0x41, 0xdf, 0xc5, 0x82, 0x93, 0xdd, 0xa6, 0xef, 0xcd, 0x8d, 0xa3, 0xae, 0x7a, 0xb6,
    0x2f, 0xfd, 0xbd, 0xe5, 0x98, 0x66, 0xf3, 0x4f, 0x57, 0x88, 0x90, 0x9c, 0x0a, 0x50, 0xe7, 0x15,
    0x7b, 0x58, 0xbc, 0x07, 0x68, 0x3a, 0x5f, 0xee, 0x32, 0x9f, 0xeb, 0xcc, 0x18, 0x8b, 0xe2, 0x57,
    0xb7, 0x49, 0x37, 0xde, 0xf5, 0x99, 0x67, 0x5b, 0x3b, 0xbb, 0x3d, 0xb5, 0x2d, 0x19, 0x2e, 0x0d,
    0x93, 0xfc, 0x7e, 0x06, 0x08, 0xbe, 0x3f, 0xd9, 0x2a, 0x70, 0x9a, 0xc8, 0x7d, 0xd8, 0x46, 0x65,
    0x22, 0xf4, 0xb9, 0xa2, 0x6f, 0x12, 0x1b, 0x14, 0x45, 0xc7, 0x87, 0x31, 0x60, 0x29, 0xf7, 0x73,
    /* Second half (256 values) */
    0x2c, 0x97, 0x72, 0xcd, 0x89, 0xa6, 0x88, 0x4c, 0xe8, 0x83, 0xeb, 0x59, 0xca, 0x50, 0x3f, 0x27,
    0x4e, 0xae, 0x43, 0xd5, 0x6e, 0xd0, 0x99, 0x7b, 0x7c, 0x40, 0x0c, 0x52, 0x86, 0xc1, 0x46, 0x12,
    0x5a, 0x28, 0xa8, 0xbb, 0xcb, 0xf0, 0x11, 0x95, 0x26, 0x0d, 0x34, 0x66, 0x22, 0x18, 0x6f, 0x51,
    0x9b, 0x3b, 0xda, 0xec, 0x5e, 0x00, 0x2a, 0xf5, 0x8f, 0x61, 0xba, 0x96, 0xb3, 0xd1, 0x30, 0xdc,
    0x33, 0x75, 0xe9, 0x6d, 0xc8, 0xa1, 0x3a, 0x3e, 0x5f, 0x9d, 0xfd, 0xa9, 0x31, 0x9f, 0xaa, 0x85,
    0x2f, 0x92, 0xaf, 0x67, 0x78, 0xa5, 0xab, 0x03, 0x21, 0x4f, 0xb9, 0xad, 0xfe, 0xf3, 0x42, 0xfc,
    0x17, 0xd7, 0xee, 0xa3, 0xd8, 0x80, 0x14, 0x2e, 0xa0, 0x47, 0x55, 0xc4, 0xff, 0xe5, 0x13, 0x3f,
    0x81, 0xb6, 0x7a, 0x94, 0xd0, 0xb5, 0x54, 0xbf, 0x91, 0xa7, 0x37, 0xf1, 0x6b, 0xc9, 0x1b, 0xb1,
    0x3c, 0xb6, 0xd9, 0x32, 0x24, 0x8d, 0xf2, 0x82, 0xb4, 0xf9, 0xdb, 0x7d, 0x44, 0xfb, 0x1e, 0xd4,
    0xea, 0x5d, 0x35, 0x69, 0x23, 0x71, 0x57, 0x01, 0x06, 0xe4, 0x55, 0x9a, 0xa4, 0x58, 0x56, 0xc7,
    0x4a, 0x8c, 0x8a, 0xd6, 0x6a, 0x49, 0x70, 0xc5, 0x8e, 0x0a, 0x62, 0xdc, 0x29, 0x4b, 0x42, 0x41,
    0xcb, 0x2b, 0xb7, 0xce, 0x08, 0xa1, 0x76, 0x1d, 0x1a, 0xb8, 0xe3, 0xcc, 0x7e, 0x48, 0x20, 0xe6,
    0xf8, 0x45, 0x93, 0xde, 0xc3, 0x63, 0x0f, 0xb0, 0xac, 0x5c, 0xba, 0xdf, 0x07, 0x77, 0xe7, 0x4e,
    0x1f, 0x28, 0x10, 0x6c, 0x59, 0xd3, 0xdd, 0x2d, 0x65, 0x39, 0xb2, 0x74, 0x84, 0x3d, 0xf4, 0xbd,
    0xc7, 0x79, 0x60, 0x0b, 0x4d, 0x33, 0x36, 0x25, 0xbc, 0xe0, 0x09, 0xcf, 0x5b, 0xe2, 0x38, 0x9e,
    0xc0, 0xef, 0xd2, 0x16, 0x05, 0xbe, 0x53, 0xf7, 0xc2, 0xc6, 0xa2, 0x24, 0x98, 0x1c, 0xad, 0x04
};

u32 diffusion[SIZE] = {
    0xf26cb481, 0x16a5dc92, 0x3c5ba924, 0x79b65248, 0x2fc64b18, 0x615acd29, 0xc3b59a42, 0x976b2584,
    0x6cf281b4, 0xa51692dc, 0x5b3c24a9, 0xb6794852, 0xc62f184b, 0x5a6129cd, 0xb5c3429a, 0x6b978425,
    0xb481f26c, 0xdc9216a5, 0xa9243c5b, 0x524879b6, 0x4b182fc6, 0xcd29615a, 0x9a42c3b5, 0x2584976b,
    0x81b46cf2, 0x92dca516, 0x24a95b3c, 0x4852b679, 0x184bc62f, 0x29cd5a61, 0x429ab5c3, 0x84256b97
};

/* --- Bit-packed Matrix-Vector "Multiplication" ---*/
static  void matrix_vector_multiply_packed(const uint32_t M_packed[SIZE], const u8 c[SIZE], u8 d[SIZE]) {
    for (int j = 0; j < SIZE; j++) {
        uint32_t mask = M_packed[j];
        u8 sum = 0;
        while (mask) {
            const int k = __builtin_ctz(mask);  // index of the lowest set bit
            sum ^= c[k];
            mask &= mask - 1;  // clear the lowest set bit
        }
        d[j] = sum;
    }
}

/* --- Forward transformation (unchanged) --- */
void Forward(u8 c[SIZE], u8 d[SIZE], const u8 s[512], const u32 p[SIZE]) {
    for (u32 i = 0; i < 256; i++) {
        for (u8 j = 0; j < SIZE; j++) {
            d[j] = s[c[j]];
            c[j] = 0;
        }
        for (u8 j = 0; j < SIZE; j++) {
            for (u8 k = 0; k < SIZE; k++) {
                c[j] ^= d[k] * (p[j] >> k & 1);
            }
        }
    }
    for (u8 i = 0; i < ENCRYPTED_SIZE; i++)
        d[i] = s[c[i * 2]] ^ s[c[i * 2 + 1] + 256];
}

/* --- Diffusion matrix construction --- */
void construct_matrix(uint8_t M[SIZE][SIZE]) {
    for (int j = 0; j < SIZE; j++) {
        for (int k = 0; k < SIZE; k++) {
            M[j][k] = (diffusion[j] >> k) & 1;
        }
    }
}

/* --- Invert matrix over GF(2) using Gaussian elimination --- */
int invert_matrix(uint8_t M[SIZE][SIZE], uint8_t M_inv[SIZE][SIZE]) {
    uint8_t augmented[SIZE][SIZE * 2];
    int i, j;
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            augmented[i][j] = M[i][j];
            augmented[i][j + SIZE] = (i == j) ? 1 : 0;
        }
    }
    for (i = 0; i < SIZE; i++) {
        if (augmented[i][i] == 0) {
            int swap_row = i + 1;
            while (swap_row < SIZE && augmented[swap_row][i] == 0)
                swap_row++;
            if (swap_row == SIZE)
                return 0; // Matrix is singular
            for (j = 0; j < SIZE * 2; j++) {
                const u8 temp = augmented[i][j];
                augmented[i][j] = augmented[swap_row][j];
                augmented[swap_row][j] = temp;
            }
        }
        for (int row = 0; row < SIZE; row++) {
            if (row != i && augmented[row][i] == 1) {
                for (j = 0; j < SIZE * 2; j++) {
                    augmented[row][j] ^= augmented[i][j];
                }
            }
        }
    }
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            M_inv[i][j] = augmented[i][j + SIZE];
        }
    }
    return 1;
}

/* --- Dynamic array structures and helper functions --- */
typedef struct {
    u8 data[SIZE];
} Output;

typedef struct {
    Output *outputs;
    int count;
    int capacity;
} OutputList;

OutputList outputs;

void push_output(OutputList *list, const u8 candidate[SIZE]) {
    if (list->count >= list->capacity) {
        list->capacity *= 2;
        list->outputs = (Output*)realloc(list->outputs, list->capacity * sizeof(Output));
        if (!list->outputs) {
            fprintf(stderr, "Memory allocation error\n");
            exit(1);
        }
    }
    memcpy(list->outputs[list->count].data, candidate, SIZE);
    list->count++;
}

/* Structure for a pair (x,y) */
typedef struct {
    u8 first;
    u8 second;
} Pair;

typedef struct {
    Pair *pairs;
    int count;
    int capacity;
} PairList;

PairList pair_lookup[256];

void push_pair(PairList *list, const u8 first, const u8 second) {
    if (list->count >= list->capacity) {
        list->capacity *= 2;
        list->pairs = (Pair*)realloc(list->pairs, list->capacity * sizeof(Pair));
        if (!list->pairs) {
            fprintf(stderr, "Memory allocation error\n");
            exit(1);
        }
    }
    list->pairs[list->count].first = first;
    list->pairs[list->count].second = second;
    list->count++;
}

/* Structure for a dynamic array of u8 values */
typedef struct {
    u8 *values;
    int count;
    int capacity;
} U8List;

U8List inv_confusion_full[256];

void push_u8(U8List *list, const u8 value) {
    if (list->count >= list->capacity) {
        list->capacity *= 2;
        list->values = (u8*)realloc(list->values, list->capacity * sizeof(u8));
        if (!list->values) {
            fprintf(stderr, "Memory allocation error\n");
            exit(1);
        }
    }
    list->values[list->count] = value;
    list->count++;
}

/* --- Build the inverse confusion mapping --- */
void initialize_inv_confusion_full() {
    for (int i = 0; i < 256; i++) {
        /* Only the first 256 entries of confusion are used */
        push_u8(&inv_confusion_full[confusion[i]], i);
    }
}

/* --- Generate lookup table for the final XOR stage --- */
void generate_pair_lookup() {
    for (int x = 0; x < 256; x++) {
        for (int y = 0; y < 256; y++) {
            const u8 value = confusion[x] ^ confusion[y + 256];
            push_pair(&pair_lookup[value], x, y);
        }
    }
}

/* --- Recursive functions for decryption --- */

/* Global variable for the bit-packed inverse matrix */
uint32_t M_inv_packed[SIZE];

void Backward(u8 c[SIZE], const uint8_t M_inv[SIZE][SIZE], int round);

void BackwardConf(const u8 d[SIZE], u8 c[SIZE], const uint8_t M_inv[SIZE][SIZE], const int round, const int inner_round) {
    if (inner_round == SIZE) {
        u8 candidate[SIZE];
        memcpy(candidate, c, SIZE);
        Backward(candidate, M_inv, round);
        return;
    }
    const u8 original = c[inner_round];
    for (int i = 0; i < inv_confusion_full[d[inner_round]].count; i++) {
        const u8 inversion = inv_confusion_full[d[inner_round]].values[i];
        c[inner_round] = inversion;
        BackwardConf(d, c, M_inv, round, inner_round + 1);
    }
    c[inner_round] = original; // Backtrack
}

void Backward(u8 c[SIZE], const uint8_t M_inv[SIZE][SIZE], const int round) {
    if (round != 256) {
        u8 d[SIZE] = {0};
        /* Use the optimized bit-packed multiplication */
        matrix_vector_multiply_packed(M_inv_packed, c, d);
        BackwardConf(d, c, M_inv, round + 1, 0);
        return;
    }
    push_output(&outputs, c);
}

/* Rebuilds the initial 32-byte state from the 16-byte encrypted output */
void init_backward(const u8 input[ENCRYPTED_SIZE], const uint8_t M_inv[SIZE][SIZE], u8 trial_d_u8[SIZE], const int depth) {
    u8 trial_d[SIZE];
    memcpy(trial_d, trial_d_u8, SIZE);
    if (depth < ENCRYPTED_SIZE) {
        for (int i = 0; i < pair_lookup[input[depth]].count; i++) {
            trial_d[depth * 2]     = pair_lookup[input[depth]].pairs[i].first;
            trial_d[depth * 2 + 1] = pair_lookup[input[depth]].pairs[i].second;
            init_backward(input, M_inv, trial_d, depth + 1);
            if (outputs.count > MAX_OUTPUTS) {
                printf("Found %d solutions, stopping search.\n", MAX_OUTPUTS);
                return;
            }
        }
        return;
    }
    u8 candidate[SIZE], d[SIZE];
    memcpy(d, trial_d, SIZE);
    memcpy(candidate, d, SIZE);
    Backward(candidate, M_inv, 0);
}

/* --- Global target (encrypted output) --- */
u8 target[ENCRYPTED_SIZE] = "Hire me!!!!!!!!";

/* ---pack 32x32 matrix --- */
void pack_matrix(const uint8_t M[SIZE][SIZE], uint32_t M_packed[SIZE]) {
    for (int j = 0; j < SIZE; j++) {
        M_packed[j] = 0;
        for (int k = 0; k < SIZE; k++) {
            if (M[j][k])
                M_packed[j] |= (1U << k);
        }
    }
}

int main() {
    u8 trial_d[SIZE] = {0};

    /* Initialize outputs dynamic array */
    outputs.count = 0;
    outputs.capacity = MAX_OUTPUTS * 2;
    outputs.outputs = (Output*)malloc(outputs.capacity * sizeof(Output));
    if (!outputs.outputs) {
        fprintf(stderr, "Memory allocation error\n");
        return 1;
    }

    /* Initialize pair_lookup buckets */
    for (int i = 0; i < 256; i++) {
        pair_lookup[i].count = 0;
        pair_lookup[i].capacity = 16;
        pair_lookup[i].pairs = (Pair*)malloc(pair_lookup[i].capacity * sizeof(Pair));
        if (!pair_lookup[i].pairs) {
            fprintf(stderr, "Memory allocation error\n");
            return 1;
        }
    }

    /* Initialize inv_confusion_full buckets */
    for (int i = 0; i < 256; i++) {
        inv_confusion_full[i].count = 0;
        inv_confusion_full[i].capacity = 16;
        inv_confusion_full[i].values = (u8*)malloc(inv_confusion_full[i].capacity * sizeof(u8));
        if (!inv_confusion_full[i].values) {
            fprintf(stderr, "Memory allocation error\n");
            return 1;
        }
    }

    /* Build inverse confusion mapping and pair lookup table */
    initialize_inv_confusion_full();
    generate_pair_lookup();

    /* Construct and invert the diffusion matrix */
    uint8_t M[SIZE][SIZE] = {0};
    uint8_t M_inv[SIZE][SIZE] = {0};
    construct_matrix(M);
    if (!invert_matrix(M, M_inv)) {
        fprintf(stderr, "Error: Matrix is not invertible\n");
        return 1;
    }

    /* Pack the inverse matrix into bit-packed representation */
    pack_matrix(M_inv, M_inv_packed);

    /* Test the matrix inversion using the bit-packed multiplication */
    const u8 test_input[SIZE] = {
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
        0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F, 0x10,
        0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18,
        0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F, 0x20
    };
    u8 output_test[SIZE] = {0};
    u8 o2[SIZE] = {0};
    uint32_t M_packed[SIZE];
    pack_matrix(M, M_packed);
    matrix_vector_multiply_packed(M_inv_packed, test_input, output_test);
    matrix_vector_multiply_packed(M_packed, output_test, o2);
    if (memcmp(test_input, o2, SIZE) == 0)
        printf("Matrix inversion test passed!\n");
    else
        printf("Matrix inversion test failed!\n");

    const clock_t start_time = clock();
    init_backward(target, M_inv, trial_d, 0);
    const clock_t end_time = clock();

    printf("\nTotal solutions found: %d\n", outputs.count);
    const double time_taken = ((double)(end_time - start_time)) / CLOCKS_PER_SEC * 1000.0;
    printf("Time taken: %.0f ms\n", time_taken);
    if (outputs.count > 0)
        printf("Time per solution: %.6f ms\n", time_taken / outputs.count);

    /* Validate each output candidate by running the forward transformation. Sanity check */
    for (int i = 0; i < outputs.count; i++) {
        u8 c[SIZE], d[SIZE];
        memcpy(c, outputs.outputs[i].data, SIZE);
        Forward(c, d, confusion, diffusion);
        if (memcmp(d, target, ENCRYPTED_SIZE) != 0) {
            printf("Output is not valid!\n");
            return 1;
        }
        for (int j = i + 1; j < outputs.count; j++) {
            if (memcmp(outputs.outputs[i].data, outputs.outputs[j].data, SIZE) == 0) {
                printf("Duplicate output found!\n");
                return 1;
            }
        }
    }

    /* Free allocated memory */
    for (int i = 0; i < 256; i++) {
        free(pair_lookup[i].pairs);
    }
    for (int i = 0; i < 256; i++) {
        free(inv_confusion_full[i].values);
    }
    free(outputs.outputs);
    return 0;
}
