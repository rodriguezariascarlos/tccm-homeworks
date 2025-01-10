#include <trexio.h> // Library for trexio functions
#include <stdio.h> // Library for printing
#include <stdlib.h> // Library for array memory allocation

// Function to allocate a 4D array dynamically for storing two-electron integrals
static double ****allocate_array(int dim) {
    double ****array = malloc(dim * sizeof(double ***)); // Allocate 3D pointers
    for (int i = 0; i < dim; i++) {
        array[i] = malloc(dim * sizeof(double **));  // Allocate 2D pointers
        for (int j = 0; j < dim; j++) {
            array[i][j] = malloc(dim * sizeof(double *)); // Allocate 1D pointers
            for (int k = 0; k < dim; k++) {
                array[i][j][k] = calloc(dim, sizeof(double)); // Initialize to zero
            }
        }
    }
    return array;
}

// Function to free a dynamically allocated 4D array
static void free_array(double ****array, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                free(array[i][j][k]); // Free 1D arrays
            }
            free(array[i][j]); // Free 2D pointers
        }
        free(array[i]); // Free 3D pointers
    }
    free(array); // Free the main array
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <molecule_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Construct file path using "../tests" directory
    char filepath[256];
    snprintf(filepath, sizeof(filepath), "../tests/%s", argv[1]); // Adding the tests folder

    // Open TREXIO file for reading
    trexio_exit_code rc;
    trexio_t *file = trexio_open(filepath, 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error opening file: %s\n", trexio_string_of_error(rc));
        return EXIT_FAILURE;
    }

    // Variables to store data from the TREXIO file
    double nucleus_repulsion;
    int n_up, mo_count;
    int64_t integral_count;

    // Read nuclear repulsion energy from file and print
    trexio_read_nucleus_repulsion(file, &nucleus_repulsion);
    printf("Nuclear Repulsion Energy: %f\n", nucleus_repulsion);

    // Read the number of occupied molecular orbitals (spin-up electrons)
    trexio_read_electron_up_num(file, &n_up);
    printf("Number of Occupied Orbitals: %d\n", n_up);

    // Read the total number of molecular orbitals
    trexio_read_mo_num(file, &mo_count);

    // Allocate memory for the one-electron integrals (core Hamiltonian)
    double *one_e_integrals = malloc(mo_count * mo_count * sizeof(double));
    trexio_read_mo_1e_int_core_hamiltonian(file, one_e_integrals);
    printf("One-electron integrals read successfully.\n");

    // Read two-electron integrals (electron repulsion integrals)
    trexio_read_mo_2e_int_eri_size(file, &integral_count);
    int32_t *indices = malloc(4 * integral_count * sizeof(int32_t));
    double *values = malloc(integral_count * sizeof(double));
    trexio_read_mo_2e_int_eri(file, 0, &integral_count, indices, values);
    printf("Two-electron integrals read successfully.\n");

    // Read molecular orbital energies
    double *orbital_energies = malloc(mo_count * sizeof(double));
    trexio_read_mo_energy(file, orbital_energies);

    // Allocate a 4D array to store the two-electron integrals
    double ****two_e_integrals = allocate_array(mo_count);
    for (int64_t n = 0; n < integral_count; n++) {
        int i = indices[4 * n], j = indices[4 * n + 1];
        int k = indices[4 * n + 2], l = indices[4 * n + 3];
        double value = values[n];

        // Fill all symmetric permutations of the two-electron integral
        two_e_integrals[i][j][k][l] = value;
        two_e_integrals[k][l][i][j] = value;
        two_e_integrals[i][l][k][j] = value;
        two_e_integrals[k][j][i][l] = value;
        two_e_integrals[j][i][l][k] = value;
        two_e_integrals[l][k][j][i] = value;
        two_e_integrals[j][k][l][i] = value;
        two_e_integrals[l][i][j][k] = value;
    }

    // Compute the Hartree-Fock energy
    double one_e_sum = 0.0, two_e_sum = 0.0;

    // Calculate the one-electron term (contributions from occupied orbitals)
    for (int i = 0; i < n_up; i++) {
        one_e_sum += one_e_integrals[i * mo_count + i];
    }
    printf("One-electron contribution: %f\n", 2.0 * one_e_sum);

    // Calculate the two-electron term using occupied orbitals
    for (int i = 0; i < n_up; i++) {
        for (int j = 0; j < n_up; j++) {
            two_e_sum += 2.0 * two_e_integrals[i][j][i][j] - two_e_integrals[i][j][j][i];
        }
    }
    printf("Two-electron contribution: %f\n", two_e_sum);

    // Final Hartree-Fock energy
    double hf_energy = nucleus_repulsion + 2.0 * one_e_sum + two_e_sum;
    printf("Hartree-Fock Energy: %f\n", hf_energy);

    // Compute MP2 energy (correlation energy)
    double mp2_energy = 0.0;

    // Loop over occupied and virtual orbitals for MP2 calculation
    for (int i = 0; i < n_up; i++) {
        for (int j = 0; j < n_up; j++) {
            for (int a = n_up; a < mo_count; a++) {
                for (int b = n_up; b < mo_count; b++) {
                    double value = two_e_integrals[i][j][a][b];

                    // Compute energy denominator for MP2 correlation energy
                    double denom = orbital_energies[i] + orbital_energies[j] - orbital_energies[a] - orbital_energies[b];
                    
                    // Avoid division by zero by checking the denominator
                    if (denom != 0) {
                        // Accumulate MP2 contribution (make sure to divide by the denominator)
                        mp2_energy += (value * value) / denom;
                    }
                }
            }
        }
    }
    printf("MP2 Energy: %f\n", mp2_energy);

    // Free allocated memory to avoid memory leaks
    free(one_e_integrals);
    free(orbital_energies);
    free(indices);
    free(values);
    free_array(two_e_integrals, mo_count);

    // Close TREXIO file
    trexio_close(file);

    return EXIT_SUCCESS;
}

