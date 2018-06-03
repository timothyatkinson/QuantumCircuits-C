#ifndef PD_Q_H
#define PD_Q_H

#include "q_circuit.h"

q_state* q_zero();
q_state* q_one();
q_state* q_rand();
q_op* q_identity(int qubits);
q_op* q_hadamard();
q_op* q_pauli_X();
q_op* q_pauli_Y();
q_op* q_pauli_Z();
q_op* q_cX();
q_op* q_cY();
q_op* q_cZ();
q_op* q_rot_z(double p);
q_op* q_crot_z(double p);
q_op* q_swap(int qubits, int* map);

/**rand_double - RAND DOUBLE
  *Computes a random double between 0 and 1 using C's inbuilt RNG
*/
double rand_double();
gsl_complex e_i_pi(double p);
#endif
