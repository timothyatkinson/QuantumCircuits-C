#ifndef PD_Q_H
#define PD_Q_H
#define TOFFOLI_QUBITS 3
#define FREDKIN_QUBITS 3
#define MARGOLUS_QUBITS 3

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

q_op* q_s();
q_op* q_t();
q_op* q_td();
q_op* q_ct();
q_op* r_x(double angle);
q_op* r_y(double angle);
q_op* r_z(double angle);
q_op* q_g();
q_op* q_gd();
q_op* q_v();
q_op* q_vd();
q_op* q_cvd();
q_op* q_toffoli();
q_op* q_fredkin();
q_op* q_margolus();
#endif
