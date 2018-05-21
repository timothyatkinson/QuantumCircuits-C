#include "q_circuit.h"

/**g_state_alloc
  *Allocates a q_state struct.
    *qubits. The number of qubits of the q_state.
  Returns the empty state "state".
*/
q_state* q_state_alloc(int qubits){
  int rows = pow(2, qubits);
  q_state* state = malloc(sizeof(q_state));
  state->qubits = qubits;
  state->vector = gsl_matrix_complex_alloc(rows, 1);
  return state;
}

/**q_state_calloc
  *Allocates a q_state struct and initialises all entries to zero. Note that this is not a valid quantum state and must be changed.
    *qubits. The number of qubits of the q_state.
  Returns the generated state "state"
*/
q_state* q_state_calloc(int qubits){
  int rows = pow(2, qubits);
  q_state* state = q_state_alloc(qubits);
  gsl_matrix_complex_set_all(state->vector, GSL_COMPLEX_ZERO);
  return state;
}

/**q_state_free
  *Frees a given q_state.
    *state. The state to free.
*/
void q_state_free(q_state* state){
  gsl_matrix_complex_free(state->vector);
  free(state);
}

/**q_state_print
  *Prints a q_state.
    *state. The state to print.
*/
void q_state_print(q_state* state){

  printf("%d qubits:\n", state->qubits);
  for (int i = 0; i < state->vector->size1; i++) {
    for (int j = 0; j < state->vector->size2; j++) {
      gsl_complex q = gsl_matrix_complex_get(state->vector, i, j);
      double re = GSL_REAL(q);
      double im = GSL_IMAG(q);
      printf("(%lf, %lfj)\t", re, im);
    }
    printf("\n");
  }

}


/**q_op_alloc
  *Allocates a q_op struct.
    *qubits. The number of qubits that the q_op operates on.
  Returns the generated operator "op"
*/
q_op* q_op_alloc(int qubits){
  int rows = pow(2, qubits);
  q_op* op = malloc(sizeof(q_op));
  op->qubits = qubits;
  op->matrix = gsl_matrix_complex_alloc(rows, rows);
  return op;
}

/**q_op_calloc
  *Allocates a q_op struct and initialises all entries to zero. Note that this is not a valid quantum operator and must be changed.
    *qubits. The number of qubits that the q_op operates on.
  Returns the generated operator "op"
*/
q_op* q_op_calloc(int qubits){
  int rows = pow(2, qubits);
  q_op* op = q_op_alloc(qubits);
  gsl_matrix_complex_set_all(op->matrix, GSL_COMPLEX_ZERO);
  return op;
}

/**q_op_free
  *Frees a given q_op.
    *op. The state to op.
*/
void q_op_free(q_op* op){
  gsl_matrix_complex_free(op->matrix);
  free(op);
}

/**q_op_print
  *Prints a q_op.
    *op. The op to print.
*/
void q_op_print(q_op* op){

  printf("%d qubits:\n", op->qubits);
  for (int i = 0; i < op->matrix->size1; i++) {
    for (int j = 0; j < op->matrix->size2; j++) {
      gsl_complex q = gsl_matrix_complex_get(op->matrix, i, j);
      double re = GSL_REAL(q);
      double im = GSL_IMAG(q);
      printf("(%lf, %lfj)\t", re, im);
    }
    printf("\n");
  }

}

/**apply_qop
  *Applies a q_op operator to a state. The old state is not destroyed and must be freed by the user. This operation is effectively op * state using matrix multiplication.
    *op. The q_op to apply.
    *state. The state to apply the q_op to.
  Returns the produced state "new_state"
*/
q_state* apply_qop(q_op* op, q_state* state){
  if(op->qubits != state->qubits){
    printf("Error: size mismatch in operator application. Terminating.\n");
    exit(0);
  }
  q_state* new_state = q_state_alloc(state->qubits);

  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, op->matrix, state->vector, GSL_COMPLEX_ZERO, new_state->vector);

  return new_state;
}

/**q_state_tensor
  *Performs the matrix tensor operation on quantum states a and b. Neither a nor b is destroyed and both must be freed by the user.
    *a. The first q_state.
    *b. The second q_state.
  *Returns the tensor of a and b "new_state".
*/
q_state* q_state_tensor(q_state* a, q_state* b){
  q_state* new_state = q_state_calloc(a->qubits + b->qubits);
  for(int i = 0; i < a->vector->size1; i++){
    for(int j = 0; j < a->vector->size2; j++){
      gsl_complex q1 = gsl_matrix_complex_get(a->vector, i, j);
      if(GSL_REAL(q1) != 0.0 || GSL_IMAG(q1) != 0.0){
        for(int k = 0; k < b->vector->size1; k++){
          for(int l = 0; l < b->vector->size2; l++){
            gsl_complex q2 = gsl_matrix_complex_get(b->vector, k, l);
            if(GSL_REAL(q1) != 0.0 || GSL_IMAG(q1) != 0.0){
              gsl_matrix_complex_set(new_state->vector, (i * b->vector->size1) + k, (j * b->vector->size2) + l, gsl_complex_mul(q1, q2));
            }
          }
        }
      }
    }
  }
  return new_state;
}

/**q_op_tensor
  *Performs the matrix tensor operation on quantum operators a and b. Neither a nor b is destroyed and both must be freed by the user.
    *a. The first q_op.
    *b. The second q_op.
  *Returns the tensor of a and b "new_op".
*/
q_op* q_op_tensor(q_op* a, q_op* b){
  q_op* new_op = q_op_calloc(a->qubits + b->qubits);
  for(int i = 0; i < a->matrix->size1; i++){
    for(int j = 0; j < a->matrix->size2; j++){
      gsl_complex q1 = gsl_matrix_complex_get(a->matrix, i, j);
      if(GSL_REAL(q1) != 0.0 || GSL_IMAG(q1) != 0.0){
        for(int k = 0; k < b->matrix->size1; k++){
          for(int l = 0; l < b->matrix->size2; l++){
            gsl_complex q2 = gsl_matrix_complex_get(b->matrix, k, l);
            if(GSL_REAL(q1) != 0.0 || GSL_IMAG(q1) != 0.0){
              gsl_matrix_complex_set(new_op->matrix, (i * b->matrix->size1) + k, (j * b->matrix->size2) + l, gsl_complex_mul(q1, q2));
            }
          }
        }
      }
    }
  }
  return new_op;
}
