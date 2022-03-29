#include "predefined_q.h"


/**rand_double - RAND double
  *Computes a random double between 0 and 1 using C's inbuilt RNG
*/
double rand_double(){
  return (double)rand() / (double)RAND_MAX;
}

gsl_complex e_i_pi(double p){
  gsl_complex ret;
  GSL_SET_COMPLEX(&ret, cos(p), sin(p));
  return ret;
}

void int_to_bin_digit(unsigned int in, int count, int* out);
//Lazy Source: https://stackoverflow.com/questions/31577866/c-programming-convert-integer-to-binary-array
void int_to_bin_digit(unsigned int in, int count, int* out)
{
    /* assert: count <= sizeof(int)*CHAR_BIT */
    unsigned int mask = 1U << (count-1);
    int i;
    for (i = 0; i < count; i++) {
        out[i] = (in & mask) ? 1 : 0;
        in <<= 1;
    }
}

int bin_to_int_digit(int count, int* in);
int bin_to_int_digit(int count, int* in){
  int val = 0;

  for(int i = 0; i < count; i++){
    if(in[i] == 1){
      val += pow(2, count - (1 + i));
    }
  }

  return val;
}

q_state* q_zero(){
  q_state* q = q_state_calloc(1);
  gsl_matrix_complex_set(q->vector, 0, 0, GSL_COMPLEX_ONE);
  return q;
}

q_state* q_one(){
  q_state* q = q_state_calloc(1);
  gsl_matrix_complex_set(q->vector, 1, 0, GSL_COMPLEX_ONE);
  return q;
}

q_state* q_rand(){
  q_state* q = q_state_calloc(1);
  double r1 = rand_double() * 2 * M_PI;
  double r2 = rand_double() * 2 * M_PI;
  gsl_complex z;
  gsl_complex o;
  GSL_SET_COMPLEX(&z, cos(r1) * cos(r2), sin(r1) * cos(r2));
  GSL_SET_COMPLEX(&o, cos(r1) * sin(r2), sin(r1) * sin(r2));
  gsl_matrix_complex_set(q->vector, 0, 0, z);
  gsl_matrix_complex_set(q->vector, 1, 0, o);
  return q;
}

q_op* q_identity(int qubits){
  q_op* op = q_op_calloc(qubits);
  for(int i = 0; i < op->matrix->size1; i++){
    gsl_matrix_complex_set(op->matrix, i, i, GSL_COMPLEX_ONE);
  }
  return op;
}

q_op* q_hadamard(){
  q_op* op = q_op_alloc(1);
  gsl_complex a;
  GSL_SET_COMPLEX(&a, 1.0 / sqrt(2), 0.0);
  gsl_complex b;
  GSL_SET_COMPLEX(&b, -1.0 / sqrt(2), 0.0);
  gsl_matrix_complex_set(op->matrix, 0, 0, a);
  gsl_matrix_complex_set(op->matrix, 1, 0, a);
  gsl_matrix_complex_set(op->matrix, 0, 1, a);
  gsl_matrix_complex_set(op->matrix, 1, 1, b);
  return op;
}

q_op* q_pauli_X(){
  q_op* op = q_op_calloc(1);
  gsl_matrix_complex_set(op->matrix, 0, 1, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 1, 0, GSL_COMPLEX_ONE);
  return op;
}

q_op* q_pauli_Y(){
  q_op* op = q_op_calloc(1);
  gsl_complex i;
  GSL_SET_COMPLEX(&i, 0.0, 1.0);
  gsl_complex m_i;
  GSL_SET_COMPLEX(&m_i, 0.0, -1.0);
  gsl_matrix_complex_set(op->matrix, 0, 1, m_i);
  gsl_matrix_complex_set(op->matrix, 1, 0, i);
  return op;
}

q_op* q_pauli_Z(){
  q_op* op = q_op_calloc(1);
  gsl_complex m_o;
  GSL_SET_COMPLEX(&m_o, -1.0, 0.0);
  gsl_matrix_complex_set(op->matrix, 0, 0, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 1, 1, m_o);
  return op;
}

q_op* q_cX(){
  q_op* op = q_op_calloc(2);
  gsl_matrix_complex_set(op->matrix, 0, 0, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 1, 1, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 2, 3, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 3, 2, GSL_COMPLEX_ONE);
  return op;
}

q_op* q_cY(){
  q_op* op = q_op_calloc(2);
  gsl_complex i;
  GSL_SET_COMPLEX(&i, 0.0, 1.0);
  gsl_complex m_i;
  GSL_SET_COMPLEX(&m_i, 0.0, -1.0);
  gsl_matrix_complex_set(op->matrix, 0, 0, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 1, 1, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 2, 3, m_i);
  gsl_matrix_complex_set(op->matrix, 3, 2, i);
  return op;
}

q_op* q_cZ(){
  q_op* op = q_op_calloc(2);
  gsl_complex m_o;
  GSL_SET_COMPLEX(&m_o, -1.0, 0.0);
  gsl_matrix_complex_set(op->matrix, 0, 0, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 1, 1, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 2, 2, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 3, 3, m_o);
  return op;
}

q_op* q_rot_z(double p){
  q_op* op = q_op_calloc(1);
  gsl_matrix_complex_set(op->matrix, 0, 0, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 1, 1, e_i_pi(p * 2.0 * M_PI));
  return op;
}

q_op* q_crot_z(double p){
  q_op* op = q_op_calloc(2);
  gsl_matrix_complex_set(op->matrix, 0, 0, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 1, 1, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 2, 2, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set(op->matrix, 3, 3, e_i_pi(p *  2.0 * M_PI));
  return op;
}

q_op* q_swap(int qubits, int* map){
  q_op* op = q_op_calloc(qubits);
  for(int i = 0; i < (int)pow(2, qubits); i++){
    int bits[qubits];
    int_to_bin_digit(i, qubits, bits);
    int new_bits[qubits];
    for(int j = 0; j < qubits; j++){
      new_bits[j] = bits[map[j]];
    }
    int index = bin_to_int_digit(qubits, new_bits);
    gsl_matrix_complex_set(op->matrix, i, index, GSL_COMPLEX_ONE);
  }
  return op;
}

q_op* q_s(){
  q_op* op = q_op_calloc(1);
  gsl_complex o;
  gsl_complex i;
  GSL_SET_COMPLEX(&o, 1.0, 0.0);
  GSL_SET_COMPLEX(&i, 0.0, 1.0);
  gsl_matrix_complex_set(op->matrix, 0, 0, o);
  gsl_matrix_complex_set(op->matrix, 1, 1, i);
  return op;
}

q_op* q_t(){
  q_op* op = q_op_calloc(1);
  gsl_complex o;
  gsl_complex i;
  GSL_SET_COMPLEX(&o, 1.0, 0.0);
  GSL_SET_COMPLEX(&i, 1.0 / sqrt(2), 1.0 / sqrt(2));
  gsl_matrix_complex_set(op->matrix, 0, 0, o);
  gsl_matrix_complex_set(op->matrix, 1, 1, i);
  return op;
}

q_op* q_td() {
    q_op* op = q_op_calloc(1);
    gsl_complex x;
    gsl_complex y;
    GSL_SET_COMPLEX(&x, 1.0, 0.0);
    GSL_SET_COMPLEX(&y, 1.0 / sqrt(2), -1.0 / sqrt(2));
    gsl_matrix_complex_set(op->matrix, 0, 0, x);
    gsl_matrix_complex_set(op->matrix, 1, 1, y);
    return op;
}

q_op* q_ct(){
  q_op* op = q_op_calloc(2);
  gsl_complex o;
  gsl_complex i;
  GSL_SET_COMPLEX(&o, 1.0, 0.0);
  GSL_SET_COMPLEX(&i, 1.0 / sqrt(2), 1.0 / sqrt(2));
  gsl_matrix_complex_set(op->matrix, 0, 0, o);
  gsl_matrix_complex_set(op->matrix, 1,1, o);
  gsl_matrix_complex_set(op->matrix, 2,2, o);
  gsl_matrix_complex_set(op->matrix, 3,3, i);
  return op;
}

q_op* r_x(double angle){
  q_op* op = q_op_calloc(1);
  gsl_complex a;
  gsl_complex b;
  gsl_complex c;
  gsl_complex d;
  GSL_SET_COMPLEX(&a, cos(angle/2.0), 0.0);
  GSL_SET_COMPLEX(&b, 0.0, -sin(angle/2.0));
  GSL_SET_COMPLEX(&c, 0.0, -sin(angle/2.0));
  GSL_SET_COMPLEX(&d, cos(angle/2.0), 0.0);
  gsl_matrix_complex_set(op->matrix, 0, 0, a);
  gsl_matrix_complex_set(op->matrix, 0, 1, b);
  gsl_matrix_complex_set(op->matrix, 1, 0, c);
  gsl_matrix_complex_set(op->matrix, 1, 1, d);
  return op;
}

q_op* r_y(double angle){
  q_op* op = q_op_calloc(1);
  gsl_complex a;
  gsl_complex b;
  gsl_complex c;
  gsl_complex d;
  GSL_SET_COMPLEX(&a, cos(angle/2.0), 0.0);
  GSL_SET_COMPLEX(&b, -sin(angle/2.0), 0.0);
  GSL_SET_COMPLEX(&c, sin(angle/2.0), 0.0);
  GSL_SET_COMPLEX(&d, cos(angle/2.0), 0.0);
  gsl_matrix_complex_set(op->matrix, 0, 0, a);
  gsl_matrix_complex_set(op->matrix, 0, 1, b);
  gsl_matrix_complex_set(op->matrix, 1, 0, c);
  gsl_matrix_complex_set(op->matrix, 1, 1, d);
  return op;
}

q_op* r_z(double angle){
  q_op* op = q_op_calloc(1);
  gsl_matrix_complex_set(op->matrix, 0, 0, e_i_pi(-angle/2.0));
  gsl_matrix_complex_set(op->matrix, 1, 1, e_i_pi(angle/2.0));
  return op;
}

static void _q_g_util(gsl_complex* x, gsl_complex* y, gsl_complex* z) {
    GSL_SET_COMPLEX(x, cos(M_PI/8), 0.0);
    GSL_SET_COMPLEX(y, sin(M_PI/8), 0.0);
    GSL_SET_COMPLEX(z, -sin(M_PI/8), 0.0);
    return;
}

q_op* q_g() {
    q_op* op = q_op_calloc(1);
    gsl_complex x, y, z;
    _q_g_util(&x, &y, &z);
    gsl_matrix_complex_set(op->matrix, 0, 0, x);
    gsl_matrix_complex_set(op->matrix, 0, 1, z);
    gsl_matrix_complex_set(op->matrix, 1, 0, y);
    gsl_matrix_complex_set(op->matrix, 1, 1, x);
    return op;
}

q_op* q_gd() {
    q_op* op = q_op_calloc(1);
    gsl_complex x, y, z;
    _q_g_util(&x, &y, &z);
    gsl_matrix_complex_set(op->matrix, 0, 0, x);
    gsl_matrix_complex_set(op->matrix, 0, 1, y);
    gsl_matrix_complex_set(op->matrix, 1, 0, z);
    gsl_matrix_complex_set(op->matrix, 1, 1, x);
    return op;
}

static void _q_v_util(gsl_complex* x, gsl_complex* y) {
    GSL_SET_COMPLEX(x, 1.0/2.0, 1.0/2.0);
    GSL_SET_COMPLEX(y, 1.0/2.0, -1.0/2.0);
    return;
}

q_op* q_v() {
    q_op* op = q_op_calloc(1);
    gsl_complex x, y;
    _q_v_util(&x, &y);
    gsl_matrix_complex_set(op->matrix, 0, 0, x);
    gsl_matrix_complex_set(op->matrix, 0, 1, y);
    gsl_matrix_complex_set(op->matrix, 1, 0, y);
    gsl_matrix_complex_set(op->matrix, 1, 1, x);
    return op;
}

q_op* q_vd() {
    q_op* op = q_op_calloc(1);
    gsl_complex x, y;
    _q_v_util(&x, &y);
    gsl_matrix_complex_set(op->matrix, 0, 0, y);
    gsl_matrix_complex_set(op->matrix, 0, 1, x);
    gsl_matrix_complex_set(op->matrix, 1, 0, x);
    gsl_matrix_complex_set(op->matrix, 1, 1, y);
    return op;
}

q_op* q_toffoli() {
    q_op* op = q_op_calloc(TOFFOLI_QUBITS);
    int n2 = pow(2, TOFFOLI_QUBITS);
    for (int i = 0; i < n2; i++) {
        if (i == 3)
            gsl_matrix_complex_set(op->matrix, i, n2 - 1, GSL_COMPLEX_ONE);
        else if (i == 7)
            gsl_matrix_complex_set(op->matrix, i, 3, GSL_COMPLEX_ONE);
        else
            gsl_matrix_complex_set(op->matrix, i, i, GSL_COMPLEX_ONE);
    }
    return op;
}
