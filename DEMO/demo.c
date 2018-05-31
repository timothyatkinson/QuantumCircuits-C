#include "../q_circuit.h"
#include "../predefined_q.h"
#include "time.h"

int main (void)
{
  srand(time(NULL));
  q_op* id = q_identity(2);
  q_op_print(id);

  q_state* z = q_rand();
  q_state_print(z);
  q_state* o = q_one();
  q_state_print(o);
  q_state* s = q_state_tensor(z, o);
  q_state_print(s);
  q_state* s2 = apply_qop(id, s);
  q_state_print(s2);

  q_op* h = q_hadamard();
  q_op* w = q_identity(1);
  q_op* wh = q_op_tensor(w, h);
  q_op* hw = q_op_tensor(h, w);
  q_state* s3 = apply_qop(wh, s2);
  q_op_print(h);
  q_op_print(wh);
  q_op_print(hw);
  q_op_free(h);
  q_op_free(w);
  q_op_free(wh);
  q_op_free(hw);
  q_state_print(s3);
  q_state_free(s3);

  q_op_free(id);
  q_state_free(z);
  q_state_free(o);
  q_state_free(s);
  q_state_free(s2);

  q_op* x = q_pauli_X();
  printf("X\n");
  q_op_print(x);
  q_op_free(x);

  q_op* y = q_pauli_Y();
  printf("Y\n");
  q_op_print(y);
  q_op_free(y);

  q_op* pz = q_pauli_Z();
  printf("Z\n");
  q_op_print(pz);
  q_op_free(pz);

  q_op* cx = q_cX();
  printf("cX\n");
  q_op_print(cx);
  q_op_free(cx);

  q_op* cy = q_cY();
  printf("cY\n");
  q_op_print(cy);
  q_op_free(cy);

  q_op* cz = q_cZ();
  printf("cZ\n");
  q_op_print(cz);
  q_op_free(cz);

  q_op* r0 = q_crot_z(0.0);
  printf("cR0\n");
  q_op_print(r0);
  q_op_free(r0);

  q_op* r1 = q_crot_z(1.0);
  printf("cR1\n");
  q_op_print(r1);
  q_op_free(r1);

  q_op* r05 = q_crot_z(0.5);
  printf("cR05\n");
  q_op_print(r05);
  q_op_free(r05);

  q_op* r025 = q_crot_z(0.25);
  printf("cR025\n");
  q_op_print(r025);
  q_op_free(r025);

  int map[3];
  map[0] = 2;
  map[1] = 0;
  map[2] = 1;
  q_op* swap = q_swap(3, map);
  printf("Swap3\n");
  q_op_print(swap);
  q_op_free(swap);

    q_state* z2 = q_rand();
    q_state_print(z2);
      q_state* z3 = q_rand();
      q_state_print(z3);
    q_state* o2 = q_one();
  q_state* z4 = q_zero();
    q_state_print(o2);
    s2 = q_state_tensor(z2, o2);
    q_state_print(s2);
    s3 = q_state_tensor(z3, o2);
    q_state_print(s3);
    printf("s2, s3: %lf\n", fidelity(s2, s3));
    printf("s2, s3: %lf\n", fidelity(s3, s2));
    printf("s2, s2: %lf\n", fidelity(s2, s2));
    printf("s3, s3: %lf\n", fidelity(s3, s3));
    q_state* s4 = q_state_tensor(o2, o2);
    q_state* s5 = q_state_tensor(o2, z4);
    printf("s4, s5: %lf\n", fidelity(s4, s5));


}
