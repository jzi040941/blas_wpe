//#include <iostream>
//#include <cstring>
//#include <cmath>
//#include <stdlib>
#include "mother.h"
int Nch = 2;
int tap = 2;
int Nfreq = 4;

CMAT* mat;
CMAT* sub;
CMAT *g, *AA, *phi, *identity, *X_delay;
MAT *idntest, *gamma_wpe;
int pregamma;
void ident_3dcmat(CMAT *mat) {
  ITER i, j, k;
  UINT d0, d2;

  ASSERT(mat->d0 == mat->d1, "This function requires SQUARE MATRIX.\n")
  if (mat->ndim == 1) {
    d0 = mat->d0;

#pragma omp parallel for schedule(dynamic,CHUNK_SIZE) shared(mat) private(i, j)
    for (i = 0; i < d0; i++) {
      for (j = 0; j < d0; j++) {
        if (i == j) {
          mat->data[i * d0 + j].re = 1;
          mat->data[i * d0 + j].im = 0;
        } else {
          mat->data[i * d0 + j].re = 0;
          mat->data[i * d0 + j].im = 0;
        }
      }
    }
  } else if (mat->ndim == 2) {
    d0 = mat->d0;
    d2 = mat->d2;
#pragma omp parallel for schedule(dynamic,CHUNK_SIZE) shared(mat) private(i, j, k)
    for (i = 0; i < d2; i++) {
      for (j = 0; j < d0; j++) {
        for (k = 0; k < d0; k++) {
          if (i == j) {
            mat->data[i + j * d0 + k * d0 * d0].re = 1;
            mat->data[i + j * d0 + k * d0 * d0].im = 0;
          } else {
            mat->data[i + j * d0 + k * d0 * d0].re = 0;
            mat->data[i + j * d0 + k * d0 * d0].im = 0;
          }
        }
      }
    }
  }
}
int main(){
    int NT = Nch*tap;
    pregamma = 0.1;
    //mat = czeros(Nch,NT,Nfreq);
    //sub =  czeros(Nch,NT);
    //csubmat(mat,sub,0,Nch,0,NT,0,1);
    //idntest = zeros(NT,NT);
    //ident_mat(idntest);
    //print_mat(idntest);
    //print_cmat(sub);

    //g,공분산 행렬 phi 초기화
    g = czeros(Nch,NT,Nfreq);
    AA = czeros(NT,NT,Nfreq);
    identity = czeros(NT,NT,Nfreq);
    ident_3dcmat(identity);
    phi = czeros(NT,NT,Nfreq);
    CTYPE epsilon= {2.220446049250313 * 1E-16,0.0}; 
    axpy_cmat(epsilon, identity, phi);
    //print_cmat(identity);
    
    //X_delay 초기화
    X_delay = czeros(NT, 1, Nfreq);
    //gamma 초기화
    gamma_wpe = alloc_mat(Nfreq);
    fill(gamma_wpe, pregamma);
    return 0;
}
