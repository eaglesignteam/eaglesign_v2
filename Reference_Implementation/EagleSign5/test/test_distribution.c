/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "../params.h"
#include "../sign.h"
#include "../packing.h"
#include "../polymatrix.h"
#include "../polyvec.h"
#include "../poly.h"
#include "../randombytes.h"
#include "../symmetric.h"
#include "../fips202.h"

#define MLEN 100
#define NTESTS 1000

int main(void)
{
  uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES];
  uint8_t tr[SEEDBYTES], *rhoprime1, *rhoprime3;
  const uint8_t *rhoprime2;
  poly g, g_inv;
  polyvecl A[K], F[L], D[K], F_INV[L], Fg[L], DF[K];
  polyveck tmp2[L];
  char fn_req[35], fn_req2[35];
  FILE *fp_req, *fp_req2;
  int siz = ETAF * TG;
  int siz2 = L * ETAF * TD;
  float stat[2 * siz + 1];
  float mean = 0;
  float vari = 0;

  float stat2[2 * siz + 1];
  float mean2 = 0;
  float vari2 = 0;

  for (int i = 0; i <= 2 * siz; i++)
    stat[i] = 0;

  for (int i = 0; i <= 2 * siz2; i++)
    stat2[i] = 0;

  for (int ite = 0; ite < NTESTS; ite++)
  {
    /* Get randomness for rho and rhoprime */
    randombytes(seedbuf, SEEDBYTES);
    shake256(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES);
    rhoprime1 = seedbuf;
    rhoprime2 = seedbuf + SEEDBYTES;
    rhoprime3 = seedbuf + 2 * SEEDBYTES;

    /* Small and Sparse polynomials based matrix D*/
    polymatrix_k_l_expand_d(D, rhoprime2);

    /* Small and Sparse polynomial g*/
  retry_g:
    shake256(rhoprime1, SEEDBYTES, rhoprime1, SEEDBYTES);
    if (poly_expand_g_invertible(&g_inv, &g, rhoprime1))
      goto retry_g;

    /* Small and Uniform polynomials based invertible matrix F and its inverse F_INV */
  retry_F:
    shake256(rhoprime3, CRHBYTES, rhoprime3, CRHBYTES);
    if (polymatrix_l_expand_f_invertible(F_INV, F, rhoprime3))
      goto retry_F;

    /* Compute Fg = F*g */
    poly_pointwise_matrix_product_l_l(Fg, g, F);
    polymatrix_invntt_tomont_l_l(Fg);

    for (int i = 0; i < L; i++)
      for (int j = 0; j < L; j++)
        for (int k = 0; k < N; k++)
        {
          stat[Fg[i].vec[j].coeffs[k] + siz] += (float)1 / (NTESTS * N * L * L);
          mean += (float)Fg[i].vec[j].coeffs[k] / (NTESTS * N * L * L);
          vari += (float)Fg[i].vec[j].coeffs[k] * Fg[i].vec[j].coeffs[k] / (NTESTS * N * L * L);
        }

    /* Compute DF*/
    polymatrix_pointwise_product(tmp2, D, F);
    polymatrix_reformat(DF, tmp2);
    polymatrix_invntt_tomont_k_l(DF);

    for (int i = 0; i < K; i++)
      for (int j = 0; j < L; j++)
        for (int k = 0; k < N; k++)
        {
          stat2[DF[i].vec[j].coeffs[k] + siz2] += (float)1 / (NTESTS * N * K * L);
          mean2 += (float)DF[i].vec[j].coeffs[k] / (NTESTS * N * K * L);
          vari2 += (float)DF[i].vec[j].coeffs[k] * DF[i].vec[j].coeffs[k] / (NTESTS * N * K * L);
        }
  }

  // Fg
  vari -= mean * mean;

  sprintf(fn_req, "test/distributionFg_%.16s.csv", CRYPTO_ALGNAME);
  if ((fp_req = fopen(fn_req, "w")) == NULL)
  {
    printf("Couldn't open <%s> for write\n", fn_req);
    return -1;
  }

  fprintf(fp_req, "coef\tcount\n");
  fprintf(fp_req, "%.10lf\t%.10lf\n", mean, vari);
  for (int i = 0; i <= 2 * siz; i++)
    fprintf(fp_req, "%d\t%.10lf\n", i - siz, stat[i]);

  fclose(fp_req);

  // DF
  vari2 -= mean2 * mean2;

  sprintf(fn_req2, "test/distributionDF_%.16s.csv", CRYPTO_ALGNAME);
  if ((fp_req2 = fopen(fn_req2, "w")) == NULL)
  {
    printf("Couldn't open <%s> for write\n", fn_req2);
    return -1;
  }

  fprintf(fp_req2, "coef\tcount\n");
  fprintf(fp_req2, "%.10lf\t%.10lf\n", mean2, vari2);
  for (int i = 0; i <= 2 * siz2; i++)
    fprintf(fp_req2, "%d\t%.10lf\n", i - siz2, stat2[i]);

  fclose(fp_req2);
}
