/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#ifndef PACKING_H
#define PACKING_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"

#define pack_pk EAGLESIGN_NAMESPACE(pack_pk)
void pack_pk(uint8_t pk[CRYPTO_EAGLESIGN_PUBLICKEYBYTES],
             const uint8_t rho[SEEDBYTES],
             const polyvecl E[K]);

#define unpack_pk EAGLESIGN_NAMESPACE(unpack_pk)
void unpack_pk(
    uint8_t rho[SEEDBYTES],
    polyvecl E[K],
    const uint8_t pk[CRYPTO_EAGLESIGN_PUBLICKEYBYTES]);

#define pack_sk EAGLESIGN_NAMESPACE(pack_sk)
void pack_sk(uint8_t sk[CRYPTO_EAGLESIGN_SECRETKEYBYTES],
             const uint8_t rho[SEEDBYTES],
             const uint8_t tr[SEEDBYTES],
             const poly g,
             const polyvecl D[K],
             const polyvecl F[L],
             const polyvecl Es[K]);

#define unpack_sk EAGLESIGN_NAMESPACE(unpack_sk)
void unpack_sk(
    uint8_t rho[SEEDBYTES],
    uint8_t tr[SEEDBYTES],
    poly *g,
    polyvecl D[K],
    polyvecl F[L],
    polyvecl Es[K],
    const uint8_t sk[CRYPTO_EAGLESIGN_SECRETKEYBYTES]);

#define pack_sig EAGLESIGN_NAMESPACE(pack_sig)
void pack_sig(uint8_t sig[CRYPTO_EAGLESIGN_BYTES],
              const uint8_t r[SEEDBYTES],
              const polyvecl *z);
#define unpack_sig EAGLESIGN_NAMESPACE(unpack_sig)
int unpack_sig(
    uint8_t r[SEEDBYTES],
    polyvecl *z,
    const uint8_t sig[CRYPTO_EAGLESIGN_BYTES]);

#endif
