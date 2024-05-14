#include "lrsinterface.hpp"
#include <iostream>
#include <algorithm>

lrs_restart_dat* lrs_alloc_restart()
{
  int i;

  lrs_restart_dat *R;

  R = (lrs_restart_dat *) malloc (sizeof (lrs_restart_dat));
  if (R == NULL)
    return R;  
  
  R->overide=0;     /* do not overide Q */
  R->restart=0;     /* do not do a restart */
  R->facet=NULL;    /* this will be allocated later when we know its size */
  R->d=0;
  R->maxcobases=0;
  R->maxdepth=-1;  /* will be set to MAXD in lrs*_main */
  R->mindepth=0;
  R->maxcobases=0;
  for(i=0;i<10;i++)
    R->count[i]=0;
  R->depth=0;
  R->lrs=1;
  R->redund=0;
  R->verifyredund=0;
  R->redineq = NULL;

  return R;
}


void toLRSNumeratorAndDenominator(std::size_t iRow, const Vrep& vrep, lrs_mp_vector num, lrs_mp_vector den) {
    // const zType one_t(1L);
    // const zType zero_t(0L);
    // if (vrep.isVertex[iRow]) {
    //     mpz_set(num[0], one_t.get_mpz_t());
    // } else {
    //     mpz_set(num[0], zero_t.get_mpz_t());
    // }
    // mpz_set(den[0], one_t.get_mpz_t());
    for (auto i = 0; i <= vrep.nDim; i++) {
        mpz_set(num[i], vrep.rows[iRow](i).get_num_mpz_t());
        mpz_set(den[i], vrep.rows[iRow](i).get_den_mpz_t());
    }
}

void toLRSNumeratorAndDenominator(std::size_t iRow, const Hrep& hrep, lrs_mp_vector num, lrs_mp_vector den) {
    qType curVal;
    mpz_set(num[0], hrep.rows[iRow](0).get_num_mpz_t());
    mpz_set(den[0], hrep.rows[iRow](0).get_den_mpz_t());
    for (auto i = 1; i < hrep.rows[iRow].cols(); i++) {
        /*
        In LRS polyhedra are given as bl+Al*x>=0 hence the regular A*x<=b has to be transformed 
        by Al = -A.
        */
        curVal = -hrep.rows[iRow](i);
        mpz_set(num[i], curVal.get_num_mpz_t());
        mpz_set(den[i], curVal.get_den_mpz_t());
    }
}


void Polytope::vertexEnumerate(void) {
    
    lrs_dic *Pv;
    lrs_dat *Qv;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    lrs_mp_init (ZERO, stdin, stdout);
    
    Qv = lrs_alloc_dat ("LRS globals");
    assert( Qv!= NULL );

    Qv->m = hRep.rows.size();
    Qv->n = hRep.nDim + 1;
    Qv->getvolume = 1L;
    
    output = lrs_alloc_mp_vector (Qv->n);

    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(Qv->n);
    den = lrs_alloc_mp_vector(Qv->n);
    
    Pv = lrs_alloc_dic (Qv);
    assert( Pv != NULL );
    

    for (auto i = 1; i <= hRep.rows.size(); ++i){
        toLRSNumeratorAndDenominator(i-1, hRep, num, den);
        auto type = std::find(hRep.linearities.begin(), hRep.linearities.end(), i-1) != hRep.linearities.end() ? EQ : GE;
        lrs_set_row_mp(Pv,Qv,i,num,den,type);
    }

    if (lrs_getfirstbasis(&Pv, Qv, &Lin, TRUE))
    {

        for (auto col = 0L; col < Qv->nredundcol; col++)
        {
            vRep.push_back(Lin[col], true);
        }

        do
        {
            for (auto col = 0L; col <= Pv->d; col++)
            {
                if (lrs_getsolution(Pv, Qv, output, col))
                {
                    vRep.push_back(output, false);
                }
            }
        } while (lrs_getnextbasis(&Pv, Qv, FALSE));
    }

    rescalevolume (Pv, Qv, Qv->Nvolume, Qv->Dvolume);
    rattodouble (Qv->Nvolume, Qv->Dvolume, &volume);

    lrs_clear_mp_vector (output, Qv->n);
    lrs_clear_mp_vector (num, Qv->n);
    lrs_clear_mp_vector (den, Qv->n);
    lrs_free_dic (Pv,Qv);
    lrs_free_dat (Qv);

}

void Polytope::facetEnumerate(void) {
    
    lrs_dic *P;
    lrs_dat *Q;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    lrs_mp_init (ZERO, stdin, stdout);
    Q = lrs_alloc_dat ("LRS globals");
    assert ( Q != NULL );
    Q->m = vRep.rows.size();
    Q->n = vRep.nDim + 1;
    Q->hull = TRUE;
    Q->polytope = TRUE;

    if ( std::find(vRep.isVertex.begin(), vRep.isVertex.end(), false) != vRep.isVertex.end() ){
        Q->polytope = FALSE;
        Q->homogeneous = FALSE;
    }

    output = lrs_alloc_mp_vector (Q->n); // b, -a1, -a2, -a3, ... , -an
    
    lrs_mp_vector num, den; // same size as output
    num = lrs_alloc_mp_vector(Q->n);
    den = lrs_alloc_mp_vector(Q->n);

    zType myOne(1L);
    for (auto i = 0; i<Q->n; i++){
        mpz_set(num[i],myOne.get_mpz_t());
        mpz_set(den[i],myOne.get_mpz_t());
    }
    P = lrs_alloc_dic (Q);
    assert ( P != NULL );
    lrs_set_row_mp(P ,Q ,0 ,num ,den , GE);
    mpz_set( P->A[0][0], myOne.get_mpz_t() );

    
    for (auto i = 0; i < vRep.rows.size(); i++) {
        toLRSNumeratorAndDenominator(i, vRep, num, den);
        auto type = std::find(vRep.linearities.begin(), vRep.linearities.end(), i+1) != vRep.linearities.end() ? EQ : GE;
        lrs_set_row_mp(P ,Q ,i+1 ,num ,den , type);
    }

    if ( lrs_getfirstbasis (&P, Q, &Lin, TRUE) ){

        for (auto col = 0L; col < Q->nredundcol; col++)
        {
            hRep.push_back(Lin[col], true);
        }

        do
        {
        for (auto col = 0; col <= P->d; col++)
            if (lrs_getsolution (P, Q, output, col))
                hRep.push_back(output, false);
        }
        while (lrs_getnextbasis (&P, Q, FALSE));
    }
    lrs_clear_mp_vector ( output, Q->n);
    lrs_clear_mp_vector ( num, Q->n);
    lrs_clear_mp_vector ( den, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );
    
}