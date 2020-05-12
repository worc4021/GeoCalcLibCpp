#include "lrsinterface.hpp"
#include <iostream>

RowXq classify(size_t nCols, const lrs_mp_vector output) {
    RowXq retVal(nCols);
    if (zero(output[0])) {
        // ray 
        RowXz ray(nCols);
        for (auto i = 0; i < ray.cols(); i++)
            ray(i) = zType(output[i]);

        zType norm(ray.squaredNorm());

        for (auto i = 0; i < retVal.cols(); i++)
            retVal(i) = qType(zType(output[i]),norm);
    } else {
        // handle hyperplanes with negative sign
        qType signOfValue = (mpz_cmp_si(output[0],0L)+1) ? 1: -1;
        zType out0(output[0]);
        for (auto i = 0; i<retVal.cols(); i++)
            retVal(i) = signOfValue * qType(zType(output[i]), out0);
    }

    return retVal;
}


void declassify(const RowXq& row, lrs_mp_vector num, lrs_mp_vector den) {
    for (auto i = 0; i < row.cols(); i++) {
    
        mpz_set(num[i], row(i).get_num_mpz_t());
        mpz_set(den[i], row(i).get_den_mpz_t());
    }
}

void Polytope::vertexEnumerate(void) {
    assert(inHrep);

    lrs_dic *Pv;
    lrs_dat *Qv;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    lrs_mp_init (ZERO, stdin, stdout);
    
    Qv = lrs_alloc_dat ("LRS globals");
    assert( Qv!= NULL );

    Qv->m = HrepIneq.rows() + HrepEq.rows();
    Qv->n = HrepIneq.cols();
    Qv->getvolume = 1L;
    
    output = lrs_alloc_mp_vector (Qv->n);

    std::list<RowXq> bufferV;
    std::list<RowXq> bufferR;


    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(HrepIneq.cols());
    den = lrs_alloc_mp_vector(HrepIneq.cols());
    
    Pv = lrs_alloc_dic (Qv);
    assert( Pv != NULL );
    
    
    for (auto i = 1; i <= HrepIneq.rows(); ++i){
        declassify(HrepIneq.row(i-1), num, den);
        lrs_set_row_mp(Pv,Qv,i,num,den,GE);
    }

    for (auto i = 1; i <= HrepEq.rows(); ++i){
        declassify(HrepIneq.row(i-1), num, den);
        lrs_set_row_mp(Pv,Qv,i+HrepIneq.rows()-1,num,den,EQ);
    }

    if (lrs_getfirstbasis (&Pv, Qv, &Lin, TRUE)) {

        for (auto col = 0L; col < Qv->nredundcol; col++)
            lrs_printoutput (Qv, Lin[col]); 

        do
            {
                for (auto col = 0L; col <= Pv->d; col++)
                    if (lrs_getsolution (Pv, Qv, output, col)) {
                        mpz_cmp_si(output[0],0L)                                    ? 
                            bufferV.push_back( classify(HrepIneq.cols(), output) )  : 
                            bufferR.push_back( classify(HrepIneq.cols(), output) )  ;
                    }
            }
        while (lrs_getnextbasis (&Pv, Qv, FALSE));
    }


    rescalevolume (Pv, Qv, Qv->Nvolume, Qv->Dvolume);
    rattodouble (Qv->Nvolume, Qv->Dvolume, &volume);

    lrs_clear_mp_vector (output, Qv->n);
    lrs_clear_mp_vector (num, Qv->n);
    lrs_clear_mp_vector (den, Qv->n);
    lrs_free_dic (Pv,Qv);
    lrs_free_dat (Qv);

    VrepV = MatrixXq(bufferV.size(), HrepIneq.cols());
    VrepR = MatrixXq(bufferR.size(), HrepIneq.cols());

    auto rowCount  = 0;

    for (auto row : bufferV){
        setRow(VrepV, rowCount, row);
        rowCount++;
    }
    
    rowCount  = 0;
    for (auto row : bufferR){
        setRow(VrepR, rowCount, row);
        rowCount++;
    }

    inVrep = true;

}

void Polytope::facetEnumerate(void) {
    assert( inVrep );
    
    lrs_dic *P;
    lrs_dat *Q;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    lrs_mp_init (ZERO, stdin, stdout);
    Q = lrs_alloc_dat ("LRS globals");
    assert ( Q != NULL );
    Q->m = VrepV.rows() + VrepR.rows();
    Q->n = VrepV.cols();
    Q->hull = TRUE;
    Q->polytope = TRUE;

    if ( VrepR.rows()>0 ) {
        Q->polytope = FALSE;
        Q->homogeneous = FALSE;
    }

    output = lrs_alloc_mp_vector (Q->n);
    std::list<RowXq> buffer;
    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(VrepV.cols());
    den = lrs_alloc_mp_vector(VrepV.cols());

    zType myOne(1L);
    for (auto i = 0; i<VrepV.cols(); i++){
        mpz_set(num[i],myOne.get_mpz_t());
        mpz_set(den[i],myOne.get_mpz_t());
    }
    P = lrs_alloc_dic (Q);
    assert ( P != NULL );
    lrs_set_row_mp(P ,Q ,0 ,num ,den , GE);
    mpz_set( P->A[0][0], myOne.get_mpz_t() );

    for (auto i = 0; i < VrepV.rows(); i++) {
        declassify(VrepV.row(i), num, den);
        lrs_set_row_mp(P ,Q ,i+1 ,num ,den , GE);
    }

    for (auto i = 0; i < VrepR.rows(); i++) {
        declassify(VrepR.row(i), num, den);
        lrs_set_row_mp(P,Q, i+1+VrepV.rows(), num, den, GE);
    }

    if ( lrs_getfirstbasis (&P, Q, &Lin, TRUE) ){

        for (auto col = 0L; col < Q->nredundcol; col++)
        lrs_printoutput (Q, Lin[col]);

        do
        {
        for (auto col = 0; col <= P->d; col++)
            if (lrs_getsolution (P, Q, output, col))
                buffer.push_back(classify(VrepV.cols(), output));
        }
        while (lrs_getnextbasis (&P, Q, FALSE));
    }
    lrs_clear_mp_vector ( output, Q->n);
    lrs_clear_mp_vector ( num, Q->n);
    lrs_clear_mp_vector ( den, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );
    HrepIneq = MatrixXq(buffer.size(), VrepV.cols());
    HrepEq = MatrixXq(0,VrepV.cols());
    auto rowCount  = 0;
    for (auto row : buffer){
        setRow(HrepIneq, rowCount, row);
        rowCount++;
    }
    inHrep = true;

}