/***/
/* shufDQ.M: shuffling SOFAST-H(Z/D)QC input to extract DQ spectrum 
/*
/* Example: | nmrPipe -fn MAC -all -noRd -noWr -macro shufDQ.M \
/***/

/***/
/* Initialization steps:
/*  ARGS: extract summation or difference mode.
/*  INIT: update header to reflect reduced number of 1D vectors.
/***/

if (sliceCode == CODE_ARGS)
   {
   };

if (sliceCode == CODE_INIT)
   {
    tdSize   = getParm( fdata, NDTDSIZE, CUR_YDIM );
    apodSize = getParm( fdata, NDAPOD,   CUR_YDIM );

    (void) setParm( fdata, NDSIZE,   integer( specnum/2 ),  CUR_YDIM ); 
    (void) setParm( fdata, NDTDSIZE, integer( tdSize/2 ),   CUR_YDIM ); 
    (void) setParm( fdata, NDAPOD,   integer( apodSize/2 ), CUR_YDIM ); 
   };

if (sliceCode < 1)
   {
    exit( 0 ); 
   };

/***/
/* Process only every group of 4 complex 1D slices:
/***/

if (sliceCode % 4)
   {
    exit( 0 );
   };

float vA[size],  vB[size],  vC[size],  vD[size],  vE[size],  vF[size],  vG[size],  vH[size];
float vAA[size], vBB[size], vCC[size], vDD[size];

(void) dReadB( inUnit, vA, wordLen*size );
(void) dReadB( inUnit, vB, wordLen*size );
(void) dReadB( inUnit, vC, wordLen*size );
(void) dReadB( inUnit, vD, wordLen*size );
(void) dReadB( inUnit, vE, wordLen*size );
(void) dReadB( inUnit, vF, wordLen*size );
(void) dReadB( inUnit, vG, wordLen*size );
(void) dReadB( inUnit, vH, wordLen*size );

/***/
/* Desired output:
/*    cos component = (fid 1)  - i(fid 2)  + (fid 3)  + i(fid 4)
/*                  = (A + iB) - i(C + iD) + (E + iF) + i(G + iH)
/*
/*    sin component = -i(fid 1)  - (fid 2)  + i(fid 3)  - (fid 4)
/*                  = -i(A + iB) - (C + iD) + i(E + iF) - (G + iH)
/*
/* AA = (t1 = 1 real, t2 = real) = A + D + E - H
/* BB = (t1 = 1 real, t2 = imag) = B - C + F + G
/* 
/* CC = (t1 = 2 imag, t2 = real) = B - C - F - G
/* DD = (t1 = 2 imag, t2 = imag) = (-A) - D + E - H
/***/

(void) vvCopy( vAA, vA, size );
(void) vvAdd( vAA, vD, size );
(void) vvAdd( vAA, vE, size );
(void) vvSub( vAA, vH, size );

(void) vvCopy( vBB, vB, size );
(void) vvSub( vBB, vC, size );
(void) vvAdd( vBB, vF, size );
(void) vvAdd( vBB, vG, size );

(void) vvCopy( vCC, vB, size );
(void) vvSub( vCC, vC, size );
(void) vvSub( vCC, vF, size );
(void) vvSub( vCC, vG, size );

(void) vvCopy( vDD, vA, size );
(void) vNeg( vDD, size );
(void) vvSub( vDD, vD, size );
(void) vvAdd( vDD, vE, size );
(void) vvSub( vDD, vH, size );

(void) dWrite( outUnit, vAA, wordLen*size );
(void) dWrite( outUnit, vBB, wordLen*size );
(void) dWrite( outUnit, vCC, wordLen*size );
(void) dWrite( outUnit, vDD, wordLen*size );
