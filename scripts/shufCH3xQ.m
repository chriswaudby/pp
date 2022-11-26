/***/
/* shufCH3xQ.M: TODO shuffling SOFAST-H(Z/D)QC input to extract DQ spectrum 
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

    (void) setParm( fdata, NDSIZE,   integer( specnum/4 ),  CUR_YDIM ); 
    (void) setParm( fdata, NDTDSIZE, integer( tdSize/4 ),   CUR_YDIM ); 
    (void) setParm( fdata, NDAPOD,   integer( apodSize/4 ), CUR_YDIM ); 
   };

if (sliceCode < 1)
   {
    exit( 0 ); 
   };

/***/
/* Process only every group of 8 complex 1D slices:
/***/

if (sliceCode % 8)
   {
    exit( 0 );
   };

float v1r[size],  v1i[size],  v2r[size],  v2i[size],  v3r[size],  v3i[size],  v4r[size],  v4i[size];
float v5r[size],  v5i[size],  v6r[size],  v6i[size],  v7r[size],  v7i[size],  v8r[size],  v8i[size];
float vA[size],  vB[size],  vC[size],  vD[size],  vE[size],  vF[size],  vG[size],  vH[size];
float vAA[size], vBB[size], vCC[size], vDD[size];

(void) dReadB( inUnit, v1r, wordLen*size );
(void) dReadB( inUnit, v1i, wordLen*size );
(void) dReadB( inUnit, v2r, wordLen*size );
(void) dReadB( inUnit, v2i, wordLen*size );
(void) dReadB( inUnit, v3r, wordLen*size );
(void) dReadB( inUnit, v3i, wordLen*size );
(void) dReadB( inUnit, v4r, wordLen*size );
(void) dReadB( inUnit, v4i, wordLen*size );
(void) dReadB( inUnit, v5r, wordLen*size );
(void) dReadB( inUnit, v5i, wordLen*size );
(void) dReadB( inUnit, v6r, wordLen*size );
(void) dReadB( inUnit, v6i, wordLen*size );
(void) dReadB( inUnit, v7r, wordLen*size );
(void) dReadB( inUnit, v7i, wordLen*size );
(void) dReadB( inUnit, v8r, wordLen*size );
(void) dReadB( inUnit, v8i, wordLen*size );

/* combine neighbouring FIDs */
(void) vvCopy( vA, v1r, size );
(void) vvAdd( vA, v2r, size );
(void) vvCopy( vB, v1i, size );
(void) vvAdd( vB, v2i, size );

(void) vvCopy( vC, v3r, size );
(void) vvAdd( vC, v4r, size );
(void) vvCopy( vD, v3i, size );
(void) vvAdd( vD, v4i, size );

(void) vvCopy( vE, v5r, size );
(void) vvAdd( vE, v6r, size );
(void) vvCopy( vF, v5i, size );
(void) vvAdd( vF, v6i, size );

(void) vvCopy( vG, v7r, size );
(void) vvAdd( vG, v8r, size );
(void) vvCopy( vH, v7i, size );
(void) vvAdd( vH, v8i, size );


/***/
/* Desired output:
/*    cos component = i(fid 12)  + i(fid 34)  + i(fid 56[-90])  + i(fid 78[+90])
/*                  = i(A + iB) + i(C + iD) + i(F - iE) + i(-H + iG)
/*                  = iA - B + iC - D + iF + E - iH - G
/*                  = (-B - D + E - G) + i(A + C + F - H)
/*
/*    sin component = (fid 12)  - (fid 34)  + (fid 56[-90])  - (fid 78[+90])
/*                  = (A + iB) - (C + iD) + (F - iE) - (-H + iG)
/*                  = A + iB - C - iD + F -iE + H - iG
/*                  = (A - C + F + H) + i(B - D - E - G)
/*
/* AA = (t1 = 1 real, t2 = real) = (-B) - D + E - G
/* BB = (t1 = 1 real, t2 = imag) = A + C + F - H
/* 
/* CC = (t1 = 2 imag, t2 = real) = A - C + F + H
/* DD = (t1 = 2 imag, t2 = imag) = B - D - E - G
/***/

(void) vvCopy( vAA, vB, size );
(void) vNeg( vAA, size );
(void) vvSub( vAA, vD, size );
(void) vvAdd( vAA, vE, size );
(void) vvSub( vAA, vG, size );

(void) vvCopy( vBB, vA, size );
(void) vvAdd( vBB, vC, size );
(void) vvAdd( vBB, vF, size );
(void) vvSub( vBB, vH, size );

(void) vvCopy( vCC, vA, size );
(void) vvSub( vCC, vC, size );
(void) vvAdd( vCC, vF, size );
(void) vvAdd( vCC, vH, size );

(void) vvCopy( vDD, vB, size );
(void) vvSub( vDD, vD, size );
(void) vvSub( vDD, vE, size );
(void) vvSub( vDD, vG, size );

(void) dWrite( outUnit, vAA, wordLen*size );
(void) dWrite( outUnit, vBB, wordLen*size );
(void) dWrite( outUnit, vCC, wordLen*size );
(void) dWrite( outUnit, vDD, wordLen*size );
