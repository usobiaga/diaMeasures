#include <R.h>
#include <Rinternals.h>

#define STRCMPR(str1, i, str2, j) (strcmp(CHAR(STRING_ELT(str1, i)), CHAR(STRING_ELT(str2, j))))
#define MIN2(nbr1, nbr2) (nbr1 > nbr2 ? nbr2 : nbr1)
#define MAX2(nbr1, nbr2) (nbr1 < nbr2 ? nbr2 : nbr1)
#define SWAPSEXP(x, y, t) do { t=x; x=y; y=t; } while (0)

enum {IRD = 1, IPI, LEVENSHTEIN};
enum {DICE = 1, JACCARD};

/* match first vector str*/
int matchFirstStr(SEXP str, SEXP vecStr){
  SEXP vecElem;
  R_len_t strl, vecStrl;
  int i, j, z, flag;

  strl = length(str);
  vecStrl = length(vecStr);

  for (i = 0; i < vecStrl; i++){
    vecElem = VECTOR_ELT(vecStr, i);
    if (length(vecElem) != length(str)) continue;
    flag = 0;
    for (j = 0; j < strl; j++){
      if (STRCMPR(vecElem, j, str, j) != 0){ flag = 1; break; }
    }
    if (flag == 0){
      return i;
    }
  }
  error ("matchFirstStr has not found any result");
}

int matchFirstChar(const char *s, SEXP str){
  R_len_t len;
  int i;
  len = length(str);
  for (i = 0; i < len; i++){
    if (strcmp(s, CHAR(STRING_ELT(str, i))) == 0) return i;
  }
  error ("matchFirstChar has not found any results");
}


/* optimize these functions */
SEXP call_unique(SEXP x){
  SEXP s, t;
  t = s = PROTECT(allocList(2));
  SET_TYPEOF(s, LANGSXP);
  SETCAR(t, install("unique")); t = CDR(t);
  SETCAR(t,  x);
  UNPROTECT(1);
  return eval(s, R_GlobalEnv);
}

SEXP call_nchar(SEXP x){
  SEXP s, t;
  t = s = PROTECT(allocList(2));
  SET_TYPEOF(s, LANGSXP);
  SETCAR(t, install("nchar")); t = CDR(t);
  SETCAR(t,  x);
  UNPROTECT(1);
  return eval(s, R_GlobalEnv);
}

SEXP call_unlist(SEXP x){
  SEXP s, t;
  t = s = PROTECT(allocList(2));
  SET_TYPEOF(s, LANGSXP);
  SETCAR(t, install("unlist")); t = CDR(t);
  SETCAR(t, x);
  UNPROTECT(1);
  return (eval(s, R_GlobalEnv));
}

SEXP call_adist(SEXP x){
  SEXP s, t;
  t = s = PROTECT(allocList(2));
  SET_TYPEOF(s, LANGSXP);
  SETCAR(t, install("adist")); t = CDR(t);
  SETCAR(t, x);
  UNPROTECT(1);
  return (eval(s, R_GlobalEnv));
}

double do_binary(SEXP elem1, SEXP elem2, int byn){
  R_len_t  nelem1, nelem2;
  double vals[3] = {0.0, 0.0, 0.0};  /* vals = [a,b,c] */
  double aux;
  int i, j, flag;
  const char *s1, *s2;
  
  flag = 0; aux = 0.0;
  for (i = 0; i < length(elem1); i++){
    for (j = 0; j < length(elem2); j++){
      if (STRCMPR(elem1, i, elem2, j) == 0){
	vals[0] = vals[0] + 1.0;
	flag = 1;
	break;
      }
      if (flag != 1) aux = aux + 1.0;
      flag = 0;
    }
  }
  
  vals[1] = (double)length(elem1) - vals[0];
  vals[2] = (double)length(elem2) - vals[0];
  
  switch (byn){

  case DICE:
    return ((2.0 * vals[0]) / (2.0 * vals[0] + vals[1] + vals[2]));

  case JACCARD:
    return (vals[0] / (vals[0] + vals[1] + vals[2]));
    
  }
}

double do_binary_levenshtein(SEXP elem1, SEXP elem2, SEXP distMat, SEXP strs){
  SEXP t;
  R_len_t  nelem1, nelem2, multiplier, nstrs;
  double  val, val2, distance;
  int i, j, pos1, pos2, totalVal, counter;
  const char *s1, *s2, *s3, *s4;

  /* elem1 must be the long one */
  if (length(elem2) > length(elem1)){
    SWAPSEXP(elem1, elem2, t);
  }

  nstrs = length(strs);
  multiplier = length(elem2);
  distance = 0.0;

  /* Rprintf(" printing binary levenshtein \n");  */
  for (i = 0; i < length(elem1); i++){
    for (j = 0; j < length(elem2); j++){
      if (STRCMPR(elem1, i, elem2, j) == 0){
	s1 = CHAR(STRING_ELT(elem1, i));
	s2 = CHAR(STRING_ELT(elem2, j));
	val = 0.0;
	break;
      }
      pos1 = matchFirstChar(CHAR(STRING_ELT(elem1, i)), strs);
      pos2 = matchFirstChar(CHAR(STRING_ELT(elem2, j)), strs);
      if (i == 0){
	val = REAL(distMat)[nstrs * pos1 + pos2];
	s1 = CHAR(STRING_ELT(elem1, i));
	s2 = CHAR(STRING_ELT(elem2, j));
      } else {
	val2 =  REAL(distMat)[nstrs * pos1 + pos2];
	if (val2 < val){
	  val = val2;
	  s1 = CHAR(STRING_ELT(elem1, i));
	  s2 = CHAR(STRING_ELT(elem2, j));
	}
      }
    }
    distance = distance + 2.0 * val / (double) MAX2(strlen(s1), strlen(s2));
  }

  distance = distance / (double)(length(elem1) + length(elem2));
  Rprintf("distance is %f \n", distance);
  return distance;
}

/* for easier posterior manipulation a matrix is returned rather
   than a vector even even at the cost of hogging more memory */
void do_ird_variable(double *answer, int *obs, SEXP strList, int binary_measure){
  SEXP first, second;
  R_len_t len;
  double binAns;
  int i, j, ij;
  
  len = length(strList);
  ij = 0;
  for (i = 0; i < len; i++){
    for (j = 0; j < i; j++){

      first = VECTOR_ELT(strList, i);
      second = VECTOR_ELT(strList, j);
      
      if (length(first) == 0 || length(second) == 0) continue;
      obs[ij]++;
      
      if (length(first) == 1 && length(second) == 1){  /* simple case */
	if (STRCMPR(first, 0, second, 0) == 0){
	  answer[ij] = answer[ij] + 1.0;
	  ij++;
	  continue;
	}
	ij++;
	continue;
      }
      
      answer[ij] = answer[ij] + do_binary(first, second, binary_measure); /* binary measure */
      ij++;
    }
  }
}

void do_ipi_variable(double *answer, int *obs, SEXP strList, int binary_measure){
  SEXP count, uniqueLemmas, valueMatrix, coidentityMatrix, first, second;
  int i, j, ij, index, index2, aux;
  R_len_t nlem, nlemu;
  double *valueMatrix_ptrs, *coidentityMatrix_ptrs, val, totalCount, *count_ptrs;

  /* count lemmas unique lemmas */

  Rprintf("calling unique list \n");
  PROTECT(uniqueLemmas = call_unique(strList)); /* 1 */
  Rprintf("out of calling unique list \n");
  PROTECT(count = allocVector(REALSXP, length(uniqueLemmas))); /* 2 */
  memset(REAL(count), 0.0, length(uniqueLemmas) * sizeof(double));
  count_ptrs = REAL(count);
  totalCount = 0.0;
  nlem = length(strList);
  for (i = 0; i < nlem; i++){
    if (length(VECTOR_ELT(strList, i)) == 0) continue;
    index = matchFirstStr(VECTOR_ELT(strList, i), uniqueLemmas);
    count_ptrs[index] = count_ptrs[index] + 1.0;
    totalCount = totalCount + 1.0;
  }
  
  /* compute the TAX values */
  
  nlemu = length(uniqueLemmas);
  PROTECT(valueMatrix = allocMatrix(REALSXP, nlemu, nlemu));  /* 3 */
  PROTECT(coidentityMatrix = allocMatrix(REALSXP, nlemu, nlemu)); /* 4 */
  valueMatrix_ptrs = REAL(valueMatrix);
  coidentityMatrix_ptrs = REAL(coidentityMatrix);
  
  for (i = 0; i < nlemu; i++){
    /* fill the diagonal */
    valueMatrix_ptrs[i + i*nlemu] = count_ptrs[i] / totalCount; 
    coidentityMatrix_ptrs[i + i*nlemu] = 0.0;
    
    for (j = 0; j < i; j++){

      first = VECTOR_ELT(uniqueLemmas, i);
      second = VECTOR_ELT(uniqueLemmas, j);
      
      if (length(first) == 1 && length(second) == 1){ /* simple case */
	valueMatrix_ptrs[j + i*nlemu] = valueMatrix_ptrs[i + j*nlemu] = 0.0;
	coidentityMatrix_ptrs[j + i*nlemu] = coidentityMatrix_ptrs[i + j*nlemu] = 1.0;
	continue;
      }
      
      val = do_binary(first, second, binary_measure);
      if (val == 0.0){ /* no matching binary */
	valueMatrix_ptrs[j + i*nlemu] = valueMatrix_ptrs[i + j*nlemu] = 0.0;
	coidentityMatrix_ptrs[j + i*nlemu] = coidentityMatrix_ptrs[i + j*nlemu] = 1.0;
	continue;
      }
      
      valueMatrix_ptrs[j + i*nlemu] = valueMatrix_ptrs[i + j*nlemu] =
	val * (count_ptrs[i] + count_ptrs[j]) / (2.0 * totalCount);

      coidentityMatrix_ptrs[j + i*nlemu] = coidentityMatrix_ptrs[i + j*nlemu] =
	1.0 - valueMatrix_ptrs[j + i*nlemu];
    }
  }

  /* apply identity values */
  
  ij = 0;
  for (i = 0; i < nlem; i++){
    for (j = 0; j < i; j++){
      first = VECTOR_ELT(strList, i);
      second = VECTOR_ELT(strList, j);
      
      if (length(first) == 0 || length(second) == 0){ /* if any na pass */
	ij++;
	continue;
      }
      obs[ij]++;
      index = matchFirstStr(first, uniqueLemmas);
      index2 = matchFirstStr(second, uniqueLemmas);
      val = 1.0 - (valueMatrix_ptrs[index + index2*nlemu] - 1.0) / totalCount;
      answer[ij] = answer[ij] + val / (val + coidentityMatrix_ptrs[index + index2*nlemu]);
      ij++;
    }
  }

  UNPROTECT(4);
}

void do_levenshtein_variable(double *answer, int *obs, SEXP strList, int binary_measure){
  SEXP uniqueStr, strs, strDists, ustrList, binaryValues, maxValues, first, second;
  int i, j, ij, index, index2;
  R_len_t nustrl, nulemmas, nlemmas;
  double *binaryValues_ptrs;

  PROTECT(uniqueStr = call_unique(call_unlist(strList)));  /* 1 */
  PROTECT(strs = call_unlist(strList)); /* 2 */
  PROTECT(strDists = call_adist(strs));  /* 3 */
  PROTECT(ustrList = call_unique(strList)); /* 4 */

  nustrl = length(ustrList);
  PROTECT(binaryValues = allocMatrix(REALSXP, nustrl, nustrl)); /* 5 */
  PROTECT(maxValues = allocMatrix(REALSXP, nustrl, nustrl)); /* 6 */
  binaryValues_ptrs = REAL(binaryValues);
  
  /* ij = 0; */
  for (i = 0; i < nustrl; i++){
    first = VECTOR_ELT(ustrList, i);
    for (j = 0; j < i; j++){
      second = VECTOR_ELT(ustrList, j);
      if (length(first) == 0 || length(second) == 0){
	Rprintf("length 0 stuff \n");
	continue;
      }
      Rprintf("length not zero stuff \n");
      binaryValues_ptrs[j + i*nustrl] = binaryValues_ptrs[i + j*nustrl] =
	do_binary_levenshtein(first, second, strDists, strs);
      /* Rprintf(" binary Values first value is %f \n", REAL(binaryValues)[ij]); */
      /* ij++; */
    }
  }

  Rprintf("not here I guess \n");

  /* TODO the binaryValues are calculated next calculate the results */
  nlemmas = length(strList);
  ij = 0;
  for (i = 0; i < nlemmas; i++){
    for (j = 0; j < i; j++){
      first = VECTOR_ELT(strList, i);
      second = VECTOR_ELT(strList, j);
      if (length(first) == 0 || length(second) == 0){ /* if any na pass */
	ij++;
	continue;
      }
      obs[ij]++;
      index = matchFirstStr(first, ustrList);
      index2 = matchFirstStr(second, ustrList);
      Rprintf("the adding value is %f \n", REAL(binaryValues)[index + index2 * nustrl]);
      answer[ij] = answer[ij] + REAL(binaryValues)[index + index2 * nustrl];
      ij++;
    }
  }
  UNPROTECT(6);
}

SEXP do_ird(SEXP inputList, int binary_measure){
  SEXP answer, obsCount;
  R_len_t nvars, nids;
  int i;

  nids = length(VECTOR_ELT(inputList, 0));
  nvars = length(inputList);
  PROTECT(answer = allocVector(REALSXP, nids * (nids - 1) / 2));
  PROTECT(obsCount = allocVector(INTSXP, nids * (nids - 1) / 2));

  memset(REAL(answer), 0.0, length(answer) * sizeof(double));
  memset(INTEGER(obsCount), 0, length(obsCount) * sizeof(int));
  
  for (i = 0; i < nvars; i++){
    do_ird_variable(REAL(answer), INTEGER(obsCount), VECTOR_ELT(inputList, i), binary_measure);
  }
  
  for (i = 0; i < length(answer); i++){
    REAL(answer)[i] = REAL(answer)[i] * 100.0 / (double) INTEGER(obsCount)[i];
  }
  
  UNPROTECT(2);
  return answer;
}

SEXP do_ipi(SEXP inputList, int binary_measure){
  SEXP answer, nobs, lemmaList;
  R_len_t nvars, nids;
  int i, ansLen;

  nids = length(VECTOR_ELT(inputList, 0));
  nvars = length(inputList);
  ansLen = nids * (nids - 1) / 2;
  
  PROTECT(answer = allocVector(REALSXP, ansLen));
  PROTECT(nobs = allocVector(INTSXP, ansLen));
  
  memset(REAL(answer), 0.0, ansLen * sizeof(double));
  memset(INTEGER(nobs), 0, ansLen * sizeof(int));

  for (i = 0; i < nvars; i++){
    Rprintf("input list has length  %d \n", length(VECTOR_ELT(inputList, i)));
    do_ipi_variable(REAL(answer), INTEGER(nobs), VECTOR_ELT(inputList, i), binary_measure);
  }

  for (i = 0; i < ansLen; i++){
    REAL(answer)[i] = 100.0 * REAL(answer)[i] / (double)(INTEGER(nobs)[i]);
  }
  
  UNPROTECT(2);
  return answer;
}

SEXP do_levenshtein(SEXP inputList, int binary_measure){
  SEXP answer, nobs;
  R_len_t nvars, nids;
  int i, ansLen;

  nids = length(VECTOR_ELT(inputList, 0));
  nvars = length(inputList);
  ansLen = nids * (nids - 1) / 2;

  PROTECT(answer = allocVector(REALSXP, ansLen));
  PROTECT(nobs = allocVector(INTSXP, ansLen));

  memset(REAL(answer), 0.0, ansLen * sizeof(double));
  memset(INTEGER(nobs), 0, ansLen * sizeof(int));

  for (i = 0; i < nvars; i++){
    do_levenshtein_variable(REAL(answer), INTEGER(nobs), VECTOR_ELT(inputList, i), binary_measure);
  }

  for (i = 0; i < ansLen; i++){
    REAL(answer)[i] = 100.0 * REAL(answer)[i] / (double)((INTEGER(nobs)[i]));
  }
  
  UNPROTECT(2);
  return answer;
  
}

SEXP diaMeasure_C(SEXP inputList, SEXP measureR, SEXP binary_measureR, SEXP attrs){
  int measure, binary_measure, i;
  SEXP answer, diagVal, names;

  measure = INTEGER(measureR)[0];
  binary_measure = INTEGER(binary_measureR)[0];
  PROTECT(diagVal = allocVector(INTSXP, 1)); /* 1 */
  
  switch (measure){
    
  case IRD:
    PROTECT(answer = do_ird(inputList, binary_measure)); /* 2 */
    INTEGER(diagVal)[0] = 100;
    break;
    
  case IPI:
    PROTECT(answer = do_ipi(inputList, binary_measure)); /* 2 */
    INTEGER(diagVal)[0] = 100;
    break;

  case LEVENSHTEIN:
    PROTECT(answer = do_levenshtein(inputList, binary_measure)); /* 2 */
    INTEGER(diagVal)[0] = 0;
  }

  PROTECT(names = getAttrib(attrs, R_NamesSymbol)); /* 3 */
  for (i = 0; i < length(attrs); i++)
    setAttrib(answer, install(translateChar(STRING_ELT(names, i))), 
	      VECTOR_ELT(attrs, i));
  setAttrib(answer, install("diagv"), diagVal);

  UNPROTECT(3);
  return answer;
}


