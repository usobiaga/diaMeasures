#include <R.h>
#include <Rinternals.h>

#define STRCMPR(str1, i, str2, j) (strcmp(CHAR(STRING_ELT(str1, i)), CHAR(STRING_ELT(str2, j))))

enum {IRD = 1, IPI};
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

/* make unique */
SEXP call_unique(SEXP x){
  SEXP s, t;
  t = s = PROTECT(allocList(2));
  SET_TYPEOF(s, LANGSXP);
  SETCAR(t, install("unique")); t = CDR(t);
  SETCAR(t,  x);
  UNPROTECT(1);
  return eval(s, R_GlobalEnv);
}


double do_binary(SEXP elem1, SEXP elem2, int byn){
  SEXP longElem, shortElem;
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
  int i, j, ij, index, index2;
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
  }

  PROTECT(names = getAttrib(attrs, R_NamesSymbol)); /* 3 */
  for (i = 0; i < length(attrs); i++)
    setAttrib(answer, install(translateChar(STRING_ELT(names, i))), 
	      VECTOR_ELT(attrs, i));
  setAttrib(answer, install("diagv"), diagVal);

  UNPROTECT(3);
  return answer;
}


