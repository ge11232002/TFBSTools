#include "TFBSTools.h"

void reverseMatrix(float **matrix1, float **matrix2, int width){
// reverse the matrix1 and put the results in matrix2
  int i,j;
  for(i=1; i<=width; i++){
    for(j=0; j<=3; j++){
      matrix2[width-i+1][3-j] = matrix1[i][j];
    }
  }
}

void printMatrix(float **matrix, int width){
// print a matrix to R console
  int i, j;
  Rprintf("Printing an matrix \n");
  for(i=0; i<=width; i++){
    for(j=0; j<=3; j++){
      Rprintf("%f\t", matrix[i][j]);
    }
    Rprintf("\n");
  }
}

void printArray(float *array, int width){
// print a array to R console
  int i;
  Rprintf("Printing an arrary\n");
  for(i=0; i<=width; i++){
    Rprintf("%f\t", array[i]);
  }
  Rprintf("\n");
}

struct alignment score(int width1, int width2, float **matrix1, float **matrix2, double open_penalty, double ext_penalty){
// scoring function, the modified Needleman algorithm
  int i,j; 
  struct entry **F, **I, **B, **E;
  F = (struct entry **) R_alloc(width1+1, sizeof(struct entry *));
  I = (struct entry **) R_alloc(width1+1, sizeof(struct entry *));
  B = (struct entry **) R_alloc(width1+1, sizeof(struct entry *));
  E = (struct entry **) R_alloc(width1+1, sizeof(struct entry *));
  for(i=0; i<=width1; i++){
    F[i] = (struct entry *) R_alloc(width2+1, sizeof(struct entry));
    I[i] = (struct entry *) R_alloc(width2+1, sizeof(struct entry));
    B[i] = (struct entry *) R_alloc(width2+1, sizeof(struct entry));
    E[i] = (struct entry *) R_alloc(width2+1, sizeof(struct entry));
  }
  float nogap_score;   // score in a position without gaps
  float start_insert; // variables for cmparing scores, basically
  float extend_insert;
  float start_deletion;
  float extend_deletion;
  float max_score=0;
  float end_insert;
  float end_deletion;
  float end_continue;
  float sum_i; // counter for sums in a position in profile1
  float sum_j; // counter for sums in a position in profile2
  int align_i[40]; // keeping alignment for printing
  int align_j[40];
  int align_length; // length of alignment
  int counter;   // another counter variable
  int number_of_gaps=0; // number of gaps
  int nucleotide; // nucleotide, 0-3 =ACGT
   
  for(i=0; i<=width1; i++){
    for(j=0; j<=width2; j++){
      F[i][j].score = 0;
      F[i][j].cell_score = 0;
      F[i][j].insert = 0;
      F[i][j].deletion = 0;
      F[i][j].father = NULL;
      F[i][j].align[0] = 0;
      F[i][j].align[1] = 0;
      F[i][j].aln_length = 0;
      F[i][j].kind = 'F';

      I[i][j].score = 0;
      I[i][j].cell_score = 0;
      I[i][j].insert = 0;
      I[i][j].deletion = 0;
      I[i][j].father = NULL;
      I[i][j].align[0] = 0;
      I[i][j].align[1] = 0;
      I[i][j].aln_length = 0;
      I[i][j].kind = 'I';

      B[i][j].score = 0;
      B[i][j].cell_score = 0;
      B[i][j].insert = 0;
      B[i][j].deletion = 0;
      B[i][j].father = NULL;
      B[i][j].align[0] = 0;
      B[i][j].align[1] = 0;
      B[i][j].aln_length = 0;
      B[i][j].kind = 'B';

      E[i][j].score = 0;
      E[i][j].cell_score = 0;
      E[i][j].insert = 0;
      E[i][j].deletion = 0;
      E[i][j].father = NULL;
      E[i][j].align[0] = 0;
      E[i][j].align[1] = 0;
      E[i][j].aln_length = 0;
      E[i][j].kind = 'E';
    }
  }
  struct entry *best_pntr; // pointer to the best entry so far
  
  /*------------scoring engine------------*/
  for(i=1; i<=width1; i++){ //go over all pos vs all pos
    for(j=1; j<=width2; j++){
      nogap_score=0; // initialized
      sum_i = 0; // calculate the sum of the position
      sum_j = 0;
      for (nucleotide=0; nucleotide<=3; ++nucleotide){
        // go through nucleotides in position
        nogap_score += pow((matrix1[i][nucleotide] - matrix2[j][nucleotide]), 2);
        sum_i += matrix1[i][nucleotide];
        sum_j += matrix2[j][nucleotide];
      }
      nogap_score = 2 - nogap_score;
      // define the three different scores
      F[i][j].score = nogap_score + F[i-1][j-1].score; // define non-gapped alignment score
      F[i][j].father= &F[i-1][j-1];
      F[i][j].cell_score = nogap_score;
      F[i][j].align[0] = i;
      F[i][j].align[1] = j;

      if(F[i][j].score >= max_score){ // check if best score yet
        max_score = F[i][j].score;
        best_pntr = &F[i][j];
      }

      // inserts in  profile1 (i) profile
      start_insert = F[i-1][j].score - open_penalty; // define cost to open gap-insertion from here
      extend_insert = F[i-1][j].score - ext_penalty; // cost of extending gap-insertion from here
      if(start_insert >= extend_insert){        // take the best one
        I[i][j].score = start_insert;
        I[i][j].father= &F[i-1][j];
        I[i][j].cell_score = - open_penalty;
        I[i][j].insert = 1;
        I[i][j].align[0] = i;
        I[i][j].align[1] = 0;
      }else{
        I[i][j].score = extend_insert;
        I[i][j].father = &I[i-1][j];
        I[i][j].cell_score = - ext_penalty;
        I[i][j].insert = 1;
        I[i][j].align[0] = i;
        I[i][j].align[1] = 0;
      }

      if(I[i][j].score >= max_score){  // update if best score yet
        max_score = I[i][j].score;
        best_pntr = &I[i][j];
      }

      // deletions in profile1 (i) profile
      start_deletion  = F[i][j-1].score - open_penalty; // open deletion gap
      extend_deletion = B[i][j-1].score - ext_penalty;  // extend deletion gap
      if(start_deletion >= extend_deletion){    // check which one is highest
        B[i][j].score = start_deletion;
        B[i][j].father = &F[i][j-1];
        B[i][j].cell_score = - open_penalty;
        B[i][j].deletion = 1;
        B[i][j].align[0] = 0;
        B[i][j].align[1] = j;
      }else{
        B[i][j].score = extend_deletion;
        B[i][j].father = &B[i][j-1];
        B[i][j].cell_score = - ext_penalty;
        B[i][j].deletion = 1;
        B[i][j].align[0] = 0;
        B[i][j].align[1] = j;
      }
      if(B[i][j].score >= max_score){ // update if best sofar
        max_score = B[i][j].score;
        best_pntr = &B[i][j];
      }
      // end alignment after gap
      end_insert = I[i-1][j-1].score + nogap_score;  // start real alignment after insertion-gap
      end_deletion = B[i-1][j-1].score + nogap_score;  // start real alignent after deletion-gap
      end_continue = E[i-1][j-1].score + nogap_score;  // continue real alignment

      if(end_insert >= end_deletion && end_insert >= end_continue ){ // check which score is highest
        E[i][j].score = end_insert;
        E[i][j].father = &I[i-1][j-1];
        E[i][j].cell_score = nogap_score;
        E[i][j].align[0] = i;
        E[i][j].align[1] = j;
      }else if(end_deletion >= end_insert && end_deletion >= end_continue){
        E[i][j].score = end_deletion;
        E[i][j].father = &B[i-1][j-1];
        E[i][j].cell_score = nogap_score;
        E[i][j].align[0] = i;
        E[i][j].align[1] = j;
      }else{
        E[i][j].score = end_continue;
        E[i][j].father = &E[i-1][j-1];
        E[i][j].cell_score = nogap_score;
        E[i][j].align[0] = i;
        E[i][j].align[1] = j;
      }
      if(E[i][j].score > max_score){ // update if highest yet
        max_score = E[i][j].score;
        best_pntr = &E[i][j];
      }
    }
  }
  /*-----------------function for walking through the alignment------------*/
  // starting with the best scoring cell, going back throgh the father-pointers
  struct alignment align;
  align.best_score = max_score;
  counter = 0;
  align_length = 0;
  struct entry *current_pntr = best_pntr; // for walking, start with the best score
  while (current_pntr->father != NULL){ // while the father of the current pointer exists, walk through the best posible alignment
    align_i[counter] = current_pntr->align[0];
    align_j[counter]=  current_pntr->align[1];
    align_length++;
    current_pntr = current_pntr->father;
    counter ++;
  }
  align.length= align_length;
  return align; 
  // Below is not necessary so far.. Do not run them
  /*for(i=align_length-1; i>=0; --i){ // walk through alignment for printing, first profile
    align->over_string[i] = align_i[align_length-i-1]; // fill in alignment in alignment-object
    if(align_i[i] == 0){ // count the number of gaps
      number_of_gaps = number_of_gaps + 1;
    }
  }
  for(i=align_length-1; i>=0; --i){ // walk through alignment for printing, second profile
    align->under_string[i] = align_j[align_length-i-1]; // fill in alignment in alignment-object
    if(align_j[i] == 0){ // count the number of gaps
      number_of_gaps =number_of_gaps +1;
    }
  }
  align->gaps = number_of_gaps; // fill in number of gaps in alignment-object
  */
  //return align;  
}



/* ----------------.Call() Entry points: the main matrixAligner function ------------- */
SEXP matrixAligner(SEXP matrixQuery, SEXP matrixSubject, SEXP open_penalty, SEXP ext_penalty){
  // matrixQuery and matrixSubject are matrix of integers.
  // open_penalty and ext_penalty are numerics.
  int vidd1; // column number of matrixQuery
  int vidd2; // column number of matrixSubject
  int i, j;

  vidd1 = INTEGER(GET_DIM(matrixQuery))[1];
  vidd2 = INTEGER(GET_DIM(matrixSubject))[1];
  //Rprintf("the matrixQuery dim is %d\n", vidd1);
  //Rprintf("the matrixSubject dim is %d\n", vidd2);

  // for simplicity with old code, but at the cost of assignment
  // make another matrix of profile with additional column so make pos 1 is index 1 in the matrix.
  float **matris1, **matris2, **matris3;
  matris1 = (float **) R_alloc(vidd1+1, sizeof(float *));
  matris2 = (float **) R_alloc(vidd2+1, sizeof(float *));
  matris3 = (float **) R_alloc(vidd2+1, sizeof(float *));
  for(i=0; i<=vidd1; i++){
    matris1[i] = (float *) R_alloc(3, sizeof(float));
  }
  for(i=0; i<=vidd2; i++){
    matris2[i] = (float *) R_alloc(3, sizeof(float));
    matris3[i] = (float *) R_alloc(3, sizeof(float));
  }
  
  float *position_weights, *position_weights2;
  position_weights = (float *) R_alloc(vidd1+1, sizeof(float));
  position_weights2 = (float *) R_alloc(vidd2+1, sizeof(float));

  // fill all the elements of matris1 with 0 and position_weights with 0
  for(j=0; j<=3; j++){
    for(i=0; i<=vidd1; i++){
      matris1[i][j] = 0;
    }
    for(i=0; i<=vidd2; i++){
      matris2[i][j] = 0;
      matris3[i][j] = 0;
    }
  }
  for(i=0; i<=vidd1; i++){
    position_weights[i] = 0;
  }
  for(i=0; i<=vidd2; i++){
    position_weights2[i] = 0;
  }

  // fill in the matrix with these data:
  for(i=1; i<=vidd1; i++){
    for(j=0; j<=3; j++){
      matris1[i][j] = (float)INTEGER(matrixQuery)[(i-1)*4+j];
      position_weights[i] += matris1[i][j];
    }
  }
  for(i=1; i<=vidd2; i++){
    for(j=0; j<=3; j++){
      matris2[i][j] = (float)INTEGER(matrixSubject)[(i-1)*4+j];
      position_weights2[i] += matris2[i][j];
    }
  }

  // normalize the data
  for(j=0; j<=3; j++){
    for(i=1; i<=vidd1; i++){
      matris1[i][j] = matris1[i][j] / position_weights[i];
    }
    for(i=1; i<=vidd2; i++){
      matris2[i][j] = matris2[i][j] / position_weights2[i];
    }
  }
  reverseMatrix(matris2, matris3, vidd2);// reverse the second profile for +- scoring

  //printMatrix(matris1, vidd1);
  //printArray(position_weights, vidd1);
  //printMatrix(matris2, vidd2);
  //printArray(position_weights2, vidd2);
  //printMatrix(matris3, vidd2);

  struct alignment score1, score2;
  score1 = score(vidd1, vidd2, matris1, matris2, REAL(open_penalty)[0], REAL(ext_penalty)[0]);
  score2 = score(vidd1, vidd2, matris1, matris3, REAL(open_penalty)[0], REAL(ext_penalty)[0]);
  SEXP maxScore;
  PROTECT(maxScore = NEW_NUMERIC(1));
  if(score1.best_score > score2.best_score){  // what final score is highest?
  //  Rprintf("The best score is %f\n", score1->best_score);
    REAL(maxScore)[0] = (double) score1.best_score;
  }else{
  //  Rprintf("The best score2 is %f\n", score2->best_score);
    REAL(maxScore)[0] = (double) score2.best_score;
  }
  UNPROTECT(1);
  return(maxScore);
  //return R_NilValue;
}

