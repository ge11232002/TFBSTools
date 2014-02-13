#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include <string.h>
#include <R_ext/Rdynload.h>

struct entry
// box in score-table
{
  float score;          // the dynamically best-score-yet
  float cell_score;     // real score in this position
  struct entry *father; //  pointer to father box
  int insert;           // insertion, 1=yes, 0 =no
  int deletion;         // deletion,  1=yes, 0 =no
  int align[2];         // first is matrix_1 pos(i) ,second is second matrix (j)
  int aln_length;       // dynamically extended length of alignment so far
  char kind;
};

struct alignment
{
  float best_score; // the final score
  int length;       // the alignment length
  int gaps;         // number of gaps
  int over_string[30];  // string matrix 1 in alignment (gap represented with -1)
  int under_string[30]; // string matrix2 in alignment
};

/* matrixAlignerDynamic.c */
SEXP matrixAligner(SEXP matrixQuery, SEXP matrixSubject, 
    SEXP open_penalty, SEXP ext_penalty);



