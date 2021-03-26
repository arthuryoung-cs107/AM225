#include "square_specs.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

square_specs::square_specs(): n_sqrs(22), N(111),
coords(imatrix(0, n_sqrs, 0, 1)), grid_mat(imatrix(0, N-1, 0, N-1)), grid_vis(imatrix(0, N-1, 0, N-1)),
ascii_num(new int[n_sqrs]), n_vec(new int[n_sqrs])
{
  int i, j,  jj;

  double *** A_mats = (double ***)malloc(n_sqrs*(sizeof(double **)));

  char ascii_name[22][2] = {".", "A", "B", "C", "D", "E", "F", "G", "H", "I",
  "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U"};
  int n_local[22] = {0, 50, 35, 27, 8, 19, 15, 17, 11, 6, 24, 29,
  25, 9, 2, 7, 16, 18, 4, 37, 42, 33};
  int col_local[22] = {0, 0, 50, 85, 85, 93, 50, 65, 82, 82, 88, 0, 29, 54, 63, 63, 54, 70, 29, 33, 70, 0};
  int row_local[22] = {0, 0, 0, 0, 27, 27, 35, 35, 35, 46, 46, 50, 50, 50, 50, 52, 59, 52, 75, 75, 70, 79};

  for ( i = 0; i < n_sqrs; i++)
  {
    ascii_num[i] = (int) *(ascii_name[i]);
    n_vec[i] = n_local[i] - 1;
    coords[i][0] = row_local[i];
    coords[i][1] = col_local[i];
    if (i > 0)
    {
      A_mats[i] = dmatrix(0, n_vec[i]-1, 0, n_vec[i]-1);
    }
  }

  grid_read();
  A_init();
}

square_specs::~square_specs()
{
  free_imatrix(coords, 0, n_sqrs, 0, 1);
  free_imatrix(grid_mat, 0, N-1, 0, N-1);
  delete [] ascii_num;
  delete [] n_vec;

}

void square_specs::grid_read()
{
  int n = N;
  int row_count, col_count, i, j, jj;
  char extract[1000];

  FILE * read_file = fopen("./ascii_data/grid_layout.dat", "r");

  row_count = 0;
  while (fgets(extract, sizeof (extract), read_file))
  {
    for (i = 0; i < n; i++)
    {
      grid_vis[row_count][i] = (int) extract[i] ;
    }
    row_count++;
  }
  for ( j = 0; j < n; j++)
  {
    jj = n-j-1;
    for ( i = 0; i < n; i++)
    {
      grid_mat[j][i] = grid_vis[jj][i];
    }
  }

  fclose(read_file);
}
void square_specs::grid_mat2vec(double ** mat, double * vec)
{
  int i, j;

  for ( i = 0; i < N; i++)
  {
  
  }
}
// void square_specs::grid_vec2mat(double ** mat, double * vec)
// {
//
// }
