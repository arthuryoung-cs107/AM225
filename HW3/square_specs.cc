#include "square_specs.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

square_specs::square_specs(): n_sqrs(22), N(111),
coords(imatrix(0, n_sqrs-1, 0, 1)), grid_mat(imatrix(0, N-1, 0, N-1)),
grid_vis(imatrix(0, N-1, 0, N-1)), vec_indices(imatrix(0, N*N-1, 0, 1)),
mat_indices(imatrix(0, N-1, 0, N-1)),
ascii_num(new int[n_sqrs]), n_vec(new int[n_sqrs]), nn_vec(new int[n_sqrs]),
sqr_indices((int ***)malloc(n_sqrs*(sizeof(int **)))),
A_mats((double ***)malloc(n_sqrs*(sizeof(double **)))),
A_glue((double ***)malloc(n_sqrs*(sizeof(double **)))),
A_glueT((double ***)malloc(n_sqrs*(sizeof(double **))))
{
  int i, j,  jj;

  char ascii_name[22][2] = {"A", "B", "C", "D", "E", "F", "G", "H", "I",
  "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "."};
  int n_local[22] = {50, 35, 27, 8, 19, 15, 17, 11, 6, 24, 29,
  25, 9, 2, 7, 16, 18, 4, 37, 42, 33};
  int col_local[22] = {0, 50, 85, 85, 93, 50, 65, 82, 82, 88, 0, 29, 54, 63, 63, 54, 70, 29, 33, 70, 0, 0};
  int row_local[22] = {0, 0, 0, 27, 27, 35, 35, 35, 46, 46, 50, 50, 50, 50, 52, 59, 52, 75, 75, 70, 79, 0};
  grid_read();
  for ( i = 0; i < n_sqrs-1; i++)
  {
    ascii_num[i] = (int) *(ascii_name[i]);
    n_vec[i] = n_local[i] - 1;
    nn_vec[i] = n_vec[i]*n_vec[i];
    coords[i][0] = row_local[i];
    coords[i][1] = col_local[i];
    A_mats[i] = dmatrix(0, nn_vec[i]-1, 0, nn_vec[i]-1);
    A_glue[i] = dmatrix(0, nn_vec[i]-1, 0, n_glue-1); // dependencies on glue terms, will be sparse
    A_glueT[i] = dmatrix(0, n_glue-1, 0, nn_vec[i]-1); // dependencies on glue terms, will be sparse
    A_gen(A_mats[i], n_vec[i]);
  }
  ascii_num[n_sqrs-1] = (int) *(ascii_name[i]);
  n_vec[n_sqrs-1] = N; // only for optimized searching
  nn_vec[n_sqrs-1] = n_glue;
  coords[n_sqrs-1][0] = 0;
  coords[n_sqrs-1][1] = 0;
  A_mats[n_sqrs-1] = dmatrix(0, n_glue-1, 0, n_glue-1); // gonna be all over the place
  A_glue[n_sqrs-1] = A_mats[n_sqrs-1];
  A_glueT[n_sqrs-1] = dmatrix(0, n_glue-1, 0, n_glue-1);
  vec_assign();
  glue_assign();
}

square_specs::~square_specs()
{
  free_imatrix(coords, 0, n_sqrs, 0, 1);
  free_imatrix(grid_mat, 0, N-1, 0, N-1);
  free_imatrix(grid_vis, 0, N-1, 0, N-1);
  free_imatrix(vec_indices, 0, N*N-1, 0, 1);
  delete [] ascii_num;
  delete [] n_vec;
  delete [] nn_vec;

}
void square_specs::grid_read()
{
  int n = N;
  int row_count, col_count, i, j, jj;
  char extract[1000];

  FILE * read_file = fopen("./ascii_data/grid_layout.dat", "r");

  row_count = 0;
  n_glue = 0;
  while (fgets(extract, sizeof (extract), read_file))
  {
    for (i = 0; i < n; i++)
    {
      grid_vis[row_count][i] = (int) extract[i];
      if ( extract[i] ==  *"." )
      {
        n_glue++;
      }
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
  nn_glue = n_glue * n_glue;
}
void square_specs::vec_assign()
{
  int i, j, k, check;
  int count = 0;

  for ( k = 0; k < n_sqrs; k++) // brute force search, row stacking
  {
    check = ascii_num[k];
    sqr_indices[k] = vec_indices + count;
    for ( j = coords[k][0]; j < coords[k][0]+n_vec[k]; j++)
    {
      for ( i = coords[k][1]; i < coords[k][1]+n_vec[k]; i++)
      {
        if (check == grid_mat[j][i] )
        {
          vec_indices[count][0] = j; // assign row
          vec_indices[count][1] = i; // asign column
          mat_indices[j][i] = count;
          count++;
        }
      }
    }
  }

  // int check_old = 0;
  // for ( i = 0; i < N*N; i++)
  // {
  //   check = grid_mat[vec_indices[i][0]][vec_indices[i][1]];
  //   if (check != check_old)
  //   {
  //     printf("%d: %c\n", i, check);
  //     check_old = check;
  //   }
  // }

  glue_indices = sqr_indices[n_sqrs-1];
}
void square_specs::grid_mat2vec(double ** mat, double * vec)
{
  int i, j;
  for ( i = 0; i < N*N; i++)
  {
    vec[i] = mat[vec_indices[i][0]][vec_indices[i][1]];
  }
}
void square_specs::grid_vec2mat(double ** mat, double * vec)
{
  int i, j;
  for ( i = 0; i < N*N; i++)
  {
    mat[vec_indices[i][0]][vec_indices[i][1]] = vec[i];
  }
}
void square_specs::A_gen(double ** mat, int n )
{
  int i, j;
  int nn = n*n;
  zerom_init(mat, 0, nn-1, 0, nn-1);

  for ( i = 0; i < nn; i++) mat[i][i] = 4.0;
  for ( i = n; i < nn; i++) mat[i-n][i] = mat[i][i-n] = -1.0;
  for ( i = 1; i < nn; i++)
  {
    if (i%n != 0)
    {
      mat[i-1][i] = -1.0;
      mat[i][i-1] = -1.0;
    }
  }
}
void square_specs::glue_assign()
{
  int i, j, k, row, col, n, i_glue, j_val, nn, row_check, col_check;

  for ( k = 0; k < n_sqrs-1; k++)
  {
    row = coords[k][0];
    col = coords[k][1];
    n = n_vec[k];
    nn = nn_vec[k];
    if (row > 0) // search for glue dependencies of top row
    {
      j_val = 0;
      for ( i = col; i < col+n; i++) // traversing columns
      {
        i_glue = mat_indices[row-1][i] - (N*N - n_glue);
        A_glue[k][j_val][i_glue] = -1.0;
        j_val++;
      }
    }
    if (col > 0) // search for glue dependencies of left column
    {
      j_val = 0;
      for ( j = row; j < row+n; j++) // traversing rows
      {
        i_glue = mat_indices[j][col-1] - (N*N - n_glue);
        A_glue[k][j_val][i_glue] = -1.0;
        j_val+=n;
      }
    }
    if (row+n < N) // search for glue dependencies of bottom row
    {
      j_val = nn - n;
      for ( i = col; i < col+n; i++) // traversing columns
      {
        i_glue = mat_indices[row+1][i] - (N*N - n_glue);
        A_glue[k][j_val][i_glue] = -1.0;
        j_val++;
      }
    }
    if (col+n < N) // search for glue dependencies of right column
    {
      j_val = n-1;
      for ( j = row; j < row+n; j++) // traversing rows
      {
        i_glue = mat_indices[j][col+1] - (N*N - n_glue);
        A_glue[k][j_val][i_glue] = -1.0;
        j_val+=n;
      }
    }
  }

  // now, search for glue-to-glue dependencies
  for ( k = 0; k < n_glue; k++)
  {
    row = glue_indices[k][0];
    col = glue_indices[k][1];
    j_val = mat_indices[row][col] - (N*N - n_glue);
    A_glue[n_sqrs-1][j_val][j_val] = 4.0;

    for ( i = 0; i < n_glue; i++)
    {
      row_check = glue_indices[i][0];
      col_check = glue_indices[i][1];

      if (col_check == col)
      {
        if (row_check == row+1 || row_check == row-1)
        {
          i_glue = mat_indices[row_check][col_check] - (N*N - n_glue);
          A_glue[n_sqrs-1][j_val][i_glue] = -1.0;
        }
      }
      else
      {
        if (row_check == row)
        {
          if (col_check == col+1 || col_check == col-1)
          {
            i_glue = mat_indices[row_check][col_check] - (N*N - n_glue);
            A_glue[n_sqrs-1][j_val][i_glue] = -1.0;
          }
        }
      }
    }
  }
  for ( k = 0; i < n_sqrs; i++)
  {
    for ( j = 0; j < nn_vec[k]; j++)
    {
      for ( i = 0; i < n_glue; i++)
      {
        A_glueT[k][i][j] = A_glue[k][j][i];
      }
    }
  }
}
