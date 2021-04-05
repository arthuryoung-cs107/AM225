#ifndef SQUARE_SPECS
#define SQUARE_SPECS

class square_specs {
    public:
        const int n_sqrs, N;
        int n_glue;
        int nn_glue;
        int **coords, **grid_mat, **grid_vis,
        ** vec_indices, // query vector index, get matrix index
        ** glue_indices, // query glue index (0-623), get matrix index
        *** sqr_indices, // query square number (0-20), get stack of vector order
        ** mat_indices, // query matrix index, get vector index
        *ascii_num, *n_vec, *nn_vec;
        double *** A_mats;
        double *** A_glue;
        double *** A_glueT;

        square_specs();
        ~square_specs();

        void grid_read();
        void vec_assign();
        void grid_mat2vec(double ** mat, double * vec);
        void grid_vec2mat(double ** mat, double * vec);
        void A_gen(double ** mat, int n);
        void glue_assign();
        void print_adjacents();

    private:

};


#endif
