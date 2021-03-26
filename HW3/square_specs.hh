#ifndef SQUARE_SPECS
#define SQUARE_SPECS

class square_specs {
    public:
        const int n_sqrs, N;
        int **coords, **grid_mat, **grid_vis;
        int * ascii_num;
        int * n_vec;
        double *** A_mats;

        square_specs();
        ~square_specs();

        void grid_read();
        void init();
        void grid_mat2vec(double ** mat, double * vec);
        void grid_vec2mat(double ** mat, double * vec);
    private:

};


#endif
