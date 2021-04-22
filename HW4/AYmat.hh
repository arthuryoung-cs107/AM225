#ifndef AYMAT_HH
#define AYMAT_HH

class AYmat {
    public:
      int M, N;
      double * A_ptr;

      AYmat(int M_, int N_);
      ~AYmat();

      void print_mat();
      void init_123();
      void init_0();
      void init_randuni();
      void gen_transpose();
      void set(int i, int j, double val);
      double get(int i, int j);

      void mult_set(AYmat * A, AYmat * B, double alpha, double beta);
      void mult_put(AYmat * B, AYmat * C, double alpha, double beta);

    private:
      double ** AT;
};


#endif
