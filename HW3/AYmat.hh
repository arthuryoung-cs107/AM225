#ifndef AYMAT_HH
#define AYMAT_HH

class AYmat {
    public:
      int M, N;
      double ** A;
      double ** AT;

      AYmat(int M_, int N_);
      ~AYmat();

      void print_mat();
      void init_123();
      void init_randuni();
      void gen_transpose();

    private:

};


#endif
