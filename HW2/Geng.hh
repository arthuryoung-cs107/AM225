#ifndef GENG_HH
#define GENG_HH

#include <cstdio>

class Geng
{
    public:
        char prefix[100];
        char specfile[200];
        FILE * out_file_ptr;

        int dof;
        int eval_count;
        int step_count;
        int writeout_width;
        double t_it;

        double * y_init;
        double * y_it;
        double * y_final;

        double * y_work1;

        Geng(int dof_in);
        ~Geng();

        virtual void init() = 0;
        virtual void eval(double t_in, double * y_in, double * y_out) = 0;
        virtual void write_out() = 0;

        void solve_fixed(double t_start, double t_end, int max_steps);
        void step(double del_t);

    private:

        double *k1;
        double *k2;
        double *k3;
        double *k1b;
        double *k2b;
        double *k3b;
};

#endif
