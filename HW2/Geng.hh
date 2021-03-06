#ifndef SOL_HAMMER_H_HH
#define SOL_HAMMER_H_HH

#include "sol_base.hh"

class Geng
{
    public:
        int dof;
        int eval_count;
        double t_it;

        double * y_init;
        double * y_final;

        Geng(int dof_in);
        virtual ~Geng();

        virtual void solve_fixed(double t_start, double t_end, double del_t);
        virtual void step(double del_t) = 0;
        
        virtual void init(double t_init, ) = 0;
        virtual void eval(double t_in, double * y_in, double * y_out) = 0;

    private:

        double *k1;
        double *k2;
        double *k3;
        double *k1b;
        double *k2b;
        double *k3b;
        double *k1c;
        double *k2c;
        double *k3c;
};

#endif
