/** Example Stan model for fitting an SIR model with gamma-distributed 
 * infectious period to incidence and serial interval data.
 * This version uses the "smoothed" hypoxponential approximation.
 */


functions {
    /** vector field.
     * j is the shape parameter of the gamma distribution we want to approximate.
     * let n = min(2, ceil(j)). let s be the stage number. When s <= n-2, 
     * the rate through the state is a constant j/tau, cf. Erlang.
     */
    vector ode(real t, vector y, real beta, real j, real tau, real h) {
        int n = num_elements(y)-2;
        vector[n+2] dy = rep_vector(0.0, n+2);
        // compute some rates
        real lambda = beta * sum(y[2:n+1]); // FOI: beta * (I_1 + I_2 + ... + I_n)
        real gamma = j / tau; // rate for first n-2 stages
        real frj = 1+j-n; // {j}
        real x1 = 2*j / tau / (1 + frj + sqrt(1-frj*frj + h*h) - h);
        real x2 = 2*j / tau / (1 + frj - sqrt(1-frj*frj + h*h) + h);
        // fill derivative
        dy[1] = -lambda * y[1]; // S, rate is -beta * S * I
        dy[2] = lambda * y[1]; // I_1, in-rate is beta * S * I 
        // first n-2 stages
        for ( s in 1:n-2 ) {
            int i = s+1; // index in y and dy: (S,I_1,I_2,...,I_n,R)
            dy[i] += -gamma * y[i];
            dy[i+1] += gamma * y[i];
        }
        // final 2 stages
        dy[n] += -x1 * y[n]; // I_{n-1}
        dy[n+1] = x1 * y[n] - x2 * y[n+1]; // I_n
        dy[n+2] = x2 * y[n+1]; // R
        return dy;
    }
    /** the distribution of a generation interval T_G,
     * given a Gamma(alpha, beta) distributed duration T_I of the infectious period.
     * The density function of T_G is equal to
     * f_{T_G}(t) = 1/E[T_I] \int_t^{\infty} f_{T_I}(s) ds
     */
    real generation_interval_lpdf(real T, real alpha, real beta) {
        return gamma_lccdf(T | alpha, beta) - log(alpha/beta);
    }
    real generation_interval_lcdf(real T, real alpha, real beta) {
        real x = gamma_lcdf(T | alpha+1, beta);
        real y = gamma_lccdf(T | alpha, beta) + log(T*beta/alpha);
        return log_sum_exp(x, y);
    }
    real generation_interval_lccdf(real T, real alpha, real beta) {
        return log1m_exp(generation_interval_lcdf(T | alpha, beta));
    }
}

data {
    int N; // number of time points
    real<lower=0> T[N]; // time points
    int<lower=0> Y[N]; // observations
    int<lower=1> K; // max number of stages
    real<lower=0> h; // avoid singularities (required for smoothed hypoexp approx)
    real<lower=0> M; // population size
    int L; // number of observed generation/serial intervals
    vector<lower=0>[L] GenInt; // observed generation/serial intervals
    int<lower=2> NumGridPts; // for plotting survival function
}

parameters {
    real<lower=0> beta; // infection rate
    real<lower=1, upper=K> j; // shape parameter
    real<lower=0> tau; // mean duration infectious period
    real<lower=0, upper=1> eps; // initial fraction infected
}

model {    
    // incidence data
    for ( n in 1:K ) { 
        if ( j < n && n < j+1 ) {
            // initial condition
            vector[n+2] u0 = rep_vector(0.0, n+2);
            u0[1] = 1-eps;
            u0[2] = eps;
            // solve ivp
            vector[n+2] sol[N] = ode_rk45(ode, u0, 0.0, T, beta, j, tau, h);
            // compute likelihood of data
            Y[1] ~ poisson((u0[1] - sol[1,1]) * M);
            for ( i in 2:N ) {
                Y[i] ~ poisson((sol[i-1,1] - sol[i,1]) * M);
            }
        } // else do nothing
    }
        
    // additional data to inform the distribution of the infectious period
    for ( ell in 1:L ) {
        GenInt[ell] ~ generation_interval(j, j/tau);
    }
}

generated quantities {
    // incidence
    real Yhat[N];
    int Ysim[N];
    matrix[3,N] sol; // S, I, R
    // gen intervals
    real gen_int_surv[NumGridPts];
    // simulate incidence data
    for ( n in 1:K ) { 
        if ( j < n && n < j+1 ) {
            // initial condition
            vector[n+2] u0 = rep_vector(0.0, n+2);
            u0[1] = 1-eps;
            u0[2] = eps;
            // solve ivp
            vector[n+2] fullsol[N] = ode_rk45(ode, u0, 0.0, T, beta, j, tau, h);
            // simulate data
            for ( i in 1:N ) {
                real yhat;
                if ( i == 1 ) {
                    yhat = (u0[1] - fullsol[i,1]) * M;
                } else {
                    yhat = (fullsol[i-1,1] - fullsol[i,1]) * M;
                }
                Yhat[i] = yhat;
                Ysim[i] = poisson_rng(yhat);
                sol[1,i] = fullsol[i][1]; // S
                sol[2,i] = sum(fullsol[i][2:n+1]); // I
                sol[3,i] = fullsol[i][n+2]; // R
            }    
        } // else do nothing
    } // for loop to determine n
    // compute survival function for generation interval
    for ( i in 1:NumGridPts ) {
        real t = (i-1) * max(GenInt) / (NumGridPts-1);
        gen_int_surv[i] = generation_interval_lccdf(t | j, j/tau);
    }
}
