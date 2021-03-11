#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <dirent.h>
#include <time.h>

typedef struct trellis_entry_t{
    double loglike;
    int backptr_f0;
    int backptr_fdot0;
} trellis_entry_t;


// 2-D Gaussian distribution
double fokker_planck_pdf(double x, double y, double sigma, double z)
{
    // det(V) = z**4*sigma**2/4
    double c = 1/(2*M_PI*pow(sigma,2)*pow(z,2)/sqrt(12.0));
    double rho = sqrt(3)/2;
    double x_variance = pow(z,3)/3*pow(sigma,2);
    double y_variance = z*pow(sigma,2);
    double exp_term = -1/(2*(1-pow(rho,2)))*(pow(x,2)/x_variance + pow(y,2)/y_variance - 2*rho*x*y/(sqrt(x_variance)*sqrt(y_variance)));
    return c*exp(exp_term);
}

double fokker_planck(double sigma, double z, double *x, int x_len, double *y, int y_len, double df, double dfdot, double *pdf)
{
    double x_mean = x[(x_len+1)/2-1];
    double y_mean = y[(y_len+1)/2-1];
    double running_sum = 0;
    double max = -INFINITY;
    for(int i = 0; i < x_len; i++)
    {
        for(int j = 0; j < y_len; j++)
        {
            pdf[j*x_len + i] = fokker_planck_pdf(x[i]-x_mean, y[j]-y_mean, sigma, z)*df*dfdot;
            running_sum += pdf[j*x_len + i];
            if(pdf[j*x_len + i] > max)
                max = pdf[j*x_len + i];
        }
    }
    if(running_sum == 0)
    {
        pdf[((y_len+1)/2)*x_len + ((x_len+1)/2)] = 1;
        return 1;
    }
    else
    {
        for(int i = 0; i < x_len; i++)
        {
            for(int j = 0; j < y_len; j++)
            {
                pdf[j*x_len + i] /= running_sum;
            }
        }
        return max / running_sum;
    }
}

void shift_pdf(double *pdf, int pdf_f_len, int pdf_fdot_len, int n, double *shifted_pdf)
{
    int mid = (pdf_fdot_len+1)/2;
    memset(shifted_pdf, 0, sizeof(double)*pdf_fdot_len*pdf_f_len);
    if(n < mid+1)
    {
        memcpy(shifted_pdf, &(pdf[(mid-n)*pdf_f_len]), sizeof(double)*pdf_f_len*(pdf_fdot_len - (mid - n)));
    }
    else
    {
        memcpy(&(shifted_pdf[(n-mid)*pdf_f_len]), &(pdf[0]), sizeof(double)*pdf_f_len*(pdf_fdot_len-(n-mid)));
    }
}

int main()
{
    opterr = 0;

    static struct option long_options[] = {
        {"f0", required_argument, 0, 0},
        {"fdot0", required_argument, 0, 0},
        {"z", required_argument, 0, 0},
        {"sigma", required_argument, 0, 0},
        {0, 0, 0, 0}
    };

    double *f0, *fdot0, *z, sigma, df0, dfdot0;
    int obs_len, f0_len, fdot0_len;

    sigma = 1e-16;

    f0_len = 143;
    df0 = 1e-08;
    f0 = malloc(sizeof(double)*f0_len);
    FILE *f0_file = fopen("sofiaf0.dat", "r");
    for(int i = 0; i < f0_len; i++)
    {
        fscanf(f0_file, "%lf", &(f0[i]));
    }
    fclose(f0_file);

    fdot0_len = 23;
    dfdot0 = 1e-14;
    fdot0 = malloc(sizeof(double)*fdot0_len);
    FILE *fdot0_file = fopen("sofiafdot0.dat", "r");
    for(int i = 0; i < fdot0_len; i++)
    {
        fscanf(fdot0_file, "%lf", &(fdot0[i]));
    }
    fclose(fdot0_file);

    obs_len = 10;
    z = malloc(sizeof(double)*obs_len);
    FILE *z_file = fopen("sofiaz.dat", "r");
    for(int i = 0; i < obs_len; i++)
    {
        fscanf(z_file, "%lf", &(z[i]));
    }
    fclose(z_file);
    printf("%.16lf\n", z[0]);
    printf("%.16lf\n", z[1]);
    trellis_entry_t *trellis = malloc(sizeof(trellis_entry_t)*f0_len*fdot0_len*obs_len);

    // Initialise first trellis timestep
    for(int i = 0; i < fdot0_len; i++)
    {
        for(int j = 0; j < f0_len; j++)
        {
            trellis[i*f0_len + j].loglike = cos(2*M_PI*(z[0]*f0[j] + 0.5*z[0]*z[0]*fdot0[i]));
            trellis[i*f0_len + j].backptr_f0 = -1;
            trellis[i*f0_len + j].backptr_fdot0 = -1;
        }
    }
    for(int t = 1; t < obs_len; t++)
    {
        for(int i = 0; i < fdot0_len; i++)
        {
            for(int j = 0; j < f0_len; j++)
            {
                trellis[(t*f0_len*fdot0_len) + i*f0_len + j].loglike = -INFINITY;
                trellis[(t*f0_len*fdot0_len) + i*f0_len + j].backptr_f0 = -1;
                trellis[(t*f0_len*fdot0_len) + i*f0_len + j].backptr_fdot0 = -1;
            }
        }
    }

    for(int t = 1; t < obs_len; t++)
    {
        // This expression is from line 37 of viterbi_Pulse_3dcol.m
        int f0_spread = ceil(sigma*pow(z[t], 1.5)/sqrt(3)*13/df0);
        // Sanitize f0_spread
        if(f0_len < f0_spread)
        {
            f0_spread = f0_len;
        }
        if(f0_spread < 1)
        {
            f0_spread = 1;
        }
        if(f0_spread % 2 == 0)
        {
            f0_spread -= 1;
        }
        printf("%d at %d\n", f0_spread, t);
        double *pdf = malloc(sizeof(double)*f0_spread*fdot0_len);
        double max = fokker_planck(sigma, z[t], f0, f0_spread, fdot0, fdot0_len, df0, dfdot0, pdf);
        printf("Max: %.16lf\n", max);

        FILE *pdf_dump = fopen("pdf_dump.dat", "w");
        for(int i = 0; i < f0_spread*fdot0_len; i++)
            fprintf(pdf_dump, "%.16lf\n", pdf[i]);
        fclose(pdf_dump);
        
        trellis_entry_t *prev = &(trellis[(t-1)*f0_len*fdot0_len]);
        trellis_entry_t *curr = &(trellis[t*f0_len*fdot0_len]);

        double *obslikes = malloc(sizeof(double)*f0_len*fdot0_len);
        for(int fdot = 0; fdot < fdot0_len; fdot++)
        {
            for(int f = 0; f < f0_len; f++)
            {
                obslikes[fdot*f0_len + f] = cos(2*M_PI*(z[t]*f0[f] + 0.5*z[t]*z[t]*fdot0[fdot]));
            }
        }
        FILE *mun_dump = fopen("mun_dump.dat", "w");
        for(int prev_fdot = 0; prev_fdot < fdot0_len; prev_fdot++)
        {
            for(int prev_f = 0; prev_f < f0_len; prev_f++)
            {
                int central_freq_bin = round(((f0[prev_f] - z[t]*fdot0[prev_fdot]) - f0[0])/df0);
                fprintf(mun_dump, "%d\n", central_freq_bin);
                for(int candidate_fdot = 0; candidate_fdot < fdot0_len; candidate_fdot++)
                {
                    for(int candidate_f = central_freq_bin-(f0_spread-1)/2; candidate_f <= central_freq_bin+(f0_spread-1)/2; candidate_f++)
                    {
                        int f_offset = candidate_f - central_freq_bin; // f offset into pdf from origin
                        int fdot_offset = prev_fdot - candidate_fdot; // fdot offset into pdf from origin
                        int pdf_origin_f = (f0_spread+1)/2-1;
                        int pdf_origin_fdot = (fdot0_len+1)/2-1;
                        double trans_prob;
                        if(pdf_origin_f + f_offset < 0 || pdf_origin_f + f_offset >= f0_spread || pdf_origin_fdot + fdot_offset < 0 || pdf_origin_fdot + fdot_offset >= fdot0_len || candidate_f < 0 || candidate_f >= f0_len)
                        {
                            continue;
                        }
                        else
                        {
                            trans_prob = pdf[(pdf_origin_fdot + fdot_offset)*f0_spread + (pdf_origin_f + f_offset)];
                        }
                        double loglike = log(trans_prob) + obslikes[candidate_fdot*f0_len + candidate_f] + prev[prev_fdot*f0_len + prev_f].loglike - log(max);
                        if(loglike > curr[candidate_fdot*f0_len + candidate_f].loglike)
                        {
                            curr[candidate_fdot*f0_len + candidate_f].loglike = loglike;
                            curr[candidate_fdot*f0_len + candidate_f].backptr_f0 = prev_f;
                            curr[candidate_fdot*f0_len + candidate_f].backptr_fdot0 = prev_fdot;
                        }
                    }
                }
            }
        }
        fclose(mun_dump);
        free(pdf);
        free(obslikes);
    }
    FILE *delta_dump = fopen("delta_dump.dat", "w");
    double best_loglike = -INFINITY;
    int best_backptr_f = -1;
    int best_backptr_fdot = -1;
    int best_f = -1;
    int best_fdot = -1;
    trellis_entry_t *last = &trellis[(obs_len-1)*f0_len*fdot0_len];
    for(int fdot = 0; fdot < fdot0_len; fdot++)
    {
        for(int f = 0; f < f0_len; f++)
        {
            if(last[fdot*f0_len + f].loglike >= best_loglike)
            {
                best_loglike = last[fdot*f0_len + f].loglike;
                best_backptr_f = last[fdot*f0_len + f].backptr_f0;
                best_backptr_fdot = last[fdot*f0_len + f].backptr_fdot0;
                best_f = f;
                best_fdot = fdot;
            }
            fprintf(delta_dump, "%.16lf\n", last[fdot*f0_len + f].loglike);
        }
    }
    int *f_path = malloc(sizeof(int)*obs_len);
    int *fdot_path = malloc(sizeof(int)*obs_len);
    f_path[obs_len-1] = best_f;
    fdot_path[obs_len-1] = best_fdot;
    for(int t = obs_len-2; t >= 0; t--)
    {
        f_path[t] = trellis[(t+1)*fdot0_len*f0_len + fdot_path[t+1]*f0_len + f_path[t+1]].backptr_f0;
        fdot_path[t] = trellis[(t+1)*fdot0_len*f0_len + fdot_path[t+1]*f0_len + f_path[t+1]].backptr_fdot0;
    }
    for(int t = 0; t < obs_len; t++)
    {
        printf("%d\t%d\n", fdot_path[t]+1, f_path[t]+1);
    }
    printf("Best loglike: %.16lf\n", best_loglike);
}
