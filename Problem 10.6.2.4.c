#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define tMax 2
#define dx (M_PI/20)
#define xf M_PI
#define xi 0
#define D 1.0
#define F 0.0
#define n 1.0

#define lambda .4
#define dt (lambda*dx*dx/D)

double u0(double x);
double exact(double x,double t);

void main()
{
    int iterations;
    double tolerance;
    printf("How many iterations (at most)?");
    scanf("%d", &iterations);
//    printf("Tolerance?");
//    scanf("%lf",&tolerance);
    int N = (int)((xf-xi)/dx), M = (int)(tMax/dt);
    double u[M+1][N+1];
    //double lambda = D*dt/dx/dx;
    double a = 1 + 2*lambda, b = -lambda, c = -lambda;
    double alpha[N],g[N];

    int count = 0, i, j, t;
    double k, l;

    for(i=0;i<=M;i++)
    {
        for(j=0;j<N;j++)
        {
            u[i][j]=u0(j*dx);
        }
    }
    for(i=0;i<=M;i++)
    {
        //u[i][0] = 0;
        u[i][N] = 0;
    }

    double temp;
    for(i=1;i<=M;i++)
    {
        for(k=0;k<1000; k++)
        {
            u[i][0]=(u[i-1][0]-b*u[i][1]-c*u[i][1])/a;
            for(j=1;j<=N-1;j++)
            {
                temp=u[i][j];
                u[i][j]=(u[i-1][j]-b*u[i][j-1]-c*u[i][j+1])/a;
//                if(fabs(temp-u[i][j])<tolerance)
//                    break;
            }
        }
    }

    FILE *fptr;
    fptr=fopen("GSdata4p2.csv", "w");
    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
            {printf("%lf ",u[i][j]);
            fprintf(fptr, "%lf, ", u[i][j]);}
        printf("\n");
        fprintf(fptr, "\n");
    }
    printf("\n\n");
    fclose(fptr);

    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
        {
            printf("%lf ",exact(j*dx,i*dt));
        }
        printf("\n");
    }
    fptr=fopen("GSError4p2.csv","w");
    for(i=0; i<M; i++)
    {
        for(j=0; j<N; j++)
        {

            if(u[i][j]==0)
                fprintf(fptr,"%lf, ",0);
            else
            fprintf(fptr,"%lf, ",fabs(u[i][j]-exact(j*dx,i*dt))/exact(j*dx,i*dt));
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

double u0(double x)
{
    return cos((n-.5)*x);
}
double exact(double x, double t)
{
    return exp(-(n-.5)*(n-.5)*t)*cos((n-.5)*x);
}

