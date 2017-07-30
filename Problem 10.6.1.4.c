#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define tMax 2.0
#define dx (M_PI/20)
#define xf M_PI
#define xi 0
#define D 1.0
#define F 0.0
#define n 1.0

#define lambda .8
#define dt (lambda*dx*dx/D)

double u0(double x);
double exact(double x,double t);

void main()
{
    double tolerance;
    printf("How many iterations (at most)?");
    scanf("%d", &iterations);
//    printf("Tolerance?");
//    scanf("%lf",&tolerance);
    int N = (int)((xf-xi)/dx), M = (int)(tMax/dt);
    double u[M+1][N+1];
    double a = 1+2*lambda;
    double  b = -lambda;
    double  c = -lambda;
    double alpha[N+1],g[N+1];
    double exacta[M+1][N+1];

    int count = 0, i, j, t;
    double k, l;
    FILE *fptr;
    fptr=fopen("GSdata8.csv", "w");

    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
        {
            printf("%lf ",exact(j*dx,i*dt));
            exacta[i][j]=exact(j*dx,i*dt);
        }
        printf("\n");
    }
printf("\n");
printf("\n");


    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
        {
            u[i][j]=u0(j*dx);
        }
    }

    for(i=0;i<=M;i++)
    {
        u[i][0] = 0;
        u[i][N] = 0;
    }

    double temp;

    for(i=1;i<=M;i++)
    {
        for(k=0;k<iterations;k++)
        {
            for(j=1;j<=N-1;j++)
            {
                temp=u[i][j];
                u[i][j]=(u[i-1][j]-b*u[i][j-1]-c*u[i][j+1])/a;
            }
        }
    }




    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
            {
            printf("%lf ",u[i][j]);
            fprintf(fptr, "%lf, ", u[i][j]);
            }
        printf("\n");
        fprintf(fptr, "\n");
    }
    fclose(fptr);


    fptr=fopen("GSError8.csv","w");
    for(i=0; i<M; i++)
    {
        for(j=0; j<N; j++)
        {

            if(u[i][j]==0)
                fprintf(fptr,"%lf, ",0);
            else
            fprintf(fptr,"%lf, ",100*fabs(u[i][j]-exact(j*dx,i*dt))/exact(j*dx,i*dt));
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

double u0(double x)
{
    return sin(n*x);
}

double exact(double x,double t)
{
    return exp(-n*n*t)*sin(n*x);
}
