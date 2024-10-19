#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define OK 0
#define ERROR 1
#define pi atan(1)*4

void fill_matrices(FILE *file, double *X, double *Y,const int nj,const int ni);
void print_martix(const double *a,const int index,const int ni,const int nj);
int read_size(FILE *file,int *ni,int *nj);
int read_flight_input(double *Mach,double *alpha,int *i_TEL,int *i_TEU,double *rho,double *P,double *gamma,float *temp,
    float *temp1,float *temp2,float *temp3,float *temp4);
int offset2d(const int i,const int j,const int ni);
int offset3d(const int i,const int j,const int k,const int ni,const int nj);
int intial_Q(const int ni,const int nj,double *q,const double rho,const double Mach,const double gamma,const double P,
    const double alpha);
void output_matrix_3D(const double *a,const int ni,const int nj);
void output_sheet(const double *a,const int ni,const int nj,const int layer);
void output_matrix_2D(double *a,const int col,const int row);
int jacobian(const int ni,const int nj,double *X,double *Y,double *J,double *xi_x,double *xi_y,double *eta_x,double *eta_y);
int RHS(const int ni,const int nj,double *S,double *W, double *q, double *J, double *xi_x, double *xi_y, double *eta_x,
    double *eta_y, const double gamma,const double dt);
int smooth(double *q, double *s, double *jac, double *xx, double *xy,double *yx, double *yy, int id, int jd, double *s2,
           double *rspec, double *qv, double *dd,double epse, double gamma, double fsmach, double dt);
int update_Q(const int ni,const int nj,double *S,double *q,double *J);
int BC(const int ni,const int nj, double *q,const int i_TEL,const int i_TEU,const double gamma, double *eta_x, double *eta_y,
    double *J, double *xi_x, double *xi_y);
int LHSX(const int j,const int ni,const int nj,double *A,double *B,double *C,double *q, double gamma,double *xi_x,
    double *xi_y, double dt);
int smoothx(double *q, double *xx, double *xy, int id, int jd, double *a,double *b, double *c, int j,double *jac, double *drr,
    double *drp, double *rspec, double *qv, double *dd,double epsi, double gamma, double fsmach, double dt);
int LHSY(const int i,const int ni,const int nj,double *A,double *B,double *C,double *q, double gamma,double *eta_x,
    double *eta_y, double dt);
int smoothy(double *q, double *yx, double *yy, int id, int jd, double *a,double *b, double *c, int i,double *jac, double *drr,
    double *drp, double *rspec, double *qv, double *dd,double epsi, double gamma, double fsmach, double dt);
int btri4s(double *a, double *b, double *c, double *f, int kd, int ks, int ke);
int STEP(const int ni,const int nj,double *S,double *W, double *q, double *J, double *xi_x, double *xi_y, double *eta_x,
    double *eta_y,const double gamma,const double dt,double *s2,double *rspec,double *qv,double *dd,const double epse
    , const double Mach,double *A, double *B, double *C, double *D, double *drr, double *drp,const double epsi
    ,const int i_TEL,const int i_TEU);

int main() {
    int nj, ni,i_TEL,i_TEU;
    double Mach,alpha,rho,P,gamma;
    float temp,temp1,temp2,temp3,temp4;
    const double dt = pow(10,-4),epse = 0.06,epsi=0.5*epse;

    FILE *grid_file = fopen("roi_mesh.txt", "r");
    if(read_size(grid_file,&ni,&nj) != OK) {
        printf("Something is worng with the input\n");
        return ERROR;
    }

    double *X = (double *)malloc((unsigned) nj * ni * sizeof(double));
    double *Y = (double *)malloc((unsigned) nj * ni * sizeof(double));
    double *J = (double *)malloc((unsigned) nj * ni * sizeof(double)); //Jacobian matrix
    double *xi_y = (double *)malloc((unsigned) nj * ni * sizeof(double));
    double *eta_y = (double *)malloc((unsigned) nj * ni * sizeof(double));
    double *xi_x = (double *)malloc((unsigned) nj * ni * sizeof(double));
    double *eta_x = (double *)malloc((unsigned) nj * ni * sizeof(double));
    double *E = (double *)malloc((unsigned) nj * ni * 4 * sizeof(double));
    double *F = (double *)malloc((unsigned) nj * ni * 4 * sizeof(double));
    double *S = (double *)malloc((unsigned) nj * ni * 4 * sizeof(double));
    double *q = (double *)malloc((unsigned) nj * ni * 4 * sizeof(double));
    double *W = (double *)malloc((unsigned)  fmax(ni,nj) * 4 * sizeof(double));
    double *s2 = (double *)malloc((unsigned)  fmax(ni,nj) *  sizeof(double));
    double *rspec = (double *)malloc((unsigned)  fmax(ni,nj) *  sizeof(double));
    double *qv = (double *)malloc((unsigned)  fmax(ni,nj) *  sizeof(double));
    double *dd = (double *)malloc((unsigned)  fmax(ni,nj) *  sizeof(double));
    double *A = (double *)malloc((unsigned)  fmax(ni,nj) * 16 * sizeof(double));
    double *B = (double *)malloc((unsigned)  fmax(ni,nj) * 16 * sizeof(double));
    double *C = (double *)malloc((unsigned)  fmax(ni,nj) * 16 * sizeof(double));
    double *D = (double *)malloc((unsigned)  fmax(ni,nj) * 4 * sizeof(double));
    double *drr = (double *)malloc((unsigned)  fmax(ni,nj) *  sizeof(double));
    double *drp = (double *)malloc((unsigned)  fmax(ni,nj) *  sizeof(double));

    if (X == NULL || Y == NULL) {
        perror("Memory allocation failed");
        fclose(grid_file);
        return ERROR;
    }

    /* Reading the input */
    if(read_flight_input(&Mach,&alpha,&i_TEL,&i_TEU,&rho,&P,&gamma,&temp,&temp1,&temp2,&temp3,&temp4) != OK){
        printf("Something is worng with the input\n");
        return ERROR;
    }
    alpha = alpha*pi/180;

    for(int i=0;i<fmax(ni,nj);i++) {
        s2[i]=0;
        rspec[i]=0;
        qv[i]=0;
        dd[i]=0;
    }

    fill_matrices(grid_file, X, Y, nj, ni);
    fclose(grid_file);
    intial_Q(ni,nj,q,rho,Mach,gamma,P,alpha);
    jacobian(ni,nj,X,Y,J,xi_x,xi_y,eta_x,eta_y);
    const int iter = STEP(ni,nj,S,W,q,J,xi_x,xi_y,eta_x,eta_y,gamma,dt,s2,rspec,qv,dd,epse,Mach,A,B,C,D,drr,drp,epsi,
        i_TEL,i_TEU);
    printf("number of iterations are: %d\n",iter);


    // output_sheet(q,ni,nj,0);
    // output_sheet(q,ni,nj,1);
    // output_sheet(q,ni,nj,2);
    // output_sheet(q,ni,nj,3);
    // output_matrix_3D(q,ni,nj);
    output_matrix_2D(J,ni,nj);

    // print_martix(X,1,ni,nj);
    // print_martix(Y,2,ni,nj);
    printf("ni = %d, nj= %d \n",ni,nj);
    printf("Mach=%g, Attack Angle=%g [rad], i_TEL=%d, i_TEU=%d, Density=%g, Pressure=%g , gamma=%g\n",Mach,alpha,
        i_TEL,i_TEU,rho,P,gamma);

    free(X);
    free(Y);
    free(J);
    free(xi_y);
    free(eta_y);
    free(xi_x);
    free(eta_x);
    free(E);
    free(F);
    free(S);
    free(q);
    free(W);
    free(s2);
    free(rspec);
    free(qv);
    free(dd);
    free(A);
    free(B);
    free(C);
    free(D);
    free(drr);
    free(drp);

    return OK;
}

int read_size(FILE *file,int *ni,int *nj) {
    if (file == NULL) {
        perror("Error opening file");
        return ERROR;
    }

    fscanf(file, "%d %d", nj, ni);
    return OK;
}

// Function to read and fill matrices X and Y from a file
void fill_matrices(FILE *file, double *X, double *Y,const int nj,const int ni) {
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            fscanf(file, "%lf", &X[offset2d(i,j,ni)]);
        }
    }

    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            fscanf(file, "%lf", &Y[offset2d(i,j,ni)]);
        }
    }
}

void print_martix(const double *a,const int index,const int ni,const int nj) {
    if(index==1) {
        printf("Matrix X:\n");
        for (int j = 0; j < nj; j++) {
            for (int i = 0; i < ni; i++) {
                printf("%lf ", a[offset2d(i,j,ni)]);
            }
            printf("\n");
        }
    }
    else {
        printf("Matrix Y:\n");
        for (int j = 0; j < nj; j++) {
            for (int i = 0; i < ni; i++) {
                printf("%lf ", a[offset2d(i,j,ni)]);
            }
            printf("\n");
        }
    }
}

void output_matrix_3D(const double *a,const int ni,const int nj) {
    FILE *fpo = NULL;
    fpo = fopen("q.txt","wt");
    for (int i = 0; i < ni; i++) {
        for(int j=0;j<nj;j++) {
            for(int k=0;k<4;k++) {
                if(k==3)
                        fprintf(fpo,"q:%d = %lf  \n", k,a[offset3d(i,j,k,ni,nj)]);
                else
                        fprintf(fpo,"q:%d = %lf, ", k,a[offset3d(i,j,k,ni,nj)]);
            }
        }
    }
}

void output_sheet(const double *a,const int ni,const int nj,const int layer) {
    if(layer==0) {
        FILE *fpo = NULL;
        fpo = fopen("Q0_mach.txt","wt");
        for (int j = 0; j < nj; j++) {
            for(int i=0;i<ni;i++) {
                if(i==ni-1)
                    // fprintf(fpo,"%lf cell:%d \n", a[offset3d(i,j,k,ni,nj)],offset3d(i,j,k,ni,nj));
                        fprintf(fpo,"%g\n",a[offset3d(i,j,layer,ni,nj)]);
                else
                    // fprintf(fpo,"%lf cell:%d ", a[offset3d(i,j,k,ni,nj)],offset3d(i,j,k,ni,nj));
                        fprintf(fpo,"%g ",a[offset3d(i,j,layer,ni,nj)]);

            }
        }
        fclose(fpo);
    }
    else if(layer==1) {
        FILE *fpo = NULL;
        fpo = fopen("Q1_mach.txt","wt");
        for (int j = 0; j < nj; j++) {
            for(int i=0;i<ni;i++) {
                if(i==ni-1)
                        fprintf(fpo,"%g\n",a[offset3d(i,j,layer,ni,nj)]);
                else
                        fprintf(fpo,"%g ",a[offset3d(i,j,layer,ni,nj)]);

            }
        }
        fclose(fpo);
    }
    else if(layer==2){
        FILE *fpo = NULL;
        fpo = fopen("Q2_mach.txt","wt");
        for (int j = 0; j < nj; j++) {
            for(int i=0;i<ni;i++) {
                if(i==ni-1)
                        fprintf(fpo,"%g\n",a[offset3d(i,j,layer,ni,nj)]);
                else
                        fprintf(fpo,"%g ",a[offset3d(i,j,layer,ni,nj)]);

            }
        }
        fclose(fpo);
    }
    else{
        FILE *fpo = NULL;
        fpo = fopen("Q3_mach.txt","wt");
        for (int j = 0; j < nj; j++) {
            for(int i=0;i<ni;i++) {
                if(i==ni-1)
                        fprintf(fpo,"%g\n",a[offset3d(i,j,layer,ni,nj)]);
                else
                        fprintf(fpo,"%g ",a[offset3d(i,j,layer,ni,nj)]);

            }
        }
        fclose(fpo);
    }
}

void output_matrix_2D(double *a,const int col,const int row){
    FILE *fpo=NULL;
    fpo = fopen("Jacobian.txt","wt");
    for (int i = 0; i < row; i++) {
        for(int j=0;j<col;j++) {
            if(j==col-1)
                fprintf(fpo,"%g\n", a[offset2d(j, i, col)]);
            else
                fprintf(fpo,"%g ", a[offset2d(j, i, col)]);
            }
        }
    fclose(fpo);
}

int read_flight_input(double *Mach,double *alpha,int *i_TEL,int *i_TEU,double *rho,double *P,double *gamma,float *temp,float *temp1,float *temp2,float *temp3,float *temp4) {
    FILE *fpi = NULL;
    if (NULL == (fpi = fopen("flight.txt","rt"))) {
        printf("Something is worng\n");
        fclose(fpi);
        return ERROR;
    }
    fscanf(fpi,"%g %g %d %d %g %g %g",temp,temp1,i_TEL,i_TEU,temp2,temp3,temp4);
    *Mach = (double)(*temp);
    *alpha = (double)(*temp1);
    *rho = (double)(*temp2);
    *P = (double)(*temp3);
    *gamma = (double)(*temp4);

    fclose(fpi);
    return OK;
}

int offset2d(const int i,const int j,const int ni)
{
    return j * ni + i;
}

int offset3d(const int i,const int j,const int k,const int ni,const int nj)
{
    return (k*nj+j) * ni + i;
}

int intial_Q(const int ni,const int nj,double *q,const double rho,const double Mach,const double gamma,const double P,const double alpha) {
    const double a = sqrt(gamma*P/rho); // speed of sound
    for(int i=0;i<ni;i++) {
        for(int j=0;j<nj;j++) {
            q[offset3d(i,j,0,ni,nj)] = rho; //density
            q[offset3d(i,j,1,ni,nj)] = rho*Mach*a*cos(alpha); //momentom x
            q[offset3d(i,j,2,ni,nj)] = rho*Mach*a*sin(alpha); //momentom y
            q[offset3d(i,j,3,ni,nj)] = P/(gamma-1) + 0.5*(pow(q[offset3d(i,j,1,ni,nj)],2)+pow(q[offset3d(i,j,2,ni,nj)],2))/q[offset3d(i,j,0,ni,nj)];
        }
    }
    return OK;
}

int jacobian(const int ni,const int nj,double *X,double *Y,double *J,double *xi_x,double *xi_y,double *eta_x,double *eta_y) { //maybe change it to second order
    double x_i ,x_j ,y_i ,y_j;
    for(int i=0;i<ni;i++) {
        for(int j=0;j<nj;j++) {
            const int index = offset2d(i,j,ni);
            if(i==0) { //forword difference for i
                x_i = X[offset2d(i+1,j,ni)]-X[offset2d(i,j,ni)]; //x_xi
                y_i = Y[offset2d(i+1,j,ni)]-Y[offset2d(i,j,ni)]; //y_xi
            } else if(i==ni-1) { //backword difference for i
                x_i = X[offset2d(i,j,ni)]-X[offset2d(i-1,j,ni)]; //x_xi
                y_i = Y[offset2d(i,j,ni)]-Y[offset2d(i-1,j,ni)]; //y_xi
            } else { // central diffrence
                x_i = 0.5*(X[offset2d(i+1,j,ni)]-X[offset2d(i-1,j,ni)]);
                y_i = 0.5*(Y[offset2d(i+1,j,ni)]-Y[offset2d(i-1,j,ni)]);
            }
            if(j==0) { //forword difference for j
                x_j = X[offset2d(i,j+1,ni)]-X[offset2d(i,j,ni)]; //x_eta
                y_j = Y[offset2d(i,j+1,ni)]-Y[offset2d(i,j,ni)]; //y_eta
            } else if(j==nj-1) { //backword difference for j
                x_j = X[offset2d(i,j,ni)]-X[offset2d(i,j-1,ni)]; //x_eta
                y_j = Y[offset2d(i,j,ni)]-Y[offset2d(i,j-1,ni)]; //y_eta
            } else { // central diffrence
                x_j = 0.5*(X[offset2d(i,j+1,ni)]-X[offset2d(i,j-1,ni)]);
                y_j = 0.5*(Y[offset2d(i,j+1,ni)]-Y[offset2d(i,j-1,ni)]);
            }
            J[index] = 1.0/(x_i*y_j-x_j*y_i);
            xi_x[index] = J[index]*y_j;
            xi_y[index] = -J[index]*x_j;
            eta_x[index] = -J[index]*y_i;
            eta_y[index] = J[index]*x_i;
        }
    }
    return OK;
}

int RHS(const int ni,const int nj,double *S,double *W, double *q, double *J, double *xi_x, double *xi_y, double *eta_x,
    double *eta_y, const double gamma,const double dt) {

    for(int i=0;i<ni;i++)
        for(int j=0;j<nj;j++)
            for(int k=0;k<4;k++)
                S[offset3d(i,j,k,ni,nj)]=0;
    /*xi direction*/
    for(int j=1;j<nj-1;j++) {
        for(int i=0;i<ni;i++) {
            const int index = offset2d(i,j,ni);
            const double q0 = q[offset3d(i,j,0,ni,nj)];
            const double q1 = q[offset3d(i,j,1,ni,nj)];
            const double q2 = q[offset3d(i,j,2,ni,nj)];
            const double q3 = q[offset3d(i,j,3,ni,nj)];
            const double p = (gamma-1)*(q3-0.5*(pow(q1,2)+pow(q2,2))/q0);
            const double U = xi_x[index]*(q1/q0) + xi_y[index]*(q2/q0);

            W[offset2d(i,0,ni)]=q0*U/J[index]; /*continuity equation*/
            W[offset2d(i,1,ni)]=(q1*U+xi_x[index]*p)/J[index];/*X momentum equation*/
            W[offset2d(i,2,ni)]=(q2*U+xi_y[index]*p)/J[index];/*Y momentum equation*/
            W[offset2d(i,3,ni)]=(q3+p)*U/J[index];/*Energy equation*/
        }
        for(int i=1;i<ni-1;i++) {
            for(int k=0;k<4;k++) {
                S[offset3d(i,j,k,ni,nj)] += -0.5*dt*(W[offset2d(i+1,k,ni)]-W[offset2d(i-1,k,ni)]);
            }
        }
    }
    /*eta direction*/
    for(int i=1;i<ni-1;i++) {
        for(int j=0;j<nj;j++) {
            const int index = offset2d(i,j,ni);
            const double q0 = q[offset3d(i,j,0,ni,nj)];
            const double q1 = q[offset3d(i,j,1,ni,nj)];
            const double q2 = q[offset3d(i,j,2,ni,nj)];
            const double q3 = q[offset3d(i,j,3,ni,nj)];
            const double p = (gamma-1)*(q3-0.5*(pow(q1,2)+pow(q2,2))/q0);
            const double V = eta_x[index]*(q1/q0) + eta_y[index]*(q2/q0);

            W[offset2d(j,0,nj)]=q0*V/J[index];/*continuity equation*/
            W[offset2d(j,1,nj)]=(q1*V+eta_x[index]*p)/J[index];/*X momentum equation*/
            W[offset2d(j,2,nj)]=(q2*V+eta_y[index]*p)/J[index];/*Y momentum equation*/
            W[offset2d(j,3,nj)]=(q3+p)*V/J[index];/*Energy equation*/
        }
        for(int j=1;j<nj-1;j++) {
            for(int k=0;k<4;k++) {
                S[offset3d(i,j,k,ni,nj)] += -0.5*dt*(W[offset2d(j+1,k,nj)]-W[offset2d(j-1,k,nj)]);
            }
        }
    }
    return OK;
}

int smooth(double *q, double *s, double *jac, double *xx, double *xy,
           double *yx, double *yy, int id, int jd, double *s2,
           double *rspec, double *qv, double *dd,
           double epse, double gamma, double fsmach, double dt)
{

    double *rho, *u_vel, *v_vel, *t_e;

    double eratio, smool, gm1, ggm1, cx, cy, eps, ra, u, v, qq, ss, st,
          qav, qxx, ssfs, qyy;
    int ib, ie, jb, je, i, j, offset, offsetp1, offsetm1, ip, ir, n,
        jp, jr;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001,gamma);
    smool = 1.0;
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    ib = 1;
    ie = id - 1;
    jb = 1;
    je = jd - 1;

    cx = 2.;
    cy = 1.;

    rho = q;
    u_vel = &q[id*jd];
    v_vel = &q[2*id*jd];
    t_e = &q[3*id*jd];

/*     smoothing in xi direction */

    for (j = jb; j < je; j++) {
	for (i = 0; i < id; i++) {
	    offset = id * j + i;
	    eps = epse / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq);
            rspec[i] = eps * (fabs(xx[offset] * u + xy[offset] * v) + sqrt((xx[offset] * xx[offset] + xy[offset] * xy[offset]) * ss + 0.01));
	    qv[i] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (i = ib; i < ie; i++) {
	    ip = i + 1;
	    ir = i - 1;
	    qxx = qv[ip] - qv[i] * 2. + qv[ir];
	    qav = (qv[ip] + qv[i] * 2. + qv[ir]) * .25;
	    dd[i] = eratio * fabs(qxx / qav);
	}

	dd[0] = dd[1];
	dd[id - 1] = dd[id - 2];

	for (n = 0; n < 4; n++) {
	    for (i = ib; i < ie; i++) {
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j) * id + i + 1;
		offsetm1 = (jd * n + j) * id + i - 1;
		s2[i] = q[offsetp1] - 2.0 * q[offset] + q[offsetm1];
	    }

	    s2[0] = s2[1] * -1.;
	    s2[id - 1] = s2[id - 2] * -1.;

	    for (i = ib; i < ie; i++) {
		ip = i + 1;
		ir = i - 1;
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j) * id + i + 1;
		offsetm1 = (jd * n + j) * id + i - 1;
		st = ((dd[ip] + dd[i]) * .5 * (q[offsetp1] - q[offset]) - cx / (cx + dd[ip] + dd[i]) * (s2[ip] - s2[i])) * (rspec[ip] + rspec[i]) + ((dd[ir] + dd[i]) * .5 * (q[offsetm1] - q[offset]) - cx / (cx + dd[ir] + dd[i]) * (s2[ir] - s2[i])) * (rspec[ir] + rspec[i]);
		s[offset] += st * .5 * dt;
	    }
	}
    }

/*     smoothing in eta direction */

    ssfs = 1. / (0.001 + fsmach * fsmach);
    for (i = ib; i < ie; i++) {
	for (j = 0; j < jd; j++) {
	    offset = id * j + i;
	    eps = epse / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq) * (1.0 - smool) + smool * qq * ssfs;
            rspec[j] = eps * (fabs(yx[offset] * u + yy[offset] * v) + sqrt((yx[offset] * yx[offset] + yy[offset] * yy[offset]) * ss + 0.01));
	    qv[j] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (j = jb; j < je; j++) {
	    jp = j + 1;
	    jr = j - 1;
	    qyy = qv[jp] - qv[j] * 2. + qv[jr];
	    qav = (qv[jp] + qv[j] * 2. + qv[jr]) * .25;
	    dd[j] = eratio * fabs(qyy / qav);
	}

	dd[0] = dd[1];
	dd[jd - 1] = dd[jd - 2];

	for (n = 0; n < 4; n++) {
	    for (j = jb; j < je; j++) {
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j + 1) * id + i;
		offsetm1 = (jd * n + j - 1) * id + i;
		s2[j] = q[offsetp1] - 2.0 * q[offset] + q[offsetm1];
	    }

	    s2[0] = s2[1] * -1.;
	    s2[jd - 1] = s2[jd - 2] * -1.;

	    for (j = jb; j < je; j++) {
		jp = j + 1;
		jr = j - 1;
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j + 1) * id + i;
		offsetm1 = (jd * n + j - 1) * id + i;
		st = ((dd[jp] + dd[j]) * .5 * (q[offsetp1] - q[offset]) - cy / (cy + dd[jp] + dd[j]) * (s2[jp] - s2[j])) * (rspec[jp] + rspec[j]) + ((dd[jr] + dd[j]) * .5 * (q[offsetm1] - q[offset]) - cy / (cy + dd[jr] + dd[j]) * (s2[jr] - s2[j])) * (rspec[jr] + rspec[j]);
		s[offset] += st * .5 * dt;
	    }
	}
    }

    return 0;
} /* smooth */

int update_Q(const int ni,const int nj,double *S,double *q,double *J) {
    for(int i=1;i<ni-1;i++)
        for(int j=1;j<nj-1;j++)
            for(int k=0;k<4;k++) {
                q[offset3d(i,j,k,ni,nj)] +=S[offset3d(i,j,k,ni,nj)]*J[offset2d(i,j,ni)];
                // if(j==1&&k==0) {
                //     printf("i=%2d,j=%2d,%g\n",i,j,q[offset3d(i,j,k,ni,nj)]);
                // }
            }
    return OK;
}

int BC(const int ni,const int nj, double *q,const int i_TEL,const int i_TEU,const double gamma, double *eta_x, double *eta_y,
            double *J, double *xi_x, double *xi_y) {
    /*wall - no penetration*/
    for(int i=i_TEL-1;i<i_TEU;i++) { //including TEL,TEU i.e parmaters for trailing edge
        const double u1 = q[offset3d(i,1,1,ni,nj)]/q[offset3d(i,1,0,ni,nj)];
        const double v1 = q[offset3d(i,1,2,ni,nj)]/q[offset3d(i,1,0,ni,nj)];
        const double e1 = q[offset3d(i,1,3,ni,nj)];
        const double rho1 = q[offset3d(i,1,0,ni,nj)];
        const double p0 = (gamma-1)*(e1-0.5*rho1*(u1*u1+v1*v1));
        const double U1 = xi_x[offset2d(i,1,ni)]*u1+xi_y[offset2d(i,1,ni)]*v1;
        const double u0 = (eta_y[offset2d(i,0,ni)]*U1)/J[offset2d(i,0,ni)];
        const double v0 = (-eta_x[offset2d(i,0,ni)]*U1)/J[offset2d(i,0,ni)];
        const double rho0 = rho1;
        q[offset3d(i,0,0,ni,nj)] = rho0; //rho
        q[offset3d(i,0,1,ni,nj)] = rho0*u0; //momentum x
        q[offset3d(i,0,2,ni,nj)] = rho0*v0; //momentum y
        q[offset3d(i,0,3,ni,nj)] = p0/(gamma-1)+0.5*rho0*(u0*u0+v0*v0); //energy
        // printf("%g , index - %d\n",p0,offset2d(i,0,ni));
    }
    /*Trailing edge*/
    const double u_teu = q[offset3d(i_TEU-1,0,1,ni,nj)]/q[offset3d(i_TEU-1,0,0,ni,nj)];
    const double v_teu = q[offset3d(i_TEU-1,0,2,ni,nj)]/q[offset3d(i_TEU-1,0,0,ni,nj)];
    const double p_teu = (gamma-1)*(q[offset3d(i_TEU-1,0,3,ni,nj)]-0.5*q[offset3d(i_TEU-1,0,0,ni,nj)]*(u_teu*u_teu + v_teu*v_teu));

    const double u_tel = q[offset3d(i_TEL-1,0,1,ni,nj)]/q[offset3d(i_TEL-1,0,0,ni,nj)];
    const double v_tel = q[offset3d(i_TEL-1,0,2,ni,nj)]/q[offset3d(i_TEL-1,0,0,ni,nj)];
    const double p_tel = (gamma-1)*(q[offset3d(i_TEL-1,0,3,ni,nj)]-0.5*q[offset3d(i_TEL-1,0,0,ni,nj)]*(u_tel*u_tel + v_tel*v_tel));

    const double rho_te = 0.5*(q[offset3d(i_TEU-1,0,0,ni,nj)]+q[offset3d(i_TEL-1,0,0,ni,nj)]); // average density
    const double u_te = 0.5*(u_teu+u_tel); // average u
    const double v_te = 0.5*(v_teu+v_tel); // average v
    const double p_te = 0.5*(p_teu+p_tel); // average p
    /*setting te at teu and tel*/
    q[offset3d(i_TEU-1,0,0,ni,nj)] = rho_te;
    q[offset3d(i_TEU-1,0,1,ni,nj)] = rho_te*u_te;
    q[offset3d(i_TEU-1,0,2,ni,nj)] = rho_te*v_te;
    q[offset3d(i_TEU-1,0,3,ni,nj)] = p_te/(gamma-1)+0.5*rho_te*(u_te*u_te+v_te*v_te);

    q[offset3d(i_TEL-1,0,0,ni,nj)] = q[offset3d(i_TEU-1,0,0,ni,nj)];
    q[offset3d(i_TEL-1,0,1,ni,nj)] = q[offset3d(i_TEU-1,0,1,ni,nj)];
    q[offset3d(i_TEL-1,0,2,ni,nj)] = q[offset3d(i_TEU-1,0,2,ni,nj)];
    q[offset3d(i_TEL-1,0,3,ni,nj)] = q[offset3d(i_TEU-1,0,3,ni,nj)];

    /*wake BC*/
    for(int i=i_TEU;i<ni;i++)
        for(int k=0;k<4;k++) {
            q[offset3d(i,0,k,ni,nj)] = 0.5*(q[offset3d(i,1,k,ni,nj)]+q[offset3d(ni-i-1,1,k,ni,nj)]);
            q[offset3d(ni-i-1,0,k,ni,nj)]=q[offset3d(i,0,k,ni,nj)];
        }

    /*outflow BC*/
    for(int j=0;j<nj;j++)
        for(int k=0;k<4;k++) {
            q[offset3d(0,j,k,ni,nj)]=q[offset3d(1,j,k,ni,nj)];
            q[offset3d(ni-1,j,k,ni,nj)]=q[offset3d(ni-2,j,k,ni,nj)];
        }
    return OK;
}

int LHSX(const int j,const int ni,const int nj,double *A,double *B,double *C,double *q, double gamma,double *xi_x,
    double *xi_y, double dt) {
    // xi direaction
    const int n_max = (int)fmax(ni,nj);
    for(int i=0;i<ni;i++) {
        const double u = q[offset3d(i,j,1,ni,nj)]/q[offset3d(i,j,0,ni,nj)];
        const double v = q[offset3d(i,j,2,ni,nj)]/q[offset3d(i,j,0,ni,nj)];
        const double psi_sqrd = 0.5*(gamma-1)*(u*u + v*v);
        const double k_x = xi_x[offset2d(i,j,ni)];
        const double k_y = xi_y[offset2d(i,j,ni)];
        const double theta = k_x*u+k_y*v;
        const double gamma1 = gamma-1;
        const double gamma2 = gamma-2;
        const double e = q[offset3d(i,j,3,ni,nj)];
        const double rho = q[offset3d(i,j,0,ni,nj)];
        const double beta = gamma*e/rho-psi_sqrd;

        B[offset3d(i,0,0,n_max,4)]=0; // k_t
        B[offset3d(i,0,1,n_max,4)]=k_x;
        B[offset3d(i,0,2,n_max,4)]=k_y;
        B[offset3d(i,0,3,n_max,4)]=0;

        B[offset3d(i,1,0,n_max,4)]=k_x*psi_sqrd-u*theta;
        B[offset3d(i,1,1,n_max,4)]=theta-k_x*gamma2*u;
        B[offset3d(i,1,2,n_max,4)]=k_y*u-gamma1*k_x*v;
        B[offset3d(i,1,3,n_max,4)]=k_x*gamma1;

        B[offset3d(i,2,0,n_max,4)]=k_y*psi_sqrd-v*theta;
        B[offset3d(i,2,1,n_max,4)]=k_x*v-k_y*gamma1*u;
        B[offset3d(i,2,2,n_max,4)]=theta-k_y*gamma2*v;
        B[offset3d(i,2,3,n_max,4)]=k_y*gamma1;

        B[offset3d(i,3,0,n_max,4)]=theta*(2*psi_sqrd-gamma*e/rho);
        B[offset3d(i,3,1,n_max,4)]=k_x*beta-gamma1*u*theta;
        B[offset3d(i,3,2,n_max,4)]=k_y*beta-gamma1*v*theta;
        B[offset3d(i,3,3,n_max,4)]=gamma*theta;
    }
    for(int n=0;n<4;n++)
        for(int m=0;m<4;m++)
            for(int i=1;i<ni-1;i++) {
                A[offset3d(i,m,n,n_max,4)]=-0.5*B[offset3d(i-1,m,n,n_max,4)]*dt;
                C[offset3d(i,m,n,n_max,4)]=0.5*B[offset3d(i-1,m,n,n_max,4)]*dt;
            }
    for(int i=1;i<ni-1;i++)
        for(int n=0;n<4;n++)
            for(int m=0;m<4;m++) {
                if(n==m)
                    B[offset3d(i,n,m,n_max,4)]=1;
                else
                    B[offset3d(i,n,m,n_max,4)]=0;
            }
    return OK;
}

int smoothx(double *q, double *xx, double *xy, int id, int jd, double *a,double *b, double *c, int j,double *jac, double *drr,
    double *drp,double *rspec, double *qv, double *dd,double epsi, double gamma, double fsmach, double dt)
{

    double *rho, *u_vel, *v_vel, *t_e;

    double eratio, gm1, ggm1, eps, ra, u, v, qq, ss,
          qav, qxx, rr, rp;
    int ib, ie, i, offset, offsetp1, offsetm1, ip, ir, n;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001,gamma);
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    ib = 1;
    ie = id - 1;

    rho = q;
    u_vel = &q[id*jd];
    v_vel = &q[2*id*jd];
    t_e = &q[3*id*jd];

/*     smoothing in xi direction */

        for (i = 0; i < id; i++) {
	    offset = id * j + i;
	    eps = epsi / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq);
            rspec[i] = eps * (fabs(xx[offset] * u + xy[offset] * v) + sqrt((xx[offset] * xx[offset] + xy[offset] * xy[offset]) * ss + 0.01));
	    qv[i] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (i = ib; i < ie; i++) {
	    ip = i + 1;
	    ir = i - 1;
	    qxx = qv[ip] - qv[i] * 2. + qv[ir];
	    qav = (qv[ip] + qv[i] * 2. + qv[ir]) * .25;
	    dd[i] = eratio * fabs(qxx / qav);
	}

	dd[0] = dd[1];
	dd[id - 1] = dd[id - 2];

        for (i = ib; i < ie; i++) {
	    ip = i + 1;
	    ir = i - 1;
	    offset = j * id + i;
	    offsetp1 = j * id + i + 1;
	    offsetm1 = j * id + i - 1;
	    rp = (0.5 * (dd[ip] + dd[i]) + 2.5) * dt * 0.5 * (rspec[ip] + rspec[i]);
	    rr = (0.5 * (dd[ir] + dd[i]) + 2.5) * dt * 0.5 * (rspec[ir] + rspec[i]);
	    qv[i] = (rr + rp) * jac[offset];
	    drr[i] = rr * jac[offsetm1];
	    drp[i] = rp * jac[offsetp1];
	}

	for (n = 0; n < 4; n++) {
	    for (i = ib; i < ie; i++) {
	        offset = (n * 4 + n) * id + i;
		a[offset] -= drr[i];
		b[offset] += qv[i];
		c[offset] -= drp[i];
	    }
        }
    return 0;
} /* smoothx */

int LHSY(const int i,const int ni,const int nj,double *A,double *B,double *C,double *q, double gamma,double *eta_x,
    double *eta_y, double dt) {
    // eta direaction
    const int n_max = (int)fmax(ni,nj);
    for(int j=0;j<nj;j++) {
        const double u = q[offset3d(i,j,1,ni,nj)]/q[offset3d(i,j,0,ni,nj)];
        const double v = q[offset3d(i,j,2,ni,nj)]/q[offset3d(i,j,0,ni,nj)];
        const double psi_sqrd = 0.5*(gamma-1)*(u*u + v*v);
        const double k_x = eta_x[offset2d(i,j,ni)];
        const double k_y = eta_y[offset2d(i,j,ni)];
        const double theta = k_x*u+k_y*v;
        const double gamma1 = gamma-1;
        const double gamma2 = gamma-2;
        const double e = q[offset3d(i,j,3,ni,nj)];
        const double rho = q[offset3d(i,j,0,ni,nj)];
        const double beta = gamma*e/rho-psi_sqrd;

        B[offset3d(j,0,0,n_max,4)]=0; // k_t
        B[offset3d(j,0,1,n_max,4)]=k_x;
        B[offset3d(j,0,2,n_max,4)]=k_y;
        B[offset3d(j,0,3,n_max,4)]=0;

        B[offset3d(j,1,0,n_max,4)]=k_x*psi_sqrd-u*theta;
        B[offset3d(j,1,1,n_max,4)]=theta-k_x*gamma2*u;
        B[offset3d(j,1,2,n_max,4)]=k_y*u-gamma1*k_x*v;
        B[offset3d(j,1,3,n_max,4)]=k_x*gamma1;

        B[offset3d(j,2,0,n_max,4)]=k_y*psi_sqrd-v*theta;
        B[offset3d(j,2,1,n_max,4)]=k_x*v-k_y*gamma1*u;
        B[offset3d(j,2,2,n_max,4)]=theta-k_y*gamma2*v;
        B[offset3d(j,2,3,n_max,4)]=k_y*gamma1;

        B[offset3d(j,3,0,n_max,4)]=theta*(2*psi_sqrd-gamma*e/rho);
        B[offset3d(j,3,1,n_max,4)]=k_x*beta-gamma1*u*theta;
        B[offset3d(j,3,2,n_max,4)]=k_y*beta-gamma1*v*theta;
        B[offset3d(j,3,3,n_max,4)]=gamma*theta;
    }
    for(int n=0;n<4;n++)
        for(int m=0;m<4;m++)
            for(int j=1;j<nj-1;j++) {
                A[offset3d(j,m,n,n_max,4)]=-0.5*B[offset3d(i-1,m,n,n_max,4)]*dt;
                C[offset3d(j,m,n,n_max,4)]=0.5*B[offset3d(i-1,m,n,n_max,4)]*dt;
            }
    for(int j=1;j<nj-1;j++)
        for(int n=0;n<4;n++)
            for(int m=0;m<4;m++) {
                if(n==m)
                    B[offset3d(j,n,m,n_max,4)]=1;
                else
                    B[offset3d(j,n,m,n_max,4)]=0;
            }
    return OK;
}

int smoothy(double *q, double *yx, double *yy, int id, int jd, double *a,double *b, double *c, int i,double *jac, double *drr,
    double *drp,double *rspec, double *qv, double *dd,double epsi, double gamma, double fsmach, double dt)
{

    double *rho, *u_vel, *v_vel, *t_e;

    double eratio, smool, gm1, ggm1, eps, ra, u, v, qq, ss,
          qav, ssfs, qyy, rp, rr;
    int jb, je, j, offset, offsetp1, offsetm1, n,
        jp, jr;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001,gamma);
    smool = 1.0;
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    jb = 1;
    je = jd - 1;

    rho = q;
    u_vel = &q[id*jd];
    v_vel = &q[2*id*jd];
    t_e = &q[3*id*jd];

/*     smoothing in eta direction */

        ssfs = 1. / (0.001 + fsmach * fsmach);
        for (j = 0; j < jd; j++) {
	    offset = id * j + i;
	    eps = epsi / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq) * (1.0 - smool) + smool * qq * ssfs;
            rspec[j] = eps * (fabs(yx[offset] * u + yy[offset] * v) + sqrt((yx[offset] * yx[offset] + yy[offset] * yy[offset]) * ss + 0.01));
	    qv[j] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (j = jb; j < je; j++) {
	    jp = j + 1;
	    jr = j - 1;
	    qyy = qv[jp] - qv[j] * 2. + qv[jr];
	    qav = (qv[jp] + qv[j] * 2. + qv[jr]) * .25;
	    dd[j] = eratio * fabs(qyy / qav);
	}

	dd[0] = dd[1];
	dd[jd - 1] = dd[jd - 2];

        for (j = jb; j < je; j++) {
	    jp = j + 1;
	    jr = j - 1;
	    offset = j * id + i;
	    offsetp1 = (j + 1) * id + i;
	    offsetm1 = (j - 1) * id + i;
	    rp = (0.5 * (dd[jp] + dd[j]) + 2.5) * dt * 0.5 * (rspec[jp] + rspec[j]);
	    rr = (0.5 * (dd[jr] + dd[j]) + 2.5) * dt * 0.5 * (rspec[jr] + rspec[j]);
	    qv[j] = (rr + rp) * jac[offset];
	    drr[j] = rr * jac[offsetm1];
	    drp[j] = rp * jac[offsetp1];
	}

	for (n = 0; n < 4; n++) {
	    for (j = jb; j < je; j++) {
	        offset = (n * 4 + n) * jd + j;
		a[offset] -= drr[j];
		b[offset] += qv[j];
		c[offset] -= drp[j];
	    }
        }
    return 0;
} /* smoothy */

int btri4s(double *a, double *b, double *c, double *f, int kd, int ks, int ke)
{
  /* Local variables */
  int k, m, n, nd, md;

  double c1, d1, d2, d3, d4, c2, c3, c4, b11, b21, b22, b31, b32, b33,
    b41, b42, b43, b44, u12, u13, u14, u23, u24, u34;


  /*   (A,B,C)F = F, F and B are overloaded, solution in F */

  md = 4;
  nd = 4;

  /*   Part 1. Forward block sweep */

  for (k = ks; k <= ke; k++)
    {

      /*      Step 1. Construct L in B */

      if (k != ks)
	{
	  for (m = 0; m < md; m++)
	    {
	      for (n = 0; n < nd; n++)
		{
		  b[k + kd * (m + md * n)] = b[k + kd * (m + md * n)]
		    - a[k + kd * (m + md * 0)] * b[k - 1 + kd * (0 + md * n)]
		    - a[k + kd * (m + md * 1)] * b[k - 1 + kd * (1 + md * n)]
		    - a[k + kd * (m + md * 2)] * b[k - 1 + kd * (2 + md * n)]
		    - a[k + kd * (m + md * 3)] * b[k - 1 + kd * (3 + md * n)] ;
		}
	    }
	}

      /*      Step 2. Compute L inverse (block matrix) */

      /*          A. Decompose L into L and U */

      b11 = 1. / b[k + kd * (0 + md * 0)];
      u12 = b[k + kd * (0 + md * 1)] * b11;
      u13 = b[k + kd * (0 + md * 2)] * b11;
      u14 = b[k + kd * (0 + md * 3)] * b11;
      b21 = b[k + kd * (1 + md * 0)];
      b22 = 1. / (b[k + kd * (1 + md * 1)] - b21 * u12);
      u23 = (b[k + kd * (1 + md * 2)] - b21 * u13) * b22;
      u24 = (b[k + kd * (1 + md * 3)] - b21 * u14) * b22;
      b31 = b[k + kd * (2 + md * 0)];
      b32 = b[k + kd * (2 + md * 1)] - b31 * u12;
      b33 = 1. / (b[k + kd * (2 + md * 2)] - b31 * u13 - b32 * u23);
      u34 = (b[k + kd * (2 + md * 3)] - b31 * u14 - b32 * u24) * b33;
      b41 = b[k + kd * (3 + md * 0)];
      b42 = b[k + kd * (3 + md * 1)] - b41 * u12;
      b43 = b[k + kd * (3 + md * 2)] - b41 * u13 - b42 * u23;
      b44 = 1. / (b[k + kd * (3 + md * 3)] - b41 * u14 - b42 * u24
		  - b43 * u34);

      /*      Step 3. Solve for intermediate vector */

      /*          A. Construct RHS */
      if (k != ks)
	{
	  for (m = 0; m < md; m++)
	    {
	      f[k + kd * m] = f[k + kd * m]
		- a[k + kd * (m + md * 0)] * f[k - 1 + kd * 0]
		- a[k + kd * (m + md * 1)] * f[k - 1 + kd * 1]
		- a[k + kd * (m + md * 2)] * f[k - 1 + kd * 2]
		- a[k + kd * (m + md * 3)] * f[k - 1 + kd * 3];
	    }
	}

      /*          B. Intermediate vector */

      /*          Forward substitution */

      d1 = f[k + kd * 0] * b11;
      d2 = (f[k + kd * 1] - b21 * d1) * b22;
      d3 = (f[k + kd * 2] - b31 * d1 - b32 * d2) * b33;
      d4 = (f[k + kd * 3] - b41 * d1 - b42 * d2 - b43 * d3) * b44;

      /*          Backward substitution */

      f[k + kd * 3] = d4;
      f[k + kd * 2] = d3 - u34 * d4;
      f[k + kd * 1] = d2 - u23 * f[k + kd * 2] - u24 * d4;
      f[k + kd * 0] = d1 - u12 * f[k + kd * 1] - u13 * f[k + kd * 2] - u14 * d4;

      /*      Step 4. Construct U = L ** (-1) * C */
      /*              by columns and store in B */

      if (k != ke)
	{
	  for (n = 0; n < nd; n++)
	    {

	      /*          Forward substitution */

	      c1 = c[k + kd * (0 + md * n)] * b11;
	      c2 = (c[k + kd * (1 + md * n)] - b21 * c1) * b22;
	      c3 = (c[k + kd * (2 + md * n)] - b31 * c1 - b32 * c2) *
		b33;
	      c4 = (c[k + kd * (3 + md * n)] - b41 * c1 - b42 * c2 -
		    b43 * c3) * b44;

	      /*          Backward substitution */

	      b[k + kd * (3 + md * n)] = c4;
	      b[k + kd * (2 + md * n)] = c3 - u34 * c4;
	      b[k + kd * (1 + md * n)] = c2 - u23 * b[k + kd * (2 + md * n)] - u24 * c4;
	      b[k + kd * (0 + md * n)] = c1 - u12 * b[k + kd * (1 + md * n)]
		- u13 * b[k + kd * (2 + md * n)] - u14 * c4;
	    }
	}
    }

  /*   Part 2. Backward block sweep */

  if (ke == ks)
    {
      return 0;
    }

  for (k = ke - 1; k >= ks; --k)
    {
      for (m = 0; m < md; m++)
	{
	  f[k + kd * m] = f[k + kd * m]
	    - b[k + kd * (m + md * 0)] * f[k + 1 + kd * 0]
	    - b[k + kd * (m + md * 1)] * f[k + 1 + kd * 1]
	    - b[k + kd * (m + md * 2)] * f[k + 1 + kd * 2]
	    - b[k + kd * (m + md * 3)] * f[k + 1 + kd * 3];
	}
    }

  return 0;

} /* btri4s_ */

int STEP(const int ni,const int nj,double *S,double *W, double *q, double *J, double *xi_x, double *xi_y, double *eta_x,
    double *eta_y,const double gamma,const double dt,double *s2,double *rspec,double *qv,double *dd,const double epse
    , const double Mach,double *A, double *B, double *C, double *D, double *drr, double *drp,const double epsi
    ,const int i_TEL,const int i_TEU) {

    double L2_max_1=0,L2_max_end=0.1;
    int runs = 0;
    FILE *fpo=NULL;
    fpo = fopen("L2.txt","wt");
    fclose(fpo);
    fpo = fopen("L2.txt","a");

    while(L2_max_end>0.0001*L2_max_1) {
        double S_max=0;
        // BC(ni,nj,q,i_TEL,i_TEU,gamma,eta_x,eta_y,J,xi_x,xi_y);
        RHS(ni,nj,S,W,q,J,xi_x,xi_y,eta_x,eta_y,gamma,dt);
        smooth(q,S,J,xi_x,xi_y,eta_x,eta_y,ni,nj,s2,rspec,qv,dd,epse,gamma,Mach,dt);

        // //xi inversions
        // for(int j=1;j<nj-1;j++) {
        //     LHSX(j,ni,nj,A,B,C,q,gamma,xi_x,xi_y,dt);
        //     smoothx(q,xi_x,xi_y,ni,nj,A,B,C,j,J,drr,drp,rspec,qv,dd,epsi,gamma,Mach,dt);
        //     for(int i=0;i<ni;i++)
        //         for(int k=0;k<4;k++)
        //             D[offset2d(i,k,ni)]=S[offset3d(i,j,k,ni,nj)];
        //     btri4s(A,B,C,D,(int)fmax(ni,nj),1,ni-2);
        //     for(int i=0;i<ni;i++)
        //         for(int k=0;k<4;k++)
        //             S[offset3d(i,j,k,ni,nj)]=D[offset2d(i,k,ni)];
        // }
        // //eta inversions
        // for(int i=1;i<ni-1;i++) {
        //     LHSY(i,ni,nj,A,B,C,q,gamma,xi_x,xi_y,dt);
        //     smoothy(q,xi_x,xi_y,ni,nj,A,B,C,i,J,drr,drp,rspec,qv,dd,epsi,gamma,Mach,dt);
        //     for(int j=0;j<nj;j++)
        //         for(int k=0;k<4;k++)
        //             D[offset2d(j,k,ni)]=S[offset3d(i,j,k,ni,nj)];
        //     btri4s(A,B,C,D,(int)fmax(ni,nj),1,nj-2);
        //     for(int j=0;j<nj;j++)
        //         for(int k=0;k<4;k++)
        //             S[offset3d(i,j,k,ni,nj)]=D[offset2d(j,k,ni)];
        // }
        update_Q(ni,nj,S,q,J);
        for(int i=0;i<ni;i++)
            for(int j=0;j<nj;j++)
                for(int k=0;k<4;k++) {
                    if(isnan(q[offset3d(i,j,k,ni,nj)])) {
                        printf("iter = %d , loc i=%d,j=%d,k=%d \n",runs+1,i,j,k);
                        output_sheet(q,ni,nj,0);
                        return ERROR;
                    }
                }

        for(int i=0;i<ni;i++)
            for(int j=0;j<nj;j++)
                for(int k=0;k<4;k++) {
                    const double value = sqrt(S[offset3d(i,j,k,ni,nj)]*S[offset3d(i,j,k,ni,nj)]);
                    if(S_max<value)
                        S_max=value;
                }
        // printf("L2 max iter %d iteration = %g \n",runs,L2_max_1);
        if(runs==0) {
            L2_max_1 = S_max;
            fprintf(fpo, "%0.15f\n",L2_max_1);
        }else {
            L2_max_end = S_max;
            fprintf(fpo, "%0.15f\n",L2_max_end);
        }
        if(runs%1000==0)
            printf("iter = %d\n",runs);
        if(runs==100000)
            return ERROR;
        runs++;
    }
    printf("L2 max first iteration = %g \n",L2_max_1);
    printf("L2 max last iteration = %g \n",L2_max_end);
    fclose(fpo);
    return  runs;
}

