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
int intial_EF(const int ni,const int nj,double *q,double *J,double *xi_x,double *xi_y,double *eta_x,double *eta_y,double *E,
    double *F,const double gamma);
int RHS(const int ni,const int nj,double *S,double *W,double *E,double *F,const double dt);
int smooth(double *q, double *s, double *jac, double *xx, double *xy,double *yx, double *yy, int id, int jd, double *s2,
    double *rspec, double *qv, double *dd,double epse, double gamma, double fsmach, double dt);
int update_Q(const int ni,const int nj,double *S,double *q,double *J);
int BC(const int ni,const int nj, double *q,const int i_TEL,const int i_TEU,const double gamma, double *eta_x, double *eta_y,
    double *J, double *xi_x, double *xi_y);

int main() {
    int nj, ni,i_TEL,i_TEU;
    double Mach,alpha,rho,P,gamma;
    float temp,temp1,temp2,temp3,temp4;
    const double dt = pow(10,-5),epse = 0.06;

    FILE *grid_file = fopen("grid.txt", "r");
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

    for(int i=0;i<5;i++) {
        // if((i+1)%100==0)
            printf("itaration %d \n",i+1);
        intial_EF(ni,nj,q,J,xi_x,xi_y,eta_x,eta_y,E,F,gamma);
        RHS(ni,nj,S,W,E,F,dt);
        smooth(q,S,J,xi_x,xi_y,eta_x,eta_y,ni,nj,s2,rspec,qv,dd,epse,gamma,Mach,dt);
        BC(ni,nj,q,i_TEL,i_TEU,gamma,eta_x,eta_y,J,xi_x,xi_y);
        update_Q(ni,nj,S,q,J);
    }
    output_matrix_3D(q,ni,nj);
    output_sheet(q,ni,nj,3);
    // output_matrix_2D(xi_x,ni,nj);


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
                    // fprintf(fpo,"%lf cell:%d \n", a[offset3d(i,j,k,ni,nj)],offset3d(i,j,k,ni,nj));
                        fprintf(fpo,"q:%d = %lf  \n", k,a[offset3d(i,j,k,ni,nj)]);
                else
                    // fprintf(fpo,"%lf cell:%d ", a[offset3d(i,j,k,ni,nj)],offset3d(i,j,k,ni,nj));
                        fprintf(fpo,"q:%d = %lf, ", k,a[offset3d(i,j,k,ni,nj)]);
            }
        }
    }
}

void output_sheet(const double *a,const int ni,const int nj,const int layer) {
    FILE *fpo = NULL;
    fpo = fopen("sheet.txt","wt");
    for (int j = 0; j < nj; j++) {
        for(int i=0;i<ni;i++) {
            if(i==ni-1)
                // fprintf(fpo,"%lf cell:%d \n", a[offset3d(i,j,k,ni,nj)],offset3d(i,j,k,ni,nj));
                fprintf(fpo,"%lf\n",a[offset3d(i,j,layer,ni,nj)]);
            else
                // fprintf(fpo,"%lf cell:%d ", a[offset3d(i,j,k,ni,nj)],offset3d(i,j,k,ni,nj));
                fprintf(fpo,"%lf ",a[offset3d(i,j,layer,ni,nj)]);

        }
    }
}

void output_matrix_2D(double *a,const int col,const int row){
    FILE *fpo=NULL;
    fpo = fopen("F_almog.txt","wt");
    for (int i = 0; i < row; i++) {
        for(int j=0;j<col;j++) {
            if(j==col-1)
                fprintf(fpo,"%lf\n", a[offset2d(j, i, col)]);
            else
                fprintf(fpo,"%lf ", a[offset2d(j, i, col)]);
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
            // printf("%g \n",xi_x[index]);
        }
    }
    return OK;
}

int intial_EF(const int ni, const int nj, double *q, double *J, double *xi_x, double *xi_y, double *eta_x,double *eta_y,
    double *E, double *F, const double gamma) {
    for(int i=0;i<ni;i++) {
        for(int j=0;j<nj;j++) {
            const int index = offset2d(i,j,ni);
            const double q0 = q[offset3d(i,j,0,ni,nj)];
            const double q1 = q[offset3d(i,j,1,ni,nj)];
            const double q2 = q[offset3d(i,j,2,ni,nj)];
            const double q3 = q[offset3d(i,j,3,ni,nj)];
            const double p = (gamma-1)*(q3-0.5*(pow(q1,2)+pow(q2,2))/q0);
            const double U = xi_x[index]*(q1/q0) + xi_y[index]*(q2/q0);
            const double V = eta_x[index]*(q1/q0) + eta_y[index]*(q2/q0);

            E[offset3d(i,j,0,ni,nj)] = q0*U/J[index];
            E[offset3d(i,j,1,ni,nj)] = (q1*U+xi_x[index]*p)/J[index];
            E[offset3d(i,j,2,ni,nj)] = (q2*U+xi_y[index]*p)/J[index];
            E[offset3d(i,j,3,ni,nj)] = (q3+p)*U/J[index];

            F[offset3d(i,j,0,ni,nj)] = q0*V/J[index];
            F[offset3d(i,j,1,ni,nj)] = (q1*V+eta_x[index]*p)/J[index];
            F[offset3d(i,j,2,ni,nj)] = (q2*V+eta_y[index]*p)/J[index];
            F[offset3d(i,j,3,ni,nj)] = (q3+p)*V/J[index];
        }
    }
    return OK;
}

int RHS(const int ni,const int nj,double *S,double *W,double *E,double *F,const double dt) {
    const int index = fmax(ni,nj);
    for(int i=0;i<ni;i++)
        for(int j=0;j<nj;j++)
            for(int k=0;k<4;k++)
                S[offset3d(i,j,k,ni,nj)]=0;
    /*xi direction*/
    for(int j=1;j<nj-1;j++) {
        for(int i=0;i<ni;i++) {
            W[offset2d(i,0,index)]=E[offset3d(i,j,0,ni,nj)];/*continuity equation*/
            W[offset2d(i,1,index)]=E[offset3d(i,j,1,ni,nj)];/*X momentum equation*/
            W[offset2d(i,2,index)]=E[offset3d(i,j,2,ni,nj)];/*Y momentum equation*/
            W[offset2d(i,3,index)]=E[offset3d(i,j,3,ni,nj)];/*Energy equation*/
        }
        for(int i=1;i<ni-1;i++) {
            for(int k=0;k<4;k++) {
                S[offset3d(i,j,k,ni,nj)] += -0.5*dt*(W[offset2d(i+1,k,index)]-W[offset2d(i-1,k,index)]);
            }
        }
    }
    /*eta direction*/
    for(int i=1;i<ni-1;i++) {
        for(int j=0;j<nj;j++) {
            W[offset2d(j,0,index)]=F[offset3d(i,j,0,ni,nj)];/*continuity equation*/
            W[offset2d(j,1,index)]=F[offset3d(i,j,1,ni,nj)];/*X momentum equation*/
            W[offset2d(j,2,index)]=F[offset3d(i,j,2,ni,nj)];/*Y momentum equation*/
            W[offset2d(j,3,index)]=F[offset3d(i,j,3,ni,nj)];/*Energy equation*/
        }
        for(int j=1;j<nj-1;j++) {
            for(int k=0;k<4;k++) {
                S[offset3d(i,j,k,ni,nj)] += -0.5*dt*(W[offset2d(j+1,k,index)]-W[offset2d(j-1,k,index)]);
            }
        }
    }
    return OK;
}


int smooth(double *q, double *s, double *jac, double *xx, double *xy,double *yx, double *yy, int id, int jd, double *s2,
    double *rspec, double *qv, double *dd,double epse, double gamma, double fsmach, double dt)
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
    for(int i=0;i<ni;i++)
        for(int j=0;j<nj;j++)
            for(int k=0;k<4;k++)
                q[offset3d(i,j,k,ni,nj)] +=S[offset3d(i,j,k,ni,nj)]*J[offset2d(i,j,ni)];
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
        const double u0 = eta_y[offset2d(i,0,ni)]*U1/J[offset2d(i,0,ni)];
        const double v0 = -eta_x[offset2d(i,0,ni)]*U1/J[offset2d(i,0,ni)];
        const double rho0 = rho1;
        q[offset3d(i,0,0,ni,nj)] = rho0; //rho
        q[offset3d(i,0,1,ni,nj)] = rho0*u0; //momentum x
        q[offset3d(i,0,2,ni,nj)] = rho0*v0; //momentum y
        q[offset3d(i,0,3,ni,nj)] = p0/(gamma-1)+0.5*rho0*(u0*u0+v0*v0); //energy
        printf("%g\n",p0);
    }
    /*Trailing edge*/
    const double u_teu = q[offset3d(i_TEU-1,0,1,ni,nj)]/q[offset3d(i_TEU-1,0,0,ni,nj)];
    const double v_teu = q[offset3d(i_TEU-1,0,2,ni,nj)]/q[offset3d(i_TEU-1,0,0,ni,nj)];
    const double p_teu = (gamma-1)*(q[offset3d(i_TEU-1,0,3,ni,nj)]-0.5*q[offset3d(i_TEU-1,0,0,ni,nj)]*(u_teu*u_teu+v_teu*v_teu));

    const double u_tel = q[offset3d(i_TEL-1,0,1,ni,nj)]/q[offset3d(i_TEL-1,0,0,ni,nj)];
    const double v_tel = q[offset3d(i_TEL-1,0,2,ni,nj)]/q[offset3d(i_TEL-1,0,0,ni,nj)];
    const double p_tel = (gamma-1)*(q[offset3d(i_TEL-1,0,3,ni,nj)]-0.5*q[offset3d(i_TEL-1,0,0,ni,nj)]*(u_tel*u_tel+v_tel*v_tel));

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
