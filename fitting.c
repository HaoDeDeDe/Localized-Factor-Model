//
//  fitting.c
//  MyProjectC
//
//  Created by 杨昊的 on 16/3/10.
//  Copyright © 2016年 杨昊的. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "fitting.h"
#include "matrixOperation.h"

//a single cell in a response matrix
struct response{
    int indicator;  //-1 for null，0 for predict, 1 for train
    double value;
};

//response matrix or feature matrix
struct resMatrix{
    int user_num;
    int item_num;
    double min_value;
    double max_value;
    struct response *matrix;
};

struct vector{
    int length;
    double *v;
};

struct matrix{
    int row_num;
    int col_num;
    double *m;
};

struct indexCell{
    int num;
    int *index;
};

//get uniform random variable
double var_uniform(double min, double max){
    //srand(time(0));
    return min+(max-min)*rand()/(RAND_MAX+1.0);
}

//density for Gaussian
double normal_pdf(double x, double mean, double std){
    if(0==std){
        std=0.000000001;
    }
    return exp(-(x-mean)*(x-mean)/(2*std*std))/(sqrt(2*3.141593)*std);
}

//generate Gaussian random variable
double var_normal(double min, double max, double mean, double std){
    double x, y, dScope;
    do{
        x=var_uniform(min, max);
        y=normal_pdf(x, mean, std);
        dScope=var_uniform(0, normal_pdf(mean, mean, std));
    }
    while(dScope>y);
    return x;
}

void paraInitialize(int context_num, struct resMatrix *context, double *alpha, double *sigma2_z, double *sigma2_y, double *sigma2_u, struct vector *beta, struct matrix *w, struct matrix *A, int *q, int *r){
    int temp, i, j, k;
    
    //initialization
    r[0]=3;
    for (k=0; k<context_num; k++) {
        q[k]=3;
    }
    
    temp=0;
    for (k=0; k<context_num; k++) {
        temp=temp+context[k].item_num;
    }
    
    for (i=0; i<temp; i++) {
        alpha[i]=var_normal(-0.3, 0.3, 0, 0.1);
    }
    
    for (k=0; k<context_num; k++) {
        sigma2_z[k]=var_normal(0.01, 0.3, 0.15, 0.05);
    }
    
    for (k=0; k<context_num; k++) {
        sigma2_y[k]=var_normal(0.01, 0.3, 0.15, 0.05);
    }
    
    sigma2_u[0]=var_normal(0.01, 0.3, 0.15, 0.05);

    temp=0;
    for (k=0; k<context_num; k++) {
        for (j=0; j<context[k].item_num; j++) {
            beta[temp+j].length=q[k];
            beta[temp+j].v=(double *)malloc(q[k]*sizeof(double));
            for (i=0; i<q[k]; i++) {
                beta[temp+j].v[i]=var_normal(-0.3, 0.3, 0, 0.1);
            }
        }
        temp=temp+context[k].item_num;
    }
    
    for (k=0; k<context_num; k++) {
        w[k].row_num=context[k].user_num;
        w[k].col_num=context[k].item_num;
        w[k].m=(double *)malloc(w[k].row_num*w[k].col_num*sizeof(double));
        for (i=0; i<w[k].row_num; i++) {
            for (j=0; j<w[k].col_num; j++) {
                w[k].m[i*w[k].col_num+j]=1;    
            }
        }
    }
    
    for (k=0; k<context_num; k++) {
        A[k].row_num=q[k];
        A[k].col_num=r[0];
        A[k].m=(double *)malloc(A[k].row_num*A[k].col_num*sizeof(double));
        for (i=0; i<A[k].row_num; i++) {
            for (j=0; j<A[k].col_num; j++) {
                A[k].m[i*A[k].col_num+j]=var_normal(-0.3, 0.3, 0, 0.1);
            }
        }
    }
    
    return;
}

void EStep(int context_num, struct resMatrix *context, double *alpha, double *sigma2_z, double *sigma2_y, double *sigma2_u, struct vector *beta, struct vector *Eu, struct vector *Ez, struct matrix *w, struct matrix *A, struct matrix *Vu, struct matrix *Vz, struct matrix *Covzu, int *q, int *r, struct indexCell *J){
    
    int i, j, k, t, temp;
    double * matrix1, * matrix2, * matrix3, * matrix4;
    
    //每个user一棵树，各树之间无关系，所以按user为单位进行循环计算
    for (i=0; i<context[0].user_num; i++) {
        
        //d is an array of vector with potentially different length
        //C, D, S, E, H and R is each an array of matrices with potentially different shapes
        struct vector * d=(struct vector *)malloc(context_num*sizeof(struct vector));
        struct matrix * C=(struct matrix *)malloc(context_num*sizeof(struct matrix));
        struct matrix * D=(struct matrix *)malloc(context_num*sizeof(struct matrix));
        struct matrix * S=(struct matrix *)malloc(context_num*sizeof(struct matrix));
        struct matrix * E=(struct matrix *)malloc(context_num*sizeof(struct matrix));
        struct matrix * H=(struct matrix *)malloc(context_num*sizeof(struct matrix));
        struct matrix * R=(struct matrix *)malloc(context_num*sizeof(struct matrix));
        
        //construct dik
        temp=0;
        for (k=0; k<context_num; k++) {
            if (J[i*context_num+k].num==0) {
                d[k].length=0;
                d[k].v=NULL;
            }
            else{
                d[k].length=J[i*context_num+k].num;
                d[k].v=(double *)malloc(d[k].length*sizeof(double));
                for (j=0; j<d[k].length; j++) {
                    d[k].v[j]=context[k].matrix[i*context[k].item_num+J[i*context_num+k].index[j]].value-alpha[temp+J[i*context_num+k].index[j]];
                }
            }
            temp=temp+context[k].item_num;
        }
        
        //construct Cik
        temp=0;
        for (k=0; k<context_num; k++) {
            if (J[i*context_num+k].num==0) {
                C[k].row_num=0;
                C[k].col_num=0;
                C[k].m=NULL;
            }
            else{
                C[k].row_num=J[i*context_num+k].num;
                C[k].col_num=q[k];
                C[k].m=(double *)malloc(C[k].row_num*C[k].col_num*sizeof(double));
                for (j=0; j<J[i*context_num+k].num; j++) {
                    for (t=0; t<q[k]; t++) {
                        C[k].m[j*q[k]+t]=beta[temp+J[i*context_num+k].index[j]].v[t];
                    }
                }
            }
            temp=temp+context[k].item_num;
        }
        
        //construct Dik
        for (k=0; k<context_num; k++) {
            if (J[i*context_num+k].num==0) {
                D[k].row_num=0;
                D[k].col_num=0;
                D[k].m=NULL;
            }
            else{
                D[k].row_num=J[i*context_num+k].num;
                D[k].col_num=J[i*context_num+k].num;
                D[k].m=(double *)malloc(D[k].row_num*D[k].col_num*sizeof(double));
                for (j=0; j<D[k].row_num; j++) {
                    for (t=0; t<D[k].col_num; t++) {
                        if (j==t) {
                            D[k].m[j*D[k].col_num+t]=sigma2_y[k]*w[k].m[i*w[k].col_num+J[i*context_num+k].index[j]];
                        }
                        else{
                            D[k].m[j*D[k].col_num+t]=0;
                        }
                    }
                }
            }
        }
        
        //construct Sik
        for (k=0; k<context_num; k++) {
            matrix1=matrixTransit(A[k].m, A[k].row_num, A[k].col_num);
            
            if (matrix1!=NULL) {
                matrix2=matrixMultiply(A[k].m, matrix1, A[k].row_num, A[k].col_num, A[k].col_num, A[k].row_num);
                free(matrix1);
            }
            else{
                printf("fail calculating S %d %d\n", i, k);
                S[k].m=NULL;
                S[k].row_num=0;
                S[k].col_num=0;
                continue;
            }
            
            if (matrix2!=NULL) {
                matrix1=constMatrixMultiply(sigma2_u[0], matrix2, A[k].row_num, A[k].row_num);
                free(matrix2);
            }
            else{
                printf("fail calculating S %d %d\n", i, k);
                S[k].m=NULL;
                S[k].row_num=0;
                S[k].col_num=0;
                continue;
            }
            
            if (matrix1!=NULL) {
                matrix2=constIdentityMatrix(sigma2_z[k], A[k].row_num);
            }
            else{
                printf("fail calculating S %d %d\n", i, k);
                S[k].m=NULL;
                S[k].row_num=0;
                S[k].col_num=0;
                continue;
            }
            
            S[k].m=matrixAddition(matrix1, matrix2, A[k].row_num, A[k].row_num, A[k].row_num, A[k].row_num);
            S[k].row_num=A[k].row_num;
            S[k].col_num=A[k].row_num;
            
            free(matrix1);
            free(matrix2);
        }
        
        //construct Eik
        for (k=0; k<context_num; k++) {
            matrix1=matrixTransit(C[k].m, C[k].row_num, C[k].col_num);
            
            if (matrix1!=NULL) {
                matrix2=matrixMultiply(S[k].m, matrix1, S[k].row_num, S[k].col_num, C[k].col_num, C[k].row_num);
                free(matrix1);
            }
            else{
                //printf("fail calculating E %d %d\n", i, k);
                E[k].m=NULL;
                E[k].row_num=0;
                E[k].col_num=0;
                continue;
            }
            
            if (matrix2!=NULL) {
                matrix1=matrixMultiply(C[k].m, matrix2, C[k].row_num, C[k].col_num, S[k].row_num, C[k].row_num);
                free(matrix2);
            }
            else{
                //printf("fail calculating E %d %d\n", i, k);
                E[k].m=NULL;
                E[k].row_num=0;
                E[k].col_num=0;
                continue;
            }
            
            if (matrix1!=NULL) {
                E[k].m=matrixAddition(matrix1, D[k].m, C[k].row_num, C[k].row_num, D[k].row_num, D[k].col_num);
                E[k].row_num=D[k].row_num;
                E[k].col_num=D[k].col_num;
                free(matrix1);
            }
            else{
                //printf("fail calculating E %d %d\n", i, k);
                E[k].m=NULL;
                E[k].row_num=0;
                E[k].col_num=0;
                continue;
            }
        
        }
        
        //construct Hik;
        for (k=0; k<context_num; k++) {
            matrix1=matrixTransit(A[k].m, A[k].row_num, A[k].col_num);
            
            if (matrix1!=NULL) {
                matrix2=constMatrixMultiply(sigma2_u[0], matrix1, A[k].col_num, A[k].row_num);
                free(matrix1);
            }
            else{
                printf("fail calculating H %d %d\n", i, k);
                H[k].m=NULL;
                H[k].row_num=0;
                H[k].col_num=0;
                continue;
            }
            
            if (matrix2!=NULL) {
                matrix1=matrixInverse(S[k].m, S[k].row_num, S[k].col_num);
            }
            else{
                printf("fail calculating H %d %d\n", i, k);
                H[k].m=NULL;
                H[k].row_num=0;
                H[k].col_num=0;
                continue;
            }
            
            if (matrix1!=NULL) {
                H[k].m=matrixMultiply(matrix2, matrix1, A[k].col_num, A[k].row_num, S[k].row_num, S[k].col_num);
                H[k].row_num=A[k].col_num;
                H[k].col_num=S[k].col_num;
                free(matrix2);
                free(matrix1);
            }
            else{
                printf("fail calculating H %d %d\n", i, k);
                H[k].m=NULL;
                H[k].row_num=0;
                H[k].col_num=0;
                continue;
            }
        }
        
        //construct Rik
        for (k=0; k<context_num; k++) {
            matrix1=matrixInverse(S[k].m, S[k].row_num, S[k].col_num);
            
            if (matrix1!=NULL) {
                matrix2=matrixMultiply(matrix1, A[k].m, S[k].row_num, S[k].col_num, A[k].row_num, A[k].col_num);
                free(matrix1);
                matrix1=matrixTransit(A[k].m, A[k].row_num, A[k].col_num);
            }
            else{
                printf("fail calculating R %d %d\n", i, k);
                R[k].m=NULL;
                R[k].row_num=0;
                R[k].col_num=0;
                continue;
            }
            
            if ((matrix1!=NULL)&&(matrix2!=NULL)) {
                matrix3=matrixMultiply(matrix1, matrix2, A[k].col_num, A[k].row_num, S[k].row_num, A[k].col_num);
                free(matrix1);
                free(matrix2);
            }
            else{
                printf("fail calculating R %d %d\n", i, k);
                R[k].m=NULL;
                R[k].row_num=0;
                R[k].col_num=0;
                continue;
            }
            
            if (matrix3!=NULL) {
                matrix1=constIdentityMatrix(sigma2_u[0], A[k].col_num);
                matrix2=constMatrixMultiply((-1)*sigma2_u[0]*sigma2_u[0], matrix3, A[k].col_num, A[k].col_num);
                free(matrix3);
            }
            else{
                printf("fail calculating R %d %d\n", i, k);
                R[k].m=NULL;
                R[k].row_num=0;
                R[k].col_num=0;
                continue;
            }
            
            if ((matrix1!=NULL)&&(matrix2!=NULL)) {
                R[k].m=matrixAddition(matrix1, matrix2, A[k].col_num, A[k].col_num, A[k].col_num, A[k].col_num);
                R[k].row_num=A[k].col_num;
                R[k].col_num=A[k].col_num;
                free(matrix1);
                free(matrix2);
            }
            else{
                printf("fail calculating R %d %d\n", i, k);
                R[k].m=NULL;
                R[k].row_num=0;
                R[k].col_num=0;
                continue;
            }
        }
        
        //calculate zik|ik
        struct vector *Ez_cd=(struct vector *)malloc(context_num*sizeof(struct vector));
        struct matrix *Vz_cd=(struct matrix *)malloc(context_num*sizeof(struct matrix));
        
        for (k=0; k<context_num; k++) {
            Ez_cd[k].length=q[k];
            Vz_cd[k].row_num=q[k];
            Vz_cd[k].col_num=q[k];
            
            if (J[i*context_num+k].num==0) {
                Ez_cd[k].v=(double *)malloc(q[k]*sizeof(double));
                for (t=0; t<Ez_cd[k].length; t++) {
                    Ez_cd[k].v[t]=0;
                }
                
                Vz_cd[k].m=(double *)malloc(q[k]*q[k]*sizeof(double));
                for (j=0; j<q[k]; j++) {
                    for (t=0; t<q[k]; t++) {
                        Vz_cd[k].m[j*q[k]+t]=S[k].m[j*q[k]+t];
                    }
                }
            }
            else{
                matrix1=matrixTransit(C[k].m, C[k].row_num, C[k].col_num);
                matrix2=matrixInverse(E[k].m, E[k].row_num, E[k].col_num);
                if ((matrix1!=NULL)&&(matrix2!=NULL)) {
                    matrix3=matrixMultiply(matrix1, matrix2, C[k].col_num, C[k].row_num, E[k].row_num, E[k].col_num);
                    free(matrix1);
                    free(matrix2);
                }
                else{
                    printf("fail calculating Ez_cd %d %d\n", i, k);
                    continue;
                }
                
                if (matrix3!=NULL) {
                    matrix2=matrixMultiply(S[k].m, matrix3, S[k].row_num, S[k].col_num, C[k].col_num, E[k].col_num);
                    free(matrix3);
                }
                else{
                    printf("fail calculating Ez_cd %d %d\n", i, k);
                    continue;
                }
                
                if (matrix2!=NULL) {
                    Ez_cd[k].v=matrixMultiply(matrix2, d[k].v, S[k].row_num, E[k].col_num, d[k].length, 1);
                    matrix1=matrixMultiply(matrix2, C[k].m, S[k].row_num, E[k].col_num, C[k].row_num, C[k].col_num);
                    free(matrix2);
                }
                else{
                    printf("fail calculating Ez_cd %d %d\n", i, k);
                    continue;
                }
                
                if (matrix1!=NULL) {
                    matrix2=matrixMultiply(matrix1, S[k].m, S[k].row_num, C[k].col_num, S[k].row_num, S[k].col_num);
                    free(matrix1);
                }
                else{
                    printf("fail calculating Vz_cd %d %d\n", i, k);
                    continue;
                }
                
                if (matrix2!=NULL) {
                    matrix1=constMatrixMultiply(-1, matrix2, S[k].row_num, S[k].col_num);
                    free(matrix2);
                }
                else{
                    printf("fail calculating Vz_cd %d %d\n", i, k);
                    continue;
                }
                
                if (matrix1!=NULL) {
                    Vz_cd[k].m=matrixAddition(S[k].m, matrix1, S[k].row_num, S[k].col_num, S[k].row_num, S[k].col_num);
                    free(matrix1);
                }
                else{
                    printf("fail calculating Vz_cd %d %d\n", i, k);
                    continue;
                }
            }

        }
            
        //calculate Vu_cd ik
        struct matrix * Vu_cd=(struct matrix *)malloc(context_num*sizeof(struct matrix));
        for (k=0; k<context_num; k++) {
            matrix1=matrixTransit(H[k].m, H[k].row_num, H[k].col_num);
            if (matrix1!=NULL) {
                matrix2=matrixMultiply(Vz_cd[k].m, matrix1, Vz_cd[k].row_num, Vz_cd[k].col_num, H[k].col_num, H[k].row_num);
                free(matrix1);
            }
            else{
                printf("fail calculating Vu_cd %d %d\n", i, k);
                continue;
            }
                
            if (matrix2!=NULL) {
                matrix1=matrixMultiply(H[k].m, matrix2, H[k].row_num, H[k].col_num, Vz_cd[k].row_num, H[k].row_num);
                free(matrix2);
            }
            else{
                printf("fail calculating Vu_cd %d %d\n", i, k);
            continue;
            }
                
            if (matrix1!=NULL) {
                Vu_cd[k].m=matrixAddition(matrix1, R[k].m, H[k].row_num, H[k].row_num, R[k].row_num, R[k].col_num);
                Vu_cd[k].row_num=H[k].row_num;
                Vu_cd[k].col_num=H[k].row_num;
                free(matrix1);
            }
            else{
                printf("fail calculating Vu_cd %d %d\n", i, k);
                continue;
            }
        }
        
        //calculate Vui
        if (Vu[i].m!=NULL) {
            free(Vu[i].m);
        }
        matrix1=(double *)malloc(r[0]*r[0]*sizeof(double));
        for (j=0; j<r[0]; j++) {
            for (t=0; t<r[0]; t++) {
                matrix1[j*r[0]+t]=0;
            }
        }
        for (k=0; k<context_num; k++) {
            matrix2=matrixInverse(Vu_cd[k].m, Vu_cd[k].row_num, Vu_cd[k].col_num);
            matrix3=constIdentityMatrix(-1/sigma2_u[0], r[0]);
            if ((matrix2!=NULL)&&(matrix3!=NULL)) {
                matrix4=matrixAddition(matrix2, matrix3, Vu_cd[k].row_num, Vu_cd[k].col_num, r[0], r[0]);
                free(matrix2);
                free(matrix3);
            }
            else{
                printf("fail calculating Vu %d\n", i);
                continue;
            }
            if (matrix4!=NULL) {
                matrix2=matrixAddition(matrix1, matrix4, r[0], r[0], r[0], r[0]);
                free(matrix1);
                free(matrix4);
            }
            else{
                printf("fail calculating Vu %d\n", i);
                continue;
            }
            if (matrix2!=NULL) {
                matrix1=matrix2;
                matrix2=NULL;
            }
            else{
                printf("fail calculating Vu %d\n", i);
                continue;
            }
        }
        matrix2=constIdentityMatrix(1/sigma2_u[0], r[0]);
        if ((matrix1!=NULL)&&(matrix2!=NULL)) {
            matrix3=matrixAddition(matrix2, matrix1, r[0], r[0], r[0], r[0]);
            free(matrix1);
            free(matrix2);
        }
        else{
            printf("fail calculating Vu %d\n", i);
            continue;
        }
        if (matrix3!=NULL) {
            Vu[i].m=matrixInverse(matrix3, r[0], r[0]);
            Vu[i].row_num=r[0];
            Vu[i].col_num=r[0];
            free(matrix3);
        }
        else{
            printf("fail calculating Vu %d\n", i);
            continue;
        }
        
        //calculate Eu
        if (Eu[i].v!=NULL) {
            free(Eu[i].v);
        }
        matrix1=(double *)malloc(r[0]*sizeof(double));
        for (j=0; j<r[0]; j++) {
            matrix1[j]=0;
        }
        for (k=0; k<context_num; k++) {
            matrix2=matrixInverse(Vu_cd[k].m, Vu_cd[k].row_num, Vu_cd[k].col_num);
            if (matrix2!=NULL) {
                matrix3=matrixMultiply(matrix2, H[k].m, Vu_cd[k].row_num, Vu_cd[k].col_num, H[k].row_num, H[k].col_num);
                free(matrix2);
            }
            else{
                printf("fail calculating Eu %d\n", i);
                continue;
            }
            if (matrix3!=NULL) {
                matrix2=matrixMultiply(matrix3, Ez_cd[k].v, Vu_cd[k].row_num, H[k].col_num, Ez_cd[k].length, 1);
                free(matrix3);
            }
            else{
                printf("fail calculating Eu %d\n", i);
                continue;
            }
            if (matrix2!=NULL) {
                matrix3=matrixAddition(matrix1, matrix2, r[0], 1, Vu_cd[k].row_num, 1);
                free(matrix1);
                free(matrix2);
            }
            else{
                printf("fail calculating Eu %d\n", i);
                continue;
            }
            if (matrix3!=NULL) {
                matrix1=matrix3;
                matrix3=NULL;
            }
            else{
                printf("fail calculating Eu %d\n", i);
                continue;
            }
        }
        Eu[i].v=matrixMultiply(Vu[i].m, matrix1, Vu[i].row_num, Vu[i].col_num, r[0], 1);
        Eu[i].length=Vu[i].row_num;
        free(matrix1);
        
        //calculate Ez ik
        for (k=0; k<context_num; k++) {
            if (Ez[k*context[0].user_num+i].v!=NULL) {
                free(Ez[k*context[0].user_num+i].v);
            }
            matrix1=matrixMultiply(H[k].m, Ez_cd[k].v, H[k].row_num, H[k].col_num, Ez_cd[k].length, 1);
            if (matrix1!=NULL) {
                matrix2=constMatrixMultiply(-1, matrix1, H[k].row_num, 1);
                free(matrix1);
            }
            else{
                printf("fail calculating Ez %d %d\n", i, k);
                continue;
            }
            if (matrix2!=NULL) {
                matrix1=matrixAddition(Eu[i].v, matrix2, Eu[i].length, 1, H[k].row_num, 1);
                matrix3=matrixInverse(Vu_cd[k].m, Vu_cd[k].row_num, Vu_cd[k].col_num);
                free(matrix2);
            }
            else{
                printf("fail calculating Ez %d %d\n", i, k);
                continue;
            }
            if ((matrix1!=NULL)&&(matrix3!=NULL)) {
                matrix2=matrixMultiply(matrix3, matrix1, Vu_cd[k].row_num, Vu_cd[k].col_num, H[k].row_num, 1);
                free(matrix1);
                free(matrix3);
                matrix3=matrixTransit(H[k].m, H[k].row_num, H[k].col_num);
            }
            else{
                printf("fail calculating Ez %d %d\n", i, k);
                continue;
            }
            if ((matrix2!=NULL)&&(matrix3!=NULL)) {
                matrix1=matrixMultiply(matrix3, matrix2, H[k].col_num, H[k].row_num, Vu_cd[k].row_num, 1);
                free(matrix2);
                free(matrix3);
            }
            else{
                printf("fail calculating Ez %d %d\n", i, k);
                continue;
            }
            if (matrix1!=NULL) {
                matrix2=matrixMultiply(Vz_cd[k].m, matrix1, Vz_cd[k].row_num, Vz_cd[k].col_num, H[k].col_num, 1);
                free(matrix1);
            }
            else{
                printf("fail calculating Ez %d %d\n", i, k);
                continue;
            }
            if (matrix2!=NULL) {
                Ez[k*context[0].user_num+i].v=matrixAddition(Ez_cd[k].v, matrix2, Ez_cd[k].length, 1, Vz_cd[k].row_num, 1);
                Ez[k*context[0].user_num+i].length=Ez_cd[k].length;
                free(matrix2);
            }
            else{
                printf("fail calculating Ez %d %d\n", i, k);
                continue;
            }
        }
        
        //calculate Vz ik
        for (k=0; k<context_num; k++) {
            if (Vz[k*context[0].user_num+i].m!=NULL) {
                free(Vz[k*context[0].user_num+i].m);
            }
            matrix1=constMatrixMultiply(-1, Vu_cd[k].m, Vu_cd[k].row_num, Vu_cd[k].col_num);
            matrix2=matrixInverse(Vu_cd[k].m, Vu_cd[k].row_num, Vu_cd[k].col_num);
            if ((matrix1!=NULL)&&(matrix2!=NULL)) {
                matrix3=matrixAddition(Vu[i].m, matrix1, Vu[i].row_num, Vu[i].col_num, Vu_cd[k].row_num, Vu_cd[k].col_num);
                free(matrix1);
                matrix1=matrixMultiply(H[k].m, Vz_cd[k].m, H[k].row_num, H[k].col_num, Vz_cd[k].row_num, Vz_cd[k].col_num);
            }
            else{
                printf("fail calculating Vz %d %d\n", i, k);
                continue;
            }
            if ((matrix3!=NULL)&&(matrix1!=NULL)) {
                matrix4=matrixMultiply(matrix2, matrix3, Vu_cd[k].row_num, Vu_cd[k].col_num, Vu[i].row_num, Vu[i].col_num);
                free(matrix3);
            }
            else{
                printf("fail calculating Vz %d %d\n", i, k);
                continue;
            }
            if (matrix4!=NULL) {
                matrix3=matrixMultiply(matrix4, matrix2, Vu_cd[k].row_num, Vu[i].col_num, Vu_cd[k].row_num, Vu_cd[k].col_num);
                free(matrix2);
                free(matrix4);
                matrix2=matrixTransit(matrix1, H[k].row_num, Vz_cd[k].col_num);
            }
            else{
                printf("fail calculating Vz %d %d\n", i, k);
                continue;
            }
            if ((matrix3!=NULL)&&(matrix2!=NULL)) {
                matrix4=matrixMultiply(matrix2, matrix3, Vz_cd[k].col_num, H[k].row_num, Vu_cd[k].row_num, Vu_cd[k].col_num);
                free(matrix2);
                free(matrix3);
            }
            else{
                printf("fail calculating Vz %d %d\n", i, k);
                continue;
            }
            if (matrix4!=NULL) {
                matrix2=matrixMultiply(matrix4, matrix1, Vz_cd[k].col_num, Vu_cd[k].col_num, H[k].row_num, Vz_cd[k].col_num);
                free(matrix4);
                free(matrix1);
            }
            else{
                printf("fail calculating Vz %d %d\n", i, k);
                continue;
            }
            if (matrix2!=NULL) {
                Vz[k*context[0].user_num+i].m=matrixAddition(Vz_cd[k].m, matrix2, Vz_cd[k].row_num, Vz_cd[k].col_num, Vz_cd[k].col_num, Vz_cd[k].col_num);
                free(matrix2);
                Vz[k*context[0].user_num+i].row_num=Vz_cd[k].row_num;
                Vz[k*context[0].user_num+i].col_num=Vz_cd[k].col_num;
            }
            else{
                printf("fail calculating Vz %d %d\n", i, k);
                continue;
            }
        }
        
        //calculate Covzu
        for (k=0; k<context_num; k++) {
            if (Covzu[k*context[0].user_num+i].m!=NULL) {
                free(Covzu[k*context[0].user_num+i].m);
            }
            matrix1=matrixTransit(H[k].m, H[k].row_num, H[k].col_num);
            matrix2=matrixInverse(Vu_cd[k].m, Vu_cd[k].row_num, Vu_cd[k].col_num);
            if ((matrix1!=NULL)&&(matrix2!=NULL)) {
                matrix3=matrixMultiply(Vz_cd[k].m, matrix1, Vz_cd[k].row_num, Vz_cd[k].col_num, H[k].col_num, H[k].row_num);
                free(matrix1);
                matrix4=matrixMultiply(matrix2, Vu[i].m, Vu_cd[k].row_num, Vu_cd[k].col_num, Vu[i].row_num, Vu[i].col_num);
                free(matrix2);
            }
            else{
                printf("fail calculating Covzu %d %d\n", i, k);
                continue;
            }
            if ((matrix3!=NULL)&&(matrix4!=NULL)) {
                Covzu[k*context[0].user_num+i].m=matrixMultiply(matrix3, matrix4, Vz_cd[k].row_num, H[k].row_num, Vu_cd[k].row_num, Vu[i].col_num);
                free(matrix3);
                free(matrix4);
                Covzu[k*context[0].user_num+i].row_num=Vz_cd[k].row_num;
                Covzu[k*context[0].user_num+i].col_num=Vu[i].col_num;
            }
            else{
                printf("fail calculating Covzu %d %d\n", i, k);
                continue;
            }
        }
        
        //set space free
        for (k=0; k<context_num; k++) {
            if (d[k].v!=NULL) {
                free(d[k].v);
            }
        }
        free(d);
        for (k=0; k<context_num; k++) {
            if (C[k].m!=NULL) {
                free(C[k].m);
            }
        }
        free(C);
        for (k=0; k<context_num; k++) {
            if (D[k].m!=NULL) {
                free(D[k].m);
            }
        }
        free(D);
        for (k=0; k<context_num; k++) {
            if (S[k].m!=NULL) {
                free(S[k].m);
            }
        }
        free(S);
        for (k=0; k<context_num; k++) {
            if (E[k].m!=NULL) {
                free(E[k].m);
            }
        }
        free(E);
        for (k=0; k<context_num; k++) {
            if (H[k].m!=NULL) {
                free(H[k].m);
            }
        }
        free(H);
        for (k=0; k<context_num; k++) {
            if (R[k].m!=NULL) {
                free(R[k].m);
            }
        }
        free(R);
        for (k=0; k<context_num; k++) {
            if (Ez_cd[k].v!=NULL) {
                free(Ez_cd[k].v);
            }
        }
        free(Ez_cd);
        for (k=0; k<context_num; k++) {
            if (Vz_cd[k].m!=NULL) {
                free(Vz_cd[k].m);
            }
        }
        free(Vz_cd);
        for (k=0; k<context_num; k++) {
            if (Vu_cd[k].m!=NULL) {
                free(Vu_cd[k].m);
            }
        }
        free(Vu_cd);
    }
    
}

void MStep(int context_num, struct resMatrix *context, double *alpha, double *sigma2_z, double *sigma2_y, double *sigma2_u, struct vector *beta, struct vector *Eu, struct vector *Ez, struct matrix *w, struct matrix *A, struct matrix *Vu, struct matrix *Vz, struct matrix *Covzu, int *q, int *r, struct indexCell *J){
    
    int i,j,k,t,l,temp,sum;
    double *matrix1, *matrix2, *matrix3, *matrix4, *matrix5;
    struct indexCell *I;
    
    //get I={i: k belongs to Ki}
    I=(struct indexCell *)malloc(context_num*sizeof(struct indexCell));
    for (k=0; k<context_num; k++) {
        I[k].num=0;
        for (i=0; i<context[k].user_num; i++) {
            if (J[i*context_num+k].num>0) {
                I[k].num++;
            }
        }
        if (I[k].num==0) {
            I[k].index=NULL;
        }
        else{
            I[k].index=(int *)malloc(I[k].num*sizeof(int));
            t=0;
            for (i=0; i<context[k].user_num; i++) {
                if (J[i*context_num+k].num>0) {
                    I[k].index[t]=i;
                    t++;
                }
            }
        }
    }
    
    //estimate Ak,l and sigma2_z[k]
    struct vector * A_rows=(struct vector *)malloc(A[0].row_num*sizeof(struct vector));
    double *loss=(double *)malloc(A[k].row_num*sizeof(double));
    
    for (k=0; k<context_num; k++) {
       
        //calculate Ak,l and loss l
        for (l=0; l<A[0].row_num; l++) {
            
            //calculate Ak,l
            matrix1=(double *)malloc(r[0]*r[0]*sizeof(double));
            for (i=0; i<r[0]; i++) {
                for (t=0; t<r[0]; t++) {
                    matrix1[i*r[0]+t]=0;
                }
            }
            for (t=0; t<I[k].num; t++) {
                matrix2=matrixTransit(Eu[I[k].index[t]].v, Eu[I[k].index[t]].length, 1);
                if (matrix2!=NULL) {
                    matrix3=matrixMultiply(Eu[I[k].index[t]].v, matrix2, Eu[I[k].index[t]].length, 1, 1, Eu[I[k].index[t]].length);
                    free(matrix2);
                }
                else{
                    printf("fail calculating A %d %d\n", k, l);
                    continue;
                }
                if (matrix3!=NULL) {
                    matrix2=matrixAddition(Vu[I[k].index[t]].m, matrix3, Vu[I[k].index[t]].row_num, Vu[I[k].index[t]].col_num, Eu[I[k].index[t]].length, Eu[I[k].index[t]].length);
                    free(matrix3);
                }
                else{
                    printf("fail calculating A %d %d\n", k, l);
                    continue;
                }
                if (matrix2!=NULL) {
                    matrix3=matrixAddition(matrix1, matrix2, r[0], r[0], Vu[I[k].index[t]].row_num, Eu[I[k].index[t]].length);
                    free(matrix2);
                    free(matrix1);
                }
                else{
                    printf("fail calculating A %d %d\n", k, l);
                    continue;
                }
                if (matrix3!=NULL) {
                    matrix1=matrix3;
                    matrix3=NULL;
                }
                else{
                    printf("fail calculating A %d %d\n", k, l);
                    continue;
                }
            }
            matrix2=(double *)malloc(r[0]*sizeof(double));
            for (t=0; t<r[0]; t++) {
                matrix2[t]=0;
            }
            for (t=0; t<I[k].num; t++) {
                matrix3=matrixMultiply(Eu[I[k].index[t]].v, Ez[k*context[0].user_num+I[k].index[t]].v+l, Eu[I[k].index[t]].length, 1, 1, 1);
                matrix4=(double *)malloc(r[0]*sizeof(double));
                for (i=0; i<r[0]; i++) {
                    matrix4[i]=Covzu[k*context[0].user_num+I[k].index[t]].m[l*r[0]+i];
                }
                if (matrix3!=NULL) {
                    matrix5=matrixAddition(matrix3, matrix4, r[0], 1, r[0], 1);
                    free(matrix3);
                    free(matrix4);
                }
                else{
                    printf("fail calculating A %d %d\n", k, l);
                    continue;
                }
                if (matrix5!=NULL) {
                    matrix3=matrixAddition(matrix2, matrix5, r[0], 1, r[0], 1);
                    free(matrix2);
                    free(matrix5);
                }
                else{
                    printf("fail calculating A %d %d\n", k, l);
                    continue;
                }
                if (matrix3!=NULL) {
                    matrix2=matrix3;
                    matrix3=NULL;
                }
                else{
                    printf("fail calculating A %d %d\n", k, l);
                    continue;
                }
            }
            if (matrix1!=NULL) {
                matrix3=matrixInverse(matrix1, r[0], r[0]);
                //free(matrix1);    
            }
            else{
                printf("fail calculating A %d %d\n", k, l);
                continue;
            }
            if ((matrix3!=NULL)&&(matrix2!=NULL)) {
                A_rows[l].length=r[0];
                A_rows[l].v=matrixMultiply(matrix3, matrix2, r[0], r[0], r[0], 1);
                free(matrix3);
                //free(matrix2); 
            }
            else{
                printf("fail calculating A %d %d\n", k, l);
                continue;
            }
            
            //calculate loss l
            loss[l]=0;
            temp=0;
            for (t=0; t<I[k].num; t++) {
                temp=temp+Ez[k*context[0].user_num+I[k].index[t]].v[l]*Ez[k*context[0].user_num+I[k].index[t]].v[l]+Vz[k*context[0].user_num+I[k].index[t]].m[l*Vz[k*context[0].user_num+I[k].index[t]].col_num+l];
            }
            loss[l]=loss[l]+temp;
            
            matrix3=matrixMultiply(A_rows[l].v, matrix1, 1, r[0], r[0], r[0]);
            free(matrix1);
            matrix4=matrixMultiply(matrix2, A_rows[l].v, 1, r[0], r[0], 1);
            free(matrix2);
            if ((matrix3!=NULL)&&(matrix4!=NULL)) {
                loss[l]=loss[l]+(-2)*matrix4[0];
                free(matrix4);
                matrix4=matrixMultiply(matrix3, A_rows[l].v, 1, r[0], r[0], 1);
            }
            else{
                printf("fail calculating loss %d %d\n", k, l);
                continue;
            }
            if (matrix4!=NULL) {
                loss[l]=loss[l]+matrix4[0];
                free(matrix4);
            }
            else{
                printf("fail calculating loss %d %d\n", k, l);
                continue;
            }
        }
        
        //construct Ak
        for (l=0; l<A[k].row_num; l++) {
            for (t=0; t<A[k].col_num; t++) {
                A[k].m[l*A[k].col_num+t]=A_rows[l].v[t];
            }
        }
        
        //calculate sigma2_z
        temp=0;
        for (l=0; l<A[k].row_num; l++) {
            temp=temp+loss[l];
        }
        sigma2_z[k]=temp/(q[k]*I[k].num);
        
        for (l=0; l<A[k].row_num; l++) {
            if (A_rows[l].v!=NULL) {
                free(A_rows[l].v);
            }
        }
    }
    free(A_rows);
    free(loss);
    
    //get I={i: k belongs to Ki & j belongs to Jik}
    temp=0;
    for (k=0; k<context_num; k++) {
        temp=temp+context[k].item_num;
        if(I[k].index!=NULL){
            free(I[k].index);
        }
    }
    free(I);
    I=(struct indexCell *)malloc(temp*sizeof(struct indexCell));
    temp=0;
    for (k=0; k<context_num; k++) {
        for (j=0; j<context[k].item_num; j++) {
            I[temp+j].num=0;
            for (i=0; i<context[k].user_num; i++) {
                if (context[k].matrix[i*context[k].item_num+j].indicator==1) {
                    I[temp+j].num++;
                }
            }
            if (I[temp+j].num==0) {
                I[temp+j].index=NULL;
            }
            else{
                I[temp+j].index=(int *)malloc(I[temp+j].num*sizeof(int));
                t=0;
                for (i=0; i<context[k].user_num; i++) {
                    if (context[k].matrix[i*context[k].item_num+j].indicator==1) {
                        I[temp+j].index[t]=i;
                        t++;
                    }
                }
            }
        }
        temp=temp+context[k].item_num;
    }
    
    //estimate alpha jk   and  beta jk  and sigma2_y
    struct vector albeta;
    albeta.length=1+q[0];
    albeta.v=NULL;
    
    temp=0;
    for (k=0; k<context_num; k++) {
        loss=(double *)malloc(context[k].item_num*sizeof(double));
        for (j=0; j<context[k].item_num; j++) {
            //estimate alpha jk and beta jk
            if (I[temp+j].num==0) {
                loss[j]=0;
                continue;
            }
            matrix1=(double *)malloc((1+q[k])*(1+q[k])*sizeof(double));
            for (i=0; i<(1+q[k]); i++) {
                for (t=0; t<(1+q[k]); t++) {
                    matrix1[i*(1+q[k])+t]=0;
                }
            }
            matrix2=(double *)malloc((1+q[k])*sizeof(double));
            for (i=0; i<(1+q[k]); i++) {
                matrix2[i]=0;
            }
            for (t=0; t<I[temp+j].num; t++) {
                matrix3=V_one(Ez[k*context[0].user_num+I[temp+j].index[t]].v, Ez[k*context[0].user_num+I[temp+j].index[t]].length);
                if (matrix3!=NULL) {
                    matrix4=constMatrixMultiply(context[k].matrix[I[temp+j].index[t]*context[k].item_num+j].value/w[k].m[I[temp+j].index[t]*w[k].col_num+j], matrix3, 1+q[k], 1);
                    matrix5=matrixMultiply(matrix3, matrix3, 1+q[k], 1, 1, 1+q[k]);
                    free(matrix3);
                }
                else{
                    printf("fail calculating alpha beta %d %d\n", j, k);
                    continue;
                }
                if (matrix4!=NULL) {
                    matrix3=matrixAddition(matrix2, matrix4, 1+q[k], 1, 1+q[k], 1);
                    free(matrix2);
                    free(matrix4);
                }
                else{
                    printf("fail calculating alpha beta %d %d\n", j, k);
                    continue;
                }
                if (matrix3!=NULL) {
                    matrix2=matrix3;
                    matrix3=NULL;
                }
                else{
                    printf("fail calculating alpha beta %d %d\n", j, k);
                    continue;
                }
                matrix3=I_zero(Vz[k*context[0].user_num+I[temp+j].index[t]].m, Vz[k*context[0].user_num+I[temp+j].index[t]].row_num, Vz[k*context[0].user_num+I[temp+j].index[t]].col_num);
                if ((matrix3!=NULL)&&(matrix5!=NULL)) {
                    matrix4=matrixAddition(matrix3, matrix5, 1+q[k], 1+q[k], 1+q[k], 1+q[k]);
                    free(matrix3);
                    free(matrix5);
                }
                else{
                    printf("fail calculating alpha beta %d %d\n", j, k);
                    continue;
                }
                if (matrix4!=NULL) {
                    matrix3=constMatrixMultiply(1/w[k].m[I[temp+j].index[t]*w[k].col_num+j], matrix4, 1+q[k], 1+q[k]);
                    free(matrix4);
                }
                else{
                    printf("fail calculating alpha beta %d %d\n", j, k);
                    continue;
                }
                if (matrix3!=NULL) {
                    matrix4=matrixAddition(matrix1, matrix3, 1+q[k], 1+q[k], 1+q[k], 1+q[k]);
                    free(matrix1);
                    free(matrix3);
                }
                else{
                    printf("fail calculating alpha beta %d %d\n", j, k);
                    continue;
                }
                if (matrix4!=NULL) {
                    matrix1=matrix4;
                    matrix4=NULL;
                }
                else{
                    printf("fail calculating alpha beta %d %d\n", j, k);
                    continue;
                }
            }
            matrix3=matrixInverse(matrix1, 1+q[k], 1+q[k]);
            if (matrix3!=NULL) {
                albeta.v=matrixMultiply(matrix3, matrix2, 1+q[k], 1+q[k], 1+q[k], 1);
                free(matrix3);
            }
            else{
                printf("fail calculating alpha beta %d %d\n", j, k);
                continue;
            }
            
            alpha[temp+j]=albeta.v[0];
            for (t=0; t<q[k]; t++) {
                beta[temp+j].v[t]=albeta.v[t+1];
            }
            
            //calculate loss jk
            loss[j]=0;
            for (t=0; t<I[temp+j].num; t++) {
                loss[j]=loss[j]+context[k].matrix[I[temp+j].index[t]*context[k].item_num+j].value*context[k].matrix[I[temp+j].index[t]*context[k].item_num+j].value/w[k].m[I[temp+j].index[t]*w[k].col_num+j];
            }
            if ((matrix1!=NULL)&&(matrix2!=NULL)) {
                matrix3=matrixMultiply(matrix1, albeta.v, 1+q[k], 1+q[k], 1+q[k], 1);
                matrix4=matrixMultiply(matrix2, albeta.v, 1, 1+q[k], 1+q[k], 1);
                free(matrix1);
                free(matrix2);
            }
            else{
                printf("fail calculating loss %d %d\n", j, k);
                continue;
            }
            if ((matrix3!=NULL)&&(matrix4!=NULL)) {
                matrix1=matrixMultiply(albeta.v, matrix3, 1, 1+q[k], 1+q[k], 1);
                free(matrix3);
                loss[j]=loss[j]+(-2)*matrix4[0];
                free(matrix4);
            }
            else{
                printf("fail calculating alpha beta %d %d\n", j, k);
                continue;
            }
            if (matrix1!=NULL) {
                loss[j]=loss[j]+matrix1[0];
                free(matrix1);
            }
            else{
                printf("fail calculating alpha beta %d %d\n", j, k);
                continue;
            }
        }
        
        //estimate sigma2_y
        sum=0;
        for (j=0; j<context[k].item_num; j++) {
            sum=sum+loss[j];
        }
        sigma2_y[k]=sum;
        sum=0;
        for (j=0; j<context[k].item_num; j++) {
            sum=sum+I[temp+j].num;
        }
        sigma2_y[k]=sigma2_y[k]/sum;
        
        free(loss);
        temp=temp+context[k].item_num;
    }
    
    
    temp=0;
    for (k=0; k<context_num; k++) {
        for (j=0; j<context[k].item_num; j++) {
            if (I[temp+j].index!=NULL) {
                free(I[temp+j].index);
            }
        }
        temp=temp+context[k].item_num;
    }
    free(I);
    
}

double RMSE_training(int context_num, struct resMatrix *context, double *alpha, struct vector *beta, struct vector *Ez, int *q, int *r){
    int i,j,k,num,temp;
    double RMSE, predicted, sum, *matrix;
    
    num=0;
    sum=0;
    temp=0;
    for (k=0; k<context_num; k++) {
        for (i=0; i<context[k].user_num; i++) {
            for (j=0; j<context[k].item_num; j++) {
                if (context[k].matrix[i*context[k].item_num+j].indicator==1) {
                    matrix=matrixMultiply(beta[temp+j].v, Ez[k*context[0].user_num+i].v, 1, beta[temp+j].length, Ez[k*context[0].user_num+i].length, 1);
                    if (matrix!=NULL) {
                        if ((matrix[0]+alpha[temp+j])<context[k].min_value) {
                            predicted=context[k].min_value;
                        }
                        else if((matrix[0]+alpha[temp+j])>context[k].max_value){
                            predicted=context[k].max_value;
                        }
                        else
                        predicted=matrix[0]+alpha[temp+j];
                        free(matrix);
                        sum=sum+(predicted-context[k].matrix[i*context[k].item_num+j].value)*(predicted-context[k].matrix[i*context[k].item_num+j].value);
                        num++;
                    }
                    else{
                        printf("fail calculating y %d %d %d\n", i, j, k);
                        continue;
                    }
                }
            }
        }
        temp=temp+context[k].item_num;
    }
    RMSE=sqrt(sum/num);
    
    return RMSE;
}

double RMSE_predict(int context_num, struct resMatrix *context, double *alpha, struct vector *beta, struct vector *Ez, int *q, int *r){
    int i,j,k,num,temp;
    double RMSE, predicted, sum, *matrix;
    
    num=0;
    sum=0;
    temp=0;
    for (k=0; k<context_num; k++) {
        for (i=0; i<context[k].user_num; i++) {
            for (j=0; j<context[k].item_num; j++) {
                if (context[k].matrix[i*context[k].item_num+j].indicator==0) {
                    matrix=matrixMultiply(beta[temp+j].v, Ez[k*context[0].user_num+i].v, 1, beta[temp+j].length, Ez[k*context[0].user_num+i].length, 1);
                    if (matrix!=NULL) {
                        if ((matrix[0]+alpha[temp+j])<context[k].min_value) {
                            predicted=context[k].min_value;
                        }
                        else if((matrix[0]+alpha[temp+j])>context[k].max_value){
                            predicted=context[k].max_value;
                        }
                        else
                            predicted=matrix[0]+alpha[temp+j];
                        free(matrix);
                        sum=sum+(predicted-context[k].matrix[i*context[k].item_num+j].value)*(predicted-context[k].matrix[i*context[k].item_num+j].value);
                        num++;
                    }
                    else{
                        printf("fail calculating y %d %d %d\n", i, j, k);
                        continue;
                    }
                }
            }
        }
        temp=temp+context[k].item_num;
    }
    RMSE=sqrt(sum/num);
    
    return RMSE;
}

void fitting(int context_num, struct resMatrix *context, double *alpha, double *sigma2_z, double *sigma2_y, double *sigma2_u, struct vector *beta, struct vector *Eu, struct vector *Ez, struct matrix *w, struct matrix *A, struct matrix *Vu, struct matrix *Vz, struct matrix *Covzu, int *q, int *r){
    
    int i,j,k,temp;
    
    //initialize parameters
    paraInitialize(context_num, context, alpha, sigma2_z, sigma2_y, sigma2_u, beta, w, A, q, r);
    
    struct indexCell *J=(struct indexCell *)malloc(context_num*context[0].user_num*sizeof(struct indexCell));
    
    //construct Jik
    for (i=0; i<context[0].user_num; i++) {
        for (k=0; k<context_num; k++) {
            temp=0;
            for (j=0; j<context[k].item_num; j++) {
                if (context[k].matrix[i*context[k].item_num+j].indicator==1) {
                    temp++;
                }
            }
            J[i*context_num+k].num=temp; 
            if(temp>0)  J[i*context_num+k].index=(int *)malloc(temp*sizeof(int));
            else    J[i*context_num+k].index=NULL;
            temp=0;
            for (j=0; j<context[k].item_num; j++) {
                if (context[k].matrix[i*context[k].item_num+j].indicator==1) {
                    J[i*context_num+k].index[temp]=j;
                    temp++;
                }
            }
        }
    }
    
    FILE *fp1=fopen("/Users/DaisyYang/Desktop/MyProjectC/MyProjectC/MyProjectC/EM_results/RSME.txt", "w");
    FILE *fp2=fopen("/Users/DaisyYang/Desktop/MyProjectC/MyProjectC/MyProjectC/EM_results/time.txt", "w");
    fprintf(fp1, "RMSE_training RSME_predict\n");
    double previous=0, present=0, predict=0;
    i=0;
    int indicator=0;
    
    time_t start_time;
    time(&start_time);
    
    while (indicator<=30) {
        i++;
        
        EStep(context_num, context, alpha, sigma2_z, sigma2_y, sigma2_u, beta, Eu, Ez, w, A, Vu, Vz, Covzu, q, r, J);
        
        MStep(context_num, context, alpha, sigma2_z, sigma2_y, sigma2_u, beta, Eu, Ez, w, A, Vu, Vz, Covzu, q, r, J);
        
        if(i>0){
            previous=present;
        }
        
        present=RMSE_training(context_num, context, alpha, beta, Ez, q, r);
        fprintf(fp1, "%f ",present);
        printf("%d: %f ",i,present);
        
        predict=RMSE_predict(context_num, context, alpha, beta, Ez, q, r);
        fprintf(fp1, "%f\n",predict);
        printf("%f\n",predict);
        
        if (((present-previous)<=0.000001)&&((present-previous)>=-0.000001)) {
            indicator++;
        }
        else{
            indicator=0;
        }
    }
    fclose(fp1);
    
    time_t end_time;
    time(&end_time);
    
    long used_time=end_time-start_time;
    long average_time=used_time/i;
    
    printf("start time: %ld\n",start_time);
    printf("end time: %ld\n",end_time);
    printf("used time: %ld\n",used_time);
    printf("iteration: %d\n",i);
    printf("average time: %ld\n",average_time);
    
    fprintf(fp2, "start time: %ld\n",start_time);
    fprintf(fp2, "end time: %ld\n",end_time);
    fprintf(fp2, "used time: %ld\n",used_time);
    fprintf(fp2, "iteration: %d\n",i);
    fprintf(fp2, "average time: %ld\n",average_time);
    
    fclose(fp2);

}






