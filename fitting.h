//
//  fitting.h
//  MyProjectC
//
//  Created by 杨昊的 on 16/3/10.
//  Copyright © 2016年 杨昊的. All rights reserved.
//

#ifndef fitting_h
#define fitting_h

#include <stdio.h>

//a single cell in a response matrix
struct response;

//response matrix or feature matrix
struct resMatrix;

struct vector;

struct matrix;

struct indexCell;

double var_uniform(double min, double max);

double normal_pdf(double x, double mean, double std);

double var_normal(double min, double max, double mean, double std);

void paraInitialize(int context_num, struct resMatrix *context,double *alpha, double *sigma2_z, double *sigma2_y, double *sigma2_u, struct vector *beta, struct matrix *w, struct matrix *A, int *q, int *r);

void EStep(int context_num, struct resMatrix *context, double *alpha, double *sigma2_z, double *sigma2_y, double *sigma2_u, struct vector *beta, struct vector *Eu, struct vector *Ez, struct matrix *w, struct matrix *A, struct matrix *Vu, struct matrix *Vz, struct matrix *Covzu, int *q, int *r, struct indexCell *J);

void MStep(int context_num, struct resMatrix *context, double *alpha, double *sigma2_z, double *sigma2_y, double *sigma2_u, struct vector *beta, struct vector *Eu, struct vector *Ez, struct matrix *w, struct matrix *A, struct matrix *Vu, struct matrix *Vz, struct matrix *Covzu, int *q, int *r, struct indexCell *J);

double RMSE_training(int context_num, struct resMatrix *context, double *alpha, struct vector *beta, struct vector *Ez, int *q, int *r);

double RMSE_predict(int context_num, struct resMatrix *context, double *alpha, struct vector *beta, struct vector *Ez, int *q, int *r);

void fitting(int context_num, struct resMatrix *context,double *alpha, double *sigma2_z, double *sigma2_y, double *sigma2_u, struct vector *beta, struct vector *Eu, struct vector *Ez, struct matrix *w, struct matrix *A, struct matrix *Vu, struct matrix *Vz, struct matrix *Covzu, int *q, int *r);




#endif /* fitting_h */
