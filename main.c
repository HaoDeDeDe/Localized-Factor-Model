//
//  main.c
//  MyProjectC
//
//  Created by 杨昊的 on 16/3/2.
//  Copyright © 2016年 杨昊的. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matrixOperation.h"
#include "quickSort.h"
#include "dataInitialize.h"
#include "fitting.h"


//a single cell in a response matrix
struct response{
    int indicator;  //-1代表没有值，0代表predict, 1代表train
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


//* 100000 ratings (1-5) from 943 users on 1682 movies.
int main(int argc, const char * argv[]) {

    int i,j,temp,context_num=2;
    struct resMatrix *context;
    double *alpha, *sigma2_z, *sigma2_y,*sigma2_u;
    struct vector *beta, *Eu, *Ez;
    struct matrix *w, *A, *Vu, *Vz, *Covzu;
    int *q, *r;
    
    
    context=(struct resMatrix *)malloc(context_num*sizeof(struct resMatrix));
    dataInitialize2(context, context_num);
    
    
    //检查是否每个context的response matrix拥有相同数量的行
    if (context_num>1) {
        for (i=0; i<context_num-1; i++) {
            if (context[i].user_num!=context[i+1].user_num) {
                printf("false construction of response matrix - unequal number of users");
                return -1;
            }
        }
    }
    
    //为各变量和参数申请空间
    temp=0;
    for (i=0; i<context_num; i++) {
        temp=temp+context[i].item_num;
    }
    alpha=(double *)malloc(temp*sizeof(double));
    beta=(struct vector *)malloc(temp*sizeof(struct vector));
    sigma2_y=(double *)malloc(context_num*sizeof(double));  //此处为sigma squared
    sigma2_z=(double *)malloc(context_num*sizeof(double));  //同上
    sigma2_u=(double *)malloc(sizeof(double));              //同上
    w=(struct matrix *)malloc(context_num*sizeof(struct matrix));
    A=(struct matrix *)malloc(context_num*sizeof(struct matrix));
    q=(int *)malloc(context_num*sizeof(int));
    r=(int *)malloc(sizeof(int));
    
    Eu=(struct vector *)malloc(context[0].user_num*sizeof(struct vector));
    Vu=(struct matrix *)malloc(context[0].user_num*sizeof(struct matrix));
    Ez=(struct vector *)malloc(context_num*context[0].user_num*sizeof(struct vector));
    Vz=(struct matrix *)malloc(context_num*context[0].user_num*sizeof(struct matrix));
    Covzu=(struct matrix *)malloc(context_num*context[0].user_num*sizeof(struct matrix));
    
    //fitting algorithm
    fitting(context_num, context, alpha, sigma2_z, sigma2_y, sigma2_u, beta, Eu, Ez, w, A, Vu, Vz, Covzu, q, r);
    
    
    //printf("%f\n",context[0].matrix[195*1682+241].value);
    
    
}
 
