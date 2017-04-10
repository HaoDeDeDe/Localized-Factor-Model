//
//  matrixOperation.h
//  MyProjectC
//
//  Created by 杨昊的 on 16/3/3.
//  Copyright © 2016年 杨昊的. All rights reserved.
//

#ifndef matrixOperation_h
#define matrixOperation_h

#include <stdio.h>

//矩阵加法
//int matrixAddition(double *A, double *B, int Am, int An, int Bm, int Bn, double *C);
double * matrixAddition(double *A, double *B, int Am, int An, int Bm, int Bn);

//矩阵减法
//int matrixSubtraction(double *A, double *B, int Am, int An, int Bm, int Bn, double *C);
double * matrixSubtraction(double *A, double *B, int Am, int An, int Bm, int Bn);

//矩阵乘法
//int matrixMultiply(double *A, double *B, int Am, int An, int Bm, int Bn, double *C);
double * matrixMultiply(double *A, double *B, int Am, int An, int Bm, int Bn);

//常数乘以矩阵
//int constMatrixMultiply(double c, double * A, int m, int n);
double * constMatrixMultiply(double c, double * A, int m, int n);

//构造const倍的单位矩阵
//int constIdentityMatrix(double c, double *I, int n);
double * constIdentityMatrix(double c, int n);

//矩阵转置
//int matrixTransit(double *A, int m, int n, double *B);
double * matrixTransit(double *matrix, int m, int n);

//矩阵加一行和一列全0
double * I_zero(double * original, int m, int n);

//向量加一个元素1
double * V_one(double * original, int n);

//矩阵求逆
//int matrixInverse(double *original, double *inverse, int m, int n);
double * matrixInverse(double *original, int m, int n);
    //矩阵行互换
    void swapRows(double *matrix, int x, int y, int n);
    //矩阵列呼唤
    void swapColumns(double *matrix, int x, int y, int n);





#endif /* matrixOperation_h */
