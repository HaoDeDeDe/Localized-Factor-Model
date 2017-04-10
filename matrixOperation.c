//
//  matrixOperation.c
//  MyProjectC
//
//  Created by 杨昊的 on 16/3/3.
//  Copyright © 2016年 杨昊的. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "matrixOperation.h"


/*int matrixAddition(double *A, double *B, int Am, int An, int Bm, int Bn, double *C){
    if((Am!=Bm)||(An!=Bn)) return 0;
    if((A==NULL)||(B==NULL)) return 0;
    
    int i,j,m=Am,n=An;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            C[i*n+j]=A[i*n+j]+B[i*n+j];
        }
    }
    return 1;
}*/

double * matrixAddition(double *A, double *B, int Am, int An, int Bm, int Bn){
    if((Am!=Bm)||(An!=Bn)) return NULL;
    if((A==NULL)||(B==NULL)) return NULL;
    
    int i,j,m=Am,n=An;
    double *C=(double *)malloc(m*n*sizeof(double));
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            C[i*n+j]=A[i*n+j]+B[i*n+j];
        }
    }
    return C;
}

/*int matrixSubtraction(double *A, double *B, int Am, int An, int Bm, int Bn, double *C){
    if((Am!=Bm)||(An!=Bn)) return 0;
    if((A==NULL)||(B==NULL)) return 0;
 
    int i,j,m=Am,n=An;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            C[i*n+j]=A[i*n+j]-B[i*n+j];
        }
    }
    return 1;
}*/


double * matrixSubtraction(double *A, double *B, int Am, int An, int Bm, int Bn){
    if((Am!=Bm)||(An!=Bn)) return NULL;
    if((A==NULL)||(B==NULL)) return NULL;
    
    int i,j,m=Am,n=An;
    double *C=(double *)malloc(m*n*sizeof(double));
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            C[i*n+j]=A[i*n+j]-B[i*n+j];
        }
    }
    return C;
}

/*int matrixMultiply(double *A, double *B, int Am, int An, int Bm, int Bn, double *C){
    if(An!=Bm) return 0;
    if((A==NULL)||(B==NULL)) return 0;
    
    int i,j,k;
    double sum=0;
    
    for(i=0;i<Am;i++){
        for(j=0;j<Bn;j++){
            for(k=0;k<An;k++){
                sum=sum+A[i*An+k]*B[k*Bn+j];
            }
            C[i*Bn+j]=sum;
            sum=0;
        }
    }
   */

double * matrixMultiply(double *A, double *B, int Am, int An, int Bm, int Bn){
    if(An!=Bm) return NULL;
    if((A==NULL)||(B==NULL)) return NULL;
    
    int i,j,k;
    double sum=0;
    double *C=(double *)malloc(Am*Bn*sizeof(double));
    
    for(i=0;i<Am;i++){
        for(j=0;j<Bn;j++){
            for(k=0;k<An;k++){
                sum=sum+A[i*An+k]*B[k*Bn+j];
            }
            C[i*Bn+j]=sum;
            sum=0;
        }
    }
    return C;
}

/*int constMatrixMultiply(double c, double * A, int m, int n){   //常数乘以矩阵
    if (A==NULL) {
        return 0;
    }
    
    int i, j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            A[i*n+j]=A[i*n+j]*c;
        }
    }
    
    return 1;
}*/

double * constMatrixMultiply(double c, double * A, int m, int n){  //常数乘以矩阵
    if (A==NULL) {
        return NULL;
    }
    
    int i, j;
    double * matrix=(double *)malloc(m*n*sizeof(double));
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            matrix[i*n+j]=A[i*n+j]*c;
        }
    }
    
    return matrix;
}

//构造const倍的单位矩阵
/*int constIdentityMatrix(double c, double *I, int n){
    if (n==0) {
        return 0;
    }
    
    int i,j;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if (i==j) {
                I[i*n+j]=c;
            }
            else{
                I[i*n+j]=0;
            }
        }
    }
    
    return 1;
}*/

double * constIdentityMatrix(double c, int n){
    if (n==0) {
        return NULL;
    }
    
    int i,j;
    double *I=(double *)malloc(n*n*sizeof(double));
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if (i==j) {
                I[i*n+j]=c;
            }
            else{
                I[i*n+j]=0;
            }
        }
    }
    
    return I;
}

//矩阵转置
/*int matrixTransit(double *A, int m, int n, double *B){
    if ((A==NULL)||(m==0)||(n==0)) {
        return 0;
    }
    
    int i,j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            B[j*m+i]=A[i*n+j];
        }
    }
    
    return 1;
}*/

double * matrixTransit(double *matrix, int m, int n){
    if ((matrix==NULL)||(m==0)||(n==0)) {
        return NULL;
    }
    
    int i,j;
    double *M=(double *)malloc(m*n*sizeof(double));
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            M[j*m+i]=matrix[i*n+j];
        }
    }
    
    return M;
}

//矩阵加一行和一列全0
double * I_zero(double * original, int m, int n){
    double *new=(double *)malloc((m+1)*(n+1)*sizeof(double));
    int i,j;
    
    for (j=0; j<(n+1); j++) {
        new[j]=0;
    }
    for (i=1; i<(m+1); i++) {
        new[i*(n+1)]=0;
    }
    
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            new[(i+1)*(n+1)+(j+1)]=original[i*n+j];
        }
    }
    
    return new;
}

//向量加一个元素1
double * V_one(double * original, int n){
    double *new=(double *)malloc((n+1)*sizeof(double));
    new[0]=1;
    
    for (int i=0; i<n; i++) {
        new[i+1]=original[i];
    }
    
    return new;
}

/*int matrixInverse(double *original, double *inverse, int m, int n){   //矩阵求逆
    double max;
    int *is,*js,i,j,k,mi=-1,mj=-1;
    
    //检查是否是正方形矩阵
    if(n!=m) return 0;
    if(original==NULL) return 0;
    
    is=(int *)malloc(n*sizeof(int));
    js=(int *)malloc(n*sizeof(int));
    
    //复制矩阵
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            inverse[i*n+j]=original[i*n+j];
        }
    }
    
    for(k=0;k<n;k++){
        
        //找到右下子阶中最大的元素
        max=inverse[k*n+k];
        for(i=k;i<n;i++){
            for(j=k;j<n;j++){
                if(inverse[i*n+j]>=max){
                    max=inverse[i*n+j];
                    mi=i;	//记录最大元素所在行列
                    mj=j;
                }
            }
        }
        if((max+1.0)==1.0){
            free(is);
            free(js);
            return 0;
        }
        
        //移动主元素
        if(mi==k)	is[k]=-1;
        else{
            swapRows(inverse,mi,k,n);
            is[k]=mi;
        }
        if(mj==k)	js[k]=-1;
        else{
            swapColumns(inverse,mj,k,n);
            js[k]=mj;
        }
        
        //各位置元素变换
        inverse[k*n+k]=1/inverse[k*n+k];	//step 1
        for(j=0;j<n;j++){	//step 2
            if(j!=k){
                inverse[k*n+j]=inverse[k*n+j]*inverse[k*n+k];
            }
        }
        for(i=0;i<n;i++){	//step 3
            for(j=0;j<n;j++){
                if((i!=k)&&(j!=k)){
                    inverse[i*n+j]=inverse[i*n+j]-inverse[i*n+k]*inverse[k*n+j];
                }
            }
        }
        for(i=0;i<n;i++){	//step 4
            if(i!=k){
                inverse[i*n+k]=-inverse[i*n+k]*inverse[k*n+k];
            }
        }
        
    }
    
    //恢复主元
    for(k=(n-1);k>=0;k--){
        if(js[k]!=(-1))	swapRows(inverse,k,js[k],n);
        if(is[k]!=(-1))	swapColumns(inverse,k,is[k],n);
    }
    
    free(is);
    free(js);
    
    return 1;
}*/

double * matrixInverse(double *original, int m, int n){   //矩阵求逆
    double max;
    int *is,*js,i,j,k,mi=-1,mj=-1;
    
    //检查是否是正方形矩阵
    if(n!=m) return NULL;
    if(original==NULL) return NULL;
    
    is=(int *)malloc(n*sizeof(int));
    js=(int *)malloc(n*sizeof(int));
    
    double *inverse=(double *)malloc(n*n*sizeof(double));
    
    //复制矩阵
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            inverse[i*n+j]=original[i*n+j];
        }
    }
    
    for(k=0;k<n;k++){
        
        //找到右下子阶中最大的元素
        max=inverse[k*n+k];
        for(i=k;i<n;i++){
            for(j=k;j<n;j++){
                if(inverse[i*n+j]>=max){
                    max=inverse[i*n+j];
                    mi=i;	//记录最大元素所在行列
                    mj=j;
                }
            }
        }
        if((max+1.0)==1.0){
            free(is);
            free(js);
            return NULL;
        }
        
        //移动主元素
        if(mi==k)	is[k]=-1;
        else{
            swapRows(inverse,mi,k,n);
            is[k]=mi;
        }
        if(mj==k)	js[k]=-1;
        else{
            swapColumns(inverse,mj,k,n);
            js[k]=mj;
        }
        
        //各位置元素变换
        inverse[k*n+k]=1/inverse[k*n+k];	//step 1
        for(j=0;j<n;j++){	//step 2
            if(j!=k){
                inverse[k*n+j]=inverse[k*n+j]*inverse[k*n+k];
            }
        }
        for(i=0;i<n;i++){	//step 3
            for(j=0;j<n;j++){
                if((i!=k)&&(j!=k)){
                    inverse[i*n+j]=inverse[i*n+j]-inverse[i*n+k]*inverse[k*n+j];
                }
            }
        }
        for(i=0;i<n;i++){	//step 4
            if(i!=k){
                inverse[i*n+k]=-inverse[i*n+k]*inverse[k*n+k];
            }
        }
        
    }
    
    //恢复主元
    for(k=(n-1);k>=0;k--){
        if(js[k]!=(-1))	swapRows(inverse,k,js[k],n);
        if(is[k]!=(-1))	swapColumns(inverse,k,is[k],n);
    }
    
    free(is);
    free(js);
    
    return inverse;
}

void swapRows(double *matrix, int x, int y, int n){   //矩阵行互换
    if(matrix==NULL) return;
    
    double temp;
    int j;
    for(j=0;j<n;j++){
        temp=matrix[x*n+j];
        matrix[x*n+j]=matrix[y*n+j];
        matrix[y*n+j]=temp;
    }
    return;
}

void swapColumns(double *matrix, int x, int y, int n){   //矩阵列互换
    if(matrix==NULL) return;
    
    double temp;
    int i;
    for(i=0;i<n;i++){
        temp=matrix[i*n+x];
        matrix[i*n+x]=matrix[i*n+y];
        matrix[i*n+y]=temp;
    }
    return;
}




