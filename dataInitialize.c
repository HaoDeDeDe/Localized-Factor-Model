//
//  dataInitialize.c
//  MyProjectC
//
//  Created by 杨昊的 on 16/3/9.
//  Copyright © 2016年 杨昊的. All rights reserved.
//



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "quickSort.h"
#include "dataInitialize.h"


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

struct indexCell{
    int num;
    int *index;
};


//convert a two-char string to a number
int strToInt(char *str){
    int length=2;
    int sum=0;
    for (int i=length-1; i>=0; i--) {
        if((str[i]>=48)&&(str[i]<=57)) sum=sum+(str[i]-48)*pow(10, length-i-1);
        else return -1;   
    }
    return sum;
}

//生成1-n的随机序列
void randomSequence(int *random, int n){
    int *flag=(int *)malloc(n*sizeof(int));
    int temp,k=0,i;
    
    for (i=0; i<n; i++) {
        flag[i]=0;
    }
    
    //srand(time(NULL));
    while(k<n){
        temp=rand()%n;
        if (flag[temp]==0) {
            random[k]=temp+1;
            flag[temp]=1;
            k++;
        }
    }
}


void dataInitialize1(struct resMatrix *context, int context_num) {
    
    //* 100000 ratings (1-5) from 943 users on 1682 movies.
    int user_num=943;
    int movie_num=1682;
    int rating_num=100000;
    int *ratings=(int *)malloc(rating_num*4*sizeof(int));
    int i,j,k;
    int *user_id=(int *)malloc(user_num*sizeof(int));
    int *user_age=(int *)malloc(user_num*sizeof(int));
    char *user_gender=(char *)malloc(user_num*sizeof(char));
    char user_job[943][30];
    char *user_zip[943];
    char *temp;
    
    FILE *fp1=fopen("/Users/DaisyYang/Desktop/MyProjectC/MyProjectC/MyProjectC/data1/u.data", "r");
    FILE *fp2=fopen("/Users/DaisyYang/Desktop/MyProjectC/MyProjectC/MyProjectC/data1/u.user", "r" );
    
    //read ratings，user id, movie id, rating, timestamp
    if(fp1!=NULL){
        for(i=0;i<rating_num;i++){
            for(j=0;j<4;j++){
                fscanf(fp1, "%d",ratings+i*4+j);
            }
        }
    }
    
    //read user features
    if(fp2!=NULL){
        for(i=0;i<943;i++){
            fscanf(fp2, "%d|", user_id+i);
            fscanf(fp2, "%d|", user_age+i);
            fscanf(fp2, "%c|", user_gender+i);
            fscanf(fp2, "%[a-zA-Z]", user_job[i]);
            fgetc(fp2);
            temp=(char *)malloc(5*sizeof(char));
            user_zip[i]=(char *)malloc(3*sizeof(char));
            fscanf(fp2, "%5c", temp);
            strncpy(user_zip[i], temp, 2);
            user_zip[i][3]='\0';
            free(temp);
            fgetc(fp2);
        }
    }
    
    fclose(fp1);
    fclose(fp2);
    
    //only movies with top 100 ratings will be used for training
    int *rating_num_ori=(int *)malloc(movie_num*sizeof(int));
    int *rating_num_new=(int *)malloc(movie_num*sizeof(int));
    //统计每个movie被rate的次数
    for (i=0; i<rating_num; i++) {
        rating_num_ori[ratings[i*4+1]-1]++;
        rating_num_new[ratings[i*4+1]-1]++;
    }
    
    //按sort ratings
    quickSort(rating_num_new, 0, movie_num-1);
    
    int rate_point=rating_num_new[movie_num-100];
    int *movie_place=(int *)malloc(movie_num*sizeof(int));
    struct indexCell chosen_movies;
    chosen_movies.num=0;
    for (i=0; i<movie_num; i++) {
        if (rating_num_ori[i]>=rate_point) {
            chosen_movies.num++;
        }
    }
    chosen_movies.index=(int *)malloc(chosen_movies.num*sizeof(int));
    j=0;
    for (i=0; i<movie_num; i++) {
        if (rating_num_ori[i]>=rate_point) {
            chosen_movies.index[j]=i;
            movie_place[i]=j;
            j++;
        }
    }
    
    
    //prune ratings
    int pruned_rating_num=0;
    for (i=0; i<rating_num; i++) {
        if (rating_num_ori[ratings[i*4+1]-1]>=rate_point) {
            pruned_rating_num++;
        }
    }
    int *pruned_ratings=(int *)malloc(pruned_rating_num*4*sizeof(int));
    j=0;
    for (i=0; i<rating_num; i++) {
        if (rating_num_ori[ratings[i*4+1]-1]>=rate_point) {
            for (k=0; k<4; k++) {
                pruned_ratings[j*4+k]=ratings[i*4+k];
                //printf("%d ",pruned_ratings[j*4+k]);
            }
            //printf("\n");
            j++;
        }
    }
    
    
    //initialize the response matrix of context 1
    context[0].user_num=user_num;
    context[0].item_num=chosen_movies.num;
    context[0].min_value=1;
    context[0].max_value=5;
    context[0].matrix=(struct response *)malloc(context[0].user_num*context[0].item_num*sizeof(struct response));
    for(i=0;i<context[0].user_num;i++){
        for(j=0;j<context[0].item_num;j++){
            context[0].matrix[i*context[0].item_num+j].indicator=-1;
        }
    }
    
    
    int *time_stamp=(int *)malloc(pruned_rating_num*sizeof(int));
    int time_point;
    
    for (i=0; i<pruned_rating_num; i++) {
        time_stamp[i]=pruned_ratings[i*4+3];
    }
    
    quickSort(time_stamp,0,pruned_rating_num-1);
    time_point=time_stamp[(int)(pruned_rating_num*0.75)];
    
    //get training set and test set for context 1
    for (i=0; i<pruned_rating_num; i++) {
        if(pruned_ratings[i*4+3]<time_point){
            context[0].matrix[(pruned_ratings[i*4]-1)*context[0].item_num+movie_place[pruned_ratings[i*4+1]-1]].indicator=1;
            context[0].matrix[(pruned_ratings[i*4]-1)*context[0].item_num+movie_place[pruned_ratings[i*4+1]-1]].value=pruned_ratings[i*4+2];
        }
        else{
            context[0].matrix[(pruned_ratings[i*4]-1)*context[0].item_num+movie_place[pruned_ratings[i*4+1]-1]].indicator=0;
            context[0].matrix[(pruned_ratings[i*4]-1)*context[0].item_num+movie_place[pruned_ratings[i*4+1]-1]].value=pruned_ratings[i*4+2];
        }
    }
    
    //printf("user:716 movie:204 rating:%f\n",context[0].matrix[715*context[0].item_num+movie_place[203]].value);

    
    //initialize matrix for context 2
    context[1].user_num=user_num;
    context[1].item_num=130;
    context[1].min_value=0;
    context[1].max_value=1;
    context[1].matrix=(struct response *)malloc(context[1].user_num*context[1].item_num*sizeof(struct response));
    for(i=0;i<context[1].user_num;i++){
        for(j=0;j<context[1].item_num;j++){
            context[1].matrix[i*context[1].item_num+j].indicator=0;
            context[1].matrix[i*context[1].item_num+j].value=0;
        }
    }
    
    //read in occupations
    FILE *fp3=fopen("/Users/DaisyYang/Desktop/MyProjectC/MyProjectC/MyProjectC/data1/u.occupation", "r");
    char jobs[21][30];
    for(i=0;i<21;i++){
        fscanf(fp3, "%[a-zA-Z]",jobs[i]);
        fgetc(fp3);
    }
    
    fclose(fp3);   
    
    
    int *random=(int *)malloc(context[1].item_num*sizeof(int));
    double f=0.4;
    
    for (i=0; i<context[1].user_num; i++) {
        randomSequence(random, context[1].item_num);
        
        for (j=0; j<context[1].item_num; j++) {
            if (random[j]<=(context[1].item_num*f)) {
                context[1].matrix[i*context[1].item_num+j].indicator=1;
            }
        }
        
        if (user_age[i]>=0&&user_age[i]<=9) {
            context[1].matrix[i*context[1].item_num+0].value=1;
        }
        else if (user_age[i]>=10&&user_age[i]<=19){
            context[1].matrix[i*context[1].item_num+1].value=1;
        }
        else if (user_age[i]>=20&&user_age[i]<=29){
            context[1].matrix[i*context[1].item_num+2].value=1;
        }
        else if (user_age[i]>=30&&user_age[i]<=39){
            context[1].matrix[i*context[1].item_num+3].value=1;
        }
        else if (user_age[i]>=40&&user_age[i]<=49){
            context[1].matrix[i*context[1].item_num+4].value=1;
        }
        else {
            context[1].matrix[i*context[1].item_num+5].value=1;
        }
        
        if (user_gender[i]=='M') {
            context[1].matrix[i*context[1].item_num+6].value=1;
        }
        else{
            context[1].matrix[i*context[1].item_num+7].value=1;
        }
        
        for (j=0; j<21; j++) {
            if (strcmp(user_job[i], jobs[j])==0) {
                context[1].matrix[i*context[1].item_num+8+j].value=1;
                break;
            }
        }
        
        if (strToInt(user_zip[i])==-1) {
            context[1].matrix[i*context[1].item_num+129].value=1;
        }
        else{
            for (j=0; j<100; j++) {
                if (strToInt(user_zip[i])==j) {
                    context[1].matrix[i*context[1].item_num+29+j].value=1;
                }
            }
        }
    }
    
    
}


void dataInitialize2(struct resMatrix *context, int context_num) {
    
    //* 1,000,209 ratings (1-5) from 6040 users on 3952 movies.
    int user_num=6040;
    int movie_num=3952;
    int rating_num=1000209;
    
    int i,j,k;
    int *ratings=(int *)malloc(rating_num*4*sizeof(int));
    int *user_id=(int *)malloc(user_num*sizeof(int));
    int *user_age=(int *)malloc(user_num*sizeof(int));
    int *user_job=(int *)malloc(user_num*sizeof(int));
    char *user_gender=(char *)malloc(user_num*sizeof(char));
    char *user_zip[6040];
    char *temp;
    
    FILE *fp1=fopen("/Users/DaisyYang/Desktop/MyProjectC/MyProjectC/MyProjectC/data2/ratings.dat", "r");
    FILE *fp2=fopen("/Users/DaisyYang/Desktop/MyProjectC/MyProjectC/MyProjectC/data2/users.dat", "r" );
    
    //read movie ratings，user id, movie id, rating, timestamp
    if(fp1!=NULL){
        for(i=0;i<rating_num;i++){
            for(j=0;j<3;j++){
                fscanf(fp1, "%d::",ratings+i*4+j);
            }
            fscanf(fp1, "%d",ratings+i*4+3);
        }
    }
    
    //read user features
    if(fp2!=NULL){
        for(i=0;i<user_num;i++){
            fscanf(fp2, "%d::", user_id+i);
            //printf("%d ",user_id[i]);
            
            fscanf(fp2, "%c::", user_gender+i);
            //printf("%c ",user_gender[i]);
            
            fscanf(fp2, "%d::", user_age+i);
            //printf("%d ",user_age[i]);
            
            fscanf(fp2, "%d::", user_job+i);
            //printf("%d ",user_job[i]);
            
            temp=(char *)malloc(20*sizeof(char));
            user_zip[i]=(char *)malloc(3*sizeof(char));
            fgets(temp, 20, fp2);
            strncpy(user_zip[i], temp, 2);
            user_zip[i][3]='\0';
            //printf("%s\n",user_zip[i]);
            
            free(temp);
        }
    }
    
    fclose(fp1);
    fclose(fp2);
    
    int *rating_num_ori=(int *)malloc(movie_num*sizeof(int));
    int *rating_num_new=(int *)malloc(movie_num*sizeof(int));

    for (i=0; i<rating_num; i++) {
        rating_num_ori[ratings[i*4+1]-1]++;
        rating_num_new[ratings[i*4+1]-1]++;
    }
    

    quickSort(rating_num_new, 0, movie_num-1);
    

    int rate_point=rating_num_new[movie_num-100];
    int *movie_place=(int *)malloc(movie_num*sizeof(int));
    struct indexCell chosen_movies;
    chosen_movies.num=0;
    for (i=0; i<movie_num; i++) {
        if (rating_num_ori[i]>=rate_point) {
            chosen_movies.num++;
        }
    }
    chosen_movies.index=(int *)malloc(chosen_movies.num*sizeof(int));
    j=0;
    for (i=0; i<movie_num; i++) {
        if (rating_num_ori[i]>=rate_point) {
            chosen_movies.index[j]=i;
            movie_place[i]=j;
            j++;
        }
    }
    
    
    //prune ratings
    int pruned_rating_num=0;
    for (i=0; i<rating_num; i++) {
        if (rating_num_ori[ratings[i*4+1]-1]>=rate_point) {
            pruned_rating_num++;
        }
    }
    int *pruned_ratings=(int *)malloc(pruned_rating_num*4*sizeof(int));
    j=0;
    for (i=0; i<rating_num; i++) {
        if (rating_num_ori[ratings[i*4+1]-1]>=rate_point) {
            for (k=0; k<4; k++) {
                pruned_ratings[j*4+k]=ratings[i*4+k];
                //printf("%d ",pruned_ratings[j*4+k]);
            }
            //printf("\n");
            j++;
        }
    }
    
    //printf("movie:1097 place:%d\n",movie_place[1096]);
    
    //initialize the response matrix of context 1
    context[0].user_num=user_num;
    context[0].item_num=chosen_movies.num;
    context[0].min_value=1;
    context[0].max_value=5;
    context[0].matrix=(struct response *)malloc(context[0].user_num*context[0].item_num*sizeof(struct response));
    for(i=0;i<context[0].user_num;i++){
        for(j=0;j<context[0].item_num;j++){
            context[0].matrix[i*context[0].item_num+j].indicator=-1;
        }
    }
    
    int *time_stamp=(int *)malloc(pruned_rating_num*sizeof(int));
    int time_point;
    
    for (i=0; i<pruned_rating_num; i++) {
        time_stamp[i]=pruned_ratings[i*4+3];
    }
    

    quickSort(time_stamp,0,pruned_rating_num-1);
    time_point=time_stamp[(int)(pruned_rating_num*0.75)];
    
    //get training set and test set for context 1
    for (i=0; i<pruned_rating_num; i++) {
        if(pruned_ratings[i*4+3]<time_point){
            context[0].matrix[(pruned_ratings[i*4]-1)*context[0].item_num+movie_place[pruned_ratings[i*4+1]-1]].indicator=1;
            context[0].matrix[(pruned_ratings[i*4]-1)*context[0].item_num+movie_place[pruned_ratings[i*4+1]-1]].value=pruned_ratings[i*4+2];
        }
        else{
            context[0].matrix[(pruned_ratings[i*4]-1)*context[0].item_num+movie_place[pruned_ratings[i*4+1]-1]].indicator=0;
            context[0].matrix[(pruned_ratings[i*4]-1)*context[0].item_num+movie_place[pruned_ratings[i*4+1]-1]].value=pruned_ratings[i*4+2];
        }
    }
    
    //printf("user:6036 movie:2997 rating:%f\n",context[0].matrix[6035*context[0].item_num+movie_place[2996]].value);
    
    
    //initialize matrix for context 2
    context[1].user_num=user_num;
    context[1].item_num=130;
    context[1].min_value=0;
    context[1].max_value=1;
    context[1].matrix=(struct response *)malloc(context[1].user_num*context[1].item_num*sizeof(struct response));
    for(i=0;i<context[1].user_num;i++){
        for(j=0;j<context[1].item_num;j++){
            context[1].matrix[i*context[1].item_num+j].indicator=-1;
            context[1].matrix[i*context[1].item_num+j].value=0;
        }
    }
    
    
    int *random=(int *)malloc(user_num*sizeof(int));
    randomSequence(random, user_num);
    double f=0.05;
    
    for (i=0; i<context[1].user_num; i++) {
        if (random[i]<=(int)(user_num*f)) {  
            for (j=0; j<130; j++) {
                context[1].matrix[i*context[1].item_num+j].indicator=1;
            }
        }
        
        if (user_gender[i]=='M') {
            context[1].matrix[i*context[1].item_num+0].value=1;
        }
        else{
            context[1].matrix[i*context[1].item_num+1].value=1;
        }
        
        if (user_age[i]==1) {
            context[1].matrix[i*context[1].item_num+2].value=1;
        }
        else if (user_age[i]==18){
            context[1].matrix[i*context[1].item_num+3].value=1;
        }
        else if (user_age[i]==25){
            context[1].matrix[i*context[1].item_num+4].value=1;
        }
        else if (user_age[i]==35){
            context[1].matrix[i*context[1].item_num+5].value=1;
        }
        else if (user_age[i]==45){
            context[1].matrix[i*context[1].item_num+6].value=1;
        }
        else if (user_age[i]==50){
            context[1].matrix[i*context[1].item_num+7].value=1;
        }
        else {
            context[1].matrix[i*context[1].item_num+8].value=1;
        }
        
        for (j=0; j<=20; j++) {
            if (user_job[i]==j) {
                context[1].matrix[i*context[1].item_num+9+j].value=1;
            }
        }
        
        for (j=0; j<100; j++) {
            if (strToInt(user_zip[i])==j) {
                context[1].matrix[i*context[1].item_num+30+j].value=1;
            }
        }
    }
    
    
}
