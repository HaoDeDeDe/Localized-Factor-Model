//
//  dataInitialize.h
//  MyProjectC
//
//  Created by 杨昊的 on 16/3/9.
//  Copyright © 2016年 杨昊的. All rights reserved.
//

#ifndef dataInitialize_h
#define dataInitialize_h

#include <stdio.h>

//a single cell in a response matrix
struct response;

//response matrix or feature matrix
struct resMatrix;

struct indexCell;

int strToInt(char *str);

void randomSequence(int *random, int n);

void dataInitialize1(struct resMatrix *context, int context_num);

void dataInitialize2(struct resMatrix *context, int context_num);

#endif /* dataInitialize_h */
