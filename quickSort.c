//
//  quickSort.c
//  MyProjectC
//
//  Created by 杨昊的 on 16/3/8.
//  Copyright © 2016年 杨昊的. All rights reserved.
//

#include "quickSort.h"
#include <stdio.h>

int partition(int ary[],int left,int right)
{
    int z,y,k;
    z=left;
    y=right;
    k=ary[left];
    while(0==0)
    {
        if(y==z)
        {
            ary[y]=k;
            return y;
        };
        while(y>z)
        {
            if(ary[y]<k)
            {
                ary[z]=ary[y];
                z=z+1;
                break;
            }
            if(ary[y]>=k) y=y-1;
        };
        while(z<y)
        {
            if(ary[z]>k)
            {
                ary[y]=ary[z];
                y=y-1;
                break;
            };
            if(ary[z]<=k) z=z+1;
        };
    };
}

void quickSort(int ary[],int left,int right)
{
    int m;
    if(left>=right) return;
    m=partition(ary,left,right);
    quickSort(ary,left,m-1);
    quickSort(ary,m+1,right);
    return;
}

/*int main()
{
    int n,i=0,a[100];
    scanf("%d",&n);
    while(i<n)
    {
        scanf("%d",&a[i]);
        i=i+1;
    };
    quicksort(a,0,n-1);
    i=0;
    while(i<n)
    {
        printf("%d ",a[i]);
        i=i+1;
    };
    return 0;
    
}*/
