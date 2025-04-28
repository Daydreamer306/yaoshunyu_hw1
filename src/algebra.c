#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows=row;
    m.cols=col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if(a.rows==b.rows&&a.cols==b.cols)
    {
        Matrix c=create_matrix(a.rows, a.cols);
        for(int i=0;i<a.rows;i++)
        {
            for(int j=0;j<a.cols;j++)
            {
                c.data[i][j]=a.data[i][j]+b.data[i][j];
            }
        }
        return c;
    }
    // 如果行数或列数不相等，给出错误提示并返回一个空矩阵
    printf("Error: Matrix a and b must have the same rows and cols.\n");
    return create_matrix(0, 0);
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if(a.rows==b.rows&&a.cols==b.cols)
    {
        Matrix c=create_matrix(a.rows, a.cols);
        for(int i=0;i<a.rows;i++)
        {
            for(int j=0;j<a.cols;j++)
            {
                c.data[i][j]=a.data[i][j]-b.data[i][j];
            }
        }
        return c;
    }
    // 如果行数或列数不相等，给出错误提示并返回一个空矩阵
    printf("Error: Matrix a and b must have the same rows and cols.\n");
    return create_matrix(0, 0);
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if(a.cols==b.rows)
    {
        Matrix c=create_matrix(a.rows, b.cols);
        for(int i=0;i<a.rows;i++)
        {
            for(int j=0;j<b.cols;j++)
            {
                c.data[i][j]=0;
                for(int k=0;k<a.cols;k++)
                {
                    c.data[i][j]+=a.data[i][k]*b.data[k][j];
                }
            }
        }
        return c;
    }
    // 如果a的列数不等于b的行数，给出错误提示并返回一个空矩阵
    printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
    return create_matrix(0, 0);
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix c=create_matrix(a.rows, a.cols);
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<a.cols;j++)
        {
            c.data[i][j]=a.data[i][j]*k;
        }
    }
    return c;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix c=create_matrix(a.cols, a.rows);
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<a.cols;j++)
        {
            c.data[j][i]=a.data[i][j];
        }
    }
    return c;
}

double det_matrix(Matrix a)
{
    if(a.rows==a.cols)
    {
        if(a.rows==1)
        {
            return a.data[0][0];
        }
        else if(a.rows==2)
        {
            return a.data[0][0]*a.data[1][1]-a.data[0][1]*a.data[1][0];
        }
        else
        {
            double det=0;
            for(int i=0;i<a.cols;i++)
            {
                Matrix b=create_matrix(a.rows-1, a.cols-1);
                for(int j=1;j<a.rows;j++)
                {
                    for(int k=0;k<a.cols;k++)
                    {
                        if(k<i)
                        {
                            b.data[j-1][k]=a.data[j][k];
                        }
                        else if(k>i)
                        {
                            b.data[j-1][k-1]=a.data[j][k];
                        }
                    }
                }
                det+=pow(-1, i)*a.data[0][i]*det_matrix(b);
            }
            return det;
        }
    }
    printf("Error: The matrix must be a square matrix.\n");
    return 0;
}

Matrix inv_matrix(Matrix a)
{
    double det=det_matrix(a);
    if(a.rows==a.cols&&det!=0)
    {
        Matrix c=create_matrix(a.rows, a.cols);
        for(int i=0;i<a.rows;i++)
        {
            for(int j=0;j<a.cols;j++)
            {
                Matrix b=create_matrix(a.rows-1, a.cols-1);
                for(int k=0;k<a.rows;k++)
                {
                    for(int l=0;l<a.cols;l++)
                    {
                        if(k<i&&l<j)
                        {
                            b.data[k][l]=a.data[k][l];
                        }
                        else if(k<i&&l>j)
                        {
                            b.data[k][l-1]=a.data[k][l];
                        }
                        else if(k>i&&l<j)
                        {
                            b.data[k-1][l]=a.data[k][l];
                        }
                        else if(k>i&&l>j)
                        {
                            b.data[k-1][l-1]=a.data[k][l];
                        }
                    }
                }
                c.data[j][i]=pow(-1, i+j)*det_matrix(b)/det;
            }
        }
        return c;
    }
    return create_matrix(0, 0);
}

int rank_matrix(Matrix A)
{
    Matrix a=A;
    int rank=0;
    double EPSILON=1e-10;//定义一个很小的数，用于判断接近零的情况
    for(int i=0;i<a.rows;i++)
    {
       int max_row=i;
        for(int k=i+1;k<a.rows;k++)
        {
            if(fabs(a.data[k][i])>fabs(a.data[max_row][i]))
            {
                max_row=k;
            }
        }
        
        //如果当前列所有元素都为零，则秩不增加
        if(fabs(a.data[max_row][i])<EPSILON)
            continue;
            
        //交换行
        double* temp=a.data[i];
        for (int j=0;j<a.cols;j++) 
        {
            double temp = a.data[i][j];
            a.data[i][j] = a.data[max_row][j];
            a.data[max_row][j] = temp;
        }      
        rank++;
            
        //消去下面的行
        for(int k=i + 1;k<a.rows;k++)
        {
            double factor=a.data[k][i]/a.data[i][i];
            for(int j=i;j<a.cols;j++){
                a.data[k][j]-=factor*a.data[i][j];
            }
        }
    }
        
    return rank;
}

double trace_matrix(Matrix a)
{
    if(a.rows==a.cols)
    {
        double trace=0;
        for(int i=0;i<a.rows;i++)
        {
            trace+=a.data[i][i];
        }
        return trace;
    }
    return 0;
}

void print_matrix(Matrix a)
{
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<a.cols;j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}