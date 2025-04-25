#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if(a.rows==b.rows&&a.cols==b.cols)
    {
        Matrix c=create_matrix(a.rows, a.cols);
        for (int i=0;i<a.rows;i++)
        {
            for (int j=0;j<a.cols;j++)
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
        for (int i=0;i<a.rows;i++)
        {
            for (int j=0;j<a.cols;j++)
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
        for (int i=0;i<a.rows;i++)
        {
            for (int j=0;j<b.cols;j++)
            {
                c.data[i][j]=0;
                for (int k=0;k<a.cols;k++)
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
    for (int i=0;i<a.rows;i++)
    {
        for (int j=0;j<a.cols;j++)
        {
            c.data[i][j]=a.data[i][j]*k;
        }
    }
    return c;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix c=create_matrix(a.cols, a.rows);
    for (int i=0;i<a.rows;i++)
    {
        for (int j=0;j<a.cols;j++)
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
        //TODO
    }
    printf("Error: The matrix must be a square matrix.\n");
    return 0;
}

Matrix inv_matrix(Matrix a)
{
    if(a.rows==a.cols&&det_matrix(a)!=0)
    {
        Matrix c=create_matrix(a.rows, a.cols);
        // ToDo
        return c;
    }
    return create_matrix(0, 0);
}

int rank_matrix(Matrix a)
{
    // ToDo
    return 0;
}

double trace_matrix(Matrix a)
{
    if(a.rows==a.cols)
    {
        double trace=0;
        for (int i=0;i<a.rows;i++)
        {
            trace+=a.data[i][i];
        }
        return trace;
    }
    return 0;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}