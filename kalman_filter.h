/**
  ******************************************************************************
  * @file    kalman_filter.h
  * @author  Hongxi Wong
  * @version V1.0.2
  * @date    2020/3/1
  * @brief   
  ******************************************************************************
  * @attention 
  *
  ******************************************************************************
  */
#ifndef __KALMAN_FILTER_H
#define __KALMAN_FILTER_H

//stm32 DSP lib
#ifndef ARM_MATH_CM4
#define ARM_MATH_CM4
#endif
#ifndef __CC_ARM
#define __CC_ARM
#endif
#ifndef ARM_MATH_MATRIX_CHECK
#define ARM_MATH_MATRIX_CHECK
#endif
#ifndef ARM_MATH_ROUNDING
#define ARM_MATH_ROUNDING
#endif

#include <math.h>
#include <stdlib.h>
#include "stm32f4xx.h"
#ifdef ARM_MATH_CM4
#include "arm_math.h"

#ifdef _CMSIS_OS_H
#define user_malloc pvPortMalloc
#else
#define user_malloc malloc
#endif

#define mat arm_matrix_instance_f32
#define Matrix_Init arm_mat_init_f32
#define Matrix_Add arm_mat_add_f32
#define Matrix_Subtract arm_mat_sub_f32
#define Matrix_Multiply arm_mat_mult_f32
#define Matrix_Transpose arm_mat_trans_f32
#define Matrix_Inverse arm_mat_inverse_f32

typedef struct
{
    uint8_t xhat_size;
    float *Raw_Value;
    float *Filtered_Value;

    mat xhat;      //状态向量
    mat xhatminus; //状态向量先验估计
    mat u;         //控制向量
    mat z;         //测量数据向量
    mat P;         //协方差矩阵
    mat Pminus;    //协方差矩阵先验估计
    mat A, AT;     //状态转移矩阵
    mat B;         //控制矩阵
    mat H, HT;     //观测矩阵
    mat Q;         //模型噪声矩阵
    mat R;         //传感器噪声矩阵
    mat K;         //卡尔曼增益
    mat temp_matrix, temp_matrix1, temp_vector, temp_vector1;

    float *xhat_data, *xhatminus_data;
    float *u_data;
    float *z_data;
    float *P_data, *Pminus_data;
    float *A_data, *AT_data;
    float *B_data;
    float *H_data, *HT_data;
    float *Q_data;
    float *R_data;
    float *K_data;
    float *temp_matrix_data, *temp_matrix_data1, *temp_vector_data, *temp_vector_data1;
} kalman_filter_t;

void Kalman_Filter_Init(kalman_filter_t *KF, uint8_t xhat_size);
float *Kalman_Filter_Calculate(kalman_filter_t *KF, float *input);

#endif //ARM_MATH_CM4

#endif //__KALMAN_FILTER_H
