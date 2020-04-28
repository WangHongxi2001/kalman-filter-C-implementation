/**
  ******************************************************************************
  * @file    kalman_filter.h
  * @author  Hongxi Wong
  * @version V1.0.3
  * @date    2020/4/27
  * @brief   
  ******************************************************************************
  * @attention 
  *
  ******************************************************************************
  */
#ifndef __KALMAN_FILTER_H
#define __KALMAN_FILTER_H

//cortex-m4 DSP lib
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
    float *Raw_Value;
    float *Filtered_Value;

    uint8_t xhat_size;
    uint8_t u_size;
    uint8_t z_size;
    uint8_t Use_Auto_Adjustment;
    uint8_t *Measurement_Reference; //how measurement relates to the state vector
    float *Measurement_Degree;      //how much directly do sensors measure states
    float *Mat_R_Diagonal_Elements; //variance for each measurement

    mat xhat;      //x(k|k)
    mat xhatminus; //x(k|k-1)
    mat u;         //control vector u
    mat z;         //measurement vector z
    mat z_buf;     //measurement vector z for update
    mat P;         //covariance matrix P(k|k)
    mat Pminus;    //covariance matrix P(k|k-1)
    mat A, AT;     //state transition matrix A AT
    mat B;         //control matrix B
    mat H, HT;     // measurement matrix H
    mat Q;         //process noise covariance matrix Q
    mat R;         //measurement noise covariance matrix R
    mat K;         //kalman gain  K
    mat temp_matrix, temp_matrix1, temp_vector, temp_vector1;

    float *xhat_data, *xhatminus_data;
    float *u_data;
    float *z_data;
    float *z_buf_data;
    float *P_data, *Pminus_data;
    float *A_data, *AT_data;
    float *B_data;
    float *H_data, *HT_data;
    float *Q_data;
    float *R_data;
    float *K_data;
    float *temp_matrix_data, *temp_matrix_data1, *temp_vector_data, *temp_vector_data1;
} kalman_filter_t;

void Kalman_Filter_Init(kalman_filter_t *KF, uint8_t xhat_size, uint8_t u_size, uint8_t z_size);
float *Kalman_Filter_Update(kalman_filter_t *KF);

#endif //ARM_MATH_CM4

#endif //__KALMAN_FILTER_H
