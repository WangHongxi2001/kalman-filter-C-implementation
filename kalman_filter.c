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
#include "kalman_filter.h"
#ifdef ARM_MATH_CM4
#include "arm_math.h"

void Kalman_Filter_Init(kalman_filter_t *KF, uint8_t xhat_size)
{
    KF->xhat_size = xhat_size;

    //滤波数据分配空间
    KF->Raw_Value = (float *)user_malloc(sizeof(float) * xhat_size);
    memset(KF->Raw_Value, 0, xhat_size);
    KF->Filtered_Value = (float *)user_malloc(sizeof(float) * xhat_size);
    memset(KF->Filtered_Value, 0, xhat_size);

    //状态向量初始化
    KF->xhat_data = (float *)user_malloc(sizeof(float) * xhat_size);
    memset(KF->xhat_data, 0, xhat_size);
    Matrix_Init(&KF->xhat, KF->xhat_size, 1, (float *)KF->xhat_data);

    //状态先验估计向量初始化
    KF->xhatminus_data = (float *)user_malloc(sizeof(float) * xhat_size);
    memset(KF->xhatminus_data, 0, xhat_size);
    Matrix_Init(&KF->xhatminus, KF->xhat_size, 1, (float *)KF->xhatminus_data);

    //控制向量初始化
    KF->u_data = (float *)user_malloc(sizeof(float) * xhat_size);
    memset(KF->u_data, 0, xhat_size);
    Matrix_Init(&KF->u, KF->xhat_size, 1, (float *)KF->u_data);

    //测量数据向量初始化
    KF->z_data = (float *)user_malloc(sizeof(float) * xhat_size);
    memset(KF->z_data, 0, xhat_size);
    Matrix_Init(&KF->z, KF->xhat_size, 1, (float *)KF->z_data);

    //协方差矩阵初始化
    KF->P_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    memset(KF->P_data, 0, xhat_size * xhat_size);
    Matrix_Init(&KF->P, KF->xhat_size, KF->xhat_size, (float *)KF->P_data);

    //先验估计协方差矩阵初始化
    KF->Pminus_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    memset(KF->Pminus_data, 0, xhat_size * xhat_size);
    Matrix_Init(&KF->Pminus, KF->xhat_size, KF->xhat_size, (float *)KF->Pminus_data);

    //状态转移矩阵初始化
    KF->A_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    KF->AT_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    memset(KF->A_data, 0, xhat_size * xhat_size);
    memset(KF->AT_data, 0, xhat_size * xhat_size);
    Matrix_Init(&KF->A, KF->xhat_size, KF->xhat_size, (float *)KF->A_data);
    Matrix_Init(&KF->AT, KF->xhat_size, KF->xhat_size, (float *)KF->AT_data);
    //Matrix_Transpose(&KF->A, &KF->AT);

    //控制矩阵初始化
    KF->B_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    memset(KF->B_data, 0, xhat_size * xhat_size);
    Matrix_Init(&KF->B, KF->xhat_size, KF->xhat_size, (float *)KF->B_data);

    //观测矩阵初始化
    KF->H_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    KF->HT_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    memset(KF->H_data, 0, xhat_size * xhat_size);
    memset(KF->HT_data, 0, xhat_size * xhat_size);
    Matrix_Init(&KF->H, KF->xhat_size, KF->xhat_size, (float *)KF->H_data);
    Matrix_Init(&KF->HT, KF->xhat_size, KF->xhat_size, (float *)KF->HT_data);
    //Matrix_Transpose(&KF->H, &KF->HT);

    //模型噪声矩阵初始化
    KF->Q_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    memset(KF->Q_data, 0, xhat_size * xhat_size);
    Matrix_Init(&KF->Q, KF->xhat_size, KF->xhat_size, (float *)KF->Q_data);

    //测量噪声矩阵初始化
    KF->R_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    memset(KF->R_data, 0, xhat_size * xhat_size);
    Matrix_Init(&KF->R, KF->xhat_size, KF->xhat_size, (float *)KF->R_data);

    //卡尔曼增益矩阵初始化
    KF->K_data = (float *)user_malloc(sizeof(float) * xhat_size * xhat_size);
    memset(KF->K_data, 0, xhat_size * xhat_size);
    Matrix_Init(&KF->K, KF->xhat_size, KF->xhat_size, (float *)KF->K_data);

    KF->temp_matrix_data = (float *)user_malloc(sizeof(float) * KF->xhat_size * KF->xhat_size);
    KF->temp_matrix_data1 = (float *)user_malloc(sizeof(float) * KF->xhat_size * KF->xhat_size);
    KF->temp_vector_data = (float *)user_malloc(sizeof(float) * KF->xhat_size);
    KF->temp_vector_data1 = (float *)user_malloc(sizeof(float) * KF->xhat_size);

    Matrix_Init(&KF->temp_matrix, KF->xhat_size, KF->xhat_size, (float *)KF->temp_matrix_data);
    Matrix_Init(&KF->temp_matrix1, KF->xhat_size, KF->xhat_size, (float *)KF->temp_matrix_data);
    Matrix_Init(&KF->temp_vector, KF->xhat_size, 1, (float *)KF->temp_vector_data);
    Matrix_Init(&KF->temp_vector1, KF->xhat_size, 1, (float *)KF->temp_vector_data1);
}

float *Kalman_Filter_Calculate(kalman_filter_t *KF, float *input)
{
    //更新测量数据
    for (uint8_t i = 0; i < KF->xhat_size; i++)
    {
        KF->Raw_Value[i] = input[i];
        KF->z_data[i] = input[i];
    }

    //1. xhat'(k)= A xhat(k-1) + B u
    Matrix_Multiply(&KF->A, &KF->xhat, &KF->temp_vector);
    Matrix_Multiply(&KF->B, &KF->u, &KF->temp_vector1);
    Matrix_Add(&KF->temp_vector, &KF->temp_vector1, &KF->xhatminus);

    //2. P'(k) = A P(k-1) AT + Q
    Matrix_Transpose(&KF->A, &KF->AT);
    Matrix_Multiply(&KF->A, &KF->P, &KF->Pminus);
    Matrix_Multiply(&KF->Pminus, &KF->AT, &KF->temp_matrix);
    Matrix_Add(&KF->temp_matrix, &KF->Q, &KF->Pminus);

    //3. K(k) = P'(k) HT / (H P'(k) HT + R)
    Matrix_Transpose(&KF->H, &KF->HT);
    Matrix_Multiply(&KF->H, &KF->Pminus, &KF->K);
    Matrix_Multiply(&KF->K, &KF->HT, &KF->temp_matrix);
    Matrix_Add(&KF->temp_matrix, &KF->R, &KF->K);
    Matrix_Inverse(&KF->K, &KF->P);
    Matrix_Multiply(&KF->Pminus, &KF->HT, &KF->temp_matrix);
    Matrix_Multiply(&KF->temp_matrix, &KF->P, &KF->K);

    //4. xhat(k) = xhat'(k) + K(k) (z(k) - H xhat'(k))
    Matrix_Multiply(&KF->H, &KF->xhatminus, &KF->temp_vector);
    Matrix_Subtract(&KF->z, &KF->temp_vector, &KF->xhat);
    Matrix_Multiply(&KF->K, &KF->xhat, &KF->temp_vector);
    Matrix_Add(&KF->xhatminus, &KF->temp_vector, &KF->xhat);

    //5. P(k) = (1-K(k)H)P'(k)
    //参考的资料这里计算的是 Q - K(k)H 我还没弄明白为什么
    //就暂且将括号乘开了
    // Matrix_Multiply(&KF->K, &KF->H, &KF->P);
    // Matrix_Subtract(&KF->Q, &KF->P, &KF->temp_matrix);
    // Matrix_Multiply(&KF->temp_matrix, &KF->Pminus, &KF->P);

    //5. P(k) = P'(k)-K(k)HP'(k)
    Matrix_Multiply(&KF->K, &KF->H, &KF->temp_matrix);
    Matrix_Multiply(&KF->temp_matrix, &KF->Pminus, &KF->temp_matrix1);
    Matrix_Subtract(&KF->Pminus, &KF->temp_matrix1, &KF->P);

    //滤波数据输出数据
    for (uint8_t i = 0; i < KF->xhat_size; i++)
        KF->Filtered_Value[i] = KF->xhat_data[i];

    return KF->Filtered_Value;
}

#endif //ARM_MATH_CM4
