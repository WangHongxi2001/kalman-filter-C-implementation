/**
  ******************************************************************************
  * @file    kalman_filter.c
  * @author  Hongxi Wong
  * @version V1.0.3
  * @date    2020/4/27
  * @brief   
  ******************************************************************************
  * @attention 
  * This kalman filter tool can dynamically adjust dimension and value of matrix 
  * H R and K according to measurement validity when the measuring frequency of 
  * sensors are not the same. 
  * 
  * So there will be some differences between matrix P A and Q between matrix H R 
  * when initialization. Additionally, "how measurement relates to the state vector"
  * and "how much directly do sensors measure states" should be given when initializing.
  * Please see the example for the detailed information. 
  * 
  * If you don't need the automatic adjustment. You can just simply remove the
  * H_K_R_Adjustment function in Kalman_Filter_Update and initialize z H R in general 
  * way like matrix P.
  * 
  * It is required that update z_data and u_data in sensor callback function. And 
  * integer 0 in measurement vector z means the current measurement is invalid. So 
  * z_data and u_data will be zeroed after each update. 
  * 
  * @example:
  * xhat = 
  *   |   height   |
  *   |  velocity  |
  *   |acceleration|
  * 
  * kalman_filter_t Height_KF;
  * 
  * void INS_Task_Init(void)
  * {
  *     static float P_Init[9] =
  *     {
  *         10, 0, 0, 
  *         0, 30, 0, 
  *         0, 0, 10, 
  *     };
  *     static float A_Init[9] =
  *     {
  *         1, dt, 0.5*dt*dt, 
  *         0, 1, dt, 
  *         0, 0, 1, 
  *     };
  *     static float Q_Init[9] =
  *     {
  *         0.25*dt*dt*dt*dt, 0.5*dt*dt*dt, 0.5*dt*dt, 
  *         0.5*dt*dt*dt,        dt*dt,         dt, 
  *         0.5*dt*dt,              dt,         1, 
  *     };
  * 
  *     //baro for height  GPS for height  IMU for acc
  *     static uint8_t measurement_reference[3] = {1, 1, 3}
  *     //barometer measures height indirectly
  *     static float measurement_degree[3] = {0.8, 1, 1}     
  *     //according to measurement_reference and measurement_degree
  *     //matrix H should be like this:
  *       |0.8   0   0|
  *       |  1   0   0|
  *       |  0   0   1|
  * 
  *     static float mat_R_diagonal_elements = {30, 25, 35}
  *     //according to mat_R_diagonal_elements 
  *     //matrix R should be like this:
  *       |30   0   0|
  *       | 0  25   0|
  *       | 0   0  35|
  * 
  *     Kalman_Filter_Init(&Height_KF, 3, 0, 3);
  *     memcpy(Height_KF.P_data, P_Init, sizeof(P_Init));
  *     memcpy(Height_KF.A_data, A_Init, sizeof(A_Init));
  *     memcpy(Height_KF.Q_data, Q_Init, sizeof(Q_Init));
  *     memcpy(Height_KF.Measurement_Reference, measurement_reference, sizeof(measurement_reference));
  *     memcpy(Height_KF.Measurement_Degree, measurement_degree, sizeof(measurement_degree));
  *     memcpy(Height_KF.Mat_R_Diagonal_Elements, mat_R_diagonal_elements, sizeof(mat_R_diagonal_elements));
  * }
  * 
  * void INS_Task(void const *pvParameters)
  * {
  *     //infinite loop
  *     Kalman_Filter_Update(&Height_KF);
  *     vTaskDelay(ts);
  * }
  * 
  * void Barometer_Read_Over(void)
  * {
  *     ......
  *     INS_KF.z_data[0] = baro_height;
  * }
  * void GPS_Read_Over(void)
  * {
  *     ......
  *     INS_KF.z_data[1] = GPS_height;
  * }
  * void Acc_Data_Process(void)
  * {
  *     ......
  *     INS_KF.z_data[2] = acc.z;
  * }
  ******************************************************************************
  */

#include "kalman_filter.h"
#ifdef ARM_MATH_CM4
#include "arm_math.h"

uint16_t sizeof_float, sizeof_double;

static uint8_t H_K_R_Adjustment(kalman_filter_t *KF);

void Kalman_Filter_Init(kalman_filter_t *KF, uint8_t xhat_size, uint8_t u_size, uint8_t z_size)
{
    sizeof_float = sizeof(float);
    sizeof_double = sizeof(double);

    KF->xhat_size = xhat_size;

    //allocate space for measurement flags
    KF->Measurement_Reference = (uint8_t *)user_malloc(sizeof(uint8_t) * z_size);
    memset(KF->Measurement_Reference, 0, sizeof(uint8_t) * z_size);
    KF->Measurement_Degree = (float *)user_malloc(sizeof_float * z_size);
    memset(KF->Measurement_Degree, 0, sizeof_float * z_size);
    KF->Mat_R_Diagonal_Elements = (float *)user_malloc(sizeof_float * z_size);
    memset(KF->Mat_R_Diagonal_Elements, 0, sizeof_float * z_size);

    //allocate space for filter data
    KF->Raw_Value = (float *)user_malloc(sizeof_float * xhat_size);
    memset(KF->Raw_Value, 0, sizeof_float * xhat_size);
    KF->Filtered_Value = (float *)user_malloc(sizeof_float * xhat_size);
    memset(KF->Filtered_Value, 0, sizeof_float * xhat_size);

    //create xhat x(k|k)
    KF->xhat_data = (float *)user_malloc(sizeof_float * xhat_size);
    memset(KF->xhat_data, 0, sizeof_float * xhat_size);
    Matrix_Init(&KF->xhat, KF->xhat_size, 1, (float *)KF->xhat_data);

    //create xhatminus x(k|k-1)
    KF->xhatminus_data = (float *)user_malloc(sizeof_float * xhat_size);
    memset(KF->xhatminus_data, 0, sizeof_float * xhat_size);
    Matrix_Init(&KF->xhatminus, KF->xhat_size, 1, (float *)KF->xhatminus_data);

    if (u_size != 0)
    {
        //create control vector u
        KF->u_data = (float *)user_malloc(sizeof_float * u_size);
        memset(KF->u_data, 0, sizeof_float * u_size);
        Matrix_Init(&KF->u, KF->u_size, 1, (float *)KF->u_data);
    }

    //create measurement vector z and z_buf
    KF->z_data = (float *)user_malloc(sizeof_float * z_size);
    memset(KF->z_data, 0, sizeof_float * z_size);
    Matrix_Init(&KF->z, KF->z_size, 1, (float *)KF->z_data);
    KF->z_buf_data = (float *)user_malloc(sizeof_float * z_size);
    memset(KF->z_buf_data, 0, sizeof_float * z_size);
    Matrix_Init(&KF->z_buf, KF->z_size, 1, (float *)KF->z_data);

    //create covariance matrix P(k|k)
    KF->P_data = (float *)user_malloc(sizeof_float * xhat_size * xhat_size);
    memset(KF->P_data, 0, sizeof_float * xhat_size * xhat_size);
    Matrix_Init(&KF->P, KF->xhat_size, KF->xhat_size, (float *)KF->P_data);

    //create covariance matrix P(k|k-1)
    KF->Pminus_data = (float *)user_malloc(sizeof_float * xhat_size * xhat_size);
    memset(KF->Pminus_data, 0, sizeof_float * xhat_size * xhat_size);
    Matrix_Init(&KF->Pminus, KF->xhat_size, KF->xhat_size, (float *)KF->Pminus_data);

    //create state transition matrix A AT
    KF->A_data = (float *)user_malloc(sizeof_float * xhat_size * xhat_size);
    KF->AT_data = (float *)user_malloc(sizeof_float * xhat_size * xhat_size);
    memset(KF->A_data, 0, sizeof_float * xhat_size * xhat_size);
    memset(KF->AT_data, 0, sizeof_float * xhat_size * xhat_size);
    Matrix_Init(&KF->A, KF->xhat_size, KF->xhat_size, (float *)KF->A_data);
    Matrix_Init(&KF->AT, KF->xhat_size, KF->xhat_size, (float *)KF->AT_data);

    if (u_size != 0)
    {
        //create control matrix B
        KF->B_data = (float *)user_malloc(sizeof_float * xhat_size * u_size);
        memset(KF->B_data, 0, sizeof_float * xhat_size * u_size);
        Matrix_Init(&KF->B, KF->xhat_size, KF->u_size, (float *)KF->B_data);
    }

    //create measurement matrix H
    KF->H_data = (float *)user_malloc(sizeof_float * z_size * xhat_size);
    KF->HT_data = (float *)user_malloc(sizeof_float * xhat_size * z_size);
    memset(KF->H_data, 0, sizeof_float * z_size * xhat_size);
    memset(KF->HT_data, 0, sizeof_float * xhat_size * z_size);
    Matrix_Init(&KF->H, KF->z_size, KF->xhat_size, (float *)KF->H_data);
    Matrix_Init(&KF->HT, KF->xhat_size, KF->z_size, (float *)KF->HT_data);
    //Matrix_Transpose(&KF->H, &KF->HT);

    //create process noise covariance matrix Q
    KF->Q_data = (float *)user_malloc(sizeof_float * xhat_size * xhat_size);
    memset(KF->Q_data, 0, sizeof_float * xhat_size * xhat_size);
    Matrix_Init(&KF->Q, KF->xhat_size, KF->xhat_size, (float *)KF->Q_data);

    //create measurement noise covariance matrix R
    KF->R_data = (float *)user_malloc(sizeof_float * xhat_size * xhat_size);
    memset(KF->R_data, 0, sizeof_float * xhat_size * xhat_size);
    Matrix_Init(&KF->R, KF->xhat_size, KF->xhat_size, (float *)KF->R_data);

    //creat kalman gain K
    KF->K_data = (float *)user_malloc(sizeof_float * xhat_size * z_size);
    memset(KF->K_data, 0, sizeof_float * xhat_size * z_size);
    Matrix_Init(&KF->K, KF->xhat_size, KF->z_size, (float *)KF->K_data);

    KF->temp_matrix_data = (float *)user_malloc(sizeof_float * KF->xhat_size * KF->xhat_size);
    KF->temp_matrix_data1 = (float *)user_malloc(sizeof_float * KF->xhat_size * KF->xhat_size);
    KF->temp_vector_data = (float *)user_malloc(sizeof_float * KF->xhat_size);
    KF->temp_vector_data1 = (float *)user_malloc(sizeof_float * KF->xhat_size);
    Matrix_Init(&KF->temp_matrix, KF->xhat_size, KF->xhat_size, (float *)KF->temp_matrix_data);
    Matrix_Init(&KF->temp_matrix1, KF->xhat_size, KF->xhat_size, (float *)KF->temp_matrix_data1);
    Matrix_Init(&KF->temp_vector, KF->xhat_size, 1, (float *)KF->temp_vector_data);
    Matrix_Init(&KF->temp_vector1, KF->xhat_size, 1, (float *)KF->temp_vector_data1);
}

float *Kalman_Filter_Update(kalman_filter_t *KF)
{
    static uint8_t valid_num = -1;

    //1. xhat'(k)= A xhat(k-1) + B u
    if (KF->u_size > 0)
    {
        Matrix_Multiply(&KF->A, &KF->xhat, &KF->temp_vector);
        Matrix_Multiply(&KF->B, &KF->u, &KF->temp_vector1);
        Matrix_Add(&KF->temp_vector, &KF->temp_vector1, &KF->xhatminus);

        memset(KF->u_data, 0, sizeof_float * KF->u_size); //turn u_data to invalid
    }
    else
    {
        Matrix_Multiply(&KF->A, &KF->xhat, &KF->xhatminus);
    }

    //2. P'(k) = A P(k-1) AT + Q
    Matrix_Transpose(&KF->A, &KF->AT);
    Matrix_Multiply(&KF->A, &KF->P, &KF->Pminus);
    Matrix_Multiply(&KF->Pminus, &KF->AT, &KF->temp_matrix);
    Matrix_Add(&KF->temp_matrix, &KF->Q, &KF->Pminus);

    valid_num = H_K_R_Adjustment(KF);

    if (valid_num != 0)
    {
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
        Matrix_Subtract(&KF->z_buf, &KF->temp_vector, &KF->xhat);
        Matrix_Multiply(&KF->K, &KF->xhat, &KF->temp_vector);
        Matrix_Add(&KF->xhatminus, &KF->temp_vector, &KF->xhat);

        //5. P(k) = (1-K(k)H)P'(k) ==> P(k) = P'(k)-K(k)HP'(k)
        Matrix_Multiply(&KF->K, &KF->H, &KF->temp_matrix);
        Matrix_Multiply(&KF->temp_matrix, &KF->Pminus, &KF->temp_matrix1);
        Matrix_Subtract(&KF->Pminus, &KF->temp_matrix1, &KF->P);
    }

    memcpy(KF->Filtered_Value, KF->xhat_data, sizeof_float * KF->xhat_size);

    return KF->Filtered_Value;
}

static uint8_t H_K_R_Adjustment(kalman_filter_t *KF)
{
    static uint8_t valid_num = 0;
    valid_num = 0;

    //if we don't use z_buf here, the latest data may be skipped through
    //when it comes in through DMA or interrupt during the update process
    memcpy(KF->z_buf_data, KF->z_data, sizeof_float * KF->z_size);
    memset(KF->z_data, 0, sizeof_float * KF->z_size);

    //recognize measurement validity and adjust matrix H R K
    for (uint8_t i = 0; i < KF->xhat_size; i++)
    {
        if (KF->z_buf_data[i] != 0)
        {
            KF->z_buf_data[valid_num] = KF->z_buf_data[i];
            KF->H_data[KF->xhat_size * valid_num + KF->Measurement_Reference[i] - 1] = KF->Measurement_Degree[i];
            KF->R_data[KF->xhat_size * valid_num + valid_num] = KF->Mat_R_Diagonal_Elements[i];
            valid_num++;
        }
    }

    Matrix_Init(&KF->H, valid_num, KF->xhat_size, (float *)KF->H_data);
    Matrix_Init(&KF->HT, KF->xhat_size, valid_num, (float *)KF->HT_data);
    Matrix_Init(&KF->R, valid_num, valid_num, (float *)KF->R_data);
    Matrix_Init(&KF->K, KF->xhat_size, valid_num, (float *)KF->K_data);

    return valid_num;
}

#endif //ARM_MATH_CM4
