/**
  ******************************************************************************
  * @file    kalman filter.c
  * @author  Hongxi Wong
  * @version V1.1.0
  * @date    2020/12/17
  * @brief   C implementation of kalman filter
  ******************************************************************************
  * @attention 
  * 该卡尔曼滤波器可以在传感器采样频率不同的情况下，动态调整矩阵H R和K的维数与数值。
  * This implementation of kalman filter can dynamically adjust dimension and  
  * value of matrix H R and K according to the measurement validity under any 
  * circumstance that the sampling rate of component sensors are different.
  * 
  * 因此矩阵H和R的初始化会与矩阵P A和Q有所不同。另外的，在初始化量测向量z时需要额外写
  * 入传感器量测所对应的状态与这个量测的方式，详情请见例程
  * Therefore, the initialization of matrix P, A, and Q is sometimes different 
  * from that of matrices H R. when initialization. Additionally, the corresponding 
  * state and the method of the measurement should be provided when initializing 
  * measurement vector z. For more details, please see the example. 
  * 
  * 若不需要动态调整量测向量z，可简单将结构体中的Use_Auto_Adjustment初始化为0，并像初
  * 始化矩阵P那样用常规方式初始化z H R即可。
  * If automatic adjustment is not required, assign zero to the UseAutoAdjustment 
  * and initialize z H R in the normal way as matrix P.
  * 
  * 要求量测向量z与控制向量u在传感器回调函数中更新。整数0意味着量测无效，即自上次卡尔曼
  * 滤波更新后无传感器数据更新。因此量测向量z与控制向量u会在卡尔曼滤波更新过程中被清零
  * MeasuredVector and ControlVector are required to be updated in the sensor 
  * callback function. Integer 0 in measurement vector z indicates the invalidity 
  * of current measurement, so MeasuredVector and ControlVector will be reset 
  * (to 0) during each update. 
  * 
  * 此外，矩阵P过度收敛后滤波器将难以再适应状态的缓慢变化，从而产生滤波估计偏差。该算法
  * 通过限制矩阵P最小值的方法，可有效抑制滤波器的过度收敛，详情请见例程。
  * Additionally, the excessive convergence of matrix P will make filter incapable
  * of adopting the slowly changing state. This implementation can effectively
  * suppress filter excessive convergence through boundary limiting for matrix P.
  * For more details, please see the example.
  * 
  * @example:
  * x = 
  *   |   height   |
  *   |  velocity  |
  *   |acceleration|
  * 
  * KalmanFilter_t Height_KF;
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
  *     // 设置最小方差
  *     static float state_min_variance[3] = {0.03, 0.005, 0.1};
  *     
  *     // 开启自动调整
  *     Height_KF.UseAutoAdjustment = 1;
  * 
  *     // 气压测得高度 GPS测得高度 加速度计测得z轴运动加速度
  *     static uint8_t measurement_reference[3] = {1, 1, 3}
  * 
  *     static float measurement_degree[3] = {1, 1, 1}     
  *     // 根据measurement_reference与measurement_degree生成H矩阵如下（在当前周期全部测量数据有效情况下）
  *       |1   0   0|
  *       |1   0   0|
  *       |0   0   1|
  * 
  *     static float mat_R_diagonal_elements = {30, 25, 35}
  *     //根据mat_R_diagonal_elements生成R矩阵如下（在当前周期全部测量数据有效情况下）
  *       |30   0   0|
  *       | 0  25   0|
  *       | 0   0  35|
  * 
  *     Kalman_Filter_Init(&Height_KF, 3, 0, 3);
  * 
  *     // 设置矩阵值
  *     memcpy(Height_KF.P_data, P_Init, sizeof(P_Init));
  *     memcpy(Height_KF.A_data, A_Init, sizeof(A_Init));
  *     memcpy(Height_KF.Q_data, Q_Init, sizeof(Q_Init));
  *     memcpy(Height_KF.MeasurementMap, measurement_reference, sizeof(measurement_reference));
  *     memcpy(Height_KF.MeasurementDegree, measurement_degree, sizeof(measurement_degree));
  *     memcpy(Height_KF.MatR_DiagonalElements, mat_R_diagonal_elements, sizeof(mat_R_diagonal_elements));
  *     memcpy(Height_KF.StateMinVariance, state_min_variance, sizeof(state_min_variance));
  * }
  * 
  * void INS_Task(void const *pvParameters)
  * {
  *     // 循环更新
  *     Kalman_Filter_Update(&Height_KF);
  *     vTaskDelay(ts);
  * }
  * 
  * // 测量数据更新应按照以下形式 即向MeasuredVector赋值
  * void Barometer_Read_Over(void)
  * {
  *     ......
  *     INS_KF.MeasuredVector[0] = baro_height;
  * }
  * void GPS_Read_Over(void)
  * {
  *     ......
  *     INS_KF.MeasuredVector[1] = GPS_height;
  * }
  * void Acc_Data_Process(void)
  * {
  *     ......
  *     INS_KF.MeasuredVector[2] = acc.z;
  * }
  ******************************************************************************
  */

#include "kalman filter.h"

uint16_t sizeof_float, sizeof_double;

static void H_K_R_Adjustment(KalmanFilter_t *KF);

void Kalman_Filter_Init(KalmanFilter_t *KF, uint8_t xhatSize, uint8_t uSize, uint8_t zSize)
{
    sizeof_float = sizeof(float);
    sizeof_double = sizeof(double);

    KF->xhatSize = xhatSize;
    KF->uSize = uSize;
    KF->zSize = zSize;

    KF->MeasurementValidNum = 0;

    // measurement flags
    KF->MeasurementMap = (uint8_t *)user_malloc(sizeof(uint8_t) * zSize);
    memset(KF->MeasurementMap, 0, sizeof(uint8_t) * zSize);
    KF->MeasurementDegree = (float *)user_malloc(sizeof_float * zSize);
    memset(KF->MeasurementDegree, 0, sizeof_float * zSize);
    KF->MatR_DiagonalElements = (float *)user_malloc(sizeof_float * zSize);
    memset(KF->MatR_DiagonalElements, 0, sizeof_float * zSize);
    KF->StateMinVariance = (float *)user_malloc(sizeof_float * xhatSize);
    memset(KF->StateMinVariance, 0, sizeof_float * xhatSize);
    KF->temp = (uint8_t *)user_malloc(sizeof(uint8_t) * zSize);
    memset(KF->temp, 0, sizeof(uint8_t) * zSize);

    // filter data
    KF->RawValue = (float *)user_malloc(sizeof_float * xhatSize);
    memset(KF->RawValue, 0, sizeof_float * xhatSize);
    KF->FilteredValue = (float *)user_malloc(sizeof_float * xhatSize);
    memset(KF->FilteredValue, 0, sizeof_float * xhatSize);
    KF->MeasuredVector = (float *)user_malloc(sizeof_float * zSize);
    memset(KF->MeasuredVector, 0, sizeof_float * zSize);
    KF->ControlVector = (float *)user_malloc(sizeof_float * uSize);
    memset(KF->ControlVector, 0, sizeof_float * uSize);

    // xhat x(k|k)
    KF->xhat_data = (float *)user_malloc(sizeof_float * xhatSize);
    memset(KF->xhat_data, 0, sizeof_float * xhatSize);
    Matrix_Init(&KF->xhat, KF->xhatSize, 1, (float *)KF->xhat_data);

    // xhatminus x(k|k-1)
    KF->xhatminus_data = (float *)user_malloc(sizeof_float * xhatSize);
    memset(KF->xhatminus_data, 0, sizeof_float * xhatSize);
    Matrix_Init(&KF->xhatminus, KF->xhatSize, 1, (float *)KF->xhatminus_data);

    if (uSize != 0)
    {
        // control vector u
        KF->u_data = (float *)user_malloc(sizeof_float * uSize);
        memset(KF->u_data, 0, sizeof_float * uSize);
        Matrix_Init(&KF->u, KF->uSize, 1, (float *)KF->u_data);
    }

    // measurement vector z and z_buf
    KF->z_data = (float *)user_malloc(sizeof_float * zSize);
    memset(KF->z_data, 0, sizeof_float * zSize);
    Matrix_Init(&KF->z, KF->zSize, 1, (float *)KF->z_data);

    // covariance matrix P(k|k)
    KF->P_data = (float *)user_malloc(sizeof_float * xhatSize * xhatSize);
    memset(KF->P_data, 0, sizeof_float * xhatSize * xhatSize);
    Matrix_Init(&KF->P, KF->xhatSize, KF->xhatSize, (float *)KF->P_data);

    //create covariance matrix P(k|k-1)
    KF->Pminus_data = (float *)user_malloc(sizeof_float * xhatSize * xhatSize);
    memset(KF->Pminus_data, 0, sizeof_float * xhatSize * xhatSize);
    Matrix_Init(&KF->Pminus, KF->xhatSize, KF->xhatSize, (float *)KF->Pminus_data);

    // state transition matrix A AT
    KF->A_data = (float *)user_malloc(sizeof_float * xhatSize * xhatSize);
    KF->AT_data = (float *)user_malloc(sizeof_float * xhatSize * xhatSize);
    memset(KF->A_data, 0, sizeof_float * xhatSize * xhatSize);
    memset(KF->AT_data, 0, sizeof_float * xhatSize * xhatSize);
    Matrix_Init(&KF->A, KF->xhatSize, KF->xhatSize, (float *)KF->A_data);
    Matrix_Init(&KF->AT, KF->xhatSize, KF->xhatSize, (float *)KF->AT_data);

    if (uSize != 0)
    {
        // control matrix B
        KF->B_data = (float *)user_malloc(sizeof_float * xhatSize * uSize);
        memset(KF->B_data, 0, sizeof_float * xhatSize * uSize);
        Matrix_Init(&KF->B, KF->xhatSize, KF->uSize, (float *)KF->B_data);
    }

    // measurement matrix H
    KF->H_data = (float *)user_malloc(sizeof_float * zSize * xhatSize);
    KF->HT_data = (float *)user_malloc(sizeof_float * xhatSize * zSize);
    memset(KF->H_data, 0, sizeof_float * zSize * xhatSize);
    memset(KF->HT_data, 0, sizeof_float * xhatSize * zSize);
    Matrix_Init(&KF->H, KF->zSize, KF->xhatSize, (float *)KF->H_data);
    Matrix_Init(&KF->HT, KF->xhatSize, KF->zSize, (float *)KF->HT_data);

    // process noise covariance matrix Q
    KF->Q_data = (float *)user_malloc(sizeof_float * xhatSize * xhatSize);
    memset(KF->Q_data, 0, sizeof_float * xhatSize * xhatSize);
    Matrix_Init(&KF->Q, KF->xhatSize, KF->xhatSize, (float *)KF->Q_data);

    // measurement noise covariance matrix R
    KF->R_data = (float *)user_malloc(sizeof_float * zSize * zSize);
    memset(KF->R_data, 0, sizeof_float * zSize * zSize);
    Matrix_Init(&KF->R, KF->zSize, KF->zSize, (float *)KF->R_data);

    // kalman gain K
    KF->K_data = (float *)user_malloc(sizeof_float * xhatSize * zSize);
    memset(KF->K_data, 0, sizeof_float * xhatSize * zSize);
    Matrix_Init(&KF->K, KF->xhatSize, KF->zSize, (float *)KF->K_data);

    KF->S_data = (float *)user_malloc(sizeof_float * KF->xhatSize * KF->xhatSize);
    KF->temp_matrix_data = (float *)user_malloc(sizeof_float * KF->xhatSize * KF->xhatSize);
    KF->temp_matrix_data1 = (float *)user_malloc(sizeof_float * KF->xhatSize * KF->xhatSize);
    KF->temp_vector_data = (float *)user_malloc(sizeof_float * KF->xhatSize);
    KF->temp_vector_data1 = (float *)user_malloc(sizeof_float * KF->xhatSize);
    Matrix_Init(&KF->S, KF->xhatSize, KF->xhatSize, (float *)KF->S_data);
    Matrix_Init(&KF->temp_matrix, KF->xhatSize, KF->xhatSize, (float *)KF->temp_matrix_data);
    Matrix_Init(&KF->temp_matrix1, KF->xhatSize, KF->xhatSize, (float *)KF->temp_matrix_data1);
    Matrix_Init(&KF->temp_vector, KF->xhatSize, 1, (float *)KF->temp_vector_data);
    Matrix_Init(&KF->temp_vector1, KF->xhatSize, 1, (float *)KF->temp_vector_data1);
}

float *Kalman_Filter_Update(KalmanFilter_t *KF)
{
    // 矩阵H K R根据量测情况自动调整
    // matrix H K R auto adjustment
    if (KF->UseAutoAdjustment != 0)
        H_K_R_Adjustment(KF);
    else
    {
        memcpy(KF->z_data, KF->MeasuredVector, sizeof_float * KF->zSize);
        memset(KF->MeasuredVector, 0, sizeof_float * KF->zSize);
    }

    memcpy(KF->u_data, KF->ControlVector, sizeof_float * KF->uSize);

    // 1. xhat'(k)= A·xhat(k-1) + B·u
    if (KF->uSize > 0)
    {
        Matrix_Multiply(&KF->A, &KF->xhat, &KF->temp_vector);
        KF->temp_vector1.numRows = KF->xhatSize;
        KF->temp_vector1.numCols = 1;
        Matrix_Multiply(&KF->B, &KF->u, &KF->temp_vector1);
        Matrix_Add(&KF->temp_vector, &KF->temp_vector1, &KF->xhatminus);
    }
    else
    {
        Matrix_Multiply(&KF->A, &KF->xhat, &KF->xhatminus);
    }

    // 预测更新
    // 2. P'(k) = A·P(k-1)·AT + Q
    Matrix_Transpose(&KF->A, &KF->AT);
    Matrix_Multiply(&KF->A, &KF->P, &KF->Pminus);
    KF->temp_matrix.numRows = KF->Pminus.numRows;
    KF->temp_matrix.numCols = KF->AT.numCols;
    Matrix_Multiply(&KF->Pminus, &KF->AT, &KF->temp_matrix); //temp_matrix = A P(k-1) AT
    Matrix_Add(&KF->temp_matrix, &KF->Q, &KF->Pminus);

    if (KF->MeasurementValidNum != 0 || KF->UseAutoAdjustment == 0)
    {
        // 量测更新
        // 3. K(k) = P'(k)·HT / (H·P'(k)·HT + R)
        Matrix_Transpose(&KF->H, &KF->HT); //z|x => x|z
        KF->temp_matrix.numRows = KF->H.numRows;
        KF->temp_matrix.numCols = KF->Pminus.numCols;
        Matrix_Multiply(&KF->H, &KF->Pminus, &KF->temp_matrix); //temp_matrix = H·P'(k)
        KF->temp_matrix1.numRows = KF->temp_matrix.numRows;
        KF->temp_matrix1.numCols = KF->HT.numCols;
        Matrix_Multiply(&KF->temp_matrix, &KF->HT, &KF->temp_matrix1); //temp_matrix1 = H·P'(k)·HT
        KF->S.numRows = KF->R.numRows;
        KF->S.numCols = KF->R.numCols;
        Matrix_Add(&KF->temp_matrix1, &KF->R, &KF->S); //S = H P'(k) HT + R
        Matrix_Inverse(&KF->S, &KF->temp_matrix1);     //temp_matrix1 = inv(H·P'(k)·HT + R)
        KF->temp_matrix.numRows = KF->Pminus.numRows;
        KF->temp_matrix.numCols = KF->HT.numCols;
        Matrix_Multiply(&KF->Pminus, &KF->HT, &KF->temp_matrix); //temp_matrix = P'(k)·HT
        Matrix_Multiply(&KF->temp_matrix, &KF->temp_matrix1, &KF->K);

        // 4. xhat(k) = xhat'(k) + K(k)·(z(k) - H·xhat'(k))
        KF->temp_vector.numRows = KF->H.numRows;
        KF->temp_vector.numCols = 1;
        Matrix_Multiply(&KF->H, &KF->xhatminus, &KF->temp_vector); //temp_vector = H xhat'(k)
        KF->temp_vector1.numRows = KF->z.numRows;
        KF->temp_vector1.numCols = 1;
        Matrix_Subtract(&KF->z, &KF->temp_vector, &KF->temp_vector1); //temp_vector1 = z(k) - H·xhat'(k)
        KF->temp_vector.numRows = KF->K.numRows;
        KF->temp_vector.numCols = 1;
        Matrix_Multiply(&KF->K, &KF->temp_vector1, &KF->temp_vector); //temp_vector = K(k)·(z(k) - H·xhat'(k))
        Matrix_Add(&KF->xhatminus, &KF->temp_vector, &KF->xhat);

        // 5. P(k) = (1-K(k)·H)·P'(k) ==> P(k) = P'(k)-K(k)·H·P'(k)
        KF->temp_matrix.numRows = KF->K.numRows;
        KF->temp_matrix.numCols = KF->H.numCols;
        KF->temp_matrix1.numRows = KF->temp_matrix.numRows;
        KF->temp_matrix1.numCols = KF->Pminus.numCols;
        Matrix_Multiply(&KF->K, &KF->H, &KF->temp_matrix);                 //temp_matrix = K(k)·H
        Matrix_Multiply(&KF->temp_matrix, &KF->Pminus, &KF->temp_matrix1); //temp_matrix1 = K(k)·H·P'(k)
        Matrix_Subtract(&KF->Pminus, &KF->temp_matrix1, &KF->P);
    }
    else
    {
        // 无有效量测 仅预测
        // xhat(k) = xhat'(k)
        // P(k) = P'(k)
        memcpy(KF->xhat_data, KF->xhatminus_data, sizeof_float * KF->xhatSize);
        memcpy(KF->P_data, KF->Pminus_data, sizeof_float * KF->xhatSize * KF->xhatSize);
    }

    // 避免滤波器过度收敛
    // suppress filter excessive convergence
    for (uint8_t i = 0; i < KF->xhatSize; i++)
    {
        if (KF->P_data[i * KF->xhatSize + i] < KF->StateMinVariance[i])
            KF->P_data[i * KF->xhatSize + i] = KF->StateMinVariance[i];
    }

    if (KF->UseAutoAdjustment != 0)
    {
        memset(KF->R_data, 0, sizeof_float * KF->zSize * KF->zSize);
        memset(KF->H_data, 0, sizeof_float * KF->xhatSize * KF->zSize);
    }

    memcpy(KF->FilteredValue, KF->xhat_data, sizeof_float * KF->xhatSize);

    return KF->FilteredValue;
}

static void H_K_R_Adjustment(KalmanFilter_t *KF)
{
    KF->MeasurementValidNum = 0;

    memcpy(KF->z_data, KF->MeasuredVector, sizeof_float * KF->zSize);
    memset(KF->MeasuredVector, 0, sizeof_float * KF->zSize);

    // 识别量测数据有效性并调整矩阵H R K
    // recognize measurement validity and adjust matrices H R K
    for (uint8_t i = 0; i < KF->zSize; i++)
    {
        if (KF->z_data[i] != 0)
        {
            // 重构向量z
            // rebuild vector z
            KF->z_data[KF->MeasurementValidNum] = KF->z_data[i];
            KF->temp[KF->MeasurementValidNum] = i;
            // 重构矩阵H
            // rebuild matrix H
            KF->H_data[KF->xhatSize * KF->MeasurementValidNum + KF->MeasurementMap[i] - 1] = KF->MeasurementDegree[i];
            KF->MeasurementValidNum++;
        }
    }
    for (uint8_t i = 0; i < KF->MeasurementValidNum; i++)
    {
        // 重构矩阵R
        // rebuild matrix R
        KF->R_data[i * KF->MeasurementValidNum + i] = KF->MatR_DiagonalElements[KF->temp[i]];
    }

    // 调整矩阵维数
    // adjust the dimensions of system matrices
    KF->H.numRows = KF->MeasurementValidNum;
    KF->H.numCols = KF->xhatSize;
    KF->HT.numRows = KF->xhatSize;
    KF->HT.numCols = KF->MeasurementValidNum;
    KF->R.numRows = KF->MeasurementValidNum;
    KF->R.numCols = KF->MeasurementValidNum;
    KF->K.numRows = KF->xhatSize;
    KF->K.numCols = KF->MeasurementValidNum;
}
