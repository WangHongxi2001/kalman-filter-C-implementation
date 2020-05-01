/**
******************************************************************************
  * @name    kalman_filter implementation
  * @author  Hongxi Wong <hongxiw0x7d1@foxmail.com>
  * @version V1.0.5
  * @date    2020/5/1
  * @brief   
******************************************************************************
  * @attention 
  * 该卡尔曼滤波器可以在传感器采样频率不同的情况下，动态调整矩阵H R和K的维数与数值。
  * This implementation of kalman filter can dynamically adjust dimension and  
  * value of matrix H R and K according to the measurement validity under any 
  * circumstance that the sampling rate of component sensors are different.
  * 
  * 因此矩阵H和R的初始化会与矩阵P A和Q有所不同。另外的，在初始化量测向量z时需要额外写
  * 入传感器量测所对应的状态与这个量测的方式，详情请见c文件内例程
  * Therefore, the initialization of matrix P, A, and Q is sometimes different 
  * from that of matrix H R. when initialization. Additionally, the corresponding 
  * state and the method of the measurement should be provided when initializing 
  * measurement vector z. For more details, please see the example in kalman filter.c 
  * 
  * 若不需要动态调整量测向量z，可简单将结构体中的Use_Auto_Adjustment初始化为0，并像初
  * 始化矩阵P那样用常规方式初始化z H R即可。
  * If automatic adjustment is not required, assign zero to the Use_Auto_Adjustment 
  * and initialize z H R in the normal way as matrix P.
  * 
  * 要求量测向量z与控制向量u在传感器回调函数中更新。整数0意味着量测无效，即自上次卡尔曼
  * 滤波更新后无传感器数据更新。因此量测向量z与控制向量u会在卡尔曼滤波更新过程中被清零
  * Measured_Vector and Control_Vector are required to be updated in the sensor 
  * callback function. Integer 0 in measurement vector z indicates the invalidity 
  * of current measurement, so Measured_Vector and Control_Vector will be reset 
  * (to 0) during each update. 
******************************************************************************
  */
