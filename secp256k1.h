#ifndef SECP256K1_H  
#define SECP256K1_H  
  
#include <miracl.h>  
  
// secp256k1曲线的参数  
extern const big secp256k1_p;  
extern const big secp256k1_a;  
extern const big secp256k1_b;  
extern const big secp256k1_n;  
extern const big secp256k1_gx;  
extern const big secp256k1_gy;  
  
// 初始化secp256k1曲线参数的函数声明  
void init_secp256k1_params();  
  
#endif // SECP256K1_H
