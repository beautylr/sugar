#include <iostream>

// 扩展欧几里得算法，返回 gcd 和逆元
int extendedGCD(int a, int b, int &x, int &y) {
    if (b == 0) {
        x = 1; 
        y = 0; 
        return a;
    }
    
    int x1, y1;
    int gcd = extendedGCD(b, a % b, x1, y1);
    
    // 更新 x 和 y
    x = y1;
    y = x1 - (a / b) * y1;
    
    return gcd;
}

// 计算 a mod m 的乘法逆元
int modInverse(int a, int m) {
    int x, y;
    int gcd = extendedGCD(a, m, x, y);
    
    // 如果 gcd 不为 1，则逆元不存在
    if (gcd != 1) {
        std::cout << "Inverse doesn't exist";
        return -1; // 或者其他错误值
    } else {
        // x 可能为负值，因此将其转换为正值
        return (x % m + m) % m;
    }
}

int main() {
    int x = 47; // 这里是你想要找逆元的数
    int y = 17; // 这里是模数

    int inverse = modInverse(x, y);
    
    if (inverse != -1) {
        std::cout << "The multiplicative inverse of " << x << " mod " << y << " is " << inverse << std::endl;
    }

    return 0;
}
