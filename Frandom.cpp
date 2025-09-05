#include <iostream>
#include "crt.h"  // 假设你有一个包含 Crt 类实现的 crt.h 头文件
#include "big.h"  // MIRACL 库的 big.h 头文件
#include "miracl.h" // MIRACL 主头文件
#include <ctime>
#include <sys/timeb.h>

using namespace std;

// 生成 n 个大素数作为模数，确保它们两两互素,p0,p1,p2...
bool generate_moduli(int n, int* bit_lengths, big* moduli, Big* moduli_new) {
    big gcd_val = mirvar(0);
    big tmp = mirvar(0);

    for (int i = 0; i < n; ++i) {
        int current_bit_length = bit_lengths[i];

        do {
            bigbits(current_bit_length, moduli[i]);  // 使用不同位数生成随机大数

            // 确保最高位为1，以防止随机数位数过小
            expb2(current_bit_length - 1, tmp);  // 生成 2^(current_bit_length-1)
            if (mr_compare(moduli[i], tmp) < 0) {  // 如果生成的随机数小于 2^(bit_length-1)
                add(moduli[i], tmp, moduli[i]);    // 将其增加到该范围内
            }

            if (!isprime(moduli[i])) continue;  // 确保该数为素数

            // 检查与之前的模数是否互素
            bool co_prime = true;
            for (int j = 0; j < i; ++j) {
                egcd(moduli[i], moduli[j], gcd_val);
                if (size(gcd_val) != 1) {  // 不互素
                    co_prime = false;
                    break;
                }
            }
            if (co_prime) break;  // 互素，跳出循环
        } while (true);
        moduli_new[i] = moduli[i];
        cout << "模数 " << i << ": " << moduli_new[i] << endl;
    }
    mirkill(gcd_val);
    return true;
}
int  main() {
    miracl *mip = mirsys(5000, 16);  // 初始化MIRACL库
    mip->IOBASE = 16; //十进制也可以

    // 定义参与者数量和阈值
    int num_participants = 6;  //实际参与者只有num_participants-1，多了一个p0，为了保证p0也和各个参与者的模数保持互素
    int threshold = 3;


    // 清理
    mirexit();

    return 0;
}


