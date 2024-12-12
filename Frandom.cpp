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

Big generate_ul(Big L){
    Big UL;
    do{
        UL=rand(L); // 随机生成 UL, 0 <= UL < L
    }while (UL==0);

    return UL;
}

// Lift：生成 S = s + p0 * UL
Big liftSecret(Big secret, Big p0, Big UL) {
    Big S;
    S= p0 * UL;          // S = p0 * UL
    S= secret + p0 * UL; // S = s + p0 * UL

    return S;
}
//密钥的份额分发
void distributesk(Big S, int num_moduli, Big shares[], Big moduli_new[]) {
    for (int i = 1; i < num_moduli; i++) {
        // 检查当前模数是否为零
        if (moduli_new[i] == 0) {
            cout << "错误：模数 moduli[" << i << "] 为零，无法进行分发" << endl;
            continue;
        }
        // 进行模运算
        shares[i] = S % moduli_new[i];  // shares[i] = S % pi
        // divide(S, moduli[i], shares[i]);  
    }
    cout << "分发的份额: " << endl;
    for (int i = 0; i < num_moduli; i++) {
        cout << "ski份额 " << i  << ": " << shares[i] << endl;
    }
}

// 随机值的份额分发
void distributeSecret(Big S, int rowi, int num_moduli, Big** Xij, Big moduli_new[]) {
    cout << "rowi:" <<rowi << endl;
    for (int i = 1; i < num_moduli; i++) {
        // 检查当前模数是否为零
        if (moduli_new[i] == 0) {
            cout << "错误：模数 moduli[" << i << "] 为零，无法进行分发" << endl;
            continue;
        }
        // 进行模运算
        Xij[rowi][i] = S % moduli_new[i];  // Xij[i] = S % pi
        // cout << "Xij:" << Xij[rowi][i] << endl; 
    }
}

void add_share(int num_participants, Big** Xij, Big moduli_new[], Big* Xi){
        //int ri[4]={0}; // 初始化使用
    for (int j = 1; j < num_participants; ++j) {
        for (int i = 1; i < num_participants; ++i) {
            Xi[j] += Xij[i][j]; // 累加份额值
        }
        // 取模
        Xi[j] %= moduli_new[j];
        cout << "Participant " << j << " secret share [X]_" << j << " = " << Xi[j] << endl;
    }
}

// 恢复秘密：从授权集合 A 中恢复秘密
Big recoverSecret(Big shares[], int* chosen_indexes, int num_chosen, Big moduli[]) {
    Big M=1;
    Big result=0;  // M = P（模数乘积）
    Big inv=1, temp;
    // 计算模数乘积 M
    for (int i = 0; i < num_chosen; i++) {
        M= M * moduli[chosen_indexes[i]];
        // multiply(M, moduli[chosen_indexes[i]], M);
    }
    big Mb=M.getbig();
    cout << "M:" << M << endl;
    //使用中国剩余定理恢复 S
    for (int i = 0; i < num_chosen; i++) {
        Big Mi;
        // big modulib=moduli[chosen_indexes[i]].getbig();
        Mi= M / moduli[chosen_indexes[i]];
        // divide(Mb, modulib, Mi);  // Mi = M / pi
        // cout << "Mi:" << Mi << endl;
        if (moduli[chosen_indexes[i]] == 0) {
            cout << "错误：模数为零" << endl;
            // return mirvar(0);
        }

        inv = inverse(Mi , moduli[chosen_indexes[i]]); // 计算逆元 inv = Mi^(-1) mod pi
        // cout << "inv:" << inv << endl;
        temp = 1;
        temp  = inv * Mi; // temp = Mi * inv
        temp = temp * shares[chosen_indexes[i]]; // temp = temp * si
        result = result + temp;  // result += temp
        // cout << "temp:" << temp << endl;
    }
    
    // result % M
    result=result % M;
    cout << "result_S:" << result << endl;

    return result;
}

// 从 S 中恢复秘密 s = S % p0
Big recoverOriginalSecret(Big S, Big p0) {
    Big secret;
    secret = S % p0;
    // divide(S, p0, secret);  // s = S % p0，余数保存在 secret 中
    return secret;
}

//生成随机值
void randomValue(int n, Big p0, Big* randV){
    for(int i = 1; i < n; i++){
        randV[i] =rand(p0);
        cout << "randV" << i << ":" << randV[i] << endl;
    }
}

void Frandom(int num_participants, Big p0, Big UL, Big** Xij, Big moduli_new[], Big* randV){
    for(int i = 1; i < num_participants; i++){
        // 生成 S = s + p0 * UL
        Big S = liftSecret(randV[i], p0, UL);
        cout << "S:"  << S << endl;

        distributeSecret(S, i, num_participants, Xij, moduli_new);

        cout << "分发的份额: " << endl;
        for (int j = 1; j < num_participants; j++) {

            cout << "X" <<i << "份额" << j  << ": " << Xij[i][j] << endl;
        
        }
    }
}

void Fdeg(int num_participants, Big p0, Big moduli_new[], Big* ki0, Big* ki1){
    //一个随机数在不同的范围中采样
    Big L1=5;
    Big L2=5000000;
    Big UL1 = generate_ul(L1);
    Big UL2 = generate_ul(L2);
    //初始化
    Big** kij0 = new Big*[num_participants];
    Big** kij1 = new Big*[num_participants];
    for (int i = 0; i < num_participants; ++i) {
        kij0[i] = new Big[num_participants];
        kij1[i] = new Big[num_participants];
    }
    cout << "UL1:" << UL1 << endl;
    cout << "UL2:" << UL2 << endl;

    Big* secret = new Big[num_participants];
    randomValue(num_participants, p0, secret);
    Frandom(num_participants, p0, UL1, kij0, moduli_new, secret);
    Frandom(num_participants, p0, UL2, kij1, moduli_new, secret);

    add_share(num_participants, kij0, moduli_new, ki0);
    add_share(num_participants, kij1, moduli_new, ki1);

    for (int j = 1; j < num_participants; ++j) {

        cout << "Participant " << j << " secret share [k]" << j << "0 = " << ki0[j] << endl;
        cout << "Participant " << j << " secret share [k]" << j << "1 = " << ki1[j] << endl;
    }
}

void Fmult(int num_participants, Big p0, Big* ri, Big* li, Big* ki0, Big* ki1, Big moduli_new[], Big* rlk_sub){
    Big rlk[num_participants];
    // Big rlk_sum;
    Big rlk_S,rlk_s;
    int num_chosen;
    int chosen_indexes[num_participants-1];
    for (int i = 1; i < num_participants; ++i) {

        rlk[i] = ri[i] * li[i];
        rlk[i] += ki1[i];
        rlk[i] %= moduli_new[i];
        // rlk_sum+=rlk[i];
    }
    num_chosen = num_participants-1;
    //选择固定的恢复秘密值的参与方（所有）
    for (int i = 0; i < num_participants; ++i) {
        chosen_indexes[i] = i+1;
        // cout << "chosen_indexes:" << chosen_indexes[i] << endl;
    }
    cout << "chosen_indexes:" << chosen_indexes << endl;

    //重构 r*l+k
    rlk_S = recoverSecret(rlk, chosen_indexes, num_chosen, moduli_new);
    // rlk_s = recoverOriginalSecret(rlk_S, p0);

    for (int i = 1; i < num_participants; ++i) {

        rlk_S %= moduli_new[i];
        rlk_sub[i] = ki0[i] - rlk_S;
        cout << "rlk_sub" << i <<":" << rlk_sub[i] << endl;

    }
}

void Fopen(int num_participants, Big p0, Big UL, Big* outi, Big moduli_new[], Big &zerout_S){
    Big* zero = new Big[num_participants];
    for (int i = 0; i < num_participants; ++i) {
        zero[i] = 0;
    }
    Big* zeroi = new Big[num_participants]; 
    Big zerouti[num_participants];
    int chosen_indexes[num_participants-1];
    int num_chosen;
    Big** zeroij = new Big*[num_participants];
    for (int i = 0; i < num_participants; ++i) {
        zeroij[i] = new Big[num_participants];
    }

    //生成0的份额--Frandom
    Frandom(num_participants, p0, UL, zeroij, moduli_new, zero);
    add_share(num_participants, zeroij, moduli_new, zeroi);

    //([0]i+[out]i) mod pi as the secret shares of 0+out
    for (int i = 1; i < num_participants; ++i) {

        zerouti[i] = zeroi[i]+ outi[i];
        zerouti[i] %= moduli_new[i];

    }

    for (int i = 0; i < num_participants; ++i) {
        chosen_indexes[i] = i+1;
        // cout << "chosen_indexes:" << chosen_indexes[i] << endl;
    }
    cout << "chosen_indexes:" << chosen_indexes << endl;

    num_chosen = num_participants-1;
    //重构 0+out
    zerout_S = recoverSecret(zerouti, chosen_indexes, num_chosen, moduli_new);
    // zerout_s = recoverOriginalSecret(zerout_S, p0);

    cout << "zerout_S:" << zerout_S << endl;
    // return zerout_S;
}

//计算 l*G、预签名的ECDSA两个部分
void preSign(int num_participants, Big p0, Big* randV_li, Big theta, Big* ri, Big* deltai, Big moduli_new[], Big &l, Big* sigmai0, Big* sigmai1){
    Big theta_inv[num_participants];
    for (int i = 1; i < num_participants; ++i) {
        l += randV_li[i];
        theta_inv[i] = inverse(theta, p0);
        sigmai0[i] = theta_inv[i] * ri[i];
        sigmai0[i] %= moduli_new[i];
        // sigmai1[i] = l_x * theta_inv[i] * deltai[i];
        sigmai1[i] = theta_inv[i] * deltai[i];
        sigmai1[i] %= moduli_new[i];
        // cout << "P" << i << "sigmai0:" << sigmai0[i] << endl;
        // cout << "P" << i << "sigmai1:" << sigmai1[i] << endl;
    }
    cout << "l:" << l << endl;
}

int  main() {
    miracl *mip = mirsys(5000, 16);  // 初始化MIRACL库
    mip->IOBASE = 16; //十进制也可以

    // 定义参与者数量和阈值
    int num_participants = 6;  //实际参与者只有num_participants-1，多了一个p0，为了保证p0也和各个参与者的模数保持互素
    int threshold = 3;

    Big L = 5;  // L 的范围
    // int bit_lengths[num_participants]={400,146,147,148,149,150,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171};
    int bit_lengths[num_participants]={400,158,159,160,161,162};

    // 分配模数的内存
    big* moduli = new big[num_participants];
    for (int i = 0; i < num_participants; ++i) {
        moduli[i] = mirvar(0);
    }
    Big* moduli_new = new Big[num_participants];
    // 生成模数,必须是互质的
    if (!generate_moduli(num_participants, bit_lengths, moduli, moduli_new)) {
        cout << "模数生成失败。" << endl;
        return 1;
    }

    Big p0=moduli[0];

    Big sk = "A56DFC67958970C946F36D4A30BA04F0178F2E1C4D9CD274" ;
    
    Big* randV_ri = new Big[num_participants];
    Big* randV_li = new Big[num_participants];

    // 定义动态二维数组 rij
    Big** rij = new Big*[num_participants];
    for (int i = 0; i < num_participants; ++i) {
        rij[i] = new Big[num_participants];
    }

    Big* ri = new Big[num_participants];

    Big** lij = new Big*[num_participants];
    for (int i = 0; i < num_participants; ++i) {
        lij[i] = new Big[num_participants];
    }

    Big* li = new Big[num_participants];

    Big UL = generate_ul(L);

    Big* ki0 = new Big[num_participants];
    Big* ki1 = new Big[num_participants];

    Big* thetai = new Big[num_participants];
    Big* deltai = new Big[num_participants];

    Big theta, l=0;

    Big* sigmai0 = new Big[num_participants];
    Big* sigmai1 = new Big[num_participants];

    //签名密钥sk份额生成
    Big S = liftSecret(sk, p0, UL);
    cout << "UL:" << UL << endl;
    cout << "p0:" << p0 << endl;
    cout << "S:"  << S << endl;
    Big ski[num_participants];
    distributesk(S, num_participants, ski, moduli_new);

    //各参与方需要生成一些随机数ri、li
    randomValue(num_participants, p0, randV_ri);
    randomValue(num_participants, p0, randV_li);

    //[r]i <-- Frandom(ri);  [l]i <-- Frandom(li)
    Frandom(num_participants, p0, UL, rij, moduli_new, randV_ri); 
    add_share(num_participants, rij, moduli_new, ri);
    Frandom(num_participants, p0, UL, lij, moduli_new, randV_li); 
    add_share(num_participants, lij, moduli_new, li);

    Fdeg(num_participants, p0, moduli_new, ki0, ki1);

    // [theta_i]=Fmult(ri*li), [delta_i]=Fmult(ri*ski)
    Fmult(num_participants, p0, ri, li, ki0, ki1, moduli_new, thetai);
    Fmult(num_participants, p0, ri, ski, ki0, ki1, moduli_new, deltai);

    Fopen(num_participants, p0, UL, thetai, moduli_new, theta);
    // cout << "theta:" << theta << endl;

    //l=sum_li, sigmai0=theta^-1 * ri, sigmai1=lx*theta^-1*deltai
    preSign(num_participants, p0, randV_li, theta, ri, deltai, moduli_new, l, sigmai0, sigmai1);
    cout << "l:" << l << endl;
    for (int i = 1; i < num_participants; ++i) {
        cout << "P" << i << "sigmai0:" << sigmai0[i] << endl;
        cout << "P" << i << "sigmai1:" << sigmai1[i] << endl; 
    }

    //计算l^-1的份额
    Big l_inv=inverse(l,p0);
    cout << "l_inv:" << l_inv << endl;
    for (int i = 1; i < num_participants; ++i) {
        Big ver = l_inv % moduli_new[i];
        cout << ver << endl;
    }
    Big l_S = liftSecret(l_inv, p0, UL);
    Big linv[num_participants];
    distributesk(S, num_participants, linv, moduli_new);


    // 释放动态分配的内存
    for (int i = 0; i < num_participants; ++i) {
        delete[] rij[i];
    }
    delete[] rij;
    delete[] ri;
    delete[] lij;
    delete[] li;
    delete[] ki0;
    delete[] ki1;
    delete[] thetai;
    delete[] deltai;
    delete[] randV_ri;
    delete[] randV_li;
    delete[] sigmai0;
    delete[] sigmai1;

    // 释放 moduli 内存
    delete[] moduli;

    // 清理
    mirexit();

    return 0;
}


