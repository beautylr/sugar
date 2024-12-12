#include <iostream>
#include <cstring>
#include <fstream>
#include "ecn.h"
#include "crt.h"  
#include "big.h"  
#include "miracl.h" 
#include <ctime>
#include <sys/timeb.h>
#include <time.h>
#include <chrono>
#include <string>
#include <flash.h>
#include <zzn.h> 
#include <openssl/hmac.h>  
#include <openssl/evp.h>  
#include <iomanip>
#include <sstream>
#include <set>
#include <cstdint>

using namespace std;
const static int POINTXY_BYTESIZE=100;
const static int SMALLVECTOR_BITSIZE=10;

/*p0=q, guarantee that q and all modulus pi are mutually prime*/
/*get current timestamp */
Big getCurTimestamp() {
    chrono::time_point<chrono::system_clock, chrono::milliseconds> tp =
            chrono::time_point_cast<chrono::milliseconds>(
                    chrono::system_clock::now()); //get current time point
    time_t timestamp = tp.time_since_epoch().count(); //compute time length until 1970-1-1,00:00
    return timestamp;
}

void strip(char *name)
{ /* strip off filename extension */
    int i;
    for (i=0;name[i]!='\0';i++)
    {
        if (name[i]!='.') continue;
        name[i]='\0';
        break;
    }
}

/*h0:{0,1}*->Z_q* ps:this hash function as basic */
static Big h0(const char *input, Big order) {
    char s[32];
    sha256 sh;
    shs256_init(&sh);
    while (*input!=0) shs256_process(&sh,*input++);
    shs256_hash(&sh,s);
    Big output=from_binary(32,s);
    // return output;
    return output % order;
}


/*h1:{0,1}* \plus G \plus G  \plus{0,1}* --> Z_q* */
Big h1(const char * join, ECn VPKm, ECn Xm, Big Tm1, Big q) {
    Big vpkmx,vpkmy,xmx,xmy;
    VPKm.getxy(vpkmx,vpkmy);
    Xm.getxy(xmx,xmy);
    char *vpkmx_char=new char[POINTXY_BYTESIZE],*vpkmy_char=new char[POINTXY_BYTESIZE],
            *xmx_char=new char[POINTXY_BYTESIZE],*xmy_char=new char[POINTXY_BYTESIZE],
            *tm1_char=new char[POINTXY_BYTESIZE];
    vpkmx_char<<vpkmx;
    vpkmy_char<<vpkmy;
    xmx_char<<xmx;
    xmy_char<<xmy;
    tm1_char<<Tm1;
    const string &istring=string(vpkmx_char)+string(vpkmy_char)+string(xmx_char)+
            string(xmy_char)+string(tm1_char)+string(join);
    const char *ichar=istring.c_str();
    return h0(ichar,q);
}

/* h2:{0,1}* \plus G \plus G \plus Z_q* \plus G \plus{0,1}* --> Z_q* */
Big h2(const char * mi, ECn Ai, ECn VPKi, Big sigma, ECn L, Big Ti, Big q) {
    Big aix,aiy,vpkix,vpkiy,lx,ly;
    Ai.getxy(aix,aiy);
    VPKi.getxy(vpkix,vpkiy);
    L.getxy(lx,ly);
    char *aix_char=new char[POINTXY_BYTESIZE],*aiy_char=new char[POINTXY_BYTESIZE],
            *vpkix_char=new char[POINTXY_BYTESIZE],*vpkiy_char=new char[POINTXY_BYTESIZE],
            *lx_char=new char[POINTXY_BYTESIZE],*ly_char=new char[POINTXY_BYTESIZE],
            *sigma_char=new char[POINTXY_BYTESIZE],*ti_char=new char[POINTXY_BYTESIZE];
    aix_char<<aix;
    aiy_char<<aiy;
    vpkix_char<<vpkix;
    vpkiy_char<<vpkiy;
    lx_char<<lx;
    ly_char<<ly;
    sigma_char<<sigma;
    ti_char<<Ti;
    const string &istring=string(aix_char)+string(aiy_char)+string(vpkix_char)+string(vpkiy_char)+
            string(lx_char)+string(ly_char)+string(ti_char)+string(sigma_char)+string(mi);
    const char *ichar=istring.c_str();
    delete[] aix_char;
    delete[] aiy_char;
    delete[] vpkix_char;
    delete[] vpkiy_char;
    delete[] lx_char;
    delete[] ly_char;
    delete[] sigma_char;
    delete[] ti_char;
    return h0(ichar,q);
}

/* h1:{0,1}* \plus G \plus{0,1}* --> Z_q* */
Big h_signf(char * encm, ECn pkfj, Big Tm2, Big q) {
    Big pkfjx,pkfjy;
    pkfj.getxy(pkfjx,pkfjy);
    char *pkfjx_char=new char[POINTXY_BYTESIZE],*pkfjy_char=new char[POINTXY_BYTESIZE],
            *tm2_char=new char[POINTXY_BYTESIZE];
    pkfjx_char<<pkfjx;
    pkfjy_char<<pkfjy;
    tm2_char<<Tm2;
    const string &istring=string(encm)+string(pkfjx_char)+string(pkfjy_char)+string(tm2_char);
    const char *ichar=istring.c_str();
    return h0(ichar,q);
}

/*get message hash */
Big getmsghash(const char* message, Big q) {
    Big h = h0(message,q);
    return h;
}

/*Since weight generation relies on continuous iterative optimization of the algorithm and can be generated in advance based on a 
priori data, only the very simple code to get the values is given here for convenience so that it can be calculated later on */
void weight_select(int count, int bit_lengths[]){
    int temp = (count-1) / 2;
    int min = 220 - temp -2; int max = 220 + temp +2;
    bit_lengths[0] = 240; int i=1;
    if (count > (max - min + 1)) cout << "Unable to generate unique random numbers, numbers in scope insufficient." << endl;
    set<int> unique_numbers;
    srand(static_cast<unsigned int>(time(nullptr))); // initialize random seed 
    while (unique_numbers.size() < count-1) {
        int num = rand() % (max - min + 1) + min; // generate random numbers within a specified range 生成指定范围内的随机数
        unique_numbers.insert(num);
    }
    // output unique random value
    for (int num : unique_numbers) {
        bit_lengths[i++] = num; 
        // cout << "weight size:" << num <<endl;
    }
}

/* generate n big prime as modulus, make sure they are both prime to each other,p0(q),p1,p2... */
bool generate_moduli(int n, int* bit_lengths, big* moduli,Big q) {
    big gcd_val = mirvar(0);
    big tmp = mirvar(0);

    for (int i = 0; i < n; ++i) {
        int current_bit_length = bit_lengths[i];
        do {
            bigbits(current_bit_length, moduli[i]);  //use different bits to generate random big numbers
            //ensure that the highest bit is 1 to preventing the random number of bits too small
            expb2(current_bit_length - 1, tmp);  //generate 2^(current_bit_length-1)
            if (mr_compare(moduli[i], tmp) < 0) {  //if the random number genration less than 2^(bit_length-1)
                add(moduli[i], tmp, moduli[i]);    //add it to this scope
            }
            if (!isprime(moduli[i])) continue;  //ensure this number is prime
            //check whether the generated modulus is commutative with q 
            egcd(moduli[i], q.getbig(), gcd_val);
            if (size(gcd_val) != 1) {  
                continue;
            }
            //check whether the generated modulus is commutative with previous modulus
            bool co_prime = true;
            for (int j = 0; j < i; ++j) {
                egcd(moduli[i], moduli[j], gcd_val);
                if (size(gcd_val) != 1) {  
                    co_prime = false;
                    break;
                }
            }
            if (co_prime) break; 
        } while (true);
        // cout << "Moduli " << i + 1 << ": " << moduli[i] << endl;
    }
    mirkill(gcd_val);
    return true;
}

Big generate_ul(Big L){
    Big UL;
    do{
        UL=rand(L); //randomly generate UL, 0 <= UL < L
    }while (UL==0);
    return UL;
}

/*generate random number*/
void randomValue(int n, Big p0, Big* randV){
    for(int i = 0; i < n; i++) randV[i] =rand(p0);
}

/*Lift: compute S = s + p0 * UL*/ 
Big liftSecret(Big secret, Big p0, Big UL) {
    Big S;
    S= p0 * UL;         
    S= secret + p0 * UL; 
    return S;
}

/* distributing shares of secret key sk */
void distributesk(Big S, int num_moduli, Big shares[], Big moduli_new[]) {
    for (int i = 0; i < num_moduli; i++) {
        //check whether current modulus is zero 
        if (moduli_new[i] == 0) {
            cout << "Error：moduli[" << i << "] is zero, unable to distribute secret shares " << endl;
            continue;
        }
        // modulo operation
        shares[i] = S % moduli_new[i];  // shares[i] = S % pi
        // divide(S, moduli[i], shares[i]);  
    }
}
/* distributing shares of random numbers  */
void distributeRandomV(Big S, int rowi, int num_moduli, Big** Xij, Big moduli_new[]) {
    // cout << "rowi:" <<rowi << endl;
    for (int i = 0; i < num_moduli; i++) {
        if (moduli_new[i] == 0) {
            cout << "Error：moduli[" << i << "] is zero, unable to distribute secret shares" << endl;
            continue;
        }
        Xij[rowi][i] = S % moduli_new[i];  // Xij[i] = S % pi
        // cout << "Xij:" << Xij[rowi][i] << endl; 
    }
}

/* compute P_A, P_A_Ba*/
void compute_PA_PABa(Big* moduli_n, int size, int maxt, int mint, Big& P_A, Big& PA_Ba) {
    int size_new = size-1;
    Big* moduli_new = new Big[size_new];
    for (int i = 0; i < size_new; ++i) moduli_new[i] = moduli_n[i+1];
    // sorting arrays
    for (int i = 0; i < size_new; ++i) {
        for (int j = i + 1; j < size_new; ++j) {
            if (mr_compare(moduli_new[i].getbig(), moduli_new[j].getbig()) < 0) {
                Big temp = moduli_new[i];
                moduli_new[i] = moduli_new[j];
                moduli_new[j] = temp;
            }
        }
    }
    P_A = 1;
    PA_Ba = 1;
    for (int i = 0; i < maxt && i < size_new; ++i)  P_A *= moduli_new[i];
    for (int i = size_new - 1; i >= size_new - mint && i >= 0; --i)  PA_Ba *= moduli_new[i];
}

void generate_random_in_range(big lower_bound, big upper_bound, big random_result) {
    big range = mirvar(0);  // storage range
    big random_offset = mirvar(0);  
    // compute range: range = upper_bound - lower_bound
    subtract(upper_bound, lower_bound, range);
    // generate random number within range
    bigrand(range, random_offset);
    add(lower_bound, random_offset, random_result);
    mirkill(range);
    mirkill(random_offset);
}

/* reconstruct secrets from authorization set A */
Big recoverSecret(Big shares[], int* chosen_indexes, int num_chosen, Big moduli[]) {
    Big M=1; // M = P(modulus product)
    Big result;  
    Big inv=1, temp;
    for (int i = 0; i < num_chosen; i++) {
        M= M * moduli[chosen_indexes[i]];
        // multiply(M, moduli[chosen_indexes[i]], M);
    }
    // cout << "M:" << M << endl;
    // recovering S by using CRT
    for (int i = 0; i < num_chosen; i++) {
        Big Mi;
        Mi= M / moduli[chosen_indexes[i]];
        // divide(Mb, modulib, Mi);  // Mi = M / pi
        if (moduli[chosen_indexes[i]] == 0) cout << "Error: moduli is zero" << endl;
        inv = inverse(Mi , moduli[chosen_indexes[i]]); // compute inverse element:inv = Mi^(-1) mod pi
        temp = 1;
        temp  = inv * Mi; 
        temp = temp * shares[chosen_indexes[i]]; 
        result = result + temp;  
    }
    // result % M
    result %= M;
    // cout << "result_S:" << result << endl;
    return result;
}

/* recovering secret s from S: s = S % p0 */
Big recoverOriginalSecret(Big S, Big p0) {
    Big secret;
    secret = S % p0;
    return secret;
}

void add_share(int num_participants, Big** Xij, Big moduli_new[], Big* Xi){
    for (int j = 0; j < num_participants; ++j) {
        for (int i = 0; i < num_participants; ++i) {
            Xi[j] += Xij[i][j]; // cumulative share values
        }
        Xi[j] %= moduli_new[j];
        // cout << "Participant " << j << " secret share [X]_" << j << " = " << Xi[j] << endl;
    }
}

void Frandom(int num_participants, Big p0, Big UL, Big** Xij, Big moduli_new[], Big* randV){
    for(int i = 0; i < num_participants; i++){
        Big S = liftSecret(randV[i], p0, UL);  // compute S = s + p0 * UL
        // cout << "S:"  << S << endl;
        distributeRandomV(S, i, num_participants, Xij, moduli_new);
    }
}

void Fdeg(int num_participants, Big p0, Big moduli_new[], Big* ki0, Big* ki1, int threshold, int namda, Big UL1){
    //sampling in different range of random value k 
    Big L1=1;
    Big L2=1;
    L1 <<= (threshold+namda);
    L2 <<= (2*threshold+5*namda);
    UL1 = generate_ul(L1);
    Big UL2 = generate_ul(L2);
    Big** kij0 = new Big*[num_participants];
    Big** kij1 = new Big*[num_participants];
    for (int i = 0; i < num_participants; ++i) {
        kij0[i] = new Big[num_participants];
        kij1[i] = new Big[num_participants];
    }
    Big* randomk = new Big[num_participants];
    // randomValue(num_participants, p0, randomk); 
    Big onlyk = rand(p0);
    for (int i = 0; i < num_participants; ++i) {
        randomk[i] = onlyk;
    }
    Frandom(num_participants, p0, UL1, kij0, moduli_new, randomk);
    Frandom(num_participants, p0, UL2, kij1, moduli_new, randomk);
    add_share(num_participants, kij0, moduli_new, ki0);
    add_share(num_participants, kij1, moduli_new, ki1);
}

void Fmult(int num_participants, Big* Xi, Big* Yi, Big moduli_new[], Big* pointopi){
    for(int i=0;i<num_participants;i++){
        pointopi[i] = Xi[i] * Yi[i];
        pointopi[i] %= moduli_new[i];
    }
}

void Fopen(int num_participants, Big p0, Big UL, Big* zeroi, Big* outi, Big moduli_new[], Big &zerout_s,
    int* chosen_indexes,int num_chosen){
    Big zerouti[num_participants];
    //([0]i+[out]i) mod pi as the secret shares of 0+out
    for (int i = 0; i < num_chosen; ++i) {
        zerouti[chosen_indexes[i]] = zeroi[chosen_indexes[i]]+ outi[chosen_indexes[i]];
        zerouti[chosen_indexes[i]] %= moduli_new[chosen_indexes[i]];
        // cout << "zerouti[i]" << zerouti[chosen_indexes[i]] << endl;
    }
    // reconstruct 0+out
    Big zerout_S = recoverSecret(zerouti, chosen_indexes, num_chosen, moduli_new);
    zerout_s = recoverOriginalSecret(zerout_S, p0);
    // cout << "zerout_s:" << zerout_s << endl;
}

/* compute sigmai0 and sigmai1 of ECDSA in presign phase  */
void preSign(int num_participants, Big q, Big theta, Big* ri, Big* deltai, Big moduli_new[], Big &l_x, Big* sigmai0, Big* sigmai1){
    Big theta_inv;
    theta_inv = inverse(theta, q); 
    for (int i = 0; i < num_participants; ++i) {
        sigmai0[i] = theta_inv;
        sigmai0[i] *= ri[i];
        // sigmai0[i] %= moduli_new[i];
        sigmai1[i] = l_x;
        sigmai1[i] *= theta_inv;
        sigmai1[i] *= deltai[i];
        // sigmai1[i] %= moduli_new[i];
    }
}

/* imitate randomly participates  */
void choose_party(int threshold, int* chosen_indexes){
    bool selected[threshold] = { false };  // mark the selected participants
    // random number generator
    irand(time(0));
    int count = 0; // randomly selected threshold participants index
    while (count < threshold) {
        big random_value = mirvar(0);
        bigrand(mirvar(threshold), random_value); // generate random value within range [0, threshold)
        int index = size(random_value);  // translating random Big to int
        if (!selected[index ]) {         // check whether selected
            selected[index] = true;     // mark it selected
            chosen_indexes[count] = index; // storage selected index
            count++;
        }
    }
}

void sign_ecdsa(Big k, Big q,Big h, Big d, Big r, Big &s){
    Big v;
    Big k_inv = inverse(k,q);
    s = ((h+d*r)*k_inv)%q;    
}

void gen_sigma(const char * mi, int num_participants, int threshold, int namda, Big &q, Big &UL, Big &L_x,Big sigmai0[], Big sigmai1[], ECn &G, ECn &PK, 
    Big moduli_new[], Big *zeroi,ECn &Lg,int maxt,int mint,int* chosen_indexes){
    bool selected[num_participants] = { false };    
    int bit_lengths[num_participants];
    weight_select(num_participants,bit_lengths);

    ifstream common("common.ecs");    /* construct file I/O streams */
    ifstream private_key("private.ecs");
    ifstream public_key("public.ecs");
    const char* message;
    char ifname[50],ofname[50];
    ECn Pub, R;
    Big a,b,p,x,y,h,r,s,secret,k,u1,u2,v;
    long seed;
    int bits,ep;

    // cout << "Enter 9 digit random number seed  = ";
    // cin >> seed;
    irand(369266);

// get common data 
    common >> bits;
    common >> p >> a >> b >> q >> x >> y;
    ecurve(a,b,p,MR_PROJECTIVE);
    G=ECn(x,y);

// get public key of signer 
    public_key >> ep >> x;
    PK=ECn(x,ep);         // decompress
    // PK=ECn(xx,yy);

// get private key of recipient 
    private_key >> secret;
    // cout << "secret:" << secret << "\n";

// allocate modulus memory
    big* moduli = new big[num_participants];
    for (int i = 0; i < num_participants; ++i) moduli[i] = mirvar(0);

// generate modulus
    if (!generate_moduli(num_participants, bit_lengths, moduli, q)) cout << "Modulus generation failure." << endl;
    for (int i = 0; i < num_participants; ++i) moduli_new[i]=moduli[i];

    Big P_A, PA_Ba;
    Big L=1; 
    compute_PA_PABa(moduli_new, num_participants, maxt, mint, P_A, PA_Ba);
    L <<= (threshold+namda);
    UL = generate_ul(L);

    Big S = liftSecret(secret, q, UL); // generate S = s + p0 * UL

    // distribute shares of secret
    Big ski[num_participants];
    distributesk(S, num_participants, ski, moduli_new);

    Big* randV_ri = new Big[num_participants];
    Big** rij = new Big*[num_participants];
    for (int i = 0; i < num_participants; ++i) {
        rij[i] = new Big[num_participants];
    }
    Big* ri = new Big[num_participants];
    Big* randV_li = new Big[num_participants];
    Big** lij = new Big*[num_participants];
    for (int i = 0; i < num_participants; ++i) {
        lij[i] = new Big[num_participants];
    }
    Big* li = new Big[num_participants];
    Big l_sum = 0;

    Big* ki0 = new Big[num_participants]; Big* ki1 = new Big[num_participants];
    Big* thetai = new Big[num_participants]; Big* deltai = new Big[num_participants];

// generate ri,li
    randomValue(num_participants, q, randV_ri);
    randomValue(num_participants, q, randV_li);
    for(int i=0;i<num_participants;i++) l_sum += randV_li[i];

// [r]i <-- Frandom(ri);  [l]i <-- Frandom(li)
    Frandom(num_participants, q, UL, rij, moduli_new, randV_ri); 
    add_share(num_participants, rij, moduli_new, ri);
    Frandom(num_participants, q, UL, lij, moduli_new, randV_li); 
    add_share(num_participants, lij, moduli_new, li);

    Fdeg(num_participants, q, moduli_new, ki0, ki1, threshold,namda,UL);

// [theta_i]=Fmult(ri*li), [delta_i]=Fmult(ri*ski)
    Fmult(num_participants, ri, li, moduli_new, thetai);
    Fmult(num_participants, ri, ski, moduli_new, deltai);
    for(int i=0;i<num_participants;i++){
        // cout << "thetai:" << thetai[i] << endl;
        // cout << "deltai:" << deltai[i] << endl;
    }
// generate shares of 0--Frandom
    Big Lzero=1;
    Lzero <<= (threshold+3*namda);
    Big ULzero = generate_ul(Lzero);
    for(int i =0;i<num_participants;i++){ 
        Big zero_S = liftSecret((Big)0, q, ULzero);
        zeroi[i]=zero_S % moduli_new[i];
        // cout << "0 shares：" << zeroi[i] << endl;
    }

// select threshold shares of authorized set
    choose_party(threshold, chosen_indexes);
// simulate multiple randomized participates
    for (int i = 0; i < threshold; ++i) cout << chosen_indexes[i] << " ";
    cout << endl;     
// compute theta
    Big theta;
    Fopen(num_participants,q,UL,zeroi,thetai,moduli_new,theta,chosen_indexes,threshold);
    cout << "theta:" << theta << endl;

// calculate L_x - this can be done off-line,and hence amortized to almost nothing
    Lg = l_sum * G;       
    Lg.get(L_x); //compute x abscissa of point l_G
    L_x%=q; 

    preSign(num_participants, q, theta, ri, deltai, moduli_new, L_x, sigmai0, sigmai1);    
}

/* generate secret keys of vehicle */
struct gen_sk_pk {
    Big sk;
    ECn pk;
};
gen_sk_pk getValues(ECn G, Big q) {
    Big vski = rand(q);
    ECn VPKi = vski * G;
    return {vski, VPKi};
}

/* vehicle group formation */
void req_vehgroup(Big q, ECn G, Big vskm, ECn &VPKm, const char* join, ECn &Xm, Big &Tm1, Big &pdm){
    Big xm = rand(q);
    Tm1 = getCurTimestamp();
    Xm = xm * G;
    Big hm = h1(join,VPKm,Xm,Tm1,q);
    pdm = (xm + vskm* hm) % q;
}

void aes_enc_cbc(char* ichar, Big Wf_x, char* ciphertext, int& ciphertext_len){
    aes a;
    char iv[32]={0};   // initiation vector
    char key[32];
    to_binary(Wf_x,32,key);
    aes_init(&a,MR_CBC,32,key,iv);
    int len = strlen(ichar); // compute the padded plaintext length 
    int padded_len = len + (16 - (len % 16)); // pad to 填充到16的倍数
    char* padded_plaintext = new char[padded_len + 1];
    // copy plaintext and padding it 
    memcpy(padded_plaintext, ichar, len);
    memset(padded_plaintext + len, 16 - (len % 16), padded_len - len); // PKCS7 padding
    padded_plaintext[padded_len] = '\0'; // null terminator for safety
    for (int i = 0; i < padded_len; i += 16) {
        aes_encrypt(&a, padded_plaintext + i); // encrypt each block
    }
    memcpy(ciphertext, padded_plaintext, padded_len); //copy encrypt result
    ciphertext_len = padded_len;

    delete[] padded_plaintext;
}

bool aes_dec_cbc(char* ciphertext, int ciphertext_len, Big Wf_x, char* plaintext) {
    aes a;
    char iv[32] = {0}; 
    char key[32];
    to_binary(Wf_x, 32, key); // translate key to char arrays
    aes_init(&a, MR_CBC, 32, key, iv);
    for (int i = 0; i < ciphertext_len; i += 16) {
        aes_decrypt(&a, const_cast<char*>(ciphertext) + i); // decrypt each block
    }
    int padding = ciphertext[ciphertext_len - 1];
    memcpy(plaintext, ciphertext, ciphertext_len - padding);
    plaintext[ciphertext_len - padding] = '\0'; // null terminator for safety
    return 0;
}

/* fog node verifies the join group request of vehicle and distributes parameters */
void Fog_resp(Big q, ECn G,Big randvalue_k, ECn k_R, ECn VPKm, const char* join, ECn Xm, Big Tm1, Big pdm, Big fskj, ECn FPKj,
    char* encm,int &ciplen,Big &Tm2,Big &Signf,Big &gskfj){
    Big hm = h1(join,VPKm,Xm,Tm1,q);
    ECn left = pdm * G;
    ECn right;
    right = mul(Big(1),Xm,hm,VPKm);
    Big sm = rand(q),pm = rand(q);
    gskfj = rand(q);
    // cout << "gskfj:" << gskfj << endl;
    char sm_char[POINTXY_BYTESIZE],pm_char[POINTXY_BYTESIZE],
            gskfj_char[POINTXY_BYTESIZE];
    sm_char<<sm;
    pm_char<<pm;
    gskfj_char<<gskfj;
    string istring=string(sm_char)+string(pm_char)+string(gskfj_char);
    char plaintext[200];
    strncpy(plaintext, istring.c_str(), sizeof(plaintext) - 1); // Copy string
    plaintext[sizeof(plaintext) - 1] = '\0';

    if (left == right){
        cout << "Message M1 verification, fog agree request" << endl;
        ECn Wf = fskj * VPKm;
        Big Wf_x, k_r;
        Wf.get(Wf_x); 
        aes_enc_cbc(plaintext,Wf_x,encm,ciplen);
        Tm2 = getCurTimestamp();
        k_R.get(k_r);
        k_r%=q;
        Big hsignf = h_signf(encm, FPKj, Tm2, q);
        sign_ecdsa(randvalue_k, q, hsignf, fskj, k_r, Signf);
    } 
    else cout << "Message M1 not verification, fog not agree" << endl;
}


bool verify(Big r, Big s, Big q, Big h, ECn Pub, ECn G){
    Big u1,u2,v;
    ECn RR;
    if (r>=q || s>=q) cout << "Signature NOT verified" << endl;
    Big s_inv=inverse(s,q);
    u1=(h*s_inv)%q;
    u2=(r*s_inv)%q;
    RR=mul(u2,Pub,u1,G);
    // G=mul(u2,Pub,u1,G);
    RR.get(v);
    v%=q;

    if (v==r){
        cout << "Signature verification for message M2" << endl; return 1;}
    else{cout << "Signature not verified for message M2"<< endl; return 0;}
}

string to_hex(const char* data, size_t len) {
    std::ostringstream oss;
    for (size_t i = 0; i < len; ++i) {
        oss << std::hex << std::setw(2) << std::setfill('0') << (int)(unsigned char)data[i];
    }
    return oss.str();
}

const char* joint_message(ECn VPKi, const char *mi, Big Ti2, size_t &total_message_len){
    Big VPKix,VPKiy;
    VPKi.getxy(VPKix,VPKiy);
    char *VPKix_char=new char[POINTXY_BYTESIZE],*VPKiy_char=new char[POINTXY_BYTESIZE],
            *ti2_char=new char[POINTXY_BYTESIZE];
    VPKix_char<<VPKix;
    VPKiy_char<<VPKiy;
    ti2_char<<Ti2;
    total_message_len = sizeof(VPKix_char) + sizeof(VPKiy_char) + sizeof(ti2_char) + strlen(mi);   
    const string str=string(VPKix_char)+string(VPKiy_char)+string(ti2_char)+string(mi);
    const char *total_char=str.c_str();
    // cout << "total_char in joint_message:" << total_char << endl;
    return total_char;
}

/* compute HMAC-SHA256  */
void calculate_hmac_sha256(Big gskfj, const unsigned char* message, size_t message_len, unsigned char* output) {  
    unsigned int output_len = EVP_MAX_MD_SIZE;  
    char *key=new char[POINTXY_BYTESIZE];
    key<<gskfj;
    int key_len = sizeof(key);
    string hex_message = to_hex((const char*)message, message_len);
    string hex_key = to_hex((const char*)key, key_len);
    HMAC_CTX* ctx = HMAC_CTX_new();  
    HMAC_Init_ex(ctx, (const unsigned char*)key, key_len, EVP_sha256(), nullptr);
    HMAC_Update(ctx, message, message_len);
    HMAC_Final(ctx, output, &output_len);
    HMAC_CTX_free(ctx);  
}  

/* reported message request assist signature */
void msg_signature(ECn VPKi, const char *mi, Big gskfj, Big Ti2, unsigned char* hmac){
    // compute the length of total message 
    size_t total_message_len;  
    //joint all message
    const char* total_message = joint_message(VPKi, mi, Ti2, total_message_len);
    string hex_message = to_hex(total_message, total_message_len);
    // compute HMAC   
    calculate_hmac_sha256(gskfj, (const unsigned char *)total_message, total_message_len, hmac); 
    // cout << "hmac:" << hmac << endl;
}

/* Function to convert hex string to binary */
void hex_to_binary(const char* hex, char* binary, int binary_len) {
    for (int i = 0; i < binary_len; ++i) {
        sscanf(hex + i * 2, "%2hhx", &binary[i]);
    }
}

/* Function to convert binary back to Big */
Big hex_string_to_big(const char* hexString) {
    // Calculate the number of bytes needed for the binary representation
    int hexLength = strlen(hexString);
    int binaryLength = hexLength / 2;  // Each byte is represented by two hex characters
    char* binaryArray = new char[binaryLength];
    // Convert the hex string to a binary array
    hex_to_binary(hexString, binaryArray, binaryLength);
    // Convert binary array to Big
    Big bigValue = from_binary(binaryLength, binaryArray);
    // Clean up memory
    delete[] binaryArray;
    return bigValue;
}

/* verify HMAC   gskfj,Ti2,VPKi,mi,hmac  */
bool verify_hmac_sha256(Big gskfj, Big Ti2, ECn VPKi, const char* mi, unsigned char* expected_hmac) {
    unsigned char hmac[EVP_MAX_MD_SIZE]; 
    size_t total_message_len;  
    //joint all message
    const char* total_message = joint_message(VPKi, mi, Ti2, total_message_len);
    string hex_message = to_hex(total_message, total_message_len);
    // compute HMAC   
    calculate_hmac_sha256(gskfj, (const unsigned char *)total_message, total_message_len, hmac); 
    // compare calculated HMAC with expected HMAC 
    return (memcmp(hmac, expected_hmac, EVP_MD_size(EVP_sha256())) == 0);

}

/*Vh choose whether to assist signature for message mi  */
void compute_thre_sig(const char* mi,Big q,int num_participants,Big* sigmai0,Big* sigmai1,Big moduli_new[], Big *sigmai){
    // get message 
    Big h = getmsghash(mi,q); 
    // calculate signature--threshold 
    cout << "Threshold signature of each participant:" << endl;
    for (int i = 0; i < num_participants; ++i) {
        sigmai[i] = h * sigmai0[i] + sigmai1[i];
        sigmai[i] %= moduli_new[i];
        // cout << "sigmai" << i << ":" << sigmai[i] << endl;
    }
}

const char* joint_message_h(Big q,const char *mi, char* enc_sigmah, Big Th, int &total_message_len){
    char *th_char=new char[POINTXY_BYTESIZE];
    th_char<<Th;
    total_message_len = sizeof(enc_sigmah) + sizeof(th_char) + strlen(mi);
    const string &istring=string(th_char)+string(mi)+string(enc_sigmah);
    const char *total_message=istring.c_str();
    return total_message;
}

struct return_sigi{
    ECn Ai;
    ECn VPKi2;
    Big sigma;
    ECn Lg;
    Big Ti;
    Big sigi;
};

/* collect threshold signatures  */
return_sigi collect_sign(const char* mi,int num_participants, Big q, Big UL, ECn G, Big *zeroi, Big moduli_new[], Big* sigmai, Big &sigma,
    ECn Lg, int sign_thre,int* chosen_indexes,Big gskfj,Big all_Th[],char** ciphertext,unsigned char** all_Tagh){
    //verify Tag of all vehicles h
    unsigned char Tagh[EVP_MAX_MD_SIZE]; int i;
    for(i=0;i<num_participants;i++){
        int total_message_len;  
        //joint all message
        const char *total_message = joint_message_h(q,mi, ciphertext[i], all_Th[i], total_message_len);
        string hex_message = to_hex(total_message, total_message_len);
        // compute HMAC   
        calculate_hmac_sha256(gskfj, (const unsigned char *)total_message, total_message_len, Tagh); 
        // cout << "Tagh:" << Tagh << endl;
        // cout << "all_Tagh:" << all_Tagh[i] << endl;
        if(memcmp(Tagh, all_Tagh[i], EVP_MD_size(EVP_sha256()))==0) continue;
        else break;
    }
    if (i==num_participants) cout << "All hmac for vehicle_h verification--message M4." << endl;
    else cout << "Hmac for vehicle_h not verified--message M4." << endl;

    // aggregate threshold signature shares 
    Fopen(num_participants,q,UL,zeroi,sigmai,moduli_new,sigma,chosen_indexes,sign_thre);
    cout << "sigma:" << sigma << endl;
    gen_sk_pk Ai = getValues(G,q);
    gen_sk_pk vehi2 = getValues(G,q);
    Big vski2 = vehi2.sk;
    ECn VPKi2 = vehi2.pk;
    Big Ti =getCurTimestamp(); 
    Big hsi = h2(mi, Ai.pk, VPKi2, sigma, Lg, Ti, q);
    Big sigi = Ai.sk + (vski2 * hsi) % q;
    return {Ai.pk,VPKi2,sigma,Lg,Ti,sigi};
}

/* verify threshold signature, message M5 */
void msg_verify(const char* mi,ECn Ai,ECn VPKi2,Big sigma,ECn Lg,Big Ti,Big sigi,Big q,ECn G,ECn pkfj){
    Big hsi = h2(mi, Ai, VPKi2, sigma, Lg, Ti, q);
    ECn left = G;
    left *= sigi;
    ECn right;
    right = mul(Big(1),Ai,hsi,VPKi2);
    // cout << "vehicle left:" << left << endl;
    // cout << "vehicle right:" << right << endl;
    if (left == right) cout << "Signture verification for message M5." << endl;
    else cout << "Signature not verified for message M5." << endl;

    Big hash_mi = getmsghash(mi,q);
    Big sigma_inv = inverse(sigma,q);
    Big l_x,thre_right_x;
    Lg.get(l_x);
    l_x %= q;

    ECn thre_right1 = G;
    thre_right1 *= sigma_inv;
    thre_right1 *= hash_mi;
    ECn thre_right2 = pkfj;
    thre_right2 *= sigma_inv;
    thre_right2 *= l_x;
    ECn thre_right;
    thre_right += thre_right1;
    thre_right += thre_right2;
    thre_right.get(thre_right_x);
    thre_right_x %= q;

    if (l_x == thre_right_x) cout << "Threshold signature verification." << endl;
    else cout << "Threshold signature not verified." << endl;
}

return_sigi veh_get_param(Big q,Big UL,Big* zeroi,ECn G,const char* mi,int num_participants,Big sigmai0[],
    Big sigmai1[],Big moduli_new[],ECn Lg,ECn pkfj,int sign_thre,int* chosen_indexes,ECn k_R,char *encm, 
    int ciphertext_len, Big Tm2,Big Signf,ECn FPKj,Big vskm,gen_sk_pk VPKi_sk_pk,Big vskh[],ECn VPKh[]){
    char decm[160]; unsigned char hmac[EVP_MAX_MD_SIZE];
    ECn Wf = vskm * FPKj;
    Big Wf_x,k_r,gskfj,Ti2,Wistar_x,Th,sigma;
    Wf.get(Wf_x);
    k_R.get(k_r);
    k_r%=q;
    Big hsignf = h_signf(encm, FPKj, Tm2, q);
    if(verify(k_r, Signf, q, hsignf, FPKj, G)){
        aes_dec_cbc(encm,ciphertext_len,Wf_x,decm);
        // split ciphertext chars, recover field 
        char dec_gskfj[160];
        int decsize = sizeof(decm);
        int segmentsize = decsize / 3;
        memcpy(dec_gskfj, decm + 96, 48);
        gskfj = hex_string_to_big(dec_gskfj);
        // gskfj = from_binary(32,dec_gskfj);
    }
    ECn VPKi = VPKi_sk_pk.pk; Big vski = VPKi_sk_pk.sk;
    Ti2 = getCurTimestamp();
    msg_signature(VPKi,mi,gskfj,Ti2,hmac);
    // cout << "Ti2:" << Ti2 << endl;
    // cout << "hmac returned:" << hmac << endl;
    // cout << "gskfj:" << gskfj <<endl;
    bool is_valid = verify_hmac_sha256(gskfj,Ti2,VPKi,mi,hmac);
    if(is_valid) cout << "HMAC is valid for message M3." << endl;
    else cout << "HMAC is invalid for message M3." << endl;
    Big* sigmai = new Big[num_participants];char** ciphertext = new char*[num_participants];
    int* ciph_len = new int[num_participants]; unsigned char* all_Tagh[EVP_MAX_MD_SIZE];
    unsigned char Tagh[EVP_MAX_MD_SIZE]; Big all_Th[num_participants];
    compute_thre_sig(mi,q,num_participants,sigmai0,sigmai1,moduli_new,sigmai);
    for(int i=0;i<num_participants;i++){
        ECn Wistar = vskh[i] * VPKi;
        char sigmai_char[32]; char cip_temp[200];
        Wistar.get(Wistar_x);
        to_binary(sigmai[i],32,sigmai_char);
        aes_enc_cbc(sigmai_char, Wistar_x, cip_temp, ciph_len[i]);
        ciphertext[i]=cip_temp;    
        int total_message_len;  
        Th = getCurTimestamp();
        const char *total_message = joint_message_h(q,mi, cip_temp, Th, total_message_len);
        string hex_message = to_hex(total_message, total_message_len); 
        calculate_hmac_sha256(gskfj, (const unsigned char *)total_message, total_message_len, Tagh); 
        all_Tagh[i] = Tagh;
        all_Th[i] = Th;
    }
    return_sigi returnSig= collect_sign(mi,num_participants,q,UL,G,zeroi,moduli_new,sigmai,sigma,Lg,sign_thre,chosen_indexes,gskfj,all_Th,ciphertext,all_Tagh);
    // cout << "sigi_returnSig:" << returnSig.Ai << endl;
    // cout << "sigi_returnSig:" << returnSig.VPKi2 << endl;
    // cout << "sigi_returnSig:" << returnSig.sigma << endl;
    // cout << "sigi_returnSig:" << returnSig.Lg << endl;
    // cout << "sigi_returnSig:" << returnSig.Ti << endl;
    // cout << "sigi_returnSig:" << returnSig.sigi << endl;
    msg_verify(mi,returnSig.Ai,returnSig.VPKi2,returnSig.sigma,returnSig.Lg,returnSig.Ti,returnSig.sigi,q,G,pkfj);
    return {returnSig.Ai,returnSig.VPKi2,returnSig.sigma,returnSig.Lg,returnSig.Ti,returnSig.sigi};
}

void batchverifytest(int num,ECn G,Big q,const char* mi, ECn* Ai_set, ECn* VPKi2_set, Big* sigma_set,ECn* Lg_set,Big* Ti_set,Big* sigi_set){
    Big visigi_sum = 0;
    ECn vA_sum = q*G;
    ECn vshVPK_sum = q*G;
    for(int i=0;i<num;i++){
        Big hsi = h2(mi,Ai_set[i],VPKi2_set[i],sigma_set[i],Lg_set[i],Ti_set[i],q);
        Big vi=rand(SMALLVECTOR_BITSIZE,2);
        Big vi_sigi = vi * sigi_set[i];   
        visigi_sum += vi_sigi;
        ECn tmp1 = vi * Ai_set[i];
        vA_sum += tmp1;
        Big vihsi = vi * hsi;
        ECn tmp2 = vihsi * VPKi2_set[i];
        vshVPK_sum += tmp2;
    }

    ECn left = G;
    left *= visigi_sum;
    ECn right = vA_sum;
    right += vshVPK_sum;
    // cout << "fog batchverify:" << left << endl;
    // cout << "fog batchverify:" << right << endl;
    if (left == right) cout << "Batch verification successful." << endl;
    else cout << "Batch verification failure." << endl;
}

/* test CPU cycle  */
inline uint64_t readTSC() {
    uint32_t lo, hi;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)hi << 32) | lo; // return the 64-bit number of the combination of the high 32 bits and the low 32 bits
}

/* test CO2 emissions  */ 
double estimateCO2Emissions(double cpuTimeSec, double cpuPowerW) {
// In 2023, the average CO2 emissions per kwh of electricity in all provinces in China will be about 0.608 kg
    const double CO2_PER_KWH = 0.608; 
    // Convert power to kilowatts
    double cpuPowerKw = cpuPowerW / 1000.0;
    // Calculate Power consumption (KWH)
    double energyKWh = cpuPowerKw * cpuTimeSec / 3600.0;
    // compute CO2 emission
    return energyKWh * CO2_PER_KWH;
}

int  main() {
    miracl *mip = mirsys(5000, 16);  // Initiation MIRACL library 
    mip->IOBASE = 16; 
    double cpuPowerW = 600.0; // estimated value
    uint64_t start_cpu = readTSC(); // start timing
    auto start = std::chrono::high_resolution_clock::now();

    int num = 20; // number of message reported
    ECn Ai_set[num];
    ECn VPKi2_set[num];
    Big sigma_set[num];
    ECn Lg_set[num];
    Big Ti_set[num];
    Big sigi_set[num];
    ECn G; Big q;
    const char* mi;

// for convenience, use the same message mi in batchverify phase
    for (int j=0;j<1;j++){
        ECn PK,VPKm,VPKi,Xm,Lg;
        Big h,UL,L_x,vskm,vski,sigma,Tm1,pdm,Tm2,Signf,gskfj;

    // define pariticipates numbers and threshold
        int num_participants = 30; 
        int threshold = 20;int t = 400; int namda = 128; 
    // for convenience, to set reconstruction and privacy threshold , put the number as input here
        int maxt=15;int mint=3;  
        Big* moduli_new = new Big[num_participants]; Big* zeroi = new Big[num_participants]; 
        Big* sigmai0 = new Big[num_participants]; Big* sigmai1 = new Big[num_participants];
        Big* sigmai = new Big[num_participants];
        Big* vehh_sk = new Big[num_participants]; ECn* vehh_pk = new ECn[num_participants];
        int chosen_indexes[num_participants];  
        mi = "这是一条交通道路用度提示消息安徽省合肥市蜀山区九龙路短拥堵500米预计通行时间约30分钟";
        const char *join = "request joining";
        char encm[200]; int ciphertext_len;
        gen_sigma(mi,num_participants,threshold,namda,q,UL,L_x,sigmai0,sigmai1,G,PK,moduli_new,zeroi,Lg,maxt,mint,chosen_indexes);

    // generate communication private and public key of vehicle
        gen_sk_pk vehm = getValues(G,q);
        gen_sk_pk vehi = getValues(G,q);
        for(int i=0;i<num_participants;i++){
            gen_sk_pk vehh = getValues(G,q);
            vehh_sk[i] = vehh.sk;
            vehh_pk[i] = vehh.pk;
        }
    // generate communication private and public key of fog node
        gen_sk_pk fogj = getValues(G,q);
    // generate one-time random value and public key
        gen_sk_pk k_R = getValues(G,q);
    // vehicle group formation
        req_vehgroup(q,G,vehm.sk,vehm.pk,join,Xm,Tm1,pdm);
        Fog_resp(q,G,k_R.sk,k_R.pk,vehm.pk,join,Xm,Tm1,pdm,fogj.sk,fogj.pk,encm,ciphertext_len,Tm2,Signf,gskfj);
        return_sigi returnParam=veh_get_param(q,UL,zeroi,G,mi,num_participants,sigmai0,sigmai1,moduli_new,
            Lg,PK,threshold,chosen_indexes,k_R.pk,encm,ciphertext_len,Tm2,Signf,fogj.pk,vehm.sk,vehi,vehh_sk,vehh_pk);
        
        Ai_set[j] = returnParam.Ai;
        VPKi2_set[j] = returnParam.VPKi2;
        sigma_set[j] = returnParam.sigma;
        Lg_set[j] = returnParam.Lg;
        Ti_set[j] = returnParam.Ti;
        sigi_set[j] = returnParam.sigi;
    }

    batchverifytest(num,G,q,mi,Ai_set,VPKi2_set,sigma_set,Lg_set,Ti_set,sigi_set);

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration = end - start; 
    uint64_t end_cpu = readTSC(); // end timing
    uint64_t cycles = end_cpu - start_cpu; // calculation of CPU cycles consumed
    double co2Emissions = estimateCO2Emissions(duration.count(), cpuPowerW);

    cout << "Algorithm execution time: " << duration.count() << " ms" << endl;
    cout << "Algorithm executed in " << cycles << " CPU cycles." << endl;
    cout << "Estimated CO2 emissions: " << co2Emissions << " kg." << endl;

    // 清理
    mirexit();

    return 0;
}