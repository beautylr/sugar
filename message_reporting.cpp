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

    return true;
}

Big generate_ul(Big L){
    Big UL;
    return UL;
}

/*generate random number*/
void randomValue(int n, Big p0, Big* randV){
    for(int i = 0; i < n; i++) randV[i] =rand(p0);
}


/* distributing shares of secret key sk */
void distributesk(Big S, int num_moduli, Big shares[], Big moduli_new[]) {
    for (int i = 0; i < num_moduli; i++) {
        //check whether current modulus is zero 
        if (moduli_new[i] == 0) {
            cout << "Error：moduli[" << i << "] is zero, unable to distribute secret shares " << endl;
            continue;
        }
    }
}
/* distributing shares of random numbers  */
void distributeRandomV(Big S, int rowi, int num_moduli, Big** Xij, Big moduli_new[]) {
    // cout << "rowi:" <<rowi << endl;
    for (int i = 0; i < num_moduli; i++) {

    }
}

/* compute P_A, P_A_Ba*/
void compute_PA_PABa(Big* moduli_n, int size, int maxt, int mint, Big& P_A, Big& PA_Ba) {
    int size_new = size-1;
    Big* moduli_new = new Big[size_new];
    for (int i = 0; i < size_new; ++i) moduli_new[i] = moduli_n[i+1];
    // sorting arrays

}

void generate_random_in_range(big lower_bound, big upper_bound, big random_result) {
    big range = mirvar(0);  // storage range
    big random_offset = mirvar(0);  
    mirkill(range);
    mirkill(random_offset);
}

/* recovering secret s from S: s = S % p0 */
Big recoverOriginalSecret(Big S, Big p0) {
    Big secret;
    secret = S % p0;
    return secret;
}

/* compute sigmai0 and sigmai1 of ECDSA in presign phase  */
void preSign(int num_participants){
    Big theta_inv;
    theta_inv = inverse(theta, q); 
}

/* imitate randomly participates  */
void choose_party(int threshold, int* chosen_indexes){
    bool selected[threshold] = { false };  // mark the selected participants
    // random number generator
    irand(time(0));
    int count = 0; // randomly selected threshold participants index
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


    delete[] padded_plaintext;
}

bool aes_dec_cbc(char* ciphertext, int ciphertext_len, Big Wf_x, char* plaintext) {
    aes a;
    char iv[32] = {0}; 
    char key[32];
    to_binary(Wf_x, 32, key); // translate key to char arrays
    aes_init(&a, MR_CBC, 32, key, iv);
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
}


bool verify(Big r, Big s, Big q, Big h, ECn Pub, ECn G){
    Big u1,u2,v;
    ECn RR;
    if (r>=q || s>=q) cout << "Signature NOT verified" << endl;
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
    delete[] binaryArray;
    return bigValue;
}

/* verify HMAC   gskfj,Ti2,VPKi,mi,hmac  */
bool verify_hmac_sha256(Big gskfj, Big Ti2, ECn VPKi, const char* mi, unsigned char* expected_hmac) {
    unsigned char hmac[EVP_MAX_MD_SIZE]; 
    size_t total_message_len;  
    //joint all message
    const char* total_message = joint_message(VPKi, mi, Ti2, total_message_len);
    return (memcmp(hmac, expected_hmac, EVP_MD_size(EVP_sha256())) == 0);

}

/*Vh choose whether to assist signature for message mi  */
void compute_thre_sig(const char* mi,Big q,int num_participants,Big* sigmai0,Big* sigmai1,Big moduli_new[], Big *sigmai){
    // get message 
    Big h = getmsghash(mi,q); 
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

    if (l_x == thre_right_x) cout << "Threshold signature verification." << endl;
    else cout << "Threshold signature not verified." << endl;
}

return_sigi veh_get_param(Big q,Big UL){
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
