#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define DEBUG 0


double log2(double x){
    return  log(x) * M_LOG2E;
}

int computeNumberOfBits(int n){
    int i=0, j=1;
    while(n>=j){
        j=j<<1;
        i++;
    }
    return i;
    
}


typedef struct rational Rational;

struct rational{
    int numerator;
    int denominator;
};

typedef struct nodelist Node;

struct nodelist{
    int i;
    Node* next;
};
/******************************************************************************
*
*                           COMPLEX NUMBER OPERATIONS
*
*
******************************************************************************/
typedef struct complex Complex;

struct complex{
  double real;
  double imaginary;
};

double complexNorm(Complex c){
    return sqrt( c.real*c.real + c.imaginary*c.imaginary  );
    //return sqrt( c.real*c.real + c.imaginary*c.imaginary  );
}

Complex* complexConjugate(Complex* c){
    Complex* c1 = (Complex*) malloc( sizeof(Complex) );
    c1->real=c1->real;
    c1->imaginary = -c->imaginary;
    
    return c1;
}

Complex* complexDivision(Complex* c1, Complex* c2){
    Complex* res = (Complex*) malloc( sizeof(Complex) );
    double den = (c2->real*c2->real + c2->imaginary*c2->imaginary); //can also be norm*norm
    
    res->real = (c1->real*c2->real + c1->imaginary*c2->imaginary)/den;
    res->imaginary = (c1->imaginary*c2->real - c1->real*c2->imaginary)/den;
    
    return res;
}

Complex* complexMultiplication(Complex* c1, Complex* c2){
    Complex* res = (Complex*) malloc( sizeof(Complex) );
    
    
    res->real = (c1->real*c2->real - c1->imaginary*c2->imaginary);
    res->imaginary = (c1->imaginary*c2->real + c1->real*c2->imaginary);
    
    return res;
}




/******************************************************************************
*
*                           QUANTUM REGISTER OPERATIONS
*
*
******************************************************************************/

typedef struct quantum_register QReg;
struct quantum_register{
  //int* values;
  int size;
  Complex* amplitudes;
};

void printRegister(QReg* reg){
    int i;
    for(i=0; i<reg->size; i++){
        printf("i: %d Cr: %f Ci: %f\n", i, reg->amplitudes[i].real, reg->amplitudes[i].imaginary);
    }
}

QReg* createRegister(int n){
    
    QReg* reg = (QReg*) malloc( sizeof(QReg) );
    //reg->values = (int*) malloc(sizeof(int) * (2<<(n-1)));
    reg->size=2<<(n-1);
    reg->amplitudes = (Complex*) malloc(sizeof(Complex) * (2<<(n-1)));
    reg->amplitudes[0].real=1;
    reg->amplitudes[0].imaginary=0;
    
    return reg;
}


int measureRegister(QReg* reg){
    int i, res=-1;
    double randn = (double) rand()/ (double) RAND_MAX; //random number between 0 and 1
    double acc=0;
    
    //printf("Here's the random number: %f.\n", randn);
    //Measure a random value
    for(i=0; i<reg->size; i++){
        acc+=(complexNorm(reg->amplitudes[i]))*(complexNorm(reg->amplitudes[i]));
        //printf("Acc: %f.\n", acc);
        if (acc>=randn){
            res=i;
            break;
        }
    }
    
    //Didn't measure anything -- the state probably wasn't normalized
    if(res==-1){
        printf("Error. NULL measurement.\n");
        return -1;
    }
   
    //Collapse the state so that further readings will always give the same value
    for(i=0; i<reg->size; i++){
        if(res==i) reg->amplitudes[i].real=1;
        else reg->amplitudes[i].real=0;
        reg->amplitudes[i].imaginary=0;
    }
    
    return res;
}

void normalizeRegister(QReg* reg){
    int i;
    double acc_r=0, acc_i=0;
    for(i=0;i<reg->size;i++){
        acc_r+=reg->amplitudes[i].real;
        acc_i+=reg->amplitudes[i].imaginary;
    }
}

//This one assumes the imaginary part is already 0 for simplicity
void normalizeRegisterREAL(QReg* reg){
    int i;
    double acc=0;
    for(i=0;i<reg->size;i++){
        acc+=(reg->amplitudes[i].real)*(reg->amplitudes[i].real);
    }
    for(i=0;i<reg->size;i++){
        reg->amplitudes[i].real/=sqrt(acc);
    }
}

void collapseRegister1(QReg* reg, int xj, int x, int n){
    int i;
    for(i=0;i<reg->size;i++){
        if(mod_pow(x, i, n)==xj)
            reg->amplitudes[i].real=1;
        else
            reg->amplitudes[i].real=0;
        reg->amplitudes[i].imaginary=0;
    }
    
    normalizeRegisterREAL(reg);
}



/*******************************************************************************
*
*                               FOURIER TRANSFORM
*
*
*******************************************************************************/


int cheatGetR(QReg* reg){
    int period=0;
    int i;
    for(i=0; i<reg->size; i++){
        if (reg->amplitudes[i].real!=0){
            if(period==-1)  period = i;
            else return i-period;
        }
    }
        
    printf("Error. Couldn't find two values != 0.\n");    
    return -1;
}

int cheatGetB(QReg* reg){
    int i;
    for(i=0; i<reg->size; i++){
        if (reg->amplitudes[i].real!=0) return i;
    }
        
    printf("Error. Couldn't find any value != 0.\n");    
    return -1;
}

void quantumDFT(QReg* reg){
    //Complex temp; //This will accumulate the values of a
    Complex* c_array = (Complex*) malloc( sizeof(Complex)*reg->size );
    
    int a, j;
    for(a=0;a<reg->size;a++){
        for(j=0; j<reg->size; j++){
            if(reg->amplitudes[a].real == 0) continue; //reg1 collapsed and no longer has has a in its superposition
            
            c_array[j].real += reg->amplitudes[a].real * pow(reg->size,-0.5) * cos(2*M_PI*a*j/reg->size);
            c_array[j].imaginary+= reg->amplitudes[a].real * pow(reg->size,-0.5) * sin(2*M_PI*a*j/reg->size);
            
        } 
    }
    
    free(reg->amplitudes);
    reg->amplitudes=c_array;
    
}



/*******************************************************************************
*
*                              CLASSICAL FUNCTIONS
*
*
*******************************************************************************/

/* returns (a to the power of b) mod n */
// Needs to be optimized to compute the exponentiation in powers of two (to take log time)
int mod_pow(a, b, n){
    int r=1;
    while(b>0){
        r*=a;
        r%=n;
        b--;
    }
    return r;
}

int find_order(int x, int N){
    int r=1;
    int res = x;
    while(res!=1)
        res=mod_pow(x, r++, N);
    return r-1; //because of the r++ above
}

//Needs to return an array. Do I know the length of the array? Or do I make it a list?
Node* continued_fraction_list(int a, int b){
    
    Node * node = (Node*) malloc(sizeof(node));
    node->i=-1;
    node->next=NULL;
    
    //printf("a: %d b: %d \n", a, b);
    if(a==1){ 
        //printf("%d \n", b); 
        node->i=b; 
        return node;
    }
    if(a==0) return NULL;
    
    int i = floor(b/a);
    //printf("%d \n", i);
    node->i=i;
    node->next=continued_fraction_list( b%a, a);

    return node;
}

int* continued_fraction(int a, int b){
    //Get the values in a list
    Node * n = continued_fraction_list(a, b);
    
    //Find the list length
    Node * n_aux=n;
    int size=1;
    while(n_aux->next!=NULL){
        size++;
        n_aux=n_aux->next;
    }
    n_aux=n;
    
    //Convert list to array
    int i;
    int* array = (int*) malloc( sizeof(int) * (size+1) );
    array[0]=size;
    for(i=1; i<size+1; i++){
        array[i] = n->i;
        n=n->next;
    }
    
    n=n_aux;
    //Free list
    while(n_aux->next != NULL){
        n_aux=n->next;
        free(n);
        n=n_aux;
    }
    
    return array;
 
}

Rational* qth_convergent(int* convergents, int q){
    int n = convergents[0];
    int i;
    Rational* r = malloc( sizeof(Rational) );
    
    Rational aux;
    
    aux.numerator=1;
    aux.denominator=convergents[q];
    //printf("Current convergent is: %d/%d.\n", aux.numerator, aux.denominator);
    
    int num;
    int den;
    for(i=q-1; i>0; i--){
        den = convergents[i] * aux.denominator + aux.numerator;
        num = aux.denominator;
        
        aux.denominator=den;
        aux.numerator=num;
        
        //printf("Current convergent is: %d/%d.\n", aux.numerator, aux.denominator);
    }
    r->denominator = aux.denominator;
    r->numerator = aux.numerator;
    
    //printf("Q-th Convergent is: %d/%d.\n", r->numerator, r->denominator);
    return r;
    
}

int euclidesAlgorithm(int a, int b){
    if(b==0) return a;
    else return euclidesAlgorithm(b, a%b);
}

int getRandomCoprime(int n){
    //The coprime can be larger than n
    //But some numbers cause trouble... mostly for even numbers? 16 didn't work, always yielded 0 on the second measurement even though it is coprime with n
    int r;
    do{ r=rand()%100;   } while (euclidesAlgorithm(r, n)!=1);
    return r;
}

int repeated_squaring(int x, int n){
    if (n < 0) return repeated_squaring(1 / x, -n);
    else if (n == 0) return  1;
    else if (n == 1) return  x;
    else if (n %2==0)  return repeated_squaring(x * x,  n / 2);
    else return x * repeated_squaring(x * x, (n - 1) / 2);
}
/*******************************************************************************
*
*                              SHOR'S ALGORITHM
*
*
*******************************************************************************/
//Basically apply Hadamard 0
void averageStates(QReg* reg){
    int i;
    double temp=1/sqrt(reg->size);
    for(i=0; i<reg->size; i++){
        reg->amplitudes[i].real=temp;
        reg->amplitudes[i].imaginary=0;
    }
    
    normalizeRegister(reg);
}

void applyVx(QReg* reg1, QReg* reg2, int x, int n){
    int i, expmod;
    double temp=1/sqrt(reg1->size);
    
    //Well, this isn't good. But since I initialize it to 1, I guess I need to reset it too...
    reg2->amplitudes[0].real=0; 
    
    for(i=0; i<reg1->size; i++){
        expmod = mod_pow(x, i, n);
        reg2->amplitudes[expmod].real+=temp;
        reg2->amplitudes[expmod].imaginary=0;
    }
    
    normalizeRegisterREAL(reg2);
}

int quantumShor(int N, int q, int x){
    if(N==0){
        printf("N is zero, jokester.\n");
        return -1;
    }
    
    //int t = (int) floor(log2(q)+1)-1;
    //int n = (int) floor(log2(N)+1);
    int t = computeNumberOfBits(q)-1;
    int n = computeNumberOfBits(N);
    
    //printf("Here's t: %d.\n", t);
    
    QReg * reg1 = createRegister(t);
    QReg * reg2 = createRegister(n);
    
    
    if(DEBUG){
    printf("==============Initialized registers. Should be defaulted to 0.=====================\n");
    printf("\nRegister1.\n");
    printRegister(reg1);
    printf("\nRegister2.\n");
    printRegister(reg2);
    }
    
    averageStates(reg1);
    if(DEBUG){
    printf("==============Register1 in equal superposition over all states.====================\n");
    printf("\nRegister1.\n");
    printRegister(reg1);
    printf("\nRegister2.\n");
    printRegister(reg2);
    }
    
    printf("Starting powering.\n");
    applyVx(reg1, reg2, x, N);
    printf("Finished powering.\n");
    if(DEBUG){
    printf("==============Register2 in equal superposition over powers of x mod N.====================\n");
    printf("\nRegister1.\n");
    printRegister(reg1);
    printf("\nRegister2.\n");
    printRegister(reg2);
    }
    
    int measure2 = measureRegister(reg2);
    collapseRegister1(reg1, measure2, x, N);
    if(DEBUG){
    printf("==============Measured Register2 & Collapsed Register1====================\n");
    printf("Measured: %d.\n", measure2);
    printf("\nRegister1.\n");
    printRegister(reg1);
    printf("\nRegister2.\n");
    printRegister(reg2);
    }
    
    printf("Starting DFT.\n");
    quantumDFT(reg1);
    printf("Finished DFT.\n");
    //quantumDFTFAKE(reg1);
    if(DEBUG){
    printf("==============Applied Quantum DFT to Register 1====================\n");
    printf("\nRegister1.\n");
    printRegister(reg1);
    printf("\nRegister2.\n");
    printRegister(reg2);
    }
    int measure1 = measureRegister(reg1);
    if(DEBUG){
    printf("==============Measured Register 1====================\n");
    printf("Measured: %d.\n", measure1);
    printf("\nRegister1.\n");
    printRegister(reg1);
    }
    printf("measurement2: %d; measurement1: %d.\n", measure2, measure1);

    return measure1;
    
    //if(measure1==0) printf("Shit happened, need to repeat this whole thing (can keep the x, though).\n");
}



int shorX(int N, int x){
    //pick n^2<q<2*n^2
    //q is the smallest power of 2 larger than n^2
    int q=1;
    int N2 = N*N;
    while(q<=N2) //<=, because I'll decrease t next
        q*=2;
        
    //int n = (int) floor(log2(q)+1)-1; 
    int n = computeNumberOfBits(q)-1;
    
    printf("Here's q: %d and N: %d and N-2: %d and 2N-2: %d.\n", q, N, N2, 2*N2);
    
    int size = 2<<(n-1);
    int cand_r;
    
    int factor=0;
    
    do{
        //Measuring 0 doesn't help, but can happen with 16.(6)% chance
        do {
           cand_r = quantumShor(N, q, x);
           printf("Candidate period is: %d.\n", cand_r);
        } while( cand_r==0 );
        
        int i;
        int* tempar=continued_fraction(cand_r, size);
        
        double approx;
        double threshold = 1.0/(2*q);
        double target = (0.0+cand_r)/q;
        
        Rational* r;
        //Go through all the convergents with denominator smaller than n
        //Only accept very good approximations (which works only in peaks)
        //the approximations are bounded by a certain value so that we only accept
        //cand_r that correspond to the peak probabilities
        for(i=1; i<=tempar[0]; i++){    
            r=qth_convergent(tempar, i);   
            
            if(r->denominator > N) break; //denominators are in crescent order
            approx = (0.0 + r->numerator) / r->denominator;
            
            if (fabs(approx-target)<threshold){
                factor=r->denominator;
                break;
            }
        }
    
    } while( factor<=0 );
    
    printf("Here's the found factor: %d.\n", factor);
    if(mod_pow(x, factor, N) == 1){
        printf("Found the period of %d mod %d: %d\n", x, N, factor); 
        return factor;
    }
    else{
        //Try to find the other factor 
        factor *= shorX(N, mod_pow(x, factor, N));
        //factor *= shorX(N, repeated_squaring(x, factor));
        printf("The period of %d is: %d\n", x, factor);
        return factor;
    }
    
    
}




/*
int shor(int N){
    //pick n^2<q<2*n^2
    //q is the smallest power of 2 larger than n^2
    int q=1;
    int N2 = N*N;
    while(q<=N2) //<=, because I'll decrease t next
        q*=2;
        
    //Choose random x, coprime with N
    int x=2;
    int n = (int) floor(log2(q)+1)-1; 
    
    
    printf("Here's q: %d and N: %d and N-2: %d and 2N-2: %d.\n", q, N, N2, 2*N2);
    
    int size = 2<<(n-1);
    int cand_r;
    
    int factor=0;
    
    do{
        //Measuring 0 doesn't help, but can happens with 16.(6)% chance
        do {
           cand_r = quantumShor(N, q, x);
        } while( cand_r==0 );
        
        int i;
        int* tempar=continued_fraction(cand_r, size);
        
        double approx;
        double threshold = 1.0/(2*q);
        double target = (0.0+cand_r)/q;
        
        Rational* r;
        //Go through all the convergents with denominator smaller than n
        //Only accept very good approximations (which works only in peaks)
        //the approximations are bounded by a certain value so that we only accept
        //cand_r that correspond to the peak probabilities
        for(i=1; i<=tempar[0]; i++){    
            r=qth_convergent(tempar, i);   
            
            if(r->denominator > N) break; //denominators are in crescent order
            approx = (0.0 + r->numerator) / r->denominator;
            
            if (fabs(approx-target)<threshold){
                factor=r->denominator;
                break;
            }
        }
    
    } while( factor<=0 );
    
    printf("Here's the found factor: %d.\n", factor);
    if(mod_pow(x, factor, N) == 1){
        printf("Found the period of %d mod %d: %d\n", x, N, factor); 
        return factor;
    }
    else{
        //Try to find the other factor 
        factor *= quantumShor(N, q, pow(x, factor));
        printf("New factor is: %d\n", factor);
    }
    
    
}
*/



int main(){
    srand(time(NULL));

    int x = 2;
    int N = 15*13;
    printf("Input the value N to be factored:\n");
    printf("N must be a product of two different primes larger than 2.\n");
    scanf("%d", &N);
    
    while(N%2==0){
        N=N/2;
        printf("One factor is 2. Retrying with N=%d.\n", N);
    }
    
    if(N%2==0) x=getRandomCoprime(N);
    
    int p = shorX(N, x);
    
    int repeat = 1;
    
    while( repeat ){
        
        while(p%2!=0){
            //I'll assume I don't need to keep track of older xs and it will choose different ones by chance alone
            x=getRandomCoprime(N);
            printf("The period (classically measured) is: %d.\n", p);
            printf("Because it was odd, we're trying again with a different coprime. The number is: %d.\n", x);
            p = shorX(N, x);
        }
        
        //This gives me two factors -- they can be N or 1, and they can be composite factors
        printf("gcd(%d, %d): %d; gcd(%d, %d): %d;\n", mod_pow(x, p/2, N) - 1, N, euclidesAlgorithm(mod_pow(x, p/2, N) - 1, N), mod_pow(x, p/2, N) + 1, N, euclidesAlgorithm(mod_pow(x, p/2, N) + 1, N));
        repeat = ((mod_pow(x, p/2, N) - 1) == 1) || ((mod_pow(x, p/2, N) - 1) == N) || ((mod_pow(x, p/2, N) + 1) == 1) || ((mod_pow(x, p/2, N) + 1) == N);
        if(repeat){ 
            printf("Found trivial factors -- retrying with a different x!");
            p=1;
        }
    }
    
    printf("These are the two factors of %d: %d and %d.\n", N, euclidesAlgorithm(mod_pow(x, p/2, N) - 1, N), euclidesAlgorithm(mod_pow(x, p/2, N) + 1, N));
    
    return 0;
}
