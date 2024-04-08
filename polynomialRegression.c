#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define __max(a,b) (((a) > (b)) ? (a) : (b))
#define __min(a,b) (((a) < (b)) ? (a) : (b))

int maxGradient = 1000;
struct polyRegression{
    int nPow;
    float stepSize;
    int *coefficients;
};
typedef struct polyRegression polynomialRegression;


float * initCoefficients(float x[], float y[],int nPow){
    // todo: smartly generate coefficients based of of x,y vals
    float *array = (float *)malloc((nPow+1) * sizeof(float));

    if (array == NULL) {
        printf("Memory allocation failed\n");
        return (float *)1;
    }

    // Initialize the array using a loop
    for (int i = 0 ; i < nPow+1; i++) {
        array[i] = rand()/(RAND_MAX*1.0f)*nPow;
    }
    return array;
}

float passThroughFunction(int numPow,float *coefficients,float xValue){
    // using current coefficients calculate the current estimated function value
    float total = 0;
    for(int i = 0;i<numPow+1;i++){
        float current = coefficients[i] * pow(xValue,i);
        total+= current;

    }
    return total;
}


float computeForEachDerivative(float x[],float y[], int valLength,int currentN,int nPow,float *coefficients){
    int chainRulePower = 0;
    if(currentN != 0){
        chainRulePower = currentN-1;
    }
    // printf("currentN:%d chainCoeff%f chainPower %d\n",currentN,chainRuleCoefficient,chainRulePower);
    float sum = 0;
    for(int i = 0; i<valLength;i++){
        float current = (y[i] - passThroughFunction(nPow,coefficients,x[i]))*pow(x[i],chainRulePower);
        sum+=current;
    }
    return sum;
}

float* generatePartialDerivatives(float x[], float y[],int nPow, float *coefficients,int valLength, float * partials){
    for(int n = 0;n<nPow+1;n++){
        float mult = coefficients[n]*n;
        float partialDerivative = (-2.0f/valLength) * mult * computeForEachDerivative(x,y,valLength,n,nPow,coefficients);
        // printf("d: %f\n",partialDerivative);
        partials[n] = partialDerivative;
    }
    return partials;

}

float computeError(float x[],float y[], int valLength,int nPow,float *coefficients){
    
    float error = 0;
    for(int i = 0; i<valLength;i++){
        float current = pow(y[i] - passThroughFunction(nPow,coefficients,x[i]),2);
        // printf("Current: %f\n",current);
        error+=current;
    }
    // printf("error: %f\n",error);

    return error * (1.0f/valLength);
}

float * applyGradientDescent(float* coefficients,float* partialDerivatives,float stepSize,int valLength){
    for(int i = 0;i<valLength;i++){
        float gradient = partialDerivatives[i] * stepSize;
        if(gradient > 0){
            gradient = __min(gradient,maxGradient);
        }
        else{
            gradient = __max(gradient,-1*maxGradient);

        }
        coefficients[i] = coefficients[i] - gradient;
    } 
    return coefficients;
}


int main(){
    srand(time(NULL));
    polynomialRegression poly = (polynomialRegression){.nPow = 1,.stepSize = .2};
    int nPow = 2;
    int arrLength = 6;
    float x[] = {11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
    float y[] = {121,144,169,196,225,256,289,324,361,400,441,484,529,576,625,676,729,784,841};
    float *coeffs = initCoefficients(x,y,nPow);
    for(int i = 0;i<nPow+1;i++){
        printf("Coefficient#%d value: %f\n",i,coeffs[i]);
    }
    float initialerror = computeError(x,y,arrLength,nPow,coeffs);
    float *partials = (float *)malloc((nPow+1) * sizeof(float));

    float *derivs = generatePartialDerivatives(x,y,nPow,coeffs,arrLength,partials);

    printf("initial error: %f\n",initialerror);
    for(int i = 0;i<2;i++){
        for(int i = 0;i<nPow+1;i++){
            printf("Coefficient#%d derivative: %f\n",i,derivs[i]);
        }
        coeffs = applyGradientDescent(coeffs,derivs,poly.stepSize,arrLength);
        derivs = generatePartialDerivatives(x,y,nPow,coeffs,arrLength,derivs);
        for(int i = 0;i<nPow+1;i++){
            printf("Coefficient#%d value: %f\n",i,coeffs[i]);
        }

    }
    
    float finalError = computeError(x,y,arrLength,nPow,coeffs);
    printf("final error: %f\n",finalError);
    for(int i = 0;i<nPow+1;i++){
        printf("Coefficient#%d value: %f\n",i,coeffs[i]);
    }
    float prediction = passThroughFunction(nPow,coeffs,20);
    printf("Prediction: %f",prediction);
    // init coefficients
    
}

