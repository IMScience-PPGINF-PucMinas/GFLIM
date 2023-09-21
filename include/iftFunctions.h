//
// Created by peixinho on 6/23/15.
//

#ifndef IFT_IFTFUNCTIONS_H
#define IFT_IFTFUNCTIONS_H

/**
 * The Sphere function has n local minima except for the global one.
 * It is continuous, convex and unimodal.
 *
 * Global Minimum f(x*) = 0 at x* = (0,0,...,0)
 *
 * Commonly evaluated in [-5.12, 5.12]
 */
float iftSphereFunction(float *x, int n);

/**
 * The Rosenbrock function, also referred to as the Valley or Banana function.
 * It is unimodal, and the global minimum lies in a narrow, parabolic valley.
 * However, even though this valley is easy to find, convergence to the minimum is difficult (Picheny et al., 2012).
 *
 * Global Minimum f(x*) = 0 at x* = (1,1,...,1)
 *
 * Commonly evaluated in [-5,10]
 */
float iftRosembrockFunction(float* x, int n);


/**
 *  The Ackley function is widely used for testing optimization algorithms. It is characterized by a nearly flat outer region, and a large hole at the centre.
 *  The function poses a risk for optimization algorithms, particularly hillclimbing algorithms, to be trapped in one of its many local minima.
 *
 *  Global Minimum f(x*) = 0 at x* = (0,0,...,0)
 *
 *  Commonly evaluated in [-32.768, 32.768]
 */
float iftAckleyFunction(float* x, int n);


/**
 * The Griewank function has many widespread local minima, which are regularly distributed.
 *
 * Global minimum f(x*) = 0 at x* = (0,0,...,0)
 *
 * Commonly evaluated in [-600, 600]
 */
float iftGrienwankFunction(float* x, int n);


/**
 * The Rastrigin function has several local minima.
 * It is highly multimodal, but locations of the minima are regularly distributed.
 *
 * Global Minimum f(x*) = 0 at x* = (0,0,...,0)
 *
 * Commonly evaluated in [-5.15, 5.12]
 */
float iftRastringinFunction(float *x, int n);


/**
 * @param x
 * @details apporixmate exp(x) with (1 + x / 256) ^ 256
 *          accuracy can be controlled by halving or doubling the constant (256)
 * @return approximation of exp(x)
 */
static inline float iftApproxExpf(float x)
{
    float eaprox = (1.0f + x / 256.0f);
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    return eaprox;
}


static inline double iftApproxExp(double x)
{
    double eaprox = (1.0 + x / 256.0);
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    eaprox *= eaprox;
    return eaprox;
}


#endif //IFT_IFTFUNCTIONS_H
