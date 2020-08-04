#ifndef __MATRICES_H__
#define __MATRICES_H__

#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <stdio.h>
// #include <stdlib.h>
#include <math.h>
#include <iostream>


// just so it's tidy
static Eigen::MatrixXcd
gamma13(int i)
{
    Eigen::MatrixXcd gamma(4, 4);

    gamma = Eigen::MatrixXcd::Zero(4, 4);

    if (i == 0) {
        gamma << 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
    } else if (i == 1)      {
        //	gamma << 0, 0, 0, I, 0, 0, I, 0, 0, -I, 0, 0, -I, 0, 0, 0;
        gamma.real() << 0, 0, 0, -1, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, 0;
    } else if (i == 2)      {
        gamma.imag() << 0, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, 1, 0, 0, 0;
    } else if (i == 3)      {
        //	gamma<<0, 0, I, 0, 0, 0, 0, -I, -I, 0, 0, 0, 0, I, 0, 0;
        gamma.real() << 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, -1, 0, 0;
    } else   {
        printf("I'm sorry there are only four gamma matrices to be had ");
    }

    return gamma;
}

// just so it's tidy
static Eigen::MatrixXcd
gamma10(int i)
{
    Eigen::MatrixXcd gamma(1, 1);

    gamma = Eigen::MatrixXcd::Zero(1, 1);

    /// no need for this either
    gamma << 1;
    return gamma;
}

// just so it's tidy
static Eigen::MatrixXcd
gamma01(int i)
{
    Eigen::MatrixXcd gamma(1, 1);

    gamma = Eigen::MatrixXcd::Zero(1, 1);

    /// no need for this either
    gamma << 1;
    return gamma;
}

// just so it's tidy
static Eigen::MatrixXcd
gamma20(int i)
{
    Eigen::MatrixXcd gamma(2, 2);

    gamma = Eigen::MatrixXcd::Zero(2, 2);

    if (i == 1) {
        gamma << 1, 0, 0, -1;
    } else if (i == 2)      {
        gamma << 0, 1, 1, 0;
    } else   {
        printf("Gamma failure!\n");
    }

    return gamma;
}

// just so it's tidy
static Eigen::MatrixXcd
gamma11(int i)
{
    Eigen::MatrixXcd gamma(2, 2);

    gamma = Eigen::MatrixXcd::Zero(2, 2);

    if (i == 1) {
        gamma << 1, 0, 0, -1;
    } else if (i == 2)      {
        gamma << 0, 1, -1, 0;
    } else   {
        printf("Gamma failure!\n");
    }

    return gamma;
}

// just so it's tidy
static Eigen::MatrixXcd
gamma02(int i)
{
    Eigen::MatrixXcd gamma(2, 2);

    gamma = Eigen::MatrixXcd::Zero(2, 2);
    if (i == 1) {
        gamma.imag() << 1, 0, 0, -1;
    } else if (i == 2)      {
        gamma << 0, 1, -1, 0;
    } else   {
        printf("Gamma failure!\n");
    }

    return gamma;
}

// just so it's tidy
static Eigen::MatrixXcd
gamma03(int i)
{
    Eigen::MatrixXcd gamma(2, 2);

    gamma = Eigen::MatrixXcd::Zero(2, 2);

    if (i == 1) {
        gamma.imag() << 1, 0, 0, -1;
    } else if (i == 2)      {
        gamma << 0, 1, -1, 0;
    } else if (i == 3)      {
        gamma.imag() << 0, 1, 1, 0;
    } else   {
        printf("Gamma failure!\n");
    }

    return gamma;
}

// I need to write a tensor product between matrices
// This assumes square matrices
static Eigen::MatrixXcd
tens(Eigen::MatrixXcd a, Eigen::MatrixXcd b)
{
    int i, j;

    i = a.cols();
    j = b.cols();

    int mps = i * j;


    Eigen::MatrixXcd mp(mps, mps);

    for (int k = 0; k < mps; k++) {
        for (int k2 = 0; k2 < mps; k2++) {
            // I am using that k, k2, i, j are integers!
            mp(k, k2) = a(k / j, k2 / j) * b(k % j, k2 % j);
        }
    }

    return mp;
}

// function to do the commutator
static Eigen::MatrixXcd
com(Eigen::MatrixXcd m)
{
    int s = m.cols();
    Eigen::MatrixXcd com(s * s, s * s);

    Eigen::MatrixXcd ident(s, s);

    ident.setIdentity();

    com = tens(ident, m) - tens(m.transpose(), ident);

    return com;
}

// function to do the anti-commutator
static Eigen::MatrixXcd
acom(Eigen::MatrixXcd m)
{
    int s = m.cols();
    Eigen::MatrixXcd acom(s * s, s * s);

    Eigen::MatrixXcd ident(s, s);

    ident.setIdentity();

    acom = tens(ident, m) + tens(m.transpose(), ident);

    return acom;
}

// simple function to give me the represenations
static Eigen::MatrixXcd
Lxy(int n, int direct)
{
    std::complex<double> ii(0., 1.); // This is the imaginary unit
    Eigen::MatrixXcd m(n, n);

    m = Eigen::MatrixXcd::Zero(n, n);

    if (direct == 1) {
        for (int i = 0; i < n - 1; i++) {
            /// n= 2s+1 -> s=(n-1)/2
            m(i, i + 1) = ii * sqrt(((n + 1) / 2.) * (2 * i + 2) - (i + 2) * (i + 1) ) / 2.;
            m(i + 1, i) = ii * sqrt(((n + 1) / 2.) * (2 * i + 2) - (i + 2) * (i + 1) ) / 2.;
        }
    } else if (direct == 2)     {
        for (int i = 0; i < n - 1; i++) {
            m(i, i + 1) = sqrt(((n + 1) / 2.) * (2 * i + 2) - (i + 2) * (i + 1) ) / 2;
            m(i + 1, i) = -sqrt(((n + 1) / 2.) * (2 * i + 2) - (i + 2) * (i + 1) ) / 2;
        }
    } else if (direct == 3)     {
        for (int i = 0; i <= n - 1; i++) {
            m(i, i) = ii * ((n + 1.) / 2 - (i + 1) );
        }
    } else   {
        printf("I'm sorry there are only three directions in SU(2) ");
    }

    return m;
}

#endif // ifndef __MATRICES_H__
