#include <iostream>
#include <Eigen/Dense>
#include <complex>
using namespace Eigen;
using namespace std;

int main()
{
    MatrixXd m = MatrixXd::Random(3, 3);

    m = (m + MatrixXd::Constant(3, 3, 1.2)) * 50;
    cout << "m =" << endl << m << endl;
    VectorXd v(3);

    v << 1, 2, 3;
    cout << "m * v =" << endl << m * v << endl;

    MatrixXcd test(3, 3);

    test = MatrixXcd::Zero(3, 3);
    std::complex<double> c(0, 1);

    // test(0,0) = c;

    cout << "test=" << endl << test << endl;

    cout << "adjoint=" << endl << test.adjoint() << endl;
}
