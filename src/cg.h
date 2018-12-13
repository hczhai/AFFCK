#pragma once
#include "cluster.h"

namespace affck {

	class MVector { // vector
	public:
		virtual double& operator[](int i) = 0;
		virtual int Length() const = 0;
		virtual vector<double> ArrayForm() const = 0;
	};

	class MVectorFMat :public MVector { // vector from a matrix column
		vector<vector<double> >& a;
		const int j;
	public:
		MVectorFMat(vector<vector<double> >& a, int j);
		double& operator[](int i);
		int Length() const;
		vector<double> ArrayForm() const;
	};

	class MVectorFVec :public MVector { // vector from a vector
		vector<double>& a;
	public:
		MVectorFVec(vector<double>& a);
		double& operator[](int i);
		int Length() const;
		vector<double> ArrayForm() const;
	};

	typedef double (ParamFitting::*QF) (const vector<double>&) const; // pointers to derivs functions
	typedef const vector<double>(ParamFitting::*QFp) (const vector<double>&) const;
	typedef const vector<vector<double> >(ParamFitting::*QFpp) (const vector<double>&) const;

	class CG {
	public:
		static int cg_imax;
		static double cg_epsilon;
		static int cg_jmax;
		static double cg_epsilon_n;
		static const vector<double> SolveNonLinear(QF f, QFp fp, QFpp fpp, const vector<double>& xx, int *step, ParamFitting *pf);
		static const vector<double> SolveLinearHH(const vector<vector<double> >& a, const vector<double>& b);
	};

}
