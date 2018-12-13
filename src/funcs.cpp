
#include "funcs.h"

namespace affck {
    
    GXL::~GXL() { }
	GFunction::~GFunction() { }
	GFunction::GFunction(const int index) : Index(index) { }

	GP2C::GP2C(const double l, const int index, const vector<int>& pi)
		: GFunction(index), l(l), PointIndex(pi) {
		Num = 3;
	}
	const string GP2C::Formula(const vector<double>& x) const {
		return "GP2C = " + ToString(x[0])
			+ " * l^2 - 2*l*" + ToString(x[1])
			+ " + " + ToString(x[2]);
	}
	double GP2C::F(const vector<double>& x) const {
		return x[0] * l*l - 2 * l*x[1] + x[2];
	}
	const vector<double> GP2C::Fp(const vector<double>& x) const {
		vector<double> r = vector<double>(Num);
		r[0] = l*l;
		r[1] = -2 * l;
		r[2] = 1;
		return r;
	}
	const vector<vector<double> > GP2C::Fpp(const vector<double>& x) const {
		vector<vector<double> > r = vector<vector<double> >(Num);
		for (int i = 0; i < Num; i++) {
			r[i] = vector<double>(Num);
			for (int j = 0; j < Num; j++)
				r[i][j] = 0;
		}
		return r;
	}
	
	const vector<int>& GP2C::GetPointIndex() const { return PointIndex; }
	double GP2C::Fl() const { return Param[0] * l*l - 2 * l*Param[1] + Param[2]; }
	double GP2C::Fpl() const { return 2 * Param[0] * l - 2 * Param[1]; }
	double GP2C::Fppl() const { return 2 * Param[0]; }

	const vector<int> GConst::SPointIndex = vector<int>();
	GConst::GConst(const int index) : GFunction(index) { Num = 1; }
	const string GConst::Formula(const vector<double>& x) const {
		return "GConst = " + ToString(x[0]);
	}

	double GConst::F(const vector<double>& x) const {
		return x[0];
	}
	const vector<double> GConst::Fp(const vector<double>& x) const {
		vector<double> r(1); r[0] = 1;
		return r;
	}
	const vector<vector<double> > GConst::Fpp(const vector<double>& x) const {
		vector<vector<double> > r(1);
		vector<double> rr(1); rr[0] = 0; r[0] = rr;
		return r;
	}

	const vector<int>& GConst::GetPointIndex() const { return SPointIndex; }
	double GConst::Fl() const { return Param[0]; }
	double GConst::Fpl() const { return 0; }
	double GConst::Fppl() const { return 0; }

	GLJ::GLJ(const double l, const int index, const vector<int>& pi)
		: GFunction(index), l(l), PointIndex(pi) {
		Num = 2;
	}
	const string GLJ::Formula(const vector<double>& x) const {
		return "GLJ = " + ToString(x[0]) + " / l^" + IntToString(lj_power)
			+ " - " + ToString(x[1]) + " / l^6";
	}

	double GLJ::F(const vector<double>& x) const {
		return x[0] / pow(l, lj_power) - x[1] / pow(l, 6);
	}
	const vector<double> GLJ::Fp(const vector<double>& x) const {
		vector<double> r = vector<double>(Num);
		r[0] = 1 / pow(l, lj_power);
		r[1] = -1 / pow(l, 6);
		return r;
	}
	const vector<vector<double> > GLJ::Fpp(const vector<double>& x) const {
		vector<vector<double> > r(2);
		vector<double> rr(2); rr[0] = 0, rr[1] = 0;
		r[0] = rr, r[1] = rr;
		return r;
	}

	const vector<int>& GLJ::GetPointIndex() const { return PointIndex; }
	double GLJ::Fl() const { return Param[0] / pow(l, lj_power) - Param[1] / pow(l, 6); }
	double GLJ::Fpl() const {
		return (-lj_power)*Param[0] / pow(l, lj_power + 1)
			+ 6 * Param[1] / pow(l, 7);
	}
	double GLJ::Fppl() const {
		return lj_power*(1 + lj_power)*Param[0] / pow(l, lj_power + 2)
			- 42 * Param[1] / pow(l, 8);
	}
}
