#pragma once
#include "base.h"

namespace affck {
	
	class GXL { // dF/dl
	public:
		virtual const vector<int>& GetPointIndex() const = 0; // the index of each atom in LA included in the interaction
		virtual double Fl() const = 0;
		virtual double Fpl() const = 0;
		virtual double Fppl() const = 0;
		virtual ~GXL();
	};

	class GFunction :public GXL { // dF/dp
	public:
		int Num; // The number of parameters need for the function
		const int Index;
		vector<double> Param;

		virtual double F(const vector<double>& x) const = 0;
		virtual const vector<double> Fp(const vector<double>& x) const = 0;
		virtual const vector<vector<double> > Fpp(const vector<double>& x) const = 0;
		virtual const string Formula(const vector<double>& x) const = 0;

		GFunction(const int index);
		virtual ~GFunction();
	};

	class GP2C : public GFunction { // quadratic form term
		const double l;
		const vector<int> PointIndex;
	public:
		GP2C(const double l, const int index, const vector<int>& pi);
		const string Formula(const vector<double>& x) const;
		double F(const vector<double>& x) const;
		const vector<double> Fp(const vector<double>& x) const;
		const vector<vector<double> > Fpp(const vector<double>& x) const;
		const vector<int>& GetPointIndex() const;
		double Fl() const;
		double Fpl() const;
		double Fppl() const;
	};

	class GConst : public GFunction { // const term
		const static vector<int> SPointIndex;
	public:
		GConst(const int index);
		const string Formula(const vector<double>& x) const;
		double F(const vector<double>& x) const;
		const vector<double> Fp(const vector<double>& x) const;
		const vector<vector<double> > Fpp(const vector<double>& x) const;
		const vector<int>& GetPointIndex() const;
		double Fl() const;
		double Fpl() const;
		double Fppl() const;
	};

	class GLJ : public GFunction { // L-J form term
		const double l;
		const vector<int> PointIndex;
	public:
		static int lj_power;
		static int lj_distance;
		GLJ(const double l, const int index, const vector<int>& pi);
		const string Formula(const vector<double>& x) const;
		double F(const vector<double>& x) const;
		const vector<double> Fp(const vector<double>& x) const;
		const vector<vector<double> > Fpp(const vector<double>& x) const;
		const vector<int>& GetPointIndex() const;
		double Fl() const;
		double Fpl() const;
		double Fppl() const;
	};

}
