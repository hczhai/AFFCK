#pragma once
#include <cmath>
#include "base.h"

namespace affck {

	class Derivs { // Derivative formulas
	public:
		class DFuntion {
		public:
			virtual double F(const vector<double>& x) const = 0;
			virtual const vector<double> Fp(const vector<double>& x) const = 0;
			virtual const vector<vector<double> > Fpp(const vector<double>& x) const = 0;
		};

		class DTheta :public DFuntion {
			static const int Num;
		private:
			void rpl_lpx(const vector<double>& x, vector<double>& rpl, vector<vector<double> >& lpx) const;
		public:
			double F(const vector<double>& x) const;
			const vector<double> Fp(const vector<double>& x) const;
			const vector<vector<double> > Fpp(const vector<double>& x) const;

		};

		class DL : public DFuntion {
			static const int Num;
		public:
			double F(const vector<double>& x) const;
			const vector<double> Fp(const vector<double>& x) const;
			const vector<vector<double> > Fpp(const vector<double>& x) const;
		};
	};

}