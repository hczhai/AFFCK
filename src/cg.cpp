
#include "cg.h"

namespace affck {

	MVectorFMat::MVectorFMat(vector<vector<double> >& a, int j) :a(a), j(j) {}
	double& MVectorFMat::operator[](int i) { return a[i][j]; }
	int MVectorFMat::Length() const { return a.size(); }
	vector<double> MVectorFMat::ArrayForm() const {
		vector<double> x(Length());
		for (int i = 0; i < Length(); i++)
			x[i] = a[i][j];
		return x;
	}

	MVectorFVec::MVectorFVec(vector<double>& a) : a(a) {}
	double& MVectorFVec::operator[](int i) { return a[i]; }
	int MVectorFVec::Length() const { return a.size(); }
	vector<double> MVectorFVec::ArrayForm() const { return a; }

	template<class T> const vector<T> LineCut(const vector<T>& x, int n) {
		vector<T> h(n);
		for (int i = 0; i < n; i++) h[i] = x[i];
		return h;
	}

	const vector<double> BackSubstitution(const vector<vector<double> >& u,
		const vector<double>& b) {
		int n = b.size();
		vector<double> x(n);
		for (int i = n - 1; i >= 0; i--) {
			x[i] = b[i];
			for (int j = n - 1; j > i; j--)
				x[i] -= x[j] * u[i][j];
			x[i] /= u[i][i];
		}
		return x;
	}

	double InnerProduct(const vector<double>& a, MVector *b) {
		double r = 0;
		for (size_t i = 0; i < a.size(); i++)
			r += a[i] * (*b)[i];
		return r;
	}

	double InnerProduct(const vector<double>& a, const vector<double>& b) {
		double r = 0;
		for (size_t i = 0; i < a.size(); i++)
			r += a[i] * b[i];
		return r;
	}

	const vector<double> VAct(const vector<vector<double> >& x, const vector<double>& y) {
		const int n = y.size();
		vector<double> vv = y;
		for (size_t i = 0; i < x.size(); i++) {
			vector<double> v = x[i];
			double beta_k = InnerProduct(v, v);
			if (beta_k == 0)continue;
			double gamma = InnerProduct(v, vv);
			for (size_t k = 0; k < vv.size(); k++)
				vv[k] += v[k] * (-2 * gamma / beta_k);
		}
		return vv;
	}

	const vector<vector<vector<double> > > HouseHolder(const vector<vector<double> >& a) {
		vector<vector<double> > g = a;
		int n = g[0].size(), m = g.size();
		vector<vector<double> > vs(n);
		for (int k = 0; k < n; k++) {
			double sq_sum = 0;
			vector<double> v(m);
			for (int i = k; i < m; i++) {
				v[i] = g[i][k];
				sq_sum += v[i] * v[i];
			}
			double alpha_k = sqrt(sq_sum);
			if (g[k][k] > 0) alpha_k = -alpha_k;
			v[k] += -alpha_k;
			double beta_k = InnerProduct(v, v);
			if (beta_k == 0) continue; // pass if current column = 0
			vs[k] = v;
			for (int j = k; j < n; j++) { // trans sub-matrix
				MVectorFMat gj(g, j);
				double gamma = InnerProduct(v, &gj);
				for (int i = 0; i < gj.Length(); i++)
					gj[i] += v[i] * (-2 * gamma / beta_k);
			}
		}
		vector<vector<vector<double> > > r(2);
		r[0] = g, r[1] = vs;
		return r;
	}

	const vector<double> CG::SolveLinearHH(const vector<vector<double> >& a,
		const vector<double>& b) {
		const int n = a[0].size();
		const vector<vector<vector<double> > >& g = HouseHolder(a);
		const vector<double>& c = LineCut(VAct(g[1], b), n);
		const vector<vector<double> >& rr = LineCut(g[0], n);
		const vector<double>& x = BackSubstitution(rr, c);
		return x;
	}

	vector<double> Times(const vector<double>& x, const double k) {
		vector<double> r(x.size());
		for (size_t i = 0; i < r.size(); i++)r[i] = x[i] * k;
		return r;
	}

	vector<double> Add(const vector<double>& x, const vector<double>& y) {
		vector<double> r(x.size());
		for (size_t i = 0; i < r.size(); i++)r[i] = x[i] + y[i];
		return r;
	}

	vector<double> Times(const vector<vector<double> >& x, const vector<double>& y) {
		vector<double> r(x.size());
		int n = x.size();
		for (int i = 0; i < n; i++) {
			r[i] = 0;
			for (size_t j = 0; j < x[i].size(); j++)
				r[i] += x[i][j] * y[j];
		}
		return r;
	}

	int cg_imax = 100;
	double cg_epsilon = 1e-7;
	int cg_jmax = 10;
	double cg_epsilon_n = 1e-6;

	const vector<double> CG::SolveNonLinear(QF f, QFp fp, QFpp fpp, const vector<double>& xx, int *step, ParamFitting *pf) {
		const int n = xx.size();
		int i = 0;
		int k = 0;
		vector<double> r = Times((pf->*fp)(xx), -1);
		vector<double> d = r;
		double delta_new = InnerProduct(r, r);
		double delta_0 = delta_new;
		double delta_mid = 0.0;
		double dff = 0.0;
		double ff = 0.0;
		vector<double> x = xx;
		while (i<cg_imax && delta_new > cg_epsilon*cg_epsilon*delta_0) {
			int j = 0;
			const double delta_d = InnerProduct(d, d);
			double alpha = 0.0;
			do {
				alpha = -InnerProduct((pf->*fp)(x), d)
					/ InnerProduct(d, Times((pf->*fpp)(x), d));
				vector<double> s = x;
				x = Add(s, Times(d, alpha));
				while ((dff = (pf->*f)(x)) > (pf->*f)(s)) {
					alpha = abs(alpha) / 2;
					x = Add(s, Times(d, alpha));
				}
				j++;
			} while (j<cg_jmax && alpha*alpha*delta_d>cg_epsilon_n*cg_epsilon);
			if (ff - dff < cg_epsilon) { ff = dff; break; }
			vector<double> s2 = r;
			r = Times((pf->*fp)(x), -1);
			double delta_old = delta_new;
			delta_mid = InnerProduct(r, s2);
			delta_new = InnerProduct(r, r);
			double beta = (delta_new - delta_mid) / delta_old;
			d = Add(r, Times(d, beta));
			k++;
			if (k == n || InnerProduct(r, d) <= 0) {
				d = r;
				k = 0;
			}
			ff = dff;
			i++;
		}
		*step = i;
		return x;
	}
}