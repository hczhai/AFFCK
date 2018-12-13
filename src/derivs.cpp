
#include "derivs.h"

namespace affck {

	const int Derivs::DTheta::Num = 3;
	const int Derivs::DL::Num = 2;

	void Derivs::DTheta::rpl_lpx(const vector<double>& x, vector<double>& rpl, vector<vector<double> >& lpx) const {
		double la = pow(x[3] - x[0], 2) + pow(x[4] - x[1], 2) + pow(x[5] - x[2], 2);
		double lb = pow(x[6] - x[3], 2) + pow(x[7] - x[4], 2) + pow(x[8] - x[5], 2);
		double lc = pow(x[6] - x[0], 2) + pow(x[7] - x[1], 2) + pow(x[8] - x[2], 2);

		if (la < 0) la = 0;
		if (lb < 0) lb = 0;
		double tmp;

		rpl[0] = (-la + lb - lc) / (2 * pow(la, 1.5) * sqrt(lb) *
			sqrt((tmp = -((la * la - 2 * la * (lb + lc) + pow(lb - lc, 2)) / (la * lb))) < 0 ? 0 : tmp));
		rpl[1] = (la - lb - lc) / (2 * sqrt(la) * pow(lb, 1.5) *
			sqrt((tmp = -((la * la - 2 * la * (lb + lc) + pow(lb - lc, 2)) / (la * lb))) < 0 ? 0 : tmp));
		rpl[2] = 1 / (2 * sqrt(la) * sqrt(lb) *
			sqrt((tmp = 1 - pow(la + lb - lc, 2) / (4 * la * lb)) < 0 ? 0 : tmp));

		for (int i = 0; i < 3; i++) if (!isfinite(rpl[i]) || isnan(rpl[i])) rpl[i] = 0;

		for (int i = 0; i < 3; i++) lpx[i] = vector<double>(Num * 3);

		lpx[0][0] = -2 * (x[3] - x[0]);
		lpx[0][1] = -2 * (x[4] - x[1]);
		lpx[0][2] = -2 * (x[5] - x[2]);
		lpx[0][3] = 2 * (x[3] - x[0]);
		lpx[0][4] = 2 * (x[4] - x[1]);
		lpx[0][5] = 2 * (x[5] - x[2]);

		lpx[1][3] = -2 * (x[6] - x[3]);
		lpx[1][4] = -2 * (x[7] - x[4]);
		lpx[1][5] = -2 * (x[8] - x[5]);
		lpx[1][6] = 2 * (x[6] - x[3]);
		lpx[1][7] = 2 * (x[7] - x[4]);
		lpx[1][8] = 2 * (x[8] - x[5]);

		lpx[2][0] = -2 * (x[6] - x[0]);
		lpx[2][1] = -2 * (x[7] - x[1]);
		lpx[2][2] = -2 * (x[8] - x[2]);
		lpx[2][6] = 2 * (x[6] - x[0]);
		lpx[2][7] = 2 * (x[7] - x[1]);
		lpx[2][8] = 2 * (x[8] - x[2]);
	}

	double Derivs::DTheta::F(const vector<double>& x) const {
		double la = pow(x[3] - x[0], 2) + pow(x[4] - x[1], 2) + pow(x[5] - x[2], 2);
		double lb = pow(x[6] - x[3], 2) + pow(x[7] - x[4], 2) + pow(x[8] - x[5], 2);
		double lc = pow(x[6] - x[0], 2) + pow(x[7] - x[1], 2) + pow(x[8] - x[2], 2);
		return acos((la + lb - lc) / (2 * sqrt(la) * sqrt(lb)));
	}

	const vector<double> Derivs::DTheta::Fp(const vector<double>& x) const {
		vector<double> r(Num * 3);
		vector<double> rpl(3);
		vector<vector<double> > lpx(3);

		rpl_lpx(x, rpl, lpx);
		for (size_t i = 0; i < r.size(); i++) {
			r[i] = 0;
			for (size_t j = 0; j < lpx.size(); j++)
				r[i] += rpl[j] * lpx[j][i];
		}
		return r;
	}

	const vector<vector<double> > Derivs::DTheta::Fpp(const vector<double>& x) const {
		vector<vector<double> > r(Num * 3);
		vector<double> rpl(3);
		vector<vector<double> > lpx(3);

		rpl_lpx(x, rpl, lpx);
		vector<vector<double> > rppl(3);
		vector<vector<vector<double> > > lppx(3);
		for (int i = 0; i < 3; i++) {
			rppl[i] = vector<double>(3);
			lppx[i] = vector<vector<double> >(Num * 3);
			for (int j = 0; j < 3 * Num; j++)
				lppx[i][j] = vector<double>(Num * 3);
		}

		for (int i = 0; i < 3 * Num; i++) r[i] = vector<double>(Num * 3);

		double la = pow(x[3] - x[0], 2) + pow(x[4] - x[1], 2) + pow(x[5] - x[2], 2);
		double lb = pow(x[6] - x[3], 2) + pow(x[7] - x[4], 2) + pow(x[8] - x[5], 2);
		double lc = pow(x[6] - x[0], 2) + pow(x[7] - x[1], 2) + pow(x[8] - x[2], 2);

		rppl[0][0] = -((lc * (la * la + 3 * lb * lb) - 3 * lc * lc * (la + lb) + pow(la - lb, 3) +
			lc * lc * lc) / (2 * pow(la, 7.0 / 2.0) * pow(lb, 3.0 / 2.0) * pow(-((la *
			la - 2 * la * (lb + lc) + pow(lb - lc, 2)) / (la * lb)), 3.0 / 2.0)));
		rppl[0][1] = rppl[1][0] = (2 * lc) / (pow(la, 3.0 / 2.0) * pow(lb, 3.0 / 2.0) *
			pow(-((la * la - 2 * la * (lb + lc) + pow(lb - lc, 2)) / (la * lb)), 3.0 / 2.0));
		rppl[0][2] = rppl[2][0] = -((-la + lb + lc) / (pow(la, 3.0 / 2.0) * pow(lb, 3.0 / 2.0)
			* pow(-((la * la - 2 * la * (lb + lc) + pow(lb - lc, 2)) / (la * lb)), 3.0 / 2.0)));

		rppl[1][1] = -((lc * (3 * la * la + lb * lb) - 3 * lc * lc * (la + lb) - pow(la - lb, 3) +
			lc * lc * lc) / (2 * pow(la, 3.0 / 2.0) * pow(lb, 7.0 / 2.0) * pow(-((la *
			la - 2 * la * (lb + lc) + pow(lb - lc, 2)) / (la * lb)), 3.0 / 2.0)));
		rppl[1][2] = rppl[2][1] = -((la - lb + lc) / (pow(la, 3.0 / 2.0) * pow(lb, 3.0 / 2.0) *
			pow(-((la * la - 2 * la * (lb + lc) + pow(lb - lc, 2)) / (la * lb)), 3.0 / 2.0)));

		rppl[2][2] = -((la + lb - lc) / (8 * pow(la, 3.0 / 2.0) * pow(lb, 3.0 / 2.0) *
			pow(1 - pow(la + lb - lc, 2) / (4 * la * lb), 3.0 / 2.0)));

		for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		if (!isfinite(rppl[i][j]) || isnan(rppl[i][j])) rppl[i][j] = 0;

		lppx[0][0][0] = 2;
		lppx[0][0][3] = lppx[0][3][0] = -2;
		lppx[0][1][1] = 2;
		lppx[0][1][4] = lppx[0][4][1] = -2;
		lppx[0][2][2] = 2;
		lppx[0][2][5] = lppx[0][5][2] = -2;
		lppx[0][3][3] = 2;
		lppx[0][4][4] = 2;
		lppx[0][5][5] = 2;

		lppx[1][3][3] = 2;
		lppx[1][3][6] = lppx[2][6][3] = -2;
		lppx[1][4][4] = 2;
		lppx[1][4][7] = lppx[2][7][4] = -2;
		lppx[1][5][5] = 2;
		lppx[1][5][8] = lppx[2][8][5] = -2;
		lppx[1][6][6] = 2;
		lppx[1][7][7] = 2;
		lppx[1][8][8] = 2;

		lppx[2][0][0] = 2;
		lppx[2][0][6] = lppx[1][6][0] = -2;
		lppx[2][1][1] = 2;
		lppx[2][1][7] = lppx[1][7][1] = -2;
		lppx[2][2][2] = 2;
		lppx[2][2][8] = lppx[1][8][2] = -2;
		lppx[2][6][6] = 2;
		lppx[2][7][7] = 2;
		lppx[2][8][8] = 2;

		for (int k = 0; k < 3 * Num; k++)
		for (int l = 0; l < 3 * Num; l++) {
			r[k][l] = 0;
			for (int i = 0; i < 3; i++) {
				double sj = 0;
				for (int j = 0; j < 3; j++)
					sj += lpx[j][k] * rppl[i][j];
				r[k][l] += sj * lpx[i][l];
				r[k][l] += rpl[i] * lppx[i][k][l];
			}
		}

		return r;
	}

	double Derivs::DL::F(const vector<double>& x) const {
		return sqrt(pow(x[0] - x[3], 2) + pow(x[1] - x[4], 2) + pow(x[2] - x[5], 2));
	}

	const vector<double> Derivs::DL::Fp(const vector<double>& x) const {
		vector<double> r(Num * 3);
		double rr = sqrt(pow(x[0] - x[3], 2) +
			pow(x[1] - x[4], 2) + pow(x[2] - x[5], 2));
		r[0] = (x[0] - x[3]) / rr;
		r[1] = (x[1] - x[4]) / rr;
		r[2] = (x[2] - x[5]) / rr;
		r[3] = (x[3] - x[0]) / rr;
		r[4] = (x[4] - x[1]) / rr;
		r[5] = (x[5] - x[2]) / rr;
		return r;
	}

	const vector<vector<double> > Derivs::DL::Fpp(const vector<double>& x) const {
		vector<vector<double> > r(Num * 3);
		for (int i = 0; i < 3 * Num; i++) r[i] = vector<double>(Num * 3);

		double rrsq = pow(x[0] - x[3], 2) + pow(x[1] - x[4], 2) + pow(x[2] - x[5], 2);
		double rqs = pow(rrsq, 3.0 / 2.0);

		r[0][0] = (pow(x[1] - x[4], 2) + pow(x[2] - x[5], 2)) / rqs;
		r[0][1] = r[1][0] = -(((x[0] - x[3]) * (x[1] - x[4])) / rqs);
		r[0][2] = r[2][0] = -(((x[0] - x[3]) * (x[2] - x[5])) / rqs);
		r[0][3] = r[3][0] = (-pow(x[1] - x[4], 2) - pow(x[2] - x[5], 2)) / rqs;
		r[0][4] = r[4][0] = ((x[0] - x[3]) * (x[1] - x[4])) / rqs;
		r[0][5] = r[5][0] = ((x[0] - x[3]) * (x[2] - x[5])) / rqs;

		r[1][1] = (pow(x[0] - x[3], 2) + pow(x[2] - x[5], 2)) / rqs;
		r[1][2] = r[2][1] = -(((x[1] - x[4]) * (x[2] - x[5])) / rqs);
		r[1][3] = r[3][1] = ((x[0] - x[3]) * (x[1] - x[4])) / rqs;
		r[1][4] = r[4][1] = (-pow(x[0] - x[3], 2) - pow(x[2] - x[5], 2)) / rqs;
		r[1][5] = r[5][1] = ((x[1] - x[4]) * (x[2] - x[5])) / rqs;

		r[2][2] = (pow(x[0] - x[3], 2) + pow(x[1] - x[4], 2)) / rqs;
		r[2][3] = r[3][2] = ((x[0] - x[3]) * (x[2] - x[5])) / rqs;
		r[2][4] = r[4][2] = ((x[1] - x[4]) * (x[2] - x[5])) / rqs;
		r[2][5] = r[5][2] = (-pow(x[0] - x[3], 2) - pow(x[1] - x[4], 2)) / rqs;

		r[3][3] = (pow(x[1] - x[4], 2) + pow(x[2] - x[5], 2)) / rqs;
		r[3][4] = r[4][3] = -(((x[0] - x[3]) * (x[1] - x[4])) / rqs);
		r[3][5] = r[5][3] = -(((x[0] - x[3]) * (x[2] - x[5])) / rqs);

		r[4][4] = (pow(x[0] - x[3], 2) + pow(x[2] - x[5], 2)) / rqs;
		r[4][5] = r[5][4] = -(((x[1] - x[4]) * (x[2] - x[5])) / rqs);

		r[5][5] = (pow(x[0] - x[3], 2) + pow(x[1] - x[4], 2)) / rqs;

		return r;
	}
}