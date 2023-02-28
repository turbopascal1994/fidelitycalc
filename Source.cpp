#include <algorithm>
#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "linalg.h"

using namespace std;
static const double PI = acos(-1.0);
static const double h = 1.054e-34;
static const double C1 = 1e-12;
static const double F0 = 2.06e-15;
struct FinalState { double Level0 = 0; double Level1 = 0; double Level2 = 0; vector<complex<double>> WF; };
struct Results { double Leakage; double Fidelity; };

class Kernel {
	double tstep, w01, w12, wt, w, T;
	vector<complex<double>> Id, H0, EigVec, EigVal, WF0, WF1, WF2, Hmatrix, InitStates, FinalStates;
		void precalcFidelityAndRotateCheck(double tstep, double w01, double w12) {

		Id = { {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0} };

		H0 = { {0, 0}, {0, 0}, {0, 0}, {0, 0}, {h * w01, 0}, {0, 0}, {0, 0}, {0, 0}, {h * (w01 + w12), 0} };
		EigVec.resize(9);
		EigVal.resize(3);
		eig(H0, EigVec, EigVal, 3);

		vector<int> indices = { 0, 1, 2 };
		sort(indices.begin(), indices.end(), [&](int a, int b) {
			return EigVal[a].real() < EigVal[b].real();
			});
		WF0 = { EigVec[indices[0]], EigVec[indices[0] + 3], EigVec[indices[0] + 6] };
		WF1 = { EigVec[indices[1]], EigVec[indices[1] + 3], EigVec[indices[1] + 6] };
		WF2 = { EigVec[indices[2]], EigVec[indices[2] + 3], EigVec[indices[2] + 6] };

		Hmatrix = { {0, 0}, {-1, 0}, {0, 0}, {1, 0}, {0, 0}, {-sqrt(2), 0}, {0, 0}, {sqrt(2), 0}, {0, 0} };
		InitStates = {
			{1, 0}, {0, 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0},
			{0, 0}, {1, 0}, {1.0 / sqrt(2), 0}, {-1.0 / sqrt(2), 0}, {0, 1.0 / sqrt(2)}, {0, -1.0 / sqrt(2)},
			{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
		};
		// конечные состояния в случае идеальной операции
		FinalStates = {
			{1.0 / sqrt(2), 0}, {-1.0 / sqrt(2), 0}, {0, 0}, {1, 0}, {0.5, -0.5}, {0.5, 0.5},
			{1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1, 0}, {0, 0}, {0.5, 0.5}, {0.5, -0.5},
			{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
		};
	}

	FinalState calcProbability(
		int N,
		const vector<int>& InputSequence,
		const int CyclePlusMinusSteps,
		const int CycleZeroSteps,
		const vector<complex<double>>& UPlus,
		const vector<complex<double>>& UMinus,
		const vector<complex<double>>& UPlus2,
		const vector<complex<double>>& UMinus2,
		const vector<complex<double>>& UZero,
		vector<complex<double>>& WF,
		const vector<complex<double>>& WF0,
		const vector<complex<double>>& WF1,
		const vector<complex<double>>& WF2
	) {
		auto UPlus80 = pow_matrix(UPlus, CyclePlusMinusSteps, N);
		auto UMinus80 = pow_matrix(UMinus, CyclePlusMinusSteps, N);
		auto UPlus280 = pow_matrix(UPlus2, CyclePlusMinusSteps, N);
		auto UMinus280 = pow_matrix(UMinus2, CyclePlusMinusSteps, N);
		auto UZero80 = pow_matrix(UZero, CyclePlusMinusSteps, N);
		auto UZero720 = pow_matrix(UZero, CycleZeroSteps, N);

		auto UPlusZero = mult(UPlus80, UZero720, N, N, N, N, N, N);
		auto UMinusZero = mult(UMinus80, UZero720, N, N, N, N, N, N);
		auto UPlus2Zero = mult(UPlus280, UZero720, N, N, N, N, N, N);
		auto UMinus2Zero = mult(UMinus280, UZero720, N, N, N, N, N, N);
		auto UZeroZero = mult(UZero80, UZero720, N, N, N, N, N, N);

		vector<complex<double>> resU(N * N);
		for (int i = 0; i < N; i++) resU[i * N + i] = { 1, 0 };
		for (int i = 0; i < InputSequence.size(); i++) {
			vector<complex<double>> U = UZeroZero;
			if (InputSequence[i] == -1) {
				U = UMinusZero;
			}
			else if (InputSequence[i] == 1) {
				U = UPlusZero;
			}
			else if (InputSequence[i] == 2) {
				U = UPlus2Zero;
			}
			else if (InputSequence[i] == -2) {
				U = UMinus2Zero;
			}
			WF = mult(U, WF, N, 1, N, N, 1, 1);
		}

		// населенности всех уровней кубита
		complex<double> res0 = {0,0}, res1 = {0,0}, res2 = {0,0};
		for (int i = 0; i < N; i++) {
			res0 += WF[i] * conj(WF0[i]);
			res1 += WF[i] * conj(WF1[i]);
			res2 += WF[i] * conj(WF2[i]);
		}
		FinalState out{ norm(res0), norm(res1), norm(res2), WF };
		return out;
	}

public:
	Kernel(
		double tstep,
		double w01,
		double w12,
		double wt,
		double w,
		double T) : tstep(tstep), w01(w01), w12(w12), wt(wt), w(w), T(T) {
		precalcFidelityAndRotateCheck(tstep, w01, w12);
	}

	Results Output(const vector<int>& SignalString, double Theta) {
		int CellsNumber = SignalString.size();
		int N = 3;
		double V = F0 / w;
		double Cc = Theta / (F0 * sqrt(2.0 * w01 / (h * C1)));
		double Amp = Cc * V * sqrt(h * w01 / (2.0 * C1));

		//	число тактов с импульсом и без него на одном периоде тактовой частоты генератора
		int CycleSteps = floor(T / tstep + 0.5);
		int CyclePlusMinusSteps = floor(CycleSteps * w / T + 0.5);
		int CycleZeroSteps = CycleSteps - CyclePlusMinusSteps;

		FinalState res;
		double leak = 0, F = 0;
		complex<double> unF;
		vector<complex<double>> HrPlus(H0.size());
		vector<complex<double>> HrMinus(H0.size());
		vector<complex<double>> HrPlus2(H0.size());
		vector<complex<double>> HrMinus2(H0.size());
		vector<complex<double>> HrZero(H0.size());
		for (size_t i = 0; i < H0.size(); ++i) {
			HrPlus[i] = { H0[i].real(), Amp * Hmatrix[i].real() };
			HrMinus[i] = { H0[i].real(), -Amp * Hmatrix[i].real() };
			HrPlus2[i] = { H0[i].real(), Amp * 2 * Hmatrix[i].real() };
			HrMinus2[i] = { H0[i].real(), -Amp * 2 * Hmatrix[i].real() };
			HrZero[i] = H0[i];
		}

		auto UPlus = getUMatrix(Id, HrPlus, tstep, h, N);
		auto UMinus = getUMatrix(Id, HrMinus, tstep, h, N);
		auto UPlus2 = getUMatrix(Id, HrPlus2, tstep, h, N);
		auto UMinus2 = getUMatrix(Id, HrMinus2, tstep, h, N);
		auto UZero = getUMatrix(Id, HrZero, tstep, h, N);

		cout << fixed;
		cout << setprecision(6);
		for (size_t IS = 0; IS < 6; ++IS) {
			cout << "N" << IS+1 << "\n";
			vector<complex<double>> IWF = { InitStates[IS], InitStates[IS + 6], InitStates[IS + 12] };
			vector<complex<double>> FWF = { FinalStates[IS], FinalStates[IS + 6], FinalStates[IS + 12] };
			cout << "initial: " << IWF[0] << " " << IWF[1] << " " << IWF[2] << "\n";
			cout << "final: " << FWF[0] << " " << FWF[1] << " " << FWF[2] << "\n";
			res = calcProbability(N, SignalString, CyclePlusMinusSteps, CycleZeroSteps, UPlus, UMinus, UPlus2, UMinus2, UZero, IWF, WF0, WF1, WF2);
			cout << "WF" << IS+1 << ": " << res.WF[0] << " " << res.WF[1] << " " << res.WF[2] << "\n";
			cout << "P0: " << res.Level0 << ", P1: " << res.Level1 << ", P2: " << res.Level2 << "\n";
			leak += res.Level2;
			for (int i = 0; i < N; i++) {
				unF += res.WF[i] * conj(FWF[i]);
			}
			F = norm(unF);
			//cout << "F = " << F << "\n";
		}
		Results out{ leak/6, F/6 };
		return out;
	}
};

// функция перевода последовательности в виде строки в массив
vector<int> convert(string InputString, int NumberOfCycles) {
	vector<int> NewSequence;
	int lenseq = InputString.length(), j = 0, minus = 0;
	char symb;
	for (int i = 0; i < lenseq; i++) {
		symb = InputString.at(i);
		if (symb == '2') {
			if (minus == 1)
			{
				NewSequence.push_back(-2);
				j++;
				minus = 0;
			}
			else {
				NewSequence.push_back(2);
				j++;
			}
		}
		else if (symb == '1') {
			if (minus == 1)
			{
				NewSequence.push_back(-1);
				j++;
				minus = 0;
			}
			else {
				NewSequence.push_back(1);
				j++;
			}
		}
		else if (symb == '0') {
			NewSequence.push_back(0);
			j++;
		}
		else if (symb == '-') minus = 1;
	}
	vector<int> CycledSequence;
	int k = 0;
	for (int i = 0; i < NumberOfCycles* NewSequence.size(); i++) {
		CycledSequence.push_back(NewSequence[k]);
		if (k == NewSequence.size() - 1) k = 0;
		else k++;
	}
	return CycledSequence;
}

int main() {
	const double w01 = 4.73047 * 2 * PI * 1e9;
	const double w12 = w01 - 0.25 * 2 * PI * 1e9;
	const double wt = 25 * 2 * PI * 1e9;
	const double w = 4e-12;
	const double T = PI / wt * 2;
	const double Theta = 0.032;
	const double tstep = 1e-15;
	const int NumberOfCycles = 8;
	Kernel kernel(tstep, w01, w12, wt, w, T);
	string SeqStr = "1110111100011000100001000110001111011";
	vector<int> Sequence = convert(SeqStr, NumberOfCycles);
	Results res = kernel.Output(Sequence, Theta);
	/*cout << "w01 = " << w01 / (2 * PI * 1e9) << "GHz" << "\n";
	cout << "w12 = " << w12/(2 * PI * 1e9) << "GHz" << "\n";
	cout << "wt = " << wt/(2 * PI * 1e9) << "GHz" << "\n";
	cout << "Angle = " << Theta << "\n";
	cout << "Sequence: ";
	for (int i = 0; i < NumberOfCycles; i++)
		cout << SeqStr;*/
	cout << "\n" << "Overall:" << "\n";
	cout << "W2 = " << res.Leakage << "\n";
	cout << "F = " << res.Fidelity << "\n";
}